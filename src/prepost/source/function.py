import numpy as np
import os
import time
import h5py
import scipy.interpolate as interp
from scipy.special import gamma
import multiprocessing
import matplotlib.pyplot as plt
import logging

from .fresnelCS import fresnelCS
from .WP_LEE import WP_LEE
from .WP_Goody import WP_Goody
from  . import BEM as bem
from prepost.les import Les
from prepost.source.wind_turbine import WindTurbine
from prepost.source.mesh import Mesh


from .utils import interp_atmos_data, interp_3D_atmos_data, cotte_conv, c_round, find_index, R1_func,Aweight
from .parallel import domain_facture, output_construct,reshuffle_output

def SPL_SWL_function(index,input_var,return_dict,BEM=False):
    wind_turbine  = input_var['wind_turbine']
    atmos = input_var['atmos']
    coord = input_var['coord']
    freq = input_var['freq']

    # receivers = cotte_conv(coord)
    x_rec = coord['x_coord']
    z_rec = coord['z_coord']
    tau = coord['tau_coord']

    # Wind turbine input
    #---------------------------------------------------------------------------
    H_ref = wind_turbine.href
    Omega = wind_turbine.omega
    Nseg = wind_turbine.Nseg
    Nblade = wind_turbine.Nblade
    seg = wind_turbine.seg.reshape(1, Nseg)
    Lspan = wind_turbine.Lspan.reshape(1, Nseg)
    c = wind_turbine.chord.reshape(1, Nseg)
    twist = wind_turbine.twist.reshape(1, Nseg)
    airfoil_data = wind_turbine.airfoil_data
    delta_beta = wind_turbine.delta_beta
    blade_length = seg[0, -1] + 0.5 * Lspan[0, -1]

    Reynolds = airfoil_data['reynolds'].reshape(-1)
    AoA = airfoil_data['AoA'].reshape(-1)
    delta_star_bot = airfoil_data['delta_star_bot']
    theta_momen_bot = airfoil_data['theta_momen_bot']
    dpdx_bot = airfoil_data['dpdx_bot']
    UeUinf_bot = airfoil_data['UeUinf_bot']
    cf_bot = airfoil_data['cf_bot']

    delta_star_top = airfoil_data['delta_star_top']
    theta_momen_top = airfoil_data['theta_momen_top']
    dpdx_top = airfoil_data['dpdx_top']
    UeUinf_top = airfoil_data['UeUinf_top']
    cf_top = airfoil_data['cf_top']
    if BEM == True :
        Cl_tck = airfoil_data['Cl_tck']
        Cd_tck = airfoil_data['Cd_tck']
        # quit()

    #Corresponds to the blade count
    dbeta = 2 * np.pi / Nblade
    # beta in clockwise directions looking in the wind direction
    # rearrange nbeta and nseg
    beta = np.linspace(0, 2 * np.pi - dbeta, Nblade).reshape(Nblade, 1) + delta_beta
    seg = seg.reshape(1,Nseg)

    # Atmospheric inputs
    #---------------------------------------------------------------------------
    # old version
    # epsilon_Kol = atmos.epsilon_Kol
    # U_inf = atmos.U_inf
    # z_atmos = atmos.z_coord

    # coordinate of the segment
    x_ref = (np.sin(beta) * seg).reshape(len(beta), seg.shape[1])                  #x_ref
    h_temp = (np.cos(beta) * seg + H_ref).reshape(len(beta), seg.shape[1])         #hight of current segment

    # interpolate atmospheric field on segment position (cubic interpolation)
    # interpolate from 1D profile
    # U_inf = interp_atmos_data(atmos.z_coord,atmos.U_inf,h_temp)
    # epsil = interp_atmos_data(z_atmos,epsilon_Kol,h_temp)
    # epsil = epsilon_Kol + np.zeros((beta.shape[0], seg.shape[1]))

    xnew = np.zeros_like(h_temp)
    xnew[...] = wind_turbine.absolute_pos[0]
    # constant value 
    #ynew = x_ref*0 + wind_turbine.absolute_pos[1]
    # znew = h_temp*0 + H_ref

    # interpolate 1D profile
    # ynew = x_ref*0 + wind_turbine.absolute_pos[1]
    # znew = h_temp

    # interpolate from 3D fields
    ynew = x_ref + wind_turbine.absolute_pos[1]
    znew = h_temp

    U_inf, epsil = interp_3D_atmos_data(atmos, xnew, ynew, znew)

    # generic atmospheric data 
    Tref = 10 + 273.15                                                        # reference temperature at height zref (10°C)
    zref = 2
    gamma0 = 1.4
    r = 287
    c0 = np.sqrt(gamma0 * r * Tref)
    rho0 = 1.225
    nu0 = 1.45e-5                                                             #kinetic viscosity of air
    Cp0 = 1004

    freq = freq.reshape(1,len(freq))
    Nfreq = freq.shape[1]
    lambda0 = c0 / freq
    k0 = 2 * np.pi / lambda0
    M = U_inf / c0

    # parameters for Amiet's model
    b = c/2                                                                        # half chord, c=2b
    d = Lspan/2;
    e = c * 0.15                                                                   # thickness of each segment (for NACA63415 thickness is 15%)
    U_rot = seg * Omega
    M_rot = U_rot / c0

    omega = 2 * np.pi * freq                                                       # frequency of incident pressure, also the frequency of induced acoustic waves
    k0 = omega / c0                                                                # acoustic wave number
    k_bar = k0.T * b
    kappa = 0.41

    #The Kolmogorov term is extracted as per the theory  from  Buck et al 2018
    A = 1.65;
    const1 =  (np.pi ** (1/2)) * (gamma(5/6) / gamma(1/3))
    Kol_term= (A * 9/55 * np.pi / const1) * epsil ** (2/3)

    AoA_seg = np.zeros((beta.shape[0], seg.shape[1]))
    rey_seg = np.zeros((beta.shape[0], seg.shape[1]))
    a_array = np.zeros((beta.shape[0], seg.shape[1]))
    adash_array = np.zeros((beta.shape[0], seg.shape[1]))
    epsilon_array = np.zeros((beta.shape[0], seg.shape[1]))
    U_rel = np.zeros((beta.shape[0], seg.shape[1]))
    Uc = np.zeros((beta.shape[0], seg.shape[1]))
    M_rel = np.zeros((beta.shape[0], seg.shape[1]))
    beta_sq = np.zeros((beta.shape[0],seg.shape[1]))
    Kx = np.zeros((freq.shape[1], beta.shape[0], seg.shape[1]))
    Kx_bar = np.zeros((freq.shape[1], beta.shape[0], seg.shape[1]))
    Kc = np.zeros((freq.shape[1], beta.shape[0], seg.shape[1]))
    Kc_bar = np.zeros((freq.shape[1], beta.shape[0], seg.shape[1]))
    mu_bar_TEN = np.zeros((freq.shape[1], beta.shape[0], seg.shape[1]))
    mu_bar_TIN = np.zeros((freq.shape[1], beta.shape[0], seg.shape[1]))
    CoefB = 1.47                                #for the calculation of the correlation length
    CoefC = 0.7                                 #for the calculation of the convection velocity

    if BEM is True :
        print(BEM)
        for ibeta in range(beta.shape[0]):
            for iseg in range(seg.shape[1]):   
                U_rel[ibeta, iseg] = np.sqrt(U_inf[ibeta, iseg] ** 2 + U_rot[0, iseg] ** 2)     
                Re = U_rel[ibeta, iseg] * c[0, iseg] / nu0
                J= Omega * blade_length / U_inf[ibeta, iseg]

                theta = 3*c[0,iseg]/(2*np.pi*seg[0,iseg])            
                # F = 2 / np.pi * np.arccos(np.exp(-(3 / 2) * (1 - seg[0, iseg] / blade_length) * (1 + J ** 2) ** 0.5))

                # AoA_seg[ibeta,iseg],a,adash,F,epsilon = bem.simple(Cl_tck[iseg],Cd_tck[iseg], Re,J,theta,twist[0,iseg],blade_length,seg[0,iseg],F)
                AoA_seg[ibeta,iseg],a,adash,F,epsilon,a_conv,adash_conv = bem.hemant1(Cl_tck[iseg],Cd_tck[iseg],Re,J,blade_length,seg[0,iseg],twist[0,iseg],c[0,iseg])
                # AoA_seg[ibeta,iseg],a,adash,F,epsilon = bem.hemant2(Cl_tck[iseg],Cd_tck[iseg],Re,J,blade_length,seg[0,iseg],twist[0,iseg],c[0,iseg])
                # AoA_seg[ibeta,iseg],a,adash,F,epsilon = bem.noBEM(Cl_tck[iseg],Cd_tck[iseg],Re,J,blade_length,seg[0,iseg],twist[0,iseg],c[0,iseg])

                a_array[ibeta,iseg] = a
                adash_array[ibeta,iseg] = adash
                epsilon_array[ibeta,iseg] = epsilon
                # U_rel[ibeta,iseg] = ((U_inf[ibeta, iseg]*(1-a))**2 +(U_rot[0,iseg]*(1+adash))**2)**0.5
                U_rel[ibeta, iseg] = np.sqrt((U_inf[ibeta, iseg] * (1 - F * a)) ** 2 + (U_rot[0, iseg] * (1 + adash)) ** 2)
                rey_seg[ibeta, iseg] = Re

                Uc[ibeta, iseg] = U_rel[ibeta, iseg] * CoefC
                M_rel[ibeta, iseg] = U_rel[ibeta, iseg] / c0
                beta_sq[ibeta, iseg] = 1 - M_rel[ibeta, iseg] ** 2
                Kx[:, ibeta, iseg] = omega / U_rel[ibeta, iseg]
                Kx_bar[:, ibeta, iseg] = Kx[:, ibeta, iseg] * b[0, iseg]
                Kc[:, ibeta, iseg] = omega / Uc[ibeta, iseg]
                Kc_bar[:, ibeta, iseg] = Kc[:, ibeta, iseg] * b[0, iseg]
                mu_bar_TEN[:, ibeta, iseg] = Kc_bar[:, ibeta, iseg] * M_rel[ibeta, iseg] / beta_sq[ibeta, iseg]
                mu_bar_TIN[:, ibeta, iseg] = Kx_bar[:, ibeta, iseg] * M_rel[ibeta, iseg] / beta_sq[ibeta,iseg]
        # print('AoA')
        # print(AoA_seg)
        # print('epsilon')
        # print(epsilon_array)


    else: 
        for ibeta in range(beta.shape[0]):
            for iseg in range(seg.shape[1]):
                U_rel[ibeta, iseg] = np.sqrt(U_inf[ibeta, iseg] ** 2 + U_rot[0, iseg] ** 2)
                Uc[ibeta, iseg] = U_rel[ibeta, iseg] * CoefC
                M_rel[ibeta, iseg] = U_rel[ibeta, iseg] / c0
                beta_sq[ibeta, iseg] = 1 - M_rel[ibeta, iseg] ** 2
                Kx[:, ibeta, iseg] = omega / U_rel[ibeta, iseg]
                Kx_bar[:, ibeta, iseg] = Kx[:, ibeta, iseg] * b[0, iseg]
                Kc[:, ibeta, iseg] = omega / Uc[ibeta, iseg]
                Kc_bar[:, ibeta, iseg] = Kc[:, ibeta, iseg] * b[0, iseg]
                mu_bar_TEN[:, ibeta, iseg] = Kc_bar[:, ibeta, iseg] * M_rel[ibeta, iseg] / beta_sq[ibeta, iseg]
                mu_bar_TIN[:, ibeta, iseg] = Kx_bar[:, ibeta, iseg] * M_rel[ibeta, iseg] / beta_sq[ibeta,iseg]
                # different convention for AoA
                # AoA_seg[ibeta, iseg] = np.arctan(U_inf[ibeta, iseg] / U_rot[0, iseg]) - twist[0, iseg]
                rey_seg[ibeta, iseg] = U_rel[ibeta, iseg] * c[0, iseg] / nu0
                AoA_seg[ibeta, iseg] = twist[0, iseg] - np.arctan(U_rot[0, iseg] / U_inf[ibeta, iseg])


    AoA_seg1 = AoA_seg * 180 / np.pi

    maxbetac_t = 0                              #corresponding to maximum Clauser's parameter calculated
    Dopp = np.zeros((x_rec.shape[0], seg.shape[1], beta.shape[0]))
    omegae = np.zeros((x_rec.shape[0], freq.shape[0], freq.shape[1]))
    Ky = np.zeros((x_rec.shape[0], freq.shape[0], freq.shape[1]))
    Ky_bar = np.zeros((x_rec.shape[0], freq.shape[0], freq.shape[1]))
    Kx_hat = np.zeros((x_rec.shape[0], freq.shape[0], freq.shape[1]))
    Ky_hat = np.zeros((x_rec.shape[0], freq.shape[0], freq.shape[1]))
    kappa_bar = np.zeros((x_rec.shape[0], freq.shape[0], freq.shape[1]))

    B = np.zeros((x_rec.shape[0], freq.shape[1], freq.shape[0]))
    C = np.zeros((x_rec.shape[0], freq.shape[1], freq.shape[0]))
    C1 = np.zeros((x_rec.shape[0], freq.shape[1], freq.shape[0]))
    S1 = np.zeros((x_rec.shape[0], freq.shape[1], freq.shape[0]))
    C2 = np.zeros((x_rec.shape[0], freq.shape[1], freq.shape[0]))
    S2 = np.zeros((x_rec.shape[0], freq.shape[1], freq.shape[0]))

    LL_TE = np.zeros((x_rec.shape[0], freq.shape[1], freq.shape[0]))
    ly = np.zeros((x_rec.shape[0], freq.shape[0], freq.shape[1]))
    Phi_pp_pres1 = np.zeros((x_rec.shape[0], freq.shape[0], freq.shape[1]))
    Phi_pp_pres2 = np.zeros((x_rec.shape[0], freq.shape[0], freq.shape[1]))
    Phi_pp_suct1 = np.zeros((x_rec.shape[0], freq.shape[0], freq.shape[1]))
    Phi_pp_suct2 = np.zeros((x_rec.shape[0], freq.shape[0], freq.shape[1]))
    Spp_suct = np.zeros((x_rec.shape[0], seg.shape[1], beta.shape[0], freq.shape[1]))
    Spp_pres = np.zeros((x_rec.shape[0], seg.shape[1], beta.shape[0], freq.shape[1]))
    Spp_TEN = np.zeros((x_rec.shape[0], seg.shape[1], beta.shape[0], freq.shape[1]))

    theta1 = np.zeros((x_rec.shape[0], freq.shape[1], freq.shape[0]))
    theta2 = np.zeros((x_rec.shape[0], freq.shape[1], freq.shape[0]))
    C4 = np.zeros((x_rec.shape[0], freq.shape[1], freq.shape[0]))
    S4 = np.zeros((x_rec.shape[0], freq.shape[1], freq.shape[0]))
    LL_TI = np.zeros((x_rec.shape[0], freq.shape[1], freq.shape[0]))
    Phi_ww = np.zeros((x_rec.shape[0], freq.shape[0], freq.shape[1]))
    Spp_TIN = np.zeros((x_rec.shape[0], seg.shape[1], beta.shape[0], freq.shape[1]))
    Spp_tot = np.zeros((x_rec.shape[0], seg.shape[1], beta.shape[0], freq.shape[1]))

    SPL_tot = np.zeros((x_rec.shape[0], Nseg, Nblade, Nfreq))
    SPL_TEN = np.zeros((x_rec.shape[0], Nseg, Nblade, Nfreq))
    SPL_TIN = np.zeros((x_rec.shape[0], Nseg, Nblade, Nfreq))
    SWL_tot = np.zeros((x_rec.shape[0], Nseg, Nblade, Nfreq))

    beta0 = []

    R0 = np.sqrt(x_rec ** 2 + (H_ref-z_rec) ** 2)

    # calculate angles teta and phi according to Schlinker and Amiet (1981) and Rozenberg (2007)

    # receiver coordinates in the (X,Y,Z) coordinate system :
    # teta = np.arccos(x_rec * np.cos(tau) / R0)
    teta = np.arccos(np.divide(x_rec * np.cos(tau) ,R0,out=np.zeros_like(R0),where=R0!=0))
    ZO = R0 * np.cos(teta)
    XO = R0 * np.sin(teta)
    YO = np.zeros(XO.shape[0])                                      # by definition

    # loop over the azimutal blade position
    for ibeta in range(beta.shape[0]):
        beta0.append(dbeta * ibeta)

        # calculation of angle psi with beta =  0 corresponding to blade pointing up (y-axis)
        numer = -(x_rec * np.sin(beta[ibeta, 0]) * np.sin(tau) + (H_ref - z_rec) * np.cos(beta[ibeta, 0]))
        denom = np.sqrt(x_rec ** 2 * np.sin(tau) ** 2 + (H_ref - z_rec) ** 2)
        cospsi = np.divide(numer,denom,out=np.zeros_like(numer),where=denom!=0)

        numer = (x_rec * np.cos(beta[ibeta, 0]) * np.sin(tau) - (H_ref - z_rec) * np.sin(beta[ibeta, 0]))
        denom = np.sqrt(x_rec ** 2 * np.sin(tau) ** 2 + (H_ref - z_rec) ** 2)
        sinpsi = np.divide(numer,denom,out=np.zeros_like(numer),where=denom!=0)

        # cospsi = -(x_rec * np.sin(beta[ibeta, 0]) * np.sin(tau) + (H_ref - z_rec) * np.cos(beta[ibeta, 0])) / np.sqrt(x_rec ** 2 * np.sin(tau) ** 2 + (H_ref - z_rec) ** 2)
        # sinpsi = (x_rec * np.cos(beta[ibeta, 0]) * np.sin(tau) - (H_ref - z_rec) * np.sin(beta[ibeta, 0])) / np.sqrt(x_rec ** 2 * np.sin(tau) ** 2 + (H_ref - z_rec) ** 2)

        # loop over the segments
        for iseg in range(seg.shape[1]):
            if np.isnan(AoA_seg[ibeta,iseg]):
                continue
            # receiver coordinates in the hub coordinate system (C,x2,y2,z2)
            #--------------------------------------------------------------------------------------------
            # rotation of pi/2-psi to align y-axis with span
            MZ = np.zeros((x_rec.shape[0], 3, 3))
            MZ[:, 0, 0] = sinpsi.reshape(x_rec.shape[0], )
            MZ[:, 0, 1] = -cospsi.reshape(x_rec.shape[0], )
            MZ[:, 1, 0] = cospsi.reshape(x_rec.shape[0], )
            MZ[:, 1, 1] = sinpsi.reshape(x_rec.shape[0], )
            MZ[:, 2, 2] = np.ones(x_rec.shape[0])
            # MZ = np.asarray([sinpsi, -cospsi, 0, cospsi, sinpsi, 0, 0, 0, 1]).reshape(3, 3)
            xyz_global = np.asarray([XO.reshape(x_rec.shape[0], ), YO.reshape(x_rec.shape[0], ), ZO.reshape(x_rec.shape[0], )]).T.reshape(x_rec.shape[0], 3, 1)
            xyzO = np.matmul(MZ, xyz_global)

            # receiver coordinates in the segment coordinate system (S,x3,y3,z3),
            MU = np.asarray([np.cos(np.pi/2 - twist[0, iseg]), 0, np.sin(np.pi/2 - twist[0, iseg]), 0, 1, 0, -np.sin(np.pi / 2 - twist[0, iseg]), 0, np.cos(np.pi / 2 - twist[0, iseg])]).reshape(3, 3)
            MU = np.repeat(MU[np.newaxis, :, :], x_rec.shape[0], axis = 0)
            xyz = np.matmul(MU, xyzO)
            M_trans = np.asarray([0, -seg[0, iseg], 0]).reshape(3, 1)
            xyz = xyz + np.repeat(M_trans[np.newaxis, :, :], x_rec.shape[0], axis = 0)

            S0 = np.sqrt(xyz[:, 0, 0] ** 2 + (1 - M_rel[ibeta, iseg] ** 2) * ((xyz[:, 1, 0]) ** 2 + xyz[:, 2, 0] ** 2))

            #calculate Doppler factor and deduce emission frequency omegae
            #--------------------------------------------------------------------------------------------
            Dopp[:, iseg, ibeta] = (1 + M_rot[0, iseg] * sinpsi * np.sin(teta) / np.sqrt(1 - M[ibeta, iseg] ** 2 * np.sin(teta) ** 2)).reshape(x_rec.shape[0], ) # omegae/omega - equation (4.4) of Rozenberg (2007)
            omegae = np.matmul(Dopp[:, iseg, ibeta].reshape(x_rec.shape[0], 1), omega)

            # variables useful for TEN and TIN model at emitted frequency omegae
            #--------------------------------------------------------------------------------------------
            Ky = np.matmul((Dopp[:, iseg, ibeta] * xyz[:, 1, 0] / S0).reshape(x_rec.shape[0], 1), k0)
            Ky_bar = Ky * b[0, iseg]

            #TEN calculation at emitted frequency omegae
            #--------------------------------------------------------------------------------------------
            kappa_bar = np.sqrt((np.matmul(Dopp[:, iseg, ibeta].reshape(x_rec.shape[0], 1), mu_bar_TEN[:, ibeta, iseg].reshape(1, Nfreq))) ** 2 - Ky_bar ** 2 / beta_sq[ibeta,iseg])


            B = np.matmul(Dopp[:, iseg, ibeta].reshape(x_rec.shape[0], 1), Kc_bar[:, ibeta, iseg].reshape(1, Nfreq)) + np.matmul(Dopp[:, iseg, ibeta].reshape(x_rec.shape[0], 1), (M_rel[ibeta, iseg] * mu_bar_TEN[:, ibeta, iseg]).reshape(1, Nfreq)) + kappa_bar
            C = Kc_bar[:, ibeta, iseg].reshape(1, Nfreq) * Dopp[:, iseg, ibeta].reshape(x_rec.shape[0], 1) - np.matmul(Dopp[:, iseg, ibeta].reshape(x_rec.shape[0], 1), mu_bar_TEN[:, ibeta, iseg].reshape(1, Nfreq)) * (xyz[:, 0, 0] / S0 - M_rel[ibeta, iseg]).reshape(x_rec.shape[0], 1)
            C1, S1 = fresnelCS(np.sqrt(2 / np.pi * 2 * np.abs(B-C)))
            C2, S2 = fresnelCS(np.sqrt(2 / np.pi * 2 * B))

            # Aeroacoustic transfer function
            LL_TE =  np.abs(-np.exp( 2 * 1j * C) / (1j * C) * ((1 + 1j) * np.exp(-2 * 1j * C) * np.sqrt(B / (B - C)) * (C1 - 1j * S1) - (1 + 1j) * (C2 - 1j * S2) + 1))

            # span-wise correlation length calculated from the Corcos model
            ly = (omegae / (CoefB * Uc[ibeta, iseg])) / (Ky ** 2 + (omegae / (CoefB * Uc[ibeta,iseg])) ** 2)

            # linear interpolation for AoA and Reynolds of the data base 
            # Warning no interpolation for the segments hence must be exactly the same 
            ind_AoA = find_index(AoA_seg[ibeta,iseg] * 180 / np.pi, AoA)
            if AoA_seg1[ibeta,iseg] > max(AoA) : 
                logging.warning('AoA outside of bound set to %s instead of %s' % (AoA[-1], AoA_seg1[ibeta,iseg]))
                ind_AoA_left = len(AoA)-1
                ind_AoA_right = 0
                ratio_AoA = 0
            elif AoA_seg1[ibeta,iseg] < min(AoA) : 
                ind_AoA_left = 0
                ind_AoA_right = 0
                ratio_AoA = 0
                logging.warning('AoA outside of bound set to %s instead of %s' % (AoA[0], AoA_seg1[ibeta,iseg]))
            else : 
                ind_AoA_left = ind_AoA - 1
                ind_AoA_right = ind_AoA
                ratio_AoA = (AoA_seg[ibeta,iseg] * 180 / np.pi - AoA[ind_AoA_left]) / (AoA[ind_AoA_right] - AoA[ind_AoA_left])

            delta_star_top1 = delta_star_top[:, ind_AoA_left, iseg] + ratio_AoA * (delta_star_top[:, ind_AoA_right, iseg] - delta_star_top[:, ind_AoA_left, iseg])
            delta_star_bot1 = delta_star_bot[:, ind_AoA_left, iseg] + ratio_AoA * (delta_star_bot[:, ind_AoA_right, iseg] - delta_star_bot[:, ind_AoA_left, iseg])


            theta_momen_top1 = theta_momen_top[:, ind_AoA_left, iseg] + ratio_AoA * (theta_momen_top[:, ind_AoA_right, iseg] - theta_momen_top[:, ind_AoA_left, iseg])
            theta_momen_bot1 = theta_momen_bot[:, ind_AoA_left, iseg] + ratio_AoA * (theta_momen_bot[:, ind_AoA_right, iseg] - theta_momen_bot[:, ind_AoA_left, iseg])

            dpdx_top1 = dpdx_top[:, ind_AoA_left, iseg] + ratio_AoA * (dpdx_top[:, ind_AoA_right, iseg] - dpdx_top[:, ind_AoA_left, iseg])
            dpdx_bot1 = dpdx_bot[:, ind_AoA_left, iseg] + ratio_AoA * (dpdx_bot[:, ind_AoA_right, iseg] - dpdx_bot[:, ind_AoA_left, iseg])

            cf_top1 = cf_top[:, ind_AoA_left, iseg] + ratio_AoA * (cf_top[:, ind_AoA_right, iseg] - cf_top[:, ind_AoA_left, iseg])
            cf_bot1 = cf_bot[:, ind_AoA_left, iseg] + ratio_AoA * (cf_bot[:, ind_AoA_right, iseg] - cf_bot[:, ind_AoA_left, iseg])


            UeUinf_top1 = UeUinf_top[:, ind_AoA_left, iseg] + ratio_AoA * (UeUinf_top[:, ind_AoA_right, iseg] - UeUinf_top[:, ind_AoA_left, iseg])
            UeUinf_bot1 = np.abs(UeUinf_bot[:, ind_AoA_left, iseg] + ratio_AoA * (UeUinf_bot[:, ind_AoA_right, iseg] - UeUinf_bot[:, ind_AoA_left, iseg]))

            Reyn = U_rel[ibeta,iseg] * c[0, iseg] / nu0
            ind_re = find_index(Reyn, Reynolds)
            ind_left = ind_re - 1
            ind_right = ind_re
            ratio = (Reyn-Reynolds[ind_left]) / (Reynolds[ind_right] - Reynolds[ind_left])

            delta_star_t = delta_star_top1[ind_left] + ratio * (delta_star_top1[ind_right] - delta_star_top1[ind_left])
            delta_star_b = delta_star_bot1[ind_left] + ratio * (delta_star_bot1[ind_right] - delta_star_bot1[ind_left])

            theta_momen_t = theta_momen_top1[ind_left] + ratio * (theta_momen_top1[ind_right] - theta_momen_top1[ind_left])
            theta_momen_b = theta_momen_bot1[ind_left] + ratio * (theta_momen_bot1[ind_right] - theta_momen_bot1[ind_left])

            dpdx_t = dpdx_top1[ind_left] + ratio * (dpdx_top1[ind_right] - dpdx_top1[ind_left])
            dpdx_b = dpdx_bot1[ind_left] + ratio * (dpdx_bot1[ind_right] - dpdx_bot1[ind_left])

            cf_t = cf_top1[ind_left] + ratio * (cf_top1[ind_right] - cf_top1[ind_left])
            cf_b = cf_bot1[ind_left] + ratio * (cf_bot1[ind_right] - cf_bot1[ind_left])

            UeUinf_t = UeUinf_top1[ind_left] + ratio * (UeUinf_top1[ind_right] - UeUinf_top1[ind_left])
            UeUinf_b = np.abs(UeUinf_bot1[ind_left] + ratio * (UeUinf_bot1[ind_right] - UeUinf_bot1[ind_left]))

            #Calculation of Clauser parameter in the pressure and suction sides
            # to determine which WPS model to use
            Ue_t = U_rel[ibeta,iseg] * np.abs(UeUinf_t)                    # external velocity (m)
            tauw_t = cf_t * 0.5 * rho0 * Ue_t ** 2                         # wall shear stress (Pa)
            betac_t = theta_momen_t / tauw_t * dpdx_t                      # Clauser parameter
            maxbetac_t = max(maxbetac_t, betac_t)

            # Wall pressure spectrum calculation 
            # Suction side
            if cf_t < 0:  
                logging.warning('top : BL separated for seg %s and angle %s, cf= %s' %  (iseg,beta[ibeta, 0], cf_t ))                                                     # no model used because BL is separated
                Phi_pp_suct1 = 0
                Phi_pp_suct2 = 0
            elif betac_t > 60:       
                logging.warning('top : AP very strong BL may be separated') 
                dpdx_t = 60 / theta_momen_t * tauw_t
                Phi_pp_suct1, Phi_pp_suct2 = WP_LEE(np.matmul(Dopp[:, iseg, ibeta].reshape(x_rec.shape[0], 1), freq.reshape(1, Nfreq)), delta_star_t, theta_momen_t, dpdx_t, cf_t, UeUinf_t, U_rel[ibeta,iseg])
            elif betac_t > 0 and betac_t < 60:                                              # Lee's model
                Phi_pp_suct1, Phi_pp_suct2 = WP_LEE(np.matmul(Dopp[:, iseg, ibeta].reshape(x_rec.shape[0], 1), freq.reshape(1, Nfreq)), delta_star_t, theta_momen_t, dpdx_t, cf_t, UeUinf_t ,U_rel[ibeta, iseg])
            elif betac_t < 0:     
                logging.warning('top : BL separated -- > should never happen')                                         #should never happen on the suction side
                Phi_pp_suct1 = 0
                Phi_pp_suct2 = 0

            # Pressure side
            Ue_b = U_rel[ibeta,iseg] * np.abs(UeUinf_b)                    # external velocity (m)
            tauw_b = cf_b * 0.5 * rho0 * Ue_b ** 2                         # wall shear stress (Pa)
            betac_b = theta_momen_b / tauw_b * dpdx_b                      # Clauser parameter
            if dpdx_b < 0:                                                  # use Goody's model if favorable pressure gradient
                Phi_pp_pres1, Phi_pp_pres2 = WP_Goody(np.matmul(Dopp[:, iseg, ibeta].reshape(x_rec.shape[0], 1), freq.reshape(1, Nfreq)), delta_star_b , theta_momen_b , cf_b , UeUinf_b , U_rel[ibeta,iseg])
            else:
                Phi_pp_pres1, Phi_pp_pres2 = WP_LEE(np.matmul(Dopp[:, iseg, ibeta].reshape(x_rec.shape[0], 1), freq.reshape(1, Nfreq)), delta_star_b, theta_momen_b, dpdx_b, cf_b, UeUinf_b, U_rel[ibeta,iseg])


            # PSD of acoustic pressure on the pressure and suction sides
            Spp_suct[:, iseg, ibeta, :] = (omegae * (xyz[:, 2, 0] * b[0, iseg] / (2 * np.pi * c0 * S0 ** 2)).reshape(x_rec.shape[0], 1)) ** 2 * Lspan[0, iseg] * 2 * LL_TE ** 2 * Phi_pp_suct1 * ly
            Spp_pres[:, iseg, ibeta, :] = (omegae * (xyz[:, 2, 0] * b[0, iseg] / (2 * np.pi * c0 * S0 ** 2)).reshape(x_rec.shape[0], 1)) ** 2 * Lspan[0, iseg] * 2 * LL_TE ** 2 * Phi_pp_pres1 * ly
            # PSD of acoustic pressure at far-field
            Spp_TEN[:, iseg, ibeta, :] = Spp_suct[:, iseg, ibeta, :] + Spp_pres[:, iseg, ibeta, :]


            # Turbulence Inflow Noise calculation at emitted frequency omegae
            #--------------------------------------------------------------------------------------------
            kappa_bar = np.sqrt(np.matmul(Dopp[:, iseg, ibeta].reshape(x_rec.shape[0], 1), mu_bar_TIN[:, ibeta, iseg].reshape(1, Nfreq)) ** 2 - Ky_bar ** 2 / beta_sq[ibeta, iseg])


            theta1 = kappa_bar - np.matmul(Dopp[:, iseg, ibeta].reshape(x_rec.shape[0], 1), mu_bar_TIN[:, ibeta, iseg].reshape(1, Nfreq)) * (xyz[:, 0, 0] / S0).reshape(x_rec.shape[0], 1)
            theta2 = np.matmul(Dopp[:, iseg, ibeta].reshape(x_rec.shape[0], 1), mu_bar_TIN[:, ibeta, iseg].reshape(1, Nfreq)) * (M_rel[ibeta, iseg] - (xyz[0, 0] / S0).reshape(x_rec.shape[0], 1)) - np.pi / 4
            C4, S4 = fresnelCS(np.sqrt(2 / np.pi * 2 * theta1))
            LL_TI =  np.abs(np.sqrt(2 / (np.matmul(Dopp[:, iseg, ibeta].reshape(x_rec.shape[0], 1), Kx_bar[:, ibeta, iseg].reshape(1, Nfreq)) + kappa_bar * beta_sq[ibeta, iseg]) / theta1) * (C4 - 1j * S4) * np.exp(1j * theta2) / np.pi)
            
            kx_kol = np.matmul(Dopp[:, iseg, ibeta].reshape(x_rec.shape[0], 1), Kx[:, ibeta, iseg].reshape(1, Nfreq))              # normalized chordwise turbulent wavenumber
            ky_kol =  Ky                                                   # normalized spanwise turbulent wavenumber

            Phi_ww =  (4 / 9 / np.pi) * ((kx_kol ** 2 + ky_kol ** 2) ** (-4/3)) * Kol_term[ibeta, iseg]
            Spp_TIN[:, iseg, ibeta, :] = 2 * (rho0 * np.matmul((Dopp[:, iseg, ibeta] * xyz[:, 2, 0] / S0 ** 2).reshape(x_rec.shape[0], 1), k0) * c[0, iseg] / 2) ** 2 * np.pi * U_rel[ibeta, iseg] * d[0, iseg] * Phi_ww * LL_TI ** 2
            

            # total contribution
            #------------------------------------------------------------------------------------------
            Spp_tot[:, iseg, ibeta, :] = Spp_TEN[:, iseg, ibeta, :] + Spp_TIN[:, iseg, ibeta, :]

    #for the calculation of SPL by individual blade segments at the given reciever location
    #for the net calculation of the SPL at the receiver location due to all the different blade segments
    for ifreq in range(Nfreq):
        # SPL tot for all segment, beta position and receivers position 
        SPL_tot[:, :, :, ifreq] = 10 * np.log10(Spp_tot[:, :, :, ifreq] * Dopp ** 2 * 2 * np.pi / 4e-10)
        SPL_TIN[:, :, :, ifreq] = 10 * np.log10(Spp_TIN[:, :, :, ifreq] * Dopp ** 2 * 2 * np.pi / 4e-10)
        SPL_TEN[:, :, :, ifreq] = 10 * np.log10(Spp_TEN[:, :, :, ifreq] * Dopp ** 2 * 2 * np.pi / 4e-10)

    R1_array = R1_func(coord, H_ref)
    SWL_tot = SPL_tot + 10 * np.log10(4 * np.pi * R1_array.reshape(-1,1,1,1) ** 2)


    output = {}
    output['AoA'] = AoA_seg1
    output['U_inf'] = U_inf
    output['epsilon'] = epsil
    output['U_rot'] = np.repeat(U_rot.reshape(1,-1),beta.shape[0],axis=1)
    output['U_rel'] = U_rel
    output['a'] = a_array
    output['adash'] = adash_array
    # output['epsilon'] = epsilon_array
    output['SPL_tot'] = SPL_tot
    output['SPL_TIN'] = SPL_TIN
    output['SPL_TEN'] = SPL_TEN
    output['SWL_tot'] = SWL_tot
    # output['OASPL'] = OASPL
    # output['OASPL_beta'] = OASPL_beta
    # output['AM_rec'] = AM_rec
    # output['SWL_freq'] = SWL_freq
    # output['SPL_freq'] = SPL_freq
    output['recon_index'] = input_var['recon_index']
    # output['SPL_tot'] = SPL_tot
    # output['Spp_tot'] = Spp_tot
    # output['Spp_TIN'] = Spp_TIN
    # output['Spp_TEN'] = Spp_TEN
    # output['beta'] = beta

    return_dict[index] = output

    # return OASPL, OASPL_beta, AM_rec, SPL_freq, SWL_freq, SPL_tot, beta

def SPL_SWL_parallel(wt:WindTurbine,atmos:Les,mesh:Mesh,freq:list,Ncore:int,BEM:bool=False, fc:np.ndarray = None, Nfc:np.ndarray = None):
    """Compute SPL in free field and SWL for a given wind turbine, atmospheric condition, mesh and set of frequency.
        The calculations are done in parrallel for different part of the mesh.

    Args:
        wt (WindTurbine): defines the wind turbine geometric property
        atmos (Les): Set the amospheric condition around the turbine
        mesh (Mesh): Defines the mesh on which to compute the resuls
        freq (list): list of frequency for which the calculation are done
        Ncore (int): Number of threads used for the parallel (minimum 2)
        BEM (bool, optional): If True uses BEM theory and drag and lift coefficient of each segment to compute the effective AoA. Defaults to False.
        fc (np.ndarray, optional): Array of third octave band. Must corresponds to freq list and is used to compute SPL and SWL in third octave band. Defaults to None.
        Nfc (np.ndarray, optional): Number of frequency per band. Must correspond with fc and freq. Defaults to None.

    Returns:
        _type_: _description_
    """

    final_coord = domain_facture(mesh, Ncore)
    nx1 = mesh.x_coord.shape[0]
    nx2 = mesh.x_coord.shape[1]

    # calculate only the first blade has reached the position of the second
    beta_array = np.linspace(0,2*np.pi*(wt.Nbeta-1)/wt.Nblade/wt.Nbeta,wt.Nbeta)#+np.pi/3

    # final results matrix to be field
    final_SPL_tot = np.zeros((nx1 * nx2,wt.Nseg*wt.Nblade, len(freq), wt.Nbeta))
    final_SPL_TIN = np.zeros((nx1 * nx2,wt.Nseg*wt.Nblade, len(freq), wt.Nbeta))
    final_SPL_TEN = np.zeros((nx1 * nx2,wt.Nseg*wt.Nblade, len(freq), wt.Nbeta))
    final_SWL_tot = np.zeros((nx1 * nx2,wt.Nseg*wt.Nblade, len(freq), wt.Nbeta))
    final_AoA = np.zeros((wt.Nseg*wt.Nblade, wt.Nbeta))
    final_Uinf = np.zeros((wt.Nseg*wt.Nblade, wt.Nbeta))
    final_epsilon = np.zeros((wt.Nseg*wt.Nblade, wt.Nbeta))
    final_Urot = np.zeros((wt.Nseg*wt.Nblade, wt.Nbeta))
    final_Urel = np.zeros((wt.Nseg*wt.Nblade, wt.Nbeta))
    final_a = np.zeros((wt.Nseg*wt.Nblade, wt.Nbeta))
    final_adash = np.zeros((wt.Nseg*wt.Nblade, wt.Nbeta))

    input_array = []
    start0 = time.time()
    for ibeta in range(wt.Nbeta):
        # offset angle
        wt.delta_beta = beta_array[ibeta]
        # loop over the sub domain
        for temp in range(Ncore):
            temp_dict = {}
            temp_dict['coord'] = final_coord[temp]
            temp_dict['wind_turbine'] = wt
            temp_dict['atmos'] = atmos
            temp_dict['freq'] = freq
            temp_dict['recon_index'] = temp
            input_array.append(temp_dict)

        # initialization of multiprocessing
        start = time.time()
        manager = multiprocessing.Manager()
        return_dict = manager.dict()
        jobs = []
        index_array = np.linspace(1, Ncore, Ncore)

        # launch of the threads
        for i in range(Ncore):
            p = multiprocessing.Process(target = SPL_SWL_function, args=(index_array[i], input_array[i], return_dict,BEM))
            jobs.append(p)
            p.start()

        for proc in jobs:
            proc.join()

        end = time.time()
        # end of the mulmtiprocdessing
        print(str(c_round(end - start)) + ' seconds')

        output_total = return_dict.values()
        # rearange output to the given shape
        final_output = reshuffle_output(output_total, Ncore)
        # construct the needed quantities
        # OASPL, OASPL_beta, AM_rec, SPL_freq, SWL_freq, SPL_tot_ff,Spp_tot,Spp_TIN,Spp_TEN, beta = output_construct(final_output, mesh, Ncore, len(freq), wt.Nblade, wt.Nseg)
        SPL_tot, SPL_TIN,SPL_TEN,SWL_tot,AoA,U_inf,epsilon,U_rot,U_rel,a,adash= output_construct(final_output, mesh, Ncore, len(freq), wt.Nblade, wt.Nseg)
        
        # recorded output over the change in blade position
        final_SPL_tot[:, :,:, ibeta] = SPL_tot
        final_SPL_TIN[:, :,:, ibeta] = SPL_TIN
        final_SPL_TEN[:, :,:, ibeta] = SPL_TEN
        final_SWL_tot[:, :,:, ibeta] = SWL_tot
        final_AoA[:, ibeta] = AoA
        final_Uinf[:,ibeta] = U_inf
        final_epsilon[:,ibeta] = epsilon
        final_Urot[:, ibeta] = U_rot
        final_Urel[:,ibeta] = U_rel
        final_a[:,ibeta] = a
        final_adash[:,ibeta] = adash
        
        final_beta_SPL_tot = 10*np.log10(np.sum(10**(final_SPL_tot/10), (1,3))/wt.Nbeta)

        print('Computations completed: ' + str(ibeta) + '/' + str(wt.Nbeta))
    end0 = time.time()
    print('total time: '+ str(c_round(end0-start0))+' seconds')

    # computation of third octave results and overall SPL results 
    if fc:
        SPL_third_tot = np.zeros((nx1 * nx2,wt.Nseg*wt.Nblade, len(fc), wt.Nbeta))
        for ifc in range(len(fc)): 
            SPL_third_tot[:,:,ifc,:] = 10*np.log10(fc[ifc]*0.232/Nfc[ifc]* np.sum(10 ** (final_SPL_tot[:,:, int(np.sum(Nfc[0:ifc])):int(np.sum(Nfc[0:ifc+1])),:]/10), axis=2))
        OASPL = 10*np.log10(np.sum(10**(SPL_third_tot/10),2))
        dBa = Aweight(np.array(fc)).reshape(1,1,-1,1)
        OASPL_A = 10*np.log10(np.sum(10**((SPL_third_tot+dBa)/10),2))
        return final_SPL_tot,final_SPL_TIN, final_SPL_TEN, final_SWL_tot,SPL_third_tot,OASPL,OASPL_A,final_AoA,final_Uinf

    return final_SPL_tot,final_SPL_TIN, final_SPL_TEN, final_SWL_tot, final_AoA,final_Uinf,final_epsilon,final_Urot,final_Urel,final_a,final_adash

