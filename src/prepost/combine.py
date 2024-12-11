import numpy as np
import time

from .utils import (computeThirdOctaveFrequencies, 
                    Aweight, 
                    integrateThirdOctave,
                    atm_absorption,
                    integrateThirdOctave)


def combine(src,field,fc,Nfc):
    # 2D OASPL computation
    #---------------------------------------------------------------------------
    print('combining ...')
    proximity = src.wt.proximity(field.height)
    # src.wt.plotProximity(proximity)
    # print(proximity[:,6])
    shape_Spp = src.Spp[...,0,0,0].shape
    shape_deltaL = field.deltaL[...,0,0].shape
    if shape_deltaL != shape_Spp :
        print('Spp and delta L must be of same size')
        print(shape_deltaL)
        print(shape_Spp)
        quit()
    if len(shape_Spp)==2 :
        (nx,nz) = shape_Spp
        adDim = 1
    elif len(shape_Spp)==1 :
        (nx,) = shape_Spp
        nz = 1
        adDim = 0
    else :
        print('space dimension not handle bny the actual code ')
        quit()
   
    final_OASPL_time = np.zeros(np.shape(src.Spp[...,0,0,:]))
    final_OASPL_A_time = np.zeros(np.shape(src.Spp[...,0,0,:]))
    final_SPL_freq_time =  np.zeros(np.shape(src.Spp[...,0,:,:]))
    final_SPL_seg_freq_time =  np.zeros(np.shape(src.Spp))
    final_SPL_freq_third_time = np.zeros(np.shape(src.Spp[...,0,0,0]) + (len(fc), src.wt.Nbeta))

    for ibeta in np.arange(0,src.wt.Nbeta):
        print('ibeta = ' + str(ibeta))
        # create coresponding matrix for Spp and Delta L fo onre angle 
        #---------------------------------------------------------------------------------------------------
        SPL_tot_ff = src.Spp[..., ibeta]
        delta_L = np.zeros_like(SPL_tot_ff)
        for iseg in range(src.wt.Nblade*src.wt.Nseg):
            print('iseg = ' + str(iseg) )
            print('delta L index : ' + str(proximity[iseg,ibeta]))
            print(delta_L.shape)
            delta_L[...,iseg, :] = field.deltaL[..., proximity[iseg,ibeta]]
        t0 = time.time()

        #   sum Spp and Delta L
        #--------------------------------------------------------------------------------------------------
        SPL_seg_freq = SPL_tot_ff + delta_L 
        final_SPL_seg_freq_time[...,ibeta] = SPL_seg_freq
        t_combine = time.time()
        pp_seg_freq = 10 ** (SPL_seg_freq / 10) * 4e-10
        print('combine in : ' + str(t_combine - t0) + 's.')

        # Compute third octave band
        #------------------------------------------------------------------------------------------------
        SPL_seg_freq_third = np.zeros(np.shape(src.Spp[...,0,0]) + (len(fc),))
        for ifreq in range(len(Nfc)):
            temp_array = 10*np.log10(fc[ifreq]*0.232/Nfc[ifreq]* np.sum(10 ** (SPL_seg_freq[..., int(np.sum(Nfc[0:ifreq])):int(np.sum(Nfc[0:ifreq+1]))]/10), axis=len(SPL_seg_freq.shape)-1))
            # temp_array[temp_array == -inf] = 0
            SPL_seg_freq_third[..., ifreq] = temp_array
        pp_seg_freq_third = 10 ** (SPL_seg_freq_third / 10) * 4e-10

        # compute sum over segments
        #--------------------------------------------------------------------------------------------------
        pp_freq = np.sum(pp_seg_freq.reshape(nx*nz,src.Nseg*src.Nblade,src.Nfreq), axis = 1)
        SPL_freq = np.squeeze(10 * np.log10(pp_freq / 4e-10).reshape(nx,nz,src.Nfreq))

        pp_freq_third = np.sum(pp_seg_freq_third.reshape(nx*nz,src.Nseg*src.Nblade,len(Nfc)), axis = 1)
        SPL_freq_third = np.squeeze(10 * np.log10(pp_freq_third/ 4e-10))
        t_sseg = time.time()
        print('sum over segments in ' + str(t_sseg - t_combine) + 's.')

        # absorpiton atmospherique 
        #-----------------------------------------------------------------------------------------------
        R = np.sqrt((np.squeeze(field.x_grid))**2 + (np.squeeze(field.z_grid)- src.wt.href)**2)
        c0 =343 
        rGP = 287.06;
        rho0 = 1.2# masse volumique de l'air au sol
        gamma = 1.4# constante adiabatique

        rh = 80# humidite relative (%)
        p0 = rho0*c0**2/gamma# pression atmospherique
        T0 = c0**2/(gamma*rGP)# t
        alpha = atm_absorption(T0,p0,rh
        ,fc).reshape(1,-1)*R.reshape(nx*nz,-1)
        SPL_freq_third = SPL_freq_third - alpha 


        # compute OASPL and OASPL_A
        #----------------------------------------------------------------------------------------------
        freq = computeThirdOctaveFrequencies(fc,Nfc)
        aweight = Aweight(fc).reshape(1,-1)

        OASPL_A = 10*np.log10(np.squeeze(np.sum(10**((SPL_freq_third+aweight)/10),axis = 1)))
        OASPL = 10*np.log10(np.squeeze(np.sum(10**(SPL_freq_third/10),axis = 1)))
        # if adDim == 0 : 
        OASPL_A = np.squeeze(OASPL_A.reshape(nx,nz))
        OASPL = np.squeeze(OASPL.reshape(nx,nz))
        SPL_freq_third = np.squeeze(SPL_freq_third.reshape(nx,nz,len(fc)))
        SPL_freq = np.squeeze(SPL_freq.reshape(nx,nz,len(freq)))

        t_sfreq = time.time()
        print('sum over frequency in ' + str(t_sfreq - t_sseg) + 's.')

        # save value 
        #-------------------------------------------------------------------------------------------------
        final_OASPL_time[..., ibeta] = OASPL
        final_OASPL_A_time[..., ibeta] = OASPL_A
        final_SPL_freq_time[..., ibeta] = SPL_freq
        final_SPL_freq_third_time[...,ibeta] = SPL_freq_third
        t_copy = time.time()
        print('copying in  ' + str(t_copy - t_sfreq) + 's.')
    final_SPL_seg_freq_time = final_SPL_seg_freq_time.reshape(nx,nz,src.Nseg,src.Nblade,len(freq),src.Nbeta)
    return final_SPL_freq_time, final_OASPL_time, final_OASPL_A_time

def combine_linear(src,field,fc,Nfc):
    # 2D OASPL computation
    #---------------------------------------------------------------------------
    print('combining ...')
    proximity = src.wt.proximityLinear(field.height)
    # src.wt.plotProximity(proximity)
    # print(proximity[:,6])
    shape_Spp = src.SppInterpolated[...,0,0,0].shape
    shape_deltaL = field.deltaLInterpolated[...,0,0].shape
    if shape_deltaL != shape_Spp :
        print('Spp and delta L must be of same size')
        print(shape_deltaL)
        print(shape_Spp)
        quit()
    if len(shape_Spp)==2 :
        (nx,nz) = shape_Spp
        adDim = 1
    elif len(shape_Spp)==1 :
        (nx,) = shape_Spp
        nz = 1
        adDim = 0
    else :
        print('space dimension not handle bny the actual code ')
        quit()
   
    final_OASPL_time = np.zeros(np.shape(src.SppInterpolated[...,0,0,:]))
    final_OASPL_A_time = np.zeros(np.shape(src.SppInterpolated[...,0,0,:]))
    final_SPL_freq_time =  np.zeros(np.shape(src.SppInterpolated[...,0,:,:]))
    final_SPL_seg_freq_time =  np.zeros(np.shape(src.SppInterpolated))
    final_SPL_freq_third_time = np.zeros(np.shape(src.SppInterpolated[...,0,0,0]) + (len(fc), src.wt.Nbeta))

    for ibeta in np.arange(0,src.wt.Nbeta):
        print('ibeta = ' + str(ibeta))
        # create coresponding matrix for Spp and Delta L fo one angle 
        #---------------------------------------------------------------------------------------------------
        SPL_tot_ff = src.SppInterpolated[..., ibeta]
        delta_L = np.zeros_like(SPL_tot_ff)
        for iseg in range(src.wt.Nblade*src.wt.Nseg):
            print('iseg = ' + str(iseg) )
            print('delta L index : ' + str(proximity[iseg,ibeta,0]) + ',' + str(proximity[iseg,ibeta,1])+','+str(proximity[iseg,ibeta,2]))

            delta_L[...,iseg, :] = 10*np.log10(proximity[iseg,ibeta,2]*10**(field.deltaLInterpolated[..., proximity[iseg,ibeta,0]]/10) + 
                                        (1-proximity[iseg,ibeta,2])*10**(field.deltaLInterpolated[..., proximity[iseg,ibeta,1]]/10))
        t0 = time.time()

        #   sum Spp and Delta L
        #--------------------------------------------------------------------------------------------------
        SPL_seg_freq = SPL_tot_ff + delta_L 
        final_SPL_seg_freq_time[...,ibeta] = SPL_seg_freq
        t_combine = time.time()
        pp_seg_freq = 10 ** (SPL_seg_freq / 10) * 4e-10
        print('combine in : ' + str(t_combine - t0) + 's.')

        # Compute third octave band
        #------------------------------------------------------------------------------------------------
        SPL_seg_freq_third = np.zeros(np.shape(src.Spp[...,0,0]) + (len(fc),))
        for ifreq in range(len(Nfc)):
            temp_array = 10*np.log10(fc[ifreq]*0.232/Nfc[ifreq]* np.sum(10 ** (SPL_seg_freq[..., int(np.sum(Nfc[0:ifreq])):int(np.sum(Nfc[0:ifreq+1]))]/10), axis=len(SPL_seg_freq.shape)-1))
            # temp_array[temp_array == -inf] = 0
            SPL_seg_freq_third[..., ifreq] = temp_array
        pp_seg_freq_third = 10 ** (SPL_seg_freq_third / 10) * 4e-10

        # compute sum over segments
        #--------------------------------------------------------------------------------------------------
        pp_freq = np.sum(pp_seg_freq.reshape(nx*nz,src.Nseg*src.Nblade,src.Nfreq), axis = 1)
        SPL_freq = np.squeeze(10 * np.log10(pp_freq / 4e-10).reshape(nx,nz,src.Nfreq))

        pp_freq_third = np.sum(pp_seg_freq_third.reshape(nx*nz,src.Nseg*src.Nblade,len(Nfc)), axis = 1)
        SPL_freq_third = np.squeeze(10 * np.log10(pp_freq_third/ 4e-10))
        t_sseg = time.time()
        print('sum over segments in ' + str(t_sseg - t_combine) + 's.')

        # absorpiton atmospherique 
        #-----------------------------------------------------------------------------------------------
        R = np.sqrt((np.squeeze(field.x_grid))**2 + (np.squeeze(field.z_grid)- src.wt.href)**2)
        c0 =343 
        rGP = 287.06;
        rho0 = 1.2# masse volumique de l'air au sol
        gamma = 1.4# constante adiabatique

        rh = 80# humidite relative (%)
        p0 = rho0*c0**2/gamma# pression atmospherique
        T0 = c0**2/(gamma*rGP)# t
        alpha = atm_absorption(T0,p0,rh,fc).reshape(1,-1)*R.reshape(nx*nz,-1)
        SPL_freq_third = SPL_freq_third - alpha 


        # compute OASPL and OASPL_A
        #----------------------------------------------------------------------------------------------
        freq = computeThirdOctaveFrequencies(fc,Nfc)
        aweight = Aweight(fc).reshape(1,-1)

        OASPL_A = 10*np.log10(np.squeeze(np.sum(10**((SPL_freq_third+aweight)/10),axis = 1)))
        OASPL = 10*np.log10(np.squeeze(np.sum(10**(SPL_freq_third/10),axis = 1)))
        # if adDim == 0 : 
        OASPL_A = np.squeeze(OASPL_A.reshape(nx,nz))
        OASPL = np.squeeze(OASPL.reshape(nx,nz))
        SPL_freq_third = np.squeeze(SPL_freq_third.reshape(nx,nz,len(fc)))
        SPL_freq = np.squeeze(SPL_freq.reshape(nx,nz,len(freq)))

        t_sfreq = time.time()
        print('sum over frequency in ' + str(t_sfreq - t_sseg) + 's.')

        # save value 
        #-------------------------------------------------------------------------------------------------
        final_OASPL_time[..., ibeta] = OASPL
        final_OASPL_A_time[..., ibeta] = OASPL_A
        final_SPL_freq_time[..., ibeta] = SPL_freq
        final_SPL_freq_third_time[...,ibeta] = SPL_freq_third
        t_copy = time.time()
        print('copying in  ' + str(t_copy - t_sfreq) + 's.')
    final_SPL_seg_freq_time = final_SPL_seg_freq_time.reshape(nx,nz,src.Nseg,src.Nblade,len(freq),src.Nbeta)
    return final_SPL_freq_time, final_OASPL_time, final_OASPL_A_time


def combine_linear_broadband(src,field,freq):
    # 2D OASPL computation
    #---------------------------------------------------------------------------
    print('combining ...')
    proximity = src.wt.proximityLinear(field.height)
    # src.wt.plotProximity(proximity)
    # print(proximity[:,6])
    shape_Spp = src.SppInterpolated[...,0,0,0].shape
    shape_deltaL = field.deltaLInterpolated[...,0,0].shape
    if shape_deltaL != shape_Spp :
        print('Spp and delta L must be of same size')
        print(shape_deltaL)
        print(shape_Spp)
        quit()
    if len(shape_Spp)==2 :
        (nx,nz) = shape_Spp
        ny=1
        adDim = 1
    if len(shape_Spp)==3 :
        print('test')
        (nx,ny,nz) = shape_Spp
        adDim = 1
    elif len(shape_Spp)==1 :
        (nx,) = shape_Spp
        nz=1
        ny=1
        adDim = 0
    else :
        print('space dimension not handle bny the actual code ')
        quit()

    final_OASPL_time = np.zeros(np.shape(src.SppInterpolated[...,0,0,:]))
    final_OASPL_A_time = np.zeros(np.shape(src.SppInterpolated[...,0,0,:]))
    final_SPL_freq_time =  np.zeros(np.shape(src.SppInterpolated[...,0,:,:]))
    final_SPL_seg_freq_time =  np.zeros(np.shape(src.SppInterpolated))
    final_SPL_freq_third_time = np.zeros(np.shape(src.SppInterpolated[...,0,0,0]) + (14, src.wt.Nbeta))

    for ibeta in np.arange(0,src.wt.Nbeta):
        print('ibeta = ' + str(ibeta))
        # create coresponding matrix for Spp and Delta L fo one angle 
        #---------------------------------------------------------------------------------------------------
        SPL_tot_ff = src.SppInterpolated[..., ibeta]
        delta_L = np.zeros_like(SPL_tot_ff)
        for iseg in range(src.wt.Nblade*src.wt.Nseg):
            print('iseg = ' + str(iseg) )
            print('delta L index : ' + str(proximity[iseg,ibeta,0]) + ',' + str(proximity[iseg,ibeta,1])+','+str(proximity[iseg,ibeta,2]))

            delta_L[...,iseg, :] = 10*np.log10(proximity[iseg,ibeta,2]*10**(field.deltaLInterpolated[..., proximity[iseg,ibeta,0]]/10) + 
                                        (1-proximity[iseg,ibeta,2])*10**(field.deltaLInterpolated[..., proximity[iseg,ibeta,1]]/10))
        t0 = time.time()

        #   sum Spp and Delta L
        #--------------------------------------------------------------------------------------------------
        SPL_seg_freq = SPL_tot_ff + delta_L 
        final_SPL_seg_freq_time[...,ibeta] = SPL_seg_freq
        t_combine = time.time()
        pp_seg_freq = 10 ** (SPL_seg_freq / 10) * 4e-10
        print('combine in : ' + str(t_combine - t0) + 's.')

        # Compute third octave band
        #------------------------------------------------------------------------------------------------
        # SPL_seg_freq_third = np.zeros(np.shape(src.Spp[...,0,0]) + (len(fc),))
        fc,pp_seg_freq_third= integrateThirdOctave(freq,10**(SPL_seg_freq/10))
        # SPL_seg_freq_third = 10*np.log10(pp_seg_freq_third)
        pp_seg_freq_third = pp_seg_freq_third * 4e-10
        # for ifreq in range(len(Nfc)):
        #     temp_array = 10*np.log10(fc[ifreq]*0.232/Nfc[ifreq]* np.sum(10 ** (SPL_seg_freq[..., int(np.sum(Nfc[0:ifreq])):int(np.sum(Nfc[0:ifreq+1]))]/10), axis=len(SPL_seg_freq.shape)-1))
        #     # temp_array[temp_array == -inf] = 0
        #     SPL_seg_freq_third[..., ifreq] = temp_array
        # pp_seg_freq_third = 10 ** (SPL_seg_freq_third / 10) * 4e-10

        # compute sum over segments
        #--------------------------------------------------------------------------------------------------
        pp_freq = np.sum(pp_seg_freq.reshape(nx*nz*ny,src.Nseg,src.Nfreq), axis = 1)
        SPL_freq = np.squeeze(10 * np.log10(pp_freq / 4e-10).reshape(nx,ny,nz,src.Nfreq))

        pp_freq_third = np.sum(pp_seg_freq_third.reshape(nx*nz*ny,src.Nseg,len(fc)), axis = 1)
        SPL_freq_third = np.squeeze(10 * np.log10(pp_freq_third/ 4e-10))
        t_sseg = time.time()
        print('sum over segments in ' + str(t_sseg - t_combine) + 's.')

        # absorpiton atmospherique 
        #-----------------------------------------------------------------------------------------------
        R = np.sqrt((np.squeeze(field.x_interpolate))**2 + (np.squeeze(field.z_interpolate)- src.wt.href)**2 + np.squeeze(field.y_interpolate)**2)
        c0 =343 
        rGP = 287.06;
        rho0 = 1.2# masse volumique de l'air au sol
        gamma = 1.4# constante adiabatique

        rh = 80# humidite relative (%)
        p0 = rho0*c0**2/gamma# pression atmospherique
        T0 = c0**2/(gamma*rGP)# t
        alpha = atm_absorption(T0,p0,rh,fc).reshape(1,-1)*R.reshape(nx*nz*ny,-1)
        SPL_freq_third = SPL_freq_third - alpha 


        # compute OASPL and OASPL_A
        #----------------------------------------------------------------------------------------------
        aweight = Aweight(fc).reshape(1,-1)

        OASPL_A = 10*np.log10(np.squeeze(np.sum(10**((SPL_freq_third+aweight)/10),axis = 1)))
        OASPL = 10*np.log10(np.squeeze(np.sum(10**(SPL_freq_third/10),axis = 1)))
        # if adDim == 0 : 
        OASPL_A = np.squeeze(OASPL_A.reshape(nx,ny,nz))
        OASPL = np.squeeze(OASPL.reshape(nx,ny,nz))
        SPL_freq_third = np.squeeze(SPL_freq_third.reshape(nx,ny,nz,len(fc)))
        SPL_freq = np.squeeze(SPL_freq.reshape(nx,ny,nz,len(freq)))

        t_sfreq = time.time()
        print('sum over frequency in ' + str(t_sfreq - t_sseg) + 's.')

        # save value 
        #-------------------------------------------------------------------------------------------------
        final_OASPL_time[..., ibeta] = OASPL
        final_OASPL_A_time[..., ibeta] = OASPL_A
        final_SPL_freq_time[..., ibeta] = SPL_freq
        final_SPL_freq_third_time[...,ibeta] = SPL_freq_third
        t_copy = time.time()
        print('copying in  ' + str(t_copy - t_sfreq) + 's.')
    final_SPL_seg_freq_time = final_SPL_seg_freq_time.reshape(nx,ny,nz,src.wt.Nseg,src.wt.Nblade,len(freq),src.Nbeta)
    return final_SPL_freq_time, final_SPL_freq_third_time, final_OASPL_time, final_OASPL_A_time






