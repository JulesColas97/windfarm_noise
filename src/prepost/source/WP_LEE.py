from cmath import isnan, nan
import numpy as np
import scipy as sci
import logging

def WP_LEE(freq, delta_star, theta_momen, dpdx, cf, UeUinf, Uinf):
    # global rho0 nu0
    rho0 = 1.225
    nu0 = 1.45e-5

    # % WPS model of Lee et al. (2018)
    #calculate input parameters for the LEE model
    Ue = Uinf * np.abs(UeUinf)                                                 # external velocity (m)
    tau_w = np.abs(cf) * 0.5 * rho0 * Ue ** 2                                  # wall shear stress (Pa)
    tau_max = tau_w                                                            # maximum shear stress (Pa)
    H = delta_star / theta_momen                                               # shape factor
    u_star = np.sqrt(tau_w / rho0)                                             # friction velocity (m/s)
    delta = theta_momen * (3.15 + 1.72 / (H-1)) + delta_star                   # BL thickness (Drela 1986 paper)
    beta_c = theta_momen / tau_w * np.abs(dpdx)                                # Clauser parameter
    
    Delta = delta / delta_star  
    if Delta<0:
        logging.warning('negative value of Delta')   

    PI_Lee = 0.8 * (beta_c + 0.5) ** (3/4)
    e_Lee=  3.7 + 1.5 * beta_c
    d_Lee = 4.76 * ((1.4 / Delta) ** 0.75) * (0.375 * e_Lee - 1)
    
    R_T = (delta / Ue) / (nu0 / u_star ** 2)                                   # ratio of outer to inner boundary layer time scales
    a_Lee = 2.82 * Delta ** 2 * ((6.13 * (Delta ** (-0.75))+ d_Lee) ** e_Lee) * (4.2 * PI_Lee / Delta + 1)
    hstar_Lee= min(3, (0.139 + 3.1043 * beta_c)) + 7
    if beta_c < 0.5:
        dstar_Lee = max(1.0, 1.5 * d_Lee)
    else:
        dstar_Lee = d_Lee

    # non-dimensional angular frequency
    omega_adim = 2 * np.pi * freq * delta_star / Ue                           

    #wall pressure spectrun in dimensional and non-dimensional forms LEE et al
    Phipp_LEE_adim = (max(a_Lee, (0.25 * beta_c - 0.52) * a_Lee) * omega_adim ** 2) / ((4.76 * omega_adim ** 0.75 + dstar_Lee) ** e_Lee +(omega_adim * 8.8 * R_T ** (-0.57)) ** hstar_Lee)
    Phipp_LEE = tau_max ** 2 * delta_star / Ue * Phipp_LEE_adim

    return Phipp_LEE, Phipp_LEE_adim
