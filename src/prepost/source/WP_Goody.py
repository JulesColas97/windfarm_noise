import numpy as np


def WP_Goody(frequ, delta_star, theta_momen, cf, UeUinf, U_rota):
    # Phi_pp1 is the wall pressure spectra, Phi_pp2 is the normalize WPS by
    # (U^3*delta_star*rho^2)
    
    nu0 = 1.45e-5
    rho0 = 1.225
    
    Ue = U_rota * np.abs(UeUinf)
    Const1 = 0.5
    Const2 = 3
    
    # Re_c = U*c(index)/nu
    # delta = c(index)*0.382/(U*c(index)/nu)^0.2
    # delta_star = c(index)*0.047*Re_c^-0.2
    # tau_w = 0.0233*rho*Uc^2*(nu/Uc/delta)^0.25
    # u_star = sqrt(tau_w/rho)
    
    
    H = delta_star / theta_momen
    delta = theta_momen * (3.15 + 1.72 / (H - 1)) + delta_star                     #Boundary layer thickness (also Drela 1986 paper )
    # tau_w = abs(cf)*0.5* 1.225*U_rota^2;%Pa wall share stress
    # use Ue instead of U_rota for wall shear stress
    tau_w = np.abs(cf) * 0.5 * 1.225 * Ue ** 2                                     # wall share stress
    u_star = np.sqrt(tau_w / rho0)
    R_T = (delta / Ue)/(nu0 / u_star ** 2)
    Const3 = 1.1 * R_T ** (-0.57)
    omega_tilta = 2 * np.pi * frequ * delta / Ue
    Phi_pp1 = tau_w ** 2 * delta * Const2 * omega_tilta ** 2 / ((omega_tilta ** 0.75 + Const1) ** 3.7 + (Const3 * omega_tilta) ** 7) / Ue
    # Phi_pp2 = Phi_pp1/rho^2/delta/U^3;
    Phi_pp2 = Phi_pp1 * Ue / tau_w ** 2 / delta

    return Phi_pp1, Phi_pp2
