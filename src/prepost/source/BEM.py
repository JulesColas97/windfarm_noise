"""
This functions allows to use Blade Element Momentum Theory to recursively compute the angle of attack on an airfoil. 
It takes into accound drag and lift coefficient previously computed using `prepost.source.wind_turbine.WindTurbine.createBLData` 

TODO add explanation here !!
"""
import scipy.interpolate as interp
import numpy as np

def simple(Cl_tck,Cd_tck,Re,J,theta,twist,blade_length,seg,F):
    a = 0
    adash = 0
    conv_limit = 1e-6
    imax = 100
    it = 0
    F = 1
    x = seg/blade_length
                
    a_conv = []
    adash_conv = []
    alpha_conv = []

    for it in range(imax):
        # phi = np.arctan((1-a)/(J*x*(1+adash)))
        # if a == 1: 
        #     phi = 0
        #     print('test')
        # else:
        #     phi = np.arctan((1-a)/(J*x*(1+adash)))
        phi = np.arctan2((1-a*F)/(J*x*(1+adash)),1)
        # phi = np.arctan(np.divide((1-a),(J*x*(1+adash)),out=np.zeros_like(1-a),where=(J*x*(1+adash)!=0)))
        alpha = (phi - np.pi/2  + twist)
        F = 2/np.pi * np.arccos(np.exp(-3*(blade_length - seg)/(2*seg*np.sin(phi))))
        # F = 1
        Cl = interp.bisplev(np.array([Re]),np.array([alpha])* 180 / np.pi,Cl_tck)
        Cd = interp.bisplev(np.array([Re]),np.array([alpha])* 180 / np.pi,Cd_tck)


        Cn = Cl*np.cos(phi) + Cd*np.sin(phi)
        Ct = Cl*np.sin(phi) - Cd*np.cos(phi)
        a_new = 1/(4*F*np.sin(phi)**2/(theta*Cn)+1)
        adash_new = 1/(4*F*np.sin(phi)*np.cos(phi)/(theta*Ct)-1)
        epsilon = ((a - a_new) ** 2 + (adash - adash_new) ** 2) ** 0.5
        if epsilon < conv_limit:
            break
        a = a_new
        adash = adash_new
        it += 1
   
    return alpha,a_new,adash_new,F,epsilon,a_conv,adash_conv,alpha_conv


# def hemant(Cl_tck,Re,J,F,Ltot,seg,twist,chord):
#     a = 0
#     a_dash = 0
#     conv_limit = 1e-6
#     imax = 100
#     Z = 3
#     a = 0; a_dash = 0
#     for it in range(imax):
#         prev_a = a; prev_a_dash = a_dash
#         phi = np.arctan2(Ltot/(seg*J) * (1 - a * F) / (1 + a_dash), 1)

#         # F = 2/np.pi * np.arccos(np.exp(-3*(Ltot - seg)/(2*seg*np.sin(phi))))
#         # alpha = phi - twist
#         alpha = (phi - np.pi/2 + twist)
#         Cl = interp.bisplev(np.array([Re]),np.array([alpha])* 180 / np.pi,Cl_tck)
#         lambda_coef = Z * chord* Cl / (8 * np.pi * seg)
#         a = (1 + F / lambda_coef * np.tan(phi) * np.sin(phi)) ** (-1)
#         a_dash = (F / lambda_coef * np.cos(phi) - 1) ** (-1)
#         epsilon = ((a - prev_a) ** 2 + (a_dash - prev_a_dash) ** 2) ** 0.5
#         if epsilon < conv_limit:
#             break
#     return alpha,a,a_dash



# use only Cl bby considering Cd/Cl<<<1
def hemant1(Cl_tck,Cd_tck,Re,J,R,r,twist,c):
    conv_limit = 1e-6
    imax = 100
    Z = 3
    a = 0; a_dash = 0
    F = 2 / np.pi * np.arccos(np.exp(-(3 / 2) * (1 - r/ R) * (1 + J ** 2) ** 0.5))

    a_conv = []
    adash_conv = []
    alpha_conv = []

    for it in range(imax):
        prev_a = a; prev_a_dash = a_dash
        phi = np.arctan2(R/(r*J) * (1 - a * F) / (1 + a_dash), 1)
        alpha = (phi - np.pi/2 + twist)

        Cl = interp.bisplev(np.array([Re]),np.array([alpha])* 180/np.pi,Cl_tck)
        lambda_coef = Z * c* Cl / (8 * np.pi * r)

        a = (1 + F / lambda_coef * np.tan(phi) * np.sin(phi)) ** (-1)
        a_dash = (F / lambda_coef * np.cos(phi) - 1) ** (-1)
        
        a_conv.append(a)
        adash_conv.append(a_dash)
        alpha_conv.append(alpha)

        epsilon = ((a - prev_a) ** 2 + (a_dash - prev_a_dash) ** 2) ** 0.5
        if epsilon < conv_limit:
            break
    eps = 0
    if it == imax-1 : 
        a = np.nan
        a_dash = np.nan
        alpha = np.nan 

    return alpha,a,a_dash,F,eps,a_conv,adash_conv,alpha_conv


# use of Cd and Cl 
def hemant2(Cl_tck,Cd_tck,Re,J,R,r,twist,c):
    a = 0
    a_dash = 0
    conv_limit = 1e-6
    imax = 100
    Z = 3
    a = 0; a_dash = 0
    F = 2 / np.pi * np.arccos(np.exp(-(3 / 2) * (1 - r/ R) * (1 + J ** 2) ** 0.5))
    a_conv = []
    adash_conv = []
    alpha_conv = []
    for it in range(imax):
        prev_a = a; prev_a_dash = a_dash
         
        phi = np.arctan2(R/(r*J) * (1 - a * F) / (1 + a_dash), 1)
        alpha = (phi - np.pi/2 + twist)
        Cl = interp.bisplev(np.array([Re]),np.array([alpha])* 180 / np.pi,Cl_tck)
        Cd = interp.bisplev(np.array([Re]),np.array([alpha])* 180 / np.pi,Cd_tck)
        # F = 2/np.pi * np.arccos(np.exp(-3*(1 - r/R)/(2*np.sin(phi))))
        lambda_coef = Z * c* Cl / (8 * np.pi * r)
        eps = Cd/Cl
        # eps = 0

        a = (1 + F*np.sin(phi)**2/(lambda_coef*(np.cos(phi) + eps*np.sin(phi))))**(-1)
        a_dash = (-1 + F*np.sin(phi)*np.cos(phi)/(lambda_coef*(np.sin(phi) - eps*np.cos(phi))))**(-1)

        a_conv.append(a)
        adash_conv.append(a_dash)
        alpha_conv.append(alpha)


        epsilon = ((a - prev_a) ** 2 + (a_dash - prev_a_dash) ** 2) ** 0.5
        if epsilon < conv_limit:
            break
    if it ==imax:
        print('not converged')

    return alpha,a,a_dash,F,eps, a_conv,adash_conv,alpha_conv



#  modification of the reynold number at each step 
def hemant3(Cl_tck,Cd_tck,Re,J,R,r,twist,c,Uinf,Urot):
    nu0 = 1.45e-5  
    a = 0
    a_dash = 0
    conv_limit = 1e-6
    imax = 99
    Z = 3
    a = 0; a_dash = 0
    F = 2 / np.pi * np.arccos(np.exp(-(3 / 2) * (1 - r/ R) * (1 + J ** 2) ** 0.5))
    a_conv = []
    adash_conv = []
    alpha_conv = []
    for it in range(imax):
        prev_a = a; prev_a_dash = a_dash
         
        phi = np.arctan2(R/(r*J) * (1 - a * F) / (1 + a_dash), 1)
        alpha = (phi - np.pi/2 + twist)
        Cl = interp.bisplev(np.array([Re]),np.array([alpha])* 180 / np.pi,Cl_tck)
        Cd = interp.bisplev(np.array([Re]),np.array([alpha])* 180 / np.pi,Cd_tck)
        # F = 2/np.pi * np.arccos(np.exp(-3*(1 - r/R)/(2*np.sin(phi))))
        lambda_coef = Z * c* Cl / (8 * np.pi * r)
        eps = Cd/Cl
        # eps = 0

        a = (1 + F*np.sin(phi)**2/(lambda_coef*(np.cos(phi) + eps*np.sin(phi))))**(-1)
        a_dash = (-1 + F*np.sin(phi)*np.cos(phi)/(lambda_coef*(np.sin(phi) - eps*np.cos(phi))))**(-1)

        a_conv.append(a)
        adash_conv.append(a_dash)
        alpha_conv.append(alpha)
        
        U_rel= np.sqrt((Uinf* (1 - F * a)) ** 2 + (Urot * (1 + a_dash)) ** 2)
        Re = U_rel* c / nu0


        epsilon = ((a - prev_a) ** 2 + (a_dash - prev_a_dash) ** 2) ** 0.5
        if epsilon < conv_limit:
            break
    if it ==imax:
        print('not converged')

    return alpha,a,a_dash,F,eps, a_conv,adash_conv, alpha_conv

def noBEM(Cl_tck,Cd_tck,Re,J,R,r,twist,c):
    a = 0; a_dash = 0
    phi = np.arctan2(R/(r*J), 1)
    alpha = (phi - np.pi/2 + twist)
    eps = 0
    F = 1
    return alpha,a,a_dash,F,eps
