import numpy as np
import os
# from scipy.io import loadmat
import h5py

class Atmosphere():
    """Class to read and used ABL data from LES. 
        The class provides methods to read velocity and temperature profile from LES. to define generic velocity profile. 
        The class also providess function to compute the turbulence dissipation rate $\epsilon$ from the LES data
        !!! warning 
            this is legacy code. 
            It is not used anymore in the `Source` class and is replaced by the `Les` class
            instead. 
            This was used when only a vertical profile was used as input for the source model. 
            Now an interpolation from the 3D les data is perfomed for each segment and blade position. 
            Maybe a different approach could be implememented where the flow data are given already interpolated 
            to the `Source` class in order to make the two class less intricated.  
            

    """
    U_inf : np.ndarray
    epsilon_Kol : np.ndarray
    z_coord : np.ndarray
    H_mos : float                                                                                                 #Sensible heat flux
    L_star_mos : float

    def __init__(self):
        return

    def read_profile_from_les(self,path,fname,ratio,xpos,ypos):
        print('reading velocity profile from : ' + path+fname )
        u = ratio*np.array(h5py.File(path+fname, 'r+')['var'])
        print(u.shape)
        #Calculate the size of the grid
        nx = np.size(u[1,1,:])
        ny = np.size(u[1,:,1])
        nz = np.size(u[:,1,1])
        nz_tot = nz+1
        z_i = 625
        alpha=0.25
        #Domain length
        L_x = 10*z_i
        L_y = 1*z_i
        L_z = 1*z_i

        #Roughness height at the ground
        z0lo = 1.67e-4*z_i
        kappa = 0.4
        u_star = 0.6

        #Define the grid
        self.x = np.linspace(0.000,L_x,nx)
        self.y = np.linspace(0.000,L_y,ny)
        z = np.linspace(0.00,L_z,nz)
        dz = z_i/(nz_tot-1)
        z_w2 = np.linspace(0,z_i,nz_tot) #Location of w
        z_w2 = z_w2[0:nz]
        z_uv2 = z_w2+z_w2[1]/2 ## Location of u,v

        self.U_inf = np.sum(u[:,:,xpos],1)/ny
        self.z_coord = z_uv2

    def constant_epsilon(self,epsilon):
            self.epsilon_Kol = epsilon*np.ones((len(self.z_coord)))
