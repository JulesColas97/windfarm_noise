import os
import numpy as np
import h5py
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


class Les:

    def __init__(self,path):
        self.path = path

    def read(self,dirname="output", ratio=None, adim=False):
        fname = self.path + '/%s/gridz.h5'%(dirname)
        vel = h5py.File(fname,"r")
        z = np.array(vel['var'])
        self.z = np.asarray(z)
        fname =  self.path + '/%s/gridzw.h5'%(dirname)
        vel = h5py.File(fname,"r")
        zw = np.array(vel['var'])
        self.zw = np.asarray(zw)
        i = 0
        tstart = 0
        lines = []
        in_file = open (self.path+'/input.conf', 'rt')
        self.u_star = 1
        for line in in_file:
            lines.append(line)
            dummy = lines[i].split()
            if (len(dummy)> 0):
                if(dummy[0] == 'Nx'):
                    self.nx = int(dummy[2])
                if(dummy[0] == 'Ny'):
                    self.ny = int(dummy[2])
                if(dummy[0] == 'Nz'):
                    self.nz = int(dummy[2])
                if(dummy[0] == 'z_i'):
                    self.z_i = float(dummy[2])
                if(dummy[0] == 'grid_size'):
                    self.dz = float(dummy[2])
                if(dummy[0] == 'Lx'):
                    self.lx = float(dummy[2])
                if(dummy[0] == 'Ly'):
                    self.ly = float(dummy[2])
                if(dummy[0] == 'Lz'):
                    self.lz = float(dummy[2])
                if(dummy[0] == 'ug'):
                    self.ug = float(dummy[2])
                if (dummy[0] == 'u_star'):
                    self.u_star = float(dummy[2])
                if (dummy[0] == 'zo'):
                    self.z0 = float(dummy[2])
                    self.z0 = self.z0*self.z_i
                if (dummy[0] == 'use_mean_p_force'):
                    if dummy[2] == '.false.':
                        self.use_mean_p_force = False
                    else:
                        self.use_mean_p_force = True
                if(dummy[0] == 'T_scale'):
                    self.T_scale = float(dummy[2])
                if(dummy[0] == 'T_init'):
                    self.T_init = float(dummy[2])
                if(dummy[0] == 'nloc'):
                    self.nturb = int(dummy[2])

            i = i + 1

        self.dx  = (self.lx / np.double(self.nx)) * self.z_i
        self.dy  = (self.ly / np.double(self.ny)) * self.z_i

        self.L_x = self.lx * self.z_i
        self.L_y = self.ly * self.z_i
        self.L_z = self.lz * self.z_i


        #Define the grid 
        self.dx = self.L_x/self.nx
        self.dy = self.L_y/self.ny
        self.x = np.linspace(0.000,self.L_x-self.dx,self.nx)
        self.y = np.linspace(0.000,self.L_y-self.dy,self.ny)
        self.z = self.z * self.z_i
        self.zw = self.zw * self.z_i

        if self.use_mean_p_force:
            self.u_dim = 1
            # print("ustar =%s"%self.u_star)
        else:
            self.u_dim = self.ug
            # print("ug =%s"%self.ug)

        self.adim = adim
        # dimension of flow
        if ratio is not None:
            self.u_dim = ratio
        if self.adim:
            self.u_dim = 1/self.u_star

        fullpath = self.path +'/' + dirname + '/'
        if os.path.isfile(fullpath+'tavg_u.h5'):
            vel = h5py.File(fullpath+'tavg_u.h5',"r")
            u = self.u_dim*np.array(vel['var'])
            self.u = np.asarray(u)
        if os.path.isfile(fullpath+'tavg_u.h5'):  
            vel = h5py.File(fullpath+'tavg_v.h5',"r")
            v = self.u_dim*np.array(vel['var'])
            self.v = np.asarray(v)
        if os.path.isfile(fullpath+'tavg_w.h5'):
            vel = h5py.File(fullpath+'tavg_w.h5',"r")
            w = self.u_dim*np.array(vel['var'])
            self.w = np.asarray(w)
        if os.path.isfile(fullpath+'tavg_theta.h5'):
            vel = h5py.File(fullpath+'tavg_theta.h5',"r")
            theta = self.T_scale*np.array(vel['var'])
            self.theta = np.asarray(theta)


        if os.path.isfile(fullpath+'tavg_pres.h5'):
            vel = h5py.File(fullpath+'tavg_pres.h5',"r")
            self.pres = np.array(vel['var'])

        else:
            self.theta = None
        if os.path.isfile(fullpath+'initial_theta.h5'):
            vel = h5py.File(fullpath+'initial_theta.h5',"r")
            theta = self.T_scale*np.array(vel['var'])
            self.initial_theta = np.asarray(theta)



        if self.theta is not None:
            # print(np.max(self.pres))
            # print(np.min(self.pres))
            # convert potential temperature to absolute temperature 
            R_over_cp = 0.286
            P0 = 1e5    # pression atmospherique de reference
            P =  101325 # reference pressure
            # Mass molar of air  
            M = 29e-3 # Kg/mol
            g = 9.81 # m/s2
            # universal gaz constant 
            R = 8.314 # J/(mol;K)
            T0 = 293 # reference temperature 

            p = P0 * np.exp(- self.z[:,None,None] * M*g/(R*T0))
            self.T = self.theta*(p/P0)**(R_over_cp)
        

            # self.T = self.theta*(P/P0)**(R_over_cp)

            gamma = 1.4 
            R = 286
            self.c = np.sqrt(gamma*R*self.T) 

        turbines = []
        if os.path.isfile(self.path + '/turbine_input.dat'):
            with open(self.path + '/turbine_input.dat') as f:
                for l in f :
                    sp = l.split(',')
                    turbines.append(sp)

            turbines.pop(0)
        self.turbines = np.array(turbines, dtype=float)
        self.epsilon = None
        self.ti = None

    # def add_first_point():


    def plotQuantityXZ(self,fname,y,**kwargs):
        iy = np.argmin(np.abs(self.y-y))

        if fname == 'u' :
            vel =  self.u
        elif fname == 'v' :
            vel = self.v
        elif fname == 'w':
            vel = self.w
        elif fname == 'c':
            vel = self.c
        elif fname == 'theta':
            vel = self.theta
        elif fname == 'T':
            vel = self.T
        elif fname == 'epsilon':
            if self.epsilon is None:
                self.dissipation()
            vel = self.epsilon
        else:
            vel = h5py.File(self.path + fname,"r")
            vel = np.array(vel['var'])

        cax = plt.pcolormesh(self.x,self.z,vel[:,iy,:],**kwargs)
        # plt.gca().set_aspect('equal', adjustable='box')
        # plt.xlabel('x (m)')
        # plt.ylabel('z (m)')
        # ax=plt.gca()
        # # plt.rcParams["contour.negative_linestyle"] = "solid"
        # divider = make_axes_locatable(ax)
        # cb = divider.append_axes("right", size="5%", pad=0.05)
        # cb.set_title(fname)
        # plt.colorbar(cax, cax=cb)
        # plt.tight_layout()
        return cax

    def plotQuantityXY(self,fname,z,**kwargs):
        iz = np.argmin(np.abs(self.z-z))

        if fname == 'u' :
            vel =  self.u
        elif fname == 'v' :
            vel = self.v
        elif fname == 'w':
            vel = self.w
        elif fname == 'c':
            vel = self.c
        elif fname == 'theta':
            vel = self.theta
        elif fname == 'T':
            vel = self.T
        elif fname == 'epsilon':
            if self.epsilon is None:
                self.dissipation()
            vel = self.epsilon
        elif fname == 'ti':
            if self.ti is None:
                self.turbulent_intensity()
            vel = self.ti
        else:
            vel = h5py.File(self.path + fname,"r")
            vel = np.array(vel['var'])

        cax = plt.pcolormesh(self.x,self.y,vel[iz,:,:],**kwargs)
        # plt.gca().set_aspect('equal', adjustable='box')
        plt.xlabel('x (m)')
        plt.ylabel('y (m)')
        # ax=plt.gca()
        # # plt.rcParams["contour.negative_linestyle"] = "solid"
        # divider = make_axes_locatable(ax)
        # cb = divider.append_axes("right", size="5%", pad=0.05)
        # cb.set_title(fname)
        # plt.colorbar(cax, cax=cb)
        # plt.tight_layout()
        # plt.rcParams["contour.negative_linestyle"] = "solid"
        #divider = make_axes_locatable(ax)
        #cb = divider.append_axes("right", size="5%", pad=0.05)
        #cb.set_title(fname)
        #plt.colorbar(cax, cax=cb)
        plt.tight_layout()
        return cax


    def plotQuantityYZ(self,fname,x,**kwargs):
        ix = np.argmin(np.abs(self.x-x))

        if fname == 'u' :
            vel =  self.u
        elif fname == 'v' :
            vel = self.v
        elif fname == 'w':
            vel = self.w
        elif fname == 'c':
            vel = self.c
        elif fname == 'theta':
            vel = self.theta
        elif fname == 'T':
            vel = self.T
        elif fname == 'epsilon':
            if self.epsilon is None:
                self.dissipation()
            vel = self.epsilon
        elif fname == 'ti':
            if self.ti is None:
                self.turbulent_intensity()
            vel = self.ti
        else:
            vel = h5py.File(self.path + fname,"r")
            vel = np.array(vel['var'])

        cax = plt.pcolormesh(self.y,self.z,vel[:,:,ix],**kwargs)
        # plt.gca().set_aspect('equal', adjustable='box')
        plt.xlabel('y (m)')
        plt.ylabel('z (m)')
        ax=plt.gca()
        # # plt.rcParams["contour.negative_linestyle"] = "solid"
        # divider = make_axes_locatable(ax)
        # cb = divider.append_axes("right", size="5%", pad=0.05)
        # cb.set_title(fname)
        # plt.colorbar(cax, cax=cb)
        plt.tight_layout()
        return ax

    def plotProfile(self, fname, x=None, y=None, scaling=1, **kwargs):
        if fname == 'u':
            vel = self.u
            unit = "ms"
        elif fname == 'v':
            vel = self.v
            unit = "ms"
        elif fname == "wind_angle":
            vel = np.arctan(self.v, self.u)*180/np.pi
        elif fname == 'w':
            vel = self.w
            unit = "ms"
        elif fname == 'c':
            vel = self.c
            unit = "ms"
        elif fname == 'theta':
            vel = self.theta
            unit = "K"
        elif fname == 'T':
            vel = self.T
            unit = "K"
        elif fname == 'epsilon':
            if self.epsilon is None:
                self.dissipation()
            vel = self.epsilon
            unit = "m2/s3"
        elif fname == 'ti':
            if self.ti is None:
                self.turbulent_intensity()
            vel = self.ti
            unit = "(-)"
        elif fname == 'initial_theta':
            vel = self.initial_theta
        else:
            vel = h5py.File(self.path + fname, "r")
            if ('_u' in fname) or ('_v' in fname) or ('_w' in fname):
                vel = self.u_dim * np.array(vel['var'])
                unit = 'm/s'
            elif 'theta' in fname:
                vel = self.T_scale * np.array(vel['var'])
                unit = 'K'
            else:
                vel = np.array(vel['var'])
                unit = 'm/s'

        if len(vel.shape) == 3:
            if x is not None:
                ix = np.argmin(np.abs(x - self.x))
                vel_profile = vel[:, :, ix]
            else:
                vel_profile = np.mean(vel, 2)

            if y is not None:
                iy = np.argmin(np.abs(y - self.y))
                vel_profile = vel_profile[:, iy]
            else:
                vel_profile = np.mean(vel_profile, 1)
        else:
            vel_profile = vel
        vel_profile = vel_profile[0:self.z.size]

        if fname == "epsilon" or fname == "ti":
            plt.plot(scaling*vel_profile[1:], self.zw[1:], **kwargs)
        else:
            plt.plot(scaling*vel_profile, self.z, **kwargs)

        # plt.plot(les.theta[:,10,10],les.z_coord)
        # plt.ylabel('$z$ (m)')
        # if self.adim:
        #     plt.xlabel(fname + ' (-)')
        # else:
        #     plt.xlabel(fname + ' ('+unit+')')
        return vel_profile

    def readConvergence(self):
        self.t = []
        self.ke = []
        self.friction = []
        self.temp = []

        in_file = open(self.path+'/output/check_ke.dat', 'rt')
        for line in in_file:
            dummy = line.split()
            self.t.append(float(dummy[0]))
            self.ke.append(float(dummy[1]))
            self.friction.append(float(dummy[2]))
            if len(dummy) == 4:
                self.temp.append(float(dummy[3]))
        self.t = np.array(self.t)
        self.ke = np.array(self.ke)

        if self.adim:
            self.friction = np.array(self.friction)*self.u_dim
        else:
            self.friction = np.array(self.friction)

        self.temp = np.array(self.temp)

    def plotConvergence(self):
        t = []
        ke = []
        friction = []
        temp = []

        in_file = open (self.path+'/output/check_ke.dat', 'rt')
        for line in in_file:
            dummy = line.split()
            t.append(float(dummy[0]))
            ke.append(float(dummy[1]))
            friction.append(float(dummy[2]))
            if len(dummy)==4:
                temp.append(float(dummy[3]))
        plt.figure(figsize=(10,4))
        plt.subplot(131)
        plt.plot(t,ke)
        plt.ylabel('ke')
        plt.xlabel('t')
        plt.subplot(132)
        plt.plot(t,friction)
        plt.ylabel('u*')
        plt.xlabel('t')
        plt.subplot(133)
        if len(dummy)==4:
            plt.plot(t,temp)
        plt.ylabel('Tsurf')
        plt.xlabel('t')
        plt.tight_layout()

        files_in_directory = os.listdir(self.path +'/output/')
        file_name_prefix = "angle_"
        matching_files = [file for file in files_in_directory if file.startswith(file_name_prefix)]
        if matching_files:
            iter = []
            angles = []
            for file in matching_files :
                iter.append(int(file[6:-3]))
                angles.append(list(h5py.File(self.path + '/output/' + file,'r')['var']))

            plt.figure()
            plt.plot(angles[-1])

            angles_array = np.array(angles)
            print(angles_array.shape)
        
            plt.figure()
            plt.plot(iter,angles_array[:,7])

            plt.figure()
            plt.plot(angles_array[-1,:])

    def takeProfile(self,x,y,ratio=1,epsilon=0.01):
        ix = np.argmin(np.abs(self.x-x))
        iy = np.argmin(np.abs(self.y-y))

        print('x=' + str(self.x[ix]) + ' instead of x_target=' + str(x))
        print('y=' + str(self.y[iy]) + ' instead of y_target=' + str(y))

        self.U_inf = self.u[:,iy,ix]*ratio
        self.z_coord = self.z

        if epsilon is None:
            print("took real eps  profile")
            self.epsilon_Kol = self.epsilon[:,iy,ix]
        else:
            print("set eps to 0.01")
            self.epsilon_Kol = epsilon*np.ones((len(self.z_coord)))


    def constant_epsilon(self,epsilon):
        self.epsilon_Kol = epsilon*np.ones((len(self.z_coord)))

    # def dissipation(self):
    #     #Read all the velocities and temperature
    #     filnam = self.path + '/output/tavg_txxs11.h5'
    #     if not os.path.isfile(filnam):
    #         print('no epsilon computed')
    #         return

    #     vel = h5py.File(filnam,"r")
    #     txxs11 = np.array(vel['var'])
    #     txxs11 = txxs11[0:self.nz,:,:]
    #     txxs11 = np.asarray(txxs11)
    #     txxs11[1:self.nz-1,:,:] = 0.5 * (txxs11[0:self.nz-2,:,:] + txxs11[1:self.nz-1,:,:])
    #     txxs11[0,:,:] = 0.
    #     #Read all the velocities and temperature 
    #     filnam = self.path + '/output/tavg_txys12.h5'
    #     vel = h5py.File(filnam,"r")
    #     txys12 = np.array(vel['var'])
    #     txys12 = txys12[0:self.nz,:,:]
    #     txys12 = np.asarray(txys12)
    #     txys12[1:self.nz-1,:,:] = 0.5 * (txys12[0:self.nz-2,:,:] + txys12[1:self.nz-1,:,:])
    #     txys12[0,:,:] = 0.
    #     # txys12 = txys12
    #     #Read all the velocities and temperature 
    #     filnam = self.path + '/output/tavg_txzs13.h5'
    #     vel = h5py.File(filnam,"r")
    #     txzs13 = np.array(vel['var'])
    #     txzs13 = txzs13[0:self.nz,:,:]
    #     txzs13 = np.asarray(txzs13)
    #     # txzs13 = txzs13
    #     #Read all the velocities and temperature 
    #     filnam = self.path + '/output/tavg_tyys22.h5'
    #     vel = h5py.File(filnam,"r")
    #     tyys22  = np.array(vel['var'])
    #     tyys22  = tyys22[0:self.nz,:,:]
    #     tyys22  = np.asarray(tyys22)
    #     tyys22[1:self.nz-1,:,:] = 0.5 * (tyys22[0:self.nz-2,:,:] + tyys22[1:self.nz-1,:,:])
    #     tyys22[0,:,:] = 0.
    #     # tyys22 = tyys22
    #     #Read all the velocities and temperature 
    #     filnam = self.path + '/output/tavg_tyzs23.h5'
    #     vel = h5py.File(filnam,"r")
    #     dummy = np.array(vel['var'])
    #     dummy = dummy[0:self.nz,:,:]
    #     tyzs23 = np.asarray(dummy)
    #     # tyzs23 = tyzs23
    #     #Read all the velocities and temperature 
    #     filnam = self.path + '/output/tavg_tzzs33.h5'
    #     vel = h5py.File(filnam,"r")
    #     tzzs33 = np.array(vel['var'])
    #     tzzs33 = tzzs33[0:self.nz,:,:]
    #     tzzs33 = np.asarray(tzzs33)
    #     # tzzs33 = tzzs33

    #     # mu = 15e-6
    #     # rho = 1.293 
    #     # self.epsilon = -(mu/rho)*(txxs11 + tyys22 +  tzzs33 + 2.0 * txzs13 + 2.0 * tyzs23 + 2.0 * txys12)*1000
    #     #self.epsilon = - self.u_dim**3 / self.z_i * (txxs11 + tyys22 +  tzzs33 + 2.0 * txzs13 + 2.0 * tyzs23 + 2.0 * txys12)
    #     self.epsilon = - self.u_dim**3 / self.z_i * (txxs11 + tyys22 +  tzzs33 + 2.0 * txzs13 + 2.0 * tyzs23 + 2.0 * txys12)
    #     # if self.adim :
    #     #     self.epsilon = - 1 / self.u_star**3 / self.z_i * (txxs11 + tyys22 +  tzzs33 + 2.0 * txzs13 + 2.0 * tyzs23 + 2.0 * txys12)
    #     # else: 
    #     #     self.epsilon = - 1/ self.z_i * (txxs11 + tyys22 +  tzzs33 + 2.0 * txzs13 + 2.0 * tyzs23 + 2.0 * txys12)
    #     #self.epsilon = - 1 / self.z_i * (txxs11 + tyys22 +  tzzs33 + 2.0 * txzs13 + 2.0 * tyzs23 + 2.0 * txys12)
    # 

    def dissipation(self):
        # Read all the velocities and temperature
        filnam = self.path + '/output/tavg_txxs11.h5'
        if not os.path.isfile(filnam):
            print('no epsilon computed')
            return
        vel = h5py.File(filnam,"r")
        txxs11 = np.array(vel['var'])
        txxs11 = txxs11[0:self.nz,:,:]
        txxs11 = np.asarray(txxs11)
        
        filnam = self.path + '/output/tavg_txys12.h5'
        vel = h5py.File(filnam,"r")
        txys12 = np.array(vel['var'])
        txys12 = txys12[0:self.nz, :, :]
        txys12 = np.asarray(txys12)

        filnam = self.path + '/output/tavg_txzs13.h5'
        vel = h5py.File(filnam, "r")
        txzs13 = np.array(vel['var'])
        txzs13 = txzs13[0:self.nz, :, :]
        txzs13 = np.asarray(txzs13)
        txzs13[0:self.nz-1,:,:] = 0.5 * (txzs13[0:self.nz-1,:,:] + txzs13[1:self.nz,:,:])

        filnam = self.path + '/output/tavg_tyys22.h5'
        vel = h5py.File(filnam,"r")
        tyys22  = np.array(vel['var'])
        tyys22  = tyys22[0:self.nz,:,:]
        tyys22  = np.asarray(tyys22)

        filnam = self.path + '/output/tavg_tyzs23.h5'
        vel = h5py.File(filnam,"r")
        dummy = np.array(vel['var'])
        dummy = dummy[0:self.nz,:,:]
        tyzs23 = np.asarray(dummy)
        tyzs23[0:self.nz-1,:,:] = 0.5 * (tyzs23[0:self.nz-1,:,:] + tyzs23[1:self.nz,:,:])

        filnam = self.path + '/output/tavg_tzzs33.h5'
        vel = h5py.File(filnam,"r")
        tzzs33 = np.array(vel['var'])
        tzzs33 = tzzs33[0:self.nz,:,:]
        tzzs33[0:self.nz-1,:,:] = 0.5 * (tzzs33[0:self.nz-1,:,:] + tzzs33[1:self.nz,:,:])
        tzzs33 = np.asarray(tzzs33)

        # mu = 15e-6
        # rho = 1.293 
        # self.epsilon = -(mu/rho)*(txxs11 + tyys22 +  tzzs33 + 2.0 * txzs13 + 2.0 * tyzs23 + 2.0 * txys12)*1000
        #self.epsilon = - self.u_dim**3 / self.z_i * (txxs11 + tyys22 +  tzzs33 + 2.0 * txzs13 + 2.0 * tyzs23 + 2.0 * txys12)
        self.epsilon = - self.u_dim**3 / self.z_i * (txxs11 + tyys22 +  tzzs33 + 2.0 * txzs13 + 2.0 * tyzs23 + 2.0 * txys12)
        # if self.adim :
        #     self.epsilon = - 1 / self.u_star**3 / self.z_i * (txxs11 + tyys22 +  tzzs33 + 2.0 * txzs13 + 2.0 * tyzs23 + 2.0 * txys12)
        # else: 
        #     self.epsilon = - 1/ self.z_i * (txxs11 + tyys22 +  tzzs33 + 2.0 * txzs13 + 2.0 * tyzs23 + 2.0 * txys12)
        #self.epsilon = - 1 / self.z_i * (txxs11 + tyys22 +  tzzs33 + 2.0 * txzs13 + 2.0 * tyzs23 + 2.0 * txys12)
        # self.epsilon = np.insert(self.epsilon, 0, 0, axis=0)
    
    def turbulent_intensity(self):
        filnam = self.path + '/output/tavg_u2.h5'
        vel = h5py.File(filnam, "r")
        u2 = np.array(vel['var'])
        #u2[1:self.nz, :, :] = 0.5 * (u2[0:self.nz-1, :, :] + u2[1:self.nz, :, :])
        # u2[0, :, :] = 0.

        filnam = self.path + '/output/tavg_v2.h5'
        vel = h5py.File(filnam, "r")
        v2 = np.array(vel['var'])
        # v2[1:self.nz, :, :] = 0.5 * (v2[0:self.nz-1, :, :] + v2[1:self.nz, :, :])
        # v2[0, :, :] = 0.

        filnam = self.path + '/output/tavg_w2.h5'
        vel = h5py.File(filnam, "r")
        w2 = np.array(vel['var'])
        w2[0:self.nz-1,:,:] = 0.5 * (w2[0:self.nz-1,:,:] + w2[1:self.nz,:,:])
        #w2[0,:,:] = 0.

        # u2 = u2*self.u_dim**2
        # v2 = v2*self.u_dim**2
        # w2 = w2*self.u_dim**2
        filnam = self.path + '/output/tavg_u.h5'
        vel = h5py.File(filnam, "r")
        u = np.array(vel['var'])
        filnam = self.path + '/output/tavg_v.h5'
        vel = h5py.File(filnam, "r")
        v = np.array(vel['var'])
        filnam = self.path + '/output/tavg_w.h5'
        vel = h5py.File(filnam, "r")
        w = np.array(vel['var'])


        # u = np.zeros((self.nz, self.ny, self.nx))
        # v = np.zeros((self.nz, self.ny, self.nx))
        # w = self.w

        #u[1:self.nz,:,:] = 0.5 * (u[0:self.nz-1,:,:] + u[1:self.nz,:,:])
        #v[1:self.nz,:,:] = 0.5 * (v[0:self.nz-1,:,:] + v[1:self.nz,:,:])
        w[0:self.nz-1,:,:] = 0.5 * (self.w[0:self.nz-1,:,:] + self.w[1:self.nz,:,:])

        TKE = 0.5*(u2 + v2 + w2 - (u**2 + v**2 + w**2))
        self.ti =  np.sqrt(2.*TKE/3) / np.sqrt(u**2 + v**2 + w**2)
        # print(self.ti.shape)
        # print(np.mean(self.ti,(1,2)))

    def plotTurbineAngle(self):
        i = 0
        lines = []
        data = []
        in_file = open (self.path+'turbine/turbine0001.dat', 'rt')
        for line in in_file:
            lines.append(line)
            dummy = lines[i].split()
            data.append(dummy)
            i+=1
        self.data = np.array(data).astype(np.float)
        
    def plotXZ(self,y,u=True,eps=False,**kwargs):
        iy = np.argmin(np.abs(self.y-y))
        if u:
            cax = plt.pcolormesh(self.x,self.z,self.u[:,iy,:],shading='gouraud',**kwargs)
            plt.gca().set_aspect('equal', adjustable='box')
            plt.xlabel('x (m)')
            plt.ylabel('z (m)')
            ax=plt.gca()
            plt.rcParams["contour.negative_linestyle"] = "solid"
            contours = ax.contour(self.x,self.z,self.u[:,iy,:], [0,2,4,6,8,10,12,14,16],colors="b",
        linewidths=0.8,
        alpha=0.3,)
            ax.clabel(contours, contours.levels, inline=True, fontsize=10, colors="black")
            divider = make_axes_locatable(ax)
            cb = divider.append_axes("right", size="5%", pad=0.05)
            cb.set_title('u (m/s)')
            plt.colorbar(cax, cax=cb)
            plt.tight_layout()
        if eps:
            cax = plt.pcolormesh(self.x,self.z,self.epsilon[:,iy,:],shading='gouraud',**kwargs)
            plt.gca().set_aspect('equal', adjustable='box')
            plt.xlabel('x (m)')
            plt.ylabel('z (m)')
            ax=plt.gca()
            divider = make_axes_locatable(ax)
            cb = divider.append_axes("right", size="5%", pad=0.05)
            cb.set_title('eps (m2/s3)')
            plt.colorbar(cax, cax=cb)
            plt.tight_layout()
        return cax,ax

    def plotXY(self,z,u=True,eps=False,**kwargs):
        iz = np.argmin(np.abs(self.z-z))
        
        if u:
            cax = plt.pcolormesh(self.x,self.y,self.u[iz,:,:],shading='gouraud',**kwargs)
            plt.gca().set_aspect('equal', adjustable='box')
            plt.xlabel('x (m)')
            plt.ylabel('y (m)')
            ax=plt.gca()
            plt.rcParams["contour.negative_linestyle"] = "solid"
            ax.contour(self.x,self.y,self.u[iz,:,:], [0,2,4,6,8,10,12,14,16],colors="b",
        linewidths=0.8,
        alpha=0.3,)
            divider = make_axes_locatable(ax)
            cb = divider.append_axes("right", size="5%", pad=0.05)
            cb.set_title('u (m/s)')
            plt.colorbar(cax, cax=cb)
            plt.tight_layout()
        if eps:
            cax = plt.pcolormesh(self.x,self.y,self.epsilon[iz,:,:],shading='gouraud',**kwargs)
            plt.gca().set_aspect('equal', adjustable='box')
            plt.xlabel('x (m)')
            plt.ylabel('y (m)')
            ax=plt.gca()
            divider = make_axes_locatable(ax)
            cb = divider.append_axes("right", size="5%", pad=0.05)
            cb.set_title('$\\epsilon$ (m2/s3)')
            plt.colorbar(cax, cax=cb)
            plt.tight_layout()
        return cax,ax
           
    def plotYZ(self,x,**kwargs):
        ix = np.argmin(np.abs(self.x-x))
        plt.pcolormesh(self.y,self.z,self.u[:,:,ix],cmap='RdBu_r',shading='gouraud',**kwargs)
        plt.colorbar()
        plt.gca().set_aspect('equal', adjustable='box')
