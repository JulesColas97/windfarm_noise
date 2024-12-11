import numpy as np
import os
from os import path
from scipy.io import loadmat
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import scipy.interpolate as interp
import pickle
import logging
from matplotlib.widgets import Slider
import time 

from .WP_LEE import WP_LEE
from .func_repo_XFOIL import XFOIL, run_Xvfb,kill_Xvfb
from .utils import computeThirdOctaveFrequencies, interp_atmos_data
from . import BEM as bem
 

class WindTurbine:
    """A class to create a wind turbine, with geometry and boundary layer quantities 
    """
    href : float                                                                # hub height
    Nseg : int                                                                  # number of segments per blade
    Nblade : int                                                                # number of blades
    Nbeta : int                                                                 # number of angular position
    omega : float                                                               # rotation speed
    delta_beta : float                                                          # angular position
    seg : np.ndarray                                                            # segment position
    chord : np.ndarray                                                          # segment chord length
    Lspan : np.ndarray                                                          # sgement span length
    twist : np.ndarray                                                          # segement twist
    airfoilIndex : np.ndarray                                                   # segment airfoil type
    pathToAirfoilGeom : str
    airfoilDataName : np.ndarray                                                # xfoil data path
    airfoil_data : dict                                                         # xfoil data
    tau: float                                                                  # orientation of the turbine (from x axis CCW)
    absolute_pos: tuple                                                         # position in LES grid
    def __init__(self):
        return


    def setOptimalTwist(self,U_ref:np.array,AoA_opt:float):
        """Create optimal twist for a given angle of attack (used to create Cotte wind turbines)

        Args:
            U_ref (np.ndarray): wind profile
            AoA_opt (float): optimal angle of attack wanted 
        """
        # twist and AoA of the wind turbine blades
        U_rot = self.seg*self.omega
        self.twist = np.arctan(U_rot / U_ref) + AoA_opt * np.pi / 180


    def default(self):
        """Create default wind turbine 
        same as Cotte but scaled up for D = 100m
        """
        self.tau = 0
        self.absolute_pos = (0.,0.)
        self.href = 100
        self.Nseg = 8
        self.Nblade = 3
        self.Nbeta = 12
        self.omega =  12.1 * 2 * np.pi / 60
        self.delta_beta = 0
        self.seg = 1.07*np.asarray([6.0454, 14.8958, 22.7309, 29.1557, 34.4240, 38.7440, 42.2865, 45.1913])  #Spnwise location of segment center point
        self.Lspan = 1.07* np.asarray([9.0908, 8.6100, 7.0602, 5.7894, 4.7473, 3.8928, 3.1921, 2.6175])       #Spanwise lengths of blade segments
        self.chord = np.asarray([3.2273, 2.6963, 2.2261, 1.8407, 1.5246, 1.2654, 1.0528, 0.8785])
        self.setOptimalTwist(8,4)
        self.airfoilIndex = np.asarray([0,0,0,0,0,0,0,0])
        path2data =  path.join(path.dirname(__file__), 'BL_data/NACA63415')
        self.airfoilDataPath = np.asarray([path2data])
        self.airfoil_data = {}
        self.airfoil_data = pickle.load(open(path2data,'rb'))
        # self.airfoil_data = loadmat(pathToAirfoilData+'NACA63415')
    
    def computeBeta(self):
        """compute beta angle according to number of blade and number of discretization per 1/Nblade turn 
        the beta array is of shape (Nblade, Nbeta). The total number of angular position is Nblade*Nbeta and cover a complete 360 rotation.
        """
        dbeta = 2 * np.pi / self.Nblade
        self.beta = (np.linspace(0, 2 * np.pi - dbeta, self.Nblade).reshape(self.Nblade, 1) + np.linspace(0,2*np.pi*(self.Nbeta-1)/self.Nblade/self.Nbeta,self.Nbeta).reshape(1,self.Nbeta)) #+ np.pi/3


    def proximity(self,height_array:np.ndarray)->np.ndarray:
        """Compute the proximity matrix to combine delta L and Spp. 

        Args:
            height_array (np.ndarray): the source heights array corresponding to the Delta L computations 

        Returns:
            np.ndarray: array of shape (Nseg*Nblade,Nbeta) of the index corsponding to the closest Delta L source 
        """

        proximity = np.zeros((self.Nseg*self.Nblade, self.Nbeta))
        self.computeBeta()
        # print(self.beta.shape)
        for ibeta in range(self.Nbeta):
            h_temp = (np.cos(self.beta[:,ibeta]).reshape(1,self.Nblade) * self.seg.reshape(self.Nseg,1) + self.href).reshape(self.Nblade * self.Nseg)
            # print(h_temp.shape)
            # print( np.abs(h_temp - height_array).argmin(axis = 0))
            proximity[:, ibeta] = np.abs(h_temp.reshape(1,-1) - height_array.reshape(-1,1)).argmin(axis = 0)
        return proximity.astype(int)


    def proximityLinear(self,height_array: np.ndarray)->np.ndarray:
        """compute th proximity matrix to combine delta L and Spp, with linear interpolation between the Delta L

        Args:
            height_array (np.ndarray): the source heights array corresponding to the Delta L computations 

        Returns:
            np.ndarray: array of shape (Nseg*Nblade,Nbeta,3) of the two index corresponding to the closest DeltaL computations, and the ratio between them. 
        """
        # height_array is the heigth of the delta L computations 
        # third dimension to store the  2 delta L index for each source and the ratio between the two delta L 
        proximity = np.zeros((self.Nseg*self.Nblade, self.Nbeta,3),dtype=object)
        # compute beta matrix 
        self.computeBeta()
        # loop on the blades position
        for ibeta in range(self.Nbeta):
            # calculate heigth of each segment at each blade position 
            h_temp = (np.cos(self.beta[:,ibeta]).reshape(1,self.Nblade) * self.seg.reshape(self.Nseg,1) + self.href).reshape(self.Nblade * self.Nseg)
            for iseg in range (self.Nseg*self.Nblade) : 
                prox0 = np.abs(h_temp[iseg] - height_array).argmin()
                if h_temp[iseg]<height_array[prox0] and prox0>0:
                    prox1 = prox0-1 
                    prox2 = prox0   
                    ratio = (height_array[prox2] - h_temp[iseg])/(height_array[prox2]-height_array[prox1])
                elif h_temp[iseg]>height_array[prox0] and prox0<len(height_array)-1:
                    prox1 = prox0
                    prox2 = prox0+1
                    ratio = (height_array[prox2] - h_temp[iseg])/(height_array[prox2]-height_array[prox1])
                else : 
                    prox1 = prox0
                    prox2 = prox0
                    ratio = 1
                proximity[iseg,ibeta,0] = int(prox1)
                proximity[iseg,ibeta,1] = int(prox2)
                proximity[iseg,ibeta,2] = float(ratio) 
        return proximity


    def plotProximity(self,height_array:np.ndarray,ax1=None,ibeta :int=None):
        """plot the corespondance between segement and Delta L (using colors)

        Args:
            height_array (np.ndarray): the source heights array corresponding to the Delta L computations 
            ax1 (_type_, optional): ax object on which to plot. Defaults to None.
            ibeta (int, optional): index of the angular position to plot, if None cycle through all beta. Defaults to None.
        """
        proximity = self.proximity(height_array)
        Ncolor = np.max(proximity)+1
        icol = np.linspace(0,1,Ncolor)
        cmap = cm.get_cmap('coolwarm')
        color = []

        for ii in icol:
            color.append(cmap(ii))

        beta = self.beta.reshape(self.Nblade,self.Nbeta)
        if ibeta is not None:
            # fig, ax = plt.subplots()
            for i, height in  enumerate(height_array):
                if ax1 is None : 
                    plt.plot([-100,100], [height,height],'--', color= color[i])
                else : 
                    ax1.plot([-100,100], [height,height],'--', color= color[i])

            proximity_temp = proximity[:,ibeta].reshape(self.Nseg, self.Nblade)
            for iblade in range(self.Nblade):
                for iseg in range(self.Nseg) :
                    # TODO probelm with angle 
                    # x = np.sin(beta[iblade,ibeta]) * self.seg[iseg] + np.sin(beta[iblade,ibeta])*0.5*self.Lspan[iseg]  - np.cos(beta[iblade,ibeta])*0.5*self.chord[iseg]
                    # y = self.href + np.cos(beta[iblade,ibeta]) * self.seg[iseg] - np.cos(beta[iblade,ibeta])*0.5*self.Lspan[iseg]  - np.sin(beta[iblade,ibeta])*0.5*self.chord[iseg]
                    # height = self.Lspan[iseg]
                    # width = self.chord[iseg]
                    # angle = beta[iblade,ibeta]
                    # if ax1 is None : 
                    #     plt.gca().add_patch(Rectangle((x, y), width,height, -angle*180/np.pi,
                    #     facecolor =color[proximity_temp[iseg,iblade]],
                    #     fill=True))
                    # else:
                    #     ax1.add_patch(Rectangle((x, y), width,height, -angle*180/np.pi,
                    #     facecolor =color[proximity_temp[iseg,iblade]],
                    #     fill=True))

                    #segment 
                    if iseg == 0 :
                            x0 = 0
                            z0 = self.href
                    else :
                        x0 =np.sin(beta[iblade,ibeta]) * self.seg[iseg-1]
                        z0 = np.cos(beta[iblade,ibeta]) * self.seg[iseg-1] + self.href
                    x1 =np.sin(beta[iblade,ibeta]) * self.seg[iseg]
                    z1 = np.cos(beta[iblade,ibeta]) * self.seg[iseg] + self.href
                    if ax1 is None : 
                        plt.plot([x0,x1],[z0,z1],'-',linewidth = 5, color = color[proximity_temp[iseg,iblade]])
                    else :
                        ax1.plot([x0,x1],[z0,z1],'-',linewidth = 5, color = color[proximity_temp[iseg,iblade]])


            if ax1 is None :         
                plt.gca().set_aspect('equal', adjustable='box')
                plt.ylim(self.href -50,self.href + 50)
                plt.xlim(-50,50) 
            else :
                ax1.set_aspect('equal', adjustable='box')   
                ax1.set_ylim(self.href -50,self.href + 50)
                ax1.set_xlim(-50,50)

        else :
            for ibeta in range(self.Nbeta):
                proximity_temp = proximity[:,ibeta].reshape(self.Nseg, self.Nblade)
                plt.figure(ibeta)
                for iblade in range(self.Nblade):
                    for iseg in range(self.Nseg) :
                        if iseg == 0 :
                            x0 = 0
                            z0 = self.href
                        else :
                            x0 =np.sin(beta[iblade,ibeta]) * self.seg[iseg-1]
                            z0 = np.cos(beta[iblade,ibeta]) * self.seg[iseg-1] + self.href
                        x1 =np.sin(beta[iblade,ibeta]) * self.seg[iseg]
                        z1 = np.cos(beta[iblade,ibeta]) * self.seg[iseg] + self.href


                        plt.plot([x0,x1],[z0,z1],'-',linewidth = 5, color = color[proximity_temp[iseg,iblade]])
                plt.gca().set_aspect('equal', adjustable='box')
                plt.xlim(-50,50)
                plt.ylim(self.href -50,self.href + 50)
                plt.show()


    def createGeometry(self,seg:np.array,Lspan:np.array,c:np.array,twist:np.array,airfoilIndex:np.array, target_l_c=3.0, D=None):
        """Create the geometry for the source model from aero geometry of a blade. The final blade respect the given target_lc_c ratio

        Args:
            seg (np.array): array of the center position of the segments
            Lspan (np.array): array of teh span length of the segments
            c (np.array): array of the chord length of the segments
            twist (np.array): array of the twist of the segments
            airfoilIndex (np.array): index of the airfoil to be used (coresponding to self.airfoilDataPath )
            target_l_c (float, optional): span/chord target ratio (to respect the incoherent source assumption). Defaults to 3.0.
            D (float, optional): new Diameter to resize the wind turbine. Defaults to None.
        """
        
        # change diameter size 
        if D is not None : 
            originalDiameter = seg[-1] + 0.5*Lspan[-1]
            seg = (0.5*D/originalDiameter) * seg
            Lspan =  (0.5*D/originalDiameter) * Lspan
            c = (0.5*D/originalDiameter) * c
            logging.info('WT radius set to %s m' % (seg[-1] + 0.5*Lspan[-1]))
        if target_l_c is None : 
            self.seg = seg
            self.Lspan = Lspan 
            self.twist = twist
            self.chord = c
            self.airfoilIndex = airfoilIndex
            return 
        # for interpolation of chord and twist 
        tck_c = interp.splrep(seg, c, k = 2)
        tck_twist = interp.splrep(seg, twist, k = 2)

        # init list for new geometry 
        new_seg = []
        new_Lspan = []
        new_twist = []
        new_c = []
        new_airfoilIndex = []

        # init first segment 
        total_span = seg[0]- 0.5*Lspan[0]
        current_span = 0 
        current_seg = seg[0]- 0.5*Lspan[0]
        current_c = float(interp.splev(current_seg, tck_c, der = 0))
        delta_span = 0.1
        # increase length until the end of the blade 
        while total_span <= seg[-1]+0.5*Lspan[-1]:
            # increase length of current segment until target is reach 
            if current_span / current_c < target_l_c:
                current_span += 2*delta_span
                total_span += 2*delta_span
                current_seg += delta_span 
                current_c = float(interp.splev(current_seg, tck_c, der = 0))
            # interpolate and store new values for the segment 
            else:
                new_seg.append(current_seg)
                new_Lspan.append(current_span)
                new_c.append(current_c)
                new_twist.append(float(interp.splev(current_seg, tck_twist, der = 0)))
                print(np.where(current_seg-seg-0.5*Lspan<0)[0][0])
                new_airfoilIndex.append(airfoilIndex[np.where(current_seg-seg-0.5*Lspan<0)[0][0]])
                # new_airfoilIndex.append(airfoilIndex[(np.abs(seg-current_seg)).argmin()])

                # move the current segment position to the end of the segment 
                current_seg += 0.5*current_span
                # reset curent span to 0
                current_span = 0
                # set current to twist according to segment position 
                current_c = float(interp.splev(current_seg, tck_c, der = 0))

        # add something to reach final    D size 
        new_seg[-1] = new_seg[-1] + 0.5*current_span
        new_Lspan[-1] = new_Lspan[-1] + current_span  
        new_twist[-1] = float(interp.splev(new_seg[-1], tck_twist, der = 0))
        new_c[-1] = float(interp.splev(new_seg[-1], tck_c, der = 0))

        new_airfoilIndex[-1] = airfoilIndex[(np.abs(seg-current_seg)).argmin()]

        self.seg = np.array(new_seg)
        self.Lspan = np.array(new_Lspan)
        self.twist = np.array(new_twist)
        self.chord = np.array(new_c)
        self.airfoilIndex = np.array(new_airfoilIndex)
        return


    def createBLData(self,AoA:np.ndarray,reynolds:np.ndarray,fname:str):
        """create the BL quantites data base for the given blade, AoA and re

        Args:
            AoA (np.ndarray): Angle of attack array to compute for the DB
            reynolds (np.ndarray): Re numbers to compute for the DB
            fname (str): file name to save the BL quantities
        """

        # Init BL array
        delta_star_bot = np.zeros((reynolds.shape[0], AoA.shape[0], self.chord.shape[0]))
        theta_momen_bot = np.zeros((reynolds.shape[0], AoA.shape[0], self.chord.shape[0]))
        cf_bot = np.zeros((reynolds.shape[0], AoA.shape[0], self.chord.shape[0]))
        dpdx_bot = np.zeros((reynolds.shape[0], AoA.shape[0], self.chord.shape[0]))
        UeUinf_bot = np.zeros((reynolds.shape[0], AoA.shape[0], self.chord.shape[0]))

        delta_star_top = np.zeros((reynolds.shape[0], AoA.shape[0], self.chord.shape[0]))
        theta_momen_top = np.zeros((reynolds.shape[0], AoA.shape[0], self.chord.shape[0]))
        cf_top = np.zeros((reynolds.shape[0], AoA.shape[0], self.chord.shape[0]))
        dpdx_top = np.zeros((reynolds.shape[0], AoA.shape[0], self.chord.shape[0]))
        UeUinf_top = np.zeros((reynolds.shape[0], AoA.shape[0], self.chord.shape[0]))

        foil = {}
        # airfoil geometry
        
        # no artificial  trip
        foil['trip_top'] = 1.0
        foil['trip_bottom'] = 1.0

        # percent of chord length
        # should be 95 :'(  but does the flow detached 
        target_frac = 0.85

        # out_dict = XFOIL(path, target_frac, inf, foil)
        # run virtual buffer server before xfoil calculation
        run_Xvfb()

        inf = {}
        inf['rho'] = 1.225
        inf['nu'] = 1.48e-5
        inf['T'] = 300
        # loop over the parameters of the dB
        for temp3 in range(self.chord.shape[0]):
            for temp2 in range(AoA.shape[0]):
                for temp1 in range(reynolds.shape[0]):
                    print(self.airfoilIndex[temp3])
                    print(self.airfoilDataName)
                    print(self.airfoilDataName[self.airfoilIndex[temp3]])
                    # create the input dictionnary for Xfoil
                    foil['name'] = self.pathToAirfoilGeom + self.airfoilDataName[self.airfoilIndex[temp3]] + '.dat'
                    foil['chord'] = self.chord[temp3]
                    foil['span'] = self.Lspan[temp3]
                    foil['AoA'] = AoA[temp2]
                    foil['re'] = reynolds[temp1]
                    inf_velocity = foil['re'] * inf['nu'] / foil['chord']
                    foil['mach'] = inf_velocity / (1.4 * 287 * inf['T']) ** 0.5
                    inf['mach'] = foil['mach']
                    # run Xfoil to obtain the BL quantities 
                    out_dict = XFOIL('./metadata/' , target_frac, inf, foil)

                    delta_star_top[temp1, temp2, temp3] = out_dict['Dstar_suction']
                    theta_momen_top[temp1, temp2, temp3] = out_dict['theta_suction']
                    cf_top[temp1, temp2, temp3] = out_dict['Cf_suction']
                    dpdx_top[temp1, temp2, temp3] = out_dict['dpdx_suction']
                    UeUinf_top[temp1, temp2, temp3] = out_dict['Ue_suction']

                    delta_star_bot[temp1, temp2, temp3] = out_dict['Dstar_pressure']
                    theta_momen_bot[temp1, temp2, temp3] = out_dict['theta_pressure']
                    cf_bot[temp1, temp2, temp3] = out_dict['Cf_pressure']
                    dpdx_bot[temp1, temp2, temp3] = out_dict['dpdx_pressure']
                    UeUinf_bot[temp1, temp2, temp3] = out_dict['Ue_pressure']
            file_name = './metadata/' + str(temp3) + '.dat'
            outfile = open(file_name, mode = 'w')
            outfile.write('Data generation has started!!!\n')
            outfile.write(str(self.chord[temp3]) + 'm has been completed')
            outfile.close()

        # kill the virtual frame buffer
        kill_Xvfb()

        # save the data base
        print("saving the data Base")
        print(AoA.shape)
        print(reynolds.shape)
        print(delta_star_bot.shape)

        airfoil_data = {}
        airfoil_data['AoA'] = AoA
        airfoil_data['reynolds'] = reynolds
        airfoil_data['chord'] = self.chord 
        airfoil_data['Lspan'] = self.Lspan 
        airfoil_data['seg'] = self.seg 
        airfoil_data['twist'] = self.twist 
        airfoil_data['airfoildIndex'] = self.airfoilIndex
        airfoil_data['airfoilDataName'] = self.airfoilDataName


        airfoil_data['delta_star_bot'] = delta_star_bot
        airfoil_data['delta_star_top'] = delta_star_top
        airfoil_data['cf_bot'] = cf_bot
        airfoil_data['cf_top'] = cf_top
        airfoil_data['dpdx_bot'] = dpdx_bot
        airfoil_data['dpdx_top'] = dpdx_top
        airfoil_data['theta_momen_bot'] = theta_momen_bot
        airfoil_data['theta_momen_top'] = theta_momen_top
        airfoil_data['UeUinf_bot'] = UeUinf_bot
        airfoil_data['UeUinf_top'] = UeUinf_top

        pickle.dump(airfoil_data,  open(fname, 'wb'))


    def createBLData2(self,AoA:np.ndarray,reynolds:np.ndarray,fname:str):
        """create the BL quantites data base for the given blade, AoA and re
        also compute Cl and Cd (for BEM implementation)
        NOT CHECKED YET 

        Args:
            AoA (np.ndarray): Angle of attack array to compute for the DB
            reynolds (np.ndarray): Re numbers to compute for the DB
            fname (str): file name to save the BL quantities
        """
        delta_star_bot = np.zeros((reynolds.shape[0], AoA.shape[0], self.chord.shape[0]))
        theta_momen_bot = np.zeros((reynolds.shape[0], AoA.shape[0], self.chord.shape[0]))
        cf_bot = np.zeros((reynolds.shape[0], AoA.shape[0], self.chord.shape[0]))
        dpdx_bot = np.zeros((reynolds.shape[0], AoA.shape[0], self.chord.shape[0]))
        UeUinf_bot = np.zeros((reynolds.shape[0], AoA.shape[0], self.chord.shape[0]))

        delta_star_top = np.zeros((reynolds.shape[0], AoA.shape[0], self.chord.shape[0]))
        theta_momen_top = np.zeros((reynolds.shape[0], AoA.shape[0], self.chord.shape[0]))
        cf_top = np.zeros((reynolds.shape[0], AoA.shape[0], self.chord.shape[0]))
        dpdx_top = np.zeros((reynolds.shape[0], AoA.shape[0], self.chord.shape[0]))
        UeUinf_top = np.zeros((reynolds.shape[0], AoA.shape[0], self.chord.shape[0]))
        Cl = np.zeros((reynolds.shape[0], AoA.shape[0], self.chord.shape[0]))
        Cd = np.zeros((reynolds.shape[0], AoA.shape[0], self.chord.shape[0]))
        Cl_tck = []
        Cd_tck = []
        foil = {}
        # no artificial  trip
        foil['trip_top'] = 0.05 #1.0
        foil['trip_bottom'] = 0.05 #1.0
        # percent of chord length
        target_frac = 0.95

        # run virtual buffer server before xfoil calculation
        run_Xvfb()

        inf = {}
        inf['rho'] = 1.2
        inf['nu'] = 1.48e-5
        # inf['nu'] = 1.62e-5
        inf['T'] = 273+20 
        # loop over the parameters of the dB (segments, reynolds, AoA) 
        for temp3 in range(self.chord.shape[0]):
            # change input 
            foil['name'] = self.pathToAirfoilGeom + self.airfoilDataName[self.airfoilIndex[temp3]] + '.dat'
            foil['chord'] = self.chord[temp3]
            foil['span'] = self.Lspan[temp3]
            for temp1 in range(reynolds.shape[0]):
                # foil['re'] = reynolds[temp1]
                # modif after B cotté reunion Re différent de Rc
                foil['re'] = reynolds[temp1]/foil['chord']
                inf_velocity = reynolds[temp1] * inf['nu'] / foil['chord']
                foil['mach'] = inf_velocity / (1.4 * 287 * inf['T']) ** 0.5
                inf['mach'] = foil['mach']

                for temp2 in range(AoA.shape[0]):
                    print(self.airfoilIndex[temp3])
                    print(self.airfoilDataName)
                    print(self.airfoilDataName[self.airfoilIndex[temp3]])
                    foil['AoA'] = AoA[temp2]
                    out_dict,polar = XFOIL('./metadata/' , target_frac, inf, foil)
                    # this is a hard coded shit in case XFOIl does converge for the given AoA 
                    # it whil decrease slightly until it converges 
                    cont = 0
                    print(polar.size)
                    
                    if polar.size ==0:
                        while polar.size == 0 and cont<10:
                            print(foil['AoA'])
                            foil['AoA'] = foil['AoA'] - 0.05
                            out_dict,polar = XFOIL('./metadata/' , target_frac, inf, foil)
                            cont+=1
                            
                        print(AoA[temp2])
                        print(cont)                        
                        print(foil['AoA'])
                        time.sleep(5)
                    if polar.size == 0:
                        print("seg No %s"%temp3)
                        print("Re = %s"%reynolds[temp1])
                        print("AoA = %s"%AoA[temp2])
                        print("error to big cant converge problem ")
                        quit()
                    delta_star_top[temp1, temp2, temp3] = out_dict['Dstar_suction']
                    theta_momen_top[temp1, temp2, temp3] = out_dict['theta_suction']
                    cf_top[temp1, temp2, temp3] = out_dict['Cf_suction']
                    dpdx_top[temp1, temp2, temp3] = out_dict['dpdx_suction']
                    UeUinf_top[temp1, temp2, temp3] = out_dict['Ue_suction']

                    delta_star_bot[temp1, temp2, temp3] = out_dict['Dstar_pressure']
                    theta_momen_bot[temp1, temp2, temp3] = out_dict['theta_pressure']
                    cf_bot[temp1, temp2, temp3] = out_dict['Cf_pressure']
                    dpdx_bot[temp1, temp2, temp3] = out_dict['dpdx_pressure']
                    UeUinf_bot[temp1, temp2, temp3] = out_dict['Ue_pressure'] 

                    Cl[temp1,temp2,temp3] = polar[1]
                    Cd[temp1,temp2,temp3] = polar[2]
            print(polar[1])
            print(polar[2])
            print(Cl.shape)
            print(Cd.shape)


            # time.sleep(5.5) 
            R,A = np.meshgrid(reynolds,AoA,indexing='ij')
            # print(R)
            # print(A)

            # print(Cd)
            Cl_tck.append(interp.bisplrep(R,A,Cl[:,:,temp3]))
            Cd_tck.append(interp.bisplrep(R,A,Cd[:,:,temp3]))


            file_name = './metadata/' + str(temp3) + '.dat'
            outfile = open(file_name, mode = 'w')
            outfile.write('Data generation has started!!!\n')
            outfile.write(str(self.chord[temp3]) + 'm has been completed')
            outfile.close()
        
        # kill the virtual frame buffer
        kill_Xvfb()

        # save the data base
        print("saving the data Base")
        print(AoA.shape)
        print(reynolds.shape)
        print(delta_star_bot.shape)


        self.airfoil_data = {}
        self.airfoil_data['AoA'] = AoA
        self.airfoil_data['reynolds'] = reynolds
        self.airfoil_data['chord'] = self.chord 
        self.airfoil_data['Lspan'] = self.Lspan 
        self.airfoil_data['seg'] = self.seg 
        self.airfoil_data['twist'] = self.twist 
        self.airfoil_data['airfoildIndex'] = self.airfoilIndex
        self.airfoil_data['airfoilDataName'] = self.airfoilDataName


        self.airfoil_data['delta_star_bot'] = delta_star_bot
        self.airfoil_data['delta_star_top'] = delta_star_top
        self.airfoil_data['cf_bot'] = cf_bot
        self.airfoil_data['cf_top'] = cf_top
        self.airfoil_data['Cl'] = Cl
        self.airfoil_data['Cd'] = Cd
        self.airfoil_data['Cl_tck'] = Cl_tck
        self.airfoil_data['Cd_tck'] = Cd_tck
        self.airfoil_data['dpdx_bot'] = dpdx_bot
        self.airfoil_data['dpdx_top'] = dpdx_top
        self.airfoil_data['theta_momen_bot'] = theta_momen_bot
        self.airfoil_data['theta_momen_top'] = theta_momen_top
        self.airfoil_data['UeUinf_bot'] = UeUinf_bot
        self.airfoil_data['UeUinf_top'] = UeUinf_top

        pickle.dump(self.airfoil_data,  open(fname, 'wb'))


    def readBladeData(self,fname:str):
        """read the BL data from file, set the seg, twist, chord and span accordingly

        Args:
            fname (str): path to the BL file
        """
        self.airfoil_data = pickle.load(open(fname,'rb'))
        self.chord = self.airfoil_data['chord'] 
        self.Lspan = self.airfoil_data['Lspan'] 
        self.seg = self.airfoil_data['seg']
        self.twist = self.airfoil_data['twist'] 
        self.airfoilIndex = self.airfoil_data['airfoildIndex'] 
        self.airfoilDataname = self.airfoil_data['airfoilDataName'] 
        self.Nseg = len(self.chord)

    # set omega according to Uref (from NREL documentation)
    def controlRotSpeed(self,Uref:float,omega:float=None):
        """set omega according to Uref from NREL documentation

        Args:
            Uref (float): reference wind speed (usually taken at hub height)
            omega (float, optional): if set, forces omega to the given value. Defaults to None.
        """
        # from Nrel decsription 
        cutin=3; rated=11.4; cutOut=25
        cutin_ratio = 6.9; rated_ratio = 12.1 

        U =  np.array([3.0760869565217392,
        4.315217391304348,
            5.016304347826087,
        6.010869565217392,
        7.478260869565218,
        8.994565217391305,
        9.711956521739129,
        10.413043478260871,
        11.032608695652172,
        11.4,
        12.23913043478261,
        13.021739130434781,
        14.815217391304348,
        16.934782608695652,
        19.46195652173913,
        24.95652173913043])

        omega_rpm = np.array([ 6.96319018404909,
        7.208588957055213,
        7.453987730061357,
        7.944785276073617,
        8.680981595092028,
        10.235173824130882,
        11.053169734151332,
        11.625766871165645,
        11.871165644171782,
        12.1,
        12.1,
        12.1,
        12.1,
        12.1,
        12.1,
        12.1])

        omega_tck = interp.splrep(U,omega_rpm,k=2)

        if omega is None:
            # self.omega = (omega_rpm[np.argmin(np.abs(Uref-U))]) * 2*np.pi/60
            self.omega = interp.splev(Uref, omega_tck) * 2*np.pi/60

            print('Omega set to :' +str( self.omega *60/(2*np.pi))+' rpm')
        else : 
           self.omega = omega

        # # linear increase of rotational speed with wind velocity 
        # if Uref > cutin and Uref < rated:
        #     ratio = (rated-Uref)/(rated-cutin)
        #    self.omega =  ((ratio)* cutin_ratio +(1-ratio)*rated_ratio) * 2 * np.pi / 60
        #     logging.info('Omega set to :' +str( self.omega *60/(2*np.pi))+' rpm')
        #elif Uref > rated and Uref < cutOut:
        #     self.omega = rated_ratio * 2 * np.pi / 60
        #    logging.info('Omega set to :' +str( self.omega *60/(2*np.pi))+' rpm')
        #else:
        #     self.omega = 0

    # set pitch of the blades according to Uref (from Nrel documentation)
    def controlPitchAngle(self,Uref:float,pitch:float=None):
        """set pitcch angle with respect to Uref according to Nrel documentation

        Args:
            Uref (float): reference wind speed (usually taken at hub height)
            pitch (float, optional): if set, forces th pitch angle to the given value. Defaults to None.
        """
        # from NREL calculation for a unifrom wind speed 
        # pitch_array = np.array([0,3.8,6.6,8.7,10.4,12,13.5,14.9,16.2,17.4,18.7,19.9,21.18,22.35,23.47])
        # U = np.array([11.4,12,13,14,15,16,17,18,19,20,21,22,23,24,25])

        U = np.array([3.0271739130434785,
            9.01086956521739,
            11.309782608695652,
            11.456521739130434,
            12.02717391304348,
            13.005434782608695,
            14,
            15.98913043478261,
            17,
            17.994565217391305,
            19.005434782608695,
            19.983695652173914,
            20.994565217391305,
            21.989130434782606,
            23.016304347826086,
            23.9945652173913,
            24.97282608695652])
        pitch_array = np.array([ 0,
                        0,
                        0,
                        0.9918200408998032,
                        3.7730061349693287,
                        6.635991820040907,
                        8.517382413087937,
                        12.034764826175874,
                        13.425357873210636,
                        14.897750511247445,
                        16.20654396728017,
                        17.43353783231084,
                        18.742331288343564,
                        19.969325153374236,
                        21.11451942740286,
                        22.259713701431494,
                        23.48670756646217])
        if pitch is None:
            alpha = 1
            self.twist = self.twist - alpha * np.pi/180 * (pitch_array[np.argmin(np.abs(Uref-U))])
            print(' wind turbine blade pitched : ' + str(pitch_array[np.argmin(np.abs(Uref-U))]) + 'deg')
        else : 
            self.twist = self.twist - np.pi/180 * pitch

    # find index in array (need for exploreDataBase)
    @staticmethod
    def find_index(value, array):
        for temp in range(len(array)):
            if value <= array[temp]:
                break
        return temp

    # calculate Wall pressure spectrum from LEE model (needed for explore Data base)
    @staticmethod
    def _WP_LEE(freq, delta_star, theta_momen, dpdx, cf, UeUinf, Uinf):
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

        return Phipp_LEE, Phipp_LEE_adim, Delta, beta_c


    def computeBEM(self,Nbeta,atmos,BEM=True):
        nu0 = 1.45e-5       

        Cl_tck = self.airfoil_data['Cl_tck']
        Cd_tck = self.airfoil_data['Cd_tck']
        
        # reshape vectors
        beta = np.linspace(0, 2 * np.pi, Nbeta).reshape(-1,1)
        seg = self.seg.reshape(1,-1)
        c = self.chord.reshape(1,-1)
        twist = self.twist.reshape(1,-1)
        blade_length = self.seg[-1] + 0.5 * self.Lspan[-1]

        # coordinate of the segment
        h_temp = (np.cos(beta) * seg + self.href).reshape(len(beta), seg.shape[1])         #hight of current segment

        # interpolate atmospheric field on segment position (cubic interpolation)
        U_inf = interp_atmos_data(atmos.z_coord,atmos.U_inf,h_temp)
        U_rot = seg * self.omega
        
        # initialize matrices
        AoA_seg = np.zeros((beta.shape[0], seg.shape[1]))
        rey_seg = np.zeros((beta.shape[0], seg.shape[1]))
        a_array = np.zeros((beta.shape[0], seg.shape[1]))
        adash_array = np.zeros((beta.shape[0], seg.shape[1]))
        epsilon_array = np.zeros((beta.shape[0], seg.shape[1]))
        U_rel = np.zeros((beta.shape[0], seg.shape[1]))
        F_array = np.zeros((beta.shape[0], seg.shape[1]))
        colors=mcp.gen_color(cmap="copper",n=self.Nseg)
        #colors = ['k','b','r','g','orange','brown','steelblue', 'plum','gray']
        for ibeta in range(beta.shape[0]):
            for iseg in range(seg.shape[1]):   
                U_rel[ibeta, iseg] = np.sqrt(U_inf[ibeta, iseg] ** 2 + U_rot[0, iseg] ** 2)     
                Re = U_rel[ibeta, iseg] * c[0, iseg] / nu0
                J= self.omega * blade_length / U_inf[ibeta, iseg]
                theta = 3*c[0,iseg]/(2*np.pi*seg[0,iseg])            
                F = 2 / np.pi * np.arccos(np.exp(-(3 / 2) * (1 - seg[0, iseg] / blade_length) * (1 + J ** 2) ** 0.5))
        
                if BEM : 
                    # AoA_seg[ibeta,iseg],a,adash,F,epsilon, a_conv,adash_conv,alpha_conv = bem.simple(Cl_tck[iseg],
                    #                                                                                  Cd_tck[iseg], 
                    #                                                                                  Re,J,theta,
                    #                                                                                  twist[0,iseg],
                    #                                                                                  blade_length,
                    #                                                                                  seg[0,iseg],F)

                    # AoA_seg[ibeta,iseg],a,adash,F,epsilon,a_conv,adash_conv,alpha_conv = bem.hemant1(Cl_tck[iseg],Cd_tck[iseg],Re,J,blade_length,seg[0,iseg],twist[0,iseg],c[0,iseg])
                    AoA_seg[ibeta,iseg],a,adash,F,epsilon, a_conv,adash_conv,alpha_conv = bem.hemant3(Cl_tck[iseg],Cd_tck[iseg],
                                                                                            Re,J,blade_length,
                                                                                            seg[0,iseg],twist[0,iseg],c[0,iseg],
                                                                                            U_inf[ibeta,iseg],U_rot[0,iseg])

                    plt.figure(10)
                    plt.plot(a_conv,colors[iseg])
                    plt.figure(11)
                    plt.plot(adash_conv,colors[iseg])
                    plt.figure(12)
                    plt.plot(np.array(alpha_conv)*180/np.pi,colors[iseg])

                else :
                     AoA_seg[ibeta,iseg],a,adash,F,epsilon = bem.noBEM(Cl_tck[iseg],Cd_tck[iseg],Re,J,blade_length,seg[0,iseg],twist[0,iseg],c[0,iseg])


                a_array[ibeta,iseg] = a
                adash_array[ibeta,iseg] = adash
                epsilon_array[ibeta,iseg] = epsilon
                F_array[ibeta,iseg] = F
                # U_rel[ibeta,iseg] = ((U_inf[ibeta, iseg]*(1-a))**2 +(U_rot[0,iseg]*(1+adash))**2)**0.5
                U_rel[ibeta, iseg] = np.sqrt((U_inf[ibeta, iseg] * (1 - F * a)) ** 2 + (U_rot[0, iseg] * (1 + adash)) ** 2)
                rey_seg[ibeta, iseg] = Re

        U_rot = np.repeat(U_rot.reshape(1,-1),beta.shape[0],axis=0)

        B,R = np.meshgrid(np.squeeze(beta),self.seg,indexing='ij')

        print(adash_array.shape)
        print(R.shape)

        fig = plt.figure(0)
        ax = fig.add_subplot(231, polar = True)
        cax = ax.pcolormesh(B + np.pi/2,R,AoA_seg*180/np.pi,shading='auto')#,edgecolors='k')   
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        cbar = plt.colorbar(cax)
        cbar.ax.set_title("AoA")

        fig = plt.figure(0)
        ax = fig.add_subplot(232, polar = True)
        cax = ax.pcolormesh(B + np.pi/2,R,a_array,shading='auto')#,edgecolors='k')   
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        cbar = plt.colorbar(cax)
        cbar.ax.set_title("a")
        
        fig = plt.figure(0)
        ax = fig.add_subplot(233, polar = True)
        cax = ax.pcolormesh(B + np.pi/2,R,adash_array,shading='auto')#,edgecolors='k')   
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        cbar = plt.colorbar(cax)
        cbar.ax.set_title("a'")

        ax = fig.add_subplot(234, polar = True)
        cax = ax.pcolormesh(B + np.pi/2,R,U_rel,shading='auto')#,edgecolors='k')   
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        cbar = plt.colorbar(cax)
        cbar.ax.set_title("U_rel")

        ax = fig.add_subplot(235, polar = True)
        cax = ax.pcolormesh(B + np.pi/2,R,U_inf,shading='auto')#,edgecolors='k')   
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        cbar = plt.colorbar(cax)
        cbar.ax.set_title("U_inf")

        ax = fig.add_subplot(236, polar = True)
        cax = ax.pcolormesh(B + np.pi/2,R,U_rot,shading='auto')#,edgecolors='k')   
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        cbar = plt.colorbar(cax)
        cbar.ax.set_title("U_rot")

        plt.tight_layout()


        fig = plt.figure(100)
        ax = fig.add_subplot(231, polar = True)
        cax = ax.pcolormesh(B + np.pi/2,R,F_array,cmap='RdBu_r',shading='auto')#,edgecolors='k')   
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        plt.colorbar(cax)
        cbar = cbar.ax.set_title("F")

        plt.tight_layout()

    def plotPolar(self,iseg,Re,AoA):

        Re = np.squeeze(np.array([Re]))
        AoA = np.squeeze(np.array([AoA]))

        Cl = interp.bisplev(np.array(Re),AoA,self.airfoil_data['Cl_tck'][iseg])
        Cd = interp.bisplev(np.array(Re),AoA,self.airfoil_data['Cd_tck'][iseg])

        if Re.size>1 and AoA.size >1: 
            fig = plt.figure(1)
            ax = fig.add_subplot(211)
            ax.pcolormesh(AoA, Re, Cl)
            ax.set_xlabel('AoA')
            ax.set_ylabel('Re')
            ax.set_title('Cl')

            ax = fig.add_subplot(212)
            ax.pcolormesh(AoA, Re, Cd)
            ax.set_xlabel('AoA')
            ax.set_ylabel('Re')
            ax.set_title('Cd')
            plt.tight_layout()
            plt.plot(AoA,Cl)
            plt.plot(self.airfoil_data['AoA'],self.airfoil_data['Cl'][iR,:,iseg],'+')
            plt.xlabel('AoA')
            plt.ylabel('Cl')

            plt.subplot(212)
            plt.plot(AoA,Cd)
            plt.plot(self.airfoil_data['AoA'],self.airfoil_data['Cd'][iR,:,iseg],'+')
            plt.xlabel('AoA')
            plt.ylabel('Cd')
            ax = fig.add_subplot(212)
            ax.pcolormesh(self.airfoil_data['AoA'],self.airfoil_data['reynolds'],self.airfoil_data['Cd'][:,:,iseg])
            ax.set_xlabel('AoA')
            ax.set_ylabel('Re')
            ax.set_title('Cd')
            plt.tight_layout()

        if Re.size==1 and AoA.size>1: 
            
            plt.figure(1)
            plt.subplot(211)
            iR = np.argmin(np.abs(Re - self.airfoil_data['reynolds']))
            ax = plt.gca()
            color = next(ax._get_lines.prop_cycler)['color']
            plt.plot(AoA,Cl,color=color)
            plt.plot(self.airfoil_data['AoA'],self.airfoil_data['Cl'][iR,:,iseg],'+',color=color)
            plt.xlabel('AoA')
            plt.ylabel('Cl')

            plt.subplot(212)
            plt.plot(AoA,Cd,color=color,label= self.airfoil_data['airfoilDataName'][self.airfoil_data['airfoildIndex'][iseg]])
            plt.plot(self.airfoil_data['AoA'],self.airfoil_data['Cd'][iR,:,iseg],'+',color=color)
            plt.xlabel('AoA')
            plt.ylabel('Cd')

            plt.legend()
            plt.tight_layout()


    # calculate WPS acording to the blade segment, Angle of attack and wind speed (needed for exploreDataBase)
    def Phi_pp(self,freq,iseg,AoA,U_inf):
        rho0 = 1.225
        nu0 = 1.45e-5 

        Nfreq=len(freq)

        AoA_array = self.airfoil_data['AoA'].reshape(-1)
        iAoA = WindTurbine.find_index(AoA, AoA_array)

        delta_star_top = self.airfoil_data['delta_star_top'][:,iAoA,iseg]
        theta_momen_top = self.airfoil_data['theta_momen_top'][:,iAoA,iseg]
        dpdx_top = self.airfoil_data['dpdx_top'][:,iAoA,iseg]
        UeUinf_top = self.airfoil_data['UeUinf_top'][:,iAoA,iseg]
        cf_top = self.airfoil_data['cf_top'][:,iAoA,iseg]
        seg = self.airfoil_data['seg'][iseg]
        chord = self.airfoil_data["chord"][iseg]
        Reynolds = self.airfoil_data['reynolds'].reshape(-1)

        self.controlRotSpeed(U_inf)
        
        # Omega = 12.1 * 2 * np.pi / 60
        U_rot = seg * self.omega
        U_rel = np.sqrt(U_inf ** 2 + U_rot ** 2)
        Reyn = U_rel* chord / nu0

        ind_re = WindTurbine.find_index(Reyn, Reynolds)
        ind_left = ind_re - 1
        ind_right = ind_re
        ratio = (Reyn-Reynolds[ind_left]) / (Reynolds[ind_right] - Reynolds[ind_left]) 
        print(ratio)
        delta_star_t = delta_star_top[ind_left] + ratio * (delta_star_top[ind_right] - delta_star_top[ind_left])
        theta_momen_t = theta_momen_top[ind_left] + ratio * (theta_momen_top[ind_right] - theta_momen_top[ind_left])
        dpdx_t = dpdx_top[ind_left] + ratio * (dpdx_top[ind_right] - dpdx_top[ind_left])
        cf_t = cf_top[ind_left] + ratio * (cf_top[ind_right] - cf_top[ind_left])
        UeUinf_t = UeUinf_top[ind_left] + ratio * (UeUinf_top[ind_right] - UeUinf_top[ind_left])

        Ue_t = U_rel * np.abs(UeUinf_t)                    # external velocity (m)
        tauw_t = cf_t * 0.5 * rho0 * Ue_t ** 2                         # wall shear stress (Pa)
        betac_t = theta_momen_t / tauw_t * dpdx_t                      # Clauser parameter
        print('test')
        print(betac_t)
        print(dpdx_t)
        print(delta_star_t)
        print(theta_momen_t)
        print(dpdx_t)
        print(cf_t)
        print(UeUinf_t)
        print(U_rel)
        if cf_t<0 :
            Phi_pp_adim = 0
            Delta = 0 
            beta_c = betac_t 
        elif betac_t > 60:       
            dpdx_t = 60 / theta_momen_t * tauw_t
            Phi_pp_suct1,phisucc2, Delta,beta_c = WindTurbine._WP_LEE(freq.reshape(1, Nfreq), delta_star_t, theta_momen_t, dpdx_t, cf_t, UeUinf_t, U_rel)
            Phi_pp_adim = np.squeeze(10*np.log10(np.squeeze(Phi_pp_suct1)*Ue_t/(tauw_t **2*delta_star_t)))
        else:
            Phi_pp_suct1,phisucc2, Delta,beta_c = WindTurbine._WP_LEE(freq.reshape(1, Nfreq), delta_star_t, theta_momen_t, dpdx_t, cf_t, UeUinf_t, U_rel)
            Phi_pp_adim = np.squeeze(10*np.log10(np.squeeze(Phi_pp_suct1)*Ue_t/(tauw_t **2*delta_star_t)))

        freq_adim = freq*delta_star_t/Ue_t 
        airfoilType = self.airfoil_data['airfoilDataName'][self.airfoil_data["airfoildIndex"][iseg]]
        chord = self.airfoil_data["chord"][iseg]
        return freq_adim, Phi_pp_adim, airfoilType, chord, Delta, beta_c


    # allow to compute and the wall pressure spectrum for all the blade segment and possible velocity and AoA 
    def exploreDataBase(self):
        # init frequency for the spectrum 
        fc = [5, 6.3, 8, 10, 12.5, 16, 20, 25,  31.5, 40, 50, 63, 80,100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000]
        Nfc = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,1,1,2,2,3,4,4,4,5,5,5,6,6,6,7,7,10]
        freq = computeThirdOctaveFrequencies(fc,Nfc)
        AoA_array = self.airfoil_data['AoA'].reshape(-1)
        seg_array = self.airfoil_data['seg']
        # init the figure 
        fig, (ax1,ax2) = plt.subplots(1,2)
        line, = ax1.semilogx(self.Phi_pp(freq, 1,0,8)[0], self.Phi_pp(freq, 1,0,8)[1], lw=2)
        line.set_label(self.Phi_pp(freq, 1,0,8)[2] + '; chord=' + str(round(self.Phi_pp(freq, 1,0,8)[3],0)) + 'm')
        ax1.set_xlabel('$\omega \delta^* /U_e$')
        ax1.set_ylabel('$\Phi_{pp} U_e / \\tau \delta^*$')
        ax1.set_ylim(-50,40)
        ax1.set_xlim(0.001,10)

        line2, = ax2.plot(self.Phi_pp(freq, 1,0,8)[4], self.Phi_pp(freq, 1,0,8)[5],'*' )
        ax2.set_ylim(0,70)
        ax2.set_xlim(0,10)
        ax2.set_xlabel('$\Delta$')
        ax2.set_ylabel('$\\beta_c$')

        # create de sliders to modify the plot 
        fig.subplots_adjust( bottom=0.3)

        axAoA = fig.add_axes([0.25, 0.15, 0.65, 0.03])
        AoA_slider = Slider(
            ax=axAoA,
            label='AoA',
            valmin = AoA_array[0],
            valmax = AoA_array[-1],
            valstep = AoA_array,
            valinit=0,
        )

        axSeg = fig.add_axes([0.25, 0.05, 0.65, 0.03])
        seg_slider = Slider(
            ax=axSeg,
            label="Seg index",
            valmin=0,
            valmax=len(seg_array )-1,
            valinit=0,
            valstep=1,
        )

        axU = fig.add_axes([0.25, 0.1, 0.65, 0.03])
        U_slider = Slider(
            ax=axU,
            label="U ",
            valmin=1,
            valmax=25,
            valinit=8,
            valstep=1,
        )

        # The function to be called anytime a slider's value changes
        def update(val):
            Phi = self.Phi_pp(freq,seg_slider.val, AoA_slider.val,U_slider.val)
            line.set_xdata(Phi[0])
            line.set_ydata(Phi[1])
            line2.set_xdata(Phi[4])
            line2.set_ydata(Phi[5])
            line.set_label(Phi[2] + '; chord=' + str(round(Phi[3],1)) + 'm')
            ax1.legend()
            # ax.set_xlim(min(Phi[0]),max(Phi[0]))
            # ax.set_ylim(min(Phi[1]),max(Phi[1]))
            fig.canvas.draw_idle()

        # register the update function with each slider
        AoA_slider.on_changed(update)
        seg_slider.on_changed(update)
        U_slider.on_changed(update)

        ax1.legend()
        ax1.grid(True)
        ax2.grid(True)
        plt.show()

    def save(self,fname: str):
        """save class as self.name.dat"""
        with open(fname,'wb') as file:
            pickle.dump(self.__dict__,file)

    def load(self,fname: str):
        """try load self.name.dat"""
        with open(fname,'rb') as file:
            self.__dict__ = pickle.load(file)

if __name__ == "__main__":
    wind_turbine = WindTurbine()

    seg = np.asarray([11.75, 15.85, 19.95, 24.05, 28.15, 32.25, 36.35, 40.45, 44.55, 48.65, 52.75, 56.1667, 58.9, 61.6333])
    Lspan = np.asarray([ 4.1, 4.1, 4.1, 4.1, 4.1, 4.1, 4.1, 4.1, 4.1, 4.1, 4.1, 2.7333, 2.7333, 2.7333])
    c = np.asarray([4.557, 4.652, 4.458, 4.249, 4.007, 3.748, 3.502, 3.256, 3.010, 2.764, 2.518, 2.313, 2.086, 1.419])
    twist = (90 - np.asarray([ 13.308, 11.48, 10.162, 9.011, 7.795, 6.544, 5.361, 4.188, 3.125, 2.319, 1.526, 0.863, 0.37, 0.106])) *np.pi/180
    airfoilType = np.asarray([0,1,1,2,3,3,4,4,5,5,5,5,5,5])
    wind_turbine.airfoilDataName = np.asarray(["DU40", "DU35", "DU30", "DU25", "DU21", "NACA64"])

    wind_turbine.pathToAirfoilGeom  = './DU_geom/'
    wind_turbine.createGeometry(seg,Lspan,c,twist,airfoilType)
    # reynolds = np.array([7000000]) #np.linspace(7000000, 13000000, 2)
    # AoA = np.array([1])   #np.linspace(1, 13, 2)

    reynolds = np.linspace(7000000, 13000000, 7)
    AoA = np.linspace(1, 13, 13)
    print(AoA)
    # quit()
    fname = 'NREL_airfoilt_data.dat'
    wind_turbine.createBLData2(AoA,reynolds,fname)
    wind_turbine.save('NREL_wt.dat')


    # orginal 
