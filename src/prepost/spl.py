import numpy as np
from math import ceil
from tqdm import tqdm
import pickle
import matplotlib.pyplot as plt 
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.animation import FuncAnimation, FFMpegWriter
import numexpr as ne
import h5py
import logging
from time import time, sleep
from .les import Les
from .utils import (integrateThirdOctave, cos_window,
                    uneven_tile,uneven_loop, interp_weights,interpolate )
from .source.main import Source
from .deltaLfield import DeltaLField
from .wape import  PeResults
from .pre import Simu


class SplField():
    """
    A class to handle and process SPL (Sound Pressure Level) fields generated from PE (Parabolic Equation) simulation results.
    The SplField class is designed to handle and process SPL (Sound Pressure Level) fields generated from PE (Parabolic Equation)simulation results.
    It provides fonctionaly to read and combine deltaLField with source power levels. 
    It can also read, store, manipulate, and visualize SPL data. 
    Several post processing like atmospheric absorption, Aweighing, third octave band integration can be performed. 
    It can also compute the SPL in time domain, compute the time averaged over a rotation or the Amplitude modulation. 
    Finally functions for visualisation are provided. 

    Attributes:
        src (Source): The source object containing the Sound power level data.
        deltaL (DeltaLField): The deltaL field object containing the deltaL data.
        SPL_time (np.ndarray): The SPL data in the time domain. Shape: (x, y, z, freq, t).
        SPL_seg (np.ndarray): The SPL data in the angular domain. Shape: (x, y, z, seg, blade, freq, beta).
        OASPL_time (np.ndarray): The overall SPL data in the time domain. Shape: (x, y, z, t).
        OASPL_seg (np.ndarray): The overall SPL data in the angular domain. Shape: (x, y, z, seg, blade, beta).
        OASPL_seg_tot (np.ndarray): The total overall SPL data in the angular domain. Shape: (x, y, z, beta).
        am_seg (np.ndarray): The amplitude modulation (AM) data in the angular domain.
        am_time (np.ndarray): The AM data in the time domain.
        OAam_seg (np.ndarray): The overall AM data in the angular domain.
        OAam_time (np.ndarray): The overall AM data in the time domain.
        mean_seg (np.ndarray): The mean data in the angular domain.
        mean_time (np.ndarray): The mean data in the time domain.
        OAmean_seg (np.ndarray): The overall mean data in the angular domain.
        OAmean_time (np.ndarray): The overall mean data in the time domain.
        Nx (int): The number of points in the x-direction.
        Ny (int): The number of points in the y-direction.
        Nz (int): The number of points in the z-direction.
        Nfreq (int): The number of frequency bands.
        Nbeta (int): The number of beta angles.
        x (np.ndarray): The x-coordinates.
        y (np.ndarray): The y-coordinates.
        third (bool): Flag indicating if the data is in third-octave bands. Default is False.
        oaspl (bool): Flag indicating if the overall SPL data is computed. Default is False.
        FULL_ROTATION (bool): Flag indicating if the full rotation is applied. Default is False.
        AWEIGHT (bool): Flag indicating if A-weighting is applied. Default is False.
        ATM_ABS (bool): Flag indicating if atmospheric absorption is applied. Default is False.
        POLAR (bool): Flag indicating if the data is in polar coordinates. Default is False.
    """

    src: Source
    deltaL: DeltaLField
    # shape: (x,y,z,freq,t)
    SPL_time: np.ndarray = None
    # shape: (x,y,z,seg,blade,freq,beta)
    SPL_seg: np.ndarray = None
    # shape: (x,y,z,t)
    OASPL_time: np.ndarray = None
    # shape: (x,y,z,seg,blade,beta)
    OASPL_seg: np.ndarray = None
    # shape: (x,y,z,beta)
    OASPL_seg_tot: np.ndarray = None

    am_seg: np.ndarray = None
    am_time: np.ndarray = None
    OAam_seg: np.ndarray = None
    OAam_time: np.ndarray = None
    mean_seg: np.ndarray = None
    mean_time: np.ndarray = None
    OAmean_seg: np.ndarray = None
    OAmean_time: np.ndarray = None

    Nx: int
    Ny: int
    Nz: int
    Nfreq: int
    Nbeta: int
    x: np.ndarray
    y: np.ndarray
    third: bool = False
    oaspl: bool = False

    FULL_ROTATION: bool = False
    AWEIGHT: bool = False
    ATM_ABS: bool = False
    POLAR = False

    def __init__(self, src: Source = None, deltaL: DeltaLField = None):
        """
        Initializes the SplField class.

        Args:
            src (Source, optional): The source object containing the sound pressure level data. Default is None.
            deltaL (DeltaLField, optional): The deltaL field object containing the deltaL data. Default is None.
        """
        if src is not None:
            self.src = src
            self.wt = src.wt
            self.Nbeta = self.wt.Nbeta

        if deltaL is not None:
            self.deltaL = deltaL
            self.xS = deltaL.xS
            self.yS = deltaL.yS

    def info(self):
        """
        Prints information about the SPL field.
        Show the different quantities loaded and the FLAG to assess what post processing were already run.
        """
        print('MESH:')

        print('quantity loaded: ')
        if self.SPL_seg is not None:
            print('SPL_seg')
        if self.OASPL_seg is not None:
            print('OASPL_seg')

        if self.SPL_time is not None:
            print('SPL_time')
        if self.OASPL_time is not None:
            print('OASPL_time')

        print('flags')

        print('oaspl', self.oaspl)
        print('third', self.third)
        print('FULL_ROTATION', self.FULL_ROTATION)
        print('AWEIGHT', self.AWEIGHT)
        print('ATM_ABS', self.ATM_ABS)

    def check_compatibility_polar(self):
        """
        Checks the compatibility between the `Source` data and `DelteLField` data in polar coordinates.
        This is usually not used and we prefer to convert $Spp$ and $\Delta L$ in cartesian coordinates and then 
        use check_compatibility_cartesian.

        Returns:
            int: Returns -1 if the data is not compatible.
        """
        shape_Spp = self.src.SppInterpolated[..., 0, 0, 0].shape
        if self.deltaL.deltaL_polar is None:
            logging.warning(" no deltaL polar in the data")
            return -1
        shape_deltaL = self.deltaL.deltaL_polar[..., 0, 0].shape

        # check space shape
        if shape_deltaL != shape_Spp:
            print('Spp and delta L must be of same size')
            print(shape_deltaL)
            print(shape_Spp)

        if not np.all(self.src.x_interpolate[:, :, 0] == self.deltaL.x_polar):
            print('x grids are not the same')
            return -1
        else:
            print('x grid ok')
        if not np.all(self.src.y_interpolate[:, :, 0] == self.deltaL.y_polar):
            print('y grids are not the same')
            # print(self.src.y_interpolate)
            # print(self.deltaL.y_polar)
            return -1
        else:
            print('y grid ok')

        if not np.all(self.src.z_interpolate[0, 0, :] == self.deltaL.z_polar):
            print('z grids are not the same')
            return -1
        else:
            print('z grid ok')

        # check frequency
        # Nfreq_Spp = self.src.SppInterpolated.shape[-2]
        # Nfreq_deltaL = self.deltaL.deltaLInterpolated.shape[-2]
        if not np.all(self.src.frequencies == self.deltaL.frequencies):
            print('frequencies are not the same')
            print(self.src.frequencies)
            print(self.deltaL.frequencies)
            return -1
        else:
            self.frequencies = self.src.frequencies
            print('frequencies ok')
        self.Nx = self.src.x_interpolate.shape[0]
        self.Ny = self.src.x_interpolate.shape[1]
        self.Nz = self.src.x_interpolate.shape[2]
        self.Nfreq = self.frequencies.size

    def check_compatibility_cart(self):
        """
        Checks the compatibility between the `Source` data and `DelteLField` data in cartesian coordinates.

        Returns:
            int: Returns -1 if the data is not compatible.
        """

        shape_Spp = self.src.SppInterpolated[..., 0, 0, 0].shape
        if self.deltaL.deltaL_cart is None:
            logging.warning(" no delta_L cart in the ddat")
            return -1
        shape_deltaL = self.deltaL.deltaL_cart[..., 0, 0].shape

        # check space shape
        if shape_deltaL != shape_Spp:
            print('Spp and delta L must be of same size')
            print(shape_deltaL)
            print(shape_Spp)

        if not np.all(self.src.x_interpolate == self.deltaL.x_cart):
            print('x grids are not the same')
            return -1
        else:
            print('x grid ok')
        if not np.all(self.src.y_interpolate == self.deltaL.y_cart):
            print('y grids are not the same')
            print(self.src.y_interpolate)
            print(self.deltaL.y_cart)
            return -1
        else:
            print('y grid ok')
        if not np.all(self.src.z_interpolate == self.deltaL.z_cart):
            print('z grids are not the same')
            return -1
        else:
            print('z grid ok')

        # check frequency
        # Nfreq_Spp = self.src.SppInterpolated.shape[-2]
        # Nfreq_deltaL = self.deltaL.deltaLInterpolated.shape[-2]
        if not np.all(self.src.frequencies == self.deltaL.frequencies):
            print('frequencies are not the same')
            print(self.src.frequencies)
            print(self.deltaL.frequencies)
            return -1
        else:
            self.frequencies = self.src.frequencies
            print('frequencies ok')
        self.Nx = self.src.x_interpolate.shape[0]
        self.Ny = self.src.x_interpolate.shape[1]
        self.Nz = self.src.x_interpolate.shape[2]
        self.Nfreq = self.frequencies.size

    # interpolate source and deltaL on the same grid
    def interpolate(self, x: np.array, y: np.array):
        """
        Interpolates the source and deltaL on the same grid.

        Args:
            x (np.ndarray): The x-coordinates for interpolation.
            y (np.ndarray): The y-coordinates for interpolation.
        """
        self.src.interpolate_xy(x, y)
        self.deltaL.interpolate_xy(x, y)
        self.x = x
        self.y = y

    def combine_linear_broadband(self,free_field=False):
        """
        Combines the $S_{pp}$ and $\Delta L$ fields with linear interpolation between to source heights.
        First the function check if the two fields are compatible then the sum of $S_{pp}$ and $\Delta L$ is done for each angular position, segment and frequency.


        Args:
            free_field (bool): Flag indicating if only the free field result must be coomputed, disregarding the Delta L. Default is False.
        """

        # 2D OASPL computation
        # ---------------------------------------------------------------------
        print('combining ...')
        proximity = self.src.wt.proximityLinear(self.deltaL.height)
        shape_Spp = self.src.SppInterpolated[..., 0, 0, 0].shape
        shape_deltaL = self.deltaL.deltaL_cart[..., 0, 0].shape
        if shape_deltaL != shape_Spp:
            print('Spp and delta L must be of same size')
            print(shape_deltaL)
            print(shape_Spp)
            quit()
        if len(shape_Spp) == 2:
            (nx, nz) = shape_Spp
            ny = 1
        if len(shape_Spp) == 3:
            (nx, ny, nz) = shape_Spp
        elif len(shape_Spp) == 1:
            (nx,) = shape_Spp
            nz = 1
            ny = 1
        else:
            print('space dimension not handle by the actual code ')
            quit()

        final_SPL_seg_freq_time = np.zeros(np.shape(self.src.SppInterpolated))

        for ibeta in np.arange(0, self.src.wt.Nbeta):
            print('ibeta = ' + str(ibeta))
            # create coresponding matrix for Spp and Delta L fo one angle 
            # ---------------------------------------------------------------
            SPL_tot_ff = self.src.SppInterpolated[..., ibeta]
            delta_L = np.zeros_like(SPL_tot_ff)
            for iseg in range(self.src.wt.Nblade*self.src.wt.Nseg):

                delta_L[..., iseg, :] = 10*np.log10(proximity[iseg, ibeta, 2] *
                                        10**(self.deltaL.deltaL_cart[..., proximity[iseg, ibeta, 0]]/10) + 
                                            (1-proximity[iseg, ibeta, 2])*
                                        10**(self.deltaL.deltaL_cart[..., proximity[iseg, ibeta, 1]]/10))

            #   sum Spp and Delta L
            #------------------------------------------------------------------
            if free_field:
                final_SPL_seg_freq_time[..., ibeta] = SPL_tot_ff
            else:
                final_SPL_seg_freq_time[..., ibeta] = SPL_tot_ff + delta_L

        self.SPL_seg = final_SPL_seg_freq_time.reshape(nx, ny, nz,
                                                       self.src.wt.Nseg,
                                                       self.src.wt.Nblade,
                                                       len(self.frequencies),
                                                       self.src.Nbeta)
        print('done combining')
        self.x_grid = self.deltaL.x_cart
        self.y_grid = self.deltaL.y_cart
        self.z_grid = self.deltaL.z_cart

    def interpolate_from_polar(self, x: np.ndarray, y: np.ndarray):
        """
        Interpolates the SPL field from polar coordinates to Cartesian coordinates.
        This is done if the the Delta L and Spp where combined in polar coordinates in the first place.

        Args:
            x (np.ndarray): The x-coordinates for interpolation.
            y (np.ndarray): The y-coordinates for interpolation.
        """ 

        xy_polar = np.concatenate(
                (self.x_grid[:, :, 0].reshape((-1, 1)),
                 self.y_grid[:, :, 0].reshape((-1, 1))), 1)
        print("start interpolation ...")

        x_cart, y_cart = np.meshgrid(x, y, indexing="ij")
        xy_cart = np.zeros([x_cart.shape[0] * x_cart.shape[1], 2])
        xy_cart[:, 0] = x_cart.flatten()
        xy_cart[:, 1] = y_cart.flatten()
        # create knots for the interpolation

        vtx, wts = interp_weights(xy_polar, xy_cart)
        print("finished creating knot ...")

        # loop of z, freq, heights and interpolate
        # using the previously computes knots
        print("starting loop on height and frequency band ...")
        tag = False
        # ---------------------------------------------------------------------
        if self.OASPL_time is not None:
            print("interpolate OASPL_time")
            OASPL_cart = np.zeros((x_cart.shape[0],
                                   x_cart.shape[1],
                                   1, self.Nt))

            for it in range(self.Nt):
                values = self.OASPL_time[:, :, 0, it].flatten()
                valuesInterpolated = interpolate(values, vtx, wts)
                OASPL_cart[:, :, 0, it] = valuesInterpolated.reshape(
                        x_cart.shape[0], x_cart.shape[1])
            self.OASPL_time = OASPL_cart
            tag = True
        # ---------------------------------------------------------------------
        if self.OASPL_seg_tot is not None:
            print("interpolate OASPL_seg_tot")
            OASPL_cart = np.zeros((x_cart.shape[0],
                                   x_cart.shape[1],
                                   1, self.Nbeta))
            for it in range(self.Nbeta):
                values = self.OASPL_seg_tot[:, :, 0, it].flatten()
                valuesInterpolated = interpolate(values, vtx, wts)
                OASPL_cart[:, :, 0, it] = valuesInterpolated.reshape(
                        x_cart.shape[0], x_cart.shape[1])
            self.OASPL_seg_tot = OASPL_cart
            tag = True

        # ---------------------------------------------------------------------
        if self.OASPL_seg is not None:
            print(self.OASPL_seg.shape)
            print("interpolate OASPL_seg")
            OASPL_cart = np.zeros((x_cart.shape[0],
                                   x_cart.shape[1],
                                   1,
                                   self.wt.Nseg, self.wt.Nblade,
                                   self.wt.Nbeta))
            for iseg in range(self.wt.Nseg):
                for iblade in range(self.wt.Nblade):
                    for it in range(self.wt.Nbeta):
                        values = self.OASPL_seg[:, :, 0, iseg,
                                                iblade, it].flatten()
                        valuesInterpolated = interpolate(values, vtx, wts)
                        OASPL_cart[:, :, 0, iseg,
                                   iblade,
                                    it] = valuesInterpolated.reshape(
                                           x_cart.shape[0], x_cart.shape[1])
            self.OASPL_seg = OASPL_cart
            tag = True

        # ---------------------------------------------------------------------
        if self.SPL_seg is not None:
            print("interpolate SPL_seg")
            SPL_cart = np.zeros((x_cart.shape[0],
                                x_cart.shape[1],
                                1,
                                self.wt.Nseg, self.wt.Nblade,
                                self.Nfreq, self.wt.Nbeta))
            print(SPL_cart.shape)

            for iseg in range(self.wt.Nseg):
                for iblade in range(self.wt.Nblade):
                    for ifreq in range(self.Nfreq):
                        for it in range(self.wt.Nbeta):
                            values = self.SPL_seg[:, :, 0, iseg,
                                                iblade, ifreq, it].flatten()
                            valuesInterpolated = interpolate(values, vtx, wts)
                            SPL_cart[:, :, 0, iseg,
                                   iblade,
                                   ifreq, it] = valuesInterpolated.reshape(
                                           x_cart.shape[0], x_cart.shape[1])
            self.SPL_seg = SPL_cart
            tag = True

        if tag:
            [self.x_grid, self.y_grid, self.z_grid] = np.meshgrid(
                    x, y, self.z_grid[0, 0, :], indexing="ij")
            self.Nx = x.size
            self.Ny = y.size
        print("done")


    def create_full_rotation(self):
        """
        The computation are usual done for a third of rotation which is sufficient to have the complete rotation because there are 3 blades.
        This function loops the results three time over the Blade angle to obtain a complete rotation of the rotor.
        """
        if self.FULL_ROTATION:
            logging.warning("solution is already convert to full rotation")
            return 
        if self.SPL_seg is not None:
            self.SPL_seg = np.concatenate((self.SPL_seg,
                                           self.SPL_seg[:, :, :, :,
                                                        [1, 2, 0], :, :],
                                           self.SPL_seg[:, :, :, :,
                                                        [2, 0, 1], :, :]),
                                          axis=6)
        if self.OASPL_seg is not None:
            self.OASPL_seg = np.concatenate((self.OASPL_seg,
                                           self.OASPL_seg[:, :, :, :,
                                                          [1, 2, 0], :],
                                           self.OASPL_seg[:, :, :, :,
                                                          [2, 0, 1], :]),
                                            axis=5)
        if self.OASPL_seg_tot is not None:
            self.OASPL_seg_tot = np.concatenate((self.OASPL_seg_tot,
                                           self.OASPL_seg_tot,
                                           self.OASPL_seg_tot),
                                            axis=3)
        self.Nbeta = self.Nbeta*3
        self.FULL_ROTATION = True
        print('create full rotation')

    def clean_full_rotation(self):
        """
        Cleans the full rotation of the SPL field.
        Check if full rotation was computed before removing the additional angle.
        """

        if not self.FULL_ROTATION:
            logging.warning("solution was not convert to full rotation")
            return
        if self.SPL_seg is not None:
            self.SPL_seg = self.SPL_seg[:, :, :, :, :, :, 0:self.Nbeta//3]
            logging.info('SPL_seg shape : ' + str(self.SPL_seg.shape))
        self.Nbeta = self.Nbeta//3
        self.FULL_ROTATION = False
        print('cleaned full rotation')

    def compute_receiver_time(self):
        """
        Computes the receiver time for each receiver/source pair.
        The distance between each blade segment and receiver at each angular position is computed.
        Then the propagation time is computed with a reference wind speed of c0=343ms.
        The modulo of the propagation time with respect to a complete rotation is then computed. 

        Hence each receiver time is comprise between 0 and $2\pi / \Omega$.
        """

        print('compute receiver time ...')
        c0 = 343
        # (nx,ny,nz,src.wt.Nseg,src.wt.Nblade,len(freq),src.Nbeta)
        x = self.x_grid[:, 0, 0].reshape(-1, 1, 1, 1, 1, 1)
        y = self.y_grid[0, :, 0].reshape(1, -1, 1, 1, 1, 1)
        z = self.z_grid[0, 0, :].reshape(1, 1, -1, 1, 1, 1)

        self.wt.computeBeta()
        beta = self.wt.beta.reshape(1, 1, 1, 1, self.wt.Nblade, self.wt.Nbeta)
        beta = np.concatenate((beta, beta+2*np.pi/3, beta+4*np.pi/3), axis=5)
        seg = self.wt.seg.reshape(1, 1, 1, self.wt.Nseg, 1, 1)

        # compute segment location
        # modif source :
        zS = (np.cos(beta) * seg + self.wt.href)
        xS = -np.sin(beta) * seg * np.sin(self.wt.tau) + self.xS
        yS = -np.sin(beta) * seg * np.cos(self.wt.tau) + self.yS

        # test
        # zS = (np.cos(beta) * seg *0. + self.wt.href)
        # xS = -np.sin(beta) * seg * np.sin(self.wt.tau) * 0. + self.xS
        # yS = -np.sin(beta) * seg * np.cos(self.wt.tau) * 0. + self.yS

        # compute distance between source and receiver
        R = np.sqrt((x - xS)**2 + (y - yS)**2 + (z - zS)**2)

        # absolute time of signal reception
        self.t = R/c0 + beta[:, :, :, :, 0, :].reshape(1, 1, 1,
                                                       1, 1, -1)/self.wt.omega

        # copy first angle at the end
        self.t = np.concatenate((self.t, self.t[..., 0:1]), 5)

        # modulo one rotation
        self.t = (self.t) % (2*np.pi/(self.wt.omega))
        print('done.')
        print(self.t.shape)

    def compute_real_receiver_time(self,loop=True,last=True):
        """
        Computes the receiver time for each receiver/source pair.
        The distance between each blade segment and receiver at each angular position is computed.
        Then the propagation time is computed with a reference wind speed of c0=343ms.

       Args:
       last (bool): Flag indicating if the last time is considered. Default is True.
       loop (bool): Flag indicating if the last time is set equal the first or not. Default is True.
        
        """
        print('compute receiver time ...')
        c0 = 343
        # (nx,ny,nz,src.wt.Nseg,src.wt.Nblade,len(freq),src.Nbeta)
        x = self.x_grid[:, 0, 0].reshape(-1, 1, 1, 1, 1, 1)
        y = self.y_grid[0, :, 0].reshape(1, -1, 1, 1, 1, 1)
        z = self.z_grid[0, 0, :].reshape(1, 1, -1, 1, 1, 1)

        self.wt.computeBeta()
        beta = self.wt.beta.reshape(1, 1, 1, 1, self.wt.Nblade, self.wt.Nbeta)
        if self.FULL_ROTATION:
            beta = np.concatenate((beta,
                                   beta+2*np.pi/3, beta+4*np.pi/3), axis=5)
        seg = self.wt.seg.reshape(1, 1, 1, self.wt.Nseg, 1, 1)

        # compute segment location
        # modif source :
        zS = (np.cos(beta) * seg + self.wt.href)
        xS = -np.sin(beta) * seg * np.sin(self.wt.tau) + self.xS
        yS = -np.sin(beta) * seg * np.cos(self.wt.tau) + self.yS

        # compute distance between source and receiver
        # zS = (np.cos(beta) * seg * 0 + self.wt.href)
        # xS = -np.sin(beta) * seg * np.sin(self.wt.tau) + self.xS
        # yS = -np.sin(beta) * seg * np.cos(self.wt.tau)  + self.yS

        R = np.sqrt((x - xS)**2 + (y - yS)**2 + (z - zS)**2)

        # absolute time of signal reception
        self.t = R/c0 + beta[:, :, :, :, 0, :].reshape(1, 1, 1,
                                                       1, 1, -1)/self.wt.omega

        if last:
            if loop:
                # for freq to time 
                self.t = np.concatenate((self.t,
                                 self.t[..., 0:1]), 5)
            else:
                # for receiver to time
                self.t = np.concatenate((self.t,
                                 self.t[..., 0:1]+2*np.pi/self.wt.omega), 5)
        print('done.')

    def angle_to_time_2(self, dt):
        """
        Converts the SPL as a function of the rotor angle (`SPL_seg`) to the SPL 
        as a function of receiver time (`SPL_time`).
        There is several method to do this transformation. 
        In this function a loop is done over all receiver positions. 
        For each ti a check is perfomed on each segment and angular position to check 
        if the receiver_time T computed satisfies ti>T and ti+1 > T and ti>ti+1.
        then the contribution of all segement is added logarythmically. 
        This is not the prefered version anymore. 

        Args:
            dt (float): The time step.
        """
        print('start converting angle to time ...')
        # define the time array according to max time and dt

        tmin = dt
        tmax = np.max(self.t)
        t_array = np.arange(tmin, tmax+dt, dt)

        t_array = np.arange(0, 2 * np.pi / (self.wt.omega), dt)

        print(t_array[-1])
        print(tmax)
        print(2*np.pi/(self.wt.omega))

        Nt = len(t_array)
        (nx, ny, nz, Nseg, Nblade, Nfreq, Nbeta) = self.SPL_seg.shape
        print('time array :' + str(Nt))

        # intialize output matrix
        self.SPL_time = np.zeros((nx, ny, nz, Nfreq, Nt))
        pp_time = np.zeros((nx, ny, nz, Nfreq, Nt))
        print('start loop over t ...')
        time0 = time()
        t_mask = 0
        t_sum = 0
        t_copy = 0
        print('compute pp_seg ...')
        ne.set_num_threads(8)
        pp_seg = ne.evaluate('10**(spl/10)', local_dict={'spl': self.SPL_seg})
        print('done.')
        # loop over T
        for ii, t in enumerate(t_array):
            print(t)
            t0 = time()
            for ix in range(self.t.shape[0]):
                for iy in range(self.t.shape[1]):
                    for iz in range(self.t.shape[2]):

                        t_loop = time()
                        # create masking array
                        # all values where ti >T and ti+1 > T and ti>ti+1
                        # the last condition is important: because of the
                        # modulo omega the t(beta) is not allways increasing
                        ti_m = self.t[ix, iy, iz, :, :, 0:-1]
                        ti_p = self.t[ix, iy, iz, :, :, 1:]

                        m1 = np.asarray(((ti_m <= t) & (ti_p >= t)
                                         & (ti_m < ti_p)), dtype=bool)

                        m2 = np.asarray(((
                            (((ti_m - 2*np.pi/(self.wt.omega)) <= t)
                             & (ti_p >= t)) |
                            ((ti_m <= t)
                             & ((ti_p + 2*np.pi / (self.wt.omega)) >= t)))
                            & (ti_m > ti_p)), dtype=bool)

                        mask = np.logical_or(m1, m2)
                        t1 = time()
                        t_mask += t1 - t_loop

                        mm = np.any(mask, (0, 1, 2))
                        if mm:
                            pp_time[ix, iy, iz, :, ii] = np.sum(
                                    pp_seg[ix, iy, iz, :, :, :, :],
                                    axis=(0, 1, 3),
                                    where=mask[:, :, None, :])
                        t2 = time()
                        t_sum += t2 - t1

            tend = time()
            print('create mask time = ' + str((t_mask)/(tend-t0)*100)+'%')
            print('sum time = ' + str((t_sum)/(tend-t0)*100)+'%')
            print('total time =' + str(tend - t0) + 's.')
            t_mask = 0
            t_sum = 0

        print('compute SPL ...')
        self.SPL_time = ne.evaluate('10*log10(pp)', local_dict={'pp': pp_time})
        # self.SPL_time = 10*np.log10(pp_time)
        self.time = t_array
        self.Nt = self.time.size
        print("done.")

    def compute_one_position_to_time(self, ix, iy, iz, pp_seg, pp_time):
        """
        Converts the SPL as a function of the rotor angle (`SPL_seg`) to the SPL 
        as a function of receiver time (`SPL_time`).
        There is several method to do this transformation.
        In this function a time signal is created from `pp_seg` for each blade at each angular position. 
        Each of this grain is shifted in time according to the `compute_receiver_time` function. 
        All grain are then added together to create the final `pp_time` signal. 
        In this method there is no overlap between grain. 

        Args:
            ix (int): The x-index.
            iy (int): The y-index.
            iz (int): The z-index.
            pp_seg (np.ndarray): The SPL segment data.
            pp_time (np.ndarray): The SPL time data.
        """
        # for iblade in [2]:
        for iblade in range(pp_seg.shape[1]):
            # for iblade in range(1):
            for iseg in range(pp_seg.shape[0]):
                # for iseg in [pp_seg.shape[0]-1]:
                t0 = self.t[ix, iy, iz, iseg, iblade, 0]
                i_first = int((t0 % self.T) * self.sampling)
                it = i_first
                # without overlap
                for ibeta in range(0, pp_seg.shape[3]):
                    tend = self.t[ix, iy, iz, iseg, iblade, ibeta+1]
                    iend = int((tend % self.T) * self.sampling)
                    if iend > it:
                        indices = np.arange(it, iend, 1)
                    else:
                        # print('ibeta loop = %s'%(ibeta))
                        indices = np.concatenate([np.arange(it, self.Nt, 1),
                                                  np.arange(0, iend, 1)],
                                                 axis=0)
                    # indices = np.arange(it,iend,1) % self.Nt
                    pp_time[:, indices] = pp_time[:, indices] + \
                            pp_seg[iseg, iblade, :, ibeta, None]
                    t0 = tend
                    it = iend
                    # plt.plot(pp_time[0,:])
                    # plt.show()

    def compute_one_position_to_time_overlap(self, ix, iy, iz, pp_seg, pp_time):
        """
        Converts the SPL as a function of the rotor angle (`SPL_seg`) to the SPL 
        as a function of receiver time (`SPL_time`).
        There is several method to do this transformation.
        In this function a time signal is created from `pp_seg` for each blade at each angular position. 
        Each of this grain is shifted in time according to the `compute_receiver_time` function. 
        All grain are then added together to create the final `pp_time` signal. 
        In this method there is an overlap between grain. 
        The number of sample where the overlap occurs is set by `self.overlap`.

        Args:
            ix (int): The x-index.
            iy (int): The y-index.
            iz (int): The z-index.
            pp_seg (np.ndarray): The SPL segment data.
            pp_time (np.ndarray): The SPL time output data.
        """
        for iblade in range(pp_seg.shape[1]):
            for iseg in range(pp_seg.shape[0]):
                t0 = self.t[ix, iy, iz, iseg, iblade, 0]
                i_first = int((t0 % self.T) * self.sampling)
                it = i_first
                # overlap
                for ibeta in range(0, pp_seg.shape[3]):
                    # deltaT = (self.t[ix, iy, iz, iseg, iblade, ibeta+1] - t0)
                    tend = self.t[ix, iy, iz, iseg, iblade, ibeta+1]
                    iend = int((tend % self.T) * self.sampling)
                    i0 = (it - self.overlap) % self.Nt
                    i1 = (iend + self.overlap) % self.Nt
                    # print(i0)
                    # print(i1)
                    if i1 > i0:
                        indices = np.arange(i0, i1, 1)
                    else:
                        # print('ibeta loop = %s'%(ibeta))
                        indices = np.concatenate([np.arange(i0, self.Nt, 1),
                                                  np.arange(0, i1, 1)],
                                                 axis=0)
                    # indices = np.arange(it,iend,1) % self.Nt

                    window = cos_window(len(indices), 2*self.overlap)

                    pp_time[:, indices] = pp_time[:, indices]
                    + window * pp_seg[iseg, iblade, :, ibeta, None]
                    t0 = tend
                    it = iend
                    # plt.plot(pp_time[0,:])
                    # plt.show()

    def angle_to_time_3(self,x:float,y:float,z:float,dt:float=0.1):
        """
        Converts the SPL as a function of the rotor angle (`SPL_seg`) to the SPL 
        as a function of receiver time (`SPL_time`) for a given position without overlap.
        This uses the function `compute_one_position_to_time`.

        Args:
            x (float): The x-coordinate.
            y (float): The y-coordinate.
            z (float): The z-coordinate.
            dt (float, optional): The time step. Default is 0.1.
        """
        # convert SPL(theta) tp spl(t) for a fiven position
        # this should be a better way than the one coded in angle_to_time_2
        ix = np.argmin(np.abs(self.x_grid[:, 0, 0] - x))
        iy = np.argmin(np.abs(self.y_grid[0, :, 0] - y))
        iz = np.argmin(np.abs(self.z_grid[0, :, 0] - z))
        # ix = 493
        # ix = 503
        # ix=504
        # ix = 531
        self.sampling = 1/dt
        pp_seg = (10**(self.SPL_seg[ix, iy, iz, :, :, :, :]/10))
        # t = self.t[ix, iy, iz, :, :, :]
        self.T = 2*np.pi/(self.wt.omega)
        self.Nt = int(self.T*self.sampling)

        pp_time = np.zeros((1, 1, 1, self.Nfreq, self.Nt))
        self.time = np.linspace(0, self.T, self.Nt)
        start_time = time()

        self.compute_one_position_to_time(ix, iy, iz, pp_seg,
                                          pp_time[0, 0, 0, ...])

        end_time = time()
        print(end_time - start_time)
        plt.plot(pp_time[0, 0, 0, 0, :])
        plt.show()

        # self.SPL_time = np.roll(10*np.log10(pp_time),offset,axis=4)
        self.SPL_time = 10*np.log10(pp_time)

    def angle_to_time_3_full(self, dt: float = 0.1, overlap=0):
        """
        Converts the SPL as a function of the rotor angle (`SPL_seg`) to the SPL 
        as a function of receiver time (`SPL_time`) for the entire domain.
        This uses the function `compute_one_position_to_time` if `overlap` is set to 0
        and `compute_one_position_to_time_overlap` if `overlap>0`.

        Args:
            dt (float, optional): The time step. Default is 0.1.
            overlap (int, optional): The overlap. Default is 0.
        """
        self.overlap = overlap
        self.sampling = 1/dt

        if self.SPL_seg is not None:
            pp_seg = (10**(self.SPL_seg/10))
            logging.info("converting spl(f) into time")
            nfreq_tmp = self.Nfreq
        elif self.OASPL_seg is not None:
            pp_seg = (10**(self.OASPL_seg/10))[:, :, :, :, :, None, :]
            logging.info("converting oaspl into time")
            nfreq_tmp = 1
        else:
            print("Error: SPL_seg or OASPL_seg mus not None")
            return -1

        t = self.t
        self.T = 2 * np.pi / (self.wt.omega)
        self.Nt = int(self.T * self.sampling)
        self.time = np.arange(0, self.T, dt)
        self.T = self.time[-1] + dt

        print('Nt = ')
        print(self.Nt)
        # copy first angle at the end
        self.t = np.concatenate((self.t, self.t[..., 0:1] + self.T), 5)
        self.t = self.t - np.min(self.t)
        pp_time = np.zeros((self.Nx, self.Ny, self.Nz, nfreq_tmp, self.Nt))

        t0 = time()
        logging.info("starting loop on receiver position...")
        for ix in tqdm(range(self.Nx)):
            for iy in range(self.Ny):
                for iz in range(self.Nz):
                    # print(ix,iy,iz)
                    if overlap == 0:
                        self.compute_one_position_to_time(ix, iy, iz,
                                                          pp_seg[ix, iy,
                                                                 iz, ...],
                                                          pp_time[ix, iy,
                                                                  iz, ...])
                    else:
                        self.compute_one_position_to_time_overlap(ix, iy, iz,
                                                                pp_seg[ix, iy,
                                                             iz, ...],
                                                      pp_time[ix, iy,
                                                              iz, ...])

        t1 = time()
        print(t1-t0)
        # plt.plot(pp_time[0,0,0,0,:])
        # plt.show()

        self.time = np.arange(0, self.T-dt, dt)
        # self.SPL_time = np.roll(10*np.log10(pp_time),offset,axis=4)
        if self.SPL_seg is not None:
            self.SPL_time = 10*np.log10(pp_time)
        elif self.OASPL_seg is not None:
            self.OASPL_time = 10*np.log10(pp_time)[:, :, :, 0, :]


    def save(self,fname: str):
        """
        Save data in a pickle .dat file. B
        Big matrix are store in h5 files in order to open only part of the data.
        """
        print('saving SPL ...')
        self.src = None
        self.deltaL = None

        with h5py.File(fname[:-3]+'h5', "w") as f:
            if self.SPL_time is not None:
                f.create_dataset('SPL_time', data=self.SPL_time)
            if self.SPL_seg is not None:
                f.create_dataset('SPL_seg', data=self.SPL_seg)
            if self.OASPL_time is not None:
                f.create_dataset('OASPL_time', data=self.OASPL_time)
            if self.OASPL_seg is not None:
                f.create_dataset('OASPL_seg', data=self.OASPL_seg)
            if self.OASPL_seg_tot is not None:
                f.create_dataset('OASPL_seg_tot', data=self.OASPL_seg_tot)

            if self.am_seg is not None:
                f.create_dataset('am_seg', data=self.am_seg)
            if self.am_time is not None:
                f.create_dataset('am_time', data=self.am_time)
            if self.OAam_seg is not None:
                f.create_dataset('OAam_seg', data=self.OAam_seg)
            if self.OAam_time is not None:
                f.create_dataset('OAam_time', data=self.OAam_time)

            if self.mean_seg is not None:
                f.create_dataset('mean_seg', data=self.mean_seg)
            if self.mean_time is not None:
                f.create_dataset('mean_time', data=self.mean_time)
            if self.OAmean_seg is not None:
                f.create_dataset('OAmean_seg', data=self.OAmean_seg)
            if self.OAmean_time is not None:
                f.create_dataset('OAmean_time', data=self.OAmean_time)

        self.SPL_time = None
        self.SPL_seg = None
        self.OASPL_time = None
        self.OASPL_seg = None
        self.OASPL_seg_tot = None

        self.am_seg = None
        self.am_time = None
        self.OAam_seg = None
        self.OAam_time = None

        self.mean_seg = None
        self.mean_time = None
        self.OAmean_seg = None
        self.OAmean_time = None

        """save class as self.name.dat"""
        with open(fname,'wb') as file:
            pickle.dump(self.__dict__,file)
        print('done.')

    def load(self, fname: str, seg: bool = True, time: bool = True, oaspl: bool = True, am: bool = True,
             mean: bool = True, x: float = None, y: float = None, z: float = None, freq: float = None, Nt: int = None):
        """
        Loads the SPL field data from a file.
        Flag are used to chose which data to load.
        WARNING: the name of the flag is very poorly chosen. Should be change at some point

        Args:
            fname (str): The file name to load the data.
            seg (bool, optional): Flag indicating if the frequency data is loaded. Default is True.
            time (bool, optional): Flag indicating if the time data is loaded. If false the data as function of beta is loaded. Default is True.
            oaspl (bool, optional): Flag indicating if the OASPL data is loaded. Default is True.
            am (bool, optional): Flag indicating if the amplitude modulation data is loaded. Default is True.
            mean (bool, optional): Flag indicating if the time (or beta) averaged data is loaded. Default is True.
            x (float, optional): The x-coordinate to load the data if None load all $x$ positions are loaded. Default is None.
            y (float, optional): The y-coordinate to load the data if None load all $y$ positions are loaded. Default is None.
            z (float, optional): The z-coordinate to load the data if None load all $z$ positions are loaded. Default is None.
            freq (float, optional): The frequency to load the data if None load all frequencies are loaded. Default is None.
            Nt (int, optional): The number of time steps to load the data. Default is None.
        """
        logging.info('loading SPL ...')
        with open(fname,'rb') as file:
            self.__dict__ = pickle.load(file)
        logging.info('done loading pickle')

        if x is None:
            ix = slice(0,self.Nx,1)
        else:
            ix = np.argmin(np.abs(self.x_grid[:,0,0]-x))
            self.Nx = 1
            self.x_grid = self.x_grid[ix:ix+1,:,:]
            self.y_grid = self.y_grid[ix:ix+1,:,:]
            self.z_grid = self.z_grid[ix:ix+1,:,:]
        if y is None:
            iy = slice(0,self.Ny,1)
        else:
            iy= np.argmin(np.abs(self.y_grid[0,:,0]-y))
            self.Ny = 1
            self.x_grid = self.x_grid[:,iy:iy+1,:]
            self.y_grid = self.y_grid[:,iy:iy+1,:]
            self.z_grid = self.z_grid[:,iy:iy+1,:]
        if z is None:
            iz = slice(0,self.Nz,1)
        else:
            iz = np.argmin(np.abs(self.z_grid[0,0,:]-z))
            self.Nz = 1
            self.x_grid = self.x_grid[:,:,iz:iz+1]
            self.y_grid = self.y_grid[:,:,iz:iz+1]
            self.z_grid = self.z_grid[:,:,iz:iz+1]
        if freq is None:
            ifreq = slice(0,self.Nfreq,1)
            nfreq = self.Nfreq
        else:
            ifreq = np.argmin(np.abs(self.frequencies-freq))
            self.Nfreq = 1
            self.frequencies = self.frequencies[ifreq:ifreq+1]

        if Nt is not None:
            if Nt>self.Nt:
                print('error new Nt to large')
                return
            self.Nt = Nt
            self.time = self.time[:self.Nt]



        with h5py.File(fname[:-3]+'h5', "r") as f:
            if seg and not(time):
                if 'SPL_seg' in f:
                    logging.info('loading SPL_seg')
                    self.SPL_seg = np.array(f["SPL_seg"][ix,iy,iz,:,:,ifreq,:]).reshape((self.Nx,self.Ny,self.Nz,
                                                                                            self.wt.Nseg,self.wt.Nblade,self.Nfreq,
                                                                                            int(self.Nbeta)))
            if seg and time:
                if 'SPL_time' in f:
                    logging.info('loading SPL_time')
                    self.SPL_time = np.array(f.get('SPL_time')[ix,iy,iz,ifreq,:self.Nt]).reshape((self.Nx,self.Ny,self.Nz,
                                                                                        self.Nfreq,self.Nt))

            if oaspl and not(time):
                if 'OASPL_seg' in f:
                    self.OASPL_seg = np.array(f.get('OASPL_seg')[ix,iy,iz,:,:,:self.Nbeta])
                    logging.info('loading OASPL_seg')

                    self.OASPL_seg = np.array(f.get('OASPL_seg')[ix,iy,iz,
                                                                 :,:,:self.Nbeta]).reshape((self.Nx,
                                                                self.Ny,self.Nz,
                                                                self.wt.Nseg,self.wt.Nblade,
                                                                self.Nbeta))

                if 'OASPL_seg_tot' in f:
                    logging.info('loading OASPL_seg_tot')
                    self.OASPL_seg_tot = np.array(f.get('OASPL_seg_tot')[ix,iy,iz,
                                                                 :self.Nbeta]).reshape((self.Nx,
                                                                self.Ny,self.Nz,
                                                             self.Nbeta))
            if oaspl and time:
                if 'OASPL_time' in f:
                    logging.info('loading SPL_time')
                    self.OASPL_time = np.array(f.get('OASPL_time')[ix,iy,iz,:self.Nt]).reshape((self.Nx,self.Ny,self.Nz,self.Nt))

            if am:
                if 'am_time' in f:
                    logging.info('loading am_tim')
                    self.am_time = np.array(f.get('am_time')[ix,iy,iz,ifreq]).reshape((self.Nx,self.Ny,self.Nz,self.Nfreq))
                if 'am_seg' in f:
                    logging.info('loading am_seg')
                    self.am_seg = np.array(f.get('am_seg')[ix,iy,iz,ifreq]).reshape((self.Nx,self.Ny,self.Nz,self.Nfreq))
                if 'OAam_time' in f:
                    logging.info('loading OAam_time')
                    self.OAam_time = np.array(f.get('OAam_time')[ix,iy,iz]).reshape((self.Nx,self.Ny,self.Nz))
                if 'OAam_seg' in f:
                    logging.info('loading OAam_seg')
                    self.OAam_seg = np.array(f.get('OAam_seg')[ix,iy,iz]).reshape((self.Nx,self.Ny,self.Nz))

            if mean:
                if 'mean_time' in f:
                    self.mean_time = np.array(f.get('mean_time')[ix,iy,iz,ifreq]).reshape((self.Nx,self.Ny,self.Nz,self.Nfreq))
                if 'mean_seg' in f:
                    self.mean_seg = np.array(f.get('mean_seg')[ix,iy,iz,ifreq]).reshape((self.Nx,self.Ny,self.Nz,self.Nfreq))
                if 'OAmean_time' in f:
                    self.OAmean_time = np.array(f.get('OAmean_time')[ix,iy,iz]).reshape((self.Nx,self.Ny,self.Nz))
                if 'OAmean_seg' in f:
                    self.OAmean_seg = np.array(f.get('OAmean_seg')[ix,iy,iz]).reshape((self.Nx,self.Ny,self.Nz))
                        
        print('done.')


    def shift(self, xshift: float, yshift: float):
        """
        Shifts the $xy$ grid by the given x and y shifts.

        Args:
            xshift (float): The x-shift.
            yshift (float): The y-shift.
        """
        self.x_grid = self.x_grid + xshift
        self.y_grid = self.y_grid + yshift

    def atm_absorption(self,c0: float = 343,rho0: float = 1.2,rh: float = 80):
        """
        Applies atmospheric absorption due to the thermo-viscous effects
        and relaxation of oxygen and nitrogen to the SPL field.

        alpha - attenuation of sound for input parameters in dB/m

        T0 - temperature in K

        p0 - static pressure in pascal

        rh - relative humidity

        f - frequency of sound (may be a vector)

        The attenuation can only be applied to `SPL_time` or `SPL_seg` and not to frequency integrated data.

        References:   Salomons p.109-111
        Args:
            c0 (float, optional): The speed of sound. Default is 343.
            rho0 (float, optional): The density of air. Default is 1.2.
            rh (float, optional): The relative humidity. Default is 80.
        """
        print('compute atmospheric absorption')
        if self.ATM_ABS:
            logging.warning("atmospheric absorption already applied")
            return
        if not self.third:
            logging.warning("results are not int third octave band")

        rGP = 287.06;
        gamma = 1.4# constante adiabatique
        T0 = c0**2/(gamma*rGP)# temperature
        p0 = rho0*c0**2/gamma# pression atmospherique

        p0_ref = 1.01325e+05 # reference static pressure (pa)

        T0_tpw = 273.15 # triple point in K
        T0_ref = 293.15 # ref temp in K

        rho = p0/p0_ref
        tau = T0/T0_ref

        # calculate saturation pressure
        Csat = -6.8346*(T0_tpw/T0)**1.261 + 4.6151
        p0_sat = p0_ref*10**Csat
        h = rh*p0_sat/p0 # absolute humidity

        # Scaled relaxation frequency for Nitrogen
        frN = rho*tau**(-1/2)*(9 + 280*h*np.exp(-4.17*(tau**(-1/3)-1)))

        # scaled relaxation frequency for Oxygen
        frO = rho*(24 + 40400*h*(0.02+h)/(0.391+h))

        # attenuation coefficient in dB/m
        b1 = 0.1068*np.exp(-3352/T0)/(frN + self.frequencies**2/frN)
        b2 = 0.01275*np.exp(-2239.1/T0)/(frO + self.frequencies**2/frO)
        alpha = 8.686*self.frequencies**2*tau**(1/2)*(1.84e-11/rho + tau**(-3)*(b1 + b2))

        R = np.sqrt((self.x_grid - self.xS)**2 + (self.y_grid - self.yS)**2 + (self.z_grid - self.wt.href)**2)
        if self.SPL_seg is not None:
            self.SPL_seg = self.SPL_seg - alpha.reshape(1,1,1,1,1,-1,1)*R.reshape(self.Nx,self.Ny,self.Nz,1,1,1,1)
        if self.SPL_time is not None:
            self.SPL_time = self.SPL_time - alpha.reshape(1,1,1,-1,1)*R.reshape(self.Nx,self.Ny,self.Nz,1,1)
        self.ATM_ABS = True

    def Aweight(self):
        """
        Applies A-weighting to the SPL field.
        The weighting can only be applied to `SPL_time` or `SPL_seg` and not to frequency integrated data.
        """
        if self.AWEIGHT:
            logging.warning("A weighting already applied")
            return
        if not self.third:
            logging.warning("results are not int third octave band")
        print('compute Aweight')
        Af = 12200**2*self.frequencies**4./(self.frequencies**2+20.6**2)/(self.frequencies**2+12200**2)/(self.frequencies**2+107.7**2)**0.5/(self.frequencies**2+737.9**2)**0.5
        dBA = 20*np.log10(Af/0.7943)
        if self.SPL_seg is not None:
            self.SPL_seg = self.SPL_seg + dBA.reshape(1,1,1,1,1,-1,1)
        if self.SPL_time is not None:
            self.SPL_time = self.SPL_time + dBA.reshape(1,1,1,-1,1)
        self.AWEIGHT = True
        return dBA

    def combine_2_turbines(self, spl2, shift: int = 0, tmax: float = None):
        """
        Combines the SPL fields of two turbines.
        The function first tries to combine `SPL_time`, then `OASPL_time` then `OASPL_seg` 
        if the array are loaded. 
        If the two signals dont have the same length the signal are looped until reaching a same value. 
        If a tmax is set this loop is cut when it reaches tmax. 

        The two time signal can also be shifted with an offset added to `spl2`. 

        Args:
            spl2 (SplField): The second SPL field object.
            shift (int, optional): The shift between the two turbines. Default is 0.
            tmax (float, optional): The maximum time. Default is None.

        Returns:
            int: Returns -1 if the combination is not possible.
        """
        seg_flag = False
        if self.SPL_time is not None:
            s1 = self.SPL_time
            if spl2.SPL_time is not None:
                s2 = spl2.SPL_time
                spl_flag = 1
            else:
                print("error spl2 has no SPL_time")
                return -1
        elif self.OASPL_time is not None: 
            s1 = self.OASPL_time 
            if spl2.OASPL_time is not None:
                s2 = spl2.OASPL_time
                spl_flag = 0
            else:
                print("error spl2 has no OASPL_time")
                return -1
        elif self.OASPL_seg is not None:

            s1 = self.OASPL_seg_tot
            if spl2.OASPL_seg_tot is not None:
                s2 = spl2.OASPL_seg_tot
                seg_flag = True
                spl_flag = 0
            print("combine turbines for OASPL_seg_tot")
        else:
            print("Error: no OASPL_time or SPL_time loaded")

        dim = len(s1.shape)-1
        print(s1.shape)
        if (s1.shape[-1] == s2.shape[-1]):
            print('signal same size, periode unchanged')
            spl_tot = 10*np.log10(
                    10**(s1/10)
                    + 10**(np.roll(s2,shift, axis=dim)/10))
            if seg_flag:
                Nt = self.Nbeta
            else:
                Nt = self.Nt
        else:
            Nt = np.lcm(self.Nt, spl2.Nt)
            if tmax is not None:
                Nt = min(Nt, int(tmax * self.sampling))
            print('signal of different length, periode set to :')
            print(Nt)
            spl_tot = 10*np.log10(10**(uneven_loop(s1, Nt)/10) +
                                  10**(uneven_loop(s2, Nt)/10))

            print('done.')

        if seg_flag:
            self.OASPL_seg_tot = spl_tot
        elif spl_flag:
            self.SPL_time = spl_tot
        else:
            self.OASPL_time = spl_tot
        if not seg_flag:
            self.Nt = Nt
            dt = self.time[1]-self.time[0]
            self.time = np.arange(0, self.Nt*dt, dt)


        # if (self.SPL_time.shape[4] == spl2.SPL_time.shape[4]):
        #     print('signal same size, periode unchanged')
        #     spl_tot = 10*np.log10(
        #             10**(self.SPL_time[:,:,:,:,:]/10) 
        #             + 10**(np.roll(spl2.SPL_time[:,:,:,:,:],shift,axis=4)/10))
        #     self.SPL_time[:,:,:,:,:] = spl_tot
        # else:
        #     Nt = np.lcm(self.Nt, spl2.Nt)
        #     if tmax is not None:
        #         Nt = min(Nt, int(tmax * self.sampling))
        #     print('signal of different length, periode set to :')
        #     print(Nt)
        #     # spl_tot = 10*np.log10(np.tile(10**(self.SPL_time/10),(1,1,1,1,int(Nt/self.Nt))) +
        #     #                      np.tile(10**(spl2.SPL_time/10),(1,1,1,1,int(Nt/spl2.Nt))))

        #     spl_tot = 10*np.log10(uneven_loop(10**(self.SPL_time/10), Nt) +
        #                           uneven_loop(10**(spl2.SPL_time/10), Nt))

        #     # spl_tot = 10*np.log10(uneven_tile(10**(self.SPL_time/10),(1,1,1,1,Nt/self.Nt)) +
        #     #                       uneven_tile(10**(spl2.SPL_time/10),(1,1,1,1,Nt/spl2.Nt)))
        #     print('done.')
        #     self.SPL_time = spl_tot
        #     self.Nt = Nt
        #     dt = self.time[1]-self.time[0]
        #     self.time = np.arange(0, self.Nt*dt, dt)


    def compute_third_octave(self, fc: no.ndarray = None, Nfc: np.ndarray = None):
        """
        Compute third octave band spectrum from previously computed frequencies.
        Be sure that fc anf Nfc corresponds to the `self.frequencies` computed. 
        If fc or Nfc set to None a default value are set (corresponding to frequency compute in Colas et al. (2023)).
        
        Modifies:

        - `SPL_time`  of shape (Nx, Ny, Nz, Nfreq, Nt) with new Nfreq
        - `SPL_seg` of shape (Nx, Ny, Nz, Nseg, Nblade, Nfreq, Nbeta) with new Nfreq

        Args:
            fc (np.ndarray): central frequency of each band. Default to None.
            Nfc (np.ndarray): number of frequency computed per band. default to None.
        """
        if (fc is None) or (Nfc is None):
            print("set to default frequency bands")
            fc = [50, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000]
            Nfc = [1,  1,  1,   1,   1,  1,   2,   2,   3,   4,   4,   4,   5,  5]
        print('computing third octave ...')
        if self.third:
            logging.warning('skipping function: results are already in 3rd octave band.')
            return
        if self.SPL_time is not None:
            SPL_time_third = np.zeros(np.shape(self.SPL_time[...,0,0]) + (len(fc),self.Nt))
        if self.SPL_seg is not None:
            SPL_seg_third =  np.zeros(np.shape(self.SPL_seg[...,0,0]) + (len(fc),self.Nbeta))
        for ifreq in range(len(Nfc)):
            # temp_array = 10*np.log10(fc[ifreq]*0.232/Nfc[ifreq]* np.sum(10 ** (self.SPL_time[..., int(np.sum(Nfc[0:ifreq])):int(np.sum(Nfc[0:ifreq+1])),:]/10), axis=3))
            # temp_array[temp_array == -inf] = 0
            if self.SPL_time is not None:
                SPL_time_third[..., ifreq,:] = 10*np.log10(fc[ifreq]*0.232/Nfc[ifreq]* np.sum(10 ** (self.SPL_time[..., int(np.sum(Nfc[0:ifreq])):int(np.sum(Nfc[0:ifreq+1])),:]/10), axis=3))
            if self.SPL_seg is not None:
                SPL_seg_third[..., ifreq,:] = 10*np.log10(fc[ifreq]*0.232/Nfc[ifreq]* np.sum(10 ** (self.SPL_seg[..., int(np.sum(Nfc[0:ifreq])):int(np.sum(Nfc[0:ifreq+1])),:]/10), axis=5))

        if self.SPL_time is not None:
            self.SPL_time = SPL_time_third
        if self.SPL_seg is not None:
            self.SPL_seg = SPL_seg_third

        self.frequencies = np.array(fc)
        self.Nfreq = len(fc)
        self.third = True
        print('done.')

    def compute_oaspl(self):
        """
        Compute OASPL from previously computed thrid octave band spectrum. 
        
        Computes:

        - `OASPL_time`  of shape (Nx, Ny, Nz,Nt)
        - `OASPL_seg` of shape (Nx, Ny, Nz, Nseg, Nblade, Nbeta)
        - `OASPL_seg_tot` of shape (Nx, Ny, Nz, Nbeta)

        """
        print('computing OASPL ...')
        if self.third is False:
            logging.warning('need to convert to third octave band spectra first.')
            return -1
        if not self.AWEIGHT:
            logging.warning("computed OASPL without A weighting")

        if not self.ATM_ABS:
            logging.warning("computed OASPL without atmospheric abs")

        if self.SPL_time is not None:
            self.OASPL_time = 10*np.log10(np.sum(10**(self.SPL_time/10), 3))
        if self.SPL_seg is not None:
            self.OASPL_seg_tot = 10*np.log10(np.sum(10**(self.SPL_seg/10), (3,4,5)))
            self.OASPL_seg = 10*np.log10(np.sum(10**(self.SPL_seg/10), (5)))
        self.oaspl = True
        print('done.')

    def compute_am(self):
        """
        Compute AM from previoulsy computed spectrum or OASPL.

        Computes:

        - `am_time`  of shape (Nx, Ny, Nz, Nfreq)
        - `am_seg` of shape (Nx, Ny, Nz, Nfreq)
        - `OAam_time` of shape (Nx, Ny, Nz)
        - `OAam_seg` of shape (Nx, Ny, Nz)
        """
        print('Computing AM ...')
        if self.SPL_time is not None:
            self.am_time = np.max(self.SPL_time,4)-np.min(self.SPL_time,4)
        if self.SPL_seg is not None:
            print(self.SPL_seg.shape)
            spl = 10*np.log10(np.sum(10**(self.SPL_seg/10),(3,4)))
            print(spl.shape)
            self.am_seg = np.max(spl,4)-np.min(spl,4)
            print(self.am_seg.shape)
        if self.OASPL_time is not None:
            self.OAam_time = np.max(self.OASPL_time,3)-np.min(self.OASPL_time,3)
        if self.OASPL_seg_tot is not None:
            self.OAam_seg = np.max(self.OASPL_seg_tot,3)-np.min(self.OASPL_seg_tot,3)
        print('done.')

    def compute_time_average(self):
        """
        Compute AM from previoulsy computed spectrum or OASPL.

        Computes:

        - `mean_time`  of shape (Nx, Ny, Nz, Nfreq)
        - `mean_seg` of shape (Nx, Ny, Nz, Nfreq)
        - `OAmean_time` of shape (Nx, Ny, Nz)
        - `OAmean_seg` of shape (Nx, Ny, Nz)
        """
        print('computing time/beta average ...')
        if self.SPL_time is not None:
            self.mean_time = np.mean(self.SPL_time,4)
        if self.SPL_seg is not None:
            self.mean_seg = 10*np.log10(np.sum(10**(self.SPL_seg/10),(3,4,6))/self.Nbeta)
            # self.mean_seg = np.mean(spl,4)
        if self.OASPL_time is not None:
            self.OAmean_time = 10*np.log10(np.mean(10**(self.OASPL_time/10),3))
        if self.OASPL_seg is not None:
            self.OAmean_seg = np.mean(self.OASPL_seg,3)
        if self.OASPL_seg_tot is not None:
            self.OAmean_seg = 10*np.log10(np.mean(10**(self.OASPL_seg_tot/10),3))
        print('done.')

    # Plot 
    #------------------------------------------------------------------------------------------------------
    def plot_spl(self, time=True, OA=False,
                 x: float = None, y: float = None, z: float = None,
                 freq: float = None, it: int = None,
                 roll=0, angle=False, **kwargs):

        """
        Plots the SPL field. 
        if `x`, `y`, `z`, `freq`, `it` is set to `None` the plot is done in this direction. 
        Plots can be 1D (SPL$(x)$ or SPL$(it)$ for example) or 2D (SPL$(x,y)$ or SPL$(freq,it)$ for example) 

        Args:
            time (bool, optional): Flag indicating if time data (True) or angular data (False) is plotted. Default is True.
            OA (bool, optional): Flag indicating if the OASPL data (True) or frequency data (False) is plotted. Default is False.
            x (float, optional): The x-coordinate to plot the data. Default is None.
            y (float, optional): The y-coordinate to plot the data. Default is None.
            z (float, optional): The z-coordinate to plot the data. Default is None.
            freq (float, optional): The frequency to plot the data. Default is None.
            it (int, optional): The time index to plot the data. Default is None.
            roll (int, optional): Allow to shift the SPL in time domain by several indices. Default is 0.
            angle (bool, optional): Flag indicating if the angle data is plotted. Default is False.
            **kwargs: Additional keyword arguments for plotting.
        """

        # choose Over all data or frequency data
        # choose between beta data or time dat
        test_dim = [x is None, y is None, z is None, freq is None, it is None]
        beta = self.wt.beta.reshape(1, 1, 1, 1, self.wt.Nblade, self.wt.Nbeta)
        self.time_src = beta[:, :, :, :, 0, :].reshape(-1)/self.wt.omega

        if time:
            absice = [self.x_grid[:, 0, 0], self.y_grid[0, :, 0],
                      self.z_grid[0, 0, :], self.frequencies, self.time]
            labels = ['x (m)', 'y (m)', 'z (m)', 'f (Hz)', 't (s)']
            time_src = self.time
        else:
            test_dim = [x is None, y is None, z is None, freq is None, it is None]
            beta = self.wt.beta.reshape(1, 1, 1, 1, self.wt.Nblade, self.wt.Nbeta)
            beta = np.concatenate((beta, beta+2*np.pi/3, beta+4*np.pi/3), axis=5)
            if angle:
                time_src = beta[:, :, :, :, 0, :].reshape(-1)*180/np.pi
            else:
                time_src = beta[:, :, :, :, 0, :].reshape(-1)/self.wt.omega

            absice = [self.x_grid[:, 0, 0],
                      self.y_grid[0, :, 0], self.z_grid[0, 0, :],
                      self.frequencies, time_src]
            if angle:
                labels = ['x (m)', 'y (m)', 'z (m)', 'f (Hz)', '$\\beta$ (deg)']
            else:
                labels = ['x (m)', 'y (m)', 'z (m)', 'f (Hz)', '$t$ (s)']
        if OA:
            test_dim[3] = False
            if time:
                # spl = np.roll(self.OASPL_time, roll,axis=-1)
                spl = self.OASPL_time
            else:
                spl = self.OASPL_seg_tot
        else:
            if time:
                # spl = np.roll(self.SPL_time,roll,axis=-1)
                spl = self.SPL_time
            else:
                spl = 10*np.log10(np.sum(10**(self.SPL_seg/10),(3,4)))

        # find coresponding freq, ix, iy, iz 
        if freq is not None:
            ifreq = np.argmin(np.abs(self.frequencies-freq))
            spl = spl[...,ifreq:ifreq+1,:]
        if x is not None:
            ix = np.argmin(np.abs(self.x_grid[:,0,0]-x))
            print(ix)
            spl = spl[ix:ix+1,...]
        if y is not None :
            iy = np.argmin(np.abs(self.y_grid[0,:,0]-y))
            print(iy)
            spl = spl[:,iy:iy+1,...]
        if z is not None :
            iz = np.nonzero(self.z_grid[0,0,:]==z)[0][0]
            print(iz)
            spl = spl[:,:,iz:iz+1,...]
        if it is not None :
            spl = spl[...,it]
        spl = np.squeeze(spl)

        # plot either pcolor or lineplot 
        if len(spl.shape) == 1:
            print("plotting 1D spl")
            plt.plot(absice[np.where(test_dim)[0][0]], spl,**kwargs)
            # plt.plot(spl,**kwargs)
            plt.xlabel(labels[np.where(test_dim)[0][0]])
            plt.ylabel('SPL (dB)')

        if len(spl.shape) == 2:
            print("plotting 2D spl")
            # if z is not None:
            # if test dim == [1,1,0,0,0]:
            #     absice = [self.x_grid[:, :, 0],
            #           self.y_grid[:, :, 0], self.z_grid[0, 0, :],
            #           self.frequencies,time_src]
            # # if y is not None:
            # if test dim == [1,0,1,0,0]:
            #     absice = [self.x_grid[:, 0, :],
            #               self.y_grid[:, 0, :], self.z_grid[:, 0, :],
            #           self.frequencies,time_src]
            # # if x is not None:
            # if test dim == [0,1,1,0,0]:
            #     absice = [self.x_grid[0, :, :],
            #               self.y_grid[0, :, :], self.z_grid[0, :, :],
            #           self.frequencies,time_src]
            #absice = [self.x_grid[:, :, 0],
            #          self.y_grid[:, :, 0], self.z_grid[0, 0, :],
            #          self.frequencies, time_src]
            #plt.pcolormesh(absice[np.where(test_dim)[0][0]],absice[np.where(test_dim)[0][1]],spl,**kwargs)
            
            
            plt.pcolormesh(absice[np.where(test_dim)[0][0]],absice[np.where(test_dim)[0][1]],spl.T,**kwargs)
            plt.xlabel(labels[np.where(test_dim)[0][0]])
            plt.ylabel(labels[np.where(test_dim)[0][1]])
        if len(spl.shape) == 3:
            return

    def plot_mean(self, time: bool = True, OA: bool = False, x: float = None, y: float = None, z: float = None, freq: float = None, **kwargs):
        """
        Plots the time averaged SPL field. 
        if `x`, `y`, `z`, `freq` is set to `None` the plot is done in this direction. 
        Plots can be 1D (SPL$(x)$ or SPL$(it)$ for example) or 2D (SPL$(x,y)$ or SPL$(freq,it)$ for example) 

        Args:
            time (bool, optional): Flag indicating if time data (True) or angular data (False) is plotted. Default is True.
            OA (bool, optional): Flag indicating if the OASPL data (True) or frequency data (False) is plotted. Default is False.
            x (float, optional): The x-coordinate to plot the data. Default is None.
            y (float, optional): The y-coordinate to plot the data. Default is None.
            z (float, optional): The z-coordinate to plot the data. Default is None.
            freq (float, optional): The frequency to plot the data. Default is None.
            **kwargs: Additional keyword arguments for plotting.
        """
        # choose Over all data or frequency data 
        # choose between beta data or time dat 
        test_dim = [x is None, y is None, z is None, freq is None] 
        absice = [self.x_grid[:,0,0],self.y_grid[0,:,0],self.z_grid[0,0,:],self.frequencies]
        labels = ['x (m)', 'y (m)', 'z (m)', 'f (Hz)']
        if OA:
            test_dim[3] = False
            if time:
                mean = self.OAmean_time
            else:
                mean = self.OAmean_seg
        else:
            if time:
                mean = self.mean_time
            else:
                mean = self.mean_seg

        # find coresponding freq, ix, iy, iz 
        if freq is not None:
            ifreq = np.argmin(np.abs(self.frequencies-freq))
            mean = mean[...,ifreq:ifreq+1]
        if x is not None:
            ix = np.argmin(np.abs(self.x_grid[:,0,0]-x))
            mean = mean[ix:ix+1,...]
        if y is not None :
            iy = np.argmin(np.abs(self.y_grid[0,:,0]-y))
            mean = mean[:,iy:iy+1,...]
        if z is not None :
            iz = np.nonzero(self.z_grid[0,0,:]==z)[0][0]
            mean = mean[:,:,iz:iz+1,...]
        mean = np.squeeze(mean)

        # plot either pcolor or lineplot 
        if len(mean.shape)==1:
            plt.plot(absice[np.where(test_dim)[0][0]], mean,**kwargs)
            plt.xlabel(labels[np.where(test_dim)[0][0]])
            plt.ylabel('OASPL(dB)')

            # if x is  None:
            #     plt.plot(self.x_grid[:,0,0],mean,**kwargs)
            # if y is  None:
            #     plt.plot(self.y_grid[0,:,0],mean,**kwargs)
            # if z is  None:
            #     plt.plot(self.z_grid[0,0,:],mean,**kwargs)
            # if freq is None:
            #     plt.plot(self.frequencies,mean,**kwargs)

        if len(mean.shape)==2:
            # if z is not None:
            #     absice = [self.x_grid[:, :, 0],
            #           self.y_grid[:, :, 0], self.z_grid[0, 0, :],
            #           self.frequencies]
            # if y is not None:
            #     absice = [self.x_grid[:, 0, :],
            #               self.y_grid[:, 0, :], self.z_grid[:, 0, :],
            #           self.frequencies]
            # if x is not None:
            #     absice = [self.x_grid[0, :, :],
            #               self.y_grid[0, :, :], self.z_grid[0, :, :],
            #           self.frequencies]
            plt.pcolormesh(absice[np.where(test_dim)[0][0]],absice[np.where(test_dim)[0][1]],mean.T,**kwargs)
            # plt.pcolormesh(absice[np.where(test_dim)[0][0]],absice[np.where(test_dim)[0][1]],mean,**kwargs)
            plt.xlabel(labels[np.where(test_dim)[0][0]])
            plt.ylabel(labels[np.where(test_dim)[0][1]])

            #plt.gca().set_aspect('equal', adjustable='box')


        if len(mean.shape)==3:
            print("use anim function to plot 3D data field")
            return

    def plot_am(self,time=True,OA=False,x:float=None,y:float=None,z:float=None,freq:float=None,**kwargs):
        """
        Plots the amplitude modulation field. 
        if `x`, `y`, `z`, `freq` is set to `None` the plot is done in this direction. 
        Plots can be 1D (AM$(x)$ or AM$(freq)$ for example) or 2D (AM$(x,y)$ or AM$(x,freq)$ for example) 

        Args:
            time (bool, optional): Flag indicating if time data (True) or angular data (False) is plotted. Default is True.
            OA (bool, optional): Flag indicating if the OASPL data (True) or frequency data (False) is plotted. Default is False.
            x (float, optional): The x-coordinate to plot the data. Default is None.
            y (float, optional): The y-coordinate to plot the data. Default is None.
            z (float, optional): The z-coordinate to plot the data. Default is None.
            freq (float, optional): The frequency to plot the data. Default is None.
            **kwargs: Additional keyword arguments for plotting.
        """
        test_dim = [x is None, y is None, z is None, freq is None] 
        absice = [self.x_grid[:,0,0],self.y_grid[0,:,0],self.z_grid[0,0,:],self.frequencies]
        labels = ['x (m)', 'y (m)', 'z (m)', 'f (Hz)']
        if OA:
            test_dim[3] = False
            if time:
                am = self.OAam_time
            else:
                am = self.OAam_seg
        else :
            if time:
                am = self.am_time
            else:
                am = self.am_seg

        if freq is not None:
            ifreq = np.argmin(np.abs(self.frequencies-freq))
            am = am[...,ifreq:ifreq+1]
        if x is not None:
            ix = np.argmin(np.abs(self.x_grid[:,0,0]-x))
            am = am[ix:ix+1,...]
        if y is not None :
            iy = np.argmin(np.abs(self.y_grid[0,:,0]-y))
            am = am[:,iy:iy+1,...]
        if z is not None :
            iz = np.nonzero(self.z_grid[0,0,:]==z)[0][0]
            am = am[:,:,iz:iz+1,...]
        am = np.squeeze(am)

        if len(am.shape)==1:
            plt.plot(absice[np.where(test_dim)[0][0]], am,**kwargs)
            plt.xlabel(labels[np.where(test_dim)[0][0]])
            plt.ylabel('am (dB)')

        if len(am.shape)==2:
            # absice = [self.x_grid[:, :, 0],
            #           self.y_grid[:, :, 0], self.z_grid[0, 0, :],
            #           self.frequencies]
            plt.pcolormesh(absice[np.where(test_dim)[0][0]],absice[np.where(test_dim)[0][1]],am.T,**kwargs)
            plt.xlabel(labels[np.where(test_dim)[0][0]])
            plt.ylabel(labels[np.where(test_dim)[0][1]])
            #plt.pcolormesh(absice[np.where(test_dim)[0][1]],absice[np.where(test_dim)[0][0]],am,**kwargs)
            #plt.xlabel(labels[np.where(test_dim)[0][1]])
            #plt.ylabel(labels[np.where(test_dim)[0][0]])
            # plt.gca().set_aspect('equal', adjustable='box')

        if len(am.shape)==3:
            return 

    def anim_spl(self,time=True,OA=False,x:float=None,y:float=None,z:float=None,freq:float=None,**kwargs):
        """
        Create an animation of the SPL field.
        if `x`, `y`, `z`, `freq`, `freq` is set to `None` the plot is done in this direction.
        Animation can be 1D (SPL$(x,it)$ or SPL$(freq,it)$ for example) or 2D (SPL$(x,y,it)$ or SPL$(x,freq,it)$ for example)

        Args:
            time (bool, optional): Flag indicating if time data (True) or angular data (False) is plotted. Default is True.
            OA (bool, optional): Flag indicating if the OASPL data (True) or frequency data (False) is plotted. Default is False.
            x (float, optional): The x-coordinate to plot the data. Default is None.
            y (float, optional): The y-coordinate to plot the data. Default is None.
            z (float, optional): The z-coordinate to plot the data. Default is None.
            freq (float, optional): The frequency to plot the data. Default is None.
            **kwargs: Additional keyword arguments for plotting.
        """
        # choose Over all data or frequency data 
        # choose between beta data or time dat
        test_dim = [x is None, y is None, z is None, freq is None] 
        absice = [self.x_grid[:,0,0],self.y_grid[0,:,0],self.z_grid[0,0,:],self.frequencies]
        labels = ['x (m)', 'y (m)', 'z (m)', 'f (Hz)', 't (s)']

        if OA:
            test_dim[3] = False
            if time:
                spl = self.OASPL_time
            else:
                spl = self.OASPL_seg_tot
        else:
            if time:
                spl = self.SPL_time
            else:
                spl = np.sum(10**(self.SPL_seg/10),(3,4))

        if time:
            nframes = self.Nt
        else:
            nframes = self.Nbeta

        # find coresponding freq, ix, iy, iz 
        if freq is not None:
            ifreq = np.argmin(np.abs(self.frequencies-freq))
            spl = spl[...,ifreq:ifreq+1,:]
        if x is not None:
            ix = np.argmin(np.abs(self.x_grid[:,0,0]-x))
            print(ix)
            spl = spl[ix:ix+1,...]
        if y is not None :
            iy = np.argmin(np.abs(self.y_grid[0,:,0]-y))
            print(iy)
            spl = spl[:,iy:iy+1,...]
        if z is not None :
            iz = np.nonzero(self.z_grid[0,0,:]==z)[0][0]
            print(iz)
            spl = spl[:,:,iz:iz+1,...]
        spl = np.squeeze(spl)

        # plot either pcolor or lineplot 
        if len(spl.shape) == 2:
            fig, ax = plt.subplots()
            line, = ax.plot(absice[np.where(test_dim)[0][0]],spl[...,0].T,**kwargs)
            plt.xlabel(labels[np.where(test_dim)[0][0]])
            plt.ylabel('SPL (dB)')
            def animate(ii):
                line.set_ydata(spl[...,ii%nframes])
                return line,
            anim = FuncAnimation(fig, animate, frames=160, interval=30, repeat=True)
            return anim

        if len(spl.shape) == 3:
            logging.info("2D animation ..")
            fig, ax = plt.subplots()
            if z is not None:
                absice2D = [self.x_grid[:,:,0], self.y_grid[:,:,0], self.z_grid[:,:,0], self.frequencies]
            if x is not None:
                absice2D = [self.x_grid[0,:,:], self.y_grid[0,:,:], self.z_grid[0,:,:], self.frequencies]
            if y is not None:
                absice2D = [self.x_grid[:,0,:], self.y_grid[:,0,:], self.z_grid[:,0,:], self.frequencies]
            # absice2D = [self.x_grid[:,0,:], self.y_grid[:,0,:], self.z_grid[:,0,:], self.frequencies]

            cax = ax.pcolormesh(absice2D[np.where(test_dim)[0][0]],absice2D[np.where(test_dim)[0][1]],spl[...,0],**kwargs)
            # plt.gca().set_aspect('equal', adjustable='box')
            plt.xlabel(labels[np.where(test_dim)[0][0]])
            plt.ylabel(labels[np.where(test_dim)[0][1]])
            # cax = ax.pcolormesh(absice[np.where(test_dim)[0][1]],absice[np.where(test_dim)[0][0]],spl[...,0],**kwargs)
            # plt.gca().set_aspect('equal', adjustable='box')
            # plt.xlabel(labels[np.where(test_dim)[0][1]])
            # plt.ylabel(labels[np.where(test_dim)[0][0]])

            # plt.gca().invert_yaxis()
            
            divider = make_axes_locatable(ax)
            cb = divider.append_axes("right", size="5%", pad=0.05)
            cb.set_title('OASPL (dBA)')
            plt.colorbar(cax, cax=cb)
            plt.tight_layout()
            def animate(ii):
                # ax.clear()
                cax.set_array(spl[...,ii%nframes])
                plt.tight_layout()
            anim = FuncAnimation(fig, animate, frames=nframes, interval=100, repeat=True)
            return anim
        
        if len(spl.shape)==4:
            logging.warning('splwith 4 dimensions not allowed')
            return 
