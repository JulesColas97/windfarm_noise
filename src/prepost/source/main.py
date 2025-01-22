import numpy as np
import matplotlib.pyplot as plt
import pickle
import logging 
from tqdm import tqdm 
from scipy.integrate import simpson

from .wind_turbine import WindTurbine
from .atmos import Atmosphere
from .mesh import Mesh
from .utils import computeThirdOctaveFrequencies
from .mio import delta_L_ip
from .function import SPL_SWL_parallel
from ..les import  Les 
from ..utils import interp_weights, interpolate, integrateThirdOctave,Aweight

class Source:
    """
    A class to represent an acoustic source model for wind turbines.

    This class provides methods to manage wind turbine acoustics data, perform spatial
    interpolations, compute sound power levels, and visualize results.

    Attributes:
        atmos (Les): Atmospheric conditions.
        wt (WindTurbine): Wind turbine object.
        mesh (Mesh): Computational mesh.
        Spp (np.array): Spectral power distribution data.
    """
    atmos : Les
    wt : WindTurbine
    mesh : Mesh
    Spp : np.array

    def __init__(self,wt: WindTurbine=None,atmos: Les=None,mesh: Mesh=None,Spp: np.array=None):
        """
        Initialize the Source object with wind turbine, atmosphere, mesh, and spectral data.

        Args:
            wt (WindTurbine, optional): Wind turbine object. Defaults to None.
            atmos (Les, optional): Atmospheric conditions. Defaults to None.
            mesh (Mesh, optional): Computational mesh. Defaults to None.
            Spp (np.array, optional): Spectral power distribution data. Defaults to None.
        """
        self.wt = wt
        self.atmos = atmos
        self.mesh = mesh
        self.Spp = Spp

        if self.Spp is not None:
            self.nx = self.Spp.shape[0]
            self.nz = self.Spp.shape[1]
            self.ntau = self.Spp.shape[2]
            self.Nfreq = self.Spp.shape[4]
            self.Nbeta = self.Spp.shape[5]
        if self.wt is not None :
            self.Nseg = self.wt.Nseg
            self.Nblade = self.wt.Nblade


    def read_pickle(self,filename: str):
        """
        Load the source data from a pickle file.

        Args:
            filename (str): Path to the pickle file.
        """
        repo_out = pickle.load(open(filename,'rb'))
        self.wt = repo_out['wind_turbine']
        self.atmos = repo_out['atmos']
        self.mesh = repo_out['mesh']
        self.Spp = repo_out['Spp']

        self.nx = self.Spp.shape[0]
        self.nz = self.Spp.shape[1]
        self.ntau = self.Spp.shape[2]
        self.Nseg = self.wt.Nseg
        self.Nblade = self.wt.Nblade
        self.Nfreq = self.Spp.shape[4]
        self.Nbeta = self.Spp.shape[5]

    def save(self, fname: str, atmos_data: bool = True) -> None:
        """
        Save the current Source object to a file.

        Args:
            fname (str): The filename to save data.
            atmos_data (bool, optional): Whether to include atmospheric data. Defaults to True.
        """

        if atmos_data == False:
            self.atmos = None
        with open(fname,'wb') as file:
            pickle.dump(self.__dict__,file)

    def load(self, fname: str) -> None:
        """
        Load the Source object from a file.

        Args:
            fname (str): The filename to load data from.
        """
        print('loading source ...')
        with open(fname,'rb') as file:
            self.__dict__ = pickle.load(file)
        print('done.')

    def computeSpp(self, freq: np.array, Ncore: int = 16, BEM: bool = False) -> None:
        """
        Compute the spectral power distribution (Spp) for given frequencies.

        Args:
            freq (np.array): Array of frequencies.
            Ncore (int, optional): Number of CPU cores to use. Defaults to 16.
            BEM (bool, optional): Use BEM (Blade Element Momentum) method. Defaults to False.
        """
        self.frequencies = freq
        self.Nfreq = len(freq)
        # rotate mesh : 
        mesh = self.mesh
        # mesh.tau_array = mesh.tau_array - self.wt.tau
        mesh.tau_coord = mesh.tau_coord - self.wt.tau
        (self.Spp, self.Spp_tin, self.Spp_ten, self.SWL,
        self.AoA, self.U_inf,self.epsilon,self.U_rot,self.U_rel, self.a, self.adash  
                    )  = SPL_SWL_parallel(self.wt,self.atmos,mesh,freq,Ncore,BEM)
        # print(self.Spp.shape)
        # print(self.Spp_tin.shape)
        if self.mesh.polar:
            self.Spp = self.Spp.reshape(self.mesh.nx,self.mesh.nz,self.mesh.ntau,self.wt.Nblade*self.wt.Nseg,self.Nfreq,self.wt.Nbeta)
            self.Spp_tin = self.Spp_tin.reshape(self.mesh.nx,self.mesh.nz,self.mesh.ntau,self.wt.Nblade*self.wt.Nseg,self.Nfreq,self.wt.Nbeta)
            self.Spp_ten = self.Spp_ten.reshape(self.mesh.nx,self.mesh.nz,self.mesh.ntau,self.wt.Nblade*self.wt.Nseg,self.Nfreq,self.wt.Nbeta)
            self.SWL = self.SWL.reshape(self.mesh.nx,self.mesh.nz,self.mesh.ntau,self.wt.Nblade*self.wt.Nseg,self.Nfreq,self.wt.Nbeta)
            self.nx = self.Spp.shape[0]
            self.nz = self.Spp.shape[1]
            self.ntau = self.Spp.shape[2]
        else : 
            self.Spp = self.Spp.reshape(self.mesh.nx,self.mesh.nz,self.mesh.ny,self.wt.Nblade*self.wt.Nseg,self.Nfreq,self.wt.Nbeta)
            self.Spp_tin = self.Spp_tin.reshape(self.mesh.nx,self.mesh.nz,self.mesh.ny,self.wt.Nblade*self.wt.Nseg,self.Nfreq,self.wt.Nbeta)
            self.Spp_ten = self.Spp_ten.reshape(self.mesh.nx,self.mesh.nz,self.mesh.ny,self.wt.Nblade*self.wt.Nseg,self.Nfreq,self.wt.Nbeta)
            self.SWL = self.SWL.reshape(self.mesh.nx,self.mesh.nz,self.mesh.ny,self.wt.Nblade*self.wt.Nseg,self.Nfreq,self.wt.Nbeta)
            self.nx = self.Spp.shape[0]
            self.nz = self.Spp.shape[1]
            self.ny = self.Spp.shape[2]
        

        self.Nseg = self.Spp.shape[3]
        self.Nfreq = self.Spp.shape[4]
        self.Nbeta = self.Spp.shape[5]

    def interpolate_xy(self, x: np.array, y: np.array) -> None:
        """
        Perform 2D interpolation over the XY plane. 
        Moves form self.Spp(mesh.x_array, mesh.y_array, mesh.z_array) to 
        self.SppInterpolated(self.x_interpolated, self.y_interpolated, self.z_interpolated)

        Args:
            x (np.array): Array of x-coordinates.
            y (np.array): Array of y-coordinates.
        """
        logging.info('xy 2D interpolation ...')
        x_grid,y_grid = np.meshgrid(self.mesh.x_array,self.mesh.y_array,indexing='ij')

        xy_old = np.concatenate((x_grid.reshape((-1,1)),y_grid.reshape((-1,1))),1)

        logging.info('start interpolation ...')
        if (len(x.shape) == 1) and (len(y.shape==1)):
            x_cart, y_cart = np.meshgrid(x, y, indexing='ij')
            nx = x_cart.shape[0]
            ny = x_cart.shape[1]
            xy_new=np.zeros([nx* ny,2])
            xy_new[:, 0] = x_cart .flatten()
            xy_new[:, 1] = y_cart.flatten()
        elif (len(x.shape) == 2) and (len(y.shape) == 2) and np.all(x.shape==y.shape):
            nx = x.shape[0]
            ny = x.shape[1]
            print(nx)
            print(ny)
            xy_new=np.zeros([nx*ny,2])
            xy_new[:, 0]=x.flatten()
            xy_new[:, 1]=y.flatten()
            
            
        vtx, wts = interp_weights(xy_old, xy_new)
        logging.info('finished creating knot ...')

        self.SppInterpolated = np.zeros((nx,ny,self.nz,self.Nseg,self.Nfreq,self.Nbeta))

        logging.info('starting loop on height and frequency band ...')
        print(self.nz)
        print(self.Nseg)
        print(self.Nfreq)
        for iz in tqdm(range(self.nz)):
            for iseg in range(self.Nseg):
                print(iseg)
                for ifreq in range(self.Nfreq):
                    for ibeta in range(self.Nbeta):
                        values = self.Spp[:,iz,:,iseg,ifreq,ibeta].flatten()
                        valuesInterpolated = interpolate(values, vtx, wts)
                        self.SppInterpolated[:,:,iz,iseg,ifreq,ibeta] = valuesInterpolated.reshape(nx,ny)
        logging.info('done')

        # create new 3D mesh 
        self.x_interpolate = np.meshgrid(xy_new[:,0],
                                         self.mesh.z_array,
                                         indexing='ij')[0].reshape(nx,ny,self.nz)
        self.y_interpolate = np.meshgrid(xy_new[:,1],
                                         self.mesh.z_array,
                                         indexing='ij')[0].reshape(nx,ny,self.nz)
        self.z_interpolate = np.meshgrid(xy_new[:,1],
                                         self.mesh.z_array,
                                         indexing='ij')[1].reshape(nx,ny,self.nz)
        
        # [self.x_interpolate,self.y_interpolate,self.z_interpolate] = np.meshgrid(x,y,self.mesh.z_array,indexing='ij')


    def interpolate_xz(self, x: np.array, z: np.array) -> None:
        """
        Perform 2D interpolation over the XZ plane.
        Moves form self.Spp(mesh.x_array, mesh.y_array, mesh.z_array) to 
        self.SppInterpolated(self.x_interpolated, self.y_interpolated, self.z_interpolated)

        Args:
            x (np.array): Array of x-coordinates.
            z (np.array): Array of z-coordinates.
        """
        logging.info('xz 2D interpolation ...')
        x_grid,z_grid = np.meshgrid(self.mesh.x_array,self.mesh.z_array,indexing='ij')

        xz_old = np.concatenate((x_grid.reshape((-1,1)),z_grid.reshape((-1,1))),1)

        logging.info('start interpolation ...')
        x_new,z_new = np.meshgrid(x,z,indexing='ij')
        xz_cart=np.zeros([x_new.shape[0]*z_new .shape[1],2])
        xz_cart[:,0]=x_new .flatten()
        xz_cart[:,1]=z_new.flatten()
        vtx, wts = interp_weights(xz_old, xz_cart)
        logging.info('finished creating knot ...')

        self.SppInterpolated = np.zeros((x_new.shape[0],self.ny,x_new.shape[1],self.Nseg,self.Nfreq,self.Nbeta))

        logging.info('starting loop on height and frequency band ...')
        # for iy in tqdm(range(self.ny)):
        for iy in range(self.ny):
            logging.debug('iy = %s'%(iy))
            for iseg in range(self.Nseg):
                logging.debug('iseg = %s'%(iseg))
                for ifreq in range(self.Nfreq):
                    logging.debug('ifreq = %s'%(ifreq))
                    for ibeta in range(self.Nbeta):
                        values = self.Spp[:,:,iy,iseg,ifreq,ibeta].flatten()
                        valuesInterpolated = interpolate(values, vtx, wts)
                        self.SppInterpolated[:,iy,:,iseg,ifreq,ibeta] = valuesInterpolated.reshape(x_new.shape[0],x_new.shape[1])
        logging.info('done')
        [self.x_interpolate,self.y_interpolate,self.z_interpolate] = np.meshgrid(x,self.mesh.y_array,z,indexing='ij')


    def interpolate_yz(self, y: np.array, z: np.array) -> None:
        """
        Perform 2D interpolation over the YZ plane.
        Moves form self.Spp(mesh.x_array, mesh.y_array, mesh.z_array) to 
        self.SppInterpolated(self.x_interpolated, self.y_interpolated, self.z_interpolated)

        Args:
            y (np.array): Array of y-coordinates.
            z (np.array): Array of z-coordinates.
        """
        y_grid,z_grid = np.meshgrid(self.mesh.y_array,self.mesh.z_array,indexing='ij')

        yz_old = np.concatenate((y_grid.reshape((-1,1)),z_grid.reshape((-1,1))),1)

        print('start interpolation ...')
        y_new,z_new = np.meshgrid(y,z,indexing='ij')
        yz_cart=np.zeros([y_new.shape[0]*z_new.shape[1],2])
        yz_cart[:,0]=y_new .flatten()
        yz_cart[:,1]=z_new.flatten()
        vtx, wts = interp_weights(yz_old, yz_cart)
        print('finished creating knot ...')

        self.SppInterpolated = np.zeros((self.nx,y_new.shape[0],y_new.shape[1],self.Nseg,self.Nfreq,self.Nbeta))

        print('starting loop on height and frequency band ...')

        # transpose y /z 
        self.Spp = np.transpose(self.Spp, axes=(0,2,1,3,4,5))
        for ix in tqdm(range(self.nx)):
            for iseg in range(self.Nseg):
                for ifreq in range(self.Nfreq):
                    for ibeta in range(self.Nbeta):
                        values = self.Spp[ix,:,:,iseg,ifreq,ibeta].flatten()
                        valuesInterpolated = interpolate(values, vtx, wts)
                        self.SppInterpolated[ix,:,:,iseg,ifreq,ibeta] = valuesInterpolated.reshape(y_new.shape[0],y_new .shape[1])
        print('done')
        [self.x_interpolate, self.y_interpolate, self.z_interpolate] = np.meshgrid(self.mesh.x_array,y,z,indexing='ij')


    def convert_to_position(self):
        """
        Convert Spp and SWL data into position-based arrays.
        Useful for polar plot.
        """

        spp = np.transpose(self.SWL.reshape(
            self.nx, self.ny, self.nz, self.Nseg//self.Nblade,
            self.Nblade, self.Nfreq, self.Nbeta),
            axes=(0, 1, 2, 3, 4, 6, 5))

        self.swlPosition = np.zeros((self.nx, self.ny, self.nz,
                                self.Nseg//self.Nblade,self.Nbeta*self.Nblade, self.Nfreq))

        self.swlPosition[:, :, :, :, 0:12, :] = spp[:, :, :, :, 0, :, :]
        self.swlPosition[:, :, :, :, 12:24, :] = spp[:, :, :, :, 1, :, :]
        self.swlPosition[:, :, :, :, 24:, :] = spp[:, :, :, :, 2, :, :]

        spp = np.transpose(self.Spp_tin.reshape(
            self.nx, self.ny, self.nz, self.Nseg//self.Nblade,
            self.Nblade, self.Nfreq, self.Nbeta),
            axes=(0, 1, 2, 3, 4, 6, 5))

        self.tinPosition = np.zeros((self.nx, self.ny, self.nz,
                                self.Nseg//self.Nblade,self.Nbeta*self.Nblade, self.Nfreq))

        self.tinPosition[:, :, :, :, 0:12, :] = spp[:, :, :, :, 0, :, :]
        self.tinPosition[:, :, :, :, 12:24, :] = spp[:, :, :, :, 1, :, :]
        self.tinPosition[:, :, :, :, 24:, :] = spp[:, :, :, :, 2, :, :]

        spp = np.transpose(self.Spp_ten.reshape(
            self.nx, self.ny, self.nz, self.Nseg//self.Nblade,
            self.Nblade, self.Nfreq, self.Nbeta),
            axes=(0, 1, 2, 3, 4, 6, 5))

        self.tenPosition = np.zeros((self.nx, self.ny, self.nz,
                                self.Nseg//self.Nblade,self.Nbeta*self.Nblade, self.Nfreq))

        self.tenPosition[:, :, :, :, 0:12, :] = spp[:, :, :, :, 0, :, :]
        self.tenPosition[:, :, :, :, 12:24, :] = spp[:, :, :, :, 1, :, :]
        self.tenPosition[:, :, :, :, 24:, :] = spp[:, :, :, :, 2, :, :]

        spp = np.transpose(self.Spp.reshape(
            self.nx, self.ny, self.nz, self.Nseg//self.Nblade,
            self.Nblade, self.Nfreq, self.Nbeta),
            axes=(0, 1, 2, 3, 4, 6, 5))

        self.sppPosition = np.zeros((self.nx, self.ny, self.nz,
                                self.Nseg//self.Nblade,self.Nbeta*self.Nblade, self.Nfreq))

        self.sppPosition[:, :, :, :, 0:12, :] = spp[:, :, :, :, 0, :, :]
        self.sppPosition[:, :, :, :, 12:24, :] = spp[:, :, :, :, 1, :, :]
        self.sppPosition[:, :, :, :, 24:, :] = spp[:, :, :, :, 2, :, :]
        return

    def convert_to_one_blade(self):
        """
        Aggregate data to represent a single blade's contribution to sound power levels.
        """

        spp = np.transpose(self.SWL.reshape(
            self.ny, self.ny, self.nz, self.Nseg//self.Nblade,
            self.Nblade, self.Nfreq, self.Nbeta),
            axes=(0, 1, 2, 3, 4, 6, 5))

        self.sppOneBlade = np.zeros((self.nx, self.ny, self.nz,
                                self.Nbeta*self.Nblade, self.Nfreq))

        self.sppOneBlade[:, :, :, 0:12, :] = 10 * \
            np.log10(np.sum(10**(spp[:, :, :, :, 0, :, :]/10), (3)))
        self.sppOneBlade[:, :, :, 12:24, :] = 10 * \
            np.log10(np.sum(10**(spp[:, :, :, :,  1, :, :]/10), (3)))
        self.sppOneBlade[:, :, :, 24:, :] = 10 * \
            np.log10(np.sum(10**(spp[:, :,:, :, 2, :, :]/10), (3)))
        return

    def compute_oaswl(self) -> int:
        """
        Compute the overall A-weighted sound power levels.

        Returns:
            int: Returns -1 if an error occurs.
        """

        if not hasattr(self, "swlPosition"):
            print('error')
            return -1

        swl = self.swlPosition
        self.fc, swl_3rd = integrateThirdOctave(self.frequencies,10**(swl/10))
        self.oaswl_db = 10*np.log10(np.sum(swl_3rd,5))
        self.oaswl_dbA = 10*np.log10(np.sum(swl_3rd * 10**(Aweight(self.fc).reshape(1,1,1,1,1,-1)/10),5))

        spp = self.tinPosition
        self.fc, swl_3rd = integrateThirdOctave(self.frequencies,10**(spp/10))
        self.tin_dbA = 10*np.log10(np.sum(swl_3rd * 10**(Aweight(self.fc).reshape(1,1,1,1,1,-1)/10),5))

        spp = self.tenPosition
        self.fc, swl_3rd = integrateThirdOctave(self.frequencies,10**(spp/10))
        self.ten_dbA = 10*np.log10(np.sum(swl_3rd * 10**(Aweight(self.fc).reshape(1,1,1,1,1,-1)/10),5))

        spp = self.sppPosition
        self.fc, swl_3rd = integrateThirdOctave(self.frequencies,10**(spp/10))
        self.spp_dbA = 10*np.log10(np.sum(swl_3rd * 10**(Aweight(self.fc).reshape(1,1,1,1,1,-1)/10),5))
        return

    def print_info(self, x: float, y: float, z: float) -> None:
        """
        Print information about the sound power level at a given position.

        Args:
            x (float): X-coordinate.
            y (float): Y-coordinate.
            z (float): Z-coordinate.
        """
        ix = np.abs(self.mesh.x_array - x).argmin()
        iy = np.abs(self.mesh.y_array - y).argmin()
        iz = np.abs(self.mesh.z_array - z).argmin()

        swl_dbA = 10*np.log10(np.sum(10**(self.SWL[ix, iz, iy, ...]/10),
                  (0, 2))/self.Nbeta)+Aweight(self.frequencies)

        swl_db = 10*np.log10(np.sum(10**(self.SWL[ix, iz, iy, ...]/10),
                  (0, 2))/self.Nbeta)

        print('SWL (dBA)=')
        print(10*np.log10(simpson(y=10**(swl_dbA/10), x=self.frequencies)))

        print('SWL (dB)=')
        print(10*np.log10(simpson(y=10**(swl_db/10), x=self.frequencies)))

        print('omega (rpm)=')
        print(self.wt.omega * 60 / (2*np.pi))


    def plot_atmos(self) -> None:
        """
        Plot atmospheric data such as wind speed and turbulence.
        """
        plt.figure(0)
        plt.plot(self.atmos.U_inf,self.atmos.z_coord)
        plt.xlabel('$u_{inf}$ (m/s)')
        plt.ylabel('$z$ (m)')

        plt.figure(1)
        plt.plot(self.atmos.epsilon_Kol,self.atmos.z_coord)
        plt.xlabel('$\epsilon$ (m2/s3)')
        plt.ylabel('$z$ (m)')


    def plot_OASWL(self, ix: int, iy: int, iz: int, **kwargs) -> None:
        """
        Plot the overall A-weighted sound power level in a polar plot.

        Args:
            ix (int): Index along the x-axis.
            iy (int): Index along the y-axis.
            iz (int): Index along the z-axis.
        """
        self.convert_to_position()
        self.compute_oaswl()
        r = self.wt.seg
        self.wt.computeBeta()
        beta = self.wt.beta.reshape(self.wt.Nblade*self.wt.Nbeta)
        R,B = np.meshgrid(r,beta)

        fig = plt.gcf() 
        ax = fig.add_subplot(polar = True)
        cax = ax.pcolormesh(B.T + np.pi/2,R.T,self.oaswl_dbA[ix,iy,iz],**kwargs)
        plt.gca().set_yticklabels([])
        plt.gca().set_xticklabels([])
        # plt.colorbar(cax)

    def plot_tin_ten(self, ix: int, iy: int, iz: int, **kwargs) -> None:
        """
        Plot tonal and broadband noise contributions.

        Args:
            ix (int): Index along the x-axis.
            iy (int): Index along the y-axis.
            iz (int): Index along the z-axis.
        """
        self.convert_to_position()
        self.compute_oaswl()

        r = self.wt.seg
        self.wt.computeBeta()
        beta = self.wt.beta.reshape(self.wt.Nblade*self.wt.Nbeta)
        R,B = np.meshgrid(r,beta)

        R_tot = np.sqrt((self.mesh.x_array[ix]-self.wt.absolute_pos[0])**2 
                        + (self.mesh.y_array[iy]-self.wt.absolute_pos[1])**2)

        fig = plt.gcf()
        ax = fig.add_subplot(1, 3, 1, polar=True)
        cax = ax.pcolormesh(B.T + np.pi/2,R.T,
                            self.tin_dbA[ix,iy,iz]+10*np.log10(2*np.pi*R_tot**2),
                            **kwargs)
        plt.gca().set_yticklabels([])
        plt.gca().set_xticklabels([])
        plt.colorbar(cax)
        ax.set_title('SWL TIN')

        #fig = plt.gcf()
        ax = fig.add_subplot(1, 3, 2, polar=True)
        cax = ax.pcolormesh(B.T + np.pi/2,R.T,
                            self.ten_dbA[ix,iy,iz]+10*np.log10(2*np.pi*R_tot**2),
                            **kwargs)
        plt.gca().set_yticklabels([])
        plt.gca().set_xticklabels([])
        plt.colorbar(cax)
        ax.set_title('SWL TEN')


        # fig = plt.figure(1)
        ax = fig.add_subplot(1, 3, 3, polar=True)
        cax = ax.pcolormesh(B.T + np.pi/2,R.T,
                            self.spp_dbA[ix,iy,iz]+10*np.log10(2*np.pi*R_tot**2),
                            **kwargs)
        plt.gca().set_yticklabels([])
        plt.gca().set_xticklabels([])
        plt.colorbar(cax)
        ax.set_title('SWL tot')


    def plot_AoA(self, **kwargs) -> plt.Axes:
        """
        Plot the angle of attack distribution.

        Returns:
            plt.Axes: Matplotlib axis with the plot.
        """
        r = self.wt.seg
        self.wt.computeBeta()
        beta = self.wt.beta.reshape(self.wt.Nblade*self.wt.Nbeta)
        R,B = np.meshgrid(r,beta,indexing='ij')
        AoA = self.AoA.reshape(self.wt.Nblade,self.wt.Nseg,self.wt.Nbeta)
        AoA = AoA.transpose((1,0,2))
        AoA = AoA.reshape(self.wt.Nseg,self.wt.Nblade*self.wt.Nbeta)

        fig = plt.gcf()
        ax = fig.add_subplot(polar=True)
        
        #B = np.concatenate((B,B[:,0:1]),1)
        #R = np.concatenate((R,R[:,0:1]),1)
        #AoA = np.concatenate((AoA,AoA[:,0:1]),1)
        cax = ax.pcolormesh(B + np.pi/2,R,AoA,vmin=2,vmax=5,shading='auto',**kwargs)#,edgecolors='k')   
        plt.gca().set_yticklabels([])
        plt.gca().set_xticklabels([])
        # plt.colorbar()
        return cax


    def plot_all_u(self, **kwargs) -> None:
        """
        Plot various wind velocity components.
        """
        r = self.wt.seg
        self.wt.computeBeta()
        beta = self.wt.beta.reshape(self.wt.Nblade*self.wt.Nbeta)
        R,B = np.meshgrid(r,beta,indexing='ij')

        # fig = plt.figure()
        fig = plt.gcf()

        #----------------------------------------------------------------------
        eps = self.epsilon.reshape(self.wt.Nblade,self.wt.Nseg,self.wt.Nbeta)
        eps = eps.transpose((1,0,2))
        eps = eps.reshape(self.wt.Nseg,self.wt.Nblade*self.wt.Nbeta)

        ax = fig.add_subplot(1,4,1,polar=True)
        cax = ax.pcolormesh(B.T + np.pi/2, R.T, eps.T,cmap="Blues",
                            vmin=0.01, vmax=0.09,
                            shading='auto',**kwargs)#,edgecolors='k')   
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        plt.colorbar(cax)
        ax.set_title('eps')

        #----------------------------------------------------------------------
        U_inf = self.U_inf.reshape(self.wt.Nblade,self.wt.Nseg,self.wt.Nbeta)
        U_inf = U_inf.transpose((1,0,2))
        U_inf = U_inf.reshape(self.wt.Nseg,self.wt.Nblade*self.wt.Nbeta)

        ax = fig.add_subplot(1,4,2,polar=True)
        cax = ax.pcolormesh(B.T + np.pi/2, R.T, U_inf.T,shading='auto',
                            vmin=4,vmax=11,
                            cmap="RdYlBu_r", **kwargs)
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        plt.colorbar(cax)
        ax.set_title('U_inf')

        # ----------------------------------------------------------------------
        U_rel = self.U_rel.reshape(self.wt.Nblade,self.wt.Nseg,self.wt.Nbeta)
        U_rel = U_rel.transpose((1,0,2))
        U_rel = U_rel.reshape(self.wt.Nseg,self.wt.Nblade*self.wt.Nbeta)

        ax = fig.add_subplot(1,4,3,polar=True)
        cax = ax.pcolormesh(B.T + np.pi/2, R.T, U_rel.T,
                            vmin=10,vmax=60,**kwargs)
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        plt.colorbar(cax)
        ax.set_title('U_rel')

        # ----------------------------------------------------------------------
        AoA = self.AoA.reshape(self.wt.Nblade,self.wt.Nseg,self.wt.Nbeta)
        AoA = AoA.transpose((1,0,2))
        AoA = AoA.reshape(self.wt.Nseg,self.wt.Nblade*self.wt.Nbeta)

        ax = fig.add_subplot(1,4,4,polar=True)
        cax = ax.pcolormesh(B.T + np.pi/2,R.T,AoA.T,
                            vmin=2, vmax=10,
                            shading='auto',**kwargs)#,edgecolors='k')   
        plt.gca().set_yticklabels([])
        plt.gca().set_xticklabels([])
        plt.colorbar(cax)
        ax.set_title('AoA')


    def plot_epsilon(self, shading: str = 'auto', **kwargs) -> None:
        """
        Plot turbulence dissipation rate.

        Args:
            shading (str, optional): Shading style for the plot. Defaults to 'auto'.
        """
        r = self.wt.seg
        self.wt.computeBeta()
        beta = self.wt.beta.reshape(self.wt.Nblade*self.wt.Nbeta)
        R,B = np.meshgrid(r,beta)

        eps = self.epsilon.reshape(self.wt.Nblade,self.wt.Nseg,self.wt.Nbeta)
        eps = eps.transpose((1,0,2))
        eps = eps.reshape(self.wt.Nseg,self.wt.Nblade*self.wt.Nbeta)

        if shading=="gouraud":
            R = np.concatenate((R,R[0:1,:]),axis=0)
            B = np.concatenate((B,B[0:1,:]),axis=0)
            eps = np.concatenate((eps,eps[:,0:1]),axis=1)

        fig = plt.figure()
        ax = fig.add_subplot(polar=True)
        cax = ax.pcolormesh(B.T + np.pi/2, R.T, eps,shading=shading,**kwargs)#,edgecolors='k')   
        ax.set_yticklabels([])
        ax.set_xticklabels([])

    def plot_Uinf(self, shading: str = 'auto', **kwargs) -> None:
        """
        Plot free-stream wind velocity.

        Args:
            shading (str, optional): Shading style for the plot. Defaults to 'auto'.
        """
        r = self.wt.seg
        self.wt.computeBeta()
        beta = self.wt.beta.reshape(self.wt.Nblade*self.wt.Nbeta)
        R,B = np.meshgrid(r,beta)

        U_inf = self.U_inf.reshape(self.wt.Nblade,self.wt.Nseg,self.wt.Nbeta)
        U_inf = U_inf.transpose((1,0,2))
        U_inf = U_inf.reshape(self.wt.Nseg,self.wt.Nblade*self.wt.Nbeta)
        if shading=="gouraud":
            R = np.concatenate((R,R[0:1,:]),axis=0)
            B = np.concatenate((B,B[0:1,:]),axis=0)
            U_inf = np.concatenate((U_inf,U_inf[:,0:1]),axis=1)

        fig = plt.figure()
        # ax = fig.add_subplot(131, polar = True)
        # cax = ax.pcolormesh(B.T + np.pi/2,R.T,U_rel,cmap='RdBu_r',shading='auto')#,edgecolors='k')   
        # ax.set_yticklabels([])
        # ax.set_xticklabels([])
        #plt.colorbar(cax)

        ax = fig.add_subplot(polar=True)
        cax = ax.pcolormesh(B.T + np.pi/2, R.T, U_inf,shading=shading,**kwargs)
        ax.set_yticklabels([])
        ax.set_xticklabels([])

    def plot_Urel(self, shading: str = 'auto', **kwargs) -> None:
        """
        Plot relative wind velocity.

        Args:
            shading (str, optional): Shading style for the plot. Defaults to 'auto'.
        """
        r = self.wt.seg
        self.wt.computeBeta()
        beta = self.wt.beta.reshape(self.wt.Nblade*self.wt.Nbeta)
        R,B = np.meshgrid(r,beta)

        U_rel = self.U_rel.reshape(self.wt.Nblade,self.wt.Nseg,self.wt.Nbeta)
        U_rel = U_rel.transpose((1,0,2))
        U_rel = U_rel.reshape(self.wt.Nseg,self.wt.Nblade*self.wt.Nbeta)

        if shading=="gouraud":
            R = np.concatenate((R,R[0:1,:]),axis=0)
            B = np.concatenate((B,B[0:1,:]),axis=0)
            U_rel = np.concatenate((U_rel,U_rel[:,0:1]),axis=1)

        fig = plt.figure()
        ax = fig.add_subplot(polar=True)
        cax = ax.pcolormesh(B.T + np.pi/2, R.T, U_rel,shading=shading,**kwargs)
        ax.set_yticklabels([])
        ax.set_xticklabels([])

    def plot_Urot(self, shading: str = 'auto') -> None:
        """
        Plot the rotational velocity.

        Args:
            shading (str, optional): Shading style for the plot. Defaults to 'auto'.
        """
        r = self.wt.seg
        self.wt.computeBeta()
        beta = self.wt.beta.reshape(self.wt.Nblade*self.wt.Nbeta)
        R,B = np.meshgrid(r,beta)

        U_rot = self.U_rot.reshape(self.wt.Nblade,self.wt.Nseg,self.wt.Nbeta)
        U_rot = U_rot.transpose((1,0,2))
        U_rot = U_rot.reshape(self.wt.Nseg,self.wt.Nblade*self.wt.Nbeta)

        if shading=="gouraud":
            R = np.concatenate((R,R[0:1,:]),axis=0)
            B = np.concatenate((B,B[0:1,:]),axis=0)
            U_rot = np.concatenate((U_rot,U_rot[:,0:1]),axis=1)

        fig = plt.figure()
        ax = fig.add_subplot(polar=True)
        cax = ax.pcolormesh(B.T + np.pi/2, R.T, U_rot,shading=shading,**kwargs)   
        ax.set_yticklabels([])
        ax.set_xticklabels([])

    def plot_induction(self) -> None:
        """
        Plot the induction factors (axial and tangential).
        """
        r = self.wt.seg
        self.wt.computeBeta()
        beta = self.wt.beta.reshape(self.wt.Nblade*self.wt.Nbeta)
        R,B = np.meshgrid(r,beta)
        a = self.a.reshape(self.wt.Nblade,self.wt.Nseg,self.wt.Nbeta)
        a = a.transpose((1,0,2))
        a = a.reshape(self.wt.Nseg,self.wt.Nblade*self.wt.Nbeta)

        adash = self.adash.reshape(self.wt.Nblade,self.wt.Nseg,self.wt.Nbeta)
        adash = adash.transpose((1,0,2))
        adash = adash.reshape(self.wt.Nseg,self.wt.Nblade*self.wt.Nbeta)


        fig = plt.figure()
        ax = fig.add_subplot(121, polar = True)
        cax = ax.pcolormesh(B.T + np.pi/2,R.T,a,cmap='RdBu_r',shading='auto')#,edgecolors='k')   
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        plt.colorbar(cax)

        ax = fig.add_subplot(122, polar = True)
        cax = ax.pcolormesh(B.T + np.pi/2,R.T,adash,cmap='RdBu_r',shading='auto')#,edgecolors='k')   
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        plt.colorbar(cax)



if __name__=='__main__':
    # INPUTS
    #-------------------------------------------------------------------------------
    # define wind turbine geometry
    wind_turbine = WindTurbine()
    wind_turbine.default()
    wind_turbine.href = 200
    wind_turbine.Nblade = 3
    wind_turbine.Nbeta = 12

    # define atmosphere
    atmos = Atmosphere()
    path = '/home/lmfa/jcolas/Documents/DEV/ECOULEMENT/2Dhill/caseC_top/blue/output/'
    fname = 'tavg_u.h5'
    ratio = 0.451
    atmos.read_profile_from_les(path,fname,ratio,10,32)
    atmos.epsilon_Kol  = 0.01 * np.ones((len(atmos.z_coord)))
    # define frequencies
    fc = [5, 6.3, 8, 10, 12.5, 16, 20, 25,  31.5, 40, 50, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000]
    Nfc = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 4, 4, 4, 5,5]
    freq = computeThirdOctaveFrequencies(fc,Nfc)
    Nfreq = len(freq)
    # Define mesh grid
    xmin = 0; xmax = 5000;nx = 1001
    # ymin = 0; ymax = 0 ; ny = 1
    # zmin = 0; zmax = 300 ; nz = 61
    zmin = 2; zmax = 2 ; nz = 1
    taumin = 0; taumax = 0; ntau = 1
    mesh = Mesh(polar=True)
    mesh.create_polar_mesh(xmin,xmax,nx,zmin,zmax,nz,taumin,taumax,ntau)
    # mesh.polar_2_curvilinear(260,260,100)
    hpos = 0;hlength = 260;hheight = 100
    print(mesh.z_coord.shape)
    mesh.add_topography(hpos,hlength,hheight)
    print(mesh.z_coord.shape)
    Ncore = 16

    # Parallel initialization
    #-------------------------------------------------------------------------------
    Spp = SPL_SWL_parallel(wind_turbine,atmos,mesh,True,freq,Ncore)

    Spp = Spp.reshape(mesh.nx,mesh.nz,mesh.ntau,wind_turbine.Nseg*wind_turbine.Nblade,Nfreq,wind_turbine.Nbeta)

    pickle_file_name = 'sourceDataLineC'
    repo_out = {}
    repo_out['wind_turbine'] = wind_turbine
    repo_out['atmos'] = atmos
    repo_out['mesh'] = mesh
    repo_out['Spp'] = Spp
    pickle.dump(repo_out, open(pickle_file_name, "wb"), protocol=4)

    quit()
