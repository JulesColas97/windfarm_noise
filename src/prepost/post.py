import numpy as np
import h5py
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import (
    bisplrep,
    bisplev,
    RectBivariateSpline,
    RegularGridInterpolator,
)
import pickle
import os.path
import time 
import logging     
import logging 

from .utils import computeThirdOctaveFrequencies, chkList, interp_weights, interpolate
from .pre import Simu


class PeResults:
    def __init__(self, casename, iTurb, height, tau, dirname="./", distribute_tau=None):
        self.casename = casename
        self.iTurb = iTurb
        self.tau = tau
        if int(tau) == tau:
            self.tau_str = format(tau, "04d")
        else:
            self.tau_str = str(tau)

        self.height = height
        self.dirname = dirname
        self.distribute_tau = distribute_tau
        self.deltaL = None
        self.receiver = None
        self.THIRD = False

    def read_carto(self, fname):
        file = h5py.File(fname, "r")
        solutions = file.get("solution")
        mesh = file.get("mesh")
        frequencies = list(solutions.keys())
        self.frequencies = np.array([float(i) for i in frequencies])
        self.x = np.array(mesh["x"])
        self.z = np.array(mesh["y"])
        self.deltaL = np.zeros((len(self.x), len(self.z), len(self.frequencies)))
        # solution matrix (Nz,Nx,Nf)
        for ii in range(len(frequencies)):
            sol = solutions.get(frequencies[ii])
            self.deltaL[:, :, ii] = np.transpose(np.array(sol["deltaL"]))

    def read_receiver(self, fname):
        # if self.distribute_tau is not None:
        #     for kk in range(self.distribute_tau):
        #         path = self.dirname + 't'+str(self.iTurb)+'/'+str(self.height)+'/tau'+str(kk)+'/'
        #         if os.path.isfile(path+self.casename+'_'+self.tau_str+'.h5'):
        #             break
        # else :
        #     path = self.dirname + 't'+str(self.iTurb)+'/'+str(self.height)+'/'

        # file = h5py.File(path+self.casename+'_'+self.tau_str+'.h5', "r")
        file = h5py.File(fname, "r")
        mesh = file.get("mesh")
        self.x = np.array(mesh["x"])

        # read delta L saved at receivers
        receivers_all = file.get("receiver")
        self.heights = np.array(receivers_all.get("heights"))
        frequencies = list(receivers_all.keys())
        id_height = frequencies.index("heights")
        frequencies.pop(id_height)
        self.frequencies = np.array([float(i) for i in frequencies])
        self.z = self.heights

        self.receiver = np.zeros(
            (len(self.x), len(self.heights), len(self.frequencies))
        )
        for ii in range(len(self.frequencies)):
            rec = receivers_all.get(frequencies[ii])
            self.receiver[:, :, ii] = np.transpose(np.array(rec["deltaL"]))

    def read_planes(self, fname):
        t0 = time.time()
        file = h5py.File(fname, "r")
        t1 = time.time()
        mesh = file.get("mesh")
        self.z = np.array(mesh["y"])
        self.x = np.array(mesh["x"])

        # read delta L saved at receivers
        planes = file.get("planes")
        frequencies = list(planes.keys())
        self.frequencies = np.array([float(i) for i in frequencies])

        rec = planes.get(frequencies[0])

        self.xycoord = np.array(rec["xycoord"])
        self.xxcoord = np.array(rec["xxcoord"])
        self.xcount = np.array(rec["xcount"])
        self.nxplanes = self.xycoord.size

        self.yxcoord = np.array(rec["yxcoord"])
        self.yycoord = np.array(rec["yycoord"])
        self.ycount = np.array(rec["ycount"])
        # print(self.xcount)
        # print(self.xycoord)
        # print(self.xxcoord)

        # print(self.ycount)
        # print(self.yxcoord)
        # print(self.yycoord)
        self.nyplanes = self.yxcoord.size
        self.xplanes = np.zeros((self.nxplanes, len(self.z), len(self.frequencies)))
        self.yplanes = np.zeros((self.nyplanes, len(self.z), len(self.frequencies)))
        t2 = time.time()
        for ii in range(len(self.frequencies)):
            rec = planes.get(frequencies[ii])
            self.xplanes[:, :, ii] = np.transpose(np.array(rec["xplane"]))
            self.yplanes[:, :, ii] = np.transpose(np.array(rec["yplane"]))

        t_end = time.time()

        t_read =  t1 - t0
        t_info = t2 - t1
        t_get_dat = t_end - t2
        t_tot = t_end - t0 

        print('time -----')
        print(t_read/t_tot, t_info/t_tot, t_get_dat/t_tot) 
        print('----------')

    def plotSide(self, freq, cmap="RdBu_r"):
        ifreq = np.nonzero(self.frequencies == freq)[0][0]
        if ifreq is None:
            print("frequence was not calculated.")
            return
        plt.pcolormesh(
            self.x, self.z, self.deltaL[:, :, ifreq].T, cmap=cmap, shading="gouraud"
        )
        #plt.clim(-10, 10)
        # plt.gca().set_aspect("equal", adjustable="box")
        #plt.ylim(0, 300)

    def plotSide3octave(self, freq, cmap="RdBu_r"):
        ifreq = np.nonzero(self.fc == freq)[0][0]
        if ifreq is None:
            print("frequence was not calculated.")
            return
        plt.pcolormesh(
            self.x,
            self.z,
            self.deltaL3octave[:, :, ifreq].T,
            cmap=cmap,
            shading="gouraud",
        )
        plt.clim(-10, 10)
        # plt.gca().set_aspect("equal", adjustable="box")
        plt.ylim(0, 300)

    def plotLineTotalDeltaL(self, height):
        iheight = np.nonzero(self.heights == height)[0][0]
        plt.plot(
            self.x,
            10 * np.log10(np.sum(10 ** (self.deltaL3octave[:, iheight, :] / 10), 1)),
        )

    def plotSideTotalDeltaL(self, cmap="RdBu_r"):
        plt.pcolormesh(
            self.x,
            self.z,
            10 * np.log10(np.sum(10 ** (self.deltaL3octave[:, :, :] / 10), 2)).T,
            cmap=cmap,
            shading="gouraud",
        )
        plt.clim(-10, 10)
        plt.gca().set_aspect("equal", adjustable="box")
        plt.ylim(0, 300)

    def plotLine(self, freq, z,**kwargs):
        ifreq = np.nonzero(self.frequencies == freq)[0][0]
        iz = np.nonzero(self.z == z)[0][0]
        if ifreq is None:
            print("frequence was not calculated.")
            return
        if iz is None:
            print("receiver's height was not calculated.")
            return
        if self.receiver is not None:
            plt.plot(self.x, self.receiver[:, iz, ifreq],**kwargs)
        else:
            plt.plot(self.x, self.deltaL[:, iz, ifreq],**kwargs)

    def plotLine3octave(self, freq, height,**kwargs):
        ifreq = np.nonzero(self.fc == freq)[0][0]
        iheight = np.nonzero(self.heights == height)[0][0]
        if ifreq is None:
            print("frequence was not calculated.")
            return
        if iheight is None:
            print("receiver's height was not calculated.")
            return
        plt.plot(self.x, self.receiver3octave[:, iheight, ifreq],**kwargs)

    def compute3octave(self, fc=None, Nfc=None):
        if (fc is None) or (Nfc is None):
            fc = [50, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000]
            Nfc = [1,  1,  1,   1,   1,  1,   2,   2,   3,   4,   4,   4,   5,  5]
        freq = computeThirdOctaveFrequencies(fc, Nfc)
        if not np.all(freq == self.frequencies):
            print("central frequencies are not the same")
            return
        self.fc = np.array(fc)
        self.Nfc = np.array(Nfc)
        if self.deltaL is not None:
            self.deltaL3octave = np.zeros((len(self.x), self.deltaL.shape[1], len(fc)))
        if self.receiver is not None:
            self.receiver3octave = np.zeros((len(self.x), len(self.heights), len(fc)))
        for ifc in range(len(fc)):
            if self.deltaL is not None:
                self.deltaL3octave[:, :, ifc] = 10 * np.log10(
                np.sum(10** (self.deltaL[:,:,
                            int(np.sum(self.Nfc[0:ifc])) : np.sum(
                                self.Nfc[0 : ifc + 1]),]/10),2,)/Nfc[ifc])
            if self.receiver is not None:
                self.receiver3octave[:, :, ifc] = 10 * np.log10(
                        np.sum(10** (self.receiver[:, :,
                                                   int(np.sum(self.Nfc[0:ifc])) : np.sum(
                                                       self.Nfc[0 : ifc + 1]),]/10),2,)/Nfc[ifc])

    def plotLineOaspl(self, height):
        iheight = np.argmin(abs(self.heights - height))
        SPL, SPL_A, OASPL, OASPL_A = computeSPLLine(
            np.squeeze(self.receiver3octave[iheight, :, :]),
            np.squeeze(self.x),
            np.squeeze(self.heights[iheight] + self.h),
            self.fc,
            self.tau,
            self.height,
            c0=343,
        )
        plt.plot(self.x, OASPL_A)

    def plotFlowSide(self, fname, cmap="RdYlBu_r", xstep=1, zstep=1):
        f = h5py.File(fname, "r")

        # read mesh
        x = np.array(f["x"])
        z = np.array(f["z"])

        # read flow
        u = np.array(f["u"])
        # v = np.array(f["v"])
        c = np.array(f["c"])
        # print(c.shape)
        print(u.shape)
        fig, ax = plt.subplots(figsize=(10, 1.8))
        cax = ax.pcolormesh(
            x[::xstep, ::zstep],
            z[::xstep, ::zstep],
            u[::xstep, ::zstep],
            cmap=cmap,
            shading="gouraud",
        )
        ax.set_aspect("equal", adjustable="box")
        divider = make_axes_locatable(ax)
        cbax = divider.append_axes("right", size="2%", pad=0.05)
        ax.set_xlabel("$x$ (m)")
        ax.set_ylabel("$z$ (m)")
        fig.colorbar(cax, cax=cbax, label="$u$ (m/s)")
        plt.tight_layout()


        fig, ax = plt.subplots(figsize=(10, 1.8))
        cax = ax.pcolormesh(
            x[::xstep, ::zstep],
            z[::xstep, ::zstep],
            c[::xstep, ::zstep],
            cmap=cmap,
            shading="gouraud",
        )
        ax.set_aspect("equal", adjustable="box")
        divider = make_axes_locatable(ax)
        cbax = divider.append_axes("right", size="2%", pad=0.05)
        ax.set_xlabel("$x$ (m)")
        ax.set_ylabel("$z$ (m)")
        fig.colorbar(cax, cax=cbax, label="$c$ (m/s)")
        plt.tight_layout()
        return u,x,z

# define a deltaL field to use for the source model
class DeltaLField:
    """A class to read, and store a 2D deltaL field. created either from PE or LEE simulation results."""

    def __init__(self, dirname: str = None, casename: str = None):
        self.dirname = dirname
        self.casename = casename
        self.height = []

        # variable for top view concatenation
        self.tau = []
        self.frequencies = []
        self.deltaLlist = []
        self.x_polar = None
        self.y_polar = None
        self.deltaL_polar = None

        # variable for plane concatenation
        self.slice_list = []
        self.plane_list = []
        self.plane_xlist = []
        self.plane_ylist = []
        self.plane_zlist = []

        # variable for cartesian interpolation
        self.x_cart = None
        self.y_cart = None
        self.z_cart = None
        self.x_array = None
        self.y_array = None
        self.z_array = None
        self.deltaL_cart = None

        self.hindex = 0
        self.tauindex = 0
        if (dirname is not None) and (casename is not None):
            self.simu = Simu(casename)
            self.simu.load(dirname + casename + ".dat")
    
    def save(self,fname: str):
        """save data in a pick le .dat file. bif matrix are store in h5 files in order to open only part of the data
        """
        print('saving Delta L ...')
        self.deltaLlist = None
        self.simu = None
        with  h5py.File(fname[:-3]+'h5', "w") as f:
            if self.x_polar is not None:
                f.create_dataset('x_polar', data=self.x_polar)
                self.x_polar = None
            if self.y_polar is not None:
                f.create_dataset('y_polar', data=self.y_polar)
                self.y_polar = None
            if self.deltaL_polar is not None:
                f.create_dataset('deltaL_polar', data=self.deltaL_polar)
                self.deltaL_polar = None

            if self.x_cart is not None:
                f.create_dataset('x_cart', data=self.x_cart)
                self.x_cart = None

            if self.y_cart is not None:
                f.create_dataset('y_cart', data=self.y_cart)
                self.y_cart = None

            if self.z_cart is not None:
                f.create_dataset('z_cart', data=self.z_cart)
                self.z_cart = None

            if self.deltaL_cart is not None:
                f.create_dataset('deltaL_cart', data=self.deltaL_cart)
                self.deltaL_cart = None

            if self.x_array is not None:
                f.create_dataset('x_array', data=self.x_array)
                self.x_array = None

            if self.y_array is not None:
                f.create_dataset('y_array', data=self.y_array)
                self.y_array= None

            if self.z_array is not None:
                f.create_dataset('z_array', data=self.z_array)
                self.z_array= None

        """save class as self.name.dat"""
        with open(fname,'wb') as file:
            pickle.dump(self.__dict__,file)
        print('done.')

    def load(self,fname: str,polar=True,cart=True):
        """try load self.name.dat"""
        print('loading delta L ...')
        with open(fname,'rb') as file:
            self.__dict__ = pickle.load(file)
        print('done loading pickle')

        with h5py.File(fname[:-3]+'h5', "r") as f:
            if ('deltaL_polar' in f) and (polar):
                self.deltaL_polar = np.array(f['deltaL_polar'])
                self.x_polar = np.array(f['x_polar'])
                self.y_polar = np.array(f['y_polar'])

            if (('deltaL_cart' in f) and cart):
                self.deltaL_cart = np.array(f['deltaL_cart'])

                self.x_cart = np.array(f['x_cart'])
                self.y_cart = np.array(f['y_cart'])
                self.z_cart = np.array(f['z_cart'])

                self.x_array = np.array(f['x_array'])
                self.y_array = np.array(f['y_array'])
                self.z_array = np.array(f['z_array'])
        print('done.')

    def compute3octave(self, fc=None, Nfc=None):
        if fc is None or Nfc is None:
            fc = [50, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000]
            Nfc = [1,  1,  1,   1,   1,  1,   2,   2,   3,   4,   4,   4,   5,  5]
        freq = computeThirdOctaveFrequencies(fc, Nfc)
        if not np.all(freq == self.frequencies):
            print("central frequencies are not the same")
            return
        self.fc = np.array(fc)
        self.Nfc = np.array(Nfc)
        shape = self.deltaL_cart.shape
        deltaL3octave = np.zeros((shape[0],shape[1], shape[2], len(fc), self.nheight))

        for ifc in range(len(fc)):
            if0 = int(np.sum(self.Nfc[0:ifc]))
            if1 = int(np.sum(self.Nfc[0:ifc+1]))
            dl = self.deltaL_cart[:, :, :, if0:if1, :]
            deltaL3octave[:, :, :, ifc, :] = 10*np.log10(np.mean(10**(dl/10), 3))

        self.THIRD = True
        self.deltaL_cart = deltaL3octave
        self.frequencies = np.array(fc)
        self.Nfreq = len(fc)

    def save_old(self, fname: str):
        """save class as self.name.dat"""
        with open(fname, "wb") as file:
            pickle.dump(self.__dict__, file)

    def load_old(self, fname: str):
        print("loading Delta L ... ")
        """try load self.name.dat"""
        with open(fname, "rb") as file:
            self.__dict__ = pickle.load(file)
        print("done.")

    def fname(self, iTurb: int, tau: float, height: float, distribute_tau: int) -> str:
        """construct the name of th h5 file storing the pe results.
            the data stucture is such that if tau is not distibuted:
              <self.dirname>/t<iTurb>/<height>/<self.casename>_<tau>.h5
            if tau is distributed
              <self.dirname>/t<iTurb>/<height>/tau<itau>/<self.casename>_<tau>.h5

        Returns:
            str: return the string for the complete path to the h5 file
        """

        # convert angle to str format
        if int(tau) == tau:
            tau_str = format(tau, "04d")
        else:
            tau_str = str(tau)

        # if tau distributed find the right folder
        if distribute_tau is not None:
            for kk in range(distribute_tau):
                path = (
                    self.dirname
                    + "t"
                    + str(iTurb)
                    + "/"
                    + str(height)
                    + "/tau"
                    + str(kk)
                    + "/"
                )
                if os.path.isfile(path + self.casename + "_" + tau_str + ".h5"):
                    break
        # define directly the folder
        else:
            path = self.dirname + "t" + str(iTurb) + "/" + str(height) + "/"
        return path + self.casename + "_" + tau_str + ".h5"

    def read_carto(self, iTurb: int, height: float, tau: float, distribute_tau: bool):
        """

        Args:
            iTurb (int): _description_
            height (float): _description_
            tau (float): _description_
            distribute_tau (bool): _description_
        """
        res = PeResults(
            self.casename,
            iTurb,
            height,
            tau,
            dirname=self.dirname,
            distribute_tau=distribute_tau,
        )
        fname = self.fname(iTurb, tau, height, distribute_tau)
        res.read_carto(fname)
        # find height and tau index of the simulation
        [x_grid, z_grid] = np.meshgrid(res.x, res.z, indexing="ij")

        # add height and tau to list if not already read
        if height not in self.height:
            self.height.append(height)
        if tau not in self.tau:
            self.tau.append(tau)

        deltaL = {}
        deltaL["height"] = height
        deltaL["tau"] = tau
        deltaL["freq"] = res.frequencies
        deltaL["x"] = res.x
        deltaL["z"] = res.z
        deltaL["x_grid"] = x_grid
        deltaL["z_grid"] = z_grid
        deltaL["val"] = res.deltaL
        self.deltaLlist.append(deltaL)

    def read_receiver(
        self, iTurb: int, height: float, tau: float, distribute_tau: bool
    ):
        """read lines of receiver from PE results. stored in self.deltaL.list

        Args:
            iTurb (int): turbine index
            height (float): source height
            tau (float): propagation angle
            distribute_tau (bool): wether simulation were parallized over tau
        """
        res = PeResults(
            self.casename,
            iTurb,
            height,
            tau,
            dirname=self.dirname,
            distribute_tau=distribute_tau,
        )
        fname = self.fname(iTurb, tau, height, distribute_tau)
        res.read_receiver(fname)
        # find height and tau index of the simulation
        [x_grid, z_grid] = np.meshgrid(res.x, res.z, indexing="ij")

        # add height and tau to list if not already read
        if height not in self.height:
            self.height.append(height)
        if tau not in self.tau:
            self.tau.append(tau)

        deltaL = {}
        deltaL["height"] = height
        deltaL["tau"] = tau
        deltaL["freq"] = res.frequencies
        deltaL["x"] = res.x
        deltaL["z"] = res.heights
        deltaL["x_grid"] = x_grid
        deltaL["z_grid"] = z_grid
        deltaL["val"] = res.receiver
        self.deltaLlist.append(deltaL)

    def check_compatibility(self):
        """check compatibility between all solution recorded in self.deltaLlist"""
        if not self.deltaLlist:
            logging.warning("you need to read some PE results first ...")
            return -1

        # check frequencies
        frequencies = [deltaL["freq"] for deltaL in self.deltaLlist]
        flag,error_index = chkList(frequencies)
        if flag:
            logging.info("frequencies OK")
            logging.info( frequencies[0] )
        else:
            logging.warning("frequency are not the same")
            logging.warning("error for (h,tau)=(%s,%s)" % (self.deltaLlist[error_index]["height"], 
                                                           self.deltaLlist[error_index]["tau"]))
            return -1

        # check receievr Heights
        heights = [deltaL["z"] for deltaL in self.deltaLlist]
        if chkList(heights)[0]:
            logging.info("receiver heights OK")
            logging.info(heights[0])
        else:
            logging.warning("receiver heights are not the same")
            return -1
        # TODO add other check on grid and dx

    def concatenate_angles(self):
        """concatenate all loaded angles in a same array.
        Find minimum and maximum values for x.
        Create the array of shape (Nx,Nz, Ntau, Nfreq, Nheight) and fiel it with the resultst from deltaLlist.
        """

        # find xmin, xmax
        # ---------------------------------------------------------------------------------------
        # initialize min max values
        xmin = self.deltaLlist[0]["x"][0]
        xmax = self.deltaLlist[0]["x"][-1]

        # browse all deltaL to find min max values
        for ii in range(1, len(self.deltaLlist)):
            if self.deltaLlist[ii]["x"][0] < xmin:
                xmin = self.deltaLlist[ii]["x"][0]
            if self.deltaLlist[ii]["x"][-1] > xmax:
                xmax = self.deltaLlist[ii]["x"][-1]

        # rearange height and tau
        # ---------------------------------------------------------------------------------------
        self.height.sort()
        self.tau.sort()
        self.height = np.array(self.height)
        self.tau = np.array(self.tau)
        self.frequencies = self.deltaLlist[0]["freq"]

        self.ntau = len(self.tau)
        self.nheight = len(self.height)

        # reshaping
        # ---------------------------------------------------------------------------------------
        # assuming all dx are the same
        dx = self.deltaLlist[-1]["x"][1] - self.deltaLlist[-1]["x"][0]
        # creating grid for x
        self.x = np.arange(xmin, xmax + dx, dx)
        self.nx = len(self.x)
        # creating grid for z
        self.z = self.deltaLlist[-1]["z"]
        self.nz = len(self.z)

        # create the complete matrix of deltaL
        self.Nfreq = len(self.deltaLlist[-1]["freq"])
        self.deltaL_polar = np.zeros(
            (self.nx, self.nz, self.ntau, self.Nfreq, self.nheight)
        )
        # [self.x_grid,self.z_grid] = np.meshgrid(self.x, self.z,indexing='ij')

        for deltaL in self.deltaLlist:
            iheight = np.argmin(abs(self.height - deltaL["height"]))
            itau = np.argmin(abs(self.tau - deltaL["tau"]))
            # find corresponding index
            ixmin = np.argmin(abs(self.x - deltaL["x"][0]))
            ixmax = np.argmin(abs(self.x - deltaL["x"][-1]))
            logging.info("height = " + str(deltaL["height"]) + ", tau = " + str(deltaL["tau"]))
            self.deltaL_polar[ixmin : ixmax + 1, :, itau, :, iheight] = deltaL["val"][
                :, :, :
            ]

        # create polar mesh
        angles = np.reshape(self.tau * np.pi / 180, (1, -1))
        r = np.reshape(self.x, (-1, 1))

        self.x_polar = r * np.cos(angles)
        self.y_polar = r * np.sin(angles)

    def create_side_view(self, iTurb, tau=[0,180]):
        if len(self.deltaLlist) != 0:
            print('WARNING simulations were already loaded')
            return -1

        for h in self.simu.heights:
            for t in tau:
                print(h, t)
                self.read_carto(iTurb, h, t, self.simu.distribute_tau)
                print(self.deltaLlist[-1]['val'].shape)

        self.check_compatibility()
        # rearange height and tau
        # ---------------------------------------------------------------------
        self.height.sort()
        self.tau.sort()
        self.height = np.array(self.height)
        self.tau = np.array(self.tau)
        self.frequencies = self.deltaLlist[0]["freq"]

        self.ntau = len(self.tau)
        self.nheight = len(self.height)
        if self.ntau == 2:
            # reshaping
            # -----------------------------------------------------------------
            # creating grid for x
            x_0 = self.deltaLlist[0]["x"]
            nx_0 = len(x_0)
            x_180 = self.deltaLlist[1]["x"]
            nx_180 = len(x_180)

            self.x = np.concatenate([-x_180[-1:0:-1], x_0], axis=0)
            self.nx = len(self.x)

            # creating grid for z
            z_0 = self.deltaLlist[0]["z"]
            z_180 = self.deltaLlist[1]["z"]

            if (z_0 != z_180).any():
                print("WARNING z array are not the same ")
                return -1
        else:
            # if nx_180 = 0 put equal to 1 to avaoid negative nx
            nx_180 = 1
            self.x = self.deltaLlist[0]["x"]
            self.nx = len(self.x)

        self.z = self.deltaLlist[0]["z"]
        self.nz = len(self.z)

        # create the complete matrix of deltaL
        self.Nfreq = len(self.deltaLlist[0]["freq"])
        self.deltaL_cart = np.zeros(
            (self.nx, 1, self.nz,
                self.Nfreq, self.nheight)
        )

        print(self.deltaL_cart.shape)
        for deltaL in self.deltaLlist:
            iheight = np.argmin(abs(self.height - deltaL["height"]))
            # itau = np.argmin(abs(self.tau - deltaL["tau"]))

            if deltaL["tau"] == 180:
                print("height = " +
                      str(deltaL["height"]) +
                      ", tau = " +
                      str(deltaL["tau"]))
                self.deltaL_cart[:nx_180-1, 0,
                                 :, :, iheight] = deltaL["val"][-1:0:-1, :, :]

            if deltaL["tau"] == 0:

                print("height = " +
                      str(deltaL["height"]) +
                      ", tau = " +
                      str(deltaL["tau"]))

                print(nx_180)
                print(deltaL["val"].shape)

                self.deltaL_cart[nx_180-1:, 0, :, :, iheight] = deltaL["val"]
                # change fi nc_180 = 0
                #self.deltaL_cart[nx_180:, 0, :, :, iheight] = deltaL["val"]

        [self.x_cart, self.y_cart, self.z_cart] = np.meshgrid(
            self.x, 0, self.z, indexing="ij"
        )
        self.x_array = self.x
        self.y_array = np.array([0])
        self.z_array = self.z
        self.ny = 1


    # linear interpolate Pe results on a cartesian grid x,y
    def interpolate_from_polar(self, x: np.ndarray, y: np.ndarray):
        """Define a cartesian mesh from x and y. Interpolate from deltaL_polar.
            create the deltaL_cartesian field.

        Args:
            x (np.ndarray): 1D array of the x position on which to interpolate
            y (np.ndarray): 1D array of the y position on which to interpolate
        """
        # create vector of coordinates of the original grid
        xy_polar = np.concatenate(
            (self.x_polar.reshape((-1, 1)), self.y_polar.reshape((-1, 1))), 1
        )

        print("start interpolation ...")
        # create vector of coordinates of the new grid
        self.x_cart, self.y_cart = np.meshgrid(x, y, indexing="ij")
        xy_cart = np.zeros([self.x_cart.shape[0] * self.x_cart.shape[1], 2])
        xy_cart[:, 0] = self.x_cart.flatten()
        xy_cart[:, 1] = self.y_cart.flatten()
        # create knots for the interpolation
        vtx, wts = interp_weights(xy_polar, xy_cart)
        print("finished creating knot ...")

        self.deltaL_cart = np.zeros(
            (
                self.x_cart.shape[0],
                self.x_cart.shape[1],
                self.nz,
                self.Nfreq,
                self.nheight,
            )
        )
        # loop of z, freq, heights and interpolate using the previously computes knots
        print("starting loop on height and frequency band ...")
        for iz in range(self.nz):
            for ifreq in range(self.Nfreq):
                for iheight in range(self.nheight):
                    values = self.deltaL_polar[:, iz, :, ifreq, iheight].flatten()
                    valuesInterpolated = interpolate(values, vtx, wts)
                    self.deltaL_cart[
                        :, :, iz, ifreq, iheight
                    ] = valuesInterpolated.reshape(
                        self.x_cart.shape[0], self.x_cart.shape[1]
                    )
        print("done")
        # save the interpolation grid
        [self.x_cart, self.y_cart, self.z_cart] = np.meshgrid(
            x, y, self.z, indexing="ij"
        )
        self.x_array = x
        self.y_array = y
        self.z_array = self.z
        self.nx = len(x)
        self.ny = len(y)
        self.nz = len(self.z)



    # def read_slice(
    #     self,
    #     iTurb: int,
    #     height: float,
    #     tau: float,
    #     distribute_tau: bool,
    #     xplane: float = None,
    #     yplane: float = None,
    # ):
    #     """read the slice recorded for plane reconstruction (only implemented for PE results)
    #         for a given height, angle turbine and a given plane position, check if the pe simulation has recorded on slice.

    #     Args:
    #         iTurb (int): index of the turbine
    #         height (float): height of the source
    #         tau (float): propagation angle
    #         distribute_tau (bool): check if tau distributed during the simulations
    #         xplane (float, optional): recorded x constant planes. Defaults to None.
    #         yplane (float, optional): recorder y constant planes. Defaults to None.
    #     """
    #     # check that at least of the input was given
    #     if (xplane is None) and (yplane is None):
    #         print("must choose a position for xplane or yplane")
    #         return -1

    #     # read Pe results
    #     res = PeResults(
    #         self.casename,
    #         iTurb,
    #         height,
    #         tau,
    #         dirname=self.dirname,
    #         distribute_tau=distribute_tau,
    #     )
    #     fname = self.fname(iTurb, tau, height, distribute_tau)

    #     # read from side output in case x==0 and angle coresponds
    #     if xplane == 0:
    #         if tau == 90 or tau == -90 or tau == 270:
    #             res.read_carto(fname)
    #             if height not in self.height:
    #                 self.height.append(height)
    #             # create a slice for each x position --> this allow to kee the same post processing rootine for all planes
    #             for ii, y in enumerate(res.x):
    #                 deltaL = {}
    #                 deltaL["height"] = height
    #                 deltaL["tau"] = tau
    #                 deltaL["freq"] = res.frequencies
    #                 deltaL["z"] = res.z
    #                 deltaL["x"] = 0
    #                 deltaL["y"] = y * np.sin(tau * np.pi / 180)
    #                 deltaL["val"] = res.deltaL[ii, :, :]
    #                 self.slice_list.append(deltaL)
    #         return
    #     # read from side output in case y==0
    #     elif yplane == 0:
    #         if tau == 0 or tau == 180 or tau == -180:
    #             res.read_carto(fname)
    #             if height not in self.height:
    #                 self.height.append(height)
    #             for ii, x in enumerate(res.x):
    #                 deltaL = {}
    #                 deltaL["height"] = height
    #                 deltaL["tau"] = tau
    #                 deltaL["freq"] = res.frequencies
    #                 deltaL["z"] = res.z
    #                 deltaL["x"] = x * np.cos(tau * np.pi / 180)
    #                 deltaL["y"] = 0
    #                 deltaL["val"] = res.deltaL[ii, :, :]
    #                 self.slice_list.append(deltaL)
    #         return

    #     # for all other cases read slice of plane
    #     else:
    #         # read slice of planes
    #         res.read_planes(fname)

    #         if height not in self.height:
    #             self.height.append(height)

    #         if yplane is None:
    #             # compar recorded slice with coresponding plane
    #             for ii, x in enumerate(res.xxcoord):
    #                 if x == xplane:
    #                     # check if a slice was recorded
    #                     if res.xcount[ii] == 1:
    #                         deltaL = {}
    #                         deltaL["height"] = height
    #                         deltaL["tau"] = tau
    #                         deltaL["freq"] = res.frequencies
    #                         deltaL["z"] = res.z
    #                         deltaL["x"] = x
    #                         deltaL["y"] = res.xycoord[ii]
    #                         deltaL["val"] = res.xplanes[ii, :]
    #                         self.slice_list.append(deltaL)
    #                         break
    #         elif xplane is None:
    #             for ii, y in enumerate(res.yycoord):
    #                 if y == yplane:
    #                     # check if slice was recorded
    #                     if res.ycount[ii] == 1:
    #                         deltaL = {}
    #                         deltaL["height"] = height
    #                         deltaL["tau"] = tau
    #                         deltaL["freq"] = res.frequencies
    #                         deltaL["z"] = res.z
    #                         deltaL["x"] = res.yxcoord[ii]
    #                         deltaL["y"] = y
    #                         deltaL["val"] = res.yplanes[ii, :]
    #                         self.slice_list.append(deltaL)
    #                         break

    # def concatenate_slice(self):
    #     # test
    #     ####################################################################
    #     if len(self.slice_list) == 0:
    #         print("no plane loaded")
    #         return -1

    #     # concatenate all the position recorded at the first height
    #     x = np.array([d["x"] for d in self.slice_list if d["height"] == self.height[0]])
    #     y = np.array([d["y"] for d in self.slice_list if d["height"] == self.height[0]])
    #     z = self.slice_list[0]["z"]
    #     frequencies = self.slice_list[0]["freq"]

    #     for ii, h in enumerate(self.height[1:]):
    #         # take all the slice position recorded for iith height, check that all position are the same
    #         x_new = np.array([d["x"] for d in self.slice_list if d["height"] == h])
    #         y_new = np.array([d["y"] for d in self.slice_list if d["height"] == h])
    #         z = self.slice_list[ii]["z"]
    #         frequencies = self.slice_list[0]["freq"]
    #         if np.any(x != x_new):
    #             print("planes are not the same for all heights")
    #             return -1
    #         if np.any(y != y_new):
    #             print("plane are not the same  for all heights")
    #             return -1

    #     # sort heights
    #     self.height.sort()
    #     self.height = np.array(self.height)
    #     self.nheight = len(self.height)
    #     # remove all duplicate of x position (coming from several hgeight sbeeing recorded )
    #     x = np.unique(x)
    #     x.sort()

    #     # remove all duplicate of y position
    #     y = np.unique(y)
    #     y.sort()

    #     # take the first recorded z position on the slice
    #     # TODO add check on z for all slice and heigths
    #     z = self.slice_list[0]["z"]

    #     # take first frequency recorded
    #     # TODO add check pn frequency
    #     self.frequencies = self.slice_list[0]["freq"]
    #     self.Nfreq = len(self.frequencies)

    #     # create deltaL matrix for one plane
    #     deltaL = np.zeros(
    #         (len(x), len(y), len(z), len(self.frequencies), len(self.height))
    #     )
    #     for deltaL_slice in self.slice_list:
    #         iheight = np.argmin(abs(self.height - deltaL_slice["height"]))
    #         ix = np.argmin(abs(x - deltaL_slice["x"]))
    #         iy = np.argmin(abs(y - deltaL_slice["y"]))
    #         deltaL[ix, iy, :, :, iheight] = deltaL_slice["val"]

    #     self.plane_list.append(deltaL)
    #     self.plane_xlist.append(x)
    #     self.plane_ylist.append(y)
    #     self.plane_zlist.append(z)

    def read_all_slices(self,
                        iTurb: int,
                        height_list: list,
                        tau_list: list,
                        distribute_tau: bool,
                        xplanes: list = [],
                        yplanes: list = []):
        

        self.xslice_list = [[] for k in range(len(xplanes))]
        self.yslice_list = [[] for k in range(len(yplanes))]
        t_tot = 0 
        t_read = 0
        t_plane = 0 
        t_save = 0 

        t_init = time.time() 
        # loop over all angles and height solution 
        for height in height_list:
            for tau in tau_list:
                print(height,tau)
                # read Pe results

                # read slice results 
                res = PeResults(
                    self.casename,
                    iTurb,
                    height,
                    tau,
                    dirname=self.dirname,
                    distribute_tau=distribute_tau,
                )
                fname = self.fname(iTurb, tau, height, distribute_tau)
                t0 = time.time()
                res.read_planes(fname)
                t_read += time.time() - t0
                # add height to vector of heights
                if height not in self.height:
                    self.height.append(height)


                t0 = time.time()
                # loop over the x constant planes 
                for ii,xplane in enumerate(xplanes):
                    # read from side output in case x==0 and angle coresponds
                    if xplane == 0:
                        if tau == 90 or tau == -90 or tau == 270:
                            res.read_carto(fname)
                            if height not in self.height:
                                self.height.append(height)
                            # create a slice for each x position --> this allow to kee the same post processing rootine for all planes
                            for jj, y in enumerate(res.x):
                                deltaL = {}
                                deltaL["height"] = height
                                deltaL["tau"] = tau
                                deltaL["freq"] = res.frequencies
                                deltaL["z"] = res.z
                                deltaL["x"] = 0
                                deltaL["y"] = y * np.sin(tau * np.pi / 180)
                                deltaL["val"] = res.deltaL[jj, :, :]
                                self.xslice_list[ii].append(deltaL)

                    # get the coresponding slice         
                    else:
                         # compar recorded slice with coresponding plane
                        for jj, x in enumerate(res.xxcoord):
                            if x == xplane:
                                # check if a slice was recorded
                                if res.xcount[jj] == 1:
                                    deltaL = {}
                                    deltaL["height"] = height
                                    deltaL["tau"] = tau
                                    deltaL["freq"] = res.frequencies
                                    deltaL["z"] = res.z
                                    deltaL["x"] = x
                                    deltaL["y"] = res.xycoord[jj]
                                    deltaL["val"] = res.xplanes[jj, :]
                                    self.xslice_list[ii].append(deltaL)
                t_save += time.time() - t0 

                # loop over the y constant planes 
                for ii,yplane in enumerate(yplanes):
                    # read from side output in case y==0
                    if yplane == 0:
                        if tau == 0 or tau == 180 or tau == -180:
                            res.read_carto(fname)
                            for jj, x in enumerate(res.x):
                                deltaL = {}
                                deltaL["height"] = height
                                deltaL["tau"] = tau
                                deltaL["freq"] = res.frequencies
                                deltaL["z"] = res.z
                                deltaL["x"] = x * np.cos(tau * np.pi / 180)
                                deltaL["y"] = 0
                                deltaL["val"] = res.deltaL[jj, :, :]
                                self.yslice_list[ii].append(deltaL)
                    else : 
                        for jj, y in enumerate(res.yycoord):
                            if y == yplane:
                                # check if slice was recorded
                                if res.ycount[jj] == 1:
                                    deltaL = {}
                                    deltaL["height"] = height
                                    deltaL["tau"] = tau
                                    deltaL["freq"] = res.frequencies
                                    deltaL["z"] = res.z
                                    deltaL["x"] = res.yxcoord[jj]
                                    deltaL["y"] = y
                                    deltaL["val"] = res.yplanes[jj, :]
                                    self.yslice_list[ii].append(deltaL)
            t_tot = time.time() - t_init
            print('t_read = ' + str (t_read/t_tot * 100))
            print('t_save = ' + str(t_save/t_tot * 100))

    def concatenate_xslices(self):
        # test
        ####################################################################
        if len(self.xslice_list) == 0:
            print("no plane loaded")
            return -1

        # loop over the planes (remember xslice_list = [[dict0,dict1,...], ... [dict0,dict1]])
        for slice_list in self.xslice_list :
            # concatenate all the position recorded at the first height
            x = np.array([d["x"] for d in slice_list if d["height"] == self.height[0]])
            y = np.array([d["y"] for d in slice_list if d["height"] == self.height[0]])
            z = slice_list[0]["z"]
            frequencies = slice_list[0]["freq"]

            for ii, h in enumerate(self.height[1:]):
                # take all the slice position recorded for iith height, check that all position are the same
                x_new = np.array([d["x"] for d in slice_list if d["height"] == h])
                y_new = np.array([d["y"] for d in slice_list if d["height"] == h])
                z = slice_list[ii]["z"]
                frequencies = slice_list[0]["freq"]
                if np.any(x != x_new):
                    print("planes are not the same for all heights")
                    return -1
                if np.any(y != y_new):
                    print("plane are not the same  for all heights")
                    return -1

            # sort heights
            self.height.sort()
            self.height = np.array(self.height)
            self.nheight = len(self.height)
            # remove all duplicate of x position (coming from several hgeight sbeeing recorded )
            x = np.unique(x)
            x.sort()

            # remove all duplicate of y position
            y = np.unique(y)
            y.sort()

            # take the first recorded z position on the slice
            # TODO add check on z for all slice and heigths
            z = slice_list[0]["z"]

            # take first frequency recorded
            # TODO add check pn frequency
            self.frequencies = slice_list[0]["freq"]
            self.Nfreq = len(self.frequencies)

            # create deltaL matrix for one plane
            deltaL = np.zeros(
                (len(x), len(y), len(z), len(self.frequencies), len(self.height))
            )
            for deltaL_slice in slice_list:
                iheight = np.argmin(abs(self.height - deltaL_slice["height"]))
                ix = np.argmin(abs(x - deltaL_slice["x"]))
                iy = np.argmin(abs(y - deltaL_slice["y"]))
                deltaL[ix, iy, :, :, iheight] = deltaL_slice["val"]

            self.plane_list.append(deltaL)
            self.plane_xlist.append(x)
            self.plane_ylist.append(y)
            self.plane_zlist.append(z)


    def concatenate_yslices(self):
        # test
        ####################################################################
        if len(self.yslice_list) == 0:
            print("no plane loaded")
            return -1

        # loop over the planes (remember xslice_list = [[dict0,dict1,...], ... [dict0,dict1]])
        for slice_list in self.yslice_list :
            # concatenate all the position recorded at the first height
            x = np.array([d["x"] for d in slice_list if d["height"] == self.height[0]])
            y = np.array([d["y"] for d in slice_list if d["height"] == self.height[0]])
            z = slice_list[0]["z"]
            frequencies = slice_list[0]["freq"]

            for ii, h in enumerate(self.height[1:]):
                # take all the slice position recorded for iith height, check that all position are the same
                x_new = np.array([d["x"] for d in slice_list if d["height"] == h])
                y_new = np.array([d["y"] for d in slice_list if d["height"] == h])
                z = slice_list[ii]["z"]
                frequencies = slice_list[0]["freq"]
                if np.any(x != x_new):
                    print("planes are not the same for all heights")
                    return -1
                if np.any(y != y_new):
                    print("plane are not the same  for all heights")
                    return -1

            # sort heights
            self.height.sort()
            self.height = np.array(self.height)
            self.nheight = len(self.height)
            # remove all duplicate of x position (coming from several hgeight sbeeing recorded )
            x = np.unique(x)
            x.sort()

            # remove all duplicate of y position
            y = np.unique(y)
            y.sort()

            # take the first recorded z position on the slice
            # TODO add check on z for all slice and heigths
            z = slice_list[0]["z"]

            # take first frequency recorded
            # TODO add check pn frequency
            self.frequencies = slice_list[0]["freq"]
            self.Nfreq = len(self.frequencies)

            # create deltaL matrix for one plane
            deltaL = np.zeros(
                (len(x), len(y), len(z), len(self.frequencies), len(self.height))
            )
            for deltaL_slice in slice_list:
                iheight = np.argmin(abs(self.height - deltaL_slice["height"]))
                ix = np.argmin(abs(x - deltaL_slice["x"]))
                iy = np.argmin(abs(y - deltaL_slice["y"]))
                deltaL[ix, iy, :, :, iheight] = deltaL_slice["val"]

            self.plane_list.append(deltaL)
            self.plane_xlist.append(x)
            self.plane_ylist.append(y)
            self.plane_zlist.append(z)


    def interpolate_planes(self, x: np.ndarray = None, y: np.ndarray = None):
        """Interpolate from previously recorded slice. 

        Args:
            x (np.ndarray, optional): _description_. Defaults to None.
            y (np.ndarray, optional): _description_. Defaults to None.

        Returns:
            _type_: _description_
        """
        if x is not None:
            # check
            if len(self.plane_list) == 0:
                print("not planes loaded")
                return -1
            if self.plane_list[0].shape[1] != 1:
                print("probably loaded x constant planes and not y constant")
                return -1

            z = self.plane_zlist[0]
            self.deltaL_cart = np.zeros(
                (len(x), len(self.plane_ylist), len(z), self.Nfreq, self.nheight)
            )
            y_coord = []
            X, Z = np.meshgrid(x, z, indexing="ij")
            for ii in range(len(self.plane_list)):
                x_grid, z_grid = np.meshgrid(
                    self.plane_xlist[ii], self.plane_zlist[ii], indexing="ij"
                )
                y_coord.append(self.plane_ylist[ii][0])
                for ifreq in range(self.Nfreq):
                    for iheight in range(self.nheight):
                        # if not a lot of point spline interpolation
                        # if x_grid.size<200000 :
                        if x_grid.size < 2:
                            tck = bisplrep(
                                x_grid,
                                z_grid,
                                np.squeeze(
                                    self.plane_list[ii][:, :, :, ifreq, iheight]
                                ),
                            )
                            self.deltaL_cart[:, ii, :, ifreq, iheight] = bisplev(
                                x, z, tck
                            )
                        # else linear interpolation
                        else:
                            tck = RegularGridInterpolator(
                                (self.plane_xlist[ii], self.plane_zlist[ii]),
                                np.squeeze(
                                    self.plane_list[ii][:, :, :, ifreq, iheight]
                                ),
                                bounds_error=False,
                                fill_value=None,
                                method="linear",
                            )
                            self.deltaL_cart[:, ii, :, ifreq, iheight] = tck((X, Z))

            y_coord = np.array(y_coord)
            [self.x_cart, self.y_cart, self.z_cart] = np.meshgrid(
                x, y_coord, z, indexing="ij"
            )
            self.x_array = x
            self.y_array = y_coord
            self.z_array = z
            self.nx = len(x)
            self.ny = len(y_coord)
            self.nz = len(z)

        if y is not None:
            # check
            if len(self.plane_list) == 0:
                print("not planes loaded")
                return -1
            if self.plane_list[0].shape[0] != 1:
                print("probably loaded y constant planes and not x constant")
                return -1

            z = self.plane_zlist[0]
            self.deltaL_cart = np.zeros(
                (len(self.plane_xlist), len(y), len(z), self.Nfreq, self.nheight)
            )

            x_coord = []
            Y, Z = np.meshgrid(y, z, indexing="ij")
            for ii in range(len(self.plane_list)):
                y_grid, z_grid = np.meshgrid(
                    self.plane_ylist[ii], self.plane_zlist[ii], indexing="ij"
                )
                x_coord.append(self.plane_xlist[ii][0])
                for ifreq in range(self.Nfreq):
                    for iheight in range(self.nheight):
                        # if y_grid.size<200000 :
                        if y_grid.size < 2:
                            tck = bisplrep(
                                y_grid,
                                z_grid,
                                np.squeeze(
                                    self.plane_list[ii][:, :, :, ifreq, iheight]
                                ),
                            )
                            self.deltaL_cart[ii, :, :, ifreq, iheight] = bisplev(
                                y, z, tck
                            )
                        else:
                            tck = RegularGridInterpolator(
                                (self.plane_ylist[ii], self.plane_zlist[ii]),
                                np.squeeze(
                                    self.plane_list[ii][:, :, :, ifreq, iheight]
                                ),
                                bounds_error=False,
                                fill_value=None,
                                method="linear",
                            )
                            self.deltaL_cart[ii, :, :, ifreq, iheight] = tck((Y, Z))
            x_coord = np.array(x_coord)

            [self.x_cart, self.y_cart, self.z_cart] = np.meshgrid(
                x_coord, y, z, indexing="ij"
            )
            self.x_array = x_coord
            self.y_array = y
            self.z_array = z
            self.nx = len(x_coord)
            self.ny = len(y)
            self.nz = len(z)
        return

    # duplicate matrix size and devide grid step by alpha
    # TODO warning the refine method for tau is not the same than for x and z
    def refine(self, alpha: float, axis: int):
        # x axis
        if axis == 0:
            # temp value
            deltaL = np.zeros(
                (
                    (self.nx - 1) * alpha + 1,
                    self.nz,
                    self.ntau,
                    self.Nfreq,
                    self.nheight,
                )
            )
            x_grid = np.zeros(((self.nx - 1) * alpha + 1, self.nz, self.ntau))
            z_grid = np.zeros(((self.nx - 1) * alpha + 1, self.nz, self.ntau))
            # tau_grid = np.zeros(((self.nx-1)*alpha +1,self.nz,self.ntau))

            deltaL[::alpha, ...] = self.deltaL
            x_grid[::alpha, ...] = self.x_grid
            z_grid[::alpha, ...] = self.z_grid
            # tau_grid[::alpha,...] = self.tau_grid
            for ii in range(alpha - 1):
                deltaL[ii + 1 :: alpha, ...] = self.deltaL[:-1, ...]
                x_grid[ii + 1 :: alpha, ...] = self.x_grid[:-1, ...]
                z_grid[ii + 1 :: alpha, ...] = self.z_grid[:-1, ...]
                # tau_grid[ii+1::alpha,...] = self.tau_grid[:-1,...]

            self.deltaL = deltaL
            self.x_grid = x_grid
            self.z_grid = z_grid
            # self.tau_grid = tau_grid
            self.nx = (self.nx - 1) * alpha + 1
            self.x = np.linspace(self.x[0], self.x[-1], self.nx)
        # z axis
        if axis == 1:
            deltaL = np.zeros(
                (
                    self.nx,
                    (self.nz - 1) * alpha + 1,
                    self.ntau,
                    self.Nfreq,
                    self.nheight,
                )
            )
            x_grid = np.zeros((self.nx, (self.nz - 1) * alpha + 1, self.ntau))
            z_grid = np.zeros((self.nx, (self.nz - 1) * alpha + 1, self.ntau))
            # tau_grid = np.zeros((self.nx,(self.nz-1)*alpha +1,self.ntau))

            deltaL[:, ::alpha, ...] = self.deltaL
            x_grid[:, ::alpha, ...] = self.x_grid
            z_grid[:, ::alpha, ...] = self.z_grid
            # tau_grid[:,::alpha,...] = self.tau_grid
            for ii in range(alpha - 1):
                deltaL[:, ii + 1 :: alpha, ...] = self.deltaL[:, :-1, ...]
                x_grid[:, ii + 1 :: alpha, ...] = self.x_grid[:, :-1, ...]
                z_grid[:, ii + 1 :: alpha, ...] = self.z_grid[:, :-1, ...]
                # tau_grid[:,ii+1::alpha,...] = self.tau_grid[:,:-1,...]

            self.nz = (self.nz - 1) * alpha + 1
            self.z = np.linspace(self.z[0], self.z[-1], self.nz)
            self.deltaL = deltaL
            self.x_grid = x_grid
            self.z_grid = z_grid
            # self.tau_grid = tau_grid

        # tau axis
        # TODO finish modify grid for tau axis
        if axis == 2:
            deltaL = np.zeros(
                (self.nx, self.nz, 2 * self.ntau - 1, self.Nfreq, self.nheight)
            )
            tau = np.zeros((2 * self.ntau - 1))

            for ii in range(self.ntau - 1):
                deltaL[:, :, 2 * ii, :] = self.deltaL[:, :, ii, :, :]
                deltaL[:, :, 2 * ii + 1, :] = 0.5 * (
                    self.deltaL[:, :, ii, :, :] + self.deltaL[:, :, ii + 1, :, :]
                )
                tau[2 * ii] = self.tau[ii]
                tau[2 * ii + 1] = 0.5 * (self.tau[ii] + self.tau[ii + 1])

            deltaL[:, :, -1, :] = self.deltaL[:, :, -1, :]
            tau[-1] = self.tau[-1]

            self.deltaL = deltaL
            self.tau = tau
            self.ntau = len(self.tau)

    def refine_angle(self):
        deltaL = np.zeros(
            (self.nx, self.nz, 2 * self.ntau - 1, self.Nfreq, self.nheight)
        )
        tau = np.zeros((2 * self.ntau - 1))

        for ii in range(self.ntau - 1):
            deltaL[:, :, 2 * ii, :] = self.deltaL_polar[:, :, ii, :, :]
            deltaL[:, :, 2 * ii + 1, :] = 0.5 * (
                self.deltaL_polar[:, :, ii, :, :]
                + self.deltaL_polar[:, :, ii + 1, :, :]
            )
            tau[2 * ii] = self.tau[ii]
            tau[2 * ii + 1] = 0.5 * (self.tau[ii] + self.tau[ii + 1])

        deltaL[:, :, -1, :] = self.deltaL_polar[:, :, -1, :]
        tau[-1] = self.tau[-1]

        self.deltaL_polar = deltaL
        self.tau = tau
        self.ntau = len(self.tau)

        # create polar mesh
        angles = np.reshape(self.tau * np.pi / 180, (1, -1))
        r = np.reshape(self.x, (-1, 1))

        self.x_polar = r * np.cos(angles)
        self.y_polar = r * np.sin(angles)

    # create a results for tau = 360 from result tau=0
    def loop_angle(self):
        """used to create tau = 360 from tau =0. This is done before using self.interpolate_from_polar. 

        """
        # check if the first angle is indeed 0 
        if self.tau[0] != 0:
            print("first angle is not 0")
            return -1

        # create the new deltaL 
        deltaL = np.zeros((self.nx, self.nz, self.ntau + 1, self.Nfreq, self.nheight))
        tau = np.zeros((self.ntau + 1,))
        deltaL[:, :, 0:-1, :, :] = self.deltaL_polar
        deltaL[:, :, -1, :, :] = self.deltaL_polar[:, :, 0, :, :]

        tau[0:-1] = self.tau
        tau[-1] = 360
        self.tau = tau
        self.deltaL_polar = deltaL
        self.ntau = len(self.tau)

    def interpolate_xz(self, x: np.ndarray, z: np.ndarray):
        x_grid, z_grid = np.meshgrid(self.x_array, self.z_array, indexing="ij")

        xz_old = np.concatenate((x_grid.reshape((-1, 1)), z_grid.reshape((-1, 1))), 1)

        print("start interpolation ...")
        self.x_new, self.z_new = np.meshgrid(x, z, indexing="ij")
        xz_new = np.zeros([self.x_new.shape[0] * self.z_new.shape[1], 2])
        xz_new[:, 0] = self.x_new.flatten()
        xz_new[:, 1] = self.z_new.flatten()
        vtx, wts = interp_weights(xz_old, xz_new)
        print("finished creating knot ...")

        deltaLInterpolatedNew = np.zeros(
            (
                self.x_new.shape[0],
                self.ny,
                self.x_new.shape[1],
                self.Nfreq,
                self.nheight,
            )
        )

        print("starting loop on height and frequency band ...")
        for iy in range(self.ny):
            for ifreq in range(self.Nfreq):
                for iheight in range(self.nheight):
                    deltaLInterpolatedNew[:, iy, :, ifreq, iheight] = interpolate(
                        self.deltaL_cart[:, iy, :, ifreq, iheight].flatten(),
                        vtx,
                        wts,
                    ).reshape(self.x_new.shape[0], self.x_new.shape[1])
        print("done")
        self.deltaL_cart = deltaLInterpolatedNew
        [self.x_cart, self.y_cart, self.z_cart] = np.meshgrid(
            x, self.y_array, z, indexing="ij"
        )
        self.x_array = x
        self.z_array = z
        self.nx = len(x)
        self.nz = len(z)

    def interpolate_xy(self, x: np.ndarray, y: np.ndarray):
        x_grid, y_grid = np.meshgrid(self.x_array, self.y_array, indexing="ij")

        xy_old = np.concatenate((x_grid.reshape((-1, 1)), y_grid.reshape((-1, 1))), 1)

        print("start interpolation ...")
        self.x_new, self.y_new = np.meshgrid(x, y, indexing="ij")
        xy_new = np.zeros([self.x_new.shape[0] * self.y_new.shape[1], 2])
        xy_new[:, 0] = self.x_new.flatten()
        xy_new[:, 1] = self.y_new.flatten()
        vtx, wts = interp_weights(xy_old, xy_new)
        print("finished creating knot ...")

        deltaLInterpolatedNew = np.zeros(
            (
                self.x_new.shape[0],
                self.x_new.shape[1],
                self.nz,
                self.Nfreq,
                self.nheight,
            )
        )

        print("starting loop on height and frequency band ...")
        for iz in range(self.nz):
            for ifreq in range(self.Nfreq):
                for iheight in range(self.nheight):
                    deltaLInterpolatedNew[:, :, iz, ifreq, iheight] = interpolate(
                        self.deltaL_cart[:, :, iz, ifreq, iheight].flatten(),
                        vtx,
                        wts,
                    ).reshape(self.x_new.shape[0], self.x_new.shape[1])
        print("done")
        self.deltaL_cart = deltaLInterpolatedNew
        [self.x_cart, self.y_cart, self.z_cart] = np.meshgrid(
            x, y, self.z_array, indexing="ij"
        )
        self.x_array = x
        self.y_array = y
        self.nx = len(x)
        self.ny = len(y)

    def interpolate_yz(self, y: np.ndarray, z: np.ndarray):
        y_grid, z_grid = np.meshgrid(self.y_array, self.z_array, indexing="ij")

        yz_old = np.concatenate((y_grid.reshape((-1, 1)), z_grid.reshape((-1, 1))), 1)

        print("start interpolation ...")
        self.y_new, self.z_new = np.meshgrid(y, z, indexing="ij")
        yz_new = np.zeros([self.y_new.shape[0] * self.z_new.shape[1], 2])
        yz_new[:, 0] = self.y_new.flatten()
        yz_new[:, 1] = self.z_new.flatten()
        vtx, wts = interp_weights(yz_old, yz_new)
        print("finished creating knot ...")

        deltaLInterpolatedNew = np.zeros(
            (
                self.nx,
                self.y_new.shape[0],
                self.x_new.shape[1],
                self.Nfreq,
                self.nheight,
            )
        )

        print("starting loop on height and frequency band ...")
        for ix in range(self.nx):
            for ifreq in range(self.Nfreq):
                for iheight in range(self.nheight):
                    deltaLInterpolatedNew[ix, :, :, ifreq, iheight] = interpolate(
                        self.deltaL_cart[ix, :, :, ifreq, iheight].flatten(),
                        vtx,
                        wts,
                    ).reshape(self.x_new.shape[0], self.x_new.shape[1])
        print("done")
        self.deltaL_cart = deltaLInterpolatedNew
        [self.x_cart, self.y_cart, self.z_cart] = np.meshgrid(
            self.x_array, y, z, indexing="ij"
        )
        self.y_array = y
        self.z_array = z
        self.ny = len(y)
        self.nz = len(z)

    def shift_domain(self, xshift:float, yshift:float):
        if self.x_polar is not None:
            self.x_polar = self.x_polar + xshift
            self.y_polar = self.y_polar + yshift

        if self.x_cart is not None:
            self.x_cart = self.x_cart + xshift
            self.y_cart = self.y_cart + yshift
            self.x_array = self.x_array + xshift
            self.y_array = self.y_array + yshift

        if self.plane_xlist is not None:
            for ii in range(len(self.plane_xlist)):
                self.plane_xlist[ii] = self.plane_xlist[ii] + xshift  
            # self.plane_xlist = np.array(self.plane_xlist, dtype=float)
            # self.plane_xlist = self.plane_xlist + xshift

        if self.plane_ylist is not None:
            for ii in range(len(self.plane_ylist)):
                self.plane_ylist[ii] = self.plane_ylist[ii] + yshift  
            # self.plane_ylist = np.array(self.plane_ylist, dtype=float)
            # self.plane_ylist = self.plane_ylist + yshift


    # def reduce(self, alpha, axis):
    #     # x axis
    #     if axis == 0:
    #         self.deltaL = self.deltaL[::alpha, ...]
    #         self.x = self.x[::alpha]
    #         self.x_grid = self.x_grid[::alpha, ...]
    #         self.z_grid = self.z_grid[::alpha, ...]
    #         # self.tau_grid = self.tau_grid [::alpha,...]
    #         self.nx = len(self.x)
    #     if axis == 1:
    #         self.deltaL = self.deltaL[:, ::alpha, ...]
    #         self.z = self.z[::alpha]
    #         self.x_grid = self.x_grid[:, ::alpha, ...]
    #         self.z_grid = self.z_grid[:, ::alpha, ...]
    #         # self.tau_grid = self.tau_grid [:,::alpha,...]
    #         self.nz = len(self.z)

    #     # TODO finish this for tau
    #     if axis == 2:
    #         self.deltaL = self.deltaL[:, :, ::alpha, ...]
    #         self.tau = self.tau[::alpha]
    #         self.ntau = len(self.tau)

    # def truncate(self, max, axis, min=0):
    #     if axis == 0:
    #         Nxmax = np.argmin(abs(self.x - max))
    #         Nmin = np.argmin(abs(self.x - min))
    #         self.deltaL = self.deltaL[Nmin : Nxmax + 1, ...]
    #         self.x = self.x[Nmin : Nxmax + 1]
    #         self.x_grid = self.x_grid[Nmin : Nxmax + 1, ...]
    #         self.z_grid = self.z_grid[Nmin : Nxmax + 1, ...]
    #     if axis == 1:
    #         Nzmax = np.argmin(abs(self.z - max))
    #         Nmin = np.argmin(abs(self.z - min))
    #         self.deltaL = self.deltaL[:, Nmin : Nzmax + 1, ...]
    #         self.z = self.z[Nmin : Nzmax + 1]
    #         self.x_grid = self.x_grid[:, Nmin : Nzmax + 1, ...]
    #         self.z_grid = self.z_grid[:, Nmin : Nzmax + 1, ...]

    # # axis is the ax orthogonal to the plane
    # def takeOnePlane(self, pos1, axis):
    #     if axis == 0:
    #         ipos1 = np.argmin(abs(self.x - pos1))
    #         self.deltaL = self.deltaL[ipos1, ...]
    #     if axis == 1:
    #         ipos1 = np.argmin(abs(self.z - pos1))
    #         self.deltaL = self.deltaL[:, ipos1, ...]
    #     if axis == 2:
    #         ipos1 = np.argmin(abs(self.tau - pos1))
    #         self.deltaL = self.deltaL[:, :, ipos1, ...]

    # # here axis is the ax of the line
    # def takeOneLine(self, pos1, pos2, axis):
    #     if axis == 0:
    #         ipos1 = np.argmin(abs(self.z - pos1))
    #         ipos2 = np.argmin(abs(self.tau - pos2))
    #         self.deltaL = self.deltaL[:, ipos1, ipos2, ...]
    #         self.x_grid = self.x_grid[:, ipos1, ...]
    #         self.z_grid = self.z_grid[:, ipos1, ...]
    #     if axis == 1:
    #         ipos1 = np.argmin(abs(self.x - pos1))
    #         ipos2 = np.argmin(abs(self.tau - pos2))
    #         self.deltaL = self.deltaL[ipos1, :, ipos2, ...]
    #         self.x_grid = self.x_grid[ipos1, :, ipos2, ...]
    #         self.z_grid = self.z_grid[ipos1, :, ipos2, ...]
    #     if axis == 2:
    #         ipos1 = np.argmin(abs(self.x - pos1))
    #         ipos2 = np.argmin(abs(self.z - pos2))
    #         self.deltaL = self.deltaL[ipos1, ipos2, ...]
    #         self.x_grid = self.x_grid[ipos1, :, ipos2, ...]
    #         self.z_grid = self.z_grid[ipos1, :, ipos2, ...]

    # def reshape(self, x_array, z_array=None):
    #     print("delta L shape : " + str(self.deltaL.shape))
    #     print("reshaping Delta L ...")
    #     # Change delta X to fit
    #     if (x_array[1] - x_array[0]) > (self.x[1] - self.x[0]):
    #         alpha = (x_array[1] - x_array[0]) / (self.x[1] - self.x[0])
    #         if int(alpha) != alpha:
    #             print("mesh incompatible")
    #             quit()
    #         self.reduce(int(alpha), 0)
    #     if (x_array[1] - x_array[0]) < (self.x[1] - self.x[0]):
    #         alpha = (self.x[1] - self.x[0]) / (x_array[1] - x_array[0])
    #         if int(alpha) != alpha:
    #             print("mesh incompatible")
    #             quit()
    #         self.refine(int(alpha), 0)
    #     if z_array is not None:
    #         # Change delta z to fit
    #         if (z_array[1] - z_array[0]) > (self.z[1] - self.z[0]):
    #             alpha = (z_array[1] - z_array[0]) / (self.z[1] - self.z[0])
    #             if int(alpha) != alpha:
    #                 print("mesh incompatible")
    #                 quit()
    #             self.reduce(int(alpha), 1)
    #         if (z_array[1] - z_array[0]) < (self.z[1] - self.z[0]):
    #             alpha = (self.z[1] - self.z[0]) / (z_array[1] - z_array[0])
    #             if int(alpha) != alpha:
    #                 print("mesh incompatible")
    #                 quit()
    #             self.refine(int(alpha), 1)

    #     # reside domain
    #     if x_array[-1] <= self.x[-1]:
    #         self.truncate(x_array[-1], 0)
    #     else:
    #         print("delta L to small")
    #     if z_array is not None:
    #         if z_array[-1] <= self.z[-1]:
    #             self.truncate(z_array[-1], 1)
    #         else:
    #             print("delta L to small ")
    #     print("delta L shape : " + str(self.deltaL.shape))

    def plot_top_polar(self, z: float, freq: float, height: float, cmap="RdBu_r", **kwargs):
        """plot top view of deltaL field at a given z, frequency, and for a given source

        Args:
            z (float): receiver height
            freq (float): frequency
            height (float): source height
            cmap (str, optional): colormap. Defaults to 'RdBu_r'.
        """
        ifreq = np.nonzero(self.frequencies == freq)[0][0]
        iz = np.nonzero(self.z == z)[0][0]
        iheight = np.nonzero(self.height == height)[0][0]

        cax = plt.pcolormesh(
            self.x_polar,
            self.y_polar,
            self.deltaL_polar[:, iz, :, ifreq, iheight],
            cmap=cmap,
            shading="auto",
            **kwargs
        )
        plt.gca().set_aspect("equal", adjustable="box")
        plt.xlabel("x (m)")
        plt.ylabel("y (m)")
        ax = plt.gca()
        divider = make_axes_locatable(ax)
        cb = divider.append_axes("right", size="5%", pad=0.05)
        cb.set_title("$\Delta L$ (dB)")
        plt.colorbar(cax, cax=cb)
        plt.tight_layout()
        return ax

    def plotLineRaw(self, tau, z, freq, height, **kwargs):
        itau = np.nonzero(self.tau == tau)[0][0]
        ifreq = np.nonzero(self.frequencies == freq)[0][0]
        iz = np.nonzero(self.z == z)[0][0]
        iheight = np.nonzero(self.height == height)[0][0]
        plt.plot(self.x_polar[:, itau], self.deltaL_polar[:, iz, itau, ifreq, iheight], **kwargs)

    def plotDirectivityRaw(self, r, z, freq, height, **kwargs):
        ir = np.nonzero(self.x == r)[0][0]
        ifreq = np.nonzero(self.frequencies == freq)[0][0]
        iz = np.nonzero(self.z == z)[0][0]
        iheight = np.nonzero(self.height == height)[0][0]

        # fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
        # plt.axes(projection='polar')
        # plt.gca().projection = 'polar'
        plt.polar(
            self.tau * np.pi / 180, self.deltaL[ir, iz, :, ifreq, iheight], **kwargs
        )
        # ax.plot(self.tau*np.pi/180,self.deltaL[ir,iz,:,ifreq,iheight])

    def plot_top_cart(self, z, freq, height,ax=None, **kwargs):
        ifreq = np.nonzero(self.frequencies == freq)[0][0]
        iz = np.nonzero(self.z == z)[0][0]
        iheight = np.nonzero(self.height == height)[0][0]

        if ax is not(None):
            cax = ax.pcolormesh(
                self.x_cart[:,:,iz],
                self.y_cart[:,:,iz],
                self.deltaL_cart[:, :, iz, ifreq, iheight],
                **kwargs
            )
        else:
            ax = plt.pcolormesh(
                self.x_cart[:,:,iz],
                self.y_cart[:,:,iz],
                self.deltaL_cart[:, :, iz, ifreq, iheight],
                **kwargs
            )
            cax = None
        # plt.gca().set_aspect("equal", adjustable="box")
        return cax

    def plot_yz(self, x, freq, height, **kwargs):
        ifreq = np.nonzero(self.frequencies == freq)[0][0]
        ix = np.nonzero(self.x_cart[:,0,0] == x)[0][0]
        iheight = np.nonzero(self.height == height)[0][0]

        plt.pcolormesh(
            self.y_cart[ix,:,:],
            self.z_cart[ix,:,:],
            self.deltaL_cart[ix, :, :, ifreq, iheight],
            **kwargs
        )
        plt.gca().set_aspect("equal", adjustable="box")

    def plot_xz(self, y, freq, height, **kwargs):
        ifreq = np.nonzero(self.frequencies == freq)[0][0]
        iy = np.nonzero(self.y_cart[0,:,0] == y)[0][0]
        iheight = np.nonzero(self.height == height)[0][0]

        plt.pcolormesh(
            self.x_cart[:,iy,:],
            self.z_cart[:,iy,:],
            self.deltaL_cart[:, iy, :, ifreq, iheight],
            **kwargs
        )
        # plt.gca().set_aspect("equal", adjustable="box")

    def plot_x(self, y, z, freq, height, **kwargs):
        ifreq = np.nonzero(self.frequencies == freq)[0][0]
        iy = np.nonzero(self.y_cart[0, :, 0] == y)[0][0]
        iz = np.nonzero(self.z_cart[0, 0, :] == z)[0][0]
        iheight = np.nonzero(self.height == height)[0][0]
        plt.plot(self.x_cart[:, iy, iz],
                 self.deltaL_cart[:, iy, iz, ifreq, iheight],
                 **kwargs)

    def plot_tau(self, tau, z, freq, height, Nmax=None, **kwargs):
        ifreq = np.nonzero(self.frequencies == freq)[0][0]
        iheight = np.nonzero(self.height == height)[0][0]
        iz = np.nonzero(self.z_cart[0, 0, :] == z)[0][0]
        itau = np.nonzero(self.tau == tau)[0][0]
        # this is used in order to not plot element that are
        # outside the computational domain
        if Nmax is not None:
            self.deltaL_cart[Nmax:, itau, iz, ifreq, iheight] = np.nan
        plt.plot(self.x_cart[:, itau, iz],
                 self.deltaL_cart[:, itau, iz, ifreq, iheight],
                 **kwargs)
