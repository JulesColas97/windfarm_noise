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
    """Class to handle PE results from an H5 file.

    Args:
        casename (str): Name of the case.
        iTurb (int): Index of the turbine.
        height (float): Source height.
        tau (float): Propagation angle.
        dirname (str, optional): Directory name. Defaults to "./".
        distribute_tau (int, optional): Number of tau distributions. Defaults to None.
    """
    def __init__(self, casename: str, iTurb: int, height: float, tau: float, dirname: str = "./", distribute_tau: int = None):
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

    def read_carto(self, fname: str):
        """Reads the cartography data from the H5 file.

        Args:
            fname (str): File name of the H5 file.
        """
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

    def read_receiver(self, fname: str):
        """Reads the receiver data from the H5 file.

        Args:
            fname (str): File name of the H5 file.
        """
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

    def read_planes(self, fname: str):
        """Reads the plane data from the H5 file.

        Args:
            fname (str): File name of the H5 file.
        """
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

    def plotSide(self, freq: float, cmap: str = "RdBu_r"):
        """Plots the side view of the deltaL for a given frequency.

        Args:
            freq (float): Frequency to plot.
            cmap (str, optional): Colormap to use. Defaults to "RdBu_r".
        """
        ifreq = np.nonzero(self.frequencies == freq)[0][0]
        if ifreq is None:
            print("frequence was not calculated.")
            return
        plt.pcolormesh(
            self.x, self.z, self.deltaL[:, :, ifreq].T, cmap=cmap, shading="gouraud"
        )

    def plotSide3octave(self, freq: float, cmap: str = "RdBu_r"):
        """Plots the side view of the deltaL for a given third-octave frequency.

        Args:
            freq (float): Third-octave frequency to plot.
            cmap (str, optional): Colormap to use. Defaults to "RdBu_r".
        """
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
        plt.ylim(0, 300)

    def plotLineTotalDeltaL(self, height: float):
        """Plots the total deltaL along a line at a given height.

        Args:
            height (float): Height to plot.
        """
        iheight = np.nonzero(self.heights == height)[0][0]
        plt.plot(
            self.x,
            10 * np.log10(np.sum(10 ** (self.deltaL3octave[:, iheight, :] / 10), 1)),
        )

    def plotSideTotalDeltaL(self, cmap: str = "RdBu_r"):
        """Plots the side view of the total deltaL.

        Args:
            cmap (str, optional): Colormap to use. Defaults to "RdBu_r".
        """
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

    def plotLine(self, freq: float, z: float, **kwargs):
        """Plots the deltaL along a line at a given frequency and height.

        Args:
            freq (float): Frequency to plot.
            z (float): Height to plot.
            **kwargs: Additional keyword arguments for plotting.
        """
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

    def plotLine3octave(self, freq: float, height: float, **kwargs):
        """Plots the deltaL along a line at a given third-octave frequency and height.

        Args:
            freq (float): Third-octave frequency to plot.
            height (float): Height to plot.
            **kwargs: Additional keyword arguments for plotting.
        """
        ifreq = np.nonzero(self.fc == freq)[0][0]
        iheight = np.nonzero(self.heights == height)[0][0]
        if ifreq is None:
            print("frequence was not calculated.")
            return
        if iheight is None:
            print("receiver's height was not calculated.")
            return
        plt.plot(self.x, self.receiver3octave[:, iheight, ifreq],**kwargs)

    def compute3octave(self, fc: list = None, Nfc: list = None):
        """Computes the third-octave frequencies.

        Args:
            fc (list, optional): List of center frequencies. Defaults to None.
            Nfc (list, optional): List of the number of frequencies per band. Defaults to None.
        """
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

    def plotLineOaspl(self, height: float):
        """Plots the OASPL along a line at a given height.

        Args:
            height (float): Height to plot.
        """
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

    def plotFlowSide(self, fname: str, cmap: str = "RdYlBu_r", xstep: int = 1, zstep: int = 1):
        """Plots the side view of the flow field from a given H5 file.

        Args:
            fname (str): File name of the H5 file containing the flow data.
            cmap (str, optional): Colormap to use. Defaults to "RdYlBu_r".
            xstep (int, optional): Step size in the x-direction. Defaults to 1.
            zstep (int, optional): Step size in the z-direction. Defaults to 1.
        """
        f = h5py.File(fname, "r")

        # read mesh
        x = np.array(f["x"])
        z = np.array(f["z"])

        # read flow
        u = np.array(f["u"])
        c = np.array(f["c"])
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
