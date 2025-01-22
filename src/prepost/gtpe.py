import numpy as np
import h5py
import matplotlib.pyplot as plt
from .utils import computeThirdOctaveFrequencies, chkList, interp_weights, interpolate
from mpl_toolkits.axes_grid1 import make_axes_locatable


class GtpeResults():
    """Class to handle GTPE results from an H5 file.

    Args:
        fname (str): File name of the H5 file.
        tau (float): Propagation angle.
        height (float): Source height.
    """
    def __init__(self, fname, tau, height):
        self.fname = fname
        self.tau = tau
        self.height = height

    # read solution from h5 file
    # read carto and receivers
    def read(self):
        """Reads the solution from the H5 file, including mesh and receivers.
        """
        fname = self.fname
        file = h5py.File(fname, "r")
        solutions = file.get('solution')
        # read mesh and create mesh grid
        mesh = file.get('mesh')
        frequencies = list(solutions.keys())
        self.frequencies = np.array([float(i) for i in frequencies])
        self.x = np.array(mesh['x'])
        self.z = np.array(mesh['y'])
        self.h = np.array(mesh['H'])
        [self.X, self.Z] = np.meshgrid(self.x, self.z)
        self.Z = self.Z*(self.z[len(self.z)-1]-self.h) / \
            self.z[len(self.z)-1] + self.h

        # read carto of delta L
        self.deltaL = np.zeros(
            (len(self.x), len(self.z), len(self.frequencies)))
        # solution matrix (Nz,Nx,Nf)
        for ii in range(len(frequencies)):
            sol = solutions.get(frequencies[ii])
            self.deltaL[:, :, ii] = np.transpose(np.array(sol['deltaL']))
        # read delta L saved at receivers
        receivers_all = file.get('receiver')
        self.heights = np.array(receivers_all.get('heights'))
        self.receiver = np.zeros(
            (len(self.x), len(self.heights), len(self.frequencies)))
        for ii in range(len(self.frequencies)):
            rec = receivers_all.get(frequencies[ii])
            self.receiver[:, :, ii] = np.transpose(np.array(rec['deltaL']))

    def plotSide(self, freq: float, cmap: str = 'RdYlBu_r'):
        """Plots the side view of the deltaL for a given frequency.

        Args:
            freq (float): Frequency to plot.
            cmap (str, optional): Colormap to use. Defaults to 'RdYlBu_r'.
        """
        ifreq = np.nonzero(self.frequencies == freq)[0][0]
        if ifreq is None:
            print('frequence was not calculated.')
            return
        plt.pcolormesh(
            self.X, self.Z, self.deltaL[:, :, ifreq].T, cmap=cmap, shading='gouraud')
        plt.clim(-10, 10)
        plt.gca().set_aspect('equal', adjustable='box')
        plt.ylim(0, 300)

    def plotSide3octave(self, freq: float, cmap: str = 'RdYlBu_r'):
        """Plots the side view of the deltaL for a given third-octave frequency.

        Args:
            freq (float): Third-octave frequency to plot.
            cmap (str, optional): Colormap to use. Defaults to 'RdYlBu_r'.
        """
        ifreq = np.nonzero(self.fc == freq)[0][0]
        if ifreq is None:
            print('frequence was not calculated.')
            return
        plt.pcolormesh(
            self.X, self.Z, self.deltaL3octave[:, :, ifreq].T, cmap=cmap, shading='gouraud')
        plt.clim(-10, 10)
        plt.gca().set_aspect('equal', adjustable='box')
        plt.ylim(0, 300)

    def plotLine(self, freq: float, height: float):
        """Plots the deltaL along a line at a given frequency and height.

        Args:
            freq (float): Frequency to plot.
            height (float): Height to plot.
        """
        ifreq = np.nonzero(self.frequencies == freq)[0][0]
        iheight = np.nonzero(self.heights == height)[0][0]
        print(self.heights[iheight])
        if ifreq is None:
            print('frequence was not calculated.')
            return
        plt.plot(self.x, self.receiver[:, iheight, ifreq])

    def plotLine3octave(self, freq: float, height: float):
        """Plots the deltaL along a line at a given third-octave frequency and height.

        Args:
            freq (float): Third-octave frequency to plot.
            height (float): Height to plot.
        """
        ifreq = np.nonzero(self.fc == freq)[0][0]
        iheight = np.nonzero(self.heights == height)[0][0]
        if ifreq is None:
            print('frequence was not calculated.')
            return
        plt.plot(self.x, self.receiver3octave[:, iheight, ifreq])

    def compute3octave(self, fc: list, Nfc: list):
        """Computes the third-octave frequencies.

        Args:
            fc (list): List of center frequencies.
            Nfc (list): List of the number of frequencies per band.
        """
        freq = computeThirdOctaveFrequencies(fc, Nfc)
        if not np.all(freq == self.frequencies):
            print('central frequencies are not the same')
            print(freq)
            print(self.frequencies)
            return
        self.fc = np.array(fc)
        self.Nfc = np.array(Nfc)
        self.deltaL3octave = np.zeros((len(self.x), len(self.z), len(fc)))
        self.receiver3octave = np.zeros(
            (len(self.x), len(self.heights), len(fc)))
        for ifc in range(len(fc)):
            self.deltaL3octave[:, :, ifc] = 10*np.log10(np.sum(10**(self.deltaL[:, :, int(
                np.sum(self.Nfc[0:ifc])):np.sum(self.Nfc[0:ifc+1])]/10), 2)/Nfc[ifc])
            self.receiver3octave[:, :, ifc] = 10*np.log10(np.sum(10**(self.receiver[:, :, int(
                np.sum(self.Nfc[0:ifc])):np.sum(self.Nfc[0:ifc+1])]/10), 2)/Nfc[ifc])

    def plotLineOaspl(self, height: float):
        """Plots the OASPL along a line at a given height.

        Args:
            height (float): Height to plot.
        """
        iheight = np.argmin(abs(self.heights - height))
        SPL, SPL_A, OASPL, OASPL_A = computeSPLLine(np.squeeze(self.receiver3octave[iheight, :, :]), np.squeeze(self.x),
                                                    np.squeeze(self.heights[iheight]+self.h), self.fc, self.tau, self.height, c0=343)
        plt.plot(self.x, OASPL_A)

    def plotFlowSide(self, fname: str, cmap: str = 'RdYlBu_r', xstep: int = 1, zstep: int = 1):
        """Plots the side view of the flow field from a given H5 file.

        Args:
            fname (str): File name of the H5 file containing the flow data.
            cmap (str, optional): Colormap to use. Defaults to 'RdYlBu_r'.
            xstep (int, optional): Step size in the x-direction. Defaults to 1.
            zstep (int, optional): Step size in the z-direction. Defaults to 1.
        """
        f = h5py.File(fname, "r")

        # read mesh
        x = np.array(f['x'])
        z = np.array(f['z'])

        # read flow
        u = np.array(f['u'])
        v = np.array(f['v'])

        print(x.shape)
        print(z.shape)
        fig, ax = plt.subplots(figsize=(10, 1.8))

        cax = ax.pcolormesh(x[::xstep, ::zstep], z[::xstep, ::zstep],
                            u[::xstep, ::zstep], cmap=cmap, shading='gouraud')
        ax.set_aspect('equal', adjustable='box')
        divider = make_axes_locatable(ax)
        cbax = divider.append_axes("right", size="2%", pad=0.05)
        ax.set_xlabel('$x$ (m)')
        ax.set_ylabel('$z$ (m)')
        fig.colorbar(cax, cax=cbax, label='$u$ (m/s)')
        plt.tight_layout()
