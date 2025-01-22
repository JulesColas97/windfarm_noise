import numpy as np


class Mesh():
    """
    A class to create, manipulate, and transform spatial meshes in Cartesian
    and polar coordinate systems.

    The Mesh class provides methods to generate Cartesian and polar meshes,
    transform between coordinate systems, add topography
    or curvilinear features, and apply shifts.
    """

    polar: bool
    xmin: float
    xmax: float
    zmin: float
    zmax: float
    taumin: float
    taumax: float
    ymin: float
    ymax: float

    dx: float
    dy: float
    dz: float
    dtau: float

    # 1D array that define the regular grid
    x_array: np.ndarray
    y_array: np.ndarray
    z_array: np.ndarray
    tau_array: np.ndarray

    # 1D array of size (nx*ny*nz) or (nx*nz*tau)
    # used afterward in class Source() to compute the Spp
    x_coord: np.ndarray
    y_coord: np.ndarray
    z_coord: np.ndarray
    tau_coord: np.ndarray

    def __init__(self, polar):
        self.polar = polar

    def create_polar_mesh(self, xmin, xmax, dx, zmin, zmax, dz, taumin, taumax, dtau):
        """
        Create a polar mesh with specified spatial and angular ranges.

        Args:
            xmin (float): Minimum x value.
            xmax (float): Maximum x value.
            dx (float): Step size for x.
            zmin (float): Minimum z value.
            zmax (float): Maximum z value.
            dz (float): Step size for z.
            taumin (float): Minimum tau value (angle in radians).
            taumax (float): Maximum tau value (angle in radians).
            dtau (float): Step size for tau.
        """
        print('creating polar mesh')
        self.xmin = xmin
        self.xmax = xmax

        self.zmin = zmin
        self.zmax = zmax

        self.taumin = taumin
        self.taumax = taumax

        self.dx = dx
        self.dz = dz
        self.dtau = dtau

        self.x_array = np.arange(xmin, xmax, dx)
        self.tau_array = np.arange(taumin, taumax, dtau)
        self.z_array = np.arange(zmin, zmax, dz)

        self.nx = len(self.x_array)
        self.nz = len(self.z_array)
        self.ntau = len(self.tau_array)

        self.x_coord, self.z_coord, self.tau_coord,  = np.meshgrid(
            self.x_array, self.z_array, self.tau_array, indexing='ij')

        self.x_coord = self.x_coord.reshape(self.nx*self.ntau*self.nz, 1)
        self.z_coord = self.z_coord.reshape(self.nx*self.ntau*self.nz, 1)
        self.tau_coord = self.tau_coord.reshape(self.nx*self.ntau*self.nz, 1)

        self.x_coord[self.x_coord == 0] = 0.00001
        self.z_coord[self.z_coord == 0] = 0.00001

    def create_cartesian_mesh(self, xmin, xmax, nx, ymin, ymax, ny, zmin, zmax, nz):
        """
        Create a Cartesian mesh with specified ranges and dimensions.

        Args:
            xmin (float): Minimum x value.
            xmax (float): Maximum x value.
            nx (int): Number of points along x-axis.
            ymin (float): Minimum y value.
            ymax (float): Maximum y value.
            ny (int): Number of points along y-axis.
            zmin (float): Minimum z value.
            zmax (float): Maximum z value.
            nz (int): Number of points along z-axis.
        """
        print('creating cartesian mesh')
        self.xmin = xmin
        self.xmax = xmax
        self.nx = nx
        self.ymin = ymin
        self.ymax = ymax
        self.ny = ny
        self.zmin = zmin
        self.zmax = zmax
        self.nz = nz

        self.x_array = np.linspace(xmin, xmax, nx)
        self.y_array = np.linspace(ymin, ymax, ny)
        self.z_array = np.linspace(zmin, zmax, nz)

        self.x_coord, self.y_coord, self.z_coord, = np.meshgrid(
            self.x_array, self.y_array, self.z_array, indexing='ij')

        self.x_coord = self.x_coord.reshape(nx*ny*nz, 1)
        self.y_coord = self.y_coord.reshape(nx*ny*nz, 1)
        self.z_coord = self.z_coord.reshape(nx*ny*nz, 1)

        self.x_coord[self.x_coord == 0] = 0.0001
        self.y_coord[self.y_coord == 0] = 0.0001
        self.z_coord[self.z_coord == 0] = 0.0001
        return

    def set_polar_mesh(self, x, tau, z):
        """
        Set a polar mesh using predefined coordinate arrays.

        Args:
            x (np.ndarray): Array of x coordinates.
            tau (np.ndarray): Array of tau coordinates (angles).
            z (np.ndarray): Array of z coordinates.
        """
        self.x_array = x
        self.tau_array = tau
        self.z_array = z

        self.xmin = x[0]
        self.xmax = x[-1]
        self.nx = len(x)
        self.taumin = tau[0]
        self.taumax = tau[-1]
        self.ntau = len(tau)
        self.zmin = z[0]
        self.zmax = z[-1]
        self.nz = len(z)

        self.x_coord, self.z_coord, self.tau_coord, = np.meshgrid(
            self.x_array, self.z_array, self.tau_array, indexing='ij')
        self.y_coord = None
        self.x_coord = self.x_coord.reshape(self.nx*self.ntau*self.nz, 1)
        self.tau_coord = self.tau_coord.reshape(self.nx*self.ntau*self.nz, 1)
        self.z_coord = self.z_coord.reshape(self.nx*self.ntau*self.nz, 1)

        self.x_coord[self.x_coord == 0] = 0.0001
        self.y_coord[self.y_coord == 0] = 0.0001
        self.z_coord[self.z_coord == 0] = 0.0001

    def set_cartesian_mesh(self, x, y, z):
        """
        Set a Cartesian mesh using predefined coordinate arrays.

        Args:
            x (np.ndarray): Array of x coordinates.
            y (np.ndarray): Array of y coordinates.
            z (np.ndarray): Array of z coordinates.
        """
        self.x_array = x
        self.y_array = y
        self.z_array = z

        self.xmin = x[0]
        self.xmax = x[-1]
        self.nx = len(x)
        self.ymin = y[0]
        self.ymax = y[-1]
        self.ny = len(y)
        self.zmin = z[0]
        self.zmax = z[-1]
        self.nz = len(z)

        self.x_coord, self.z_coord, self.y_coord, = np.meshgrid(
            self.x_array, self.z_array, self.y_array, indexing='ij')
        self.tau_coord = None
        self.x_coord = self.x_coord.reshape(self.nx*self.ny*self.nz, 1)
        self.y_coord = self.y_coord.reshape(self.nx*self.ny*self.nz, 1)
        self.z_coord = self.z_coord.reshape(self.nx*self.ny*self.nz, 1)

        self.x_coord[self.x_coord == 0] = 0.0001
        self.y_coord[self.y_coord == 0] = 0.0001
        self.z_coord[self.z_coord == 0] = 0.0001

    def polar_2_curvilinear(self, hpos, hlength, hheight):
        """
        Transform polar mesh to curvilinear coordinates with a 2D hill.
        Usully used to create the mesh for a 2D xz plane

        Args:
            hpos (float): Position of the hilltop.
            hlength (float): Length of the hill.
            hheight (float): Height of the hill.
        """
        self.hpos = hpos
        self.hlength = hlength
        self.hheight = hheight
        self.z_coord = self.z_coord.reshape(self.nx, self.nz, self.ntau)
        for itau in range(self.ntau):
            m1 = (self.x_array*np.cos(self.tau_array[itau]) > (hpos-hlength)) * (
                self.x_array*np.cos(self.tau_array[itau]) < (hpos+hlength))
            H = np.zeros(self.nx)
            H[m1] = hheight * np.cos(np.pi * (self.x_array[m1]
                                     * np.cos(self.tau_array[itau])-hpos)/(2*hlength))**2
            self.z_coord[:, :, itau] = H.reshape(-1, 1) + self.z_coord[:, :, itau]*(
                self.zmax - H.reshape(-1, 1))/self.zmax

        self.z_coord = self.z_coord.reshape(self.nx*self.nz*self.ntau, 1)

    def add_topography(self, hpos, hlength, hheight):
        """
        Add topography to the mesh. All z coordinates are modified similarly in function of the 2D hill.
        Usually used to create the mesh for the 2D xy plane, where receiver are at 2m height.

        Args:
            hpos (float): Position of the topography peak.
            hlength (float): Length of the topography.
            hheight (float): Height of the topography.
        """
        self.z_coord = self.z_coord.reshape(self.nx, self.nz, self.ntau)
        for itau in range(self.ntau):
            m1 = ((self.x_array*np.cos(self.tau_array[itau]) > (hpos-hlength))
                  * (self.x_array*np.cos(self.tau_array[itau]) < (hpos+hlength)))
            H = np.zeros(self.nx)
            H[m1] = hheight * np.cos(np.pi * (self.x_array[m1]
                                              * np.cos(self.tau_array[itau])-hpos)/(2*hlength))**2
            self.z_coord[:, :, itau] = (H.reshape(-1, 1)
                                        + self.z_coord[:, :, itau])

        self.z_coord = self.z_coord.reshape(self.nx*self.nz*self.ntau, 1)

    def cartesian_2_polar(self):
        """
        Convert the Cartesian mesh to polar coordinates.
        """
        self.tau_coord = (np.arctan2(self.y_coord, self.x_coord))
        self.x_coord = np.sqrt(self.x_coord ** 2 + self.y_coord ** 2)
        self.y_coord = None

    def polar_2_cartesian(self):
        """
        Convert the polar mesh to Cartesian coordinates.
        """
        self.y_coord = np.sin(self.tau_coord)*self.x_coord
        self.x_coord = np.cos(self.tau_coord)*self.x_coord
        self.tau_coord = None

    def shift(self, x, y):
        """
        Shift the mesh in Cartesian coordinates.
        This is usefull to create first a general cartesian mesh.
        Then shift it to put the wind turbine at position (0,0).
        Compute the Spp using Source class.
        And then shift it again.

        Args:
            x (float): Shift along x-axis.
            y (float): Shift along y-axis.
        """
        if self.tau_coord is not None:
            print('warning not in cartesian coordinate')
            return

        self.x_array = self.x_array + x
        self.y_array = self.y_array + y
        self.z_array = self.z_array

        self.xmin = self.x_array[0]
        self.xmax = self.x_array[-1]
        self.nx = len(self.x_array)
        self.ymin = self.y_array[0]
        self.ymax = self.y_array[-1]
        self.ny = len(self.y_array)

        self.x_coord, self.z_coord, self.y_coord, = np.meshgrid(
            self.x_array, self.z_array, self.y_array, indexing='ij')

        self.x_coord = self.x_coord.reshape(self.nx*self.ny*self.nz, 1)
        self.y_coord = self.y_coord.reshape(self.nx*self.ny*self.nz, 1)
        self.z_coord = self.z_coord.reshape(self.nx*self.ny*self.nz, 1)

        self.x_coord[self.x_coord == 0] = 0.0001
        self.y_coord[self.y_coord == 0] = 0.0001
        self.z_coord[self.z_coord == 0] = 0.0001
