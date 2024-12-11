import numpy as np

class Mesh():
    polar:bool
    def __init__(self,polar):
        self.polar = polar

    def create_polar_mesh(self,xmin,xmax,dx,zmin,zmax,dz,taumin,taumax,dtau):
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


        self.x_array =np.arange(xmin,xmax,dx)
        self.tau_array =np.arange(taumin,taumax,dtau)
        self.z_array =np.arange(zmin,zmax,dz)

        self.nx = len(self.x_array)
        self.nz =  len(self.z_array)
        self.ntau =  len(self.tau_array)

        self.x_coord, self.z_coord, self.tau_coord,  = np.meshgrid(self.x_array, self.z_array, self.tau_array,indexing='ij')

        self.x_coord = self.x_coord.reshape(self.nx*self.ntau*self.nz, 1)
        self.z_coord = self.z_coord.reshape(self.nx*self.ntau*self.nz, 1)
        self.tau_coord = self.tau_coord.reshape(self.nx*self.ntau*self.nz, 1)

        self.x_coord[self.x_coord == 0] = 0.00001
        self.z_coord[self.z_coord == 0] = 0.00001


    def create_cartesian_mesh(self,xmin,xmax,nx,ymin,ymax,ny,zmin,zmax,nz):
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

        self.x_array =np.linspace(xmin,xmax,nx)
        self.y_array =np.linspace(ymin,ymax,ny)
        self.z_array =np.linspace(zmin,zmax,nz)

        self.x_coord, self.y_coord, self.z_coord, = np.meshgrid(self.x_array, self.y_array, self.z_array,indexing='ij')

        self.x_coord = self.x_coord.reshape(nx*ny*nz, 1)
        self.y_coord = self.y_coord.reshape(nx*ny*nz, 1)
        self.z_coord = self.z_coord.reshape(nx*ny*nz, 1)

        self.x_coord[self.x_coord == 0] = 0.0001
        self.y_coord[self.y_coord == 0] = 0.0001
        self.z_coord[self.z_coord == 0] = 0.0001
        return

    def set_polar_mesh(self,x,tau,z):
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

        self.x_coord, self.z_coord, self.tau_coord, = np.meshgrid(self.x_array, self.z_array, self.tau_array, indexing='ij')
        self.y_coord = None
        self.x_coord = self.x_coord.reshape(self.nx*self.ntau*self.nz, 1)
        self.tau_coord = self.tau_coord.reshape(self.nx*self.ntau*self.nz, 1)
        self.z_coord = self.z_coord.reshape(self.nx*self.ntau*self.nz, 1)

        self.x_coord[self.x_coord == 0] = 0.0001
        self.y_coord[self.y_coord == 0] = 0.0001
        self.z_coord[self.z_coord == 0] = 0.0001

    def set_cartesian_mesh(self,x,y,z):
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

        self.x_coord, self.z_coord, self.y_coord, = np.meshgrid(self.x_array, self.z_array, self.y_array, indexing='ij')
        self.tau_coord = None
        self.x_coord = self.x_coord.reshape(self.nx*self.ny*self.nz, 1)
        self.y_coord = self.y_coord.reshape(self.nx*self.ny*self.nz, 1)
        self.z_coord = self.z_coord.reshape(self.nx*self.ny*self.nz, 1)

        self.x_coord[self.x_coord == 0] = 0.0001
        self.y_coord[self.y_coord == 0] = 0.0001
        self.z_coord[self.z_coord == 0] = 0.0001

    def polar_2_curvilinear(self,hpos,hlength,hheight):
        self.hpos = hpos
        self.hlength = hlength
        self.hheight = hheight
        self.z_coord = self.z_coord.reshape(self.nx,self.nz, self.ntau)
        for itau in range(self.ntau) :
            m1 = (self.x_array*np.cos(self.tau_array[itau])>(hpos-hlength)) * (self.x_array*np.cos(self.tau_array[itau])<(hpos+hlength))
            H = np.zeros(self.nx)
            H[m1] = hheight * np.cos (np.pi *(self.x_array[m1]*np.cos(self.tau_array[itau])-hpos)/(2*hlength))**2
            self.z_coord[:,:,itau] = H.reshape(-1,1) + self.z_coord[:,:,itau]*(self.zmax - H.reshape(-1,1))/self.zmax

        self.z_coord = self.z_coord.reshape(self.nx*self.nz*self.ntau,1)


    def add_topography(self,hpos,hlength,hheight):
        self.z_coord = self.z_coord.reshape(self.nx,self.nz, self.ntau)
        for itau in range(self.ntau) :
            m1 = (self.x_array*np.cos(self.tau_array[itau])>(hpos-hlength)) * (self.x_array*np.cos(self.tau_array[itau])<(hpos+hlength))
            H = np.zeros(self.nx)
            H[m1] = hheight * np.cos (np.pi *(self.x_array[m1]*np.cos(self.tau_array[itau])-hpos)/(2*hlength))**2
            self.z_coord[:,:,itau] = H.reshape(-1,1) + self.z_coord[:,:,itau]

        self.z_coord = self.z_coord.reshape(self.nx*self.nz*self.ntau,1)

    def cartesian_2_polar(self):
        self.tau_coord = (np.arctan2(self.y_coord, self.x_coord))
        self.x_coord = np.sqrt(self.x_coord ** 2 + self.y_coord ** 2)
        self.y_coord = None

    def polar_2_cartesian(self):
        self.y_coord = np.sin(self.tau_coord)*self.x_coord
        self.x_coord = np.cos(self.tau_coord)*self.x_coord
        self.tau_coord = None

    def shift(self,x,y):
        if self.tau_coord is not None:
            print('warning not in cartesian coordinate')
            return

        self.x_array = self.x_array + x
        self.y_array = self.y_array + y
        self.z_array = self.z_array

        self.xmin = self.x_array [0]
        self.xmax = self.x_array [-1]
        self.nx = len(self.x_array )
        self.ymin = self.y_array[0]
        self.ymax = self.y_array[-1]
        self.ny = len(self.y_array)

        self.x_coord, self.z_coord, self.y_coord, = np.meshgrid(self.x_array, self.z_array, self.y_array, indexing='ij')

        self.x_coord = self.x_coord.reshape(self.nx*self.ny*self.nz, 1)
        self.y_coord = self.y_coord.reshape(self.nx*self.ny*self.nz, 1)
        self.z_coord = self.z_coord.reshape(self.nx*self.ny*self.nz, 1)

        self.x_coord[self.x_coord == 0] = 0.0001
        self.y_coord[self.y_coord == 0] = 0.0001
        self.z_coord[self.z_coord == 0] = 0.0001

