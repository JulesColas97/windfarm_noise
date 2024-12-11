import numpy as np
import h5py
import matplotlib.pyplot as plt
import os
import shutil
import pickle

from .les import Les
from .utils import computeThirdOctaveFrequencies, chkList, interp_weights, interpolate

# class to init and launch PE simulation for wind farm noise propagation 
class Simu():
    
    # initialize the input dictionnary for the PE simulation 
    def __init__(self,casename: str):
        self.casename=casename
        self.inputs = {
            "case_name"               : str(casename),
            "var0%c"                  : 343,
            "var0%rho"                : 1.2,
            "var0%gamma"              : 1.4,

            "nb_freq"                 : 3,
            "frequencies(1)"          : [100,500,1000],

            "nb_theta"                : 1,
            "theta(1)"                : 0,

            "Lx1"                     : 3000.,
            "Lx2"                     : 1000.,
            "Ly"                      : 1000.,
            "Lz"                      : 300.,
            "dx"                      : 0.5,
            "cfl" 	                  : 0.1,

            "size"                    : 30,
            "param"                   : 5e4,
            "n"                       : 2.5,

            "imp%sigmae"              : 50e3,
            "imp%alphae"              : 100.,
            "imp%rigid"               : False,

            "src%pos_x"               : 0.,
            "src%pos_z"               : 100.,

            "external_flow"           : True,
            "uniform"                 : False,
            "logarithmic"             : False,
            "u0"                      : 0,
            "arbitrary"                : True,
            "arbitrary_new"           : False,
            "interpolation"           : True,

            "fdir"                    : '/2Dhill/caseA_flat/blue/output/',
            "fname"                   : 'tavg',

            "tinput%ratio"            : 1,
            "tinput%z_i"              : 625.,
            "tinput%delta"            : 0.25,
            "tinput%Lx"               : 10.,
            "tinput%Ly"               : 1.,
            "tinput%Lz"               : 1.,

            "tinput%posx"             : 1.875,
            "tinput%posy"             : 0.5,
            "tinput%posz"             : 0.04,


            "dout"                    : 1,
            "nb_receiver"             : 3,
            "heights(1)"              : [2.,10,100],
            "side"                    : True,
            "top"                     : False,
            "continuation"            : False,

        }

    # read data from LES 
    def readLes(self,path):
        self.les = Les(path)
        self.les.read()
        print('les domain dimension:')
        print(self.les.L_x,self.les.L_y,self.les.L_z)
        print(str(self.les.turbines.shape[0])+' wind turbines :')
        print(self.les.turbines)

    # create input.nml files from inputs dictionnary
    def createInputFile(self,directory='./'):
        f= open(directory + '/input.nml','w')
        f.write("$input\n")
        for key,value in self.inputs.items():
            if type(value)==str:
                f.write("%s='%s'\n"%(key,value))
            elif type(value)==bool :
                if value :
                    f.write("%s=%s\n"%(key,'.true.'))
                else :
                    f.write("%s=%s\n"%(key,'.false.'))
            elif isinstance(value, list):
                f.write("%s=%s\n"%(key,str(value)[1:-1]))
            else :
                f.write("%s=%s\n"%(key,str(value)))
        f.write('$end input')
        f.close()
        return

    # set an input in inputs dictionnary
    def set_input(self,key,value):
        self.inputs[key] = value

    # set frequencies 
    # if Fc and Nfc fc are the central frequency and Nfc are the number of frequencies to be computed er band
    # if only fc is given only these frequencies are computed
    def set_frequencies(self,fc,Nfc=None):
        if Nfc is None:
            self.inputs['nb_freq'] = len(fc)
            self.inputs['frequencies(1)'] = list(fc)
        else : 
            self.inputs['nb_freq'] = sum(Nfc)
            self.inputs['frequencies(1)'] = list(computeThirdOctaveFrequencies(fc,Nfc))

    # define the case to be simulate with number of source heights, angles and flow data 
    def defineCases(self,heights,angles,src_path='/home/lmfa/jcolas/Documents/DEV/wf_phd/src/kernel/New_PE/PE_2D_WAPE',
                                        flow_path='/home/lmfa/jcolas/Documents/DEV/LES/'):
        self.heights=heights
        self.angles = angles
        self.src_path = src_path

        self.inputs['fdir'] = flow_path+'output/'
        print(self.inputs['fdir'])
        self.les = Les(flow_path)
        self.les.read()
        self.tx = self.les.turbines[:,1]
        self.ty = self.les.turbines[:,2]

        self.inputs['tinput%z_i'] = self.les.z_i
        self.inputs['tinput%Lx'] = self.les.lx 
        self.inputs['tinput%Ly'] = self.les.ly 
        self.inputs['tinput%Lz'] = self.les.lz 

    # define the total domain
    def defineDomain(self,x1,x2,y1,y2,h):
        self.x1 = x1
        self.x2 = x2
        self.y1 = y1
        self.y2 = y2
        self.h = h 

    # set the specific domain for each turbine 
    def set_domain(self,iturbine):
        Lx1 = self.x2-self.tx[iturbine]*self.les.z_i
        Lx2 = self.tx[iturbine]*self.les.z_i - self.x1 
        Ly1 = self.y2 - self.ty[iturbine]*self.les.z_i
        Ly2 = self.ty[iturbine]*self.les.z_i - self.y1
        Ly = max(Ly1,Ly2)
        Lz = self.h

        self.set_input('Lx1',Lx1)
        self.set_input('Lx2',Lx2)
        self.set_input('Ly',Ly)
        self.set_input('Lz',Lz)

    # create launch file for slurm 
    def createLaunch(self,dirname,jname,mem=5000,time="2:00:00") :
        f= open(dirname+'/launch.sh','w')
        f.write("#!/bin/bash\n")
        f.write("#SBATCH --job-name="+jname+"\n")
        f.write("#SBATCH --output=out.out # output messages go here\n")
        f.write("#SBATCH --error=err.err    # error messages go here\n")
        f.write("#SBATCH --mail-user=jules.colas@ecl17.ec-lyon.fr\n")
        f.write("#SBATCH --mail-type=ALL\n")
        f.write("#SBATCH --partition=haswell # partition name\n")
        f.write("#SBATCH --nodes=1\n")
        f.write("#SBATCH --cpus-per-task=1\n")
        f.write("#SBATCH --mem=" + str(mem) + "\n")
        f.write("#SBATCH --time=" + str(time) + "\n")


        f.write("module purge\n")
        f.write("module load  HDF5/1.10.1-intel-2018a\n")
        f.write("export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK\n")
        f.write("export KMP_AFFINITY=granularity=fine,compact,1,0\n")
        f.write("export OMP_STACKSIZE=1g\n")
        f.write("ulimit -s unlimited\n")

        f.write("time ./PE_2D\n")

    # create the directories inpu.nml, launch.sh files for running the complete simulation 
    # this wont launch the jobs 
    def distributeCases(self,distribute_tau=False,mem=5000,time="2:00:00"):
        dir = os.getcwd()
        self.distribute_tau = distribute_tau
        if distribute_tau==False: 
            self.inputs['nb_theta'] = len(self.angles)
            self.inputs['theta(1)'] = list(self.angles)

        for ii in range(len(self.tx)):
            self.inputs["tinput%posx"] = self.tx[ii]
            self.inputs["tinput%posy"] = self.ty[ii]
            self.set_domain(ii)
            for jj in range(len(self.heights)):
                self.inputs['src%pos_z'] = self.heights[jj]
                if distribute_tau:
                    for kk in range(len(self.angles)):
                        self.inputs['nb_theta'] = 1
                        self.inputs['theta(1)'] = self.angles[kk]

                        path = os.path.join(dir, 't'+str(ii)+'/'+str(self.heights[jj])+'/'+str(self.angles[kk])+'/')
                        os.makedirs(path,exist_ok=True)
                        self.createInputFile(directory=path)
                        self.createLaunch(path,self.casename+'t'+str(ii)+'h'+str(self.heights[jj])+'a'+str(self.angles[kk]),mem=mem,time=time)
                        shutil.copy(self.src_path,path+'/PE_2D')
                else : 
                    path = os.path.join(dir, 't'+str(ii)+'/'+str(self.heights[jj]))
                    os.makedirs(path,exist_ok=True)
                    self.createInputFile(directory=path)
                    self.createLaunch(path,self.casename+'t'+str(ii)+'h'+str(self.heights[jj]),mem=mem,time=time)
                    shutil.copy(self.src_path,path+'/PE_2D')

    # launch the simulation 
    # must be used after distribute cases 
    def launchCases(self):
        dir = os.getcwd()
        for ii in range(len(self.tx)):
            for jj in range(len(self.heights)):
                if self.distribute_tau:
                    for kk in range(len(self.angles)):
                        path = os.path.join(dir, 't'+str(ii)+'/'+str(self.heights[jj])+'/'+str(self.angles[kk])+'/')
                        os.system("cd "+path+"; sbatch launch.sh")
                else : 
                    path = os.path.join(dir, 't'+str(ii)+'/'+str(self.heights[jj]))  
                    os.system("cd "+path+"; sbatch launch.sh")

    
    def save(self):
        """save class as self.name.dat"""
        with open(self.casename+'.dat','wb') as file:
            pickle.dump(self.__dict__,file)

    def load(self):
        """try load self.name.dat"""
        with open(self.casename+'.dat','rb') as file:
            self.__dict__ = pickle.load(file)


class PeResults():
    def __init__(self,casename,iTurb,height,tau,dirname='./',distribute_tau=False):
        self.casename=casename
        self.iTurb = iTurb
        self.tau = tau
        if int(tau)==tau:
            self.tau_str = format(tau,'04d')
        else :
            self.tau_str = str(tau)

        self.height = height
        self.dirname = dirname 
        self.distribute_tau = distribute_tau
        
    def read(self):
        if self.distribute_tau:
            path = self.dirname + 't'+str(self.iTurb)+'/'+str(self.height)+'/'+self.tau_str+'/'
        else : 
            path = self.dirname + 't'+str(self.iTurb)+'/'+str(self.height)+'/'

        file = h5py.File(path+self.casename+'_'+self.tau_str+'.h5', "r")
        solutions = file.get('solution')
        mesh = file.get('mesh')
        frequencies = list(solutions.keys())
        self.frequencies= np.array([float(i) for i in frequencies])
        self.x =  np.array(mesh['x'])
        self.z =  np.array(mesh['y'])
        self.deltaL = np.zeros((len(self.x),len(self.z),len(self.frequencies)))
        # solution matrix (Nz,Nx,Nf)
        for ii in range(len(frequencies)):
            sol = solutions.get(frequencies[ii])
            self.deltaL[:,:,ii] = np.transpose(np.array(sol['deltaL']))

        # read delta L saved at receivers
        receivers_all= file.get('receiver')
        self.heights = np.array(receivers_all.get('heights'))
        self.receiver = np.zeros((len(self.x),len(self.heights),len(self.frequencies)))
        for ii in range(len(self.frequencies)) :
            rec = receivers_all.get(frequencies[ii])
            self.receiver[:,:,ii] = np.transpose(np.array(rec['deltaL']))

    def plotSide(self,freq,cmap='RdBu_r'):
        ifreq = np.nonzero(self.frequencies==freq)[0][0]
        if ifreq is None:
            print('frequence was not calculated.')
            return
        plt.pcolormesh(self.x,self.z,self.deltaL[:,:,ifreq].T,cmap=cmap,shading='gouraud')
        plt.clim(-10,10)
        plt.gca().set_aspect('equal', adjustable='box')
        plt.ylim(0,300)

    def plotSide3octave(self,freq,cmap='RdBu_r'):
        ifreq = np.nonzero(self.fc==freq)[0][0]
        if ifreq is None:
            print('frequence was not calculated.')
            return
        plt.pcolormesh(self.x,self.z,self.deltaL3octave[:,:,ifreq].T,cmap=cmap,shading='gouraud')
        plt.clim(-10,10)
        plt.gca().set_aspect('equal', adjustable='box')
        plt.ylim(0,300)

    def plotLineTotalDeltaL(self,height):
        iheight = np.nonzero(self.heights==height)[0][0]
        plt.plot(self.x,10*np.log10(np.sum(10**(self.deltaL3octave[:,iheight,:]/10),1)))

    def plotSideTotalDeltaL(self,cmap='RdBu_r'):
        plt.pcolormesh(self.x,self.z,10*np.log10(np.sum(10**(self.deltaL3octave[:,:,:]/10),2)).T,cmap=cmap,shading='gouraud')
        plt.clim(-10,10)
        plt.gca().set_aspect('equal', adjustable='box')
        plt.ylim(0,300)

    def plotLine(self,freq,height):
        ifreq = np.nonzero(self.frequencies==freq)[0][0]
        iheight = np.nonzero(self.heights==height)[0][0]
        if ifreq is None:
            print('frequence was not calculated.')
            return
        if iheight is None:
            print("receiver's height was not calculated.")
            return
        plt.plot(self.x,self.receiver[:,iheight,ifreq])

    def plotLine3octave(self,freq,height):
        ifreq = np.nonzero(self.fc==freq)[0][0]
        iheight = np.nonzero(self.heights==height)[0][0]
        if ifreq is None:
            print('frequence was not calculated.')
            return
        if iheight is None:
            print("receiver's height was not calculated.")
            return
        plt.plot(self.x,self.receiver3octave[:,iheight,ifreq])

    def compute3octave(self,fc,Nfc):
        freq =  computeThirdOctaveFrequencies(fc,Nfc)
        if not np.all(freq==self.frequencies):
            print('central frequencies are not the same')
            return
        self.fc = np.array(fc)
        self.Nfc = np.array(Nfc)
        self.deltaL3octave = np.zeros((len(self.x),len(self.z),len(fc)))
        self.receiver3octave = np.zeros((len(self.x),len(self.heights),len(fc)))
        for ifc in range (len(fc)) :
            self.deltaL3octave[:,:,ifc]   = 10*np.log10(np.sum(10**(self.deltaL  [:,:,int(np.sum(self.Nfc[0:ifc])):np.sum(self.Nfc[0:ifc+1])]/10),2)/Nfc[ifc])
            self.receiver3octave[:,:,ifc] = 10*np.log10(np.sum(10**(self.receiver[:,:,int(np.sum(self.Nfc[0:ifc])):np.sum(self.Nfc[0:ifc+1])]/10),2)/Nfc[ifc])

    def plotLineOaspl(self,height):
        iheight = np.argmin(abs(self.heights - height))
        SPL, SPL_A,OASPL,OASPL_A = computeSPLLine(np.squeeze(self.receiver3octave[iheight,:,:]),np.squeeze(self.x),
                                                  np.squeeze(self.heights[iheight]+self.h),self.fc,self.tau,self.height,c0=343)
        plt.plot(self.x,OASPL_A)
   
    def plotFlowSide(self,fname,cmap='RdYlBu_r',xstep=1,zstep=1):
        f = h5py.File(fname, "r")

        # read mesh
        x =  np.array(f['x'])
        z =  np.array(f['z'])

        #read flow
        u =  np.array(f['u'])
        v =  np.array(f['v'])

        fig,ax = plt.subplots(figsize=(10, 1.8))
        cax = ax.pcolormesh(x[::xstep,::zstep],z[::xstep,::zstep],u[::xstep,::zstep],cmap=cmap,shading='gouraud')
        ax.set_aspect('equal', adjustable='box')
        divider = make_axes_locatable(ax)
        cbax = divider.append_axes("right", size="2%", pad=0.05)
        ax.set_xlabel('$x$ (m)')
        ax.set_ylabel('$z$ (m)')
        fig.colorbar(cax, cax=cbax,label='$u$ (m/s)')
        plt.tight_layout()


# define a deltaL field to use for the source model 
class DeltaLField():
    def __init__(self,dirname, casename):
        self.dirname=dirname
        self.casename=casename
        self.height = []
        self.tau = []
        self.frequencies = []
        self.deltaL = None
        self.deltaLlist = []
        self.hindex = 0
        self.tauindex = 0

    # read solution from file
    # one file corespond to one plane of propagation with one source altitude
    def read_carto(self, iTurb,height,tau,distribute_tau):
        res = PeResults(self.casename,iTurb,height,tau,dirname=self.dirname,distribute_tau=distribute_tau)
        res.read()
        # find height and tau index of the simulation 
        [x_grid,z_grid] = np.meshgrid(res.x, res.z,indexing='ij')

        # add height and tau to list if not already read 
        if height not in self.height:
            self.height.append(height)
        if tau not in self.tau:
            self.tau.append(tau)

        deltaL = {}
        deltaL['height'] = height
        deltaL['tau'] = tau
        deltaL['freq'] = res.frequencies
        deltaL['x'] = res.x
        deltaL['z'] = res.z
        deltaL['x_grid'] = x_grid
        deltaL['z_grid'] = z_grid
        deltaL['val'] = res.deltaL 
        self.deltaLlist.append(deltaL)

    # read receiver from file 
    # one file corespond to one plane of propagation with one source altitude
    def read_receiver(self,iTurb,height,tau,distribute_tau):
        res = PeResults(self.casename,iTurb,height,tau,dirname=self.dirname,distribute_tau=distribute_tau)
        res.read()
        # find height and tau index of the simulation 
        [x_grid,z_grid] = np.meshgrid(res.x, res.z,indexing='ij')

        # add height and tau to list if not already read 
        if height not in self.height:
            self.height.append(height)
        if tau not in self.tau:
            self.tau.append(tau)

        deltaL = {}
        deltaL['height'] = height
        deltaL['tau'] = tau
        deltaL['freq'] = res.frequencies
        deltaL['x'] = res.x
        deltaL['z'] = res.heights
        deltaL['x_grid'] = x_grid
        deltaL['z_grid'] = z_grid
        deltaL['val'] = res.receiver 
        self.deltaLlist.append(deltaL)

    # check compatibility in frequency and receiver heigths of all PE calculation
    def checkCompatibility(self):
        if not self.deltaLlist:
            print('you need to read some PE results first ...')
            return -1

        # check frequencies 
        frequencies = [deltaL['freq'] for deltaL in self.deltaLlist]
        if chkList(frequencies):
            print('frequencies OK')
        else :
            print('frequency are not the same')
            return -1 

        #check receievr Heights
        heights = [deltaL['z'] for deltaL in self.deltaLlist]
        if chkList(heights):
            print('receiver heights OK')
        else :
            print('receiver heights are not the same')

        #TODO add other check on grid and dx 

    # concatenate all pe results in on matrix (nx,nz,ntau)
    # if rec == true wont try to modify the size of z (could be removed) 
    def concatenate(self, rec = True):
        # find xmin, xma, zmin, zmax
        #---------------------------------------------------------------------------------------
        # initialize min max values 
        xmin = self.deltaLlist[0]['x'][0]
        xmax = self.deltaLlist[0]['x'][-1]
        zmin = self.deltaLlist[0]['z'][0]
        zmax = self.deltaLlist[0]['z'][-1]

        # browse all deltaL to find min max values 
        for ii in range(1,len(self.deltaLlist)):
            if self.deltaLlist[ii]['x'][0]<xmin :
                xmin = self.deltaLlist[ii]['x'][0]
            if self.deltaLlist[ii]['z'][0]<zmin :
                zmin = self.deltaLlist[ii]['z'][0]

            if self.deltaLlist[ii]['x'][-1]>xmax :
                xmax = self.deltaLlist[ii]['x'][-1]
            if self.deltaLlist[ii]['z'][-1]>zmax :
                zmax = self.deltaLlist[ii]['z'][-1]

        # rearange height and tau 
        #---------------------------------------------------------------------------------------
        self.height.sort()
        self.tau.sort()
        self.height = np.array(self.height)
        self.tau = np.array(self.tau)
        self.frequencies = self.deltaLlist[0]['freq']

        self.ntau = len(self.tau)
        self.nheight = len(self.height)

        # reshaping
        #---------------------------------------------------------------------------------------
        # assuming all dx are the same 
        dx = self.deltaLlist[-1]['x'][1] - self.deltaLlist[-1]['x'][0]
        # creating grid for x 
        self.x = np.arange(xmin,xmax+dx,dx)
        self.nx = len(self.x)
        # creating grid for z 
        if rec ==False:
            self.z = np.arange(zmin,zmax+dx,dx)
        else :
            self.z  = self.deltaLlist[-1]['z']
            self.nz = len(self.z)

        # create the complete matrix of deltaL 
        self.Nfreq = len(self.deltaLlist[-1]['freq'])
        self.deltaL = np.zeros((self.nx, self.nz, self.ntau, self.Nfreq, self.nheight))
        print(self.deltaL[:,:,0,:,0].shape)
        [self.x_grid,self.z_grid] = np.meshgrid(self.x, self.z,indexing='ij')

        for deltaL in self.deltaLlist:
            iheight = np.argmin(abs(self.height - deltaL['height']))
            itau    = np.argmin(abs(self.tau - deltaL['tau']))
            # find corresponding index 
            ixmin = np.argmin(abs(self.x - deltaL['x'][0]))
            izmin = np.argmin(abs(self.z - deltaL['z'][0]))
            ixmax = np.argmin(abs(self.x - deltaL['x'][-1]))
            izmax = np.argmin(abs(self.z - deltaL['z'][-1]))
            print('height = ' + str(deltaL['height']) + ', tau = '+str(deltaL['tau']))
            self.deltaL[ixmin:ixmax+1,izmin:izmax+1,itau,:,iheight] = deltaL['val'][:,:,:]

    # create the polar coordinate mesh for the concaten,ated field
    def createMesh(self):
        angles = np.reshape(self.tau*np.pi/180,(1,-1))
        r = np.reshape(self.x,(-1,1))

        self.x_polar = r*np.cos(angles)
        self.y_polar = r*np.sin(angles)

    # duplicate matrix size and devide grid step by alpha
    # TODO warning the refine method for tau is not the same than for x and z  
    def refine(self,alpha,axis):
        # x axis
        if axis ==0 :
            # temp value
            deltaL = np.zeros(((self.nx-1)*alpha +1,self.nz,self.ntau,self.Nfreq,self.nheight))
            x_grid = np.zeros(((self.nx-1)*alpha +1,self.nz,self.ntau))
            z_grid = np.zeros(((self.nx-1)*alpha +1,self.nz,self.ntau))
            # tau_grid = np.zeros(((self.nx-1)*alpha +1,self.nz,self.ntau))

            deltaL[::alpha,...] = self.deltaL
            x_grid[::alpha,...] = self.x_grid
            z_grid[::alpha,...] = self.z_grid
            # tau_grid[::alpha,...] = self.tau_grid
            for ii in range(alpha-1):
                deltaL[ii+1::alpha,...] = self.deltaL[:-1,...]
                x_grid[ii+1::alpha,...] = self.x_grid[:-1,...]
                z_grid[ii+1::alpha,...] = self.z_grid[:-1,...]
                # tau_grid[ii+1::alpha,...] = self.tau_grid[:-1,...]
            
            self.deltaL = deltaL
            self.x_grid = x_grid
            self.z_grid = z_grid
            # self.tau_grid = tau_grid 
            self.nx = (self.nx-1)*alpha +1
            self.x = np.linspace(self.x[0],self.x[-1],self.nx)
        # z axis
        if axis == 1:
            deltaL = np.zeros((self.nx,(self.nz-1)*alpha+1,self.ntau,self.Nfreq,self.nheight))
            x_grid = np.zeros((self.nx,(self.nz-1)*alpha +1,self.ntau))
            z_grid = np.zeros((self.nx,(self.nz-1)*alpha +1,self.ntau))
            # tau_grid = np.zeros((self.nx,(self.nz-1)*alpha +1,self.ntau))

            deltaL[:,::alpha,...] = self.deltaL
            x_grid[:,::alpha,...] = self.x_grid
            z_grid[:,::alpha,...] = self.z_grid
            # tau_grid[:,::alpha,...] = self.tau_grid
            for ii in range(alpha-1):
                deltaL[:,ii+1::alpha,...] = self.deltaL[:,:-1,...]
                x_grid[:,ii+1::alpha,...] = self.x_grid[:,:-1,...]
                z_grid[:,ii+1::alpha,...] = self.z_grid[:,:-1,...]
                # tau_grid[:,ii+1::alpha,...] = self.tau_grid[:,:-1,...]

            self.nz = (self.nz-1)*alpha +1
            self.z = np.linspace(self.z[0],self.z[-1],self.nz)
            self.deltaL = deltaL
            self.x_grid = x_grid
            self.z_grid = z_grid
            # self.tau_grid = tau_grid 

        # tau axis
        # TODO finish modify grid for tau axis 
        if axis == 2:
            deltaL= np.zeros((self.nx,self.nz,2*self.ntau-1,self.Nfreq,self.nheight))
            tau = np.zeros((2*self.ntau-1))

            for ii in range(self.ntau-1):
                deltaL[:,:,2*ii,:] = self.deltaL[:,:,ii,:,:]
                deltaL[:,:,2*ii +1,:] = 0.5*(self.deltaL[:,:,ii,:,:] + self.deltaL[:,:,ii+1,:,:])
                tau[2*ii] = self.tau[ii]
                tau[2*ii+1] = 0.5*(self.tau[ii] + self.tau[ii+1])

            deltaL[:,:,-1,:] = self.deltaL[:,:,-1,:]
            tau[-1] = self.tau[-1]

            self.deltaL = deltaL
            self.tau = tau
            self.ntau = len(self.tau)

    # create a results for tau = 360 from result tau=0 
    def loopAngle(self):
        if self.tau[0] != 0:
            print('first angle is not 0')
            return -1

        deltaL = np.zeros((self.nx, self.nz, self.ntau+1, self.Nfreq, self.nheight))
        tau = np.zeros((self.ntau+1,))
        deltaL[:,:,0:-1,:,:] = self.deltaL
        deltaL[:,:,-1,:,:] =self.deltaL[:,:,0,:,:]

        tau[0:-1] = self.tau
        tau[-1] = 360
        self.tau = tau 
        self.deltaL = deltaL 
        self.ntau = len(self.tau)

    # linear interpolate Pe results on a cartesian grid x,y 
    def interpolate(self,x,y):
        # create vector of coordinates of the original grid 
        xy_polar = np.concatenate((self.x_polar.reshape((-1,1)),self.y_polar.reshape((-1,1))),1)

        print('start interpolation ...')
        # create vector of coordinates of the new grid 
        self.x_cart,self.y_cart = np.meshgrid(x,y)
        xy_cart=np.zeros([self.x_cart.shape[0]*self.x_cart .shape[1],2])
        xy_cart[:,0]=self.x_cart .flatten()
        xy_cart[:,1]=self.y_cart.flatten()
        # create knots for the interpolation 
        vtx, wts = interp_weights(xy_polar, xy_cart)
        print('finished creating knot ...')

        self.deltaLInterpolated = np.zeros((self.x_cart.shape[0],self.x_cart.shape[1],self.nz,self.Nfreq,self.nheight))
        # loop of z, freq, heights and interpolate using the previously computes knots 
        print('starting loop on height and frequency band ...')
        for iz in range(self.nz):
            for ifreq in range(self.Nfreq):
                for iheight in range(self.nheight):
                    values = self.deltaL[:,iz,:,ifreq,iheight].flatten()
                    valuesInterpolated = interpolate(values, vtx, wts)
                    self.deltaLInterpolated[:,:,iz,ifreq,iheight] = valuesInterpolated.reshape(self.x_cart.shape[0],self.x_cart .shape[1])
        print('done')
        # save the interpolation grid 
        [self.x_interpolate,self.y_interpolate,self.z_interpolate] = np.meshgrid(x,y,self.z,indexing='ij')
        

    def shiftDomain(self,xshift,yshift):
        self.x_polar = self.x_polar + xshift
        self.y_polar = self.y_polar + yshift

        self.x_cart = self.x_cart + xshift
        self.y_cart = self.y_cart + yshift

    def reduce(self,alpha,axis) :
        # x axis
        if axis == 0 :
            self.deltaL = self.deltaL[::alpha,...]
            self.x = self.x[::alpha]
            self.x_grid = self.x_grid[::alpha,...]
            self.z_grid = self.z_grid[::alpha,...]
            # self.tau_grid = self.tau_grid [::alpha,...]
            self.nx = len(self.x)
        if axis == 1 :
            self.deltaL = self.deltaL[:,::alpha,...]
            self.z = self.z[::alpha]
            self.x_grid = self.x_grid[:,::alpha,...]
            self.z_grid = self.z_grid[:,::alpha,...]
            # self.tau_grid = self.tau_grid [:,::alpha,...]
            self.nz = len(self.z)

        # TODO finish this for tau 
        if axis == 2 :
            self.deltaL = self.deltaL[:,:,::alpha,...]
            self.tau = self.tau[::alpha]
            self.ntau = len(self.tau)

    def truncate(self,max,axis, min=0):
        if axis == 0:
            Nxmax = np.argmin(abs(self.x-max))
            Nmin = np.argmin(abs(self.x-min))
            self.deltaL = self.deltaL[Nmin:Nxmax+1,...]
            self.x = self.x[Nmin:Nxmax+1]
            self.x_grid = self.x_grid[Nmin:Nxmax+1,...]
            self.z_grid = self.z_grid[Nmin:Nxmax+1,...]
        if axis == 1:
            Nzmax = np.argmin(abs(self.z-max))
            Nmin = np.argmin(abs(self.z-min))
            self.deltaL = self.deltaL[:,Nmin:Nzmax+1,...]
            self.z = self.z[Nmin:Nzmax+1]
            self.x_grid = self.x_grid[:,Nmin:Nzmax+1,...]
            self.z_grid = self.z_grid[:,Nmin:Nzmax+1,...]

    # axis is the ax orthogonal to the plane
    def takeOnePlane(self,pos1,axis):
        if axis == 0:
            ipos1 = np.argmin(abs(self.x - pos1))
            self.deltaL = self.deltaL[ipos1,...]
        if axis == 1:
            ipos1 = np.argmin(abs(self.z - pos1))
            self.deltaL = self.deltaL[:,ipos1,...]
        if axis == 2:
            ipos1 = np.argmin(abs(self.tau - pos1))
            self.deltaL = self.deltaL[:,:,ipos1,...]

    # here axis is the ax of the line
    def takeOneLine(self,pos1,pos2,axis):
        if axis == 0:
            ipos1 = np.argmin(abs(self.z - pos1))
            ipos2 = np.argmin(abs(self.tau - pos2))
            self.deltaL = self.deltaL[:,ipos1,ipos2,...]
            self.x_grid = self.x_grid[:,ipos1,...]
            self.z_grid = self.z_grid[:,ipos1,...]
        if axis == 1:
            ipos1 = np.argmin(abs(self.x - pos1))
            ipos2 = np.argmin(abs(self.tau - pos2))
            self.deltaL = self.deltaL[ipos1,:,ipos2,...]
            self.x_grid = self.x_grid[ipos1,:,ipos2,...]
            self.z_grid = self.z_grid[ipos1,:,ipos2,...]
        if axis == 2:
            ipos1 = np.argmin(abs(self.x - pos1))
            ipos2 = np.argmin(abs(self.z - pos2))
            self.deltaL = self.deltaL[ipos1,ipos2,...]
            self.x_grid = self.x_grid[ipos1,:,ipos2,...]
            self.z_grid = self.z_grid[ipos1,:,ipos2,...]

    def reshape(self,x_array,z_array=None):
        print('delta L shape : ' + str(self.deltaL.shape))
        print('reshaping Delta L ...')
        # Change delta X to fit
        if (x_array[1]-x_array[0])>(self.x[1]-self.x[0]):
            alpha = (x_array[1]-x_array[0])/(self.x[1]-self.x[0])
            if int(alpha)!=alpha:
                print('mesh incompatible')
                quit()
            self.reduce(int(alpha),0)
        if (x_array[1]-x_array[0])<(self.x[1]-self.x[0]):
            alpha = (self.x[1]-self.x[0])/(x_array[1]-x_array[0])
            if int(alpha)!=alpha:
                print('mesh incompatible')
                quit()
            self.refine(int(alpha),0)
        if z_array is not None :
            # Change delta z to fit
            if (z_array[1]-z_array[0])>(self.z[1]-self.z[0]):
                alpha = (z_array[1]-z_array[0])/(self.z[1]-self.z[0])
                if int(alpha)!=alpha:
                    print('mesh incompatible')
                    quit()
                self.reduce(int(alpha),1)
            if (z_array[1]-z_array[0])<(self.z[1]-self.z[0]):
                alpha = (self.z[1]-self.z[0])/(z_array[1]-z_array[0])
                if int(alpha)!=alpha:
                    print('mesh incompatible')
                    quit()
                self.refine(int(alpha),1)

        # reside domain
        if (x_array[-1]<=self.x[-1]):
            self.truncate(x_array[-1],0)
        else:
            print('delta L to small')
        if z_array is not None :
            if (z_array[-1]<=self.z[-1]):
                self.truncate(z_array[-1],1)
            else:
                print('delta L to small ')
        print('delta L shape : ' + str(self.deltaL.shape))

    def plotTopRaw(self,z,freq,height,cmap='RdBu_r'):
        ifreq = np.nonzero(self.frequencies==freq)[0][0]
        iz = np.nonzero(self.z==z)[0][0]
        iheight = np.nonzero(self.height==height)[0][0]

        plt.pcolormesh(self.x_polar,self.y_polar,self.deltaL[:,iz,:,ifreq,iheight],cmap=cmap,shading='auto')
        plt.gca().set_aspect('equal', adjustable='box')

    def plotTopInterpolated(self,z,freq,height,cmap='RdBu_r'):
        ifreq = np.nonzero(self.frequencies==freq)[0][0]
        iz = np.nonzero(self.z==z)[0][0]
        iheight = np.nonzero(self.height==height)[0][0]

        plt.pcolormesh(self.x_cart,self.y_cart,self.deltaLInterpolated[:,:,iz,ifreq,iheight],cmap=cmap,shading='auto')
        plt.gca().set_aspect('equal', adjustable='box')
