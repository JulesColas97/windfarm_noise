import os
import shutil
import pickle
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
from .les import Les
from .source.mesh import Mesh
from .source.main import Source
from .source.utils import interp_3D_atmos_data
from .source.wind_turbine import WindTurbine
from .utils import computeThirdOctaveFrequencies


class Simu():
    ##
    # @brief Initializes the Simu class with a dictionary of inputs for a PE simulation.
    #
    # @param casename (str): Name of the case. 
    # initialize the input dictionnary for the PE simulation 
    def __init__(self,casename: str):
        """initialize Instance Simu, with a dictionnary of input for a PE simulation

        Args:
            casename (str): name of the case
        """
        self.casename=casename
        self.xplanes = None
        self.yplanes = None
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
            "Ly1"                     : 1000.,
            "Ly2"                     : 1000.,
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
            "tinput%T_scale_K"        : 265,
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

    ##
    # @brief Reads and adds a LES flow simulation and stores it as a Les object.
    #
    # @param path (str): Path to the LES data.
    def readLes(self,path: str):
        """read and add a LES flow simulation and store it as a :class:`les.LES`

        Args:
            path (_type_): _description_
        """
        self.les = Les(path)
        self.les.read()
        print('les domain dimension:')
        print(self.les.L_x,self.les.L_y,self.les.L_z)
        print(str(self.les.turbines.shape[0])+' wind turbines :')
        print(self.les.turbines)

    ##
    # @brief Creates an input.nml file from the inputs dictionary.
    #
    # @param directory (str, optional): Directory to save the input file. Defaults to './'.
    def createInputFile(self,directory: str ='./'):
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

    ##
    # @brief Sets an input in the inputs dictionary.
    #
    # @param key (str): Key of the input.
    # @param value (Any): Value of the input.
    def set_input(self,key: str, value: any):
        self.inputs[key] = value

    ##
    # @brief Sets the frequencies for the simulation. If only fc is given only these frequencies are computed
    #
    # @param fc (list): List of center frequencies.
    # @param Nfc (list, optional): List of the number of frequencies per band. Defaults to None.
    def set_frequencies(self,fc: list, Nfc: list = None):
        if Nfc is None:
            self.inputs['nb_freq'] = len(fc)
            self.inputs['frequencies(1)'] = list(fc)
            self.frequencies = np.array(fc)
        else : 
            self.inputs['nb_freq'] = sum(Nfc)
            self.inputs['frequencies(1)'] = list(computeThirdOctaveFrequencies(fc,Nfc))
            self.frequencies = np.array(computeThirdOctaveFrequencies(fc,Nfc))

    # define the case to be simulate with number of source heights, angles and flow data 
    ##
    # @brief Defines the cases to be simulated with the number of source heights, angles, and flow data.
    #
    # @param heights (list): List of source heights.
    # @param angles (list): List of propagation angles.
    # @param src_path (str, optional): Path to the source code. Defaults to '/home/lmfa/jcolas/Documents/DEV/wf_phd/src/kernel/New_PE/PE_2D_WAPE'.
    # @param flow_path (str, optional): Path to the flow data. Defaults to '/home/lmfa/jcolas/Documents/DEV/LES/'.
    # @param ratio (float, optional): Ratio for the LES data. Defaults to None.
    def defineCases(self, heights: list, angles: list, 
                    src_path: str = '/home/lmfa/jcolas/Documents/DEV/wf_phd/src/kernel/New_PE/PE_2D_WAPE',
                    flow_path: str = '/home/lmfa/jcolas/Documents/DEV/LES/', ratio: float = None):
        self.heights=heights
        self.tau = angles
        self.src_path = src_path

        self.inputs['fdir'] = flow_path+'output/'
        print(self.inputs['fdir'])
        self.les = Les(flow_path)

        self.les.read(ratio=ratio)
        self.tx = self.les.turbines[:,1]
        self.ty = self.les.turbines[:,2]
        if ratio is not None:
            self.inputs['tinput%ratio'] = ratio
        else :
            self.inputs['tinput%ratio'] = self.les.ug
            self.inputs['tinput%ratio'] = 1
        self.inputs['tinput%z_i'] = self.les.z_i
        self.inputs['tinput%T_scale_K'] = self.les.T_scale
        self.inputs['tinput%Lx'] = self.les.lx 
        self.inputs['tinput%Ly'] = self.les.ly 
        self.inputs['tinput%Lz'] = self.les.lz 

    ##
    # @brief Defines the total domain for the simulation.
    #
    # @param x1 (float): Start of the x-domain.
    # @param x2 (float): End of the x-domain.
    # @param y1 (float): Start of the y-domain.
    # @param y2 (float): End of the y-domain.
    # @param h (float): Height of the domain.
    def defineDomain(self,x1: float,x2: float,y1: float,y2: float,h: float):
        self.x1 = x1
        self.x2 = x2
        self.y1 = y1
        self.y2 = y2
        self.h = h 


    ##
    # @brief Sets the x constant planes for the simulation.
    #
    # @param xplanes (list): List of x-planes.
    def set_xplanes(self,xplanes: list):
        # x constant planes
        self.xplanes = np.array(xplanes)

    ##
    # @brief Sets the y constant planes for the simulation.
    #
    # @param yplanes (list): List of y-planes
    def set_yplanes(self,yplanes: list):
        # y constant planes
        self.yplanes = np.array(yplanes)

    ##
    # @brief Sets the specific domain for each turbine.
    #
    # @param iturbine (int): Index of the turbine.
    def set_domain(self,iturbine: int):
        Lx1 = self.x2-self.tx[iturbine]*self.les.z_i
        Lx2 = self.tx[iturbine]*self.les.z_i - self.x1 
        Ly1 = self.y2 - self.ty[iturbine]*self.les.z_i
        Ly2 = self.ty[iturbine]*self.les.z_i - self.y1
        Ly = max(Ly1,Ly2)
        Lz = self.h

        self.set_input('Lx1',Lx1)
        self.set_input('Lx2',Lx2)
        self.set_input('Ly1',Ly1)
        self.set_input('Ly2',Ly2)
        self.set_input('Lz',Lz)

        if self.xplanes is not None:
            self.set_input('nb_xplane',len(self.xplanes))
            self.set_input('xplane(1)',list(self.xplanes-self.tx[iturbine]*self.les.z_i))
        else:
            self.set_input('nb_xplane',0)

        if self.yplanes is not None:
            self.set_input('nb_yplane',len(self.yplanes))
            self.set_input('yplane(1)',list(self.yplanes-self.ty[iturbine]*self.les.z_i))
        else:
            self.set_input('nb_yplane',0)

    ##
    # @brief Computes an estimation of the simulation time.
    #
    # @return (float): Total simulation time.
    def computeSimulationTime(self):
        TIME_1km_45freq_s = 195

        Lx1 = self.inputs['Lx1']
        Lx2 = self.inputs['Lx2']
        Ly1 = self.inputs['Ly1']
        Ly2 = self.inputs['Ly2']
        print(Lx1)
        print(Lx2)

        area_rectangle = (Lx1 + Lx2) * (Ly1 + Ly2)
        area_circle = np.pi*1000**2
        t_tot = TIME_1km_45freq_s * self.tau.size * area_rectangle / area_circle
        return t_tot


    # create launch file for slurm 
    ##
    # @brief Creates a launch file for SLURM.
    #
    # @param dirname (str): Directory name.
    # @param jname (str): Job name.
    # @param mem (int, optional): Memory allocation. Defaults to 5000.
    # @param time (str, optional): Time limit. Defaults to "2:00:00".
    def createLaunch(self,dirname: str,jname: str,mem: int = 5000,time: str = "2:00:00") :
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

    ##
    # @brief Creates a parallel launch file for SLURM with distribution over propagation angles
    #
    # @param dirname (str): Directory name.
    # @param jname (str): Job name.
    # @param distribute_tau (int): Number of tau distributions.
    # @param mem (int, optional): Memory allocation. Defaults to 5000.
    # @param time (str, optional): Time limit. Defaults to "2:00:00".
    def createParallelLaunch(self,dirname: str,jname: str,distribute_tau: int,mem: int = 5000,time: str = "2:00:00") :
        f= open(dirname+'/launch.sh','w')
        f.write("#!/bin/bash\n")
        f.write("#SBATCH --job-name="+jname+"\n")
        f.write("#SBATCH --output=out.out # output messages go here\n")
        f.write("#SBATCH --error=err.err    # error messages go here\n")
        f.write("#SBATCH --mail-user=jules.colas@ecl17.ec-lyon.fr\n")
        f.write("#SBATCH --mail-type=ALL\n")
        f.write("#SBATCH --partition=haswell # partition name\n")
        f.write("#SBATCH --nodes=1\n")
        f.write("#SBATCH --cpus-per-task="+str(distribute_tau)+"\n")
        f.write("#SBATCH --mem=" + str(mem*distribute_tau) + "\n")
        f.write("#SBATCH --time=" + str(time) + "\n")


        f.write("module purge\n")
        f.write("module load  HDF5/1.10.1-intel-2018a\n")
        f.write("export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK\n")
        f.write("export KMP_AFFINITY=granularity=fine,compact,1,0\n")
        f.write("export OMP_STACKSIZE=1g\n")
        f.write("ulimit -s unlimited\n")
        for kk in range(distribute_tau):
            f.write("cd ./tau" +str(kk)+"\n")
            f.write("time ./PE_2D &\n")
            f.write("cd .. " +"\n")
        f.write("wait")


    ##
    # @brief Distributes cases for the simulation. create the directorise, input.nml files, launch.sh files for running the complete simulation.I does not launch the jobs
    #
    # @param distribute_tau (int, optional): Number of tau distributions. Defaults to None.
    # @param mem (int, optional): Memory allocation. Defaults to 5000.
    # @param time (str, optional): Time limit. Defaults to "2:00:00".
    # @param turbine_index (list, optional): List of turbine indices. Defaults to None.
    def distributeCases(self,distribute_tau: int = None,mem: int = 5000,
                        time: str = "2:00:00",turbine_index: list = None):
        dir = os.getcwd()
        self.distribute_tau = distribute_tau
        if distribute_tau==False: 
            self.inputs['nb_theta'] = len(self.tau)
            self.inputs['theta(1)'] = list(self.tau)

        if turbine_index is None:
            turbine_index = np.arange(0,len(self.tx))
        
        print(turbine_index)
        # loop over the wind turbine 
        for ii in turbine_index:
            print('t%s :'%(ii), self.tx[ii], self.ty[ii])
            # set wind turbines location in input.nml
            # the PE simulation wil always 
            # be perfomed with the source at (0,0,zS)
            self.inputs["tinput%posx"] = self.tx[ii]
            self.inputs["tinput%posy"] = self.ty[ii]
            # set the PE domain according to source position and total domain
            self.set_domain(ii)
            # loop over the source heights
            for jj in range(len(self.heights)):
                # set source heights in input.nml
                self.inputs['src%pos_z'] = self.heights[jj]
                # loop over the propagation angles
                if distribute_tau is not None:
                    # if tau distributed launch a job for each angle
                    ntau = self.tau.size//distribute_tau
                    rtau = self.tau.size % distribute_tau
                    for kk in range(distribute_tau-1):
                        self.inputs['nb_theta'] = ntau
                        self.inputs['theta(1)'] = list(
                                self.tau[kk*ntau:(kk+1)*ntau])
                        path = os.path.join(dir,
                                            't' + str(ii) + '/' +
                                            str(self.heights[jj]) +
                                            '/tau' + str(kk) + '/')
                        os.makedirs(path, exist_ok=True)
                        self.createInputFile(directory=path)
                        shutil.copy(self.src_path, path + '/PE_2D')

                    self.inputs['nb_theta'] = ntau + rtau
                    self.inputs['theta(1)'] = list(self.tau[(distribute_tau - 1) * ntau:])
                    path = os.path.join(dir, 't' + str(ii) + '/' + str(self.heights[jj])
                                        + '/tau' + str(distribute_tau-1) + '/')
                    os.makedirs(path,exist_ok=True)
                    self.createInputFile(directory=path)
                    shutil.copy(self.src_path,path+'/PE_2D')


                    self.createParallelLaunch(os.path.join(dir, 't'+str(ii)+'/'+str(self.heights[jj])+'/'),
                                      self.casename+'t'+str(ii)+'h'+str(self.heights[jj]),
                                      distribute_tau,
                                      mem=mem,time=time)
                # if tau is not distributed one job is launched for all propagation 
                else:
                    self.inputs['nb_theta'] = self.tau.size
                    self.inputs['theta(1)'] = list(self.tau)
                    path = os.path.join(dir, 't'+str(ii)+'/'+str(self.heights[jj]))
                    os.makedirs(path,exist_ok=True)
                    self.createInputFile(directory=path)
                    self.createLaunch(path,self.casename+'t'+str(ii)+'h'+str(self.heights[jj]),mem=mem,time=time)
                    shutil.copy(self.src_path,path+'/PE_2D')

    ##
    # @brief Creates launch files for the simulation. This was coded specifically for haswell partition on Newton HPC.
    # This done to launch several PE cases with the same job in order to avoid launchind hundreds of jobs in the case of a wind farm
    #
    # @param mem (int): Memory allocation.
    # @param time (str): Time limit.
    # @param turbine_index (list, optional): List of turbine indices. Defaults to None.
    def createLaunchFiles(self, mem: int, time: str, turbine_index: list = None):
        dirname = os.getcwd()
        haswell_mem_size = 64000
        haswell_nb_core = 16
        count = 0
        countjobs = 0
        nlaunchPerJob = min(haswell_nb_core, int(haswell_mem_size / mem))

        if turbine_index is None:
            turbine_index = np.arange(0, len(self.tx))

        self.nJob = int((len(turbine_index)*len(self.heights))
                        // nlaunchPerJob)
        self.rest_cases = int((len(turbine_index)*len(self.heights))
                        % nlaunchPerJob)

        for ii in turbine_index:
            for jj in range(len(self.heights)):
                if count % nlaunchPerJob == 0:
                    print(countjobs)

                    f = open(dirname + '/launch' + str(count//nlaunchPerJob) + '.sh','w')
                    f.write("#!/bin/bash\n")
                    f.write("#SBATCH --job-name=PE"+str(count//nlaunchPerJob)+"\n")
                    f.write("#SBATCH --output=out.out # output messages go here\n")
                    f.write("#SBATCH --error=err.err    # error messages go here\n")
                    f.write("#SBATCH --mail-user=jules.colas@ecl17.ec-lyon.fr\n")
                    f.write("#SBATCH --mail-type=ALL\n")
                    f.write("#SBATCH --partition=cascade # partition name\n")
                    f.write("#SBATCH --nodes=1\n")
                    if countjobs == (self.nJob):
                        f.write("#SBATCH --ntasks="+str(self.rest_cases)+"\n")
                    else:
                        f.write("#SBATCH --ntasks="+str(nlaunchPerJob)+"\n")
                    f.write("#SBATCH --cpus-per-task=1\n")
                    #f.write("#SBATCH --mem=" + str(mem) + "\n")
                    f.write("#SBATCH --mem=64000\n")
                    f.write("#SBATCH --time=" + str(time) + "\n")

                    f.write("module purge\n")
                    f.write("module load  HDF5/1.10.1-intel-2018a\n")
                    f.write("export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK\n")
                    f.write("export KMP_AFFINITY=granularity=fine,compact,1,0\n")
                    f.write("export OMP_STACKSIZE=1g\n")
                    f.write("ulimit -s unlimited\n")

                f.write("cd t%s/%s/ \n"%(ii,self.heights[jj]))
                f.write("time ./PE_2D >out.out & \n")
                f.write("cd ../../ \n")
                count += 1

                if (count) % nlaunchPerJob == 0:
                    f.write("wait")
                    f.close()
                    countjobs += 1
        f.write("wait")
        f.close()

    def createLocalLaunchFiles(self, mem: int, time: str, turbine_index: list = None):
        dirname = os.getcwd()
        haswell_mem_size = 64000
        haswell_nb_core = 16
        count = 0
        countjobs = 0
        nlaunchPerJob = min(haswell_nb_core, int(haswell_mem_size / mem))

        if turbine_index is None:
            turbine_index = np.arange(0, len(self.tx))
        
        f = open(dirname + '/launch.sh','w')
        for ii in turbine_index:
            for jj in range(len(self.heights)):
                print(jj)

                f.write("#!/bin/bash\n")
                f.write("cd t%s/%s/ \n"%(ii,self.heights[jj]))
                f.write("time ./PE_2D >out.out\n")
                f.write("cd ../../ \n")

        f.close()

    ##
    # @brief Launches the simulation cases. When one job corresponds to one Height and one turbine 
    #
    # @param turbine_index (list, optional): List of turbine indices. Defaults to None.
    def launchCases(self,turbine_index: list = None):
        dir = os.getcwd()
        if turbine_index is None:
            turbine_index = np.arange(0, len(self.tx))
        for ii in turbine_index:
            print('t%s'%(ii))
            for jj in range(len(self.heights)):
                path = os.path.join(dir, 't'+str(ii)+'/'+str(self.heights[jj]))  
                os.system("cd "+path+"; sbatch launch.sh")

    ##
    # @brief Launches the simulation cases. When the function  createLaunchFiles() has been used. 
    def launchCases2(self):
        dir = os.getcwd()
        for ii in range(self.nJob):
            os.system("sbatch launch%s.sh" % (ii))

    ##
    # @brief Constructs the name of the h5 file storing the PE results.
    # the data stucture is such that if tau is not distibuted: <self.dirname>/t<iTurb>/<height>/<self.casename>_<tau>.h5
    #if tau is distributed:
    #<self.dirname>/t<iTurb>/<height>/tau<itau>/<self.casename>_<tau>.h5
    #
    # @param iTurb (int): Index of the turbine.
    # @param tau (float): Propagation angle.
    # @param height (float): Source height.
    # @return (str): Complete path to the h5 file.
    def fname(self, iTurb: int, 
              tau: float, height: float) -> str:
        # convert angle to str format
        if int(tau) == tau:
            tau_str = format(tau, "04d")
        else:
            tau_str = str(tau)

        dirname = os.getcwd()
        # if tau distributed find the right folder
        if self.distribute_tau is not None:
            for kk in range(self.distribute_tau):
                path = (
                    dirname
                    + "/t"
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
            path = dirname + "/t" + str(iTurb) + "/" + str(height) + "/"
        return path + self.casename + "_" + tau_str + ".h5"

    ##
    # @brief Checks if the simulation cases have run. This only check if the file has been created.
    def check_run_cases(self):
        for ii in range(len(self.tx)):
            for jj in range(len(self.heights)):
                for tt in range(len(self.tau)):
                    fname = self.fname(ii, self.tau[tt], self.heights[jj])
                    if not os.path.isfile(fname):
                        print(fname)
    ##
    # @brief Makes the source for the simulation.
    #
    # @param wt (WindTurbine): Wind turbine object.
    # @param mesh (Mesh): Mesh object.
    # @param offset (float, optional): Offset for the wind speed profile. Defaults to 100.
    # @param U_inf (list, optional): Wind speed profile. if None self.les is used. Defaults to None.
    # @param z_coord (list, optional): Vertical coordinates. Defaults to None.
    # @param fname (str, optional): File name to save the source. Defaults to None.
    # @param iTurb (list, optional): List of turbine indices. Defaults to None.
    # @param Ncore (int, optional): Number of cores. Defaults to 16.
    # @param plot (bool, optional): Whether to plot the wind speed profile. Defaults to False.
    # @param ratio (float, optional): Ratio for the LES data. Defaults to 1.
    # @param epsilon (float, optional): To set a constant epsilon instead of the LES data. Defaults to None.
    # @param omega (float, optional): To set a constant Omega for the wind turbine instead of using NREL table. Defaults to None.
    def makeSource(self, wt: WindTurbine, mesh: Mesh, offset: float = 100,
                   U_inf: np.ndarray = None, z_coord: np.ndarray = None,
                   fname: str = None, iTurb: list = None,
                   Ncore: int=16, plot: bool = False, ratio: float=1,
                   epsilon: float=None, omega: float=None):
        # check
        try: self.tx
        except:
            print('turbine not define')
            return

        if mesh.polar:
            print('mesh polar is not handled')
            return
        if iTurb is None:
            iTurb = np.arange(0, len(self.tx))

        for ii in iTurb:
            print('start computation for Turbine ', ii)
            mesh1 = deepcopy(mesh)
            mesh1.shift(-self.tx[ii]*self.les.z_i,
                        -self.ty[ii]*self.les.z_i)
            mesh1.cartesian_2_polar()

            # set wind speed profile according to les data or input
            if (U_inf is not None) and (z_coord is not None):
                self.les.U_inf = U_inf
                self.les.z_coord = z_coord
            else:
                self.les.dissipation()
                self.les.takeProfile(self.tx[ii]*self.les.z_i-offset,
                                     self.ty[ii]*self.les.z_i,
                                     ratio=ratio, epsilon=epsilon)

            xnew = np.array([[self.tx[ii]*self.les.z_i-offset]])
            ynew = np.array([[self.ty[ii]*self.les.z_i]])
            znew = np.array([[wt.href]])
            U_hub,epsilon_hub = interp_3D_atmos_data(self.les,xnew,ynew,znew)


            #U_hub = self.les.U_inf[np.argmin(np.abs(self.les.z_coord
            #                                        - wt.href))]
            #epsilon_hub = self.les.epsilon_Kol[np.argmin(np.abs(self.les.z_coord
            #                                                    - wt.href))]
            wt.controlRotSpeed(U_hub, omega=omega)
            wt.setOptimalTwist(U_hub, 4)
            wt.absolute_pos = (self.tx[ii]*self.les.z_i,
                               self.ty[ii]*self.les.z_i)
            self.U_hub = U_hub
            print('U_hub', U_hub)
            print('epsilon_hub', epsilon_hub)
            print('Omega', wt.omega*60/(2*np.pi))
            # self.les.constant_epsilon(0.01)

            if plot:
                plt.plot(self.les.U_inf, self.les.z_coord)
                plt.figure()
                plt.plot(self.les.epsilon_Kol, self.les.z_coord)
                plt.show()

            src = Source(wt, self.les, mesh1)
            src.computeSpp(self.frequencies, Ncore)
            src.mesh.polar_2_cartesian()
            src.mesh.shift(self.tx[ii]*self.les.z_i,
                           self.ty[ii]*self.les.z_i)

            if fname is not None:
                src.save(fname+str(ii)+'.dat')
            else:
                src.save('Spp'+str(ii)+'.dat')
            src = None
            mesh1 = None
    
    ##
    # @brief Saves the class as a .dat file using pickle format.
    #
    # @param fname (str, optional): File name to save the class. Defaults to None.
    def save(self, fname: str = None):
        if fname is None:
            fname = self.casename
        with open(fname+'.dat', 'wb') as file:
            pickle.dump(self.__dict__, file)

    ##
    # @brief Loads the class from a .dat file using pickle format.
    #
    # @param fname (str, optional): File name to save the class. Defaults to None.
    def load(self, fname: str):
        print('loading simulation ...')
        with open(fname, 'rb') as file:
            self.__dict__ = pickle.load(file)
        print('done.')
