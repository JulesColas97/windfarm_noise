from prepost import Les
from prepost import Simu
from prepost import source as src
import numpy as np
import matplotlib.pyplot as plt

# create simu object
# -----------------------------------------------------------------------------
simu = Simu('c0')

# define frequency content
# -----------------------------------------------------------------------------
fc = [50, 100, 200]
simu.set_frequencies(fc)

#  define output parameters 
# -----------------------------------------------------------------------------
simu.set_input('side', True)
#set plane to save solution 
# same referentiel as LES sim 
simu.set_xplanes([1000])
simu.set_yplanes([1720])

# define domain 
# -----------------------------------------------------------------------------
xmin = 0
xmax = 1000 + 2000
ymin = 1720 - 1000
ymax = 1720 + 1000
zmax = 300
simu.defineDomain(xmin, xmax, ymin, ymax, zmax)

# discretize the rotor 
R = 60
hub = 90
Nh = 3
dh = 2*R/Nh
height = np.arange(hub-R+dh/2, hub+R, dh)
heights = np.rint(height).astype(int)
# tau = np.array([0, 90, 180, 270])
tau = np.array([0, 20,40,69,80, 90,100,120,140,160, 180,200,220,240,260,270])



# define cases with height, angles and flow 
# les_path = '../../les_data/2TSBL/blue/'
les_path = '/store/lmfa-2/acoustique-nl/simu_jules/LES/2T/S2/blue/'
# src_path = '../../../src/kernel/New_PE_c/PE_2D_WAPE'
src_path = '/home/lmfa/jcolas/Documents/DEV/windfarm_noise/src/kernel/New_PE_c/PE_2D_WAPE'
simu.defineCases(heights, tau, flow_path=les_path,
                 src_path=src_path, ratio=1.0)

# distribute tau --> parellized th simulation over tau 
# (divide the number of angle done by each process by distribute_tau value
# mem --> set the ram needed for one process,
# for 1kHz simulation up to 5km set to 6Gb
simu.distributeCases(distribute_tau=None, mem=1000,
                     time="01:00:00", turbine_index=None)

# to create launch files that can directly run on computer 
# be aware that for very large cases this is not adivsable and could take a very long time 
simu.createLocalLaunchFiles(mem=6000, time="01:00:00", turbine_index=None)
# to create launch case that can  be run on a super computer 
# Here the launch files are using SLURM managmenet system and were set up for the 
# Newton HPC
simu.createLaunchFiles(mem=6000, time="01:00:00", turbine_index=None)
# simu.launchCases2()
# simu.launchCases(turbine_index=None)
simu.save()
