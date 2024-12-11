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
simu.set_input('side', False)
#set plane to save solution 
# same referentiel as LES sim 
simu.set_xplanes([1000])
simu.set_yplanes([1720])

# define domain 
# -----------------------------------------------------------------------------
xmin = 0
xmax = 5500
ymin = 0
ymax = 3400
zmax = 300
simu.defineDomain(xmin, xmax, ymin, ymax, zmax)

# discretize the rotor 
R = 60
hub = 90
Nh = 7
dh = 2*R/Nh
height = np.arange(hub-R+dh/2, hub+R, dh)
heights = np.rint(height).astype(int)
tau = np.array([0, 90, 180, 270])


# define cases with height, angles and flow 
path = '../../les_data/2TSBL/blue/'
src_path = '../../../src/kernel/New_PE_c/PE_2D_WAPE'
simu.defineCases(heights, tau, flow_path=path,
                 src_path=src_path, ratio=1.0)

# distribute tau --> parellized th simulation over tau 
# (divide the number of angle done by each process by distribute_tau value
# mem --> set the ram needed for one process,
# for 1kHz simulation up to 5km set to 6Gb
simu.distributeCases(distribute_tau=None, mem=6000,
                     time="01:00:00", turbine_index=None)
simu.createLocalLaunchFiles(mem=6000, time="01:00:00", turbine_index=None)
# simu.launchCases2()
# simu.launchCases(turbine_index=None)
simu.save()
quit()
