from prepost import Les
from prepost import Simu
from prepost import source as src
import sys
import numpy as np
import matplotlib.pyplot as plt

# Input for the simulation
# -----------------------------------------------------------------------------
simu = Simu('c0')
xmin = 500
xmax = 4000
ymin = 700
ymax = 2740
zmax = 300

les_path = '/store/lmfa-2/acoustique-nl/simu_jules/LES/2T/C2/blue/'

# vcentral frequency band 
fc = [50, 63, 80, 100]
# number of simulation per frequency band
Nfc = [1,  1,  1,   1]
simu.set_frequencies(fc, Nfc)
# if set to true save the complete solution in the xz plane 
simu.set_input('side', False)

# define Acoustic domain inside the LES domain 
# if acoustic domain larger the LES fiels are copied outside of the domain box
simu.defineDomain(xmin, xmax, ymin, ymax, zmax)

# define sources height to compute 
# discretize a rotor of D = 120m with 3 source height in this case
Nh = 3
dh = 120/Nh
heights = np.arange(90-60+dh/2, 90+60, dh)
print(heights)
heights = np.rint(heights).astype(int)
print(heights)

# set planes where solution must be saved 
# ususlly ususefull to set the planes corresponding to tau=0 and tau = 180
# here we set the planes corresponding to tau =0,180,90,270
simu.set_xplanes([1000])
simu.set_yplanes([1720])

# define computational angle 
tau = np.array([0,30,60,90,120,150,180,210,240,270,300,330])


# Define the simulation in the Simu() class
# -----------------------------------------------------------------------------
# define the cases for the objects simu. 
# number of height, angles, turbines (as described in the LES simulations)
# ratio allows to scale the LES fields.
simu.defineCases(heights, tau, flow_path=les_path,src_path=src_path,ratio=1.04)

# Creates the files for the simulation
# distribute tau allows to compute angles in parallels, not very usefull for large wind farm. 
# the is an inherent parallelization on the turbines and height computed. 
# each case has uts own launch.sh file 
simu.distributeCases(distribute_tau=None, mem=6000,
                     time="0-01:00:00", turbine_index=None)
# show simulation time in days
print("Simu time = %s days"%(simu.computeSimulationTime()/60/60/24))

# create launch files to run several simulation in parallel in the same job
simu.createLaunchFiles(mem=6000,
                       time="0-01:00:00", turbine_index=None)

# Launch simulation from python
# Usually I launch the simulation from the terminal an check that every thing is 
# running fine
# ------------------------------------------------------------------------------
# launch cases using the files created in simu.createLaunchFiles
# simu.launchCases2()

# launch each cases in a different job using the launch.sh created in simu.distributeCases
# simu.launchCases(turbine_index=None)

# -----------------------------------------------------------------------------
simu.save()
