import prepost as pp
from prepost import Simu
from prepost import concatenate_all_dl
from prepost import concatenate_side_dl
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
import logging

# workdir = '/home/lmfa/jcolas/Documents/DEV/wf_phd/simu/windfarm/2T/C2/'
workdir = './'
savdir = workdir + '/xy/figure/'
case = 'c0'
path2Pe = workdir
# to choose for which turbine the concatenation occurs
iTurb = [0]

simu = Simu(workdir + case)
simu.load(workdir + 'c0.dat')
print(simu.frequencies.shape)
simu.check_run_cases()

# post process from PE h5 files to create top view
# -----------------------------------------------------------------------------
print('reading delta L ...')
print('------------------------------------------------------------------')
concatenate_all_dl(case, path2Pe, refine=5, iTurb=iTurb,
                      z=2, stepx=5,
                      dl_fname=workdir+'/xy/DL_polar',
                      spl_fname=workdir+'/xy/spl_polar',
                      spp_fname=workdir+'/xy/spp')
