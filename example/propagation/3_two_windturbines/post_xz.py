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
case = 'c0'
path2Pe = workdir
iTurb = [0, 1]

simu = Simu(workdir + case)
simu.load(workdir + 'c0.dat')
print(simu.frequencies.shape)

simu.check_run_cases()


# create side view fields
# -----------------------------------------------------------------------------
pp.concatenate_side_dl(case, path2Pe, iTurb=iTurb,
            tau=[0,180], nx=900, nz=151, dl_fname=workdir+'/xz/DL')

