from prepost import PeResults
import numpy as np
import matplotlib.pyplot as plt

# Inputs
# ------------------------------------------------------------------------------
dirname = '/home/lmfa/jcolas/Documents/DEV/wf_phd/simu/windfarm/16T/'
pe = PeResults(casename='c0', iTurb=0, height=90, tau=0)
# pe.read_carto('./c0_0000.h5')
pe.read_receiver('./c0_0000.h5')

plt.figure()
pe.plotLine(freq=1000, z=2)
# plt.clim(-10, 10)
plt.show()
