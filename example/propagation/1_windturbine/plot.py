from prepost import PeResults
import numpy as np
import matplotlib.pyplot as plt

# Inputs
# ------------------------------------------------------------------------------
pe = PeResults(casename='c0', iTurb=0, height=90, tau=0)

pe.plotFlowSide(fname='./flow.h5', cmap='RdYlBu_r')
pe.read_carto('./c0_0000.h5')
# pe.read_receiver('./c0_0000.h5')

plt.figure()
pe.plotSide(freq=50)
# pe.plotLine(freq=500, z=2)
plt.clim(-10, 10)
plt.show()
