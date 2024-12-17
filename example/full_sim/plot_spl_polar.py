from prepost import Les
from prepost import Simu
from prepost import SplField
from prepost import source as src
import sys
import numpy as np
import matplotlib.pyplot as plt

workdir = './'

spl = SplField()
print('toto')
spl.load(workdir +'/xy/spl_s_polar_0.dat',seg=True,time=False,oaspl=False)
# spl.plot_spl(z=2, freq=100, )
spl.info()
plt.figure()
spl.compute_time_average()
plt.pcolormesh(spl.x_grid[:,:,0],
               spl.y_grid[:,:,0],
               spl.mean_seg[:,:,0,0])
plt.xlabel('x (m)')
plt.ylabel('y (m)'
plt.show()
