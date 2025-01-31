from prepost import Les
from prepost import Simu
from prepost import SplField
from prepost import source as src
import sys
import numpy as np
import matplotlib.pyplot as plt

workdir = './'

spl = SplField()
spl.load(workdir +'/xy/spl_s_cart_0.dat',seg=True,time=False,oaspl=False)
# spl.plot_spl(z=2, freq=100, )
spl.info()
print(spl.frequencies)
plt.figure()
spl.compute_time_average()
plt.subplot(211)
spl.plot_mean(time=True,OA=False,z=2,freq=100)
plt.clim(8,44)
plt.gca().set_aspect('equal',"box")



spl = SplField()
spl.load(workdir +'/xy/spl_s_cart_1.dat',seg=True,time=False,oaspl=False)
# spl.plot_spl(z=2, freq=100, )
spl.info()
spl.compute_time_average()
plt.subplot(212)
spl.plot_mean(time=False,OA=False,z=2,freq=50)
plt.clim(8,44)
plt.gca().set_aspect('equal',"box")
plt.show()
