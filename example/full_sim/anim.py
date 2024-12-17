from prepost import Les
from prepost import Simu
from prepost import SplField
from prepost import source as src
import sys
import numpy as np
import matplotlib.pyplot as plt

workdir = './'

spl = SplField()
spl.load(workdir +'/xy/spl_t_cart_0.dat',seg=True,time=True,oaspl=False)
# spl.plot_spl(z=2, freq=100, )
spl.info()
print(spl.frequencies)
spl.compute_time_average()
anim = spl.anim_spl(time=True, OA=False, z=2, freq=100, vmin=8, vmax=44)
# plt.gca().set_aspect('equal',"box")

spl = SplField()
spl.load(workdir +'/xy/spl_t_cart_1.dat',seg=True,time=True,oaspl=False)
# spl.plot_spl(z=2, freq=100, )
spl.info()
spl.compute_time_average()
anim2 = spl.anim_spl(time=True, OA=False, z=2, freq=100, vmin=8, vmax=44)


plt.show()
