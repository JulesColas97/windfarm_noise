from prepost import Les
from prepost import Simu
from prepost import source as src
import sys
import numpy as np
import matplotlib.pyplot as plt

path = '/home/lmfa/jcolas/Documents/DEV/LES/windfarm/2T/C/blue/'
path = '/store/lmfa-2/acoustique-nl/simu_jules/LES/2T/C2/blue/'
les = Les(path)
les.read(ratio=1.04)
les.plotProfile('u', 500, 1720, color="b", linestyle=':')
les.plotProfile('u', 800, 1720, color="b", linestyle='--')
les.plotProfile('u', 1200, 1720, color="b", linestyle='-.')
plt.plot([11.4], [90], '+r')
les.dissipation()
plt.figure(2)
les.plotProfile('epsilon', 500, 1000, color="b", linestyle=':')
les.plotProfile('epsilon', 1200, 1000, color="b", linestyle='--')
les.plotProfile('epsilon', 2000, 1000, color="b", linestyle='-.')
plt.figure()
les.plotQuantityXY('u',90,cmap='coolwarm')

path = '/store/lmfa-2/acoustique-nl/simu_jules/LES/2T/C2/red/'
les = Les(path)
les.read(ratio=1.04)
plt.figure(1)
les.plotProfile('u', 500, 1720, color="r", linestyle=':')
les.plotProfile('u', 800, 1720, color="r", linestyle='--')
les.plotProfile('u', 1200, 1720, color="r", linestyle='-.')

plt.show()
quit()
