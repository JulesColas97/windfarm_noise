from prepost import Source 
from prepost.utils import Aweight
import numpy as np
import matplotlib.pyplot as plt

src = Source()
src.load("./swl.dat")

x = 200
y = 0
ix = np.abs(src.mesh.x_array - x).argmin()
iy = np.abs(src.mesh.y_array - y).argmin()

x = src.mesh.x_array[ix]
y = src.mesh.y_array[iy]
R = np.sqrt((x-src.wt.absolute_pos[0])**2 + (y-src.wt.absolute_pos[1])**2)

# print(iy)
# print(src.SWL.shape)
# print(src.Spp_tin.shape)
SWL_true = 10*np.log10(np.sum(10**(src.SWL[ix, 0, iy, ...]/10),
                  (0, 2))/src.Nbeta) + Aweight(src.frequencies)

SWL = 10*np.log10(np.sum(10**(src.Spp[ix, 0, iy, ...]/10),
                         (0,2))/src.Nbeta)+Aweight(src.frequencies)+10*np.log10(4*np.pi*R*R)
TIN = 10*np.log10(np.sum(10**(src.Spp_tin[ix, 0, iy, ...]/10),
                         (0,2))/src.Nbeta)+Aweight(src.frequencies)+10*np.log10(4*np.pi*R*R)
TEN = 10*np.log10(np.sum(10**(src.Spp_ten[ix, 0, iy, ...]/10),
                         (0,2))/src.Nbeta)+Aweight(src.frequencies)+10*np.log10(4*np.pi*R*R)


plt.plot(src.frequencies, TEN, 'r:')
plt.plot(src.frequencies, TIN, 'b:')
plt.plot(src.frequencies, SWL, 'k')
plt.plot(src.frequencies, SWL_true, 'g--')
plt.legend(['TEN','TIN','SWL','SWL_true'])
plt.xscale('log')

plt.show()
