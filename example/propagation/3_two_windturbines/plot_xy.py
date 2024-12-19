from prepost import DeltaLField
import matplotlib.pyplot as plt


dl = DeltaLField()
dl.load('./xy/DL_polar0.dat')

plt.subplot(211)
dl.plot_top_cart(z=2, freq=50, height=90, cmap='RdBu_r')
plt.clim(-8, 8)
plt.xlim(0, 3000)
plt.ylim(720, 2720)
plt.gca().set_aspect('equal', adjustable='box')
# les_path = '../../les_data/2TSBL/blue/'

dl.load('./xy/DL_polar1.dat')
plt.subplot(212)
dl.plot_top_cart(z=2, freq=50, height=90, cmap='RdBu_r')
plt.clim(-8, 8)
plt.xlim(0, 3000)
plt.ylim(720, 2720)
plt.gca().set_aspect('equal', adjustable='box')
plt.show()
