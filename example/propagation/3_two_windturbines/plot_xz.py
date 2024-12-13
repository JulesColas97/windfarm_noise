from prepost import DeltaLField
import matplotlib.pyplot as plt

dl = DeltaLField()
dl.load('./xz/DL0.dat')
print(dl.y_array)

plt.subplot(211)
dl.plot_xz(y=1480,freq=100,height=90,cmap="RdBu_r")
plt.clim(-10,10)


dl = DeltaLField()
dl.load('./xz/DL1.dat')
print(dl.y_array)

plt.subplot(212)
dl.plot_xz(y=1960,freq=100,height=90,cmap="RdBu_r")
plt.clim(-10,10)

plt.show()
