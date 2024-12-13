from prepost import DeltaLField
import matplotlib.pyplot as plt

# First turbine
#------------------------------------------------------------------------------
plt.subplot(211)
dl = DeltaLField()
dl.load('./xy/DL_polar0.dat')
dl.plot_tau(tau=0,z=2,height=90,freq=50,color='k')
dl.plot_tau(tau=180,z=2,height=90,freq=50,color='k')


dl = DeltaLField()
dl.load('./xz/DL0.dat')
dl.plot_x(y=1480,z=2,height=90,freq=50,color='b',linestyle='--')
plt.ylim(-8,8)
plt.xlim(500,3000)

# second turbine
#------------------------------------------------------------------------------
plt.subplot(212)
dl = DeltaLField()
dl.load('./xy/DL_polar1.dat')
dl.plot_tau(tau=0,z=2,height=90,freq=50,color='k')
dl.plot_tau(tau=180,z=2,height=90,freq=50,color='k')

dl = DeltaLField()
dl.load('./xz/DL1.dat')
dl.plot_x(y=1960,z=2,height=90,freq=50,color='b',linestyle='--')
plt.ylim(-8,8)
plt.xlim(500,3000)
plt.show()
