import prepost as pp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
import logging

# logging.basicConfig(level=logging.DEBUG, force=True)
logging.basicConfig(format='%(levelname)s: %(message)s', force=True)

iTurb = [0,1]
print('Making source ...')
print('------------------------------------------------------------------')
simu = pp.Simu('c0')
simu.load('./c0.dat')

wt = pp.source.WindTurbine()
wt.default()
wt.Nblade = 3
wt.Nbeta = 12
wt.Lspan = 1.2*wt.Lspan
wt.seg = 1.2*wt.seg
wt.href = 90

ny = 60
nx = 90

x = np.linspace(simu.x1-100, simu.x2+100, nx)
y = np.linspace(simu.y1-100, simu.y2+100, ny)
z = np.array([2])

mesh = pp.source.Mesh(polar=False)
mesh.set_cartesian_mesh(x, y, z)
print(simu.frequencies)

omega = 12.1 * np.pi * 2 / 60
omega = None
simu.makeSource(wt, mesh, offset=0, plot=False, iTurb=iTurb, omega=omega,
                Ncore=16, fname='./xy/spp_polar')
