import prepost as pp
from prepost.source.utils import interp_3D_atmos_data
import numpy as np
import matplotlib.pyplot as plt


# Create the object wind turbine
# -----------------------------------------------------------------------------
wt = pp.source.WindTurbine()
# set the parameter to the default turbine 
# the default corresponds to the turbine define in Cotte et al . 
wt.default()
# set the number of blades
wt.Nblade = 3
# set the number of angles to descretize 1/Nblade of the rotor
wt.Nbeta = 12
# Resest the span to reach R=60m
wt.Lspan = 1.2*wt.Lspan
wt.seg = 1.2*wt.seg
# Set th wind turbine hub 
wt.href = 90


# define the computational domain 
# the wind turbine is assumed to be at (0,0)
# -----------------------------------------------------------------------------
ny = 20
nx = 40

x = np.linspace(-200, 1000, nx)
y = np.linspace(-100, 100, ny)
z = np.array([2])

# create the Mesh object
mesh = pp.source.Mesh(polar=False)
# set the mesh arrays 
mesh.set_cartesian_mesh(x, y, z)
mesh.cartesian_2_polar()

# Create the flow object
# read LES file
# -----------------------------------------------------------------------------
les = pp.Les('/store/lmfa-2/acoustique-nl/simu_jules/LES/2T/C1/blue/')
les.read(ratio=1)
les.dissipation()

offset = 0
xS = np.array([[les.turbines[0,1]*les.z_i-offset]])
yS = np.array([[les.turbines[0,2]*les.z_i]])
zS = np.array([[wt.href]])
U_hub,epsilon_hub = interp_3D_atmos_data(les,xS,yS,zS)

wt.controlRotSpeed(U_hub)
wt.setOptimalTwist(U_hub, 4)
wt.absolute_pos = (les.turbines[0,1]*les.z_i,
                   les.turbines[0,2]*les.z_i)
print('U_hub', U_hub)
print('epsilon_hub', epsilon_hub)
print('Omega', wt.omega*60/(2*np.pi))

# define frequencies 
# -----------------------------------------------------------------------------
frequencies = np.array([50,100,200,500,800,1000])

# Compute sound power
# -----------------------------------------------------------------------------
Ncore = 4
src = pp.Source(wt, les, mesh)
src.computeSpp(frequencies, Ncore)
src.mesh.polar_2_cartesian()

src.save("./swl.dat", atmos_data=False)
