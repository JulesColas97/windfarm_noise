import prepost as pp
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
ny = 40
nx = 20

x = np.linspace(-200, 200, nx)
y = np.linspace(-100, 100, ny)
z = np.array([2])

# create the Mesh object
mesh = pp.source.Mesh(polar=False)
# set the mesh arrays 
mesh.set_cartesian_mesh(x, y, z)
mesh.cartesian_2_polar()


# Create the flow object 
# for now this is a bit complicated because the code only manages Les type input
# so to enter a simple profile or even a constant value you still need to input 
# an entire flow field
# -----------------------------------------------------------------------------
atmos = pp.Les('./')
atmos.z = np.linspace(0,200,100)
atmos.y = np.array([-1000,1000])
atmos.x = np.array([-1000,1000])
Z,Y,X = np.meshgrid(atmos.z,atmos.y,atmos.x,indexing='ij')
print(Z.shape)

U_hub = 11.5 # ms 
alpha = 0.0
epsilon = 0.01

atmos.u = U_hub * (Z/wt.href) ** alpha
print(atmos.u.shape)

# plt.plot(atmos.u[:,0,0],atmos.z)
atmos.epsilon = np.zeros_like(atmos.u) +  epsilon



# define frequencies 
frequencies = np.array([50,100,200,500,800,1000])


Ncore = 4
# Compute sound power 
# -----------------------------------------------------------------------------

omega = 12.1 * 2 * np.pi / 60
wt.controlRotSpeed(U_hub, omega=omega)
wt.setOptimalTwist(U_hub, 4)

print('U_hub', U_hub)
print('Omega', wt.omega*60/(2*np.pi))

src = pp.Source(wt, atmos, mesh)
src.computeSpp(frequencies, Ncore)
src.mesh.polar_2_cartesian()

src.save("./swl_constant.dat")
