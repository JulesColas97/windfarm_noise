
First example allows to compute the noise generated from a wind turbine with asimple power profile as input for the flow. 


First you need to load the usefull libraries.
```python
from prepost.source import Windturbine
from prepost.source import Mesh
from prepost.source import Source
from prepost import Les
import numpy as np
import matplotlib.pyplot as plt
```


### Wind turbine geometry
Then we first define thegeometry of the wind turbine with the class `WindTurbine`.

```python
wt = WindTurbine()
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
```

The default windturbine is created. It is similar to the one used in Cotté et al (2019). 
But with aslightly wider diameter (50m radius instead of 48m). 

The number of blade is set to 3 and the number of discretization angle to 12. 
The angular step is equal to $\Delta \beta = 2\pi /$(Nbeta.Nblade).

We then reajust the diameter of the blade by scaling the span and the position of each segment to obtain a radius of 60m. 
Finaly the hub height is set a 90m. 

#### Mesh definition
After defining the wind turbine geometry the mesh on which the simulation must be performed is defined. 
Here we do not compute a sound power level (which would be indepedent of the position) but the sound pressure level in free field, which depends on the distance between the source and the receiver and the angle between each pair of receievr and sources. 

```python
ny = 40
nx = 20

x = np.linspace(-200, 200, nx)
y = np.linspace(-100, 100, ny)
z = np.array([2])

# create the Mesh object
mesh = Mesh(polar=False)
# set the mesh arrays 
mesh.set_cartesian_mesh(x, y, z)
mesh.cartesian_2_polar()
```
Here a cartesian mesh is defined from the three vectoir $x$, $y$, and $z$.
The code also allows to define a polar mesh with $x$, $\tau$ and $z$, where x is the radius and $\tau$ the propagation angle.
In order for the mesh to be used in the source model the coordinates must allways be converted to polar which is done at the last step.
To do so the wind turbine is allways considered to bee at position $(x=0,y=0)$.



### Atmospheric conditions
As previously stated the code developped in this project is tailored to work with 
Lage eddy simulation results from the PoF group at the university of twente. 
However it also possible to create your own flow data that you can input to the source model. 

The flow data must contain a cartesian mesh defined by three vector x, y, and z and the streamiwse componenent of the wind speed and turbulence dissipation rate as 2D matrix which size corresponds to (z,y,x). 

In this example we create a simple flow data to obtain a startified flow  (constant in x and y direction) with a power profile and a wind speed at hub height of 11.4 m/s. 
Here the turbulent dissipation rate is set constant equal to 0.01 m2/s3
```python
atmos = Les('./')
atmos.z = np.linspace(0,200,100)
atmos.y = np.array([-1000,1000])
atmos.x = np.array([-1000,1000])
Z,Y,X = np.meshgrid(atmos.z,atmos.y,atmos.x,indexing='ij')

U_hub = 11.5 # ms 
alpha = 0.4
epsilon = 0.01

atmos.u = U_hub * (Z/wt.href) ** alpha
print(atmos.u.shape)

# plt.plot(atmos.u[:,0,0],atmos.z)
atmos.epsilon = np.zeros_like(atmos.u) +  epsilon
```


### Source computation
Finally the source computation is perfomed using the class `Source`


```python
# define frequencies
frequencies = np.array([50,100,200,500,800,1000])

# set the rotational speed 
omega = 12.1 * 2 * np.pi / 60
wt.controlRotSpeed(U_hub, omega=omega)
#define the blade twist 
wt.setOptimalTwist(U_hub, 4)

print('U_hub', U_hub)
print('Omega', wt.omega*60/(2*np.pi))

src = Source(wt, atmos, mesh)

Ncore = 4
src.computeSpp(frequencies, Ncore)
src.mesh.polar_2_cartesian()

src.save("./swl.dat")
```
The data computed are saved in the file `swl.dat`


### Data visualisation 

To visualise the data previously computed first the data are loaded susing the `Source` class.
```python
src = Source()
src.load("./swl.dat")
```

Then a receiver position is choosen to plot the data
```python
x = 200
y = 0
ix = np.abs(src.mesh.x_array - x).argmin()
iy = np.abs(src.mesh.y_array - y).argmin()
```

The SWL (sound power  level) computed for each segment is logarithmicaly sum over all the segments and is averaged over the angular position of the rotor to obtain the total SWL of the turbine. 
```python
SWL_true = 10*np.log10(np.sum(10**(src.SWL[ix, 0, iy, ...]/10),
                  (0, 2))/src.Nbeta) + Aweight(src.frequencies)
```
!!! note
    here a A weighting is added to the spectrum with the function `Aweight()`

An other possibility is to compute the distance between the wind turbine hub and the receiver
```python 
x = src.mesh.x_array[ix]
y = src.mesh.y_array[iy]
R = np.sqrt((x-src.wt.absolute_pos[0])**2 + (y-src.wt.absolute_pos[1])**2)
```
then a pseudo Sound power level is computed from the the Sound pressure level in freefield. 
The differenrence with the previous sound power level (SWL_true) is that here the distance between the source and the receiver is the same for all segments.
In our case the computation is done for the trailing edge (TEN) and turbulent inflow (TIN) noise and for the total SPL.
```python
SWL = 10*np.log10(np.sum(10**(src.Spp[ix, 0, iy, ...]/10),
                         (0,2))/src.Nbeta)+Aweight(src.frequencies)+10*np.log10(4*np.pi*R*R)
TIN = 10*np.log10(np.sum(10**(src.Spp_tin[ix, 0, iy, ...]/10),
                         (0,2))/src.Nbeta)+Aweight(src.frequencies)+10*np.log10(4*np.pi*R*R)
TEN = 10*np.log10(np.sum(10**(src.Spp_ten[ix, 0, iy, ...]/10),
                         (0,2))/src.Nbeta)+Aweight(src.frequencies)+10*np.log10(4*np.pi*R*R)
```


Finally the spectrum at this receiver location can be plotted using mlatplotlib
```python
plt.plot(src.frequencies, TEN, 'r:')
plt.plot(src.frequencies, TIN, 'b:')
plt.plot(src.frequencies, SWL, 'k')
plt.plot(src.frequencies, SWL_true, 'g--')
plt.legend(['TEN','TIN','SWL','SWL_true'])
plt.xscale('log')

plt.show()
```
