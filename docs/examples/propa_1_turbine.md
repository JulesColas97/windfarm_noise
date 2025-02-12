

### Launching the PE fortran code


The simulation are run using the Fortran code `src/kernel/New_PE_c/`. 
See the installation section for compilation instruction. 
The code is also briefly explained in the [WAPE](../WAPE.md) section.

An example of an input file is given here

```fortran 
$input
! general parameters 
!------------------------------------------------------------------------------
case_name='c0'
var0%c=343
var0%rho=1.2
var0%gamma=1.4

! Frequencies 
nb_freq=4
frequencies(1)=50,100,200,500

! Angles 
nb_theta=1
theta(1)=0

! Domain definition 
Lx1=3500.0
Lx2=500.0
Ly1=1000.0
Ly2=1000.0
Lz=300

! Numerical parameters 
!------------------------------------------------------------------------------
dx=0.5
cfl=0.1
! PML
size=30
param=50000.0
n=2.5
! impedance
imp%sigmae=50000.0
imp%alphae=100.0
imp%rigid=.false.

! Source 
src%pos_x=0.0
src%pos_z=90

! External Flow 
!------------------------------------------------------------------------------
! if .true. read from twente LES
external_flow=.true.
! if .true. constant flow equal to u0
uniform=.false.

! if .true. log profile with z0=0.002
logarithmic=.false.
u0=11

! Set the type of WAPE used
! if .true. Ostashev WAPE, else classic WAPE
arbitrary=.true.

! if true newa phase preserving WAPE
arbitrary_new=.false.

! (i am not sure this does somethind) 
interpolation=.true.

! path to LES data
fdir='/store/lmfa-2/acoustique-nl/simu_jules/LES/2T/C1/blue/output/'

! twente paramaters needed to rescale and positionned the flow 
fname='tavg'
tinput%ratio=1
tinput%z_i=1000.0
tinput%delta=0.25
tinput%T_scale_K=293.15
tinput%Lx=6.0
tinput%Ly=3.44
tinput%Lz=1.0
tinput%posx=1.0
tinput%posy=1.72
tinput%posz=0.04

! Set output paramaters 
!------------------------------------------------------------------------------
dout=1
nb_receiver=3
heights(1)=2.0, 10, 100
side=.true.

! does not doe anything 
top=.false.

! if continuation set to true resart from the last computed angles
! not yet a continuation with respect to frequency 
continuation=.false.

! Set the plane in the complete domain where the solution must be recorded
nb_xplane=1
xplane(1)=1000
nb_yplane=1
yplane(1)=1000
$end input
```


### Data visualization

Once the code has finished running

The PE results can be loaded using the Python class `PeResults`

```python
from prepost import PeResults
import matplotlib.pyplot as plt

# Inputs
# ------------------------------------------------------------------------------
pe = PeResults(casename='c0', iTurb=0, height=90, tau=0)
```

!!! warning
    The class needs a name case, a Turbine index a source height and a propagation angle 
    due to a previous implementation where the file name was defined from these quantities. 
    In the newest version the file name is directly passed to the read function.

The fortran code saves the $\Delta L$ (sound pressure level relative to the free field) either in the entire domain (2D propagation plane) or for a set of receivers at given heights.
In the first case we can plot the 2D cartography for the propagation angle $\tau=0$ with :
```python 
pe.read_carto(fname='./c0_0000.h5')
pe.plotSide(freq=500)
plt.clim(-10,10)
plt.show()
```


In the second case we can plot the $\Delta L$ as a function of $x$ for a given receiver height $z=2$ m with:
```python
pe.read_receiver('./c0_0000.h5')
pe.plotLine(freq=500, z=2)
plt.ylim(-10,10)
plt.show()
```


