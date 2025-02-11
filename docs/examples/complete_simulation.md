### Setting up the propagation simulation 


The propagation simulation is similar to the [two turbine](../examples/propa_2_turbines.md) example.

a typicall input.py file would be the following.

```python 
from prepost import Simu
from prepost import source as src
import sys
import numpy as np
import matplotlib.pyplot as plt

# Input for the simulation
# -----------------------------------------------------------------------------
simu = Simu('c0')
xmin = 500
xmax = 4000
ymin = 700
ymax = 2740
zmax = 300

# absolute path
les_path = '/<path>/<to>/<les>/<results>/'
src_path = '<path>/<to>/<src>/<folder>/kernel/New_PE_c/PE_2D_WAPE'

# vcentral frequency band
fc = [50, 63, 80, 100]
# number of simulation per frequency band
Nfc = [1,  1,  1,   1]
simu.set_frequencies(fc, Nfc)
# if set to true save the complete solution in the xz plane 
simu.set_input('side', False)

# define Acoustic domain inside the LES domain 
# if acoustic domain larger the LES fiels are copied outside of the domain box
simu.defineDomain(xmin, xmax, ymin, ymax, zmax)

# define sources height to compute 
# discretize a rotor of D = 120m with 3 source height in this case
Nh = 3
dh = 120/Nh
heights = np.arange(90-60+dh/2, 90+60, dh)
print(heights)
heights = np.rint(heights).astype(int)
print(heights)

# set planes where solution must be saved 
# ususlly ususefull to set the planes corresponding to tau=0 and tau = 180
# here we set the planes corresponding to tau =0,180,90,270
simu.set_xplanes([1000])
simu.set_yplanes([1720])

# define computational angle 
tau = np.array([0,30,60,90,120,150,180,210,240,270,300,330])

# Define the simulation in the Simu() class
# -----------------------------------------------------------------------------
# define the cases for the objects simu. 
# number of height, angles, turbines (as described in the LES simulations)
# ratio allows to scale the LES fields.
simu.defineCases(heights, tau, flow_path=les_path,src_path=src_path,ratio=1.04)

# Creates the files for the simulation
# distribute tau allows to compute angles in parallels, not very usefull for large wind farm. 
# the is an inherent parallelization on the turbines and height computed. 
# each case has uts own launch.sh file 
simu.distributeCases(distribute_tau=None, mem=6000,
                     time="0-01:00:00", turbine_index=None)
# show simulation time in days
print("Simu time = %s days"%(simu.computeSimulationTime()/60/60/24))

# create launch files to run several simulation in parallel in the same job
simu.createLaunchFiles(mem=6000,
                       time="0-01:00:00", turbine_index=None)

# Launch simulation from python
# Usually I launch the simulation from the terminal an check that every thing is 
# running fine
# ------------------------------------------------------------------------------
# launch cases using the files created in simu.createLaunchFiles
simu.launchCases2()

# -----------------------------------------------------------------------------
simu.save()
```

<hr>
### Setting up the source simulation 

To launch source simulation that match the propagation simulation a function in class `Simu`
is used to launch the different source simulation for each turbine of the wind farm. 

First the previously saved `Simu` object is loaded. 
```python
import numpy as np
from prepost import Simu
from prepost.source import WindTurbine, Mesh


iTurb = [0,1]
print('Making source ...')
print('------------------------------------------------------------------')
simu = Simu('c0')
simu.load('./c0.dat')
```

Then a wind turbine geometry is define using the `WindTurbine` class

```python 
wt = WindTurbine()
wt.default()
wt.Nblade = 3
wt.Nbeta = 12
wt.Lspan = 1.2*wt.Lspan
wt.seg = 1.2*wt.seg
wt.href = 90
```

!!! note 
    The span and hub height of the simulation must match the geometry defined for the LES, 
    and also the source height for which the PE simulation were computed.


Then a mesh is set up for the source simulation.
```python
ny = 60
nx = 90

x = np.linspace(simu.x1-100, simu.x2+100, nx)
y = np.linspace(simu.y1-100, simu.y2+100, ny)
z = np.array([2])

mesh = Mesh(polar=False)
mesh.set_cartesian_mesh(x, y, z)
```
here the boundary of the mesh are taken slightly larger than the propagation boundary. 

!!!note 
    The source mesh is very coarse here because away from the source there is not a lot of variation.
    A way to create a more refine mesh without increasing to much the computational cost is to use a polar mesh with increased number of point close to the source and cross wind. 
    Then interpolate this mesh to the same cartesian mesh as the propagation results.



Finally the simulation are launched 
```python 
print(simu.frequencies)

omega = 12.1 * np.pi * 2 / 60
omega = None
simu.makeSource(wt, mesh, offset=0, plot=False, iTurb=iTurb, omega=omega,
                Ncore=8, fname='./xy/spp_polar')
```

For large simulation run the case using an sbatch script such as:

```bash
#!/bin/bash
#SBATCH --job-name=source
#SBATCH --output=out.out # output messages go here
#SBATCH --error=err.err    # error messages go here
#SBATCH --mail-type=ALL
#SBATCH --partition=haswell # partition name
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64000
#SBATCH --time=0-01:00:00
module purge 
module load Python/3.10.8-GCCcore-12.2.0
export KMP_AFFINITY=granularity=fine,compact,1,0
export OMP_STACKSIZE=1g
limit -s unlimited

source venv/bin/activate

python input_source.py
```



<hr>
### Post processing
Once both the propagation and source simulation are performed, the post processing 
function can be used to compute the SPL, OASPL, AM etc ...
The different post processing functions can be found in [spl_process](../reference/spl_process.md)

First a check can be performed using the `Simu`class to verify that all propagation simulation were corrrectly run. 
```python 
from prepost import Source

workdir = './'
case = 'c0'

simu = pp.Simu(workdir + case)
simu.load(workdir + 'c0.dat')

print('checking delta L ...')
print('------------------------------------------------------------------')
simu.check_run_cases()
```

Then  the different PE simulations are concatenated and saved in one `DeltaLField` object for each turbine. 
```python 
from prepost import concatenate_all_dl

# -----------------------------------------------------------------------------
workdir = './'
case = 'c0'
path2Pe = workdir
iTurb = [0, 1]

# post process from PE h5 files
# -----------------------------------------------------------------------------
print('reading delta L ...')
print('------------------------------------------------------------------')
concatenate_all_dl(case, path2Pe, refine=4, iTurb=iTurb,
                          z=2, stepx=5,
                          dl_fname=workdir+'/xy/DL_polar')
```

!!! note
    here `nx` and `ny` are note specified hence the `DeltaL` fields are saved in polar coordinates system.


The $\Delta L$ and $S_{ff}$ are then added using the `conbine_dl_function`. 
The method used is described in [Colas Phd Thesis](https://acoustique.ec-lyon.fr/publi/colas_thesis.pdf#page=71)


```python 
from prepost import combine_dl_src
print('combining ...')
print('------------------------------------------------------------------')
combine_dl_src(case, path2Pe, iTurb=iTurb,
                  dl_fname=workdir+'/xy/DL_polar',
                  spp_fname=workdir+'/xy/spp_polar',
                  spl_fname=workdir+'/xy/spl_s_polar_',
                  polar=False,third=False)
```
Details on the function can be found [here](../reference/spl_process.md#src.prepost.spl_process.combine_dl_src).








