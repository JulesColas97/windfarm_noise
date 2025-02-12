
### Setting up and launching the simulations

To set up multiple propagation simulation from the LES results we first need to create a simulation object
```python
from prepost import Simu
# create simu object
# -----------------------------------------------------------------------------
simu = Simu('c0')
```
When creating this object a default dictionary of the input needed for the Fortran code is created. 
We can then modify the different inputs for our specific case: 

```python
# define frequency content
# -----------------------------------------------------------------------------
fc = [50, 100, 200]
simu.set_frequencies(fc)

#  define output parameters
# -----------------------------------------------------------------------------
simu.set_input('side', True)
#set plane to save solution
# same referentiel as LES sim
simu.set_xplanes([1000])
simu.set_yplanes([1720])

# define domain
# -----------------------------------------------------------------------------
xmin = 0
xmax = 1000 + 2000
ymin = 1720 - 1000
ymax = 1720 + 1000
zmax = 300
simu.defineDomain(xmin, xmax, ymin, ymax, zmax)
```

!!! note 
    In the `Simu` class the domain is rectangular with the same coordinates system as the LES results
    For the propagation simulation the domain is rectangular, and the wind turbine position is at the origin. This translation is performed when creating the input.nml files for each Fortran simulation.


We can then define the source heights used to discretize the wind turbine rotor and the propagation angles:

```python
# discretize the rotor with equivalent source height
R = 60
hub = 90
Nh = 3
dh = 2*R/Nh
height = np.arange(hub-R+dh/2, hub+R, dh)
heights = np.rint(height).astype(int)

tau = np.array([0,20,40,60,80, 90,100,120,140,160,180,200,220,240,260,270])
```

Finally, the path to the LES results used as input for the simulation and the absolute path to the Fortran executable must be given to the `defineCases` function
```python 
# define cases with height, angles and flow 
# needs to be absolute path
les_path = '<absolute>/<path>/LES/2T/S2/blue/'
src_path = '<absolute>/<path>/windfarm_noise/src/kernel/New_PE_c/PE_2D_WAPE'
simu.defineCases(heights, tau, flow_path=les_path,
                 src_path=src_path, ratio=1.0)
```

Once the cases are defined in the `Simu` object we can create all the input.nml files need to run the cases. 
```python
simu.distributeCases(distribute_tau=None, mem=1000,
                     time="01:00:00", turbine_index=None)
```

The cases are organized in a tree with a folder for each turbine simulation, a subfolder for each source height and finally each propagation angle is store in the different file which contains all the frequencies computed. 

```
├── t0
│   ├──<height_0>
│   │       c0_0000.dat
│   │       c0_0020.dat
│   │     
│   ├──<height_1>
│   │       c0_0000.dat
│   │       c0_0020.dat
│   │     
│   ├──<height_1>
│   │       c0_0000.dat
│   │       c0_0020.dat
│   │     
├── t1
│   ├──<height_0>
│   │       c0_0000.dat
│   │       c0_0020.dat
│   │     
│   ├──<height_1>
│   │       c0_0000.dat
│   │       c0_0020.dat
│   │     
│   ├──<height_1>
│   │       c0_0000.dat
│   │       c0_0020.dat
│   │     
```

In each subfolder a `launch.sh` and `input.nml` file is present to run the simulation.
To cases can then be submitted manually with command
```bash
sbacth launch.sh
```
or all at once with the python function 

```python 
simu.launchCases(turbine_index=None)
```

!!! warning
    In the case of large wind farm this can result in a very large number of jobs being submitted at the same which can creates problem in SLURM queue. 
    To avoid this, the following method is preferable. 

After creating all the input files the function can be used 

```python 
simu.createLocalLaunchFiles(mem=6000, time="01:00:00", turbine_index=None)
```
This regroups several cases in one node with respect to the maximum memory allocated per cases (here 6Gb).

!!! warning
    This method was tailored to run on the haswell partition of the PMCS2I cluster, the maximum size of a haswell node (64Gb) is hard coded in the function

This function creates several launch<index>.sh files in the parent directory which should reduce largely the number of SLURM jobs. 

The job can then be launched manually or with the command 
```python
simu.launchCases2()
```

Finally, the `Simu` can be saved to be use in the post-processing routines. 

```python
simu.save()
```

<hr>
### Post processing for top view

To post process the simulations we can either use the `PeResults` to read individual propagation planes or we can use the `DeltaLfield` class which allows concatenating all propagation angle and source height into one matrix.

```python 
from prepost import concatenate_all_dl
```

First we can check that all simulation were run properly 
```python
from prepost import Simu
workdir = './'
case = 'c0'
simu = Simu(workdir + case)
simu.load(workdir + 'c0.dat')
simu.check_run_cases()
```
!!! warning
    The `check_run_cases` only check that the output files exists. It does not check if all frequency were computed or if the files are not corrupted. 


A wrapper for the class is also written which allows to read all the simulation as defined in the object `Simu` and save it in a file as a `DeltaLField` object. 

```python
from prepost import concatenate_all_dl
path2Pe = './'
case = 'c0'
iTurb = [0]
concatenate_all_dl(case, path2Pe, refine=5, iTurb=iTurb,
                      z=2, stepx=5,
                      dl_fname=workdir+'/xy/DL_polar')
```

This function reads the `Simu` object, then according to the case defined it reads all the PE receiver results and concatenate the propagation fields.
Details of the function can be found [here](../reference/spl_process.md#src.prepost.spl_process.concatenate_all_dl).


<hr>

To visualize the data first load the [deltaLField](../reference/deltaLfield.md) object:
Then we can use the `plot_top_cartesian` function to plot a top view of the results. 


```python
from prepost import DeltaLField
import matplotlib.pyplot as plt

dl = DeltaLField()
dl.load('./xy/DL_polar0.dat')

plt.subplot(211)
dl.plot_top_cart(z=2, freq=200, height=130, cmap='RdBu_r')
plt.clim(-8, 8)
plt.xlim(0, 3000)
plt.ylim(720, 2720)
plt.gca().set_aspect('equal', adjustable='box')

dl.load('./xy/DL_polar1.dat')
plt.subplot(212)
dl.plot_top_cart(z=2, freq=50, height=90, cmap='RdBu_r')
plt.clim(-8, 8)
plt.xlim(0, 3000)
plt.ylim(720, 2720)
plt.gca().set_aspect('equal', adjustable='box')
plt.show()
```


<br><hr>
### Post-processing for side view 

Create a side view of the delta L field a wrapper is also implemented. 
The `concatenate_side_dl` function allows to concatenate PE results to obtain a 2D field in 
a plane of propagation. 
The fields are interpolated on a grid `(nx,nz)`.
Note that in this example two opposite propagation planes are concatenated.
Another set of angle could be `[20 ,200]` or just one angle. 
Note that order is important the second propagation angle will be reverse so that the increasing x direction is in the first angles given direction. 


```python 
from prepost import concatenate_side_dl
# create side view fields
# -----------------------------------------------------------------------------
pp.concatenate_side_dl(case, path2Pe, iTurb=iTurb,
            tau=[0,180], nx=900, nz=151, dl_fname=workdir+'/xz/DL')

```
After this post-processing, the Delta is store in the same variable as for the top view.
The field is a ND matrix of size `(nx, ny, nz, Nfreq, Nheight)` with `ny=1`. 
The same set of plotting function or interpolating function of the `DeltaLField`class can then be used to visualize the data. 

```python 
from prepost import DeltaLField
import matplotlib.pyplot as plt

dl = DeltaLField()
dl.load('./xz/DL0.dat')

plt.subplot(211)
dl.plot_xz(y=1480,freq=100,height=130,cmap="RdBu_r")
plt.clim(-10,10)


dl = DeltaLField()
dl.load('./xz/DL1.dat')
print(dl.y_array)

plt.subplot(212)
dl.plot_xz(y=1960,freq=100,height=130,cmap="RdBu_r")
plt.clim(-10,10)

plt.show()
```


