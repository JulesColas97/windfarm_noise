"""
`wind_turbine.Windtrubine`
`prepost.source.func_repo_XFOIL.run_Xvfb`
`prepost.source.wind_turbine.Windturbine`

`source.func_repo_XFOIL.run_Xvfb`
# Contains:

- `prepost.source.wind_turbine.WindTurbine`: a class able to define a wind turbine geometry, and BL quantities.
- `prepost.source.mesh.Mesh`: a class to define the mesh used to comput Spp.
- `prepost.source.main.Source`: a class to create $Spp_{ff}$ field accroding to ABL, wind turbine géometry and mesh definition.

"""
from .wind_turbine import WindTurbine
from .atmos import Atmosphere
from .mesh import Mesh
from .main import Source



# from .func_repo_XFOIL import XFOIL, run_Xvfb,kill_Xvfb
# from .utils import *
# from .function import *
# from .main import *
