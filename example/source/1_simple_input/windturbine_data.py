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

wt.exploreDataBase()
