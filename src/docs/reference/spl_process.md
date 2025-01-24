# Reference for `prepost/spl_process.py`


These functions are designed to process and analyze simulation
data from wind turbine and especially wind farm noise. 
They facilitate the concatenation, refinement, and interpolation of delta L fields
for multiple turbines.
Functions like `concatenate_angles_dl` and `concatenate_all_dl` are used to combine angle-specific
data from multiple turbines, refine the angular resolution, and interpolate the results onto a common mesh.
Other functions, such as `concatenate_side_dl` and `concatenate_planes_dl`, 
are used for creating side views and plane-specific analyses.
Additionally, `combine_dl_src` and `convert_to_receiver_time` integrate delta L fields with source fields
and converts SPL as function of the rotor position into SPL as a function of receiver time. 



<br>

## ::: src.prepost.spl_process.concatenate_angles_dl

<br><br><hr><br>

## ::: src.prepost.spl_process.init_deltaL_field

<br><br><hr><br>

## ::: src.prepost.spl_process.add_1Pe_to_deltaL_field

<br><br><hr><br>

## ::: src.prepost.spl_process.concatenate_angles

<br><br><hr><br>

## ::: src.prepost.spl_process.concatenate_all_dl

<br><br><hr><br>

## ::: src.prepost.spl_process.concatenate_side_dl

<br><br><hr><br>

## ::: src.prepost.spl_process.combine_dl_src

<br><br><hr><br>

## ::: src.prepost.spl_process.convert_to_receiver_time

<br><br><hr><br>

## ::: src.prepost.spl_process.concatenate_planes_dl

<br><br><hr><br>

## ::: src.prepost.spl_process.postprocessingPlanesCombine

<br><br>
