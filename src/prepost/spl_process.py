"""
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
"""

from .spl import SplField
from .pre import Simu
from .post import Simu, DeltaLField
from .wape import PeResults
import numpy as np

import matplotlib.pyplot as plt
import logging
from .source.main import Source

def concatenate_angles_dl(casename: str, path2Pe: str, refine: int = 2, nx: int = 350, ny: int = 320, iTurb: np.ndarray = None, plot: bool = False, dl_fname: str = None, spp_fname: str = None, spl_fname: str = None) -> None:
    """
    Concatenates angles for delta L field from simulation results.

    This function reads simulation parameters, initializes the delta L field,
    and processes the data to concatenate angles. It optionally plots the results
    and saves the concatenated delta L field to a file.
    WARNING: this function is not used that much and `concatenate_all_dl` is the prefered method.  

    Args:
        casename (str): The name of the simulation case.
        path2Pe (str): The path to the PE results.
        refine (int, optional): The number of refinement steps (artificially increase the number of propagation angle). Defaults to 2.
        nx (int, optional): The number of points in the x-direction for cartesian interpolation. Defaults to 350.
        ny (int, optional): The number of points in the y-direction for cartesian interpolation. Defaults to 320.
        iTurb (np.ndarray, optional): The indices of the turbines to process. Defaults to None.
        plot (bool, optional): Whether to plot the results. Defaults to False.
        dl_fname (str, optional): The filename to save  the delta L field. Defaults to None.
        spp_fname (str, optional): The filename to save the Spp field. Defaults to None.
        spl_fname (str, optional): The filename to save the Spl field. Defaults to None.

    Returns:
        None
    """
    print('Start concatenating angles ...')
    if dl_fname is None :
        dl_fname = 'DL'
    if spp_fname is None:
        spp_fname = 'Spp'
    if spl_fname is None:
        spl_fname = 'Spl'

    # load simulation parameters
    simu = Simu(casename)
    simu.load(path2Pe+casename+'.dat')

    if iTurb is None:
        iTurb = np.arange(0,len(simu.tx))
    # if not trubine index given loop over all turbines 
    for ii in iTurb:
        #--------------------------------------------------------------
        # create delta L field (read from Pe results)
        deltaL = DeltaLField(path2Pe,casename)
        for h in simu.heights:
            for t in simu.tau:
                logging.info("h,t= %s,%s" % (h, t))
                deltaL.read_receiver(ii,h,t,simu.distribute_tau)
        flag = deltaL.check_compatibility()
        if flag == -1:
            logging.warning('error in check compatibility')
            quit()

        deltaL.concatenate_angles()
        deltaL.loop_angle()
        for kk in range(refine):
            deltaL.refine_angle()

        deltaL.shift_domain(simu.tx[ii]*simu.les.z_i,simu.ty[ii]*simu.les.z_i)
        if plot:
            deltaL.plotTopRaw(2,100,100)
            plt.xlim(simu.x1,simu.x2)
            plt.ylim(simu.y1,simu.y2)
            plt.clim(-10,10)
            plt.show()
        deltaL.xS = simu.tx[ii]*simu.les.z_i
        deltaL.yS = simu.ty[ii]*simu.les.z_i

       # define common mesh 
        x = np.linspace(simu.x1,simu.x2,nx)
        y = np.linspace(simu.y1,simu.y2,ny)

        deltaL.interpolate_from_polar(x,y)

        logging.info("saving delta L as :" + dl_fname+str(ii)+'.dat')
        deltaL.save(dl_fname+str(ii)+'.dat')

        # reset fields
        deltaL = None
        spl = None
        src = None
    print('done.')



def init_deltaL_field(
    path2Pe: str,
    casename: str,
    simu: Simu,
    iTurb: int
) -> DeltaLField:
    """
    Initializes the delta L field for a given turbine.

    This function sets up the delta L field with the necessary parameters and
    dimensions based on the simulation data `Simu`.

    Args:
        path2Pe (str): The path to the Pe results.
        casename (str): The name of the simulation case.
        simu (Simu): The simulation object containing parameters.
        iTurb (int): The index of the turbine.

    Returns:
        DeltaLField: The initialized delta L field.
    """
    deltaL = DeltaLField(path2Pe,casename)

    Lx1 = simu.x2-simu.tx[iTurb]*simu.les.z_i
    Lx2 = simu.tx[iTurb]*simu.les.z_i - simu.x1 
    Ly1 = simu.y2 - simu.ty[iTurb]*simu.les.z_i
    Ly2 = simu.ty[iTurb]*simu.les.z_i - simu.y1

    xmax = np.sqrt(max(Lx1,Lx2)**2 + max(Ly1,Ly2)**2) + 50
    dx = simu.inputs["dout"]
    deltaL.x = np.arange(0, xmax + dx , dx)
    deltaL.nx = np.size(deltaL.x)

    deltaL.ntau = simu.tau.size
    deltaL.tau = simu.tau
    deltaL.height = np.array(simu.heights)
    deltaL.nheight = len(simu.heights)

    deltaL.nz = simu.inputs["nb_receiver"]
    deltaL.z =  np.array(simu.inputs["heights(1)"])

    deltaL.Nfreq = simu.frequencies.size
    deltaL.frequencies = simu.frequencies

    # create matrice for deltaL polar 
    deltaL.deltaL_polar = np.zeros(
        (deltaL.nx, deltaL.nz, deltaL.ntau, deltaL.Nfreq, deltaL.nheight)
    )
    return deltaL


def add_1Pe_to_deltaL_field(
    deltaL: DeltaLField,
    iTurb: int,
    height: float,
    tau: float
) -> None:
    """
    Adds Pe results to the delta L field for a given turbine, height, and tau.

    This function reads the receiver data from the Pe results and updates the
    `deltaL_polar` field accordingly.

    Args:
        deltaL (DeltaLField): The delta L field object.
        iTurb (int): The index of the turbine.
        height (float): The height at which to add the Pe results.
        tau (float): The tau value at which to add the Pe results.

    Returns:
        None
    """
    res = PeResults(
                deltaL.casename,
                iTurb,
                height,
                tau,
                dirname=deltaL.dirname,
                distribute_tau=deltaL.simu.distribute_tau,
            )

    fname = deltaL.fname(iTurb, tau, height, deltaL.simu.distribute_tau)
    res.read_receiver(fname)
    # find height and tau index of the simulation
    if (res.frequencies!=deltaL.simu.frequencies).any():
        logging.warning("frequencies not the same")
    if (res.heights!=deltaL.z).any():
        logging.warning("heights not the same")

    iheight = np.argmin(abs(deltaL.height - height))
    itau = np.argmin(abs(deltaL.tau - tau))
    # find corresponding index
    ixmin = np.argmin(abs(deltaL.x - res.x[0]))
    ixmax = np.argmin(abs(deltaL.x - res.x[-1]))
    logging.info("height = " + str(height) + ", tau = " + str(tau))
    deltaL.deltaL_polar[ixmin : ixmax + 1, :, itau, :, iheight] = res.receiver
    # copy last values
    deltaL.deltaL_polar[ixmax+1:, :, itau, :, iheight] = deltaL.deltaL_polar[ixmax, :, itau, :, iheight][None,:,:,]

def concatenate_angles(
    casename: str,
    path2Pe: str,
    simu: Simu,
    iTurb: int,
    refine: int,
    z: float = None,
    stepx: int = None
) -> DeltaLField:
    """
    Concatenates angles for the delta L field for a given turbine.

    This function initializes the delta L field, adds Pe results, and processes
    the data to concatenate angles. It optionally refines the angles and selects
    a specific height from the receiver heights. It can also increase the grid size in the radius direction.

    Args:
        casename (str): The name of the simulation case.
        path2Pe (str): The path to the Pe results.
        simu (Simu): The simulation object containing parameters.
        iTurb (int): The index of the turbine.
        refine (int): The number of refinement steps.
        z (float, optional): The height to select. Defaults to None.
        stepx (int, optional): The step size for refinement. Defaults to None.

    Returns:
        DeltaLField: The concatenated delta L field.
    """
    # create delta L field (read from Pe results)
    deltaL = init_deltaL_field(path2Pe, casename, simu, iTurb)

    for h in simu.heights:
        for t in simu.tau:
            add_1Pe_to_deltaL_field(deltaL, iTurb, h, t)

    deltaL.loop_angle()
    # create polar mesh
    angles = np.reshape(deltaL.tau * np.pi / 180, (1, -1))
    r = np.reshape(deltaL.x, (-1, 1))

    deltaL.x_polar = r * np.cos(angles)
    deltaL.y_polar = r * np.sin(angles)
    deltaL.z_polar = deltaL.z

    # select one heights 
    if z is not None:
        if z in deltaL.z:
            iz = np.where(z==deltaL.z)[0][0]
        else:
            logging.warning("z is not in deltaL.z")

        deltaL.z_polar = deltaL.z_polar[iz:iz+1]
        deltaL.deltaL_polar = deltaL.deltaL_polar[:,iz:iz+1,...]

    if stepx is not None:
        deltaL.deltaL_polar = deltaL.deltaL_polar[::stepx, ...]
        deltaL.x_polar = deltaL.x_polar[::stepx, ...]
        deltaL.x = deltaL.x[::stepx]
        deltaL.y_polar = deltaL.y_polar[::stepx, ...]

    deltaL.nx = deltaL.x_polar.shape[0]
    deltaL.ntau = deltaL.x_polar.shape[1]
    deltaL.nz = deltaL.z_polar.shape[0]

    # refine angles (linera interpolation)
    for kk in range(refine):
        deltaL.refine_angle()

    deltaL.nx = deltaL.x_polar.shape[0]
    deltaL.ntau = deltaL.x_polar.shape[1]
    deltaL.nz = deltaL.z_polar.shape[0]


    deltaL.shift_domain(simu.tx[iTurb]*simu.les.z_i,simu.ty[iTurb]*simu.les.z_i)
    deltaL.xS = simu.tx[iTurb]*simu.les.z_i
    deltaL.yS = simu.ty[iTurb]*simu.les.z_i
    return deltaL 


def concatenate_all_dl(
    casename: str,
    path2Pe: str,
    refine: int = 2,
    z: float = None,
    nx: int = None,
    ny: int = None,
    stepx: int = None,
    iTurb: np.ndarray = None,
    plot: bool = False,
    dl_fname: str = None,
    spp_fname: str = None,
    spl_fname: str = None
) -> None:
    """
    Concatenates all delta L fields for all turbines.

    This function processes the delta L fields for all turbines, optionally
    refines the angles, and saves the concatenated delta L fields into different files.

    Args:
        casename (str): The name of the simulation case.
        path2Pe (str): The path to the Pe results.
        refine (int, optional): The number of refinement steps. Defaults to 2.
        z (float, optional): The receiver height to select. Defaults to None.
        nx (int, optional): The number of points in the x-direction for cartesian interpolation. Defaults to None.
        ny (int, optional): The number of points in the y-direction for cartesian interpolation. Defaults to None.
        stepx (int, optional): The step size for a coarser grid in the radius direction. Defaults to None.
        iTurb (np.ndarray, optional): The indices of the turbines to process. Defaults to None.
        plot (bool, optional): Whether to plot the results. Defaults to False.
        dl_fname (str, optional): The filename to save the delta L field. Defaults to None.
        spp_fname (str, optional): The filename to save the Spp field. Defaults to None.
        spl_fname (str, optional): The filename to save the Spl field. Defaults to None.

    Returns:
        None
    """
    print('Start concatenating angles ...')
    if dl_fname is None :
        dl_fname = 'DL'
    if spp_fname is None:
        spp_fname = 'Spp'
    if spl_fname is None:
        spl_fname = 'Spl'

    # load simulation parameters
    simu = Simu(casename)
    simu.load(path2Pe+casename+'.dat')

    if iTurb is None:
        iTurb = np.arange(0,len(simu.tx))
    # if not trubine index given loop over all turbines 
    for ii in iTurb:
        deltaL = concatenate_angles(simu.casename,path2Pe,simu, ii, refine, z, stepx)
        if nx is None or ny is None:
            logging.info("saving delta L polar shape as :" + dl_fname+str(ii)+'.dat')
            print(deltaL.x_polar.shape)
            print(deltaL.z_polar.shape)

            deltaL.x_cart = np.meshgrid(deltaL.x_polar.flatten(),
                                         deltaL.z_polar,
                                         indexing='ij')[0].reshape(deltaL.nx,deltaL.ntau,deltaL.nz)
            deltaL.y_cart = np.meshgrid(deltaL.y_polar.flatten(),
                                         deltaL.z_polar,
                                         indexing='ij')[0].reshape(deltaL.nx,deltaL.ntau,deltaL.nz)
            deltaL.z_cart = np.meshgrid(deltaL.y_polar.flatten(),
                                         deltaL.z_polar,
                                         indexing='ij')[1].reshape(deltaL.nx, deltaL.ntau,deltaL.nz)
            deltaL.deltaL_cart = np.transpose(deltaL.deltaL_polar,(0,2,1,3,4))
            deltaL.x_array = deltaL.x_cart[:,0]
            deltaL.y_array = deltaL.y_cart[0,:]
            deltaL.z_array = deltaL.z_polar
            deltaL.x_polar = None
            deltaL.y_polar = None
            deltaL.deltaL_polar = None
            deltaL.deltaLlist = []
            deltaL.save(dl_fname+str(ii)+'.dat')
        else:
            # define common mesh  for interpolation from simu
            x = np.linspace(simu.x1,simu.x2,nx)
            y = np.linspace(simu.y1,simu.y2,ny)

            deltaL.interpolate_from_polar(x,y)

            logging.info("saving delta L cartesian shape as :" + dl_fname+str(ii)+'.dat')
            deltaL.x_polar = None
            deltaL.y_polar = None
            deltaL.deltaL_polar = None
            deltaL.deltaLlist = []
            deltaL.save(dl_fname+str(ii)+'.dat')

        # reset fields
        deltaL = None
        spl = None
        src = None
    print('done.')


def concatenate_side_dl(
    casename: str,
    path2Pe: str,
    iTurb: np.ndarray = None,
    tau: list = [0, 180],
    nx: int = None,
    nz: int = None,
    dl_fname: str = None,
    spp_fname: str = None,
    spl_fname: str = None
) -> None:
    """
    Concatenates side views of the delta L field for all turbines.

    This function processes the delta L fields to create side views and saves
    the results to files.

    Args:
        casename (str): The name of the simulation case.
        path2Pe (str): The path to the Pe results.
        iTurb (np.ndarray, optional): The indices of the turbines to process. Defaults to None.
        tau (list, optional): The tau values to process. Defaults to [0, 180].
        nx (int, optional): The number of points in the x-direction. Defaults to None.
        nz (int, optional): The number of points in the z-direction. Defaults to None.
        dl_fname (str, optional): The filename for the delta L field. Defaults to None.
        spp_fname (str, optional): The filename for the Spp field. Defaults to None.
        spl_fname (str, optional): The filename for the Spl field. Defaults to None.

    Returns:
        None
    """
    simu = Simu(casename)
    print('Start concatenating side view...')
    if dl_fname is None :
        dl_fname = './xz/DL'
    if spp_fname is None:
        spp_fname = './xz/Spp'
    if spl_fname is None:
        spl_fname = './xz/Spl'
    if iTurb is None:
        iTurb = np.arange(0,len(simu.tx))

    # load simulation parameters
    simu = Simu(casename)
    simu.load(path2Pe+casename+'.dat')
    for ii in iTurb:
        #--------------------------------------------------------------
        # create delta L field (read from Pe results)
        deltaL = DeltaLField(path2Pe, casename)
        deltaL.create_side_view(ii, tau)

        deltaL.shift_domain(simu.tx[ii]*simu.les.z_i,simu.ty[ii]*simu.les.z_i)

        deltaL.xS = simu.tx[ii]*simu.les.z_i
        deltaL.yS = simu.ty[ii]*simu.les.z_i

        # interpolate on coarse grid
        if (nx is not None) and (nz is not None):
            x = np.linspace(simu.x1, simu.x2, nx)
            z = np.linspace(0, simu.h, nz)
            deltaL.interpolate_xz(x, z)

        deltaL.save(dl_fname+str(ii)+'.dat')

        # reset fields
        deltaL = None
    print('done.')


def combine_dl_src(
    casename: str,
    path2Pe: str,
    iTurb: np.ndarray = None,
    dl_fname: str = None,
    spp_fname: str = None,
    spl_fname: str = None,
    polar: bool = False,
    third: bool = True,
    free_field: bool = False
) -> None:
    """
    Combines the delta L field with the source field for all turbines.

    This function processes the delta L and source fields, combines them, and
    saves the results to files.
    The field must be in cartesian coordinate system bu the matrix shape can be polar.
    If the size of the grid dont match the source results are interpolated on the delta L results grid.


    Args:
        casename (str): The name of the simulation case.
        path2Pe (str): The path to the Pe results.
        iTurb (np.ndarray, optional): The indices of the turbines to process. Defaults to None.
        dl_fname (str, optional): The filename for the delta L field. Defaults to None.
        spp_fname (str, optional): The filename for the Spp field. Defaults to None.
        spl_fname (str, optional): The filename for the Spl field. Defaults to None.
        polar (bool, optional): Whether to use polar coordinates. Defaults to False.
        third (bool, optional): Whether to compute SPL in third octave bands before saving (reduce the size of files). Defaults to True.
        free_field (bool, optional): Whether to consider free field conditions. Defaults to False.

    Returns:
        None
    """

    print('Start combining DeltaL and Spp...')
    if dl_fname is None :
        dl_fname = './xz/DL'
    if spp_fname is None:
        spp_fname = './xz/Spp'
    if spl_fname is None:
        spl_fname = './xz/Spl'

    # load simulation parameters 
    simu = Simu(casename)
    simu.load(path2Pe+casename+'.dat')

    # compute spl for gievn wind turbin
    if iTurb is None:
        iTurb = np.arange(0,len(simu.tx))

    for ii in iTurb:
        #-----------------------------------------------------------------
        # create delta L field (read from Pe results)
        deltaL = DeltaLField(path2Pe,casename)

        # loadDelta L field
        deltaL.load(dl_fname+str(ii)+'.dat')

        # read source file
        #-----------------------------------------------------------------
        src = Source()
        src.load(spp_fname+str(ii)+'.dat')

        # deltaL.deltaL_cart = np.transpose(deltaL.deltaL_polar, (0,2,1,3,4))

        # if deltaL.x_cart is None:
        #     print("interpolating Spp to polar ")
        #     if src.mesh.x_coord.size != deltaL.x_polar[:,:].size:
        #         src.interpolate_xy(deltaL.x_polar[:, : ], deltaL.y_polar[:, :])
        if src.mesh.y_array.size != deltaL.x_cart.shape[1]:
            src.interpolate_xy(deltaL.x_cart[:, :, 0], deltaL.y_cart[:, :, 0])
        elif src.mesh.z_array.size != deltaL.z_cart.shape[2]:
            src.interpolate_xz(deltaL.x_cart[:, 0, 0], deltaL.z_cart[0, 0, :])

        # create spl field
        #-----------------------------------------------------------------
        spl = SplField(src,deltaL)
        if polar:
            flag = spl.check_compatibility_polar()
        else:
            flag = spl.check_compatibility_cart()
        if flag == -1:
            return flag
        spl.src.wt.tau = 0

        # combine Delta L spp
        print(deltaL.height)
        spl.combine_linear_broadband(free_field=free_field)

        if third:
            fc = [50, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000]
            Nfc = [1,  1,  1,   1,   1,  1,   2,   2,   3,   4,   4,   4,   5,  5]
            spl.compute_third_octave(fc, Nfc)
        spl.atm_absorption()
        spl.save(spl_fname+str(ii)+'.dat')

        # reset fields
        deltaL = None
        spl = None
        src = None


def convert_to_receiver_time(
    casename: str,
    path2Pe: str,
    iTurb: np.ndarray = None,
    spl_fname: str = None,
    oaspl: bool = False,
    spl_fname_out: str = None,
    **kwargs
) -> None:
    """
    Converts the SPL field to receiver time for all turbines.

    Args:
        casename (str): The name of the simulation case.
        path2Pe (str): The path to the Pe results.
        iTurb (np.ndarray, optional): The indices of the turbines to process. Defaults to None.
        spl_fname (str, optional): The filename for the SPL field. Defaults to None.
        oaspl (bool, optional): Whether to load OASPL data. Defaults to False.
        spl_fname_out (str, optional): The output filename for the SPL field. Defaults to None.
        **kwargs: Additional keyword arguments.

    Returns:
        None
    """
    if spl_fname is None:
        spl_fname = 'Spl'

    # load simulation parameters
    simu = Simu(casename)
    simu.load(path2Pe+casename+'.dat')
    if iTurb is None:
        iTurb = np.arange(0, len(simu.tx))

    for ii in iTurb:
        # -----------------------------------------------------------------
        # create spl field
        # -----------------------------------------------------------------
        spl = SplField()
        if oaspl:
            spl.load(spl_fname+str(ii)+'.dat',seg=False,time=False,oaspl=True)
        else:
            spl.load(spl_fname+str(ii)+'.dat',seg=True,time=False,oaspl=False)

        spl.info()

        # Receiver Time computation
        # -----------------------------------------------------------------
        # compute Spl at receiver time
        spl.create_full_rotation()
        # compute source/receiver time matrix
        spl.compute_real_receiver_time(last=True, loop=True)
        # create time solution from angle
        spl.angle_to_time_3_full(dt=0.1)
        # save spl  in file
        spl.clean_full_rotation()
        if spl_fname_out is not None:
            spl.SPL_seg = None
            spl.OASPL_seg = None
            spl.save(spl_fname_out+str(ii)+'.dat')
        else:
            spl.save(spl_fname+str(ii)+'.dat')
        # reset fields
        spl = None


def concatenate_planes_dl(
    casename: str,
    path2Pe: str,
    nx: int = 350,
    ny: int = 320,
    iTurb: int = None,
    xplane: np.ndarray = None,
    yplane: np.ndarray = None,
    dl_fname: str = None
) -> None:
    """
    This function processes the delta L field to concatenate planes and saves
    the results to files.
    WARNING: this was not used a lot and it may be a bit broken. 

    Args:
        casename (str): The name of the simulation case.
        path2Pe (str): The path to the Pe results.
        nx (int, optional): The number of points in the x-direction. Defaults to 350.
        ny (int, optional): The number of points in the y-direction. Defaults to 320.
        iTurb (int, optional): The index of the turbine to process. Defaults to None.
        xplane (np.ndarray, optional): The x-planes to process. Defaults to None.
        yplane (np.ndarray, optional): The y-planes to process. Defaults to None.
        dl_fname (str, optional): The filename for the delta L field. Defaults to None.

    Returns:
        None
    """
    # load simulation parameters 
    simu = Simu(casename) 
    simu.load(path2Pe+casename+'.dat')

    # compute spl for given wind turbine
    if iTurb is not None:
        #-----------------------------------------------------------------
        # create delta L field (read from Pe results)
        
        # retrieve planes to be computed 
        if xplane is not None:
            deltaL = DeltaLField(path2Pe,casename)
            # check that planes was set to be recorded
            if not np.all(np.isin(xplane,simu.xplanes)) :
                print('planes were not calculated')
                return -1
            
            # set plane in the turbine coordinate system      
            xplane = xplane-simu.tx[iTurb]*simu.les.z_i
            # read all slice for all heights and angles  
            deltaL.read_all_slices(iTurb,simu.heights,simu.tau,simu.distribute_tau,xplanes=xplane)
            
            # concatenate slice x constant
            deltaL.concatenate_xslices()
            # # shift domain to absolute coordinate
            deltaL.shift_domain(simu.tx[iTurb]*simu.les.z_i,simu.ty[iTurb]*simu.les.z_i)
            # # create interpolation mesh 
            x = np.linspace(simu.x1,simu.x2,nx)
            y = np.linspace(simu.y1,simu.y2,ny)
            # set sour receiver (usefull for later)
            deltaL.xS = simu.tx[iTurb]*simu.les.z_i
            deltaL.yS = simu.ty[iTurb]*simu.les.z_i
            # interpolate mesh on regular y values 
            print('interpolating ...')
            # deltaL.interpolate_planes(y=y)
            deltaL.interpolate_planes(y=y)

            # save mesh
            if dl_fname is None :
                dl_fname = './yz/DL'
            print('saving ...')
            deltaL.save(dl_fname+str(iTurb)+'.dat')
            deltaL = None
            quit()

        # retrieve planes to be computed 
        if yplane is not None:
            deltaL = DeltaLField(path2Pe,casename)
            # check that planes were recorded
            if not np.all(np.isin(yplane,simu.yplanes)) :
                print('planes were not calculated')

            yplane = yplane-simu.ty[iTurb]*simu.les.z_i
            # read planes from all angles all receiver
            deltaL.read_all_slices(iTurb,simu.heights,simu.tau,simu.distribute_tau,yplanes=yplane)

            # concatenate slice x constant
            deltaL.concatenate_yslices()
            # shift domain to absolute coordinate
            deltaL.shift_domain(simu.tx[iTurb]*simu.les.z_i,simu.ty[iTurb]*simu.les.z_i)
            # create interpolation mesh 
            x = np.linspace(simu.x1,simu.x2,nx)
            y = np.linspace(simu.y1,simu.y2,ny)
            # set sour receiver (usefull for later)
            deltaL.xS = simu.tx[iTurb]*simu.les.z_i
            deltaL.yS = simu.ty[iTurb]*simu.les.z_i
            # interpolate mesh on regular y values 
            print('interpolating ...')
            deltaL.interpolate_planes(x=x)

            # save mesh
            if dl_fname is None :
                dl_fname = './xz/DL'
            print('saving ...')
            deltaL.save(dl_fname+str(iTurb)+'.dat')
            deltaL = None
    # if not trubine index given loop over all turbines 
    else:
        for iTurb in range(len(simu.tx)):
            #-----------------------------------------------------------------
            # create delta L field (read from Pe results)
            deltaL = DeltaLField(path2Pe,casename)

            # retrieve planes to be computed 
            if xplane is not None:
                # check that planes were recorded
                if not np.all(np.isin(xplane,simu.xplanes)) :
                    print('planes were not calculated')
                else:
                    xplane = xplane-simu.tx[iTurb]*simu.les.z_i
                # read planes from all angles all receiver
                for plane in xplane:
                    for t in simu.tau:
                        for h in simu.heights:
                            deltaL.read_planes(iTurb,h,t,simu.distribute_tau,xplane=plane)
                    deltaL.concatenate_plane()
                    deltaL.plane_list=[]
                # shift domain to absolute coordinate
                deltaL.shiftDomain(simu.tx[iTurb]*simu.les.z_i,simu.ty[iTurb]*simu.les.z_i)
                # create interpolation mesh 
                x = np.linspace(simu.x1,simu.x2,nx)
                y = np.linspace(simu.y1,simu.y2,ny)
                # set sour receiver (usefull for later)
                deltaL.xS = simu.tx[iTurb]*simu.les.z_i
                deltaL.yS = simu.ty[iTurb]*simu.les.z_i
                # interpolate mesh on regular y values 
                deltaL.interpolate_planes(y=y)
                # save mesh
                if dl_fname is None :
                    dl_fname = 'DLx'
                deltaL.save(dl_fname+str(iTurb)+'.dat')
                deltaL = None
            
            deltaL = DeltaLField(path2Pe,casename)
            # retrieve planes to be computed 
            if yplane is not None:
                # check that planes were recorded
                if not np.all(np.isin(yplane,simu.yplanes)) :
                    print('planes were not calculated')
                else:
                    yplane = yplane-simu.ty[iTurb]*simu.les.z_i
                # read planes from all angles all receiver
                for plane in yplane:
                    for t in simu.tau:
                        for h in simu.heights:
                            deltaL.read_planes(iTurb,h,t,simu.distribute_tau,yplane=plane)
                    deltaL.concatenate_plane()
                    deltaL.plane_list=[]
                # shift domain to absolute coordinate
                deltaL.shiftDomain(simu.tx[iTurb]*simu.les.z_i,simu.ty[iTurb]*simu.les.z_i)
                # create interpolation mesh 
                x = np.linspace(simu.x1,simu.x2,nx)
                y = np.linspace(simu.y1,simu.y2,ny)
                # set sour receiver (usefull for later)
                deltaL.xS = simu.tx[iTurb]*simu.les.z_i
                deltaL.yS = simu.ty[iTurb]*simu.les.z_i
                # interpolate mesh on regular y values 
                deltaL.interpolate_planes(x=x)
                # save mesh
                if dl_fname is None :
                    dl_fname = 'DLy'
                deltaL.save(dl_fname+str(iTurb)+'.dat')
                deltaL = None


def postprocessingPlanesCombine(
    casename: str,
    path2Pe: str,
    iTurb: int = None,
    dl_fname: str = None,
    spp_fname: str = None,
    spl_fname: str = None
) -> None:
    """
    Combines the delta L field with the source field for a given turbine and plane.

    Args:
        casename (str): The name of the simulation case.
        path2Pe (str): The path to the Pe results.
        iTurb (int, optional): The index of the turbine to process. Defaults to None.
        dl_fname (str, optional): The filename for the delta L field. Defaults to None.
        spp_fname (str, optional): The filename for the Spp field. Defaults to None.
        spl_fname (str, optional): The filename for the Spl field. Defaults to None.

    Returns:
        None
    """
    if dl_fname is None :
        dl_fname = 'DLy'
    if spp_fname is None:
        spp_fname = 'SppY_'
    if spl_fname is None:
        spl_fname = 'SplY_'
    
    # load simulation parameters 
    simu = Simu(casename) 
    simu.load(path2Pe+casename+'.dat')
    
    #-----------------------------------------------------------------
    # create delta L field (read from Pe results)
    deltaL = DeltaLField(path2Pe,casename)

    # loadDelta L field 
    deltaL.load(dl_fname+str(iTurb)+'.dat')

    # read source file 
    #-----------------------------------------------------------------
    src = Source()
    src.load(spp_fname+str(iTurb)+'.dat')

    print(src.Spp.shape)

    print(deltaL.deltaLInterpolated.shape)


    # create spl field 
    #-----------------------------------------------------------------
    z = np.linspace(0,300,100)
    x = np.linspace(500,4000,100)


    # src.interpolateXY(x,np.array([1600]))
    src.interpolateXZ(x,z)
    deltaL.interpolateXZ(x,z)

    print(src.y_interpolate.shape)
    print(deltaL.y_interpolate.shape)

    print(src.y_interpolate[0,:,0])
    print(deltaL.y_interpolate[0,:,0])
    spl = SplField(src,deltaL)

    spl.check_compatibility()
    spl.src.wt.tau = 0

    # combine Delta L spp 
    spl.combine_linear_broadband()
    # spl.src.save(spp_fname+str(iTurb)+'.dat')
    # spl.deltaL.save(dl_fname+str(iTurb)+'.dat')
    spl.save(spl_fname+str(iTurb)+'.dat')
