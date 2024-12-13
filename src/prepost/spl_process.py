from .spl import SplField
from .post import Simu, DeltaLField
from .post import PeResults
import numpy as np

import matplotlib.pyplot as plt
import logging
from .source.main import Source



def concatenate_angles_dl(casename,path2Pe,refine=2,nx=350,ny=320,iTurb=None,plot=False,
                        dl_fname=None,spp_fname=None,spl_fname=None):
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



def init_deltaL_field(path2Pe,casename,simu,iTurb):
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


def add_1Pe_to_deltaL_field(deltaL, iTurb, height, tau):
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

def concatenate_angles(casename, path2Pe, simu, iTurb, refine, z, stepx):
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
    print(deltaL.x_polar.shape)
    print(deltaL.y_polar.shape)
    print(deltaL.nx)
    print(deltaL.nz)

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


def concatenate_all_dl(casename,path2Pe,refine=2,z=None,nx=None,ny=None,stepx=None,iTurb=None,plot=False,
                        dl_fname=None,spp_fname=None,spl_fname=None):
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
            # define common mesh 
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


def concatenate_side_dl(casename,path2Pe, iTurb=None, tau=[0,180], nx=None, nz=None,
                        dl_fname=None, spp_fname=None, spl_fname=None):

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


def combine_dl_src(casename,path2Pe,iTurb=None,
                            dl_fname=None,spp_fname=None,spl_fname=None,
                   polar=False,third=True,free_field=False):

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
        # if polar:
        #     deltaL.x_cart = deltaL.x_polar[::10,:]
        #     deltaL.y_cart = deltaL.y_polar[::10,:]
        #     deltaL.deltaL_polar = deltaL.deltaL_polar[::10,:,...]

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
        # spl.src.save(spp_fname+str(ii)+'.dat')
        # spl.deltaL.save(dl_fname+str(ii)+'.dat')

        if third:
            fc = [50, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000]
            Nfc = [1,  1,  1,   1,   1,  1,   2,   2,   3,   4,   4,   4,   5,  5]
            spl.compute_third_octave(fc, Nfc)
        spl.atm_absorption()
        # spl.Aweight()

        spl.save(spl_fname+str(ii)+'.dat')

        # remove interpolated field from source
        # src.SppInterpolated = None
        # src.x_interpolate = None
        # src.y_interpolate = None
        # src.z_interpolate = None
        # src.save(spp_fname+str(ii)+'.dat')
        # reset fields
        deltaL = None
        spl = None
        src = None


def convert_to_receiver_time(casename, path2Pe, iTurb=None, spl_fname=None,
                             oaspl=False, spl_fname_out=None, **kwargs):
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
        # print(spl.SPL_seg.shape)
        # print(spl.OASPL_seg.shape)

        # Receiver Time computation
        # -----------------------------------------------------------------
        # compute Spl at receiver time
        spl.create_full_rotation()
        # compute source/receiver time matrix
        spl.compute_real_receiver_time(last=True, loop=True)
        # create time solution from angle
        spl.angle_to_time_3_full(dt=0.1)
        # save spl  in file
        # spl.save(spl_fname+str(ii)+'.dat')
        spl.clean_full_rotation()
        if spl_fname_out is not None:
            spl.SPL_seg = None
            spl.OASPL_seg = None
            spl.save(spl_fname_out+str(ii)+'.dat')
        else:
            spl.save(spl_fname+str(ii)+'.dat')
        # reset fields
        spl = None


def concatenate_planes_dl(casename, path2Pe,nx=350,ny=320,iTurb=None,xplane=None,yplane=None,dl_fname=None):
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


def postprocessingPlanesCombine(casename,path2Pe,iTurb=None,
                            dl_fname=None,spp_fname=None,spl_fname=None):
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
