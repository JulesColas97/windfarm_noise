import prepost as pp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FFMpegWriter
import logging

# logging.basicConfig(level=logging.DEBUG, force=True)
logging.basicConfig(format='%(levelname)s: %(message)s', force=True)

# define which post processing are run 
# -----------------------------------------------------------------------------
check = False
read_dl = False
combine = True

interpolate = False
oaspl = False
time = False
combine_turbines = False

# -----------------------------------------------------------------------------
workdir = './'
case = 'c0'
path2Pe = workdir
iTurb = [0, 1]

simu = pp.Simu(workdir + case)
simu.load(workdir + 'c0.dat')
print(simu.frequencies.shape)

fc = [50, 63, 80, 100]
Nfc = [1,  1,  1,   1]

# post process from PE h5 files
# -----------------------------------------------------------------------------
if check:
    print('checking delta L ...')
    print('------------------------------------------------------------------')
    simu.check_run_cases()

# post process from PE h5 files
# -----------------------------------------------------------------------------
if read_dl:
    print('reading delta L ...')
    print('------------------------------------------------------------------')
    pp.concatenate_all_dl(case, path2Pe, refine=4, iTurb=iTurb,
                          z=2, stepx=5,
                          dl_fname=workdir+'/xy/DL_polar',
                          spl_fname=workdir+'/xy/spl_polar',
                          spp_fname=workdir+'/xy/spp')


# combine
# -----------------------------------------------------------------------------
if combine:
    print('combining ...')
    print('------------------------------------------------------------------')
    pp.combine_dl_src(case, path2Pe, iTurb=iTurb,
                      dl_fname=workdir+'/xy/DL_polar',
                      spp_fname=workdir+'/xy/spp_polar',
                      spl_fname=workdir+'/xy/spl_s_polar_',
                      polar=False,third=False)

# interpolate sourtce
if interpolate:
    for ii in iTurb:
        spl = pp.SplField()
        spl.load(workdir + '/xy/spl_s_polar_' + str(ii) + '.dat', seg=True,
                 time=False, oaspl=False, am=False, mean=False, z=2)

        x = np.linspace(simu.x1, simu.x2, 350)
        y = np.linspace(simu.y1, simu.y2, 220)
        spl.interpolate_from_polar(x, y)
        spl.save(workdir + '/xy/spl_s_cart_'+str(ii) + '.dat')

# OASPL
# -----------------------------------------------------------------------------
if oaspl:
    for ii in iTurb:
        spl = pp.SplField()
        spl.load(workdir+'/xy/spl_s_cart_%s.dat' % ii, seg=True, time=False,
                 am=False, oaspl=False, mean=False)
        spl.compute_third_octave(fc, Nfc)
        spl.atm_absorption()
        spl.Aweight()
        spl.compute_oaspl()
        # spl.compute_am()
        # spl.compute_time_average()
        spl.SPL_seg = None
        spl.save(workdir+'/xy/oaspl_dbA_s_cart_%s.dat' % ii)
        spl = None


# time
# -----------------------------------------------------------------------------
if time:
    print('compute time reciver ...')
    print('------------------------------------------------------------------')
    pp.convert_to_receiver_time(
        case, path2Pe, iTurb=iTurb, spl_fname=workdir+'/xy/spl_s_cart_',
        spl_fname_out=workdir+'/xy/spl_t_cart_', oaspl=False)


if combine_turbines:
    spl0 = pp.SplField()
    spl0.load(workdir + '/xy/spl_t_cart_0.dat', seg=True,
              time=True, am=False, mean=False, z=2)
    print(spl0.SPL_time.shape)
    for ii in iTurb[1:]:
        spl1 = pp.SplField()
        spl1.load(workdir + "/xy/spl_t_cart_%s.dat"%ii, seg=True,
                  time=True, am=False, mean=False, z=2)
        spl0.combine_2_turbines(spl1, tmax=250)
        spl1 = None

    spl0.save(workdir + '/xy/oaspl_t_cart_sum.dat')
