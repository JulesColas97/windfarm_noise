!*******************************************************************************
subroutine initialize
!*******************************************************************************
!
! Driver routine for initialization
!
use sim_param, only: u_dim 
use param, only: z_i,rp,coord,dt,jt_total,path,clock,mean_p_force,use_mean_p_force,grid_stretch, grid_size 
use param, only: DYN_init,sgs_model,nz_tot,nx,nx2,ny,ny2,nz,nproc,use_cfl_dt
use param, only: dx,dy,dz,tadv1,tadv2,tadv3,tfac1,tfac2,L_x,L_y,L_z,ld,lh,ld_big,lh_big
use param, only: idx,idy,idz,inxny,inx2ny2,u_star,vonk
#ifdef SCALAR 
use param, only: Kiw, Kpw, dmpfr, g
#endif SCALAR 
use sgs_param, only: sgs_param_init
use mpi, only: mpi_wtime
use sim_param,only: sim_param_init, u_dim 
use hdf5
#ifdef CORIOLIS
use param, only: coriol,coriol_x,coriol_y,z_i,adjust_wind_angle
use param, only: ug, vg, Kiw, Kpw
use sim_param, only: vgMag
#ifdef CPS
use param, only: fringe_factor
#endif CPS
#endif CORIOLIS
#ifdef BAROCLINIC
use param, only: ug_delta, vg_delta, bar_start, bar_end, geo_force_x, geo_force_y
use grid_defs, only: gridz
#endif BAROCLINIC
#ifdef SCALAR
use scalars_param 
#endif SCALAR
#ifdef LVLSET
use level_set_mod, only: zwall
#endif LVLSET
#ifdef ATM
use atm_lesgo_interface, only: atm_lesgo_initialize
#endif ATM

implicit none
integer:: hdf_error
character(len=100) :: make_output_dir
character(len=100) :: make_output_dir2
#ifdef TURBINES
character(len=100) :: make_turbine_dir
#endif TURBINES
integer :: ival_read

make_output_dir = 'mkdir -p ' // trim(path) // 'output'
make_output_dir2 = 'mkdir -p ' // trim(path) // 'continuation'
#ifdef TURBINES
make_turbine_dir = 'mkdir -p ' // trim(path) // 'turbine'
#endif TURBINES

! Initialize mpi
call initialize_mpi

if( coord == 0 ) call system( make_output_dir )
if( coord == 0 ) call system( make_output_dir2 )
#ifdef TURBINES
if( coord == 0 ) call system( make_turbine_dir )
#endif TURBINES

! Read input file
call read_input_conf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Set parameters based on read input values
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
vonk = 0.4_rp   ! von Karman constant

! Set the processor owned vertical grid spacing
nz = floor ( real( nz_tot, rp ) / nproc ) + 1
if(nproc.eq.1) nz=nz-1

! Recompute nz_tot to be compliant with computed nz
ival_read = nz_tot
nz_tot = ( nz - 1 ) * nproc + 1
if(coord==0 .AND. ival_read /= nz_tot ) write(*,*)  'Reseting Nz (total) to: ',nz_tot

! Grid size for dealiasing
nx2 = 3 * nx / 2
ny2 = 3 * ny / 2
! Grid size for FFT's
lh = nx / 2 + 1
ld = 2 * lh
lh_big = nx2 / 2 + 1
ld_big = 2 * lh_big

! Grid spacing (x direction)
dx = L_x / nx

! Grid spacing (y direction)
dy = L_y / ny

! Grid spacing (z direction) 
if(grid_stretch) then 
  !Set grid size in the uniform portion of the domain 
  dz = (grid_size)/z_i 
else 
  !Calculate dz value 
  dz = L_z / ( nz_tot - 1 )
endif 

idx    =1.0_rp/dx
idy    =1.0_rp/dy
idz    =1.0_rp/dz
inxny  =1.0_rp/(nx*ny)
inx2ny2=1.0_rp/(nx2*ny2)

! Constants time integration
if( .not. use_cfl_dt ) then
! Set AB3 integration coefficients
tadv1 = 23.0_rp / 12.0_rp
tadv2 = -16.0_rp / 12.0_rp
tadv3 = 5.0_rp / 12.0_rp
tfac1 = -1.0_rp
tfac2 = 2.0_rp
endif

#ifdef CORIOLIS
u_star = 1._rp
#else 
u_dim = 1.0_rp
#endif CORIOLIS

if(use_mean_p_force) then
  ! Evaluate the mean pressure force
  ! length normalized by z_i
#ifdef LVLSET
  mean_p_force = u_star**2 / (L_z - zwall)
  if (zwall.ne.0.0_rp) then
  ! Note no coord.eq.0 to make sure it is printed many times
  write(*,*) 'zwall not equal to zero. Please double check this case'
  write(*,*) 'Used for Liu and Stevens, Computers and Fluids 2020 paper'
  ! https://doi.org/10.1016/j.compfluid.2020.104604
  endif
#else
  mean_p_force = u_star**2 / L_z
#endif LVLSET
else
  mean_p_force = 0.0_rp
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

clock(1) = mpi_wtime()

! Start HDF5
call h5open_f(hdf_error)

! Define simulation parameters
call sim_param_init

#ifdef CORIOLIS
u_dim = sqrt(ug**2 + vg**2)  
! Srinidhi : Feb 20, 2022 
! Somehow the code likes small values (normalized values around 1) for Coriolis cases 
! Therefore the velocities for Coriolis cases are being normalized w.r.t. geostrophic wind velocity. 
! Non-dimensionalize Coriolis parameter with z_i and u=1 (m/s)
coriol   = coriol   *  (z_i/u_dim)
coriol_x = coriol_x *  (z_i/u_dim)
coriol_y = coriol_y *  (z_i/u_dim)
ug = ug/u_dim; vg=vg/u_dim;
vgMag = 1._rp;
#ifdef SCALAR
dmpfr = dmpfr / u_dim
#endif SCALAR
Kiw = Kiw*z_i**2/u_dim**2
Kpw = Kpw*z_i/u_dim

!#ifdef CPS
!fringe_factor = fringe_factor /u_dim
!#endif CPS
#else 
u_dim = 1.0_rp
#endif CORIOLIS

#ifdef BAROCLINIC
ug_delta = ug_delta/u_dim; vg_delta = vg_delta/u_dim
bar_start = bar_start/z_i; bar_end = bar_end/z_i

allocate( geo_force_x(1:nz) ); geo_force_x = 0.0_rp
allocate( geo_force_y(1:nz) ); geo_force_y = 0.0_rp
#endif BAROCLINIC

#ifdef SCALAR
! Non-dimensionalize gravity
g = g*z_i/u_dim**2
#endif SCALAR

! Initialize sgs variables
call sgs_param_init

! Initialize uv grid (calculate x,y,z vectors)
call grid_build

! Initialize all the scalar terms 
#ifdef SCALAR
call scalars_init 
#endif SCALAR

!Initialize wind angle controller
#ifdef CORIOLIS
if(adjust_wind_angle) call wind_angle_controller_init
#endif

! Initialize turbines
#ifdef TURBINES
call turbine_init  !must occur before initial is called
#endif TURBINES

! Initialize windbreaks
#ifdef WINDBREAKS
call windbreak_init
#endif WINDBREAKS

#ifdef ATM
call atm_lesgo_initialize
#endif ATM

! Formulate the fft plans--may want to use FFTW_USE_WISDOM
! Initialize the kx,ky arrays
call init_fft

! Initialize test filter(s)
! this is used for lower BC, even if no dynamic model
call test_filter_init

#ifdef LVLSET
! Initialize level set method
call level_set_init
#endif LVLSET

! Initialize velocity field
call initial 

! Initialize concurrent precursor stuff
#ifdef CPS
call initialize_cps
#endif CPS

! Initialize pressure solver (should be after init_fft)
call press_stag_init

if(sgs_model.eq.1 .or. jt_total<DYN_init) call smagorinsky_model_init 
if(sgs_model.eq.5 .and. coord.eq.0) print *,'running dynamic sgs_model = ',sgs_model

! Print out warnings

#ifdef LVLSET
#ifdef SCALAR
! update surface temperature and heat flux
if(hbc.eq.0)
! No coord.eq.0. to make sure it is printed many times
write(*,*) 'Check note 3 in the “to watch list in the documentation”'
endif

if(hbc.eq.1 .and. gabls_test) then
! No coord.eq.0. to make sure it is printed many times
write(*,*) 'Check note 3 in the “to watch list in the documentation”'
endif

#endif SCALAR
#endif LVLSET

end subroutine initialize
