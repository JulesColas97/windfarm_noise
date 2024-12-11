module param
use mpi, only: mpi_STATUS_SIZE,MPI_COMM_WORLD
implicit none
integer, parameter :: sp = kind(1.0)
integer, parameter :: rp = selected_real_kind(2*precision(1.0_sp))
integer, parameter :: qp = selected_real_kind(2*precision(1.0_rp))
integer, parameter :: lprec=rp
!---------------------------------------------------
! GLOBAL PARAMETERS
!---------------------------------------------------
character(len=32) :: PATH = './'
real(rp),dimension(2) :: clock,clock_total,clock_loc

!---------------------------------------------------
! MPI PARAMETERS
!---------------------------------------------------
integer :: status(mpi_STATUS_SIZE)

integer :: nproc = 1 !--this must be 1 if no MPI
integer :: rank = 0   !--init to 0 (so its defined, even if no MPI)
integer :: coord = 0 

integer :: ierr,up,down,global_rank
integer :: mpi_rp,mpi_cp
integer, allocatable, dimension(:) ::  rank_of_coord
integer :: kmin,kmin2,kmax,kmax2  ! levels that "belong" to this processor
integer :: localComm

!---------------------------------------------------
! COMPUTATIONAL DOMAIN PARAMETERS
!---------------------------------------------------
real(rp),parameter::pi=3.1415926535897932384626433_rp

integer :: Nx=64, Ny=64, Nz=64
integer :: nz_tot = 64
integer :: nx2, ny2
integer :: lh, ld, lh_big, ld_big

! this value is dimensional [m]:
real(rp) :: z_i = 1000.0_rp

! these values should be non-dimensionalized by z_i: 
! set as multiple of BL height (z_i) then non-dimensionalized by z_i
real(rp) :: L_x = 2.0*pi, L_y=2.0*pi, L_z=1.0_rp

!Switch on grid_stretch 
logical :: grid_stretch = .false.
real(rp) :: grid_size = 10.0_rp 

! these values are also non-dimensionalized by z_i:
real(rp) :: dx,dy,dz

! Normalization factors
real(rp) :: idx,idy,idz,inxny,inx2ny2

!---------------------------------------------------
! MODEL PARAMETERS
!---------------------------------------------------
! Model type: 1->Smagorinsky; 5-> Lagragian scale-dep
integer :: sgs_model=5, wall_damp_exp=2
real(rp) :: clip=0.2

! timesteps between dynamic Cs updates
integer :: cs_count = 9

! When to start dynamic Cs calculations
integer :: DYN_init = 1000
  
! Cs is the Smagorinsky Constant
! Co and wall_damp_exp are used in the mason model for Smagorinsky coeff
real(rp) :: Co = 0.16_rp
  
! von Karman constant
real(rp) :: vonk = 0.4_rp
  
#ifdef CORIOLIS
!---------------------------------------------------
! CORIOLIS  PARAMETERS
!---------------------------------------------------
!Turbine hub height 
real(rp) :: z_wind = 100.0_rp !(in m) this is the hub-height at which the wind angle is controlled
real(rp) :: u_hub = 1.0_rp, v_hub=1.0_rp

logical :: neutral_ekman_test = .false.
!Non-dimensional coriolis parameter, geostrophic velocities, baroclinic parameters
logical :: coriolis_forcing = .true.
real(rp) :: coriol = 1.0e-4,coriol_x = 0.0,coriol_y = 0.0, ug=1.0_rp, vg=0.0_rp
real(rp) :: wind_angle = 0.0_rp

!PID controller for wind angle controller
logical :: adjust_wind_angle = .false. !Activates wind angle controller
real(rp) :: Kpw = 1.0_rp !Proportional constant of PID controller
real(rp) :: Kiw = 1.0_rp !Integral constant of PID controller
real(rp) :: Kdw = 0.0_rp !Derivative constant of PID controller
integer  :: height_index = 1
integer  :: rank_hub
logical  :: fixed_height = .false.
#endif

#ifdef BAROCLINIC
!---------------------------------------------------
! BAROCLINIC  PARAMETERS
!---------------------------------------------------
real(rp) :: ug_delta = 0.0_rp, vg_delta = 0.0_rp
real(rp) :: bar_start = 0.0_rp, bar_end = 100.0_rp
logical :: temp_adv = .false.
real(rp), dimension(:), allocatable :: geo_force_x, geo_force_y
#endif

#ifdef SCALAR
!---------------------------------------------------
! SCALAR  PARAMETERS
!---------------------------------------------------
!Prandtl number : This is only used with Smagorinsky model
real(rp) :: Pr = 0.3_rp 

!Roughness for the heat (non-dimensionalized)
real(rp) :: zo_scalar = 1.e-4
 
!Gabls_test: Stable BL, GABLS1 test case, see Beare et al. (2006) 
!Sullivan_test: Unstable BL test case, Sullivan et al. (2011) or Abkar and Moin (2017)
logical :: gabls_test = .false., sullivan_test = .false., neutral_test = .false., manual_test = .false.

!Parameters for manual case
real(rp) :: cap_height=1000.0_rp, cap_strength=0.0_rp, strat_height=1200.0_rp, strat_strength=0.0_rp
integer :: surface_cond = 0
real(rp) :: surface_rate=0.0_rp, surface_flux=0.0_rp

logical :: initsc = .false. 
!g: acceleration due to gravity (m/s2), inv_strength: Capping inversion strength K/Km
real(rp) :: g = 9.81_rp, inv_strength = 0.012_rp 
!T_scale: Temperature for non-dimensionalizing, wt_s: Surface heat flux, T_init:
!Initial temperature
real(rp) :: cooling_rate = 0.25 !K/h, cooling rate used in GABLS test
real(rp) :: theta_s1 = 300._rp,T_scale = 300.0_rp, wt_s = 0.2_rp, T_init = 283._rp
real(rp) :: z01 = 0.1_rp !Roughness height for temperature in (m)

! heat/thermal boundary condition specified with a constant surface
! hbc=0: potential temperature; hbc=1: heat flux
integer :: hbc=1

!ubc=0 : No damping layer at the top; ubc=1 : Activates damping layer
!damping_method=1 => Rayleigh damping layer
!damping_method=2 => Klemp and Lilly damping layer (1978)
!dmpfr: Damping factor in non-dimensional seconds
!ra_damp_factor: Rayleigh damping factor generally fixed at 2.0
integer :: ubc=0, damping_method=1, dmpfr = 1.0, ra_damp_exponent = 2.0
real(rp):: damp_hgt = 1.0 !(default) damping layer starts at the 3/4th of the L_z domain 

!PID controller to limit the BL growth
logical :: damping_x = .false. 
logical :: temp_sink = .false. !.true. activates the PID controller to limit the BL growth
!The layer in which the PID controller is activated, BL does not grow beyond z_heat_min
real(rp) :: z_heat_min, z_heat_max
!Time at which the PID controller should be activated (non-dimensional time)
real(rp) :: time_heat = 10._rp  
!To tune the controller just set Kd and Ki to zero and set Cp to very high value
!(say 10 or 100), if the controller starts doing what it's supposed to be doing
!then reduce the P gain and increase the I and D gains for optimum value
real(rp) :: Cp = 10.2_rp  !Proportional constant for the PID controller
real(rp) :: Ci = 0.5_rp   !Integral constant for the PID controller
real(rp) :: Cd = 0.0_rp   !Derivative constant of the PID controller
#endif 

!---------------------------------------------------
! TIMESTEP PARAMETERS
!---------------------------------------------------
integer :: nsteps = 50000
! -- Maximum runtime in seconds. Simulation will exit if exceeded. (disabled by default)
integer :: runtime = -1 

logical :: use_cfl_dt = .false.
logical :: cont_reduce = .false.
real(rp) :: cfl = 0.08
real(rp) :: dt_f=2.0e-4
real(rp) :: dt_ff=2.0e-4
real(rp) :: dt = 2.0e-4
  
! time advance parameters (Adams-Bashforth, 3rd order accurate)
real(rp) :: tadv1,tadv2,tadv3,tfac1, tfac2
 
integer :: jt_total=0 ! Global time counter

! Used to track the computational time in the code
real(rp) :: total_time,time_sgs,time_convec,time_press
real(rp) :: time_project,time_deriv,time_divstress
real(rp) :: time_wall,time_output,time_pdf,time_spec,time_RHS,time_cfl
#ifdef CORIOLIS
real(rp) :: time_coriolis
#endif 
#ifdef SCALAR 
real(rp) :: time_scalar,time_scalar_deriv,time_temperature
#endif 
#ifdef TURBINES
real(rp) :: time_turbine
real(rp) :: R1_ind ! Precalculated factor indicator function
#endif
#ifdef LVLSET
real(rp) :: time_level_set,time_level_set2
#endif  
#ifdef ATM
real(rp) :: time_atm
#endif ATM
#ifdef WINDBREAKS
real(rp) :: u_d_filter
#endif WINDBREAKS

!---------------------------------------------------
! BOUNDARY/INITIAL CONDITION PARAMETERS
!---------------------------------------------------  

! initu = true to read from a file; false to create with random noise
logical :: initu = .false.
! only used if initu=.false.
logical :: initu_interp = .false.

! lbc: lower boundary condition:  0 - stress free, 1 - wall 
! NOTE: the upper boundary condition is stress free
integer :: lbc_mom = 1
  
! lower boundary condition, roughness length
real(rp) :: zo = 0.0001_rp ! nondimensional

! prescribed inflow:
logical :: inflow = .false.
logical :: symmetric_fringe = .false.
real(rp) :: fringe_factor = 30.0 
! if inflow is true the following should be set:
! position of right end of fringe region, as a fraction of L_x
real(rp) :: fringe_region_end  = 1.0_rp
! length of fringe region as a fraction of L_x
real(rp) :: fringe_region_len = 0.125_rp
! position of right end of fringe region, as a fraction of L_y
real(rp) :: fringe_region_end_y  = 1.0_rp
! length of fringe region as a fraction of L_y
real(rp) :: fringe_region_len_y = 0.125_rp

! Use uniform inflow instead of concurrent precursor inflow
logical :: uniform_inflow = .false.
real(rp) :: inflow_velocity = 1.0_rp

! if true, imposes a pressure gradient in the x-direction to force the flow
logical :: use_mean_p_force = .true.
real(rp) :: mean_p_force = 1.0_rp
real(rp) :: u_star = 1.0_rp
  
!---------------------------------------------------
! DATA OUTPUT PARAMETERS
!---------------------------------------------------
! how often to display output 
integer :: wbase = 100

! Flags for controling checkpointing data
logical :: checkpoint_data = .false.
integer :: checkpoint_nskip = 10000

! records time-averaged data to files ./output/*_avg.dat
logical :: tavg_calc = .false.
integer :: tavg_nstart = 1,  tavg_nskip = 100
real(rp) :: tavg_start_time = 1.0_rp, tavg_end_time = 2.0_rp

! domain instantaneous output
logical :: domain_calc=.false.
integer :: itype = 3 
integer :: domain_nstart=10000,domain_nskip=10000

! pdf-spec output:
real(lprec) :: pdf_total_time,pdf_dt
real(lprec) :: spec_total_time,spec_dt

! pdf bins:
integer :: nbins = 50
real(rp) :: bin_min_u = 10
real(rp) :: bin_max_u = 40
real(rp) :: bin_min_v = -30
real(rp) :: bin_max_v = 30
real(rp) :: bin_min_w = -30
real(rp) :: bin_max_w = 30
end module param


module interp
use param, only: rp
! Parameters for interpolation
real(rp):: xdiff,ydiff,zdiff
integer:: istart,jstart,kstart,istart1,jstart1,kstart1
end module interp

module filters 
use param, only: rp
real(rp),dimension(:,:),allocatable :: G_test,G_test_test
real(rp):: delta
end module filters

module fft
use param, only: rp
use iso_c_binding
implicit none
include'fftw3.f'
real(rp), allocatable, dimension(:,:) :: kx,ky
integer*8 :: forw,back,forw_big,back_big
#ifdef TAVG_SPECTRUM
real(rp), dimension(:), allocatable :: dummy_in,dummy2_in
complex(rp), dimension(:), allocatable :: dummy_out,dummy2_out
integer*8 :: forw_1d
#endif TAVG_SPECTRUM
end module fft

module mpi_defs
implicit none
integer, parameter :: mpi_SYNC_DOWN=1
integer, parameter :: mpi_SYNC_UP=2
integer, parameter :: mpi_SYNC_DOWNUP=3
end module mpi_defs

module grid_defs
use param, only: rp
implicit none
real(rp),allocatable,dimension(:) :: gridx,gridy,gridz,gridzw,awrap_i,awrap_j,z_tot
end module grid_defs

module io
use param, only : rp,lprec 
implicit none
real(lprec) :: tavg_total_time,tavg_dt
logical :: tavg_initialized = .false.

! Windfarm
type tavg_t
#ifdef STATSTAU
 real(lprec) :: txx,tyy,tzz,txy,txz,tyz
#endif
 real(lprec) :: pres,u,u2,v,v2,w,w2,uv,uw,vw
#ifdef SCALAR
 real(lprec) :: csopt2,theta,T2,sgst3,dsopt2
 !Horizontal and vertical fluxes of temperature and pressure energy
 real(lprec) :: buoy, wT 
#endif
#ifdef TURBINES 
 real(lprec) :: fx, fy, fu, fv 
#endif TURBINES
#ifdef ATM 
 real(lprec) :: fxa, fya,fza, fua, fva,fwa
#endif ATM 
 
#ifdef STATSKE
 !Terms in turbulent kinetic energy equation
 real(lprec) :: u3, v2u, w2u, u2v, v3, w2v, u2w, v2w, w3
 real(lprec) :: txxs11, txys12, txzs13, tyys22, tyzs23, tzzs33
 real(lprec) :: s11, s12, s13, s22, s23, s33 
 real(lprec) :: utxx, vtyx, wtzx, utxy, vtyy, wtzy, utxz, vtyz, wtzz 
 real(lprec) :: pu, pv, pw 
#endif STATSKE 
 
end type tavg_t

type(tavg_t),allocatable,dimension(:,:,:) :: tavg
#ifdef SCALAR
real(rp),allocatable,dimension(:) :: tavg_ustar
#endif SCALAR  

#if defined(TAVG_PDF)
logical :: pdf_initialized = .false.
real(rp),allocatable,dimension(:,:) :: tavg_u_histogram,tavg_v_histogram,tavg_w_histogram
real(rp),allocatable,dimension(:) :: tavg_u_xy,tavg_u2_xy,tavg_v_xy,tavg_v2_xy,tavg_w_xy, tavg_w2_xy
real(rp) :: binsize_u, binsize_v, binsize_w
integer  :: startbin_u,startbin_v,startbin_w
#endif 

#ifdef TAVG_SPECTRUM
logical :: spec_initialized = .false.
real(rp),allocatable,dimension(:,:) :: tavg_u_stream_spec, tavg_v_stream_spec, tavg_w_stream_spec
real(rp),allocatable,dimension(:,:) :: tavg_u_span_spec, tavg_v_span_spec, tavg_w_span_spec
real(rp) :: dkx,u_mean,real_energy
#endif TAVG_SPECTRUM 
end module io 

module convec_mod
use param, only : rp
implicit none
real(rp),allocatable,dimension(:,:,:) :: u1_big,u2_big,u3_big
real(rp),allocatable,dimension(:,:,:) :: vort1_big,vort2_big,vort3_big
real(rp),allocatable,dimension(:,:) :: cc_big
end module convec_mod

#ifdef CPS
module cps_mod
use param, only : rp
implicit none
integer, parameter :: RED=0 ! Upstream domain (producer)
integer, parameter :: BLUE=1 ! Downstream domain (consumer)
integer :: interComm, color
integer :: nx_conp
#ifdef CPS_Y 
integer :: ny_conp
#endif CPS_Y


type vel_sample_type
   integer, allocatable, dimension(:) :: iwrap
   real(rp), allocatable, dimension(:,:,:) :: u, v, w
#ifdef SCALAR
   real(rp), allocatable, dimension(:,:,:) :: theta 
#endif SCALAR

#ifdef CPS_Y 
   integer, allocatable, dimension(:) :: iwrap_y 
   real(rp), allocatable, dimension(:,:,:) :: u_y, v_y, w_y
#ifdef SCALAR
   real(rp), allocatable, dimension(:,:,:) :: theta_y
#endif SCALAR
#endif CPS_Y

end type vel_sample_type

type(vel_sample_type), target :: vel_sample_t

! Weights used in fringe region
real(kind=16), allocatable, dimension(:) :: alpha, beta
#ifdef CPS_Y 
real(kind=16), allocatable, dimension(:) :: alpha_y, beta_y
#endif CPS_Y 

end module cps_mod
#endif CPS

#ifdef LVLSET
module level_set_mod
use param, only : rp
implicit none
real(rp) :: zwall = 0._rp
real(rp) :: zo_level_set = 0.0001_rp ! roughness of IB surface
real(rp) :: coeff_phi_cutoff = 1._rp
real(rp) :: coeff_phi_band = 1._rp
real(rp) :: phi_cutoff, phi_band
real(rp), allocatable, dimension(:,:,:) :: phi_u, normx_u, normy_u, normz_u
real(rp), allocatable, dimension(:,:,:) :: phi_w, normx_w, normy_w, normz_w
end module level_set_mod
#endif LVLSET
