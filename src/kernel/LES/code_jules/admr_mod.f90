! Information on this routine can be found on:
! https://palm.muk.uni-hannover.de/trac/wiki/doc/tec/wtm#wu
!
! For now, this routine is hardcoded for NREL-5MW turbines
! Necessary input-files:
! - admr_input.dat (number and position of turbines can be changed here)
! - NREL_DATA.dat  (airfoil data for NREL-5MW)
!
! routine separated in:
! admr_read_input.f90
! admr.f90
! admr_nodes.f90 (I still need to check if this was a good idea)
!
#ifdef ADMR
module admr_mod
use param, only : rp,nproc,nx,ny,nz,dx,dy,dz,idx,idy,idz,dt,nz_tot,pi,z_i
use param, only : mpi_rp, rank_of_coord, coord, status, ierr
implicit none

integer ::  nturb =1                   !< number of turbines
integer ::  admr_tbase = 20           !< number of timesteps to collect data for output file
integer ::  admr_nodes = 0
logical ::  admr_nodes_allocated = .false.
logical ::  admr_in_proc = .false.
real(rp)::  T_avg_admr
real(rp)::  eps_kernel
real(rp)::  pro_radius

! Variables for one turbine 
type admr_t
 real(rp)  :: rr,dia,delta_r
 real(rp)  :: torque_total !< [m^5/s^2]
 real(rp)  :: thrust_rotor !< [m^4/s^2]
 real(rp)  :: omega_rot,phi_yaw
 real(rp),dimension(3) :: loc              ! location vector
 real(rp),dimension(2) :: nhat            ! unit normal for each turbine
 real(rp),dimension(3)  :: rotx,roty,rotz
 real(rp),allocatable,dimension(:)  :: smear,cur_r,alpha_attack_in,chord
 integer :: num_nodes,nrings,nsegs_total                     ! number of nodes associated with each turbine
 integer :: i_hub,j_hub,k_hub,i_smear,j_smear,k_smear,min_k,max_k
 integer,allocatable,dimension(:,:) :: nodes ! (i,j,k) of each included node
 integer,allocatable,dimension(:,:) :: ringseg ! (i,j,k) of each included node
 integer,allocatable,dimension(:)   :: nsegs
end type admr_t

type(admr_t), pointer, dimension(:) :: admr

integer ::  nairfoils = 8             !< number of airfoils of the used turbine model (for ADM-R and ALM)
logical ::  adjust_yaw   = .FALSE.   !< switch for use of yaw controller

real(rp) :: segment_length  = 1.0_rp            !< length of the segments, the rotor area is divided into
                                                    !< (in tangential direction, as factor of MIN(dx,dy,dz))
real(rp) :: segment_width   = 0.5_rp            !< width of the segments, the rotor area is divided into
                                                    !< (in radial direction, as factor of MIN(dx,dy,dz))
real(rp) :: tilt            = 0.0_rp           !< vertical tilt of the rotor [degree] ( positive = backwards )
real(rp) :: delta_t_factor


!-- Variables for initialization of the turbine model

integer ::  inot         !< turbine loop index (turbine id)
integer ::  nsegs_max    !< maximum number of segments (all turbines, required for allocation of arrays)
integer ::  nrings_max   !< maximum number of rings (all turbines, required for allocation of arrays)
integer ::  ring         !< ring loop index (ring number)
integer ::  rr_int       !<
integer ::  upper_end    !<

integer, dimension(1) ::  lct   !<

INTEGER, DIMENSION(:,:), ALLOCATABLE ::  nsegs   !< number of segments per ring and turbine

!-- Variables for the calculation of lift and drag coefficients
REAL(rp), DIMENSION(:), ALLOCATABLE  ::  ard     !< [deg]
REAL(rp), DIMENSION(:), ALLOCATABLE  ::  crd     !< [m]
REAL(rp), DIMENSION(:), ALLOCATABLE  ::  delta_r !< radial segment length
REAL(rp), DIMENSION(:), ALLOCATABLE  ::  lrd     !< [m]

REAL(rp) ::  accu_cl_cd_tab = 0.1_rp  !< Accuracy of the interpolation of 
                                      !< the lift and drag coeff [deg] 

REAL(rp), DIMENSION(:,:), ALLOCATABLE :: turb_cd_tab   !< table of the blade drag coefficient [-]
REAL(rp), DIMENSION(:,:), ALLOCATABLE :: turb_cl_tab   !< table of the blade lift coefficient [-]

!-- Variables for the calculation of the forces
     
REAL(rp) ::  delta_t                     !<  tangential segment length
REAL(rp) ::  phi_rotor                   !< 
REAL(rp) ::  pre_factor                  !< [m^4/s^2]  
REAL(rp) ::  torque_seg                  !< [m^4/s^2]
REAL(rp) ::  u_int_l                     !< [m/s]
REAL(rp) ::  u_int_u                     !< [m/s]
REAL(rp) ::  u_rot                       !< [m/s]
REAL(rp) ::  v_int_l                     !< [m/s]
REAL(rp) ::  v_int_u                     !< [m/s]
REAL(rp) ::  w_int_l                     !< [m/s]
REAL(rp) ::  w_int_u                     !< [m/s]

REAL(rp), DIMENSION(:), ALLOCATABLE ::  omega_gen    !< curr. generator speed

REAL(rp), DIMENSION(:,:), ALLOCATABLE ::  rbx, rby, rbz     !< coordinates of the blade elements [m]


!-  Fields for the interpolation of velocities on the rotor grid
REAL(rp), DIMENSION(:,:,:), ALLOCATABLE ::  u_int       !< 
REAL(rp), DIMENSION(:,:,:), ALLOCATABLE ::  u_int_1_l   !< 
REAL(rp), DIMENSION(:,:,:), ALLOCATABLE ::  v_int       !< 
REAL(rp), DIMENSION(:,:,:), ALLOCATABLE ::  v_int_1_l   !< 
REAL(rp), DIMENSION(:,:,:), ALLOCATABLE ::  w_int       !< 
REAL(rp), DIMENSION(:,:,:), ALLOCATABLE ::  w_int_1_l   !< 

!
!-  rotor tendencies on the segments 
REAL(rp), DIMENSION(:), ALLOCATABLE :: thrust_seg   !< [m^4/s^2]
REAL(rp), DIMENSION(:), ALLOCATABLE :: torque_seg_y !< [m^4/s^2]
REAL(rp), DIMENSION(:), ALLOCATABLE :: torque_seg_z !< [m^4/s^2]

!
!-  rotor tendencies on the rings
REAL(rp), DIMENSION(:,:), ALLOCATABLE ::  thrust_ring       !< [m^4/s^2] 
REAL(rp), DIMENSION(:,:), ALLOCATABLE ::  torque_ring_y     !< [m^4/s^2]
REAL(rp), DIMENSION(:,:), ALLOCATABLE ::  torque_ring_z     !< [m^4/s^2]
    
!
!-  rotor tendencies on rotor grids for all turbines
REAL(rp), DIMENSION(:,:,:), ALLOCATABLE ::  thrust      !< [m/s^2]
REAL(rp), DIMENSION(:,:,:), ALLOCATABLE ::  torque_y    !< [m/s^2]
REAL(rp), DIMENSION(:,:,:), ALLOCATABLE ::  torque_z    !< [m/s^2]

!
!-  rotor tendencies on coordinate grid
REAL(rp), DIMENSION(:,:,:), ALLOCATABLE ::  rot_tend_x  !< [m/s^2] 
REAL(rp), DIMENSION(:,:,:), ALLOCATABLE ::  rot_tend_y  !< [m/s^2]
REAL(rp), DIMENSION(:,:,:), ALLOCATABLE ::  rot_tend_z  !< [m/s^2]
!    
!-  variables for the rotation of the rotor coordinates        
real(rp), DIMENSION(1:100,1:3,1:3) ::  rot_coord_trans  !< matrix for rotation of rotor coordinates
real(rp), DIMENSION(1:3) ::  rot_eigen_rad   !< 
real(rp), DIMENSION(1:3) ::  rot_eigen_azi   !< 
real(rp), DIMENSION(1:3) ::  rot_eigen_nor   !<

! dimensional dx etc. [m] or idx_dim [1/m]    
real(rp) :: dx_dim, dy_dim, dz_dim 
real(rp) :: idx_dim, idy_dim, idz_dim

! Parameters used for buffering the writing of turbine data
integer,parameter :: admr_max_buf=10 ! Number of timesteps to keep in memory
integer :: admr_count_buf
real(rp),allocatable,dimension(:,:,:) :: admr_data_buf

SAVE
contains
!------------------------------------------------------------------------------!
! Description:
! ------------
!> The projection matrix for the coordinate system of the rotor disc in respect
!> to the yaw and tilt angle of the rotor is calculated
!------------------------------------------------------------------------------!
    SUBROUTINE wtm_rotate_rotor( inot )

      IMPLICIT NONE
      integer :: inot

!--    Calculation of the rotation matrix for the application of the tilt to
!--    the rotors
       rot_eigen_rad(1) = SIN( admr(inot)%phi_yaw )    ! x-component of the radial eigenvector
       rot_eigen_rad(2) = COS( admr(inot)%phi_yaw )    ! y-component of the radial eigenvector 
       rot_eigen_rad(3) = 0.0_rp                  ! z-component of the radial eigenvector

       rot_eigen_azi(1) = 0.0_rp                  ! x-component of the azimuth eigenvector
       rot_eigen_azi(2) = 0.0_rp                  ! y-component of the azimuth eigenvector
       rot_eigen_azi(3) = 1.0_rp                  ! z-component of the azimuth eigenvector

       rot_eigen_nor(1) =  COS( admr(inot)%phi_yaw )   ! x-component of the normal eigenvector
       rot_eigen_nor(2) = -SIN( admr(inot)%phi_yaw )   ! y-component of the normal eigenvector
       rot_eigen_nor(3) = 0.0_rp                  ! z-component of the normal eigenvector

!!
!!--    Calculation of the coordinate transformation matrix to apply a tilt to
!!--    the rotor. If tilt = 0, rot_coord_trans is a unit matrix.
!

!!
!!--    Calculation of the coordinate transformation matrix to apply a tilt to
!!--    the rotor. If tilt = 0, rot_coord_trans is a unit matrix.
!
       rot_coord_trans(inot,1,1) = rot_eigen_rad(1)**2                   *     &
                                   ( 1.0_rp - COS( tilt ) ) + COS( tilt )
       rot_coord_trans(inot,1,2) = rot_eigen_rad(1) * rot_eigen_rad(2)   *     &
                                   ( 1.0_rp - COS( tilt ) )              -     &
                                   rot_eigen_rad(3) * SIN( tilt )
       rot_coord_trans(inot,1,3) = rot_eigen_rad(1) * rot_eigen_rad(3)   *     &
                                   ( 1.0_rp - COS( tilt ) )              +     &
                                   rot_eigen_rad(2) * SIN( tilt )
       rot_coord_trans(inot,2,1) = rot_eigen_rad(2) * rot_eigen_rad(1)   *     &
                                   ( 1.0_rp - COS( tilt ) )              +     &
                                   rot_eigen_rad(3) * SIN( tilt )
       rot_coord_trans(inot,2,2) = rot_eigen_rad(2)**2                   *     &
                                   ( 1.0_rp - COS( tilt ) ) + COS( tilt )
       rot_coord_trans(inot,2,3) = rot_eigen_rad(2) * rot_eigen_rad(3)   *     &
                                   ( 1.0_rp - COS( tilt ) )              -     &
                                   rot_eigen_rad(1) * SIN( tilt )
       rot_coord_trans(inot,3,1) = rot_eigen_rad(3) * rot_eigen_rad(1)   *     &
                                   ( 1.0_rp - COS( tilt ) )              -     &
                                   rot_eigen_rad(2) * SIN( tilt )
       rot_coord_trans(inot,3,2) = rot_eigen_rad(3) * rot_eigen_rad(2)   *     &
                                   ( 1.0_rp - COS( tilt ) )              +     &
                                   rot_eigen_rad(1) * SIN( tilt )
       rot_coord_trans(inot,3,3) = rot_eigen_rad(3)**2                   *     &
                                   ( 1.0_rp - COS( tilt ) ) + COS( tilt )


!--    Vectors for the Transformation of forces from the rotor's spheric
!--    coordinate system to the cartesian coordinate system
       admr(inot)%rotx(:) = MATMUL( rot_coord_trans(inot,:,:), rot_eigen_nor )
       admr(inot)%roty(:) = MATMUL( rot_coord_trans(inot,:,:), rot_eigen_rad )
       admr(inot)%rotz(:) = MATMUL( rot_coord_trans(inot,:,:), rot_eigen_azi )

    END SUBROUTINE wtm_rotate_rotor
end module

!*******************************************************************************
subroutine wtm_init_arrays
!*******************************************************************************
use param, only : rp,nproc,nx,ny,nz,dx,dy,dz,idx,idy,idz,dt,nz_tot,pi,z_i
use admr_mod
implicit none

!--    Allocate 1D arrays (dimension = number of rotor segments)
ALLOCATE( thrust_seg(   1:nsegs_max) )
ALLOCATE( torque_seg_y( 1:nsegs_max) )
ALLOCATE( torque_seg_z( 1:nsegs_max) )

!--    Allocate 2D arrays (dimension = number of rotor rings and segments)
ALLOCATE( rbx(1:nrings_max,1:nsegs_max) )
ALLOCATE( rby(1:nrings_max,1:nsegs_max) )
ALLOCATE( rbz(1:nrings_max,1:nsegs_max) )
ALLOCATE( thrust_ring(   1:nrings_max,1:nsegs_max) )
ALLOCATE( torque_ring_y( 1:nrings_max,1:nsegs_max) )
ALLOCATE( torque_ring_z( 1:nrings_max,1:nsegs_max) )

!--    Allocate 3D arrays (dimension = number of grid points)
ALLOCATE( rot_tend_x( 1:nx,1:ny,1:nz) )
ALLOCATE( rot_tend_y( 1:nx,1:ny,1:nz) )
ALLOCATE( rot_tend_z( 1:nx,1:ny,1:nz) )
ALLOCATE( thrust(     1:nx,1:ny,1:nz) )
ALLOCATE( torque_y(   1:nx,1:ny,1:nz) )
ALLOCATE( torque_z(   1:nx,1:ny,1:nz) )

ALLOCATE( u_int(     1:nturb,1:nrings_max,1:nsegs_max) )
ALLOCATE( u_int_1_l( 1:nturb,1:nrings_max,1:nsegs_max) )
ALLOCATE( v_int(     1:nturb,1:nrings_max,1:nsegs_max) )
ALLOCATE( v_int_1_l( 1:nturb,1:nrings_max,1:nsegs_max) )
ALLOCATE( w_int(     1:nturb,1:nrings_max,1:nsegs_max) )
ALLOCATE( w_int_1_l( 1:nturb,1:nrings_max,1:nsegs_max) )


thrust_seg(:)            = 0.0_rp
torque_seg_y(:)          = 0.0_rp
torque_seg_z(:)          = 0.0_rp

rbx(:,:)                 = 0.0_rp
rby(:,:)                 = 0.0_rp
rbz(:,:)                 = 0.0_rp
thrust_ring(:,:)         = 0.0_rp
torque_ring_y(:,:)       = 0.0_rp
torque_ring_z(:,:)       = 0.0_rp


rot_tend_x(:,:,:)        = 0.0_rp
rot_tend_y(:,:,:)        = 0.0_rp
rot_tend_z(:,:,:)        = 0.0_rp
thrust(:,:,:)            = 0.0_rp
torque_y(:,:,:)          = 0.0_rp
torque_z(:,:,:)          = 0.0_rp

end subroutine wtm_init_arrays

!*******************************************************************************
subroutine admr_init
!*******************************************************************************
use param, only: z_i,dx,dy,dz,nz,pi,coord,rp,path,L_z,nx,ny,nz,pi
use param, only : nproc,idx,idy,idz,dt,nz_tot
use admr_mod
implicit none
integer, parameter :: lun=42
real :: b,c,d,e,f
real(rp) :: x_dist,y_dist
integer :: a,j,k_start,k_end
integer :: min_i,max_i,min_j,max_j
integer ::  i_rseg,i_ring  !<
real(rp) ::  delta_r_factor   !< 
real(rp) ::  delta_r_init     !< 


! Define the turbine array
nullify(admr)
allocate(admr(nturb))

!Global z indices 
k_start= 1   +coord*(nz-1)
k_end  = nz-1+coord*(nz-1)

admr(:)%omega_rot = 0.96_rp
segment_length    = 1.0_rp
segment_width     = 0.5_rp
nairfoils         = 8

dx_dim  = dx * z_i
dy_dim  = dy * z_i
dz_dim  = dz * z_i
idx_dim = 1._rp / dx_dim
idy_dim = 1._rp / dy_dim
idz_dim = 1._rp / dz_dim

!
!--    To be able to allocate arrays with dimension of rotor rings and segments,
!--    the maximum possible numbers of rings and segments have to be calculated:

admr(:)%nrings  = 0
admr(:)%delta_r = 0.0_rp

!
!--    Thickness (radial) of each ring and length (tangential) of each segment:
delta_r_factor = segment_width
delta_t_factor = segment_length
delta_r_init   = delta_r_factor * MIN( dx_dim, dy_dim, dz_dim)
delta_t        = delta_t_factor * MIN( dx_dim, dy_dim, dz_dim)


! Read the data from the input table
if (coord.eq.0) write(*,*) '   x_loc    y_loc     Diameter  Height    Thick Vol       Theta     Ctprime'

! Read turbine locations
open(lun,file=trim(path)//'admr_input.dat',status='unknown',form='formatted',action='read',position='rewind')
read(lun,*)
do inot = 1,nturb
  read(lun,*) a,b,c,d,e,f
  admr(inot)%loc(1)  = b
  admr(inot)%loc(2)  = c
  admr(inot)%rr      = d /2._rp
  admr(inot)%loc(3)  = e/z_i
  admr(inot)%phi_yaw = f * pi / 180.0_rp !read in degrees,transformed to rad
  tilt = tilt * pi / 180.0_rp

  if (coord.eq.0) write(*,"(8F10.4)") admr(inot)%loc(1),admr(inot)%loc(2),admr(inot)%rr,admr(inot)%loc(3),admr(inot)%phi_yaw

  ! Index of turbines center:
  admr(inot)%i_hub  = nint(admr(inot)%loc(1)/dx)
  admr(inot)%j_hub  = nint(admr(inot)%loc(2)/dy)
  admr(inot)%k_hub  = nint(admr(inot)%loc(3)/dz)

!-- Determining the area to which the smearing of the forces is applied.
!-- As smearing now is effectively applied only for distances smaller than
!-- eps_min, the smearing area can be further limited and regarded as a
!-- function of eps_min:
  admr(inot)%i_smear = CEILING( (admr(inot)%rr / dx_dim) + 2._rp )
  admr(inot)%j_smear = CEILING( (admr(inot)%rr / dy_dim) + 2._rp )
  admr(inot)%k_smear = CEILING( (admr(inot)%rr / dz_dim) + 7._rp )

  ! Save index for loop in turbine_nodes_func (see below)
  min_i = admr(inot)%i_hub-admr(inot)%i_smear
  max_i = admr(inot)%i_hub+admr(inot)%i_smear
  min_j = admr(inot)%j_hub-admr(inot)%j_smear
  max_j = admr(inot)%j_hub+admr(inot)%j_smear
  admr(inot)%min_k = max(admr(inot)%k_hub-admr(inot)%k_smear,k_start)
  admr(inot)%max_k = min(admr(inot)%k_hub+admr(inot)%k_smear,k_end)
!
!--       Determine number of rings:
 admr(inot)%nrings  = NINT(admr(inot)%rr / delta_r_init )
 admr(inot)%delta_r = admr(inot)%rr / float(admr(inot)%nrings)
enddo

nrings_max = MAXVAL(admr(:)%nrings)
admr(:)%nsegs_total = 0
nsegs_max = 0

do inot = 1, nturb
allocate( admr(inot)%alpha_attack_in( 1:nrings_max) )
allocate( admr(inot)%chord(           1:nrings_max) )
allocate( admr(inot)%nsegs(1:nrings_max) )
allocate( admr(inot)%cur_r(1:nrings_max) )
admr(inot)%nsegs(:)     = 0
admr(inot)%cur_r(:)     = 0

 do i_ring = 1, admr(inot)%nrings

!--          Determine number of segments for each ring:
  admr(inot)%nsegs(i_ring) = MAX( 8, CEILING( delta_r_factor * pi *         &
                                      ( 2.0_rp * i_ring - 1.0_rp ) /  &
                                                delta_t_factor ) )
  admr(inot)%cur_r(i_ring) = (float(i_ring) - 0.5_rp) * admr(inot)%delta_r
 enddo

!--       Total sum of all rotor segments:
 admr(inot)%nsegs_total = SUM( admr(inot)%nsegs(:) )

 if (MAXVAL(admr(inot)%nsegs(:)) > nsegs_max)  nsegs_max = MAXVAL(admr(inot)%nsegs(:))
enddo

do inot = 1,nturb
  !Check whether the turbine is in this processor
  if(admr(inot)%k_hub+admr(inot)%k_smear.ge.k_start .and. &
     admr(inot)%k_hub-admr(inot)%k_smear.le.k_end) then
   write(*,*) "admr in proc:", coord
   admr_in_proc = .true.
   admr_nodes   = 4*nsegs_max*nrings_max*(admr(inot)%i_smear+1)*(admr(inot)%j_smear+1)*(nz-1)
   allocate(admr(inot)%nodes(3,admr_nodes))
   allocate(admr(inot)%smear(admr_nodes))
   allocate(admr(inot)%ringseg(2,admr_nodes))
   admr(inot)%num_nodes = 0
   admr(inot)%nodes     = 0._rp
   admr(inot)%smear     = 0._rp
  end if

  ! Check if turbines are in domain / not too close to border of domain:
  if( admr(inot)%i_hub <= admr(inot)%i_smear .or. &
      admr(inot)%i_hub >= nx - admr(inot)%i_smear ) stop 'Turbine out of x domain.'
  if( admr(inot)%j_hub <= admr(inot)%j_smear .or. &
      admr(inot)%j_hub >= nx - admr(inot)%j_smear ) stop 'Turbine out of y domain.'
  !if( turbine(k)%loc(3) + 0.5_rp*turbine(k)%dia >  0.5_rp * L_z ) stop 'Turbine
  !too high.'

enddo
close(lun)
! End read the data from the input table 

! Check if spacing of turbines is at least 2.5 times turbine diameter:
do inot = 1,nturb
  do j = 1,nturb
    if (inot .ne. j) then
      x_dist = abs(admr(inot)%loc(1) - admr(j)%loc(1))
      y_dist = abs(admr(inot)%loc(2) - admr(j)%loc(2))
      if(sqrt(x_dist**2 + y_dist**2) < 5._rp * admr(inot)%rr) stop 'Turbine spacing too small.'
    endif
  enddo
enddo

! Check to see whether the grid spacing is smaller than the turbine radius
do inot = 1,nturb
if (max(dx_dim,dy_dim,dz_dim).ge.admr(inot)%rr) then
 write(*,'(5F12.4)') inot,dx_dim,dy_dim,dz_dim,admr(inot)%rr
 stop 'Grid way too coarse to accurate resolve turbine'
endif
enddo

call wtm_init_arrays

!! Initialize buffer feature for turbine data
admr_count_buf=0
allocate(admr_data_buf(6,nturb,admr_max_buf)); admr_data_buf=0.0_rp

do inot = 1, nturb
 call wtm_rotate_rotor( inot )
enddo
!- read blade geometry and airfoil tables
call admr_read_blade_tables

do inot = 1, nturb
do i_ring = 1,admr(inot)%nrings
 do i_rseg = 1, admr(inot)%nsegs(i_ring)
    phi_rotor = float(i_rseg) * 2.0_rp * pi / float(admr(inot)%nsegs(i_ring))
    rbx(i_ring,i_rseg) = admr(inot)%loc(1) * z_i + admr(inot)%cur_r(i_ring) * COS( phi_rotor ) * SIN(admr(inot)%phi_yaw)
    rby(i_ring,i_rseg) = admr(inot)%loc(2) * z_i + admr(inot)%cur_r(i_ring) * COS( phi_rotor ) * COS(admr(inot)%phi_yaw)
    rbz(i_ring,i_rseg) = admr(inot)%loc(3) * z_i + admr(inot)%cur_r(i_ring) *SIN( phi_rotor )
 enddo
enddo
end do
call admr_nodes_func
end subroutine admr_init

#endif ADMR
