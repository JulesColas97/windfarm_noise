#ifdef ATM
!*******************************************************************************
module atm_lesgo_interface
!*******************************************************************************
!
! This module interfaces actuator turbine module with lesgo. It is a lesgo 
! specific module, unlike the atm module. The MPI management is done only in 
! this section of the code. This is very code dependent and will have to be 
! modified according to the code being used. In this case LESGO has its own MPI 
! details. Look into mpi_defs.f90 for the details
!
use param, only: rp
implicit none
save
private
public :: atm_lesgo_initialize
public :: atm_lesgo_forcing
public :: atm_lesgo_finalize

! Variable for interpolating the velocity in w onto the uv grid
real(rp), allocatable, dimension(:,:,:) :: w_uv

! This is a list that stores all the points in the domain with a body force due
! to the turbines.
type bodyForce_t
  integer :: c                           ! Number of cells
  integer, allocatable :: ijk(:,:)       ! index for the point in the domain
  real(rp), allocatable :: force(:,:)    ! Force vector on uv grid
  real(rp), allocatable :: location(:,:) ! Position vector on uv grid
end type bodyForce_t

! Body force field
type(bodyForce_t), allocatable, target, dimension(:) :: forceFieldUV,forceFieldW

contains

!*******************************************************************************
subroutine atm_lesgo_initialize 
!*******************************************************************************
!
! Initialize the actuator turbine model (ATM)
!
use param, only: ld, ny, nz, coord
use atm_input_util, only: numberOfTurbines
use actuator_turbine_model, only: atm_initialize
use actuator_turbine_model, only: atm_initialize_output
implicit none
integer :: i

! Initialize the atm 
call atm_initialize 

! Allocate space for the w_uv variable
allocate(w_uv(ld,ny,0:nz))

! Allocate the body force variables. It is an array with per turbine.
allocate(forceFieldUV(numberOfTurbines))
allocate(forceFieldW(numberOfTurbines))

! find all the cells that surround the turbines
do i=1,numberOfTurbines
  call atm_lesgo_findCells(i)
enddo

! This will create the output files and write initialization to the screen
if(coord == 0) call atm_initialize_output

end subroutine atm_lesgo_initialize


!*******************************************************************************
subroutine atm_lesgo_finalize 
!*******************************************************************************
!
! Initialize the actuator turbine model
!
use param, only: coord
use atm_input_util, only: turbineArray
use atm_input_util, only: numberOfTurbines
use actuator_turbine_model, only: atm_write_restart
implicit none
integer :: i

if(coord == 0) write(*,*) 'Finalizing ATM...'

! Loop through all turbines and finalize
do i = 1, numberOfTurbines
  if(coord == turbineArray(i) % master) then
    ! Write the restart file
    call atm_write_restart(i) 
  endif
end do

if(coord == 0) write(*,*) 'Done finalizing ATM'

end subroutine atm_lesgo_finalize


!*******************************************************************************
subroutine atm_lesgo_findCells(m)
!*******************************************************************************
!
! This subroutine finds all the cells that surround the turbines
!
use param, only: rp, coord, ierr, nz_tot, z_i
use param, only: nx, ny, nz, localComm
use grid_defs, only: gridx, gridy, gridz, gridzw 
use atm_input_util, only: turbineArray
use mpi, only: MPI_INTEGER, MPI_MIN, MPI_SUM
implicit none
integer, intent(in) :: m  ! The turbine number
integer :: cUV, cW        ! Counters for points affected on UV and W grids
integer :: i, j, k
real(rp), dimension(3) :: vector_point ! Vector used to store x, y, z locations
integer :: base_group      ! The base group from comm --> MPI_COMM_WORLD 
integer :: local_group     ! The local group of processors
integer :: member          ! (1 or 0) yes or no
integer :: num_of_members  ! total number of members
! List of all the cores that belong to this turbine
integer, allocatable, dimension(:) :: ls_of_cores
real(rp) :: dist, a(3), b(3)

! Initialize internal counter to zero
forceFieldUV(m) % c = 0

! Find all locations that are influenced by each turbine
cUV = 0
cW = 0
do i=1,nx
do j=1,ny
do k=1,nz
  vector_point(1) = gridx(i)*z_i
  vector_point(2) = gridy(j)*z_i

  ! Take into account the UV grid
  vector_point(3) = gridz(k)*z_i

  a = vector_point
  b = turbineArray(m) % towerShaftIntersect
  dist = sqrt(dot_product(a-b,a-b))
  if(dist .le. turbineArray(m) % sphereRadius) then
    cUV = cUV+1
  endif

  ! Take into account the W grid
  vector_point(3) = gridzw(k)*z_i

  a = vector_point
  b = turbineArray(m) % towerShaftIntersect
  dist = sqrt(dot_product(a-b,a-b))
  if(dist .le. turbineArray(m) % sphereRadius) then
    cW = cW+1
  endif
enddo
enddo
enddo

! Allocate space for the force fields in UV and W grids
forceFieldUV(m) % c = cUV ! Counter
allocate(forceFieldUV(m) % force(3,cUV))
allocate(forceFieldUV(m) % location(3,cUV))
allocate(forceFieldUV(m) % ijk(3,cUV))

forceFieldW(m) % c = cW  ! Counter
allocate(forceFieldW(m) % force(3,cW))
allocate(forceFieldW(m) % location(3,cW))
allocate(forceFieldW(m) % ijk(3,cW))

print*, 'Number of cells affected by turbine', m, 'cUV, cW =', cUV, cW

! Run the same loop and save all variables
! The forceField arrays include all the forces which affect the domain
cUV = 0
cW = 0
do i=1,nx
do j=1,ny
do k=1,nz
  vector_point(1) = gridx(i)*z_i
  vector_point(2) = gridy(j)*z_i

  ! Take into account the UV grid
  vector_point(3) = gridz(k)*z_i

  a = vector_point
  b = turbineArray(m) % towerShaftIntersect
  dist = sqrt(dot_product(a-b,a-b))
  if(dist .le. turbineArray(m) % sphereRadius) then
    cUV = cUV+1
    forceFieldUV(m) % ijk(1,cUV) = i
    forceFieldUV(m) % ijk(2,cUV) = j
    forceFieldUV(m) % ijk(3,cUV) = k
    forceFieldUV(m) % location(1:3,cUV) = vector_point(1:3)
    forceFieldUV(m) % force(1:3,cUV) = 0._rp
  endif

  ! Take into account the W grid
  vector_point(3) = gridzw(k)*z_i

  a = vector_point
  b = turbineArray(m) % towerShaftIntersect
  dist = sqrt(dot_product(a-b,a-b))
  if(dist .le. turbineArray(m) % sphereRadius) then
    cW = cW+1
    forceFieldW(m) % ijk(1,cW) = i
    forceFieldW(m) % ijk(2,cW) = j
    forceFieldW(m) % ijk(3,cW) = k
    forceFieldW(m) % location(1:3,cW) = vector_point(1:3)
    forceFieldW(m) % force(:,cW) = 0._rp
  endif
enddo
enddo
enddo

! Store the base group from the global communicator mpi_comm_world
call MPI_COMM_GROUP(localComm, base_group, ierr)

! Assign member
member = 0
! Flag to know if this turbine is operating or not
turbineArray(m) % operate = .FALSE.

! Assign proper values if turbine affects processors in this region
if(cUV>0 .or. cW>0) then
  member = 1
  turbineArray(m) % operate = .TRUE.
endif

! Find the total number of processors for each turbine
call mpi_allreduce(member, num_of_members, 1, MPI_INTEGER, MPI_SUM,            &
                   localComm, ierr)

if (turbineArray(m) % operate) then
  ! Find the master processor for each turbine
  call mpi_allreduce(coord, turbineArray(m) % master, 1, MPI_INTEGER, MPI_MIN, &
                     localComm, ierr)
else
  ! This is bogus since nz will always be less than number of processors
  call mpi_allreduce(nz_tot, turbineArray(m) % master, 1, MPI_INTEGER, MPI_MIN,&
                     localComm, ierr)
endif

allocate(ls_of_cores(num_of_members))
ls_of_cores(1) = turbineArray(m) % master

! Notice this list is valid only for decomposition in 1 direction
do i=2,num_of_members
  ls_of_cores(i) = ls_of_cores(i-1) + 1
enddo

! Write if this processor is the master
if(coord == turbineArray(m) % master) then
  write(*,*) 'Master for turbine', m, 'is processor', turbineArray(m) % master
endif

! Create the new communicator and group for this turbine
call MPI_GROUP_INCL(base_group, num_of_members, ls_of_cores, local_group, ierr)
call MPI_COMM_CREATE(localComm, local_group, turbineArray(m) % TURBINE_COMM_WORLD, ierr)

if(turbineArray(m) % operate) then
  write(*,*) 'Processor', coord, 'has elements in turbine   ', m
else
  write(*,*) 'Processor', coord, 'has NO elements in turbine', m
endif

call mpi_barrier(localComm, ierr)

end subroutine atm_lesgo_findCells


!*******************************************************************************
subroutine atm_lesgo_forcing
!*******************************************************************************
!
! This subroutines calls the update function from the ATM Library and calculates 
! the body forces needed in the domain
!
use param, only: localComm
use mpi_defs, only: MPI_SYNC_DOWNUP
use param, only: rp, ierr, total_time, dt, z_i, coord
use param, only: nz, jt_total, coord, nproc
use sim_param, only: w
use atm_input_util, only: turbineArray
use atm_input_util, only: updateInterval
use atm_input_util, only: numberOfTurbines
use actuator_turbine_model, only: atm_computeRotorSpeed
use actuator_turbine_model, only: atm_rotateBlades
use actuator_turbine_model, only: atm_control_yaw
use actuator_turbine_model, only: atm_output
implicit none
integer :: i

! Get the velocity from w onto the uv grid
w_uv(:,:,1:nz-1) = 0.5_rp*(w(:,:,2:nz) + w(:,:,1:nz-1))
! Take care of top "physical" boundary
if(coord == nproc-1) w_uv(:,:,nz) = w_uv(:,:,nz-1)
call mpi_sync_real_array(w_uv, MPI_SYNC_DOWNUP)

! Loop through all turbines and rotate the blades
do i=1,numberOfTurbines
  if(turbineArray(i) % operate) then
    ! Compute the deltaAzimuth
    call atm_computeRotorSpeed(i, dt*z_i)
#ifndef ADMR
    ! Rotate the blades
    call atm_rotateBlades(i)
#endif
    ! yaw the nacelle (if applicable)
    call atm_control_yaw(i, total_time*z_i)
  endif
enddo

! Only calculate new forces if interval is correct
if(mod(jt_total-1, updateInterval) == 0) then

  ! Establish all turbine properties as zero, essential for lesgo_mpi
  do i=1,numberOfTurbines
    turbineArray(i) % bladeForces = 0._rp
    turbineArray(i) % torqueRotor = 0._rp
    turbineArray(i) % thrust = 0._rp
    turbineArray(i) % alpha = 0._rp
    turbineArray(i) % Cd = 0._rp
    turbineArray(i) % Cl = 0._rp
    turbineArray(i) % lift = 0._rp
    turbineArray(i) % drag = 0._rp
    turbineArray(i) % Vmag = 0._rp
    turbineArray(i) % windVectors = 0._rp
    turbineArray(i) % nacelleForce = 0._rp
    turbineArray(i) % induction_a = 0._rp
    turbineArray(i) % u_infinity = 0._rp
    turbineArray(i) % bladeAlignedVectors = 0._rp
    turbineArray(i) % axialForce = 0._rp
    turbineArray(i) % tangentialForce = 0._rp

    if(turbineArray(i) % operate) then
      ! Set body forces to zero
      forceFieldUV(i) % force = 0._rp
      forceFieldW(i) % force = 0._rp
      ! Calculate forces for all turbines
      call atm_lesgo_force(i)
    endif
  enddo

  ! This will gather all the blade forces from all processors
  call atm_lesgo_mpi_gather

  do i=1,numberOfTurbines
    ! Only perform if turbine is active in this processor
    if(turbineArray(i) % operate) then
      ! Convolute force onto the domain
      call atm_lesgo_convolute_force(i)
    endif
  enddo
endif

call atm_lesgo_apply_force

do i=1,numberOfTurbines
  if(coord == turbineArray(i) % master) then
    ! dt has the dimension [s/m]
    call atm_output(i, jt_total, total_time*z_i)
  endif
enddo

! Make sure all processors stop wait for the output to be completed
call mpi_barrier(localComm, ierr)

end subroutine atm_lesgo_forcing


!*******************************************************************************
subroutine atm_lesgo_mpi_gather
!*******************************************************************************
!
! This subroutine will gather the necessary outputs from the turbine models so
! all processors have acces to it. This is done by means of allreduce sum since
! the value is either zero or unique correct value.
use param, only: mpi_rp, ierr, rp
use atm_input_util, only: turbineArray
use atm_input_util, only: numberOfTurbines
#ifdef ADMR
use atm_input_util, only: numberOfActuatorLines
#endif ADMR
use mpi, only: mpi_sum
implicit none
integer :: i
real(rp) :: torqueRotor, thrust
real(rp), dimension(3) :: nacelleForce
integer, pointer :: TURBINE_COMMUNICATOR ! Pointer for MPI communicator

do i=1,numberOfTurbines
  ! Only do MPI sums if processors are operating in this turbine
  if(turbineArray(i) % operate) then
    TURBINE_COMMUNICATOR => turbineArray(i) % TURBINE_COMM_WORLD

    ! Sync all the blade forces
    turbineArray(i) % bladeVectorDummy = turbineArray(i) % bladeForces
    call mpi_allreduce(turbineArray(i) % bladeVectorDummy,                     &
                       turbineArray(i) % bladeForces,                          &
                       size(turbineArray(i) % bladeVectorDummy),               &
                       mpi_rp, mpi_sum, TURBINE_COMMUNICATOR, ierr)
#ifdef ADMR
  turbineArray(i) % bladeForces = turbineArray(i) % bladeForces /&
      real(numberOfActuatorLines)
#endif ADMR

    ! Sync bladeAlignedVectors
    turbineArray(i) % bladeVectorDummy =                                       &
    turbineArray(i) % bladeAlignedVectors(:,:,1,:)
    call mpi_allreduce(turbineArray(i) % bladeVectorDummy,                     &
                       turbineArray(i) % bladeAlignedVectors(:,:,1,:),         &
                       size(turbineArray(i) % bladeVectorDummy),               &
                       mpi_rp, mpi_sum, TURBINE_COMMUNICATOR, ierr)

    turbineArray(i) % bladeVectorDummy =                                       &
    turbineArray(i) % bladeAlignedVectors(:,:,2,:)
    call mpi_allreduce(turbineArray(i) % bladeVectorDummy,                     &
                       turbineArray(i) % bladeAlignedVectors(:,:,2,:),         &
                       size(turbineArray(i) % bladeVectorDummy),               &
                       mpi_rp, mpi_sum, TURBINE_COMMUNICATOR, ierr)

    turbineArray(i) % bladeVectorDummy =                                       &
    turbineArray(i) % bladeAlignedVectors(:,:,3,:)
    call mpi_allreduce(turbineArray(i) % bladeVectorDummy,                     &
                       turbineArray(i) % bladeAlignedVectors(:,:,3,:),         &
                       size(turbineArray(i) % bladeVectorDummy),               &
                       mpi_rp, mpi_sum, TURBINE_COMMUNICATOR, ierr)

    ! Sync alpha
    turbineArray(i) % bladeScalarDummy = turbineArray(i) % alpha
    call mpi_allreduce(turbineArray(i) % bladeScalarDummy,                     &
                       turbineArray(i) % alpha,                                &
                       size(turbineArray(i) % bladeScalarDummy),               &
                       mpi_rp, mpi_sum, TURBINE_COMMUNICATOR, ierr)

    ! Sync lift
    turbineArray(i) % bladeScalarDummy = turbineArray(i) % lift
    call mpi_allreduce(turbineArray(i) % bladeScalarDummy,                     &
                       turbineArray(i) % lift,                                 &
                       size(turbineArray(i) % bladeScalarDummy),               &
                       mpi_rp, mpi_sum, TURBINE_COMMUNICATOR, ierr)

    ! Sync drag
    turbineArray(i) % bladeScalarDummy = turbineArray(i) % drag
    call mpi_allreduce(turbineArray(i) % bladeScalarDummy,                     &
                       turbineArray(i) % drag,                                 &
                       size(turbineArray(i) % bladeScalarDummy),               &
                       mpi_rp, mpi_sum, TURBINE_COMMUNICATOR, ierr)

    ! Sync Cl
    turbineArray(i) % bladeScalarDummy = turbineArray(i) % Cl
    call mpi_allreduce(turbineArray(i) % bladeScalarDummy,                     &
                       turbineArray(i) % Cl,                                   &
                       size(turbineArray(i) % bladeScalarDummy),               &
                       mpi_rp, mpi_sum, TURBINE_COMMUNICATOR, ierr)

    ! Sync Cd
    turbineArray(i) % bladeScalarDummy = turbineArray(i) % Cd
    call mpi_allreduce(turbineArray(i) % bladeScalarDummy,                     &
                       turbineArray(i) % Cd,                                   &
                       size(turbineArray(i) % bladeScalarDummy),               &
                       mpi_rp, mpi_sum, TURBINE_COMMUNICATOR, ierr)

    ! Sync Vmag
    turbineArray(i) % bladeScalarDummy = turbineArray(i) % Vmag
    call mpi_allreduce(turbineArray(i) % bladeScalarDummy,                     &
                       turbineArray(i) % Vmag,                                 &
                       size(turbineArray(i) % bladeScalarDummy),               &
                       mpi_rp, mpi_sum, TURBINE_COMMUNICATOR, ierr)

    ! Sync axialForce
    turbineArray(i) % bladeScalarDummy = turbineArray(i) % axialForce
    call mpi_allreduce(turbineArray(i) % bladeScalarDummy,                     &
                       turbineArray(i) % axialForce,                           &
                       size(turbineArray(i) % bladeScalarDummy),               &
                       mpi_rp, mpi_sum, TURBINE_COMMUNICATOR, ierr)

    ! Sync tangentialForce
    turbineArray(i) % bladeScalarDummy = turbineArray(i) % tangentialForce
    call mpi_allreduce(turbineArray(i) % bladeScalarDummy,                     &
                       turbineArray(i) % tangentialForce,                      &
                       size(turbineArray(i) % bladeScalarDummy),               &
                       mpi_rp, mpi_sum, TURBINE_COMMUNICATOR, ierr)

    ! Sync wind Vectors (Vaxial, Vtangential, Vradial)
    turbineArray(i) % bladeVectorDummy = turbineArray(i) % windVectors(:,:,1:3)
    call mpi_allreduce(turbineArray(i) % bladeVectorDummy,                     &
                       turbineArray(i) % windVectors(:,:,1:3),                 &
                       size(turbineArray(i) % bladeVectorDummy),               &
                       mpi_rp, mpi_sum, TURBINE_COMMUNICATOR, ierr)

    ! Sync induction factor
    turbineArray(i) % bladeScalarDummy = turbineArray(i) % induction_a(:,:)
    call mpi_allreduce(turbineArray(i) % bladeScalarDummy,                     &
                       turbineArray(i) % induction_a(:,:),                     &
                       size(turbineArray(i) % bladeScalarDummy),               &
                       mpi_rp, mpi_sum, TURBINE_COMMUNICATOR, ierr)

    ! Sync u infinity
    turbineArray(i) % bladeScalarDummy = turbineArray(i) % u_infinity(:,:)
    call mpi_allreduce(turbineArray(i) % bladeScalarDummy,                     &
                       turbineArray(i) % u_infinity(:,:),                      &
                       size(turbineArray(i) % bladeScalarDummy),               &
                       mpi_rp, mpi_sum, TURBINE_COMMUNICATOR, ierr)

    ! Store the torqueRotor.
    ! Needs to be a different variable in order to do MPI Sum
    torqueRotor = turbineArray(i) % torqueRotor
    thrust = turbineArray(i) % thrust
    nacelleForce = turbineArray(i) % nacelleForce

    ! Sum all the individual torqueRotor from different blade points
    call mpi_allreduce(torqueRotor, turbineArray(i) % torqueRotor,             &
                       1, mpi_rp, mpi_sum, TURBINE_COMMUNICATOR, ierr)
#ifdef ADMR
  turbineArray(i) % torqueRotor = turbineArray(i) % torqueRotor /&
      real(numberOfActuatorLines)
#endif ADMR

    ! Sum all the individual thrust from different blade points
    call mpi_allreduce(thrust, turbineArray(i) % thrust,                       &
                       1, mpi_rp, mpi_sum, TURBINE_COMMUNICATOR, ierr)
#ifdef ADMR
  turbineArray(i) % thrust = turbineArray(i) % thrust / &
      real(numberOfActuatorLines)
#endif ADMR

    ! Sync the nacelle force
    call mpi_allreduce(nacelleForce, turbineArray(i) % nacelleForce,           &
                       size(turbineArray(i) % nacelleForce), mpi_rp,           &
                       mpi_sum, TURBINE_COMMUNICATOR, ierr)

  endif
enddo

end subroutine atm_lesgo_mpi_gather


!*******************************************************************************
subroutine atm_lesgo_force(i)
!*******************************************************************************
!
! This will feed the velocity at all the actuator points into the atm. This is 
! done by using trilinear interpolation from lesgo. Force will be calculated 
! based on the velocities and stored on forceField
!
use param, only: rp, mpi_rp, z_i
use param, only: ld, ny, nz, ierr
use grid_defs, only: gridz 
use atm_input_util, only: turbineArray, turbineModel
use sim_param, only: u, v
use actuator_turbine_model, only: atm_computeBladeForce
use actuator_turbine_model, only: atm_computeNacelleForce
use mpi, only: mpi_sum
implicit none
integer, intent(in) :: i               ! The turbine number
integer :: m,q,j
real(rp), dimension(3) :: velocity
real(rp), dimension(3) :: mpi_velocity ! only used for Spalart method
real(rp), dimension(3) :: xyz          ! Point on which to interpolate velocity
integer, pointer :: TURBINE_COMM       ! The MPI turbine communcator

TURBINE_COMM => turbineArray(i) % TURBINE_COMM_WORLD

! The turbine type ID
j = turbineArray(i) % turbineTypeID 

if(turbineArray(i) % sampling == 'Spalart') then
  ! This loop goes through all the blade points and calculates the respective
  ! body forces then imposes it onto the force field
  ! Ref. Churchfield et al. 2017, AIAA 2017-1998
  do q=1,turbineArray(i) % numBladePoints
  do m=1,turbineModel(j) % numBl
    ! Actuator point onto which to interpolate the velocity
    xyz = turbineArray(i) % bladePoints(m,q,1:3)
    velocity = 0._rp
    mpi_velocity = 0._rp

    ! Smooth velocity
    ! Ref. Churchfield et al. 2017, AIAA 2017-1998
    call atm_lesgo_compute_spalart_u(i, xyz, velocity)

    mpi_velocity = velocity

    ! Sync all the blade forces
    call mpi_allreduce(mpi_velocity, velocity, size(velocity),                 &
                       mpi_rp, mpi_sum, TURBINE_COMM, ierr)

    ! This will compute the blade force for the specific point
    if((gridz(1) <= xyz(3)/z_i) .and. (xyz(3)/z_i < gridz(nz))) then
      call atm_computeBladeForce(i, m, q, velocity)
    else
      velocity = 0._rp
    endif
  enddo
  enddo

elseif(turbineArray(i) % sampling == 'atPoint') then
  ! This loop goes through all the blade points and calculates the respective
  ! body forces then imposes it onto the force field
  do q=1,turbineArray(i) % numBladePoints
  do m=1,turbineModel(j) % numBl
    ! Actuator point onto which to interpolate the velocity
    xyz = turbineArray(i) % bladePoints(m,q,1:3)
    ! Non-dimensionalizes the point location
    xyz = xyz/z_i

    ! Interpolate velocities if inside the domain, dimension [m/s]
    if((gridz(1) <= xyz(3)) .and. (xyz(3) < gridz(nz))) then
      ! Interpolate: flag_grid: 1 indicates uv grid
      call trilinear_interp_init(xyz, 1)
      call interpolation(u   (1:ld,1:ny,0:nz), velocity(1))
      call interpolation(v   (1:ld,1:ny,0:nz), velocity(2))
      call interpolation(w_uv(1:ld,1:ny,0:nz), velocity(3))

      ! This will compute the blade force for the specific point
      ! BladeForce, dimension [m^4/s^2] or [N] if density=1 [g/m^3]
      call atm_computeBladeForce(i, m, q, velocity)
    endif
  enddo
  enddo
endif

! Calculate Nacelle force
if(turbineArray(i) % nacelle) then
  xyz = turbineArray(i) % nacelleLocation
  xyz = xyz/z_i

  if((gridz(1) <= xyz(3)) .and. (xyz(3) < gridz(nz))) then
    ! Interpolate: flag_grid: 1 indicates uv grid
    call trilinear_interp_init(xyz, 1)
    call interpolation(u   (1:ld,1:ny,0:nz), velocity(1))
    call interpolation(v   (1:ld,1:ny,0:nz), velocity(2))
    call interpolation(w_uv(1:ld,1:ny,0:nz), velocity(3))

    call atm_computeNacelleForce(i, velocity)
  endif
endif

end subroutine atm_lesgo_force


!*******************************************************************************
subroutine atm_lesgo_compute_Spalart_u(i, xyz, velocity)
!*******************************************************************************
!
! This will calculate the sampling velocity using the proposed method
! from Spalart
! Ref. Churchfield et al. 2017, AIAA 2017-1998
! n turbine number
! xyz actuator point position vector
! velocity reference velocity for computing lift and drag
!
use param, only: rp, z_i, dx, dy, dz, pi, nz
use grid_defs, only: gridz 
use sim_param, only: u, v, w
use atm_input_util, only: turbineArray
implicit none
integer, intent(in) :: i
real(rp), intent(in) :: xyz(3)
real(rp), intent(inout) :: velocity(3)
integer :: c, m, n, q
real(rp) :: dist, a(3), projectradius, epsilon

! Value of epsilon
epsilon = turbineArray(i) % epsilon

! Projection radius
projectradius = turbineArray(i) % projectionRadius

! Set the velocity to zero
velocity = 0._rp

! uv-grid
do c=1,forceFieldUV(i) % c
  a = forceFieldUV(i) % location(1:3, c)
  m = forceFieldUV(i) % ijk(1, c)
  n = forceFieldUV(i) % ijk(2, c)
  q = forceFieldUV(i) % ijk(3, c)

  dist = ((a(1)-xyz(1))**2+(a(2)-xyz(2))**2+(a(3)-xyz(3))**2)**0.5_rp

  if (dist .le. (projectradius*z_i)) then
    if ((gridz(1) <= a(3)/z_i) .and. (a(3)/z_i < gridz(nz))) then

      ! The value of the kernel. This is the actual smoothing function
      velocity(1) = velocity(1) + u(m,n,q) * exp(-(dist/epsilon)**2)           &
                                /  ((epsilon**3._rp)*(pi**1.5_rp))
      velocity(2) = velocity(2) + v(m,n,q) * exp(-(dist/epsilon)**2)           &
                                /  ((epsilon**3._rp)*(pi**1.5_rp))
    endif
  endif
enddo

! W-grid (still on uv-grid)
do c=1,forceFieldW(i) % c
  a = forceFieldW(i) % location(1:3, c)
  m = forceFieldW(i) % ijk(1, c)
  n = forceFieldW(i) % ijk(2, c)
  q = forceFieldW(i) % ijk(3, c)

  dist = ((a(1)-xyz(1))**2+(a(2)-xyz(2))**2+(a(3)-xyz(3))**2)**0.5_rp

  if (dist .le. projectradius) then
    if ((gridz(1) <= a(3)/z_i) .and. (a(3)/z_i < gridz(nz))) then

      ! The value of the kernel. This is the actual smoothing function
      velocity(3) = velocity(3) + w(m,n,q) * exp(-(dist/epsilon)**2)           &
                                /  ((epsilon**3._rp)*(pi**1.5_rp))
    endif
  endif
enddo

velocity = velocity * z_i * dx * z_i * dy *z_i * dz

end subroutine atm_lesgo_compute_Spalart_u


!*******************************************************************************
subroutine atm_lesgo_convolute_force(i)
!*******************************************************************************
!
! This will convolute the forces for each turbine
! Ref. Churchfield et al. 2017 AIAA 2017-1998
!
use param, only: rp, pi, z_i
use atm_input_util, only: turbineArray, turbineModel
use sim_param, only: u, v, w
implicit none
integer, intent(in) :: i
integer :: j, m, q, c, mmend, qqend
integer :: ii, jj, kk 
real(rp) :: dist, a(3), b(3)
real(rp) :: nacelleProjectRiep2, projectRiep2
real(rp) :: nacelleEpsilon, invNacelleEpsilon2, invNacelleEp3pi1p5zi, epsilon, invepsilon2, invep3pi1p5zi
real(rp) :: kernel
real(rp), pointer, dimension(:,:,:) :: bladeForces, bladePoints
real(rp), pointer, dimension(:,:) :: bodyForceUV, bodyForceW

bladeForces => turbineArray(i) % bladeForces
bladePoints => turbineArray(i) % bladePoints
bodyForceUV => forceFieldUV(i) % force
bodyForceW => forceFieldW(i) % force

! The turbine type ID
j = turbineArray(i) % turbineTypeID 

! This will convolute the blade force onto the grid points
! affected by the turbines on both grids
! Only if the distance is less than specified value
mmend = turbineModel(j) % numBl
qqend = turbineArray(i) % numBladePoints
nacelleEpsilon = turbineArray(i) % nacelleEpsilon
epsilon = turbineArray(i) % epsilon
! For code performance
invepsilon2=1.0_rp/(epsilon**2)
invep3pi1p5zi=1.0_rp/(epsilon**3 * pi**1.5_rp)*z_i
projectRiep2 = turbineArray(i) % projectionRadius **2 *invepsilon2
invNacelleEpsilon2=1.0_rp/(nacelleEpsilon**2)
invNacelleEp3pi1p5zi=1.0_rp/(nacelleEpsilon**3 * pi**1.5_rp)*z_i
nacelleProjectRiep2 = turbineArray(i) % projectionRadius **2 *invNacelleEpsilon2

! Body Force implementation using velocity sampling at the actuator point
if (turbineArray(i) % sampling == 'atPoint') then
  do c=1, forceFieldUV(i) % c
    a = forceFieldUV(i) % location(1:3,c)
    bodyForceUV(1:2,c) = 0._rp

    ! Blade forces
    do m=1, mmend
    do q=1, qqend
      b = bladePoints(m,q,:)
      dist = (dot_product(a-b,a-b))*invepsilon2
      if(dist .le. projectRiep2) then
        ! The value of the kernel. This is the actual smoothing function
        ! Basically epsilon=2.5*dx with dimension [m]
        ! BladeForces, dimension [m^4/s^2] or [N] if density=1 [g/m^3]
        kernel = exp(-dist) * invep3pi1p5zi
        bodyForceUV(1:2,c) = bodyForceUV(1:2,c) + bladeForces(m,q,1:2) * kernel
      endif
    enddo
    enddo

    ! Nacelle force
    if(turbineArray(i) % nacelle) then
      b = turbineArray(i) % nacelleLocation
      dist = (dot_product(a-b,a-b))*invNacelleEpsilon2
      if(dist .le. nacelleProjectRiep2) then
        ! The value of the kernel. This is the actual smoothing function
        kernel = exp(-dist) * invNacelleEp3pi1p5zi
        bodyForceUV(1:2,c) = bodyForceUV(1:2,c) + turbineArray(i) % nacelleForce(1:2) * kernel
      endif
    endif

  enddo

  do c=1, forceFieldW(i) % c
    a = forceFieldW(i) %  location(1:3,c)
    bodyForceW(3,c) = 0._rp

    ! Blade forces
    do m=1,mmend
    do q=1,qqend
      b = bladePoints(m,q,:)
      dist = (dot_product(a-b,a-b))*invepsilon2
      if(dist .le. projectRiep2) then
        ! The value of the kernel. This is the actual smoothing function
        kernel = exp(-dist) * invep3pi1p5zi
        bodyForceW(3,c) = bodyForceW(3,c) + bladeForces(m,q,3) * kernel
      endif
    enddo
    enddo

    ! Nacelle force
    if(turbineArray(i) % nacelle) then
      b = turbineArray(i) % nacelleLocation
      dist = (dot_product(a-b,a-b))*invNacelleEpsilon2
      if(dist .le. nacelleProjectRiep2) then
        ! The value of the kernel. This is the actual smoothing function
        kernel = exp(-dist) * invNacelleEp3pi1p5zi
        bodyForceW(3,c) = bodyForceW(3,c) + turbineArray(i) % nacelleForce(3) * kernel
      endif
    endif

  enddo

! The Spalart method uses the local velocity field.
! For this reason it needs to be done explicitly in this module
! and cannot be generally coded from the actuator_turbine_model module
elseif (turbineArray(i) % sampling == 'Spalart') then
  do c=1,forceFieldUV(i) % c
    a = forceFieldUV(i) % location(1:3,c)
    bodyForceUV(1:2,c) = 0._rp

    ! Indices for velocity field
    ii = forceFieldUV(i) % ijk(1,c)
    jj = forceFieldUV(i) % ijk(2,c)
    kk = forceFieldUV(i) % ijk(3,c)

    ! Blade forces
    do m=1, mmend
    do q=1, qqend
      b = bladePoints(m,q,:)
      dist = (dot_product(a-b,a-b))*invepsilon2
      if(dist .le. projectRiep2) then
        ! The value of the kernel.
        ! This is the actual smoothing function
        ! Divide by velocity magnitude
        kernel = exp(-dist) * invep3pi1p5zi
        bodyForceUV(1,c) = bodyForceUV(1,c) + bladeForces(m,q,1) * kernel                      &
                 / turbineArray(i) % Vmag(m,q) * (u(ii,jj,kk) +                &
                   turbineArray(i) % rotSpeed *                                &
                   turbineArray(i) % bladeRadius(m,q) *                        &
                   cos(turbineModel(j) % PreCone))

        bodyForceUV(2,c) = bodyForceUV(2,c) + bladeForces(m,q,2) * kernel                      &
                 / turbineArray(i) % Vmag(m,q) * (v(ii,jj,kk) +                &
                   turbineArray(i) % bladeAlignedVectors(m,q,2,2) *            &
                   turbineArray(i) % rotSpeed *                                &
                   turbineArray(i) % bladeRadius(m,q) *                        &
                   cos(turbineModel(j) % PreCone))
      endif
    enddo
    enddo

    ! Nacelle force
    if(turbineArray(i) % nacelle) then
      b = turbineArray(i) % nacelleLocation
      dist = (dot_product(a-b,a-b))*invNacelleEpsilon2
      if(dist .le. nacelleProjectRiep2) then
        ! The value of the kernel. This is the actual smoothing function
        kernel = exp(-dist) * invNacelleEp3pi1p5zi
        bodyForceUV(1:2,c) = bodyForceUV(1:2,c) + turbineArray(i) % nacelleForce(1:2) * kernel
      endif
    endif

  enddo

  do c=1, forceFieldW(i) % c
    a = forceFieldW(i) %  location(1:3,c)
    bodyForceW(3,c) = 0._rp

    ! Indices for velocity field
    ii = forceFieldW(i) % ijk(1,c)
    jj = forceFieldW(i) % ijk(2,c)
    kk = forceFieldW(i) % ijk(3,c)

    ! Blade forces
    do m=1,mmend
    do q=1,qqend
      b= bladePoints(m,q,:)
      dist = (dot_product(a-b,a-b))*invepsilon2
      if(dist .le. projectRiep2) then
        ! The value of the kernel. This is the actual smoothing function
        kernel = exp(-dist) * invep3pi1p5zi
        bodyForceW(3,c) = bodyForceW(3,c) + bladeForces(m,q,3) * kernel                      &
                 / turbineArray(i) % Vmag(m,q) * (w(ii,jj,kk) +                &
                   turbineArray(i) % bladeAlignedVectors(m,q,2,3) *            &
                   turbineArray(i) % rotSpeed *                                &
                   turbineArray(i) % bladeRadius(m,q) *                        &
                   cos(turbineModel(j) % PreCone))
      endif
    enddo
    enddo

    ! Nacelle force
    if(turbineArray(i) % nacelle) then
      b = turbineArray(i) % nacelleLocation
      dist = (dot_product(a-b,a-b))*invNacelleEpsilon2
      if(dist .le. nacelleProjectRiep2) then
        ! The value of the kernel. This is the actual smoothing function
        kernel = exp(-dist) * invNacelleEp3pi1p5zi
        bodyForceW(3,c) = bodyForceW(3,c) + turbineArray(i) % nacelleForce(3) * kernel
      endif
    endif

  enddo
endif

end subroutine atm_lesgo_convolute_force


!*******************************************************************************
subroutine atm_lesgo_apply_force
!*******************************************************************************
!
! This will apply the blade force onto the CFD grid by using the convolution
! function in the ATM library
!
use atm_input_util, only: turbineArray
use atm_input_util, only: numberOfTurbines
use sim_param, only: fxa, fya, fza 
implicit none
integer :: c,m
integer :: i,j,k

do m=1,numberOfTurbines
  if(turbineArray(m) % operate) then
    ! Impose force field onto the flow field variables
    ! The forces are non-dimensionalized here as well
    do c=1, forceFieldUV(m) % c
      i = forceFieldUV(m) % ijk(1,c)
      j = forceFieldUV(m) % ijk(2,c)
      k = forceFieldUV(m) % ijk(3,c)
      fxa(i,j,k) = fxa(i,j,k) + forceFieldUV(m) % force(1,c)
      fya(i,j,k) = fya(i,j,k) + forceFieldUV(m) % force(2,c)
    enddo

    do c=1, forceFieldW(m) % c
      i = forceFieldW(m) % ijk(1,c)
      j = forceFieldW(m) % ijk(2,c)
      k = forceFieldW(m) % ijk(3,c)
      fza(i,j,k) = fza(i,j,k) + forceFieldW(m) % force(3,c)
      enddo
    endif
enddo

end subroutine atm_lesgo_apply_force

end module atm_lesgo_interface
#endif ATM
