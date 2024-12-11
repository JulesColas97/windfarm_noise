#ifdef ATM
!*******************************************************************************
module actuator_turbine_model
!*******************************************************************************
!
! This module has the subroutines to provide all calculations for use in the 
! actuator turbine model (ATM)
!
use param, only: rp, pi
implicit none
save
private
public :: atm_computeBladeForce
public :: atm_computeNacelleForce
public :: atm_computeRotorSpeed
public :: atm_control_yaw
public :: atm_initialize
public :: atm_initialize_output
public :: atm_output
public :: atm_rotateBlades
public :: atm_write_restart

! Establishes if we are at the first time step
logical :: pastFirstTimeStep 
! Degrees to radians conversion
real(rp) :: degRad = pi/180._rp
! Set the revolutions/min to radians/s 
real(rp) :: rpmRadSec = pi/30._rp

contains

!*******************************************************************************
subroutine atm_initialize
!*******************************************************************************
!
! This subroutine initializes the ATM. It calls the subroutines in
! atm_input_util to read the input data and creates the initial geometry
!
use param, only: coord, initu, path
use atm_input_util, only: atm_read_input_conf
use atm_input_util, only: numberOfTurbines
use atm_input_util, only: turbineArray
implicit none
integer :: i
logical :: file_exists

! The first time step not reached yet
pastFirstTimeStep = .false. 

if(coord.eq.0) write(*,*) 'Reading Actuator Turbine Model Input...'
! Read input data
call atm_read_input_conf
if(coord.eq.0) write(*,*) 'Done Reading Actuator Turbine Model Input'

do i=1,numberOfTurbines
  ! Creates the ATM points defining the geometry
  call atm_create_points(i) 

  ! This will create the first yaw alignment
  turbineArray(i) % deltaNacYaw = turbineArray(i) % nacYaw
  call atm_yawNacelle(i)

  ! Calculates variables depending on input
  call atm_calculate_variables(i) 

  if(initu) then
    inquire(file = trim(path)//"turbineOutput/"//trim(turbineArray(i) % turbineName)//   &
                   "/actuatorPoints", exist=file_exists)
    if(file_exists) then
      if(coord.eq.0) print*, 'Reading bladePoints from Previous Simulation'
      call atm_read_actuator_points(i)
    else
      if(coord.eq.0) print*, 'BladePoints are not present: spinup simulation, starting from beginning'
    endif

    inquire(file = trim(path)//"turbineOutput/"//trim(turbineArray(i) % turbineName)//   &
                   "/restart", exist=file_exists)
    if(file_exists) then
      if(coord.eq.0) print*, 'Reading Turbine Properties from Previous Simulation'
      call atm_read_restart(i)
    else
      if(coord.eq.0) print*, 'Turbine properties are not present: spinup simulation, starting from beginning'
    endif
  endif
end do

! Past the first time step
pastFirstTimeStep = .true.

end subroutine atm_initialize


!*******************************************************************************
subroutine atm_read_actuator_points(i)
!*******************************************************************************
!
! This subroutine reads the location of the actuator points
! It is used if the simulation wants to start from a previous simulation
! without having to start the turbine from the original location
!
use param, only: path 
use atm_input_util, only: turbineArray
use atm_input_util, only: turbineModel
integer, intent(in) :: i ! Indicates the turbine number
integer :: j, m, q

! The turbine type ID
j = turbineArray(i) % turbineTypeID

open(unit=1, file=trim(path)//"turbineOutput/"//trim(turbineArray(i) % turbineName)//    &
                  "/actuatorPoints", action='read')

do m=1, turbineModel(j) % numBl
do q=1, turbineArray(i) % numBladePoints
  read(1,*) turbineArray(i) % bladePoints(m,q,:)
enddo
enddo

close(1)

end subroutine atm_read_actuator_points


!*******************************************************************************
subroutine atm_read_restart(i)
!*******************************************************************************
!
! This subroutine reads the rotor speed
! It is used if the simulation wants to start from a previous simulation
! without having to start the turbine from the original omega
!
use param, only: coord, path
use atm_input_util, only: turbineArray
integer, intent(in) :: i  ! Indicates the turbine number

! Open the file at the last line (append)
open(unit=1, file=trim(path)//"turbineOutput/"//trim(turbineArray(i) % turbineName)//    &
                  "/restart", action='read')

! Read the restart variables
read(1,*)
read(1,*) turbineArray(i) % rotSpeed
read(1,*)
read(1,*) turbineArray(i) % torqueGen
read(1,*)
read(1,*) turbineArray(i) % torqueRotor
read(1,*)
read(1,*) turbineArray(i) % u_infinity
read(1,*)
read(1,*) turbineArray(i) % induction_a
read(1,*)
read(1,*) turbineArray(i) % PitchControlAngle
read(1,*)
read(1,*) turbineArray(i) % IntSpeedError
read(1,*)
read(1,*) turbineArray(i) % nacYaw
read(1,*)
read(1,*) turbineArray(i) % rotorApex
read(1,*)
read(1,*) turbineArray(i) % uvShaft 

close(1)

if(coord .eq. 0) then
  write(*,*) ' RotSpeed Value from previous simulation is ',                   &
               turbineArray(i) % rotSpeed
  write(*,*) ' torqueGen Value from previous simulation is ',                  &
               turbineArray(i) % torqueGen
  write(*,*) ' torqueRotor Value from previous simulation is ',                &
               turbineArray(i) % torqueRotor
  write(*,*) ' PitchControlAngle Value from previous simulation is ',          &
               turbineArray(i) % PitchControlAngle
  write(*,*) ' IntSpeedError Value from previous simulation is ',              &
               turbineArray(i) % IntSpeedError
  write(*,*) ' Yaw Value from previous simulation is ',                        &
               turbineArray(i) % nacYaw
  write(*,*) ' Rotor Apex Value from previous simulation is ',                 &
               turbineArray(i) % rotorApex
  write(*,*) ' uvShaft Value from previous simulation is ',                    &
               turbineArray(i) % uvShaft
endif

end subroutine atm_read_restart


!*******************************************************************************
subroutine atm_write_restart(i)
!*******************************************************************************
!
! This subroutine writes the rotor speed
! It is used if the simulation wants to start from a previous simulation
! without having to start the turbine from the original omega
!
use param, only : path 
use atm_input_util, only: turbineArray
use atm_input_util, only: turbineModel
implicit none
integer, intent(in) :: i  ! Indicates the turbine number
integer :: pointsFile=787 ! File to write the actuator points
integer :: restartFile=21 ! File to write restart data
integer j,m,q             ! counters

! Open the file 
open(unit=restartFile, file=trim(path)//"turbineOutput/"//                               &
     trim(turbineArray(i) % turbineName)//"/restart", status="replace")

! Store the rotSpeed value 
write(restartFile,*) 'RotSpeed'
write(restartFile,*) turbineArray(i) % rotSpeed
write(restartFile,*) 'torqueGen'
write(restartFile,*) turbineArray(i) % torqueGen
write(restartFile,*) 'torqueRotor'
write(restartFile,*) turbineArray(i) % torqueRotor
write(restartFile,*) 'u_infinity'
write(restartFile,*) turbineArray(i) % u_infinity
write(restartFile,*) 'induction_a'
write(restartFile,*) turbineArray(i) % induction_a
write(restartFile,*) 'PitchControlAngle'
write(restartFile,*) turbineArray(i) % PitchControlAngle
write(restartFile,*) 'IntSpeedError'
write(restartFile,*) turbineArray(i) % IntSpeedError
write(restartFile,*) 'nacYaw'
write(restartFile,*) turbineArray(i) % nacYaw
write(restartFile,*) 'rotorApex'
write(restartFile,*) turbineArray(i) % rotorApex
write(restartFile,*) 'uvShaft'
write(restartFile,*) turbineArray(i) % uvShaft

close(restartFile)

! Write the actuator points at every time-step regardless
! The turbine type ID
j = turbineArray(i) % turbineTypeID

open(unit=pointsFile, status="replace", file=trim(path)//"turbineOutput/"//              &
     trim(turbineArray(i) % turbineName)//"/actuatorPoints")

do m=1, turbineModel(j) % numBl
do q=1, turbineArray(i) % numBladePoints
  ! A new file will be created each time-step with the proper
  ! location of the blades
  write(pointsFile,*) turbineArray(i) % bladePoints(m,q,:)
enddo
enddo

close(pointsFile)

end subroutine atm_write_restart


!*******************************************************************************
subroutine atm_initialize_output
!*******************************************************************************
!
! This subroutine initializes the output files for the ATM
!
use param, only: path 
use atm_input_util, only: numberOfTurbines
use atm_input_util, only: turbineArray
implicit none
integer :: i
logical :: file_exists

! Write to the screen output start
write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
write(*,*) 'Actuator Turbine Model has been implemented with',                 &
              numberOfTurbines ,'turbines'
write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

do i=1,numberOfTurbines

  inquire(file = trim(path)//"turbineOutput/"//trim(turbineArray(i) % turbineName)//   &
                   "/actuatorPoints", exist=file_exists)

  if(.not. file_exists) then
    ! Create turbineOutput directory    
    call system('mkdir -p '//trim(path)//'turbineOutput/'//trim(turbineArray(i) % turbineName))

    open(unit=1, file=trim(path)//"turbineOutput/"//trim(turbineArray(i) % turbineName)//"/power") 
    write(1,*) 'time PowerRotor powerGen'
    close(1)

    open(unit=1, file=trim(path)//"turbineOutput/"//trim(turbineArray(i) % turbineName)//"/thrust") 
    write(1,*) 'time thrust'
    close(1)

    open(unit=1, file=trim(path)//"turbineOutput/"//trim(turbineArray(i) % turbineName)//"/RotSpeed") 
    write(1,*) 'time RotSpeed'
    close(1)
        
    open(unit=1, file=trim(path)//"turbineOutput/"//trim(turbineArray(i) % turbineName)//"/Yaw") 
    write(1,*) 'time deltaNacYaw NacYaw'
    close(1)
        
    open(unit=1, file=trim(path)//"turbineOutput/"//trim(turbineArray(i) % turbineName)//"/lift")
    write(1,*) 'turbineNumber bladeNumber lift'
    close(1)

    open(unit=1, file=trim(path)//"turbineOutput/"//trim(turbineArray(i) % turbineName)//"/drag")
    write(1,*) 'turbineNumber bladeNumber drag'
    close(1)
        
    open(unit=1, file=trim(path)//"turbineOutput/"//trim(turbineArray(i) % turbineName)//"/Cl")
    write(1,*) 'turbineNumber bladeNumber Cl'
    close(1)

    open(unit=1, file=trim(path)//"turbineOutput/"//trim(turbineArray(i) % turbineName)//"/Cd")
    write(1,*) 'turbineNumber bladeNumber Cd'
    close(1)

    open(unit=1, file=trim(path)//"turbineOutput/"//trim(turbineArray(i) % turbineName)//"/alpha")
    write(1,*) 'turbineNumber bladeNumber alpha'
    close(1)

    open(unit=1, file=trim(path)//"turbineOutput/"//trim(turbineArray(i) % turbineName)//"/Vrel")
    write(1,*) 'turbineNumber bladeNumber Vrel'
    close(1)

    open(unit=1, file=trim(path)//"turbineOutput/"//trim(turbineArray(i) % turbineName)//"/Vaxial")
    write(1,*) 'turbineNumber bladeNumber Vaxial'
    close(1)

    open(unit=1, file=trim(path)//"turbineOutput/"//trim(turbineArray(i) % turbineName)//"/Vtangential")
    write(1,*) 'turbineNumber bladeNumber Vtangential'
    close(1)

    open(unit=1, file=trim(path)//"turbineOutput/"//trim(turbineArray(i) % turbineName)//"/tangentialForce")
    write(1,*) 'turbineNumber bladeNumber tangentialForce'
    close(1)

    open(unit=1, file=trim(path)//"turbineOutput/"//trim(turbineArray(i) % turbineName)//"/axialForce")
    write(1,*) 'turbineNumber bladeNumber axialForce'
    close(1)

    open(unit=1, file=trim(path)//"turbineOutput/"//trim(turbineArray(i) % turbineName)//"/nacelle")
    write(1,*) 'time Velocity-no-correction Velocity-w-correction'
    close(1)

    ! Writes the ATM points defining the geometry
    call atm_write_restart(i) 
  endif
enddo

end subroutine atm_initialize_output


!*******************************************************************************
subroutine atm_create_points(i)
!*******************************************************************************
!
! This subroutine generate the set of blade points for each turbine
!
use sim_param, only: u_dim
use param, only: rp
use atm_input_util, only: turbineArray
use atm_input_util, only: turbineModel
use atm_base, only: rotatePoint
implicit none
integer, intent(in) :: i       ! Indicates the turbine number
integer :: j                   ! Indicates the turbine type
integer :: m                   ! Indicates the blade point number
integer :: k
real(rp), dimension(3) :: root ! Location of rotor apex
real(rp) :: beta               ! angle between vertical direction and blade
real(rp) :: dist               ! Distance from each actuator point
integer,  pointer :: numBladePoints
integer,  pointer :: numBl
real(rp), pointer :: NacYaw             
real(rp), pointer :: db(:)     ! Distance between two actuator points
real(rp), pointer :: bladePoints(:,:,:)
real(rp), pointer :: bladeRadius(:,:)
real(rp), pointer :: azimuth
real(rp), pointer :: rotSpeed
real(rp), pointer :: ShftTilt
real(rp), pointer :: PreCone
real(rp), pointer :: towerShaftIntersect(:)
real(rp), pointer :: baseLocation(:)
real(rp), pointer :: TowerHt
real(rp), pointer :: Twr2Shft
real(rp), pointer :: rotorApex(:)
real(rp), pointer :: OverHang
real(rp), pointer :: UndSling
real(rp), pointer :: uvShaftDir
real(rp), pointer :: uvShaft(:)
real(rp), pointer :: uvTower(:)
real(rp), pointer :: TipRad
real(rp), pointer :: HubRad

! Identifies the turbineModel being used
! The type of turbine (eg. NREL5MW)
j = turbineArray(i) % turbineTypeID 

! Variables to be used locally. They are stored in local variables within the 
! subroutine for easier code following. The values are then passed to the 
! proper type
numBladePoints => turbineArray(i) % numBladePoints
numBl => turbineModel(j) % numBl

! Allocate variables depending on specific turbine properties and general
! turbine model properties
allocate(turbineArray(i) % db(numBladePoints))
allocate(turbineArray(i) % bladePoints(numBl,              &
         numBladePoints,3))
allocate(turbineArray(i) % bladeRadius(numBl,numBladePoints))  

! Assign Pointers turbineArray denpendent (i)
db => turbineArray(i) % db
bladePoints => turbineArray(i) % bladePoints
bladeRadius => turbineArray(i) % bladeRadius
azimuth => turbineArray(i) % azimuth
rotSpeed => turbineArray(i) % rotSpeed
towerShaftIntersect => turbineArray(i) % towerShaftIntersect
baseLocation => turbineArray(i) % baseLocation
uvShaft => turbineArray(i) % uvShaft
uvTower => turbineArray(i) % uvTower
rotorApex => turbineArray(i) % rotorApex
uvShaftDir => turbineArray(i) % uvShaftDir
nacYaw => turbineArray(i) % nacYaw

! Assign Pointers turbineModel (j)
ShftTilt => turbineModel(j) % ShftTilt
preCone  => turbineModel(j) % preCone
TowerHt  => turbineModel(j) % TowerHt
Twr2Shft => turbineModel(j) % Twr2Shft
OverHang => turbineModel(j) % OverHang
UndSling => turbineModel(j) % UndSling
TipRad   => turbineModel(j) % TipRad
HubRad   => turbineModel(j) % HubRad
PreCone  => turbineModel(j) % PreCone

! Turbine specific
azimuth  = degRad * azimuth
rotSpeed = rpmRadSec * rotSpeed / u_dim
nacYaw   = degRad * nacYaw

! Turbine model specific
shftTilt = degRad * shftTilt 
preCone  = degRad * preCone

! Calculate tower shaft intersection and rotor apex locations. (The i-index is 
! at the turbine array level for each turbine and the j-index is for each type 
! of turbine--if all turbines are the same, j- is always 0.)  The rotor apex is
! not yet rotated for initial yaw that is done below.
towerShaftIntersect = turbineArray(i) % baseLocation
towerShaftIntersect(3) = towerShaftIntersect(3) + TowerHt + Twr2Shft
rotorApex = towerShaftIntersect
rotorApex(1) = rotorApex(1) + (OverHang + UndSling) * cos(ShftTilt)
rotorApex(3) = rotorApex(3) + (OverHang + UndSling) * sin(ShftTilt)

! Create Nacelle Point
turbineArray(i) % nacelleLocation = rotorApex

! Create the first set of actuator points
! Define the vector along the shaft pointing in the direction of the wind
uvShaftDir = OverHang / abs(OverHang)

! Define the vector along the shaft pointing in the direction of the wind
uvShaft = rotorApex - towerShaftIntersect
uvShaft = uvShaft / sqrt(dot_product(uvShaft,uvShaft))
uvShaft = uvShaft * uvShaftDir

! Define vector aligned with the tower pointing from the ground to the nacelle
uvTower = towerShaftIntersect - baseLocation
uvTower = uvTower / sqrt(dot_product(uvTower,uvTower))

! Define thickness of each blade section
do k=1,numBladePoints
  db(k) = (TipRad - HubRad) / real(numBladePoints)
enddo

! This creates the first set of points
do k=1,numBl
  root = rotorApex
  beta = PreCone - ShftTilt
  root(1) = root(1) + HubRad*sin(beta)
  root(3) = root(3) + HubRad*cos(beta)
  dist = 0._rp
    
  ! Number of blade points
  do m=1,numBladePoints
    dist = dist + 0.5_rp*db(m)

    bladePoints(k,m,1) = root(1) + dist*sin(beta)
    bladePoints(k,m,2) = root(2)
    bladePoints(k,m,3) = root(3) + dist*cos(beta)
    bladeRadius(k,m) = dist + HubRad

    dist = dist + 0.5_rp*db(m)
  enddo

  ! If there are more than one blade create the points of other blades by
  ! rotating the points of the first blade
  if(k > 1) then
    do m=1,numBladePoints
      call rotatePoint(bladePoints(k,m,:), rotorApex, uvShaft,               &
           360._rp/real(NumBl)*(k-1._rp)*degRad)
    enddo
  endif

enddo

! Apply the first rotation
turbineArray(i) % deltaAzimuth = azimuth

call atm_rotateBlades(i)

end subroutine atm_create_points


!*******************************************************************************
subroutine atm_control_yaw(i, time)
!*******************************************************************************
!
! This subroutine updates the model each time-step
!
use param, only: rp
use atm_input_util, only: turbineArray
use atm_input_util, only: turbineModel
use atm_base, only: interpolate
implicit none
integer, intent(in) :: i     ! Turbine number
real(rp), intent(in) :: time ! Simulation time, dimension [second]
integer :: j                 ! Turbine Type ID

! Identifies the turbineModel being used
! The type of turbine (eg. NREL5MW)
j = turbineArray(i) % turbineTypeID 

! Will calculate the yaw angle and yaw the nacelle (from degrees to radians)
if (turbineModel(j) % YawControllerType == "timeYawTable" ) then
  turbineArray(i) % deltaNacYaw = interpolate(time,                            &
  turbineModel(j) % yaw_time(:), turbineModel(j) % yaw_angle(:)) * degRad -    &
  turbineArray(i) % NacYaw

  ! Yaw only if angle is greater than given tolerance
  if (abs(turbineArray(i) % deltaNacYaw) > 0.00000001_rp) then
    call atm_yawNacelle(i)
  endif
endif

end subroutine atm_control_yaw


!*******************************************************************************
subroutine atm_computeRotorSpeed(i, dt)
!*******************************************************************************
!
! This subroutine computes the deltaAzimuth used to rotate blades
!
use param, only: rp
use atm_input_util, only: turbineArray
use atm_input_util, only: turbineModel
implicit none
integer, intent(in) :: i      ! Turbine number
real(rp), intent(in) :: dt    ! time step, dimension [second]
integer :: j                  ! Turbine type
! Pointers to turbineArray(i)
real(rp), pointer :: rotSpeed, torqueGen, torqueRotor
! Pointers to turbineModel(j)
real(rp), pointer :: GBRatio, CutInGenSpeed, RatedGenSpeed
real(rp), pointer :: Region2StartGenSpeed, Region2EndGenSpeed
real(rp), pointer :: CutInGenTorque, RateLimitGenTorque, RatedGenTorque
real(rp), pointer :: KGen,TorqueControllerRelax, DriveTrainIner
! Other variables to be used
real(rp) :: torqueGenOld, genSpeed, dGenSpeed, Region2StartGenTorque
real(rp) :: torqueSlope, Region2EndGenTorque
! Pitch Controller values
real(rp) :: GK, KP, KI, SpeedError
real(rp), pointer :: IntSpeedError, PitchControlAngle

j = turbineArray(i) % turbineTypeID

rotSpeed     => turbineArray(i) % rotSpeed
torqueGen    => turbineArray(i) % torqueGen
torqueRotor  => turbineArray(i) % torqueRotor

GBRatio => turbineModel(j) % GBRatio
CutInGenSpeed => turbineModel(j) % CutInGenSpeed
CutInGenTorque => turbineModel(j) % CutInGenTorque
Region2StartGenSpeed => turbineModel(j) % Region2StartGenSpeed
KGen => turbineModel(j) % KGen
RatedGenSpeed => turbineModel(j) % RatedGenSpeed
Region2EndGenSpeed => turbineModel(j) % Region2EndGenSpeed
RatedGenTorque => turbineModel(j) % RatedGenTorque
RateLimitGenTorque => turbineModel(j) % RateLimitGenTorque
TorqueControllerRelax => turbineModel(j) % TorqueControllerRelax
DriveTrainIner => turbineModel(j) % DriveTrainIner

IntSpeedError => turbineArray(i) % IntSpeedError
PitchControlAngle => turbineArray(i) % PitchControlAngle

! Torque controllers (If there's no torque controller, then don't do anything)
! Ref. Jonkman et al. (2009, NREL)
if(turbineModel(j) % TorqueControllerType == "fiveRegion") then

  ! Get the generator speed.
  genSpeed = (rotSpeed/rpmRadSec)*GBRatio

  ! Save the generator torque from the last time step.
  torqueGenOld = torqueGen
            
  ! Region 1.
  if(genSpeed < CutInGenSpeed) then
    torqueGen = CutInGenTorque

  ! Region 1-1/2.
  elseif((genSpeed >= CutInGenSpeed) .and.                                     &
         (genSpeed < Region2StartGenSpeed)) then
    dGenSpeed = genSpeed - CutInGenSpeed
    Region2StartGenTorque = KGen * Region2StartGenSpeed *                      &
                                   Region2StartGenSpeed
    torqueSlope = (Region2StartGenTorque - CutInGenTorque) /                   &
                  (Region2StartGenSpeed  - CutInGenSpeed )
    torqueGen = CutInGenTorque + torqueSlope*dGenSpeed

  ! Region 2.
  elseif((genSpeed >= Region2StartGenSpeed) .and.                              &
         (genSpeed < Region2EndGenSpeed)) then
    torqueGen = KGen * genSpeed * genSpeed

  ! Region 2-1/2.
  elseif((genSpeed >= Region2EndGenSpeed) .and.                                &
          (genSpeed < RatedGenSpeed)) then
    dGenSpeed = genSpeed - Region2EndGenSpeed
    Region2EndGenTorque = KGen * Region2EndGenSpeed *                          &
                                 Region2EndGenSpeed
    torqueSlope = (RatedGenTorque - Region2EndGenTorque) /                     &
                  (RatedGenSpeed  - Region2EndGenSpeed )
    torqueGen = Region2EndGenTorque + torqueSlope*dGenSpeed

  ! Region 3.
  elseif(genSpeed >= RatedGenSpeed) then
    torqueGen = RatedGenTorque
  endif

  ! Limit the change in generator torque if after first time step
  ! (otherwise it slowly ramps up from its zero initialized value--we
  ! want it to instantly be at its desired value on the first time
  ! step, but smoothly vary from there).
  if((abs((torqueGen - torqueGenOld)/dt) > RateLimitGenTorque)                 &
      .and. pastFirstTimeStep) then
    if(torqueGen > torqueGenOld) then
      torqueGen = torqueGenOld + (RateLimitGenTorque * dt)
    elseif(torqueGen <= torqueGenOld) then
      torqueGen = torqueGenOld - (RateLimitGenTorque * dt)
    endif
  endif

  ! Update the rotor speed.
  rotSpeed = rotSpeed + TorqueControllerRelax * (dt/DriveTrainIner) *          &
            (torqueRotor - GBRatio*torqueGen)

  if(turbineModel(j) % PitchControllerType == "none") then
    ! Limit the rotor speed to be positive and such that the generator 
    ! does not turn faster than rated.
    rotSpeed = max(0._rp, rotSpeed)
    rotSpeed = min(rotSpeed,(RatedGenSpeed*rpmRadSec)/GBRatio)
  endif

! Torque control for fixed tip speed ratio
! Note that this current method does NOT support Coning in the rotor
elseif(turbineModel(j) % TorqueControllerType == "fixedTSR") then

  if(pastFirstTimeStep) then
    ! Integrate the velocity along all actuator points
    call atm_integrate_u(i)   
    
    ! Match the rotor speed to a given TSR
    rotSpeed = turbineArray(i) % u_infinity_mean *                             &
               turbineArray(i) % TSR / turbineModel(j) % tipRad

    ! Important to get rid of negative values
    rotSpeed = max(0._rp, rotSpeed)
  endif
endif

! calculate u_infinity_mean, used in tip-loss correction
if(turbineArray(i) % tipALMCorrection) then
  if(turbineModel(j) % TorqueControllerType .ne. "fixedTSR") then
    call atm_integrate_u(i)   
  endif
endif

! Pitch controllers (If there's no pitch controller, then don't do anything)
! update PitchControlAngle
! Ref. Jonkman et al. (2009, NREL)
if(turbineModel(j) % PitchControllerType == "gainScheduledPI") then

  ! Get the generator speed.
  genSpeed = (rotSpeed/rpmRadSec)*GBRatio

  ! Calculate the gain
  GK =  1._rp/(1._rp + PitchControlAngle/turbineModel(j) % PitchControlAngleK)

  ! Calculate the Proportional and Integral terms
  KP = GK*turbineModel(j) % PitchControlKP0
  KI = GK*turbineModel(j) % PitchControlKI0

  ! Get speed error (generator in rpm) and update integral
  ! Integral is saturated to not push the angle beyond its limits
  SpeedError = genSpeed - RatedGenSpeed

  IntSpeedError = IntSpeedError + SpeedError*dt
  IntSpeedError = min( max(IntSpeedError,                                      &
                  turbineModel(j) % PitchControlAngleMin/KI),                  &
                  turbineModel(j) % PitchControlAngleMax/KI)

  ! Apply PI controller and saturate
  PitchControlAngle = KP*SpeedError + KI*IntSpeedError
  PitchControlAngle = min( max( PitchControlAngle,                             &
                      turbinemodel(j) % PitchControlAngleMin),                 &
                      turbineModel(j) % PitchControlAngleMax)
endif

! Compute the change in blade position at new rotor speed.
turbineArray(i) % deltaAzimuth = rotSpeed * dt

end subroutine atm_computeRotorSpeed


!*******************************************************************************
subroutine atm_rotateBlades(i)
!*******************************************************************************
!
! This subroutine rotates the turbine blades 
!
use param, only: rp, pi
use atm_input_util, only: turbineArray
use atm_input_util, only: turbineModel
use atm_base, only: rotatePoint
implicit none
integer, intent(in) :: i                  ! Turbine number
integer :: j                              ! Turbine type
integer :: m, q                           ! Counters to be used in do loops
real(rp) :: deltaAzimuth, deltaAzimuthI   ! Angle of rotation
real(rp), pointer :: rotorApex(:)
real(rp), pointer :: rotSpeed
real(rp), pointer :: uvShaft(:)
real(rp), pointer :: azimuth

j = turbineArray(i) % turbineTypeID

! Variables which are used by pointers
rotorApex => turbineArray(i) % rotorApex
rotSpeed  => turbineArray(i) % rotSpeed
uvShaft   => turbineArray(i) % uvShaft
azimuth   => turbineArray(i) % azimuth

! Angle of rotation
deltaAzimuth = turbineArray(i) % deltaAzimuth

! Check the rotation direction first and set the local delta azimuth
! variable accordingly.
if(turbineArray(i) % rotationDir == "cw") then ! clockwise
  deltaAzimuthI = deltaAzimuth
elseif(turbineArray(i) % rotationDir == "ccw") then ! counter-clockwisee
  deltaAzimuthI = -deltaAzimuth
endif

! Loop through all the points and rotate them accordingly
do q=1,turbineArray(i) % numBladePoints
do m=1,turbineModel(j) % numBl
  call rotatePoint(turbineArray(i) % bladePoints(m,q,:), rotorApex, uvShaft, &
       deltaAzimuthI)
enddo
enddo

! Compute the new azimuth angle and make sure it isn't bigger than 2*pi
if(pastFirstTimeStep) then
  azimuth = azimuth + deltaAzimuth
  azimuth = modulo(azimuth, 2._rp*pi)
endif

end subroutine atm_rotateBlades


!*******************************************************************************
subroutine atm_compute_cl_correction(i,m,q)
!*******************************************************************************
!
! This subroutine compute the tip-loss correction
! Ref. Shen et al. J. Sol. Energy Eng. 127, 209 (2005)
!
use param, only: rp, pi
use atm_input_util, only: turbineArray
use atm_input_util, only: turbineModel
implicit none
integer, intent(in) :: i,m,q
integer :: j
real(rp) :: F_tip, g, tmp
real(rp), parameter :: c1=-0.125_rp, c2=21._rp, c3=0.1_rp
real(rp), parameter :: eps=0.001_rp
real(rp), pointer :: cl, cd
real(rp) :: alpha, bladeRadius, rotSpeed, TSR
real(rp) :: TipRad, u_infinity_mean
integer :: NumBl

! variables
j = turbineArray(i) % turbineTypeID
TipRad = turbineModel(j) % TipRad
NumBl = turbineModel(j) % NumBl
rotSpeed = turbineArray(i) % rotSpeed
alpha = turbineArray(i) % alpha(m,q)
bladeRadius = turbineArray(i) % bladeRadius(m,q)
TSR = turbineArray(i) % TSR
u_infinity_mean = turbineArray(i) % u_infinity_mean

! pointers
cl => turbineArray(i) % cl(m,q)
cd => turbineArray(i) % cd(m,q)

! update TSR, tip speed ratio
if(turbineModel(j) % TorqueControllerType .ne. "fixedTSR") then
  if(pastFirstTimeStep) then
    TSR = rotSpeed*TipRad/u_infinity_mean
  endif
endif

! alpha: local angle of attack, degrees
tmp = 2._rp*bladeRadius*sin(abs(alpha)*degRad)
g = exp(-c1*(real(NumBl)*TSR-c2)) + c3

if(tmp .lt. eps) then
  F_tip = 1._rp
else
  F_tip = 2._rp/pi * acos( exp(-g*real(NumBl)*(TipRad-bladeRadius)/tmp) )
endif

! correct the lift and drag coefficients
cl = cl*F_tip
cd = cd*F_tip

end subroutine atm_compute_cl_correction


!*******************************************************************************
subroutine atm_calculate_variables(i)
!*******************************************************************************
!
! Calculates the variables of the model that need information from the input
! files. It runs after reading input information.
!
use param, only: rp, pi
use atm_input_util, only: turbineArray
use atm_input_util, only: turbineModel
use atm_base, only: interpolate
use atm_base, only: interpolate_i
implicit none
integer, intent(in) :: i    ! Indicates the turbine number
integer :: j                ! Indicates the turbine type
integer :: m, q             ! Looping indices
integer,  pointer :: NumSec ! Number of sections in lookup table
real(rp), pointer :: bladeRadius(:,:)
real(rp), pointer :: projectionRadius, projectionRadiusNacelle
real(rp), pointer :: sphereRadius
real(rp), pointer :: OverHang
real(rp), pointer :: UndSling
real(rp), pointer :: TipRad
real(rp), pointer :: PreCone

! Identifies the turbineModel being used
! The type of turbine (eg. NREL5MW)
j = turbineArray(i) % turbineTypeID 

! Number of sections in lookup table
NumSec => turbineModel(j) % NumSec

! Pointers dependent on turbineArray (i)
bladeRadius             => turbineArray(i) % bladeRadius
projectionRadius        => turbineArray(i) % projectionRadius
projectionRadiusNacelle => turbineArray(i) % projectionRadiusNacelle
sphereRadius            => turbineArray(i) % sphereRadius

! Pointers dependent on turbineType (j)
OverHang => turbineModel(j) % OverHang
UndSling => turbineModel(j) % UndSling
TipRad   => turbineModel(j) % TipRad
PreCone  => turbineModel(j) % PreCone

! First compute the radius of the force projection (to the radius where the 
! projection is only 0.001 its maximum value - this seems to recover 99.9% of 
! the total forces when integrated
! Basically epsilon=2.5*dx with dimension [m]
projectionRadius = turbineArray(i) % epsilon * sqrt(-log(.001_rp))

projectionRadiusNacelle = turbineArray(i) % nacelleEpsilon * sqrt(-log(.001_rp))

sphereRadius = sqrt(((OverHang + UndSling) + TipRad*sin(PreCone))**2 +         &
               (TipRad*cos(PreCone))**2) + projectionRadius

! Compute the lift coefficient slope without dimension
do m=1,turbineModel(j) % numBl
do q=1,turbineArray(i) % numBladePoints
  ! Interpolate quantities through section
  turbineArray(i) % twistAng(m,q) = interpolate(bladeRadius(m,q),          &
                                      turbineModel(j) % radius(1:NumSec),      &
                                      turbineModel(j) % twist(1:NumSec) )
  turbineArray(i) % chord(m,q)    = interpolate(bladeRadius(m,q),          &
                                      turbineModel(j) % radius(1:NumSec),      &
                                      turbineModel(j) % chord(1:NumSec) )
  turbineArray(i) % sectionType(m,q) = interpolate_i(bladeRadius(m,q),     &
                                      turbineModel(j) % radius(1:NumSec),      &
                                      turbineModel(j) % sectionType(1:NumSec))
enddo
enddo

end subroutine atm_calculate_variables


!*******************************************************************************
subroutine atm_computeBladeForce(i, m, q, U_local)
!*******************************************************************************
!
! This subroutine will compute the wind vectors by projecting the velocity 
! onto the transformed coordinates system
!
use param, only: rp
use atm_input_util, only: turbineArray
use atm_input_util, only: turbineModel
use atm_base, only: cross_product
use atm_base, only: interpolate
implicit none
integer, intent(in) :: i,m,q
real(rp), intent(in) :: U_local(3) ! local velocity, dimension [m/s]
integer :: j,k
integer :: sectionType
real(rp) :: twistAng, chord, windAng, db, sigma
real(rp), dimension(3) :: dragVector, liftVector
real(rp), pointer :: rotorApex(:)
real(rp), pointer :: bladeAlignedVectors(:,:,:,:)
real(rp), pointer :: windVectors(:,:,:)
real(rp), pointer :: bladePoints(:,:,:)
real(rp), pointer :: rotSpeed
real(rp), pointer :: bladeRadius(:,:)
real(rp), pointer :: PreCone
real(rp), pointer :: cl(:,:), cd(:,:), alpha(:,:)
real(rp), pointer :: Vmag(:,:)
integer,  pointer :: numSec     

! Identifier for the turbine type
j = turbineArray(i) % turbineTypeID

! Pointers to trubineArray (i)
rotorApex => turbineArray(i) % rotorApex
bladeAlignedVectors => turbineArray(i) % bladeAlignedVectors
windVectors => turbineArray(i) % windVectors
bladePoints => turbineArray(i) % bladePoints
rotSpeed => turbineArray(i) % rotSpeed
bladeRadius => turbineArray(i) % bladeRadius
cd => turbineArray(i) % cd       ! Drag coefficient
cl => turbineArray(i) % cl       ! Lift coefficient
alpha => turbineArray(i) % alpha ! Angle of attack
Vmag => turbineArray(i) % Vmag   ! Velocity magnitude

! Pointers for turbineModel (j)
NumSec => turbineModel(j) % NumSec
PreCone => turbineModel(j) % PreCone

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This will compute the vectors defining the local coordinate 
! system of the actuator point
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Define vector in z'
! If  cw rotating, this vector points along the blade towards the tip.
! If ccw rotating, this vector points along the blade towards the root.
if(turbineArray(i) % rotationDir == "cw")  then
  bladeAlignedVectors(m,q,3,:) =  bladePoints(m,q,:) - rotorApex
elseif (turbineArray(i) % rotationDir == "ccw") then
  bladeAlignedVectors(m,q,3,:) = -bladePoints(m,q,:) + rotorApex
endif
bladeAlignedVectors(m,q,3,:) = bladeAlignedVectors(m,q,3,:) /                  &
sqrt(dot_product(bladeAlignedVectors(m,q,3,:),bladeAlignedVectors(m,q,3,:)))

! Define vector in y'
bladeAlignedVectors(m,q,2,:) = cross_product(bladeAlignedVectors(m,q,3,:),     &
                               turbineArray(i) % uvShaft)
bladeAlignedVectors(m,q,2,:) = bladeAlignedVectors(m,q,2,:) /                  &
sqrt(dot_product(bladeAlignedVectors(m,q,2,:),bladeAlignedVectors(m,q,2,:)))

! Define vector in x'
bladeAlignedVectors(m,q,1,:) = cross_product(bladeAlignedVectors(m,q,2,:),     &
                               bladeAlignedVectors(m,q,3,:))
bladeAlignedVectors(m,q,1,:) = bladeAlignedVectors(m,q,1,:) /                  &
sqrt(dot_product(bladeAlignedVectors(m,q,1,:),bladeAlignedVectors(m,q,1,:)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This concludes the definition of the local coordinate system
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Now put the velocity in that cell into blade-oriented coordinates and add on 
! the velocity due to blade rotation.
windVectors(m,q,1) = dot_product(bladeAlignedVectors(m,q,1,:), U_local)
windVectors(m,q,2) = dot_product(bladeAlignedVectors(m,q,2,:), U_local) +      &
                                 rotSpeed*bladeRadius(m,q)*cos(PreCone)
windVectors(m,q,3) = dot_product(bladeAlignedVectors(m,q,3,:), U_local)

! Interpolated quantities through section
twistAng = turbineArray(i) % twistAng(m,q)
chord = turbineArray(i) % chord(m,q)
sectionType = turbineArray(i) % sectionType(m,q)

! Velocity magnitude
Vmag(m,q) = sqrt(windVectors(m,q,1)**2 + windVectors(m,q,2)**2)

! Angle between wind vector components
windAng = atan2(windVectors(m,q,1), windVectors(m,q,2)) / degRad

! Local angle of attack
alpha(m,q) = windAng - twistAng - turbineArray(i) % Pitch -                    &
             turbineArray(i) % PitchControlAngle

! Total number of entries in lists of AOA, cl and cd
k = turbineModel(j) % airfoilType(sectionType) % n

! Lift coefficient
cl(m,q) = interpolate(alpha(m,q),                                              &
          turbineModel(j) % airfoilType(sectionType) % AOA(1:k),               &
          turbineModel(j) % airfoilType(sectionType) % cl(1:k))

! Drag coefficient
cd(m,q) = interpolate(alpha(m,q),                                              &
          turbineModel(j) % airfoilType(sectionType) % AOA(1:k),               &
          turbineModel(j) % airfoilType(sectionType) % cd(1:k))

! tip-loss correction
! Ref. Shen et al. J. Sol. Energy Eng. 127, 209-213 (2005)
if(turbineArray(i) % tipALMCorrection) then
  if(pastFirstTimeStep) then
    call atm_compute_cl_correction(i,m,q)
  endif
endif

! The blade section width, dimension [m]
db = turbineArray(i) % db(q)

! Lift force, dimension [m^4/s^2]
turbineArray(i) % lift(m,q) = 0.5_rp*cl(m,q) * Vmag(m,q)**2 * chord*db 

! Drag force, dimension [m^4/s^2]
turbineArray(i) % drag(m,q) = 0.5_rp*cd(m,q) * Vmag(m,q)**2 * chord*db 

! This vector projects the drag onto the local coordinate system
dragVector = bladeAlignedVectors(m,q,1,:)*windVectors(m,q,1) +                 &
             bladeAlignedVectors(m,q,2,:)*windVectors(m,q,2)
dragVector = dragVector / sqrt(dot_product(dragVector, dragVector))

! Lift vector
liftVector = cross_product(dragVector, bladeAlignedVectors(m,q,3,:))
liftVector = liftVector / sqrt(dot_product(liftVector, liftVector))

! Apply the lift and drag as vectors
liftVector = -turbineArray(i) % lift(m,q) * liftVector
dragVector = -turbineArray(i) % drag(m,q) * dragVector

! The blade force is the total lift and drag vectors, dimension [m^4/s^2]
turbineArray(i) % bladeForces(m,q,:) = liftVector + dragVector

! Find the component of the blade element force/density in the axial 
! (along the shaft) direction.
turbineArray(i) % axialForce(m,q) = dot_product(                               &
turbineArray(i) % bladeForces(m,q,:), -turbineArray(i) % uvShaft)

! Find the component of the blade element force/density in the tangential 
! (torque-creating) direction.
turbineArray(i) % tangentialForce(m,q) = dot_product(                          &
turbineArray(i) % bladeForces(m,q,:), bladeAlignedVectors(m,q,2,:))

! Add this blade element's contribution to thrust to the total turbine thrust
! Total thrust, dimension [m^4/s^2] or [N] if density=1 [g/m^3]
turbineArray(i) % thrust = turbineArray(i) % thrust +                          &
                           turbineArray(i) % axialForce(m,q)

! Add this blade element's contribution to torque to the total turbine torque
! Total torque, dimension [m^5/s^2] or [N*m] if density=1 [g/m^3]
turbineArray(i) % torqueRotor = turbineArray(i) % torqueRotor +                &
                                turbineArray(i) % tangentialForce(m,q) *       &
                                bladeRadius(m,q)*cos(PreCone)

! Change this back to radians
windAng = windAng * degRad

! The solidity
sigma = chord * turbineModel(j) % NumBl / (2._rp*pi * bladeRadius(m,q))

! Calculate the induction factor
! Ref. Burton et al. 2011 Wind Energy Handbook (2nd ed.) John Wiley & Sons
turbineArray(i) % induction_a(m,q) = 1._rp / (4._rp * sin(windAng)**2 /        &
    (sigma * (Cl(m,q)*cos(windAng) + Cd(m,q)*sin(windAng))) + 1._rp)

! Calculate u infinity along the axial direction x'
turbineArray(i) % u_infinity(m,q) = windVectors(m,q,1)

end subroutine atm_computeBladeForce


!******************************************************************************
subroutine atm_computeNacelleForce(i, U_local)
!******************************************************************************
!
! This subroutine will compute the force from the nacelle
!
use param, only: rp
use atm_input_util, only: turbineArray
use atm_input_util, only: turbineModel
implicit none
integer, intent(in) :: i
real(rp), intent(in), dimension(3) :: U_local  ! Velocity input
integer :: j 
real(rp) :: V                                  ! Velocity projected 
real(rp), dimension(3) :: nacelleAlignedVector ! Nacelle vector
real(rp) :: area, drag

! Identifier for the turbine type
j = turbineArray(i) % turbineTypeID

area = pi * (turbineModel(j) % hubRad)**2 

nacelleAlignedVector = turbineArray(i) % uvShaft

! Velocity projected in the direction of the nacelle
V = dot_product(nacelleAlignedVector, U_local)

if(V .ge. 0._rp) then
  ! Drag force
  drag = 0.5_rp * turbineArray(i) % nacelleCd * (V**2) * area
  ! Drag Vector
  turbineArray(i) % nacelleForce = - drag * nacelleAlignedVector
endif

end subroutine atm_computeNacelleForce


!*******************************************************************************
subroutine atm_integrate_u(i)
!*******************************************************************************
!
! This subroutine will compute the induction factor a
! for each actuator point
!
use param, only: rp
use atm_input_util, only: turbineArray
implicit none
integer, intent(in) :: i
real(rp), pointer, dimension(:,:) :: u_infinity, induction_a
real(rp), pointer :: u_infinity_mean

induction_a => turbineArray(i) % induction_a
u_infinity_mean => turbineArray(i) % u_infinity_mean
u_infinity => turbineArray(i) % u_infinity

u_infinity_mean = sum(u_infinity) / size(u_infinity) /                         &
                  (1._rp - sum(induction_a) / size(induction_a))

end subroutine atm_integrate_u


!*******************************************************************************
subroutine atm_yawNacelle(i)
!*******************************************************************************
!
! This subroutine yaws the nacelle according to the yaw angle
!
use param, only: rp, pi
use atm_input_util, only: turbineArray
use atm_input_util, only: turbineModel
use atm_base, only: rotatePoint
implicit none
integer, intent(in) :: i
integer :: j 
integer :: m,q

! Perform rotation for the turbine.
j = turbineArray(i) % turbineTypeID 

! Rotate the rotor apex first.
call rotatePoint(turbineArray(i) % rotorApex,                                  &  
                 turbineArray(i) % towerShaftIntersect,                        &
                 turbineArray(i) % uvTower,                                    &
                 turbineArray(i) % deltaNacYaw)

! Recompute the shaft unit vector since the shaft has rotated.
turbineArray(i) % uvShaft = turbineArray(i) % rotorApex -                      &
                            turbineArray(i) % towerShaftIntersect

turbineArray(i) % uvShaft = turbineArray(i) % uvShaft / sqrt(dot_product(      &
turbineArray(i) % uvShaft,  turbineArray(i) % uvShaft))

turbineArray(i) % uvShaft =                                                    &
turbineArray(i) % uvShaft * turbineArray(i) % uvShaftDir

! Rotate turbine blades, blade by blade, point by point.
do q=1,turbineArray(i) % numBladePoints
do m=1,turbineModel(j) % numBl
  call rotatePoint(turbineArray(i) % bladePoints(m,q,:),                     &
                   turbineArray(i) % towerShaftIntersect,                      &
                   turbineArray(i) % uvTower,                                  &
                   turbineArray(i) % deltaNacYaw)                                    
enddo
enddo

! Compute the new yaw angle and make sure it isn't bigger than 2*pi
if(pastFirstTimeStep) then
  turbineArray(i) % nacYaw = turbineArray(i) % nacYaw +                        &
                             turbineArray(i) % deltaNacYaw
  turbineArray(i) % nacYaw = modulo(turbineArray(i) % nacYaw, 2._rp*pi)
endif

end subroutine atm_yawNacelle


!*******************************************************************************
subroutine atm_output (i, jt_total, time)
!*******************************************************************************
!
! This subroutine will calculate and write the output
!
use param, only: rp, path
use atm_input_util, only: outputInterval
use atm_input_util, only: turbineArray
use atm_input_util, only: turbineModel
implicit none
integer, intent(in) :: jt_total ! Number of iteration fed in from solver
real(rp), intent(in) :: time  ! time from simulation, dimensional
integer, intent(in) :: i  ! The turbine number
integer :: j, m
integer :: powerFile=11, rotSpeedFile=12, bladeFile=13, liftFile=14, dragFile=15
integer :: ClFile=16, CdFile=17, alphaFile=18, VrelFile=19, AzimuthFile = 40
integer :: VaxialFile=20, VtangentialFile=21, pitchFile=22, thrustFile=23
integer :: tangentialForceFile=24, axialForceFile=25, yawfile=26, nacelleFile=27
real(rp), pointer :: azimuth

azimuth => turbineArray(i)%azimuth

! Output only if the number of intervals is right
if(mod(jt_total-1, outputInterval) == 0) then
  write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  write(*,*) '!  Writing Actuator Turbine Model output  !'
  write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    
  ! The turbine type ID
  j = turbineArray(i) % turbineTypeID 

  ! File for power output
  open(unit=powerFile,position="append",                                       &
  file=trim(path)//"turbineOutput/"//trim(turbineArray(i) % turbineName)//"/power")

  ! File for thrust
  open(unit=thrustFile,position="append",                                      &
  file=trim(path)//"turbineOutput/"//trim(turbineArray(i) % turbineName)//"/thrust")

  ! File for rotor speed
  open(unit=RotSpeedFile,position="append",                                    &
  file=trim(path)//"turbineOutput/"//trim(turbineArray(i) % turbineName)//"/RotSpeed")
  
  ! File for azimuthal angle (in rad)
  open(unit=AzimuthFile,position="append",                                    &
  file=trim(path)//"turbineOutput/"//trim(turbineArray(i) % turbineName)//"/Azimuth")

  ! File for yaw
  open(unit=YawFile,position="append",                                         &
  file=trim(path)//"turbineOutput/"//trim(turbineArray(i) % turbineName)//"/Yaw")

  ! File for blade output
  open(unit=bladeFile,position="append",                                       &
  file=trim(path)//"turbineOutput/"//trim(turbineArray(i) % turbineName)//"/blade")

  open(unit=liftFile,position="append",                                        &
  file=trim(path)//"turbineOutput/"//trim(turbineArray(i) % turbineName)//"/lift")

  open(unit=dragFile,position="append",                                        &
  file=trim(path)//"turbineOutput/"//trim(turbineArray(i) % turbineName)//"/drag")
    
  open(unit=ClFile,position="append",                                          &
  file=trim(path)//"turbineOutput/"//trim(turbineArray(i) % turbineName)//"/Cl")

  open(unit=CdFile,position="append",                                          &
  file=trim(path)//"turbineOutput/"//trim(turbineArray(i) % turbineName)//"/Cd")

  open(unit=alphaFile,position="append",                                       &
  file=trim(path)//"turbineOutput/"//trim(turbineArray(i) % turbineName)//"/alpha")

  open(unit=VrelFile,position="append",                                        &
  file=trim(path)//"turbineOutput/"//trim(turbineArray(i) % turbineName)//"/Vrel")

  open(unit=VaxialFile,position="append",                                      &
  file=trim(path)//"turbineOutput/"//trim(turbineArray(i) % turbineName)//"/Vaxial")

  open(unit=VtangentialFile,position="append",                                 &
  file=trim(path)//"turbineOutput/"//trim(turbineArray(i) % turbineName)//"/Vtangential")

  open(unit=tangentialForceFile,position="append",                             &
  file=trim(path)//"turbineOutput/"//trim(turbineArray(i) % turbineName)//"/tangentialForce")

  open(unit=axialForceFile,position="append",                                  &
  file=trim(path)//"turbineOutput/"//trim(turbineArray(i) % turbineName)//"/axialForce")

  open(unit=pitchFile,position="append",                                       &
  file=trim(path)//"turbineOutput/"//trim(turbineArray(i) % turbineName)//"/pitch")

  open(unit=nacelleFile,position="append",                                     &
  file=trim(path)//"turbineOutput/"//trim(turbineArray(i) % turbineName)//"/nacelle")

  call atm_compute_power(i)

  write(powerFile   ,*) time, turbineArray(i) % powerRotor,                    &
                              turbineArray(i) % powerGen
  write(thrustFile  ,*) time, turbineArray(i) % thrust
  write(RotSpeedFile,*) time, turbineArray(i) % RotSpeed
  write(AzimuthFile,*) time,  turbineArray(i) % azimuth
  write(pitchFile   ,*) time, turbineArray(i) % PitchControlAngle,             &
                              turbineArray(i) % IntSpeedError
  write(YawFile     ,*) time, turbineArray(i) % deltaNacYaw,                   &
                              turbineArray(i) % NacYaw

  ! Will write only the first actuator section of the blade
  do m=1,turbineModel(j) % numBl
    write(bladeFile          ,*) i, m, turbineArray(i) % bladeRadius(m,:)
    write(liftFile           ,*) i, m, turbineArray(i) % lift(m,:) / turbineArray(i) % db(:)
    write(dragFile           ,*) i, m, turbineArray(i) % drag(m,:) / turbineArray(i) % db(:)
    write(ClFile             ,*) i, m, turbineArray(i) % cl(m,:)
    write(CdFile             ,*) i, m, turbineArray(i) % cd(m,:)
    write(alphaFile          ,*) i, m, turbineArray(i) % alpha(m,:)
    write(VrelFile           ,*) i, m, turbineArray(i) % Vmag(m,:)
    write(VaxialFile         ,*) i, m, turbineArray(i) % windVectors(m,:,1)
    write(VtangentialFile    ,*) i, m, turbineArray(i) % windVectors(m,:,2)
    write(tangentialForceFile,*) i, m, turbineArray(i) % tangentialForce(m,:)
    write(axialForceFile     ,*) i, m, turbineArray(i) % axialForce(m,:)
  enddo

  ! Close all the files 
  close(powerFile)
  close(thrustFile)
  close(rotSpeedFile)
  close(AzimuthFile)
  close(bladeFile)
  close(liftFile)
  close(dragFile)
  close(ClFile)
  close(CdFile)
  close(alphaFile)
  close(VrelFile)
  close(VaxialFile)
  close(VtangentialFile)
  close(pitchFile)
  close(tangentialForceFile)
  close(axialForceFile)
  close(yawFile)
  close(nacelleFile)

  write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  write(*,*) '!  Done Writing Actuator Turbine Model output  !'
  write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

  ! It is used if the simulation wants to start from a previous simulation
  ! without having to start the turbine from the original omega
  call atm_write_restart(i) 
endif

end subroutine atm_output


!******************************************************************************
subroutine atm_compute_power(i)
!******************************************************************************
!
! This subroutine will calculate the total power of the turbine
!
use atm_input_util, only: turbineArray
use atm_input_util, only: turbineModel
implicit none
integer, intent(in) :: i
integer :: j

j = turbineArray(i) % turbineTypeID

turbineArray(i) % powerRotor = turbineArray(i) % torqueRotor *                 &
                               turbineArray(i) % rotSpeed 

if(turbineModel(j) % TorqueControllerType == "fiveRegion") then
  turbineArray(i) % powerGen = turbineArray(i) % torqueGen *                   &
  turbineArray(i) % rotSpeed * turbineModel(j) % GBRatio
else
  turbineArray(i) % powerGen = turbineArray(i) % powerRotor
endif

write(*,*) 'Turbine ',i,' (Aerodynamic, Generator) Power is: ',                &
            turbineArray(i) % powerRotor, turbineArray(i) % powerGen

end subroutine atm_compute_power

end module actuator_turbine_model
#endif ATM
