#ifdef ATM
!*******************************************************************************
module atm_input_util
!*******************************************************************************
!
! This module reads the input files for the actuator turbine model module (ATM)
! Module for dynamic allocation variables
!
use param, only: rp
implicit none
save 
private
public :: numberOfTurbines
public :: outputInterval
public :: updateInterval
public :: turbineArray
public :: turbineModel
public :: atm_read_input_conf
#ifdef ADMR
public :: numberOfActuatorLines
#endif ADMR

! The variables for the ATM are defined here
integer :: numberOfTurbines
integer :: outputInterval
integer :: updateInterval
#ifdef ADMR
integer :: numberOfActuatorLines = 1
#endif ADMR

! This type will store the necessary variables for -each- turbine
! To declare: type(turbineArray_t), allocatable, dimension(:) :: turbineArray
! To access the variables do: turbineArray(n) % variableNeeded
type turbineArray_t
  ! Variables that are read from input files
  character(128) :: turbineName            ! Name of turbine ('turbine1')
  character(128) :: turbineType            ! Name of turbine type ('NREL5MWRef')
  real(rp), dimension(3) :: baseLocation   ! Location of the base (0 0 0)
  integer :: numBladePoints                ! Number of points along each blade
  real(rp) :: epsilon                      ! Width of the smearing Gaussian function
  character(128) :: sampling               ! Sampling method for velocity atPoint or Spalart
  character(128) :: rotationDir            ! Direction of rotation ('cw')
  real(rp) :: Azimuth                      ! Angle of rotation of the rotor
  real(rp) :: RotSpeed                     ! Speed of the rotor (rpm)
  real(rp) :: Pitch
  real(rp) :: NacYaw                       ! The yaw angle of the nacelle
  real(rp) :: TSR = 0._rp                  ! Tip speed ratio
  logical :: tipALMCorrection = .false.    ! Includes a correction for tip
  logical :: nacelle = .false.             ! Includes a nacelle yes or no
  real(rp) :: nacelleEpsilon               ! Width of the smearing Gaussian function
  real(rp) :: nacelleCd = 0._rp            ! Drag coefficient for the nacelle

  ! Not read variables
  real(rp) :: thrust                       ! Total turbine thrust
  real(rp) :: torqueRotor                  ! Rotor torque
  real(rp) :: torqueGen                    ! Generator torque
  real(rp) :: powerRotor                   ! Rotor Power
  real(rp) :: powerGen                     ! Generator Power
  real(rp) :: u_infinity_mean = 0._rp      ! Mean velocity
  real(rp) :: deltaNacYaw = 0._rp          ! Change in nacelle angle
  real(rp) :: PitchControlAngle = 0._rp
  real(rp) :: IntSpeedError = 0._rp
  real(rp) :: IntPowerError = 0._rp

  ! The MPI communicator for this turbine
  integer :: TURBINE_COMM_WORLD

  ! Master processor for turbine (will write output)
  integer :: master

  ! Flag to know if this turbine operates or not in this processor
  logical :: operate = .false.

  ! Identifies the type of turbine
  integer :: turbineTypeID

  ! Important geometry data
  ! Collection of all the actuator points (blade, point, 3)
  real(rp), allocatable, dimension(:,:,:) :: bladePoints
  ! The solidity at each actuator section
  real(rp), allocatable, dimension(:,:) :: solidity
  ! Collection of radius of each point (different because of coning)
  real(rp), allocatable, dimension(:,:) :: bladeRadius
  ! Twist angle along the blade
  real(rp), allocatable, dimension(:,:) :: twistAng
  ! Chord along the blade
  real(rp), allocatable, dimension(:,:) :: chord
  ! Section type along the blade
  integer,  allocatable, dimension(:,:) :: sectionType
  ! Forces on each actuator point (blade, point, 3)
  real(rp), allocatable, dimension(:,:,:) :: bladeForces
  ! Drag force of Nacelle
  real(rp), dimension(3) :: nacelleForce
  ! Vectors at each actuator point defining the local reference frame
  ! (blade, point, 3, 3) (three vectors)
  real(rp), allocatable, dimension(:,:,:,:) :: bladeAlignedVectors
  ! The wind U projected onto the bladeAlignedVectors plus rotational speed
  ! (blade, point, 3, 3) (three vectors)
  real(rp), allocatable, dimension(:,:,:) :: windVectors
  ! Angle of attack at each each actuator point
  real(rp), allocatable, dimension(:,:) :: alpha
  ! Velocity magnitude at each each actuator point
  real(rp), allocatable, dimension(:,:) :: Vmag
  ! Lift coefficient at each actuator point
  real(rp), allocatable, dimension(:,:) :: Cl
  ! Drag coefficient at each actuator point
  real(rp), allocatable, dimension(:,:) :: Cd
  ! Lift at each actuator point
  real(rp), allocatable, dimension(:,:) :: lift
  ! Drag at each actuator point
  real(rp), allocatable, dimension(:,:) :: drag
  ! Axial force at each actuator point
  real(rp), allocatable, dimension(:,:) :: axialForce
  ! Tangential force at each actuator point
  real(rp), allocatable, dimension(:,:) :: tangentialForce

  ! Induction factor and u infinity
  real(rp), allocatable, dimension(:,:) :: induction_a
  real(rp), allocatable, dimension(:,:) :: u_infinity

  ! These are dummies meant to be used for parallelization
  ! bladeVectorDummy store quantities along the blades which are vectors
  ! such as Force
  ! bladeScalarDummy will store scalar quantities along the blades such
  ! as lift coefficitent, angle of attack, etc
  real(rp), allocatable, dimension(:,:,:) :: bladeVectorDummy
  real(rp), allocatable, dimension(:,:) :: bladeScalarDummy

  ! An indicator of shaft direction. The convention is that when viewed
  ! from upwind, the rotor turns clockwise for positive rotation angles,
  ! regardless of if it is an upwind or downwind turbine. uvShaft is
  ! found by subtracting the rotor apex location from the tower shaft
  ! intersection point. This vector switches direciton depending on
  ! if the turbine is upwind or downwind, so this uvShaftDir multiplier
  ! makes the vector consistent no matter what kind of turbine
  real(rp) :: uvShaftDir

  ! Define the vector along the shaft pointing in the direction of the wind
  real(rp), dimension(3) :: uvShaft

  ! List of locations of the rotor apex relative to the origin (m)
  real(rp), dimension(3) :: rotorApex

  ! List of locations of the intersection of the tower axis and the shaft
  ! centerline relative to the origin (m).
  real(rp), dimension(3) :: towerShaftIntersect

  ! Unit vector pointing along the tower (axis of yaw).
  real(rp), dimension(3) :: uvTower

  ! Width of the actuator section
  real(rp), allocatable, dimension(:) :: db

  ! Sphere radius which defines a sphere from the center of the rotor and
  ! identifies the volume onto which forces are applied
  real(rp) :: projectionRadius
  real(rp) :: projectionRadiusNacelle
  real(rp) :: sphereRadius

  ! Location of Nacelle
  real(rp), dimension(3) :: nacelleLocation

  ! Change of azimuth angle every time-step
  real(rp) :: deltaAzimuth=0._rp
end type turbineArray_t

! This type stores all the airfoils and their AOA, Cd, Cl values
type airfoilType_t
  character(128) :: airfoilName     ! The type of Airfoil ('Cylinder1')
  integer :: n                      ! Number of data points
  ! The maximum number of points is chosen to be 150. If airfoil data has
  ! more than this then this number should be modified!
  real(rp), dimension(150) :: AOA   ! Angle of Attack
  real(rp), dimension(150) :: Cd    ! Drag coefficient
  real(rp), dimension(150) :: Cl    ! Lift coefficient
  real(rp), dimension(150) :: Cm    ! Moment coefficient
end type airfoilType_t

type turbineModel_t
  character(128) :: turbineType     ! The type of turbine ('NREL5MWRef')
  integer :: NumBl                  ! Number of turbine blades
  integer :: NumSec                 ! Number of sections
  real(rp) :: TipRad                ! Radius from the root to the tip of the blade
  real(rp) :: HubRad                ! Radius of the hub
  real(rp) :: UndSling              ! Undersling length [distance from teeter pin to the rotor apex]
  real(rp) :: OverHang              ! Distance from yaw axis to rotor apex
  real(rp) :: TowerHt               ! Tower height
  real(rp) :: Twr2Shft              ! Vertical distance from tower-top to rotor shaft
  real(rp) :: ShftTilt              ! Rotor shaft tilt angle
  real(rp) :: PreCone               ! Angle which the blades are coned
  real(rp) :: GBRatio               ! Gearbox Ratio
  real(rp) :: GenIner               ! Generator inertia
  real(rp) :: HubIner               ! Inertia of the hub
  real(rp) :: BladeIner             ! Inertia of the blades
  real(rp) :: DriveTrainIner        ! Inertia of the Drive train

  ! Blade section quantities (maximum number of sections 100, easy modify)
  real(rp), dimension(100) :: chord, twist, radius
  integer, dimension(100) :: sectionType

  ! The airfoil type properties (includes AOA, Cl, and Cd) Attempt 1
  type(airfoilType_t), allocatable, dimension(:) :: airfoilType

  ! Torque controller variables
  character(64) :: TorqueControllerType
  real(rp) :: CutInGenSpeed
  real(rp) :: RatedGenSpeed
  real(rp) :: Region2StartGenSpeed
  real(rp) :: Region2EndGenSpeed
  real(rp) :: CutInGenTorque
  real(rp) :: RatedGenTorque
  real(rp) :: RateLimitGenTorque
  real(rp) :: KGen
  real(rp) :: TorqueControllerRelax

  ! Pitch controller variables
  character(64) :: PitchControllerType
  real(rp) :: PitchControlAngleMax     ! Maximum pitch angle
  real(rp) :: PitchControlAngleMin     ! Minimum pitch angle
  real(rp) :: PitchControlAngleK       ! Angle at which sensitivity doubles
  real(rp) :: PitchControlKP0          ! Proportional term at angle = 0
  real(rp) :: PitchControlKI0          ! Integral term at angle = 0

  ! Yaw controller variables
  character(64) :: YawControllerType   ! Name of yaw controller type
  character(64) :: YawControllerFile   ! File that contains the yaw time info
  real(rp), allocatable, dimension(:) :: yaw_time  ! time for yaw (seconds)
  real(rp), allocatable, dimension(:) :: yaw_angle ! yaw angle (degrees)
end type turbineModel_t

! Declare turbine array variable
type(turbineArray_t), allocatable, dimension(:), target :: turbineArray

! Declare turbine model variable (stores information for turbine models)
type(turbineModel_t), allocatable, dimension(:), target :: turbineModel

! Name of the utility used
character (*), parameter :: input_conf = './inputATM/turbineArrayProperties'
character (*), parameter :: comment = '!'
character (*), parameter :: block_entry = '{' ! The start of a block
character (*), parameter :: block_exit = '}'  ! The end of a block
character (*), parameter :: equal = '='
character (*), parameter :: esyntax = 'syntax error at line'
character (*), parameter :: array_entry = '(' ! The start of an array
character (*), parameter :: array_exit  = ')' ! The end of an array

! Delimiters used for reading vectors and points
character(*), parameter :: delim_minor=','
character(*), parameter :: delim_major='//'

! Variables used to read lines in file
integer :: block_entry_pos  ! Determines if the line is the start of a block
integer :: block_exit_pos   ! Determines if the line is the end of a block
integer :: array_entry_pos  ! Determines if the line is the start of a block
integer :: array_exit_pos   ! Determines if the line is the end of a block
integer :: equal_pos        ! Determines if there is an equal sign
integer :: ios
logical :: exst             ! Used to determine existence of a file

contains

!*******************************************************************************
subroutine atm_read_input_conf
!*******************************************************************************
implicit none
character(128) :: buff ! Stored the read line
integer :: n = 0       ! Counter for the wind turbines
integer :: lun =1      ! Reference number for input file
integer :: line        ! Counts the current line in a file

! Check that the configuration file exists
inquire (file=input_conf, exist=exst)

! Open file
if (exst) then
  ! Open the input file
  open (lun, file=input_conf, action='read')
else
  ! Error for non existing file
  print*, 'file ' // input_conf // ' does not exist'
  stop
endif

! Read the file line by line *Counter starts at 0 and modified inside subroutine
line = 0
do
  ! Read line by line (lun=file number)
  call atm_readline(lun, line, buff, block_entry_pos, block_exit_pos,          &
                    array_entry_pos, array_exit_pos, equal_pos, ios)

  if(ios /= 0) exit ! Exit if reached end of file

  ! This will read the numberOfTurbines integer
  if(buff(1:16) == 'numberOfTurbines') then
    read(buff(17:), *) numberOfTurbines
    ! Allocate space for the wind turbine variables
    allocate(turbineArray(numberOfTurbines))
    cycle
  elseif(buff(1:14) == 'outputInterval') then
    read(buff(15:), *) outputInterval
    cycle
  elseif(buff(1:14) == 'updateInterval') then
    read(buff(15:), *) updateInterval
    cycle
#ifdef ADMR
  elseif(buff(1:21) == 'numberOfActuatorLines') then
    read(buff(22:), *) numberOfActuatorLines
    cycle
#endif ADMR
  endif

  ! This will start reading turbine name
  if(block_entry_pos /= 0) then 
    ! Increment turbine counter
    n = n + 1 
    if(n .gt. numberOfTurbines) exit
    ! Read the name of the turbine
    read(buff(1:index(buff, block_entry)-1), *) turbineArray(n) % turbineName
  endif

  ! This will start reading turbine block
  if(block_entry_pos == 0) then
    if(buff(1:11) == 'turbineType') then
      read(buff(12:), *) turbineArray(n) % turbineType
    endif
    if(buff(1:12) == 'baseLocation') then
      read(buff(13:), *) turbineArray(n) % baseLocation
    endif
    if(buff(1:14) == 'numBladePoints') then
      read(buff(15:), *) turbineArray(n) % numBladePoints
    endif
    if(buff(1:7) == 'epsilon') then
      read(buff(8:), *) turbineArray(n) % epsilon
    endif
    if(buff(1:8) == 'sampling') then
      read(buff(9:), *) turbineArray(n) % sampling
    endif
    if(buff(1:11) == 'rotationDir') then
      read(buff(12:), *) turbineArray(n) % rotationDir
    endif
    if(buff(1:7) == 'Azimuth') then
      read(buff(8:), *) turbineArray(n) % Azimuth
    endif
    if(buff(1:8) == 'RotSpeed') then
      read(buff(9:), *) turbineArray(n) % RotSpeed
    endif
    if(buff(1:5) == 'Pitch') then
      read(buff(6:), *) turbineArray(n) % Pitch
    endif
    if(buff(1:6) == 'NacYaw') then
      read(buff(7:), *) turbineArray(n) % NacYaw
    endif
    if(buff(1:11) == 'nacelleFlag') then
      read(buff(12:), *) turbineArray(n) % nacelle
    endif
    if(buff(1:9) == 'nacelleCd') then
      read(buff(10:), *) turbineArray(n) % nacelleCd
    endif
    if(buff(1:14) == 'nacelleEpsilon') then
      read(buff(15:), *) turbineArray(n) % nacelleEpsilon
    endif
    if(buff(1:3) == 'TSR') then
      read(buff(4:), *) turbineArray(n) % TSR
    endif
    if(buff(1:16) == 'tipALMCorrection') then
      read(buff(17:), *) turbineArray(n) % tipALMCorrection
    endif
  endif
end do

close (lun)

if(.not. allocated(turbineArray)) then
  write(*,*) 'Did not allocate memory for turbineArray'
  stop
endif

call read_turbine_model_variables

end subroutine atm_read_input_conf


!*******************************************************************************
subroutine read_turbine_model_variables
!*******************************************************************************
use param, only: rp
implicit none
integer :: i, j, c                   ! counter
integer :: numTurbinesDistinct       ! Number of different turbine types
character(128) :: currentTurbineType ! Will store turbineType in loop
character(139) :: input_turbine
integer :: lun =19                   ! Reference number for input file
integer :: line                      ! Counts the current line in a file
character (128) :: buff              ! Stored the read line
integer:: numAirfoils                ! Number of distinct airfoils
integer :: k, p                      ! Used to loop through aifoil types 
integer :: numBladePoints, numBl, numSec
! Variables for reading yaw control file
integer :: ios, N
integer :: yawfile = 29
! Name of all the airfoils types (max 20) If more than this increase the number
character(128), dimension(20) :: airfoils

! Initialize variables for the loop
! Will find the number of distincet turbines to allocate memory for turbineModel
numTurbinesDistinct = 0
currentTurbineType = turbineArray(1) % turbineType

! Counter for number of distinct turbines
c = 0
! Find how many turbine types there are
do i=1,numberOfTurbines

  ! Find if the name is repeated
  do j=1,i-1
    if(turbineArray(i) % turbineType .eq. turbineArray(j) % turbineType) then
      c = 1
      turbineArray(i) % turbineTypeID = turbineArray(j) % turbineTypeID
    endif
  enddo

  ! Assign a new turbine if c = 0
  if(c .eq. 0) then
    numTurbinesDistinct = numTurbinesDistinct + 1
    turbineArray(i) % turbineTypeID = numTurbinesDistinct
  endif

  ! Restart the counter at 0
  c = 0
enddo

! Allocate space for turbine model variables
allocate(turbineModel(numTurbinesDistinct))

do i=1,numberOfTurbines
  turbineModel(turbineArray(i) % turbineTypeID) % turbineType =                &
               turbineArray(i) % turbineType
enddo

! Read the input properties for each turbine type
do i=1,numTurbinesDistinct
  input_turbine = './inputATM/' // turbineModel(i) % turbineType

  ! Check that the configuration file exists
  inquire(file=input_turbine, exist=exst)

  if(exst) then
    ! Open the input file
    open(lun, file=input_turbine, action='read')
  else
    ! Error for non existing file
    print*, 'file ' // input_turbine // ' does not exist'
    stop
  endif

  ! Read the file line by line - starts at 0 and modified inside subroutine
  line = 0
  do
    ! Read line by line (lun=file number)
    call atm_readline(lun, line, buff, block_entry_pos, block_exit_pos,        &
                      array_entry_pos, array_exit_pos, equal_pos, ios)
    if(ios /= 0) exit

    ! This will all the input variables
    if(buff(1:5) == 'NumBl') then
      read(buff(6:), *) turbineModel(i) % NumBl
#ifdef ADMR
      turbineModel(i) % NumBl = turbineModel(i) % NumBl * numberOfActuatorLines
#endif ADMR
    endif
    if(buff(1:6) == 'TipRad') then
      read(buff(7:), *) turbineModel(i) % TipRad
    endif
    if(buff(1:6) == 'HubRad') then
      read(buff(7:), *) turbineModel(i) % HubRad
    endif
    if(buff(1:8) == 'UndSling') then
      read(buff(9:), *) turbineModel(i) % UndSling
    endif
    if(buff(1:8) == 'OverHang') then
      read(buff(9:), *) turbineModel(i) % OverHang
    endif
    if(buff(1:7) == 'TowerHt') then
      read(buff(8:), *) turbineModel(i) % TowerHt
    endif
    if(buff(1:8) == 'Twr2Shft') then
      read(buff(9:), *) turbineModel(i) % Twr2Shft
    endif
    if(buff(1:8) == 'ShftTilt') then
      read(buff(9:), *) turbineModel(i) % ShftTilt
    endif
    if(buff(1:7) == 'PreCone') then
      read(buff(8:), *) turbineModel(i) % PreCone
    endif

    ! This will read the torque controller type and its properties
    if(buff(1:20) == 'TorqueControllerType') then
      read(buff(21:), *) turbineModel(i) % TorqueControllerType
    endif
    if(buff(1:13) == 'CutInGenSpeed') then
      read(buff(14:), *) turbineModel(i) % CutInGenSpeed
    endif
    if(buff(1:13) == 'RatedGenSpeed') then
      read(buff(14:), *) turbineModel(i) % RatedGenSpeed
    endif
    if(buff(1:20) == 'Region2StartGenSpeed') then
      read(buff(21:), *) turbineModel(i) % Region2StartGenSpeed
    endif
    if(buff(1:18) == 'Region2EndGenSpeed') then
      read(buff(19:), *) turbineModel(i) % Region2EndGenSpeed
    endif
    if(buff(1:14) == 'CutInGenTorque') then
      read(buff(15:), *) turbineModel(i) % CutInGenTorque
    endif
    if(buff(1:14) == 'RatedGenTorque') then
      read(buff(15:), *) turbineModel(i) % RatedGenTorque
    endif
    if(buff(1:18) == 'RateLimitGenTorque') then
      read(buff(19:), *) turbineModel(i) % RateLimitGenTorque
    endif
    if(buff(1:4) == 'KGen') then
      read(buff(5:), *) turbineModel(i) % KGen
    endif
    if(buff(1:21) == 'TorqueControllerRelax') then
      read(buff(22:), *) turbineModel(i) % TorqueControllerRelax
    endif
    if(buff(1:7) == 'GBRatio') then
      read(buff(8:), *) turbineModel(i) % GBRatio
    endif
    if(buff(1:9) == 'BladeIner') then
      read(buff(10:), *) turbineModel(i) % BladeIner
    endif
    if(buff(1:7) == 'HubIner') then
      read(buff(8:), *) turbineModel(i) % HubIner
    endif
    if(buff(1:7) == 'GenIner') then
      read(buff(8:), *) turbineModel(i) % GenIner
    endif

    ! This will read the pitch control values
    if(buff(1:19) == 'PitchControllerType') then
      read(buff(20:), *) turbineModel(i) % PitchControllerType
    endif
    if(buff(1:20) == 'PitchControlAngleMax') then
      read(buff(21:), *) turbineModel(i) % PitchControlAngleMax
    endif
    if(buff(1:20) == 'PitchControlAngleMin') then
      read(buff(21:), *) turbineModel(i) % PitchControlAngleMin
    endif
    if(buff(1:18) == 'PitchControlAngleK') then
      read(buff(19:), *) turbineModel(i) % PitchControlAngleK
    endif
    if(buff(1:15) == 'PitchControlKP0') then
      read(buff(16:), *) turbineModel(i) % PitchControlKP0
    endif
    if(buff(1:15) == 'PitchControlKI0') then
      read(buff(16:), *) turbineModel(i) % PitchControlKI0
    endif

    ! This will read the yaw control values
    if(buff(1:17) == 'YawControllerType') then
      read(buff(18:), *) turbineModel(i) % YawControllerType
    endif
    if(buff(1:17) == 'YawControllerFile') then
      read(buff(18:), *) turbineModel(i) % YawControllerFile
    endif

    ! This will read the airfoils 
    if(buff(1:8) == 'Airfoils') then
      ! Conuter for the number of distince airfoils
      numAirfoils=0
      ! If 'Airfoils(' then make array_entry_pos zero
      array_entry_pos=0

      do while(array_entry_pos == 0)
        call atm_readline(lun, line, buff, block_entry_pos, block_exit_pos,    &
                          array_entry_pos, array_exit_pos, equal_pos, ios)
        ! exit if end of file reached
        if(ios /= 0) exit 

        if(array_entry_pos /= 0) then
          call atm_readline(lun, line, buff, block_entry_pos, block_exit_pos,  &
                            array_entry_pos, array_exit_pos, equal_pos, ios)
          ! exit if end of file reached
          if(ios /= 0) exit
        endif

        ! exit if end of file reached
        if(array_exit_pos /= 0) exit

        ! Increment airfoil counter
        numAirfoils = numAirfoils + 1
        ! Stores the name of the airfoil
        airfoils(numAirfoils) = buff
      enddo

      ! Allocate the airfoilTypes
      allocate(turbineModel(i) % airfoilType(numAirfoils))

      ! Loop through all the airfoil types and look-up lift and drag
      do k=1,numAirfoils
        ! Eliminate white space
        call eat_whitespace(airfoils(k))
        ! Length without last element
        p = len(trim(airfoils(k)))-1
        ! Airfoil type (2:p) is used to eliminate the commas ""
        turbineModel(i) % airfoilType % airfoilName = airfoils(k)(2:p)
        ! Read each airfoil accordingly
        call read_airfoil(turbineModel(i) % airfoilType(k))
      enddo
    endif

    ! This will read the blade properties
    if(buff(1:9) .eq. 'BladeData') then
      NumSec = 0

      do while(array_exit_pos .eq. 0)
        ! Read the comment line
        read(lun, '(a)', iostat=ios) buff
        if(scan(buff,'(') .ne. 0) cycle
        if(scan(buff,'!') .ne. 0) cycle
        if(scan(buff,')') .ne. 0) exit

        ! Number of sections will account for all blade sections
        NumSec = NumSec+1
        turbineModel(i) % NumSec = NumSec

        ! Read in radius, chord, twist, type
        read(buff,*) turbineModel(i) % radius(NumSec),                         &
                     turbineModel(i) % chord(NumSec),                          &
                     turbineModel(i) % twist(NumSec),                          &
                     turbineModel(i) % sectionType(NumSec)

        ! Add one to airfoil identifier. 
        ! List starts at 0, now it will start at 1
        turbineModel(i) % sectionType(NumSec) =                                &
        turbineModel(i) % sectionType(NumSec) + 1
      enddo
    endif
  enddo

  close (lun)

  ! Read the yaw control file, if applicable
  if(turbineModel(i) % YawControllerType == "timeYawTable") then
    open(unit = yawfile, file="inputATM/"//                                    &
         trim(turbineModel(i) % YawControllerFile),                            &
         form = "formatted", status = "old", action = "read")

    ! Determine number of lines
    ios = 0
    N = -1
    do while(ios == 0)
      N = N + 1
      read(unit = yawfile, fmt = *, iostat = ios)
    enddo
    rewind(unit = yawfile)

    ! Allocate variables
    allocate(turbineModel(i) % yaw_time(N))
    allocate(turbineModel(i) % yaw_angle(N))

    ! now read the variables and close file
    do j=1,N
      read(unit = yawfile, fmt = *) turbineModel(i) % yaw_time(j),             &
                                    turbineModel(i) % yaw_angle(j)
    enddo

    close(unit = yawfile)
  endif

  ! Calculate drive train inertia
  turbineModel(i) % DriveTrainIner = (real(turbineModel(i) % NumBl,rp))    *   &
               (turbineModel(i) % BladeIner) + (turbineModel(i) % HubIner) +   &
               (turbineModel(i) % GBRatio  ) * (turbineModel(i) % GBRatio) *   &
               (turbineModel(i) % GenIner  )
enddo

! Allocate variables inside turbineArray
do i=1,numberOfTurbines
  numBladePoints = turbineArray(i) % numBladePoints
  j = turbineArray(i) % turbineTypeID
  numBl = turbineModel(j) % numBl

  allocate(turbineArray(i) % bladeForces(numBl,numBladePoints,3))
  allocate(turbineArray(i) % bladeAlignedVectors(numBl,numBladePoints,3,3))
  allocate(turbineArray(i) % windVectors(numBl,numBladePoints,3))
  allocate(turbineArray(i) % alpha(numBl,numBladePoints))
  allocate(turbineArray(i) % Vmag(numBl,numBladePoints))
  allocate(turbineArray(i) % Cl(numBl,numBladePoints))
  allocate(turbineArray(i) % Cd(numBl,numBladePoints))
  allocate(turbineArray(i) % lift(numBl,numBladePoints))
  allocate(turbineArray(i) % drag(numBl,numBladePoints))
  allocate(turbineArray(i) % axialForce(numBl,numBladePoints))
  allocate(turbineArray(i) % tangentialForce(numBl,numBladePoints))
  allocate(turbineArray(i) % induction_a(numBl,numBladePoints))
  allocate(turbineArray(i) % u_infinity(numBl,numBladePoints))
  allocate(turbineArray(i) % chord(numBl,numBladePoints))
  allocate(turbineArray(i) % twistAng(numBl,numBladePoints))
  allocate(turbineArray(i) % sectionType(numBl,numBladePoints))

  ! Variables meant for parallelization
  allocate(turbineArray(i) % bladeVectorDummy(numBl,numBladePoints,3))
  allocate(turbineArray(i) % bladeScalarDummy(numBl,numBladePoints))

enddo

end subroutine read_turbine_model_variables


!*******************************************************************************
subroutine read_airfoil(airfoilType)
!*******************************************************************************
!
! This subroutine reads the angle of attack, lift and drag for a specific
! airfoil type
!
use param, only: rp
implicit none
integer :: lun=17 ! File identifier
integer :: q      ! Counter
character(128) :: input_airfoil
type(airfoilType_t), intent(inout) :: airfoilType
real(rp) :: AOA, Cd, Cl, Cm

! Name of the input file to read
input_airfoil='./inputATM/AeroData/'//trim (airfoilType % airfoilName)//'.dat'
open (lun, file=input_airfoil, action='read')

! Skip the first 14 lines of the Fortran input file
do q=1,14
  read(lun,*)
enddo

! Initialize the counter back to 1 for the first element in the list
q = 0
AOA = -181._rp
do while(AOA .lt. 180._rp)
  q = q+1
  read(lun,*) AOA, Cl , Cd, Cm
  airfoilType % AOA(q) = AOA
  airfoilType % Cd(q) = Cd
  airfoilType % Cl(q) = Cl
  airfoilType % Cm(q) = Cm
enddo
airfoilType % n = q

close(lun)

end subroutine read_airfoil


!*******************************************************************************
subroutine atm_readline(lun, line, buff, block_entry_pos, block_exit_pos,      &
                        array_entry_pos, array_exit_pos, equal_pos, ios)
!*******************************************************************************
!
! This subroutine reads the specified line and determines the attributes of the
! contents of the line.
!
implicit none
integer, intent(in) :: lun
integer, intent(inout) :: line
character(*), intent(inout) :: buff
integer, intent(out) :: block_entry_pos, block_exit_pos, equal_pos, ios,       &
                        array_entry_pos, array_exit_pos

block_entry_pos = 0
block_exit_pos = 0
equal_pos = 0
ios = -1

do
  line = line + 1
  read(lun, '(a)', iostat=ios) buff
  if(ios /= 0) exit
  call eat_whitespace(buff)

  ! drop blank lines
  if(verify (buff, ' ') == 0) cycle 
  ! drop comment lines
  if(buff (1:len (comment)) == comment) cycle 

  block_entry_pos = index(buff, block_entry)
  block_exit_pos  = index(buff, block_exit )
  array_entry_pos = index(buff, array_entry)
  array_exit_pos  = index(buff, array_exit )
  equal_pos       = index(buff, equal      )

  exit
enddo

end subroutine atm_readline


!*******************************************************************************
subroutine eat_whitespace(buff, whtspc)
!*******************************************************************************
!
! eats leading and intermediate whitespace, fill trailing space with blanks
!
implicit none
character(*), intent(inout) :: buff
character(*), intent(in), optional :: whtspc
! add more characters here if needed
character(*), parameter :: whtspc_default = achar (9) // achar (32)
character(1), parameter :: fill_char = ' '
character(1) :: tmp (len(buff))
character(1) :: fill(len(buff))

fill = fill_char
tmp = transfer(buff, tmp)

if(present(whtspc)) then
  tmp = pack(tmp, scan(tmp, whtspc) == 0, fill)
else
  tmp = pack(tmp, scan(tmp, whtspc_default) == 0, fill)
endif

buff = transfer(tmp, buff)

end subroutine eat_whitespace

end module atm_input_util
#endif ATM
