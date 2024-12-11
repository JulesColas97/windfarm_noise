#ifdef TURBINES
module turbine_mod
use param, only: rp,nproc
use turbine_indicator, only: turb_ind_func_t
implicit none

integer :: turb_nodes = 0
logical :: turbine_nodes_allocated = .false.
logical :: turbine_in_proc = .false.

! Variables for one turbine 
type turbine_t
  real(rp),dimension(3) ::loc              ! location vector
  real(rp) :: dia                          ! Diameter of turbine
  real(rp) :: thk                          ! Thickness of turbine
  real(rp) :: thk_half                     ! Half thickness of turbine
  real(rp) :: vol_c                        ! term used for volume correction
  real(rp) :: turbine_vol                  ! used for normalization of indicator function
  real(rp) :: theta1                       ! angle CCW(from above) from -x direction [degrees]
  real(rp),dimension(2) :: nhat            ! unit normal for each turbine
  integer :: num_nodes                     ! number of nodes associated with each turbine
  integer :: min_i,max_i,min_j,max_j,min_k,max_k
  integer,allocatable,dimension(:,:) :: nodes ! (i,j,k) of each included node
  real(rp),allocatable,dimension(:) :: ind ! Value indicator function
  type(turb_ind_func_t) :: turb_ind_func   ! Indicator function
  real(rp) :: Ct_prime                     ! Thrust coefficient turbine
end type turbine_t

type(turbine_t), pointer, dimension(:) :: turbine

integer :: nloc                           ! number of turbines 
logical :: dyn_theta1                     ! Are turbines allowed to turn?
real(rp) :: filter_size                   ! filter size as multiple of grid spacing (indicator func)
real(rp) :: T_avg                         ! Time filter for the rotating turbines
real(rp) :: filter_cutoff                 ! indicator function threshold
integer :: tbase                          ! Number of timesteps between the output

real(rp),allocatable,dimension(:) :: sumA,buffer_array

! Parameters used for buffering the writing of turbine data
integer,parameter :: turb_max_buf = 250 ! Number of timesteps to keep in memory
integer :: turb_count_buf
real(rp),allocatable,dimension(:,:,:) :: turb_data_buf

end module turbine_mod
#endif TURBINES
