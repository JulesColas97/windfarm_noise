#ifdef WINDBREAKS
module windbreak_mod
use param, only: rp,nproc
implicit none

integer :: windbreak_nodes = 0
logical :: windbreak_in_proc = .false.

! Variables for one windbreak
type windbreak_t
  real(rp),dimension(3) ::loc              ! location vector
  real(rp) :: width                        ! Width of windbreak
  real(rp) :: height                       ! Height of windbreak
  real(rp) :: thk                          ! Thickness of windbreak
  real(rp) :: vol_c                        ! term used for volume correction
  real(rp) :: vol                          ! used for normalization of indicator function
  real(rp) :: theta1                       ! angle CCW(from above) from -x direction [degrees]
  real(rp),dimension(2) :: nhat            ! unit normal for each windbreak
  integer :: num_nodes                     ! number of nodes associated with each windbreak
  integer :: min_i,max_i,min_j,max_j,min_k,max_k
  integer,allocatable,dimension(:,:) :: nodes ! (i,j,k) of each included node
  real(rp),allocatable,dimension(:) :: ind ! Value indicator function
  real(rp) :: Cd                           ! Drag coefficient windbreak, depends on porosity
end type windbreak_t

type(windbreak_t), pointer, dimension(:) :: windbreak

! parameters read from input.conf
integer  :: nloc_wb = 0                       ! number of windbreaks
real(rp) :: filter_size_wb = 1._rp           ! filter size as multiple of grid spacing
real(rp) :: filter_cutoff_wb = 0.01_rp       ! indicator function threshold

real(rp) :: sqrt6overdelta_wb
real(rp),allocatable,dimension(:) :: sumB,buffer_array_wb

end module windbreak_mod
#endif WINDBREAKS
