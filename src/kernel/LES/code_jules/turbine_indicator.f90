#ifdef TURBINES
module turbine_indicator 
use param, only: rp
save
private
public :: turb_ind_func_t, sqrt6overdelta

type turb_ind_func_t
  real(rp), dimension(:), allocatable :: r
  real(rp), dimension(:), allocatable :: R23
  real(rp) :: M
contains
  procedure, public :: init
  procedure, public :: val
end type turb_ind_func_t

real(rp) :: sqrt6overdelta
integer, dimension(:), allocatable :: ind
real(rp), dimension(:), allocatable :: yz
real(rp), dimension(:,:), allocatable :: g, f, h
real(rp), dimension(:), allocatable :: xi
real(rp), dimension(:), allocatable :: dummy1D
real(rp), dimension(:,:), allocatable :: dummy2D
complex(rp), dimension(:,:), allocatable :: ghat, fhat, hhat

contains

!*******************************************************************************
function val(this, r_disk) result(Rval)
!*******************************************************************************
!
! Calculate the Gaussian-filtered indicator function, which equals to 
! a product of a normal component and a radial component.
!
! Ref: 
!     Shapiro, Gayme and Meneveau (2019) Filtered actuator disks: Theory
!     and application to wind turbine models in large eddy simulation. 
!     Wind Energy 1–7.
!
use param, only: R1_ind, rp
implicit none 
class(turb_ind_func_t), intent(in) :: this
real(rp), intent(in) :: r_disk
real(rp) :: Rval

Rval = R1_ind*linear_interp(this%r, this%R23, r_disk)

end function val


!*******************************************************************************
subroutine init(this, delta2, dia, ctprime)
!*******************************************************************************
!
! Calculate the radial component of Gaussian-filtered indicator function
!
! Ref: 
!     Shapiro, Gayme and Meneveau (2019) Filtered actuator disks: Theory
!     and application to wind turbine models in large eddy simulation. 
!     Wind Energy 1–7.
!
use param, only: nx,ny,nz,pi,rp
use fft, only: FFTW_ESTIMATE
implicit none
class(turb_ind_func_t), intent(inout) :: this
real(rp), intent(in) :: delta2, dia, ctprime
integer :: N 
real(rp) :: L, d, R
real(rp) :: dr, Lr
integer :: i, j
integer*8 :: plan

N = max(nx,ny,nz) 
L = 4._rp * dia
d = L / N
R = 0.5_rp * dia

allocate(yz(N))
allocate(dummy1D(N)) 
allocate(ind(N))
allocate(g(N, N))
allocate(h(N, N))
allocate(f(N, N))
allocate(dummy2D(N, N)) 
allocate(ghat(N/2+1, N))
allocate(hhat(N/2+1, N))
allocate(fhat(N/2+1, N))

! Calculate constants
sqrt6overdelta = sqrt(6._rp) / sqrt(delta2)

! Calculate yz and indices to sort the result
do i = 1, N/2
  yz(i) = d*(i-0.5_rp)
  ind(i) = N/2+i
end do
do i = N/2+1, N
  yz(i) = -L + d*(i-0.5_rp)
  ind(i) = i-N/2
end do

! Calculate g and f
do j = 1, N
do i = 1, N
  g(i,j) = exp(-6._rp*(yz(i)**2+yz(j)**2)/delta2)
  if ((yz(i)**2 + yz(j)**2) < R**2) then
    h(i,j) = 1._rp
  else
    h(i,j) = 0._rp
  end if
end do
end do

! Do the convolution f = g*h in fourier space
call dfftw_plan_dft_r2c_2d(plan, N, N, g, ghat, FFTW_ESTIMATE)
call dfftw_execute_dft_r2c(plan, g, ghat)
call dfftw_destroy_plan(plan)

call dfftw_plan_dft_r2c_2d(plan, N, N, h, hhat, FFTW_ESTIMATE)
call dfftw_execute_dft_r2c(plan, h, hhat)
call dfftw_destroy_plan(plan)

fhat = ghat*hhat

! Compute the inverse fft of fhat
call dfftw_plan_dft_c2r_2d(plan, N, N, fhat, f, FFTW_ESTIMATE)
call dfftw_execute_dft_c2r(plan, fhat, f)
call dfftw_destroy_plan(plan)

! Normalize
f = f / N**2 * d**2

dummy2D = f 
dummy1D = yz 

!Sort the results 
do i=1,N 
  do j=1,N 
    f(j,i) =  dummy2D(ind(j),ind(i)) 
  enddo 
  yz(i) = dummy1D(ind(i)) 
enddo 

! Sort the results
!f = f(ind,ind)
!yz = yz(ind);

! Interpolate onto the lookup table
allocate(xi(N))
allocate(this%r(N))
allocate(this%R23(N))

Lr = R + 2._rp * sqrt(delta2)
dr = Lr / (N - 1)
do i = 1,N
  this%r(i) = real(i-1)*dr
  xi(i) = 0._rp
end do
this%R23 = bilinear_interp(yz, yz, f, xi, this%r)
this%R23 = this%R23 / this%R23(1)

! Calculate correction factor (see Shapiro, Gayme and Meneveau (2019) "Filtered
! actuator disks: ...", formula 26)
this%M = 1._rp / (1._rp + ctprime * sqrt(delta2) / (4._rp * sqrt(3._rp*pi) * R))

! Version in LESGO:
!this%M = sum(2._rp*pi*this%r*this%R23**2*dr)*pi (Is there a R^2 missing in
!formula ? , compared to 22 in their paper)

deallocate(xi,yz,ind,g,h,f,ghat,hhat,fhat)
deallocate(dummy2D, dummy1D)

end subroutine init


!*******************************************************************************
function bilinear_interp(x, y, v, xq, yq) result(vq)
!*******************************************************************************
use param, only: rp
implicit none
real(rp), dimension(:), intent(in) :: x, y !x=y=yz
real(rp), dimension(:,:), intent(in) :: v !v=f  
real(rp), dimension(:), intent(in) :: xq, yq !xq=xi=0,yq=this%r
real(rp), dimension(:), allocatable :: vq
real(rp) :: v1,v2
integer  :: i,j,k,N

N = size(x)
allocate(vq(N))

i = binary_search(x, xq(1)) ! = N/2?!
do k = 1,N
  j = binary_search(y, yq(k))
  if (j == 0) then
    vq(k) = v(i,1) + (xq(k) - x(i)) * (v(i+1,1)-v(i,1)) / (x(i+1) - x(i))
  else if (j == N) then
    vq(k) = v(i,N) + (xq(k) - x(i)) * (v(i+1,N)-v(i,N)) / (x(i+1) - x(i))
  else
    v1 = v(i,j)   + (xq(k) - x(i)) * (v(i+1,j)   - v(i,j)  ) / (x(i+1) - x(i))
    v2 = v(i,j+1) + (xq(k) - x(i)) * (v(i+1,j+1) - v(i,j+1)) / (x(i+1) - x(i))
    vq(k) = v1 + (yq(k) -y(j)) * (v2 - v1) / (y(j+1) - y(j))
  end if
end do

end function bilinear_interp


!*******************************************************************************
function linear_interp(x, v, xq) result(vq)
!*******************************************************************************
use param, only: rp
implicit none
real(rp), dimension(:), intent(in) :: x, v
real(rp), intent(in) :: xq
real(rp) :: vq
integer :: i, N

N = size(v)
i = binary_search(x, xq)
if (i == 0) then
  vq = v(1)
else if (i == N) then
  vq = v(N)
else
  vq = v(i) + (xq - x(i)) * (v(i+1)-v(i)) / (x(i+1) - x(i))
end if

end function linear_interp


!*******************************************************************************
function binary_search(arr,val) result(low)
!*******************************************************************************
!
!  This function performs a binary search on a sorted array. Given the
!  provided value, adjacent low and high indices of the array are found
!  such that the provided value is bracketed. Guaranteed log2(N) search.
!
!  Inputs:
!  arr          - sorted array of values to search
!  val          - value to be bracketed
!
!  Output:
!  low          - lower index of the array bracket
!                 0 if val < arr(1), N if val < arr(N))
!
implicit none
real(rp), dimension(:) :: arr
real(rp) :: val
integer :: low, mid, high, N

! Size of array
N = size(arr)

! Check if value is outside bounds
if ( val < arr(1) ) then
  low = 0
  return
end if
if ( val > arr(N) ) then
  low = N
  return
end if

! Otherwise perform bisection
low = 1
high = N
do while (high - low > 1)
  mid = (low + high) / 2
  if ( arr(mid) > val ) then
    high = mid
  elseif ( arr(mid) < val ) then
    low = mid
  else
    low = mid
    return
  endif
end do

end function binary_search

end module turbine_indicator
#endif TURBINES
