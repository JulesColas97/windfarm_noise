subroutine fringe_init( istart, iplateau, iend )
! 
! Sets the beginning, ending and plateau index within the fringe
! region. Provides a common routine to do this.
!
use param, only : rp,nx,fringe_region_end,fringe_region_len
implicit none

integer, intent(out) :: istart, iplateau, iend

iend = floor (fringe_region_end * nx )
iplateau = floor (( fringe_region_end - fringe_region_len / 4 ) * nx)
istart = floor ((fringe_region_end - fringe_region_len) * nx )

end subroutine fringe_init

#ifdef CPS_Y 
subroutine fringe_init_y( istart_y, iplateau_y, iend_y )
! 
! Sets the beginning, ending and plateau index within the fringe
! region. Provides a common routine to do this.

use param, only : rp,ny,fringe_region_end_y,fringe_region_len_y 
implicit none

integer, intent(out) :: istart_y, iplateau_y, iend_y

iend_y = floor (fringe_region_end_y * ny )
iplateau_y = floor (( fringe_region_end_y - fringe_region_len_y / 4 ) * ny)
istart_y = floor ((fringe_region_end_y - fringe_region_len_y) * ny )

end subroutine fringe_init_y 
#endif CPS_Y 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine fringe_weighting( i, istart, iend ,w ) 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! Provides the weigting for the sampled (desired) velocity field in the
! fringe region. The result from this function (w) is used to compute
! the fringe region velocity as
!
!   u_fringe = w * u_sample + (1-w)*u
!
! where u_fringe is the fringe velocity, u_sample is the sampled
! velocity and u is the unmodulated velocity in the fringe region right
! after the projection step.
!
use param, only : rp,pi, fringe_factor 
implicit none
integer, intent(in) :: i, istart, iend
real(kind=16), intent(out) :: w

! Sine profile with plateau
if ( i > iend ) then 
   w = 1.0Q0
else
   w = 0.5Q0 * ( 1.0Q0 - cos (pi * real (i - istart,kind=16)  &
        / real(iend - istart,kind=16)) )
endif

w = w * fringe_factor 

end subroutine fringe_weighting

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine fringe_weighting_symmetric(i, istart, iend ,w)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use param, only : rp,pi,fringe_factor 
implicit none
integer, intent(in) :: i, istart, iend
real(kind=16), intent(out) :: w
real(rp) :: Delta_Fr

Delta_Fr = real(iend - istart, kind=16) / 3.Q0

if (i < istart + int(Delta_Fr)) then
   w = 0.5Q0 * ( 1.0Q0 - cos(pi * real (i - istart,kind=16)  &
        / Delta_Fr) )
elseif (istart + int(Delta_Fr) <= i .and. &
        i <= istart + int(2.Q0 * Delta_Fr)) then
  w = 1.0Q0
elseif (i > istart + int(2.Q0*Delta_Fr)) then
  w = 0.5Q0 *  cos(pi * (real(i - istart, kind=16) - 2.Q0 * Delta_Fr) &
        / Delta_Fr ) + 0.5Q0
endif

w = w * fringe_factor 

end subroutine fringe_weighting_symmetric
