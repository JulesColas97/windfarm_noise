#ifdef ATM
!*******************************************************************************
module atm_base
!*******************************************************************************
!
! This module provides basic functionalities to the actuator turbine model
! Real precision variable and dynamic allocation types are stored here
!
implicit none
save
private
public :: cross_product
public :: interpolate
public :: interpolate_i
public :: rotatePoint

contains

!*******************************************************************************
function interpolate(xp, x, y)
!*******************************************************************************
!
! This function interpolates xp from x and y
!
use param, only: rp
implicit none
real(rp), dimension(:), intent(in) :: x,y
real(rp), intent(in) :: xp
real(rp) :: xa,xb,ya,yb
integer :: i,p
real(rp) :: interpolate

p = size(x)
interpolate = 0._rp

if(xp .le. x(1)) then
  xa = x(1)
  xb = x(2)
  ya = y(1)
  yb = y(2)
  interpolate = ya + (yb-ya)*(xp-xa)/(xb-xa)
elseif(xp .ge. x(p)) then
  xa = x(p)
  xb = x(p-1)
  ya = y(p)
  yb = y(p-1)
  interpolate = ya + (yb-ya)*(xp-xa)/(xb-xa)
else
  do i=2,p
    if((xp .gt. x(i-1)) .and. (xp .le. x(i))) then
      xa = x(i-1)
      xb = x(i)
      ya = y(i-1)
      yb = y(i)
      interpolate = ya + (yb-ya)*(xp-xa)/(xb-xa)
    endif
  enddo
endif

end function interpolate


!*******************************************************************************
integer function interpolate_i(xp, x, y)
!*******************************************************************************
!
! This function interpolates xp from x and y
!
use param, only: rp
implicit none
real(rp), dimension(:), intent(in) :: x
integer, dimension(:), intent(in) :: y
real(rp), intent(in) :: xp
real(rp) :: xa,xb,ya,yb
integer :: i,p

p = size(x)
interpolate_i = 0

if(xp .le. x(1)) then
  interpolate_i = y(1)
elseif(xp .ge. x(p)) then
  interpolate_i = y(p)
else
  do i=2,p
    if((xp .gt. x(i-1)) .and. (xp .le. x(i))) then
      xa = x(i-1)
      xb = x(i)
      ya = real(y(i-1), rp)
      yb = real(y(i), rp)
      interpolate_i = nint(ya + (yb-ya)*(xp-xa)/(xb-xa))
    endif
  enddo
endif

end function interpolate_i


!*******************************************************************************
subroutine rotatePoint(point_in, rotationPoint, axis, angle)
!*******************************************************************************
!
! This subroutine performs rotation of a point with respect to an axis and a
! certain angle
!
use param, only: rp
implicit none
real(rp), dimension(3), intent(inout) :: point_in
real(rp), dimension(3), intent(in) :: rotationPoint
real(rp), dimension(3), intent(in) :: axis
real(rp), intent(in) :: angle
real(rp), dimension(3,3) :: RM 
real(rp), dimension(3) :: point
real(rp) :: cosangle, sinangle

! buffer these values for efficiency
cosangle=cos(angle)
sinangle=sin(angle)

! Rotation matrix tensor
RM(1,1) = axis(1)**2 + (1._rp - axis(1)**2) * cosangle
RM(1,2) = axis(1) * axis(2) * (1._rp - cosangle) - axis(3) * sinangle
RM(1,3) = axis(1) * axis(3) * (1._rp - cosangle) + axis(2) * sinangle
RM(2,1) = axis(1) * axis(2) * (1._rp - cosangle) + axis(3) * sinangle
RM(2,2) = axis(2)**2 + (1._rp - axis(2)**2) * cosangle
RM(2,3) = axis(2) * axis(3) * (1._rp - cosangle) - axis(1) * sinangle
RM(3,1) = axis(1) * axis(3) * (1._rp - cosangle) - axis(2) * sinangle
RM(3,2) = axis(2) * axis(3) * (1._rp - cosangle) + axis(1) * sinangle
RM(3,3) = axis(3)**2 + (1._rp - axis(3)**2) * cosangle

! Rotation matrices make a rotation about the origin, so need to subtract
! rotation point off the point to be rotated
point = point_in - rotationPoint

! Perform rotation (multiplication matrix and vector)
point_in(1) = RM(1,1)*point(1) + RM(1,2)*point(2) + RM(1,3)*point(3)
point_in(2) = RM(2,1)*point(1) + RM(2,2)*point(2) + RM(2,3)*point(3)
point_in(3) = RM(3,1)*point(1) + RM(3,2)*point(2) + RM(3,3)*point(3)

! Return the rotated point to its new location relative to the rotation point
point_in = point_in + rotationPoint

end subroutine rotatePoint


!*******************************************************************************
function cross_product(u, v)
!*******************************************************************************
!
! This function calculates the cross product of 2 vectors
!
use param, only: rp
implicit none
real(rp), intent(in) :: u(3),v(3)
real(rp) :: cross_product(3)

cross_product(1) = u(2)*v(3) - u(3)*v(2)
cross_product(2) = u(3)*v(1) - u(1)*v(3)
cross_product(3) = u(1)*v(2) - u(2)*v(1)

end function cross_product

end module atm_base
#endif ATM
