subroutine muli( a, a_c ,b)
!  This function emulates the multiplication of two complex 2D array by
!  emulating the input real array (a) as a complex type. This subroutine
!  ignores the real part of a_c (e.g. would use this when real(a_c) = 0)
!
!  Input:
!  
!    a (real,size(nx_r,ny))     - input real array
!    a_c (real,size(nx_c,ny))   - input imaginary part of complex array
!
!  Output:
!
!    b (real,size(nx_r,ny))     - output real array
!
!  Note: nx_c must be nx_r/2
use param, only : rp,ld,ny,lh
implicit none
real(rp),dimension(1:ld,1:ny),intent(in) :: a
real(rp),dimension(1:lh,1:ny),intent(in) :: a_c
real(rp),dimension(1:ld,1:ny),intent(out) :: b
real(rp) :: cache !  Cached variables
integer :: i,j,ii,ir

!  Emulate complex multiplication
do j=1,ny !  Using outer loop to get contiguous memory access
 do i=1,lh
    !  Real and imaginary indicies of a
    ii = 2*i
    ir = ii-1

    !  Perform multiplication (cache data to ensure sequential access)
    cache   =   a(ir,j)*a_c(i,j) 
    b(ir,j) = - a(ii,j)*a_c(i,j)
    b(ii,j) =  cache
  enddo
enddo

end subroutine muli 


subroutine mulr( a, a_c)
!  This function emulates the multiplication of two complex 2D array by
!  emulating the input real array (a) as a complex type. This subroutine
!  ignores the imaginary part of a_c (e.g. would use this when imag(a_c)
!  = 0).
!
!  Input/output:
!    a (real,size(nx_r,ny))     - input/output real array
!    a_c (real,size(nx_c,ny))   - input real part of complex array
use param, only : rp,ld,ny,lh
implicit none
real(rp),dimension(1:ld,1:ny),intent(inout) :: a
real(rp),dimension(1:lh,1:ny),intent(in) :: a_c
integer :: i,j,ii,ir

!  Emulate complex multiplication
do j=1,ny !  Using outer loop to get contiguous memory access
  do i=1,lh
    !  Real and imaginary indicies of a
    ii = 2*i
    ir = ii-1

    !  Perform multiplication
    a(ir:ii,j) = a(ir:ii,j)*a_c(i,j)
  enddo
enddo

end subroutine mulr
