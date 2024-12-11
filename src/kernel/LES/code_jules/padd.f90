! puts arrays into larger, zero-padded arrays, automatically zeroes the oddballs
subroutine padd (u_big,u)
use param,only:rp,ld,ld_big,nx,ny,ny2,inxny
implicit none
!  u and u_big are interleaved as complex arrays
real(rp), dimension(ld,ny), intent(in) :: u
real(rp), dimension(ld_big,ny2), intent(out) :: u_big

integer :: ny_h, j_s, j_big_s

ny_h = ny/2

! make sure the big array is zeroed!
u_big(:,:) = 0._rp

! note: split access in an attempt to maintain locality
u_big(:nx,:ny_h) = u(:nx,:ny_h)*inxny

! Compute starting j locations for second transfer
j_s = ny_h + 2
j_big_s = ny2 - ny_h + 2

u_big(:nx,j_big_s:ny2) = u(:nx,j_s:ny)*inxny
end subroutine padd
