subroutine unpadd(cc,cc_big)
use param,only:rp,ld,nx,ny,ny2,ld_big,inx2ny2
implicit none
!  cc and cc_big are interleaved as complex arrays
real(rp), dimension( ld, ny ) :: cc
real(rp), dimension( ld_big, ny2 ) :: cc_big

integer :: ny_h, j_s, j_big_s

ny_h = ny/2

cc(:nx,:ny_h) = cc_big(:nx,:ny_h)*inx2ny2

! oddballs
cc(ld-1:ld,:) = 0._rp
cc(:,ny_h+1) = 0._rp

! Compute starting j locations for second transfer
j_s = ny_h + 2
j_big_s = ny2 - ny_h + 2
cc(:nx,j_s:ny) = cc_big(:nx,j_big_s:ny2)*inx2ny2

end subroutine unpadd
