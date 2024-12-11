#ifdef ADMR
!*******************************************************************************
subroutine admr_nodes_func
!*******************************************************************************
use param, only: dx,dy,dz,nz,pi,mpi_rp,ierr,rp,coord
use mpi
use admr_mod
implicit none
real(rp) :: smearing, dist_3d
!real(rp) :: rx,ry,rz,delx,dely,delz,del
integer :: count_i,i,j,k,i_ring,i_rseg
if (admr_in_proc) then
 do inot=1,nturb
  ! Calculation of the boundaries:
  admr(inot)%i_smear = CEILING( ( admr(inot)%rr * ABS( admr(inot)%roty(1) ) +       &
                                  2._rp*dx_dim ) / dx_dim )
  admr(inot)%j_smear = CEILING( ( admr(inot)%rr * ABS( admr(inot)%roty(2) ) +       &
                                  2._rp*dx_dim ) / dy_dim )

  admr(inot)%num_nodes=0
  count_i = 0    !index count - used for writing to array "nodes"
  do i = admr(inot)%i_hub - admr(inot)%i_smear,admr(inot)%i_hub + admr(inot)%i_smear
   do j = admr(inot)%j_hub - admr(inot)%j_smear,admr(inot)%j_hub + admr(inot)%j_smear
    do k = admr(inot)%min_k,admr(inot)%max_k
     do i_ring = 1, admr(inot)%nrings
      do i_rseg = 1, admr(inot)%nsegs(i_ring)

       dist_3d = ( float(i) * dx_dim - rbx(i_ring,i_rseg) )**2 + &
                 ( float(j) * dy_dim - rby(i_ring,i_rseg) )**2 + &
                 ( float(k) * dz_dim - rbz(i_ring,i_rseg) )**2

       if ( dist_3d <= (pro_radius*dx_dim)**2 ) then
        smearing = exp(-(dist_3d/(eps_kernel*dx_dim)**2))/(((eps_kernel*dx_dim)**3._rp)*(pi**1.5_rp)) 
       ! Alternative: 
       ! Smoothing with delta-function (Peskin 2002, The immersed
       ! boundary method, Cambridge University Press)
       !
       ! If this is used, comment out the previous two lines
       !
       ! distance between the current grid point and each rotor area segment:
       !rx = (float(i) * dx_dim - rbx(i_ring,i_rseg)) / dx_dim
       !ry = (float(j) * dy_dim - rby(i_ring,i_rseg)) / dy_dim
       !rz = (float(k) * dz_dim - rbz(i_ring,i_rseg)) / dz_dim
       !if (abs(rx)<= 2._rp .and. abs(ry)<= 2._rp .and. abs(rz) <= 2._rp) then
       ! delx = 0.25_rp * (1._rp + cos(pi*rx/2._rp))
       ! dely = 0.25_rp * (1._rp + cos(pi*ry/2._rp))
       ! delz = 0.25_rp * (1._rp + cos(pi*rz/2._rp))
       ! smearing = (delx*dely*delz)/(dx_dim*dy_dim*dz_dim)

        ! Store this point as admr point
        if (smearing .ne. 0._rp) then
        count_i = count_i + 1
        admr(inot)%smear(count_i) = smearing
        admr(inot)%nodes(1,count_i) = i
        admr(inot)%nodes(2,count_i) = j
        admr(inot)%nodes(3,count_i) = k-coord*(nz-1)!local
        admr(inot)%ringseg(1,count_i) = i_ring
        admr(inot)%ringseg(2,count_i) = i_rseg
        endif
       endif
      enddo ! end i_rseg
     enddo ! end i_rind
    enddo ! end k
   enddo ! end j
  enddo ! end i
 admr(inot)%num_nodes = count_i
 if ((admr(inot)%num_nodes).gt.admr_nodes) write(*,*) 'Memory problem in turbine routine'
 end do ! end inot
endif ! end if(admr_in_proc)
end subroutine admr_nodes_func

#endif ADMR
