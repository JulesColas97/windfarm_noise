!*******************************************************************************
subroutine press_stag_init
!*******************************************************************************
!
! Determines a,b,c for tridiagonal solver
!
use param, only: rp,lh,ny,nz,kmin
use sim_param, only: a,b,c
use fft, only: kx,ky
use grid_defs, only: gridzw
implicit none
real(rp) :: const3, dz1, dz2, dz3 
integer :: jx,jy,k


! JDA dissertation, eqn(2.85) a,b,c=coefficients

b(:,:,1)=-1._rp
c(1)=1._rp
a(nz+1) =-1._rp
b(:,:,nz+1) = 1._rp

do k=kmin,nz
 dz1 = 0.5_rp * (gridzw(k) - gridzw(k-2)) 
 dz2 = 0.5_rp * (gridzw(k+1) - gridzw(k-1)) 
 dz3 = (gridzw(k) - gridzw(k-1)) 
 const3 = (dz1 + dz2) / (dz1 * dz2 * dz3) 

 a(k)= 1._rp/(dz1 * dz3) 
 c(k)= 1._rp/(dz2 * dz3) 

 do jy=1,ny
  do jx=1,lh-1
   b(jx,jy,k)=-(kx(jx,jy)**2+ky(jx,jy)**2+const3)
  end do
 enddo

end do


end subroutine press_stag_init


!*******************************************************************************
subroutine press_stag
!*******************************************************************************
!
! Solves pressure correction; provides p, dpdx, dpdy, dpdz
!
use param, only: rp,ld,ny,mpi_rp,localComm,status,ierr,inxny
use param, only: nx,nz,dt,tadv1,coord,lh,up,down,kmin
use sim_param, only: pres 
use grid_defs, only: gridzw, gridz 
#ifdef LVLSET 
use sim_param, only: u,v,w,divtz,p=>tzz,dpdx,dpdy,dpdz
#else
use sim_param, only :u,v,w,divtz,p=>tzz,dpdx=>txx,dpdy=>txy,dpdz=>tyy
#endif LVLSET
use sim_param,only: rH_x=>dummy1,rH_y=>dummy2,rH_z=>dummy3,RHS_col=>dummy4
use fft, only: kx,ky,forw,back
use mpi_defs, only: mpi_sync_down,mpi_sync_up
implicit none
real(rp) :: const2, dz1 , dzfact
real(rp), dimension(2) :: aH_x,aH_y ! Used to emulate complex scalar
integer :: i,ir,ii,j,k,km

! Zero out dummy array for pressure calculations
rH_x(:,:,nz)=0._rp
rH_y(:,:,nz)=0._rp
rH_z(:,:,nz)=0._rp
RHS_col     =0._rp

! Specifiy cached constants
const2=inxny/tadv1/dt

do k=1,nz-1
  rH_x(:,:,k) = const2*u(:,:,k)
  rH_y(:,:,k) = const2*v(:,:,k)
  rH_z(:,:,k) = const2*w(:,:,k)
  call dfftw_execute_dft_r2c(forw,rH_x(:,:,k),rH_x(:,:,k))
  call dfftw_execute_dft_r2c(forw,rH_y(:,:,k),rH_y(:,:,k))
  call dfftw_execute_dft_r2c(forw,rH_z(:,:,k),rH_z(:,:,k))
enddo

if (coord.eq.0) then
  p(:,:,0)=0._rp
  divtz(:, :)=inxny*divtz(:,:) ! Reuse divtz space
  call dfftw_execute_dft_r2c(forw,divtz,divtz)
  RHS_col(:,:,0) =-0.5_rp * (gridzw(2) - gridzw(0))*divtz(:,:)
end if

! send to nz
call mpi_sync_real_array(rH_z,mpi_sync_down)

do k=kmin,nz
  km=k-1
  dzfact = 1.0_rp / (gridzw(k) - gridzw(km))
  do j=1,ny
  do i= 1,lh-1
    ii=2*i  ! imaginary index
    ir=ii-1 ! real index

    aH_x(1) = -rH_x(ii,j,k-1) * kx(i,j)
    aH_x(2) =  rH_x(ir,j,k-1) * kx(i,j)
    aH_y(1) = -rH_y(ii,j,k-1) * ky(i,j)
    aH_y(2) =  rH_y(ir,j,k-1) * ky(i,j)

    RHS_col(ir:ii,j,km)=aH_x+aH_y+(rH_z(ir:ii,j,k)-rH_z(ir:ii,j,km))*dzfact
  end do
  end do
end do

! Tridiagonal solver for non-zero wave numbers
call tridag

! Begin zero-wavenumber solution (only real part is solved; imaginary part is zero)
call mpi_recv(p(1,1,1),1,mpi_rp,down,8,localComm,status,ierr)

if (coord == 0) then
  p(1,1,0)=0._rp
  p(1,1,1)=-(gridz(1) - gridz(0))*divtz(1,1)
end if

! JDA dissertation, eqn(2.88)
do k=2,nz
  dz1 = 0.5_rp * (gridzw(k+1) - gridzw(k-1)) 
  p(1,1,k)=p(1,1,k-1)+rH_z(1,1,k)*dz1
end do

call mpi_send(p(1,1,nz),1,mpi_rp,up,8,localComm,ierr)

! zero out imaginary part 1st wavenumber
p(2,1,:)=0._rp
! End zero-wavenumber solution

! zero the nyquist freqs
p(ld-1:ld,:,:)=0._rp
p(:,ny/2+1,:) =0._rp

!Calculate horizontal pressure gradients
do k=1,nz-1
  do j=1,ny
  do i=1,lh-1
    ii = 2*i
    ir = ii-1

    dpdx(ir,j,k) = -p(ii,j,k) * kx(i,j)
    dpdx(ii,j,k) =  p(ir,j,k) * kx(i,j)
    dpdy(ir,j,k) = -p(ii,j,k) * ky(i,j)
    dpdy(ii,j,k) =  p(ir,j,k) * ky(i,j)
  end do
  end do

  call dfftw_execute_dft_c2r(back,dpdx(:,:,k), dpdx(:,:,k))
  call dfftw_execute_dft_c2r(back,dpdy(:,:,k), dpdy(:,:,k))
end do

! send to 0
call mpi_sync_real_array(p,mpi_sync_up)

do k=0,nz-1
  call dfftw_execute_dft_c2r(back,p(:,:,k),p(:,:,k))
enddo

! Calculate dpdz: p has additional level at z=-dz/2 for this derivative
do k = 1, nz-1 
 dz1 = 1._rp/(0.5_rp * (gridzw(k+1) - gridzw(k-1))) 
 dpdz(1:nx,1:ny,k)=(p(1:nx,1:ny,k)-p(1:nx,1:ny,k-1))*dz1
enddo

! Pressure copy 
pres = p 

end subroutine press_stag
