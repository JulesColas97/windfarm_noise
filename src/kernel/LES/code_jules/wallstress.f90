!*******************************************************************************
subroutine wallstress
!*******************************************************************************
!
! For use with staggered grid LES
! provides txz, tyz (w-nodes) and dudz, dvdz (w-nodes) at jz=1
!
use param, only: rp,lh,ld,nx,ny,vonk,lbc_mom,zo
use grid_defs, only: gridzw, gridz 
#ifdef SCALAR
use sim_param, only: u,v,dudz,dvdz,txz,tyz,phi_m,psi_m
#else
use sim_param, only: u,v,dudz,dvdz,txz,tyz
#endif SCALAR
use sim_param, only: ustar 
use filters, only: G_test
use fft, only: forw,back
implicit none
real(rp), dimension(ld,ny) :: u1,v1
#ifdef SCALAR
real(rp), dimension(ld,ny) :: logprof
#else
real(rp) :: logprof
#endif SCALAR
#ifdef LVLSET
real(rp), parameter :: eps = 100._rp * epsilon (0._rp)
#endif LVLSET
real(rp) :: const,u_avg,factor
integer :: i,j
real(rp) :: u_mean, v_mean, u_mean2, v_mean2, ufluc, vfluc, dz1 

select case (lbc_mom)
case (0) ! Stress free
txz(:,:,1)  = 0._rp
tyz(:,:,1)  = 0._rp
dudz(:,:,1) = 0._rp
dvdz(:,:,1) = 0._rp

! See John D. Albertson's dissertation, eqns (2.46)-(2.52)
! Bou-Zeid, Meneveau, Parlange, "A scale-dependent Lagrangian dynamic model
! for large eddy simulation of complex turbulent flows" (2005)
case (1) ! Wall

u1=u(:,:,1)
v1=v(:,:,1)

! Perform 2delta filtering on u1
call dfftw_execute_dft_r2c(forw,u1,u1)
call mulr(u1(1:ld,1:ny),G_test(1:lh,1:ny))
call dfftw_execute_dft_c2r(back,u1,u1)

! Perform 2delta filtering on v1
call dfftw_execute_dft_r2c(forw,v1,v1)
call mulr(v1(1:ld,1:ny),G_test(1:lh,1:ny))
call dfftw_execute_dft_c2r(back,v1,v1)

factor=1.0/(gridz(1)*vonK)

#ifdef SCALAR
#ifdef LVLSET
logprof = vonk/(log(gridz(1)/zo)-psi_m(:,:,1))
#else
logprof = vonk/(log(gridz(1)/zo)-psi_m(:,:))
#endif LVLSET

do j=1,ny
do i=1,nx
  u_avg     =sqrt(u1(i,j)*u1(i,j)+v1(i,j)*v1(i,j))
  ustar(i,j)=u_avg*logprof(i,j)
  const     =-(ustar(i,j)*ustar(i,j))/u_avg
  txz(i,j,1)=const*u1(i,j)
  tyz(i,j,1)=const*v1(i,j)

#ifdef LVLSET
  const = ustar(i,j)*factor/u_avg * phi_m(i,j,1)
#else
  const = ustar(i,j)*factor/u_avg * phi_m(i,j)
#endif LVLSET
  !this is as in Moeng 84
  dudz(i,j,1)=u(i,j,1)*const
  dvdz(i,j,1)=v(i,j,1)*const

#ifdef LVLSET
  ! ustar can be zero inside the body; set values appropriately
  if (u_avg .le. eps) then
    txz(i,j,1) = 0._rp
    tyz(i,j,1) = 0._rp
    dudz(i,j,1) = 0._rp
    dvdz(i,j,1) = 0._rp
  endif
#endif LVLSET
enddo
enddo

#else
logprof= vonk/(log(gridz(1)/zo))
do j=1,ny
do i=1,nx
  u_avg     =sqrt(u1(i,j)*u1(i,j)+v1(i,j)*v1(i,j))
  ustar(i,j)     =u_avg*logprof
  const     =-(ustar(i,j)*ustar(i,j))/u_avg
  txz(i,j,1)=const*u1(i,j)
  tyz(i,j,1)=const*v1(i,j)

  const      =ustar(i,j)*factor/(u_avg)
  !this is as in Moeng 84
!  dudz(i,j,1)=u(i,j,1)*const
!  dvdz(i,j,1)=v(i,j,1)*const

#ifdef LVLSET
  ! ustar can be zero inside the body; set values appropriately
  if (u_avg .le. eps) then
    txz(i,j,1) = 0._rp
    tyz(i,j,1) = 0._rp
    dudz(i,j,1) = 0._rp
    dvdz(i,j,1) = 0._rp
  endif
#endif LVLSET
enddo
enddo
#endif SCALAR

!This is done to remove log-layer mismatch 
!This is not appropriate for complex terrain cases 
!Essentially the gradient at the first grid point is calculated by the gradient
!above 

u_mean = sum(u(1:nx,1:ny,1))/(nx*ny) 
v_mean = sum(v(1:nx,1:ny,1))/(nx*ny) 
u_mean2 = sum(u(1:nx,1:ny,2))/(nx*ny) 
v_mean2 = sum(v(1:nx,1:ny,2))/(nx*ny) 

do j=1,ny 
 do i=1,nx 
    dz1 = gridzw(2) - gridzw(1) 
    ufluc = ((u(i,j,2) - u_mean2) - (u(i,j,1) - u_mean))/dz1 
    vfluc = ((v(i,j,2) - v_mean2) - (v(i,j,1) - v_mean))/dz1
    dudz(i,j,1) = ufluc 
    dvdz(i,j,1) = vfluc 

    !Log correction Port\'e-Agel (2000) 
    !dudz(i,j,2) = dudz(i,j,2) + ((u_mean2-u_mean)/dz) * (1._rp -
    !(1._rp/log(3.rp)))
    !dvdz(i,j,2) = dvdz(i,j,2) + ((v_mean2-v_mean)/dz) * (1._rp -
    !(1._rp/log(3._rp))) 
  enddo 
enddo 

end select

end subroutine wallstress
