!*******************************************************************************
subroutine inflow_cond
!*******************************************************************************
!
! Enforces prescribed inflow condition based on an uniform inflow velocity
!
use param, only: rp,inflow_velocity, nx, ny, nz, fringe_factor
use sim_param, only : u, v, w
#ifdef SCALAR 
use scalars_param, only: theta
#endif SCALAR
implicit none
integer :: i, i_w
integer :: istart, istart_w
integer :: iplateau
integer :: iend, iend_w
real(rp) :: alpha
real(kind=16) :: beta
#ifdef CPS_Y
integer :: istart_yw, iend_yw
integer :: istart_y, iplateau_y, iend_y
real(kind=16) ::  beta_y
real(rp) :: alpha_y 
#endif CPS_Y

! these may be out of 1, ..., nx
call fringe_init(istart, iplateau, iend)

! wrapped versions
iend_w = modulo(iend-1, nx) + 1
istart_w = modulo(istart-1, nx) + 1

! Set end of domain
u(iend_w,:,:) = inflow_velocity
v(iend_w,:,:) = 0._rp
w(iend_w,:,:) = 0._rp
#ifdef SCALAR 
theta(iend_w,:,:) = 1._rp 
#endif SCALAR

! skip istart since we know vel at istart, iend already
do i = istart+1, iend-1

  i_w = modulo(i-1, nx) + 1

  call fringe_weighting(i, istart, iplateau, beta)
  ! beta has to be divided by the fringe_factor, as the fringe_factor is not needed for uniform inflow cases, but only when CPS is activated.
  beta = beta / fringe_factor
  alpha = 1.0_rp - beta

  u(i_w,1:ny,1:nz) = alpha * u(i_w,1:ny,1:nz) + beta * inflow_velocity
  v(i_w,1:ny,1:nz) = alpha * v(i_w,1:ny,1:nz)
  w(i_w,1:ny,1:nz) = alpha * w(i_w,1:ny,1:nz)

#ifdef SCALAR
  theta(i_w,1:ny,1:nz) = alpha * theta(i_w,1:ny,1:nz)
#endif SCALAR
end do

#ifdef CPS_Y
! these may be out of 1, ..., nx
call fringe_init_y(istart_y, iplateau_y, iend_y)

! wrapped versions
iend_yw = modulo(iend_y-1, ny) + 1
istart_yw = modulo(istart_y-1, ny) + 1

! Set end of domain
u(:,iend_yw,:) = 0._rp  
v(:,iend_yw,:) = 0._rp
w(:,iend_yw,:) = 0._rp

! skip istart since we know vel at istart, iend already
do i = istart_y+1, iend_y-1

  i_w = modulo(i-1, ny) + 1

  call fringe_weighting(i, istart_y, iplateau_y, beta_y)
  alpha_y = 1.0_rp - beta_y

  u(1:nx,i_w,1:nz) = alpha_y * u(1:nx,i_w,1:nz)
  v(1:nx,i_w,1:nz) = alpha_y * v(1:nx,i_w,1:nz)
  w(1:nx,i_w,1:nz) = alpha_y * w(1:nx,i_w,1:nz)

#ifdef SCALAR
  theta(1:nx,i_w,1:nz) = alpha_y * theta(1:nx,i_w,1:nz)
#endif SCALAR
end do
#endif CPS_Y

end subroutine inflow_cond


!*******************************************************************************
subroutine project
!*******************************************************************************
!
! Add pressure gradients to RHS variables (for next time step)
! Update u, v, w 
!
use param, only: rp,nx,ny,nz,inflow,dt,tadv1,kmin
#ifdef LVLSET 
use sim_param, only: u,v,w,RHSx,RHSy,RHSz,dpdx,dpdy,dpdz
#else
use sim_param, only: u,v,w,RHSx,RHSy,RHSz,dpdx=>txx,dpdy=>txy,dpdz=>tyy
#endif LVLSET
use mpi_defs, only: mpi_sync_downup
#ifndef SCALAR 
#ifdef CPS
#ifndef LVLSET  
#ifndef TURBINES
use param, only: ld,jt_total
use param, only: coord, sgs_model
use sim_param, only: RHSx_f,RHSy_f,RHSz_f,dummy1
use sgs_param, only: F_LM,F_MM,F_NN,F_QN
#endif TURBINES
#endif LVLSET
#endif CPS
#endif SCALAR
implicit none
integer :: i, j, k 
real(rp) :: tconst

tconst = -tadv1 * dt

do k=1,nz-1
do j=1,ny
do i=1,nx
  RHSx(i,j,k) = RHSx(i,j,k) - dpdx(i,j,k)
  RHSy(i,j,k) = RHSy(i,j,k) - dpdy(i,j,k)
  RHSz(i,j,k) = RHSz(i,j,k) - dpdz(i,j,k)

  u(i,j,k) = u(i,j,k) + tconst*dpdx(i,j,k)
  v(i,j,k) = v(i,j,k) + tconst*dpdy(i,j,k)
enddo
enddo
enddo

! k_global=0 is zero (boundary condition)
do k= kmin,nz-1
do j= 1,ny
do i= 1,nx
  w(i,j,k) = w(i,j,k) + tconst*dpdz(i,j,k)
enddo
enddo
enddo

#ifndef CPS
if(inflow) call inflow_cond
#endif CPS

! Exchange ghost node information (since coords overlap)
call mpi_sync_real_array(u,mpi_sync_downup)
call mpi_sync_real_array(v,mpi_sync_downup)
call mpi_sync_real_array(w,mpi_sync_downup)

! Shift the solution in the precursor domain to eliminate streak effect
#ifdef CPS
#ifndef LVLSET  
#ifndef WINDBREAKS
#ifndef SCALAR 
#ifndef TURBINES
#ifndef ATM
if (modulo (jt_total, 1000) == 0) then
  if (coord == 0) then
    write(*,*) '--------------------red shift-------------------'
  endif

  dummy1(1:ld,1:ny,0:nz  )=     u(1:ld,1:ny  ,0:nz  )
  u     (1:ld,2:ny,0:nz  )=dummy1(1:ld,1:ny-1,0:nz  )
  u     (1:ld,1   ,0:nz  )=dummy1(1:ld,ny    ,0:nz  )

  dummy1(1:ld,1:ny,0:nz  )=     v(1:ld,1:ny  ,0:nz  )
  v     (1:ld,2:ny,0:nz  )=dummy1(1:ld,1:ny-1,0:nz  )
  v     (1:ld,1   ,0:nz  )=dummy1(1:ld,ny    ,0:nz  )

  dummy1(1:ld,1:ny,0:nz  )=     w(1:ld,1:ny  ,0:nz  )
  w     (1:ld,2:ny,0:nz  )=dummy1(1:ld,1:ny-1,0:nz  )
  w     (1:ld,1   ,0:nz  )=dummy1(1:ld,ny    ,0:nz  )

  dummy1(1:ld,1:ny,1:nz-1)=  RHSx(1:ld,1:ny  ,1:nz-1)
  RHSx  (1:ld,2:ny,1:nz-1)=dummy1(1:ld,1:ny-1,1:nz-1)
  RHSx  (1:ld,1   ,1:nz-1)=dummy1(1:ld,ny    ,1:nz-1)

  dummy1(1:ld,1:ny,1:nz-1)=RHSx_f(1:ld,1:ny  ,1:nz-1)
  RHSx_f(1:ld,2:ny,1:nz-1)=dummy1(1:ld,1:ny-1,1:nz-1)
  RHSx_f(1:ld,1   ,1:nz-1)=dummy1(1:ld,ny    ,1:nz-1)

  dummy1(1:ld,1:ny,1:nz-1)=  RHSy(1:ld,1:ny  ,1:nz-1)
  RHSy  (1:ld,2:ny,1:nz-1)=dummy1(1:ld,1:ny-1,1:nz-1)
  RHSy  (1:ld,1   ,1:nz-1)=dummy1(1:ld,ny    ,1:nz-1)

  dummy1(1:ld,1:ny,1:nz-1)=RHSy_f(1:ld,1:ny  ,1:nz-1)
  RHSy_f(1:ld,2:ny,1:nz-1)=dummy1(1:ld,1:ny-1,1:nz-1)
  RHSy_f(1:ld,1   ,1:nz-1)=dummy1(1:ld,ny    ,1:nz-1)

  dummy1(1:ld,1:ny,1:nz-1)=  RHSz(1:ld,1:ny  ,1:nz-1)
  RHSz  (1:ld,2:ny,1:nz-1)=dummy1(1:ld,1:ny-1,1:nz-1)
  RHSz  (1:ld,1   ,1:nz-1)=dummy1(1:ld,ny    ,1:nz-1)

  dummy1(1:ld,1:ny,1:nz-1)=RHSz_f(1:ld,1:ny  ,1:nz-1)
  RHSz_f(1:ld,2:ny,1:nz-1)=dummy1(1:ld,1:ny-1,1:nz-1)
  RHSz_f(1:ld,1   ,1:nz-1)=dummy1(1:ld,ny    ,1:nz-1)

  if(sgs_model==5) then
    dummy1(1:ld,1:ny,0:nz  )=  F_LM(1:ld,1:ny  ,0:nz  )
    F_LM  (1:ld,2:ny,0:nz  )=dummy1(1:ld,1:ny-1,0:nz  )
    F_LM  (1:ld,1   ,0:nz  )=dummy1(1:ld,ny    ,0:nz  )

    dummy1(1:ld,1:ny,0:nz  )=  F_MM(1:ld,1:ny  ,0:nz  )
    F_MM  (1:ld,2:ny,0:nz  )=dummy1(1:ld,1:ny-1,0:nz  )
    F_MM  (1:ld,1   ,0:nz  )=dummy1(1:ld,ny    ,0:nz  )

    dummy1(1:ld,1:ny,0:nz  )=  F_QN(1:ld,1:ny  ,0:nz  )
    F_QN  (1:ld,2:ny,0:nz  )=dummy1(1:ld,1:ny-1,0:nz  )
    F_QN  (1:ld,1   ,0:nz  )=dummy1(1:ld,ny    ,0:nz  )

    dummy1(1:ld,1:ny,0:nz  )=  F_NN(1:ld,1:ny  ,0:nz  )
    F_NN  (1:ld,2:ny,0:nz  )=dummy1(1:ld,1:ny-1,0:nz  )
    F_NN  (1:ld,1   ,0:nz  )=dummy1(1:ld,ny    ,0:nz  )
  endif
endif
#endif ATM
#endif TURBINES
#endif SCALAR
#endif WINDBREAKS
#endif LVLSET
#endif CPS
 
end subroutine project
