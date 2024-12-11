!*******************************************************************************
subroutine get_cfl_dt(dt)
!*******************************************************************************
!
! This functions determines the maximum allowable time step based on the CFL 
! number
!
use sim_param, only: u,v,w
use param, only: cfl,rp,dx,dy,nx,ny,nz,mpi_rp,ierr
use grid_defs, only: gridzw 
use mpi, only: mpi_max,mpi_min,mpi_comm_world
!#ifdef ATM
!#ifndef ADMR
!use param, only: dz
!use atm_input_util, only: numberOfTurbines
!use atm_input_util, only: turbineArray
!use atm_input_util, only: turbineModel
!#endif
!#endif ATM
implicit none
real(rp), intent(out) :: dt
real(rp) :: maxu,maxv,maxw,maxi,maxi_buf,dt_buf
real(rp), dimension(1:nz-1) :: dummy 
!#ifdef ATM
!#ifndef ADMR
!real(rp) :: rotSpeed ! rad/s
!real(rp) :: TipRad
!real(rp) :: u_tip,dx_tip,u_tip_buf
!integer :: j
!#endif
!#endif ATM
integer :: i

maxu = maxval( abs(u(1:nx,1:ny,1:nz-1)) )
maxv = maxval( abs(v(1:nx,1:ny,1:nz-1)) )
do i = 1,nz-1 
   dummy(i) = maxval(abs(w(1:nx,1:ny,i)/(gridzw(i+1) - gridzw(i))))
enddo 
maxw = maxval( abs(dummy) ) 

maxi = maxval( (/ maxu/dx,maxv/dy,maxw /) )
call mpi_allreduce(maxi,maxi_buf,1,mpi_rp,mpi_max,mpi_comm_world,ierr)

! dt has the dimension [s/m]
dt    =cfl/maxi_buf
dt_buf= dt

!Srinidhi: Feb 20, 2022. This has been removed. Anja might have the updated files! 

!This is just additional protection in case someone uses very high CFL 
!The blade tip should not move more than 1 grid point in a time step 
!This condition is generally satisfied with the CFL < 0.1 which we generally use 
!However if higher CFL is used, it is possible that the blade tip might move more than one grid point 
!In such a case this is is activated to limit the time step such that it doesn't move too much. 
!#ifdef ATM
!#ifndef ADMR
!u_tip_buf = 0._rp
!do i=1,numberOfTurbines
!  ! RS: right now the rotspeed is constant and known to all processors.
!  ! RS: When rotation speed is changed make sure this is syncronized to all processors.
!  rotSpeed = turbineArray(i) % rotSpeed ! rad/s
!  ! The turbine type ID
!  j = turbineArray(i) % turbineTypeID
!  TipRad   = (turbineModel(j) % TipRad)  ! m
!  ! maximum tip velocity, m/s
!  u_tip_buf = max(u_tip_buf, rotSpeed*TipRad)
!enddo
!call mpi_allreduce(u_tip_buf,u_tip,1,mpi_rp,mpi_max,mpi_comm_world,ierr)

!!Right now the turbine heads do not rotate dynamically with changing wind direction 
!!ALM model assumes the flow goes from West to East for this calculation
!if ((u_tip * dt) > 0.75_rp * min(dy,dz)) then
!   !The tip can only move less than 3/4th of the grid distance 
!   dx_tip = 0.75_rp * min(dy,dz)
!
!   !If the blade tip moves more than the grid size then recalculate the time step
!   !This will recalculate the time step such that the blade tip can only move 0.75 * dy or dz
!   !This would probably never happen if the cfl < 0.1
!   dt_buf = dx_tip / u_tip
! endif 
!#endif
!#endif

!call mpi_allreduce(dt_buf,dt,1,mpi_rp,mpi_min,mpi_comm_world,ierr)
end subroutine get_cfl_dt

!*******************************************************************************
subroutine get_dt_cfl(dt, cfl)
!*******************************************************************************
!
! This functions determines the maximum CFL based on the chosen dt
!
use sim_param, only: u,v,w
use param, only: rp,dx,dy,dz,nx,ny,nz,mpi_rp,ierr
use mpi, only: mpi_max,mpi_comm_world
implicit none
real(rp), intent(in) :: dt
real(rp), intent(inout) :: cfl
real(rp) :: cfl_u,cfl_v,cfl_w
#ifdef MPI
real(rp) :: cfl_buf
#endif

cfl_u = maxval( abs(u(1:nx,1:ny,1:nz-1)) ) / dx
cfl_v = maxval( abs(v(1:nx,1:ny,1:nz-1)) ) / dy
cfl_w = maxval( abs(w(1:nx,1:ny,1:nz-1)) ) / dz

cfl = dt * maxval( (/ cfl_u, cfl_v, cfl_w /) )

#ifdef MPI
  call mpi_allreduce(cfl, cfl_buf, 1, mpi_rp, mpi_max, mpi_comm_world, ierr)
  cfl = cfl_buf
#endif

end subroutine get_dt_cfl


