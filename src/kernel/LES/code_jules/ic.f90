!*******************************************************************************
subroutine ic
!*******************************************************************************
!
! This routine generates a new random initial velocity field
!
use mpi, only: mpi_comm_world
use param, only: rp,nx,ny,nz,path,coord,nproc,ierr
use param, only: z_i,vonk,zo
#ifndef CPS
use param, only: inflow,inflow_velocity
#endif
use sim_param, only :u,v,w
#ifdef CORIOLIS
use param, only: ug,vg,pi, adjust_wind_angle
use sim_param, only: vgMag, alpha_wind
#endif  
#ifdef BAROCLINIC
use param, only: ug_delta, vg_delta, bar_start, bar_end, geo_force_x, geo_force_y
#endif BAROCLINIC
use mpi_defs, only: mpi_sync_downup
implicit none
logical :: exst,exst2

! Check whether there data from a previous simulation. If this is the case stop
! the simulation as the simulation is probably configured incorrectly.
inquire(file=trim(path)//'continuation/total_time.dat',exist=exst)
inquire(file=trim(path)//'continuation/tavg_total_time.dat',exist=exst2)
if (exst .or. exst2) then
  if(coord == 0 ) write(*,'(a)') 'Trying to generate a new initial velocity field while data from previous simulation is found'
  if(coord == 0 ) write(*,'(a)') 'Termintate the simulation and check configuration files'
  call mpi_barrier(mpi_comm_world, ierr)
  call mpi_finalize(ierr)
  stop
endif

! Check whether there data from a previous simulation. If this is the case stop
! the simulation as the simulation is probably configured incorrectly.
#ifdef CPS
call boundary_layer_ic
#else
  
if (inflow) then  !--no turbulence
  ! Generate uniform initial condition
  u = inflow_velocity
  v = 0._rp
  w = 0._rp
else
  call boundary_layer_ic
end if
#endif 

contains

!*******************************************************************************
subroutine boundary_layer_ic
!*******************************************************************************
!
! This routine generates a new random initial velocity field assuming a 
! logarithmic velocity profile.
!
use param, only: u_star
use grid_defs, only: gridz 
implicit none
real(rp),dimension(nz) :: ubar,vbar
real(rp) :: rms,noise,z,arg,arg2
integer :: seed,i,j,k,k_abs
interface
  function ran3(idum)
    integer(4) idum
    real(8) :: ran3
  end function ran3
end interface

if(coord == 0) write(*,*) '------> Creating initial velocity field'

! Set logarithmic velocity profile
do k=1,nz

  arg2=gridz(k)/zo
  arg = (1._rp/vonk)*log(arg2)
  arg = arg*u_star

  ubar(k) = arg
  vbar(k) = 0._rp
#ifdef CORIOLIS
#ifdef BAROCLINIC
  ! Linear geostrophic shear with height (corresponding to a height and time invariant horizontal temperature gradient)
  ! as previously used by Momen (2018, Journal of Atmospheric Sciences) and Sorbjan (2004, Bound.-Layer Meteor.)
  if( gridz(k) <= bar_start) then
    geo_force_x(k) = ug
    geo_force_y(k) = vg
  elseif(bar_start < gridz(k) .and. gridz(k) <= bar_end) then
    geo_force_x(k) = ug + ( ug_delta/(bar_end-bar_start) ) * ( gridz(k)-bar_start )
    geo_force_y(k) = vg + ( vg_delta/(bar_end-bar_start) ) * ( gridz(k)-bar_start )
  else
    geo_force_x(k) = ug + ug_delta
    geo_force_y(k) = vg + vg_delta
  endif

  if(adjust_wind_angle) then
    ubar(k)=geo_force_x(k) * cosd(180.0 * alpha_wind/pi) - geo_force_y(k) * sind(180.0 * alpha_wind/pi)
    vbar(k)=geo_force_x(k) * sind(180.0 * alpha_wind/pi) + geo_force_y(k) * cosd(180.0 * alpha_wind/pi)
  else
    ubar(k)=geo_force_x(k)
    vbar(k)=geo_force_y(k)
  endif

#else
  if(adjust_wind_angle) then
    ubar(k)=vgMag * cosd(180.0 * alpha_wind/pi)
    vbar(k)=vgMag * sind(180.0 * alpha_wind/pi)
  else
    ubar(k)=ug
    vbar(k)=vg
  endif
#endif BAROCLINIC  
#endif CORIOLIS
end do

! Add random noise to the logarithmic velocity profile
rms = 3._rp
do k=1,nz
  k_abs = coord * (nz-1) + k 
  z = gridz(k) * z_i    !dimensions in meters
  seed = -80 - k_abs  ! Make consistent for MPI
  do j=1,ny
  do i=1,nx
    ! Ran3 returns uniform random variable between 0 and 1
    if (z.le.z_i) then
      noise=rms/.289_rp*(ran3(seed)-0.5_rp)
      u(i,j,k)=noise*(1._rp-z/z_i)+ubar(k)
      noise=rms/.289_rp*(ran3(seed)-0.5_rp)
      v(i,j,k)=noise*(1._rp-z/z_i)+vbar(k)
      noise=rms/.289_rp*(ran3(seed)-0.5_rp)
      w(i,j,k)=noise*(1._rp-z/z_i)
    else
      u(i,j,k)=ubar(k)
      v(i,j,k)=vbar(k)
      w(i,j,k)=0.0_rp 
    end if
  end do
  end do
end do
 
! Set boundary conditions at the bottom
if (coord == 0) then
  w(1:nx,1:ny,1) = 0._rp
end if

! Set boundary conditions at the top
if (coord == nproc-1) then
  w(1:nx,1:ny,nz) = 0._rp
  u(1:nx,1:ny,nz) = u(1:nx,1:ny,nz-1)
  v(1:nx,1:ny,nz) = v(1:nx,1:ny,nz-1)
end if

! Syncronize data to make sure ghost layers are updated
call mpi_sync_real_array(u,mpi_sync_downup)
call mpi_sync_real_array(v,mpi_sync_downup)
call mpi_sync_real_array(w,mpi_sync_downup)

end subroutine boundary_layer_ic
end subroutine ic
