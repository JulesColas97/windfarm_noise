!*******************************************************************************
subroutine divstress_uv
!*******************************************************************************
use sim_param, only: u_dim,txx,txy,txz,tyy,tyz,RHSx,RHSy
use sim_param, only: dtxdx =>dummy1,dtydy =>dummy2,dtzdz =>dummy3
use sim_param, only: dtxdx2=>dummy4,dtydy2=>dummy1,dtzdz2=>dummy2
use param, only: z_i,rp,nz,ld,ny,mean_p_force
#ifdef CORIOLIS
use sim_param, only: u, v, w, vgMag
use param, only: adjust_wind_angle, ug, vg, coriol, coriol_x, coriol_y, ny, pi
use sim_param, only: omega_eff, alpha_wind_old, alpha_wind
use grid_defs, only: gridzw,gridz
#ifdef CPS
use cps_mod, only: color,RED,BLUE,interComm
use param, only: mpi_rp, rank_of_coord, coord, status, ierr
#endif CPS
#endif CORIOLIS
#ifdef BAROCLINIC
use param, only: ug_delta, vg_delta, bar_start, bar_end, geo_force_x, geo_force_y
use sim_param, only: geo_force_x_rot, geo_force_y_rot
#endif BAROCLINIC
#ifdef SCALAR
#ifndef CORIOLIS
use sim_param, only: u, v
#endif CORIOLIS
use param, only: damping_x, ubc,nx,inxny
use scalars_param, only: sponge, sponge_x 
#endif SCALAR
implicit none
integer :: i,j,k,dz1,dz2,dzfact

call ddx  (txx(1:ld,1:ny,1:nz-1),dtxdx (1:ld,1:ny,0:nz))
call ddxy (txy(1:ld,1:ny,1:nz-1),dtxdx2(1:ld,1:ny,0:nz),dtydy(1:ld,1:ny,0:nz))
call ddz_w(txz(1:ld,1:ny,0:nz  ),dtzdz (1:ld,1:ny,0:nz))

RHSx(:,:,1:nz-1)=-RHSx(:,:,1:nz-1)-(dtxdx(:,:,1:nz-1)+dtydy(:,:,1:nz-1)+dtzdz(:,:,1:nz-1))+mean_p_force

call ddy  (tyy(1:ld,1:ny,1:nz-1),dtydy2(1:ld,1:ny,0:nz))
call ddz_w(tyz(1:ld,1:ny,0:nz)  ,dtzdz2(1:ld,1:ny,0:nz))

RHSy(:,:,1:nz-1)=-RHSy(:,:,1:nz-1)-(dtxdx2(:,:,1:nz-1)+dtydy2(:,:,1:nz-1)+dtzdz2(:,:,1:nz-1))

#ifdef CORIOLIS

#ifdef CPS 
if(color==BLUE) then
  call mpi_recv(alpha_wind,1,mpi_rp,rank_of_coord(coord),21,interComm,status,ierr)
  call mpi_recv(omega_eff,1,mpi_rp,rank_of_coord(coord),22,interComm,status,ierr)
elseif(color==RED) then
  call mpi_send(alpha_wind,1,mpi_rp,rank_of_coord(coord),21,interComm,ierr)
  call mpi_send(omega_eff,1,mpi_rp,rank_of_coord(coord),22,interComm,ierr)
endif 
#endif CPS

#ifdef BAROCLINIC
if(adjust_wind_angle) then 
  alpha_wind_old = alpha_wind 
  do k = 1,nz-1
    RHSx(:, :, k) = RHSx(:, :, k) + coriol * v(:, :, k) - coriol * geo_force_y_rot(k) + 2._rp * omega_eff * v(:, :, k) - coriol_y * 0.5_rp * (w(:,:,k) + w(:,:,k+1))
    RHSy(:, :, k) = RHSy(:, :, k) - coriol * u(:, :, k) + coriol * geo_force_x_rot(k) - 2._rp * omega_eff * u(:, :, k) + coriol_x * 0.5_rp * (w(:,:,k) + w(:,:,k+1))
  end do
else
  do k = 1,nz-1
    RHSx(:, :, k) = RHSx(:, :, k) + coriol * v(:, :, k) - coriol * geo_force_y(k) - coriol_y * 0.5_rp * (w(:,:,k) + w(:,:,k+1))
    RHSy(:, :, k) = RHSy(:, :, k) - coriol * u(:, :, k) + coriol * geo_force_x(k) + coriol_x * 0.5_rp * (w(:,:,k) + w(:,:,k+1))
  enddo
endif

!else meaning 'not baroclinic'
#else
if(adjust_wind_angle) then
  alpha_wind_old = alpha_wind
  RHSx(:, :,1:nz-1) = RHSx(:, :,1:nz-1) + coriol * v(:, :,1:nz-1) - coriol * vgMag * sind(180.0 * alpha_wind/pi) + 2._rp * omega_eff * v(:, :,1:nz-1) - coriol_y * 0.5_rp * (w(:,:,1:nz-1) + w(:,:,2:nz))
  RHSy(:, :,1:nz-1) = RHSy(:, :,1:nz-1) - coriol * u(:, :,1:nz-1) + coriol * vgMag * cosd(180.0 * alpha_wind/pi) - 2._rp * omega_eff * u(:, :,1:nz-1) + coriol_x * 0.5_rp * (w(:,:,1:nz-1) + w(:,:,2:nz))
else
  RHSx(:, :,1:nz-1) = RHSx(:, :,1:nz-1) + coriol * v(:, :,1:nz-1) - coriol * vg - coriol_y * 0.5_rp * (w(:,:,1:nz-1) + w(:,:,2:nz))
  RHSy(:, :,1:nz-1) = RHSy(:, :,1:nz-1) - coriol * u(:, :,1:nz-1) + coriol * ug + coriol_x * 0.5_rp * (w(:,:,1:nz-1) + w(:,:,2:nz))
endif
#endif BAROCLINIC

#endif CORIOLIS


#ifdef SCALAR
!If damping layer is activated
if(ubc == 1) then

    ! Srinidhi : Feb 20, 2022
    ! The velocity above the damping layer is slowly relaxed to the geostrophic wind speed by this implementation.  
    ! This is a good approximation, however sometimes (for mountain wave simulations) one may want to relax 
    ! the values to some other standard values 
#ifdef BAROCLINIC
  if(adjust_wind_angle) then 
    do k=1, nz-1
      RHSx(:,:,k) = RHSx(:,:,k) - 0.5_rp* (sponge(k)+sponge(k+1)) * (u(:,:,k)-geo_force_x_rot(k))
      RHSy(:,:,k) = RHSy(:,:,k) - 0.5_rp* (sponge(k)+sponge(k+1)) * (v(:,:,k)-geo_force_y_rot(k))
    end do
  else 
    do k=1, nz-1
      RHSx(:,:,k) = RHSx(:,:,k) - 0.5_rp* (sponge(k)+sponge(k+1)) * (u(:,:,k)-geo_force_x(k))
      RHSy(:,:,k) = RHSy(:,:,k) - 0.5_rp* (sponge(k)+sponge(k+1)) * (v(:,:,k)-geo_force_y(k))
    enddo  
  endif 

#else
  if(adjust_wind_angle) then 
    do k=1, nz-1
      RHSx(:,:,k) = RHSx(:,:,k) - 0.5_rp* (sponge(k)+sponge(k+1)) * (u(:,:,k)-vgMag * cosd(180.0 * alpha_wind/pi))
      RHSy(:,:,k) = RHSy(:,:,k) - 0.5_rp* (sponge(k)+sponge(k+1)) * (v(:,:,k)-vgMag * sind(180.0 * alpha_wind/pi))
    end do
  else 
    do k=1, nz-1
      RHSx(:,:,k) = RHSx(:,:,k) - 0.5_rp* (sponge(k)+sponge(k+1)) * (u(:,:,k)-ug)
      RHSy(:,:,k) = RHSy(:,:,k) - 0.5_rp* (sponge(k)+sponge(k+1)) * (v(:,:,k)-vg)
    end do
  endif
#endif BAROCLINIC

  !For periodic runs with mountain waves
  !Switching on damping_x essentially removes all the wave perturbations in the
  !streamwise direction and removes all the gravity waves 
  !Present implementation can only used for trial runs as this also damps out
  !turbulence 
  !Srinidhi: I am also implementing concurrent precursor version of this which
  !would preserve even the turbulence 
  !Srinidhi, Feb 20, 2022: This basically relaxes the inlet velocity to the planar averaged value. 
  !For mountain wave simulations I also implemented an upstream fringe method 
  !That still needs to be checked out 
  if(damping_x) then
  do k=1, nz-1
   do j=1,ny
    RHSx(:,j,k) = RHSx(:,j,k) - sponge_x(:) * (u(:,j,k)-sum(u(1:nx,1:ny,k))*inxny)
    RHSy(:,j,k) = RHSy(:,j,k) - sponge_x(:) * (v(:,j,k)-sum(v(1:nx,1:ny,k))*inxny)
   enddo
  end do
  endif
endif
#endif SCALAR

end subroutine divstress_uv
