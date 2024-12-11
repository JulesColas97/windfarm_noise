!*******************************************************************************
subroutine divstress_w
!*******************************************************************************
use sim_param, only: divtz,txz,tyz,tzz,RHSz
use sim_param, only: dtxdx =>dummy1,dtydy =>dummy2,dtzdz =>dummy3
use param, only: rp,nz,coord,ld,ny
use param, only: nx,inxny
#ifdef CORIOLIS
use param, only: coriol_x, coriol_y, ug, vg, adjust_wind_angle,pi,z_i
use sim_param, only: u, v, vgMag, alpha_wind
use grid_defs, only: gridzw,gridz
#endif CORIOLIS
#ifdef BAROCLINIC
use param, only: ug_delta, vg_delta, bar_start, bar_end, geo_force_x, geo_force_y
use sim_param, only: geo_force_x_rot, geo_force_y_rot
#endif BAROCLINIC
#ifdef SCALAR
use scalars_param, only: beta_scal
use param, only: ubc, damping_x  
use sim_param, only: w 
use scalars_param, only: sponge, sponge_x
#endif SCALAR
implicit none
#ifdef CORIOLIS
integer :: dz1,dz2,dzfact
real(rp) :: ug_loc, vg_loc
#endif CORIOLIS
#ifdef SCALAR
integer :: j
#endif SCALAR
integer :: k

! Removing the mean stress in the wall normal direction. This component results
! in a mean acceleration in the vertical direction which is non-physical.
! Srinidhi : Feb 20, 2022 
! This is the modification that I introduced as it is mathematically right. 
! However this is not generally implemented in general LES codes.
! So we are doing something which is not normal. 
do k=0,nz
  tzz(1:ld,1:ny,k) = tzz(1:ld,1:ny,k) - sum(tzz(1:nx,1:ny,k))*inxny
enddo
 
! compute stress gradients
call ddx   (txz(1:ld,1:ny,1:nz-1), dtxdx)
call ddy   (tyz(1:ld,1:ny,1:nz-1), dtydy)
call ddz_uv(tzz(1:ld,1:ny,0:nz  ), dtzdz)

if (coord.eq.0) dtzdz(:,:,1) = 0.0_rp ! at wall dz(tzz)=0._rp
if (coord.eq.0) divtz(:,:) = dtxdx(:,:,1)+dtydy(:,:,1)

do k=1,nz-1
  RHSz(:,:,k) = -RHSz(:,:,k)-(dtxdx(:,:,k)+dtydy(:,:,k)+dtzdz(:,:,k))
end do

#ifdef CORIOLIS 

#ifndef BAROCLINIC
if(adjust_wind_angle) then
  ug_loc = vgMag * cosd(180.0 * alpha_wind/pi)
  vg_loc = vgMag * sind(180.0 * alpha_wind/pi)
else
  ug_loc = ug
  vg_loc = vg
endif
#endif BAROCLINIC

if(coord==0) then

  do k = 2,nz-1
#ifdef BAROCLINIC    
    if(adjust_wind_angle) then
      ug_loc = geo_force_x_rot(k)
      vg_loc = geo_force_y_rot(k)
    else
      ug_loc = geo_force_x(k)
      vg_loc = geo_force_y(k)
    endif
#endif BAROCLINIC
    
    dz1 = 0.5_rp * (gridzw(k) - gridzw(k-1))
    dz2 = 0.5_rp * (gridzw(k+1) - gridzw(k))
    dzfact = 1._rp / (dz1 + dz2)
    RHSz(:, :, k) = RHSz(:, :, k) + &
                  coriol_y * (dzfact * (dz2 * u(:,:,k-1) + dz1 * u(:,:,k)) - ug_loc) - &
                  coriol_x * (dzfact * (dz2 * v(:,:,k-1) + dz1 * v(:,:,k)) - vg_loc)
  enddo

else
  do k = 1,nz-1
#ifdef BAROCLINIC 
    if(adjust_wind_angle) then
      ug_loc = geo_force_x_rot(k)
      vg_loc = geo_force_y_rot(k)
    else
      ug_loc = geo_force_x(k)
      vg_loc = geo_force_y(k)
    endif
#endif BAROCLINIC

    dz1 = 0.5_rp * (gridzw(k) - gridzw(k-1))
    dz2 = 0.5_rp * (gridzw(k+1) - gridzw(k))
    dzfact = 1._rp / (dz1 + dz2)
    RHSz(:, :, k) = RHSz(:, :, k) + &
                  coriol_y * (dzfact * (dz2 * u(:,:,k-1) + dz1 * u(:,:,k)) - ug_loc) - &
                  coriol_x * (dzfact * (dz2 * v(:,:,k-1) + dz1 * v(:,:,k)) - vg_loc)
  enddo
endif

#endif CORIOLIS

#ifdef SCALAR
RHSz(:,:,1:nz-1) = RHSz(:,:,1:nz-1) + beta_scal(:,:,1:nz-1)

! If damping layer is activated
if(ubc == 1) then
  do k=1, nz-1
    RHSz(:,:,k) = RHSz(:,:,k) - sponge(k)*w(:,:,k)
  end do

  if(damping_x) then
  do k=1, nz-1
   do j=1,ny
    RHSz(:,j,k) = RHSz(:,j,k) - sponge_x(:) * (w(:,j,k)-sum(w(1:nx,1:ny,k))*inxny)
   enddo 
  end do
  endif 
endif
#endif SCALAR

end subroutine divstress_w
