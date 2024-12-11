#ifdef SCALAR 
!*******************************************************************************
subroutine ic_scal
!*******************************************************************************
use param, only: T_init, z_i, nz, sullivan_test, hbc
use sim_param, only: u_dim 
#ifdef CORIOLIS
use param, only: ug, vg, pi, adjust_wind_angle
use sim_param, only: vgMag, alpha_wind
#else
use param, only: zo, vonk, u_star
#endif CORIOLIS
#ifdef BAROCLINIC
use param, only: ug_delta, vg_delta, bar_start, bar_end, geo_force_x, geo_force_y
#endif BAROCLINIC
use param, only: coord, ny, nx, T_scale, gabls_test, neutral_test, nproc
use param, only: manual_test, cap_height, cap_strength, strat_height, strat_strength
use param, only: dz, inv_strength, rp
use grid_defs, only: gridz, gridzw
use mpi_defs, only: mpi_sync_downup
use sim_param,only: u,v,w
use scalars_param, only: theta, ubar, vbar, wbar
#ifdef LVLSET
use level_set_mod, only: zwall
#endif LVLSET
implicit none
real(rp) :: rms,noise,ran3,z_inv,z_inv2,DeltaTheta
real(rp) :: z_rel, zw_rel
#ifndef CORIOLIS
real(rp) :: arg2, arg
#endif CORIOLIS
integer :: jx,jy,jz,seed,jz_abs
real(rp) :: dz1 

!Geostrophic velocities (ug,vg) have to input in m/s 

do jz=1,nz

#ifdef CORIOLIS
  ubar(jz)=ug
  vbar(jz)=vg
  wbar(jz)=0._rp

#ifdef BAROCLINIC
  ! Linear geostrophic shear with height (corresponding to a height and time invariant horizontal temperature gradient)
  ! as previously used by Momen (2018, Journal of Atmospheric Sciences) and Sorbjan (2004, Bound.-Layer Meteor.)
  if( gridz(jz) <= bar_start) then
    geo_force_x(jz) = ug
    geo_force_y(jz) = vg
  elseif(bar_start < gridz(jz) .and. gridz(jz) <= bar_end) then
    geo_force_x(jz) = ug + ( ug_delta/(bar_end-bar_start) ) * ( gridz(jz)-bar_start )
    geo_force_y(jz) = vg + ( vg_delta/(bar_end-bar_start) ) * ( gridz(jz)-bar_start )
  else
    geo_force_x(jz) = ug + ug_delta
    geo_force_y(jz) = vg + vg_delta
  endif

  if(adjust_wind_angle) then
    ubar(jz)=geo_force_x(jz) * cosd(180.0 * alpha_wind/pi) - geo_force_y(jz) * sind(180.0 * alpha_wind/pi)
    vbar(jz)=geo_force_x(jz) * sind(180.0 * alpha_wind/pi) + geo_force_y(jz) * cosd(180.0 * alpha_wind/pi)
  else
    ubar(jz)=geo_force_x(jz)
    vbar(jz)=geo_force_y(jz)
  endif

#else
  if(adjust_wind_angle) then
    ubar(jz)=vgMag * cosd(180.0 * alpha_wind/pi)
    vbar(jz)=vgMag * sind(180.0 * alpha_wind/pi)
  endif
#endif BAROCLINIC  

#else 
  !Initialized with non-dimensional logarithmic velocity profile 
  !if coriolis_forcing is false 
  arg2=gridz(jz)/zo
  arg=(1._rp/vonk)*log(arg2)
  arg = arg*u_star

  ubar(jz)=arg
  vbar(jz)=0._rp
  wbar(jz)=0._rp
#endif CORIOLIS

end do

rms = 3.0_rp

! hbc=1: heat flux or cooling rate heat boundary condition
if(hbc .eq. 1) then

  if(manual_test) then
    !MANUAL defined in input.conf
    do jz=1,nz
      jz_abs = coord * (nz-1) + jz
      seed = -80 - jz_abs
      z_rel  = (gridz (jz))*z_i
      zw_rel = (gridzw(jz))*z_i
   
      do jy=1,ny
      do jx=1,nx
        if (z_rel .le. 200._rp) then
          noise = rms/0.289_rp*(ran3(seed)-0.5_rp)
          u(jx,jy,jz)     = 1.e-6 * noise*(1._rp-z_rel/cap_height) + ubar(jz)
          noise = rms/0.289_rp*(ran3(seed)-0.5_rp)
          v(jx,jy,jz)     = 1.e-6 * noise*(1._rp-z_rel/cap_height) + vbar(jz)
          noise = rms/0.289_rp*(ran3(seed)-0.5_rp)
          w(jx,jy,jz)     = 1.e-6* noise*(1._rp-zw_rel/cap_height) + wbar(jz)
          noise = rms/0.289_rp*(ran3(seed)-0.5_rp)
          theta(jx,jy,jz) = 1.e-6 * noise*(1._rp-z_rel/cap_height) + (T_init/T_scale)
        elseif ((z_rel .gt. 200._rp) .and. (z_rel .lt. cap_height)) then
          u(jx,jy,jz) =   ubar(jz)
          v(jx,jy,jz) =   vbar(jz)
          w(jx,jy,jz) =   wbar(jz)
          theta(jx,jy,jz) = (T_init/T_scale)
        elseif ((z_rel .ge. cap_height) .and. (z_rel .le. strat_height)) then
          u(jx,jy,jz) =   ubar(jz)
          v(jx,jy,jz) =   vbar(jz)
          w(jx,jy,jz) =   wbar(jz)
          theta(jx,jy,jz) = (T_init/T_scale) + (z_rel-cap_height)*cap_strength/T_scale
        else
          u(jx,jy,jz) =  ubar(jz)
          v(jx,jy,jz) =  vbar(jz)
          w(jx,jy,jz) =  wbar(jz)
          theta(jx,jy,jz) = (T_init/T_scale) + (strat_height-cap_height)*cap_strength/T_scale &
                                             + (z_rel-strat_height)*strat_strength/T_scale
        end if
      end do
      end do
    end do

  elseif(sullivan_test) then
    !Initialized as in Abkar and Moin(2017), Minimum dissipation model for LES,
    !Boundary-Layer Meteorology  
    ! CBL case
    z_inv = 0.9512_rp*z_i
    do jz=1,nz
      jz_abs = coord*(nz-1) + jz
      seed = -80 - jz_abs

#ifndef LVLSET
      z_rel  = (gridz (jz))*z_i
      zw_rel = (gridzw(jz))*z_i
#else
      z_rel  = (gridz (jz) - zwall)*z_i
      zw_rel = (gridzw(jz) - zwall)*z_i
#endif LVLSET

      do jy=1,ny
      do jx=1,nx
        if(z_rel .lt. z_inv) then
          noise = rms/0.289_rp*(ran3(seed)-0.5_rp)
          u(jx,jy,jz) = 0.01_rp*noise*(1._rp-z_rel/z_inv) + ubar(jz)
          noise = rms/0.289_rp*(ran3(seed)-0.5_rp)
          v(jx,jy,jz) = 0.01_rp*noise*(1._rp-z_rel/z_inv) + vbar(jz)
          noise = rms/0.289_rp*(ran3(seed)-0.5_rp)
          w(jx,jy,jz) = 0.01_rp*noise*(1._rp-zw_rel/z_inv) + wbar(jz)
          noise = rms/0.289_rp*(ran3(seed)-0.5_rp)
          theta(jx,jy,jz) = (T_init/T_scale) + 0.01_rp*noise*(1._rp-z_rel/z_inv)
        elseif((z_rel .gt. z_inv) .and. (z_rel .lt. 1.0488_rp*z_i)) then
          u(jx,jy,jz) = ubar(jz)
          v(jx,jy,jz) = vbar(jz)
          w(jx,jy,jz) = wbar(jz)
          theta(jx,jy,jz) = (T_init/T_scale) + (z_rel-z_inv)*0.08_rp/T_scale
        else
          u(jx,jy,jz) = ubar(jz)
          v(jx,jy,jz) = vbar(jz)
          w(jx,jy,jz) = wbar(jz)
          theta(jx,jy,jz) = (T_init/T_scale) + (8.0/T_scale) + (z_rel-1.0488_rp*z_i)*(-inv_strength)/T_scale
        end if
      end do
      end do
    end do

  elseif(gabls_test) then
    !Initialized as in Abkar and Moin(2017), Minimum dissipation model for LES,
    !Boundary-Layer Meteorology  
    z_inv = 1200._rp 
    z_inv2= z_inv + 200._rp
    DeltaTheta= 0.015_rp
    do jz=1,nz
      jz_abs = coord * (nz-1) + jz
      seed = -80 - jz_abs  !--trying to make consistent init for MPI

#ifndef LVLSET
      z_rel  = (gridz (jz))*z_i
      zw_rel = (gridzw(jz))*z_i
#else
      z_rel  = (gridz (jz) - zwall)*z_i
      zw_rel = (gridzw(jz) - zwall)*z_i
#endif LVLSET
   
      do jy=1,ny
      do jx=1,nx
        if (z_rel .le. 200._rp) then
          noise = rms/0.289_rp*(ran3(seed)-0.5_rp)
          u(jx,jy,jz) = 1.e-6 * noise*(1._rp-z_rel/z_inv) + ubar(jz)
          noise = rms/0.289_rp*(ran3(seed)-0.5_rp)
          v(jx,jy,jz) = 1.e-6 * noise*(1._rp-z_rel/z_inv) + vbar(jz) !noise
          noise = rms/0.289_rp*(ran3(seed)-0.5_rp)
          w(jx,jy,jz) = 1.e-6* noise*(1._rp-zw_rel/z_inv) + wbar(jz)
          noise = rms/0.289_rp*(ran3(seed)-0.5_rp)
          theta(jx,jy,jz) = (T_init/T_scale) + (1.e-6)*noise*(1._rp-z_rel/z_inv)
        elseif ((z_rel .gt. 200._rp) .and. (z_rel .lt. z_inv)) then
          noise = rms/0.289_rp*(ran3(seed)-0.5_rp)
          u(jx,jy,jz) =   ubar(jz)
          v(jx,jy,jz) =   vbar(jz)
          w(jx,jy,jz) =   wbar(jz)
          theta(jx,jy,jz) = (T_init/T_scale) !+ (z_rel-z_inv)*DeltaTheta/T_scale
        elseif ((z_rel .ge. z_inv) .and. (z_rel .le. z_inv2)) then
          u(jx,jy,jz) =   ubar(jz)
          v(jx,jy,jz) =   vbar(jz)
          w(jx,jy,jz) =   wbar(jz)
          theta(jx,jy,jz) = (T_init/T_scale) + (z_rel-z_inv)*DeltaTheta/T_scale
        else
          u(jx,jy,jz) =  ubar(jz)
          v(jx,jy,jz) =  vbar(jz)
          w(jx,jy,jz) =  wbar(jz)
          theta(jx,jy,jz) =(T_init/T_scale) + (z_inv2-z_inv)*DeltaTheta/T_scale + (z_rel-z_inv2)*(-inv_strength)/T_scale
        end if
      end do
      end do
    end do

  elseif(neutral_test) then 
    !Initialized as in Abkar and Moin(2017), Minimum dissipation model for LES,
    !Boundary-Layer Meteorology  
    z_inv = 1200._rp 
    z_inv2= z_inv + 200._rp
    DeltaTheta= 0.015_rp
    do jz=1,nz
      jz_abs = coord * (nz-1) + jz
      seed = -80 - jz_abs  !--trying to make consistent init for MPI

#ifndef LVLSET
      z_rel  = (gridz (jz))*z_i
      zw_rel = (gridzw(jz))*z_i
#else
      z_rel  = (gridz (jz) - zwall)*z_i
      zw_rel = (gridzw(jz) - zwall)*z_i
#endif LVLSET
   
      do jy=1,ny
      do jx=1,nx
        if (z_rel .le. 200._rp) then
          noise = rms/0.289_rp*(ran3(seed)-0.5_rp)
          u(jx,jy,jz) = 1.e-6 * noise*(1._rp-z_rel/z_inv) + ubar(jz)
          noise = rms/0.289_rp*(ran3(seed)-0.5_rp)
          v(jx,jy,jz) = 1.e-6 * noise*(1._rp-z_rel/z_inv) + vbar(jz) !noise
          noise = rms/0.289_rp*(ran3(seed)-0.5_rp)
          w(jx,jy,jz) = 1.e-6* noise*(1._rp-zw_rel/z_inv) + wbar(jz)
          noise = rms/0.289_rp*(ran3(seed)-0.5_rp)
          theta(jx,jy,jz) = (T_init/T_scale) + (1.e-6)*noise*(1._rp-z_rel/z_inv)
        elseif ((z_rel .gt. 200._rp) .and. (z_rel .lt. z_inv)) then
          noise = rms/0.289_rp*(ran3(seed)-0.5_rp)
          u(jx,jy,jz) =   ubar(jz)
          v(jx,jy,jz) =   vbar(jz)
          w(jx,jy,jz) =   wbar(jz)
          theta(jx,jy,jz) = (T_init/T_scale) !+ (z_rel-z_inv)*DeltaTheta/T_scale
        elseif ((z_rel .ge. z_inv) .and. (z_rel .le. z_inv2)) then
          u(jx,jy,jz) =   ubar(jz)
          v(jx,jy,jz) =   vbar(jz)
          w(jx,jy,jz) =   wbar(jz)
          theta(jx,jy,jz) = (T_init/T_scale) + (z_rel-z_inv)*DeltaTheta/T_scale
        else
          u(jx,jy,jz) =  ubar(jz)
          v(jx,jy,jz) =  vbar(jz)
          w(jx,jy,jz) =  wbar(jz)
          theta(jx,jy,jz) =(T_init/T_scale) + (z_inv2-z_inv)*DeltaTheta/T_scale + (z_rel-z_inv2)*(-inv_strength)/T_scale
        end if
      end do
      end do
    end do
  endif 

  !Boundary conditions 
  if(coord == 0) w(:,:,1) = 0._rp
  if(coord == nproc-1) then
    w(:, :, nz) = 0._rp
    u(:, :, nz) = u(:, :, nz-1)
    v(:, :, nz) = v(:, :, nz-1)
   
    dz1 = 0.5_rp * (gridzw(nz+1) - gridzw(nz-1)) 
    !theta(:, :, nz) = (T_init/T_scale)!theta(:, :, nz-1) !- (inv_strength*z_i/T_scale)*dz1
  endif
endif

! hbc=0: constant temperature heat (top and bottom) boundary condition
if(hbc .eq. 0) then

  z_inv = 0.5_rp*z_i
  do jz=1,nz
    jz_abs = coord * (nz-1) + jz
    seed = -80 - jz_abs

#ifndef LVLSET
    z_rel  = (gridz (jz))*z_i
    zw_rel = (gridzw(jz))*z_i
#else
    z_rel  = (gridz (jz) - zwall)*z_i
    zw_rel = (gridzw(jz) - zwall)*z_i
#endif LVLSET

    do jy=1,ny
    do jx=1,nx
      if(z_rel .lt. z_inv) then
        noise = rms/0.289_rp*(ran3(seed)-0.5_rp)
        u(jx,jy,jz) = 0.01_rp*noise*(1._rp-z_rel/z_inv) + ubar(jz)
        noise = rms/0.289_rp*(ran3(seed)-0.5_rp)
        v(jx,jy,jz) = 0.01_rp*noise*(1._rp-z_rel/z_inv) + vbar(jz) 
        noise = rms/0.289_rp*(ran3(seed)-0.5_rp)
        w(jx,jy,jz) = 0.01_rp*noise*(1._rp-zw_rel/z_inv) + wbar(jz)
        noise = rms/0.289_rp*(ran3(seed)-0.5_rp)
        theta(jx,jy,jz) = T_init/T_scale + 0.01_rp*noise*(1._rp-z_rel/z_inv)
      else
        u(jx,jy,jz) = ubar(jz)
        v(jx,jy,jz) = vbar(jz)
        w(jx,jy,jz) = wbar(jz)
        theta(jx,jy,jz) = T_init/T_scale
      end if
    end do
    end do
  end do

  ! Boundary conditions 
  if(coord == 0) w(:,:,1) = 0._rp
  if(coord == nproc-1) then
    w(:,:,nz) = 0._rp
    u(:,:,nz) = u(:,:,nz-1)
    v(:,:,nz) = v(:,:,nz-1)
    theta(:,:,nz) = T_init/T_scale
  endif
endif

! Syncronize data to make sure ghost layers are updated
call mpi_sync_real_array(u    ,mpi_sync_downup)
call mpi_sync_real_array(v    ,mpi_sync_downup)
call mpi_sync_real_array(w    ,mpi_sync_downup)
call mpi_sync_real_array(theta,mpi_sync_downup)

end subroutine ic_scal
#endif
