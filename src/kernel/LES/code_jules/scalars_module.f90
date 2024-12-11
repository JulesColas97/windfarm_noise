#ifdef SCALAR
!*******************************************************************************
subroutine theta_all_in_one
!*******************************************************************************
use mpi_defs, only: mpi_sync_downup
use param, only: inflow,jt_total,rp, damping_x
use scalars_param, only: theta,RHS_T, RHS_Tf, RHS_Tff 
use param, only: tfac1, tfac2,cont_reduce 
use param, only: ubc,nx,ny,nz,inxny
use scalars_param, only: sponge, sponge_x 
implicit none
integer :: j,k,jz

RHS_Tff = RHS_Tf 
RHS_Tf=RHS_T

! Calculate the RHS of the scalar equation 
call scalar_RHS_calc  

! Add the damping term to the scalar equation
! sponge=0 except close to top boundary
if (ubc==1) then 
  do jz=1,nz-1                                      
    RHS_T(1:nx,1:ny,jz)=RHS_T(1:nx,1:ny,jz)-0.5_rp*(sponge(jz)+sponge(jz+1))*  &
                       (theta(1:nx,1:ny,jz)-sum(theta(1:nx,1:ny,jz))*inxny)
  end do

  !Damping layer in the horizontal direction
  if(damping_x) then 
  do k=1, nz-1
   do j=1,ny
    RHS_T(:,j,k) = RHS_T(:,j,k) - sponge_x(:) * (theta(:,j,k)-sum(theta(1:nx,1:ny,k))*inxny)
   enddo 
  end do
  endif 
end if

#ifdef CPS 
if(inflow) call inflow_cond_cps_scalar 
#endif CPS 

if(jt_total==1) then 
  RHS_Tf = RHS_T
  RHS_Tff= RHS_T 
endif 

if(jt_total==2 .OR. (cont_reduce)) then
  RHS_Tff= tfac1*RHS_T + tfac2*RHS_Tf
endif

! Calculates the buoyancy term which gets added to the vertical momentum equation
call calcbeta

! Step the scalar values 
call step_scalar

! Sync the new theta values 
call mpi_sync_real_array(theta,mpi_sync_downup)  

end subroutine theta_all_in_one


!*******************************************************************************
subroutine calcbeta
!*******************************************************************************
!
! This calculates the buoyancy term (beta_scal) to be added to the vertical
! momentum equation. Beta_scal is located on the w-node 
!
use param, only: kmin,nx,ny,nz,rp,g,z_i
use grid_defs, only: gridzw 
#ifndef LVLSET
use param, only: inxny
#endif LVLSET
use scalars_param, only: theta, beta_scal, scalar_bar  
use sim_param, only: u_dim 
implicit none
integer :: i, j, k
real(rp) :: above, below, dz1, dz2, dzfact 

! Calculate the mean scalar values
do k=kmin-1,nz-1
#ifndef LVLSET
  ! planar average
  scalar_bar(k) = sum(theta(1:nx,1:ny,k))*inxny
#else
  ! take as the reference temperature
  ! Ref. J. Fluid Mech. (2020), vol. 883, A39.
  !This will change the governing equations from (T-<T>) to (T-T0) when using
  !LVLSET: It should be noted. 
  scalar_bar(k) = 1._rp
#endif LVLSET
enddo

!Calculates the buoyancy term g(theta/theta0)
!theta is already non-dimensionalized so just g * theta  
do k=kmin,nz-1
do j=1,ny
do i=1,nx
  dz1 = 0.5_rp * (gridzw(k) - gridzw(k-1)) 
  dz2 = 0.5_rp * (gridzw(k+1) - gridzw(k)) 
  dzfact = (1._rp/(dz1+dz2))
  
  above=(theta(i,j,k)-scalar_bar(k))/scalar_bar(k)
  below=(theta(i,j,k-1)-scalar_bar(k-1))/scalar_bar(k-1)
  beta_scal(i,j,k)=g*dzfact*(above * dz1 + below * dz2)
end do
end do     
end do

end subroutine calcbeta


!*******************************************************************************
subroutine scalar_RHS_calc
!*******************************************************************************
use fft, only: forw,back,forw_big,back_big 
use param, only: kmin,coord,rp,ld, nx, ny, nz, nx2, ny2,time_heat
use scalars_param, only: RHS_T, dTdx,dTdy, dTdz
use sgs_param, only: txTh, tyTh, tzTh
use convec_mod, only: u_m=>u1_big,v_m=>u2_big,w_m=>u3_big
#ifndef MEMORY 
use convec_mod, only: dsdx_m=>vort1_big, dsdy_m=>vort2_big, RHS_m=>vort3_big
#else 
use scalars_param, only: dsdx_m, dsdy_m, RHS_m
#endif MEMORY
use sim_param, only : dplane=>dummy_plane
use scalars_param, only : dtemp, dsdz_m
#ifdef CORIOLIS
use param, only: Temp_sink,total_time,dt
use param, only: 
use param, only: Kpw, Kdw, Kiw, inxny
use scalars_param, only: error_int_sc, error_der_sc, error0_sc, error_sc, div_HF_var
use scalars_param, only: theta_aver_z, theta, theta_const
#endif CORIOLIS
#ifdef BAROCLINIC
use grid_defs, only: gridzw,gridz
use sim_param, only: u, v, u_dim
use param, only: coriol, bar_start, bar_end, ug_delta, vg_delta, geo_force_x, geo_force_y
use param, only: z_i, g, T_init, T_scale, ug, vg, temp_adv
#endif BAROCLINIC

implicit none
integer:: i, j, k, jz
#ifdef BAROCLINIC
real(rp) :: bar_fac, bar_grad_x, bar_grad_y
#endif BAROCLINIC

do jz=kmin,nz-1
  dplane(1:ld,1:ny) = dTdx(1:ld,1:ny,jz)  
  call dfftw_execute_dft_r2c(forw,dplane,dplane)
  call padd(dsdx_m(:,:,jz),dplane)
  call dfftw_execute_dft_c2r(back_big,dsdx_m(:,:,jz),dsdx_m(:,:,jz))

  dplane(1:ld,1:ny) = dTdy(1:ld,1:ny,jz) 
  call dfftw_execute_dft_r2c(forw,dplane,dplane)
  call padd(dsdy_m(:,:,jz),dplane)
  call dfftw_execute_dft_c2r(back_big,dsdy_m(:,:,jz),dsdy_m(:,:,jz))
enddo

do jz=kmin,nz
  dplane(1:ld,1:ny) = dTdz(1:ld,1:ny,jz)
  call dfftw_execute_dft_r2c(forw,dplane,dplane)
  call padd(dsdz_m(:,:,jz),dplane)
  call dfftw_execute_dft_c2r(back_big,dsdz_m(:,:,jz),dsdz_m(:,:,jz))
enddo

! Note that this is the advection term with the scalar as 
! the diffusion term has been thrown away. This is done step 
! by step for each of the expanded arrays from dealias1 separately
! for the first node & last node AND the rest of the nodes.

if (coord == 0)  then
  jz=1

  dplane(1:ld,1:ny) = dTdx(1:ld,1:ny,jz)
  call dfftw_execute_dft_r2c(forw,dplane,dplane)
  call padd(dsdx_m(:,:,jz),dplane)
  call dfftw_execute_dft_c2r(back_big,dsdx_m(:,:,jz),dsdx_m(:,:,jz))

  dplane(1:ld,1:ny) = dTdy(1:ld,1:ny,jz)
  call dfftw_execute_dft_r2c(forw,dplane,dplane)
  call padd(dsdy_m(:,:,jz),dplane)
  call dfftw_execute_dft_c2r(back_big,dsdy_m(:,:,jz),dsdy_m(:,:,jz))

  do j=1,ny2
  do i=1,nx2
    RHS_m(i,j,1) = u_m(i,j,1)*dsdx_m(i,j,1) + v_m(i,j,1)*dsdy_m(i,j,1) +       &
                   0.5_rp*w_m(i,j,2)*dsdz_m(i,j,2)
  end do
  end do
end if

do k=kmin,nz-1
do j=1,ny2
do i=1,nx2
  RHS_m(i,j,k) = u_m(i,j,k)*dsdx_m(i,j,k) + v_m(i,j,k)*dsdy_m(i,j,k) +         &
                 0.5_rp*(w_m(i,j,k)*dsdz_m(i,j,k)+w_m(i,j,k+1)*dsdz_m(i,j,k+1))
end do
end do
end do

do jz = 1,nz-1
  call dfftw_execute_dft_r2c(forw_big,RHS_m(:,:,jz),RHS_m(:,:,jz))
  call unpadd(RHS_T(:,:,jz),RHS_m(:,:,jz)) 
  call dfftw_execute_dft_c2r(back,RHS_T(:,:,jz),RHS_T(:,:,jz)) 
end do 

! Calculation of SGS fluxes 
call ddx(txTh,dtemp)
RHS_T(1:nx,1:ny,1:nz-1) =-RHS_T(1:nx,1:ny,1:nz-1) + dtemp(1:nx,1:ny,1:nz-1) 

call ddy(tyTh,dtemp)
RHS_T(1:nx,1:ny,1:nz-1) = RHS_T(1:nx,1:ny,1:nz-1) + dtemp(1:nx,1:ny,1:nz-1) 

call ddz_w(tzTh,dtemp) 
RHS_T(1:nx,1:ny,1:nz-1) = RHS_T(1:nx,1:ny,1:nz-1) + dtemp(1:nx,1:ny,1:nz-1) 

#ifdef BAROCLINIC
! Adding baroclinic temperature advection term
! Based on 'geostrophic advection'-approximation, as in Sorbjan (2004)
if (temp_adv) then
  bar_fac = (T_init/T_scale)*coriol / g

  do k = 1, nz-1

    ! Note this is taken constant over the entire domain, even though the shear layer
    ! may only cover part of the domain, to preserve imposed temperature inversion
    bar_grad_x = ( vg_delta/(bar_end-bar_start) ) * bar_fac
    bar_grad_y = - ( ug_delta/(bar_end-bar_start) ) * bar_fac

    RHS_T(1:nx,1:ny,k) = RHS_T(1:nx,1:ny,k) - geo_force_x(k)*bar_grad_x - geo_force_y(k)*bar_grad_y
  end do
end if  
#endif BAROCLINIC

#ifdef CORIOLIS
! PID controller for limiting the growth of the BL 
! Srinidhi : Feb 20, 2022  
! This is essentially for keeping a constant temperature profile 
! The idea was that keeping the temperature profile constant would fix the velocity profile 
! rendering the boundary layer completely stationary.
! However with Coriolis force the velocity profile oscillates about a mean profile and is never stationary. 
! Therefore, flows with Coriolis force would never reach an 'ideal steady state' 

if (Temp_sink) then
  do k = 1, nz-1

    theta_aver_z(k) = 0.0_rp 
    theta_aver_z(k) = sum(theta(1:nx,1:ny,k))*inxny

    !For a specified amount of time it steps the temperature
    !When current_time > time_heat this calculation stops &
    !the BL growth is restricted to theta_const
    !The temperature is calculated only at the top 
    ! if(total_time > (time_heat - 2._rp*dt) .and. & 
    !    total_time < (       time_heat       ) .and. & 
    !   (gridz(k)*z_i) > z_heat_min .and. (gridz(k)*z_i) < z_heat_max) then 

    if(total_time > (time_heat - 2._rp*dt) .and. total_time < time_heat) then 
      theta_const(k) = theta_aver_z(k)
    endif 

    ! if((gridz(k)*z_i) > z_heat_min .and. (gridz(k)*z_i) < z_heat_max .and. total_time > time_heat) then 
    if(total_time > time_heat) then 

      ! proportional error
      error_sc = theta_const(k) - theta_aver_z(k) 

      ! integral error  
      error_int_sc = error_int_sc + error_sc * dt 
      error_int_sc = error_int_sc/total_time
 
      ! derivative error
      error_der_sc = (error_sc - error0_sc)/dt  

      div_HF_var(k) = Kpw * error_sc + Kdw*error_der_sc + Kiw*error_int_sc
    endif 

    error0_sc = error_sc 
  enddo

  ! Add the PID output to the RHS to limit the growth of the boundary layer    
  ! if ((gridz(k)*z_i)  > z_heat_min .and. (gridz(k)*z_i) < z_heat_max .and. total_time > time_heat) then
  if(total_time > time_heat) then 
    do k=1,nz-1
      RHS_T(1:nx,1:ny,k) = RHS_T(1:nx,1:ny,k) + div_HF_var(k)
    enddo 
  end if
end if
#endif CORIOLIS

end subroutine scalar_RHS_calc


!*******************************************************************************
subroutine step_scalar
!*******************************************************************************
use param, only: T_scale, dt, rp,dz, nx, ny, nz, coord, nproc
use param, only: tadv1, tadv2,tadv3, inv_strength, z_i, hbc, T_init
use scalars_param, only: theta, RHS_T, RHS_Tf,RHS_Tff
use grid_defs, only: gridzw, gridz 
 
implicit none
real(rp) :: dz1 

theta(1:nx,1:ny,1:nz-1) = theta(1:nx,1:ny,1:nz-1) + dt*                        &
                         (tadv1*RHS_T (1:nx,1:ny,1:nz-1) +                     &
                          tadv2*RHS_Tf(1:nx,1:ny,1:nz-1) +                     &
                          tadv3*RHS_Tff(1:nx,1:ny,1:nz-1))

if(hbc .eq. 1) then
  if(coord == nproc-1) then
      dz1 = (gridz(nz-1) - gridz(nz-2)) 
      ! When inversion layer is used 
!      theta(1:nx,1:ny,nz-1) = theta(1:nx,1:ny,nz-2) - inv_strength/T_scale*z_i*dz1
    
      dz1 = (gridz(nz) - gridz(nz-1)) 
!      theta(1:nx,1:ny,nz)   = theta(1:nx,1:ny,nz-1) - inv_strength/T_scale*z_i*dz1
  endif
endif

if(hbc .eq. 0) theta(1:nx,1:ny,nz) = T_init/T_scale

end subroutine step_scalar


!*******************************************************************************
subroutine obukhov
!*******************************************************************************
use param, only: rp, ld, nx, ny, jt_total, dz, vonk
use param, only: wt_s, T_scale, zo_scalar 
use sim_param,only: u,v,psi_m,ustar
use sgs_param, only: tzTh
use grid_defs, only: gridzw, gridz 
use scalars_param, only: theta_avg
use scalars_param, only: theta, wt_avg
#ifndef LVLSET
use param, only: lh
use filters, only: G_test
use fft, only: forw, back 
#endif LVLSET
implicit none
real(rp), dimension(ld,ny) :: u1,v1 
real(rp), dimension(nx,ny) :: u_avg
real(rp) :: dummy

  !Thermally stratified calculations  
  if (jt_total .eq. 1) then 
    print *,'Thermally Stratified run !'
    tzTh(:,:,1) = wt_s/T_scale 
  endif   

#ifndef LVLSET
  u1 = u(:,:,1) 
  call dfftw_execute_dft_r2c(forw,u1,u1)
  call mulr(u1(1:ld,1:ny),G_test(1:lh,1:ny))
  call dfftw_execute_dft_c2r(back,u1,u1)

  v1 = v(:,:,1) 
  call dfftw_execute_dft_r2c(forw,v1,v1)
  call mulr(v1(1:ld,1:ny),G_test(1:lh,1:ny))
  call dfftw_execute_dft_c2r(back,v1,v1)

  ! Srinidhi : Feb 20, 2022
  ! Planar averaged values are being used here for calculating 
  ! surface heat flux and u_star values 
  ! The problem is for some cases the code was being unstable 
  ! when local variation was being considered. 
  ! Planar average ALWAYS works. 
  theta_avg = sum(theta(1:nx,1:ny,1))/(nx*ny) 
!  call dfftw_execute_dft_r2c(forw,theta_avg,theta_avg)
!  call mulr(theta_avg(1:ld,1:ny),G_test(1:lh,1:ny))
!  call dfftw_execute_dft_c2r(back,theta_avg,theta_avg)

  wt_avg(:,:) = -sum(tzTh(1:nx,1:ny,1))/(nx*ny) 
!  call dfftw_execute_dft_r2c(forw,wt_avg,wt_avg)
!  call mulr(wt_avg(1:ld,1:ny),G_test(1:lh,1:ny))
!  call dfftw_execute_dft_c2r(back,wt_avg,wt_avg)

!  theta_avg = theta(:,:,1)
!  call dfftw_execute_dft_r2c(forw,theta_avg,theta_avg)
!  call mulr(theta_avg(1:ld,1:ny),G_test(1:lh,1:ny))
!  call dfftw_execute_dft_c2r(back,theta_avg,theta_avg)
!
!  wt_avg(:,:) = -tzTh(:,:,1)
!  call dfftw_execute_dft_r2c(forw,wt_avg,wt_avg)
!  call mulr(wt_avg(1:ld,1:ny),G_test(1:lh,1:ny))
!  call dfftw_execute_dft_c2r(back,wt_avg,wt_avg)

#else
! If LVLSET use instantaneous rather than 2delta filtered values
u1 = u(:,:,1) 
v1 = v(:,:,1) 
theta_avg = theta(:,:,1)
wt_avg(:,:) = -tzTh(:,:,1)
#endif LVLSET

  ! predict friction velocity 
  u_avg(1:nx,1:ny) = sqrt(u1(1:nx,1:ny)**2 + v1(1:nx,1:ny)**2)
  dummy            = dlog(gridz(1)/zo_scalar)
#ifdef LVLSET
  ustar(1:nx,1:ny) = (vonK * u_avg(1:nx,1:ny)) / (dummy - psi_m(1:nx,1:ny,1))
#else
  ustar(1:nx,1:ny) = (vonK * u_avg(1:nx,1:ny)) / (dummy - psi_m(1:nx,1:ny))
#endif LVLSET

  ! compute Obukhov length L, stability corrections psi_m, psi_h, and
  ! non-dimensional shear and temperature gradient phi_m, phi_h
  call PSI_compute

  ! correct friction velocity 
#ifdef LVLSET
  ustar(1:nx,1:ny) = (vonK * u_avg(1:nx,1:ny)) / (dummy - psi_m(1:nx,1:ny,1))
#else
  ustar(1:nx,1:ny) = (vonK * u_avg(1:nx,1:ny)) / (dummy - psi_m(1:nx,1:ny))
#endif LVLSET

end subroutine obukhov 


!*******************************************************************************
subroutine PSI_compute
!*******************************************************************************
use param, only: rp, nx, ny
use sim_param, only: phi_m, psi_m
use scalars_param, only: L
use scalars_param, only: phi_h, psi_h
use param, only: vonk, g, z_i, dz, zo_scalar
use sim_param, only: ustar
use scalars_param, only: theta_avg, wt_avg
use grid_defs, only: gridz 
use sim_param, only: u_dim 
implicit none
real(rp), parameter :: eps = 100._rp * epsilon (0._rp)
real(rp), parameter :: inf = 10000._rp
real(rp):: y,y0,grav
integer:: i,j

grav = g * z_i/u_dim**2  

! Compute Obukhov Length
! Srinidhi : Feb 20, 2022 
! While this works for most cases. 'Sometimes' this can cause the code to crash (no idea why!). 
! In those cases when code crashes it is wiser to calculate an average Obukhov length instead of 
! locally varying Obukhov length! 
do j=1,ny
do i=1,nx
  if((abs(wt_avg(i,j)).le.eps) .or. (abs(ustar(i,j)).le.eps)) then
    L(i,j) = inf
  else    
    L(i,j) = -ustar(i,j)**3/(vonk*(grav/theta_avg(i,j))*wt_avg(i,j))
  endif
enddo 
enddo 

! Obukhov length < 0 : Unstable BL 
! Obukhov length > 0 : Stable BL 

do j=1,ny
do i=1,nx

  !For unstable BL calculations 
  !if ((L(i,j)<0._rp) .and. (wt_avg(i,j) .ne. 0._rp)) then 
  !if(sullivan_test) then 
  if ((L(i,j) .lt. -eps)) then
    !Adapted from Brutsaert text book: Evaporation into the atmosphere, Section 4.2
    !The constant 16 in the calculation of y is not a universally accepted constant
    !This value also depends on the Von Karman constant 
    !The present equations are taken from Brutsaert (1982) for a Karman constant of 0.4 
    !If you wish to change this, use the Karman constant associated with respective empirical
    !formulations 

    ! Brutsaert, W. (1982). Evaporation into the Atmosphere: Theory, History
    ! and Applications. Netherlands: Springer. Section 4.2, pp. 68-70.
    y = (1._rp - 16._rp*((gridz(1))/L(i,j)) )**.25_rp
    y0= (1._rp - 16._rp*((zo_scalar)/L(i,j)) )**.25_rp

#ifdef LVLSET
    ! empirical formula
    phi_m(i,j,1) = y**(-1._rp)
    phi_h(i,j,1) = y**(-2._rp)

    ! derived from the empirical formula
    psi_m(i,j,1) = dlog((1._rp+y )**2._rp*(1._rp+y **2._rp)) - 2._rp*datan(y) -&
                  (dlog((1._rp+y0)**2._rp*(1._rp+y0**2._rp)) - 2._rp*datan(y0))
    psi_h(i,j,1) = 2._rp*dlog((1._rp+y**2._rp)/(1._rp+y0**2._rp)) 
#else
    ! empirical formula
    phi_m(i,j) = y**(-1._rp)
    phi_h(i,j) = y**(-2._rp)

    ! derived from the empirical formula
    psi_m(i,j) = dlog((1._rp+y )**2._rp*(1._rp+y **2._rp)) - 2._rp*datan(y ) - &
                (dlog((1._rp+y0)**2._rp*(1._rp+y0**2._rp)) - 2._rp*datan(y0))
    psi_h(i,j) = 2._rp*dlog((1._rp+y**2._rp)/(1._rp+y0**2._rp)) 
#endif LVLSET

  !For stable BL calculations 
  !elseif ((L(i,j)>0._rp).and.(wt_avg(i,j).ne. 0._rp)) then
  !elseif(gabls_test) then      
  elseif ((L(i,j) .gt. eps) .and. (L(i,j) .ne. inf)) then
    ! The formulations are from GABLS1 study 
    ! For more information see Beare et al. (2006) or GABLS1 documentation

    ! Stoll, R., & Porté-Agel, F. (2008). Large-eddy simulation of the stable
    ! atmospheric boundary layer using dynamic models with different averaging
    ! schemes. Boundary-Layer Meteorology, 126, 1-28. 
    y = (gridz(1))/L(i,j)
    y0= zo_scalar/L(i,j)

#ifdef LVLSET
    ! empirical formula
    phi_m(i,j,1) = 1._rp + 4.8_rp * y
    phi_h(i,j,1) = 1._rp + 7.8_rp * y

    ! derived from the empirical formula
    psi_m(i,j,1) = -4.8_rp * (y - y0)
    psi_h(i,j,1) = -7.8_rp * (y - y0)
#else
    ! empirical formula
    phi_m(i,j) = 1._rp + 4.8_rp * y
    phi_h(i,j) = 1._rp + 7.8_rp * y

    ! derived from the empirical formula
    psi_m(i,j) = -4.8_rp * (y - y0)
    psi_h(i,j) = -7.8_rp * (y - y0)
#endif LVLSET

  !For neutral BL calculations 
  !elseif(neutral_test) then 
  else
#ifdef LVLSET
    psi_m(i,j,1) = 0.0_rp
    psi_h(i,j,1) = 0.0_rp
    phi_m(i,j,1) = 1.0_rp
    phi_h(i,j,1) = 1.0_rp
#else
    psi_m(i,j) = 0.0_rp
    psi_h(i,j) = 0.0_rp
    phi_m(i,j) = 1.0_rp
    phi_h(i,j) = 1.0_rp
#endif LVLSET
  endif

enddo 
enddo 

end subroutine PSI_compute
#endif SCALAR
