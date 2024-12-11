#ifdef SCALAR 
!*******************************************************************************
subroutine wallstress_scalar
!*******************************************************************************
use mpi_defs, only: mpi_sync_up
use param, only: lh,dt,gabls_test,sullivan_test,cooling_rate,rp,T_scale,z_i
use param, only: dz,nx,ny,nz,vonk,coord,nproc,inv_strength,wt_s,T_init
use param, only: ld, neutral_test, hbc, theta_s1, zo_scalar 
use param, only: manual_test, surface_cond, surface_rate, surface_flux
use filters, only: G_test
use fft, only: forw, back
use scalars_param, only: dTdz,theta
use scalars_param, only: T_s,phi_h,psi_h,surf_flux
use sim_param, only: ustar, u_dim 
use grid_defs, only: gridzw, gridz 
#ifdef CORIOLIS 
use param, only: ug 
#endif CORIOLIS 
implicit none
integer :: i,j
real(rp), dimension(ld,ny) :: scalar_node_1
real(rp), parameter :: eps = 100._rp * epsilon (0._rp)
#ifdef LVLSET
real(rp) :: dummy
#endif LVLSET

! hbc=1: heat flux or cooling rate heat boundary condition
if(hbc .eq. 1) then

  if (coord == 0) then
    scalar_node_1 = theta(:,:,1)
    call dfftw_execute_dft_r2c(forw,scalar_node_1,scalar_node_1)
    call mulr(scalar_node_1(1:ld,1:ny),G_test(1:lh,1:ny))
    call dfftw_execute_dft_c2r(back,scalar_node_1,scalar_node_1)

    do j=1,ny
    do i=1,nx
      
      if(manual_test) then
        if(surface_cond .eq. 0) then
          T_s(i,j) = T_init/T_scale
          surf_flux(i,j) = 0.0_rp
        elseif(surface_cond .eq. 1) then
          T_s(i,j) = T_s(i,j) + (surface_rate/(T_scale*u_dim*3600._rp))*(dt * z_i) 
          surf_flux(i,j) = (T_s(i,j) - scalar_node_1(i,j))*vonk*ustar(i,j) /     &
                           (dlog(gridz(1)/(zo_scalar)) - psi_h(i,j))
        elseif(surface_cond .eq. 2) then
          surf_flux(i,j) = surface_flux/T_scale/u_dim
          T_s(i,j) = scalar_node_1(i,j) + surf_flux(i,j)/(vonk*ustar(i,j)) *     &
                     (dlog(gridz(1)/(zo_scalar)) - psi_h(i,j))
        endif

      elseif(gabls_test) then 
        ! The dimension of dt is [second/meter]
#ifdef LVLSET
        T_s(i,j,1) = T_s(i,j,1) - (cooling_rate/(T_scale*u_dim*3600._rp))*(dt*z_i) 
        surf_flux(i,j) = (T_s(i,j,1) - scalar_node_1(i,j))*vonk*ustar(i,j) /   &
                         (dlog(gridz(1)/(zo_scalar)) - psi_h(i,j,1))
#else
        !Set T_s(i,j) = T_init/T_scale for a neutral case (CNBL/TNBL) 
        !For CNBL case : surf_flux = 0.0 
        !This has to be manually forced for the CNBL case otherwise there will be a slight nonzero flux. 
        T_s(i,j) = T_s(i,j) - (cooling_rate/(T_scale*u_dim*3600._rp))*(dt * z_i) 
        surf_flux(i,j) = (T_s(i,j) - scalar_node_1(i,j))*vonk*ustar(i,j) /     &
                         (dlog(gridz(1)/(zo_scalar)) - psi_h(i,j))
#endif LVLSET

      elseif(sullivan_test) then 
        surf_flux(i,j) = wt_s/T_scale
#ifdef LVLSET
        T_s(i,j,1) = scalar_node_1(i,j) + surf_flux(i,j)/(vonk*ustar(i,j)) *   &
                     (dlog(gridz(1)/(zo_scalar)) - psi_h(i,j,1))
        ! ustar can be zero inside the body; this prevents a zero variable in the denominator
        if(abs(ustar(i,j)) .lt. eps) T_s(i,j,1) = scalar_node_1(i,j)
#else
        T_s(i,j) = scalar_node_1(i,j) + surf_flux(i,j)/(vonk*ustar(i,j)) *     &
                   (dlog(gridz(1)/(zo_scalar)) - psi_h(i,j))
#endif LVLSET

      elseif(neutral_test) then 
        surf_flux(i,j) = 0._rp
#ifdef LVLSET
        T_s(i,j,1) = scalar_node_1(i,j)
#else
        T_s(i,j) = T_init/T_scale !scalar_node_1(i,j)
#endif LVLSET
      endif
    enddo
    enddo

    do j=1,ny
    do i=1,nx
#ifdef LVLSET
      dTdz(i,j,1) = -phi_h(i,j,1)*surf_flux(i,j)/(ustar(i,j)*vonk*dz/2._rp)
      ! ustar can be zero inside the body; this prevents a zero variable in the denominator
      if(abs(ustar(i,j)) .lt. eps) dTdz(i,j,1) = 0._rp
#else
      dTdz(i,j,1) = -phi_h(i,j)*surf_flux(i,j)/(ustar(i,j)*vonk*gridz(1))
#endif LVLSET
    end do
    end do
  endif

#ifdef LVLSET 
  ! The dimension of dt is [second/meter]
  if((coord.ne.0) .and. gabls_test) then
    dummy = cooling_rate / (T_scale*3600._rp) * (dt*z_i) 
    do j=1,ny
    do i=1,nx
      T_s(i,j,1) = T_s(i,j,1) - dummy 
    end do
    end do
  endif
#endif LVLSET
endif

! hbc=0: constant temperature heat (top and bottom) boundary condition
if(hbc .eq. 0) then

  if (coord == nproc-1) then
    dTdz(:,:,nz-1) = 0._rp
    dTdz(:,:,nz)   = 0._rp
  endif

  if (coord == 0) then
    scalar_node_1=theta(:,:,1)

    call dfftw_execute_dft_r2c(forw,scalar_node_1,scalar_node_1)
    call mulr(scalar_node_1(1:ld,1:ny),G_test(1:lh,1:ny))
    call dfftw_execute_dft_c2r(back,scalar_node_1,scalar_node_1)

    do j=1,ny
    do i=1,nx
#ifdef LVLSET
      T_s(i,j,1) = theta_s1/T_scale
      surf_flux(i,j) = (T_s(i,j,1) - scalar_node_1(i,j))*vonk*ustar(i,j) /     &
                       (dlog(gridz(1)/(zo_scalar)) - psi_h(i,j,1))
#else
      T_s(i,j) = theta_s1/T_scale
      surf_flux(i,j) = (T_s(i,j) - scalar_node_1(i,j))*vonk*ustar(i,j) /       &
                       (dlog(gridz(1)/(zo_scalar)) - psi_h(i,j))
#endif
    enddo
    enddo

    do j=1,ny
    do i=1,nx
#ifdef LVLSET
      dTdz(i,j,1) =-phi_h(i,j,1)*surf_flux(i,j)/(ustar(i,j)*vonk*dz/2._rp)
      ! ustar can be zero inside the body; this prevents a zero variable in the denominator
      if(abs(ustar(i,j)) .lt. eps) dTdz(i,j,1) = 0._rp
#else
      dTdz(i,j,1) =-phi_h(i,j)*surf_flux(i,j)/(ustar(i,j)*vonk*gridz(1))
#endif LVLSET
    end do
    end do
  endif
endif

end subroutine wallstress_scalar
#endif SCALAR
