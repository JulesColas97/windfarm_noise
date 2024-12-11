#ifdef LVLSET
!*******************************************************************************
! Description
!
! In the present immersed boundary method, level set function or signed distance  
! function is used to distinguish the fluid region and the solid region.
!
! phi, signed distance function (provided by user):
!  phi>0, fluid region; |phi|<=phi_band, band region; phi<0, solid region. 
!
! In fluid region, no special treatment is needed.
!
! In band region, the log-law wall-model is used to calculate the sgs stress 
! on the wall, which is used to reconstruct the sgs stress on the band points. 
!
! In solid region, all velocity components are enforced zero, as well as the
! Smagorinsky coeffecient and other quantities involved in sgs dynamic model. 
!
! Input files (for both uv-grid and w-grid): 
!  phi_u, normx_u, normy_u, normz_u
!  phi_w, normx_w, normy_w, normz_w
! norm is the unit normal vector of iso-surface of phi, normz >= 0. 
!
!*******************************************************************************
subroutine level_set_init
!*******************************************************************************
!
! This subroutine should be called for only once. 
!
use level_set_mod, only: phi_u,phi_w,normx_u,normy_u,normz_u,normx_w,normy_w,normz_w
use level_set_mod, only: phi_band,phi_cutoff,coeff_phi_band,coeff_phi_cutoff
use param, only: rp,dz,path,nx,ny,nz,ld,kmax
use mpi_defs, only: mpi_sync_downup
implicit none
character(64) :: folder

! Allocate all arrays to characterize the geometry
allocate (phi_u  (ld,ny,0:nz)); phi_u   = 1._rp 
allocate (normx_u(ld,ny,0:nz)); normx_u = 0._rp
allocate (normy_u(ld,ny,0:nz)); normy_u = 0._rp
allocate (normz_u(ld,ny,0:nz)); normz_u = 1._rp
  
allocate (phi_w  (ld,ny,0:nz)); phi_w   = 1._rp 
allocate (normx_w(ld,ny,0:nz)); normx_w = 0._rp
allocate (normy_w(ld,ny,0:nz)); normy_w = 0._rp
allocate (normz_w(ld,ny,0:nz)); normz_w = 1._rp

folder = path
call read_3Dfield (phi_u  (1:nx,1:ny,1:kmax), 'phi_u'  , folder, 1, kmax)
call read_3Dfield (normx_u(1:nx,1:ny,1:kmax), 'normx_u', folder, 1, kmax)
call read_3Dfield (normy_u(1:nx,1:ny,1:kmax), 'normy_u', folder, 1, kmax)
call read_3Dfield (normz_u(1:nx,1:ny,1:kmax), 'normz_u', folder, 1, kmax)
 
call read_3Dfield (phi_w  (1:nx,1:ny,1:kmax), 'phi_w'  , folder, 1, kmax)
call read_3Dfield (normx_w(1:nx,1:ny,1:kmax), 'normx_w', folder, 1, kmax)
call read_3Dfield (normy_w(1:nx,1:ny,1:kmax), 'normy_w', folder, 1, kmax)
call read_3Dfield (normz_w(1:nx,1:ny,1:kmax), 'normz_w', folder, 1, kmax)
 
call mpi_sync_real_array (phi_u  , mpi_sync_downup)
call mpi_sync_real_array (normx_u, mpi_sync_downup)
call mpi_sync_real_array (normy_u, mpi_sync_downup)
call mpi_sync_real_array (normz_u, mpi_sync_downup)

call mpi_sync_real_array (phi_w  , mpi_sync_downup)
call mpi_sync_real_array (normx_w, mpi_sync_downup)
call mpi_sync_real_array (normy_w, mpi_sync_downup)
call mpi_sync_real_array (normz_w, mpi_sync_downup)

! set the location where the velocity needs to be interpolated
phi_cutoff = dz * coeff_phi_cutoff
! detect closest grid nodes in band region with width phi_band
phi_band = dz * coeff_phi_band

end subroutine level_set_init


!*******************************************************************************
subroutine level_set_Cs_smag
!*******************************************************************************
!
! This subroutine should be called for only once. 
!
use param, only: rp,nx,ny,nz,dz,coord,Co, wall_damp_exp, vonk 
use sgs_param, only: Cs_opt2 
use filters, only: delta 
use level_set_mod,only: phi_u,phi_w,zwall
#ifdef SCALAR 
use param, only: Pr
use sgs_param, only: Ds_opt2
#endif 
implicit none
integer :: i, j, k, ktot
real(rp) :: phix
real(rp) :: dmin, ztot, lambda, lambda0

lambda0 = Co*delta

! The variable length calculated below is lambda 
! see ref. Bou-Zeid et al. (2005) eqn(2)
do k = 1, nz
do j = 1, ny
do i = 1, nx
  ktot = coord*(nz-1) + k
  if (ktot .eq. 1) then
    phix = phi_u(i,j,k)
    ztot = real (ktot - 0.5) * dz - zwall
  else
    phix = phi_w(i,j,k)
    ztot = real (ktot - 1) * dz - zwall
  end if

  if ((phix .le. 0._rp) .or. (ztot .le. 0._rp)) then
    Cs_opt2(i,j,k) = 0._rp
  else
    dmin = min (phix, ztot)
    lambda = lambda0 * (vonk*dmin) / (lambda0**wall_damp_exp +                 &
                       (vonk*dmin)**wall_damp_exp)**(1._rp/wall_damp_exp)
    Cs_opt2(i,j,k) = lambda**2
  end if
end do
end do
end do

#ifdef SCALAR 
Ds_opt2 = Cs_opt2/Pr 
if(coord.eq.0) print*, 'Pr',Pr 
#endif SCALAR

end subroutine level_set_Cs_smag


!*******************************************************************************
subroutine level_set_Cs_lag_dyn
!*******************************************************************************
!
! sets Cs2, F_LM, F_MM, F_QN, F_NN to zero inside solid
! These variables are stored on w-grid, Except jz=1 which are stored on u-grid
! Cs_opt2 are defined for k=1, 2, ..., nz
! F_LM, F_MM, F_QN, F_NN are defined for k=0, 1, ..., nz
!
use param, only: rp, nx, ny, nz, coord
use sgs_param, only: Cs_opt2, F_LM, F_MM, F_QN, F_NN
use level_set_mod,only: phi_u,phi_w
#ifdef SCALAR
use sgs_param, only: Ds_opt2, F_KX, F_XX, F_PY, F_YY 
#endif 
implicit none
integer :: i, j, k, ktot
real(rp) :: phix

do k = 1, nz
do j = 1, ny
do i = 1, nx
  ktot = coord*(nz-1) + k
  if (ktot .eq. 1) then
    phix = phi_u(i,j,k)
  else
    phix = phi_w(i,j,k)
  end if

  if (phix .le. 0._rp) then
    Cs_opt2(i,j,k) = 0._rp
    F_LM(i,j,k) = 0._rp
    F_MM(i,j,k) = 0._rp
    F_QN(i,j,k) = 0._rp
    F_NN(i,j,k) = 0._rp

#ifdef SCALAR
    Ds_opt2(i,j,k) = 0._rp
    F_KX(i,j,k) = 0._rp
    F_XX(i,j,k) = 0._rp
    F_PY(i,j,k) = 0._rp
    F_YY(i,j,k) = 0._rp
#endif 
  end if
end do
end do
end do

end subroutine level_set_Cs_lag_dyn


!*******************************************************************************
subroutine level_set_Cs_amd
!*******************************************************************************
!
! sets Cs2 to zero inside solid
! These variables are stored on w-grid, Except jz=1 which are stored on u-grid
! Cs_opt2 are defined for k=1, 2, ..., nz
!
use param, only: rp, nx, ny, nz, coord
use sgs_param, only: Cs_opt2
use level_set_mod,only: phi_u,phi_w
#ifdef SCALAR
use sgs_param, only: Ds_opt2
#endif 
implicit none
integer :: i, j, k, ktot
real(rp) :: phix

do k = 1, nz
do j = 1, ny
do i = 1, nx
  ktot = coord*(nz-1) + k
  if (ktot .eq. 1) then
    phix = phi_u(i,j,k)
  else
    phix = phi_w(i,j,k)
  end if

  if (phix .le. 0._rp) then
    Cs_opt2(i,j,k) = 0._rp
#ifdef SCALAR
    Ds_opt2(i,j,k) = 0._rp
#endif 
  end if
end do
end do
end do

end subroutine level_set_Cs_amd


!*******************************************************************************
subroutine level_set_sgs
!*******************************************************************************
!
! Calculate wall stress in the vicinity of the immersed boundary 
! at this point tij are only set for kmin:nz-1
!
! Ref: 
!     Chester, Meneveau and Parlange (2007) J. Comput. Phys. 225: 427-448
!
use param, only: rp,dz,kmin,nx,ny,nz,vonk
use sim_param, only: u, v, w
use sim_param, only: txx, txy, txz, tyy, tyz, tzz
use grid_defs,only: gridx,gridy
use level_set_mod, only: phi_u,phi_w,normx_u,normy_u,normz_u,normx_w,normy_w,normz_w
use level_set_mod, only: phi_band,phi_cutoff
use level_set_mod, only: zo_level_set
#ifdef SCALAR
use param, only: T_scale, wt_s, hbc, theta_s1
use param, only: sullivan_test, gabls_test, neutral_test
use sim_param, only: psi_m 
use sgs_param, only: txTh, tyTh, tzTh
use scalars_param, only: theta, T_s
use scalars_param, only: psi_h
#endif SCALAR
implicit none
real(rp), parameter :: eps = 100._rp * epsilon (0._rp)
real(rp), dimension(3) :: n_hat, vel, x_hat, x
real(rp) :: phix, ustar, tau, phi_c, mag_vel, invzo_level_set
#ifdef SCALAR
real(rp) :: theta_c, surf_flux
#endif SCALAR
integer :: i, j, k

invzo_level_set=1._rp/zo_level_set

! u-grid
do k = kmin, nz-1
do j = 1, ny
do i = 1, nx
  phix = phi_u(i,j,k)

  ! field band cells (in the fluid region only)
  if ((phix .ge. 0._rp) .and. (phix .le. phi_band)) then
    x = (/ gridx(i), gridy(j), (k-0.5_rp)*dz /)

    n_hat(1) = normx_u(i,j,k)
    n_hat(2) = normy_u(i,j,k)
    n_hat(3) = normz_u(i,j,k)

    ! determine velocity vector at point with phi ~ phi_c
    ! only uses information that is available to each processor to
    ! limit communication between processors 
    ! Therefore, x_hat must locate in the same process
    phi_c = max (phi_cutoff, phix)
    if (abs(n_hat(3)*(phi_c - phix)) .gt. 0.5_rp*dz) then
      phi_c = phix + 0.5_rp*dz
    endif
    x_hat = x + n_hat * (phi_c - phix)

    call interp_scal (u(:,:,0:nz), x_hat, vel(1), 0.5_rp) !0.5 is u node
    call interp_scal (v(:,:,0:nz), x_hat, vel(2), 0.5_rp) !0.5 is u node
    call interp_scal (w(:,:,0:nz), x_hat, vel(3), 1.0_rp) !1.0 is w node
#ifdef SCALAR
    call interp_scal (theta(:,:,0:nz), x_hat, theta_c, 0.5_rp) !0.5 is u node
#endif SCALAR

    vel = vel - dot_product (vel, n_hat) * n_hat  
    mag_vel=sqrt(dot_product(vel,vel))

    if (mag_vel < eps) then
      txx(i,j,k) = 0._rp
      txy(i,j,k) = 0._rp
      tyy(i,j,k) = 0._rp
      tzz(i,j,k) = 0._rp
#ifdef SCALAR
      txTh(i,j,k) = 0._rp
      tyTh(i,j,k) = 0._rp
#endif SCALAR
    else
      ! unit tangential vector
      x_hat = vel / mag_vel

#ifndef SCALAR
      ! assume a log-law layer between 0 and phi_c
      ustar = (vonk*mag_vel)/(dlog(phi_c*invzo_level_set))
#else
      ! predict friction velocity 
      ustar = (vonk*mag_vel)/(dlog(phi_c*invzo_level_set) - psi_m(i,j,k))

      ! update surface temperature and heat flux
      if(hbc.eq.0) then
        ! constant surface potential temperature
        T_s(i,j,1:k) = theta_s1/T_scale
        surf_flux = (vonk*ustar) * (T_s(i,j,1)-theta_c) /                      &
                    (dlog(phi_c*invzo_level_set)-psi_h(i,j,k))

      elseif(hbc.eq.1) then
        if(gabls_test) then
          ! T_s(i,j,1) is updated in wallstress_scalar.f90
          T_s(i,j,1:k) = T_s(i,j,1)
          surf_flux = (vonk*ustar) * (T_s(i,j,1)-theta_c) /                    &
                      (dlog(phi_c*invzo_level_set)-psi_h(i,j,k))

        elseif(sullivan_test) then
          ! constant surface heat flux
          surf_flux = wt_s/T_scale  
          T_s(i,j,1:k) = theta_c + surf_flux / (vonk*ustar) *                  &
                      (dlog(phi_c*invzo_level_set)-psi_h(i,j,k))

        elseif(neutral_test) then
          ! zero surface heat flux
          surf_flux = 0._rp
          T_s(i,j,1:k) = theta_c
        endif
      endif

      ! compute stability corrections psi_m, psi_h
      call level_set_PSI_compute (ustar, surf_flux, theta_c, phi_c, i, j, k)

      ! correct friction velocity 
      ustar = (vonk*mag_vel)/(dlog(phi_c*invzo_level_set) - psi_m(i,j,k))

      ! correct surface temperature and heat flux
      if(hbc.eq.0) then
        ! constant surface potential temperature
        surf_flux = (vonk*ustar) * (T_s(i,j,1)-theta_c) /                      &
                    (dlog(phi_c*invzo_level_set)-psi_h(i,j,k))

      elseif(hbc.eq.1) then
        if(gabls_test) then
          ! T_s(i,j,1) is updated in wallstress_scalar.f90
          T_s(i,j,1:k) = T_s(i,j,1)
          surf_flux = (vonk*ustar) * (T_s(i,j,1)-theta_c) /                    &
                      (dlog(phi_c*invzo_level_set)-psi_h(i,j,k))

        elseif(sullivan_test) then
          ! constant surface heat flux
          surf_flux = wt_s/T_scale  
          T_s(i,j,1:k) = theta_c + surf_flux / (vonk*ustar) *                  &
                      (dlog(phi_c*invzo_level_set)-psi_h(i,j,k))

        elseif(neutral_test) then
          ! zero surface heat flux
          surf_flux = 0._rp
          T_s(i,j,1:k) = theta_c
        endif
      endif
#endif SCALAR

      tau = -ustar * ustar
      ! sgs stress tensor in global coordinate system
      txx(i,j,k) = (x_hat(1)*n_hat(1) + x_hat(1)*n_hat(1)) * tau
      txy(i,j,k) = (x_hat(1)*n_hat(2) + x_hat(2)*n_hat(1)) * tau
      tyy(i,j,k) = (x_hat(2)*n_hat(2) + x_hat(2)*n_hat(2)) * tau
      tzz(i,j,k) = (x_hat(3)*n_hat(3) + x_hat(3)*n_hat(3)) * tau
#ifdef SCALAR
      ! heat flux vector in global coordinate system
      txTh(i,j,k) = -n_hat(1) * surf_flux
      tyTh(i,j,k) = -n_hat(2) * surf_flux
#endif
    end if
  end if
end do
end do
end do

! w-grid
do k = kmin, nz-1
do j = 1, ny
do i = 1, nx
  phix = phi_w(i,j,k)

  ! field band cells (in the fluid or solid region, closest to the wall)
  if (abs(phix) .le. 0.5_rp*phi_band) then
    x = (/ gridx(i), gridy(j), (k-1)*dz /)

    n_hat(1) = normx_w(i,j,k)
    n_hat(2) = normy_w(i,j,k)
    n_hat(3) = normz_w(i,j,k)

    ! determine velocity vector at point with phi ~ phi_c
    ! only uses information that is available to each processor to
    ! limit communication between processors 
    ! Therefore, x_hat must locate in the same process
    phi_c = max (phi_cutoff, phix)
    if (abs(n_hat(3)*(phi_c - phix)) .gt. dz) then
      phi_c = phix + dz
    endif
    x_hat = x + n_hat * (phi_c - phix)

    call interp_scal (u(:,:,0:nz), x_hat, vel(1), 0.5_rp) !0.5 is u node
    call interp_scal (v(:,:,0:nz), x_hat, vel(2), 0.5_rp) !0.5 is u node
    call interp_scal (w(:,:,0:nz), x_hat, vel(3), 1.0_rp) !1.0 is w node
#ifdef SCALAR
    call interp_scal (theta(:,:,0:nz), x_hat, theta_c, 0.5_rp) !0.5 is u node
#endif SCALAR

    vel = vel - dot_product (vel, n_hat) * n_hat
    mag_vel=sqrt(dot_product(vel,vel))

    if (mag_vel < eps) then
      txz(i,j,k) = 0._rp
      tyz(i,j,k) = 0._rp
#ifdef SCALAR
      tzTh(i,j,k) = 0._rp
#endif SCALAR
    else
      ! unit tangential vector
      x_hat = vel / mag_vel

#ifndef SCALAR
      ! assume a log-law layer between 0 and phi_c
      ustar = (vonk*mag_vel)/(dlog(phi_c*invzo_level_set))
#else
      ! predict friction velocity 
      ustar = (vonk*mag_vel)/(dlog(phi_c*invzo_level_set) - psi_m(i,j,k))

      ! update surface temperature and heat flux
      if(hbc.eq.0) then
        ! constant surface potential temperature
        T_s(i,j,1:k) = theta_s1/T_scale
        surf_flux = (vonk*ustar) * (T_s(i,j,1)-theta_c) /                      &
                    (dlog(phi_c*invzo_level_set)-psi_h(i,j,k))

      elseif(hbc.eq.1) then
        if(gabls_test) then
          ! T_s(i,j,1) is updated in wallstress_scalar.f90
          T_s(i,j,1:k) = T_s(i,j,1)
          surf_flux = (vonk*ustar) * (T_s(i,j,1)-theta_c) /                    &
                      (dlog(phi_c*invzo_level_set)-psi_h(i,j,k))

        elseif(sullivan_test) then
          ! constant surface heat flux
          surf_flux = wt_s/T_scale  
          T_s(i,j,1:k) = theta_c + surf_flux / (vonk*ustar) *                  &
                      (dlog(phi_c*invzo_level_set)-psi_h(i,j,k))

        elseif(neutral_test) then
          ! zero surface heat flux
          surf_flux = 0._rp
          T_s(i,j,1:k) = theta_c
        endif
      endif

      ! compute stability corrections psi_m, psi_h
      call level_set_PSI_compute (ustar, surf_flux, theta_c, phi_c, i, j, k)

      ! correct friction velocity 
      ustar = (vonk*mag_vel)/(dlog(phi_c*invzo_level_set) - psi_m(i,j,k))

      ! correct surface temperature and heat flux
      if(hbc.eq.0) then
        ! constant surface potential temperature
        T_s(i,j,1:k) = theta_s1/T_scale
        surf_flux = (vonk*ustar) * (T_s(i,j,1)-theta_c) /                      &
                    (dlog(phi_c*invzo_level_set)-psi_h(i,j,k))

      elseif(hbc.eq.1) then
        if(gabls_test) then
          ! T_s(i,j,1) is updated in wallstress_scalar.f90
          T_s(i,j,1:k) = T_s(i,j,1)
          surf_flux = (vonk*ustar) * (T_s(i,j,1)-theta_c) /                    &
                      (dlog(phi_c*invzo_level_set)-psi_h(i,j,k))

        elseif(sullivan_test) then
          ! constant surface heat flux
          surf_flux = wt_s/T_scale  
          T_s(i,j,1:k) = theta_c + surf_flux / (vonk*ustar) *                  &
                      (dlog(phi_c*invzo_level_set)-psi_h(i,j,k))

        elseif(neutral_test) then
          ! zero surface heat flux
          surf_flux = 0._rp
          T_s(i,j,1:k) = theta_c
        endif
      endif
#endif SCALAR

      tau = -ustar * ustar
      ! stress tensor in global coordinate
      txz(i,j,k) = (x_hat(1)*n_hat(3) + x_hat(3)*n_hat(1)) * tau
      tyz(i,j,k) = (x_hat(2)*n_hat(3) + x_hat(3)*n_hat(2)) * tau
#ifdef SCALAR
      ! heat flux vector in global coordinate system
      tzTh(i,j,k) = -n_hat(3) * surf_flux
#endif SCALAR
    end if
  end if
end do
end do
end do

end subroutine level_set_sgs


!*******************************************************************************
subroutine level_set_ic
!*******************************************************************************
!
! Set the initial condition inside the body: zero velocity and uniform
! potential temperature
!
use param, only: nx, ny, nz, rp
use sim_param, only: u, v, w
use level_set_mod, only: phi_u, phi_w
#ifdef SCALAR
use scalars_param, only: theta, T_s
#endif SCALAR
implicit none
integer :: i, j, k

! u-grid
do k = 1, nz
do j = 1, ny
do i = 1, nx
  if (phi_u(i,j,k) .le. 0._rp) then
    u(i,j,k) = 0._rp
    v(i,j,k) = 0._rp
#ifdef SCALAR
    ! set a uniform initial temperature inside the body
    ! which has the value as the surface temperature 
    theta(i,j,k) = T_s(i,j,k)
#endif SCALAR
  endif
enddo
enddo
enddo

! w-grid
do k = 1, nz
do j = 1, ny
do i = 1, nx
  if (phi_w(i,j,k) .le. 0._rp) then
    w(i,j,k) = 0._rp
  endif 
end do
end do
end do

end subroutine level_set_ic

 
!*******************************************************************************
subroutine level_set_forcing 
!*******************************************************************************
!
! Set the physical intermediate velocity inside solid zero 
! Actually here shows the numerical intermediate velocity
! dpdx, dpdy, dpdz are available for 1:nz-1, no mpi_sync
!
! u, v, w are defined for k=0, 1, ..., nz
! dpdx, dpdy, dpdz are defined for k=1, ..., nz-1
!
use param, only: tadv1, dt,nx,ny,nz,rp,kmin 
use sim_param, only: u, v, w, dpdx, dpdy, dpdz
use level_set_mod, only: phi_u,phi_w
#ifdef SCALAR
use scalars_param, only: theta, T_s
#endif SCALAR
implicit none
integer :: i, j, k
real(rp) :: dttadv1

dttadv1 = dt*tadv1

! u-grid
do k=1,nz-1
do j=1,ny
do i=1,nx
  if(phi_u(i,j,k).le.0._rp) then
    u(i,j,k) = dttadv1 * dpdx(i,j,k) 
    v(i,j,k) = dttadv1 * dpdy(i,j,k) 
#ifdef SCALAR
    theta(i,j,k) = T_s(i,j,k)
#endif SCALAR
  endif
enddo
enddo
enddo

! w-grid
do k=kmin,nz-1
do j=1,ny
do i=1,nx
  if(phi_w(i,j,k).le.0._rp) then
    w(i,j,k) = dttadv1 * dpdz(i,j,k) 
  endif 
end do
end do
end do

end subroutine level_set_forcing

 
!*******************************************************************************
subroutine interp_scal (a, x, a_x, s)
!*******************************************************************************
!
! performs tri-linear interpolation to obtain a at point x
!
use param, only: rp,ld,ny,dx,dy,nz,L_x,L_y,idx,idy,idz
use grid_defs, only: awrap_i, awrap_j
implicit none
real (rp), intent (in) :: a(ld,ny,0:nz)
real (rp), intent (in) :: x(3)
real (rp), intent (out) :: a_x
real (rp), intent (in) :: s

integer :: i, j, i1, j1, ks, ks1
real (rp) :: x1, x2, x3
real (rp) :: f1, f2, f3, f4, f5, f6, f7, f8
real (rp) :: w1, w2, w3, w4, w5, w6, w7, w8
real (rp) :: xmod(3)

xmod(1:3)=x(1:3)
xmod(1) = modulo(xmod(1), L_x) ! Ensures i is located in the domain
xmod(2) = modulo(xmod(2), L_y) ! Ensures j is located in the domain

! calculate indices
i = awrap_i( floor (xmod(1) * idx + 1._rp) )
j = awrap_j( floor (xmod(2) * idy + 1._rp) )

! try to handle boundaries nicely for the +1 indices
i1 = awrap_i( i + 1 )
j1 = awrap_j( j + 1 )

! This needs special treatment for u and w
ks = floor (xmod(3) *idz + s)
ks1 = ks + 1

! calculate interpolation weights
x1 = modulo (xmod(1), dx) * idx
x2 = modulo (xmod(2), dy) * idy
x3 = xmod(3) *idz - (ks - s)

w1 = (1._rp - x1) * (1._rp - x2) * (1._rp - x3)
w2 = (    x1    ) * (1._rp - x2) * (1._rp - x3)
w3 = (1._rp - x1) * (    x2    ) * (1._rp - x3)
w4 = (    x1    ) * (    x2    ) * (1._rp - x3)
w5 = (1._rp - x1) * (1._rp - x2) * (    x3    )
w6 = (    x1    ) * (1._rp - x2) * (    x3    )
w7 = (1._rp - x1) * (    x2    ) * (    x3    )
w8 = (    x1    ) * (    x2    ) * (    x3    )

f1 = a(i , j , ks)
f2 = a(i1, j , ks)
f3 = a(i , j1, ks)
f4 = a(i1, j1, ks)
f5 = a(i , j , ks1)
f6 = a(i1, j , ks1)
f7 = a(i , j1, ks1)
f8 = a(i1, j1, ks1)

a_x = w1*f1 + w2*f2 + w3*f3 + w4*f4 + w5*f5 + w6*f6 + w7*f7 + w8*f8

end subroutine interp_scal


#ifdef SCALAR
!*******************************************************************************
subroutine level_set_PSI_compute (ustar, wt_avg, theta_avg, phi_c, i, j, k)
!*******************************************************************************
!
! calculate stability corrections psi_m and psi_h
!
use param, only: rp, vonk, g, z_i
use sim_param, only: psi_m 
use scalars_param, only: psi_h
use level_set_mod, only: zo_level_set
implicit none
integer, intent(in) :: i, j, k
real(rp), intent(in) :: ustar, wt_avg, theta_avg, phi_c
real(rp), parameter :: eps = 100._rp * epsilon (0._rp)
real(rp), parameter :: inf = 10000._rp
real(rp) :: y,y0,grav
real(rp) :: L

grav = g * z_i 

! Compute Obukhov Length
if((abs(wt_avg).le.eps) .or. (abs(ustar).le.eps)) then
  L = inf
else    
  L = -ustar**3/(vonk*(grav/theta_avg)*wt_avg)
endif

! compute stability corrections psi_m, psi_h
if((L.lt.-eps)) then 
  ! Unstable BL
  ! Brutsaert, W. (1982). Evaporation into the Atmosphere: Theory, History
  ! and Applications. Netherlands: Springer. Section 4.2, pp. 68-70.
  y = (1._rp - 16._rp*(phi_c/L))**.25_rp
  y0= (1._rp - 16._rp*(zo_level_set/L))**.25_rp

  ! empirical formula
  psi_m(i,j,k) = dlog((1._rp+y )**2._rp*(1._rp+y **2._rp)) - 2._rp*datan(y ) - &
                (dlog((1._rp+y0)**2._rp*(1._rp+y0**2._rp)) - 2._rp*datan(y0))
  psi_h(i,j,k) = 2._rp*dlog((1._rp+y**2._rp)/(1._rp+y0**2._rp)) 

elseif((L.gt.eps) .and. (L.ne.inf)) then
  ! Stable BL
  ! Stoll, R., & Porté-Agel, F. (2008). Large-eddy simulation of the stable
  ! atmospheric boundary layer using dynamic models with different averaging
  ! schemes. Boundary-Layer Meteorology, 126, 1-28. 
  y = phi_c/L
  y0= zo_level_set/L

  ! empirical formula
  psi_m(i,j,k) = -4.8_rp*(y - y0)
  psi_h(i,j,k) = -7.8_rp*(y - y0)

else
  ! Neutral BL
  psi_m(i,j,k) = 0._rp
  psi_h(i,j,k) = 0._rp
end if

end subroutine level_set_PSI_compute
#endif
#endif
