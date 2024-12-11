module sim_param
use param, only : rp
implicit none

real(rp) :: u_dim = 1.0 
real(rp), dimension(:,:,:), allocatable :: u, v, w
real(rp), dimension(:,:,:), allocatable :: pres 
real(rp), dimension(:,:,:), allocatable :: dudx, dudy, dudz
real(rp), dimension(:,:,:), allocatable :: dvdx, dvdy, dvdz
real(rp), dimension(:,:,:), allocatable :: dwdx, dwdy, dwdz
real(rp), dimension(:,:,:), allocatable :: RHSx, RHSy, RHSz
real(rp), dimension(:,:,:), allocatable :: RHSx_f, RHSy_f, RHSz_f
real(rp), dimension(:,:,:), allocatable :: RHSx_ff, RHSy_ff, RHSz_ff
real(rp), dimension(:,:,:), allocatable :: txx, txy, tyy
real(rp), dimension(:,:,:), allocatable :: txz, tyz, tzz
#ifdef LVLSET
real(rp), dimension(:,:,:), allocatable :: dpdx, dpdy, dpdz
#endif LVLSET
#ifdef ATM
real(rp), dimension(:,:,:), allocatable :: fxa, fya, fza
#endif ATM
#ifdef TURBINES 
real(rp), dimension(:,:,:), allocatable :: fx, fy
#endif TURBINES

! Dummy variables used throughout the code
real(rp), dimension(:,:,:), allocatable :: dummy1,dummy2,dummy3,dummy4

! Used as buffer array in sgs_stag/lagrange_Sdep
real(rp), dimension(:,:  ), allocatable :: dummy_plane

! Used for the tridiagonal solver in the pressure solver
real(rp), dimension(:,:,:), allocatable :: b
real(rp), dimension(:    ), allocatable :: a,c

real(rp), dimension(:,:  ), allocatable :: divtz,ustar
#ifdef SCALAR
#ifdef LVLSET
real(rp), dimension(:,:,:), allocatable :: phi_m,psi_m
#else
real(rp), dimension(:,:  ), allocatable :: phi_m,psi_m
#endif LVLSET
#endif SCALAR

#ifdef CORIOLIS 
real(rp) :: error_int, error_der, error_pro, error_pre
real(rp) :: alpha_wind, alpha_wind_old, omega, omega_eff, phi_new, phi_0, time1, time2,vgMag 
#endif CORIOLIS

#ifdef BAROCLINIC
real(rp), dimension(:), allocatable :: geo_force_x_rot, geo_force_y_rot
#endif BAROCLINIC

contains

!*******************************************************************************
subroutine sim_param_init
!*******************************************************************************
!
! Initialize all variables to zero
!
use param, only: coord,nx,ld,ny,nz,lh,sgs_model,nz_tot
use filters, only: G_test,G_test_test 
use fft, only: kx,ky
use grid_defs, only: gridx,gridy,gridz,gridzw,awrap_i,awrap_j,z_tot
#ifndef MEMORY
use param, only: ld_big,ny2,kmin
use convec_mod, only: u1_big,u2_big,u3_big,vort1_big,vort2_big,vort3_big,cc_big
#else
#ifdef SCALAR
use param, only: ld_big,ny2,kmin
use convec_mod, only: u1_big, u2_big, u3_big 
#endif SCALAR
#endif MEMORY
implicit none

u_dim = 1.0 
allocate(u     (ld,ny,0:nz  )); u     =0.0_rp
allocate(v     (ld,ny,0:nz  )); v     =0.0_rp
allocate(w     (ld,ny,0:nz  )); w     =0.0_rp
allocate(pres  (ld,ny,0:nz  )); pres  =0.0_rp
allocate(dudx  (ld,ny,0:nz  )); dudx  =0.0_rp
allocate(dudy  (ld,ny,0:nz  )); dudy  =0.0_rp
allocate(dudz  (ld,ny,0:nz  )); dudz  =0.0_rp
allocate(dvdx  (ld,ny,0:nz  )); dvdx  =0.0_rp
allocate(dvdy  (ld,ny,0:nz  )); dvdy  =0.0_rp
allocate(dvdz  (ld,ny,0:nz  )); dvdz  =0.0_rp
allocate(dwdx  (ld,ny,0:nz  )); dwdx  =0.0_rp
allocate(dwdy  (ld,ny,0:nz  )); dwdy  =0.0_rp
allocate(dwdz  (ld,ny,0:nz  )); dwdz  =0.0_rp
allocate(RHSx  (ld,ny,1:nz-1)); RHSx  =0.0_rp
allocate(RHSy  (ld,ny,1:nz-1)); RHSy  =0.0_rp
allocate(RHSz  (ld,ny,1:nz-1)); RHSz  =0.0_rp
allocate(RHSx_f(ld,ny,1:nz-1)); RHSx_f=0.0_rp
allocate(RHSy_f(ld,ny,1:nz-1)); RHSy_f=0.0_rp
allocate(RHSz_f(ld,ny,1:nz-1)); RHSz_f=0.0_rp
allocate(RHSx_ff(ld,ny,1:nz-1)); RHSx_ff=0.0_rp
allocate(RHSy_ff(ld,ny,1:nz-1)); RHSy_ff=0.0_rp
allocate(RHSz_ff(ld,ny,1:nz-1)); RHSz_ff=0.0_rp
allocate(txx   (ld,ny,1:nz  )); txx   =0.0_rp
allocate(txy   (ld,ny,1:nz  )); txy   =0.0_rp
allocate(tyy   (ld,ny,1:nz  )); tyy   =0.0_rp
allocate(txz   (ld,ny,0:nz  )); txz   =0.0_rp
allocate(tyz   (ld,ny,0:nz  )); tyz   =0.0_rp
allocate(tzz   (ld,ny,0:nz  )); tzz   =0.0_rp
#ifdef LVLSET
allocate(dpdx  (ld,ny,1:nz-1)); dpdx  =0.0_rp
allocate(dpdy  (ld,ny,1:nz-1)); dpdy  =0.0_rp
allocate(dpdz  (ld,ny,1:nz-1)); dpdz  =0.0_rp
#endif LVLSET
#ifdef TURBINES 
allocate(fx    (ld,ny,0:nz-1)); fx   =0.0_rp
allocate(fy    (ld,ny,0:nz-1)); fy   =0.0_rp
#endif TURBINES
#ifdef ATM
allocate(fxa   (ld,ny,0:nz)); fxa   =0.0_rp
allocate(fya   (ld,ny,0:nz)); fya   =0.0_rp
allocate(fza   (ld,ny,0:nz)); fza   =0.0_rp
#endif ATM
#ifdef BAROCLINIC
allocate(geo_force_x_rot  (1:nz)); geo_force_x_rot =0.0_rp
allocate(geo_force_y_rot  (1:nz)); geo_force_y_rot =0.0_rp
#endif BAROCLINIC

allocate(ustar (ld,ny)); ustar =0.0_rp

if(coord.eq.0) then 
  allocate(divtz (ld,ny)); divtz =0.0_rp
endif

#ifdef SCALAR
#ifdef LVLSET
allocate(phi_m (ld,ny,0:nz)); phi_m =1.0_rp
allocate(psi_m (ld,ny,0:nz)); psi_m =0.0_rp
#else
if(coord.eq.0) then 
  allocate(phi_m (ld,ny)); phi_m =1.0_rp
  allocate(psi_m (ld,ny)); psi_m =0.0_rp
endif
#endif LVLSET
#endif SCALAR

allocate(dummy_plane(ld,ny  ));dummy_plane=0.0_rp

! Dummy variables used throughout the code
allocate(dummy1(ld,ny,0:nz  )); dummy1=0.0_rp
allocate(dummy2(ld,ny,0:nz  )); dummy2=0.0_rp
allocate(dummy3(ld,ny,0:nz  )); dummy3=0.0_rp
allocate(dummy4(ld,ny,0:nz  )); dummy4=0.0_rp

allocate(b(lh,ny,nz+1)); b=0.0_rp
allocate(a(nz+1));       a=0.0_rp
allocate(c(nz+1));       c=0.0_rp

! From filters 
allocate(G_test(lh,ny) ); G_test=0.0_rp
if (sgs_model==5) then  !--scale dependent dynamic
  allocate(G_test_test(lh,ny) ); G_test_test=0.0_rp
endif 

! From FFT 
allocate(kx(lh,ny) ); kx=0.0_rp
allocate(ky(lh,ny) ); ky=0.0_rp

! Griddefs
allocate(gridx  (nx+1)  ); gridx  =0.0_rp
allocate(gridy  (ny+1)  ); gridy  =0.0_rp
allocate(gridz  (0:nz+1)  ); gridz  =0.0_rp
allocate(gridzw (-1:nz+1)  ); gridzw =0.0_rp
allocate(awrap_i(0:nx+1)); awrap_i=0
allocate(awrap_j(0:ny+1)); awrap_j=0
allocate(z_tot(nz_tot));   z_tot =0.0_rp

! Convec
#ifndef MEMORY
allocate(u1_big   (ld_big,ny2,kmin-1:nz-1)); u1_big   =0.0_rp
allocate(u2_big   (ld_big,ny2,kmin-1:nz-1)); u2_big   =0.0_rp
allocate(u3_big   (ld_big,ny2,kmin  :nz  )); u3_big   =0.0_rp
allocate(vort1_big(ld_big,ny2,kmin  :nz  )); vort1_big=0.0_rp
allocate(vort2_big(ld_big,ny2,kmin  :nz  )); vort2_big=0.0_rp
allocate(vort3_big(ld_big,ny2,     1:nz-1)); vort3_big=0.0_rp
allocate(cc_big   (ld_big,ny2            )); cc_big   =0.0_rp
#else
#ifdef SCALAR
allocate(u1_big   (ld_big,ny2,kmin-1:nz-1)); u1_big   =0.0_rp
allocate(u2_big   (ld_big,ny2,kmin-1:nz-1)); u2_big   =0.0_rp
allocate(u3_big   (ld_big,ny2,kmin  :nz  )); u3_big   =0.0_rp
#endif SCALAR
#endif MEMORY

end subroutine sim_param_init
end module sim_param
