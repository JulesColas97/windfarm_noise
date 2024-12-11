#ifdef SCALAR
module scalars_param
use param
implicit none

save
public

real(rp),dimension(:,:,:),allocatable :: theta                  !Temperature variable
real(rp),dimension(:,:,:),allocatable :: RHS_T, RHS_Tf, RHS_Tff         !Resolved RHS terms
real(rp),dimension(:,:,:),allocatable :: dTdx, dTdy, dTdz       !Temperature gradients
real(rp),dimension(:,:,:),allocatable :: beta_scal              !Buoyancy term added to the w-momentum equation
real(rp),dimension(:,:)  ,allocatable :: L                      !Obukhov length
real(rp),dimension(:,:)  ,allocatable :: wt_avg, theta_avg      !Averaged surface flux and temperature
real(rp),dimension(:,:)  ,allocatable :: surf_flux              !Surface_flux
#ifdef LVLSET
real(rp),dimension(:,:,:),allocatable :: T_s                    !surface temperature
real(rp),dimension(:,:,:),allocatable :: phi_h, psi_h           !Stability corrections
#else
real(rp),dimension(:,:)  ,allocatable :: T_s                    !surface temperature
real(rp),dimension(:,:)  ,allocatable :: phi_h, psi_h           !Stability corrections
#endif LVLSET
real(rp),dimension(:)    ,allocatable :: ubar,vbar,wbar,Tbar    !Used in initializations, can be removed
real(rp),dimension(:)    ,allocatable :: theta_const,div_HF_var !Terms used in PID controller for temperature
real(rp),dimension(:)    ,allocatable :: sponge                 !Sponge layer strength calculated in damping_layer.f90
real(rp),dimension(:)    ,allocatable :: sponge_x               !Sponge layer strength calculated in damping_layer.f90
real(rp) :: error_int_sc, error_der_sc, error0_sc, error_sc     !Variables for PID controller in scalars_module.f90
real(rp), dimension(:,:,:), allocatable :: dtemp
#ifdef MEMORY 
real(rp), dimension(:,:,:), allocatable :: dsdx_m, dsdy_m, RHS_m
#endif 
real(rp), dimension(:,:,:), allocatable :: dsdz_m
real(rp), dimension(:), allocatable :: theta_aver_z
real(rp), dimension(:), allocatable :: scalar_bar
 
contains

!*******************************************************************************
subroutine scalars_init
!*******************************************************************************
use param, only: rp,ld,nx,ny,nz,ubc
implicit none

allocate(dTdz       (ld,ny,0:nz)); dTdz       =0.0_rp
allocate(dTdx       (ld,ny,0:nz)); dTdx       =0.0_rp
allocate(dTdy       (ld,ny,0:nz)); dTdy       =0.0_rp
allocate(theta      (ld,ny,0:nz)); theta      =1.0_rp
allocate(RHS_T      (ld,ny,0:nz)); RHS_T      =0.0_rp
allocate(beta_scal  (ld,ny,0:nz)); beta_scal  =0.0_rp
allocate(RHS_Tf     (ld,ny,0:nz)); RHS_Tf     =0.0_rp
allocate(RHS_Tff     (ld,ny,0:nz)); RHS_Tff     =0.0_rp
allocate(L          (ld,ny     )); L          =0.0_rp
allocate(theta_avg  (ld,ny     )); theta_avg  =1.0_rp
allocate(wt_avg     (ld,ny     )); wt_avg     =0.0_rp
allocate(surf_flux  (ld,ny     )); surf_flux  =0.0_rp
#ifdef LVLSET
allocate(T_s        (nx,ny,1:nz-1)); T_s      =1.0_rp
allocate(phi_h      (nx,ny,1:nz-1)); phi_h    =1.0_rp
allocate(psi_h      (nx,ny,1:nz-1)); psi_h    =0.0_rp
#else
if(coord .eq. 0) then 
 allocate(T_s        (nx,ny     )); T_s        =1.0_rp
 allocate(phi_h      (nx,ny     )); phi_h      =1.0_rp
 allocate(psi_h      (nx,ny     )); psi_h      =0.0_rp
endif 
#endif
allocate(ubar       (nz));         ubar       =0.0_rp
allocate(vbar       (nz));         vbar       =0.0_rp
allocate(wbar       (nz));         wbar       =0.0_rp
allocate(Tbar       (nz));         Tbar       =1.0_rp
allocate(theta_const(nz));         theta_const=1.0_rp
allocate(div_HF_var (nz));         div_HF_var =0.0_rp
allocate(sponge     (0:nz));       sponge     =0.0_rp
allocate(sponge_x   (ld));         sponge_x   =0.0_rp

allocate(dtemp(ld,ny,0:nz)); dtemp = 0.0_rp 
#ifdef MEMORY 
allocate(dsdx_m(ld_big,ny2,1:nz)); dsdx_m = 0.0_rp
allocate(dsdy_m(ld_big,ny2,1:nz)); dsdy_m = 0.0_rp
allocate(RHS_m(ld_big,ny2,1:nz-1)); RHS_m = 0.0_rp
#endif MEMORY
allocate(dsdz_m(ld_big,ny2,1:nz)); dsdz_m = 0.0_rp
allocate(theta_aver_z(nz)); theta_aver_z = 1.0_rp
allocate(scalar_bar(0:nz)); scalar_bar = 1.0_rp

error_int_sc=0.0_rp
error_der_sc=0.0_rp
error0_sc   =0.0_rp
error_sc    =0.0_rp

if(ubc==1) then
  call damping_layer
else 
  sponge = 0.0_rp
endif 

!Initialize surface temperature (only if not restarting)
!If restarting this will be read from a file  
if(coord.eq.0 .and. (.not.initu)) then
 T_s = theta_s1/T_scale
endif

end subroutine scalars_init

end module scalars_param
#endif SCALAR
