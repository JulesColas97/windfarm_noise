module sgs_param
use param, only : rp

! For all sgs models
real(rp),dimension(:,:,:),allocatable :: Cs_opt2 ! (C_s)^2, Dynamic Smag coeff

! For dynamic model (5)
real(rp),dimension(:,:,:),allocatable :: F_LM,F_MM,F_QN,F_NN

#ifdef SCALAR 
real(rp),dimension(:,:,:),allocatable :: Nu_s      ! eddy viscosity
real(rp),dimension(:,:,:),allocatable :: Ds_opt2
real(rp),dimension(:,:,:),allocatable :: F_KX, F_XX, F_PY, F_YY
real(rp),dimension(:,:,:),allocatable :: txTh, tyTh, tzTh 
#endif 

real(rp),dimension(:,:),allocatable :: u_bar,v_bar,w_bar
real(rp),dimension(:,:),allocatable :: u_hat,v_hat,w_hat
real(rp),dimension(:,:),allocatable :: L11,L12,L13,L22,L23,L33
real(rp),dimension(:,:),allocatable :: Q11,Q12,Q13,Q22,Q23,Q33
real(rp),dimension(:,:),allocatable :: M11,M12,M13,M22,M23,M33
real(rp),dimension(:,:),allocatable :: S11_hat,S12_hat,S13_hat
real(rp),dimension(:,:),allocatable :: S22_hat,S23_hat,S33_hat
real(rp),dimension(:,:),allocatable :: S_S11_hat,S_S12_hat,S_S13_hat
real(rp),dimension(:,:),allocatable :: S_S22_hat,S_S23_hat,S_S33_hat

#ifdef SCALAR 
real(rp),dimension(:,:,:),allocatable :: dbdx, dbdy, dbdz
real(rp),dimension(:,:),allocatable :: Ds_opt2_2d,Ds_opt2_4d
real(rp),dimension(:,:),allocatable :: HX,XX,PY,YY
real(rp),dimension(:,:),allocatable :: H1,H2,H3,X1,X2,X3
real(rp),dimension(:,:),allocatable :: P1,P2,P3,Y1,Y2,Y3
real(rp),dimension(:,:),allocatable :: dThdx_bar,dThdy_bar,dThdz_bar
real(rp),dimension(:,:),allocatable :: dThdx_hat,dThdy_hat,dThdz_hat
real(rp),dimension(:,:),allocatable :: S_dThdx_bar,S_dThdy_bar,S_dThdz_bar
real(rp),dimension(:,:),allocatable :: S_dThdx_hat,S_dThdy_hat,S_dThdz_hat
real(rp),dimension(:,:),allocatable :: theta_bar, theta_hat
real(rp),dimension(:,:),allocatable :: beta_s
#endif 

! Used as buffer array in lagrange_Sdep
real(rp):: lagran_dt = 0._rp

contains

!*******************************************************************************
subroutine sgs_param_init
!*******************************************************************************
use param, only: ld,ny,nz,sgs_model
implicit none

! For all sgs models:
allocate(Cs_opt2(ld,ny,0:nz));  Cs_opt2=0.0_rp

! For dynamic models:
if (sgs_model .eq. 5) then
allocate(F_LM(ld,ny,0:nz) ); F_LM=0.0_rp
allocate(F_MM(ld,ny,0:nz) ); F_mm=0.0_rp
allocate(F_QN(ld,ny,0:nz) ); F_QN=0.0_rp
allocate(F_NN(ld,ny,0:nz) ); F_NN=0.0_rp
endif

#ifdef SCALAR 
allocate(Nu_s   (ld,ny,nz    )); Nu_s   =0.0_rp
allocate(Ds_opt2(ld,ny,0:nz  )); Ds_opt2=0.0_rp
allocate(txTh   (ld,ny,1:nz-1)); txTh   =0.0_rp 
allocate(tyTh   (ld,ny,1:nz-1)); tyTh   =0.0_rp
allocate(tzTh   (ld,ny,0:nz  )); tzTh   =0.0_rp

! For AMD model
if (sgs_model .eq. 4) then
allocate(dbdx  (ld,ny,0:nz  )); dbdx  =0.0_rp
allocate(dbdy  (ld,ny,0:nz  )); dbdy  =0.0_rp
allocate(dbdz  (ld,ny,0:nz  )); dbdz  =0.0_rp

! For dynamic models:
elseif (sgs_model .eq. 5) then
allocate(F_KX   (ld,ny,0:nz  )); F_KX=0.0_rp
allocate(F_XX   (ld,ny,0:nz  )); F_XX=0.0_rp
allocate(F_PY   (ld,ny,0:nz  )); F_PY=0.0_rp
allocate(F_YY   (ld,ny,0:nz  )); F_YY=0.0_rp
endif
#endif SCALAR

#ifndef MEMORY
! For dynamic models:
if (sgs_model .eq. 5) then
  call sgs_param_local
endif
#endif MEMORY

end subroutine sgs_param_init


!*******************************************************************************
subroutine sgs_param_local
!*******************************************************************************
use param, only: ld,ny
implicit none

! For dynamic models:
allocate(u_bar    (ld,ny)); u_bar    =0.0_rp
allocate(v_bar    (ld,ny)); v_bar    =0.0_rp
allocate(w_bar    (ld,ny)); w_bar    =0.0_rp
allocate(u_hat    (ld,ny)); u_hat    =0.0_rp
allocate(v_hat    (ld,ny)); v_hat    =0.0_rp
allocate(w_hat    (ld,ny)); w_hat    =0.0_rp
allocate(L11      (ld,ny)); L11      =0.0_rp
allocate(L12      (ld,ny)); L12      =0.0_rp
allocate(L13      (ld,ny)); L13      =0.0_rp
allocate(L22      (ld,ny)); L22      =0.0_rp
allocate(L23      (ld,ny)); L23      =0.0_rp
allocate(L33      (ld,ny)); L33      =0.0_rp
allocate(Q11      (ld,ny)); Q11      =0.0_rp
allocate(Q12      (ld,ny)); Q12      =0.0_rp
allocate(Q13      (ld,ny)); Q13      =0.0_rp
allocate(Q22      (ld,ny)); Q22      =0.0_rp
allocate(Q23      (ld,ny)); Q23      =0.0_rp
allocate(Q33      (ld,ny)); Q33      =0.0_rp
allocate(M11      (ld,ny)); M11      =0.0_rp
allocate(M12      (ld,ny)); M12      =0.0_rp
allocate(M13      (ld,ny)); M13      =0.0_rp
allocate(M22      (ld,ny)); M22      =0.0_rp
allocate(M23      (ld,ny)); M23      =0.0_rp
allocate(M33      (ld,ny)); M33      =0.0_rp
allocate(S11_hat  (ld,ny)); S11_hat  =0.0_rp
allocate(S12_hat  (ld,ny)); S12_hat  =0.0_rp
allocate(S13_hat  (ld,ny)); S13_hat  =0.0_rp
allocate(S22_hat  (ld,ny)); S22_hat  =0.0_rp
allocate(S23_hat  (ld,ny)); S23_hat  =0.0_rp
allocate(S33_hat  (ld,ny)); S33_hat  =0.0_rp
allocate(S_S11_hat(ld,ny)); S_S11_hat=0.0_rp
allocate(S_S12_hat(ld,ny)); S_S12_hat=0.0_rp
allocate(S_S13_hat(ld,ny)); S_S13_hat=0.0_rp
allocate(S_S22_hat(ld,ny)); S_S22_hat=0.0_rp
allocate(S_S23_hat(ld,ny)); S_S23_hat=0.0_rp
allocate(S_S33_hat(ld,ny)); S_S33_hat=0.0_rp

#ifdef SCALAR 
allocate(Ds_opt2_2d (ld,ny)); Ds_opt2_2d =0.0_rp
allocate(Ds_opt2_4d (ld,ny)); Ds_opt2_4d =0.0_rp
allocate(HX         (ld,ny)); HX         =0.0_rp
allocate(XX         (ld,ny)); XX         =0.0_rp
allocate(PY         (ld,ny)); PY         =0.0_rp
allocate(YY         (ld,ny)); YY         =0.0_rp
allocate(H1         (ld,ny)); H1         =0.0_rp
allocate(H2         (ld,ny)); H2         =0.0_rp
allocate(H3         (ld,ny)); H3         =0.0_rp
allocate(X1         (ld,ny)); X1         =0.0_rp
allocate(X2         (ld,ny)); X2         =0.0_rp
allocate(X3         (ld,ny)); X3         =0.0_rp
allocate(P1         (ld,ny)); P1         =0.0_rp
allocate(P2         (ld,ny)); P2         =0.0_rp
allocate(P3         (ld,ny)); P3         =0.0_rp
allocate(Y1         (ld,ny)); Y1         =0.0_rp
allocate(Y2         (ld,ny)); Y2         =0.0_rp
allocate(Y3         (ld,ny)); Y3         =0.0_rp
allocate(dThdx_bar  (ld,ny)); dThdx_bar  =0.0_rp
allocate(dThdy_bar  (ld,ny)); dThdy_bar  =0.0_rp
allocate(dThdz_bar  (ld,ny)); dThdz_bar  =0.0_rp
allocate(dThdx_hat  (ld,ny)); dThdx_hat  =0.0_rp
allocate(dThdy_hat  (ld,ny)); dThdy_hat  =0.0_rp
allocate(dThdz_hat  (ld,ny)); dThdz_hat  =0.0_rp
allocate(S_dThdx_bar(ld,ny)); S_dThdx_bar=0.0_rp
allocate(S_dThdy_bar(ld,ny)); S_dThdy_bar=0.0_rp
allocate(S_dThdz_bar(ld,ny)); S_dThdz_bar=0.0_rp
allocate(S_dThdx_hat(ld,ny)); S_dThdx_hat=0.0_rp
allocate(S_dThdy_hat(ld,ny)); S_dThdy_hat=0.0_rp
allocate(S_dThdz_hat(ld,ny)); S_dThdz_hat=0.0_rp
allocate(theta_bar  (ld,ny)); theta_bar  =0.0_rp
allocate(theta_hat  (ld,ny)); theta_hat  =0.0_rp
allocate(beta_s     (ld,ny)); beta_s     =0.0_rp
#endif SCALAR

end subroutine sgs_param_local

end module sgs_param
