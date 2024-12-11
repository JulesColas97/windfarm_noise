subroutine lagrange_Sdep
use param
use sim_param,only:u,v,w
use sgs_param,only:sgs_param_local,F_LM,F_MM,F_QN,F_NN,Cs_opt2,lagran_dt
use grid_defs,only: gridzw
use sgs_param,only:u_bar,v_bar,w_bar
use sgs_param,only:u_hat,v_hat,w_hat
use sgs_param,only:L11,L12,L13,L22,L23,L33
use sgs_param,only:Q11,Q12,Q13,Q22,Q23,Q33
use sgs_param,only:M11,M12,M13,M22,M23,M33
use sgs_param,only:S11_hat,S12_hat,S13_hat
use sgs_param,only:S22_hat,S23_hat,S33_hat
use sgs_param,only:S_S11_hat,S_S12_hat,S_S13_hat
use sgs_param,only:S_S22_hat,S_S23_hat,S_S33_hat

#ifdef SCALAR 
use sgs_param,only:F_KX,F_XX,F_PY,F_YY,Ds_opt2

use sgs_param,only:Ds_opt2_2d,Ds_opt2_4d
use sgs_param,only:HX,XX,PY,YY
use sgs_param,only:H1,H2,H3,X1,X2,X3
use sgs_param,only:P1,P2,P3,Y1,Y2,Y3
use sgs_param,only:dThdx_bar,dThdy_bar,dThdz_bar
use sgs_param,only:dThdx_hat,dThdy_hat,dThdz_hat
use sgs_param,only:S_dThdx_bar,S_dThdy_bar,S_dThdz_bar
use sgs_param,only:S_dThdx_hat,S_dThdy_hat,S_dThdz_hat
use sgs_param,only:theta_bar, theta_hat,beta_s
use scalars_param,only:theta,dTdx,dTdy,dTdz
#endif 

! Used for memory saving and efficiency:
! Watch out:Data of variable destroyed when other variable is used
use sim_param,only:S11=>txx,S12=>txy,S13=>txz,S22=>tyy,S23=>tyz,S33=>tzz
use sim_param,only:Tn        =>dummy_plane
use sim_param,only:beta      =>dummy_plane
use sim_param,only:epsi      =>dummy_plane
use sim_param,only:S         =>dummy_plane
use sim_param,only:Cs_opt2_4d=>dummy_plane
use sgs_param,only:LM        =>u_bar
use sgs_param,only:MM        =>v_bar
use sgs_param,only:QN        =>w_bar
use sgs_param,only:NN        =>u_hat
use sgs_param,only:Cs_opt2_2d=>v_hat
use sgs_param,only:N11      =>M11  ,N12      =>M12  ,N13=>M13
use sgs_param,only:N22      =>M22  ,N23      =>M23  ,N33=>M33
use sgs_param,only:S11_bar  =>u_bar,S12_bar  =>v_bar,S13_bar=>w_bar
use sgs_param,only:S22_bar  =>u_hat,S23_bar  =>v_hat,S33_bar=>w_hat
use sgs_param,only:S_S11_bar=>M11  ,S_S12_bar=>M12  ,S_S13_bar=>M13
use sgs_param,only:S_S22_bar=>M22  ,S_S23_bar=>M23  ,S_S33_bar=>M33
use filters ,only:delta
implicit none
real(rp):: const,constconst,const2,opftdelta,fractus
real(rp),parameter:: zero=1.e-24_rp
real(rp),parameter:: opftime = 1.5_rp ! (Meneveau,Lund,Cabot; JFM 1996)
real(rp),parameter:: powcoeff = -1._rp/8._rp
real(rp),parameter:: tf1=2._rp
real(rp),parameter:: tf2=4._rp
real(rp),parameter:: tf3=1._rp/(tf1*tf2)
real(rp),parameter:: tf1_2=tf1**2
real(rp),parameter:: tf2_2=tf2**2
logical, save :: F_LM_MM_init = .false.
real(rp) :: dz1, dz2, dzfact
#ifdef CPS
integer :: k 
#else
integer :: k,istart,iend
#endif

#ifdef MEMORY
call sgs_param_local
#endif

! Set coefficients
opftdelta = opftime*delta
fractus= 1._rp/real(ny*nx,kind=rp)
const = 2*(delta**2)
constconst=const*const
const2 = delta**2  

! Update Cs_opt2 per horizontal slice
! See Bou-Zeid et al, Physics of Fluids 17, 025105 (2005)
! See ..
do k = 1,nz
  ! Calculate Lij
  ! Interp u,v,w onto w-nodes and store result as u_bar,v_bar,w_bar
  ! (except for very first level which should be on uvp-nodes)
  if(( coord.eq.0 ) .and. (k.eq.1) ) then  ! uvp-nodes
   u_bar(:,:)     = u(:,:,1)
   v_bar(:,:)     = v(:,:,1)
   w_bar(:,:)     = 0.5_rp*w(:,:,2)

#ifdef SCALAR 
   theta_bar(:, :)= theta(:, :, 1) 
#endif 
  else  ! w-nodes
   dz1 = 0.5_rp * (gridzw(k) - gridzw(k-1)) 
   dz2 = 0.5_rp * (gridzw(k+1) - gridzw(k)) 
   dzfact = (1._rp/(dz1+dz2))
   u_bar(:,:)      = dzfact*(u(:,:,k)*dz1 + u(:,:,k-1)*dz2) 
   v_bar(:,:)      = dzfact*(v(:,:,k)*dz1 + v(:,:,k-1)*dz2)  
   w_bar(:,:)      = w(:,:,k)

#ifdef SCALAR
   theta_bar(:, :) = dzfact*(theta(:,:,k)*dz1 + theta(:,:,k-1)*dz2)
#endif
  end if

  ! Set up L** based on unfiltered u_bar, etc values 
  L11=u_bar*u_bar
  L12=u_bar*v_bar
  L13=u_bar*w_bar
  L22=v_bar*v_bar
  L23=v_bar*w_bar
  L33=w_bar*w_bar
#ifdef SCALAR 
  !RS: This part needs to be tested
  H1 =theta_bar*u_bar
  H2 =theta_bar*v_bar
  H3 =theta_bar*w_bar
#endif

  ! Filtering u_bar and u_hat.
  ! Afterwards u_bar 2 delta filter, u_hat 4 delta filter
  call filtering(u_bar,u_hat)
  call filtering(v_bar,v_hat)
  call filtering(w_bar,w_hat)
#ifdef SCALAR 
  !RS: This part needs to be tested
  call filtering(theta_bar,theta_hat)
#endif

  ! L** is 2 delta filter, Q** 4 delta filter
  call filtering(L11,Q11)
  L11 = L11 - u_bar*u_bar
  Q11 = Q11 - u_hat*u_hat
  call filtering(L12,Q12)
  L12 = L12 - u_bar*v_bar
  Q12 = Q12 - u_hat*v_hat
  call filtering(L13,Q13)
  L13 = L13 - u_bar*w_bar
  Q13 = Q13 - u_hat*w_hat
  call filtering(L22,Q22)
  L22 = L22 - v_bar*v_bar
  Q22 = Q22 - v_hat*v_hat
  call filtering(L23,Q23)
  L23 = L23 - v_bar*w_bar
  Q23 = Q23 - v_hat*w_hat
  call filtering(L33,Q33)
  L33 = L33 - w_bar*w_bar
  Q33 = Q33 - w_hat*w_hat

#ifdef SCALAR
  !RS: This part needs to be tested
  call filtering(H1,P1)
  H1 = H1 - theta_bar*u_bar
  P1 = P1 - theta_hat*u_hat
  call filtering(H2,P2)
  H2 = H2 - theta_bar*v_bar
  P2 = P2 - theta_hat*v_hat
  call filtering(H3,P3)
  H3 = H3 - theta_bar*w_bar
  P3 = P3 - theta_hat*w_hat
#endif 

  ! Select Sij for this level for test-filtering, saving results as Sij_bar
  S11_bar(:,:) = S11(:,:,k)  
  S12_bar(:,:) = S12(:,:,k)  
  S13_bar(:,:) = S13(:,:,k)  
  S22_bar(:,:) = S22(:,:,k)  
  S23_bar(:,:) = S23(:,:,k)  
  S33_bar(:,:) = S33(:,:,k)

  call filtering(S11_bar,S11_hat)
  call filtering(S12_bar,S12_hat)
  call filtering(S13_bar,S13_hat)
  call filtering(S22_bar,S22_hat)
  call filtering(S23_bar,S23_hat)
  call filtering(S33_bar,S33_hat)

  ! Calculate |S|
  S(:,:) = sqrt(2*(S11(:,:,k)**2+S22(:,:,k)**2+S33(:,:,k)**2+&
                2*(S12(:,:,k)**2+S13(:,:,k)**2+S23(:,:,k)**2)))

  ! Calculate |S|Sij then test-filter this quantity 
  S_S11_bar(:,:) = S(:,:)*S11(:,:,k)
  S_S12_bar(:,:) = S(:,:)*S12(:,:,k)
  S_S13_bar(:,:) = S(:,:)*S13(:,:,k)
  S_S22_bar(:,:) = S(:,:)*S22(:,:,k)
  S_S23_bar(:,:) = S(:,:)*S23(:,:,k)
  S_S33_bar(:,:) = S(:,:)*S33(:,:,k)

  call filtering(S_S11_bar,S_S11_hat)
  call filtering(S_S12_bar,S_S12_hat)
  call filtering(S_S13_bar,S_S13_hat)
  call filtering(S_S22_bar,S_S22_hat)
  call filtering(S_S23_bar,S_S23_hat)
  call filtering(S_S33_bar,S_S33_hat)

#ifdef SCALAR
  !RS: This part needs to be tested
  if(coord.eq.0 .and. k.eq.1) then
   dThdx_bar = dTdx(:,:,1)
   dThdy_bar = dTdy(:,:,1)
  else
   dz1 = 0.5_rp * (gridzw(k) - gridzw(k-1)) 
   dz2 = 0.5_rp * (gridzw(k+1) - gridzw(k)) 
   dzfact = (1._rp/(dz1+dz2))
   
   dThdx_bar = dzfact * (dTdx(:,:,k) * dz1 + dTdx(:,:,k-1) * dz2)
   dThdy_bar = dzfact * (dTdy(:,:,k) * dz1 + dTdy(:,:,k-1) * dz2)
  endif
  dThdz_bar = dTdz(:,:,k)

  !Calculate Xi 
  S_dThdx_bar(:,:) = S(:,:) * dThdx_bar(:,:)
  S_dThdy_bar(:,:) = S(:,:) * dThdy_bar(:,:)
  S_dThdz_bar(:,:) = S(:,:) * dThdz_bar(:,:)

  call filtering(dThdx_bar,dThdx_hat)
  call filtering(dThdy_bar,dThdy_hat)
  call filtering(dThdz_bar,dThdz_hat)

  call filtering(S_dThdx_bar,S_dThdx_hat)
  call filtering(S_dThdy_bar,S_dThdy_hat)
  call filtering(S_dThdz_bar,S_dThdz_hat)
#endif 

  ! Calculate |S_bar| (2 delta filtered S) * normalization factor 
  S = tf1_2*sqrt(2*(S11_bar**2 + S22_bar**2 + S33_bar**2 +&
                 2*(S12_bar**2 + S13_bar**2 + S23_bar**2)))

  ! Calculate Mij (const directly applied at LM,MM)
  M11 = S_S11_bar - S*S11_bar
  M12 = S_S12_bar - S*S12_bar
  M13 = S_S13_bar - S*S13_bar
  M22 = S_S22_bar - S*S22_bar
  M23 = S_S23_bar - S*S23_bar
  M33 = S_S33_bar - S*S33_bar

#ifdef SCALAR
  X1  = const2*(S_dThdx_bar - S*dThdx_bar)
  X2  = const2*(S_dThdy_bar - S*dThdy_bar)
  X3  = const2*(S_dThdz_bar - S*dThdz_bar)
#endif 

  ! Calculate LijMij,MijMij
  LM = const     *(L11*M11+L22*M22+L33*M33+2*(L12*M12+L13*M13+L23*M23))
  MM = constconst*(M11*M11+M22*M22+M33*M33+2*(M12*M12+M13*M13+M23*M23))

  ! Calculate |S_hat| (4 delta filtered S) * normalization factor 
  S = tf2_2*sqrt(2*(S11_hat**2 + S22_hat**2 + S33_hat**2 +&
                 2*(S12_hat**2 + S13_hat**2 + S23_hat**2)))

  ! Calculate Nij (const directly applied at QN,NN)
  N11 = S_S11_hat - S*S11_hat
  N12 = S_S12_hat - S*S12_hat
  N13 = S_S13_hat - S*S13_hat
  N22 = S_S22_hat - S*S22_hat
  N23 = S_S23_hat - S*S23_hat
  N33 = S_S33_hat - S*S33_hat

#ifdef SCALAR
  Y1 = const2*(S_dThdx_hat - S*dThdx_hat)
  Y2 = const2*(S_dThdy_hat - S*dThdy_hat)
  Y3 = const2*(S_dThdz_hat - S*dThdz_hat)
#endif

  ! Calculate QijNij,NijNij
  QN = const     *(Q11*N11+Q22*N22+Q33*N33+2*(Q12*N12+Q13*N13+Q23*N23))
  NN = constconst*(N11*N11+N22*N22+N33*N33+2*(N12*N12+N13*N13+N23*N23))

#ifdef SCALAR
  ! Calculate KiXi, XiXi, PiYi, YiYi  
  HX = H1*X1 + H2*X2 + H3*X3
  XX = X1**2 + X2**2 + X3**2
  PY = P1*Y1 + P2*Y2 + P3*Y3
  YY = Y1**2 + Y2**2 + Y3**2
#endif 

  ! Using local time counter to reinitialize SGS quantities when restarting
  if (initu .and. (.not.F_LM_MM_init)) then
    if ((k.eq.1) .and. (coord.eq.0)) write(*,*) 'Reading F_LM, F_MM, F_QN, F_NN from file'
    F_LM_MM_init = .true.
  else
   if ((.not. F_LM_MM_init) .and. (jt_total.eq.cs_count.or.jt_total.eq.DYN_init)) then
    if ((k.eq.1) .and. (coord.eq.0)) write(*,*) 'F_MM, F_LM, F_NN and F_QN initialized'
    F_MM (:,:,k)     = MM
    F_LM (:,:,k)     = 0.03_rp*MM
    F_MM(ld-1:ld,:,k)= 1._rp !RS maybe not used
    F_LM(ld-1:ld,:,k)= 1._rp !RS maybe not used
    F_NN (:,:,k)     = NN
    F_QN (:,:,k)     = 0.03_rp*NN
    F_NN(ld-1:ld,:,k)= 1._rp !RS maybe not used
    F_QN(ld-1:ld,:,k)= 1._rp !RS maybe not used
#ifdef SCALAR 
    F_XX(:,:,k)      = XX/Pr
    F_KX(:,:,k)      = 0.03_rp*XX/Pr 
    F_XX(ld-1:ld,:,k)= 1._rp/Pr !RS maybe not used
    F_KX(ld-1:ld,:,k)= 1._rp/Pr !RS maybe not used

    F_YY(:,:,k)      = YY/Pr
    F_PY(:,:,k)      = 0.03_rp*YY/Pr 
    F_YY(ld-1:ld,:,k)= 1._rp/Pr !RS maybe not used
    F_PY(ld-1:ld,:,k)= 1._rp/Pr !RS maybe not used
#endif 
    if (k.eq.nz) F_LM_MM_init = .true.
   end if
  end if

#ifndef CPS  
  if (inflow) then
   iend   = floor (fringe_region_end * nx + 1._rp)
   iend   = min(nx,iend)
   istart = floor ((fringe_region_end - fringe_region_len) * nx + 1._rp)
   istart = modulo(istart-1,nx) ! prevents problems at start fringe region 
       
   ! Put these values 0 for laminar flow            
   LM  (istart:iend,1:ny)  = 0._rp
   MM  (istart:iend,1:ny)  = 0._rp
   F_LM(istart:iend,1:ny,k)= 0._rp ! gives Tn=zero first step, i.e. epsi=1
   F_MM(istart:iend,1:ny,k)= 0._rp
   QN  (istart:iend,1:ny)  = 0._rp
   NN  (istart:iend,1:ny)  = 0._rp  
   F_QN(istart:iend,1:ny,k)= 0._rp
   F_NN(istart:iend,1:ny,k)= 0._rp

#ifdef SCALAR
   HX  (istart:iend,1:ny)  = 0._rp
   XX  (istart:iend,1:ny)  = 0._rp
   F_KX(istart:iend,1:ny,k)= 0._rp 
   F_XX(istart:iend,1:ny,k)= 0._rp
   PY  (istart:iend,1:ny)  = 0._rp
   YY  (istart:iend,1:ny)  = 0._rp
   F_PY(istart:iend,1:ny,k)= 0._rp
   F_YY(istart:iend,1:ny,k)= 0._rp
#endif 
  endif
#endif

  ! Update running averages (F_LM, F_MM)
  ! Determine averaging timescale (for 2-delta filter)
  Tn(:,:)   = max(F_LM(:,:,k) * F_MM(:,:,k), zero)
  Tn(:,:)   = lagran_dt/(opftdelta*(Tn(:,:)**powcoeff))
  epsi(:,:) = Tn(:,:) / (1._rp+Tn(:,:))

  ! Calculate new running average = old*(1-epsi) + instantaneous*epsi
  F_LM(:,:,k)= max(0._rp,(epsi*LM + (1._rp-epsi)*F_LM(:,:,k)))
  F_MM(:,:,k)= max(zero ,(epsi*MM + (1._rp-epsi)*F_MM(:,:,k)))

#ifdef SCALAR 
  F_KX(:,:,k) = max(0._rp,(epsi*HX + (1._rp-epsi)*F_KX(:,:,k)))
  F_XX(:,:,k) = max(zero ,(epsi*XX + (1._rp-epsi)*F_XX(:,:,k)))
#endif

  ! Calculate Cs_opt2 (for 2-delta filter)
  Cs_opt2_2d(:,:) = F_LM(:,:,k)/F_MM(:,:,k)

#ifdef SCALAR 
  Ds_opt2_2d(:,:) = F_KX(:,:,k)/F_XX(:,:,k)
#endif 

  ! Update running averages (F_QN, F_NN)
  ! Determine averaging timescale (for 4-delta filter)
  Tn(:,:)   = max( F_QN(:,:,k)*F_NN(:,:,k), zero)
  Tn(:,:)   = lagran_dt/(opftdelta*(Tn(:,:)**powcoeff)) 
  epsi(:,:) = Tn(:,:) / (1._rp+Tn(:,:))

  ! Calculate new running average = old*(1-epsi) + instantaneous*epsi
  F_QN(:,:,k)=max(0._rp,(epsi*QN + (1._rp-epsi)*F_QN(:,:,k)))
  F_NN(:,:,k)=max(zero ,(epsi*NN + (1._rp-epsi)*F_NN(:,:,k)))

#ifdef SCALAR 
  F_PY(:,:,k)=max(0._rp,(epsi*PY + (1._rp-epsi)*F_PY(:,:,k)))
  F_YY(:,:,k)=max(zero ,(epsi*YY + (1._rp-epsi)*F_YY(:,:,k)))
#endif

  ! Calculate Cs_opt2 (for 4-delta filter)
  Cs_opt2_4d(:,:)=F_QN(:,:,k)/F_NN(:,:,k)

  beta(:,:)=max(Cs_opt2_4d(:,:)/Cs_opt2_2d(:,:),clip)

  if ((coord.eq.nproc-1).and.(k.eq.nz)) then
   beta(:,:)=1._rp
  endif

  ! Clip Beta and set Cs_opt2 for each point in the plane
  ! Multiply with delta**2 (const2) 
  Cs_opt2(:,:,k)=max(zero,Cs_opt2_2d(:,:)/Beta(:,:)*const2) 
#ifdef SCALAR 
  Ds_opt2_4d(:,:)    = F_PY(:,:,k)/F_YY(:,:,k)

  beta_s(:,:) = max(Ds_opt2_4d(:,:)/Ds_opt2_2d(:,:),clip)

  if ((coord.eq.nproc-1).and.(k.eq.nz)) then
   beta_s(:,:)=1._rp
  endif

  Ds_opt2(:,:,k)=max(zero,Ds_opt2_2d(:,:)/beta_s(:,:)*const2) 
#endif 

end do

#ifdef LVLSET
! Zero Cs_opt2 inside objects
call level_set_Cs_lag_dyn
#endif

! Reset variable for use during next set of cs_count timesteps
lagran_dt = 0.0_rp

#ifdef MEMORY
deallocate(u_bar,v_bar,w_bar,u_hat,v_hat,w_hat)
deallocate(L11,L12,L13,L22,L23,L33)
deallocate(Q11,Q12,Q13,Q22,Q23,Q33)
deallocate(M11,M12,M13,M22,M23,M33)
deallocate(S11_hat,S12_hat,S13_hat)
deallocate(S22_hat,S23_hat,S33_hat)
deallocate(S_S11_hat,S_S12_hat,S_S13_hat)
deallocate(S_S22_hat,S_S23_hat,S_S33_hat)
#ifdef SCALAR
deallocate(Ds_opt2_2d,Ds_opt2_4d)
deallocate(HX,XX,PY,YY,H1,H2,H3,X1,X2,X3)
deallocate(P1,P2,P3,Y1,Y2,Y3)
deallocate(dThdx_bar,dThdy_bar,dThdz_bar)
deallocate(dThdx_hat,dThdy_hat,dThdz_hat)
deallocate(S_dThdx_bar,S_dThdy_bar,S_dThdz_bar)
deallocate(S_dThdx_hat,S_dThdy_hat,S_dThdz_hat)
deallocate(theta_bar,theta_hat,beta_s)
#endif
#endif

end subroutine lagrange_Sdep
