!*******************************************************************************
subroutine amd_model
!*******************************************************************************
! SGS model: anisotropic minimum dissipation (AMD) model
! Ref. Abkar & Moin, 2017, BLM, 165, 405
use param, only: rp,dx,dy,dz,nx,ny,nz,coord,kmin
use sgs_param, only: Cs_opt2
use grid_defs, only: gridz, gridzw 
use sim_param, only: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
#ifdef SCALAR
use param, only: ld
use sgs_param, only: Ds_opt2,dbdx,dbdy,dbdz
use scalars_param, only: beta_scal,dTdx,dTdy,dTdz
#endif SCALAR

! Used for memory saving and efficiency:
! Watch out: Data of variable destroyed when other variable is used
use sim_param,only: S11=>txx,S12=>txy,S13=>txz,S22=>tyy,S23=>tyz,S33=>tzz
implicit none
real(rp) :: rhs1,den,tf1,tf2,tf3,ux,uy,uz,vx,vy,vz,wx,wy,wz
#ifdef SCALAR 
real(rp) :: bx, by, bz, theta_x, theta_y, theta_z
#endif SCALAR 
real(rp), parameter:: zero=1.e-24_rp
!real(rp) :: dz1, dz2, dzfact
integer :: i, j
integer :: k,km
real(rp) :: dz1, dz2, dzfact

! Modified Poincare constant
tf1 = -(1._rp/12._rp) * (dx**2)
tf2 = -(1._rp/12._rp) * (dy**2)
tf3 = -(1._rp/3._rp)  * (dz**2)

#ifdef SCALAR 
call filt_da (beta_scal, dbdx, dbdy,2)
call ddz_w (beta_scal(1:ld,1:ny,0:nz),dbdz(1:ld,1:ny,0:nz))
#endif SCALAR 

if(coord==0) then
  ! The lines including stretched grid modifications are commented out / 
  ! and replaced for now. For specific cases it might be necessary to change
  !dz1 = (gridzw(2) - gridzw(1))
  !tf3 = -(1._rp/3._rp) * (dz1**2)

  do j=1,ny
  do i=1,nx

    ux = dudx(i,j,1); uy = dudy(i,j,1); uz = dudz(i,j,1)
    vx = dvdx(i,j,1); vy = dvdy(i,j,1); vz = dvdz(i,j,1)
    wx = 0.5_rp*(dwdx(i,j,2))
    wy = 0.5_rp*(dwdy(i,j,2))
    wz = dwdz(i,j,1)

    den = (ux**2 + vx**2 + wx**2 + uy**2 + vy**2 + wy**2 + uz**2 + vz**2 +wz**2)

    rhs1 = tf1 * (S11(i,j,1)*ux*ux + S12(i,j,1)*2*ux*vx + S13(i,j,1)*2*ux*wx  &
               +  S22(i,j,1)*vx*vx + S23(i,j,1)*2*vx*wx + S33(i,j,1)*  wx*wx) &
         + tf2 * (S11(i,j,1)*uy*uy + S12(i,j,1)*2*uy*vy + S13(i,j,1)*2*uy*wy  &
               +  S22(i,j,1)*vy*vy + S23(i,j,1)*2*vy*wy + S33(i,j,1)*  wy*wy) &
         + tf3 * (S11(i,j,1)*uz*uz + S12(i,j,1)*2*uz*vz + S13(i,j,1)*2*uz*wz  &
               +  S22(i,j,1)*vz*vz + S23(i,j,1)*2*vz*wz + S33(i,j,1)*  wz*wz)

#ifdef SCALAR 
    bx = 0.5_rp * (dbdx(i,j,1) + dbdx(i,j,2))
    by = 0.5_rp * (dbdy(i,j,1) + dbdy(i,j,2))
    bz = dbdz(i,j,1)
    theta_x = dTdx(i,j,1)
    theta_y = dTdy(i,j,1)
    theta_z = dTdz(i,j,1)

#ifndef LVLSET
    rhs1 = rhs1+(-tf1*wx*bx - tf2*wy*by - tf3*wz*bz)
#endif LVLSET 

#endif SCALAR  

    Cs_opt2(i,j,1) = max(rhs1/(den+zero),zero)

#ifdef SCALAR 
    rhs1  = tf1 * (theta_x * (ux * theta_x + vx * theta_y + wx * theta_z)) &
          + tf2 * (theta_y * (uy * theta_x + vy * theta_y + wy * theta_z)) &
          + tf3 * (theta_z * (uz * theta_x + vz * theta_y + wz * theta_z))

    den = (theta_x**2 + theta_y**2 + theta_z**2)
    Ds_opt2(i,j,1) = max(rhs1/(den + zero),zero)
#endif SCALAR 

  enddo
  enddo
endif

do k=kmin,nz
  km=k-1
  ! Srinidhi : Feb 20, 2022 
  ! Technically speaking the value of tf3 would need to be fixed based on the local grid spacing. 
  ! However, that starts adding unnecessary shear stress in regions above the boundary layer (where the flow is nearly laminar) 
  ! Smagorinsky approximation takes into account an 'averaged delta' i.e. (dx * dy * dz)**(1/3)
  ! This basic Smagorinsky idealization would not work in a region where dz values are of the O(100) meters. 
  ! Therefore, to avoid the unecessary shear stress being added in non-turbulent regions I have used dz value in the uniform grid
  ! region. 
  ! If and when someone figures out what goes wrong, this can be modified. However, it is worth noting that present implementation 
  ! is not technically wrong. It is just an alternative way of looking at Smagorinsky approximation. 

  !dz1 = (gridzw(k+1) - gridzw(k))
  !tf3 = -(1._rp/3._rp) * (dz1**2)

  dz2 = 0.5_rp * (gridzw(k+1) - gridzw(k))
  dz1 = 0.5_rp * (gridzw(k) - gridzw(k-1))
  
  dzfact=(1._rp/(dz1+dz2))
  
  do j=1,ny
  do i=1,nx
    ux = dzfact * (dz1 * dudx(i,j,k) + dz2 * dudx(i,j,km))
    uy = dzfact * (dz1 * dudy(i,j,k) + dz2 * dudy(i,j,km))
    vx = dzfact * (dz1 * dvdx(i,j,k) + dz2 * dvdx(i,j,km))
    vy = dzfact * (dz1 * dvdy(i,j,k) + dz2 * dvdy(i,j,km))
    wz = dzfact * (dz1 * dwdz(i,j,k) + dz2 * dwdz(i,j,km))
    uz = dudz(i,j,k); vz = dvdz(i,j,k)
    wx = dwdx(i,j,k); wy = dwdy(i,j,k)

    den = (ux**2 + vx**2 + wx**2 + uy**2 + vy**2 + wy**2 + uz**2 + vz**2 +wz**2)

    rhs1 = tf1 * (S11(i,j,k)*ux*ux + S12(i,j,k)*2*ux*vx + S13(i,j,k)*2*ux*wx   &
               +  S22(i,j,k)*vx*vx + S23(i,j,k)*2*vx*wx + S33(i,j,k)*  wx*wx)  &
         + tf2 * (S11(i,j,k)*uy*uy + S12(i,j,k)*2*uy*vy + S13(i,j,k)*2*uy*wy   &
               +  S22(i,j,k)*vy*vy + S23(i,j,k)*2*vy*wy + S33(i,j,k)*  wy*wy)  &
         + tf3 * (S11(i,j,k)*uz*uz + S12(i,j,k)*2*uz*vz + S13(i,j,k)*2*uz*wz   &
               +  S22(i,j,k)*vz*vz + S23(i,j,k)*2*vz*wz + S33(i,j,k)*  wz*wz)

#ifdef SCALAR 
    bx = dbdx(i,j,k)
    by = dbdy(i,j,k)
    bz = dzfact * (dz1 * dbdz(i,j,k) + dz2 * dbdz(i,j,km))
    theta_x = dzfact * (dz1 * dTdx(i,j,k) + dz2 * dTdx(i,j,km))
    theta_y = dzfact * (dz1 * dTdy(i,j,k) + dz2 * dTdy(i,j,km))
    theta_z = dTdz(i,j,k)
#ifndef LVLSET 
! Srinidhi : Feb 20, 2022
! Level set crashes with the addition of this term (this term is not significant) anyway  
! For the time being the below term is added only if the level set is absent 
    rhs1 = rhs1 +(-tf1 * wx * bx - tf2 * wy * by - tf3 * wz * bz)
#endif LVLSET 
#endif 

    Cs_opt2(i,j,k) = max(rhs1/(den + zero),zero)

#ifdef SCALAR 
    rhs1  = tf1 * (theta_x * (ux * theta_x + vx * theta_y + wx * theta_z)) &
          + tf2 * (theta_y * (uy * theta_x + vy * theta_y + wy * theta_z)) &
          + tf3 * (theta_z * (uz * theta_x + vz * theta_y + wz * theta_z))

    den = (theta_x**2 + theta_y**2 + theta_z**2)
    Ds_opt2(i,j,k) = max(rhs1/(den + zero),zero)
#endif SCALAR 

  enddo
  enddo
enddo

#ifdef LVLSET
! Zero Cs_opt2 inside objects
call level_set_Cs_amd
#endif

end subroutine amd_model
