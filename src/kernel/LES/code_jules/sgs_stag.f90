!*******************************************************************************
subroutine sgs_stag
!*******************************************************************************
!
! Calculates turbulent (subgrid) stress for entire domain
! txx, txy, tyy, tzz (uvp-nodes) and txz, tyz (w-nodes)
! Sij values are stored on w-nodes (1:nz)
!
use param, only: rp,nproc,nx,ny,nz,dt,sgs_model,kmin
use param, only: coord,jt_total,cs_count,DYN_init
use sim_param, only: txx,txy,txz,tyy,tyz,tzz
use sim_param, only: dudx,dudy,dvdx,dvdy,dudz,dvdz,dwdx,dwdy,dwdz
use sim_param, only: Nu_t=>dummy1,dummy2,dummy3
use sim_param, only: S11=>txx,S12=>txy,S13=>txz,S22=>tyy,S23=>tyz,S33=>tzz
use sgs_param, only: Cs_opt2,lagran_dt
use mpi_defs, only: mpi_sync_down,mpi_sync_up
#ifdef SCALAR
use scalars_param, only: dTdx, dTdy, dTdz, surf_flux
use sgs_param, only: Ds_opt2,Nu_s,txTh,tyTh,tzTh
#endif SCALAR
#ifdef LVLSET
use param, only: time_level_set2
use mpi, only: mpi_wtime
#endif LVLSET
implicit none
real(rp) :: const,const2,S
#ifdef LVLSET
real(rp) :: clock_locloc(2)
#endif LVLSET
integer::i,j,k,kp

! This calculates the time between the calls to lagran_Sdep 
if ((jt_total.ge.DYN_init-cs_count+1)) then
  lagran_dt = lagran_dt + dt
endif

if(sgs_model.eq.5) then ! Lagrangian scale dependent model 
  if ((jt_total.ge.DYN_init) .and. (mod(jt_total,cs_count).eq.0)) then
    ! Track lagrangian trajectory of F_* averages
    call interpolag_Sdep
#ifdef SCALAR
    call interpolag_Sdep_Scalar
#endif SCALAR
  endif
endif 

! These are determined in wallstress and need to be buffered 
if(coord.eq.0) then 
  dummy2(:,:,1)=txz(:,:,1)
  dummy3(:,:,1)=tyz(:,:,1)
endif

! Calculate S12, S13, S23, etc.
call calc_Sij

if(sgs_model.eq.4) then ! AMD model
  if ((jt_total.ge.DYN_init) .and. (mod(jt_total,cs_count).eq.0)) then
    call amd_model
  endif
do k=1,nz
do j=1,ny
do i=1,nx
  Nu_t(i,j,k) = Cs_opt2(i,j,k)
#ifdef SCALAR
  Nu_s(i,j,k) = Ds_opt2(i,j,k)
#endif SCALAR
end do
end do
end do

elseif(sgs_model.eq.5) then ! Lagrangian scale dependent model 
  if ((jt_total.ge.DYN_init) .and. (mod(jt_total,cs_count).eq.0)) then
    call lagrange_Sdep
  endif
! Define |S| and eddy viscosity (nu_t= c_s^2 l^2 |S|) for entire domain
! stored on w-nodes (on uvp node for k=1 and 'wall' BC only) 
do k=1,nz
do j=1,ny
do i=1,nx
    S = sqrt(2._rp*(S11(i,j,k)**2 + S22(i,j,k)**2 + S33(i,j,k)**2 + 2._rp*(S12(i,j,k)**2 + S13(i,j,k)**2 + S23(i,j,k)**2)))
  Nu_t(i,j,k) = S*Cs_opt2(i,j,k)
#ifdef SCALAR
  Nu_s(i,j,k) = S*Ds_opt2(i,j,k)
#endif SCALAR
end do
end do
end do

elseif(sgs_model.eq.1) then ! Lagrangian scale dependent model 
! Define |S| and eddy viscosity (nu_t= c_s^2 l^2 |S|) for entire domain
! stored on w-nodes (on uvp node for k=1 and 'wall' BC only) 
do k=1,nz
do j=1,ny
do i=1,nx
    S = sqrt(2._rp*(S11(i,j,k)**2 + S22(i,j,k)**2 + S33(i,j,k)**2 +2._rp*(S12(i,j,k)**2 + S13(i,j,k)**2 + S23(i,j,k)**2)))
  Nu_t(i,j,k) = S*Cs_opt2(i,j,k)
#ifdef SCALAR
  Nu_s(i,j,k) = S*Ds_opt2(i,j,k)
#endif SCALAR
end do
end do
end do

endif

! Calculate txx, txy, tyy, tzz for bottom level: k=1 node (coord=0 only)
if (coord.eq.0) then
  do j=1,ny
  do i=1,nx
    const      =-2._rp * Nu_t(i,j,1)
    txx(i,j,1) = const * S11(i,j,1)  
    txy(i,j,1) = const * S12(i,j,1) 
    tyy(i,j,1) = const * S22(i,j,1)
    tzz(i,j,1) = const * S33(i,j,1)
#ifdef SCALAR 
    const2      = Nu_s(i,j,1) 
    txTh(i,j,1) = const2 * dTdx(i,j,1) 
    tyTh(i,j,1) = const2 * dTdy(i,j,1) 
#endif SCALAR
  end do
  end do
end if

! Calculate all tau for the rest of the domain
! txx, txy, tyy, tzz not needed at nz (so they aren't calculated)
! txz, tyz at nz will be done later
! txx, txy, tyy, tzz (uvp-nodes) and txz, tyz (w-nodes)
do k=kmin,nz-1
  kp=k+1
  do j=1,ny
  do i=1,nx

    const  = -(Nu_t(i,j,k) + Nu_t(i,j,kp))
    const2 = -2._rp*Nu_t(i,j,k)

    txx(i,j,k)=const * dudx(i,j,k)
    txy(i,j,k)=const * 0.5_rp* (dudy(i,j,k) + dvdx(i,j,k))
    tyy(i,j,k)=const * dvdy(i,j,k)
    tzz(i,j,k)=const * dwdz(i,j,k)
    txz(i,j,k)=const2* 0.5_rp* (dudz(i,j,k) + dwdx(i,j,k))
    tyz(i,j,k)=const2* 0.5_rp* (dvdz(i,j,k) + dwdy(i,j,k))

#ifdef SCALAR
    const  = 0.5_rp*(Nu_s(i,j,k) + Nu_s(i,j,kp))
    const2 = Nu_s(i,j,k)
    txTh(i,j,k) = const  * dTdx(i,j,k)
    tyTh(i,j,k) = const  * dTdy(i,j,k)
    tzTh(i,j,k) = const2 * dTdz(i,j,k)
#endif SCALAR
  enddo
  enddo
enddo

#ifdef SCALAR
if(coord.eq.0) then
  tzTh(1:nx,1:ny,1) = -surf_flux(1:nx,1:ny)
end if
#endif SCALAR

#ifdef LVLSET
! Calculate sgs stress in the vicinity of immersed boundary
! This must be called before copy back buffered wall stress
! at this point tij are only set for 1:nz-1
clock_locloc(1) = mpi_wtime()

call level_set_sgs

clock_locloc(2) = mpi_wtime()
time_level_set2=time_level_set2+clock_locloc(2)-clock_locloc(1)
#endif LVLSET

call mpi_sync_real_array(txz ,mpi_sync_down)
call mpi_sync_real_array(tyz ,mpi_sync_down)
call mpi_sync_real_array(tzz ,mpi_sync_up)
#ifdef SCALAR
call mpi_sync_real_array(tzTh,mpi_sync_down) 
#endif SCALAR

! Assuming stress-free lid
if (coord.eq.nproc-1) then
  txz (:,:,nz)=0.0_rp
  tyz (:,:,nz)=0.0_rp
#ifdef SCALAR
  tzTh(:,:,nz)=0.0_rp
#endif SCALAR
end if

! Copy back buffered results 
if(coord.eq.0) then
  txz(:,:,1)=dummy2(:,:,1)
  tyz(:,:,1)=dummy3(:,:,1)
endif

end subroutine sgs_stag

!*******************************************************************************
subroutine calc_Sij
!*******************************************************************************
!
! Calculate the resolved strain rate tensor, Sij = 0.5(djui - diuj)
! Stored on w-nodes for interior nodes and uv-nodes at the first gridpoint
!
use param, only: rp,nx,ny,nz,coord,kmin
use sim_param, only: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
use sim_param,only: S11=>txx,S12=>txy,S13=>txz,S22=>tyy,S23=>tyz,S33=>tzz
use grid_defs, only: gridzw 
use mpi_defs, only: mpi_sync_down
#ifdef SCALAR
use scalars_param, only: dTdz
#endif SCALAR

implicit none
real(rp) :: wx,wy,uy,vx, dz1, dz2, dzfact 
integer :: i,j,k,km

! Calculate Sij for k=1 (coord=0 only) on uvp-nodes (this level only)
! dudz and dvdz are stored on uvp-nodes for first level only, 'wall' only
! dwdx and dwdy are stored on w-nodes (always)
if (coord.eq.0) then
  do j=1,ny
  do i=1,nx
    S11(i,j,1) =         dudx(i,j,1)
    S12(i,j,1) = 0.5_rp*(dudy(i,j,1)+dvdx(i,j,1))
    wx         = 0.5_rp*(dwdx(i,j,1)+dwdx(i,j,2))
    S13(i,j,1) = 0.5_rp*(dudz(i,j,1)+wx)
    S22(i,j,1) =         dvdy(i,j,1)
    wy         = 0.5_rp*(dwdy(i,j,1)+dwdy(i,j,2))
    S23(i,j,1) = 0.5_rp*(dvdz(i,j,1)+wy)
    S33(i,j,1) =         dwdz(i,j,1)
  enddo
  enddo
endif

call mpi_sync_real_array(dwdz,mpi_sync_down)
#ifdef SCALAR 
call mpi_sync_real_array(dTdz,mpi_sync_down) !! RS: Is this required??
#endif SCALAR

! Calculate Sij for the rest of the domain
! values are stored on w-nodes
! dudz, dvdz, dwdx, dwdy are already stored on w-nodes
do k=kmin,nz
  km=k-1
  dz1 = 0.5_rp * (gridzw(k) - gridzw(k-1)) 
  dz2 = 0.5_rp * (gridzw(k+1) - gridzw(k)) 
  dzfact = (1._rp/(dz1+dz2))
 
  do j=1,ny
  do i=1,nx
    S11(i,j,k) = dzfact * (dudx(i,j,k) * dz1 + dudx(i,j,km) * dz2)
    uy         = dzfact * (dudy(i,j,k) * dz1 + dudy(i,j,km) * dz2)
    vx         = dzfact * (dvdx(i,j,k) * dz1 + dvdx(i,j,km) * dz2)
    S12(i,j,k) = 0.5_rp*(uy+vx)
    S13(i,j,k) = 0.5_rp *(dudz(i,j,k) + dwdx(i,j,k))
    S22(i,j,k) = dzfact *(dvdy(i,j,k) * dz1 + dvdy(i,j,km) * dz2)
    S23(i,j,k) = 0.5_rp *(dvdz(i,j,k) + dwdy(i,j,k))
    S33(i,j,k) = dzfact *(dwdz(i,j,k) * dz1 + dwdz(i,j,km) * dz2)
  enddo
  enddo
enddo

end subroutine calc_Sij
