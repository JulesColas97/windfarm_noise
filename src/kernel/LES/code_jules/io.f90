!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Performs io operations and calculates the statistical data.
! Contains the following subroutines:
!
! - output_loop     : Driver routine, determines when to calculate statistics or store data
! - inst_write      : Makes snapshots of the flow field
! - checkpoint      : Stores field field and statistical data
!
! The routines below are used to calculate the statistical data during the simulation
! - tavg_init       : Initilizes arrays used in computation of statistical data
! - tavg_compute    : Performs the actual computation of the statistical data
! - tavg_compute_tau: Calculates subgrid scale stress stress
!                     These arrays are re-used in the code, so calculate stats just after
!                     Determining txx,txy,tyz,tyy,tyz,tzz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Driver routine, determines when to calculate statistics or to store flow snapshots
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_loop
use param, only : checkpoint_data,checkpoint_nskip,jt_total
use param, only : tavg_calc,tavg_nstart,tavg_start_time,tavg_end_time,total_time
use param, only : domain_calc,domain_nstart,domain_nskip,tavg_nskip
implicit none

! Store data both on numered tavg* and vel* files as well as normal tavg* and vel* files
if( checkpoint_data ) then
  if ( modulo (jt_total, checkpoint_nskip) == 0) call checkpoint(2) ! Numbered files
  if ( modulo (jt_total, checkpoint_nskip) == 0) call checkpoint(1) ! Regular files
endif

! Determine whether "tavg_compute" should be activated
if (tavg_calc) then
 if (total_time.ge.tavg_start_time .and. total_time.le.tavg_end_time) then
  if(mod(jt_total-tavg_nstart,tavg_nskip).eq.0.and.jt_total-tavg_nstart.gt.0) then
   call tavg_compute
  endif
 endif  ! after nstart
endif  ! tavg_calc

!  Determine if instantaneous domain velocities are to be recorded
if(domain_calc) then
  if(jt_total.ge.domain_nstart.and.(mod(jt_total-domain_nstart,domain_nskip).eq.0)) then
    call inst_write
  endif
endif

end subroutine output_loop


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This subroutine is used to write instantaneous snapshots of the flow
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine inst_write
use param, only : itype,coord,rp,path,nx,ny,kmax
use sim_param, only : u,v,w
#ifdef SCALAR 
use scalars_param, only: theta 
#endif 
implicit none
integer :: kplane 
character (64) :: folder

!  Write instantaneous snapshot for the entire domain
if(itype==3) then

  folder=trim(path)//'output'
  call write_3Dfield(2,    u(1:nx,1:ny,1:kmax),'field_u'    ,folder,1,kmax)
  call write_3Dfield(2,    v(1:nx,1:ny,1:kmax),'field_v'    ,folder,1,kmax)
  call write_3Dfield(2,    w(1:nx,1:ny,1:kmax),'field_w'    ,folder,1,kmax)
#ifdef SCALAR
  call write_3Dfield(2,theta(1:nx,1:ny,1:kmax),'field_theta',folder,1,kmax)
#endif
!  Write instantaneous snapshot for horizontal plane
elseif(itype==2) then 
  kplane = 1 
  folder=trim(path)//'output'
  !coord==18 is case specific, this needs to be changed depending on where and what height 
  !the user wants the snapshot at 
  if(coord==18) call write_2Dplane(2,u(1:nx,1:ny,kplane),'plane_u',folder)
  if(coord==18) call write_2Dplane(2,v(1:nx,1:ny,kplane),'plane_v',folder)
endif

end subroutine inst_write


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Checkpoint routine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine checkpoint (version)
use param, only : rp,dt,jt_total,coord,path,total_time,tavg_nstart
use param, only : nx,ny,nz,kmax,sgs_model,tavg_calc
use sim_param, only : u,v,w,RHSx,RHSy,RHSz
use sgs_param, only : F_LM,F_MM,F_QN,F_NN,Cs_opt2,lagran_dt
#ifdef CORIOLIS
use param, only: adjust_wind_angle
use sim_param, only : alpha_wind,time1,error_pre,error_int,omega_eff
#endif
#ifdef LVLSET
use sim_param, only : dpdx, dpdy, dpdz
#endif LVLSET
#ifdef SCALAR
use param, only : temp_sink
use sim_param, only: psi_m
#ifdef LVLSET
use scalars_param, only: psi_h
#endif LVLSET
use sgs_param, only: tzTh
use scalars_param, only: theta,RHS_T
use scalars_param, only: theta_const
use sgs_param, only: F_KX,F_XX,F_PY,F_YY,Ds_opt2
#endif LVLSET
use io, only : tavg,tavg_total_time,tavg_initialized
#ifdef SCALAR
use io, only : tavg_ustar
use scalars_param, only: T_s
#endif SCALAR
implicit none
integer,intent(IN) :: version
character (64) :: folder
character*8 :: ipfi

! Store timestep number
if (version.eq.2) then
  82 format(i8.8)
  write(ipfi,82) jt_total
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Store continuation files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (coord == 0) then
  if (version.eq.1) then
  open(1,file=trim(path)//'continuation/total_time.dat',status='unknown',form='formatted')
  elseif (version.eq.2) then
  open(1,file=trim(path)//'continuation/total_time_'//trim(ipfi)//'.dat',status='unknown',form='formatted')
  endif
  write(1,*) jt_total,total_time,dt,tavg_nstart,lagran_dt
  close(1)
end if

#ifdef CORIOLIS
if (adjust_wind_angle) then
 if (coord == 0) then
  if (version.eq.1) then
   open(26,file=trim(path)//'continuation/alpha_wind.dat',status='unknown',form='formatted')
  elseif (version.eq.2) then
   open(26,file=trim(path)//'continuation/alpha_wind_'//trim(ipfi)//'.dat',status='unknown',form='formatted')
  endif
  write(26,*) time1,alpha_wind,omega_eff,error_pre,error_int
  close(26)
 end if
end if
#endif CORIOLIS

folder=trim(path)//'continuation'
call write_3Dfield(version,u      (1:nx,1:ny,1:kmax),'continua_u'    ,folder,1,kmax)
call write_3Dfield(version,v      (1:nx,1:ny,1:kmax),'continua_v'    ,folder,1,kmax)
call write_3Dfield(version,w      (1:nx,1:ny,1:kmax),'continua_w'    ,folder,1,kmax)
call write_3Dfield(version,RHSx   (1:nx,1:ny,1:nz-1),'continua_RHSx' ,folder,0,nz-1)
call write_3Dfield(version,RHSy   (1:nx,1:ny,1:nz-1),'continua_RHSy' ,folder,0,nz-1)
call write_3Dfield(version,RHSz   (1:nx,1:ny,1:nz-1),'continua_RHSz' ,folder,0,nz-1)
#ifdef LVLSET
call write_3Dfield(version,dpdx   (1:nx,1:ny,1:nz-1),'continua_dpdx' ,folder,0,nz-1)
call write_3Dfield(version,dpdy   (1:nx,1:ny,1:nz-1),'continua_dpdy' ,folder,0,nz-1)
call write_3Dfield(version,dpdz   (1:nx,1:ny,1:nz-1),'continua_dpdz' ,folder,0,nz-1)
#endif 
#ifdef SCALAR 
call write_3Dfield(version,theta  (1:nx,1:ny,1:kmax),'continua_theta',folder,1,kmax)
call write_3Dfield(version,RHS_T  (1:nx,1:ny,1:nz-1),'continua_RHSt' ,folder,0,nz-1)
#ifdef LVLSET
call write_3Dfield(version,T_s    (1:nx,1:ny,1:nz-1),'continua_Ts'   ,folder,1,nz-1)
call write_3Dfield(version,psi_m  (1:nx,1:ny,1:nz-1),'continua_psim' ,folder,1,nz-1)
call write_3Dfield(version,psi_h  (1:nx,1:ny,1:nz-1),'continua_psih' ,folder,1,nz-1)
#endif LVLSET
if(coord.eq.0) then
call write_2Dplane(version,tzTh(1:nx,1:ny,1)     ,'continua_sgst3',folder)
#ifndef LVLSET
call write_2Dplane(version,T_s    (1:nx,1:ny)       ,'continua_Ts'   ,folder)
call write_2Dplane(version,psi_m  (1:nx,1:ny)       ,'continua_psim' ,folder)
#endif LVLSET
endif
if(temp_sink) then
call write_1Dfield(version,theta_const(1:nz-1)      ,'theta_const'   ,folder,0,nz-1)
endif
#endif 

if(sgs_model.ne.1) then 
call write_3Dfield(version,Cs_opt2(1:nx,1:ny,1:kmax),'continua_csopt2',folder,1,kmax)
#ifdef SCALAR 
call write_3Dfield(version,Ds_opt2(1:nx,1:ny,1:kmax),'continua_dsopt2',folder,1,kmax)
#endif SCALAR 
endif 

if(sgs_model.eq.5) then
call write_3Dfield(version,F_LM   (1:nx,1:ny,1:kmax),'continua_F_LM'  ,folder,1,kmax)
call write_3Dfield(version,F_MM   (1:nx,1:ny,1:kmax),'continua_F_MM'  ,folder,1,kmax)
call write_3Dfield(version,F_QN   (1:nx,1:ny,1:kmax),'continua_F_QN'  ,folder,1,kmax)
call write_3Dfield(version,F_NN   (1:nx,1:ny,1:kmax),'continua_F_NN'  ,folder,1,kmax)
#ifdef SCALAR
call write_3Dfield(version,F_KX   (1:nx,1:ny,1:kmax),'continua_F_KX'  ,folder,1,kmax)
call write_3Dfield(version,F_XX   (1:nx,1:ny,1:kmax),'continua_F_XX'  ,folder,1,kmax)
call write_3Dfield(version,F_PY   (1:nx,1:ny,1:kmax),'continua_F_PY'  ,folder,1,kmax)
call write_3Dfield(version,F_YY   (1:nx,1:ny,1:kmax),'continua_F_YY'  ,folder,1,kmax)
#endif SCALAR 
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Make backup of the statistical data files 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if( tavg_calc .and. tavg_initialized ) then

folder=trim(path)//'output'

if (coord.eq.0) then
if (version.eq.1) then
open(72,file=trim(folder)//'/tavg_total_time.dat',status='unknown',form='formatted')
elseif (version.eq.2) then
open(72,file=trim(folder)//'/tavg_total_time_'//trim(ipfi)//'.dat',status='unknown',form='formatted')
endif
write (72,*) tavg_total_time
close(72)
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Write 3D statisatical data (when selected)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef TAVG_SP 
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%pres ) ,'tavg_pres'  ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%u )    ,'tavg_u'     ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%u2)    ,'tavg_u2'    ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%v )    ,'tavg_v'     ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%v2)    ,'tavg_v2'    ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%w )    ,'tavg_w'     ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%w2)    ,'tavg_w2'    ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%uw)    ,'tavg_uw'    ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%vw)    ,'tavg_vw'    ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%uv)    ,'tavg_uv'    ,folder,0,nz-1)

#ifdef STATSTAU
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%txx)   ,'tavg_txx'   ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%txy)   ,'tavg_txy'   ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%txz)   ,'tavg_txz'   ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%tyy)   ,'tavg_tyy'   ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%tyz)   ,'tavg_tyz'   ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%tzz)   ,'tavg_tzz'   ,folder,0,nz-1)
#endif

#ifdef TURBINES 
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%fu) ,'tavg_fu' ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%fv) ,'tavg_fv' ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%fx) ,'tavg_fx' ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%fy) ,'tavg_fy' ,folder,0,nz-1)
#endif TURBINES

#ifdef ATM
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%fua) ,'tavg_fua' ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%fva) ,'tavg_fva' ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%fwa) ,'tavg_fwa' ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%fxa) ,'tavg_fxa' ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%fya) ,'tavg_fya' ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%fza) ,'tavg_fza' ,folder,0,nz-1)
#endif ATM 

#ifdef STATSKE
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%u3  )  ,'tavg_u3'    ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%v2u )  ,'tavg_v2u'   ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%w2u )  ,'tavg_w2u'   ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%u2v )  ,'tavg_u2v'   ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%v3  )  ,'tavg_v3'    ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%w2v )  ,'tavg_w2v'   ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%u2w )  ,'tavg_u2w'   ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%v2w )  ,'tavg_v2w'   ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%w3  )  ,'tavg_w3'    ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%s11 )  ,'tavg_s11'    ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%s12 )  ,'tavg_s12'    ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%s13 )  ,'tavg_s13'    ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%s22 )  ,'tavg_s22'    ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%s23 )  ,'tavg_s23'    ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%s33 )  ,'tavg_s33'    ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%pu  )  ,'tavg_pu'    ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%pv  )  ,'tavg_pv'    ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%pw  )  ,'tavg_pw'    ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%txxs11  )  ,'tavg_txxs11'    ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%txys12  )  ,'tavg_txys12'    ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%txzs13  )  ,'tavg_txzs13'    ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%tyys22  )  ,'tavg_tyys22'    ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%tyzs23  )  ,'tavg_tyzs23'    ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%tzzs33  )  ,'tavg_tzzs33'    ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%utxx )   ,'tavg_utxx'    ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%vtyx )  ,'tavg_vtyx'    ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%wtzx )   ,'tavg_wtzx'    ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%utxy )   ,'tavg_utxy'    ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%vtyy )   ,'tavg_vtyy'    ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%wtzy )   ,'tavg_wtzy'    ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%utxz )   ,'tavg_utxz'    ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%vtyz )   ,'tavg_vtyz'    ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%wtzz )   ,'tavg_wtzz'    ,folder,0,nz-1)
#endif 

#ifdef SCALAR
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%csopt2),'tavg_csopt2',folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%theta ),'tavg_theta' ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%T2    ),'tavg_T2'    ,folder,0,nz-1)
!sgst3 is printed in double precision 
call write_3Dfield(version,(tavg(1:nx,1:ny,1:nz-1)%sgst3 ),'tavg_sgst3' ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%dsopt2),'tavg_dsopt2',folder,0,nz-1)
call write_3Dfield(version,(tavg(1:nx,1:ny,1:nz-1)%wT    ),'tavg_wT'    ,folder,0,nz-1)
call write_3Dfield_SinglePrecision(version,real(tavg(1:nx,1:ny,1:nz-1)%buoy  ),'tavg_buoy'  ,folder,0,nz-1)
call write_1Dfield(version,(tavg_ustar(1:nz-1)           ),'tavg_ustar' ,folder,0,nz-1)
#endif

#else 
!else write everything in double precision 
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%pres  ,'tavg_pres'  ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%u     ,'tavg_u'     ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%u2    ,'tavg_u2'    ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%v     ,'tavg_v'     ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%v2    ,'tavg_v2'    ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%w     ,'tavg_w'     ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%w2    ,'tavg_w2'    ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%uw    ,'tavg_uw'    ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%vw    ,'tavg_vw'    ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%uv    ,'tavg_uv'    ,folder,0,nz-1)

#ifdef STATSTAU
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%txx   ,'tavg_txx'   ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%txy   ,'tavg_txy'   ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%txz   ,'tavg_txz'   ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%tyy   ,'tavg_tyy'   ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%tyz   ,'tavg_tyz'   ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%tzz   ,'tavg_tzz'   ,folder,0,nz-1)
#endif

#ifdef TURBINES 
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%fu ,'tavg_fu' ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%fv ,'tavg_fv' ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%fx ,'tavg_fx' ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%fy ,'tavg_fy' ,folder,0,nz-1)
#endif TURBINES 

#ifdef ATM
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%fua ,'tavg_fua' ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%fva ,'tavg_fva' ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%fwa ,'tavg_fwa' ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%fxa ,'tavg_fxa' ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%fya ,'tavg_fya' ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%fza ,'tavg_fza' ,folder,0,nz-1)
#endif ATM

#ifdef STATSKE
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%u3    ,'tavg_u3'    ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%v2u   ,'tavg_v2u'   ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%w2u   ,'tavg_w2u'   ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%u2v   ,'tavg_u2v'   ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%v3    ,'tavg_v3'    ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%w2v   ,'tavg_w2v'   ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%u2w   ,'tavg_u2w'   ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%v2w   ,'tavg_v2w'   ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%w3    ,'tavg_w3'    ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%s11   ,'tavg_s11'    ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%s12   ,'tavg_s12'    ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%s13   ,'tavg_s13'    ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%s22   ,'tavg_s22'    ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%s23   ,'tavg_s23'    ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%s33   ,'tavg_s33'    ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%pu    ,'tavg_pu'    ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%pv    ,'tavg_pv'    ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%pw    ,'tavg_pw'    ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%txxs11    ,'tavg_txxs11'    ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%txys12    ,'tavg_txys12'    ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%txzs13    ,'tavg_txzs13'    ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%tyys22    ,'tavg_tyys22'    ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%tyzs23    ,'tavg_tyzs23'    ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%tzzs33    ,'tavg_tzzs33'    ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%utxx    ,'tavg_utxx'    ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%vtyx   ,'tavg_vtyx'    ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%wtzx    ,'tavg_wtzx'    ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%utxy    ,'tavg_utxy'    ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%vtyy    ,'tavg_vtyy'    ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%wtzy    ,'tavg_wtzy'    ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%utxz    ,'tavg_utxz'    ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%vtyz    ,'tavg_vtyz'    ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%wtzz    ,'tavg_wtzz'    ,folder,0,nz-1)
#endif 

#ifdef SCALAR
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%csopt2,'tavg_csopt2',folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%theta ,'tavg_theta' ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%T2    ,'tavg_T2'    ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%sgst3 ,'tavg_sgst3' ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%dsopt2,'tavg_dsopt2',folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%wT    ,'tavg_wT'    ,folder,0,nz-1)
call write_3Dfield(version,tavg(1:nx,1:ny,1:nz-1)%buoy  ,'tavg_buoy'  ,folder,0,nz-1)
call write_1Dfield(version,tavg_ustar(1:nz-1)           ,'tavg_ustar' ,folder,0,nz-1)
#endif

#endif TAVG_SP  

endif

end subroutine checkpoint


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize the statistical data variable to zero or read from file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine tavg_init
use param, only:nx,ny,nz,inxny,mpi_comm_world
use sim_param, only: u
use param, only:rp,mpi_rp,ierr,localComm,lprec,path,coord,tavg_calc
use io, only : tavg,tavg_total_time,tavg_dt,tavg_initialized
#ifdef SCALAR
use io, only: tavg_ustar
#endif
implicit none
character (64) :: folder
logical :: exst
real(rp) :: tavg_u_val, u_val
integer :: k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Allocated the tavg structure
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if( tavg_calc ) then
  allocate(tavg(nx,ny,1:nz-1))
#ifdef SCALAR
  allocate(tavg_ustar(1:nz-1)) 
#endif 
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize tavg structure
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

folder=trim(path)//'output'

!Initialize time averaging structures 
  tavg_dt           =0.0_lprec
  tavg_total_time   =0.0_lprec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize 3D arrays (when requested)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  tavg(:,:,:)%pres  =0.0_rp
  tavg(:,:,:)%u     =0.0_rp
  tavg(:,:,:)%u2    =0.0_rp
  tavg(:,:,:)%v     =0.0_rp
  tavg(:,:,:)%v2    =0.0_rp
  tavg(:,:,:)%w     =0.0_rp
  tavg(:,:,:)%w2    =0.0_rp
  tavg(:,:,:)%uv    =0.0_rp
  tavg(:,:,:)%uw    =0.0_rp
  tavg(:,:,:)%vw    =0.0_rp

#ifdef STATSTAU
  tavg(:,:,:)%txx   =0.0_rp
  tavg(:,:,:)%txy   =0.0_rp
  tavg(:,:,:)%tyy   =0.0_rp
  tavg(:,:,:)%txz   =0.0_rp
  tavg(:,:,:)%tyz   =0.0_rp
  tavg(:,:,:)%tzz   =0.0_rp
#endif

#ifdef TURBINES 
  tavg(:,:,:)%fx   =0.0_rp
  tavg(:,:,:)%fy   =0.0_rp
  tavg(:,:,:)%fu   =0.0_rp
  tavg(:,:,:)%fv   =0.0_rp
#endif TURBINES 

#ifdef ATM
  tavg(:,:,:)%fxa   =0.0_rp
  tavg(:,:,:)%fya   =0.0_rp
  tavg(:,:,:)%fza   =0.0_rp
  tavg(:,:,:)%fua   =0.0_rp
  tavg(:,:,:)%fva   =0.0_rp
  tavg(:,:,:)%fwa   =0.0_rp
#endif ATM 

#ifdef STATSKE
  tavg(:,:,:)%u3    =0.0_rp
  tavg(:,:,:)%v2u   =0.0_rp
  tavg(:,:,:)%w2u   =0.0_rp
  tavg(:,:,:)%u2v   =0.0_rp
  tavg(:,:,:)%v3    =0.0_rp
  tavg(:,:,:)%w2v   =0.0_rp
  tavg(:,:,:)%u2w   =0.0_rp
  tavg(:,:,:)%v2w   =0.0_rp
  tavg(:,:,:)%w3    =0.0_rp
  tavg(:,:,:)%pu = 0._rp
  tavg(:,:,:)%pv = 0._rp
  tavg(:,:,:)%pw = 0._rp
  tavg(:,:,:)%txxs11 = 0._rp
  tavg(:,:,:)%txys12 = 0._rp
  tavg(:,:,:)%txzs13 = 0._rp
  tavg(:,:,:)%tyys22 = 0._rp
  tavg(:,:,:)%tyzs23 = 0._rp
  tavg(:,:,:)%tzzs33 = 0._rp
  tavg(:,:,:)%s11 = 0._rp
  tavg(:,:,:)%s12 = 0._rp
  tavg(:,:,:)%s13 = 0._rp
  tavg(:,:,:)%s22 = 0._rp
  tavg(:,:,:)%s23 = 0._rp
  tavg(:,:,:)%s33 = 0._rp
  tavg(:,:,:)%utxx = 0._rp
  tavg(:,:,:)%vtyx = 0._rp
  tavg(:,:,:)%wtzx = 0._rp
  tavg(:,:,:)%utxy = 0._rp
  tavg(:,:,:)%vtyy = 0._rp
  tavg(:,:,:)%wtzy = 0._rp
  tavg(:,:,:)%utxz = 0._rp
  tavg(:,:,:)%vtyz = 0._rp
  tavg(:,:,:)%wtzz = 0._rp
#endif 

#ifdef SCALAR
  tavg(:,:,:)%csopt2=0.0_rp
  tavg(:,:,:)%theta =0.0_rp
  tavg(:,:,:)%T2    =0.0_rp
  tavg(:,:,:)%sgst3 =0.0_rp
  tavg(:,:,:)%dsopt2=0.0_rp
  tavg(:,:,:)%buoy  =0.0_rp
  tavg(:,:,:)%wT    =0.0_rp
  tavg_ustar(:)     =0.0_rp
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Check whether information should be read from files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

inquire (file=trim(folder)//'/tavg_total_time.dat',exist=exst)
if (.not. exst) then ! Initialize when no previous data are present
  if (coord == 0) then
    write(*,*)'No previous time averaged data - starting from scratch.'
  endif
else ! Read data from file
  if (coord.eq.0) then
  open(73,file=trim(folder)//'/tavg_total_time.dat',status='unknown',form='formatted')
  read (73,*) tavg_total_time
  write(*,*) 'Read in time averaged data from previous run'
  close(73)
  endif
  !Send information from master to other tasks
  call mpi_bcast(tavg_total_time,1,mpi_rp,0,localComm,ierr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read 3D data *(when requested)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%pres    ,'tavg_pres'     ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%u     ,'tavg_u'     ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%u2    ,'tavg_u2'    ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%v     ,'tavg_v'     ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%v2    ,'tavg_v2'    ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%w     ,'tavg_w'     ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%w2    ,'tavg_w2'    ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%uw    ,'tavg_uw'    ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%vw    ,'tavg_vw'    ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%uv    ,'tavg_uv'    ,folder,0,nz-1)

#ifdef STATSTAU
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%txx   ,'tavg_txx'   ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%txy   ,'tavg_txy'   ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%txz   ,'tavg_txz'   ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%tyy   ,'tavg_tyy'   ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%tyz   ,'tavg_tyz'   ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%tzz   ,'tavg_tzz'   ,folder,0,nz-1)
#endif


#ifdef TURBINES 
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%fu ,'tavg_fu' ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%fv ,'tavg_fv' ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%fx ,'tavg_fx' ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%fy ,'tavg_fy' ,folder,0,nz-1)
#endif TURBINES 

#ifdef ATM 
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%fua ,'tavg_fua' ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%fva ,'tavg_fva' ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%fwa ,'tavg_fwa' ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%fxa ,'tavg_fxa' ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%fya ,'tavg_fya' ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%fza ,'tavg_fza' ,folder,0,nz-1)
#endif ATM

#ifdef STATSKE
  
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%u3    ,'tavg_u3'    ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%v2u   ,'tavg_v2u'   ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%w2u   ,'tavg_w2u'   ,folder,0,nz-1)
  
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%u2v   ,'tavg_u2v'   ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%v3    ,'tavg_v3'    ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%w2v   ,'tavg_w2v'   ,folder,0,nz-1)
  
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%u2w   ,'tavg_u2w'   ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%v2w   ,'tavg_v2w'   ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%w3    ,'tavg_w3'    ,folder,0,nz-1)

  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%pu    ,'tavg_pu'    ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%pv    ,'tavg_pv'    ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%pw    ,'tavg_pw'    ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%txxs11 ,'tavg_txxs11' ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%txys12 ,'tavg_txys12' ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%txzs13 ,'tavg_txzs13' ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%tyys22 ,'tavg_tyys22' ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%tyzs23 ,'tavg_tyzs23' ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%tzzs33 ,'tavg_tzzs33' ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%s11 ,'tavg_s11' ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%s12 ,'tavg_s12' ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%s13 ,'tavg_s13' ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%s22 ,'tavg_s22' ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%s23 ,'tavg_s23' ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%s33 ,'tavg_s33' ,folder,0,nz-1)
  
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%utxx ,'tavg_utxx' ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%vtyx ,'tavg_vtyx' ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%wtzx ,'tavg_wtzx' ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%utxy ,'tavg_utxy' ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%vtyy ,'tavg_vtyy' ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%wtzy ,'tavg_wtzy' ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%utxz ,'tavg_utxz' ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%vtyz ,'tavg_vtyz' ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%wtzz ,'tavg_wtzz' ,folder,0,nz-1)
#endif
 
#ifdef SCALAR
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%csopt2,'tavg_csopt2',folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%theta ,'tavg_theta' ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%T2    ,'tavg_T2'    ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%sgst3 ,'tavg_sgst3' ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%dsopt2,'tavg_dsopt2',folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%buoy    ,'tavg_buoy'    ,folder,0,nz-1)
  call read_3Dfield(tavg(1:nx,1:ny,1:nz-1)%wT    ,'tavg_wT'    ,folder,0,nz-1)
  call read_1Dfield(tavg_ustar(1:nz-1)           ,'tavg_ustar' ,folder,0,nz-1)
#endif 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Check if the values in time averaged files are messed up
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do k = 1,nz-1
  tavg_u_val = sum(tavg(1:nx,1:ny,k)%u)*inxny
  u_val = sum(u(1:nx,1:ny,k))*inxny
 
  !If the tavg values are within certian limits the simulation continues 
  !Else the simulation is terminated 
  if(abs(tavg_u_val/u_val-1.0) .gt. 0.25_rp) then
     if(coord == 0 ) write(*,*) 'Something has messed up the time averaged files'
     call mpi_barrier(mpi_comm_world,ierr )
     call mpi_finalize (ierr)
     stop
  endif
enddo

endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize time counters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

tavg_dt=0.0_rp ! Initialize tavg_dt
tavg_initialized=.true. ! Set global switch that tavg has been initialized

end subroutine tavg_init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  This subroutine calculates basic statistical flow quantities 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine tavg_compute
use param, only : lprec,rp,ld,nx,ny,nz
use sim_param,only: pres,u,v,w,u_w=>dummy1,v_w=>dummy2
use io, only: tavg,tavg_total_time,tavg_dt
use grid_defs, only: gridzw
#ifdef TURBINES 
use sim_param, only: fx, fy
#endif TURBINES 
#ifdef ATM 
use sim_param, only: fxa, fya,fza
#endif ATM
#ifdef SCALAR
use io, only: tavg_ustar
use param, only: inxny
use sgs_param, only: Cs_opt2,Ds_opt2
use sim_param,only: theta_w=>dummy3,ustar
use scalars_param, only: theta, beta_scal
use sgs_param, only: tzTh 
#endif
implicit none
real(lprec) :: gg1,gg2
integer :: i,j,k
real(rp) :: dz1, dz2, dzfact 
real(rp), dimension(ld,ny) :: pr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interpolate u and v to the w grid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do k=1,nz-1
dz2 = 0.5_rp * (gridzw(k+1) - gridzw(k)) 
dz1 = 0.5_rp * (gridzw(k) - gridzw(k-1)) 
dzfact = (1._rp/(dz1+dz2))
u_w(1:ld,1:ny,k)=dzfact*(u(1:ld,1:ny,k) * dz1 + u(1:ld,1:ny,k-1) * dz2)
v_w(1:ld,1:ny,k)=dzfact*(v(1:ld,1:ny,k) * dz1 + v(1:ld,1:ny,k-1) * dz2)
#ifdef SCALAR
theta_w(1:ld,1:ny,k) = dzfact* (theta(1:ld,1:ny,k) * dz1 + theta(1:ld,1:ny,k-1) * dz2)
#endif SCALAR 
enddo 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set time counters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

gg1=tavg_total_time/(tavg_total_time+tavg_dt)
gg2=tavg_dt        /(tavg_total_time+tavg_dt)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate statistics in 3D array (when requested)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do k=1,nz-1
 dz2 = 0.5_rp * (gridzw(k+1) - gridzw(k))
 dz1 = 0.5_rp * (gridzw(k) - gridzw(k-1))
 dzfact = (1._rp/(dz1+dz2))
 pr(1:nx,1:ny) = (dzfact * (pres(1:nx,1:ny,k) * dz1 +pres(1:nx,1:ny,k-1) * dz2) - 0.5_rp * (u_w(1:nx,1:ny,k)**2 + v_w(1:nx,1:ny,k)**2 + w(1:nx,1:ny,k)**2))

 do j=1,ny
  do i=1,nx
   tavg(i,j,k)%pres = gg1 *  tavg(i,j,k)%pres +  gg2 * pr(i,j)
   tavg(i,j,k)%u     = gg1 * tavg(i,j,k)%u      +gg2 * u(i,j,k)
   tavg(i,j,k)%u2    = gg1 * tavg(i,j,k)%u2     +gg2 * u(i,j,k)**2
   tavg(i,j,k)%v     = gg1 * tavg(i,j,k)%v      +gg2 * v(i,j,k)
   tavg(i,j,k)%v2    = gg1 * tavg(i,j,k)%v2     +gg2 * v(i,j,k)**2
   tavg(i,j,k)%w     = gg1 * tavg(i,j,k)%w      +gg2 * w(i,j,k)
   tavg(i,j,k)%w2    = gg1 * tavg(i,j,k)%w2     +gg2 * w(i,j,k)**2
   
   tavg(i,j,k)%uv    = gg1 * tavg(i,j,k)%uv     +gg2 *   u(i,j,k)*v(i,j,k)
   tavg(i,j,k)%uw    = gg1 * tavg(i,j,k)%uw     +gg2 * u_w(i,j,k)*w(i,j,k)
   tavg(i,j,k)%vw    = gg1 * tavg(i,j,k)%vw     +gg2 * v_w(i,j,k)*w(i,j,k)

#ifdef STATSKE
   tavg(i,j,k)%pu = gg1 *  tavg(i,j,k)%pu +  gg2 * (u_w(i,j,k) * pr(i,j))
   tavg(i,j,k)%pv = gg1 *  tavg(i,j,k)%pv +  gg2 * (v_w(i,j,k) * pr(i,j))
   tavg(i,j,k)%pw = gg1 *  tavg(i,j,k)%pw + gg2 * (w(i,j,k) * pr(i,j))
   tavg(i,j,k)%u3    = gg1 * tavg(i,j,k)%u3     +gg2 *   u(i,j,k) *   u(i,j,k) *   u(i,j,k)
   tavg(i,j,k)%v2u   = gg1 * tavg(i,j,k)%v2u    +gg2 *   v(i,j,k) *   v(i,j,k) *   u(i,j,k)
   tavg(i,j,k)%w2u   = gg1 * tavg(i,j,k)%w2u    +gg2 *   w(i,j,k) *   w(i,j,k) * u_w(i,j,k)
   tavg(i,j,k)%u2v   = gg1 * tavg(i,j,k)%u2v    +gg2 *   u(i,j,k) *   u(i,j,k) *   v(i,j,k)
   tavg(i,j,k)%v3    = gg1 * tavg(i,j,k)%v3     +gg2 *   v(i,j,k) *   v(i,j,k) *   v(i,j,k)
   tavg(i,j,k)%w2v   = gg1 * tavg(i,j,k)%w2v    +gg2 *   w(i,j,k) *   w(i,j,k) * v_w(i,j,k)
   tavg(i,j,k)%u2w   = gg1 * tavg(i,j,k)%u2w    +gg2 * u_w(i,j,k) * u_w(i,j,k) *   w(i,j,k)
   tavg(i,j,k)%v2w   = gg1 * tavg(i,j,k)%v2w    +gg2 * v_w(i,j,k) * v_w(i,j,k) *   w(i,j,k)
   tavg(i,j,k)%w3    = gg1 * tavg(i,j,k)%w3     +gg2 *   w(i,j,k) *   w(i,j,k) *   w(i,j,k)
#endif 

#ifdef TURBINES 
   tavg(i,j,k)%fu = gg1 * tavg(i,j,k)%fu + gg2 * fx(i,j,k) * u(i,j,k)
   tavg(i,j,k)%fv = gg1 * tavg(i,j,k)%fv + gg2 * fy(i,j,k) * v(i,j,k)
   tavg(i,j,k)%fx = gg1 * tavg(i,j,k)%fx + gg2 * fx(i,j,k)
   tavg(i,j,k)%fy = gg1 * tavg(i,j,k)%fy + gg2 * fy(i,j,k)
#endif TURBINES 

#ifdef ATM 
   tavg(i,j,k)%fua = gg1 * tavg(i,j,k)%fua + gg2 * fxa(i,j,k) * u(i,j,k)
   tavg(i,j,k)%fva = gg1 * tavg(i,j,k)%fva + gg2 * fya(i,j,k) * v(i,j,k)
   tavg(i,j,k)%fza = gg1 * tavg(i,j,k)%fza + gg2 * fza(i,j,k) * w(i,j,k)
   tavg(i,j,k)%fxa = gg1 * tavg(i,j,k)%fxa + gg2 * fxa(i,j,k)
   tavg(i,j,k)%fya = gg1 * tavg(i,j,k)%fya + gg2 * fya(i,j,k)
   tavg(i,j,k)%fza = gg1 * tavg(i,j,k)%fza + gg2 * fza(i,j,k)
#endif ATM

#ifdef SCALAR
   tavg(i,j,k)%csopt2= gg1 * tavg(i,j,k)%csopt2 +gg2 * cs_opt2(i,j,k)
   tavg(i,j,k)%theta = gg1 * tavg(i,j,k)%theta  +gg2 *   theta(i,j,k)
   tavg(i,j,k)%t2    = gg1 * tavg(i,j,k)%t2     +gg2 *   theta(i,j,k) *   theta(i,j,k)
   tavg(i,j,k)%sgst3 = gg1 * tavg(i,j,k)%sgst3  +gg2 *  (-tzTh(i,j,k))
   tavg(i,j,k)%dsopt2= gg1 * tavg(i,j,k)%dsopt2 +gg2 * ds_opt2(i,j,k)
   tavg(i,j,k)%buoy    = gg1 * tavg(i,j,k)%buoy     +gg2 *   (beta_scal(i,j,k) * w(i,j,k)) 
   tavg(i,j,k)%wT    = gg1 * tavg(i,j,k)%wT     +gg2 *       w(i,j,k) * theta_w(i,j,k)
#endif 
  enddo
 enddo
#ifdef SCALAR
 tavg_ustar(k)       = gg1 * tavg_ustar(k)      +gg2 * inxny * sum(ustar(1:nx,1:ny))
#endif 
enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update time counters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

tavg_total_time= tavg_total_time+tavg_dt ! Update total averaging time
tavg_dt=0.0_rp ! To determine time interval with next call to "tavg_compute"

end subroutine tavg_compute

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  This subroutine calculates basic statistical flow quantities 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine tavg_compute_tau
use param,only: tavg_calc,ld,tavg_nstart,tavg_start_time,tavg_end_time,total_time
use param,only: rp,jt_total,dt,coord
use io, only: tavg,tavg_total_time,tavg_dt,tavg_initialized
#ifdef STATSKE 
use sim_param, only : dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
use sim_param,only: u,v,w,u_w=>dummy1,v_w=>dummy2
use sim_param,only: theta_w=>dummy3,ustar
#ifdef SCALAR
use scalars_param, only: theta
#endif SCALAR
use grid_defs, only: gridzw, gridz 
#endif STATSKE

#ifdef STATSTAU
use param,only: tavg_nskip,nx,ny,nz
use sim_param,only: txx,txy,tyy,txz,tyz,tzz
#endif
#ifdef STATSTAU 
use io, only: tavg
#endif
implicit none
real(rp) :: gg1,gg2
#ifdef STATSTAU
integer :: i,j,k
#endif
real(rp) :: dz1, dz2, dzfact

#ifdef STATSKE 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interpolate u and v to the w grid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do k=1,nz-1
dz2 = 0.5_rp * (gridzw(k+1) - gridzw(k))
dz1 = 0.5_rp * (gridzw(k) - gridzw(k-1))
dzfact = (1._rp/(dz1+dz2))
u_w(1:ld,1:ny,k)=dzfact*(u(1:ld,1:ny,k) * dz1 + u(1:ld,1:ny,k-1) * dz2)
v_w(1:ld,1:ny,k)=dzfact*(v(1:ld,1:ny,k) * dz1 + v(1:ld,1:ny,k-1) * dz2)
enddo
#endif STATSKE  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialization of the time averaged statistics
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (tavg_calc) then
 if (total_time >= tavg_start_time .and. total_time <= tavg_end_time) then
  tavg_dt = tavg_dt + dt
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Set time counters
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  gg1=tavg_total_time/(tavg_total_time+tavg_dt)
  gg2=tavg_dt        /(tavg_total_time+tavg_dt)

  if(total_time >= tavg_start_time .and. .not.tavg_initialized) then
   tavg_nstart = jt_total
    if (coord == 0) then
     write(*,*) '-------------------------------'
     write(*,"(1a,i9,1a,i9)") 'Starting statistical calculations from ', tavg_nstart
     write(*,*) '-------------------------------'
    endif  ! coord==0
   call tavg_init
  endif


if(mod(jt_total-tavg_nstart,tavg_nskip).eq.0.and.jt_total-tavg_nstart.gt.0) then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate statistics in 3D array (when requested)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   do k=1,nz-1
    do j=1,ny
     do i=1,nx
#ifdef STATSTAU
      tavg(i,j,k)%txx   = gg1*tavg(i,j,k)%txx    +gg2* txx(i,j,k)
      tavg(i,j,k)%txy   = gg1*tavg(i,j,k)%txy    +gg2* txy(i,j,k)
      tavg(i,j,k)%tyy   = gg1*tavg(i,j,k)%tyy    +gg2* tyy(i,j,k)
      tavg(i,j,k)%txz   = gg1*tavg(i,j,k)%txz    +gg2* txz(i,j,k)
      tavg(i,j,k)%tyz   = gg1*tavg(i,j,k)%tyz    +gg2* tyz(i,j,k)
      tavg(i,j,k)%tzz   = gg1*tavg(i,j,k)%tzz    +gg2* tzz(i,j,k)
#endif STATSTAU 

#ifdef STATSKE
      tavg(i,j,k)%utxx = gg1 * tavg(i,j,k)%utxx + gg2 * u(i,j,k) * txx(i,j,k)
      tavg(i,j,k)%vtyx = gg1 * tavg(i,j,k)%vtyx + gg2 * v(i,j,k) * txy(i,j,k)
      tavg(i,j,k)%wtzx = gg1 * tavg(i,j,k)%wtzx + gg2 * w(i,j,k) * txz(i,j,k)
      tavg(i,j,k)%utxy = gg1 * tavg(i,j,k)%utxy + gg2 * u(i,j,k) * txy(i,j,k)
      tavg(i,j,k)%vtyy = gg1 * tavg(i,j,k)%vtyy + gg2 * v(i,j,k) * tyy(i,j,k)
      tavg(i,j,k)%wtzy = gg1 * tavg(i,j,k)%wtzy + gg2 * w(i,j,k) * tyz(i,j,k)
      tavg(i,j,k)%utxz = gg1 * tavg(i,j,k)%utxz + gg2 * u_w(i,j,k) * txz(i,j,k)
      tavg(i,j,k)%vtyz = gg1 * tavg(i,j,k)%vtyz + gg2 * v_w(i,j,k) * tyz(i,j,k)
      tavg(i,j,k)%wtzz = gg1 * tavg(i,j,k)%wtzz + gg2 * w(i,j,k) * tzz(i,j,k)
      tavg(i,j,k)%txxs11 = gg1 * tavg(i,j,k)%txxs11 + gg2 * txx(i,j,k) * 0.5_rp * (dudx(i,j,k) + dudx(i,j,k))
      tavg(i,j,k)%txys12 = gg1 * tavg(i,j,k)%txys12 + gg2 * txy(i,j,k) * 0.5_rp * (dudy(i,j,k) + dvdx(i,j,k))
      tavg(i,j,k)%txzs13 = gg1 * tavg(i,j,k)%txzs13 + gg2 * txz(i,j,k) * 0.5_rp * (dudz(i,j,k) + dwdx(i,j,k))
      tavg(i,j,k)%tyys22 = gg1 * tavg(i,j,k)%tyys22 + gg2 * tyy(i,j,k) * 0.5_rp * (dvdy(i,j,k) + dvdy(i,j,k))
      tavg(i,j,k)%tyzs23 = gg1 * tavg(i,j,k)%tyzs23 + gg2 * tyz(i,j,k) * 0.5_rp * (dwdy(i,j,k) + dvdz(i,j,k))
      tavg(i,j,k)%tzzs33 = gg1 * tavg(i,j,k)%tzzs33 + gg2 * tzz(i,j,k) * 0.5_rp * (dwdz(i,j,k) + dwdz(i,j,k))
      tavg(i,j,k)%s11 = gg1 * tavg(i,j,k)%s11 + gg2 * 0.5_rp * (dudx(i,j,k) + dudx(i,j,k))
      tavg(i,j,k)%s12 = gg1 * tavg(i,j,k)%s12 + gg2 * 0.5_rp * (dudy(i,j,k) + dvdx(i,j,k))
      tavg(i,j,k)%s13 = gg1 * tavg(i,j,k)%s13 + gg2 * 0.5_rp * (dudz(i,j,k) + dwdx(i,j,k))
      tavg(i,j,k)%s22 = gg1 * tavg(i,j,k)%s22 + gg2 * 0.5_rp * (dvdy(i,j,k) + dvdy(i,j,k))
      tavg(i,j,k)%s23 = gg1 * tavg(i,j,k)%s23 + gg2 * 0.5_rp * (dwdy(i,j,k) + dvdz(i,j,k))
      tavg(i,j,k)%s33 = gg1 * tavg(i,j,k)%s33 + gg2 * 0.5_rp * (dwdz(i,j,k) + dwdz(i,j,k))
#endif STATSKE 
     enddo
    enddo
   enddo
  endif
 endif  ! after nstart
endif  ! tavg_calc

end subroutine tavg_compute_tau
