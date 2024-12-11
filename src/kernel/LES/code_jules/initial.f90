!*******************************************************************************
subroutine initial
!*******************************************************************************
!
! Driver routine to generate the initial condition
!
use param, only: rp,nx,ny,nz,kmax,jt_total,dt,sgs_model,tavg_nstart,coord
use param, only: use_cfl_dt,initu,path,total_time,initu_interp,cont_reduce
use param, only: total_time,ierr,inxny,nsteps
use mpi, only: mpi_comm_world
use grid_defs, only: gridz
use sim_param, only: u,v,w,RHSx,RHSy,RHSz
use sgs_param, only: F_LM,F_MM,F_QN,F_NN,Cs_opt2,lagran_dt
use mpi_defs, only: mpi_sync_downup,mpi_sync_down
#ifdef CORIOLIS
use param, only: adjust_wind_angle,ug,vg
use sim_param, only: alpha_wind_old,time1,error_pre,error_int,omega_eff,vgMag
#endif CORIOLIS
#ifdef LVLSET
use sim_param, only: dpdx, dpdy, dpdz
#endif LVLSET
#ifdef SCALAR
#ifdef CORIOLIS
use param, only: temp_sink
use scalars_param, only: theta_const
#endif CORIOLIS
use sim_param, only: psi_m
#ifdef LVLSET
use scalars_param, only: psi_h
#endif LVLSET
use sgs_param, only: tzTh 
use scalars_param, only: theta, RHS_T
use sgs_param, only: F_KX, F_XX, F_PY, F_YY, Ds_opt2
use scalars_param, only: T_s
#endif SCALAR
#ifdef BAROCLINIC
use param, only: ug, vg
use param, only: ug_delta, vg_delta, bar_start, bar_end, geo_force_x, geo_force_y
#endif BAROCLINIC
implicit none
real(rp), dimension(1:nz-1) :: theta_temp
character (64) :: folder
real(rp) :: dt_r ! Temporary values used to read time step and CFL from
logical :: exst,exst2
integer :: k, jz

if(initu) then
! Read the initial condition from the files (when available)

! Check whether the velocity files are present (in fact check total_time.dat)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  folder=trim(path)//'continuation'
  inquire(file=trim(folder)//'/total_time.dat',exist=exst)
  if (exst) then
    ! Read the data from the file
    open(1,file=trim(folder)//'/total_time.dat')
    read(1,*) jt_total,total_time,dt_r,tavg_nstart,lagran_dt
    close(1)
    ! Update dynamic time stepping info if required
    if( use_cfl_dt ) then
      dt = dt_r
    endif
    cont_reduce = .true.
    ! Stop the simulation when the number of timesteps is not set correctly
    if( jt_total >= nsteps ) then
      if(coord == 0) write(*,*) 'Full number of time steps reached'
      call mpi_barrier(mpi_comm_world, ierr)
      call mpi_finalize (ierr)
      stop
    endif
  else ! Stop the simulation when the datafiles are not present
    if(coord == 0) write(*,*) 'Data files are not present'
    call mpi_barrier(mpi_comm_world, ierr)
    call mpi_finalize (ierr)
    stop
  end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef CORIOLIS
  if (adjust_wind_angle) then
   call wind_angle_controller_init
   if(coord == 0) write(*,*) '--> Reading PID control data from file'
   inquire(file=trim(folder)//'/alpha_wind.dat',exist=exst)
   if (exst) then
    ! Read the data from the file
    open(26,file=trim(folder)//'/alpha_wind.dat')
    read(26,*) time1,alpha_wind_old,omega_eff,error_pre,error_int
    close(26)
   else ! Stop the simulation when the datafiles are not present
    if(coord == 0) write(*,*) '! Restart-files needed for PID controller not present : Resetting PID parameters!'
    alpha_wind_old=0.0; omega_eff = 0.0; error_pre = 0.0; error_int = 0.0 
    time1 = total_time  
   end if
  endif
#endif CORIOLIS

#ifdef BAROCLINIC
  ! Linear geostrophic shear with height (corresponding to a height and time invariant horizontal temperature gradient)
  ! as previously used by Momen (2018, Journal of Atmospheric Sciences) and Sorbjan (2004, Bound.-Layer Meteor.)
  do jz=1,nz
    if( gridz(jz) <= bar_start) then
      geo_force_x(jz) = ug
      geo_force_y(jz) = vg
    elseif(bar_start < gridz(jz) .and. gridz(jz) <= bar_end) then
      geo_force_x(jz) = ug + ( ug_delta/(bar_end-bar_start) ) * ( gridz(jz)-bar_start )
      geo_force_y(jz) = vg + ( vg_delta/(bar_end-bar_start) ) * ( gridz(jz)-bar_start )
    else
      geo_force_x(jz) = ug + ug_delta
      geo_force_y(jz) = vg + vg_delta
    endif
  enddo
#endif

  if(coord == 0) write(*,*) '--> Reading initial velocity field from file'

  ! Below we first read the information from the HDF5 files and subsequently we
  ! synchronize the information in the ghost layers. One could theoretically
  ! prevent the synchronization step by reading in the ghost layers directly,
  ! however, this leads to several processors that want to read the same data from
  ! the HDF5 files. The current setup prevents this data reading congestion.

  folder=trim(path)//'continuation'
  call read_3Dfield(u   (1:nx,1:ny,1:kmax),'continua_u'   ,folder,1,kmax)
  call read_3Dfield(v   (1:nx,1:ny,1:kmax),'continua_v'   ,folder,1,kmax)
  call read_3Dfield(w   (1:nx,1:ny,1:kmax),'continua_w'   ,folder,1,kmax)
  call read_3Dfield(RHSx(1:nx,1:ny,1:nz-1),'continua_RHSx',folder,0,nz-1)
  call read_3Dfield(RHSy(1:nx,1:ny,1:nz-1),'continua_RHSy',folder,0,nz-1)
  call read_3Dfield(RHSz(1:nx,1:ny,1:nz-1),'continua_RHSz',folder,0,nz-1)
#ifdef LVLSET
  ! Srinidhi : Feb 20, 2022 
  ! Sometimes we might want to run a spinup simulation and thereafter introduce 
  ! complex terrain in the generated initial fields.
  ! In that particular case the grads dpdx, dpdy, dpdz are not available for restart. 
  ! Therefore, this workaround. The files are read iff they are available. If they aren't
  ! the simulation is simply started with the spinup values. 

  folder=trim(path)//'continuation'
  inquire(file=trim(folder)//'/continua_dpdx.h5',exist=exst)

  if(exst) call read_3Dfield(dpdx(1:nx,1:ny,1:nz-1),'continua_dpdx',folder,0,nz-1)
  if(exst) call read_3Dfield(dpdy(1:nx,1:ny,1:nz-1),'continua_dpdy',folder,0,nz-1)
  if(exst) call read_3Dfield(dpdz(1:nx,1:ny,1:nz-1),'continua_dpdz',folder,0,nz-1)
#endif LVLSET

#ifdef SCALAR
  call read_3Dfield(theta (1:nx,1:ny,1:kmax),'continua_theta'   ,folder,1,kmax)
  call read_3Dfield(RHS_T(1:nx,1:ny,1:nz-1) ,'continua_RHSt'    ,folder,0,nz-1)
  if(coord==0) call read_2Dplane(tzTh(1:nx,1:ny,1),'continua_sgst3',folder)
 
#ifdef LVLSET
  if(exst) then 
        call read_3Dfield(psi_m(1:nx,1:ny,1:nz-1), 'continua_psim', folder, 1, nz-1)
        call read_3Dfield(psi_h(1:nx,1:ny,1:nz-1), 'continua_psih', folder, 1, nz-1)
        call read_3Dfield(T_s(1:nx,1:ny,1:nz-1)  , 'continua_Ts'  , folder, 1, nz-1)
  else
        if(coord==0)  call read_2Dplane(psi_m(1:nx,1:ny), 'continua_psim', folder)
        if(coord==0)  call read_2Dplane(T_s(1:nx,1:ny)  , 'continua_Ts'  , folder)
  endif 
#else
  if(coord==0)  call read_2Dplane(psi_m(1:nx,1:ny), 'continua_psim', folder)
  if(coord==0)  call read_2Dplane(T_s(1:nx,1:ny)  , 'continua_Ts'  , folder)
#endif LVLSET
#ifdef CORIOLIS
  if(temp_sink) call read_1Dfield(theta_const(1:nz-1),'theta_const'   ,folder)
#endif CORIOLIS
#endif SCALAR

  call mpi_sync_real_array(u    ,mpi_sync_downup)
  call mpi_sync_real_array(v    ,mpi_sync_downup)
  call mpi_sync_real_array(w    ,mpi_sync_downup)
#ifdef SCALAR 
  call mpi_sync_real_array(theta,mpi_sync_downup)
#endif SCALAR

  if(sgs_model .ne. 1) then  
  ! These variables are not used in the Smagorinsky model
    call read_3Dfield(Cs_opt2(1:nx,1:ny,1:kmax),'continua_csopt2',folder,1,kmax)
    call mpi_sync_real_array(Cs_opt2,mpi_sync_downup) !RS: only down is really required, can be completed prevented by reading correct part from continuation file
#ifdef SCALAR 
    call read_3Dfield(Ds_opt2(1:nx,1:ny,1:kmax),'continua_dsopt2',folder,1,kmax)
    call mpi_sync_real_array(Ds_opt2,mpi_sync_downup) !RS: only down is really required, can be completed prevented by reading correct part from continuation file
#endif SCALAR 
  endif 
  if(sgs_model.eq.5) then
    call read_3Dfield(F_LM   (1:nx,1:ny,1:kmax),'continua_F_LM'  ,folder,1,kmax)
    call read_3Dfield(F_MM   (1:nx,1:ny,1:kmax),'continua_F_MM'  ,folder,1,kmax)
    call read_3Dfield(F_QN   (1:nx,1:ny,1:kmax),'continua_F_QN'  ,folder,1,kmax)
    call read_3Dfield(F_NN   (1:nx,1:ny,1:kmax),'continua_F_NN'  ,folder,1,kmax)
#ifdef SCALAR
    call read_3Dfield(F_KX   (1:nx,1:ny,1:kmax),'continua_F_KX'  ,folder,1,kmax)
    call read_3Dfield(F_XX   (1:nx,1:ny,1:kmax),'continua_F_XX'  ,folder,1,kmax)
    call read_3Dfield(F_PY   (1:nx,1:ny,1:kmax),'continua_F_PY'  ,folder,1,kmax)
    call read_3Dfield(F_YY   (1:nx,1:ny,1:kmax),'continua_F_YY'  ,folder,1,kmax)
#endif SCALAR
    call mpi_sync_real_array(F_LM   ,mpi_sync_downup)
    call mpi_sync_real_array(F_MM   ,mpi_sync_downup)
    call mpi_sync_real_array(F_QN   ,mpi_sync_downup)
    call mpi_sync_real_array(F_NN   ,mpi_sync_downup)
#ifdef SCALAR
    call mpi_sync_real_array(F_KX   ,mpi_sync_downup)
    call mpi_sync_real_array(F_XX   ,mpi_sync_downup)
    call mpi_sync_real_array(F_PY   ,mpi_sync_downup)
    call mpi_sync_real_array(F_YY   ,mpi_sync_downup)
#endif SCALAR

  endif

else !initu else

  folder=trim(path)//'output'
  inquire(file=trim(folder)//'/check_ke.dat',exist=exst)
  inquire(file=trim(folder)//'/tavg_u.h5',exist=exst2)
  if (exst.or.exst2) then
    if(coord == 0) write(*,*) 'Previous simulation data seem to be present in output folder'
    call mpi_barrier(mpi_comm_world, ierr)
    call mpi_finalize (ierr)
    stop
  endif

  if (initu_interp) then
    ! Read the initial velocity (and temperature) condition from the files
    if (coord == 0) write(*,*) '--> Reading initial velocity field from file'
    ! Note no coord.eq.0 so it is printed many times and easier to notice
    write(*,*) '--> Please read the comments in the code'
    ! Note that other variables may also need to be interpolated. See the part of the code starting
    ! at “Reading initial velocity field from file” to see what variables are required for your case
    ! When these variables are not read in it may take longer for your simulation to reach the 
    ! statistical stationary state; when the grid resolution is changed some adjustment period is required
    ! anyhow, and variables (hopefully) adjust automatically
    
    ! Below we first read the information from the HDF5 files and subsequently we
    ! synchronize the information in the ghost layers.
    folder = trim(path) // 'continuation'
    call read_3Dfield (u(1:nx,1:ny,1:kmax), 'continua_u', folder, 1, kmax)
    call read_3Dfield (v(1:nx,1:ny,1:kmax), 'continua_v', folder, 1, kmax)
    call read_3Dfield (w(1:nx,1:ny,1:kmax), 'continua_w', folder, 1, kmax)
    call mpi_sync_real_array (u, mpi_sync_downup)
    call mpi_sync_real_array (v, mpi_sync_downup)
    call mpi_sync_real_array (w, mpi_sync_downup)
#ifdef SCALAR
    call read_3Dfield (theta(1:nx,1:ny,1:kmax), 'continua_theta', folder, 1, kmax)
    call mpi_sync_real_array (theta, mpi_sync_downup)
#endif SCALAR

  else

    ! Generate a new random initial velocity field
#ifdef SCALAR
    call ic_scal
#else
    call ic
#endif SCALAR

    ! Generate a zero velocity (and a uniform temperature) inside the body
#ifdef LVLSET
    call level_set_ic
#endif LVLSET

  end if ! initu_interp end
end if ! initu end

! Display the mean vertical profiles of the initial velocity field
#ifdef SCALAR
7780 format('k,z,ubar,vbar,wbar,thetabar:',(1x,I3,5(1x,F9.4)))
do k=1,nz
  write(6,7780) k,gridz(k),sum(u(1:nx,1:ny,k))*inxny,sum(v(1:nx,1:ny,k))*inxny,sum(w(1:nx,1:ny,k))*inxny,sum(theta(1:nx,1:ny,k))*inxny
end do
do k = 1,nz-1
  theta_temp(k) = sum(theta(1:nx,1:ny,k))*inxny
enddo
call write_1Dfield(1,theta_temp(1:nz-1),'initial_theta',folder,0,nz-1)
do k = 1,nz-1
  theta_temp(k) = sum(u(1:nx,1:ny,k))*inxny
enddo
call write_1Dfield(1,theta_temp(1:nz-1),'initial_u',folder,0,nz-1)
do k = 1,nz-1
  theta_temp(k) = sum(v(1:nx,1:ny,k))*inxny
enddo
call write_1Dfield(1,theta_temp(1:nz-1),'initial_v',folder,0,nz-1)
#else
7780 format('k,z,ubar,vbar,wbar:',(1x,I3,4(1x,F9.4)))
do k=1,nz
  write(6,7780) k,gridz(k),sum(u(1:nx,1:ny,k))*inxny,sum(v(1:nx,1:ny,k))*inxny,sum(w(1:nx,1:ny,k))*inxny
end do
#endif SCALAR

end subroutine initial
