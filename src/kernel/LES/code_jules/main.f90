!*******************************************************************************
program main
!*******************************************************************************
!
! Driver routine for flow solver and initialization
!
use mpi, only: mpi_max, mpi_comm_world, mpi_wtime
use param, only: rp,path,ny,nz,ld,clock_loc,tadv1,tadv2,tadv3,tfac1,tfac2,jt_total,cfl,cont_reduce
use param, only: wbase,cs_count,clock,clock_total,tavg_end_time
use param, only: coord,mpi_rp,ierr,runtime,initu
use param, only: inflow,dt_f,dt_ff,dt,use_cfl_dt,total_time,nsteps
use param, only: time_press,time_convec,time_sgs,time_project,time_deriv
use param, only: time_divstress,time_output,time_pdf,time_spec,time_wall,time_RHS,time_cfl
use sim_param, only: u, v, w
use sim_param, only: RHSx,RHSy,RHSz,RHSx_f,RHSy_f,RHSz_f, RHSx_ff, RHSy_ff, RHSz_ff
use sim_param, only: dudx,dvdx,dwdx,dudy,dvdy,dwdy,dudz,dvdz,dwdz
#ifdef TURBINES 
use param, only: time_turbine
#endif TURBINES
#ifdef ATM
use param, only: time_atm
use sim_param, only: fxa,fya,fza
use atm_lesgo_interface, only: atm_lesgo_finalize
use atm_lesgo_interface, only: atm_lesgo_forcing
#endif ATM
#ifdef CORIOLIS
use param, only: time_coriolis
use param, only: adjust_wind_angle
#endif CORIOLIS
#ifdef SCALAR
use param, only: time_scalar,time_scalar_deriv,time_temperature
use scalars_param, only: theta, dTdx, dTdy, dTdz
#endif SCALAR
#ifdef LVLSET
use param, only: time_level_set,time_level_set2
#endif LVLSET
use hdf5
implicit none
real(rp):: rmsdivvel,ke,dummy2,dummy3,rbuffer,tconst1,tconst2,tconst3
integer :: nstart,hdf_error
integer :: iarg = 0
character(len=32) :: arg

! Initialize time variable
jt_total          =0
time_press        =0._rp
time_convec       =0._rp
time_sgs          =0._rp
time_project      =0._rp
time_deriv        =0._rp
time_divstress    =0._rp
time_output       =0._rp
time_pdf          =0._rp
time_spec         =0._rp
time_wall         =0._rp
time_RHS          =0._rp
time_cfl          =0._rp
#ifdef SCALAR
time_scalar       =0._rp
time_scalar_deriv =0._rp
time_temperature  =0._rp
#endif SCALAR
#ifdef TURBINES
time_turbine      =0._rp
#endif TURBINES
#ifdef LVLSET
time_level_set    =0._rp
time_level_set2   =0._rp
#endif LVLSET
#ifdef ATM
time_atm      =0._rp
#endif ATM

!Get the path from command line
do
  call get_command_argument(iarg, arg)
  if(len_trim(arg) == 0) exit
  if(arg == 'red') then
    PATH =trim('./red/')
  elseif(arg == 'blue') then
    PATH =trim('./blue/')
  endif
  iarg = iarg+1
enddo

! Initialize all data
call initialize

clock(2) = mpi_wtime()
if(coord == 0) write(*,'(1a,E15.7)') 'Initialization wall time: ', clock(2)-clock(1)
clock_total(1) = mpi_wtime()

! Initialize starting loop index 
! If new simulation jt_total=0 by definition, if restarting jt_total
! provided by total_time.dat
if(initu) then 
  nstart = jt_total
else
  nstart = jt_total + 1
endif 

! Begin time loop 
time_loop: do jt_total = nstart, nsteps

  ! Get the starting time for the iteration
  clock(1) = mpi_wtime()
  clock_loc(1) = mpi_wtime()

  if( use_cfl_dt ) then
   dt_ff = dt_f
   dt_f = dt
   call get_cfl_dt(dt)


   ! 3rd order time integration implemented by Jens H. Kasper  
   tadv1 = ((dt*((3._rp*dt_f + 2._rp*dt)/(dt_ff + dt_f) + 3._rp))/dt_f + 6._rp)/6._rp
   tadv2 = -(dt*(3._rp*dt_ff + 3._rp*dt_f + 2._rp*dt))/(6._rp*dt_ff*dt_f)
   tadv3 = (dt*(3._rp*dt_f + 2._rp*dt))/(6._rp*dt_ff*(dt_ff + dt_f))
 
   if (jt_total .eq. 2 .OR. (cont_reduce)) then
     tfac1 = -dt_ff / dt_f
     tfac2 = 1._rp-tfac1
   end if
  else
  
  call get_dt_cfl(dt, cfl)
  
  endif

  ! Advance time
  total_time = total_time + dt

  clock_loc(2) = mpi_wtime()
  time_cfl=time_cfl+clock_loc(2)-clock_loc(1)

  clock_loc(1) = mpi_wtime()

  ! Save previous time's right-hand-sides for Adams-Bashforth Integration
  ! NOTE: RHS does not contain the pressure gradient
  ! NOTE: RHS_f does contain the pressure gradient
  RHSx_ff(:,:,1:nz-1) = RHSx_f(:,:,1:nz-1)
  RHSy_ff(:,:,1:nz-1) = RHSy_f(:,:,1:nz-1)
  RHSz_ff(:,:,1:nz-1) = RHSz_f(:,:,1:nz-1)

  RHSx_f(:,:,1:nz-1) = RHSx(:,:,1:nz-1)
  RHSy_f(:,:,1:nz-1) = RHSy(:,:,1:nz-1)
  RHSz_f(:,:,1:nz-1) = RHSz(:,:,1:nz-1)

  clock_loc(2) = mpi_wtime()
  time_RHS=time_RHS+clock_loc(2)-clock_loc(1)
  
#ifdef SCALAR
  !//////////////////////////////////////////////////////
  !/// OBUKHOV LENGTH                                 ///
  !//////////////////////////////////////////////////////
  !Calculate obukhov length, mean surface flux, friction velocity on the wall

  clock_loc(1) = mpi_wtime()
  if(coord==0) call obukhov
  clock_loc(2) = mpi_wtime()
  time_scalar=time_scalar+clock_loc(2)-clock_loc(1)
#endif SCALAR

#ifdef CORIOLIS
  !//////////////////////////////////////////////////////
  !/// WIND CONTROLLER                                ///
  !//////////////////////////////////////////////////////

  !Activate the wind angle controller
  clock_loc(1) = mpi_wtime()
  if(adjust_wind_angle) call wind_angle_controller 
  clock_loc(2) = mpi_wtime()
  time_coriolis=time_coriolis+clock_loc(2)-clock_loc(1)
#endif CORIOLIS

  !//////////////////////////////////////////////////////
  !/// DERIVATIVES                                    ///
  !//////////////////////////////////////////////////////

  clock_loc(1) = mpi_wtime()

  ! Calculate velocity derivatives dudx,dudy,dvdx,dvdy,dwdx,dwdy
  ! Enforces no stress condition at the top for u,v and zero vertical velocity for w 
  call filt_da (u, dudx, dudy, 1)
  call filt_da (v, dvdx, dvdy, 1)
  call filt_da (w, dwdx, dwdy, 2)   

  ! Calculate dudz,dvdz for on uv-nodes (1:nz levels) 
  call ddz_uv(u(1:ld,1:ny,0:nz), dudz(1:ld,1:ny,0:nz))
  call ddz_uv(v(1:ld,1:ny,0:nz), dvdz(1:ld,1:ny,0:nz))

  ! Calculate dwdz on w-nodes (0:nz-1 levels) 
  call ddz_w (w(1:ld,1:ny,0:nz), dwdz(1:ld,1:ny,0:nz))

  clock_loc(2) = mpi_wtime()
  time_deriv=time_deriv+clock_loc(2)-clock_loc(1)

#ifdef SCALAR
  !//////////////////////////////////////////////////////
  !/// WALLSTRESS AND DERIVATIVES SCALAR              ///
  !//////////////////////////////////////////////////////

  clock_loc(1) = mpi_wtime()
  call filt_da(theta, dTdx, dTdy, 3)
  call ddz_uv(theta(1:ld,1:ny,0:nz), dTdz(1:ld,1:ny,0:nz))
  call wallstress_scalar
  clock_loc(2) = mpi_wtime()
  time_scalar_deriv=time_scalar_deriv+clock_loc(2)-clock_loc(1)
#endif SCALAR

  !//////////////////////////////////////////////////////
  !/// WALLSTRESS                                     ///
  !//////////////////////////////////////////////////////

  clock_loc(1) = mpi_wtime()

  ! Calculate wall stress and derivatives at the wall 
  ! (txz, tyz, dudz, dvdz at jz=1) using the velocity log-law
  if (coord == 0) then
    call wallstress
  end if

  clock_loc(2) = mpi_wtime()
  time_wall=time_wall+clock_loc(2)-clock_loc(1)

  !//////////////////////////////////////////////////////
  !/// CONVEC                                         ///
  !//////////////////////////////////////////////////////
  ! Computes the rotation convective term in physical space

  clock_loc(1) = mpi_wtime()

  call convec

  clock_loc(2) = mpi_wtime()
  time_convec=time_convec+clock_loc(2)-clock_loc(1)

  !//////////////////////////////////////////////////////
  !/// SUBGRID SCALE MODEL                            ///
  !//////////////////////////////////////////////////////

  ! Calculate turbulent (subgrid) stress for entire domain
  clock_loc(1) = mpi_wtime()

  call sgs_stag

  clock_loc(2) = mpi_wtime()
  time_sgs=time_sgs+clock_loc(2)-clock_loc(1)

#ifdef CPS
  call synchronize_cps
#endif CPS

#ifdef SCALAR
  !//////////////////////////////////////////////////////
  !/// TIME STEPPING TEMPERATURE                      ///
  !//////////////////////////////////////////////////////
  clock_loc(1) = mpi_wtime()
  call theta_all_in_one 
  clock_loc(2) = mpi_wtime()
  time_temperature=time_temperature+clock_loc(2)-clock_loc(1)
#endif SCALAR

  !//////////////////////////////////////////////////////
  !/// UPDATE RHS WITH DIVERGENCE SGS SHEAR STRESSES  ///
  !//////////////////////////////////////////////////////
 
  clock_loc(1) = mpi_wtime()

  call divstress_uv
  call divstress_w

#ifdef CPS 
  if(inflow) call inflow_cond_cps 
#endif CPS 

  clock_loc(2) = mpi_wtime()
  time_divstress=time_divstress+clock_loc(2)-clock_loc(1)

  !//////////////////////////////////////////////////////
  !/// CALCULATE txx,txy,txz,tyy,tyz,tzz stats         //
  !//  later during the timestep these are overwritten //
  !//////////////////////////////////////////////////////

  clock_loc(1) = mpi_wtime()

  call tavg_compute_tau 

  clock_loc(2) = mpi_wtime()
  time_output=time_output+clock_loc(2)-clock_loc(1)

#ifdef TURBINES
  !//////////////////////////////////////////////////////
  !/// APPLIED TURBINE FORCING                        ///
  !//////////////////////////////////////////////////////

  clock_loc(1) = mpi_wtime()
        
  call turbine_forcing

  clock_loc(2) = mpi_wtime()
  time_turbine=time_turbine+clock_loc(2)-clock_loc(1)
#endif TURBINES

#ifdef WINDBREAKS
  !//////////////////////////////////////////////////////
  !/// APPLIED WINDBREAK FORCING                        ///
  !//////////////////////////////////////////////////////

  call windbreak_forcing
#endif WINDBREAKS

#ifdef ATM
  !//////////////////////////////////////////////////////
  !/// APPLIED TURBINE FORCING                        ///
  !//////////////////////////////////////////////////////

  clock_loc(1) = mpi_wtime()

  fxa = 0._rp
  fya = 0._rp
  fza = 0._rp

  call atm_lesgo_forcing

  RHSx(:,:,1:nz-1) = RHSx(:,:,1:nz-1) + fxa(:,:,1:nz-1)
  RHSy(:,:,1:nz-1) = RHSy(:,:,1:nz-1) + fya(:,:,1:nz-1)
  RHSz(:,:,1:nz-1) = RHSz(:,:,1:nz-1) + fza(:,:,1:nz-1)

  clock_loc(2) = mpi_wtime()
  time_atm=time_atm+clock_loc(2)-clock_loc(1)
#endif ATM

  !//////////////////////////////////////////////////////
  !/// INTERMEDIATE VELOCITY                          ///
  !//////////////////////////////////////////////////////
  ! Calculate intermediate velocity field at 1:nz-1 are valid

  clock_loc(1) = mpi_wtime()

  ! Set RHS*_f if necessary (first timestep) 
  if (jt_total.eq.1) then
   ! For the first step put RHS_f=RHS, i.e. take an Euler step
   RHSx_ff(:,:,1:nz-1)=RHSx(:,:,1:nz-1)
   RHSy_ff(:,:,1:nz-1)=RHSy(:,:,1:nz-1)
   RHSz_ff(:,:,1:nz-1)=RHSz(:,:,1:nz-1)
   RHSx_f(:,:,1:nz-1)=RHSx(:,:,1:nz-1)
   RHSy_f(:,:,1:nz-1)=RHSy(:,:,1:nz-1)
   RHSz_f(:,:,1:nz-1)=RHSz(:,:,1:nz-1)
   if(coord == 0) write(*,*) '--> Using Euler for first step.'
  end if

  if (jt_total .eq. 2 .OR. (cont_reduce)) then
   ! For the second step, take a 2nd order AB step
   RHSx_ff(:,:,1:nz-1) = tfac1*RHSx(:,:,1:nz-1) + tfac2*RHSx_f(:,:,1:nz-1)
   RHSy_ff(:,:,1:nz-1) = tfac1*RHSy(:,:,1:nz-1) + tfac2*RHSy_f(:,:,1:nz-1)
   RHSz_ff(:,:,1:nz-1) = tfac1*RHSz(:,:,1:nz-1) + tfac2*RHSz_f(:,:,1:nz-1)
   if(coord == 0) write(*,*) '--> Using 2nd order AB for second step.'
   cont_reduce = .false.
  end if

  tconst1=dt*tadv1
  tconst2=dt*tadv2
  tconst3=dt*tadv3
  u(:,:,1:nz-1)=u(:,:,1:nz-1)+tconst1*RHSx(:,:,1:nz-1)+tconst2*RHSx_f(:,:,1:nz-1)+tconst3*RHSx_ff(:,:,1:nz-1)
  v(:,:,1:nz-1)=v(:,:,1:nz-1)+tconst1*RHSy(:,:,1:nz-1)+tconst2*RHSy_f(:,:,1:nz-1)+tconst3*RHSy_ff(:,:,1:nz-1)
  w(:,:,1:nz-1)=w(:,:,1:nz-1)+tconst1*RHSz(:,:,1:nz-1)+tconst2*RHSz_f(:,:,1:nz-1)+tconst3*RHSz_ff(:,:,1:nz-1) 

  clock_loc(2) = mpi_wtime()
  time_RHS=time_RHS+clock_loc(2)-clock_loc(1)

#ifdef LVLSET
  !//////////////////////////////////////////////////////
  !/// APPLIED LEVEL SET FORCING                      ///
  !//////////////////////////////////////////////////////
  ! Compute the level set IBM forces
  ! Actually it's setting intermediate velocity zero in solid

  clock_loc(1) = mpi_wtime()

  call level_set_forcing

  clock_loc(2) = mpi_wtime()
  time_level_set=time_level_set+clock_loc(2)-clock_loc(1)
#endif LVLSET

  !//////////////////////////////////////////////////////
  !/// PRESSURE SOLUTION                              ///
  !//////////////////////////////////////////////////////
  ! Solve Poisson equation for pressure
  clock_loc(1) = mpi_wtime()

  call press_stag

  clock_loc(2) = mpi_wtime()
  time_press=time_press+clock_loc(2)-clock_loc(1)

  !//////////////////////////////////////////////////////
  !/// PROJECTION STEP                                ///
  !//////////////////////////////////////////////////////
  clock_loc(1) = mpi_wtime()

  call project

  clock_loc(2) = mpi_wtime()
  time_project=time_project+clock_loc(2)-clock_loc(1)

  !//////////////////////////////////////////////////////
  !/// CALCULATE STATISTICS, MAKE FLOW SNAPSHOTS, ETC ///
  !//////////////////////////////////////////////////////

  clock_loc(1) = mpi_wtime()

  call output_loop

  clock_loc(2) = mpi_wtime()
  time_output=time_output+clock_loc(2)-clock_loc(1)

  !//////////////////////////////////////////////////////
  !/// CALCULATE PDF & SPECTRUM                       ///
  !//////////////////////////////////////////////////////
#ifdef TAVG_PDF
  clock_loc(1) = mpi_wtime()

  call pdf

  clock_loc(2) = mpi_wtime()
  time_pdf=time_pdf+clock_loc(2)-clock_loc(1)
#endif TAVG_PDF

#ifdef TAVG_SPECTRUM
  clock_loc(1) = mpi_wtime()

  call spec

  clock_loc(2) = mpi_wtime()
  time_spec=time_pdf+clock_loc(2)-clock_loc(1)
#endif TAVG_SPECTRUM

  !//////////////////////////////////////////////////////
  !/// SIMULATION MONITORING                          ///
  !//////////////////////////////////////////////////////

  if (modulo (jt_total, wbase) == 0) then

    ! Calculate energy in the domain and rms divergence of velocity 
    call monitor(ke,rmsdivvel)
       
    ! Get the ending time for the iteration
    clock(2) = mpi_wtime()
    clock_total(2) = mpi_wtime()

    if (coord == 0) then
    write(*,*)
    write(*,'(a)') '========================================'
    write(*,'(a)') 'Time step information:'
    write(*,'(a,i9)') '  Iteration: ', jt_total
    write(*,'(a,E15.7)') '  Time step: ', dt
    write(*,'(a,E15.7)') '  CFL: ', cfl
    write(*,'(a,E15.7)') '  Lagrangian CFL: ', cfl*cs_count ! Estimate,does not take varition of dt into account
    write(*,'(a,3E15.7)') '  AB3 TADV1, TADV2, TADV3: ', tadv1, tadv2, tadv3
    write(*,*) 
    write(*,'(a)') 'Flow field information:'          
    write(*,'(a,E15.7)') '  Velocity divergence metric: ', rmsdivvel
    write(*,'(a,E15.7)') '  Kinetic energy: ', ke
    write(*,*)
    write(*,'(1a)') 'Simulation wall times (s): '
    write(*,'(1a,E15.7)') '  Iteration: ', clock(2)-clock(1)
    dummy2=clock_total(2)-clock_total(1)
    write(*,'(1a,E15.7)') '  Cumulative: ', dummy2
    write(*,'(1a)') 'CPU time usage in code (in percentage)'
    write(*,'(1a,f10.2)')' Sgs_stag         : ', time_sgs           /dummy2*100
    write(*,'(1a,f10.2)')' Convec           : ', time_convec        /dummy2*100
    write(*,'(1a,f10.2)')' Pressure         : ', time_press         /dummy2*100
    write(*,'(1a,f10.2)')' Derivatives      : ', time_deriv         /dummy2*100
    write(*,'(1a,f10.2)')' Divstress        : ', time_divstress     /dummy2*100
    write(*,'(1a,f10.2)')' Project          : ', time_project       /dummy2*100
    write(*,'(1a,f10.2)')' RHS_updates      : ', time_RHS           /dummy2*100
    write(*,'(1a,f10.2)')' Wallstress       : ', time_wall          /dummy2*100
    write(*,'(1a,f10.2)')' CFL              : ', time_cfl           /dummy2*100
    write(*,'(1a,f10.2)')' Output           : ', time_output        /dummy2*100
#ifdef SCALAR 
    write(*,'(1a,f10.2)')' Scalar boundary  : ', (time_scalar)      /dummy2*100
    write(*,'(1a,f10.2)')' Scalar derivative: ', (time_scalar_deriv)/dummy2*100
    write(*,'(1a,f10.2)')' theta_all_in_one : ', (time_temperature) /dummy2*100
#endif 
    dummy3=0._rp
#ifdef SCALAR 
    dummy3 = dummy3 + time_scalar_deriv + time_temperature 
    write(*,'(1a,f10.2)')' Temperature      : ',  dummy3            /dummy2*100
#endif 
#ifdef TURBINES
    dummy3=dummy3+time_turbine
    write(*,'(1a,f10.2)')' Turbine          : ', time_turbine       /dummy2*100
#endif
#ifdef ATM
    dummy3=dummy3+time_atm
    write(*,'(1a,f10.2)')' ATM              : ', time_atm       /dummy2*100
#endif ATM
#if defined(TAVG_PDF)
    dummy3=dummy3+time_pdf
    write(*,'(1a,f10.2)')' PDF              : ', time_pdf      /dummy2*100
#endif
#if defined(TAVG_SPECTRUM)
    dummy3=dummy3+time_spec
    write(*,'(1a,f10.2)')' Spectrum              : ', time_spec      /dummy2*100
#endif
#ifdef LVLSET
    dummy3=dummy3+time_level_set
    write(*,'(1a,f10.2)')' Level Set forcing: ', time_level_set     /dummy2*100
    write(*,'(1a,f10.2)')' Level Set sgsstag: ', time_level_set2    /dummy2*100    
#endif
    dummy3=dummy3+time_press+time_convec+time_sgs+time_project+time_RHS
    dummy3=dummy3+time_deriv+time_divstress+time_wall+time_output+time_cfl
    write(*,'(1a,f10.2)')' Other            : ', 100._rp-dummy3     /dummy2*100
    write(*,'(a)') '========================================'
   end if

   ! Check if we are to check the allowable runtime
   if( runtime > 0 ) then
     ! Determine the maximum simulation time among all processors and send result around
     dummy2=clock_total(2)-clock_total(1)
     call mpi_allreduce(dummy2, rbuffer, 1, mpi_rp, mpi_MAX, mpi_comm_world, ierr)
       
     ! If maximum time is surpassed go to the end of the program
     if ( rbuffer.gt.runtime) then
       if(coord.eq.0) write(*,*) 'Specified runtime exceeded. Exiting simulation.'
       exit time_loop
     endif
    endif ! End check runtime
  endif ! End check wbase

  ! End simulation based on chosen averaging window
  if(total_time > tavg_end_time) exit time_loop

end do time_loop

if(rbuffer .gt. runtime) jt_total = jt_total + 1
if(total_time > tavg_end_time) jt_total = jt_total + 1

! Write output files
call checkpoint(2) ! Numbered files (additional security)
call checkpoint(1) ! Regular files
#ifdef TURBINES
call turbine_write_buffer(1) ! makes sure to write buffered data
#endif
#ifdef ATM
call atm_lesgo_finalize ! write the restart files
#endif
! Write PDF and Spectrum
#if defined(TAVG_PDF) 
call checkpoint_pdf(2)
call checkpoint_pdf(1)
#endif
#if defined(TAVG_SPECTRUM)
call checkpoint_spec(2)
call checkpoint_spec(1)
#endif
! Close HDF5
call h5close_f(hdf_error)

! Write simulation time
clock_total(2) = mpi_wtime()
if(coord.eq.0) write(*,"(a,e15.7)") 'Simulation wall time (s) : ',clock_total(2)-clock_total(1)

! Close MPI 
call mpi_barrier(mpi_comm_world, ierr)
call mpi_finalize(ierr)

if(coord.eq.0) write(*,'(a)') 'Simulation complete'

end program main
