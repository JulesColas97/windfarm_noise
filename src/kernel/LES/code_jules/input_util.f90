!*******************************************************************************
subroutine read_input_conf
!*******************************************************************************
use param, only: path,coord
implicit none
integer, parameter :: lun = 1
character (128) :: buff
integer :: block_entry_pos,block_exit_pos,equal_pos,ios,line
logical :: exst
character(20) :: ucstr

! Check that the configuration file exists
inquire (file=trim(path)//'input.conf', exist=exst)

if (exst) then
  open (lun, file=trim(path)//'input.conf', action='read')
else
  write(*,*)  '1 file does not exist'
end if

line = 0
do

  call readline(lun,line,buff,block_entry_pos,block_exit_pos,equal_pos,ios)

  if (ios /= 0) exit

  ! for now, invalid format if no block entry found
  if (block_entry_pos == 0) then
    write(*,*)  '2 block entry not found '
  end if

  ! Find block
  ucstr=''
  call uppercase(buff(1:block_entry_pos-1),ucstr)

  select case(ucstr)

  case ('DOMAIN')
    call domain_block
  case ('MODEL')
    call model_block
#ifdef CORIOLIS
  case ('CORIOLIS')
    call coriolis_block
#endif CORIOLIS
#ifdef SCALAR
  case ('SCALAR')
    call scalar_block
#endif SCALAR
#ifdef BAROCLINIC
  case ('BAROCLINIC')
    call baroclinic_block
#endif BAROCLINIC    
  case ('TIME')
    call time_block
  case ('FLOW_COND')
    call flow_cond_block
  case ('OUTPUT')
    call output_block
#ifdef LVLSET
  case ('LEVEL_SET')
    call level_set_block
#endif LVLSET
#ifdef TURBINES
  case ('TURBINES')
    call turbines_block
#endif TURBINES
#ifdef WINDBREAKS
  case ('WINDBREAKS')
    call windbreak_block
#endif WINDBREAKS

  case default
    if(coord == 0) write(*,*) 'Found unused input block: '                     &
                   // buff(1:block_entry_pos-1)
    ! Now need to 'fast-forward' untile we reach the end of the block
    do while ( block_exit_pos == 0 )
      call readline(lun,line,buff,block_entry_pos,block_exit_pos,equal_pos,ios)
        if (ios /= 0) exit ! exit if end of file is reached
    enddo
  end select
  
end do

close (lun)

contains

!*******************************************************************************
subroutine domain_block
!*******************************************************************************
use param, only: Nx, Ny, Nz_tot, z_i, L_x, L_y, L_z, grid_stretch, grid_size
implicit none
character(*), parameter :: block_name = 'DOMAIN'
character(20) :: ucstr

do 

  call readline(lun,line,buff,block_entry_pos,block_exit_pos,equal_pos,ios)
  if (ios /= 0) write(*,*)  'Bad read in block'

  if( block_exit_pos == 0 ) then

    ! Check that the data entry conforms to correct format
    call checkentry

    ucstr=''
    call uppercase(buff(1:equal_pos-1),ucstr)
    select case(ucstr)

    case ('NX')
      read (buff(equal_pos+1:), *) Nx
    case ('NY')
      read (buff(equal_pos+1:), *) Ny
    case ('NZ') 
      read (buff(equal_pos+1:), *) Nz_tot
    case ('Z_I')
      read (buff(equal_pos+1:), *) z_i
    case ('GRID_STRETCH')
      read (buff(equal_pos+1:), *) grid_stretch
    case ('GRID_SIZE')
      read (buff(equal_pos+1:), *) grid_size
    case ('LX')
      read (buff(equal_pos+1:), *) L_x
    case ('LY')
      read (buff(equal_pos+1:), *) L_y
    case ('LZ')
      read (buff(equal_pos+1:), *) L_z

    case default
      if(coord == 0) write(*,*)  'Found unused data value in '                 &
                     // block_name // ' block: ' // buff(1:equal_pos-1)
    end select

  elseif( block_exit_pos == 1 ) then
    return
  else
    write(*,*) block_name // ' data block not formatted correctly: '           &
                          // buff(1:equal_pos-1)
    stop
  endif
     
enddo

end subroutine domain_block


!*******************************************************************************
subroutine model_block
!*******************************************************************************
use param, only: sgs_model, wall_damp_exp, cs_count, clip, co
implicit none
character(*), parameter :: block_name = 'MODEL'
character(22) :: ucstr

do 

  call readline(lun,line,buff,block_entry_pos,block_exit_pos,equal_pos,ios)
  if (ios /= 0) write(*,*) 'Bad read in block'

  if( block_exit_pos == 0 ) then

    ! Check that the data entry conforms to correct format
    call checkentry

    ucstr=''
    call uppercase(buff(1:equal_pos-1),ucstr)
    select case(ucstr)

    case ('SGS_MODEL')
      read (buff(equal_pos+1:), *) sgs_model
    case ('WALL_DAMP_EXP') 
      read (buff(equal_pos+1:), *) wall_damp_exp
    case ('CS_COUNT')
      read (buff(equal_pos+1:), *) cs_count
    case ('CLIP')
      read (buff(equal_pos+1:), *) clip
    case ('CO')
      read (buff(equal_pos+1:), *) Co

    case default
      if(coord == 0) write(*,*)  'Found unused data value in '                 &
                     // block_name // ' block: ' // buff(1:equal_pos-1)
    end select

  elseif( block_exit_pos == 1 ) then
    return
  else
    write(*,*) block_name // ' data block not formatted correctly: '           &
                          // buff(1:equal_pos-1)
    stop
  endif

enddo

end subroutine model_block


#ifdef CORIOLIS
!*******************************************************************************
subroutine coriolis_block
!*******************************************************************************
use param, only: neutral_ekman_test, coriolis_forcing, coriol, coriol_x, coriol_y, ug, vg
use param, only: z_wind, wind_angle, adjust_wind_angle, kpw, kiw, kdw
implicit none
character(*), parameter :: block_name = 'CORIOLIS'
character(20) :: ucstr

do 

  call readline(lun,line,buff,block_entry_pos,block_exit_pos,equal_pos,ios)
  if (ios /= 0) write(*,*) 'Bad read in block'

  if( block_exit_pos == 0 ) then

    ! Check that the data entry conforms to correct format
    call checkentry

    ucstr='' 
    call uppercase(buff(1:equal_pos-1),ucstr)
    select case(ucstr)

    case ('NEUTRAL_EKMAN_TEST')
      read (buff(equal_pos+1:), *) neutral_ekman_test 
    case ('CORIOLIS_FORCING')
      read (buff(equal_pos+1:), *) coriolis_forcing
    case ('CORIOL')
      read (buff(equal_pos+1:), *) coriol
    case ('CORIOL_X')
      read (buff(equal_pos+1:), *) coriol_x
    case ('CORIOL_Y')
      read (buff(equal_pos+1:), *) coriol_y
    case ('UG')
      read (buff(equal_pos+1:), *) ug
    case ('VG')
      read (buff(equal_pos+1:), *) vg
    case ('WIND_ANGLE')
      read (buff(equal_pos+1:), *) wind_angle
    case ('Z_WIND')
      read (buff(equal_pos+1:), *) z_wind 
    case ('ADJUST_WIND_ANGLE')
      read (buff(equal_pos+1:), *) adjust_wind_angle
    case ('KPW')
      read (buff(equal_pos+1:), *) kpw
    case ('KIW')
      read (buff(equal_pos+1:), *) kiw
    case ('KDW')
      read (buff(equal_pos+1:), *) kdw

    case default
      if(coord == 0) write(*,*)  'Found unused data value in '                 &
                     // block_name // ' block: ' // buff(1:equal_pos-1)
    end select

  elseif( block_exit_pos == 1 ) then
    return
  else
    write(*,*) block_name // ' data block not formatted correctly: '           &
                          // buff(1:equal_pos-1)
    stop
  endif

enddo

end subroutine coriolis_block
#endif 


#ifdef BAROCLINIC
!*******************************************************************************
subroutine baroclinic_block
!*******************************************************************************
use param, only: ug_delta, vg_delta
use param, only: bar_start, bar_end, temp_adv
implicit none
character(*), parameter :: block_name = 'BAROCLINIC'
character(20) :: ucstr

do 

  call readline(lun,line,buff,block_entry_pos,block_exit_pos,equal_pos,ios)
  if (ios /= 0) write(*,*) 'Bad read in block'

  if( block_exit_pos == 0 ) then

    ! Check that the data entry conforms to correct format
    call checkentry

    ucstr='' 
    call uppercase(buff(1:equal_pos-1),ucstr)
    select case(ucstr)

    case ('UG_DELTA')
      read (buff(equal_pos+1:), *) ug_delta
    case ('VG_DELTA')
      read (buff(equal_pos+1:), *) vg_delta
    case ('BAR_START')
      read (buff(equal_pos+1:), *) bar_start
    case ('BAR_END')
      read (buff(equal_pos+1:), *) bar_end
    case ('TEMP_ADV')
      read (buff(equal_pos+1:), *) temp_adv

    case default
      if(coord == 0) write(*,*)  'Found unused data value in '                 &
                     // block_name // ' block: ' // buff(1:equal_pos-1)
    end select

  elseif( block_exit_pos == 1 ) then
    return
  else
    write(*,*) block_name // ' data block not formatted correctly: '           &
                          // buff(1:equal_pos-1)
    stop
  endif

enddo

end subroutine baroclinic_block
#endif


#ifdef SCALAR 
!*******************************************************************************
subroutine scalar_block
!*******************************************************************************
use param, only: gabls_test, sullivan_test, neutral_test, Pr
use param, only: manual_test, cap_height, cap_strength, strat_height, strat_strength
use param, only: surface_cond, surface_rate, surface_flux
use param, only: damping_x,inv_strength, T_scale, theta_s1, wt_s, cooling_rate, zo_scalar 
use param, only: T_init, hbc, ubc, damping_method, dmpfr, ra_damp_exponent, damp_hgt 
use param, only: temp_sink, z_heat_min, z_heat_max, time_heat, Cp, Ci, Cd
implicit none
character(*), parameter :: block_name = 'SCALAR'
character(20) :: ucstr

do

  call readline(lun,line,buff,block_entry_pos,block_exit_pos,equal_pos,ios)
  if (ios /= 0) write(*,*) 'Bad read in block'

  if( block_exit_pos == 0 ) then

    ! Check that the data entry conforms to correct format
    call checkentry

    ucstr='' 
    call uppercase(buff(1:equal_pos-1),ucstr)
    select case(ucstr)

    case ('GABLS_TEST')
      read (buff(equal_pos+1:), *) gabls_test
    case ('SULLIVAN_TEST')
      read (buff(equal_pos+1:), *) sullivan_test
    case ('NEUTRAL_TEST')
      read (buff(equal_pos+1:), *) neutral_test
    case ('MANUAL_TEST')
      read (buff(equal_pos+1:), *) manual_test
    case ('CAP_HEIGHT')
      read (buff(equal_pos+1:), *) cap_height
    case ('STRAT_HEIGHT')
      read (buff(equal_pos+1:), *) strat_height
    case ('CAP_STRENGTH')
      read (buff(equal_pos+1:), *) cap_strength
    case ('STRAT_STRENGTH')
      read (buff(equal_pos+1:), *) strat_strength
    case ('SURFACE_COND')
      read (buff(equal_pos+1:), *) surface_cond
    case ('SURFACE_RATE')
      read (buff(equal_pos+1:), *) surface_rate
    case ('SURFACE_FLUX')
      read (buff(equal_pos+1:), *) surface_flux      
    case ('PR')
      read (buff(equal_pos+1:), *) Pr
    case ('ZO_SCALAR')
      read (buff(equal_pos+1:), *) zo_scalar  
    case ('INV_STRENGTH')
      read (buff(equal_pos+1:), *) inv_strength
    case ('T_SCALE')
      read (buff(equal_pos+1:), *) T_scale
    case ('THETA_S1')
      read (buff(equal_pos+1:), *) theta_s1
    case ('WT_S')
      read (buff(equal_pos+1:), *) wt_s
    case ('COOLING_RATE')
      read (buff(equal_pos+1:), *) cooling_rate
    case ('T_INIT')
      read (buff(equal_pos+1:), *) T_init
    case ('HBC')
      read (buff(equal_pos+1:), *) hbc
    case ('UBC')
      read (buff(equal_pos+1:), *) ubc
    case ('DAMPING_METHOD')
      read (buff(equal_pos+1:), *) damping_method
    case ('DAMPING_X')
      read (buff(equal_pos+1:), *) damping_x
    case ('DMPFR')
      read (buff(equal_pos+1:), *) dmpfr
    case ('DAMP_HGT')
      read (buff(equal_pos+1:), *) damp_hgt 
    case ('RA_DAMP_EXPONENT')
      read (buff(equal_pos+1:), *) ra_damp_exponent
    case ('TEMP_SINK')
      read (buff(equal_pos+1:), *) temp_sink
    case ('Z_HEAT_MIN')
      read (buff(equal_pos+1:), *) z_heat_min
    case ('Z_HEAT_MAX')
      read (buff(equal_pos+1:), *) z_heat_max
    case ('TIME_HEAT')
      read (buff(equal_pos+1:), *) time_heat
    case ('CP')
      read (buff(equal_pos+1:), *) Cp
    case ('CI')
      read (buff(equal_pos+1:), *) Ci
    case ('CD')
      read (buff(equal_pos+1:), *) Cd

    case default
      if(coord == 0) write(*,*)  'Found unused data value in '                 &
                     // block_name // ' block: ' // buff(1:equal_pos-1)
    end select

  elseif( block_exit_pos == 1 ) then
    return
  else
    write(*,*) block_name // ' data block not formatted correctly: '           &
                          // buff(1:equal_pos-1)
    stop
  endif

enddo

end subroutine scalar_block
#endif 


!*******************************************************************************
subroutine time_block
!*******************************************************************************
use param, only: nsteps, runtime, use_cfl_dt, cfl, dt
implicit none
character(*), parameter :: block_name = 'TIME'
character(20) :: ucstr

do

  call readline(lun, line,buff,block_entry_pos,block_exit_pos,equal_pos,ios)
  if (ios /= 0) write(*,*) 'Bad read in block'


  if( block_exit_pos == 0 ) then

    ! Check that the data entry conforms to correct format
    call checkentry

    ucstr='' 
    call uppercase(buff(1:equal_pos-1),ucstr)
    select case(ucstr)

    case ('NSTEPS')
      read (buff(equal_pos+1:), *) nsteps
    case ('RUNTIME')
      read (buff(equal_pos+1:), *) runtime
    case ('USE_CFL_DT') 
      read (buff(equal_pos+1:), *) use_cfl_dt
    case ('CFL')
      read (buff(equal_pos+1:), *) cfl
    case('DT')
      read (buff(equal_pos+1:), *) dt

    case default
      if(coord == 0) write(*,*)  'Found unused data value in '                 &
                     // block_name // ' block: ' // buff(1:equal_pos-1)
    end select

  elseif( block_exit_pos == 1 ) then
    return
  else
    write(*,*) block_name // ' data block not formatted correctly: '           &
                          // buff(1:equal_pos-1)
    stop
  endif

enddo

end subroutine time_block


!*******************************************************************************
subroutine flow_cond_block
!*******************************************************************************
use param, only: initu, initu_interp, lbc_mom, zo, inflow_velocity
use param, only: inflow, fringe_factor, symmetric_fringe, fringe_region_end, fringe_region_len
use param, only: fringe_region_end_y, fringe_region_len_y
use param, only: use_mean_p_force, u_star
implicit none
character(*), parameter :: block_name = 'FLOW_COND'
character(20) :: ucstr

do

  call readline( lun, line,buff,block_entry_pos,block_exit_pos,equal_pos,ios)
  if (ios /= 0) write(*,*) 'Bad read in block'

  if( block_exit_pos == 0 ) then

    ! Check that the data entry conforms to correct format
    call checkentry

    ucstr='' 
    call uppercase(buff(1:equal_pos-1),ucstr)
    select case(ucstr)

    case ('INITU')
      read (buff(equal_pos+1:), *) initu
    case ('INITU_INTERP')
      read (buff(equal_pos+1:), *) initu_interp
    case ('LBC_MOM')
      Read (buff(equal_pos+1:), *) lbc_mom
    case ('ZO')
      read (buff(equal_pos+1:), *) zo
    case ('INFLOW')
      read (buff(equal_pos+1:), *) inflow
    case ('SYMMETRIC_FRINGE')
      read (buff(equal_pos+1:), *) symmetric_fringe
    case ('FRINGE_FACTOR')
      read (buff(equal_pos+1:), *) fringe_factor 
    case ('FRINGE_REGION_END')
      read (buff(equal_pos+1:), *) fringe_region_end
    case ('FRINGE_REGION_LEN')
      read (buff(equal_pos+1:), *) fringe_region_len
    case ('FRINGE_REGION_END_Y')
      read (buff(equal_pos+1:), *) fringe_region_end_y
    case ('FRINGE_REGION_LEN_Y')
      read (buff(equal_pos+1:), *) fringe_region_len_y
    case ('INFLOW_VELOCITY')
      read (buff(equal_pos+1:), *) inflow_velocity
    case ('USE_MEAN_P_FORCE')
      read (buff(equal_pos+1:), *) use_mean_p_force
    case ('U_STAR')
      read (buff(equal_pos+1:), *) u_star

    case default
      if(coord == 0) write(*,*)  'Found unused data value in '                 &
                     // block_name // ' block: ' // buff(1:equal_pos-1)
    end select

  elseif( block_exit_pos == 1 ) then
    return
  else
    write(*,*) block_name // ' data block not formatted correctly: '           &
                          // buff(1:equal_pos-1)
    stop
  endif

enddo

end subroutine flow_cond_block


!*******************************************************************************
subroutine output_block
!*******************************************************************************
use param, only: wbase, checkpoint_data, checkpoint_nskip, tavg_calc, tavg_nskip
use param, only: tavg_start_time, tavg_end_time, domain_calc, itype
use param, only: domain_nstart, domain_nskip, nbins
use param, only: bin_min_u, bin_min_v, bin_min_w
use param, only: bin_max_u, bin_max_v, bin_max_w
implicit none
character(*), parameter :: block_name = 'OUTPUT'
character(20) :: ucstr

do

  call readline(lun,line,buff,block_entry_pos,block_exit_pos,equal_pos,ios)
  if (ios /= 0) write(*,*) 'Bad read in block'

  if( block_exit_pos == 0 ) then

    ! Check that the data entry conforms to correct format
    call checkentry

    ucstr='' 
    call uppercase(buff(1:equal_pos-1),ucstr)
    select case(ucstr)

    case ('WBASE')
      read (buff(equal_pos+1:), *) wbase
    case ('CHECKPOINT_DATA')
      read (buff(equal_pos+1:), *) checkpoint_data
    case ('CHECKPOINT_NSKIP')
      read (buff(equal_pos+1:), *) checkpoint_nskip
    case ('TAVG_CALC')
      read (buff(equal_pos+1:), *) tavg_calc
    case ('TAVG_NSKIP')
      read (buff(equal_pos+1:), *) tavg_nskip
    case ('TAVG_START_TIME')
      read (buff(equal_pos+1:), *) tavg_start_time
    case ('TAVG_END_TIME')
      read (buff(equal_pos+1:), *) tavg_end_time
    case ('DOMAIN_CALC')
      read (buff(equal_pos+1:), *) domain_calc
    case ('ITYPE')
      read (buff(equal_pos+1:), *) itype 
    case ('DOMAIN_NSTART')
      read (buff(equal_pos+1:), *) domain_nstart
    case ('DOMAIN_NSKIP')
      read (buff(equal_pos+1:), *) domain_nskip
    case ('NBINS')
      read (buff(equal_pos+1:), *) nbins
    case ('BIN_MIN_U')
      read (buff(equal_pos+1:), *) bin_min_u
    case ('BIN_MAX_U')
      read (buff(equal_pos+1:), *) bin_max_u
    case ('BIN_MIN_V')
      read (buff(equal_pos+1:), *) bin_min_v
    case ('BIN_MAX_V')
      read (buff(equal_pos+1:), *) bin_max_v
    case ('BIN_MIN_W')
      read (buff(equal_pos+1:), *) bin_min_w
    case ('BIN_MAX_W')
      read (buff(equal_pos+1:), *) bin_max_w

    case default
      if(coord == 0) write(*,*)  'Found unused data value in '                 &
                     // block_name // ' block: ' // buff(1:equal_pos-1)
    end select

  elseif( block_exit_pos == 1 ) then
    return
  else
    write(*,*) block_name // ' data block not formatted correctly: '           &
                          // buff(1:equal_pos-1)
    stop
  endif

enddo

end subroutine output_block


#ifdef WINDBREAKS
!*******************************************************************************
subroutine windbreak_block
!*******************************************************************************
use windbreak_mod, only: nloc_wb, filter_size_wb, filter_cutoff_wb
implicit none
character(*), parameter :: block_name = 'WINDBREAKS'
character(20) :: ucstr

do

  call readline(lun,line,buff,block_entry_pos,block_exit_pos,equal_pos,ios)
  if (ios /= 0) write(*,*) 'Bad read in block'

  if( block_exit_pos == 0 ) then

    ! Check that the data entry conforms to correct format
    call checkentry

    ucstr=''
    call uppercase(buff(1:equal_pos-1),ucstr)
    select case(ucstr)

    case ('NLOC_WB')
      read (buff(equal_pos+1:), *) nloc_wb
    case ('FILTER_SIZE_WB')
      read (buff(equal_pos+1:), *) filter_size_wb
    case ('FILTER_CUTOFF_WB')
      read (buff(equal_pos+1:), *) filter_cutoff_wb

    case default
      if(coord == 0) write(*,*)  'Found unused data value in '                 &
                     // block_name // ' block: ' // buff(1:equal_pos-1)
    end select

  elseif( block_exit_pos == 1 ) then
    return
  else
    write(*,*) block_name // ' data block not formatted correctly: '           &
                          // buff(1:equal_pos-1)
    stop
  endif

enddo

end subroutine windbreak_block
#endif WINDBREAKS


#ifdef LVLSET
!*******************************************************************************
subroutine level_set_block
!*******************************************************************************
use level_set_mod, only: zo_level_set, zwall, coeff_phi_cutoff, coeff_phi_band
implicit none
character(*), parameter :: block_name = 'LEVEL_SET'
character(20) :: ucstr

do

  call readline(lun,line,buff,block_entry_pos,block_exit_pos,equal_pos,ios)
  if (ios /= 0) write(*,*) 'Bad read in block'

  if( block_exit_pos == 0 ) then

    ! Check that the data entry conforms to correct format
    call checkentry

    ucstr=''
    call uppercase(buff(1:equal_pos-1),ucstr)
    select case(ucstr)

    case ('ZO_LEVEL_SET')
      read (buff(equal_pos+1:), *) zo_level_set
    case ('ZWALL')
      read (buff(equal_pos+1:), *) zwall
    case ('COEFF_PHI_CUTOFF')
      read (buff(equal_pos+1:), *) coeff_phi_cutoff
    case ('COEFF_PHI_BAND')
      read (buff(equal_pos+1:), *) coeff_phi_band

    case default
      if(coord == 0) write(*,*)  'Found unused data value in '                 &
                     // block_name // ' block: ' // buff(1:equal_pos-1)
    end select

  elseif( block_exit_pos == 1 ) then
    return
  else
    write(*,*) block_name // ' data block not formatted correctly: '           &
                          // buff(1:equal_pos-1)
    stop
  endif

enddo

end subroutine level_set_block
#endif LVLSET


#ifdef TURBINES
!*******************************************************************************
subroutine turbines_block
!*******************************************************************************
use turbine_mod, only: nloc, dyn_theta1, filter_size, t_avg
use turbine_mod, only: filter_cutoff, tbase
implicit none
character(*), parameter :: block_name = 'TURBINES'
character(20) :: ucstr

do

  call readline(lun,line,buff,block_entry_pos,block_exit_pos,equal_pos,ios)
  if (ios /= 0) write(*,*) 'Bad read in block'

  if( block_exit_pos == 0 ) then

    ! Check that the data entry conforms to correct format
    call checkentry

    ucstr='' 
    call uppercase(buff(1:equal_pos-1),ucstr)
    select case(ucstr)

    case ('NLOC')
      read (buff(equal_pos+1:), *) nloc
    case ('DYN_THETA1')
      read (buff(equal_pos+1:), *) dyn_theta1
    case ('FILTER_SIZE')
      read (buff(equal_pos+1:), *) filter_size
    case ('T_AVG')
      read (buff(equal_pos+1:), *) t_avg
    case ('FILTER_CUTOFF')
      read (buff(equal_pos+1:), *) filter_cutoff
    case ('TBASE')
      read (buff(equal_pos+1:), *) tbase
    case default
      if(coord == 0) write(*,*)  'Found unused data value in '                 &
                     // block_name // ' block: ' // buff(1:equal_pos-1)
    end select

  elseif( block_exit_pos == 1 ) then
    return
  else
    write(*,*) block_name // ' data block not formatted correctly: '           &
                          // buff(1:equal_pos-1)
    stop
  endif

enddo

end subroutine turbines_block
#endif TURBINES


!*******************************************************************************
subroutine checkentry
!*******************************************************************************
implicit none

if(equal_pos == 0) then
  write(*,*) 'Bad read in block at line', line, ': ' // trim(adjustl(buff))
endif

! invalid if nothing after equals
if(len_trim (buff) == equal_pos) then
  write(*,*) 'nothing after equals sign in line', line 
endif

end subroutine checkentry  

end subroutine read_input_conf


!*******************************************************************************
subroutine readline(lun,line,buff,block_entry_pos,block_exit_pos,equal_pos,ios)
!*******************************************************************************
!
! This subroutine reads the specified line and determines the attributes of the
! contents of the line.
!
implicit none
integer, intent(in) :: lun
integer, intent(inout) :: line
character(*), intent(inout) :: buff
integer, intent(out) :: block_entry_pos,block_exit_pos,equal_pos, ios
character (*), parameter :: comment = '!'
character (*), parameter :: block_entry = '{'
character (*), parameter :: block_exit = '}'
character (*), parameter :: equal = '='
character (1), parameter :: fill_char = ' '
character (*), parameter :: whtspc_default = achar (9) // achar (32)
character (1) :: tmp (len (buff))
character (1) :: fill (len (buff))

block_entry_pos = 0
block_exit_pos = 0
equal_pos = 0
ios = -1

do     
  line = line + 1
  read (lun, '(a)', iostat=ios) buff
  if (ios /= 0) exit

   fill = fill_char
   tmp = transfer (buff, tmp)
   tmp = pack (tmp, scan (tmp, whtspc_default) == 0, fill)
   buff = transfer (tmp, buff)

  if (verify (buff, ' ') == 0) cycle  !--drop blank lines
  
  if (buff (1:len (comment)) == comment) cycle  !--drop comment lines

  block_entry_pos = index( buff, block_entry )
  block_exit_pos  = index( buff, block_exit )
  equal_pos       = index( buff, equal )

  exit
enddo 

end subroutine readline


!*******************************************************************************
subroutine uppercase(str,ucstr)
!*******************************************************************************
!
! convert specified string to upper case
!
character (len=*):: str
character (len=len_trim(str)):: ucstr
integer :: i, ilen, iav, ioffset, iqc, iquote

ilen=len_trim(str) 
ioffset=iachar('A')-iachar('a')
iquote=0
ucstr=str
do i=1,ilen
  iav=iachar(str(i:i))
  if(iquote==0 .and. (iav==34 .or.iav==39)) then
    iquote=1
    iqc=iav
    cycle
  end if
  if(iquote==1 .and. iav==iqc) then
    iquote=0
    cycle
  end if
  if (iquote==1) cycle
  if(iav >= iachar('a') .and. iav <= iachar('z')) then
    ucstr(i:i)=achar(iav+ioffset)
  else
    ucstr(i:i)=str(i:i)
  end if
end do

end subroutine uppercase
