!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! IMPORTANT: When using this routine, time-averaged velocities have to be
! substracted from the velocities 
! (this has to be done in line 52: call tavg_spectrum(u(1:nx,1:ny,1:nz-1),... 
!  and for v and w respectively)
! Note that this can be only done when the simulation is restarted and 
! time-averaged quantities are known.
!
! Calculates and writes out time averaged velocity-spectra
! Contains the following subroutines: 
!
! - spec            : Driver routine, determines when to calculate and write out spec
! - spec_init       : Initilizes arrays used in computation of spectra
! - checkpoint_spec : Writes Spec to file
! - tavg_spectrum   : Calculates stream- and spanwise spectrum of u,v,w
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#if defined(TAVG_SPECTRUM)
subroutine spec
use sim_param, only: u,v,w
use param, only : jt_total,total_time,dt,rp,coord,nx,ny,nz
use param, only : checkpoint_nskip,tavg_nskip,tavg_nstart,tavg_start_time,tavg_end_time
use param, only : spec_dt,spec_total_time
use io, only : spec_initialized
use io, only : tavg_u_stream_spec,tavg_u_span_spec
use io, only : tavg_v_stream_spec,tavg_v_span_spec
use io, only : tavg_w_stream_spec,tavg_w_span_spec
implicit none
! Store data both on numered as well as normal spectrum files
if ( modulo (jt_total, checkpoint_nskip) == 0) call checkpoint_spec(2)
!Numbered files
if ( modulo (jt_total, checkpoint_nskip) == 0) call checkpoint_spec(1)
!Regular files

!  Determine if time summations are to be calculated
if (total_time >= tavg_start_time .and. total_time <= tavg_end_time) then
 spec_dt = spec_dt + dt

 ! Initialization:
 if(total_time >= tavg_start_time .and. .not.spec_initialized) then
  tavg_nstart = jt_total
  if (coord == 0) then
   write(*,*) '-------------------------------'
   write(*,"(1a,i9,1a,i9)") 'Starting Spectrum calculations from ', tavg_nstart
   write(*,*) '-------------------------------'
   write(*,*) ''
   write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   write(*,*) '!  ATTENTION                  !'
   write(*,*) '!  Did you adjust spec.f90    !'
   write(*,*) '!  before starting this       !'
   write(*,*) '!  simulation?                !'
   write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   write(*,*) ''
  endif  ! coord==0
  call spec_init
 endif ! end initialization

 if(mod(jt_total-tavg_nstart,tavg_nskip).eq.0.and.jt_total-tavg_nstart.gt.0) then
  call tavg_spectrum(u(1:nx,1:ny,1:nz-1),tavg_u_stream_spec(1:nx/2+1,1:nz-1),tavg_u_span_spec(1:ny/2+1,1:nz-1))
  call tavg_spectrum(v(1:nx,1:ny,1:nz-1),tavg_v_stream_spec(1:nx/2+1,1:nz-1),tavg_v_span_spec(1:ny/2+1,1:nz-1))
  call tavg_spectrum(w(1:nx,1:ny,1:nz-1),tavg_w_stream_spec(1:nx/2+1,1:nz-1),tavg_w_span_spec(1:ny/2+1,1:nz-1))
 
  spec_total_time= spec_total_time+spec_dt ! Update total averaging time
  spec_dt=0.0_rp ! To determine time interval with next call to "spec"
 endif
endif
end subroutine spec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Checkpoint routine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine checkpoint_spec (version)
use param, only : rp,path,nx,ny,nz
use io, only : spec_initialized
use io, only : tavg_u_stream_spec,tavg_u_span_spec
use io, only : tavg_v_stream_spec,tavg_v_span_spec
use io, only : tavg_w_stream_spec,tavg_w_span_spec
implicit none
integer,intent(IN) :: version
character (64) :: folder

! Make backup of the spec data files 
if( spec_initialized ) then
 folder=trim(path)//'output'
 call write_streamw_spec(version,tavg_u_stream_spec(1:nx/2+1,1:nz-1),'tavg_u_streamw_spec',folder,0,nz-1)
 call write_streamw_spec(version,tavg_v_stream_spec(1:nx/2+1,1:nz-1),'tavg_v_streamw_spec',folder,0,nz-1)
 call write_streamw_spec(version,tavg_w_stream_spec(1:nx/2+1,1:nz-1),'tavg_w_streamw_spec',folder,0,nz-1)
 call write_spanw_spec(version,tavg_u_span_spec(1:ny/2+1,1:nz-1),'tavg_u_spanw_spec',folder,0,nz-1)
 call write_spanw_spec(version,tavg_v_span_spec(1:ny/2+1,1:nz-1),'tavg_v_spanw_spec',folder,0,nz-1)
 call write_spanw_spec(version,tavg_w_span_spec(1:ny/2+1,1:nz-1),'tavg_w_spanw_spec',folder,0,nz-1)
endif
end subroutine checkpoint_spec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize the spec variables to zero or read from file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine spec_init
use param, only : rp,nx,ny,nz,pi,path,coord,mpi_rp,localComm,ierr
use param, only : spec_total_time,spec_dt
use io, only : spec_initialized
use io, only : tavg_u_stream_spec,tavg_v_stream_spec,tavg_w_stream_spec
use io, only : tavg_u_span_spec,tavg_v_span_spec,tavg_w_span_spec
implicit none
character (64) :: folder
logical :: exst

allocate(tavg_u_stream_spec(nx/2+1,1:nz-1))
allocate(tavg_v_stream_spec(nx/2+1,1:nz-1))
allocate(tavg_w_stream_spec(nx/2+1,1:nz-1))
allocate(tavg_u_span_spec(ny/2+1,1:nz-1))
allocate(tavg_v_span_spec(ny/2+1,1:nz-1))
allocate(tavg_w_span_spec(ny/2+1,1:nz-1))

folder=trim(path)//'output'

!Initialize time averaging structures 
  spec_dt           =0.0_rp
  spec_total_time   =0.0_rp

  tavg_u_stream_spec(:,:) = 0.0_rp
  tavg_v_stream_spec(:,:) = 0.0_rp
  tavg_w_stream_spec(:,:) = 0.0_rp
  tavg_u_span_spec(:,:)   = 0.0_rp
  tavg_v_span_spec(:,:)   = 0.0_rp
  tavg_w_span_spec(:,:)   = 0.0_rp

inquire (file=trim(folder)//'/tavg_u_streamw_spec.h5',exist=exst)
if (.not. exst) then ! Initialize when no previous data are present
  if (coord == 0) then
    write(*,*)'No previous Spec data - starting from scratch.'
  endif
else ! Read data from file
  if (coord == 0) then
   write(*,*) 'Read in Spec from previous run'
  endif
  !Send information from master to other tasks
  call mpi_bcast(spec_total_time,1,mpi_rp,coord,localComm,ierr)

  call read_streamw_spec(tavg_u_stream_spec(1:nx/2+1,1:nz-1),'tavg_u_streamw_spec',folder,0,nz-1)
  call read_streamw_spec(tavg_v_stream_spec(1:nx/2+1,1:nz-1),'tavg_v_streamw_spec',folder,0,nz-1)
  call read_streamw_spec(tavg_w_stream_spec(1:nx/2+1,1:nz-1),'tavg_w_streamw_spec',folder,0,nz-1)
  call read_spanw_spec(tavg_u_span_spec(1:ny/2+1,1:nz-1),'tavg_u_spanw_spec',folder,0,nz-1)
  call read_spanw_spec(tavg_v_span_spec(1:ny/2+1,1:nz-1),'tavg_v_spanw_spec',folder,0,nz-1)
  call read_spanw_spec(tavg_w_span_spec(1:ny/2+1,1:nz-1),'tavg_w_spanw_spec',folder,0,nz-1)
endif
! Set global switch that tavg of spec has been initialized:
spec_initialized=.true.
end subroutine spec_init

subroutine tavg_spectrum(f,tavg_f_stream_spec,tavg_f_span_spec)
use param, only:nx,ny,nz,rp,inxny
use param, only: spec_total_time,spec_dt
use fft, only : forw_1d,dummy_in,dummy_out,dummy2_in,dummy2_out
implicit none
real(rp),dimension(1:nx,1:ny,1:nz-1), intent(in) :: f
real(rp),dimension(1:nx/2+1,1:nz-1), intent(inout) :: tavg_f_stream_spec
real(rp),dimension(1:ny/2+1,1:nz-1), intent(inout) :: tavg_f_span_spec
real(rp),dimension(1:nx/2+1,1:nz-1) :: temp_spec
real(rp),dimension(1:ny/2+1,1:nz-1) :: temp_spec_s
real(rp) :: gg1,gg2
integer :: i,j,k

gg1=spec_total_time/(spec_total_time+spec_dt)
gg2=spec_dt        /(spec_total_time+spec_dt)

dummy_in   = 0._rp
dummy_out  = 0._rp
dummy2_in  = 0._rp
dummy2_out = 0._rp

do k=1,nz-1
 temp_spec(:,k) = 0._rp
 temp_spec_s(:,k) = 0._rp
 ! Streamwise spectrum is calculated: 1. FFT to spectral space
 do j = 1,ny
  dummy_in(1:nx) = f(1:nx,j,k)
  call dfftw_execute_dft_r2c(forw_1d, dummy_in, dummy_out)
  dummy_out = dummy_out/nx**0.5_rp
  ! 2. half is omitted due to symmetry, spectrum is averaged for y,and z
  ! direction:
  temp_spec(1,k) = temp_spec(1,k) &
               & + 0.5_rp*real(dummy_out(1)*conjg(dummy_out(1)))
  do i=2,nx/2
   temp_spec(i,k) = temp_spec(i,k) + real(dummy_out(i)*conjg(dummy_out(i)))
  end do
  temp_spec(nx/2+1,k) = temp_spec(nx/2+1,k) &
                    & + 0.5_rp*real(dummy_out(nx/2+1)*conjg(dummy_out(nx/2+1)))
 end do
! Spanwise spectrum is calculated: 1. FFT to spectral space
 do j = 1,nx
  dummy2_in(1:ny) = f(j,1:ny,k)
  call dfftw_execute_dft_r2c(forw_1d, dummy2_in, dummy2_out)
  dummy2_out = dummy2_out/ny**0.5_rp
  ! 2. half is omitted due to symmetry, spectrum is averaged for y,and z
  ! direction:
  temp_spec_s(1,k) = temp_spec_s(1,k) &
                   & + 0.5_rp*real(dummy2_out(1)*conjg(dummy2_out(1)))
  do i=2,ny/2
   temp_spec_s(i,k) = temp_spec_s(i,k) +real(dummy2_out(i)*conjg(dummy2_out(i)))
  end do
  temp_spec_s(ny/2+1,k) = temp_spec_s(ny/2+1,k) &
                      & + 0.5_rp*real(dummy2_out(ny/2+1)*conjg(dummy2_out(ny/2+1)))
 end do

 tavg_f_stream_spec(:,k) = gg1 * tavg_f_stream_spec(:,k) + gg2 * temp_spec(:,k)
 tavg_f_span_spec(:,k)   = gg1 * tavg_f_span_spec(:,k)   + gg2 * temp_spec_s(:,k)
end do
tavg_f_stream_spec  = tavg_f_stream_spec * inxny !/ dkx
tavg_f_span_spec  = tavg_f_span_spec * inxny !/ dky
! result is written to file in subroutine checkpoint using write_spectrum.f90
end subroutine tavg_spectrum
#endif

