!*******************************************************************************
subroutine init_fft
!*******************************************************************************
use param, only: rp,nx,ny,nx2,ny2,ny,ny2
#ifdef MEMORY
use param, only: ld_big
#endif MEMORY
use convec_mod, only: cc_big
use sim_param, only: data=>dummy_plane
#ifdef TAVG_SPECTRUM
use fft, only: dummy_in, dummy_out, dummy2_in, dummy2_out, forw_1d
#endif TAVG_SPECTRUM
use fft, only: forw, forw_big, back, back_big
use fft, only: FFTW_ESTIMATE, FFTW_UNALIGNED, FFTW_PATIENT
implicit none

#ifdef MEMORY
allocate(cc_big   (ld_big,ny2))       ; cc_big   =0.0_rp
#endif MEMORY

#ifdef TAVG_SPECTRUM
allocate(dummy_in(nx))
allocate(dummy_out(nx/2+1))
allocate(dummy2_in(ny))
allocate(dummy2_out(ny/2+1))
dummy_in(:)  = 0._rp
dummy_out(:) = 0._rp
dummy2_in(:)  = 0._rp
dummy2_out(:) = 0._rp
#endif TAVG_SPECTRUM

! Create the forward and backward plans for the unpadded and padded domains.
! Notice we are using FFTW_UNALIGNED since the arrays used will not be
! guaranteed to be memory aligned. 
#ifdef DEBUG
call dfftw_plan_dft_r2c_2d(forw    ,nx ,ny ,data  ,data  ,FFTW_ESTIMATE,FFTW_UNALIGNED)
call dfftw_plan_dft_c2r_2d(back    ,nx ,ny ,data  ,data  ,FFTW_ESTIMATE,FFTW_UNALIGNED)
call dfftw_plan_dft_r2c_2d(forw_big,nx2,ny2,cc_big,cc_big,FFTW_ESTIMATE,FFTW_UNALIGNED)
call dfftw_plan_dft_c2r_2d(back_big,nx2,ny2,cc_big,cc_big,FFTW_ESTIMATE,FFTW_UNALIGNED)
#ifdef TAVG_SPECTRUM
call dfftw_plan_dft_r2c_1d(forw_1d,nx,dummy_in ,dummy_out ,FFTW_ESTIMATE,FFTW_UNALIGNED)
call dfftw_plan_dft_r2c_1d(forw_1d,ny,dummy2_in,dummy2_out,FFTW_ESTIMATE,FFTW_UNALIGNED)
#endif TAVG_SPECTRUM
#else
call dfftw_plan_dft_r2c_2d(forw    ,nx ,ny ,data  ,data  ,FFTW_PATIENT,FFTW_UNALIGNED)
call dfftw_plan_dft_c2r_2d(back    ,nx ,ny ,data  ,data  ,FFTW_PATIENT,FFTW_UNALIGNED)
call dfftw_plan_dft_r2c_2d(forw_big,nx2,ny2,cc_big,cc_big,FFTW_PATIENT,FFTW_UNALIGNED)
call dfftw_plan_dft_c2r_2d(back_big,nx2,ny2,cc_big,cc_big,FFTW_PATIENT,FFTW_UNALIGNED)
#ifdef TAVG_SPECTRUM
call dfftw_plan_dft_r2c_1d(forw_1d,nx,dummy_in ,dummy_out ,FFTW_PATIENT,FFTW_UNALIGNED)
call dfftw_plan_dft_r2c_1d(forw_1d,ny,dummy2_in,dummy2_out,FFTW_PATIENT,FFTW_UNALIGNED)
#endif TAVG_SPECTRUM
#endif DEBUG

#ifdef MEMORY
deallocate(cc_big)
#endif MEMORY

call init_wavenumber

end subroutine init_fft


!*******************************************************************************
subroutine init_wavenumber
!*******************************************************************************
use param, only: rp,lh,ny,L_x,L_y,pi,lh,ny
use fft, only: kx,ky
implicit none
integer :: jx,jy

do jx=1,lh-1
  kx(jx,:) = real(jx-1,kind=rp)
end do

do jy=1,ny
  ky(:,jy) = real(modulo(jy - 1 + ny/2,ny) - ny/2,kind=rp)
end do

! Nyquist: makes doing derivatives easier
kx(lh,:)    =0.0_rp
ky(lh,:)    =0.0_rp
kx(:,ny/2+1)=0.0_rp
ky(:,ny/2+1)=0.0_rp

! for the aspect ratio change
kx = 2._rp*pi/L_x*kx
ky = 2._rp*pi/L_y*ky

end subroutine init_wavenumber
