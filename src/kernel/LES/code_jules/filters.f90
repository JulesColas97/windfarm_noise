!*******************************************************************************
subroutine test_filter_init
!*******************************************************************************
!
! Creates the kernels which will be used for filtering the field
!
use filters, only: delta,G_test,G_test_test
use param, only: rp,lh,ny,dx,dy,dz,pi,sgs_model,inxny
use fft, only: kx,ky
implicit none
real(rp) :: delta_test,kc2_test,delta_test_test,kc2_test_test

delta=(dx*dy*dz)**(1._rp/3._rp) ! nondimensional

! Include the normalization
G_test = inxny  ! normalization for the forward FFT

! Filter characteristic width
delta_test = 2.0_rp * sqrt(dx*dy)  ! "2d-delta", not full 3d one

! Calculate the kernel
kc2_test = (pi/(delta_test))**2
where (real(kx*kx + ky*ky,rp) >= kc2_test) G_test = 0._rp

! since our k2 has zero at Nyquist, we have to do this by hand
G_test(lh,:) = 0._rp
G_test(:,ny/2+1) = 0._rp

! Second test filter, if necessary
if (sgs_model==5) then  !--scale dependent dynamic

  ! Include the normalization
  G_test_test = inxny

  ! Filter characteristic width
  delta_test_test = 4.0_rp * sqrt(dx*dy)

  ! Calculate the kernel
  kc2_test_test = (pi/(delta_test_test))**2
  where (real(kx*kx + ky*ky,rp) >= kc2_test_test) G_test_test = 0._rp

  ! since our k2 has zero at Nyquist, we have to do this by hand
  G_test_test(lh,:)     = 0._rp
  G_test_test(:,ny/2+1) = 0._rp
endif

end subroutine test_filter_init


!*******************************************************************************
subroutine filtering(f,f4)
!*******************************************************************************
!
! Perform 2delta and 4delta filtering together
!
use param, only: rp,ld,ny,lh
use filters, only: G_test,G_test_test
use fft,only: forw,back
implicit none
real(rp), dimension(1:ld,1:ny), intent(inout) :: f
real(rp), dimension(1:ld,1:ny), intent(out) :: f4

!  Perform in-place FFT
call dfftw_execute_dft_r2c(forw,f,f)

! Perform f = G_test*f, emulating f as complex
! Nyquist frequency and normalization is taken care of with G_test_test
! 4 delta filter
f4(1:ld,1:ny)=f(1:ld,1:ny)
call mulr(f4(1:ld,1:ny),G_test_test(1:lh,1:ny))
call dfftw_execute_dft_c2r(back,f4,f4)

! 2 delta filter
call mulr(f (1:ld,1:ny),G_test     (1:lh,1:ny))
call dfftw_execute_dft_c2r(back,f,f)

end subroutine filtering 
