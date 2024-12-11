!*******************************************************************************
subroutine ddx(f, dfdx)
!*******************************************************************************
!
! First derivative in x direction
! Derivatives and input array on vertical levels 1:nz-1 (see divstress routines)
! Output array has vertical levels 0:nz (see dummy array divstress routines)
!
use param, only: rp,ld,ny,nz,inxny
use fft, only: back,forw,kx
implicit none
real(rp), dimension(1:ld,1:ny,1:nz-1), intent(in) :: f
real(rp), dimension(1:ld,1:ny,0:nz  ), intent(inout) :: dfdx
integer :: k

do k=1,nz-1
  ! Go to spectral space
  dfdx(:,:,k)=inxny*f(:,:,k)
  call dfftw_execute_dft_r2c(forw,dfdx(:,:,k),dfdx(:,:,k))

  ! Calculate derivatives: kx zeros out Nyquist frequencies
  call muli(dfdx(1:ld,1:ny,k),kx,dfdx(1:ld,1:ny,k))

  ! Go back to physical space
  call dfftw_execute_dft_c2r(back,dfdx(:,:,k),dfdx(:,:,k))
enddo

end subroutine ddx


!*******************************************************************************
subroutine ddy(f, dfdy)
!*******************************************************************************
!
! First derivative in y direction
! Derivatives and input array on vertical levels 1:nz-1 (see divstress routines)
! Output array has vertical levels 0:nz (see dummy array divstress routines)
!
use param, only: rp,ld,ny,nz,inxny
use fft, only: back,forw,ky
implicit none 
real(rp), dimension(1:ld,1:ny,1:nz-1), intent(in) :: f
real(rp), dimension(1:ld,1:ny,0:nz  ), intent(inout) :: dfdy
integer :: k

do k=1,nz-1
  ! Go to spectral space
  dfdy(:,:,k)=inxny*f(:,:,k)  
  call dfftw_execute_dft_r2c(forw, dfdy(:,:,k), dfdy(:,:,k))

  ! Calculate derivatives: ky zeros out Nyquist frequencies
  call muli(dfdy(1:ld,1:ny,k),ky,dfdy(1:ld,1:ny,k))

  ! Go back to physical space
  call dfftw_execute_dft_c2r(back, dfdy(:,:,k), dfdy(:,:,k))
enddo

end subroutine ddy


!*******************************************************************************
subroutine ddxy(f, dfdx, dfdy)
!*******************************************************************************
!
! First derivative in x and y direction together
! Derivatives and input array on vertical levels 1:nz-1 (see divstress routines)
! Output array has vertical levels 0:nz (see dummy array divstress routines)
!
use param, only: rp,ld,ny,lh,nz,inxny
use fft, only: back,forw,kx,ky
implicit none
real(rp), dimension(1:ld,1:ny,1:nz-1), intent(in) :: f
real(rp), dimension(1:ld,1:ny,0:nz  ), intent(inout) :: dfdx,dfdy
integer :: k

do k=1,nz-1
  ! Go to spectral space
  dfdx(:,:,k)=inxny*f(:,:,k)   
  call dfftw_execute_dft_r2c(forw,dfdx(:,:,k),dfdx(:,:,k))

  ! Calculate derivatives: kx,ky zero out Nyquist frequencies
  call muli(dfdx(1:ld,1:ny,k),ky(1:lh,1:ny),dfdy(1:ld,1:ny,k))
  call muli(dfdx(1:ld,1:ny,k),kx(1:lh,1:ny),dfdx(1:ld,1:ny,k))

  ! Go back to physical space
  call dfftw_execute_dft_c2r(back,dfdx(:,:,k),dfdx(:,:,k))
  call dfftw_execute_dft_c2r(back,dfdy(:,:,k),dfdy(:,:,k))
enddo

end subroutine ddxy


!*******************************************************************************
subroutine filt_da(f, dfdx, dfdy, intflag)
!*******************************************************************************
!
! First derivative in x and y direction together
! Also takes out Nyquist frequency in input array (u,v, and w) (see main)
! Derivatives and input array on vertical levels 0:nz (see main)
! Output array has vertical levels 0:nz (see main)
!
use param, only: rp,ld,ny,nz,lh,coord,nproc,inxny
use fft, only: back,forw,kx,ky
implicit none
real(rp), dimension(1:ld,1:ny,0:nz), intent(inout) :: f
real(rp), dimension(1:ld,1:ny,0:nz), intent(inout) :: dfdx, dfdy
integer :: k,kstart,kend
integer, intent(in) :: intflag

if(coord.eq.0) then
  if(intflag.eq.1) then
    ! For u,v velocties
    kstart=1
  elseif(intflag.eq.2) then
    ! For w velocties (zero velocity at wall)
    f(:,:,1)   =0.0_rp
    dfdx(:,:,1)=0.0_rp
    dfdy(:,:,1)=0.0_rp
    kstart=2
  endif
else
  kstart=0
endif

if(coord.eq.nproc-1) then
  kend=nz-1
else
  kend=nz
endif

! For temperature
if(intflag .eq. 3) then
  kstart = 0
  kend = nz
endif

do k=kstart,kend
  ! Go to spectral space
  f(:,:,k)=inxny*f(:,:,k)
  call dfftw_execute_dft_r2c(forw,f(:,:,k),f(:,:,k))

  ! Take out nyquist frequency 
  f(ld-1:ld,:,k)=0._rp
  f(:,ny/2+1 ,k)=0._rp

  ! Calculate derivatives: kx,ky zero out Nyquist frequencies
  call muli(f(1:ld,1:ny,k),kx(1:lh,1:ny),dfdx(1:ld,1:ny,k))
  call muli(f(1:ld,1:ny,k),ky(1:lh,1:ny),dfdy(1:ld,1:ny,k))

  ! Go back to physical space
  call dfftw_execute_dft_c2r(back, dfdx(:,:,k), dfdx(:,:,k))
  call dfftw_execute_dft_c2r(back, dfdy(:,:,k), dfdy(:,:,k))

  ! Go back to physical space 
  call dfftw_execute_dft_c2r(back,f(:,:,k),f(:,:,k))
enddo

if(coord.eq.nproc-1) then
  if(intflag.eq.1) then
    ! For u,v velocties (no stress at top)
    f(:,:,nz)    = f(:,:,nz-1)
    dfdx(:,:,nz) = dfdx(:,:,nz-1)
    dfdy(:,:,nz) = dfdy(:,:,nz-1)
  elseif(intflag.eq.2) then
    ! For w velocties (zero velocity at top)
    f(:,:,nz)    = 0.0_rp 
    dfdx(:,:,nz) = 0.0_rp
    dfdy(:,:,nz) = 0.0_rp
  endif
endif

end subroutine filt_da


!*******************************************************************************
subroutine ddz_uv(f, dfdz)
!*******************************************************************************
!
! Vertical derivative on uv-grid
! First level is not calculated as tzz(:,:,1)=0.0_rp and dudz(:,:,1) and
! dvdz(:,:,1) are set by wallstress
! At top level a zero velocity derivative is imposed (for u and v, see main)
! Derivatives  calculated on vertical levels 1:nz
! Input/output array have vertical levels 0:nz
!
use param, only: rp,ld,ny,nz,kmin
use grid_defs, only: gridzw 
implicit none
real(rp), dimension(1:ld,1:ny,0:nz), intent (in) :: f
real(rp), dimension(1:ld,1:ny,0:nz), intent (inout) :: dfdz
integer :: k
real(rp) :: idz1

do k=kmin,nz
  idz1 = 1._rp/(0.5_rp * (gridzw(k+1) - gridzw(k-1))) 
  dfdz(:,:,k) = idz1*(f(:,:,k) - f(:,:,k-1))
end do

end subroutine ddz_uv


!*******************************************************************************
subroutine ddz_w(f, dfdz)
!*******************************************************************************
!
! Vertical derivative on w-grid
! Derivatives  calculated on vertical levels 0:nz-1
! Input/output array have vertical levels 0:nz
! 
use param, only: rp,ld,ny,nz
use grid_defs, only: gridzw 
implicit none
real(rp), dimension(1:ld,1:ny,0:nz), intent (in) :: f
real(rp), dimension(1:ld,1:ny,0:nz), intent (inout) :: dfdz
integer :: k
real(rp) :: idz1 

do k=0,nz-1
  idz1 = 1._rp/(gridzw(k+1) - gridzw(k)) 
  dfdz(:,:,k) = idz1*(f(:,:,k+1) - f(:,:,k))
enddo

end subroutine ddz_w
