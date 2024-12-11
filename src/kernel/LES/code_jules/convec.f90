!*******************************************************************************
subroutine convec
!*******************************************************************************
!
! Computes the rotation convective term in physical space
!       c = - (u X vort)
! Uses 3/2-rule for dealiasing for, see Canuto 1991 Spectral Methods, chapter 7
!
use param, only: rp,nz,kmin,coord
use sim_param, only: u,v,w,dudy,dudz,dvdx,dvdz,dwdx,dwdy
use grid_defs, only: gridzw 
#ifdef MEMORY
use sim_param, only: dummy1,dummy2,dummy3,dummy4
use param, only: ld_big,ny2,ld,ny
#endif MEMORY
use sim_param, only: RHSx,RHSy,RHSz,dummy=>dummy_plane
use convec_mod, only: u1_big,u2_big,u3_big,vort1_big,vort2_big,vort3_big,cc_big
use fft, only: forw,back,forw_big,back_big
implicit none
integer :: k,km,kp
real(rp) :: dz1, dz2, dzfact

#ifdef MEMORY
deallocate(dummy1,dummy2,dummy3,dummy4)
#ifndef SCALAR
allocate(u1_big   (ld_big,ny2,kmin-1:nz-1)); !u1_big(:,:,:)=0.0_rp 
allocate(u2_big   (ld_big,ny2,kmin-1:nz-1)); !u2_big(:,:,:)=0.0_rp
allocate(u3_big   (ld_big,ny2,kmin:nz  )); !u3_big(:,:,:)=0.0_rp
#endif SCALAR
allocate(vort1_big(ld_big,ny2,kmin:nz  )); !vort1_big(:,:,:)=0.0_rp
allocate(vort2_big(ld_big,ny2,kmin:nz  )); !vort2_big(:,:,:)=0.0_rp
allocate(vort3_big(ld_big,ny2,1:nz-1)); !vort3_big(:,:,:)=0.0_rp
allocate(cc_big   (ld_big,ny2))       ; !cc_big(:,:)=0.0_rp
#endif MEMORY

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! RHSx - Determine u1 and vort1 on extended physical grid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do k=kmin-1,nz-1
  dummy(:,:)=u(:,:,k)

  ! do forward fft to spectral space 
  call dfftw_execute_dft_r2c(forw,dummy(:,:),dummy(:,:))

  ! Back to physical space on padded array
  call padd(u1_big(:,:,k),dummy(:,:))
  call dfftw_execute_dft_c2r(back_big,u1_big(:,:,k),u1_big(:,:,k))
end do

do k=kmin,nz
  dummy(:,:)=(dwdy(:,:,k)-dvdz(:,:,k))

  ! do forward fft to spectral space 
  call dfftw_execute_dft_r2c(forw,dummy(:,:),dummy(:,:))

  ! Back to physical space on padded array
  call padd(vort1_big(:,:,k),dummy(:,:))
  call dfftw_execute_dft_c2r(back_big,vort1_big(:,:,k),vort1_big(:,:,k))
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! RHSy - Determine v and vort2 on extended physical grid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do k=kmin-1,nz-1
  dummy(:,:)=v(:,:,k)

  ! do forward fft to spectral space 
  call dfftw_execute_dft_r2c(forw,dummy(:,:),dummy(:,:))

  ! Back to physical space on padded array
  call padd(u2_big(:,:,k),dummy(:,:))
  call dfftw_execute_dft_c2r(back_big,u2_big(:,:,k),u2_big(:,:,k))
end do

do k=kmin,nz
  dummy(:,:)=(dudz(:,:,k)-dwdx(:,:,k))

  ! do forward fft to spectral space 
  call dfftw_execute_dft_r2c(forw,dummy(:,:),dummy(:,:))

  ! Back to physical space on padded array
  call padd(vort2_big(:,:,k),dummy(:,:))
  call dfftw_execute_dft_c2r(back_big,vort2_big(:,:,k),vort2_big(:,:,k))
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! RHSz - Determine w and vort3 on extended physical grid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do k=kmin,nz
  dummy(:,:)=w(:,:,k)

  ! do forward fft to spectral space 
  call dfftw_execute_dft_r2c(forw,dummy(:,:),dummy(:,:))

  ! Back to physical space on padded array
  call padd(u3_big(:,:,k),dummy(:,:))
  call dfftw_execute_dft_c2r(back_big,u3_big(:,:,k),u3_big(:,:,k))
end do

do k=1,nz-1 
  dummy(:,:)=(dvdx(:,:,k)-dudy(:,:,k))

  ! do forward fft to spectral space 
  call dfftw_execute_dft_r2c(forw,dummy(:,:),dummy(:,:))

  ! Back to physical space
  call padd(vort3_big(:,:,k),dummy(:,:))
  call dfftw_execute_dft_c2r(back_big,vort3_big(:,:,k),vort3_big(:,:,k))
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! RHSx Compute c = - (u X vort) on extended grid and transform back
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! RHSx
do k=1,nz-1
  kp=k+1
  if (coord.eq.0 .and. k.eq.1) then
    cc_big(:,:) = u2_big(:,:,1 )*(-vort3_big(:,:,1 )) + 0.5_rp*                &
                  u3_big(:,:,2 )*( vort2_big(:,:,2 ))
  else
    cc_big(:,:) = u2_big(:,:,k )*(-vort3_big(:,:,k )) + 0.5_rp*                &
                 (u3_big(:,:,kp)*( vort2_big(:,:,kp)) +                        &
                  u3_big(:,:,k )*( vort2_big(:,:,k )))
  endif

  ! To spectral space
  call dfftw_execute_dft_r2c(forw_big, cc_big(:,:),cc_big(:,:))

  ! Back to physical space
  call unpadd(RHSx(:,:,k),cc_big(:,:))
  call dfftw_execute_dft_c2r(back,RHSx(:,:,k),RHSx(:,:,k))
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! RHSy Compute c = - (u X vort) on extended grid and transform back
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do k=1,nz-1
  kp=k+1
  if (coord.eq.0 .and. k.eq.1) then
    cc_big(:,:) = u1_big(:,:,1 )*( vort3_big(:,:,1 )) + 0.5_rp*                &
                  u3_big(:,:,2 )*(-vort1_big(:,:,2 ))
  else
    cc_big(:,:) = u1_big(:,:,k )*( vort3_big(:,:,k )) + 0.5_rp*                &
                 (u3_big(:,:,kp)*(-vort1_big(:,:,kp)) +                        &
                  u3_big(:,:,k )*(-vort1_big(:,:,k )))
  endif

  ! To spectral space
  call dfftw_execute_dft_r2c(forw_big, cc_big(:,:),cc_big(:,:))

  ! Back to physical space
  call unpadd(RHSy(:,:,k),cc_big(:,:))
  call dfftw_execute_dft_c2r(back,RHSy(:,:,k),RHSy(:,:,k))     
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! RHSz Compute c = - (u X vort) on extended grid and transform back
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(coord==0) RHSz(:,:,1)=0.0_rp
do k=kmin,nz-1
  km=k-1
  dz1 = 0.5_rp * (gridzw(k) - gridzw(k-1)) 
  dz2 = 0.5_rp * (gridzw(k+1) - gridzw(k)) 
  dzfact=(1._rp/(dz1+dz2))
  cc_big(:,:) = dzfact *((u1_big(:,:,k)*dz1 + u1_big(:,:,km) * dz2)*(-vort2_big(:,:,k)) + &
                        (u2_big(:,:,k) * dz1 + u2_big(:,:,km) * dz2)*( vort1_big(:,:,k)))

  ! To spectral space
  call dfftw_execute_dft_r2c(forw_big,cc_big(:,:),cc_big(:,:))

  ! Back to physical space
  call unpadd(RHSz(:,:,k),cc_big(:,:))
  call dfftw_execute_dft_c2r(back,RHSz(:,:,k),RHSz(:,:,k))
end do

#ifdef MEMORY
#ifndef SCALAR
deallocate(u1_big,u2_big,u3_big)
#endif SCALAR
deallocate(vort1_big,vort2_big,vort3_big,cc_big)
allocate(dummy1(ld,ny,0:nz)); dummy1(:,:,:)=0.0_rp
allocate(dummy2(ld,ny,0:nz)); dummy2(:,:,:)=0.0_rp
allocate(dummy3(ld,ny,0:nz)); dummy3(:,:,:)=0.0_rp
allocate(dummy4(ld,ny,0:nz)); dummy4(:,:,:)=0.0_rp
#endif MEMORY

end subroutine convec
