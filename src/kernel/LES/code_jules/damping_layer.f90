#ifdef SCALAR
!*******************************************************************************
subroutine damping_layer  
!*******************************************************************************
use param, only: nx,L_x, L_z, z_i, nz, pi, rp, coord, nproc
use param, only: damping_x,damping_method, dmpfr, ra_damp_exponent, damp_hgt 
use grid_defs, only: gridx, gridzw
use scalars_param, only: sponge, sponge_x
implicit none
real(rp) :: z_d,z_local
real(rp) :: x_d, x_local 
integer :: k

if (damping_method==2) then
  !Damping explained in Klemp and Lilly (1978) 
  !J. Atmos. Sciences 
  !Numerical Simulation of Hydrostatic Mountain Waves
  !Also see, Klemp and Dudhia (2008): An Upper Gravity-Wave Absorbing
  !Layer for NWP Applications 
  !Height where damping layer starts in (m) 
  z_d = damp_hgt * z_i  

  do k=1,nz-1
    z_local = gridzw(k)*z_i
    if ((z_local .gt. z_d) .and. (z_local .le. L_z*z_i)) then
      sponge(k) = dmpfr*0.5_rp*(1._rp-cos(pi*(z_local-z_d)/((L_z*z_i)-z_d)))
    else
      sponge(k) = 0._rp
    end if
  end do

  if (coord == nproc-1) then
    sponge(nz) = sponge(nz-1)
  end if

  if(damping_x) then 
     sponge_x = 0.0_rp 
     do k=1,nx 
       x_local = gridx(k)*z_i 
       x_d = 0.1 * L_x * z_i  
       if(x_local .gt. x_d .and. x_local .le. (L_x*z_i)) then  
        sponge_x(k) = dmpfr * 0.5_rp *(1._rp - cos(pi * (x_local - x_d)/((L_x * z_i)-x_d)))
       else 
        sponge_x(k) = 0._rp 
       endif 
     enddo 
  endif 

elseif (damping_method==1) then
  !Rayleigh Damping layer 
  z_d = damp_hgt * z_i  

  do k=1,nz-1
    z_local = gridzw(k)*z_i
    if ((z_local .gt. z_d) .and. (z_local .le. L_z*z_i)) then
      sponge(k) = dmpfr*((z_local-z_d)/((L_z*z_i)-z_d))**ra_damp_exponent
    else
      sponge(k) = 0._rp
    end if
  end do

  if (coord == nproc-1) then
    sponge(nz) = sponge(nz-1)
  end if
end if

end subroutine damping_layer 
#endif  
