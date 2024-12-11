!*******************************************************************************
subroutine smagorinsky_model_init
!*******************************************************************************
use param, only: coord,lbc_mom,nz,Co,wall_damp_exp,vonk, rp
use sgs_param, only: Cs_opt2
#ifdef SCALAR 
use param, only: Pr
use sgs_param, only: Ds_opt2
#endif SCALAR
use filters, only: delta
use grid_defs, only: gridz, gridzw
implicit none
integer :: k 
real(rp), dimension(nz) :: l

! Parameters (Co and nn) for wallfunction defined in param.f90
Cs_opt2 = Co**2  ! constant coefficient
#ifdef SCALAR 
!Technically this should be (1/Pr)*Cs2: 1/Pr = 3 => Pr = 0.33
Ds_opt2 = Cs_opt2/Pr 
if(coord.eq.0) print*, 'cs2, ds2',Co**2,Co**2/Pr 
#endif SCALAR

l = delta 

if (lbc_mom == 0) then
  l = delta
else
  ! The variable "l" calculated below is l_sgs/Co where l_sgs is
  ! from JDA eqn(2.30)
  if (coord == 0) then
    l(1) = ( Co**wall_damp_exp * (vonk*gridz(1))**(-wall_damp_exp) +           &
             delta**(-wall_damp_exp) )**(-1._rp/wall_damp_exp)
    do k = 2, nz
      l(k) = ( Co**wall_damp_exp * (vonk*gridzw(k))**(-wall_damp_exp) +        &
               delta**(-wall_damp_exp) )**(-1._rp/wall_damp_exp)
    end do
  else
    do k = 1, nz
      l(k) = ( Co**wall_damp_exp * (vonk*gridzw(k))**(-wall_damp_exp) +        &
               delta**(-wall_damp_exp) )**(-1._rp/wall_damp_exp)
    end do
  end if
end if

do k = 1,nz
  Cs_opt2(:,:,k) = Cs_opt2(:,:,k) * l(k)**2
#ifdef SCALAR 
  Ds_opt2(:,:,k) = Ds_opt2(:,:,k) * l(k)**2
#endif SCALAR
enddo 

#ifdef LVLSET
call level_set_Cs_smag
#endif LVLSET

end subroutine smagorinsky_model_init
