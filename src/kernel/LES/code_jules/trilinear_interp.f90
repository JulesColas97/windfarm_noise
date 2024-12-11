!*******************************************************************************
subroutine trilinear_interp_init(xyz,flag_grid)
!*******************************************************************************
!
! This subroutine perform trilinear interpolation for a point that
! exists in the cell with lower dimension (cell index) : istart,jstart,kstart
! for the point xyz
!  
! istart, jstart, kstart are used to determine the cell location on the
! uv-grid; these are defined in output_init
!
! trilinear_interp_init calculates the general grid information
! trilinear_interp      calculates actual interpolation
!
use grid_defs,only:gridx,gridy,gridz,gridzw,awrap_i,awrap_j
use param, only : rp,nz,idx,idy,L_x,L_y
use interp, only : xdiff,ydiff,zdiff,istart,jstart,kstart,istart1,jstart1,kstart1
implicit none
integer,intent(IN):: flag_grid
real(rp), intent(IN), dimension(3) :: xyz
real(rp) :: px,py
integer :: k

! Determine istart,jstart,kstart
px = modulo(xyz(1),L_x)
istart = floor (px * idx) + 1
py = modulo(xyz(2),L_y)
jstart = floor (py * idy) + 1

!flag_grid 1: Interpolates on uv grid 
if(flag_grid.eq.1) then
   do k=0,nz 
    if(xyz(3) .ge. gridz(k) .and. xyz(3) .lt. gridz(k+1)) then 
      kstart = k 
    endif 
   enddo 

   kstart = max(kstart,0) 
   kstart = min(kstart,nz-1) 

   !  Compute zdiff
   zdiff = xyz(3) - gridz(kstart)

!flag_grid 2: Interpolates on w-grid 
elseif(flag_grid .eq. 2) then
   do k=0,nz 
    if(xyz(3) .ge. gridzw(k) .and. xyz(3) .lt. gridzw(k+1)) then 
      kstart = k 
    endif 
   enddo 

   kstart = max(kstart,0) 
   kstart = min(kstart,nz-1) 

   !  Compute zdiff
   zdiff = xyz(3) - gridzw(kstart)
endif

! Set +1 values
istart1 = awrap_i(istart+1) ! Autowrap index
jstart1 = awrap_j(jstart+1) ! Autowrap index
kstart1 = kstart + 1

!  Compute xdiff
xdiff = px - gridx(istart)
!  Compute ydiff
ydiff = py - gridy(jstart)

end subroutine trilinear_interp_init


!*******************************************************************************
subroutine interpolation(var,answer)
!*******************************************************************************
use param, only : rp,ld,ny,nz,idx,idy
use interp, only : xdiff,ydiff,zdiff,istart,jstart,kstart,istart1,jstart1,kstart1
use grid_defs, only: gridzw 

implicit none
real(rp),dimension(1:ld,1:ny,0:nz),intent(IN):: var
real(rp),intent(OUT)::answer
real(rp) :: u1,u2,u3,u4,u5,u6

!  Perform interpolations in x-direction 
u1=var(istart,jstart ,kstart)  + xdiff*(var(istart1,jstart ,kstart) -var(istart,jstart ,kstart)) *idx
u2=var(istart,jstart1,kstart)  + xdiff*(var(istart1,jstart1,kstart) -var(istart,jstart1,kstart)) *idx
u3=var(istart,jstart ,kstart1) + xdiff*(var(istart1,jstart ,kstart1)-var(istart,jstart ,kstart1))*idx
u4=var(istart,jstart1,kstart1) + xdiff*(var(istart1,jstart1,kstart1)-var(istart,jstart1,kstart1))*idx

!  Perform interpolations in y-direction
u5=u1 + (ydiff) * (u2 - u1) * idy
u6=u3 + (ydiff) * (u4 - u3) * idy
!  Perform interpolation in z-direction
answer = u5 + (zdiff) * (u6 - u5) / (gridzw(kstart1) - gridzw(kstart))  

end subroutine interpolation 
