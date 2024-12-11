#ifdef SCALAR
! This subroutine takes the arrays F_{LM,MM,QN,NN} from the previous
! timestep and essentially moves the values around to follow the
! corresponding particles. The (x, y, z) value at the current
! timestep will be the (x-u*lagran_dt, y-v*lagran_dt, z-w*lagran_dt) value at
! the previous timestep. Since particle motion does not conform to
! the grid, an interpolation will be required.
! This subroutine assumes that the Lagrangian CFL will never exceed 1.
subroutine interpolag_Sdep_Scalar
use param,only:rp,nx,ny,nz,ld,coord,nproc,kmin
use sgs_param,only:lagran_dt
use sim_param,only:u,v,w
use sim_param,only:tempF_KX=>dummy1,tempF_XX=>dummy2,tempF_PY=>dummy3,tempF_YY=>dummy4
use sgs_param, only: F_KX, F_XX, F_PY, F_YY 
use grid_defs, only : gridx, gridy, gridz, gridzw
use mpi_defs, only:mpi_sync_down,mpi_sync_up
implicit none
real(rp), dimension(3) :: xyz_past
integer :: i,j,k
real(rp) :: dz1, dz2, dzfact

call mpi_sync_real_array(F_KX,mpi_sync_up)
call mpi_sync_real_array(F_XX,mpi_sync_up)
call mpi_sync_real_array(F_PY,mpi_sync_up)
call mpi_sync_real_array(F_YY,mpi_sync_up)

tempF_KX(1:ld,1:ny,0:nz) = F_KX(1:ld,1:ny,0:nz)
tempF_XX(1:ld,1:ny,0:nz) = F_XX(1:ld,1:ny,0:nz)
tempF_PY(1:ld,1:ny,0:nz) = F_PY(1:ld,1:ny,0:nz)
tempF_YY(1:ld,1:ny,0:nz) = F_YY(1:ld,1:ny,0:nz)

!Interpolation first uv node 
if(coord.eq.0) then
do j=1,ny
  do i=1,nx
    ! Determine position at previous timestep (uv grid)
    xyz_past(1) = gridx(i) - u(i,j,1)*lagran_dt
    xyz_past(2) = gridy(j) - v(i,j,1)*lagran_dt
    xyz_past(3) = gridz(1) - 0.5_rp * w(i,j,2)*lagran_dt

    ! Interpolate: flag_grid:1 indicates uv grid  
    call trilinear_interp_init(xyz_past,1)
    call interpolation(tempF_KX(1:ld,1:ny,0:nz),F_KX(i,j,1))
    call interpolation(tempF_XX(1:ld,1:ny,0:nz),F_XX(i,j,1))
    call interpolation(tempF_PY(1:ld,1:ny,0:nz),F_PY(i,j,1))
    call interpolation(tempF_YY(1:ld,1:ny,0:nz),F_YY(i,j,1))
  enddo 
enddo
endif

!Interior points 
do k=kmin,nz-1
 dz1 = 0.5_rp * (gridzw(k) - gridzw(k-1)) 
 dz2 = 0.5_rp * (gridzw(k+1) - gridzw(k)) 
 dzfact=(1._rp/(dz1+dz2))

 do j=1,ny
  do i=1,nx
    ! Determine position at previous timestep on w-grid
    xyz_past(1) = gridx(i) - dzfact*(u(i,j,k-1)*dz2+u(i,j,k)*dz1)*lagran_dt
    xyz_past(2) = gridy(j) - dzfact*(v(i,j,k-1)*dz2+v(i,j,k)*dz1)*lagran_dt
    xyz_past(3) = gridzw(k) - w(i,j,k)*lagran_dt

    ! Interpolate: flag_grid:2 indicates w-grid
    call trilinear_interp_init(xyz_past,2)
    call interpolation(tempF_KX(1:ld,1:ny,0:nz),F_KX(i,j,k))
    call interpolation(tempF_XX(1:ld,1:ny,0:nz),F_XX(i,j,k))
    call interpolation(tempF_PY(1:ld,1:ny,0:nz),F_PY(i,j,k))
    call interpolation(tempF_YY(1:ld,1:ny,0:nz),F_YY(i,j,k))
  enddo
 enddo
enddo

! Top-most level should not allow negative w
if (coord.eq.nproc-1) then
k = nz
dz1 = 0.5_rp * (gridzw(k) - gridzw(k-1)) 
dz2 = 0.5_rp * (gridzw(k+1) - gridzw(k)) 
dzfact = (1._rp/(dz1+dz2))

 do j=1,ny
   do i=1,nx
    ! Determine position at previous timestep on w-grid
    xyz_past(1) = gridx(i) - dzfact*(u(i,j,k-1)*dz2+u(i,j,k)*dz1)*lagran_dt
    xyz_past(2) = gridy(j) - dzfact*(v(i,j,k-1)*dz2+v(i,j,k)*dz1)*lagran_dt
    xyz_past(3) = gridzw(k) - max(0.0_rp,w(i,j,k))*lagran_dt

    ! Interpolate: flag_grid:2 indicates w-grid  
    call trilinear_interp_init(xyz_past,2)
    call interpolation(tempF_KX(1:ld,1:ny,0:nz),F_KX(i,j,k))
    call interpolation(tempF_XX(1:ld,1:ny,0:nz),F_XX(i,j,k))
    call interpolation(tempF_PY(1:ld,1:ny,0:nz),F_PY(i,j,k))
    call interpolation(tempF_YY(1:ld,1:ny,0:nz),F_YY(i,j,k))
   enddo
 enddo
endif

! Send data to k=nz nodes
call mpi_sync_real_array(F_KX,mpi_sync_down)
call mpi_sync_real_array(F_XX,mpi_sync_down)
call mpi_sync_real_array(F_PY,mpi_sync_down)
call mpi_sync_real_array(F_YY,mpi_sync_down)

end subroutine interpolag_Sdep_Scalar
#endif 
