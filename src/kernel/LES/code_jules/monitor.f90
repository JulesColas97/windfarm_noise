!*******************************************************************************
subroutine monitor(ke_global,rms_global)
!*******************************************************************************
use mpi, only: mpi_sum
use param, only: rp,nx,ny,nz,nz_tot,localComm,ierr,mpi_rp,coord,path,total_time
use sim_param, only: u,v,w,dudx,dvdy,dwdz,ustar
#ifdef SCALAR 
use scalars_param, only: T_s
#endif SCALAR
implicit none
real(rp) :: us,ke,temp_w,ke_global,rms,rms_global
integer :: i,j,k
#ifdef SCALAR
real(rp) :: temp_s
#endif SCALAR

ke=0._rp
rms=0._rp
do k=1,nz-1
do j=1,ny
do i=1,nx
  temp_w = 0.5_rp*(w(i,j,k)+w(i,j,k+1))
  ke = ke + (u(i,j,k)**2+v(i,j,k)**2+temp_w**2)
  rms = rms + abs(dudx(i,j,k)+dvdy(i,j,k)+dwdz(i,j,k))
end do
end do
end do

! Calculate energy in the domain 
ke = ke*0.5_rp/float(nx*ny*(nz_tot-1))
call mpi_reduce(ke,ke_global,1,mpi_rp,mpi_sum,0,localComm,ierr)

! Calculated velocity divergence metrix
rms = rms/float(nx*ny*(nz_tot-1))
call mpi_reduce (rms,rms_global,1,mpi_rp,mpi_sum,0,localComm,ierr)

if (coord.eq.0) then
  us = sum(ustar(1:nx,1:ny))/float(nx*ny)
#ifdef SCALAR 
#ifdef LVLSET
  temp_s = sum(T_s(1:nx,1:ny,1))/float(nx*ny) 
#else
  temp_s = sum(T_s(1:nx,1:ny))/float(nx*ny) 
#endif LVLSET
#endif SCALAR

  open(2,file=trim(path)//'output/check_ke.dat',status='unknown',form='formatted',position='append')
#ifdef SCALAR 
182 format(5(e14.7,1x))
  write(2,182) total_time,ke_global,us,temp_s
#else
182 format(4(e14.7,1x))
  write(2,182) total_time,ke_global,us
#endif SCALAR
  close(2)
endif

end subroutine monitor
