!*******************************************************************************
subroutine grid_build
!*******************************************************************************
!
! This subroutine creates the grid information for the domain
!
use grid_defs, only: gridx,gridy,gridz,gridzw,awrap_i,awrap_j,z_tot
use param, only: path,mpi_rp,localComm,ierr,rp,nx,ny,nz,z_i,kmin,kmin2,kmax,kmax2,dx,dy,dz,coord,nproc,nz_tot,grid_stretch 
implicit none
integer :: i,j,k,k_global 
character(len=90) :: filnam,folder 
real(rp), dimension(-1:nz_tot+1) :: read_grid 

! x and y go to nx+1, ny+1 respectively for adding
! the buffered points for periodicity
do j=1,ny+1
  gridy(j) = (j-1)*dy
enddo

do i=1,nx+1
  gridx(i) = (i-1)*dx
enddo

!Create uniform grid (default) 
do k=-1,nz+1 
   k_global = (coord*(nz-1) + k) 
   !Uniform grid points 
   gridzw(k) = (k_global - 1) * dz 
enddo 

!Read stretched grid from file 
if(grid_stretch) then 
   !Read grid file on the master process and distribute it to all processes  
   if(coord==0) then 
     filnam = trim(path) // trim('/') // 'grid.dat'
     open(1, file=filnam, form='formatted', action='read') 
  
     !On multi-processor nz_tot is one more than the number of cells 
     if(nproc > 1) then 
       do i=1,nz_tot 
         read(1,*) read_grid(i) 
       enddo 
     else
     !On single processor full grid needs to be read 
       do i=1,nz_tot+1 
         read(1,*) read_grid(i) 
       enddo 
     endif 

     close(1) 
     
     !Create a grid point above the domain 
     if(nproc > 1) then 
       read_grid(nz_tot+1) = read_grid(nz_tot) + (read_grid(nz_tot) - read_grid(nz_tot-1))
     endif 
     read_grid(0)  = -read_grid(2)  
     read_grid(-1) = -read_grid(3)
   endif

   !Broadcast all the grid points to processors 
   call mpi_bcast(read_grid,nz_tot+3,mpi_rp,0,localComm,ierr)

   !Store the grid values in gridzw array 
   do k=-1,nz+1
      k_global = (coord*(nz-1) + k) 
      gridzw(k) = read_grid(k_global) 
   enddo 

   !Non-dimensionalize the gridzw array 
   gridzw = gridzw/z_i 
endif     

do k=0,nz 
   gridz(k) = gridzw(k) + 0.5_rp * (gridzw(k+1) - gridzw(k))
enddo 
gridz(nz+1) = gridz(nz) + gridz(nz) - gridz(nz-1)

!Print the grids 
folder = trim(path)//'output' 
call write_1Dfield(1,gridzw(1:nz-1),'gridzw',folder,0,nz-1)
call write_1Dfield(1,gridz(1:nz-1),'gridz',folder,0,nz-1) 

!Srinidhi: This is now used by wind-break and turbine modules. I am not sure how gridstretching will affect this 
do k=1,nz_tot
  z_tot(k) = (k - 0.5_rp) * dz
enddo

! Set index autowrapping arrays
awrap_i(0) = nx
do i=1,nx
  awrap_i(i)=i
enddo
awrap_i(nx+1) = 1

awrap_j(0) = ny
do j=1,ny
  awrap_j(j)=j
enddo
awrap_j(ny+1) = 1
     
! Set kmax
if (coord == 0) then
  kmin = 2
  kmin2= 1 ! Used in tridag
  kmax = nz-1
  kmax2= nz ! Used in tridag
elseif (coord == nproc-1) then
  kmin = 1
  kmin2= 2 ! Used in tridag
  kmax = nz
  kmax2= nz+1 ! Used in tridag
else
  kmin = 1
  kmin2= 2 ! Used in tridag
  kmax = nz-1 
  kmax2= nz ! Used in tridag
endif

if (nproc.eq.1) then
  kmin = 2
  kmin2= 1 ! Used in tridag
  kmax = nz
  kmax2= nz+1 ! Used in tridag
endif

end subroutine grid_build
