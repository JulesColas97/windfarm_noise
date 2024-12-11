!*******************************************************************************
subroutine initialize_mpi
!*******************************************************************************
!
! Initializes mpi communicators
!
use param, only: ierr,up,down,mpi_rp,mpi_cp,nproc,localComm
use param, only: rank,coord,global_rank
use param, only: rank_of_coord
use mpi, only: mpi_double_precision,mpi_double_complex, mpi_comm_world
implicit none
integer :: ip,coords(1)

call mpi_init (ierr)

! Set the local communicator
#ifdef CPS
! Create the local communicator (split from mpi_comm_world)
! This also sets the globally defined intercommunicator (bridge)
call create_mpi_comms_cps
#else
localComm = mpi_comm_world
#endif CPS

call mpi_comm_size (localComm, nproc, ierr)
call mpi_comm_rank (localComm, global_rank, ierr)

!--set up a 1d cartesian topology 
call mpi_cart_create (localComm, 1, (/ nproc /), (/ .false. /), .false., localComm, ierr)

!  u-node info needs to be shifted up to proc w/ rank "up",
!  w-node info needs to be shifted down to proc w/ rank "down"
call mpi_cart_shift  (localComm, 0, 1, down, up, ierr)
call mpi_comm_rank   (localComm, rank, ierr)
call mpi_cart_coords (localComm, rank, 1, coords, ierr)

coord = coords(1)  !--use coord (NOT rank) to determine global position

!--rank->coord conversions
allocate( rank_of_coord(0:nproc-1))
do ip = 0, nproc-1
  call mpi_cart_rank (localComm, (/ ip /), rank_of_coord(ip), ierr)
end do

mpi_rp = mpi_double_precision
mpi_cp = mpi_double_complex

end subroutine initialize_mpi


!*******************************************************************************
subroutine mpi_sync_real_array(var,isync)
!*******************************************************************************
!
! This subroutine syncing arrays in vertical direction
! sync_down : Send data from k = 1 at coord+1 to k=nz at coord
! sync_up   : Send data from k = nz-1 at coord to k=0 at coord+1
!
! Arguments:
!
! var   : three dimensional array to be sync'd accross processors 
! isync : flag for determining the type of synchronization and takes on values, 
!         mpi_sync_down,mpi_sync_up,mpi_sync_downup
!
use param, only : ld,ny,rp,nz,mpi_rp,down,up,status,ierr,nz,localComm
use mpi_defs,only : mpi_sync_down,mpi_sync_up,mpi_sync_downup
implicit none
real(rp), dimension(1:ld,1:ny,0:nz), intent(INOUT) :: var
integer, intent(in) :: isync
integer:: splane 

splane = ld*ny

if(isync == mpi_sync_down) then
call mpi_sendrecv(var(:,:,1)   ,splane,mpi_rp,down,1,&
                  var(:,:,nz)  ,splane,mpi_rp,up  ,1,&
     localComm, status, ierr)

elseif( isync == mpi_sync_up) then
call mpi_sendrecv(var(:,:,nz-1),splane,mpi_rp,up  ,2,&
                  var(:,:,0)   ,splane,mpi_rp,down,2,&
     localComm, status, ierr)

elseif( isync == mpi_sync_downup) then
call mpi_sendrecv(var(:,:,1)   ,splane,mpi_rp,down,1,&
                  var(:,:,nz)  ,splane,mpi_rp,up  ,1,&
     localComm, status, ierr)
call mpi_sendrecv(var(:,:,nz-1),splane,mpi_rp,up  ,2,&
                  var(:,:,0)   ,splane,mpi_rp,down,2,&
     localComm, status, ierr)
endif

end subroutine mpi_sync_real_array
