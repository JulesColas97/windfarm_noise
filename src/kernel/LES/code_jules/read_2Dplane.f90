! subroutine read_2Dplane(var,varname,folder)
! var     : The variable that needs to be read 
! varname : String with the name of the variable (will be added to filename)
! folder  : Folder from which the datafile needs to be read
subroutine read_2Dplane(var,varname,folder)
use param,only:rp,nx,ny,coord,ierr
use mpi, only: mpi_comm_world
use hdf5
implicit none      
real(rp), dimension(1:nx,1:ny), intent(out) :: var
character(len=*), intent(in) :: varname,folder
integer(HID_T) :: file_id,filespace,dset_name
integer(HSIZE_T) :: dims(2)
integer :: ndims,hdf_error
character*50 filnam1
logical :: exst

filnam1=trim(folder) // trim('/') // trim(varname) // '.h5'

inquire (file=filnam1,exist=exst)
if(exst) then
 ndims=2
 dims(1)=nx
 dims(2)=ny
 
 call h5fopen_f(filnam1,H5F_ACC_RDONLY_F,file_id,hdf_error) 
 call h5dopen_f(file_id,'var',dset_name,hdf_error) 
 call h5screate_simple_f(ndims,dims,filespace,hdf_error)
 call h5dread_f(dset_name,H5T_NATIVE_DOUBLE,var(1:nx,1:ny),dims,hdf_error)
 call h5dclose_f(dset_name,hdf_error)
 call h5sclose_f(filespace,hdf_error)
 call h5fclose_f(file_id,hdf_error)

else
 if(coord == 0 ) write(*,*) 'Why! Why! Why would you restart a simulation without the input files?'
 call mpi_barrier(mpi_comm_world,ierr )
 call mpi_finalize (ierr)
 stop
endif


end subroutine read_2Dplane
