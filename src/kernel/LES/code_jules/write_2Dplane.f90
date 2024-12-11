! subroutine write_2Dplane(version,var,varname,folder)
! version : 1 without timestamp; 2 add timestep number to filename 
! var     : The variable that needs to be printed 
! varname : String with the name of the variable (will be added to filename)
! folder  : Folder in which the datafile needs to be saved
subroutine write_2Dplane(version,var,varname,folder)
use param,only:rp,nx,ny,jt_total,total_time
use hdf5
implicit none      
real(rp), dimension(1:nx,1:ny), intent(IN) :: var
character(len=*), intent(in) :: varname,folder
integer, intent(in) :: version
integer(HID_T) :: file_id,filespace,dset_name
integer(HSIZE_T) :: dims(2),dims2(1)
integer :: ndims,hdf_error
character*50 filnam1
character*8 ipfi

! Form the name of the file
if (version.eq.1) then
filnam1=trim(folder) // trim('/') // trim(varname) // '.h5'
elseif (version.eq.2) then
write(ipfi,82)jt_total
82 format(i8.8)
filnam1=trim(folder) // trim('/') // trim(varname) // '_' //trim(ipfi) //'.h5'
endif

ndims=2
dims(1)=nx
dims(2)=ny

! Sort out MPI definitions and open file
call h5fcreate_f(filnam1,H5F_ACC_TRUNC_F,file_id,hdf_error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Write the 2D data in the file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Create dataspace
call h5screate_simple_f(ndims,dims,filespace,hdf_error)
! Create the dataset with default properties.
call h5dcreate_f(file_id,'var',H5T_NATIVE_DOUBLE,filespace,dset_name,hdf_error)
! Select hyperslab and then write it
call h5dwrite_f(dset_name,H5T_NATIVE_DOUBLE,var(1:nx,1:ny),dims,hdf_error)
! Close properties and file
call h5dclose_f(dset_name,hdf_error)
call h5sclose_f(filespace,hdf_error)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Write timestamp in the file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Create dataspace
call h5screate_simple_f(1,dims2,filespace,hdf_error)
! Create the dataset with default properties.
call h5dcreate_f(file_id,'time', H5T_NATIVE_DOUBLE,filespace,dset_name,hdf_error)
! Select hyperslab and then write it
call h5dwrite_f(dset_name,H5T_NATIVE_DOUBLE,total_time,dims2,hdf_error)
! Close properties and file
call h5dclose_f(dset_name,hdf_error)
call h5sclose_f(filespace,hdf_error)

! Close the file
call h5fclose_f(file_id,hdf_error)

end subroutine write_2Dplane 
