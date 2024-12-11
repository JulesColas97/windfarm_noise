!*******************************************************************************
! subroutine write_3dfield(version,var,varname,folder)
!   version : 1 without timestamp; 2 add timestep number to filename 
!   var     : The variable that needs to be printed 
!   varname : String with the name of the variable (will be added to filename)
!   folder  : Folder in which the datafile needs to be saved
!   nzsize  : 1: total array size (nz_tot); 0 total array size (nz_tot-1)
!   nz_local: Define local top boundary in the file
!*******************************************************************************
subroutine write_3Dfield(version,var,varname,folder,nzsize,nz_local)
use param,only:rp,nx,ny,nz,nz_tot,jt_total,localComm,coord
use mpi, only: mpi_info_null
use hdf5
implicit none
real(rp), dimension(1:nx,1:ny,1:nz_local), intent(IN) :: var
character(len=*), intent(in) :: varname,folder
integer, intent(in) :: version
integer, intent(in) :: nzsize
integer, intent(in) :: nz_local
integer(HID_T) :: file_id,filespace,slabspace,memspace,plist_id,dset_name
integer(HSIZE_T),  dimension(3) :: data_count,dims
integer(HSSIZE_T), dimension(3) :: data_offset 
integer :: ndims,hdf_error
character*50 :: filnam1
character*8 :: ipfi
 
! Form the name of the file
if (version.eq.1) then
  filnam1=trim(folder) // trim('/') // trim(varname) // '.h5'
elseif (version.eq.2) then 
  write(ipfi,82)jt_total
82 format(i8.8)
  filnam1=trim(folder) // trim('/') // trim(varname) // '_' //trim(ipfi) // '.h5'
endif

! Set offsets and element counts
ndims = 3
dims(1)=nx 
dims(2)=ny
if (nzsize.eq.1) then
  dims(3)=nz_tot
elseif (nzsize.eq.0) then
  dims(3)=nz_tot-1
endif

call h5screate_simple_f(ndims,dims,filespace,hdf_error)

data_count(1) = nx
data_count(2) = ny
data_count(3) = nz_local 

data_offset(1) = 0
data_offset(2) = 0
data_offset(3) = coord*(nz-1)

! Write variable in HDF5 file
call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,hdf_error)
call h5pset_fapl_mpio_f(plist_id,localComm,mpi_info_null,hdf_error)
call h5fcreate_f(filnam1,H5F_ACC_TRUNC_F,file_id,hdf_error,access_prp=plist_id)
call h5pclose_f(plist_id, hdf_error)

call h5dcreate_f(file_id,'var',H5T_NATIVE_DOUBLE,filespace,dset_name,hdf_error)
call h5screate_simple_f(ndims,data_count,memspace,hdf_error) 
call h5dget_space_f(dset_name,slabspace,hdf_error)
call h5sselect_hyperslab_f(slabspace,H5S_SELECT_SET_F,data_offset,data_count,hdf_error)
call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,hdf_error) 
call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,hdf_error)

call h5dwrite_f(dset_name, H5T_NATIVE_DOUBLE, var(1:nx,1:ny,1:nz_local), dims,& 
           hdf_error,file_space_id=slabspace,mem_space_id=memspace,xfer_prp=plist_id)

call h5pclose_f(plist_id ,hdf_error)
call h5dclose_f(dset_name,hdf_error)
call h5sclose_f(memspace ,hdf_error)
call h5fclose_f(file_id  ,hdf_error)

end subroutine write_3Dfield 
