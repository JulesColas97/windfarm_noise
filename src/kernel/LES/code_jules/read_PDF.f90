!*******************************************************************************
! subroutine read_PDF(var,varname,folder)
!   var     : The variable that needs to be read 
!   varname : String with the name of the variable to selct the appropriate file
!   folder  : Folder from which the file should be read
!   nzsize  : 1: total array size (nz_tot); 0 total array size (nz_tot-1)
!   nz_local: Define local top boundary in the file
!*******************************************************************************
subroutine read_PDF(var,varname,folder,nzsize,nz_local)
use param,only:nbins,rp,lprec,nz,nz_tot,kmax,localComm,coord,ierr
use mpi, only: mpi_info_null, mpi_comm_world
use hdf5
implicit none
real(lprec), dimension(1:nbins,1:nz_local), intent(out) :: var
real(rp), dimension(1:nbins,1:nz_local) :: var_local
character(len=*), intent(in) :: varname,folder
integer, intent(in) :: nzsize
integer, intent(in) :: nz_local
integer(HID_T) :: file_id,slabspace,memspace,plist_id,dset_name
integer(HSIZE_T),  dimension(2) :: data_count,dims
integer(HSSIZE_T), dimension(2) :: data_offset 
integer :: hdf_error,ndims
character*50 :: filnam1
logical :: exst
 
! Form the name of the file
filnam1=trim(folder) // trim('/') // trim(varname) // '.h5'

inquire (file=filnam1,exist=exst)
if(exst) then 
   ! Set offsets and element counts
   ndims = 2
   dims(1)=nbins
   if (nzsize.eq.1) then
     dims(2)=nz_tot
   elseif (nzsize.eq.0) then
     dims(2)=nz_tot-1
   endif
   
   data_count(1) = nbins
   data_count(2) = nz_local 
            
   data_offset(1) = 0
   data_offset(2) = coord*(nz-1)
   
   ! Read variable in HDF5 file
   call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,hdf_error)
   call h5pset_fapl_mpio_f(plist_id,localComm,mpi_info_null,hdf_error)
   call h5fopen_f(filnam1,H5F_ACC_RDONLY_F,file_id,hdf_error,access_prp=plist_id)
   call h5pclose_f(plist_id, hdf_error)
   
   call h5dopen_f(file_id,'var',dset_name,hdf_error)
   call h5screate_simple_f(ndims, data_count, memspace, hdf_error) 
   call h5dget_space_f(dset_name,slabspace,hdf_error)
   call h5sselect_hyperslab_f(slabspace,H5S_SELECT_SET_F,data_offset,data_count,hdf_error)
   call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,hdf_error) 
   call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,hdf_error)
   
   call h5dread_f(dset_name,H5T_NATIVE_DOUBLE, var_local(1:nbins,1:kmax), &
        dims,hdf_error,file_space_id=slabspace,mem_space_id=memspace,xfer_prp=plist_id)
   
   call h5pclose_f(plist_id ,hdf_error)
   call h5dclose_f(dset_name,hdf_error)
   call h5sclose_f(memspace ,hdf_error)
   call h5fclose_f(file_id  ,hdf_error)

   var=var_local

else
 if(coord == 0 ) write(*,*) 'Why would you restart a simulation without the PDF input files?' 
 call mpi_barrier(mpi_comm_world,ierr )
 call mpi_finalize (ierr)
 stop
endif 
end subroutine read_PDF 
