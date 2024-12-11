      subroutine write_spanw_spec(version,var,varname,folder,nzsize,nz_local)
      use param,only:ny,rp,nz,nz_tot,kmax,jt_total,localComm,mpi_info_null,coord
      use hdf5
      implicit none
      real(rp), dimension(1:ny/2+1,1:nz_local), intent(IN) :: var
      character(len=*), intent(in) :: varname,folder
      integer, intent(in) :: version
      integer, intent(in) :: nzsize
      integer, intent(in) :: nz_local
      integer(HID_T) :: file_id,filespace,slabspace,memspace,plist_id,dset_name
      integer(HSIZE_T),  dimension(2) :: data_count,dims
      integer(HSSIZE_T), dimension(2) :: data_offset 
      integer :: ndims,hdf_error
      character*50 filnam1
      character*8 ipfi
 
!     Form the name of the file
      if (version.eq.1) then
      filnam1=trim(folder) // trim('/') // trim(varname) // '.h5'
      elseif (version.eq.2) then 
      write(ipfi,82)jt_total
   82 format(i8.8)
      filnam1=trim(folder) // trim('/') // trim(varname) // '_' //trim(ipfi) // '.h5'
      endif

!     Set offsets and element counts
      ndims = 2
      dims(1)=ny/2+1
      if (nzsize.eq.1) then
      dims(2)=nz_tot
      elseif (nzsize.eq.0) then
      dims(2)=nz_tot-1
      endif

      call h5screate_simple_f(ndims,dims,filespace,hdf_error)

      data_count(1) = ny/2+1
      data_count(2) = nz_local 
      
      data_offset(1) = 0
      data_offset(2) = coord*(nz-1)

!     Write variable in HDF5 file
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
      call h5dwrite_f(dset_name, H5T_NATIVE_DOUBLE,&
           var(1:ny/2+1,1:kmax),dims,& 
           hdf_error,file_space_id=slabspace,mem_space_id=memspace,xfer_prp=plist_id)

      call h5pclose_f(plist_id ,hdf_error)
      call h5dclose_f(dset_name,hdf_error)
      call h5sclose_f(memspace ,hdf_error)
      call h5fclose_f(file_id  ,hdf_error)

      end subroutine write_spanw_spec
