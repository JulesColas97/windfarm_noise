module utils

implicit none

private

public :: res_trid
public :: write_1Darray
public :: write_2Darray
public :: write_1integer
public :: write_1real
public :: init_file
public :: linspace
public :: stop_error
public :: float_to_str
public :: int_to_str

contains

!! convert float to string
character(len=20) function float_to_str(k)
    real, intent(in) :: k
    write (float_to_str, '(F5.1)') k
    float_to_str = adjustl(float_to_str)
end function float_to_str

!! convert integer to string
character(len=20) function int_to_str(k)
    integer, intent(in) :: k
    write (int_to_str, '(I5.4)') k
    int_to_str = adjustl(int_to_str)
end function int_to_str

!! @brief Tridiagonal matrix inversion using Thomas alogrithm
!!
!! @details this function return the product between phi and a tri diagonal matruix define by (a,b,c)
!!
!!
!! @param [in]      a: lower diagonal
!! @param [in]      b: middle diagonal
!! @param [in]      c: upper diagonal
!! @param[inout]    phi: a vector
!! @param[inout]    dom: dimension of the matrix and vector phi
!! @return  this function returns the product between the tridiagonal matrix and the vector phi
subroutine res_trid(phi,a,b,c,dms)

    complex, intent(in)                         :: a(:)
    complex, intent(in)                         :: b(:)
    complex, intent(in)                         :: c(:)

    complex, intent(inout)                      :: phi(:)
    integer, intent(in)                         :: dms

    complex, allocatable, dimension(:)          :: u,a2,b2,c2
    integer                                     :: kk

    allocate(u(dms))
    allocate(a2(dms-1))
    allocate(c2(dms-1))
    allocate(b2(dms-1))
    u = phi
    a2 = a
    b2 = b
    c2 = c
    a2(dms-1)=c2(dms-1)/b2(dms)
    u(dms)=u(dms)/b2(dms)

    do kk=dms-1,2,-1
        b2(kk)=b2(kk)-c2(kk)*a2(kk)
        a2(kk-1)=a2(kk-1)/b2(kk)
        u(kk)=(u(kk)-c2(kk)*u(kk+1))/b2(kk)
    enddo

    b2(1)=b2(1)-c2(1)*a2(1)
    u(1)=(u(1)-c2(1)*u(2))/b2(1)

    do kk=2,dms
        u(kk)=u(kk)-a2(kk-1)*u(kk-1)
    enddo
    phi=u
end subroutine res_trid

!! @brief Write a integer inside a h5 files
!!
!! @details this function create the data space and the data set  and write an integer inside an h5 file
!!
!!
!! @param [in]      a: the integer to be written
!! @param [in]      fname: the name of the variables (and hence name of the data set)
!! @param [in]      file_id: the file _id which must have be previously open
subroutine write_1integer(a,fname,file_id)

    use hdf5

    character(len=*),intent(in)                             :: fname
    integer,intent(in)                                      :: a
    integer(HID_T),intent(in)                               :: file_id

    integer(HSIZE_T),dimension(1)                           :: dim
    integer(HID_T)                                          :: dspace_id, dset_id
    integer                                                 :: error
    dim(1) = 1

    call h5screate_simple_f(1,dim,dspace_id,error)
    call h5dcreate_f(file_id,fname, H5T_NATIVE_INTEGER, dspace_id, dset_id,error)
    call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,a, dim,error)
    call h5dclose_f(dset_id,error)
    call h5sclose_f(dspace_id,error)

end subroutine write_1integer

!! @brief Write a real inside a h5 files
!!
!! @details this function create the data space and the data set  and write a real double inside an h5 file
!!
!!
!! @param [in]      a: the double to be written
!! @param [in]      fname: the name of the variables (and hence name of the data set)
!! @param [in]      file_id: the file _id which must have be previously open
subroutine write_1real(a,fname,file_id)

    use hdf5

    character(len=*),intent(in)                             :: fname
    real,intent(in)                                         :: a
    integer(HID_T),intent(in)                               :: file_id

    integer(HSIZE_T),dimension(1)                           :: dim
    integer(HID_T)                                          :: dspace_id, dset_id
    integer                                                 :: error
    dim(1) = 1
    call h5screate_simple_f(1,dim,dspace_id,error)
    call h5dcreate_f(file_id,fname, H5T_NATIVE_DOUBLE, dspace_id, dset_id,error)
    call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,a, dim,error)
    call h5dclose_f(dset_id,error)
    call h5sclose_f(dspace_id,error)

end subroutine write_1real

!! @brief Write a 1D double array inside a h5 files
!!
!! @details this function create the data space and the data set  and write a double 1D array inside an h5 file
!!
!!
!! @param [in]      a: the arrayto be written
!! @parma [in]      n: the size of the vector
!! @param [in]      fname: the name of the variable (and hence name of the data set)
!! @param [in]      file_id: the file _id which must have be previously open
subroutine write_1Darray(a,n,fname,file_id)

    use hdf5

    character(len=*),intent(in)                             :: fname
    integer,intent(in)                                      :: n
    real,dimension(:),intent(in)                            :: a
    integer(HID_T),intent(in)                               :: file_id

    integer(HSIZE_T),dimension(1)                           :: dim
    integer(HID_T)                                          :: dspace_id, dset_id
    integer                                                 :: error
    dim(1) = n

    call h5screate_simple_f(1,dim,dspace_id,error)
    call h5dcreate_f(file_id,fname, H5T_NATIVE_DOUBLE, dspace_id, dset_id,error)
    call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,a, dim,error)
    call h5dclose_f(dset_id,error)
    call h5sclose_f(dspace_id,error)

end subroutine write_1Darray

!! @brief Write a 2D double array inside a h5 files
!!
!! @details this function create the data space and the data set  and write a double 2D array inside an h5 file
!!
!!
!! @param [in]      a: the arrayto be written
!! @parma [in]      n: the size of the first dimension of the array
!! @param [in]      m: the size of the second dimension of the array
!! @param [in]      fname: the name of the variable (and hence name of the data set)
!! @param [in]      file_id: the file _id which must have be previously open
subroutine write_2Darray(a,n,m,fname,file_id)

    use hdf5

    character(len=*),intent(in)                            :: fname
    integer,intent(in)                                     :: n,m
    real,dimension(:,:),intent(in)                         :: a
    integer(HID_T),intent(in)                              :: file_id

    integer(HSIZE_T),dimension(2)                           :: dim
    integer(HID_T)                                          :: dspace_id, dset_id
    integer                                                 :: error
    dim(1) = n
    dim(2) = m


    call h5screate_simple_f(2,dim,dspace_id,error)
    call h5dcreate_f(file_id,fname, H5T_NATIVE_DOUBLE, dspace_id, dset_id,error)
    call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,a, dim,error)
    call h5dclose_f(dset_id,error)
    call h5sclose_f(dspace_id,error)

end subroutine write_2Darray


!! @brief Init the H5 file where the solution will be written
!!
!! @details this function initialize the h5 files with the relative path given from the execution folder
!! the solution will be written in a file ./solution/case_name_theta.h5
!! the matrix of p_real and imaginary will be stored in a group '/solution/$frequency$' and the mesh will be stored
!! in a group 'mesh'
!!
!!
!! @param [in]      case_name: the name of the simulation
!! @parma [in]      theta: the angle of propagation
subroutine init_file(case_name,theta, heights)
    use hdf5
    implicit none

    character(len=40), intent(in)                           :: case_name
    real,dimension(:),intent(in)                            :: heights
    real,intent(in)                                         :: theta 
    integer(HID_T)                                          :: file_id, group_id
    integer                                                 :: error,a

    character(len=10)                                       :: theta_str
    character(len=100)                                      :: file_name

    file_name =  trim(case_name)
    if (int(theta)==theta) then 
        theta_str = trim(adjustl(int_to_str(int(theta))))
    else 
        theta_str = trim(adjustl(float_to_str(theta)))
    end if

    file_name=trim('./'//trim(case_name)//'_'//trim(theta_str)//'.h5')
    
    call h5open_f(error)
    ! create new file
    call h5fcreate_f(file_name, H5F_ACC_TRUNC_F,file_id,error)
    call h5gcreate_f(file_id, '/solution', group_id, error)
    CALL h5gclose_f(group_id, error)
    call h5gcreate_f(file_id, '/planes', group_id, error)
    CALL h5gclose_f(group_id, error)
    call h5gcreate_f(file_id, '/receiver', group_id, error)
    call write_1Darray(heights, size(heights) ,'heights',group_id)
    CALL h5gclose_f(group_id, error)
    ! close file
    call h5fclose_f(file_id,error)
    !close fortran interface
    call h5close_f(error)
end subroutine init_file

!! @brief equivalent to the linspace function in Matlab
!!
!!  @param [in]     from: first value pf the vector
!!  @param [in]     to : last value of the vector
!!  @param [in]     n: number of element of the vector
!!  @param [inout]  array : vector of n elements from 'from' to 'to' with equal
!!                           step
subroutine linspace(from, to, n, array)
    real, intent(in)                                    :: from, to
    integer, intent(in)                                 :: n
    real,allocatable,dimension(:), intent(inout)        :: array

    real :: range
    integer :: i

    allocate(array(n))
    range = to - from

    if (n == 0) return

    if (n == 1) then
        array(1) = from
        return
    end if

    do i=1, n
        array(i) = from + range * (i - 1) / (n - 1)
    end do
end subroutine linspace

subroutine stop_error(msg)
! Aborts the program with nonzero exit code
!
! The statement "stop msg" will return 0 exit code when compiled using
! gfortran. stop_error() uses the statement "stop 1" which returns an exit code
! 1 and a print statement to print the message.
!
! Example
! -------
!
! call stop_error("Invalid argument")

character(len=*) :: msg ! Message to print on stdout
print *, msg
stop 1
end subroutine
end module utils
