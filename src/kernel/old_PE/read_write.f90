module  read_write

implicit none

public :: read_param
public :: write_solution
public :: write_2Dflow
public :: read_2Dflow
public ::  read_h5flow3D

contains
! !! convert float to string
! character(len=20) function float_to_str(k)
!     real, intent(in) :: k
!     write (str, '(F4.1)') k
!     str = adjustl(str)
! end function float_to_str

! !! convert integer to string
! character(len=20) function int_to_str(k)
!     integer, intent(in) :: k
!     write (str, '(I5.4)') k
!     str = adjustl(str)
! end function int_to_str
!! @brief Read the inpu parameter file for the simulation by uysing namelist
!!
!! @details this function read the input2.dat files that contaions the parameter for the simulation
!! the input file must be in the same folder that the executable, all vaiables are define like
!! in fortran code (you can use ! to comment line and x=12.5 to define variables )
!!
!!
!! @param [inout]     case_name: the name of the simulation
!! @param [inout]     var0: the medium general parameters
!! @param [inout]     simu: the global 3D parameters for the case
!! @param [inout]     src: position of the source
!! @param [inout]     imp: parameters for the ground impedance model
!! @param [inout]     output: parameters for the solution output
!! @param [inout]     pml: Perfectly matched layer boundary condition parameters
!! @param [inout]     flow: mean 3D (or 2D) flow parameters
subroutine read_param(case_name,var0,simu,src,imp,output,pml,flow,tinput)
    use types

    implicit none
    ! input variables
    !--------------------------------------------------------------------

    character(len=40), intent(inout)                        :: case_name
    type(variables_0), intent(inout)                        :: var0
    type(simulation), intent(inout)                         :: simu
    type(source),intent(inout)                              :: src
    type(impedance),intent(inout)                           :: imp
    type(output_parameters), intent(inout)                  :: output
    type(perfectly_matched_layer),intent(inout)             :: pml
    type(meanflow3D),intent(inout)                          :: flow
    type(twente_input),intent(inout)                        :: tinput

    !local variables
    !--------------------------------------------------------------------
    real                                                    :: ratio
    real                                                    :: pi = 3.1415
    complex                                                 :: eye = (0,1)
    integer                                                 :: ii

    !temporary simu
    integer                             :: nb_theta=1000
    integer                             :: nb_freq=1000
    real,allocatable,dimension(:)        :: theta
    real, allocatable,dimension(:)      :: frequencies

    real                                :: Lx1,Lx2,Ly1,Ly2,Lz
    real                                :: dx, cfl

    logical                             :: external_flow
    logical                             :: interpolation
    logical                             :: uniform
    logical                             :: logarithmic
    real                                :: u0
    logical                             :: arbitrary,arbitrary_new
    logical                             :: continuation ! not implemented yet, but implemented for GTPE

    !temporary output
    real                                :: dout
    integer                             :: nbout
    integer                             :: nb_receiver=100
    real, allocatable,dimension(:)      :: heights
    real, allocatable,dimension(:)      :: xplane
    real, allocatable,dimension(:)      :: yplane
    INTEGER                             :: nb_xplane
    INTEGER                             :: nb_yplane
    logical                             :: side
    logical                             :: top

    ! temporary pml
    integer                             :: size
    real                                :: param
    real                                :: n

    !temporary flow
    character(len=100)                     :: fname, fdir
    real                                   :: XS, YS,ZS

    ! defin element to be read from input.nml
    namelist/input/nb_freq, nb_theta, nb_receiver ,theta, frequencies, Lx1,Lx2,Ly1,Ly2,Lz,dx,cfl,external_flow,interpolation,&
                    uniform, logarithmic,u0, arbitrary,arbitrary_new,continuation,dout, nbout, nb_receiver, heights,xplane,yplane,nb_xplane,nb_yplane, side, top, &
                  case_name,var0, tinput, size,param,n, imp, src, fname, fdir

    allocate(theta(nb_theta))
    allocate(frequencies(nb_freq))
    allocate(heights(nb_receiver))
    allocate(xplane(10))
    allocate(yplane(10))
    nb_xplane = 0
    nb_yplane = 0

    arbitrary = .false.
    arbitrary_new = .false.
    uniform = .false.
    external_flow = .false.
    logarithmic = .false.
    
    ! read inout file
    open(unit=1,file= 'input.nml')
    read(unit=1,nml=input)
    close(1)

    ! because we are using define type some of elements must be read and the assign to the structure
    simu%nb_theta = nb_theta
    simu%nb_freq = nb_freq
    allocate(simu%theta(simu%nb_theta))
    allocate(simu%frequencies(simu%nb_freq))

    simu%theta = theta(1:simu%nb_theta)
    simu%frequencies= frequencies(1:simu%nb_freq)

    simu%Lx1 = Lx1
    simu%Lx2 = Lx2
    simu%Ly1 = Ly1
    simu%Ly2 = Ly2
    simu%Lz = Lz
    simu%dx = dx
    simu%cfl = cfl
    simu%external_flow = external_flow
    simu%interpolation = interpolation
    simu%uniform       = uniform
    simu%logarithmic   = logarithmic
    simu%u0            = u0
    simu%arbitrary     = arbitrary
    simu%arbitrary_new = arbitrary_new
    simu%continuation  = continuation 

    pml%size = size
    pml%param = param
    pml%n = n

    flow%fname = fname
    flow%fdir = fdir

    output%dout = dout
    output%nb_receiver = nb_receiver
    allocate(output%heights(output%nb_receiver))
    output%heights = heights(1:output%nb_receiver)
    output%side = side
    output%top = top
    ! if (nb_xplane>0) then
    !     allocate(output%xplane(nb_xplane))
    !     output%xplane = xplane(1:nb_xplane)
    ! endif
    ! if (nb_yplane>0) then 
    !     allocate(output%yplane(nb_yplane))
    !     output%yplane = yplane(1:nb_yplane)
    ! endif 
    allocate(output%xplane(nb_xplane))
    output%xplane = xplane(1:nb_xplane)    
    allocate(output%yplane(nb_yplane))
    output%yplane = yplane(1:nb_yplane)

end subroutine read_param


!! @briefWrite the solution cartographies  in h5 file
subroutine write_solution(case_name,theta,freq,sol)
    use hdf5
    use utils
    use types
    implicit none

    ! input variables
    character(len=40), intent(in)                           :: case_name
    type(solution), INTENT(IN)                              :: sol
    real,intent(in)                                         :: theta, freq

    ! HDF5 variables
    integer(HID_T)                                          :: file_id, group_id
    integer                                                 :: error,a
    character(len=50)                                       :: file_name, group_name
    character(len=10)                                       :: theta_str, freq_str
    integer                                                 :: nx, ny

    nx = size(sol%deltaL_int,1)
    ny = size(sol%deltaL_int,2)


    if (int(theta)==theta) then 
        theta_str = trim(adjustl(int_to_str(int(theta))))
    else 
        theta_str = trim(adjustl(float_to_str(theta)))
    end if
    file_name =  trim(case_name)
    file_name=trim('./'//trim(case_name)//'_'//trim(theta_str)//'.h5')

    !file=trim(adjustl(fname))//'_sol_.h5'
    ! Initialize fortran interface
    call h5open_f(error)
    ! create new file
    !call h5fcreate_f(fname, H5F_ACC_TRUNC_F,file_id,error)
    call h5fopen_f (file_name, H5F_ACC_RDWR_F, file_id, error)

    ! create a group (with frequency as name)
    a = int(freq)
    freq_str = adjustl(int_to_str(a))
    group_name = trim('/solution/'//freq_str)
    call h5gcreate_f(file_id, group_name, group_id, error)

    ! write p imag and real
  !  call write_2Darray(real(p_vect), ceiling((1.*nx)/nbout), ceiling((1.*ny)/nbout),'p_real',group_id)
  !  call write_2Darray(aimag(p_vect), ceiling((1.*nx)/nbout), ceiling((1.*ny)/nbout),'p_imag',group_id)
    call write_2Darray(sol%deltaL_int, nx, ny,'deltaL',group_id)


    ! close group
    call h5gclose_f(group_id, error)

    ! close file
    call h5fclose_f(file_id,error)
    !close fortran interface
    call h5close_f(error)
end subroutine write_solution



!! @briefWrite the solution cartographies  in h5 file
subroutine write_planes(case_name,theta,freq,sol,output)
    use hdf5
    use utils
    use types
    implicit none

    ! input variables
    character(len=40), intent(in)                           :: case_name
    type(solution), INTENT(IN)                              :: sol
    type(output_parameters),INTENT(IN)                      :: output
    real,intent(in)                                         :: theta, freq

    ! HDF5 variables
    integer(HID_T)                                          :: file_id, group_id
    integer                                                 :: error,a
    character(len=50)                                       :: file_name, group_name
    character(len=10)                                       :: theta_str, freq_str
    integer                                                 :: nx, ny,nz

    nx = size(sol%xplane,1)
    ny = size(sol%yplane,1)
    nz = size(sol%yplane,2)


    if (int(theta)==theta) then 
        theta_str = trim(adjustl(int_to_str(int(theta))))
    else 
        theta_str = trim(adjustl(float_to_str(theta)))
    end if
    file_name =  trim(case_name)
    file_name=trim('./'//trim(case_name)//'_'//trim(theta_str)//'.h5')

    !file=trim(adjustl(fname))//'_sol_.h5'
    ! Initialize fortran interface
    call h5open_f(error)
    ! create new file
    !call h5fcreate_f(fname, H5F_ACC_TRUNC_F,file_id,error)
    call h5fopen_f (file_name, H5F_ACC_RDWR_F, file_id, error)

    ! create a group (with frequency as name)
    a = int(freq)
    freq_str = adjustl(int_to_str(a))
    group_name = trim('/planes/'//freq_str)
    call h5gcreate_f(file_id, group_name, group_id, error)

    ! write p imag and real
  !  call write_2Darray(real(p_vect), ceiling((1.*nx)/nbout), ceiling((1.*ny)/nbout),'p_real',group_id)
  !  call write_2Darray(aimag(p_vect), ceiling((1.*nx)/nbout), ceiling((1.*ny)/nbout),'p_imag',group_id)
    call write_2Darray(sol%xplane, nx, nz,'xplane',group_id)
    call write_1Darray(sol%ycoord,nx,'xycoord',group_id)
    call write_1Darray(output%xplane,nx,'xxcoord',group_id)
    call write_1Darray(sol%xcount,nx,'xcount',group_id)

    call write_2Darray(sol%yplane, ny, nz,'yplane',group_id)
    call write_1Darray(sol%xcoord,ny,'yxcoord',group_id)
    call write_1Darray(output%yplane,ny,'yycoord',group_id)
    call write_1Darray(sol%ycount,ny,'ycount',group_id)

    ! close group
    call h5gclose_f(group_id, error)

    ! close file
    call h5fclose_f(file_id,error)
    !close fortran interface
    call h5close_f(error)
end subroutine write_planes



!! @brief Write receiver
subroutine write_receiver(case_name,theta,freq, nx,heights, receiver)
    use hdf5
    use utils
    implicit none

    ! input variables
    character(len=40), intent(in)                           :: case_name
    integer, intent(in)                                     :: nx
    real,dimension(:,:), intent(in)                          :: receiver
    real,intent(in)                                         :: theta, freq
    real,dimension(:),intent(in)                            :: heights

    ! HDF5 variables
    integer(HID_T)                                          :: file_id, group_id
    integer                                                 :: error,a,ii
    character(len=50)                                       :: file_name, group_name
    character(len=10)                                       :: theta_str, freq_str, height_str


    ! compute file name
    !-----------------------------------------------------------------------
    file_name =  trim(case_name)
    ! theta_str = adjustl(float_to_str(theta))
    if (int(theta)==theta) then 
        theta_str = trim(adjustl(int_to_str(int(theta))))
    else 
        theta_str = trim(adjustl(float_to_str(theta)))
    end if
    file_name=trim('./'//trim(case_name)//'_'//trim(theta_str)//'.h5')

    ! Initialize fortran interface
    call h5open_f(error)
    ! create new file
    call h5fopen_f (file_name, H5F_ACC_RDWR_F, file_id, error)

    ! create a group (with frequency as name)
    !-----------------------------------------------------------------------
    ! a = int(freq)
    ! freq_str = adjustl(int_to_str(a))
    if (int(freq)==freq) then 
        freq_str = trim(adjustl(int_to_str(int(freq))))
    else 
        freq_str = trim(adjustl(float_to_str(freq)))
    end if
    group_name = trim('/receiver/'//freq_str)

    call h5gcreate_f(file_id, group_name, group_id, error)

    call write_2Darray(receiver, nx, size(heights),'deltaL',group_id)
    call h5gclose_f(group_id, error)
    ! close file
    call h5fclose_f(file_id,error)
    !close fortran interface
    call h5close_f(error)
    print*, 'done.'
end subroutine write_receiver


subroutine write_mesh(case_name,theta,x, y)
    use hdf5
    use utils
    implicit none

    ! input variables
    character(len=40), intent(in)                           :: case_name
    real,dimension(:), intent(in)                           :: x, y
    real,intent(in)                                         :: theta

    ! HDF5 variables
    integer(HID_T)                                          :: file_id, group_id
    integer                                                 :: error
    character(len=50)                                       :: file_name
    character(len=10)                                       :: theta_str, freq_str



    ! theta_str = adjustl(float_to_str(theta))
    if (int(theta)==theta) then 
        theta_str = trim(adjustl(int_to_str(int(theta))))
    else 
        theta_str = trim(adjustl(float_to_str(theta)))
    end if
    file_name =  trim(case_name)
    file_name=trim('./'//trim(case_name)//'_'//trim(theta_str)//'.h5')
    

    ! Initialize fortran interface
    call h5open_f(error)
    ! create new file
    !call h5fcreate_f(fname, H5F_ACC_TRUNC_F,file_id,error)
    call h5fopen_f (file_name, H5F_ACC_RDWR_F, file_id, error)
    ! create group mesh
    call h5gcreate_f(file_id, 'mesh', group_id, error)

    !write mesh
    call  write_1Darray(x,size(x),'x',group_id)
    call  write_1Darray(y,size(y),'y',group_id)
    call  write_1real(theta,'tau',group_id)

    ! close group
    call h5gclose_f(group_id, error)

    ! close file
    call h5fclose_f(file_id,error)
    !close fortran interface
    call h5close_f(error)
end subroutine write_mesh

subroutine write_2Dflow(u,v, c, x, z)
    use hdf5
    use utils
    use types
    implicit none

    ! input variables
    real,dimension(:,:), intent(in)                           :: u, v, c
    real,dimension(:,:), intent(in)                           :: x, z

    ! HDF5 variables
    integer(HID_T)                                          :: file_id
    integer                                                 :: error

    ! Initialize fortran interface
    call h5open_f(error)
    ! create new file
    call h5fcreate_f(trim('./'//'flow.h5'), H5F_ACC_TRUNC_F,file_id,error)

    ! write mesh
    call  write_2Darray(x,size(x,1),size(x,2),'x',file_id)
    call  write_2Darray(z,size(z,1),size(z,2),'z',file_id)

    ! write flow
    print*, 'TEST TEST '
    print*, size(c)
    call write_2Darray(u, size(u,1), size(u,2), 'u',file_id)
    call write_2Darray(v, size(v,1), size(v,2), 'v',file_id)
    call write_2Darray(c, size(c,1), size(c,2), 'c',file_id)
    ! call write_2Darray(flow2D%tau, ceiling((1.*nx)/nbout), nz,'tau',file_id)
    ! call write_2Darray(flow2D%epsilon, ceiling((1.*nx)/nbout), nz,'epsilon',file_id)
    ! call write_2Darray(flow2D%gamma, ceiling((1.*nx)/nbout), nz,'gamma',file_id)

    ! close file
    call h5fclose_f(file_id,error)
    !close fortran interface
    call h5close_f(error)
end subroutine write_2Dflow

subroutine read_2Dflow(theta, u,v, epsilon, x, z, nbout,nx, nz, dx, dz)
    use hdf5
    use utils
    implicit none

    ! input variables
    real,intent(in)                                         :: theta
    integer, intent(inout)                                  :: nx, nz, nbout
    real,allocatable,dimension(:,:), intent(inout)          :: u, v, epsilon
    real,allocatable,dimension(:), intent(inout)            :: x, z
    real, intent(inout)                                     :: dx, dz

    ! HDF5 variables
    integer(HID_T)                                          :: file_id, dset_id
    integer                                                 :: error,a
    character(len=10)                                       :: theta_str
    integer(HSIZE_T), dimension(1)                          :: dims1
    integer(HSIZE_T), dimension(2)                          :: dims2

    ! convert theta to str
    ! a=int(theta)
    ! theta_str = trim(adjustl(str(a)))

    call h5open_f(error)
    call h5fopen_f (trim('./'//'flow.h5'), H5F_ACC_RDWR_F, file_id, error)

    ! read parameters
    dims1(1) = 1
    call h5dopen_f(file_id, 'nx', dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, nx, dims1, error)
    call h5dclose_f(dset_id, error)

    call h5dopen_f(file_id, 'nz', dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, nz, dims1, error)
    call h5dclose_f(dset_id, error)

    call h5dopen_f(file_id, 'nbout', dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, nbout, dims1, error)
    call h5dclose_f(dset_id, error)

    call h5dopen_f(file_id, 'dx', dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, dx, dims1, error)
    call h5dclose_f(dset_id, error)

    call h5dopen_f(file_id, 'dz', dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, dz, dims1, error)
    call h5dclose_f(dset_id, error)


    ! read mesh
    dims1(1) = nx
    if (.not. allocated(x)) allocate(x(nx))
    call h5dopen_f(file_id, 'x', dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, x, dims1, error)
    call h5dclose_f(dset_id, error)

    dims1(1) = nz
    if (.not. allocated(z)) allocate(z(nz))
    call h5dopen_f(file_id, 'z', dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, z, dims1, error)
    call h5dclose_f(dset_id, error)

    ! read flow
    dims2(1) = ceiling((1.*nx)/nbout)
    dims2(2) = nz
    if (.not. allocated(u))  allocate(u(ceiling((1.*nx)/nbout),nz))
    call h5dopen_f(file_id, 'u', dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, u, dims2, error)
    call h5dclose_f(dset_id, error)

    if (.not. allocated(v))  allocate(v(ceiling((1.*nx)/nbout),nz))
    call h5dopen_f(file_id, 'v', dset_id, error)
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, v, dims2, error)
    call h5dclose_f(dset_id, error)

    ! Close the file
    call h5fclose_f(file_id, error)
    ! Close FORTRAN interface.
    call h5close_f(error)
end subroutine read_2Dflow


subroutine read_h5flow3D(dirname,fname,flow3D, tinput, c0) !,ratio,Lx, Ly, Lz,z_i, delta, flow3D)
  use hdf5
  use utils
  use types
  !implicit none

  ! input variables
  character(len=*), intent(in)                            :: fname, dirname
  type(meanflow3D), intent(inout)                         :: flow3D
  type(twente_input), intent(in)                          :: tinput
  real, intent(in)                                        :: c0
  !local variables
  real                                                    :: Lx, Ly, Lz
  integer                                                 :: nx, ny, nz, nz_tot
  ! HDF5 variables
  integer(HID_T)                                          :: file_id, dset_id,dspace_id
  integer                                                 :: error,a
  character(len=10)                                       :: theta_str
  integer(HSIZE_T), dimension(3)                          :: dims_u, dimz_v, dimz_w, maxdims
  real, allocatable, dimension(:,:,:)                     :: u_temp
  logical                                                 :: file_exists

  ! read of the 3 component of the flow
  !-------------------------------------------------------------------------
  ! read u
  call h5open_f(error)
  error = -1
  do while (error==-1)
    call h5fopen_f (trim(trim(dirname)//trim(fname)//'_u.h5'), H5F_ACC_RDWR_F, file_id, error)
    call sleep(1)
    print*, error
  end do

  call h5dopen_f(file_id, 'var', dset_id, error)
  call h5dget_space_f(dset_id,dspace_id,error)
  call h5sget_simple_extent_dims_f(dspace_id, dims_u, maxdims, error)

  nx = int(dims_u(1))
  ny = int(dims_u(2))
  nz = int(dims_u(3))
  nz_tot = nz+1;


  if (.not. allocated(flow3D%u)) allocate(flow3D%u(nx+1, ny+1, nz))
  allocate(u_temp(nx,ny,nz))
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE,flow3D%u(1:nx,1:ny,1:nz) , dims_u, error)
  call h5dclose_f(dset_id, error)
  ! Close the file
  call h5fclose_f(file_id, error)
  ! Close FORTRAN interface.
  call h5close_f(error)
  ! read v
  call h5open_f(error)
  call h5fopen_f (trim(trim(dirname)//trim(fname)//'_v.h5'), H5F_ACC_RDWR_F, file_id, error)
  call h5dopen_f(file_id, 'var', dset_id, error)
  call h5dget_space_f(dset_id,dspace_id,error)
  call h5sget_simple_extent_dims_f(dspace_id, dims_u, maxdims, error)
  if (.not. allocated(flow3D%v)) allocate(flow3D%v(nx+1, ny+1, nz))
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, flow3D%v(1:nx,1:ny,1:nz), dims_u, error)
  call h5dclose_f(dset_id, error)
  ! Close the file
  call h5fclose_f(file_id, error)
  ! Close FORTRAN interface.
  call h5close_f(error)

  ! read w
  call h5open_f(error)
  call h5fopen_f (trim(trim(dirname)//trim(fname)//'_w.h5'), H5F_ACC_RDWR_F, file_id, error)
  call h5dopen_f(file_id, 'var', dset_id, error)
  call h5dget_space_f(dset_id,dspace_id,error)
  call h5sget_simple_extent_dims_f(dspace_id, dims_u, maxdims, error)
  if (.not. allocated(flow3D%w)) allocate(flow3D%w(nx+1, ny+1, nz))
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, flow3D%w(1:nx,1:ny,1:nz), dims_u, error)
  call h5dclose_f(dset_id, error)
  ! Close the file
  call h5fclose_f(file_id, error)
  ! Close FORTRAN interface.
  call h5close_f(error)


  ! read theta
  INQUIRE(FILE=trim(trim(dirname)//trim(fname)//'_theta.h5'), EXIST=file_exists)
  if (file_exists .eqv. .true.) then
  print*, 'read temperature profile from files'
  call h5open_f(error)
  call h5fopen_f (trim(trim(dirname)//trim(fname)//'_theta.h5'), H5F_ACC_RDWR_F, file_id, error)
  call h5dopen_f(file_id, 'var', dset_id, error)
  call h5dget_space_f(dset_id,dspace_id,error)
  call h5sget_simple_extent_dims_f(dspace_id, dims_u, maxdims, error)
  if (.not. allocated(flow3D%c)) allocate(flow3D%c(nx+1, ny+1, nz))
  call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, flow3D%c(1:nx,1:ny,1:nz), dims_u, error)
  call h5dclose_f(dset_id, error)
  ! Close the file
  call h5fclose_f(file_id, error)
  ! Close FORTRAN interface.
  call h5close_f(error)

  ! conversion from potential temp to absolute temperature 
  flow3D%c = flow3D%c * (1013.25/1000)**(0.286) * tinput%T_scale_K
  ! conversion from temperature to sound speed 
  flow3D%c = sqrt(1.4*286*flow3D%c)

  else

  if (.not. allocated(flow3D%c)) allocate(flow3D%c(nx+1, ny+1, nz))
  print*, 'set sound speed constant'
  flow3D%c = c0
  end if 
  ! print*, 'set sound speed constant'
  ! flow3D%c = c0
  ! define dimsension and source position
  !-------------------------------------------------------------------------
  ! dimensional size of the domain
  Lx = tinput%Lx*tinput%z_i
  Ly = tinput%Ly*tinput%z_i
  Lz = tinput%Lz*tinput%z_i

  !pos_x = 1.875*z_i; %for case A
  !pos_x = 0.83333*z_i; %for case B
  flow3D%Xs = tinput%posx*tinput%z_i
  flow3D%Ys = tinput%posy*tinput%z_i
  flow3D%Zs = tinput%posz*tinput%z_i/tinput%delta

  !ipos_x = floor(pos_x/L_x*nx) + 1;
  !ipos_y = floor(pos_y/L_y*ny) + 1;

  ! Copy boundary conditions
  !-------------------------------------------------------------------------
  flow3D%u(nx+1,:,:) = flow3D%u(1,:,:)
  flow3D%u(:,ny+1,:) = flow3D%u(:,1,:)
  flow3D%u(:,:,2:nz) = 0.5*(flow3D%u(:,:,2:nz) + flow3D%u(:,:,1:nz-1))
  flow3D%u(:,:,1) = 0.

  flow3D%v(nx+1,:,:) = flow3D%v(1,:,:);
  flow3D%v(:,ny+1,:) = flow3D%v(:,1,:);
  flow3D%v(:,:,2:nz) = 0.5*(flow3D%v(:,:,2:nz) + flow3D%v(:,:,1:nz-1))
  flow3D%v(:,:,1) = 0.

  flow3D%w(nx+1,:,:) = flow3D%w(1,:,:);
  flow3D%w(:,ny+1,:) = flow3D%w(:,1,:);

  ! apply ratio on mean flow
  flow3D%u= tinput%ratio*flow3D%u
  flow3D%v= tinput%ratio*flow3D%v
  flow3D%w= tinput%ratio*flow3D%w
  ! define grid mesh
  !-------------------------------------------------------------------------
  call linspace(0.,Lx,nx+1,flow3D%x)
  call linspace(0.,Ly,ny+1,flow3D%y)
  call linspace(0.,Lz,nz, flow3D%z)
  !z_w2 = z_w2(1:nz);
  !z_uv2 = z_w2+z_w2(2)/2 ; %% Location of u,v
  !z_uv3 = [0 z_uv2];
end subroutine read_h5flow3D
end module read_write
