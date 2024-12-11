module interpo

implicit none

public :: read_binflow
public :: write_binflow
public :: interpolate2
public :: create_2D_interpoland
public :: interpolate_2D
public :: interpo_vector
public :: interp1
public :: interp2
public :: interpo_1D
public :: interplog
public :: interpMean

contains

subroutine read_binflow(fname,u,v,w,x,y,z)
    implicit none

    ! input variables
    !-------------------------------------------------------------------------------

    character(len=100),   intent(in)                     :: fname     !écoulement d'entré
    real,allocatable, dimension(:,:,:), intent(inout)    :: u,v,w       ! check for temp
    real,allocatable, dimension(:), intent(inout)        :: x,y,z      ! check for temp

    ! local vraiables
    !-------------------------------------------------------------------------------
    integer                                 :: nx, ny, nz
    real                                    :: dx, dy, dz
    !                           Lecture ecoulement
    !-------------------------------------------------------------------------------
   print*, fname
    open(91,file=fname,STATUS = "OLD", ACCESS = "STREAM", FORM = "UNFORMATTED")
    read(91)
    read(91) nx,ny, nz
    read(91) dx,dy, dz

    allocate(x(nx))
    allocate(y(ny))
    allocate(z(nz))

    allocate(u(nx,ny,nz))
    allocate(v(nx,ny,nz))
    allocate(w(nx,ny,nz))

    read(91) x, y, z, u, v, w
    close(91)
    print*, " "

    print'(A65)', '***************      MEAN FLOW PARAMETERS        ***************'
    print'(A65)', '*****************************************************************'
    print'(A13,F9.5,F9.5,F9.5,A1)', '(dx,dy,dz) =(',dx, dy, dz, ')'
    print'(A13,I9,I9,I9,A1)', '(nx,ny,nz) =(',nx, ny, nz, ')'
end subroutine read_binflow


subroutine write_binflow(fname, x, y, u, v)

    implicit none
    ! input variables
    !-------------------------------------------------------------------------------
    character(len=*),      intent(in)      :: fname
    real,dimension(:,:),     intent(in)      :: x, y, u, v

    ! local variables
    !-------------------------------------------------------------------------------
    integer                                 :: nx, ny

    ! write binary file
    !-------------------------------------------------------------------------------
    nx = size(x, 1)
    ny = size(x,2)

    open(90,file=trim(adjustl(fname))//'_flow.bin',   form='unformatted')
    write(90) nx, ny
    write(90) x, y
    write(90) u
    write(90) v
    close(90)
end subroutine write_binflow


subroutine interpo_receiver(x,f,nx,x0,f0)
  use bspline_sub_module
  implicit none

  ! inout variables
  real, dimension(:), intent(inout)                         :: x
  complex, dimension(:), intent(inout)                      :: f
  real, intent(in)                                          :: x0
  integer, intent(in)                                       :: nx

  ! output
  complex, intent(inout)                                    :: f0

  ! local variables
  integer                                                   :: kx = 3     !! order in x
  integer                                                   :: iknot = 0
  real                                                      :: tx(nx+4)
  real,dimension(3*4)                                       :: w1_1d
  integer                                                   :: iflag
  integer                                                   :: inbvx, idx
  real                                                      :: u(nx)
  real                                                      :: fi, fr
  idx = 0
  inbvx = 1


  ! real part interpolation
  call db1ink(x,nx,real(f),kx,iknot,tx,u,iflag)
  call db1val(x0,idx,tx,nx,kx,u,fr,iflag,inbvx,w1_1d)

  idx = 0
  inbvx = 1
  ! imaginary part
  call db1ink(x,nx,aimag(f),kx,iknot,tx,u,iflag)
  call db1val(x0,idx,tx,nx,kx,u,fi,iflag,inbvx,w1_1d)

  f0 = cmplx(fr, fi)
end subroutine interpo_receiver


subroutine interpo_1D(x,y,xnew,ynew)
  use bspline_sub_module
  implicit none

  ! inout variables
  real, dimension(:), intent(in)                            :: x
  real, dimension(:), intent(in)                            :: y
  real, dimension(:), intent(in)                            :: xnew

  ! output
  real,dimension(:), intent(inout)                          :: ynew

  ! local variables
  integer                                                   :: nx
  integer                                                   :: kx = 3     !! order in x
  integer                                                   :: iknot = 0
  real, dimension(:),allocatable                            :: tx
  real,dimension(3*3)                                       :: w1_1d
  integer                                                   :: iflag
  integer                                                   :: inbvx, idx
  real, dimension(:), allocatable                           :: knot
  integer                                                   :: ii

  idx = 0
  inbvx = 1

  nx = size(x)
  allocate(knot(nx))
  allocate(tx(nx+3))
  ! real part interpolation
  call db1ink(x,nx,y,kx,iknot,tx,knot,iflag)
  do ii=1, size(xnew)
    call db1val(xnew(ii),idx,tx,nx,kx,knot,ynew(ii),iflag,inbvx,w1_1d)
  enddo
end subroutine interpo_1D


subroutine interpo_vector(z1,f1,nz1,z2,f2,nz2,kz)
  use bspline_sub_module
  implicit none

  ! inout variables
  real, dimension(:), intent(in)                            :: z1
  real, dimension(:), intent(in)                            :: f1
  real, dimension(:), intent(inout)                         :: z2
  real, dimension(:), intent(inout)                         :: f2
  integer, intent(in)                                       :: nz1, nz2
  integer, intent(in)                                       :: kz


  ! local variables
  integer                                                   :: iknot = 0
  real                                                      :: tx(nz1+kz)
  real,dimension(3*kz)                                       :: w1_1d
  integer                                                   :: iflag
  integer                                                   :: inbvx, idx
  real                                                      :: u(nz1)
  integer                                                   :: kk

  idx = 0
  inbvx = 1
  ! real part interpolation
  call db1ink(z1,nz1,f1,kz,iknot,tx,u,iflag)
  do kk =1, nz2
    call db1val(z2(kk),idx,tx,nz1,kz,u,f2(kk),iflag,inbvx,w1_1d)
  enddo

end subroutine interpo_vector


subroutine interp1( x, y, xNew, yNew )
  ! Inputs: x = a vector of the x-values of the data to be interpolated
  !         y = a vector of the y-values of the data to be interpolated
  !         xVal  = a vector of the x-values where interpolation should be performed
  ! Output: yNew  = a vector of the resulting interpolated values

  implicit none

  real, intent(in) :: x(:), y(:), xNew(:)
  real, intent(out) :: yNew(:)
  integer :: newIndex, dataIndex, dataIndexOld
  real :: weight
  !   size(x) == size(y)
  !   size(xNew) == size(yNew)
  ! loop over the new points on which to interpolate
  dataIndexOld = 1
  do newIndex = 1, size(xNew)
      dataIndex = dataIndexOld 
      ! find the two data points around the new x value
      do while ((x(dataIndex+1)<xNew(newIndex)) .and. ((dataIndex+2) < size(x)))
        ! exit the loop if end of the xdata is reached
        dataIndex = dataIndex + 1
      enddo
      dataIndexOld = dataIndex
      ! compute weghts for linear interpolation
      if ((dataIndex+1) <= size(x)) then
        weight = (xNew(newIndex) - x(dataIndex))/(x(dataIndex+1)-x(dataIndex))
        ! print*, y(dataIndex), y(dataIndex+1), weight
        yNew(newIndex) = (1.0-weight)*y(dataIndex) + weight*y(dataIndex+1)
      endif
  end do
end subroutine interp1


subroutine interpMean( xData, yData, xVal, yVal )
  ! Inputs: xData = a vector of the x-values of the data to be interpolated
  !         yData = a vector of the y-values of the data to be interpolated
  !         xVal  = a vector of the x-values where interpolation should be performed
  ! Output: yVal  = a vector of the resulting interpolated values

  implicit none

  real, intent(in) :: xData(:), yData(:), xVal(:)
  real, intent(out) :: yVal(:)
  integer :: inputIndex, dataIndex
  real :: weight
  !   size(xData) == size(yData)
  !   size(xVal) == size(yVal)
  do inputIndex = 1, size(xVal)-1
      dataIndex = 1
      do while ((xData(dataIndex+1)<xVal(inputIndex)))
        dataIndex = dataIndex + 1
      enddo
      ! weight = (xVal(inputIndex) - xData(dataIndex))/(xData(dataIndex+1)-xData(dataIndex))
      yVal(inputIndex) = 0.5*(yData(dataIndex) + yData(dataIndex+1))
  end do
end subroutine


subroutine interplog( x, y, xNew, yNew )
  ! Inputs: x = a vector of the x-values of the data to be interpolated
  !         y = a vector of the y-values of the data to be interpolated
  !         xVal  = a vector of the x-values where interpolation should be performed
  ! Output: yNew  = a vector of the resulting interpolated values
  implicit none
  real, intent(in) :: x(:), y(:), xNew(:)
  real, intent(out) :: yNew(:)
  integer :: newIndex, dataIndex, dataIndexOld
  real :: weight
  if ((size(x) .ne. size(y)) .or. (size(xNew) .ne. size(yNew))) then 
    print*, 'error_interpo : x and y must have same size '
    yNew = 0.
    return
  endif 
  dataIndexOld = 1
  ! loop over the new points on which to interpolate
  do newIndex = 1, size(xNew)
      if ((xNew(newIndex)>x(size(x))) .or. (xNew(newIndex)<x(1))) then 
        print*, 'error_interpo : xNew outside of xData'
        print*, xNew(newIndex), x(size(x)), xNew(newIndex), x(1)
        yNew(newIndex) = 0
      else 
        dataIndex = dataIndexOld
        ! find the two data points around the new x value
        do while ((x(dataIndex+1)<xNew(newIndex)) .and. ((dataIndex+2) <= size(x)))
          ! exit the loop if end of the xdata is reached
          dataIndex = dataIndex + 1
        enddo
        dataIndexOld = dataIndex
        ! compute weghts for linear interpolation
        if ((dataIndex+1) <= size(x)) then
          weight = (xNew(newIndex) - x(dataIndex))/(x(dataIndex+1)-x(dataIndex))
          ! print*, dataIndex,size(x), newIndex, size(xNew)
          ! print*, y(dataIndex)
          ! print*,  y(dataIndex+1)
          ! print*, weight


          yNew(newIndex) = 10*log10((1.0-weight)*10**(y(dataIndex)/10) + weight*10**(y(dataIndex+1)/10))
        endif
      endif 
  end do
end subroutine


subroutine interp2( xData, yData, xVal, yVal )
  ! Inputs: xData = a vector of the x-values of the data to be interpolated
  !         yData = a vector of the y-values of the data to be interpolated
  !         xVal  = a vector of the x-values where interpolation should be performed
  ! Output: yVal  = a vector of the resulting interpolated values

  implicit none

  real, intent(in) :: xData(:), yData(:), xVal(:)
  real, intent(out) :: yVal(:)
  integer :: inputIndex, dataIndex
  real :: weight

  ! Possible checks on inputs could go here
  ! Things you may want to check:
  !   monotonically increasing xData
  !   size(xData) == size(yData)
  !   size(xVal) == size(yVal)

  do inputIndex = 1, size(xVal)
    ! print*, inputIndex, 'over',  size(xVal)
      dataIndex = 1
      do while ((xData(dataIndex+1)<xVal(inputIndex)) .or. (dataIndex+1<size(xData)))
        dataIndex = dataIndex + 1
      enddo
      weight = (xVal(inputIndex) - xData(dataIndex))/(xData(dataIndex+1)-xData(dataIndex))
      yVal(inputIndex) = (1.0-weight)*yData(dataIndex) + weight*yData(dataIndex+1)
  end do
end subroutine


subroutine interpolate2(u1, v1, w1, x1, y1, z1, u2, v2, Xs, Ys, Zs, x2, y2, theta)
    use bspline_sub_module

    implicit none
    ! input variables
    !-------------------------------------------------------------------------------
    real, dimension (:,:,:),intent(in)   :: u1
    real, dimension (:,:,:),intent(in)   :: v1
    real, dimension (:,:,:),intent(in)   :: w1
    real, dimension (:), intent(in)      :: x1
    real, dimension (:), intent(in)      :: y1
    real, dimension (:), intent(in)      :: z1

    real,allocatable,  dimension (:,:),intent(inout)     :: u2
    real,allocatable,  dimension (:,:),intent(inout)     :: v2

    real, intent(in)                        :: Xs
    real, intent(in)                        :: Ys
    real, intent(in)                        :: Zs
    real, intent(in)                        :: theta

    real, dimension (:,:),intent(in)        :: x2
    real, dimension (:,:),intent(in)        :: y2

    real,allocatable, dimension (:,:,:)                 :: uKnot
    real,allocatable, dimension (:,:,:)                 :: vKnot
    real,allocatable, dimension (:,:,:)                 :: wKnot
    real,allocatable, dimension (:,:,:)                 :: u1_new
    real,allocatable, dimension (:,:,:)                 :: v1_new
    real,allocatable, dimension (:,:,:)                 :: w1_new
    real,allocatable,dimension(:)                       :: x1_new
    real,allocatable,dimension(:)                       :: y1_new
    real,allocatable,dimension(:)                       :: z1_new


    ! local variables
    integer                                 :: nx1,ny1, nz1, nx2, ny2
    integer                                 :: nx1_new,ny1_new, nz1_new
    real                                    :: dx1, dy1, dx2
    real                                    :: xmax, xmin, xmaxOverx,xmaxOvery,xminOverx,xminOvery
    integer                                 :: nxAddPos,nxAddNeg,nyAddPos,nyAddNeg
    real                                    :: V0ux,V0uy     !! velocity component over x and y
    integer                                 :: ixi, ieta,ii

    !         Declaration of locale variabales used for the interpolation
    !-------------------------------------------------------------------------------
    integer,parameter                       :: kx    = 3     !! order in x
    integer,parameter                       :: ky    = 3     !! order in y
    integer,parameter                       :: kz    = 3     !! order in z
    integer,parameter                       :: iknot = 0
    real, allocatable, dimension(:,:)       :: tx,ty,tz
    real,dimension(ky,kz,3)                 :: w1_3d
    real,dimension(kz,3)                    :: w2_3d
    real,dimension(3*max(kx,ky,kz),3)       :: w3_3d
    integer, dimension(3)                   :: inbvx,inbvy,inbvz,iloy,iloz, idx,idy,idz

    integer,dimension(6,3)                  :: iflag
    real                                    :: pi = 3.141592653589793

    ! initialisation
    !-------------------------------------------------------------------------------
    nx1 = size(x1)
    ny1 = size(y1)
    nz1 = size(z1)
    dx1 = x1(2) - x1(1)
    dy1 = y1(2) - y1(1)

    nx2 = size(x2, 1)
    ny2 = size(x2, 2)
    dx2 = x2(2,1) - x2(1,1)

    allocate(u2(nx2,ny2))
    allocate(v2(nx2,ny2))

    !have to set these before the first evaluate call:
    inbvx = 1
    inbvy = 1
    inbvz = 1
    iloy  = 1
    iloz  = 1
    idx = 0
    idy = 0
    idz = 0

    !-------------------------------------------------------------------------------------------
    !       calculation of  nxi max if the 2D mesh is larger than the 3D grid in x direction
    !-------------------------------------------------------------------------------------------

    xmax = x2(nx2,1)
    xmin = x2(1,1)

    ! projection of 2D plane dimension to the 3D axis
    xmaxOverx = xmax*cos(theta)
    xmaxOvery = xmax*sin(theta)
    xminOverx = xmin*cos(theta + pi)
    xminOvery = xmin*sin(theta + pi)

    ! different cases according to the size of the projection
    ! add points  at the beginning or ath end of the array
    if (xmaxOverx >= 0) then
      nxAddPos = floor(max (xmaxOverx + Xs - x1(nx1) , 0.)/dx1)
      nxAddNeg =  floor(max(xminOverx - Xs, 0.)/dx1)
    else
      nxAddNeg =  floor(max(xmaxOverx - Xs, 0.)/dx1)
      nxAddPos = floor(max (xminOverx - x1(nx1) + Xs, 0.)/dx1)
    endif

    if (xmaxOvery >= 0) then
      nyAddPos = floor(max (xmaxOvery - y1(ny1) + Ys, 0.)/dy1)
      nyAddNeg =  floor(max(xminOvery - Ys, 0.)/dy1)
    else
      nyAddNeg =  floor(max(xmaxOvery - Ys, 0.)/dy1)
      nyAddPos = floor(max (xminOvery - y1(ny1) + Ys, 0.)/dy1)
    endif
    
    print*, 'added ', nxAddNeg, ' in -x direction'
    print*, 'added ', nxAddPos, ' in +x direction'
    print*, 'added ', nyAddNeg, ' in -y direction'
    print*, 'added ', nyAddPos, ' in +y direction'

    ! set new dimensions
    nx1_new = nx1+nxAddNeg+nxAddPos
    ny1_new = ny1+nyAddNeg+nyAddPos
    nz1_new = nz1


    allocate(u1_new(nx1_new,ny1_new,nz1_new))
    allocate(v1_new(nx1_new,ny1_new,nz1_new))
    allocate(w1_new(nx1_new,ny1_new,nz1_new))
    allocate(uKnot(nx1_new,ny1_new,nz1_new))
    allocate(vKnot(nx1_new,ny1_new,nz1_new))
    allocate(wKnot(nx1_new,ny1_new,nz1_new))
    allocate(tx(nx1_new+kx,3))
    allocate(ty(ny1_new+ky,3))
    allocate(tz(nz1_new+kz,3))
    allocate(x1_new(nx1_new))
    allocate(y1_new(ny1_new))
    allocate(z1_new(nz1_new))


    ! create the x,y,z vector
    do ii=1,nx1_new
      x1_new(ii) = dx1 * (ii-1 - nxAddNeg)
    end do

    do ii=1,ny1_new
      y1_new(ii) = dy1 * (ii-1 - nyAddNeg)
    end do

    z1_new = z1

    ! copy of the orginal values in the new array
    u1_new(nxAddneg+1:nx1_new-nxAddPos,nyAddNeg+1:ny1_new-nyAddPos,:) = u1
    v1_new(nxAddneg+1:nx1_new-nxAddPos,nyAddNeg+1:ny1_new-nyAddPos,:) = v1
    w1_new(nxAddneg+1:nx1_new-nxAddPos,nyAddNeg+1:ny1_new-nyAddPos,:) = w1

    ! recopy of the first or last value if the domain extend further
    do ii = 1,nxAddNeg
      u1_new(ii,:,:) = u1_new(nxAddNeg,:,:)
      v1_new(ii,:,:) = v1_new(nxAddNeg,:,:)
      w1_new(ii,:,:) = w1_new(nxAddNeg,:,:)
    end do

    do ii = nx1_new-nxAddPos+1,nx1_new
      u1_new(ii,:,:) = u1_new(nx1_new-nxAddPos,:,:)
      v1_new(ii,:,:) = v1_new(nx1_new-nxAddPos,:,:)
      w1_new(ii,:,:) = w1_new(nx1_new-nxAddPos,:,:)
    end do

    do ii =1,NyAddNeg
      u1_new(:,ii,:) = u1_new(:,nyAddNeg,:)
      v1_new(:,ii,:) = v1_new(:,nyAddNeg,:)
      w1_new(:,ii,:) = w1_new(:,nyAddNeg,:)
    end do

    do ii = ny1_new-nyAddPos+1,ny1_new
      u1_new(:,ii,:) = u1_new(:,ny1_new-nyAddPos,:)
      v1_new(:,ii,:) = v1_new(:,ny1_new-nyAddPos,:)
      w1_new(:,ii,:) = w1_new(:,ny1_new-nyAddPos,:)
    end do

    ! create knots for interpolation
    call db3ink(x1_new-Xs,nx1_new,y1_new-Ys,ny1_new,z1_new,nz1_new,u1_new, &
            & kx,ky,kz,iknot,tx(:,1),ty(:,1),tz(:,1),uKnot,iflag(3,1))
    call db3ink(x1_new-Xs,nx1_new,y1_new-Ys,ny1_new,z1_new,nz1_new,v1_new,&
            & kx,ky,kz,iknot,tx(:,2),ty(:,2),tz(:,2),vKnot,iflag(3,2))
    call db3ink(x1_new-Xs,nx1_new,y1_new-Ys,ny1_new,z1_new,nz1_new,w1_new,&
            & kx,ky,kz,iknot,tx(:,3),ty(:,3),tz(:,3),wKnot,iflag(3,3))

    ! loop over the 2D points of the domain
    do ieta = 1,ny2
        do ixi = 1,nx2
            ! interpolation of v
           call db3val(cos(theta)*x2(ixi,ieta),sin(theta)*x2(ixi,ieta),y2(ixi,ieta),idx(3),idy(3),idz(3),&
                                               tx(:,3),ty(:,3),tz(:,3),nx1_new,ny1_new,nz1_new,kx,ky,kz,&
                                               wKnot,v2(ixi,ieta) ,iflag(3,3),&
                                               inbvx(3),inbvy(3),inbvz(3),iloy(3),iloz(3),&
                                               w1_3d(:,:,3),w2_3d(:,3),w3_3d(:,3),.true.)
            ! interpolation of u
            call db3val(cos(theta)*x2(ixi,ieta),sin(theta)*x2(ixi,ieta),y2(ixi,ieta),idx(1),idy(1),idz(1),&
                                                tx(:,1),ty(:,1),tz(:,1),nx1_new,ny1_new,nz1_new,kx,ky,kz,&
                                                uKnot,V0ux ,iflag(3,1),&
                                                inbvx(1),inbvy(1),inbvz(1),iloy(1),iloz(1),&
                                                w1_3d(:,:,1),w2_3d(:,1),w3_3d(:,1),.true.)

            call db3val(cos(theta)*x2(ixi,ieta),sin(theta)*x2(ixi,ieta),y2(ixi,ieta),idx(2),idy(2),idz(2),&
                                                tx(:,2),ty(:,2),tz(:,2),nx1_new,ny1_new,nz1_new,kx,ky,kz,&
                                                vKnot,V0uy ,iflag(3,2),&
                                                inbvx(2),inbvy(2),inbvz(2),iloy(2),iloz(2),&
                                                w1_3d(:,:,2),w2_3d(:,2),w3_3d(:,2),.true.)
            u2(ixi,ieta) = cos(theta)*V0ux + sin(theta)*V0uy
        enddo
    enddo
end subroutine interpolate2


subroutine create_2D_interpoland(u1, x1, y1, uKnot,tx,ty)
    use bspline_sub_module

    implicit none
    ! input variables
    !-------------------------------------------------------------------------------
    real, dimension (:,:),intent(in)                      :: u1
    real, dimension (:), intent(in)                       :: x1
    real, dimension (:), intent(in)                       :: y1

    real,allocatable,  dimension (:,:),intent(inout)      :: uKnot
    real,allocatable,  dimension (:),intent(inout)        :: tx,ty

    ! local variables
    integer                                 :: nx1,ny1

    !         Declaration of locale variabales used for the interpolation
    !-------------------------------------------------------------------------------
    integer,parameter                       :: kx    = 3     !! order in x
    integer,parameter                       :: ky    = 3     !! order in y
    integer,parameter                       :: iknot = 0
    INTEGER                                 :: iflag

    ! initialisation
    !-------------------------------------------------------------------------------
    nx1 = size(x1)
    ny1 = size(y1)

    allocate(uKnot(nx1,ny1))
    allocate(tx(nx1+kx))
    allocate(ty(ny1+ky))

    ! create knots for interpolation
    call db2ink(x1,nx1,y1,ny1,u1,kx,ky,iknot,tx,ty,uKnot,iflag)
end subroutine create_2D_interpoland


subroutine interpolate_2D(uKnot,tx,ty, x2, y2, uval)
    use bspline_sub_module

    implicit none
    ! input variables
    !-------------------------------------------------------------------------------

    real,allocatable,  dimension (:,:),intent(in)       :: uKnot
    real,allocatable,  dimension (:),intent(in)         :: tx,ty

    real, dimension (:,:),intent(in)                    :: x2
    real, dimension (:,:),intent(in)                    :: y2
    real,allocatable,  dimension (:,:),intent(inout)    :: uval

    ! local variables
    integer                                             :: nx1,ny1
    integer                                             :: nx2,ny2
    integer                                             :: ii,jj

    !         Declaration of locale variabales used for the interpolation
    !-------------------------------------------------------------------------------
    integer,parameter                       :: kx    = 3     !! order in x
    integer,parameter                       :: ky    = 3     !! order in y
    integer,parameter                       :: iknot = 0
    real,dimension(ky)                      :: w1
    real,dimension(3*max(kx,ky))            :: w0
    integer                                 :: inbvx,inbvy,inbvz,iloy, idx,idy,iflag

    ! initialisation
    !-------------------------------------------------------------------------------
    nx1 = size(uKnot,1)
    ny1 = size(uKnot,2)

    nx2 = size(x2,1)
    ny2 = size(x2,2)

    allocate(uval(nx2,ny2))

    !have to set these before the first evaluate call:
    inbvx = 1
    inbvy = 1
    iloy  = 1
    idx = 0
    idy = 0

    ! evalutae function 
    do ii = 1,nx2
        do jj = 1,ny2
          call db2val(x2(ii,jj),y2(ii,jj),idx,idy,tx,ty,nx1,ny1,kx,ky,uKnot,uval(ii,jj),iflag,inbvx,inbvy,iloy,w1,w0,.true.)
        enddo 
    enddo 
  end subroutine interpolate_2D



end module interpo
