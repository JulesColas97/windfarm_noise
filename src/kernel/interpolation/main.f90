program main
use interpo 


implicit none

character(len=100)                              :: fname
real,allocatable, dimension(:,:,:)              :: u1,v1,w1
real,allocatable,dimension(:)                   :: x1,y1,z1

real,allocatable, dimension(:,:)                :: u2, v2, x2, y2
real,allocatable, dimension(:)                  :: xi, eta
real                                            :: Xs, Ys, Zs, theta 

integer                                         :: ii, jj, kk
integer                                         :: nx, ny
real                                            :: dx

! USER input 

nx = 2000
ny = 400
dx = 0.5
Xs = 1171.875
Ys = 312.5
Zs = 100
theta = 0

! initilisation fine mesh 
allocate(x2(nx, ny))
allocate(y2(nx, ny))
allocate(xi(nx))
allocate(eta(ny))
do ii = 1, nx
    xi(ii) = dx*(ii - 1)
enddo

do ii = 1, ny
    eta(ii) = dx*(ii - 1)
    x2(:,ii) = xi(:)
enddo

do ii = 1, nx
    y2(ii,:) = eta(:)
enddo
fname = "flow.bin"
print*, fname
call  read_binflow(fname, u1,v1,w1,x1,y1,z1)
call    interpolate(u1, v1, w1, x1, y1, z1, u2, v2, Xs, Ys, Zs, x2,y2, theta)

fname = "interpolated_flow"
print*, fname
call  write_binflow(fname, x2, y2, u2, v2)

end program main 