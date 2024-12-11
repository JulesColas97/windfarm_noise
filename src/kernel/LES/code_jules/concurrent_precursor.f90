#ifdef CPS
!*******************************************************************************
subroutine create_mpi_comms_cps
!*******************************************************************************
!
! This subroutine does two things. It first splits the mpi_comm_world
! communicator into two communicators (localComm). The two new communicators are 
! then bridged to create an intercommunicator (interComm).
!
use param, only: ierr,localComm
use cps_mod, only: color,RED,BLUE,interComm
use mpi, only: mpi_comm_world
implicit none
integer :: world_np,world_rank,remoteLeader,memberKey

! Get number of processors in world comm
call mpi_comm_size (mpi_comm_world, world_np, ierr)
call mpi_comm_rank (mpi_comm_world, world_rank, ierr)

! Set color and remote leader for intercommunicator interComm
if(world_rank < world_np/2) then
  color = RED
  remoteLeader = world_np/2
else
  color = BLUE
  remoteLeader = 0
endif

! Generate member key
memberKey = modulo(world_rank, world_np/2)

! Split the world communicator into intracommunicators localComm
call mpi_Comm_split(mpi_comm_world, color, memberKey, localComm, ierr)

! Create intercommunicator interComm
call mpi_intercomm_create(localComm, 0, mpi_comm_world, remoteLeader, 1,interComm, ierr)

end subroutine create_mpi_comms_cps


!*******************************************************************************
subroutine initialize_cps
!*******************************************************************************
use param, only: nx,ny,nz,rp,sgs_model,symmetric_fringe
use param, only: coord,rank_of_coord,status,ierr
use cps_mod, only: color,RED,BLUE,interComm
use cps_mod, only: alpha,beta,vel_sample_t 
use mpi, only: mpi_INTEGER
use cps_mod, only: nx_conp
#ifdef CPS_Y
use cps_mod, only: ny_conp
use param, only: fringe_region_end_y
use cps_mod, only: alpha_y,beta_y
#endif CPS_Y
implicit none
integer :: i,index
integer :: istart_p, iplateau_p, iend_p
#ifdef CPS_Y
integer :: istart_yp, iplateau_yp, iend_yp 
#endif CPS_Y

if(color == BLUE) then
  call fringe_init(istart_p, iplateau_p, iend_p)

  ! Sample size same as buffer region (omitting istart from the block
  ! since velocity is already set there)
  nx_conp = iend_p - istart_p
   
  ! Send size of the sample block to upstream domain (RED)
  call mpi_send(nx_conp, 1, mpi_INTEGER, rank_of_coord(coord), 1, interComm, ierr)

  ! Now setup fringe weights 
  allocate(alpha(nx_conp), beta(nx_conp))
  index = 0
  do i = istart_p+1, iend_p
    index = index+1
    if(symmetric_fringe == .true.) then
     call fringe_weighting_symmetric(i, istart_p, iend_p, beta(index))
    else
     call fringe_weighting(i, istart_p, iplateau_p, beta(index))
    endif
  enddo
  alpha = 1.0_rp - beta

#ifdef CPS_Y
  call fringe_init_y(istart_yp, iplateau_yp, iend_yp)

  ! Sample size same as buffer region (omitting istart from the block since 
  ! velocity is already set there)
  ny_conp = iend_yp - istart_yp
   
  ! Send size of the sample block to upstream domain (RED)
  call mpi_send(ny_conp, 1, mpi_INTEGER, rank_of_coord(coord), 1, interComm, ierr)
 
  ! Now setup fringe weights 
  allocate(alpha_y(ny_conp), beta_y(ny_conp))
  index = 0
  do i = istart_yp+1, iend_yp
    index = index+1
    if(symmetric_fringe == .true.) then
     call fringe_weighting_symmetric(i, istart_yp, iend_yp, beta_y(index))
    else
     call fringe_weighting(i, istart_yp, iplateau_yp, beta_y(index))
    endif
  enddo
  alpha_y = 1.0_rp - beta_y
#endif CPS_Y

elseif(color == RED) then
  ! Receive from downstream domain (BLUE) the length of the sample block
  call mpi_recv(nx_conp, 1, mpi_INTEGER, rank_of_coord(coord), 1, interComm,status, ierr)

  ! Should end up as nx + 1 (this eventually gets wrapped)
  iend_p = nx
  ! Plateau location not used since no fringe treatment on the RED domain, but 
  ! setting so it is at least initialized.
  iplateau_p = iend_p
  ! Set istart based on the size of the sample block
  istart_p = iend_p - nx_conp

#ifdef CPS_Y
  ! Receive from downstream domain (BLUE) the length of the sample block
  call mpi_recv(ny_conp, 1, mpi_INTEGER, rank_of_coord(coord), 1, interComm,status, ierr)
  
  if(fringe_region_end_y >= 1._rp) then  
    ! Should end up as nx + 1 (this eventually gets wrapped) 
    iend_yp = ny
    ! Plateau location not used since no fringe treatment on the RED domain,
    ! but setting so it is at least initialized.
    iplateau_yp = iend_yp
    ! Set istart based on the size of the sample block
    istart_yp = iend_yp - ny_conp
  else 
    istart_yp = 1 
    iplateau_yp = iend_yp 
    iend_yp = istart_yp + ny_conp
  endif 
#endif CPS_Y
endif

! Allocate and assign wrapped index and fringe weights
allocate(vel_sample_t % iwrap(nx_conp))
index = 0
do i = istart_p+1, iend_p
  index = index+1
  vel_sample_t % iwrap(index) = modulo(i-1, nx) + 1
enddo

! Allocate the sample block
allocate( vel_sample_t % u( nx_conp, ny, nz ) )
allocate( vel_sample_t % v( nx_conp, ny, nz ) )
allocate( vel_sample_t % w( nx_conp, ny, nz ) )

#ifdef SCALAR  
allocate( vel_sample_t % theta(nx_conp, ny, nz) ) 
#endif SCALAR

#ifdef CPS_Y
! Allocate and assign wrapped index and fringe weights
allocate(vel_sample_t % iwrap_y(ny_conp))
index = 0
do i = istart_yp+1, iend_yp
  index = index+1
  vel_sample_t % iwrap_y(index) = modulo(i-1, ny) + 1
enddo
  
! Allocate the sample block
allocate( vel_sample_t % u_y( nx, ny_conp, nz ) )
allocate( vel_sample_t % v_y( nx, ny_conp, nz ) )
allocate( vel_sample_t % w_y( nx, ny_conp, nz ) )
  
#ifdef SCALAR  
allocate( vel_sample_t % theta_y(nx, ny_conp, nz) ) 
#endif SCALAR
#endif CPS_Y

end subroutine initialize_cps


!*******************************************************************************
subroutine synchronize_cps
!*******************************************************************************
use param, only: ny,nz,rp
use param, only: coord, rank_of_coord, status, ierr, mpi_rp
use sim_param, only: u,v,w
use cps_mod, only: color,RED,BLUE,interComm
use cps_mod, only: vel_sample_t
use cps_mod, only: nx_conp
#ifdef CPS_Y
use cps_mod, only: ny_conp
use param, only: nx
#endif CPS_Y
#ifdef SCALAR 
use scalars_param, only: theta
#endif SCALAR
implicit none

integer :: sendsize, recvsize
integer, pointer, dimension(:) :: iwrap_p
real(rp), pointer, dimension(:,:,:) :: u_p,v_p,w_p
#ifdef SCALAR 
real(rp), pointer, dimension(:,:,:) :: theta_p
#endif SCALAR

#ifdef CPS_Y
integer :: sendsize_y, recvsize_y
integer, pointer, dimension(:) :: iwrap_yp
real(rp), pointer, dimension(:,:,:) :: u_yp,v_yp,w_yp
#ifdef SCALAR 
real(rp), pointer, dimension(:,:,:) :: theta_yp
#endif SCALAR
#endif CPS_Y
 
nullify( u_p, v_p, w_p ,iwrap_p)

iwrap_p  => vel_sample_t % iwrap
u_p      => vel_sample_t % u
v_p      => vel_sample_t % v
w_p      => vel_sample_t % w

sendsize = nx_conp * ny * nz
recvsize = sendsize

#ifdef SCALAR
nullify(theta_p)
theta_p   => vel_sample_t % theta 
#endif SCALAR

#ifdef CPS_Y
nullify(u_yp, v_yp, w_yp,iwrap_yp)

iwrap_yp   => vel_sample_t % iwrap_y 
u_yp       => vel_sample_t % u_y
v_yp       => vel_sample_t % v_y
w_yp       => vel_sample_t % w_y

sendsize_y = nx * ny_conp * nz
recvsize_y = sendsize_y

#ifdef SCALAR
nullify(theta_yp)
theta_yp   => vel_sample_t % theta_y 
#endif SCALAR
#endif CPS_Y

if(color == BLUE) then

! Receive sampled velocities from upstream (RED)
call mpi_recv(    u_p(1,1,1), recvsize, mpi_rp, rank_of_coord(coord), 1,interComm, status, ierr)
call mpi_recv(    v_p(1,1,1), recvsize, mpi_rp, rank_of_coord(coord), 2,interComm, status, ierr)
call mpi_recv(    w_p(1,1,1), recvsize, mpi_rp, rank_of_coord(coord), 3,interComm, status, ierr)
#ifdef SCALAR 
call mpi_recv(theta_p(1,1,1), recvsize, mpi_rp, rank_of_coord(coord), 4,interComm, status, ierr)
#endif SCALAR
  
#ifdef CPS_Y
! Recieve sampled velocities from upstream (RED)
call mpi_recv(    u_yp(1,1,1), recvsize_y, mpi_rp, rank_of_coord(coord), 5,interComm, status, ierr)
call mpi_recv(    v_yp(1,1,1), recvsize_y, mpi_rp, rank_of_coord(coord), 6,interComm, status, ierr)
call mpi_recv(    w_yp(1,1,1), recvsize_y, mpi_rp, rank_of_coord(coord), 7,interComm, status, ierr)
#ifdef SCALAR 
call mpi_recv(theta_yp(1,1,1), recvsize_y, mpi_rp, rank_of_coord(coord), 8,interComm, status, ierr)
#endif SCALAR
#endif CPS_Y

elseif(color == RED) then

! Sample velocity and copy to buffers
u_p(:,:,:) = u(iwrap_p(:),1:ny,1:nz)
v_p(:,:,:) = v(iwrap_p(:),1:ny,1:nz)
w_p(:,:,:) = w(iwrap_p(:),1:ny,1:nz)

! Send sampled velocities to downstream domain (BLUE)
call mpi_send(u_p    (1,1,1), sendsize, mpi_rp, rank_of_coord(coord), 1,interComm, ierr)
call mpi_send(v_p    (1,1,1), sendsize, mpi_rp, rank_of_coord(coord), 2,interComm, ierr)
call mpi_send(w_p    (1,1,1), sendsize, mpi_rp, rank_of_coord(coord), 3,interComm, ierr)
#ifdef SCALAR 
theta_p(:,:,:) = theta(iwrap_p(:),1:ny,1:nz)
call mpi_send(theta_p(1,1,1), sendsize, mpi_rp, rank_of_coord(coord), 4,interComm, ierr)
#endif SCALAR 

#ifdef CPS_Y
  ! Sample velocity and copy to buffers
  u_yp(:,:,:) = u(1:nx,iwrap_yp(:),1:nz)
  v_yp(:,:,:) = v(1:nx,iwrap_yp(:),1:nz)
  w_yp(:,:,:) = w(1:nx,iwrap_yp(:),1:nz)

  ! Send sampled velocities to downstream domain (BLUE)
  call mpi_send(    u_yp(1,1,1), sendsize_y, mpi_rp, rank_of_coord(coord), 5,interComm, ierr)
  call mpi_send(    v_yp(1,1,1), sendsize_y, mpi_rp, rank_of_coord(coord), 6,interComm, ierr)
  call mpi_send(    w_yp(1,1,1), sendsize_y, mpi_rp, rank_of_coord(coord), 7,interComm, ierr)
#ifdef SCALAR 
  theta_yp(:,:,:) = theta(1:nx,iwrap_yp(:),1:nz)
  call mpi_send(theta_yp(1,1,1), sendsize_y, mpi_rp, rank_of_coord(coord), 8,interComm, ierr)
#endif SCALAR
#endif CPS_Y

endif

nullify(u_p, v_p, w_p, iwrap_p)
#ifdef SCALAR
nullify(theta_p)
#endif SCALAR
 
#ifdef CPS_Y
nullify(u_yp, v_yp, w_yp, iwrap_yp)
#ifdef SCALAR
nullify(theta_yp)
#endif
#endif CPS_Y
 
end subroutine synchronize_cps

!*******************************************************************************
subroutine inflow_cond_cps
!*******************************************************************************
!
! Enforces prescribed inflow condition from an inlet velocity field generated 
! from a precursor simulation. The inflow condition is enforced by direct
! modulation on the velocity in the fringe region.
!
use param, only: rp,ny,nz,rp
use sim_param, only: u, v, w, RHSx, RHSy, RHSz
use cps_mod, only: alpha,beta,vel_sample_t
#ifdef CPS_Y
use param, only: nx
use cps_mod, only: alpha_y,beta_y
#endif CPS_Y
implicit none
integer :: j,k 

integer, pointer, dimension(:) :: iwrap_p
real(rp), pointer, dimension(:,:,:) :: u_p, v_p, w_p

#ifdef CPS_Y
integer, pointer, dimension(:) :: iwrap_yp
real(rp), pointer, dimension(:,:,:) :: u_yp, v_yp, w_yp
#endif CPS_Y

nullify( u_p, v_p, w_p )
nullify( iwrap_p )

u_p     => vel_sample_t % u
v_p     => vel_sample_t % v
w_p     => vel_sample_t % w
iwrap_p => vel_sample_t % iwrap

! Srinidhi : Feb 20, 2022 
! The concurrent precursor is now divergence free. However, it now requires a relaxation factor called 'fringe_factor'. 
! The 'fringe_factor' needs to be carefully selected for each case. 
! It is mostly trial and error. Some value which has worked for Anja & Jens has been set for now.
! Run some coarse grid case on the final simulation domain and make sure the outlet velocity in the BLUE domain is the same as in
! the RED domain. (Always run coarse grid simulations before running the final one). 
do k=1,nz-1
do j=1,ny
  RHSx(iwrap_p(:),j,k) = RHSx(iwrap_p(:),j,k) - beta(:) * (u(iwrap_p(:),j,k) - u_p(:,j,k))
  RHSy(iwrap_p(:),j,k) = RHSy(iwrap_p(:),j,k) - beta(:) * (v(iwrap_p(:),j,k) - v_p(:,j,k))
  RHSz(iwrap_p(:),j,k) = RHSz(iwrap_p(:),j,k) - beta(:) * (w(iwrap_p(:),j,k) - w_p(:,j,k))
enddo
enddo

nullify( u_p, v_p, w_p )
nullify( iwrap_p )

#ifdef CPS_Y
nullify( u_yp, v_yp, w_yp )
nullify( iwrap_yp )

u_yp     => vel_sample_t % u_y
v_yp     => vel_sample_t % v_y
w_yp     => vel_sample_t % w_y
iwrap_yp => vel_sample_t % iwrap_y

do k=1,nz-1
do j=1,nx
  RHSx(j,iwrap_yp(:),k) = RHSx(j,iwrap_yp(:),k) - beta_y(:)*(u(j,iwrap_yp(:),k)-u_yp(j,:,k))
  RHSy(j,iwrap_yp(:),k) = RHSy(j,iwrap_yp(:),k) - beta_y(:)*(v(j,iwrap_yp(:),k)-v_yp(j,:,k))
  RHSz(j,iwrap_yp(:),k) = RHSz(j,iwrap_yp(:),k) - beta_y(:)*(w(j,iwrap_yp(:),k)-w_yp(j,:,k))
enddo
enddo

nullify( u_yp, v_yp, w_yp )
nullify( iwrap_yp )
#endif CPS_Y

end subroutine inflow_cond_cps

#ifdef SCALAR 
!*******************************************************************************
subroutine inflow_cond_cps_scalar 
!*******************************************************************************
!
! Enforces prescribed inflow condition from an inlet velocity field generated 
! from a precursor simulation. The inflow condition is enforced by direct
! modulation on the velocity in the fringe region.
!
use param, only: rp,ny,nz,rp
use cps_mod, only: alpha,beta,vel_sample_t
#ifdef CPS_Y
use param, only: nx
use cps_mod, only: alpha_y,beta_y
#endif CPS_Y
use scalars_param, only: theta, RHS_t 
implicit none
integer :: j,k 

integer, pointer, dimension(:) :: iwrap_p
real(rp), pointer, dimension(:,:,:) :: theta_p

#ifdef CPS_Y
integer, pointer, dimension(:) :: iwrap_yp
real(rp), pointer, dimension(:,:,:) :: theta_yp
#endif CPS_Y

nullify( iwrap_p )

iwrap_p => vel_sample_t % iwrap
theta_p => vel_sample_t % theta

do k=1,nz-1
do j=1,ny
  RHS_t(iwrap_p(:),j,k) = RHS_t(iwrap_p(:),j,k) - beta(:) * (theta(iwrap_p(:),j,k) - theta_p(:,j,k))
enddo
enddo

nullify( iwrap_p )
nullify( theta_p )

#ifdef CPS_Y
nullify( iwrap_yp )
iwrap_yp => vel_sample_t % iwrap_y
theta_yp => vel_sample_t % theta_y

do k=1,nz-1
do j=1,nx
  RHS_t(j,iwrap_yp(:),k) = RHS_t(j,iwrap_yp(:),k) - beta_y(:)*(theta(j,iwrap_yp(:),k)-theta_yp(j,:,k))
enddo
enddo

nullify( iwrap_yp )
nullify( theta_yp )
#endif CPS_Y

end subroutine inflow_cond_cps_scalar 
#endif SCALAR 

#endif CPS
