#ifdef WINDBREAKS
!*******************************************************************************
subroutine windbreak_init
!*******************************************************************************
use param, only: z_i,dx,dy,dz,coord,rp,path,nx,ny,nz
use windbreak_mod, only: windbreak,filter_size_wb,nloc_wb,sqrt6overdelta_wb
use windbreak_mod, only: windbreak_nodes,windbreak_in_proc
use windbreak_mod, only: sumB,buffer_array_wb
implicit none
integer, parameter :: lun=42
real :: a,b,c,d,e,f,g,h
real(rp) :: delta2
integer :: k,j,icp,jcp,kcp,imax,jmax,ktmax,k_start,k_end 

! Define the windbreak array
nullify(windbreak)
allocate(windbreak(nloc_wb))
allocate(sumB(nloc_wb))
allocate(buffer_array_wb(nloc_wb))

! Global z indices 
k_start= 1   +coord*(nz-1)
k_end  = nz-1+coord*(nz-1)

! Read the data from the input table
if(coord.eq.0) write(*,*) '   x_loc    y_loc     z_loc     Width  ' //         &
                          '   Height   Thick     Theta     Cd'

! Read windbreak locations
open(lun, file=trim(path)//'windbreak_input.dat', status='unknown',            &
          form='formatted', action='read', position='rewind')
read(lun,*)

do k=1,nloc_wb
  read(lun,*) j,a,b,c,d,e,f,g,h
  windbreak(k)%loc(1)  = a
  windbreak(k)%loc(2)  = b
  windbreak(k)%loc(3)  = c
  windbreak(k)%width   = max(d/z_i, dy*1.01_rp)
  windbreak(k)%height  = e/z_i
  windbreak(k)%thk     = max(f/z_i, dx*1.01_rp)
  windbreak(k)%vol     = windbreak(k)%width*windbreak(k)%height*windbreak(k)%thk
  windbreak(k)%vol_c   = dx*dy*dz/windbreak(k)%vol
  windbreak(k)%theta1  = g
  windbreak(k)%Cd      = h

  ! Index of windbreak center:
  icp = nint(windbreak(k)%loc(1)/dx)
  jcp = nint(windbreak(k)%loc(2)/dy)
  kcp = nint(windbreak(k)%loc(3)/dz)

  ! Index of maximal distance to windbreak center
  imax = nint(windbreak(k)%thk   /dx+2)
  jmax = nint(windbreak(k)%width /dy+2)
  ktmax= nint(windbreak(k)%height/dz+2)

  ! Save index for loop in windbreak_nodes_func (see below)
  windbreak(k)%min_i = max(icp-imax,1)
  windbreak(k)%max_i = min(icp+imax,nx)
  windbreak(k)%min_j = max(jcp-jmax,1)
  windbreak(k)%max_j = min(jcp+jmax,ny)
  windbreak(k)%min_k = max(kcp-ktmax,k_start)
  windbreak(k)%max_k = min(kcp+ktmax,k_end)

  ! Check whether the windbreak is in this processor
  if(kcp+ktmax.ge.k_start .and. kcp-ktmax.le.k_end) then
    windbreak_in_proc = .true.
    windbreak_nodes = 4*(imax+1)*(jmax+1)*(nz-1)
    allocate(windbreak(k)%nodes(3,windbreak_nodes))
    allocate(windbreak(k)%ind(windbreak_nodes))
    windbreak(k)%nodes = 0._rp
    windbreak(k)%ind = 0._rp
    windbreak(k)%num_nodes = 0
  endif

  if (coord.eq.0) then
    write(*,"(8F10.4)") windbreak(k)%loc(1), windbreak(k)%loc(2),              &
                        windbreak(k)%loc(3), windbreak(k)%width,               &
                        windbreak(k)%height, windbreak(k)%thk,                 &
                        windbreak(k)%theta1, windbreak(k)%Cd
  endif
enddo

close(lun)
! End read the data from the input table 

! Check to see whether the grid spacing is smaller than the windbreak height
do k = 1,nloc_wb
  if(dz .ge. 0.5_rp*windbreak(k)%height) then
    print*, 'k=', k, 'dz=', dz, 'height_half=', 0.5_rp*windbreak(k)%height
    print*, 'Grid way too coarse to accurate resolve windbreak'
    stop
  endif
enddo

! Ref. Shapiro et al. Wind Energy 1-7 (2019)
delta2 = filter_size_wb**2 * (dx**2 + dy**2 + dz**2) 
sqrt6overdelta_wb = sqrt(6._rp) / sqrt(delta2)

! Finds the grid cells that are located in the windbreak area
call windbreak_nodes_func

end subroutine windbreak_init


!*******************************************************************************
subroutine windbreak_nodes_func
!*******************************************************************************
use param, only: dx, dy, dz, ny, nz, localComm, mpi_rp, ierr, rp, coord
use grid_defs, only: gridx, gridy,z_tot
use mpi, only: mpi_sum
use windbreak_mod, only: windbreak, filter_cutoff_wb, nloc_wb
use windbreak_mod, only: windbreak_in_proc, windbreak_nodes
use windbreak_mod, only: sqrt6overdelta_wb, sumB, buffer_array_wb
implicit none
real(rp) :: rx,ry,rx2,ry2,r_norm,r_norm2,dxyz
real(rp) :: filt, R1_ind, R2_ind, R3_ind
integer :: count_i,i,j,k,s

sumB(:) = 0.0_rp
dxyz = dx*dy*dz

if(windbreak_in_proc) then
  do s=1,nloc_wb

    ! determine unit normal vector for each windbreak
    windbreak(s)%nhat(1) = -cosd(windbreak(s)%theta1)
    windbreak(s)%nhat(2) = -sind(windbreak(s)%theta1)

    windbreak(s)%num_nodes = 0
    ! index count - used for writing to array "nodes"
    count_i = 0

    ! vector from center point to this node is (rx,ry,rz) with length r
    do j=windbreak(s)%min_j,windbreak(s)%max_j
      ry  = gridy(j) - windbreak(s)%loc(2)
      ry2 = ry*ry

      do i=windbreak(s)%min_i,windbreak(s)%max_i
        rx = gridx(i) - windbreak(s)%loc(1)
        rx2= rx*rx

        ! length projected onto unit normal for this windbreak
        r_norm = abs(rx*windbreak(s)%nhat(1) + ry*windbreak(s)%nhat(2))

        ! Factor in normal plane. Calculate this in advanc before doing vertical
        ! loop. Here assume R1_ind is negligible small when r_norm>2*thk.
        if(r_norm .le. 2._rp*windbreak(s)%thk) then
          r_norm2 = r_norm*r_norm
          ! normalized normal component of Gaussian-filtered indicator function
          ! see ref. Shapiro et al. (2019) eqn(10)
          R1_ind = 0.5_rp*(erf(sqrt6overdelta_wb*(windbreak(s)%thk/2._rp+r_norm))+&
                           erf(sqrt6overdelta_wb*(windbreak(s)%thk/2._rp-r_norm)))

          ! square of length projected onto unit tangential for this windbreak
          r_norm2 = rx2 + ry2 - r_norm2
          r_norm = sqrt(r_norm2)
          R2_ind = 0.5_rp*(erf(sqrt6overdelta_wb*(windbreak(s)%width/2._rp+r_norm))+&
                           erf(sqrt6overdelta_wb*(windbreak(s)%width/2._rp-r_norm)))

          ! Two-dimensional windbreak: no smooth in this direction
          if(windbreak(s)%min_j.eq.1 .or. windbreak(s)%max_j.eq.ny) R2_ind=1._rp

          ! Vertical loop
          do k=windbreak(s)%min_k,windbreak(s)%max_k
            r_norm2 = (z_tot(k) - windbreak(s)%loc(3))**2
            r_norm = sqrt(r_norm2)
            R3_ind = 0.5_rp*(erf(sqrt6overdelta_wb*(windbreak(s)%height/2._rp+r_norm))+&
                             erf(sqrt6overdelta_wb*(windbreak(s)%height/2._rp-r_norm)))

            filt = R1_ind*R2_ind*R3_ind
            ! Store this point as windbreak point when indicator value is high &
            ! enough
            if(filt > filter_cutoff_wb) then
              count_i = count_i + 1
              sumB(s) = sumB(s) + filt*dxyz
              windbreak(s)%ind(count_i) = filt
              windbreak(s)%nodes(1,count_i) = i
              windbreak(s)%nodes(2,count_i) = j
              windbreak(s)%nodes(3,count_i) = k-coord*(nz-1) ! local
            endif
          enddo ! end vertical loop
        endif ! end rnorm if statement
      enddo
    enddo

    windbreak(s)%num_nodes = count_i
    if ((windbreak(s)%num_nodes) .gt. windbreak_nodes) then
      write(*,*) 'Memory problem in windbreak routine'
      stop
    endif

  enddo
endif

! Sum the indicator function across all processors and incorporate normalization
! in sumB
buffer_array_wb = sumB
call mpi_Allreduce(buffer_array_wb, sumB, nloc_wb, mpi_rp, mpi_sum, localComm, ierr)

if(windbreak_in_proc) then
  ! Normalize the indicator function
  do s=1,nloc_wb
  do i=1,windbreak(s)%num_nodes
    windbreak(s)%ind(i) = windbreak(s)%ind(i) * windbreak(s)%vol / sumB(s)
  enddo
  enddo 
endif

end subroutine windbreak_nodes_func


!*******************************************************************************
subroutine windbreak_forcing
!*******************************************************************************
!
! This subroutine applies the drag-disk forcing
!
use sim_param, only: u,v,RHSx,RHSy
use mpi, only: mpi_sum
use param, only: localComm,mpi_rp,ierr,rp,jt_total,dt,z_i,coord,u_d_filter
use windbreak_mod, only: nloc_wb, windbreak
implicit none
real(rp) :: u_d,f_n,ind2,windbreak_avg_vels,windbreak_avg_vels2,eps_avg
integer :: s,i2,j2,k2,l

do s=1,nloc_wb
  ! Calculate the disk-averaged velocity using the updated indicator function
  windbreak_avg_vels = 0._rp 

  do l=1,windbreak(s)%num_nodes   
    i2 = windbreak(s)%nodes(1,l)
    j2 = windbreak(s)%nodes(2,l)
    k2 = windbreak(s)%nodes(3,l)
    windbreak_avg_vels = windbreak_avg_vels + windbreak(s)%ind(l) *            &
          (windbreak(s)%nhat(1)*u(i2,j2,k2) + windbreak(s)%nhat(2)*v(i2,j2,k2))
  end do

  call mpi_allreduce(windbreak_avg_vels,windbreak_avg_vels2,1,mpi_rp,mpi_sum,  &
                     localComm,ierr)
   
  ! volume correction:
  ! since sum of ind is windbreak volume/(dx*dy*dz) (not exactly 1.)
  u_d = windbreak_avg_vels2 * windbreak(s)%vol_c 

  ! Time average (to avoid instabilities)
  if(jt_total < 10000) u_d_filter = u_d
  
  ! In general, eps_avg = (DT / T_avg) / (1 + DT / T_avg). T_avg is chosen as 1 s.
  eps_avg = (dt * z_i / 1._rp) / (1._rp + dt * z_i / 1._rp)
  u_d_filter = (1._rp - eps_avg) * u_d_filter + eps_avg * u_d
  u_d = u_d_filter
  
  ! calculate total thrust force for each windbreak (per unit mass)
  ! force is normal to the surface
  ! write force to array that will be transferred via MPI
  f_n = 0.5_rp*windbreak(s)%Cd*abs(u_d)*u_d/windbreak(s)%thk
  
  ! If using coarse grid or the code diverges, the next line could be needed. It slowly adds the force until jt = 10000.
  ! if(jt_total.lt.10000) f_n = real(jt_total)/real(10000)*f_n

  do l=1,windbreak(s)%num_nodes
    i2 = windbreak(s)%nodes(1,l)
    j2 = windbreak(s)%nodes(2,l)
    k2 = windbreak(s)%nodes(3,l)
    ind2 = -f_n*windbreak(s)%ind(l)
    ! Applied forcing directly to right-hand side 
    RHSx(i2,j2,k2) = RHSx(i2,j2,k2) + ind2*windbreak(s)%nhat(1)
    RHSy(i2,j2,k2) = RHSy(i2,j2,k2) + ind2*windbreak(s)%nhat(2)
  end do

end do 

end subroutine windbreak_forcing
#endif WINDBREAKS
