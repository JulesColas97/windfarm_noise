#ifdef TURBINES
!*******************************************************************************
subroutine turbine_init
!*******************************************************************************
use param, only: z_i,dx,dy,dz,pi,coord,rp,path,L_z,nx,ny,nz,initu
use turbine_mod, only: turbine,filter_size,nloc,turb_count_buf,turb_max_buf
use turbine_mod, only: turb_data_buf,turb_nodes,turbine_in_proc,dyn_theta1
use turbine_mod, only: sumA,buffer_array
implicit none
integer, parameter :: lun=42
real :: b,c,d,e,f,g,h,t1,t2,t3,t4
real(rp) :: delta2,x_dist,y_dist
integer :: k,j,icp,jcp,kcp,imax,jmax,ktmax,k_start,k_end,ios 
logical :: exst
character(4) :: ipfi
82 format(i4.4)
182 format(4(e16.7))

! Define the turbine array
nullify(turbine)
allocate(turbine(nloc))
allocate(sumA(nloc))
allocate(buffer_array(nloc))

! Global z indices 
k_start= 1   +coord*(nz-1)
k_end  = nz-1+coord*(nz-1)

! Read the data from the input table
if (coord.eq.0) write(*,*) '   x_loc    y_loc     Diameter  Height    Thick     Vol       Theta     Ctprime'

! Read turbine locations
open(lun,file=trim(path)//'turbine_input.dat',status='unknown',form='formatted',action='read',position='rewind')
read(lun,*)

do k=1,nloc
  read(lun,*) j,b,c,d,e,f,g,h
  turbine(k)%loc(1)      = b
  turbine(k)%loc(2)      = c
  turbine(k)%dia         = d/z_i
  turbine(k)%loc(3)      = e/z_i
  turbine(k)%thk         = max ( f/z_i, dx*1.01_rp )
  turbine(k)%theta1      = g
  turbine(k)%Ct_prime    = h
  turbine(k)%thk_half    = 0.5_rp*turbine(k)%thk
  turbine(k)%turbine_vol = pi/4._rp * turbine(k)%dia**2 * turbine(k)%thk
  turbine(k)%vol_c       = dx*dy*dz / turbine(k)%turbine_vol

  if (initu == .true. .and. dyn_theta1 == .true.) then
   ios = 0
   write(ipfi,82) k
   inquire(file=trim(path)//'/turbine/turbine'//trim(ipfi)//'.dat',exist=exst)
   if(exst) then
    open(1000+k,file=trim(path)//'/turbine/turbine'//trim(ipfi)//'.dat',status='unknown',action='read',position='rewind')
     do while(ios.eq.0)
      read(1000+k,182,IOSTAT=ios) t1,t2,t3,t4
     end do
     turbine(k)%theta1= t3
     if (coord == 0) write(*,*) "Restarted with yaw:",turbine(k)%theta1
     close(1000+k)
    else
     if (coord == 0) write(*,*) "RESTART DATA FOR TURBINE NOT PRESENT - yaw taken from turbine_input.dat instead."
    endif
  endif

  ! Index of turbines center:
  icp = nint(turbine(k)%loc(1)/dx)
  jcp = nint(turbine(k)%loc(2)/dy)
  kcp = nint(turbine(k)%loc(3)/dz)

  ! Index of maximal distance to turbine center
  imax = nint(turbine(k)%dia/dx+2)
  jmax = nint(turbine(k)%dia/dy+2)
  ktmax= nint(turbine(k)%dia/dz+2)

  ! Save index for loop in turbine_nodes_func (see below)
  turbine(k)%min_i = icp-imax
  turbine(k)%max_i = icp+imax
  turbine(k)%min_j = jcp-jmax
  turbine(k)%max_j = jcp+jmax
  turbine(k)%min_k = max(kcp-ktmax,k_start)
  turbine(k)%max_k = min(kcp+ktmax,k_end)

  ! Check whether the turbine is in this processor
  if(kcp+ktmax.ge.k_start .and. kcp-ktmax.le.k_end) then
    turbine_in_proc = .true.
    turb_nodes = 4*(imax+1)*(jmax+1)*(nz-1)*turbine(k)%thk/turbine(k)%dia
    allocate(turbine(k)%nodes(3,turb_nodes))
    allocate(turbine(k)%ind(turb_nodes))
    turbine(k)%nodes = 0
    turbine(k)%ind = 0._rp
    turbine(k)%num_nodes = 0
  endif

  ! Check if turbines are in domain / not too close to border of domain:
  if( icp <= imax .or. icp >= nx - imax ) stop 'Turbine out of x domain.'
  if( jcp <= jmax .or. jcp >= ny - jmax ) stop 'Turbine out of y domain. (If this is wanted, comment out this line and comment in line ***)'
  if( turbine(k)%loc(3) + 0.5_rp*turbine(k)%dia >  0.5_rp * L_z ) stop 'Turbine too high.'


if (coord.eq.0) write(*,"(8F10.4)") turbine(k)%loc(1), turbine(k)%loc(2),    &
                     turbine(k)%dia,  turbine(k)%loc(3), turbine(k)%thk,       &
                     turbine(k)%vol_c,turbine(k)%theta1, turbine(k)%Ct_prime
enddo

close(lun)
! End read the data from the input table 

! Check if spacing of turbines is at least 2.5 times turbine diameter:
do k=1,nloc
do j=1,nloc
  if (k .ne. j) then
    x_dist = abs(turbine(k)%loc(1) - turbine(j)%loc(1))
    y_dist = abs(turbine(k)%loc(2) - turbine(j)%loc(2))
    if(sqrt(x_dist**2 + y_dist**2) < 2.5_rp * turbine(k)%dia) then
      print*, 'Turbine spacing too small.'
      stop
    endif
  endif
enddo
enddo

! Check to see whether the grid spacing is smaller than the turbine radius
do k = 1,nloc
  if (max(dy,dz) .ge. 0.5_rp*turbine(k)%dia) then
    write(*,'(5F12.4)') real(k),dx,dy,dz,0.5_rp*turbine(k)%dia
    print*, 'Grid way too coarse to accurate resolve turbine'
    stop
  endif
enddo

! Calculate the radial component of indicator function using FFT
! Ref. Shapiro et al. Wind Energy 1-7 (2019)
delta2 = filter_size**2 * (dx**2 + dy**2 + dz**2) 
do k=1,nloc 
  call turbine(k)%turb_ind_func%init(delta2,turbine(k)%dia,turbine(k)%Ct_prime)
enddo 

! Finds the grid cells that are located in the turbine area
call turbine_nodes_func

! Initialize buffer feature for turbine data
turb_count_buf = 0
allocate(turb_data_buf(4,nloc,turb_max_buf)); turb_data_buf = 0.0_rp

end subroutine turbine_init


!*******************************************************************************
subroutine turbine_nodes_func
!*******************************************************************************
use param, only: dx,dy,dz,ny,nz,pi,localComm,mpi_rp,ierr,rp,coord,R1_ind,L_y
use grid_defs, only: gridx,gridy,z_tot
use mpi, only: mpi_sum
use turbine_mod, only: turbine,filter_cutoff,nloc
use turbine_mod, only: turbine_in_proc,turb_nodes,sumA,buffer_array
use turbine_indicator, only: sqrt6overdelta
implicit none
real(rp) :: rx,ry,rx2,ry2,r2,r_norm,r_norm2,r_disk,filt,dxyz
integer :: count_i,i,j,k,s,j2

sumA(:) = 0.0_rp
dxyz = dx*dy*dz

if(turbine_in_proc) then
  do s=1,nloc

    ! determine unit normal vector for each turbine
    turbine(s)%nhat(1) = -cosd(turbine(s)%theta1)
    turbine(s)%nhat(2) = -sind(turbine(s)%theta1)

    turbine(s)%num_nodes = 0
    ! index count - used for writing to array "nodes"
    count_i = 0

    ! vector from center point to this node is (rx,ry,rz) with length r
    do j=turbine(s)%min_j,turbine(s)%max_j
      j2 = mod(j+ny-1,ny)+1
      ry = gridy(j2) - turbine(s)%loc(2)
      !*** In case turbines touch the boundaries, uncomment the line below: 
      !if (j2 .ne. j) ry = min(abs(ry-L_y),abs(ry+L_y))
      ry2 = ry*ry

      do i=turbine(s)%min_i,turbine(s)%max_i
        rx = gridx(i) - turbine(s)%loc(1)
        rx2= rx*rx

        ! length projected onto unit normal for this turbine
        r_norm = abs(rx*turbine(s)%nhat(1) + ry*turbine(s)%nhat(2))

        ! Factor in normal plane. Calculate this in advanc before doing vertical loop.
        ! Here assume R1_ind is negligible small when r_norm>2*thk.
        if(r_norm .le. 2._rp*turbine(s)%thk) then
          r_norm2 = r_norm*r_norm
          ! normalized normal component of Gaussian-filtered indicator function
          ! see ref. Shapiro et al. (2019) eqn(10)
          R1_ind = 0.5_rp*(erf(sqrt6overdelta*(turbine(s)%thk_half+r_norm)) +  &
                           erf(sqrt6overdelta*(turbine(s)%thk_half-r_norm)))

          ! Vertical loop
          do k=turbine(s)%min_k,turbine(s)%max_k
            r2 = rx2 + ry2 + (z_tot(k) - turbine(s)%loc(3))**2
            ! (remaining) length projected onto turbine disk
            r_disk = sqrt(r2 - r_norm2)
            ! evaluate the Gaussian-filtered indicator function
            filt = turbine(s)%turb_ind_func%val(r_disk)

            ! Store this point as turbine point when indicator value is high enough
            if(filt > filter_cutoff) then
              count_i = count_i + 1
              sumA(s) = sumA(s) + filt*dxyz
              turbine(s)%ind(count_i) = filt
              turbine(s)%nodes(1,count_i) = i
              turbine(s)%nodes(2,count_i) = j2
              turbine(s)%nodes(3,count_i) = k-coord*(nz-1) ! local
            endif
          enddo ! end vertical loop
        endif ! end rnorm if statement
      enddo
    enddo

    turbine(s)%num_nodes = count_i
    if ((turbine(s)%num_nodes) .gt. turb_nodes) then
      write(*,*) 'Memory problem in turbine routine'
      stop
    endif

  enddo
endif

! Sum the indicator function across all processors and incorporate normalization in sumA
buffer_array = sumA
call mpi_Allreduce(buffer_array, sumA, nloc, mpi_rp, mpi_sum, localComm, ierr)

if(turbine_in_proc) then
  ! Normalize the indicator function
  do s=1,nloc
  do i=1,turbine(s)%num_nodes
    turbine(s)%ind(i) = turbine(s)%ind(i) * turbine(s)%turbine_vol / sumA(s)
  enddo
  enddo 
endif

! Print out some information on the turbine; used to check the code
!write(*,*) coord,turbine(1)%num_nodes,sum(turbine(1)%ind(1:turbine(1)%num_nodes)), sumA(1)

end subroutine turbine_nodes_func


!*******************************************************************************
subroutine turbine_forcing
!*******************************************************************************
!
! This subroutine applies the drag-disk forcing
!
use sim_param, only: u,v,RHSx,RHSy,fx,fy
use mpi, only: mpi_sum
use param, only: dt,z_i,pi,dx,dy,localComm,mpi_rp,ierr,jt_total,rp,total_time
use turbine_mod, only: tbase,nloc,dyn_theta1,turbine_nodes_allocated
use turbine_mod, only: T_avg,turbine,turb_count_buf,turb_data_buf
implicit none
real(rp), pointer :: p_Ct_prime => null()
real(rp) :: u_d,f_n,ind2,turbine_avg_vels,turbine_avg_vels2,phi,phi_old 
real(rp) :: turbine_u,turbine_v,turbine_u2,turbine_v2,dist_1D,eps
integer :: x_int,y_int,s,i2,j2,k2,l

! add one to the write buffer counter
if(modulo(jt_total,tbase)==0) turb_count_buf = turb_count_buf+1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! If dyn_theta = .true. turbine heads are rotated 
if(dyn_theta1 .and. turbine_nodes_allocated) then 

  do s=1,nloc

    turbine_u = 0._rp
    turbine_v = 0._rp
    turbine_u2 = 0._rp
    turbine_v2 = 0._rp

    phi_old = turbine(s)%theta1
    dist_1D = turbine(s)%dia/((dx + dy)*0.5_rp)

    ! Counts the number of indices upstream of the turbine
    x_int = -nint(dist_1D * cosd(phi_old))
    y_int = -nint(dist_1D * sind(phi_old))

    do l=1,turbine(s)%num_nodes
      i2 = turbine(s)%nodes(1,l)
      j2 = turbine(s)%nodes(2,l)
      k2 = turbine(s)%nodes(3,l)

      turbine_u = turbine_u + turbine(s)%ind(l)*u(i2+x_int,j2+y_int,k2)
      turbine_v = turbine_v + turbine(s)%ind(l)*v(i2+x_int,j2+y_int,k2)
    end do

    call mpi_allreduce(turbine_u,turbine_u2,1,mpi_rp,mpi_sum,localComm,ierr)
    call mpi_allreduce(turbine_v,turbine_v2,1,mpi_rp,mpi_sum,localComm,ierr)

    turbine_u2 = turbine_u2 * turbine(s)%vol_c
    turbine_v2 = turbine_v2 * turbine(s)%vol_c

    ! Angle is fixed between 0 and 360
    ! phi = modulo(atan2d(turbine_v2,turbine_u2),360.0) 

    phi = atan2d(turbine_v2, turbine_u2)

    ! Time filter for turbine yaw 
    eps = ((dt*z_i)/T_avg)/(1._rp + (dt*z_i)/T_avg)

    phi = (1._rp-eps)*phi_old + eps*phi 

    ! Phi has to be in degrees 
    turbine(s)%theta1 = phi
  end do

  ! Recompute the turbine position if theta1 or theta2 can change
  call turbine_nodes_func

endif 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do s=1,nloc
  ! Calculate the disk-averaged velocity using the updated indicator function
  turbine_avg_vels = 0._rp 

  do l=1,turbine(s)%num_nodes   
    i2 = turbine(s)%nodes(1,l)
    j2 = turbine(s)%nodes(2,l)
    k2 = turbine(s)%nodes(3,l)
    turbine_avg_vels = turbine_avg_vels + turbine(s)%ind(l) *                  &
        (turbine(s)%nhat(1)*u(i2,j2,k2) + turbine(s)%nhat(2)*v(i2,j2,k2))
  end do

  call mpi_allreduce(turbine_avg_vels,turbine_avg_vels2,1,mpi_rp,mpi_sum,localComm,ierr)

  ! set pointers
  p_Ct_prime => turbine(s)%Ct_prime
   
  ! volume correction:
  ! since sum of ind is turbine volume/(dx*dy*dz) (not exactly 1.)
  u_d = turbine_avg_vels2 * turbine(s)%vol_c 

  ! Multiplying correction factor (see Shapiro, Gayme and Meneveau (2019)
  ! "Filtered actuator disks: ...", formula 26)
  u_d = u_d * turbine(s)%turb_ind_func%M

  ! calculate total thrust force for each turbine (per unit mass)
  ! force is normal to the surface (calc from u_d_T, normal to surface)
  ! write force to array that will be transferred via MPI
  f_n = 0.5_rp*p_Ct_prime*abs(u_d)*u_d/turbine(s)%thk

  ! Buffer the write statements to the turbine files
  if (modulo (jt_total,tbase)==0) then
    turb_data_buf(1,s,turb_count_buf) = total_time
    turb_data_buf(2,s,turb_count_buf) = -u_d
    turb_data_buf(3,s,turb_count_buf) = turbine(s)%theta1
    turb_data_buf(4,s,turb_count_buf) = f_n
  endif

  do l=1,turbine(s)%num_nodes
    i2 = turbine(s)%nodes(1,l)
    j2 = turbine(s)%nodes(2,l)
    k2 = turbine(s)%nodes(3,l)
    ind2 = -f_n*turbine(s)%ind(l)
    ! Applied forcing directly to right-hand side 
    RHSx(i2,j2,k2) = RHSx(i2,j2,k2) + ind2*turbine(s)%nhat(1)
    RHSy(i2,j2,k2) = RHSy(i2,j2,k2) + ind2*turbine(s)%nhat(2)

#ifdef TURBINES
    fx(i2,j2,k2) = ind2 * turbine(s)%nhat(1)
    fy(i2,j2,k2) = ind2 * turbine(s)%nhat(2)
#endif TURBINES 
  end do

end do 

! Prevent rotating the turbine head in the first step already initialized by the
! turbine_nodes_func call in turbines_init
turbine_nodes_allocated = .true.

call turbine_write_buffer(0)

end subroutine turbine_forcing


!*******************************************************************************
subroutine turbine_write_buffer(endflag)
!*******************************************************************************
use param, only: coord,path
use turbine_mod, only: turb_count_buf,turb_max_buf,turb_data_buf,nloc
implicit none
integer,intent(in) :: endflag ! 1 forces the data to be written
integer :: k,s
character(4) :: ipfi

82 format(i4.4)
182 format(4(e16.7))

! Write the buffered turbine data to file
if ((turb_count_buf.eq.turb_max_buf) .or. (endflag.eq.1)) then
 if(coord.eq.0) then ! Only do the writing at master
  do s=1,nloc
   write(ipfi,82) s
   open(1000+s,file=trim(path)//'/turbine/turbine'//trim(ipfi)//'.dat',status='unknown',form='formatted',position='append')
   do k=1,turb_count_buf
    write(1000+s,182) turb_data_buf(1:4,s,k)
   enddo
   close(1000+s)
  enddo
 endif
  turb_count_buf = 0 ! Reset the counter
endif

end subroutine turbine_write_buffer
#endif TURBINES
