#ifdef ADMR
!*******************************************************************************
subroutine admr_forcing
!*******************************************************************************
use param, only : rp,nproc,nx,ny,nz,dx,dy,dz,idx,idy,idz,dt,nz_tot,pi,total_time,jt_total
use param, only : mpi_rp, rank_of_coord, coord, status,ierr,localComm
use mpi, only: mpi_sum
use sim_param, only: u,v,w,RHSx,RHSy,RHSz
use admr_mod
implicit none

integer ::  i, j, k,l         !< loop indices
integer ::  iialpha, iir     !<
integer ::  i_rseg  !<
integer ::  i_ring   !< 
integer ::  ii, jj, kk       !<
integer ::  k_start,k_end,kkp    

real(rp) :: sin_rot, cos_rot   !<
real(rp) :: sin_yaw, cos_yaw,phi_yaw_old   !<
       
real(rp) :: aa, bb, cc, dd  !< interpolation distances
real(rp) :: xd,yd,zd        !< interpolation volume var  
real(rp) :: alpha_attack,phi_rel,turb_cl,turb_cd
real(rp) :: vtheta, vrel
!
! for yaw-control
real(rp) :: admr_u,admr_v,admr_u2,admr_v2,dist_1D,eps_admr
integer :: x_int,y_int
! add one to the write buffer counter
if (modulo (jt_total,admr_tbase)==0) admr_count_buf=admr_count_buf+1

!Global z indices 
k_start= 1   +coord*(nz-1)
k_end  = nz-1+coord*(nz-1)

!
!--       Set forces to zero for each new time step:
thrust(:,:,:)         = 0.0_rp
torque_y(:,:,:)       = 0.0_rp
torque_z(:,:,:)       = 0.0_rp
rot_tend_x(:,:,:)     = 0.0_rp
rot_tend_y(:,:,:)     = 0.0_rp
rot_tend_z(:,:,:)     = 0.0_rp
admr(:)%torque_total  = 0.0_rp
admr(:)%thrust_rotor  = 0.0_rp

!--             Loop over number of turbines:
do inot = 1,nturb
! if (admr_in_proc) then
  cos_yaw = COS(admr(inot)%phi_yaw)
  sin_yaw = SIN(admr(inot)%phi_yaw)

!--             Loop over rings of each turbine:
  do i_ring = 1,admr(inot)%nrings
   thrust_seg(:)   = 0.0_rp
   torque_seg_y(:) = 0.0_rp
   torque_seg_z(:) = 0.0_rp

!--             Determine distance between each ring (center) and the hub:

!--             Loop over segments of each ring of each turbine:
   do i_rseg = 1, admr(inot)%nsegs(i_ring)

!--                !-----------------------------------------------------------!
!--                !-- Determine coordinates of the ring segments            --!
!--                !-----------------------------------------------------------!
!
!--                Determine angle of ring segment towards zero degree angle of
!--                rotor system (at zero degree rotor direction vectors aligned
!--                with y-axis):
    phi_rotor = float(i_rseg) * 2.0_rp * pi / float(admr(inot)%nsegs(i_ring))
    cos_rot   = COS( phi_rotor )
    sin_rot   = SIN( phi_rotor )

!
!--                Coordinates of the single segments (center points):
    rbx(i_ring,i_rseg) = admr(inot)%loc(1) * z_i + admr(inot)%cur_r(i_ring) * cos_rot * sin_yaw
    rby(i_ring,i_rseg) = admr(inot)%loc(2) * z_i + admr(inot)%cur_r(i_ring) * cos_rot * cos_yaw
    rbz(i_ring,i_rseg) = admr(inot)%loc(3) * z_i + admr(inot)%cur_r(i_ring) * sin_rot

!--                !-----------------------------------------------------------!
!--                !-- Interpolation of the velocity components from the     --!
!--                !-- cartesian grid point to the coordinates of each ring  --!
!--                !-- segment (trilinear interpolation) --!
!--                !-----------------------------------------------------------!

    u_int(inot,i_ring,i_rseg)     = 0.0_rp
    u_int_1_l(inot,i_ring,i_rseg) = 0.0_rp
    v_int(inot,i_ring,i_rseg)     = 0.0_rp
    v_int_1_l(inot,i_ring,i_rseg) = 0.0_rp
    w_int(inot,i_ring,i_rseg)     = 0.0_rp
    w_int_1_l(inot,i_ring,i_rseg) = 0.0_rp
!
!--                Interpolation of the u-component:
    ii = int(  rbx(i_ring,i_rseg) * idx_dim)
    jj = int(  rby(i_ring,i_rseg) * idy_dim)
    kk = int(  rbz(i_ring,i_rseg) * idz_dim)

    if ( kk >= k_start .and. kk <= k_end) then
     kkp = kk - coord*(nz-1)

    xd = (rbx(i_ring,i_rseg) - float(ii)*dx_dim) * idx_dim
    yd = (rby(i_ring,i_rseg) - float(jj)*dy_dim) * idy_dim
    zd = (rbz(i_ring,i_rseg) - ((float(kk))*dz_dim)) * idz_dim

    aa = u(ii,jj,kkp)     * (1._rp - xd) + u(ii+1,jj,kkp)     * xd
    bb = u(ii,jj,kkp+1)   * (1._rp - xd) + u(ii+1,jj,kkp+1)   * xd
    cc = u(ii,jj+1,kkp)   * (1._rp - xd) + u(ii+1,jj+1,kkp)   * xd
    dd = u(ii,jj+1,kkp+1) * (1._rp - xd) + u(ii+1,jj+1,kkp+1) * xd

    u_int_l =  aa * (1._rp-yd) + cc * yd 
    u_int_u =  bb * (1._rp-yd) + dd * yd
    u_int_1_l(inot,i_ring,i_rseg) = u_int_l * (1._rp-zd) + &
                                    u_int_u * zd

!
!--                Interpolation of the v-component:
    aa = v(ii,jj,kkp)     * (1._rp - xd) + v(ii+1,jj,kkp)     * xd
    bb = v(ii,jj,kkp+1)   * (1._rp - xd) + v(ii+1,jj,kkp+1)   * xd
    cc = v(ii,jj+1,kkp)   * (1._rp - xd) + v(ii+1,jj+1,kkp)   * xd
    dd = v(ii,jj+1,kkp+1) * (1._rp - xd) + v(ii+1,jj+1,kkp+1) * xd

    v_int_l =  aa * (1._rp-yd) + cc * yd  
    v_int_u =  bb * (1._rp-yd) + dd * yd
    v_int_1_l(inot,i_ring,i_rseg) = v_int_l * (1._rp-zd) + &
                                    v_int_u * zd
!
!--                Interpolation of the w-component:
    aa = w(ii,jj,kkp)     * (1._rp - xd) + w(ii+1,jj,kkp)     * xd
    bb = w(ii,jj,kkp+1)   * (1._rp - xd) + w(ii+1,jj,kkp+1)   * xd
    cc = w(ii,jj+1,kkp)   * (1._rp - xd) + w(ii+1,jj+1,kkp)   * xd
    dd = w(ii,jj+1,kkp+1) * (1._rp - xd) + w(ii+1,jj+1,kkp+1) * xd

    w_int_l =  aa * (1._rp-yd) + cc * yd  
    w_int_u =  bb * (1._rp-yd) + dd * yd
    w_int_1_l(inot,i_ring,i_rseg) = w_int_l * (1._rp-zd) + &
                                    w_int_u * zd
    end if 
   end do
  end do
 !endif
end do

!dummy = MAXVAL(admr(:)%nsegs)
! Summation of u_int over all turbines:
call mpi_Allreduce( u_int_1_l, u_int, nturb * nrings_max *   &
                    nsegs_max, mpi_rp, mpi_sum, localComm, ierr)
call mpi_Allreduce( v_int_1_l, v_int, nturb * nrings_max *   &
                    nsegs_max, mpi_rp, mpi_sum, localComm, ierr)
call mpi_Allreduce( w_int_1_l, w_int, nturb * nrings_max *   &
                    nsegs_max, mpi_rp, mpi_sum, localComm, ierr)
!
!--       Loop over number of turbines:
do inot = 1, nturb
  cos_yaw = COS(admr(inot)%phi_yaw)
  sin_yaw = SIN(admr(inot)%phi_yaw)

!
!--          Loop over rings of each turbine:
  do i_ring = 1, admr(inot)%nrings
!
!--             Determine distance between each ring (center) and the hub:
!
!--             Loop over segments of each ring of each turbine:
   do i_rseg = 1, admr(inot)%nsegs(i_ring)
!
!--                Determine angle of ring segment towards zero degree angle of
!--                rotor system (at zero degree rotor direction vectors aligned
!--                with y-axis):
    phi_rotor = float(i_rseg) * 2.0_rp * pi / float(admr(inot)%nsegs(i_ring))
    cos_rot   = COS(phi_rotor)
    sin_rot   = SIN(phi_rotor)
!
!--                Coordinates of the single segments (center points):
    rbx(i_ring,i_rseg) = admr(inot)%loc(1) * z_i + admr(inot)%cur_r(i_ring) * cos_rot * sin_yaw
    rby(i_ring,i_rseg) = admr(inot)%loc(2) * z_i + admr(inot)%cur_r(i_ring) * cos_rot * cos_yaw
    rbz(i_ring,i_rseg) = admr(inot)%loc(3) * z_i + admr(inot)%cur_r(i_ring) * sin_rot


!
!--                !-----------------------------------------------------------!
!--                !-- Calculation of various angles and relative velocities --!
!--                !-----------------------------------------------------------!
!
!--                In the following the 3D-velocity field is projected its
!--                components perpedicular and parallel to the rotor area
!--                The calculation of forces will be done in the rotor-
!--                coordinates y' and z.
!--                The yaw angle will be reintroduced when the force is applied
!--                on the hydrodynamic equations
!
!--                Projection of the xy-velocities relative to the rotor area
!
!--                Velocity perpendicular to the rotor area:
    u_rot = u_int(inot,i_ring,i_rseg) *     cos_yaw +                     &
            v_int(inot,i_ring,i_rseg) * ( - sin_yaw)
            !+w_int(inot,i_ring,i_rseg)*rotn(3) << rotn always 0
!
!--                Projection of the 3D-velocity vector in the azimuthal
!--                direction:
    vtheta = (- sin_rot * sin_yaw) * u_int(inot,i_ring,i_rseg) +            & 
             (- sin_rot * cos_yaw) * v_int(inot,i_ring,i_rseg) +            &
                           cos_rot * w_int(inot,i_ring,i_rseg)
!
!--                Determination of the angle phi_rel between the rotor plane
!--                and the direction of the flow relative to the rotor:
    phi_rel = ATAN( u_rot /( admr(inot)%omega_rot * admr(inot)%cur_r(i_ring) - vtheta ) )

!
!--                Substraction of the local pitch angle to obtain the local
!--                angle of attack:
    alpha_attack = (phi_rel - admr(inot)%alpha_attack_in(i_ring)) * ( 360.0_rp / (2.0_rp*pi) )

!
!--                Determination of the magnitude of the flow velocity relative
!--                to the rotor:
    vrel = SQRT( u_rot**2 +                         &
                   ( admr(inot)%omega_rot * admr(inot)%cur_r(i_ring) -              &
                     vtheta )**2 )

!
!--                Determine index of current angle of attack, needed for
!--                finding the appropriate interpolated values of the lift and
!--                drag coefficients (-180.0 degrees = 0, +180.0 degrees = 36000,
!--                so one index every 0.01 degrees):
    iialpha = CEILING( ( alpha_attack + 180.0_rp )        &
                      * ( 1.0_rp / accu_cl_cd_tab ) )
!
!--                Determine index of current radial position, needed for
!--                finding the appropriate interpolated values of the lift and
!--                drag coefficients (one index every 0.1 m):
    iir = CEILING( admr(inot)%cur_r(i_ring) * 10.0_rp ) !AS: Change if different turbine is used!!!

!
!--                Read in interpolated values of the lift and drag coefficients
!--                for the current radial position and angle of attack:
    turb_cl = turb_cl_tab(iialpha,iir)
    turb_cd = turb_cd_tab(iialpha,iir)

!
!--                !-----------------------------------------------------!
!--                !-- Calculation of the forces                       --!
!--                !-----------------------------------------------------!
!
!--                Calculate the pre_factor for the thrust and torque forces:
    pre_factor = 0.5_rp * (vrel**2) * 3.0_rp *  &
                 admr(inot)%chord(i_ring) * admr(inot)%delta_r / admr(inot)%nsegs(i_ring)

!
!--                Calculate the thrust force (x-component of the total force)
!--                for each ring segment:
    thrust_seg(i_rseg) = pre_factor * ( turb_cl * COS(phi_rel) + turb_cd * SIN(phi_rel) )

!
!--                Determination of the second of the additional forces acting
!--                on the flow in the azimuthal direction: force vector as basis
!--                for torque (torque itself would be the vector product of the
!--                radius vector and the force vector):
    torque_seg = pre_factor * ( turb_cl * SIN(phi_rel) - turb_cd * COS(phi_rel) )
!
!--                Decomposition of the force vector into two parts:
!--                One acting along the y-direction and one acting along the
!--                z-direction of the rotor coordinate system:
    torque_seg_y(i_rseg) = -torque_seg * sin_rot
    torque_seg_z(i_rseg) =  torque_seg * cos_rot
!
!--                Add the segment thrust to the thrust of the whole rotor
    admr(inot)%thrust_rotor = admr(inot)%thrust_rotor + thrust_seg(i_rseg)                   
    admr(inot)%torque_total = admr(inot)%torque_total + (torque_seg * admr(inot)%cur_r(i_ring))

   end do !-- end of loop over ring segments

!
!--             Restore the forces into arrays containing all the segments of
!--             each ring:
  thrust_ring(i_ring,:)   = thrust_seg(:)
  torque_ring_y(i_ring,:) = torque_seg_y(:)
  torque_ring_z(i_ring,:) = torque_seg_z(:)


 end do  !-- end of loop over rings

!--          !-----------------------------------------------------------------!
!--          !--                  Regularization kernel                      --!
!--          !-- Smearing of the forces and interpolation to cartesian grid  --!
!--          !-----------------------------------------------------------------!
!--
!--          The aerodynamic blade forces need to be distributed smoothly on
!--          several mesh points in order to avoid singular behaviour
!
!--          Summation over sum of weighted forces. The weighting factor
!--          (calculated in user_init) includes information on the distance
!--          between the center of the grid cell and the rotor segment under
!--          consideration
!--
 do l = 1,admr(inot)%num_nodes
  i = admr(inot)%nodes(1,l)
  j = admr(inot)%nodes(2,l)
  k = admr(inot)%nodes(3,l)
  i_ring = admr(inot)%ringseg(1,l) 
  i_rseg = admr(inot)%ringseg(2,l)
  thrust(i,j,k)   = thrust(i,j,k)   + &
                    thrust_ring(i_ring,i_rseg)   * admr(inot)%smear(l)
  torque_y(i,j,k) = torque_y(i,j,k) + &
                    torque_ring_y(i_ring,i_rseg) * admr(inot)%smear(l) 
  torque_z(i,j,k) = torque_z(i,j,k) + &
                    torque_ring_z(i_ring,i_rseg) * admr(inot)%smear(l)

 end do

 do i = admr(inot)%i_hub - admr(inot)%i_smear,admr(inot)%i_hub + admr(inot)%i_smear
  do j = admr(inot)%j_hub - admr(inot)%j_smear,admr(inot)%j_hub + admr(inot)%j_smear
   do kkp = admr(inot)%min_k,admr(inot)%max_k
    if ( kkp >= k_start .and. kkp <= k_end) then
     k = kkp - coord*(nz-1)
!
!--                   Rotation of force components:
     rot_tend_x(i,j,k) = rot_tend_x(i,j,k)                  +  & 
                         thrust(i,j,k)  *admr(inot)%rotx(1) +  &
                         torque_y(i,j,k)*admr(inot)%roty(1) +  &
                         torque_z(i,j,k)*admr(inot)%rotz(1)
                                
     rot_tend_y(i,j,k) = rot_tend_y(i,j,k)                  +  &
                         thrust(i,j,k)  *admr(inot)%rotx(2) +  &
                         torque_y(i,j,k)*admr(inot)%roty(2) +  &
                         torque_z(i,j,k)*admr(inot)%rotz(2)                              
                              
     rot_tend_z(i,j,k) = rot_tend_z(i,j,k)                  +  &
                         thrust(i,j,k)  *admr(inot)%rotx(3) +  &
                         torque_y(i,j,k)*admr(inot)%roty(3) +  &
                         torque_z(i,j,k)*admr(inot)%rotz(3)                               


! Add forces to RHS. 
! Comment on dimensions: [rot_tend_*] = m/s^2  
!                        [RHSx]       = m^2/s^2 >> multiply by z_i >> get m^2/s^2

     RHSx(i,j,k) = RHSx(i,j,k) - rot_tend_x(i,j,k)*z_i           
     RHSy(i,j,k) = RHSy(i,j,k) - rot_tend_y(i,j,k)*z_i
     RHSz(i,j,k) = RHSz(i,j,k) - rot_tend_z(i,j,k)*z_i                       
    endif
   end do 
  end do
 enddo

! Buffer the write statements to the turbine files
! Written will be: 
! 1. Time[-], 2. Omega[1/s], 3. AeroTorque[m5/s2], 4. RotThrust [m/s2]
! 5. YawOrient
  if (modulo (jt_total,admr_tbase)==0) then
    admr_data_buf(1,inot,admr_count_buf)=total_time
    admr_data_buf(2,inot,admr_count_buf)=admr(inot)%omega_rot
    admr_data_buf(3,inot,admr_count_buf)=admr(inot)%torque_total
    admr_data_buf(4,inot,admr_count_buf)=admr(inot)%thrust_rotor
    admr_data_buf(5,inot,admr_count_buf)=(admr(inot)%phi_yaw)*180.0_rp/pi
  endif
end do      !-- end of loop over turbines
call admr_write_buffer(0)

! Update yaw angle
if ( adjust_yaw )  then
 do inot = 1,nturb
  admr_u  = 0.0_rp
  admr_v  = 0.0_rp
  admr_u2 = 0.0_rp
  admr_v2 = 0.0_rp

  ! AS: CHECK DIMENSIONS
  phi_yaw_old = admr(inot)%phi_yaw
  dist_1D = admr(inot)%rr/((dx_dim + dy_dim)*0.25_rp)

!   !Counts the number of indices upstream of the turbine
   x_int = -nint(dist_1D * COS(phi_yaw_old))
   y_int = -nint(dist_1D * SIN(phi_yaw_old))

   do l=1,admr(inot)%num_nodes
    i = admr(inot)%nodes(1,l)
    j = admr(inot)%nodes(2,l)
    k = admr(inot)%nodes(3,l)

    admr_u=admr_u + u(i+x_int,j+y_int,k)
    admr_v=admr_v + v(i+x_int,j+y_int,k)
   end do
   call mpi_allreduce(admr_u,admr_u2,1,mpi_rp,mpi_sum,localComm,ierr)
   call mpi_allreduce(admr_v,admr_v2,1,mpi_rp,mpi_sum,localComm,ierr)

   admr(inot)%phi_yaw = atan2d(admr_v2,admr_u2)
!  Time filter for turbine yaw 
   eps_admr = ((dt*z_i)/T_avg_admr)/(1.0_rp + (dt*z_i)/T_avg_admr)
   admr(inot)%phi_yaw = (1.0_rp-eps_admr)*phi_yaw_old + eps_admr * admr(inot)%phi_yaw
end do
call admr_nodes_func
end if

end SUBROUTINE admr_forcing

!*******************************************************************************
subroutine admr_write_buffer(admrflag)
!*******************************************************************************
use param, only : coord,path
use admr_mod, only : admr_count_buf,admr_max_buf,admr_data_buf,nturb
implicit none
integer,intent(in) :: admrflag ! 1 forces the data to be written
integer :: k,inot
character(4) :: admr_id

82 format(i4.4)
182 format(5(e16.7))
!106 FORMAT (F9.5,2X,F7.3,3X,F9.1,2X,F9.1,2X,F7.2,1X,F7.2)
!Write the buffered turbine data to file
if ((admr_count_buf.eq.admr_max_buf).or. (admrflag.eq.1)) then
 if(coord.eq.0) then ! Only do the writing at master
  do inot = 1,nturb
   write(admr_id,82)  inot
   open(1000+inot,file=trim(path)//'admr-output/admr_'//trim(admr_id)//'.dat',status='unknown',form='formatted',position='append')
   do k=1,admr_count_buf
    write(1000+inot,182) admr_data_buf(1:5,inot,k)
   enddo
   close(1000+inot)
  enddo
 endif
 admr_count_buf=0 ! Reset the counter
endif

end subroutine admr_write_buffer
#endif ADMR 
