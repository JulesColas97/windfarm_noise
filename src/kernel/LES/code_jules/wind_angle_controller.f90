#ifdef CORIOLIS 
subroutine wind_angle_controller_init
use param, only : nz, z_i, z_wind, coord, localComm, mpi_rp, ierr, rank_hub, fixed_height,dz, height_index,initu,rp
use param, only: ug, vg
use sim_param, only : alpha_wind, alpha_wind_old, phi_new, phi_0, time1, time2
use sim_param, only : error_der, error_pre, error_int, error_pro, omega_eff
use grid_defs, only : gridz
use mpi, only: mpi_comm_world,mpi_integer
implicit none
integer:: k, buffernum,count_level_per_proc

if(initu == .false.) then
 error_int     =0.0_rp
 error_der     =0.0_rp
 error_pro     =0.0_rp
 error_pre     =0.0_rp
 alpha_wind    =atan2(vg, ug)
 alpha_wind_old=atan2(vg, ug)
 phi_0    =0.0_rp
 omega_eff=0.0_rp
 phi_new  =0.0_rp
 time1    =0.0_rp
 time2    =0.0_rp
endif

rank_hub = 0
buffernum = rank_hub
if ( coord == 0) then
 count_level_per_proc = 0
 do k=1,nz-1
  count_level_per_proc = count_level_per_proc + 1
 end do

 ! round in the same way as the below (floor / ceil )
 !rank_hub = nint(z_wind/(dz*z_i*real(...)-0.5)
 rank_hub = nint(z_wind/(dz*z_i*count_level_per_proc))-1
endif

call mpi_bcast(rank_hub, 1, mpi_integer, 0, localComm, ierr)

do k=1,nz-1
 if(z_wind > (gridz(k)*z_i) .and. z_wind <= (gridz(k+1)*z_i)) then
  fixed_height = .true.
  height_index = k
   if (coord .ne. rank_hub) then
    write(*,*) '! in wind_angle_controller.f90: Rank not defined correctly.'
    call mpi_barrier(mpi_comm_world, ierr)
    call mpi_finalize(ierr)
    stop
   end if
   
 end if
end do

end subroutine wind_angle_controller_init
#endif

#ifdef CORIOLIS 
!THIS HAS BEEN TESTED COMPLETELY
!DO NOT CHANGE ANYTHING IF YOU DO NOT UNDERSTAND IT COMPLETELY
!IF SOMETHING LOOKS WRONG, IT IS BY DESIGN NOT BY OVERSIGHT!!!!!
subroutine wind_angle_controller
use param, only : path, jt_total, nx, ny, nz, rp, Kdw, Kpw, Kiw, dt, u_hub, v_hub
use param, only : localComm, mpi_rp, ierr, wind_angle
use param, only : z_i, z_wind, total_time, inxny
use param, only : fixed_height, height_index,rank_hub
use grid_defs, only : gridz
use sim_param, only : u, v, alpha_wind, alpha_wind_old, phi_new, phi_0, time1, time2,u_dim
use sim_param, only : error_der, error_pre, error_int, error_pro, omega_eff
#ifdef BAROCLINIC
use param, only : geo_force_x, geo_force_y, pi
use sim_param, only : geo_force_x_rot, geo_force_y_rot
#endif BAROCLINIC

implicit none
integer:: k 
real(rp), dimension(nz) :: angle, u_aver_z, v_aver_z
real(rp) :: w1, w2, ts
character(len=64):: folder

character*8 :: ipfi
82 format(i8.8)
write(ipfi,82) jt_total
time2 = total_time
phi_0 = wind_angle

!Uncomment this block if you want to plot wind angle vs height 
!This can be used to check the performance of controller 
do k=1,nz-1
   v_aver_z(k) = 0._rp
   u_aver_z(k) = 0._rp
   u_aver_z(k) = sum(u(1:nx,1:ny,k))*inxny
   v_aver_z(k) = sum(v(1:nx,1:ny,k))*inxny
   angle(k) = atan2d(v_aver_z(k),u_aver_z(k))
enddo

folder = trim(path)//'output'
if(modulo(jt_total,5000) == 0) then
   call write_1Dfield(2,angle(1:nz-1),'angle',folder,0,nz-1)
endif
  
if(fixed_height) then
 k  = height_index

 u_aver_z(k)   = sum(u(1:nx,1:ny,k))*inxny
 v_aver_z(k)   = sum(v(1:nx,1:ny,k))*inxny
 u_aver_z(k+1) = sum(u(1:nx,1:ny,k+1))*inxny
 v_aver_z(k+1) = sum(v(1:nx,1:ny,k+1))*inxny

 !Calculate the mean wind angle at the hub height
  w1 = abs(z_wind - gridz(k)*z_i)
  w2 = abs(z_wind - gridz(k+1)*z_i)

  v_hub = (w2*v_aver_z(k) + w1*v_aver_z(k+1)) / (w1 + w2)
  u_hub = (w2*u_aver_z(k) + w1*u_aver_z(k+1)) / (w1 + w2)
  phi_new = atan2(v_hub,u_hub)

 if(modulo(jt_total,2500) == 0) then
  open(261,file=trim(path)//'output/alpha_'//trim(ipfi)//'.dat',status='unknown',form='formatted')
  write(261,*) total_time,phi_new,alpha_wind
  close(261)
 endif

  ts = time2 - time1

  !Proportional part of the error 
  error_pro = (phi_new - phi_0)

  !Integral part of the error 
  error_int = error_int + 0.5_rp* (error_pro+error_pre) * dt

  !Derivative part of the error 
  error_der = (error_pro - error_pre)/(dt) 

  omega_eff  = Kpw * error_pro + Kiw * error_int + Kdw * error_der
  alpha_wind = alpha_wind_old - 2._rp * omega_eff * ts

  error_pre = error_pro
endif

call mpi_bcast(omega_eff, 1, mpi_rp, rank_hub, localComm, ierr)
call mpi_bcast(alpha_wind,1, mpi_rp, rank_hub, localComm, ierr)
call mpi_bcast(error_pre, 1, mpi_rp, rank_hub, localComm, ierr)
call mpi_bcast(error_int, 1, mpi_rp, rank_hub, localComm, ierr)

! Update the forcing term for baroclinic ABL with new rotation angle. Doing this here prevents
! needlessly doing it locally multiple times in different files.
#ifdef BAROCLINIC
do k=1,nz
   geo_force_y_rot(k) = geo_force_x(k) * sind(180.0 * alpha_wind/pi) + geo_force_y(k) * cosd(180.0 * alpha_wind/pi)
   geo_force_x_rot(k) = geo_force_x(k) * cosd(180.0 * alpha_wind/pi) - geo_force_y(k) * sind(180.0 * alpha_wind/pi)
enddo
#endif BAROCLINIC

!Save previous time 
!Should I use dimensional time????
time1 = time2 

end subroutine wind_angle_controller
#endif 
