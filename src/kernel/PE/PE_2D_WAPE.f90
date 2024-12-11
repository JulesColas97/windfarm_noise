program PE_2D_WAPE

use types
use splines
use utils
use interpo
use read_write
use hdf5
use wape
use init
!$ use OMP_LIB


implicit none
character(len=40)                              :: case_name
type(simulation)                                :: simu
type(output_parameters)                         :: output
type(domain)                                    :: dom
type(frequency)                                 :: freq
type(impedance)                                 :: imp
type(variables_0)                               :: var0
type(perfectly_matched_layer)                   :: pml
type(source)                                    :: src
type(mesh_2D)                                   :: mesh
type(solve)                                     :: M
type(solution)                                  :: sol
type(meanflow3D)                                :: flow3D
type(twente_input)                              :: tinput
type(meanflow2D)                                :: flow2D

real, allocatable,dimension(:)                  :: frequencies
real, allocatable,dimension(:)                  :: angles
integer                                         :: ii,kk,jj

! read parameters from input config file
!-------------------------------------------------------------------------------
print*, 'read parameters'
call read_param(case_name,var0,simu,src,imp,output,pml,flow3D,tinput)

! read flow field from hdf5 file
!-------------------------------------------------------------------------------
if (simu%external_flow) then
  if (simu%interpolation) then
      call read_h5flow3D(flow3D%fdir,flow3D%fname,flow3D, tinput)
  endif
endif


! if continuation start loop on the first not coputed angle
!TODO

! start loop on angles
do kk = 1,size(simu%theta)
    print*, 'init propa angle : ', simu%theta(kk)
    ! init 2D domain and coarse mesh for recording
    call init_domain(simu%theta(kk), simu, mesh, dom,output%dout)

    ! initiliaze output file
    call init_file(case_name,simu%theta(kk), output%heights)
    ! write mesh in h5 file
    call write_mesh(case_name,simu%theta(kk),mesh%x_coarse,mesh%z_coarse)

    ! write flow 
    !-------------------------------------------------------------------------
    ! initialize mesh for propa at one frequency
    ! call init_mesh(1.0,var0%c,simu, mesh, dom,pml,output%dout)

    ! ! nbout for the meanflow interpolation over a coarse grid
    ! output%nbout = ceiling(output%dout/dom%dx)



    ! interpolation of the meanflow over the coarse grid
    if (simu%arbitrary .or. simu%arbitrary_new) then
      ! interpolation and calculation of value for arbitrary mach number
      call init_meanflow_arbitrary(kk,simu,mesh,flow3D,flow2D,var0%c)
    else
      ! classic effective sound speed approach
      call init_meanflow(kk,simu,mesh,flow3D,tinput,flow2D,output%nbout,dom%nx,dom%nz,var0%c)
    endif
    ! write 2D inetrpolated meanflow in h5 file
    call  write_2Dflow(flow2D%u,flow2D%v, flow2D%epsilon, mesh%x_coarse, mesh%z_coarse)
    print*, 'mean flow written'


    ! deallocate mesh
    deallocate(mesh%x)
    deallocate(mesh%z)

    !deallocate flow
    if (simu%external_flow) then
    deallocate(flow2D%u)
    deallocate(flow2D%v)
        if (simu%interpolation) then
            deallocate(flow2D%x)
            deallocate(flow2D%z)
        endif
    endif
    deallocate(flow2D%epsilon)
    if (simu%arbitrary .or. simu%arbitrary_new) then
      deallocate(flow2D%m)
      deallocate(flow2D%gamma)
      deallocate(flow2D%tau)
    endif
    !-------------------------------------------------------------------------


    ! loop over the frequencies
    !-------------------------------------------------------------------------
    do ii = 1, simu%nb_freq
        print*, 'start propa for frequency: ', simu%frequencies(ii)

        ! initialize mesh for propa at one frequency
        call init_mesh(simu%frequencies(ii),var0%c,simu, mesh, dom,pml,output%dout)

        ! nbout for the meanflow interpolation over a coarse grid
        output%nbout = ceiling(output%dout/dom%dx)

        ! interpolation of the meanflow over the coarse grid
        print*,  'init meanflow ...'
        if (simu%arbitrary .or. simu%arbitrary_new) then
          ! interpolation and calculation of value for arbitrary mach number
          call init_meanflow_arbitrary(kk,simu,mesh,flow3D,tinput,flow2D,output%nbout,dom%nx,dom%nz,var0%c)
        else
          ! classic effective sound speed approach
          call init_meanflow(kk,simu,mesh,flow3D,tinput,flow2D,output%nbout,dom%nx,dom%nz,var0%c)
        endif
        print*, 'done.'

        ! ! write 2D inetrpolated meanflow in h5 file
        ! call  write_2Dflow(flow2D%u,flow2D%v, flow2D%epsilon, flow2D%x, flow2D%z)


        ! Init vector for solution
        !-----------------------------------------------------------------------
        if (output%side) then
            allocate(sol%deltaL_vect(dom%nx,dom%nz_coarse))
            allocate(sol%deltaL_int(dom%nx_coarse,dom%nz_coarse))
            sol%deltaL_vect = 0
            sol%deltaL_int = 0
        endif
        allocate(sol%receiver(dom%nx,output%nb_receiver))
        allocate(sol%receiver_int(dom%nx_coarse,output%nb_receiver))
        sol%receiver  = 0
        sol%receiver_int = 0

        ! perform PE resolution
        if (simu%arbitrary) then
          ! zrbitrary mach number
          call  propa_1angle_1freq_arbitrary(var0,simu%frequencies(ii),src,imp,dom,pml,flow2D, mesh, sol, output)
        else if (simu%arbitrary_new) then 
          call  propa_1angle_1freq_arbitrary_new(var0,simu%frequencies(ii),src,imp,dom,pml,flow2D, mesh, sol, output)
        else
          ! classic wape
          call  propa_1angle_1freq(var0,simu%frequencies(ii),src,imp,dom,pml,flow2D, mesh, sol, output)
        end if

        ! Interpolation over X dimension for carto
        if (output%side) then
          do jj = 1,dom%nz_coarse
            call interplog(mesh%x,sol%deltaL_vect(:,jj),mesh%x_coarse,sol%deltaL_int(:,jj))
            !call interpo_1D(mesh%x,sol%deltaL_vect(:,jj),mesh%x_coarse,sol%deltaL_int(:,jj))
            !sol%deltaL_int(:,jj) = spline3(mesh%x, sol%deltaL_vect(:,jj), mesh%x_coarse)
          enddo
        endif

        ! Interpolation over X dimension for receivers
        do jj = 1,output%nb_receiver
        ! print*, 'test1'
        ! print*, size(mesh%x),size(sol%receiver(:,jj)),size(mesh%x_coarse),size(sol%receiver_int(:,jj))
          call interplog(mesh%x,sol%receiver(:,jj),mesh%x_coarse,sol%receiver_int(:,jj))
          !call interpo_1D(mesh%x,sol%receiver(:,jj),mesh%x_coarse,sol%receiver_int(:,jj))
          !sol%receiver_int(:,jj) = spline3(mesh%x, sol%receiver(:,jj), mesh%x_coarse)
        enddo

        ! writing solution
        print*, 'writing solution'
        call write_solution(case_name,simu%theta(kk),simu%frequencies(ii),sol%deltaL_int)

         call write_receiver(case_name,simu%theta(kk),simu%frequencies(ii),dom%nx_coarse,&
                     & output%heights, sol%receiver_int)
        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        !call write_receiver(case_name,simu%theta(kk),simu%frequencies(ii),dom%nx,&
              !      & output%heights, sol%receiver)

        !dealocate solution
        if (output%side) then
          deallocate(sol%deltaL_vect)
          deallocate(sol%deltaL_int)
        endif

        deallocate(sol%receiver)
        deallocate(sol%receiver_int)

        ! deallocate mesh
        deallocate(mesh%x)
        deallocate(mesh%z)

        !deallocate flow
        if (simu%external_flow) then
        deallocate(flow2D%u)
        deallocate(flow2D%v)
            if (simu%interpolation) then
                deallocate(flow2D%x)
                deallocate(flow2D%z)
            endif
        endif
        deallocate(flow2D%epsilon)
        if (simu%arbitrary .or. simu%arbitrary_new) then
          deallocate(flow2D%m)
          deallocate(flow2D%gamma)
          deallocate(flow2D%tau)
        endif

    enddo
    ! deallocate mesh
    deallocate(mesh%x_coarse)
    deallocate(mesh%z_coarse)
    ! if delta L true then compute delta L over all frequencies
enddo
end program PE_2D_WAPE
