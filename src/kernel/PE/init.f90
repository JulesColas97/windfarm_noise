module init

implicit none

public :: init_domain
public :: init_mesh
public ::  init_meanflow
public ::  init_meanflow_arbitrary

contains

!! @brief define the parameters and the mesh of the 2D simulation coreseponding to the kk angle
!!
!! @details the function first determine the dimension of the 2D mesh from the total dimension of the simulation
!! and from the angle, The it define the corresponding mesh
!!
!!
!! @param [in]      simu: parameters of the global simulation
!! @param [in]      theta: angle
!! @param [in]      dout: coarse deltax
!! @param[inout]    mesh: the 2D mesh coresponding to theta propagation
!! @param[inout]    dom: parameters of the 2D mesh coresponding to the kk angles
!! @return
subroutine init_domain(theta,simu,mesh,dom,dout)
    use types
    implicit none

    type(simulation), intent(in)                            :: simu
    real, intent(in)                                        :: theta
    type(mesh_2D), intent(inout)                            :: mesh
    type(domain), intent(inout)                             :: dom

    real, intent(in)                                        :: dout

    integer                                                 :: ii
    real                                                    :: pi = 3.141592653589793, theta_rad

    ! conversion in rad
    theta_rad = theta*pi/180
    if (theta_rad > pi) then 
        theta_rad = theta_rad - 2*pi
    end if 

    ! define size 2D domain from 3D dimension and angle
    if (theta_rad == 0) then
        dom%Lx = simu%Lx1
    elseif (theta_rad == -pi) then
        dom%Lx = simu%Lx2
    elseif (theta_rad == pi/2 .or. theta_rad == -pi/2) then
        dom%Lx = simu%Ly
    elseif (-pi/2<theta_rad .and.  theta_rad<pi/2) then
        dom%Lx = min(simu%Lx1/abs(cos(theta_rad)),simu%Ly/abs(sin(theta_rad)))
    else
        dom%Lx = min(simu%Lx2/abs(cos(theta_rad)),simu%Ly/abs(sin(theta_rad)))
    end if

    ! set dimension of the recording grid
    dom%Lx_coarse = dom%Lx
    dom%Lz_coarse = simu%Lz
    dom%nx_coarse = floor(dom%Lx_coarse/dout)
    dom%nz_coarse =  floor(dom%Lz_coarse/dout)

    ! create recording grid
    allocate(mesh%x_coarse(dom%nx_coarse))
    allocate(mesh%z_coarse(dom%nz_coarse))
    ! init coarse mesh
    do ii = 1, dom%nx_coarse
        mesh%x_coarse(ii) = dout*(ii -1)
    enddo

    do ii = 1, dom%nz_coarse
        mesh%z_coarse(ii) = dout*(ii -1)
    enddo
end subroutine init_domain


!! @brief Create the computing mesh for a given frequency and angle
!!
!! @details define a dx given the frequency and 'cfl', create the coresponding mesh
!!
!!
!! @param [in]      freq: frequency of the calculation
!! @param [in]      c0 : sound sped
!! @param [in]      simu: parameters of the global simulation
!! @param[inout]    mesh: the 2D mesh coresponding to the computed angle
!! @param[inout]    dom: parameters of the 2D mesh coresponding to the angle computed
!! @return
subroutine init_mesh(freq,c0,simu,mesh,dom,pml,dout)
    use types
    implicit none

    type(simulation), intent(in)                            :: simu
    real, intent(in)                                        :: freq,c0
    type(mesh_2D), intent(inout)                            :: mesh
    type(domain), intent(inout)                             :: dom
    type(perfectly_matched_layer), intent(in)               :: pml
    real, intent(in)                                        :: dout

    integer                                                 :: ii

    ! define solving mesh quantities
    dom%dx = simu%cfl*c0/freq
    if (dom%dx> dout) then
      dom%dx = dout
    endif
    dom%dz = dom%dx

    ! add pml on top of the domain
    dom%Lz = simu%Lz + pml%size*dom%dz
    dom%nx = floor(dom%Lx/dom%dx)+1
    dom%nz = floor(dom%Lz/dom%dz)+1

    ! init mesh
    allocate(mesh%x(dom%nx))
    allocate(mesh%z(dom%nz))
    do ii = 1, dom%nx
        mesh%x(ii) = dom%dx*(ii -1)
    enddo

    do ii = 1, dom%nz
        mesh%z(ii) = dom%dz*(ii -1)
    enddo
end subroutine init_mesh

!! @brief initilization of the 2D mean flow for the propagation
!!
!! @details The function either take the 3D mesh as input and interpolate in the 2D mesh for the propagation
!! or it read an already interpolated 2D mesh and check that the dimension c orrespond to the parameters of the simulation.
!! Note that mean flow is interpolated on a coarser grid coresponding to dout in the x direction

!! @param [in]      kk : number of the angle that is evalutated
!! @param [in]      simu: parameters of the global simulation
!! @param [in]      mesh: the 2D mesh coresponding to the kkth angles
!! @param [in]      nbout: the number of grid step bewteen each recording of the results
!! @param [in]      nx : number of point in the x direction
!! @param [in]      nz : number of point in the z direction
!! @param [in]      c0: the sound speed
!! @param[inout]    flow3D: the 3D flow before interpolation
!! @param[inout]    flow2D: the interpolated 2D flow inside the mesh
!! @param[inout]    dom: parameters of the 2D mesh coresponding to the kk angles
!! @return
subroutine init_meanflow_arbitrary(kk, simu,mesh,flow3D, flow2D, c0)
        use types
        use interpo
        use read_write
        implicit none
        ! input variables
        integer, intent(in)                                 :: kk
        type(simulation), intent(in)                        :: simu
        type(mesh_2D), intent(inout)                        :: mesh
        type(meanflow3D), intent(inout)                     :: flow3D
        type(meanflow2D), intent(inout)                     :: flow2D
        real, intent(in)                                    :: c0

        !local variables
        integer                                             :: ii
        integer                                             :: nx, nz
        real                                                :: start, finish
        real                                                :: theta_rad
        real                                                :: pi =3.141592653589793
        ! junk
        real                                                :: dx, dz

        theta_rad = simu%theta(kk)*pi/180
        ! Case where the flow is read from external h5 file and interpolated
        print*, 'init mean flow for arbitrary mach numbers'
        nx = len(mesh%x_coarse)
        nz = len(mesh%z_coarse)

        allocate(flow2D%epsilon(nx,ny))
        allocate(flow2D%m(nx,ny))
        allocate(flow2D%gamma(nx,ny))
        allocate(flow2D%tau(nx,ny))
        allocate(flow2D%x(nx,nz))
        allocate(flow2D%z(nx,nz))


        do ii = 1,nz
            flow2D%x(:,ii) = mesh%x_coarse(1:nx)
        enddo
        do ii = 1,ceiling(nz)
            flow2D%z(ii,:) = mesh%z_coarse(1:nz)
        enddo


        if (simu%external_flow) then
            print*,'Interpolation of 3D lES'
            ! interpolation of 3D mean flow on 2D  coarse meshgrid
            ! interpolation of the
            CALL CPU_TIME(start)
            call  interpolate2(flow3D%u,flow3D%v,flow3D%w,flow3D%x,flow3D%y,flow3D%z,flow2D%u,flow2D%v, &
                            & flow3D%Xs,flow3D%Ys, flow3D%Zs,dom%x_coarse,dom%z_coarse,theta_rad)

            CALL CPU_TIME(finish)
            print*, 'interpolation time :', finish -start, 's'
            flow2D%epsilon              = ((c0/c0)**2) - 1
            flow2D%m                    = flow2D%u/c0
            flow2D%gamma                = 1/sqrt(1- (flow2D%m**2))
            flow2D%tau                  = flow2D%m * (flow2D%gamma**2) * (sqrt(1+flow2D%epsilon) - flow2D%m)

        ! uniform mean flow u0
        else if (simu%uniform) then
          print*,'Uniform flow, u0 =', simu%u0

          flow2D%epsilon              = ((c0/c0)**2) - 1
          flow2D%m                    = simu%u0/c0
          flow2D%gamma                = 1/sqrt(1- (flow2D%m**2))
          flow2D%tau                  = flow2D%m * (flow2D%gamma**2) * (sqrt(1+flow2D%epsilon) - flow2D%m)
        ! Logarithmic profile  hard coded
        !TODO give parameter from input file
        else if (simu%logarithmic) then
          print*,'Logarithmic flow profile, u0 =', simu%u0
          allocate(flow2D%u(nx,ny))
          do ii = 1,ceiling((1.*nx)/nbout)
        !    flow2D%u(ii,:) = (simu%u0/6.9088)*log((mesh%z+0.1)/0.1)
           flow2D%u(ii,:) = (simu%u0/0.4)*log((mesh%z+0.1)/0.1)
          enddo
          flow2D%epsilon              = ((c0/c0)**2) - 1
          flow2D%m                    = flow2D%u/c0
          flow2D%gamma                = 1/sqrt(1- (flow2D%m**2))
          flow2D%tau                  = flow2D%m * (flow2D%gamma**2) * (sqrt(1+flow2D%epsilon) - flow2D%m)

        else
            print*,'No flow'
            ! flow equal 0
            flow2D%epsilon              = 0.
            flow2D%m                    = 0.
            flow2D%gamma                = 1.
            flow2D%tau                  = 0.
        endif
end subroutine init_meanflow_arbitrary


!! @param [in]      kk : number of the angle that is evalutated
!! @param [in]      simu: parameters of the global simulation
!! @param [in]      mesh: the 2D mesh coresponding to the kkth angles
!! @param [in]      nbout: the number of grid step bewteen each recording of the results
!! @param [in]      nx : number of point in the x direction
!! @param [in]      nz : number of point in the z direction
!! @param [in]      c0: the sound speed
!! @param[inout]    flow3D: the 3D flow before interpolation
!! @param[inout]    flow2D: the interpolated 2D flow inside the mesh
!! @param[inout]    dom: parameters of the 2D mesh coresponding to the kk angles
!! @return
subroutine init_meanflow(kk, simu,mesh,flow3D,tinput, flow2D, nbout, nx, nz, c0)
        use types
        use interpo
        use read_write
        implicit none
        ! input variables
        integer, intent(in)                                 :: kk
        type(simulation), intent(in)                        :: simu
        type(mesh_2D), intent(inout)                        :: mesh
        type(meanflow3D), intent(inout)                     :: flow3D
        type(meanflow2D), intent(inout)                     :: flow2D
        integer, intent(inout)                              :: nbout, nx, nz
        real, intent(in)                                    :: c0
        type(twente_input), intent(in)                      :: tinput

        !local variables
        integer                                             :: ii
        real                                                :: start, finish
        real                                                :: theta_rad
        real                                                :: pi =3.141592653589793
        ! junk
        real                                                :: dx, dz

        theta_rad = simu%theta(kk)*pi/180
        print*, 'init mean flow for classic wape'
        if (simu%external_flow) then
            print*,'Interpolation of 3D lES'
            if (simu%interpolation) then
                ! interpolation of 3D mean flow on 2D meshgrid
                ! definition of the 2D grid on which the flow woill be interpolated
                allocate(flow2D%x(ceiling((1.*nx)/nbout),nz))
                allocate(flow2D%z(ceiling((1.*nx)/nbout),nz))
                do ii = 1,nz
                    flow2D%x(:,ii) = mesh%x(1:nx:nbout)
                enddo
                do ii = 1,ceiling((1.*nx)/nbout)
                    flow2D%z(ii,:) = mesh%z(1:nz)
                enddo
                ! interpolation of the
                CALL CPU_TIME(start)
                call  interpolate2(flow3D%u,flow3D%v,flow3D%w,flow3D%x,flow3D%y,flow3D%z,flow2D%u,flow2D%v, &
                                & flow3D%Xs,flow3D%Ys, flow3D%Zs,flow2D%x,flow2D%z,theta_rad)
                CALL CPU_TIME(finish)
                print*, 'interpolation time :', finish -start, 's'
            else
                ! read already interpolatred flow
                call read_2Dflow(simu%theta(kk),flow2D%u,flow2D%v, flow2D%epsilon, &
                                    & mesh%x, mesh%z, nbout,nx, nz, dx, dz)
                ! check flow and mesh are the same */
            endif
            allocate(flow2D%epsilon(ceiling((1.*nx)/nbout),nz))
            flow2D%epsilon                 = (c0/(c0 + flow2D%u))**2-1


        else if (simu%uniform) then
            print*,'Uniform flow, u0 =', simu%u0
            allocate(flow2D%epsilon(ceiling((1.*nx)/nbout),nz))
            flow2D%epsilon                 = (c0/(c0 + simu%u0))**2-1


        else if (simu%logarithmic) then
            print*,'Logarithmic flow profile, u0 =', simu%u0
            allocate(flow2D%u(ceiling((1.*nx)/nbout),nz))
            allocate(flow2D%epsilon(ceiling((1.*nx)/nbout),nz))
            do ii = 1,ceiling((1.*nx)/nbout)
                ! flow2D%u(ii,:) = (simu%u0/6.9088)*log((mesh%z+0.1)/0.1)
                flow2D%u(ii,:) = (simu%u0/0.4)*log((mesh%z+0.1)/0.1)
            enddo
            flow2D%epsilon                 = (c0/(c0 + flow2D%u))**2-1
        else
            ! flow equal 0
            print*,'No flow'
            allocate(flow2D%epsilon(ceiling((1.*nx)/nbout),nz))
            flow2D%epsilon              = 0.
        endif
end subroutine init_meanflow

end module init
