
!! Types module for definition of all the customs types of the progra=
module types

    ! Simulation global parameters
    type simulation
        integer                             :: nb_theta                         ! number of angles
        integer                             :: nb_freq                          ! number of frequency
        real, allocatable,dimension (:)     :: theta                            ! angles
        real, allocatable,dimension (:)     :: frequencies                      ! frequencies
        logical                             :: arbitrary                        ! type of PE (arbuirary Mach or classic)
        logical                             :: arbitrary_new


        real                                :: Lx1,Lx2,Ly,Lz                    ! 3D domain dimension
        real                                :: dx                               ! grid step
        real                                :: cfl                              ! ratio betwen frequency and dx

        logical                             :: external_flow                    ! read from h5
        logical                             :: interpolation                    ! interpolated from 3D field
        logical                             :: uniform                          ! uniform flow
        logical                             :: logarithmic                      ! logarithmic profil
        real                                :: u0                               ! unifrom value
    end type simulation

    ! 2D mesh parameters
    type domain
        integer                             :: nx                               ! number of point in x
        integer                             :: nz                               ! number of point in z
        real                                :: dx                               ! grid step
        real                                :: dz                               ! grid step
        real                                :: Lx                               ! x dimension
        real                                :: Lz                               ! z dimension
        real                                :: x0                               ! source position
        real                                :: z0                               ! source position
        ! parameter for recording mesh
        integer                             :: nx_coarse
        integer                             :: nz_coarse
        real                                :: Lx_coarse
        real                                :: Lz_coarse
    end type domain

    ! 2D mesh vector
    type mesh_2D
        real, allocatable, dimension(:)     :: x, x_coarse                      ! x vector
        real, allocatable, dimension(:)     :: z, z_coarse                      ! z vector
    end type mesh_2D

    !impedance parameters
    type impedance
        real                                :: sigmae                           ! effective porosity
        real                                :: alphae                           ! effective resistivity
        complex                             :: Zs                               ! impedance
        complex                             :: beta                             ! admitance
        logical                             :: rigid                            ! if true beta=0
    end type impedance

    ! 3D fields for the mean flow
    type meanflow3D
        logical                             :: ext                              ! external file
        character(len=100)                  :: fdir                             ! directory path
        character(len=100)                  :: fname                            ! prefix  of hdf5
        real, allocatable,dimension(:,:,:)  :: u, v, w                          ! 3D fields
        real, allocatable,dimension(:)      :: x, y, z                          ! coordinates
        real                                :: Xs, Ys, Zs, theta                ! WT position in flow
    end type meanflow3D

    ! parameter from twente LES simulation
    type twente_input
        real                                :: delta                            ! adimensional parameter
        real                                :: z_i                              ! adimensional length parameter
        real                                :: Lx, Ly, Lz                       ! adimensional size of the domain
        real                                :: posx,posy,posz                   ! adimensional position of the source
        real                                :: ratio                            ! to redimension the flow
    end type twente_input

    ! 2D fields for the mean flow
    type meanflow2D
        logical                             :: ext
        character(len=100)                  :: fdir                             ! directory for 2D data
        character(len=100)                  :: fname                            ! namle of h5 file
        real, allocatable,dimension(:,:)    :: u, v, epsilon, gamma, tau, m     ! flow quantities
        real, allocatable,dimension(:,:)    :: x, z                             ! coarse grid
        real                                :: Xs, Ys, Zs, theta                ! WT position
    end type meanflow2D

    ! state varibales
    type variables_0
        real                                :: rho                              !masse volumique
        real                                :: c                                !vitesse du son
        real                                :: gamma                            !indice adiabatique
    end type variables_0

    ! frequency quantities
    type frequency
        real                                :: freq                             ! frequency
        real                                :: omega                            ! pulsation
        real                                :: lambda                           ! wave length
        real                                :: k0                               ! wave number
    end type frequency

    ! surce term
    type source
        real                                :: pos_x                            ! position of acoustic source
        real                                :: pos_z
        real                                :: amp                              ! amplitude
    end type source

    ! pml parameters
    type perfectly_matched_layer
        integer                             :: size                             ! nb of points of the layer
        real                                :: param                            ! strength
        real                                :: n                                ! power
        real, allocatable, dimension(:)     :: sigma                            ! coeficient of the layer (along z)
        real, allocatable, dimension(:)     :: sigma_de                         ! first derivative
    end type perfectly_matched_layer

    ! vector and matrix for the resolution
    type solve
        real, allocatable,dimension(:,:)      :: h10, h12, h20, h22             ! vector for arbitrati M

        ! vector for classic WAPE
        complex                             :: nu
        complex                             :: mu
        complex, allocatable,dimension(:)   :: nu_vect
        complex, allocatable,dimension(:)   :: mu_vect

        ! vector for tridiagonal matrix split step
        complex, allocatable,dimension(:)   :: a_vect
        complex, allocatable,dimension(:)   :: b_vect
        complex, allocatable,dimension(:)   :: c_vect

        complex, allocatable,dimension(:)   :: d_vect
        complex, allocatable,dimension(:)   :: e_vect
        complex, allocatable,dimension(:)   :: f_vect
    end type solve

    type solution
        complex, allocatable,dimension(:)   :: u
        complex, allocatable,dimension(:)   :: phi, phi_m1, phi_p1, p
        complex, allocatable,dimension(:,:) :: phi_vect
        complex, allocatable,dimension(:,:) :: p_vect
        real, allocatable,dimension(:)      :: deltaL
        real, allocatable,dimension(:,:)    :: deltaL_vect
        real, allocatable,dimension(:,:)    :: deltaL_int
        real,allocatable,dimension(:,:)     :: receiver
        real,allocatable,dimension(:,:)     :: receiver_int
        integer                             :: nbout
    end type solution

    type output_parameters
        real                                :: dout
        integer                             :: nbout
        integer                             :: nb_receiver
        real, allocatable,dimension(:)      :: heights
        logical                             :: side
        logical                             :: top
    end type output_parameters

end module types
