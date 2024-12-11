module wape
implicit none
public :: propa_1angle_1freq
public :: propa_1angle_1freq_arbitrary
public :: propa_1angle_1freq_arbitrary2
contains

subroutine propa_1angle_1freq_arbitrary(var0,f,src,imp,dom,pml,flow2D,mesh,sol,output)
    use types
    use utils
    use interpo
    use read_write
    use hdf5
    use splines
    implicit none
    ! input parameters
    type(domain), intent(inout)                             :: dom
    type(real), intent(in)                                  :: f
    type(impedance), intent(inout)                          :: imp
    type(variables_0), intent(inout)                        :: var0
    type(perfectly_matched_layer),intent(inout)             :: pml
    type(source), intent(inout)                             :: src
    type(mesh_2D), intent(inout)                            :: mesh
    type(meanflow2D), intent(in)                            :: flow2D
    type(output_parameters),intent(in)                      :: output
    integer                                                 :: nbout, nb_receiver
    ! output
    type(solution), intent(inout)                           :: sol


    type(solve)                                             :: M


    ! interpolated flow
    real,allocatable, dimension(:,:)                        :: u2, v2,x2,y2

    ! generel variables
    type(frequency)                                         :: freq
    real                                                    :: pi = 3.1415, ztest,a
    complex                                                 :: eye = (0,1), phi_interp
    complex                                                 :: kappa
    complex,allocatable,dimension(:)                        :: AlphamBeta, AlphapBeta
    integer                                                 :: ii, kk, index
    real                                                    :: start, finish
    integer, parameter                                      :: single = kind(1e0), double = kind(1d0)
    print*, 'start vector PE propagation...'
    nbout = output%nbout
    ! init frequency param
    freq%freq = f
    freq%omega = 2*pi*freq%freq
    freq%lambda = var0%c/freq%freq
    freq%k0 = 2*pi/freq%lambda

    imp%Zs = sqrt(4*imp%sigmae/(-eye*freq%omega*var0%gamma*var0%rho)) + var0%c*imp%alphae/(-eye*freq%omega*4*var0%gamma)
    imp%beta = 1/imp%Zs
    if (imp%rigid .eqv. .true.) imp%beta = 0
    !totally reflective ground hard code
    !imp%beta = 0
    ! PML definition
    !-------------------------------------------------------------------------
    allocate(pml%sigma(dom%nz))
    allocate(pml%sigma_de(dom%nz))

    pml%sigma = 0.
    pml%sigma_de = 0.

    do ii = dom%nz - pml%size, dom%nz
        pml%sigma(ii) = pml%param * (real(ii- dom%nz + pml%size)/real(pml%size))**pml%n

        pml%sigma_de(ii) = pml%param*pml%n*real((ii- dom%nz + pml%size))**(pml%n-1)/(dom%dz*(pml%size**pml%n))
    enddo

    ! init vector for matrix resolution
    !-------------------------------------------------------------------------
    allocate(m%h10(ceiling((1.*dom%nx)/nbout),dom%nz))
    allocate(m%h12(ceiling((1.*dom%nx)/nbout),dom%nz))
    allocate(m%h20(ceiling((1.*dom%nx)/nbout),dom%nz))
    allocate(m%h22(ceiling((1.*dom%nx)/nbout),dom%nz))

    kappa = eye*freq%k0*dom%dx/2

    m%h10 =  1. + 0.25*flow2D%epsilon
    m%h12 = 1/(4*(flow2D%gamma**2))
    m%h20 = ((flow2D%gamma**2)*flow2D%epsilon/2) - (1. + 0.25*flow2D%epsilon)*flow2D%tau
    m%h22 = 0.5 - flow2D%tau/(4*(flow2D%gamma**2))

   ! version wape normale
   !m%h10 =  1 + 0.25*flow2D%epsilon
   !m%h12 = 0.25
   !m%h20 = 0.5*flow2D%epsilon
   !m%h22 = 0.5

    allocate(m%a_vect(dom%nz))
    allocate(m%b_vect(dom%nz))
    allocate(m%c_vect(dom%nz))

    allocate(m%d_vect(dom%nz))
    allocate(m%e_vect(dom%nz))
    allocate(m%f_vect(dom%nz))
    allocate(AlphamBeta(dom%nz))
    allocate(AlphapBeta(dom%nz))

    AlphapBeta = ((1/(((freq%k0+eye*pml%sigma/var0%c)**2)*dom%dz**2)) - &
                & (eye*pml%sigma_de)/(2*(freq%omega+eye*pml%sigma)*((freq%k0+eye*pml%sigma/var0%c)**2)*dom%dz))
    AlphamBeta = ((1/(((freq%k0+eye*pml%sigma/var0%c)**2)*dom%dz**2)) + &
                & (eye*pml%sigma_de)/(2*(freq%omega+eye*pml%sigma)*((freq%k0+eye*pml%sigma/var0%c)**2)*dom%dz))

    m%a_vect = (m%h10(1,:) - kappa*m%h20(1,:)) - &
                2*(m%h12(1,:)- kappa*m%h22(1,:))/(((freq%k0+eye*pml%sigma/var0%c)**2)*dom%dz**2)
    m%b_vect = (m%h12(1,:)-kappa*m%h22(1,:)) * AlphamBeta
    m%c_vect = (m%h12(1,:)-kappa*m%h22(1,:)) * AlphapBeta


    m%d_vect = (m%h10(1,:) + kappa*m%h20(1,:)) - &
                2*(m%h12(1,:) + kappa*m%h22(1,:))/(((freq%k0+eye*pml%sigma/var0%c)**2)*dom%dz**2)
    m%e_vect = (m%h12(1,:) + kappa*m%h22(1,:)) * AlphamBeta
    m%f_vect = (m%h12(1,:) + kappa*m%h22(1,:)) * AlphapBeta

    m%c_vect(1) = m%c_vect(1)+m%b_vect(1)
    m%f_vect(1) = m%f_vect(1)+m%e_vect(1)

    m%a_vect(1) = m%a_vect(1) + 2*eye*freq%k0*imp%beta*dom%dz*m%b_vect(1);
    m%d_vect(1) = m%d_vect(1) + 2*eye*freq%k0*imp%beta*dom%dz*m%e_vect(1);

    ! propagation
    !----------------------------------------------------------------------------------------
    ! definition de phi(xn), phi(xn-1) and phi (xn+1) for the calculation of p (cf ostashev)

    allocate(sol%phi(dom%nz))
    allocate(sol%p(dom%nz))
    allocate(sol%phi_m1(dom%nz))
    allocate(sol%phi_p1(dom%nz))
    allocate(sol%u(dom%nz))
    allocate(sol%deltaL(dom%nz))

    ! Calcul du starter : starter de Salomons (Eqs. (G75) et (G76))
    sol%phi = sqrt(eye*freq%k0)*((1.3717-0.3701*(freq%k0*(mesh%z-src%pos_z))**2)* &
            & (exp(-((freq%k0*(mesh%z-src%pos_z))**2)/3)) + &
            & (1-imp%beta)/(1+imp%beta)*(1.3717-0.3701*(freq%k0*(mesh%z+src%pos_z))**2)* &
            & (exp(-((freq%k0*(mesh%z+src%pos_z))**2)/3)))

        sol%phi = sqrt(eye*freq%k0)*(1.3717-0.3701*(freq%k0*(mesh%z-src%pos_z))**2)* &
            & (exp(-((freq%k0*(mesh%z-src%pos_z))**2)/3))
    sol%phi_m1 = 0.
    sol%phi_p1 = sol%phi

    ! x advancement
    !---------------------------------------------------------------------------
    print*,"begin x loop ..."
    CALL CPU_TIME(start)
    do kk = 2,dom%nx
        index = (kk-1)/nbout+1
        if(mod(kk-1,nbout) == 0) then
            m%a_vect = (m%h10(index,:) - kappa*m%h20(index,:)) - &
                        & 2*(m%h12(index,:)-kappa*m%h22(index,:))/(((freq%k0+eye*pml%sigma/var0%c)**2)*dom%dz**2)
            m%b_vect = (m%h12(index,:)-kappa*m%h22(index,:)) * AlphamBeta
            m%c_vect = (m%h12(index,:)-kappa*m%h22(index,:)) * AlphapBeta


            m%d_vect = (m%h10(index,:) + kappa*m%h20(index,:)) - &
                        & 2*(m%h12(index,:)+kappa*m%h22(index,:))/(((freq%k0+eye*pml%sigma/var0%c)**2)*dom%dz**2)
            m%e_vect = (m%h12(index,:) + kappa*m%h22(index,:)) * AlphamBeta
            m%f_vect = (m%h12(index,:) + kappa*m%h22(index,:)) * AlphapBeta

            m%c_vect(1) = m%c_vect(1)+m%b_vect(1)
            m%f_vect(1) = m%f_vect(1)+m%e_vect(1)

            m%a_vect(1) = m%a_vect(1) + 2*eye*freq%k0*imp%beta*dom%dz*m%b_vect(1);
            m%d_vect(1) = m%d_vect(1) + 2*eye*freq%k0*imp%beta*dom%dz*m%e_vect(1);
        endif


        

        sol%u = sol%phi_p1;
        sol%phi_p1(1) = m%d_vect(1)*sol%u(1) + m%f_vect(1)*sol%u(2)
        sol%phi_p1(2:dom%nz-1) = m%e_vect(2:dom%nz-1)*sol%u(1:dom%nz-2) + m%d_vect(2:dom%nz-1)*sol%u(2:dom%nz-1) &
            & + m%f_vect(2:dom%nz-1)*sol%u(3:dom%nz)
        sol%phi_p1(dom%nz) = m%e_vect(dom%nz)*sol%u(dom%nz-1) + m%d_vect(dom%nz)*sol%u(dom%nz)
        call res_trid(sol%phi_p1, m%b_vect(2:dom%nz), m%a_vect, m%c_vect(1:dom%nz-1),dom%nz)
        sol%phi_p1(dom%nz) = 0;

        ! interpolation on general mesh grid
        if (kk>1) then
          do ii=1,dom%nz
            sol%p(ii) = (exp(eye*freq%k0*mesh%x(kk))/sqrt(mesh%x(kk))) * ((1- flow2D%M(index,ii))*sol%phi(ii) + &
              (eye*flow2D%M(index,ii)/(2*freq%k0*dom%dx))*(sol%phi_p1(ii)-sol%phi_m1(ii)))
            sol%deltaL(ii) = 20*log10(MAX(abs(sol%p(ii)),1e-10)*sqrt(mesh%x(kk)**2+(mesh%z(ii)-src%pos_z)**2))
          end do
          ! interpolation on coarse grid for side view
          if (output%side) then
            call interplog(mesh%z,sol%deltaL,mesh%z_coarse,sol%deltaL_vect(kk,:))
          endif
          ! call interpo_1D(mesh%z,sol%deltaL,mesh%z_coarse,sol%deltaL_vect(kk,:))
          !sol%deltaL_vect(kk,:) = spline3(mesh%z, sol%deltaL, mesh%z_coarse)

          ! interpolation for receiver
        !   print*, 'test0'
        !   print*, size(mesh%z),size(sol%deltaL),size(output%heights),size(sol%receiver(kk,:))
          call interplog(mesh%z,sol%deltaL,output%heights,sol%receiver(kk,:))
          !call interpo_1D(mesh%z,sol%deltaL,output%heights,sol%receiver(kk,:))
          !sol%receiver(kk,:) = spline3(mesh%z, sol%deltaL, output%heights)
        endif
        ! advancing of phi(xn) and phi(xn-1)
        sol%phi_m1 = sol%phi
        sol%phi = sol%phi_p1
    enddo
    CALL CPU_TIME(finish)
    print*,'done.'
    print*, 'x loop time : ', finish-start, 's'


    deallocate(pml%sigma)
    deallocate(pml%sigma_de)

    deallocate(AlphamBeta)
    deallocate(AlphapBeta)

    deallocate(sol%phi)
    deallocate(sol%phi_p1)
    deallocate(sol%phi_m1)
    deallocate(sol%u)
    deallocate(sol%p)
    deallocate(sol%deltaL)

    deallocate(m%a_vect)
    deallocate(m%b_vect)
    deallocate(m%c_vect)

    deallocate(m%d_vect)
    deallocate(m%e_vect)
    deallocate(m%f_vect)
end subroutine propa_1angle_1freq_arbitrary

subroutine propa_1angle_1freq(var0,f,src,imp,dom,pml,flow2D,mesh,sol,output)
    use types
    use utils
    use interpo
    use read_write
    use hdf5
    implicit none
    ! input parameters
    type(domain), intent(inout)                             :: dom
    type(real), intent(in)                                  :: f
    type(impedance), intent(inout)                          :: imp
    type(variables_0), intent(inout)                        :: var0
    type(perfectly_matched_layer),intent(inout)             :: pml
    type(source), intent(inout)                             :: src
    type(mesh_2D), intent(inout)                            :: mesh
    type(meanflow2D), intent(in)                            :: flow2D
    type(output_parameters),intent(in)                      :: output
    integer                                                 :: nbout, nb_receiver
    ! output
    type(solution), intent(inout)                           :: sol


    type(solve)                                             :: M


    ! interpolated flow
    real,allocatable, dimension(:,:)                        :: u2, v2,x2,y2

    ! generel variables
    type(frequency)                                         :: freq
    real                                                    :: pi = 3.1415, ztest
    complex                                                 :: eye = (0,1), phi_interp
    integer                                                 :: ii, kk
    real                                                    :: start, finish
    integer, parameter                                      :: single = kind(1e0), double = kind(1d0)

    print*, 'start scalar PE propagation...'
    nbout = output%nbout
    ! init frequency param
    freq%freq = f
    freq%omega = 2*pi*freq%freq
    freq%lambda = var0%c/freq%freq
    freq%k0 = 2*pi/freq%lambda

    imp%Zs = sqrt(4*imp%sigmae/(-eye*freq%omega*var0%gamma*var0%rho)) + var0%c*imp%alphae/(-eye*freq%omega*4*var0%gamma)
    imp%beta = 1/imp%Zs
    !totally reflective ground
    if (imp%rigid .eqv. .true.) imp%beta = 0
    ! PML definition
    !-------------------------------------------------------------------------
    allocate(pml%sigma(dom%nz))
    allocate(pml%sigma_de(dom%nz))

    pml%sigma = 0.
    pml%sigma_de = 0.

    do ii = dom%nz - pml%size, dom%nz
        pml%sigma(ii) = pml%param * (real(ii- dom%nz + pml%size)/real(pml%size))**pml%n

        pml%sigma_de(ii) = pml%param*pml%n*real((ii- dom%nz + pml%size))**(pml%n-1)/(dom%dz*(pml%size**pml%n))
    enddo

    ! init vector for matrix resolution
    !-------------------------------------------------------------------------
    allocate(m%nu_vect(dom%nz))
    allocate(m%mu_vect(dom%nz))

    allocate(m%a_vect(dom%nz))
    allocate(m%b_vect(dom%nz))
    allocate(m%c_vect(dom%nz))

    allocate(m%d_vect(dom%nz))
    allocate(m%e_vect(dom%nz))
    allocate(m%f_vect(dom%nz))

    m%nu = 0.25- eye*freq%k0*dom%dx/4
    m%mu = 0.25 + eye*freq%k0*dom%dx/4

    m%nu_vect = m%nu/(((freq%k0+eye*pml%sigma/var0%c)**2)*(dom%dz**2))
    m%mu_vect = m%mu/(((freq%k0+eye*pml%sigma/var0%c)**2)*(dom%dz**2))

    m%c_vect = m%nu_vect*(1-eye*pml%sigma_de*dom%dz/((freq%omega+eye*pml%sigma)*2))
    m%f_vect = m%mu_vect*(1-eye*pml%sigma_de*dom%dz/((freq%omega+eye*pml%sigma)*2));
    m%b_vect = m%nu_vect*(1+eye*pml%sigma_de*dom%dz/((freq%omega+eye*pml%sigma)*2));
    m%e_vect = m%mu_vect*(1+eye*pml%sigma_de*dom%dz/((freq%omega+eye*pml%sigma)*2));

    m%c_vect(1) = m%c_vect(1)+m%b_vect(1)
    m%f_vect(1) = m%f_vect(1)+m%e_vect(1)

    m%a_vect = 1 + m%nu*flow2D%epsilon(1,:) - 2*m%nu_vect;
    m%d_vect = 1 + m%mu*flow2D%epsilon(1,:) - 2*m%mu_vect;
    m%a_vect(1) = m%a_vect(1) + 2*eye*freq%k0*imp%beta*dom%dz*m%b_vect(1);
    m%d_vect(1) = m%d_vect(1) + 2*eye*freq%k0*imp%beta*dom%dz*m%e_vect(1);


    ! propagation
    !----------------------------------------------------------------------------------------
    allocate(sol%phi(dom%nz))
    allocate(sol%p(dom%nz))
    allocate(sol%deltaL(dom%nz))
    allocate(sol%u(dom%nz))

    ! Calcul du starter : starter de Salomons (Eqs. (G75) et (G76))
    sol%phi = sqrt(eye*freq%k0)*((1.3717-0.3701*(freq%k0*(mesh%z-src%pos_z))**2)* &
            & (exp(-((freq%k0*(mesh%z-src%pos_z))**2)/3)) + &
            & (1-imp%beta)/(1+imp%beta)*(1.3717-0.3701*(freq%k0*(mesh%z+src%pos_z))**2)* &
            & (exp(-((freq%k0*(mesh%z+src%pos_z))**2)/3)))
    
    sol%phi = sqrt(eye*freq%k0)*(1.3717-0.3701*(freq%k0*(mesh%z-src%pos_z))**2)* &
            & (exp(-((freq%k0*(mesh%z-src%pos_z))**2)/3))

    ! x advancement
    !----------------------------------------------------------------------------------------
    print*,"begin x loop ..."
    CALL CPU_TIME(start)
    do kk = 2,dom%nx
        if(mod(kk-1,nbout) == 0) then
            m%a_vect = 1 + m%nu*flow2D%epsilon((kk-1)/nbout+1,:) - 2*m%nu_vect
            m%d_vect = 1 + m%mu*flow2D%epsilon((kk-1)/nbout+1,:) - 2*m%mu_vect
            m%a_vect(1) = m%a_vect(1) + 2*eye*freq%k0*imp%beta*dom%dz*m%b_vect(1)
            m%d_vect(1) = m%d_vect(1) + 2*eye*freq%k0*imp%beta*dom%dz*m%e_vect(1)
        endif
        sol%u = sol%phi;
        sol%phi(1) = m%d_vect(1)*sol%u(1) + m%f_vect(1)*sol%u(2)
        sol%phi(2:dom%nz-1) = m%e_vect(2:dom%nz-1)*sol%u(1:dom%nz-2) + m%d_vect(2:dom%nz-1)*sol%u(2:dom%nz-1) &
            & + m%f_vect(2:dom%nz-1)*sol%u(3:dom%nz)
        sol%phi(dom%nz) = m%e_vect(dom%nz)*sol%u(dom%nz-1) + m%d_vect(dom%nz)*sol%u(dom%nz)
        call res_trid(sol%phi, m%b_vect(2:dom%nz), m%a_vect, m%c_vect(1:dom%nz-1),dom%nz)
        sol%phi(dom%nz) = 0

        ! interpolation on general mesh grid
        do ii=1,dom%nz
          sol%p(ii) = sol%phi(ii)/sqrt(mesh%x(kk))*exp(eye*freq%k0*mesh%x(kk))
          sol%deltaL(ii) = 20*log10(MAX(abs(sol%p(ii)),1e-10)*sqrt(mesh%x(kk)**2+(mesh%z(ii)-src%pos_z)**2))
        end do
        ! interpolation on coarse grid for side view
        call interplog(mesh%z,sol%deltaL,mesh%z_coarse,sol%deltaL_vect(kk,:))
        ! interpolation for receiver
        call interplog(mesh%z,sol%deltaL,output%heights,sol%receiver(kk,:))
    enddo
    CALL CPU_TIME(finish)
    print*,'done.'
    print*, 'x loop time : ', finish-start, 's'

    deallocate(pml%sigma)
    deallocate(pml%sigma_de)

    deallocate(m%nu_vect)
    deallocate(m%mu_vect)

    deallocate(m%a_vect)
    deallocate(m%b_vect)
    deallocate(m%c_vect)

    deallocate(m%d_vect)
    deallocate(m%e_vect)
    deallocate(m%f_vect)

    deallocate(sol%phi)
    deallocate(sol%p)
    deallocate(sol%deltaL)
    deallocate(sol%u)

end subroutine propa_1angle_1freq


subroutine propa_1angle_1freq_arbitrary_new(var0,f,src,imp,dom,pml,flow2D,mesh,sol,output)
    use types
    use utils
    use interpo
    use read_write
    use hdf5
    use splines
    implicit none
    ! input parameters
    type(domain), intent(inout)                             :: dom
    type(real), intent(in)                                  :: f
    type(impedance), intent(inout)                          :: imp
    type(variables_0), intent(inout)                        :: var0
    type(perfectly_matched_layer),intent(inout)             :: pml
    type(source), intent(inout)                             :: src
    type(mesh_2D), intent(inout)                            :: mesh
    type(meanflow2D), intent(in)                            :: flow2D
    type(output_parameters),intent(in)                      :: output
    integer                                                 :: nbout, nb_receiver
    ! output
    type(solution), intent(inout)                           :: sol


    type(solve)                                             :: M


    ! interpolated flow
    real,allocatable, dimension(:,:)                        :: u2, v2,x2,y2

    ! generel variables
    type(frequency)                                         :: freq
    real                                                    :: pi = 3.1415, ztest,a
    complex                                                 :: eye = (0,1), phi_interp
    complex                                                 :: kappa
    complex,allocatable,dimension(:)                        :: AlphamBeta, AlphapBeta
    integer                                                 :: ii, kk, index
    real                                                    :: start, finish
    integer, parameter                                      :: single = kind(1e0), double = kind(1d0)
    print*, 'start vector new PE propagation...'
    nbout = output%nbout
    ! init frequency param
    freq%freq = f
    freq%omega = 2*pi*freq%freq
    freq%lambda = var0%c/freq%freq
    freq%k0 = 2*pi/freq%lambda

    imp%Zs = sqrt(4*imp%sigmae/(-eye*freq%omega*var0%gamma*var0%rho)) + var0%c*imp%alphae/(-eye*freq%omega*4*var0%gamma)
    imp%beta = 1/imp%Zs
    if (imp%rigid .eqv. .true.) imp%beta = 0
    !totally reflective ground hard code
    !imp%beta = 0
    ! PML definition
    !-------------------------------------------------------------------------
    allocate(pml%sigma(dom%nz))
    allocate(pml%sigma_de(dom%nz))

    pml%sigma = 0.
    pml%sigma_de = 0.

    do ii = dom%nz - pml%size, dom%nz
        pml%sigma(ii) = pml%param * (real(ii- dom%nz + pml%size)/real(pml%size))**pml%n

        pml%sigma_de(ii) = pml%param*pml%n*real((ii- dom%nz + pml%size))**(pml%n-1)/(dom%dz*(pml%size**pml%n))
    enddo

    ! init vector for matrix resolution
    !-------------------------------------------------------------------------
    allocate(m%h10(ceiling((1.*dom%nx)/nbout),dom%nz))
    allocate(m%h12(ceiling((1.*dom%nx)/nbout),dom%nz))
    allocate(m%h20(ceiling((1.*dom%nx)/nbout),dom%nz))
    allocate(m%h22(ceiling((1.*dom%nx)/nbout),dom%nz))

    kappa = eye*freq%k0*dom%dx/2

    ! m%h10 =  1. + 0.25*flow2D%epsilon
    ! m%h12 = 1/(4*(flow2D%gamma**2))
    ! m%h20 = ((flow2D%gamma**2)*flow2D%epsilon/2) - (1. + 0.25*flow2D%epsilon)*flow2D%tau
    ! m%h22 = 0.5 - flow2D%tau/(4*(flow2D%gamma**2))

    ! new ostashev 
    m%h10 =  1.
    m%h12 = 0.
    m%h20 = -flow2D%m/(1+flow2D%m) ! which is C0/ceff -1 
    m%h22 = 0.5
   ! version wape normale
   !m%h10 =  1 + 0.25*flow2D%epsilon
   !m%h12 = 0.25
   !m%h20 = 0.5*flow2D%epsilon
   !m%h22 = 0.5

    allocate(m%a_vect(dom%nz))
    allocate(m%b_vect(dom%nz))
    allocate(m%c_vect(dom%nz))

    allocate(m%d_vect(dom%nz))
    allocate(m%e_vect(dom%nz))
    allocate(m%f_vect(dom%nz))
    allocate(AlphamBeta(dom%nz))
    allocate(AlphapBeta(dom%nz))

    AlphapBeta = ((1/(((freq%k0+eye*pml%sigma/var0%c)**2)*dom%dz**2)) - &
                & (eye*pml%sigma_de)/(2*(freq%omega+eye*pml%sigma)*((freq%k0+eye*pml%sigma/var0%c)**2)*dom%dz))
    AlphamBeta = ((1/(((freq%k0+eye*pml%sigma/var0%c)**2)*dom%dz**2)) + &
                & (eye*pml%sigma_de)/(2*(freq%omega+eye*pml%sigma)*((freq%k0+eye*pml%sigma/var0%c)**2)*dom%dz))

    m%a_vect = (m%h10(1,:) - kappa*m%h20(1,:)) - &
                2*(m%h12(1,:)- kappa*m%h22(1,:))/(((freq%k0+eye*pml%sigma/var0%c)**2)*dom%dz**2)
    m%b_vect = (m%h12(1,:)-kappa*m%h22(1,:)) * AlphamBeta
    m%c_vect = (m%h12(1,:)-kappa*m%h22(1,:)) * AlphapBeta


    m%d_vect = (m%h10(1,:) + kappa*m%h20(1,:)) - &
                2*(m%h12(1,:) + kappa*m%h22(1,:))/(((freq%k0+eye*pml%sigma/var0%c)**2)*dom%dz**2)
    m%e_vect = (m%h12(1,:) + kappa*m%h22(1,:)) * AlphamBeta
    m%f_vect = (m%h12(1,:) + kappa*m%h22(1,:)) * AlphapBeta

    m%c_vect(1) = m%c_vect(1)+m%b_vect(1)
    m%f_vect(1) = m%f_vect(1)+m%e_vect(1)

    m%a_vect(1) = m%a_vect(1) + 2*eye*freq%k0*imp%beta*dom%dz*m%b_vect(1);
    m%d_vect(1) = m%d_vect(1) + 2*eye*freq%k0*imp%beta*dom%dz*m%e_vect(1);

    ! propagation
    !----------------------------------------------------------------------------------------
    ! definition de phi(xn), phi(xn-1) and phi (xn+1) for the calculation of p (cf ostashev)

    allocate(sol%phi(dom%nz))
    allocate(sol%p(dom%nz))
    allocate(sol%phi_m1(dom%nz))
    allocate(sol%phi_p1(dom%nz))
    allocate(sol%u(dom%nz))
    allocate(sol%deltaL(dom%nz))

    ! Calcul du starter : starter de Salomons (Eqs. (G75) et (G76))
    sol%phi = sqrt(eye*freq%k0)*((1.3717-0.3701*(freq%k0*(mesh%z-src%pos_z))**2)* &
            & (exp(-((freq%k0*(mesh%z-src%pos_z))**2)/3)) + &
            & (1-imp%beta)/(1+imp%beta)*(1.3717-0.3701*(freq%k0*(mesh%z+src%pos_z))**2)* &
            & (exp(-((freq%k0*(mesh%z+src%pos_z))**2)/3)))
    sol%phi_m1 = 0.
    sol%phi_p1 = sol%phi

    ! x advancement
    !---------------------------------------------------------------------------
    print*,"begin x loop ..."
    CALL CPU_TIME(start)
    do kk = 2,dom%nx
        index = (kk-1)/nbout+1
        if(mod(kk-1,nbout) == 0) then
            m%a_vect = (m%h10(index,:) - kappa*m%h20(index,:)) - &
                        & 2*(m%h12(index,:)-kappa*m%h22(index,:))/(((freq%k0+eye*pml%sigma/var0%c)**2)*dom%dz**2)
            m%b_vect = (m%h12(index,:)-kappa*m%h22(index,:)) * AlphamBeta
            m%c_vect = (m%h12(index,:)-kappa*m%h22(index,:)) * AlphapBeta


            m%d_vect = (m%h10(index,:) + kappa*m%h20(index,:)) - &
                        & 2*(m%h12(index,:)+kappa*m%h22(index,:))/(((freq%k0+eye*pml%sigma/var0%c)**2)*dom%dz**2)
            m%e_vect = (m%h12(index,:) + kappa*m%h22(index,:)) * AlphamBeta
            m%f_vect = (m%h12(index,:) + kappa*m%h22(index,:)) * AlphapBeta

            m%c_vect(1) = m%c_vect(1)+m%b_vect(1)
            m%f_vect(1) = m%f_vect(1)+m%e_vect(1)

            m%a_vect(1) = m%a_vect(1) + 2*eye*freq%k0*imp%beta*dom%dz*m%b_vect(1);
            m%d_vect(1) = m%d_vect(1) + 2*eye*freq%k0*imp%beta*dom%dz*m%e_vect(1);
        endif

        sol%u = sol%phi_p1;
        sol%phi_p1(1) = m%d_vect(1)*sol%u(1) + m%f_vect(1)*sol%u(2)
        sol%phi_p1(2:dom%nz-1) = m%e_vect(2:dom%nz-1)*sol%u(1:dom%nz-2) + m%d_vect(2:dom%nz-1)*sol%u(2:dom%nz-1) &
            & + m%f_vect(2:dom%nz-1)*sol%u(3:dom%nz)
        sol%phi_p1(dom%nz) = m%e_vect(dom%nz)*sol%u(dom%nz-1) + m%d_vect(dom%nz)*sol%u(dom%nz)
        call res_trid(sol%phi_p1, m%b_vect(2:dom%nz), m%a_vect, m%c_vect(1:dom%nz-1),dom%nz)
        sol%phi_p1(dom%nz) = 0;

        ! interpolation on general mesh grid
        if (kk>1) then
          do ii=1,dom%nz
            sol%p(ii) = (exp(eye*freq%k0*mesh%x(kk))/sqrt(mesh%x(kk))) * ((1- flow2D%M(index,ii))*sol%phi(ii) + &
              (eye*flow2D%M(index,ii)/(2*freq%k0*dom%dx))*(sol%phi_p1(ii)-sol%phi_m1(ii)))
            sol%deltaL(ii) = 20*log10(MAX(abs(sol%p(ii)),1e-10)*sqrt(mesh%x(kk)**2+(mesh%z(ii)-src%pos_z)**2))
          end do
          ! interpolation on coarse grid for side view
          call interplog(mesh%z,sol%deltaL,mesh%z_coarse,sol%deltaL_vect(kk,:))
          ! call interpo_1D(mesh%z,sol%deltaL,mesh%z_coarse,sol%deltaL_vect(kk,:))
          !sol%deltaL_vect(kk,:) = spline3(mesh%z, sol%deltaL, mesh%z_coarse)

          ! interpolation for receiver
          call interplog(mesh%z,sol%deltaL,output%heights,sol%receiver(kk,:))
          !call interpo_1D(mesh%z,sol%deltaL,output%heights,sol%receiver(kk,:))
          !sol%receiver(kk,:) = spline3(mesh%z, sol%deltaL, output%heights)
        endif
        ! advancing of phi(xn) and phi(xn-1)
        sol%phi_m1 = sol%phi
        sol%phi = sol%phi_p1
    enddo
    CALL CPU_TIME(finish)
    print*,'done.'
    print*, 'x loop time : ', finish-start, 's'


    deallocate(pml%sigma)
    deallocate(pml%sigma_de)

    deallocate(AlphamBeta)
    deallocate(AlphapBeta)

    deallocate(sol%phi)
    deallocate(sol%phi_p1)
    deallocate(sol%phi_m1)
    deallocate(sol%u)
    deallocate(sol%p)
    deallocate(sol%deltaL)

    deallocate(m%a_vect)
    deallocate(m%b_vect)
    deallocate(m%c_vect)

    deallocate(m%d_vect)
    deallocate(m%e_vect)
    deallocate(m%f_vect)
end subroutine propa_1angle_1freq_arbitrary_new


subroutine propa_1angle_1freq_arbitrary2(var0,f,src,imp,dom,pml,flow2D,mesh,sol,output)
    use types
    use utils
    use interpo
    use read_write
    use hdf5
    use splines
    implicit none
    ! input parameters
    type(domain), intent(inout)                             :: dom
    type(real), intent(in)                                  :: f
    type(impedance), intent(inout)                          :: imp
    type(variables_0), intent(inout)                        :: var0
    type(perfectly_matched_layer),intent(inout)             :: pml
    type(source), intent(inout)                             :: src
    type(mesh_2D), intent(inout)                            :: mesh
    type(meanflow2D), intent(in)                            :: flow2D
    type(output_parameters),intent(in)                      :: output
    integer                                                 :: nbout, nb_receiver
    ! output
    type(solution), intent(inout)                           :: sol


    type(solve)                                             :: M


    ! interpolated flow
    real,allocatable, dimension(:,:)                        :: u2, v2,x2,y2

    ! generel variables
    type(frequency)                                         :: freq
    real                                                    :: pi = 3.1415, ztest,a
    complex                                                 :: eye = (0,1), phi_interp
    complex                                                 :: kappa
    complex,allocatable,dimension(:)                        :: AlphamBeta, AlphapBeta
    integer                                                 :: ii, kk, index
    real                                                    :: start, finish
    integer, parameter                                      :: single = kind(1e0), double = kind(1d0)

    print*, 'start new vector PE propagation...'
    nbout = output%nbout
    ! init frequency param
    freq%freq = f
    freq%omega = 2*pi*freq%freq
    freq%lambda = var0%c/freq%freq
    freq%k0 = 2*pi/freq%lambda

    imp%Zs = sqrt(4*imp%sigmae/(-eye*freq%omega*var0%gamma*var0%rho)) + var0%c*imp%alphae/(-eye*freq%omega*4*var0%gamma)
    imp%beta = 1/imp%Zs
    if (imp%rigid .eqv. .true.) imp%beta = 0
    !totally reflective ground hard code
    !imp%beta = 0
    ! PML definition
    !-------------------------------------------------------------------------
    allocate(pml%sigma(dom%nz))
    allocate(pml%sigma_de(dom%nz))

    pml%sigma = 0.
    pml%sigma_de = 0.

    do ii = dom%nz - pml%size, dom%nz
        pml%sigma(ii) = pml%param * (real(ii- dom%nz + pml%size)/real(pml%size))**pml%n

        pml%sigma_de(ii) = pml%param*pml%n*real((ii- dom%nz + pml%size))**(pml%n-1)/(dom%dz*(pml%size**pml%n))
    enddo

    ! init vector for matrix resolution
    !-------------------------------------------------------------------------
    allocate(m%h10(ceiling((1.*dom%nx)/nbout),dom%nz))
    allocate(m%h12(ceiling((1.*dom%nx)/nbout),dom%nz))
    allocate(m%h20(ceiling((1.*dom%nx)/nbout),dom%nz))
    allocate(m%h22(ceiling((1.*dom%nx)/nbout),dom%nz))

    kappa = eye*freq%k0*dom%dx/2

    m%h10 =  1. + 0.*flow2D%epsilon
    m%h12 = 0.25/((flow2D%gamma**2)*(1+flow2D%epsilon))
    m%h20 = -flow2D%tau
    m%h22 = 0.5/(1+flow2D%epsilon)**0.5 - 0.25*flow2D%tau/((flow2D%gamma**2)*(1+flow2D%epsilon))
   ! version wape normale
   !m%h10 =  1 + 0.25*flow2D%epsilon
   !m%h12 = 0.25
   !m%h20 = 0.5*flow2D%epsilon
   !m%h22 = 0.5

    allocate(m%a_vect(dom%nz))
    allocate(m%b_vect(dom%nz))
    allocate(m%c_vect(dom%nz))

    allocate(m%d_vect(dom%nz))
    allocate(m%e_vect(dom%nz))
    allocate(m%f_vect(dom%nz))
    allocate(AlphamBeta(dom%nz))
    allocate(AlphapBeta(dom%nz))

    AlphapBeta = ((1/(((freq%k0+eye*pml%sigma/var0%c)**2)*dom%dz**2)) - &
                & (eye*pml%sigma_de)/(2*(freq%omega+eye*pml%sigma)*((freq%k0+eye*pml%sigma/var0%c)**2)*dom%dz))
    AlphamBeta = ((1/(((freq%k0+eye*pml%sigma/var0%c)**2)*dom%dz**2)) + &
                & (eye*pml%sigma_de)/(2*(freq%omega+eye*pml%sigma)*((freq%k0+eye*pml%sigma/var0%c)**2)*dom%dz))

    m%a_vect = (m%h10(1,:) - kappa*m%h20(1,:)) - &
                2*(m%h12(1,:)- kappa*m%h22(1,:))/(((freq%k0+eye*pml%sigma/var0%c)**2)*dom%dz**2)
    m%b_vect = (m%h12(1,:)-kappa*m%h22(1,:)) * AlphamBeta
    m%c_vect = (m%h12(1,:)-kappa*m%h22(1,:)) * AlphapBeta


    m%d_vect = (m%h10(1,:) + kappa*m%h20(1,:)) - &
                2*(m%h12(1,:) + kappa*m%h22(1,:))/(((freq%k0+eye*pml%sigma/var0%c)**2)*dom%dz**2)
    m%e_vect = (m%h12(1,:) + kappa*m%h22(1,:)) * AlphamBeta
    m%f_vect = (m%h12(1,:) + kappa*m%h22(1,:)) * AlphapBeta

    m%c_vect(1) = m%c_vect(1)+m%b_vect(1)
    m%f_vect(1) = m%f_vect(1)+m%e_vect(1)

    m%a_vect(1) = m%a_vect(1) + 2*eye*freq%k0*imp%beta*dom%dz*m%b_vect(1);
    m%d_vect(1) = m%d_vect(1) + 2*eye*freq%k0*imp%beta*dom%dz*m%e_vect(1);

    ! propagation
    !----------------------------------------------------------------------------------------
    ! definition de phi(xn), phi(xn-1) and phi (xn+1) for the calculation of p (cf ostashev)

    allocate(sol%phi(dom%nz))
    allocate(sol%p(dom%nz))
    allocate(sol%phi_m1(dom%nz))
    allocate(sol%phi_p1(dom%nz))
    allocate(sol%u(dom%nz))
    allocate(sol%deltaL(dom%nz))

    ! Calcul du starter : starter de Salomons (Eqs. (G75) et (G76))
    sol%phi = sqrt(eye*freq%k0)*((1.3717-0.3701*(freq%k0*(mesh%z-src%pos_z))**2)* &
            & (exp(-((freq%k0*(mesh%z-src%pos_z))**2)/3)) + &
            & (1-imp%beta)/(1+imp%beta)*(1.3717-0.3701*(freq%k0*(mesh%z+src%pos_z))**2)* &
            & (exp(-((freq%k0*(mesh%z+src%pos_z))**2)/3)))
    sol%phi_m1 = 0.
    sol%phi_p1 = sol%phi

    ! x advancement
    !---------------------------------------------------------------------------
    print*,"begin x loop ..."
    CALL CPU_TIME(start)
    do kk = 2,dom%nx
        index = (kk-1)/nbout+1
        if(mod(kk-1,nbout) == 0) then
            m%a_vect = (m%h10(index,:) - kappa*m%h20(index,:)) - &
                        & 2*(m%h12(index,:)-kappa*m%h22(index,:))/(((freq%k0+eye*pml%sigma/var0%c)**2)*dom%dz**2)
            m%b_vect = (m%h12(index,:)-kappa*m%h22(index,:)) * AlphamBeta
            m%c_vect = (m%h12(index,:)-kappa*m%h22(index,:)) * AlphapBeta


            m%d_vect = (m%h10(index,:) + kappa*m%h20(index,:)) - &
                        & 2*(m%h12(index,:)+kappa*m%h22(index,:))/(((freq%k0+eye*pml%sigma/var0%c)**2)*dom%dz**2)
            m%e_vect = (m%h12(index,:) + kappa*m%h22(index,:)) * AlphamBeta
            m%f_vect = (m%h12(index,:) + kappa*m%h22(index,:)) * AlphapBeta

            m%c_vect(1) = m%c_vect(1)+m%b_vect(1)
            m%f_vect(1) = m%f_vect(1)+m%e_vect(1)

            m%a_vect(1) = m%a_vect(1) + 2*eye*freq%k0*imp%beta*dom%dz*m%b_vect(1);
            m%d_vect(1) = m%d_vect(1) + 2*eye*freq%k0*imp%beta*dom%dz*m%e_vect(1);
        endif

        sol%u = sol%phi_p1;
        sol%phi_p1(1) = m%d_vect(1)*sol%u(1) + m%f_vect(1)*sol%u(2)
        sol%phi_p1(2:dom%nz-1) = m%e_vect(2:dom%nz-1)*sol%u(1:dom%nz-2) + m%d_vect(2:dom%nz-1)*sol%u(2:dom%nz-1) &
            & + m%f_vect(2:dom%nz-1)*sol%u(3:dom%nz)
        sol%phi_p1(dom%nz) = m%e_vect(dom%nz)*sol%u(dom%nz-1) + m%d_vect(dom%nz)*sol%u(dom%nz)
        call res_trid(sol%phi_p1, m%b_vect(2:dom%nz), m%a_vect, m%c_vect(1:dom%nz-1),dom%nz)
        sol%phi_p1(dom%nz) = 0;

        ! interpolation on general mesh grid
        if (kk>1) then
          do ii=1,dom%nz
            sol%p(ii) = (exp(eye*freq%k0*mesh%x(kk))/sqrt(mesh%x(kk))) * ((1- flow2D%M(index,ii))*sol%phi(ii) + &
              (eye*flow2D%M(index,ii)/(2*freq%k0*dom%dx))*(sol%phi_p1(ii)-sol%phi_m1(ii)))
            sol%deltaL(ii) = 20*log10(MAX(abs(sol%p(ii)),1e-10)*sqrt(mesh%x(kk)**2+(mesh%z(ii)-src%pos_z)**2))
          end do
          ! interpolation on coarse grid for side view
          call interp1(mesh%z,sol%deltaL,mesh%z_coarse,sol%deltaL_vect(kk,:))
          ! interpolation for receiver
          call interp1(mesh%z,sol%deltaL,output%heights,sol%receiver(kk,:))
          ! sol%receiver(kk,:) = spline3(mesh%z, sol%deltaL, output%heights)
        endif
        ! advancing of phi(xn) and phi(xn-1)
        sol%phi_m1 = sol%phi
        sol%phi = sol%phi_p1
    enddo
    CALL CPU_TIME(finish)
    print*,'done.'
    print*, 'x loop time : ', finish-start, 's'


    deallocate(pml%sigma)
    deallocate(pml%sigma_de)

    deallocate(AlphamBeta)
    deallocate(AlphapBeta)

    deallocate(sol%phi)
    deallocate(sol%phi_p1)
    deallocate(sol%phi_m1)
    deallocate(sol%u)
    deallocate(sol%p)
    deallocate(sol%deltaL)

    deallocate(m%a_vect)
    deallocate(m%b_vect)
    deallocate(m%c_vect)

    deallocate(m%d_vect)
    deallocate(m%e_vect)
    deallocate(m%f_vect)
end subroutine propa_1angle_1freq_arbitrary2
end module wape
