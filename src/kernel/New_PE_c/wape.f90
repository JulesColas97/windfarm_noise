module wape
implicit none
public :: propa_1angle_1freq
public :: propa_1angle_1freq_arbitrary
public :: propa_1angle_1freq_arbitrary2
contains

subroutine propa_1angle_1freq_arbitrary(var0,f,theta,src,imp,dom,pml,flow2D,mesh,sol,output)
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
    type(real), intent(in)                                  :: theta
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
    complex, allocatable,dimension(:)                       :: temp_rec
    ! generel variables
    type(frequency)                                         :: freq
    REAL                                                    :: pi = 3.141592653589793, theta_rad
    real                                                    :: ztest,a
    complex                                                 :: eye = (0,1), phi_interp
    complex                                                 :: kappa
    complex,allocatable,dimension(:)                        :: AlphamBeta, AlphapBeta
    complex, allocatable, dimension(:)                      :: a_temp,b_temp,c_temp
    integer                                                 :: ii, kk, index
    real                                                    :: start, finish, t_tot
    real                                                    :: start_propa, finish_propa, t_propa
    real                                                    :: start_dl, finish_dl, t_dl
    real                                                    :: start_interp, finish_interp, t_interp
    real                                                    :: start_interpr, finish_interpr, t_interpr
    real                                                    :: start_interps, finish_interps, t_interps
    real                                                    :: start_interpp, finish_interpp, t_interpp
    integer, parameter                                      :: single = kind(1e0), double = kind(1d0)
    print*, 'start vector PE propagation...'
    nbout = output%nbout
    allocate(temp_rec(output%nb_receiver))

    theta_rad = theta * pi/180
    if (theta_rad > pi) then 
        theta_rad = theta_rad - 2*pi
    end if 

    ! init frequency param
    freq%freq = f
    freq%omega = 2*pi*freq%freq
    freq%lambda = var0%c/freq%freq
    freq%k0 = 2*pi/freq%lambda

    imp%Zs = sqrt(4*imp%sigmae/(-eye*freq%omega*var0%gamma*var0%rho)) + var0%c*imp%alphae/(-eye*freq%omega*4*var0%gamma)
    imp%beta = 1/imp%Zs
    if (imp%rigid .eqv. .true.) imp%beta = 0

    !totally reflective ground hard code
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
    allocate(m%h10(dom%nz))
    allocate(m%h12(dom%nz))
    allocate(m%h20(dom%nz))
    allocate(m%h22(dom%nz))

    kappa = eye*freq%k0*dom%dx/2

    m%h10 =  1. + 0.25*flow2D%epsilon(1,:)
    m%h12 = 1/(4*(flow2D%gamma(1,:)**2))
    m%h20 = ((flow2D%gamma(1,:)**2)*flow2D%epsilon(1,:)/2) - (1. + 0.25*flow2D%epsilon(1,:))*flow2D%tau(1,:)
    m%h22 = 0.5 - flow2D%tau(1,:)/(4*(flow2D%gamma(1,:)**2))

    allocate(m%a_vect(dom%nz))
    allocate(m%b_vect(dom%nz))
    allocate(m%c_vect(dom%nz))

    allocate(a_temp(dom%nz))
    allocate(b_temp(dom%nz-1))
    allocate(c_temp(dom%nz-1))
    

    allocate(m%d_vect(dom%nz))
    allocate(m%e_vect(dom%nz))
    allocate(m%f_vect(dom%nz))
    allocate(AlphamBeta(dom%nz))
    allocate(AlphapBeta(dom%nz))

    AlphapBeta = ((1/(((freq%k0+eye*pml%sigma/var0%c)**2)*dom%dz**2)) - &
                & (eye*pml%sigma_de)/(2*(freq%omega+eye*pml%sigma)*((freq%k0+eye*pml%sigma/var0%c)**2)*dom%dz))
    AlphamBeta = ((1/(((freq%k0+eye*pml%sigma/var0%c)**2)*dom%dz**2)) + &
                & (eye*pml%sigma_de)/(2*(freq%omega+eye*pml%sigma)*((freq%k0+eye*pml%sigma/var0%c)**2)*dom%dz))

    m%a_vect = (m%h10 - kappa*m%h20) - &
                2*(m%h12- kappa*m%h22)/(((freq%k0+eye*pml%sigma/var0%c)**2)*dom%dz**2)
    m%b_vect = (m%h12-kappa*m%h22) * AlphamBeta
    m%c_vect = (m%h12-kappa*m%h22) * AlphapBeta


    m%d_vect = (m%h10 + kappa*m%h20) - &
                2*(m%h12 + kappa*m%h22)/(((freq%k0+eye*pml%sigma/var0%c)**2)*dom%dz**2)
    m%e_vect = (m%h12 + kappa*m%h22) * AlphamBeta
    m%f_vect = (m%h12 + kappa*m%h22) * AlphapBeta

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

    sol%phi = sqrt(eye*freq%k0)*(1.3717-0.3701*(freq%k0*(mesh%z-src%pos_z))**2)* &
            & (exp(-((freq%k0*(mesh%z-src%pos_z))**2)/3))
    sol%phi_m1 = 0.
    sol%phi_p1 = sol%phi

    ! x advancement
    !---------------------------------------------------------------------------
    print*,"begin x loop ..."
    t_propa = 0 
    t_dl = 0 
    t_interps = 0
    t_interp = 0
    t_interpr = 0
    t_interpp = 0
    CALL CPU_TIME(start)
    do kk = 2,dom%nx
        index = (kk-1)/nbout+1
        if(mod(kk-1,nbout) == 0) then
            m%h10 =  1. + 0.25*flow2D%epsilon(index,:)
            m%h12 = 1/(4*(flow2D%gamma(index,:)**2))
            m%h20 = ((flow2D%gamma(index,:)**2)*flow2D%epsilon(index,:)/2) - (1. + 0.25*flow2D%epsilon(index,:))*flow2D%tau(index,:)
            m%h22 = 0.5 - flow2D%tau(index,:)/(4*(flow2D%gamma(index,:)**2))

            m%a_vect = (m%h10 - kappa*m%h20) - &
                        & 2*(m%h12-kappa*m%h22)/(((freq%k0+eye*pml%sigma/var0%c)**2)*dom%dz**2)
            m%b_vect = (m%h12-kappa*m%h22) * AlphamBeta
            m%c_vect = (m%h12-kappa*m%h22) * AlphapBeta


            m%d_vect = (m%h10 + kappa*m%h20) - &
                        & 2*(m%h12+kappa*m%h22)/(((freq%k0+eye*pml%sigma/var0%c)**2)*dom%dz**2)
            m%e_vect = (m%h12 + kappa*m%h22) * AlphamBeta
            m%f_vect = (m%h12 + kappa*m%h22) * AlphapBeta

            m%c_vect(1) = m%c_vect(1)+m%b_vect(1)
            m%f_vect(1) = m%f_vect(1)+m%e_vect(1)

            m%a_vect(1) = m%a_vect(1) + 2*eye*freq%k0*imp%beta*dom%dz*m%b_vect(1);
            m%d_vect(1) = m%d_vect(1) + 2*eye*freq%k0*imp%beta*dom%dz*m%e_vect(1);
        endif

        CALL CPU_TIME(start_propa)
        sol%u = sol%phi_p1;
        sol%phi_p1(1) = m%d_vect(1)*sol%u(1) + m%f_vect(1)*sol%u(2)
        sol%phi_p1(2:dom%nz-1) = m%e_vect(2:dom%nz-1)*sol%u(1:dom%nz-2) + m%d_vect(2:dom%nz-1)*sol%u(2:dom%nz-1) &
            & + m%f_vect(2:dom%nz-1)*sol%u(3:dom%nz)
        sol%phi_p1(dom%nz) = m%e_vect(dom%nz)*sol%u(dom%nz-1) + m%d_vect(dom%nz)*sol%u(dom%nz)
        call res_trid(sol%phi_p1, m%b_vect(2:dom%nz), m%a_vect, m%c_vect(1:dom%nz-1),dom%nz)
        
        !sol%u = sol%phi_p1;
        !a_temp = m%a_vect
        !b_temp = m%b_vect(2:dom%nz)
        !c_temp = m%c_vect(1:dom%nz-1)
        !call res_trid_2(sol%phi_p1,sol%u, b_temp, a_temp, c_temp,dom%nz)
        
        
        sol%phi_p1(dom%nz) = 0;
        CALL CPU_TIME(finish_propa)

        t_propa = t_propa + finish_propa - start_propa

        ! interpolation on general mesh grid
        if (kk>1) then
            CALL CPU_TIME(start_dl)
            !sol%p(:) = (exp(eye*freq%k0*mesh%x(kk))/sqrt(mesh%x(kk))) * ((1- flow2D%M(index,:))*sol%phi(:) + &
            !    (eye*flow2D%M(index,:)/(2*freq%k0*dom%dx))*(sol%phi_p1(:)-sol%phi_m1(:)))
            !sol%deltaL(:) = 20*log10(MAX(abs(sol%p(:)),1e-10)*sqrt(mesh%x(kk)**2+(mesh%z(:)-src%pos_z)**2))
            
            ! original version for computinf p then dl 
            !do ii=1,dom%nz
            !sol%p(ii) = (exp(eye*freq%k0*mesh%x(kk))/sqrt(mesh%x(kk))) * ((1- flow2D%M(index,ii))*sol%phi(ii) + &
            !(eye*flow2D%M(index,ii)/(2*freq%k0*dom%dx))*(sol%phi_p1(ii)-sol%phi_m1(ii)))
            !sol%deltaL(ii) = 20*log10(MAX(abs(sol%p(ii)),1e-10)*sqrt(mesh%x(kk)**2+(mesh%z(ii)-src%pos_z)**2))
            !end do

            ! optimized verion
            sol%p(:) = sqrt(mesh%x(kk) +  (mesh%z(:)-src%pos_z)**2/mesh%x(kk)) * ((1- flow2D%M(index,:))*sol%phi(:) + &
                (eye*flow2D%M(index,:)/(2*freq%k0*dom%dx))*(sol%phi_p1(:)-sol%phi_m1(:)))

            CALL CPU_TIME(finish_dl)

            CALL CPU_TIME(start_interp)
            ! interpolation on coarse grid for side view
            CALL CPU_TIME(start_interps)
            if (output%side) then
                !call interplog(mesh%z,sol%deltaL,mesh%z_coarse,sol%deltaL_vect(kk,:))
                call interpcomplextoreal(mesh%z,sol%p,mesh%z_coarse,sol%deltaL_vect(kk,:))
                sol%deltaL_vect(kk,:) = 20*log10(MAX(sol%deltaL_vect(kk,:),1e-10))
                CALL CPU_TIME(finish_interps)
            elseif (kk<dom%nx)  then
                CALL CPU_TIME(finish_interps)

                ! interpolation at given point to construct final plane
                !--------------------------------------------------------------
                CALL CPU_TIME(start_interpp)
                do ii=1,size(output%xplane)
                if (output%xplane(ii)>=0) then
                    if ((mesh%x(kk)*cos(theta_rad)<=output%xplane(ii)) &
                        .and. (mesh%x(kk+1)*cos(theta_rad)>=output%xplane(ii))) then 
                        !call interplog(mesh%z,sol%deltaL,mesh%z_coarse,sol%xplane(ii,:))
                        
                        call interpcomplextoreal(mesh%z,sol%p,mesh%z_coarse,sol%xplane(ii,:))
                        sol%xplane(ii,:) = 20*log10(MAX(sol%xplane(ii,:),1e-10))

                        sol%ycoord(ii) = mesh%x(kk)*sin(theta_rad)
                        sol%xcount(ii) = 1 
                    endif 
                else 
                    if ((mesh%x(kk)*cos(theta_rad)>=output%xplane(ii)) &
                        .and. (mesh%x(kk+1)*cos(theta_rad)<=output%xplane(ii))) then 
                        !call interplog(mesh%z,sol%deltaL,mesh%z_coarse,sol%xplane(ii,:))
                        call interpcomplextoreal(mesh%z,sol%p,mesh%z_coarse,sol%xplane(ii,:))
                        sol%xplane(ii,:) = 20*log10(MAX(sol%xplane(ii,:),1e-10))
                        sol%ycoord(ii) = mesh%x(kk)*sin(theta_rad)
                        sol%xcount(ii) = 1 
                    endif
                endif  
                enddo 

                do ii=1,size(output%yplane)
                if (output%yplane(ii)>0) then 
                    if ((mesh%x(kk)*sin(theta_rad)<=output%yplane(ii)) &
                        .and. (mesh%x(kk+1)*sin(theta_rad)>=output%yplane(ii))) then 
                        !call interplog(mesh%z,sol%deltaL,mesh%z_coarse,sol%yplane(ii,:))
                        call interpcomplextoreal(mesh%z,sol%p,mesh%z_coarse,sol%yplane(ii,:))
                        sol%yplane(ii,:) = 20*log10(MAX(sol%yplane(ii,:),1e-10))
                        sol%xcoord(ii) = mesh%x(kk)*cos(theta_rad)
                        sol%ycount(ii) =  1
                    endif 
                else 
                    if ((mesh%x(kk)*sin(theta_rad)>=output%yplane(ii)) &
                        .and. (mesh%x(kk+1)*sin(theta_rad)<=output%yplane(ii))) then 
                        !call interplog(mesh%z,sol%deltaL,mesh%z_coarse,sol%yplane(ii,:))
                        call interpcomplextoreal(mesh%z,sol%p,mesh%z_coarse,sol%yplane(ii,:))
                        sol%yplane(ii,:) = 20*log10(MAX(sol%yplane(ii,:),1e-10))
                        sol%xcoord(ii) = mesh%x(kk)*cos(theta_rad)
                        sol%ycount(ii) =  1
                    endif 
                endif 
                enddo 
                CALL CPU_TIME(finish_interpp)
            endif 

            !------------------------------------------------------------------
            CALL CPU_TIME(start_interpr)
            ! interpolation for receiver
            ! call interpcomplex(mesh%z,sol%deltaL,output%heights,sol%receiver(kk,:))
            ! call interpcomplex(mesh%z,sol%p,output%heights,temp_rec) 
            !call interplog(mesh%z,sol%deltaL,output%heights,sol%receiver(kk,:))
            ! record deltaL in the receiver 
            ! og version 
            !sol%receiver(kk,:) = 20*log10(MAX(sol%receiver(kk,:),1e-10)*sqrt(mesh%x(kk)**2+(output%heights - src%pos_z)**2))
            ! new version 
            call interpcomplextoreal(mesh%z,sol%p,output%heights,sol%receiver(kk,:)) 
            sol%receiver(kk,:) = 20*log10(MAX(sol%receiver(kk,:),1e-10))

            CALL CPU_TIME(finish_interpr)
        endif
    CALL CPU_TIME(finish_interp)

    t_interp = t_interp + finish_interp- start_interp
    t_dl = t_dl + finish_dl - start_dl
    t_interps = t_interps + finish_interps- start_interps
    t_interpr = t_interpr + finish_interpr- start_interpr
    t_interpp = t_interpp + finish_interpp- start_interpp
    ! advancing of phi(xn) and phi(xn-1)
    sol%phi_m1 = sol%phi
    sol%phi = sol%phi_p1
    enddo
    CALL CPU_TIME(finish)
    print*,'done.'
    print*, 'x loop time : ', finish-start, 's'
    print*, 'resolution: ', t_propa/(finish-start)*100, '%'
    print*, 'compute dl  : ', t_dl/(finish-start)*100, '%'
    print*, 'interpolation : ', t_interp/(finish-start)*100, '%'
    print*, 'interpolation side : ', t_interps/(finish-start)*100, '%'
    print*, 'interpolation rec : ', t_interpr/(finish-start)*100, '%'
    print*, 'interpolation plane : ', t_interpp/(finish-start)*100, '%'


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
    allocate(m%h10(dom%nz))
    allocate(m%h12(dom%nz))
    allocate(m%h20(dom%nz))
    allocate(m%h22(dom%nz))

    kappa = eye*freq%k0*dom%dx/2

    ! m%h10 =  1. + 0.25*flow2D%epsilon
    ! m%h12 = 1/(4*(flow2D%gamma**2))
    ! m%h20 = ((flow2D%gamma**2)*flow2D%epsilon/2) - (1. + 0.25*flow2D%epsilon)*flow2D%tau
    ! m%h22 = 0.5 - flow2D%tau/(4*(flow2D%gamma**2))

    ! new ostashev 
    m%h10 =  1.
    m%h12 = 0.
    m%h20 = -flow2D%m(1,:)/(1+flow2D%m(1,:)) ! which is C0/ceff -1 
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

    m%a_vect = (m%h10 - kappa*m%h20) - &
                2*(m%h12- kappa*m%h22)/(((freq%k0+eye*pml%sigma/var0%c)**2)*dom%dz**2)
    m%b_vect = (m%h12-kappa*m%h22) * AlphamBeta
    m%c_vect = (m%h12-kappa*m%h22) * AlphapBeta


    m%d_vect = (m%h10 + kappa*m%h20) - &
                2*(m%h12 + kappa*m%h22)/(((freq%k0+eye*pml%sigma/var0%c)**2)*dom%dz**2)
    m%e_vect = (m%h12 + kappa*m%h22) * AlphamBeta
    m%f_vect = (m%h12 + kappa*m%h22) * AlphapBeta

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
            m%h20 = -flow2D%m(index,:)/(1+flow2D%m(index,:))

            m%a_vect = (m%h10 - kappa*m%h20) - &
                        & 2*(m%h12-kappa*m%h22)/(((freq%k0+eye*pml%sigma/var0%c)**2)*dom%dz**2)
            m%b_vect = (m%h12-kappa*m%h22) * AlphamBeta
            m%c_vect = (m%h12-kappa*m%h22) * AlphapBeta


            m%d_vect = (m%h10 + kappa*m%h20) - &
                        & 2*(m%h12+kappa*m%h22)/(((freq%k0+eye*pml%sigma/var0%c)**2)*dom%dz**2)
            m%e_vect = (m%h12 + kappa*m%h22) * AlphamBeta
            m%f_vect = (m%h12 + kappa*m%h22) * AlphapBeta

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
    allocate(m%h10(dom%nz))
    allocate(m%h12(dom%nz))
    allocate(m%h20(dom%nz))
    allocate(m%h22(dom%nz))

    kappa = eye*freq%k0*dom%dx/2

    m%h10 =  1. + 0.*flow2D%epsilon(1,:)
    m%h12 = 0.25/((flow2D%gamma(1,:)**2)*(1+flow2D%epsilon(1,:)))
    m%h20 = -flow2D%tau(1,:)
    m%h22 = 0.5/(1+flow2D%epsilon(1,:))**0.5 - 0.25*flow2D%tau(1,:)/((flow2D%gamma(1,:)**2)*(1+flow2D%epsilon(1,:)))
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

    m%a_vect = (m%h10 - kappa*m%h20) - &
                2*(m%h12- kappa*m%h22)/(((freq%k0+eye*pml%sigma/var0%c)**2)*dom%dz**2)
    m%b_vect = (m%h12-kappa*m%h22) * AlphamBeta
    m%c_vect = (m%h12-kappa*m%h22) * AlphapBeta


    m%d_vect = (m%h10 + kappa*m%h20) - &
                2*(m%h12 + kappa*m%h22)/(((freq%k0+eye*pml%sigma/var0%c)**2)*dom%dz**2)
    m%e_vect = (m%h12 + kappa*m%h22) * AlphamBeta
    m%f_vect = (m%h12 + kappa*m%h22) * AlphapBeta

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

            m%h10 =  1. + 0.*flow2D%epsilon(index,:)
            m%h12 = 0.25/((flow2D%gamma(index,:)**2)*(1+flow2D%epsilon(index,:)))
            m%h20 = -flow2D%tau(index,:)
            m%h22 = 0.5/(1+flow2D%epsilon(index,:))**0.5 - &
                0.25*flow2D%tau(index,:)/((flow2D%gamma(index,:)**2)*(1+flow2D%epsilon(index,:)))
            
            m%a_vect = (m%h10 - kappa*m%h20) - &
                        & 2*(m%h12-kappa*m%h22)/(((freq%k0+eye*pml%sigma/var0%c)**2)*dom%dz**2)
            m%b_vect = (m%h12-kappa*m%h22) * AlphamBeta
            m%c_vect = (m%h12-kappa*m%h22) * AlphapBeta


            m%d_vect = (m%h10 + kappa*m%h20) - &
                        & 2*(m%h12+kappa*m%h22)/(((freq%k0+eye*pml%sigma/var0%c)**2)*dom%dz**2)
            m%e_vect = (m%h12 + kappa*m%h22) * AlphamBeta
            m%f_vect = (m%h12 + kappa*m%h22) * AlphapBeta

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
