#ifdef ADMR
!------------------------------------------------------------------------------!
! Description:
! ------------
!> Read in layout of the rotor blade , the lift and drag tables
!> and the distribution of lift and drag tables along the blade
!> Note: admr_input.dat is read in admr_init_arrays for now
!------------------------------------------------------------------------------!
!*******************************************************************************
subroutine admr_read_blade_tables
!*******************************************************************************
use param, only : rp,nproc,nx,ny,nz,dx,dy,dz,idx,idy,idz,dt,nz_tot,pi,z_i
use admr_mod
IMPLICIT NONE

INTEGER ::  ii,jj,i_ring,ierrn,iialpha,iir   !< running index
    
CHARACTER(200) :: chmess     !< Read in string 

INTEGER ::  dlen        !< no. rows of local table 
INTEGER ::  dlenbl      !< no. rows of cd, cl table
INTEGER ::  ialpha      !< table position of current alpha value
INTEGER ::  radres      !< radial resolution
INTEGER ::  t1          !< no. of airfoil
INTEGER ::  t2          !< no. of airfoil
INTEGER ::  trow        !< 
INTEGER ::  dlenbl_int  !< no. rows of interpolated cd, cl tables
    
REAL(rp) :: alpha_attack_i   !<
REAL(rp) :: weight_a         !< 
REAL(rp) :: weight_b         !< 
REAL(rp) :: cur_rr,rad_d

INTEGER, DIMENSION(:), ALLOCATABLE :: ttoint1    !<
INTEGER, DIMENSION(:), ALLOCATABLE :: ttoint2    !<
    
REAL(rp), DIMENSION(:), ALLOCATABLE :: turb_cd_sel1   !< 
REAL(rp), DIMENSION(:), ALLOCATABLE :: turb_cd_sel2   !<
REAL(rp), DIMENSION(:), ALLOCATABLE :: turb_cl_sel1   !< 
REAL(rp), DIMENSION(:), ALLOCATABLE :: turb_cl_sel2   !<
REAL(rp), DIMENSION(:), ALLOCATABLE :: read_cl_cd     !< read in var array
              
REAL(rp), DIMENSION(:), ALLOCATABLE    :: alpha_attack_tab   !< 
REAL(rp), DIMENSION(:), ALLOCATABLE    :: trad1              !< 
REAL(rp), DIMENSION(:), ALLOCATABLE    :: trad2              !<          
REAL(rp), DIMENSION(:,:), ALLOCATABLE  :: turb_cd_table      !<
REAL(rp), DIMENSION(:,:), ALLOCATABLE  :: turb_cl_table      !<
                                          
allocate( read_cl_cd(1:2*nairfoils+1) )

!
!--    Read in the distribution of lift and drag tables along the blade, the
!--    layout of the rotor blade and the lift and drag tables:

open( 201, FILE='NREL_DATA.dat', STATUS='OLD', FORM='FORMATTED', IOSTAT=ierrn )

if ( ierrn /= 0 )  then
 if(coord == 0) write(*,*) 'file NREL_DATA does not exist'
endif
!
!--    Read distribution table:
dlen = 0
        READ ( 201, '(3/)' )

       rloop3: DO
          READ ( 201, *, IOSTAT=ierrn ) chmess
          IF ( ierrn < 0  .OR.  chmess == '#'  .OR.  chmess == '')  EXIT rloop3
          dlen = dlen + 1
       ENDDO rloop3

       ALLOCATE( trad1(1:dlen), trad2(1:dlen), ttoint1(1:dlen), ttoint2(1:dlen))

       DO jj = 1,dlen+1
          BACKSPACE ( 201, IOSTAT=ierrn )
       ENDDO

       DO jj = 1,dlen
          READ ( 201, * ) trad1(jj), trad2(jj), ttoint1(jj), ttoint2(jj)
       ENDDO

!--    Read layout table:
       dlen = 0 

       READ ( 201, '(3/)')

       rloop1: DO
          READ ( 201, *, IOSTAT=ierrn ) chmess
          IF ( ierrn < 0  .OR.  chmess == '#'  .OR.  chmess == '')  EXIT rloop1
          dlen = dlen + 1
       ENDDO rloop1

       ALLOCATE( lrd(1:dlen), ard(1:dlen), crd(1:dlen) )
       DO jj = 1, dlen+1
          BACKSPACE ( 201, IOSTAT=ierrn )
       ENDDO             
       DO jj = 1, dlen
          READ ( 201, * ) lrd(jj), ard(jj), crd(jj) 
       ENDDO
       ! lrd = Pos[m] ard =twist[deg]  crd =chord[m]

!
!--    Read tables (turb_cl(alpha),turb_cd(alpha) for the different profiles:

       dlen = 0

       ! Skip first 3 lines of file
       READ ( 201, '(3/)' )

       rloop2: DO
          READ ( 201, *, IOSTAT=ierrn ) chmess
          IF ( ierrn < 0  .OR.  chmess == '#'  .OR.  chmess == '')  EXIT rloop2
          dlen = dlen + 1
       ENDDO rloop2 


       ALLOCATE( alpha_attack_tab(1:dlen), turb_cl_table(1:dlen,1:nairfoils),  &
                 turb_cd_table(1:dlen,1:nairfoils) )

       DO jj = 1,dlen+1
          BACKSPACE ( 201, IOSTAT=ierrn )
       ENDDO 

       DO jj = 1,dlen
          READ ( 201, * ) read_cl_cd
          alpha_attack_tab(jj) = read_cl_cd(1)
          DO ii= 1, nairfoils
             turb_cl_table(jj,ii) = read_cl_cd(ii*2)
             turb_cd_table(jj,ii) = read_cl_cd(ii*2+1)
          ENDDO
       ENDDO

       dlenbl = dlen 

       CLOSE ( 201 )

!
!--    For each possible radial position (resolution: 0.1 m --> 630 values) and
!--    each possible angle of attack (resolution: 0.01 degrees --> 36000 values!)
!--    determine the lift and drag coefficient by interpolating between the
!--    tabulated values of each table (interpolate to current angle of attack)
!--    and between the tables (interpolate to current radial position):

       ALLOCATE( turb_cl_sel1(0:dlenbl) )  
       ALLOCATE( turb_cl_sel2(0:dlenbl) )  
       ALLOCATE( turb_cd_sel1(0:dlenbl) )
       ALLOCATE( turb_cd_sel2(0:dlenbl) )

       radres     = INT( admr(1)%rr * 10.0_rp ) + 1
       dlenbl_int = INT( 360.0_rp / accu_cl_cd_tab ) + 1


       ALLOCATE( turb_cl_tab(0:dlenbl_int,1:radres) )
       ALLOCATE( turb_cd_tab(0:dlenbl_int,1:radres) )

       DO iir = 1, radres ! loop over radius
             cur_rr = ( iir - 1 ) * 0.1_rp
!
!--          Find position in table
             lct = MINLOC( ABS( trad1 - cur_rr ) )
!                lct(1) = lct(1)

             IF ( ( trad1(lct(1)) - cur_rr ) .GT. 0.0 )  THEN
                lct(1) = lct(1) - 1
             ENDIF

             trow = lct(1)


!
!--          Calculate weights for interpolation
             weight_a = ( trad2(trow) - cur_rr ) / ( trad2(trow) - trad1(trow) )
             weight_b = ( cur_rr - trad1(trow) ) / ( trad2(trow) - trad1(trow) )
             t1 = ttoint1(trow)
             t2 = ttoint2(trow)

             IF ( t1 .EQ. t2 ) THEN  ! if both are the same, the weights are NaN
                weight_a = 0.5_rp    ! then do interpolate in between same twice
                weight_b = 0.5_rp    ! using 0.5 as weight
             ENDIF

             IF ( t1 == 0 .AND. t2 == 0 ) THEN
                turb_cd_sel1 = 0.0_rp
                turb_cd_sel2 = 0.0_rp
                turb_cl_sel1 = 0.0_rp
                turb_cl_sel2 = 0.0_rp
                turb_cd_tab(1,iir) = 0.0_rp  ! For -180Â° (iialpha=1) the values  
                turb_cl_tab(1,iir) = 0.0_rp  ! for each radius has to be set explicitly  
             ELSE
                turb_cd_sel1 = turb_cd_table(:,t1)
                turb_cd_sel2 = turb_cd_table(:,t2)
                turb_cl_sel1 = turb_cl_table(:,t1)
                turb_cl_sel2 = turb_cl_table(:,t2)
                turb_cd_tab(1,iir) = ( weight_a * turb_cd_table(1,t1) + weight_b  &
                                    * turb_cd_table(1,t2) ) 
                turb_cl_tab(1,iir) = ( weight_a * turb_cl_table(1,t1) + weight_b  & 
                                    * turb_cl_table(1,t2) ) 
             ENDIF


          DO iialpha = 2, dlenbl_int  ! loop over angles

             alpha_attack_i = -180.0_rp + REAL( iialpha-1 ) * accu_cl_cd_tab
             ialpha = 1 ! AS: Can lead to segmentations fault below

             !DO WHILE ( alpha_attack_i > alpha_attack_tab(ialpha) )
             DO WHILE ( ( alpha_attack_i > alpha_attack_tab(ialpha) ) .AND. (ialpha <= dlen ) )
                ialpha = ialpha + 1
             ENDDO

!
!--          Interpolation of lift and drag coefficiencts on fine grid of radius 
!--          segments and angles of attack

             turb_cl_tab(iialpha,iir) = ( alpha_attack_tab(ialpha) -           &
                                          alpha_attack_i ) /                   &
                                        ( alpha_attack_tab(ialpha) -           &
                                          alpha_attack_tab(ialpha-1) ) *       &
                                        ( weight_a * turb_cl_sel1(ialpha-1) +  &
                                          weight_b * turb_cl_sel2(ialpha-1) ) +&
                                        ( alpha_attack_i             -         &
                                          alpha_attack_tab(ialpha-1) ) /       &
                                        ( alpha_attack_tab(ialpha) -           &
                                          alpha_attack_tab(ialpha-1) ) *       &
                                        ( weight_a * turb_cl_sel1(ialpha) +    &
                                          weight_b * turb_cl_sel2(ialpha) )
             turb_cd_tab(iialpha,iir) = ( alpha_attack_tab(ialpha) -           &
                                          alpha_attack_i ) /                   &
                                        ( alpha_attack_tab(ialpha) -           &
                                          alpha_attack_tab(ialpha-1) ) *       &
                                        ( weight_a * turb_cd_sel1(ialpha-1) +  &
                                          weight_b * turb_cd_sel2(ialpha-1) ) +&
                                        ( alpha_attack_i             -         &
                                          alpha_attack_tab(ialpha-1) ) /       &
                                        ( alpha_attack_tab(ialpha) -           &
                                          alpha_attack_tab(ialpha-1) ) *       &
                                        ( weight_a * turb_cd_sel1(ialpha) +    &
                                          weight_b * turb_cd_sel2(ialpha) )
   
          ENDDO   ! end loop over angles of attack

       ENDDO   ! end loop over radius


!
!--                Interpolation of the local pitch angle from tabulated values
!--                to the current radial position:
do inot = 1,nturb
 do i_ring = 1,admr(inot)%nrings
  lct=minloc(ABS(admr(inot)%cur_r(i_ring)-lrd))
  rad_d=admr(inot)%cur_r(i_ring)-lrd(lct(1))

  if (admr(inot)%cur_r(i_ring) == 0.0_rp) then
   admr(inot)%alpha_attack_in(i_ring) = 0.0_rp
  elseif( lct(1) >= size(ard)) then
   admr(inot)%alpha_attack_in(i_ring) = ( ard(size(ard)) + ard(size(ard)-1) ) / 2.0_rp
  else
   admr(inot)%alpha_attack_in(i_ring) = ( ard(lct(1)) *                                  &
                                        ( ( lrd(lct(1)+1) - admr(inot)%cur_r(i_ring) ) / &
                                        ( lrd(lct(1)+1) - lrd(lct(1)) )                  &
                                        ) ) + ( ard(lct(1)+1) *                          &
                                        ( ( admr(inot)%cur_r(i_ring) - lrd(lct(1)) ) /   &
                                        ( lrd(lct(1)+1) - lrd(lct(1)) ) ) )
  end if

!--                In Fortran radian instead of degree is used as unit for all
!--                angles. Therefore, a transformation from angles given in
!--                degree to angles given in radian is necessary here:
    admr(inot)%alpha_attack_in(i_ring) = admr(inot)%alpha_attack_in(i_ring) * ( (2.0_rp*pi) / 360.0_rp )


!--                Interpolation of the chord_length from tabulated values to
!--                the current radial position:

  if (admr(inot)%cur_r(i_ring) == 0.0_rp) then
   admr(inot)%chord(i_ring) = 0.0_rp
  elseif (lct(1) >= size(ard)) then
   admr(inot)%chord(i_ring) = (crd(size(crd)) + ard(size(crd)-1)) / 2.0_rp
  else
   admr(inot)%chord(i_ring) = ( crd(lct(1)) *                                  &
                              ( ( lrd(lct(1)+1) - admr(inot)%cur_r(i_ring) ) / &
                              ( lrd(lct(1)+1) - lrd(lct(1)) ) ) ) +            &
                              ( crd(lct(1)+1) *                                &
                              ( ( admr(inot)%cur_r(i_ring)-lrd(lct(1)) ) /     &
                              ( lrd(lct(1)+1) - lrd(lct(1)) ) ) )
  end if
 enddo
enddo  
end subroutine admr_read_blade_tables
#endif ADMR
