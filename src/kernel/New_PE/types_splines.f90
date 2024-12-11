module types_splines
implicit none
private
public dp, hp, ivector, dvector, zvector

!integer, parameter :: dp=kind(0.d0), &          ! double precision
integer, parameter :: dp=selected_real_kind(8), &       ! MODIF 10/04/2020
                      hp=selected_real_kind(15) ! high precision

type ivector                       ! allocatable integer vector
   integer, pointer :: vec(:) => null()
end type

type dvector                       ! allocatable real double precision vector
   real(dp), pointer :: vec(:) => null()
end type

type zvector                       ! allocatable complex double precision vector
   complex(dp), pointer :: vec(:) => null()
end type

end module





