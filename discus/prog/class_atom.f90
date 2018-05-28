MODULE atom_class
!
!  Defines a type "cl_atom".
!  This basic atom has a scattering type, a position (x,y,z) and a property
!
!  Procedures are bound to this type that do:
!     set_atom        set type, position and property
!     get_atom        get type, position and property
!  ==>more are needed to set/get iscat only etc
!
IMPLICIT none
!
PUBLIC cl_atom
!
TYPE :: cl_atom                        ! the basic type definition
   PRIVATE
   INTEGER                :: iscat     ! the scattering type
   REAL   , DIMENSION(3)  :: pos       ! fractional coordinates
   INTEGER                :: prop      ! a property flag
   INTEGER, DIMENSION(0:3):: surf      ! surface charcteristics
   INTEGER, DIMENSION(1:2):: in_mole   ! In which molecule are we
!
CONTAINS
   PROCEDURE, PUBLIC, PASS :: set_atom ! set verything
   PROCEDURE, PUBLIC, PASS :: get_atom ! get everything
END TYPE cl_atom                       ! end of basic type definition
!
CONTAINS
!
!  Methods to work on atoms
!
   SUBROUTINE set_atom(this, itype, posit, iprop, isurface , iin_mole)
!
!  Set atom type, position and property
!
!  As the procedure is defined with "PASS", the first argument is 
!  always the object itself and thus the call to this subroutine
!  must be:
!  CALL object_name%set_atom ( itype, posit, iprop ) 
!  The object is placed before the "%" symbol
!
   IMPLICIT none
!
   CLASS (cl_atom)                  :: this  ! work on "this" atoms
   INTEGER,              INTENT(IN) :: itype ! atom type
   REAL   , DIMENSION(3),INTENT(IN) :: posit ! atom position
   INTEGER,              INTENT(IN) :: iprop ! atom property
   INTEGER, DIMENSION(0:3),INTENT(IN) :: isurface ! atom surface
   INTEGER, DIMENSION(1:2),INTENT(IN) :: iin_mole   ! In which molecule are we
!
   this%iscat   = itype
   this%pos(:)  = posit(:)
   this%prop    = iprop
   this%surf(:) = isurface(:)
   this%in_mole(:) = iin_mole(:)
!
   END SUBROUTINE set_atom
!
   SUBROUTINE get_atom(this, itype, posit, iprop, isurface, iin_mole)
!
!  Get atom type, position and property
!
   IMPLICIT none
!
   CLASS (cl_atom)                   :: this
   INTEGER,              INTENT(OUT) :: itype
   REAL   , DIMENSION(3),INTENT(OUT) :: posit
   INTEGER,              INTENT(OUT) :: iprop
   INTEGER, DIMENSION(0:3),INTENT(OUT) :: isurface ! atom surface
   INTEGER, DIMENSION(1:2),INTENT(OUT) :: iin_mole   ! In which molecule are we
!
   itype       = this%iscat
   posit(:)    = this%pos(:)
   iprop       = this%prop
   isurface(:) = this%surf(:)
   iin_mole(:) = this%in_mole(:)
!
   END SUBROUTINE get_atom
!
END MODULE atom_class
