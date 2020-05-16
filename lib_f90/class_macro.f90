MODULE macro_class
!
USE precision_mod
!
TYPE cl_macro
   INTEGER                  :: macro_length
   CHARACTER (LEN=PREC_STRING),  DIMENSION(:), ALLOCATABLE :: macro_line
   LOGICAL                  :: lmacro
!  TYPE (cl_macro), POINTER :: before
!  TYPE (cl_macro), POINTER :: after
   CONTAINS
!
   PROCEDURE, PUBLIC, PASS :: alloc_arrays       ! Allocate the lines for a macro
   PROCEDURE, PUBLIC, PASS :: set_lines          ! Write the length  of a macro
   PROCEDURE, PUBLIC, PASS :: set_macro          ! Write the content of a macro
   PROCEDURE, PUBLIC, PASS :: get_lines          ! Get all macro lines
   PROCEDURE, PUBLIC, PASS :: get_single         ! Get a single line of a macro
   PROCEDURE, PUBLIC, PASS :: get_macro          ! Get   the content of a macro
   PROCEDURE, PUBLIC, PASS :: finalize_macro     ! Finalize the cl_macro class
END TYPE cl_macro
!
CONTAINS
!
!******************************************************************************
!
SUBROUTINE alloc_arrays   ( this, nlines)
!
!  Allocate the arrays for "this" macro 
!  Initialize all variables
!
IMPLICIT none
!
CLASS (cl_macro)                 :: this         ! Work on "this" macro
INTEGER,              INTENT(IN) :: nlines       ! Number of lines for this macro
!
INTEGER  :: istatus
!
IF(this%lmacro) THEN                             ! macro was allocated 
   DEALLOCATE ( this%macro_line , STAT=istatus ) ! Deallocate macro lines
ENDIF
ALLOCATE ( this%macro_line(1:nlines), STAT=istatus ) ! Allocate macro lines
!
this%macro_length = nlines                       ! Store number of macro lines
this%lmacro = .true.                             ! Flag that macro is allocated
!
END SUBROUTINE alloc_arrays
!
!******************************************************************************
!
SUBROUTINE set_macro (this, nlines, content )
!
!  Set the number of lines in "this" macro
!
IMPLICIT none
!
CLASS (cl_macro) :: this
INTEGER,                                INTENT(IN) :: nlines   ! Number of lines for this macro
CHARACTER (LEN=*), DIMENSION(1:nlines), INTENT(IN) :: content  ! Actual lines for this macro
!
INTEGER :: i
!
DO i=1, nlines
   this%macro_line(i) = content(i)
ENDDO
!
END SUBROUTINE set_macro
!
!******************************************************************************
!
SUBROUTINE set_lines (this, nlines )
!
!  Set the number of lines in "this" macro
!
IMPLICIT none
!
CLASS (cl_macro) :: this
INTEGER,              INTENT(IN) :: nlines       ! Number of lines for this macro
!
this%macro_length = nlines
!
END SUBROUTINE set_lines
!
!******************************************************************************
!
INTEGER FUNCTION get_lines (this )
!
!  Return the number of lines in "this" macro
!
IMPLICIT none
!
CLASS (cl_macro) :: this
!
get_lines = this%macro_length
!
END FUNCTION get_lines
!
!******************************************************************************
!
SUBROUTINE get_single (this, number, string )
!
!  Return the number of lines in "this" macro
!
IMPLICIT none
!
CLASS (cl_macro) :: this
INTEGER,          INTENT(IN)  :: number   ! Number of lines for this macro
CHARACTER (LEN=*), INTENT(OUT) :: string   ! Actual lines for this macro
!
string = this%macro_line(number)
!
END SUBROUTINE get_single
!
!******************************************************************************
!
SUBROUTINE get_macro (this, nlines, content )
!
!  Return the macro content
!
IMPLICIT none
!
CLASS (cl_macro) :: this
INTEGER,                               INTENT(IN)  :: nlines   ! Number of lines for this macro
CHARACTER (LEN=*), DIMENSION(1:nlines), INTENT(OUT) :: content  ! Actual lines for this macro
!
INTEGER :: i
!
DO i=1, nlines
   content(i) = this%macro_line(i)
ENDDO
!
END SUBROUTINE get_macro
!
!******************************************************************************
!
SUBROUTINE finalize_macro   ( this)
!
!  Finalize the class cl_macro
!  Deallocate the arrays for "this" macro 
!
IMPLICIT none
!
CLASS (cl_macro)                 :: this         ! Work on "this" crystal
!
INTEGER  :: istatus
!
IF(ALLOCATED(this%macro_line)) THEN
   DEALLOCATE ( this%macro_line , STAT=istatus )    ! Deallocate macro lines
ENDIF
!
this%macro_length = 0                            ! Store number of macro lines
this%lmacro       = .false.                      ! Flag that macro is allocated
!
END SUBROUTINE finalize_macro
!
!******************************************************************************
!
END MODULE macro_class
