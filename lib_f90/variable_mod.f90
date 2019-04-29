MODULE variable_mod
!
   USE constants_mod
!
   IMPLICIT NONE
   PUBLIC
   SAVE
!
   TYPE var_arrays
      INTEGER, DIMENSION(2  )              :: var_shape
      REAL   , DIMENSION(:,:), ALLOCATABLE :: var_value
   END TYPE var_arrays
!
   INTEGER, PARAMETER :: VAR_MAX       = 999
!
!  INTEGER, PARAMETER :: VAR_TYPE_INTE = 0
!  INTEGER, PARAMETER :: VAR_TYPE_REAL = 1
!  INTEGER, PARAMETER :: VAR_TYPE_CHAR = 2
!
   INTEGER, PARAMETER :: VAR_UNDEF     =  1
   INTEGER, PARAMETER :: VAR_TRUE      =  2
   INTEGER, PARAMETER :: VAR_FALSE     =  3
   INTEGER, PARAMETER :: VAR_IS_TOP    =  4
   INTEGER, PARAMETER :: VAR_IS_SECTION=  5
   INTEGER, PARAMETER :: VAR_IS_BRANCH =  6
   INTEGER, PARAMETER :: VAR_SUITE     =  7
   INTEGER, PARAMETER :: VAR_DISCUS    =  8
   INTEGER, PARAMETER :: VAR_DIFFEV    =  9
   INTEGER, PARAMETER :: VAR_KUPLOT    = 10
   INTEGER, PARAMETER :: VAR_REFINE    = 11
   INTEGER, PARAMETER :: VAR_MPI_OFF   = 12
   INTEGER, PARAMETER :: VAR_MPI_ON    = 13
   INTEGER, PARAMETER :: VAR_PROGRAM   = 14
   INTEGER, PARAMETER :: VAR_STATE     = 15
   INTEGER, PARAMETER :: VAR_MPI_FIRST = 16
   INTEGER, PARAMETER :: VAR_MPI       = 17
   INTEGER, PARAMETER :: VAR_NUM_NODES = 18
   INTEGER, PARAMETER :: VAR_FIT_VALUE = 19
   INTEGER, PARAMETER :: VAR_SLOW      = 20
   INTEGER, PARAMETER :: VAR_LOOP      = 21
!
   INTEGER                                   :: var_num   = 0   ! Number of variables
   INTEGER                                   :: var_con   = 0   ! Number of internal constants
   INTEGER                                   :: var_sys   = 0   ! Number of system variables
   INTEGER                                   :: var_ref   = 0   ! Number of first refinement variable
   INTEGER                                   :: var_n_arr = 0   ! Number of variable arrays
   CHARACTER (LEN=200), DIMENSION(VAR_MAX)   :: var_name  = ' ' ! Variable name
   CHARACTER (LEN=200), DIMENSION(VAR_MAX)   :: var_char  = ' ' ! Character string if 'char' variable
   INTEGER            , DIMENSION(VAR_MAX)   :: var_l     = 0   ! Length of variable name
   INTEGER            , DIMENSION(VAR_MAX)   :: var_type  = 0   ! INT=0, REAL=1, CHAR=2
   INTEGER            , DIMENSION(VAR_MAX)   :: var_entry = 0   ! Entry in 2D variable list var_field
   LOGICAL            , DIMENSION(VAR_MAX)   :: var_diff  = .FALSE.   ! This is a DIFFEV variable yes/no
   REAL               , DIMENSION(VAR_MAX)   :: var_val   = 0.0 ! Value for INT or REAL
!
   TYPE(var_arrays), DIMENSION(VAR_MAX) :: var_field
!
CONTAINS
!
   SUBROUTINE variable_init
!
   INTEGER :: i
   i = 0
!
   i = VAR_UNDEF
   var_name ( i) = 'UNDEF'             ! A state is undefined
   var_char ( i) = ' '
   var_l    ( i) = 5
   var_type ( i) = IS_INTE
   var_diff ( i) = .FALSE.
   var_val  ( i) = -1
!
   i = VAR_TRUE
   var_name ( i) = 'TRUE'              ! A state is TRUE 
   var_char ( i) = ' '
   var_l    ( i) = 4
   var_type ( i) = IS_INTE
   var_diff ( i) = .FALSE.
   var_val  ( i) =  1
!
   i = VAR_FALSE
   var_name ( i) = 'FALSE'             ! A state is FALSE
   var_char ( i) = ' '
   var_l    ( i) = 5
   var_type ( i) = IS_INTE
   var_diff ( i) = .FALSE.
   var_val  ( i) =  0
!
   i = VAR_IS_TOP
   var_name ( i) = 'IS_TOP'            ! A section is in STATE 'TOP', 'Section', 'BRANCH'
   var_char ( i) = ' '
   var_l    ( i) = 6
   var_type ( i) = IS_INTE
   var_diff ( i) = .FALSE.
   var_val  ( i) = 0
!
   i = VAR_IS_SECTION
   var_name ( i) = 'IS_SECTION'        ! A section is in STATE 'TOP', 'Section', 'BRANCH'
   var_char ( i) = ' '
   var_l    ( i) = 10
   var_type ( i) = IS_INTE
   var_diff ( i) = .FALSE.
   var_val  ( i) = 1
!
   i = VAR_IS_BRANCH
   var_name ( i) = 'IS_BRANCH'         ! A section is in STATE 'TOP', 'Section', 'BRANCH'
   var_char ( i) = ' '
   var_l    ( i) = 9
   var_type ( i) = IS_INTE
   var_diff ( i) = .FALSE.
   var_val  ( i) = 2
!
   i = VAR_SUITE
   var_name ( i) = 'SUITE'          ! A PROGRAM is 'SUITE', 'DISCUS', 'KUPLOT', 'DIFFEV'
   var_char ( i) = ' '
   var_l    ( i) = 5
   var_type ( i) = IS_INTE
   var_diff ( i) = .FALSE.
   var_val  ( i) = 0
!
   i = VAR_DISCUS
   var_name ( i) = 'DISCUS'         ! A PROGRAM is 'SUITE', 'DISCUS', 'KUPLOT', 'DIFFEV'
   var_char ( i) = ' '
   var_l    ( i) = 6
   var_type ( i) = IS_INTE
   var_diff ( i) = .FALSE.
   var_val  ( i) = 1
!
   i = VAR_DIFFEV
   var_name ( i) = 'DIFFEV'         ! A PROGRAM is 'SUITE', 'DISCUS', 'KUPLOT', 'DIFFEV'
   var_char ( i) = ' '
   var_l    ( i) = 6
   var_type ( i) = IS_INTE
   var_diff ( i) = .FALSE.
   var_val  ( i) = 2
!
   i = VAR_KUPLOT
   var_name ( i) = 'KUPLOT'         ! A PROGRAM is 'SUITE', 'DISCUS', 'KUPLOT', 'DIFFEV'
   var_char ( i) = ' '
   var_l    ( i) = 6
   var_type ( i) = IS_INTE
   var_diff ( i) = .FALSE.
   var_val  ( i) = 3
!
   i = VAR_REFINE
   var_name ( i) = 'REFINE'         ! A PROGRAM is 'SUITE', 'DISCUS', 'KUPLOT', 'DIFFEV'
   var_char ( i) = ' '
   var_l    ( i) = 6
   var_type ( i) = IS_INTE
   var_diff ( i) = .FALSE.
   var_val  ( i) = 3
!
   i = VAR_MPI_OFF
   var_name ( i) = 'MPI_OFF'         ! A PROGRAM is 'SUITE', 'DISCUS', 'KUPLOT', 'DIFFEV'
   var_char ( i) = ' '
   var_l    ( i) = 7
   var_type ( i) = IS_INTE
   var_diff ( i) = .FALSE.
   var_val  ( i) = 0
!
   i = VAR_MPI_ON
   var_name ( i) = 'MPI_ON'         ! A PROGRAM is 'SUITE', 'DISCUS', 'KUPLOT', 'DIFFEV'
   var_char ( i) = ' '
   var_l    ( i) = 6
   var_type ( i) = IS_INTE
   var_diff ( i) = .FALSE.
   var_val  ( i) = 1
!
   i = VAR_PROGRAM
   var_name ( i) = 'PROGRAM'        ! A PROGRAM is 'SUITE', 'DISCUS', 'KUPLOT', 'DIFFEV'
   var_char ( i) = ' '
   var_l    ( i) = 7
   var_type ( i) = IS_INTE
   var_diff ( i) = .FALSE.
   var_val  ( i) = 0                ! Default to 'SUITE'
!
   i = VAR_STATE
   var_name ( i) = 'STATE'          ! A section is in STATE 'TOP', 'Section', 'BRANCH'
   var_char ( i) = ' '
   var_l    ( i) = 5
   var_type ( i) = IS_INTE
   var_diff ( i) = .FALSE.
   var_val  ( i) = var_val(VAR_IS_TOP)  ! Default to 'TOP'
!
   i = VAR_MPI_FIRST
   var_name ( i) = 'MPI_FIRST'        ! MPI_FIRST can be UNDEFINED, TRUE, FALSE
   var_char ( i) = ' '
   var_l    ( i) = 9
   var_type ( i) = IS_INTE
   var_diff ( i) = .FALSE.
   var_val  ( i) = var_val(VAR_UNDEF)
!
   i = VAR_MPI
   var_name ( i) = 'MPI'            ! MPI can be ON/OFF
   var_char ( i) = ' '
   var_l    ( i) = 3
   var_type ( i) = IS_INTE
   var_diff ( i) = .FALSE.
   var_val  ( i) = var_val(VAR_MPI_OFF)
!
   i = VAR_NUM_NODES
   var_name ( i) = 'NUM_NODES'      ! MPI can be ON/OFF
   var_char ( i) = ' '
   var_l    ( i) = 9
   var_type ( i) = IS_INTE
   var_diff ( i) = .FALSE.
   var_val  ( i) = 0
!
   i = VAR_FIT_VALUE
   var_name ( i) = 'F_VALUE'        ! Function value of a fit
   var_char ( i) = ' '
   var_l    ( i) = 7
   var_type ( i) = IS_REAL
   var_diff ( i) = .FALSE.
   var_val  ( i) = 0.0
!
   i = VAR_SLOW
   var_name ( i) = 'SLOW'          ! Generic loop variable 0
   var_char ( i) = ' '
   var_l    ( i) = 4
   var_type ( i) = IS_INTE
   var_diff ( i) = .FALSE.
   var_val  ( i) = 0
!
   i = VAR_LOOP
   var_name ( i) = 'LOOP'           ! Generic loop variable 1
   var_char ( i) = ' '
   var_l    ( i) = 4
   var_type ( i) = IS_INTE
   var_diff ( i) = .FALSE.
   var_val  ( i) = 0
!
   var_con      = i                ! Store number of system "constants" set internally
!
   i = i + 1
   var_ref      = i                ! Store first refinement variable entry number
   var_name ( i) = 'REF_GENERATION'
   var_char ( i) = ' '
   var_l    ( i) = 14
   var_type ( i) = IS_INTE
   var_diff ( i) = .FALSE.
   var_val  ( i) = 0
!
   i = i + 1
   var_name ( i) = 'REF_MEMBER'
   var_char ( i) = ' '
   var_l    ( i) = 10
   var_type ( i) = IS_INTE
   var_diff ( i) = .FALSE.
   var_val  ( i) = 0
!
   i = i + 1
   var_name ( i) = 'REF_CHILDREN'
   var_char ( i) = ' '
   var_l    ( i) = 12
   var_type ( i) = IS_INTE
   var_diff ( i) = .FALSE.
   var_val  ( i) = 0
!
   i = i + 1
   var_name ( i) = 'REF_DIMENSION'
   var_char ( i) = ' '
   var_l    ( i) = 13
   var_type ( i) = IS_INTE
   var_diff ( i) = .FALSE.
   var_val  ( i) = 0
!
   i = i + 1
   var_name ( i) = 'REF_KID'
   var_char ( i) = ' '
   var_l    ( i) =  7
   var_type ( i) = IS_INTE
   var_diff ( i) = .FALSE.
   var_val  ( i) = 0
!
   i = i + 1
   var_name ( i) = 'REF_INDIV'
   var_char ( i) = ' '
   var_l    ( i) =  9
   var_type ( i) = IS_INTE
   var_diff ( i) = .FALSE.
   var_val  ( i) = 0
!
   i = i + 1
   var_name ( i) = 'REF_NINDIV'
   var_char ( i) = ' '
   var_l    ( i) = 10
   var_type ( i) = IS_INTE
   var_diff ( i) = .FALSE.
   var_val  ( i) = 0
!
   i = i + 1
   var_name ( i) = 'INDI_PARALLEL'
   var_char ( i) = ' '
   var_l    ( i) = 13
   var_type ( i) = IS_INTE
   var_diff ( i) = .FALSE.
   var_val  ( i) = var_val(VAR_FALSE)
!
   i = i + 1
   var_name ( i) = 'PI'
   var_char ( i) = ' '
   var_l    ( i) = 2
   var_type ( i) = IS_REAL
   var_diff ( i) = .FALSE.
   var_val  ( i) = REAL(3.1415926535897932384626433832795028841971693993751D0)
!
!
   var_num      = i
   var_sys      = i
!
   var_entry(1:i) = 0          ! all system variables are scalar, thus no entry
!
   END SUBROUTINE variable_init
!
   LOGICAL FUNCTION is_variable(string)
!
   USE param_mod
!
   CHARACTER(LEN=*), INTENT(IN) :: string
!
   INTEGER :: i,ihyp, ihyp2
   INTEGER :: i1, i2
!
   i1 = 1
   i2 = LEN_TRIM(string)
   ihyp = MAX (INDEX (string, '''') , INDEX (string, '"') )
   IF(ihyp > 0) THEN
      i1 = ihyp+1
      ihyp2 = ihyp + MAX (INDEX (string(ihyp+1:i2), '''') ,  &
                          INDEX (string(ihyp+1:i2), '"')   )
      i2 = ihyp2-1
   ENDIF
   is_variable = .FALSE.
   IF(i2 >= i1) THEN
      loop: DO i=1, var_num
         IF(string(i1:i2)==var_name(i)) THEN
            is_variable = .TRUE.
            res_para(0) = 2
            res_para(1) = 1
            res_para(2) = var_type(i)
            EXIT LOOP
         ENDIF
      ENDDO loop
   ENDIF
   END FUNCTION is_variable
!
END MODULE variable_mod
