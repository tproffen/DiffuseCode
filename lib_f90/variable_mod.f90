MODULE variable_mod
!
   IMPLICIT NONE
   PUBLIC
   SAVE
!
   INTEGER, PARAMETER :: VAR_MAX       = 150
!
   INTEGER, PARAMETER :: VAR_TYPE_INTE = 0
   INTEGER, PARAMETER :: VAR_TYPE_REAL = 1
   INTEGER, PARAMETER :: VAR_TYPE_CHAR = 2
!
   INTEGER                                 :: var_num  = 0   ! Number of variables
   INTEGER                                 :: var_sys  = 0   ! Number of system variables
   INTEGER                                 :: var_ref  = 0   ! Number of first refinement variable
   CHARACTER (LEN=200), DIMENSION(VAR_MAX) :: var_name = ' ' ! (VAR_MAX)
   CHARACTER (LEN=200), DIMENSION(VAR_MAX) :: var_char = ' ' ! (VAR_MAX)
   INTEGER            , DIMENSION(VAR_MAX) :: var_l    = 0   ! (VAR_MAX)
   INTEGER            , DIMENSION(VAR_MAX) :: var_type = 0   ! (VAR_MAX)
   REAL               , DIMENSION(VAR_MAX) :: var_val  = 0.0 ! (VAR_MAX)
!
CONTAINS
!
   SUBROUTINE variable_init
!
   INTEGER :: i
   i = 0
!
   i = i + 1
   var_ref      = i                ! Store first refinement variable entry number
   var_name( i) = 'REF_GENERATION'
   var_char( i) = ' '
   var_l   ( i) = 14
   var_type( i) = 0
   var_val ( i) = 0
!
   i = i + 1
   var_name( i) = 'REF_MEMBER'
   var_char( i) = ' '
   var_l   ( i) = 10
   var_type( i) = 0
   var_val ( i) = 0
!
   i = i + 1
   var_name( i) = 'REF_CHILDREN'
   var_char( i) = ' '
   var_l   ( i) = 12
   var_type( i) = 0
   var_val ( i) = 0
!
   i = i + 1
   var_name( i) = 'REF_DIMENSION'
   var_char( i) = ' '
   var_l   ( i) = 13
   var_type( i) = 0
   var_val ( i) = 0
!
   i = i + 1
   var_name( i) = 'REF_KID'
   var_char( i) = ' '
   var_l   ( i) =  7
   var_type( i) = 0
   var_val ( i) = 0
!
   i = i + 1
   var_name( i) = 'REF_INDIV'
   var_char( i) = ' '
   var_l   ( i) =  9
   var_type( i) = 0
   var_val ( i) = 0
!
   i = i + 1
   var_name( i) = 'REF_NINDIV'
   var_char( i) = ' '
   var_l   ( i) = 10
   var_type( i) = 0
   var_val ( i) = 0
!
   i = i + 1
   var_name( i) = 'PI'
   var_char( i) = ' '
   var_l   ( i) = 2
   var_type( i) = 1
   var_val ( i) = 3.1415926535897932384626433832795028841971693993751D0
!
!
   var_num      = i
   var_sys      = 7
!
   END SUBROUTINE variable_init
!
END MODULE variable_mod
