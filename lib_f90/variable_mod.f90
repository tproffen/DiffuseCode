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
   INTEGER                                 :: var_num  = 0
   CHARACTER (LEN=200), DIMENSION(VAR_MAX) :: var_name = ' ' ! (VAR_MAX)
   CHARACTER (LEN=200), DIMENSION(VAR_MAX) :: var_char = ' ' ! (VAR_MAX)
   INTEGER            , DIMENSION(VAR_MAX) :: var_l    = 0   ! (VAR_MAX)
   INTEGER            , DIMENSION(VAR_MAX) :: var_type = 0   ! (VAR_MAX)
   REAL               , DIMENSION(VAR_MAX) :: var_val  = 0.0 ! (VAR_MAX)
!
END MODULE variable_mod
