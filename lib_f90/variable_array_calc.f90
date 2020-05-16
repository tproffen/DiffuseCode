MODULE variable_array_calc_mod
!
CONTAINS
!
SUBROUTINE var_arr_mul(line, length)
!
USE get_params_mod
USE ber_params_mod
USE precision_mod
USE variable_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: line
INTEGER         , INTENT(INOUT) :: length
!
INTEGER, PARAMETER :: MAXP = 3
CHARACTER(LEN=MAX(PREC_STRING,LEN(line))), DIMENSION(MAXP) :: cpara
INTEGER            , DIMENSION(MAXP) :: lpara
REAL(KIND=PREC_DP) , DIMENSION(MAXP) :: werte
INTEGER  :: ianz
!
CALL get_params (line, ianz, cpara, lpara, MAXP, length)
!
END SUBROUTINE var_arr_mul
END MODULE variable_array_calc_mod
