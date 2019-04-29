MODULE refine_set_param_mod
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE refine_set_param(npara, params, k, wert)
!
USE calc_expr_mod
!
IMPLICIT NONE
!
INTEGER                           , INTENT(IN) :: npara   ! Number of parameters
CHARACTER(LEN=*)                  , INTENT(IN) :: params  ! Parameter names
INTEGER                           , INTENT(IN) :: k       ! number to be updated
REAL                              , INTENT(IN) :: wert    ! Target value
!
CHARACTER(LEN=1024) :: string   ! dumy string variable
INTEGER             :: lpname   ! Length of a parameter name
INTEGER             :: indxg    ! Location of "=" in string
INTEGER             :: length   ! Length of a string
!
lpname = LEN_TRIM(params)
WRITE(string,'(a,a,G20.12E3)') params(1:lpname), ' = ', wert
indxg = lpname + 2
length = LEN_TRIM(string)
!
CALL do_math (string, indxg, length)
!
END SUBROUTINE refine_set_param
!
!*******************************************************************************
!
END MODULE refine_set_param_mod
