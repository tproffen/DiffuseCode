MODULE refine_current_mod
!
CONTAINS
!
SUBROUTINE refine_current_all
!-
!  Calculate current parameter values refined and fixed
!+
USE refine_params_mod
!
IMPLICIT NONE
!
CALL refine_current(refine_par_n, refine_params, refine_p)
CALL refine_current(refine_fix_n, refine_fixed , refine_f)
!
END SUBROUTINE refine_current_all
!
!*******************************************************************************
!
SUBROUTINE refine_current(par_n, params, p)
!-
!  Calculate current parameter values refined or fixed
!+
USE ber_params_mod
USE precision_mod
!
IMPLICIT NONE
!
INTEGER                               , INTENT(IN)  :: par_n
CHARACTER(LEN=*)  , DIMENSION(1:par_n), INTENT(IN)  :: params
REAL(KIND=PREC_SP), DIMENSION(1:par_n), INTENT(OUT) :: p
!
CHARACTER(LEN=1024), DIMENSION(:), ALLOCATABLE :: cpara
INTEGER            , DIMENSION(:), ALLOCATABLE :: lpara
REAL(KIND=PREC_DP) , DIMENSION(:), ALLOCATABLE :: werte
INTEGER :: k
INTEGER :: ianz
INTEGER :: MAXW
!
!  Calculate current parameter values
!
ALLOCATE(cpara(1:par_n))
ALLOCATE(lpara(1:par_n))
ALLOCATE(werte(1:par_n))
DO k=1, par_n
   cpara(k) = params(k)
   lpara(k) = LEN_TRIM(params(k))
ENDDO
ianz = par_n
MAXW = par_n
CALL ber_params(ianz, cpara, lpara, werte, MAXW)
DO k=1, par_n
   p(k) = werte(k)
ENDDO
DEALLOCATE(cpara)
DEALLOCATE(lpara)
DEALLOCATE(werte)
!
END SUBROUTINE refine_current
!
!*******************************************************************************
!
END MODULE refine_current_mod
