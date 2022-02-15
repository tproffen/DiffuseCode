MODULE refine_fit_erg
!
use precision_mod
IMPLICIT NONE
!
REAL(kind=PREC_DP)  :: refine_chisqr = -1.0  ! Sum of dev^2 in previous cycle
REAL(kind=PREC_DP)  :: refine_conf   = -1.0  ! sum of dev^2
REAL(kind=PREC_DP)  :: refine_lamda  = -1.0  ! Final convergence criterion
REAL(kind=PREC_DP)  :: refine_rval   = -1.0  ! weighted R-value
REAL(kind=PREC_DP)  :: refine_rexp   = -1.0  ! unweighted R-value
!
END MODULE refine_fit_erg
