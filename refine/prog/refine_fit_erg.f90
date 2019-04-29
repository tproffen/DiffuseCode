MODULE refine_fit_erg
!
IMPLICIT NONE
!
REAL  :: refine_chisqr = -1.0  ! Sum of dev^2 in previous cycle
REAL  :: refine_conf   = -1.0  ! sum of dev^2
REAL  :: refine_lamda  = -1.0  ! Final convergence criterion
REAL  :: refine_rval   = -1.0  ! weighted R-value
REAL  :: refine_rexp   = -1.0  ! unweighted R-value
!
END MODULE refine_fit_erg
