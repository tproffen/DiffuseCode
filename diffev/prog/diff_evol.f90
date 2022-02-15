MODULE diff_evol
!-
!      Variables needed to define the differential evolution process
!+
use precision_mod
SAVE
PUBLIC
!
INTEGER, PARAMETER :: ADD_TO_RANDOM = 0
INTEGER, PARAMETER :: ADD_TO_BEST   = 1
!
INTEGER, PARAMETER :: SEL_COMP       = 0
INTEGER, PARAMETER :: SEL_BEST_ALL   = 1
INTEGER, PARAMETER :: SEL_BEST_CHILD = 2
!
INTEGER            :: diff_donor_mode
INTEGER            :: diff_sel_mode
!
REAL(kind=PREC_DP) :: diff_cr
REAL(kind=PREC_DP) :: diff_f 
REAL(kind=PREC_DP) :: diff_k 
REAL(kind=PREC_DP) :: diff_local
!
END MODULE diff_evol
