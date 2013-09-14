MODULE diff_evol
!-
!      Variables needed to define the differential evolution process
!+
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
REAL               :: diff_cr
REAL               :: diff_f 
REAL               :: diff_k 
REAL               :: diff_local
!
END MODULE diff_evol
