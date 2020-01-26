MODULE refine_random_mod
!
! Handles the random state 
!
! Needed to ensure that all callc to the macro for the derivatives use the 
! same set of random seeds at startup.
!
INTEGER                              :: refine_nseed
INTEGER, DIMENSION(0:64)             :: refine_seeds  ! Current seeds
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE refine_save_seeds
!-
!   Get and save the current seeds
!+
!
USE random_state_mod
IMPLICIT NONE
!
refine_nseed = random_nseeds()
!
CALL random_current(refine_nseed, refine_seeds)
!
END SUBROUTINE refine_save_seeds
!
!*******************************************************************************
!
SUBROUTINE refine_restore_seeds
!-
!   Get and save the current seeds
!+
!
USE random_state_mod
IMPLICIT NONE
!
CALL put_seeds(refine_nseed, refine_seeds)
!
END SUBROUTINE refine_restore_seeds
!
!*******************************************************************************
!
END MODULE refine_random_mod
