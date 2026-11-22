MODULE random_mod
!+
!
!     This file contains variables for seed for random number
!     generator and the definition for the random number generator
!     function.
!-
   IMPLICIT NONE
   PUBLIC
   SAVE
!
   INTEGER  ::    idum = -182783467
   INTEGER  ::    iset = 0
!  gasdev makes Gaussian numbers in pairs; with iset = 1 the second one is
!  waiting in gset. Here (not a SAVEd local of gasdev) so that the complete
!  generator state can be saved and restored (Python: random_state()).
   REAL(KIND=KIND(1.0D0)) :: gset = 0.0D0
!
!
END MODULE random_mod
