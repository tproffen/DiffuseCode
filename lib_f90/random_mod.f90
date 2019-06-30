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
!
!
END MODULE random_mod
