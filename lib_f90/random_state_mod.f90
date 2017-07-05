MODULE random_state_mod
!+
!
!     This file contains variables for the state of the random number
!     generator
!-
IMPLICIT NONE
PUBLIC
SAVE
!
INTEGER  ::    iff = 0
INTEGER  ::    ix1 = 0
INTEGER  ::    ix2 = 0
INTEGER  ::    ix3 = 0
!
CONTAINS
!
SUBROUTINE random_current(idum_out, iff_out, ix1_out, ix2_out, ix3_out)
!
USE random_mod
!
INTEGER, INTENT(OUT)  :: idum_out
INTEGER, INTENT(OUT)  :: iff_out
INTEGER, INTENT(OUT)  :: ix1_out
INTEGER, INTENT(OUT)  :: ix2_out
INTEGER, INTENT(OUT)  :: ix3_out
!
IF(idum < 0 .OR. iff==0) THEN
   idum_out = idum
   iff_out  = iff
   ix1_out  = 0
   ix2_out  = 0
   ix3_out  = 0
ELSE
   idum_out  = idum
   iff_out   = iff
   ix1_out   = ix1
   ix2_out   = ix2
   ix3_out   = ix3
ENDIF
!
END SUBROUTINE random_current
!
SUBROUTINE ini_ran_ix(np, werte)
!
! Initializes the random sequence or places it at a previous state
!
USE random_mod
USE times_mod
!
IMPLICIT NONE
!
INTEGER           , INTENT(IN) :: np
REAL, DIMENSION(:), INTENT(IN) :: werte
!
IF(np==1) THEN
   idum = NINT(werte(1))         ! Get user value
   IF(idum==0) THEN              ! If negative initialize randomly
      CALL  datum_intrinsic ()   !    by getting time since midnight
      idum = - midnight
   ENDIF
   iff  = 0
ELSE
   idum = 0                      ! User provided three numbers
   iff  = 1
   ix1  = NINT(werte(1))         ! Set points within sequence
   ix2  = NINT(werte(2))         !
   ix3  = NINT(werte(3))         !
ENDIF
!
END SUBROUTINE ini_ran_ix
!
END MODULE random_state_mod
