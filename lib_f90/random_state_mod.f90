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
!
CONTAINS
!
INTEGER FUNCTION random_nseeds()
!
IMPLICIT NONE
!
INTEGER :: nseeds
CALL RANDOM_SEED(SIZE=nseeds)        ! Get seed size
random_nseeds = nseeds
!
END FUNCTION random_nseeds
!
SUBROUTINE random_current(nseeds, seed_val)
!
IMPLICIT NONE
!
INTEGER, INTENT(OUT) :: nseeds
INTEGER, DIMENSION(12),INTENT(OUT) :: seed_val
!NTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: seed_val
!
CALL RANDOM_SEED(SIZE=nseeds)        ! Get seed size 
!IF(ALLOCATED(seed_val)) DEALLOCATE(seed_val)
!ALLOCATE(seed_val(1:nseeds))
CALL RANDOM_SEED(GET=seed_val)
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
INTEGER  :: i
INTEGER  :: nseeds
!INTEGER, DIMENSION(:), ALLOCATABLE :: seed_val
INTEGER, DIMENSION(12)              :: seed_val
REAL     :: r
!
CALL RANDOM_SEED()                   ! Set at default value
CALL RANDOM_SEED(SIZE=nseeds)        ! Get seed size 
!ALLOCATE(seed_val(1:nseeds))
!
! If one value == 0 initialize automatically
! Else take user values
IF(np==1 .AND. IABS(NINT(werte(1)))==0) THEN
      CALL  datum_intrinsic ()       ! get time since midnight
      idum =   midnight              ! idum is preserved for backwards compatibility
      seed_val(:) = midnight         ! Set all seeds
      CALL RANDOM_SEED(PUT=seed_val) 
      DO i=1,nseeds                  ! "randomly" populate all seeds
         CALL RANDOM_NUMBER(r)
         seed_val(i) = INT(midnight*r) + 1
      ENDDO
      CALL RANDOM_SEED(PUT=seed_val)
ELSE                                 ! more than one value or non-zero value
   DO i=1, MIN(np,nseeds)            ! Loop over all seeds or all provided values
      seed_val(i) = IABS(NINT(werte(i)))
   ENDDO
   CALL RANDOM_SEED(PUT=seed_val)
!  idum = 0                      ! User provided three numbers
ENDIF
!
END SUBROUTINE ini_ran_ix
!
END MODULE random_state_mod
