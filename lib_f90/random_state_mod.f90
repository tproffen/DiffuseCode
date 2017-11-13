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
LOGICAL                            :: random_linit = .FALSE.  ! Is not yet initialized
INTEGER                            :: nseed_comp       ! Compiler size of seeds
INTEGER, DIMENSION(:), ALLOCATABLE :: seed_vals        ! Compiler seeds
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
INTEGER, DIMENSION(:),INTENT(OUT) :: seed_val
!INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: seed_val
!
CALL alloc_random()                ! Just in case, if random had not been initialized
CALL RANDOM_SEED(SIZE=nseeds)      ! Get seed size 
!
CALL RANDOM_SEED(GET=seed_vals)    ! Use the global variable, is allocated to proper size
seed_val(:) = 0                    ! Initialize user array
nseeds = MIN(nseeds,UBOUND(seed_val,1),UBOUND(seed_vals,1))
seed_val(1:nseeds) = seed_vals(1:nseeds)
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
!INTEGER, DIMENSION(12)              :: seed_val
REAL     :: r
!
CALL RANDOM_SEED()                   ! Set at default value
random_linit = .FALSE.
CALL alloc_random()
!
! If one value == 0 initialize automatically
! Else take user values
IF(np==1 .AND. IABS(NINT(werte(1)))==0) THEN
      CALL  datum_intrinsic ()       ! get time since midnight
      idum =   midnight              ! idum is preserved for backwards compatibility
      seed_vals(:) = midnight        ! Set all seeds
      CALL RANDOM_SEED(PUT=seed_vals) 
      DO i=1,nseeds                  ! "randomly" populate all seeds
         CALL RANDOM_NUMBER(r)
         seed_vals(i) = INT(midnight*r) + 1
      ENDDO
      CALL RANDOM_SEED(PUT=seed_vals)
ELSE                                 ! more than one value or non-zero value
   DO i=1, MIN(np,nseeds)            ! Loop over all seeds or all provided values
      seed_vals(i) = IABS(NINT(werte(i)))
   ENDDO
   CALL RANDOM_SEED(PUT=seed_vals)
!  idum = 0                      ! User provided three numbers
ENDIF
!
END SUBROUTINE ini_ran_ix
!
SUBROUTINE alloc_random()
!
IMPLICIT NONE
!
nseed_comp = 1
CALL RANDOM_SEED(SIZE=nseed_comp)        ! Get seed size 
IF(.NOT.ALLOCATED(seed_vals)) THEN
   ALLOCATE(seed_vals(1:nseed_comp))
   seed_vals(:) = 0
ELSEIF(UBOUND(seed_vals,1)/=nseed_comp) THEN
   DEALLOCATE(seed_vals)
   ALLOCATE(seed_vals(1:nseed_comp))
   seed_vals(:) = 0
ENDIF
!
END SUBROUTINE alloc_random

!
END MODULE random_state_mod
