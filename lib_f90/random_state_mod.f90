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
INTEGER, PARAMETER                :: np = 1
INTEGER, DIMENSION(1)             :: werte
!INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: seed_val
!
CALL RANDOM_SEED(SIZE=nseeds)      ! Get seed size 
IF(.NOT. ALLOCATED(seed_vals)) THEN
CALL alloc_random()                ! Just in case, if random had not been initialized
   werte(:) = 0
   CALL ini_ran_ix(np,werte)
ENDIF
!
CALL RANDOM_SEED(GET=seed_vals)    ! Use the global variable, is allocated to proper size
seed_val(:) = 0                    ! Initialize user array
nseeds = MIN(nseeds,UBOUND(seed_val,1),UBOUND(seed_vals,1))
seed_val(1:nseeds) = seed_vals(1:nseeds)
!
END SUBROUTINE random_current
!
SUBROUTINE ini_ran_ix(np, iwerte)
!
! Initializes the random sequence or places it at a previous state
!
USE random_mod
USE times_mod
!
IMPLICIT NONE
!
INTEGER              , INTENT(IN) :: np
INTEGER, DIMENSION(:), INTENT(IN) :: iwerte
!
INTEGER  :: i
INTEGER  :: nseeds
!INTEGER, DIMENSION(:), ALLOCATABLE :: seed_val
!INTEGER, DIMENSION(12)              :: seed_val
REAL     :: r
!
CALL alloc_random()
nseeds = UBOUND(seed_vals,1)
!
! If one value == 0 initialize automatically
! Else take user values
IF(np==1 .AND. IABS(iwerte(1))==0) THEN
   CALL RANDOM_SEED()                   ! Set at default value
   random_linit = .FALSE.
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
      seed_vals(i) = IABS(iwerte(i))
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
!*****7***********************************************************      
SUBROUTINE do_seed (zeile, lp) 
!                                                                       
USE ber_params_mod
USE errlist_mod 
USE get_params_mod
!USE random_state_mod
USE take_param_mod
!
IMPLICIT none 
!                                                                       
INTEGER, PARAMETER :: maxw = 193 
!                                                                       
CHARACTER (LEN=*), INTENT(INOUT) :: zeile 
INTEGER          , INTENT(INOUT) :: lp
!
CHARACTER(LEN=1024), DIMENSION(MAXW) ::  cpara !(maxw) 
INTEGER            , DIMENSION(MAXW) ::  lpara !(maxw)
INTEGER  :: ianz, np
REAL               , DIMENSION(MAXW) ::  werte !(maxw) 
INTEGER            , DIMENSION(MAXW-1) :: iwerte !(maxw) 
INTEGER, PARAMETER :: NOPTIONAL = 1
CHARACTER(LEN=1024), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=1024), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
REAL               , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 1 ! Number of values to calculate 
INTEGER :: igroup
INTEGER :: i,ind, ip
!
DATA oname  / 'group'  /
DATA loname /  5       /
opara  =  (/ '1.0000'  /)   ! Always provide fresh default values
lopara =  (/  6        /)
owerte =  (/  1.0      /)
!
!                                                                       
werte(:)  = 0.0
iwerte(:) = 0
IF (zeile.ne.' ') THEN 
   CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
   CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                     oname, loname, opara, lopara, owerte)
   igroup = NINT(owerte(1))
!
   IF (ier_num == 0) THEN 
      IF (ianz <= MAXW) THEN 
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         IF (ier_num == 0) THEN 
            IF(MOD(ianz,igroup)==0) THEN
               np = ianz/igroup
               DO i = 1, ianz   ! loop over the current number of parameters
                  ind = (i-1)/igroup+1   ! New parameter index
                  ip  = igroup - MOD(i-1,igroup) - 1   ! Power for 10**ip
                  iwerte(ind) = iwerte(ind) + IABS(NINT(werte(i)))*10000**ip
               ENDDO
               DO i = 1, np
                  ind = (i-1)*igroup + 1
                  IF(werte(ind)<0) iwerte(i) = -iwerte(i)
               ENDDO
               CALL ini_ran_ix(np, iwerte)
            ELSE
               ier_num = - 6 
               ier_typ = ER_COMM 
               ier_msg(1) = 'Nonmatching modulo'
            ENDIF 
         ENDIF 
!           IF (ianz.eq.1) THEN 
!              CALL ber_params (ianz, cpara, lpara, werte, maxw) 
!              IF (ier_num.eq.0) THEN 
!                 iflag = - iabs (nint (werte (1) ) ) 
!                 CALL ini_ran (iflag) 
!              ENDIF 
!           ELSEIF (ianz.eq.3) THEN 
!              CALL ber_params (ianz, cpara, lpara, werte, maxw) 
!              IF (ier_num.eq.0) THEN 
!                 CALL ini_ran_ix(3, iwerte)
!              ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
   ENDIF 
ELSE 
   np = 1
   iwerte(1) = 0
   CALL ini_ran_ix(np, iwerte)
!        CALL ini_ran (0) 
ENDIF 
END SUBROUTINE do_seed                        
!*****7***********************************************************      
!
END MODULE random_state_mod
