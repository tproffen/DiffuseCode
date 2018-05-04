MODULE diffev_random
!
!  Handles logging of the random number status 
!
!  The status of the random number generator at the start of the 
!  currently best individuum is preserved.
!
IMPLICIT NONE
!
CHARACTER(LEN=100) :: random_macro
CHARACTER(LEN=100) :: random_prog
LOGICAL :: write_random_state = .FALSE.
LOGICAL :: l_get_random_state = .TRUE.
!INTEGER, DIMENSION(:,:), ALLOCATABLE :: random_state  ! Status for current members
INTEGER                              :: random_nseed
INTEGER, DIMENSION(0:64)             :: random_best   ! Status for best    member
!
CONTAINS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE diffev_random_on
!
IMPLICIT NONE
!
l_get_random_state = .TRUE.
!
END SUBROUTINE diffev_random_on
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE diffev_random_off
!
IMPLICIT NONE
!
l_get_random_state = .FALSE.
!
END SUBROUTINE diffev_random_off
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
LOGICAL FUNCTION diffev_random_status()
!
IMPLICIT NONE
!
diffev_random_status = l_get_random_state
!
END FUNCTION diffev_random_status
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE diffev_random_write_on(prog, prog_l, macro, macro_l)
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(IN) :: prog
CHARACTER(LEN=*), INTENT(IN) :: macro
INTEGER         , INTENT(IN) :: prog_l
INTEGER         , INTENT(IN) :: macro_l
!
write_random_state = .TRUE.     ! Turn on  documentation
random_prog  = prog(1:prog_l)
random_macro = macro(1:macro_l)
!
END SUBROUTINE diffev_random_write_on
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE diffev_random_write_off
!
IMPLICIT NONE
!
write_random_state = .FALSE.    ! Turn off documentation
!
END SUBROUTINE diffev_random_write_off
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE diffev_random_save(new)
!
IMPLICIT NONE
!
INTEGER, DIMENSION(0:64), INTENT(IN) :: new
!
random_best(:) = new(:)
random_nseed   = MIN(64, new(0))
!
END SUBROUTINE diffev_random_save
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE diffev_best_macro
!
USE population
USE run_mpi_mod
!
USE random_state_mod
!
IMPLICIT NONE
!
INTEGER, PARAMETER :: IWR = 88
!
CHARACTER(LEN=40) :: macro_file = 'diffev_best.mac'
CHARACTER(LEN=1024) :: line
INTEGER :: i, i1, ir1, ir2, ir3
INTEGER :: nseed_run    ! Actual number of seed used by compiler
!
nseed_run = random_nseeds()
random_nseed   = MIN(RUN_MPI_NSEEDS, nseed_run)  !  to be debugged depend on compiler ???
IF(write_random_state) THEN
   CALL oeffne(IWR, macro_file, 'unknown')
!
   WRITE(IWR,'(a)') random_prog(1:LEN_TRIM(random_prog))
   WRITE(IWR,'(a)') '#@ HEADER'
   WRITE(IWR,'(a)') '#@ NAME         diffev_best.mac'
   WRITE(IWR,'(a)') '#@ '
   WRITE(IWR,'(a)') '#@ KEYWORD      diffev, best member, initialize'
   WRITE(IWR,'(a)') '#@ '
   WRITE(IWR,'(a)') '#@ DESCRIPTION  This macro contains the parameters for the current best'
   WRITE(IWR,'(a)') '#@ DESCRIPTION  member. If run the best member will be recreated.'
   WRITE(IWR,'(a)') '#@ DESCRIPTION  As the random state is explicitely contained as well, the'
   WRITE(IWR,'(a)') '#@ DESCRIPTION  best member will be recreated exactly.'
   WRITE(IWR,'(a)') '#@ DESCRIPTION'
   WRITE(IWR,'(a)') '#@ DESCRIPTION  This macro uses the original macro on the run_mpi command'
   WRITE(IWR,'(a)') '#@ DESCRIPTION  line. Make sure to turn on writing of desired output files.'
   WRITE(IWR,'(a)') '#@'
   WRITE(IWR,'(a)') '#@ PARAMETER    $0, 0'
   WRITE(IWR,'(a)') '#@'
   WRITE(IWR,'(a)') '#@ USAGE        @diffev_best.mac'
   WRITE(IWR,'(a)') '#@'
   WRITE(IWR,'(a)') '#@ END'
   WRITE(IWR,'(a)') '#'
   WRITE(IWR,'(a,i12)') 'REF_GENERATION = ',pop_gen
   WRITE(IWR,'(a,i12)') 'REF_MEMBER     = ',pop_n
   WRITE(IWR,'(a,i12)') 'REF_CHILDREN   = ',pop_c
   WRITE(IWR,'(a,i12)') 'REF_DIMENSION  = ',pop_dimx
   WRITE(IWR,'(a,i12)') 'REF_KID        = ',9999
   WRITE(IWR,'(a,i12)') 'REF_INDIV      = ',9999
   DO i=1,pop_dimx
      WRITE(IWR,'(A,A)') 'variable real, ', pop_name(i)
   ENDDO
   DO i=1,pop_dimx
      WRITE(IWR,'(A,       A,E17.10)') pop_name(i),      ' = ',child(i,pop_best)
!     WRITE(IWR,'(A,I12,A,E17.10)') 'ref_para[',i,'] = ',child(i,pop_best)
   ENDDO
!
!  IF(random_nseed>0) THEN
!     line = ' '
!     line(1:5) = 'seed '
!     DO i=1, random_nseed - 1
!        i1 = 6 + (i-1)*10
!        WRITE(line(i1:i1+9),'(I8,A2)') random_best(i),', '
!     ENDDO
!     i = random_nseed
!     i1 = 6 + (i-1)*10
!     WRITE(line(i1:i1+7),'(I8)') random_best(i)
!     WRITE(IWR,'(a)') line(1:LEN_TRIM(line))
!  ENDIF
!
   IF(random_nseed>0) THEN
      line = ' '
      line(1:5) = 'seed '
      i1 = 6
      DO i=1, random_nseed 
         i1 = 6 + (i-1)*16
         ir1 =         random_best(i)/ 100000000
         ir2 =     MOD(random_best(i), 100000000)/10000
         ir3 =     MOD(random_best(i), 10000)
         WRITE(line(i1:i1+15),'(I4,A1,I4,A1,I4,A2)') ir1,',',ir2,',',ir3,', '
      ENDDO
      WRITE(line(i1+16:i1+23),'(a8)') ' group:3'
      WRITE(IWR,'(a)') line(1:LEN_TRIM(line))
   ENDIF
   WRITE(IWR,'(a)') '#'
   WRITE(IWR,'(a)') '#'
   WRITE(IWR,'(a1,a,a)') '@',random_macro(1:LEN_TRIM(random_macro)),'  ., REF_KID, REF_INDIV'
   WRITE(IWR,'(a)') '#'
   WRITE(IWR,'(a)') 'exit'
!
   CLOSE(IWR)
!
ENDIF
!
CALL diffev_random_write_off    ! Turn off documentation
!
END SUBROUTINE diffev_best_macro
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
END MODULE diffev_random
