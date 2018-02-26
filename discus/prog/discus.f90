PROGRAM discus 
!                                                                       
USE discus_setup_mod
USE discus_loop_mod
USE variable_mod
!
IMPLICIT none 
!
LOGICAL, PARAMETER :: standalone = .true.
!
EXTERNAL :: discus_sigint
!
CALL discus_setup(standalone)
CALL discus_set_sub
CALL SIGNAL(2, discus_sigint)
var_val(VAR_STATE)   = var_val(VAR_IS_TOP)
var_val(VAR_PROGRAM) = var_val(VAR_DISCUS)
CALL discus_loop
!                                                                       
END PROGRAM discus                            
!
SUBROUTINE discus_sigint
!
!     Handle DISCUS specific part of a CTRL-C interrupt.
!     CALLED within standalone DISCUS only
!     This subroutine calls all DISCUS specific emergency handlers
!     which can also be called from the SUITE
!
USE discus_exit_mod
IMPLICIT NONE
CHARACTER(LEN=1) :: dummy
WRITE(*,*)
WRITE(*,*) ' EMERGENCY Shutdown with USER CTRL-C Interrupt'
!
CALL discus_emergency_save
CALL exit_all
!
WRITE(*,*) 
WRITE(*,*) ' DISCUS closed by User Request CTRL-C '
WRITE(*,*) ' For final close down hit ENTER key'
READ(*,'(a)') dummy
STOP        ! Terminate program
!
END SUBROUTINE discus_sigint
