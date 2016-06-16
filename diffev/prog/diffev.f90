PROGRAM diffev 
!                                                                       
USE diffev_setup_mod
USE diffev_loop_mod
USE diffev_do_exit_mod
USE diffev_mpi_mod
USE run_mpi_mod
!
IMPLICIT none 
!
LOGICAL, PARAMETER :: standalone = .true.
EXTERNAL :: diffev_sigint
!
!*****7*****************************************************************
!                                                                       
!     Main program for DIFFEV                                           
!                                                                       
!     This is the main program for DIFFEV. It sets up most              
!     variables and calls the loop interpreting the commands.           
!                                                                       
!     Authors : R.B. Neder  (reinhard.neder@fau.de)      
!                                                                       
!*****7*****************************************************************
!
run_mpi_myid      = 0
!
CALL diffev_setup(standalone)
CALL run_mpi_init
CALL diffev_set_sub
CALL diffev_set_sub_cost
CALL SIGNAL(2, diffev_sigint)
CALL diffev_loop
!                                                                       
CALL run_mpi_finalize
CALL diffev_do_exit
!                                                                       
END PROGRAM diffev                            
!
SUBROUTINE diffev_sigint
!
!     Handle DIFFEV specific part of a CTRL-C interrupt.
!     CALLED within standalone DIFFEV only
!     This subroutine calls all DIFFEV specific emergency handlers
!     which can also be called from the SUITE
!
USE diffev_do_exit_mod
IMPLICIT NONE
CHARACTER(LEN=1) :: dummy
WRITE(*,*)
WRITE(*,*) ' EMERGENCY Shutdown with USER CTRL-C Interrupt'
!
CALL diffev_emergency_save
CALL diffev_emergency_mpi
CALL exit_all
!
WRITE(*,*) 
WRITE(*,*) ' DIFFEV closed by User Request CTRL-C '
WRITE(*,*) ' For final close down hit ENTER key'
READ(*,'(a)') dummy
STOP        ! Terminate program
!
END SUBROUTINE diffev_sigint
