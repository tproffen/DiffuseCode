PROGRAM diffev 
!                                                                       
USE diffev_setup_mod
USE diffev_loop_mod
USE diffev_mpi_mod
USE run_mpi_mod
!
IMPLICIT none 
!
LOGICAL, PARAMETER :: standalone = .true.
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
CALL diffev_loop
!                                                                       
CALL run_mpi_finalize
!                                                                       
END PROGRAM diffev                            
