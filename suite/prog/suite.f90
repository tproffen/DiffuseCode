PROGRAM discus_suite 
!                                                                       
USE suite_setup_mod
USE suite_loop_mod
USE suite_init_mod
USE discus_setup_mod
USE kuplot_setup_mod
USE diffev_setup_mod
USE diffev_loop_mod
USE diffev_mpi_mod
USE run_mpi_mod
!
USE prompt_mod
!                                                                       
IMPLICIT none 
!
!*****7*****************************************************************
!                                                                       
!     Main program for DISCUS_SUITE                                           
!                                                                       
!     This is the main program for DISCUS_SUITE. It sets up most              
!     variables and calls the loop interpreting the commands.           
!                                                                       
!     Authors : R.B. Neder  (reinhard.neder@fau.de)      
!                                                                       
!*****7*****************************************************************
!
INTEGER, PARAMETER :: master = 0 ! Master ID for MPI
run_mpi_myid      = 0
lstandalone       = .false.      ! No standalone for DIFFEV, DISCUS, KUPLOT
!lstandalone       = .true.      ! No standalone for DIFFEV, DISCUS, KUPLOT
!
CALL setup_suite     ! Define initial parameter, array values
CALL run_mpi_init    ! Do the initial MPI configuration for slave DIFFEV
CALL suite_set_sub   ! Point to specific subroutines
!
IF(run_mpi_myid /= master) THEN   !  "DIFFEV" slave, directly go to diffev
   CALL suite_set_hlp ()
   CALL diffev_setup   (lstandalone)
   CALL diffev_set_sub ()
   CALL suite_set_sub_cost ()
   CALL diffev_loop    ()
   pname     = 'diffev'
   pname_cap = 'DIFFEV'
   prompt    = pname
ELSE
   CALL discus_setup(lstandalone)
   CALL kuplot_setup(lstandalone)
   CALL diffev_setup(lstandalone)
   suite_discus_init = .TRUE.
   suite_kuplot_init = .TRUE.
   suite_diffev_init = .TRUE.
   pname     = 'suite'
   pname_cap = 'SUITE'
   prompt    = pname
   CALL suite_loop      ! Perform the normal main loop
ENDIF
!                                                                       
CALL run_mpi_finalize
!                                                                       
END PROGRAM discus_suite
