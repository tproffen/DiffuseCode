MODULE suite
!
!  Module to interface python with the discus_suite
!
!  reinhard.neder@fau.de
!
!  Version initial test not ready for release
!
USE suite_python_support
!
PRIVATE
PUBLIC initialize_suite
PUBLIC interactive
PUBLIC py_read_structure
!
CONTAINS
!
SUBROUTINE initialize_suite()
!
!   Initialization of the discus_suite, to be run only once
!
USE suite_setup_mod
USE suite_loop_mod
USE suite_init_mod
USE discus_setup_mod
USE kuplot_setup_mod
USE diffev_setup_mod
USE diffev_mpi_mod
USE run_mpi_mod
!
USE prompt_mod
USE envir_mod
!
IMPLICIT NONE
!
INTEGER, PARAMETER :: master = 0 ! Master ID for MPI
EXTERNAL :: suite_sigint
!
run_mpi_myid      = 0
lstandalone       = .false.      ! No standalone for DIFFEV, DISCUS, KUPLOT
!
IF( .NOT. lsetup_done ) THEN    ! ! If necessary do initial setup
   CALL run_mpi_init    ! Do the initial MPI configuration
   CALL setup_suite     ! Define initial parameter, array values
   CALL suite_set_sub
   CALL discus_setup(lstandalone)
   CALL kuplot_setup(lstandalone)
   CALL diffev_setup(lstandalone)
   suite_discus_init = .TRUE.
   suite_kuplot_init = .TRUE.
   suite_diffev_init = .TRUE.
   pname     = 'suite'
   pname_cap = 'SUITE'
   prompt    = pname
   hlpfile   = hlpdir(1:hlp_dir_l)//pname(1:LEN(TRIM(pname)))//'.hlp'
   hlpfile_l = LEN(TRIM(hlpfile))
   IF(.NOT.run_mpi_active) THEN
      CALL suite_set_sub_cost ()
   ENDIF
   lsetup_done = .TRUE.
ELSE
   CALL suite_set_sub
ENDIF
lstandalone = .false.
!
END SUBROUTINE initialize_suite
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE interactive()
!
!  Generic interface routine to start an interactive discus_suite session
!  from the python host
!
USE suite_loop_mod
USE prompt_mod
!
IMPLICIT NONE
!
IF( .NOT. lsetup_done) CALL initialize_suite
CALL suite_loop      ! Perform the normal main loop
lsetup_done = .TRUE.
!
END SUBROUTINE interactive
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE py_read_structure(line)
!
!  A first interface that allows to read a structre from python via
!  suite.read_cell( python_string )
!  where python string is any of: 
!       cell      crystal_structure.cell, nx, ny, nz
!       lcell     crystal_structure.cell, nx, ny, nz
!       structure crystal_structure.cell
!       free      [optional parameters]
!
USE structur
USE prompt_mod
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(IN) :: line
!
IF( .NOT. lsetup_done) CALL initialize_suite   ! Do we need to initialize?
!
linteractive=.FALSE.    ! Tell get_cmd to get input from input_gui
CALL discus_prae        ! Switch to discus section
input_gui = line        ! copy the input line to the automatic command line
CALL read_struc         ! Call the actual task at hand
CALL back_to_suite      ! Go back to the suite
linteractive=.TRUE.     ! Tell get_cmd to read input from standard I/O
!
END SUBROUTINE py_read_structure
!
END MODULE suite
