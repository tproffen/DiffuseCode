MODULE suite
!
!  Module to interface python with the discus_suite
!
!  reinhard.neder@fau.de
!
!  Version initial test not ready for release
!
USE suite_python_support
USE lib_macro_func
USE suite_setup_mod
USE suite_set_sub_mod
USE suite_loop_mod
USE suite_init_mod
!
USE discus_setup_mod
!
USE kuplot_setup_mod
!
USE diffev_setup_mod
USE diffev_setup_sub_mod
USE diffev_loop_mod
USE diffev_mpi_mod
!
USE refine_setup_mod
!
USE appl_env_mod
USE cmdline_args_mod
USE envir_mod
USE errlist_mod
use lib_errlist_func
USE gen_mpi_mod
USE prompt_mod
USE variable_mod
!
PRIVATE
PUBLIC initialize_suite    ! Initialize the discus_suite as if started directly
PUBLIC interactive         ! start an interactive suite session
PUBLIC execute_macro       ! Execute a macro at suite, discus, diffev, kuplot
PUBLIC execute_help        ! Execute the help
PUBLIC execute_command     ! Execute a single command
PUBLIC test_macro_param    ! Test a macro for the number of parameters required by the macro
PUBLIC gui_do_init         ! Initialize a Do/IF block
PUBLIC gui_do_insert       ! Read Commands into a Do/If block and execute once finished
PUBLIC send_i              ! Send an integer array to the suite
PUBLIC send_r              ! Send a real valued array to the suite
PUBLIC get_i               ! Get an integer valued array from the suite
PUBLIC get_r               ! Get a real valued array from the suite
PUBLIC suite_learn         ! Transfer a line for learning
PUBLIC discus_get_wave_number  ! Interface to get number of wave length symbols
PUBLIC discus_get_wave_symbol  ! Interface to get wave length symbols and values
PUBLIC discus_get_spcgr_number ! Interface to get number of wave length symbols
PUBLIC discus_get_spcgr_symbol ! Interface to get wave length symbols and values
PUBLIC discus_get_nscat_number  ! Interface to get number of wave length symbols
PUBLIC discus_get_scat_symbol  ! Interface to get wave length symbols and values
PUBLIC discus_get_natoms       ! Interface to get number of atoms in the structure
PUBLIC discus_read_structure   ! Use discus/read to read a structure or unit cell
PUBLIC discus_calc_fourier     ! Use discus/fourier to calculate a Fourier
PUBLIC discus_get_fourier      ! Interface to get Fourier menu items from DISCUS
PUBLIC discus_get_four_last    ! Interface to get last Fourier type
PUBLIC discus_set_fourier      ! Interface to set Fourier menu items to   DISCUS
PUBLIC discus_get_powder       ! Interface to get Powder  menu items from DISCUS
PUBLIC discus_calc_powder      ! Use discus/do_powder to calculate a Powder
PUBLIC discus_get_pdf          ! Interface to get PDF     menu items from DISCUS
PUBLIC discus_calc_pdf         ! Use discus/pdf to calculate a PDF
PUBLIC discus_get_save         ! Use discus/pdf to save a crystal structure
PUBLIC discus_run_save         ! Use discus/pdf to save a crystal structure
PUBLIC discus_output           ! Interface to run OUTPUT
PUBLIC discus_run_sro          ! Interface to run SRO
PUBLIC kuplot_load             ! Use kuplot/load to load a data set
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
USE suite_set_sub_mod
USE discus_setup_mod
USE kuplot_setup_mod
USE diffev_setup_mod
USE diffev_mpi_mod
USE run_mpi_mod
!
USE gen_mpi_mod
USE charact_mod
USE prompt_mod
USE envir_mod
USE terminal_mod
!
IMPLICIT NONE
!
INTEGER, PARAMETER :: master = 0 ! Master ID for MPI
EXTERNAL :: suite_sigint
!
gen_mpi_myid      = 0
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
   IF(.NOT.gen_mpi_active) THEN
      CALL suite_set_sub_cost ()
   ENDIF
   lsetup_done = .TRUE.
ELSE
   CALL suite_set_sub
ENDIF
lstandalone = .false.
WRITE(output_io,'(a5,a,a5)') COLOR_HIGH,'Control turned to GUI ...',COLOR_FG_DEFAULT
!
END SUBROUTINE initialize_suite
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE interactive(prog)
!
!  Generic interface routine to start an interactive discus_suite session
!  from the python host
!
USE suite_loop_mod
USE discus_loop_mod
USE diffev_loop_mod
USE kuplot_loop_mod
USE prompt_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(IN)    :: prog
!
IF( .NOT. lsetup_done) CALL initialize_suite
linteractive = .TRUE.
section: SELECT CASE (prog)
   CASE ('suite')
      CALL suite_prae
      CALL suite_loop      ! Perform the normal main loop
   CASE ('discus')
      CALL discus_prae
      CALL discus_loop     ! Perform the normal discus loop
   CASE ('diffev')
      CALL diffev_prae
      CALL diffev_loop     ! Perform the normal discus loop
   CASE ('kuplot')
      CALL kuplot_prae
      CALL kuplot_loop     ! Perform the normal discus loop
END SELECT section
lsetup_done = .TRUE.
WRITE(output_io,'(a)') 'Control returned to GUI ...'
CALL back_to_suite      ! Go back to the suite
!
END SUBROUTINE interactive
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE execute_macro(prog, line)
!
!  Execute the macro given on line for the section defined by prog
!  The macro name and all parameters must be specified, parameters
!  must be separated from each other by a comma.
!
USE suite_loop_mod
USE discus_loop_mod
USE diffev_loop_mod
USE kuplot_loop_mod
!
USE charact_mod
USE prompt_mod
USE terminal_mod
USE lib_macro_func
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(IN)    :: prog
CHARACTER(LEN=*), INTENT(INOUT) :: line
!
INTEGER :: length
!
linteractive = .FALSE.
section: SELECT CASE (prog)
   CASE ('suite')
      CALL suite_prae
   CASE ('discus')
      CALL discus_prae
   CASE ('diffev')
      CALL diffev_prae
   CASE ('kuplot')
      CALL kuplot_prae
END SELECT section
!
length = LEN_TRIM(line)
IF(line(1:1) == '@' ) THEN
   line = line(2:length)
   length = length - 1
ENDIF
!
WRITE(output_io,'(a5,''@''a,a5)') COLOR_INFO,line(1:length),COLOR_FG_DEFAULT
CALL file_kdo(line,length)
exec: SELECT CASE (prog)
   CASE ('suite')
      CALL suite_loop
   CASE ('discus')
      CALL discus_loop
   CASE ('diffev')
      CALL diffev_loop
   CASE ('kuplot')
      CALL kuplot_loop
END SELECT exec
!
linteractive = .TRUE.
WRITE(output_io,'(a)') 'Control returned to GUI ...'
CALL back_to_suite      ! Go back to the suite
!
END SUBROUTINE execute_macro
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE execute_help(prog)
!
!  Execute the help for the section defined by prog
!
!
USE charact_mod
USE prompt_mod
USE terminal_mod
USE lib_help
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(IN)    :: prog
!
INTEGER :: length
!
linteractive = .TRUE.
section: SELECT CASE (prog)
   CASE ('suite')
      CALL suite_prae
   CASE ('discus')
      CALL discus_prae
   CASE ('diffev')
      CALL diffev_prae
   CASE ('kuplot')
      CALL kuplot_prae
END SELECT section
!
length = LEN_TRIM(prog)
CALL do_hel(prog, length)
!
linteractive = .TRUE.
WRITE(output_io,'(a5,a,a5)') COLOR_HIGH,'Control returned to GUI ...',COLOR_FG_DEFAULT
CALL back_to_suite      ! Go back to the suite
!
END SUBROUTINE execute_help
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE execute_command(prog,line)
!
!   Executes a single command as given on line
!   Returns to the GUI immediately
!
USE discus_loop_mod
USE diffev_loop_mod
USE kuplot_loop_mod
USE suite_loop_mod
USE charact_mod
USE prompt_mod
USE terminal_mod
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(IN) :: prog
CHARACTER(LEN=*), INTENT(IN) :: line
INTEGER                      :: length
LOGICAL                      :: lend
!
IF( .NOT. lsetup_done) CALL initialize_suite   ! Do we need to initialize?
!
length = LEN_TRIM(line)
lend   = .FALSE.
linteractive = .FALSE.
input_gui = line
WRITE(output_io,'(a5,a,a5)') COLOR_INFO,line(1:length),COLOR_FG_DEFAULT
section: SELECT CASE (prog)
   CASE ('suite')
      CALL suite_prae
      CALL suite_loop
   CASE ('discus')
      CALL discus_prae
      CALL discus_loop
   CASE ('diffev')
      CALL diffev_prae
      CALL diffev_loop
   CASE ('kuplot')
      CALL kuplot_prae
      CALL kuplot_loop
END SELECT section
CALL back_to_suite      ! Go back to the suite
linteractive=.TRUE.     ! Tell get_cmd to read input from standard I/O
CALL back_to_suite      ! Go back to the suite
!
END SUBROUTINE execute_command
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION test_macro_param(line)
!
CHARACTER(LEN=*), INTENT(IN) :: line
!
INTEGER :: length
INTEGER :: numpar
!
length = LEN_TRIM(line)
CALL test_macro(line,length, numpar)
test_macro_param = numpar
!
END FUNCTION test_macro_param
!
END MODULE suite
