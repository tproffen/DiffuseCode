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
PUBLIC discus_read_structure   ! Use discus/read to read a structure or unit cell
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
USE discus_setup_mod
USE kuplot_setup_mod
USE diffev_setup_mod
USE diffev_mpi_mod
USE run_mpi_mod
!
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
WRITE(output_io,'(a)') 'Contol returned to GUI ...'
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
WRITE(output_io,'(a)') 'Contol returned to GUI ...'
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE gui_do_init(prog,line)
!
USE charact_mod
USE doact_mod
USE prompt_mod
USE terminal_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(IN) :: prog
CHARACTER(LEN=*), INTENT(IN) :: line
INTEGER                      :: length
LOGICAL                      :: lend
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
WRITE(output_io,'(a5,a,a5)') COLOR_INFO,line(1:length), COLOR_FG_DEFAULT
IF(.NOT.lblock_read) THEN   ! This is the first DO/IF statement
   CALL do_do_init (line, lend, length)
ELSE
   input_gui = line
   CALL do_insert_line
ENDIF
CALL back_to_suite      ! Go back to the suite
!
length = LEN_TRIM(prog)
END SUBROUTINE gui_do_init
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE gui_do_insert(prog,line)
!
USE charact_mod
USE doact_mod
USE doexec_mod
USE errlist_mod
USE class_macro_internal
USE prompt_mod
USE terminal_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(IN) :: prog
CHARACTER(LEN=*), INTENT(IN) :: line
LOGICAL                      :: lend
!
IF(lblock_read) THEN    ! Only if we are reading into a Do/If block
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
   insert: DO
      IF(level > -1) THEN
         WRITE(output_io,'(a5,a,a5)') COLOR_INFO,line(1:LEN_TRIM(line)), COLOR_FG_DEFAULT
         input_gui = line
         CALL do_insert_line
      ELSE
         EXIT insert
      ENDIF
      IF(.NOT.lmakro) EXIT insert
   ENDDO insert
   IF(level < 0) THEN   ! Reached last enddo/endif, execute block
      lblock_read = .FALSE.
      CALL do_execute_block(lend)
      lblock = .false.
   ENDIF
!
   IF(ier_num /= 0) THEN
      CALL errlist
   ENDIF
!
   CALL back_to_suite      ! Go back to the suite
!
ENDIF
!
!
END SUBROUTINE gui_do_insert
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
! These are generic send and get routines that are used by 
! the programs in the interface to python. They allow
! python to send / get parts of the arrays i[] and r[]
! to / from (discus, diffev, kuplot)
!
! This section used to be an independent file in lib_f90.
! As the pythion interface has moved to an explcit 
! directorty it is no longer needed in lib_f90
!
SUBROUTINE send_i (iin, lower, upper )
!
! The outer routine sends integer valued numbers for i[:]
!
USE param_mod
USE prompt_mod
IMPLICIT NONE
!
INTEGER,                         INTENT(IN) :: lower
INTEGER,                         INTENT(IN) :: upper
INTEGER, DIMENSION(lower:upper), INTENT(IN) :: iin
!
IF( .not. lsetup_done ) THEN    ! If necessary do initial setup
   CALL initialize_suite
ENDIF
!
IF(lower>0 .and. upper<UBOUND(inpara,1) .and. lower <= upper) THEN
   inpara(lower:upper) = iin(lower:upper)
ENDIF
!
END SUBROUTINE send_i
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE send_r (rin, lower, upper )
!
! The outer routine sends real valued numbers for r[:]
!
USE param_mod
USE prompt_mod
IMPLICIT NONE
!
INTEGER,                      INTENT(IN) :: lower
INTEGER,                      INTENT(IN) :: upper
REAL, DIMENSION(lower:upper), INTENT(IN) :: rin
!
IF( .not. lsetup_done ) THEN    ! If necessary do initial setup
   CALL initialize_suite
ENDIF
!
IF(lower>0 .and. upper<UBOUND(rpara,1) .and. lower <= upper) THEN
   rpara(lower:upper) = rin(lower:upper)
ENDIF
!
END SUBROUTINE send_r
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE get_i (iout, lower, upper )
!
! The outer routine gets integer valued numbers from i[:]
!
USE param_mod
USE prompt_mod
IMPLICIT NONE
!
INTEGER,                         INTENT(IN ) :: lower
INTEGER,                         INTENT(IN ) :: upper
INTEGER, DIMENSION(lower:upper), INTENT(OUT) :: iout
!
IF( .not. lsetup_done ) THEN    ! If necessary do initial setup
   CALL initialize_suite
ENDIF
!
IF(lower>0 .and. upper<UBOUND(inpara,1) .and. lower <= upper) THEN
   iout(lower:upper) = inpara(lower:upper) 
ENDIF
!
END SUBROUTINE get_i
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE get_r (rout, lower, upper )
!
! The outer routine gets real valued numbers from r[:]
!
USE param_mod
USE prompt_mod
IMPLICIT NONE
!
INTEGER,                      INTENT(IN ) :: lower
INTEGER,                      INTENT(IN ) :: upper
REAL, DIMENSION(lower:upper), INTENT(OUT) :: rout
!
IF( .not. lsetup_done ) THEN    ! If necessary do initial setup
   CALL initialize_suite
ENDIF
!
IF(lower>0 .and. upper<UBOUND(rpara,1) .and. lower <= upper) THEN
   rout(lower:upper) = rpara(lower:upper) 
ENDIF
!
END SUBROUTINE get_r
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE discus_read_structure(line)
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
END SUBROUTINE discus_read_structure
!
!________KUPLOT_________________________________________________________________
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE kuplot_load(line)
!
!  A first interface that allows to load a data file from python via
!  suite.kuplot_load( python_string )
!  where python string is any of: 
!       load xy, filename
!
USE prompt_mod
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(IN) :: line
INTEGER                      :: length
!
IF( .NOT. lsetup_done) CALL initialize_suite   ! Do we need to initialize?
!
linteractive=.FALSE.    ! Tell get_cmd to get input from input_gui
CALL kuplot_prae        ! Switch to discus section
length = LEN_TRIM(line)
CALL do_load(line, length) ! Call the actual task at hand
CALL back_to_suite      ! Go back to the suite
linteractive=.TRUE.     ! Tell get_cmd to read input from standard I/O
!
END SUBROUTINE kuplot_load
END MODULE suite
