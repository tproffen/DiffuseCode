MODULE suite
!
!  Module to interface python with the discus_suite
!
!  reinhard.neder@fau.de
!  tproffen@ornl.gov
!
!  Version initial test not ready for release
!
!
USE suite_setup_mod
USE suite_loop_mod
USE suite_init_mod
USE suite_set_sub_mod
USE discus_setup_mod
USE kuplot_setup_mod
USE diffev_setup_mod
USE experi_setup_mod
USE diffev_mpi_mod
USE run_mpi_mod
USE gen_mpi_mod
USE charact_mod
USE prompt_mod
USE envir_mod
USE terminal_mod
!
PRIVATE
PUBLIC initialize_suite    ! Initialize the discus_suite as if started directly
PUBLIC execute_macro       ! Execute macro
PUBLIC set_value           ! Sets value of DISCUS variable 
PUBLIC get_value           ! Gets value of DISCUS variable 
PUBLIC get_data            ! Gets data from DISCUS
PUBLIC get_data_3d         ! Gets 3d data from DISCUS
PUBLIC get_data_length     ! Returns length of KUPLOT data set ik
PUBLIC get_data_sets       ! Returns number of loaded KUPLOT data sets
public python_error
PUBLIC flush_output       ! Flush the Fortran output units (for Python output capture)
!  Core interface (migration plan, WP2). These take and return plain
!  integers, strings and integer arrays only, so that they can move to a
!  bind(C) layer later without changing the Python side.
PUBLIC dc_init             ! Initialize in library mode (no STOP on errors)
PUBLIC dc_reset            ! Reset DISCUS, KUPLOT, DIFFEV, REFINE (suite "reset")
PUBLIC dc_clear_error      ! Forget the last recorded error
PUBLIC dc_last_error       ! First error shown since dc_clear_error
PUBLIC dc_seed             ! Seed the random number generator (DISCUS "seed n")
PUBLIC dc_random_size      ! Number of integers in the generator state
PUBLIC dc_get_random_state ! Current generator state
PUBLIC dc_set_random_state ! Put the generator back into a saved state
PUBLIC dc_run_commands     ! Run lines of commands (a macro held in memory, no file)
!
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE initialize_suite()
!
!   Initialization of the discus_suite, to be run only once
!
IMPLICIT NONE
!
INTEGER, PARAMETER :: master = 0 ! Master ID for MPI
EXTERNAL :: suite_sigint
!
gen_mpi_myid      = 0
!
IF( .NOT. lsetup_done ) THEN    ! ! If necessary do initial setup
   CALL run_mpi_init    ! Do the initial MPI configuration
   CALL setup_suite_start     ! Define initial parameter, array values
   CALL setup_suite     ! Define initial parameter, array values
   CALL suite_set_sub
   CALL discus_setup
   CALL kuplot_setup
   CALL diffev_setup
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
linteractive = .false.
!
END SUBROUTINE initialize_suite
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION get_data_sets()
!
! Returns number of data sets in KUPLOT
!
USE kuplot_config
USE kuplot_mod
!
IMPLICIT NONE
!
get_data_sets = iz-1
!
END FUNCTION get_data_sets
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION get_data_length(ik)
!
! Returns length of KUPLOT data set ik
!
USE kuplot_config
USE kuplot_mod
USE errlist_mod
USE lib_errlist_func
!
IMPLICIT NONE
!
INTEGER :: ik
INTEGER :: l = 0
!
IF (ik.le.(iz - 1) .and. ik.ge.1) THEN
   IF (lni(ik)) THEN
      l = nx(ik)*ny(ik)
   ELSE
      l = lenc(ik)
   ENDIF
ELSE
   CALL python_error(-1)
ENDIF
!
get_data_length = l
!
END FUNCTION get_data_length
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE get_data(ik,xpy,ypy,n) 
!
!  Return data from KUPLOT x,y arrays
!
USE kuplot_config
USE kuplot_mod
USE errlist_mod
!
IMPLICIT NONE
!
INTEGER :: ik
INTEGER :: n
REAL, INTENT(OUT) :: xpy(n)
REAL, INTENT(OUT) :: ypy(n)
!
IF (ik.le.(iz - 1) .and. ik.ge.1) THEN
   xpy(1:n) = x(offxy(ik-1)+1:offxy(ik-1)+n)
   ypy(1:n) = y(offxy(ik-1)+1:offxy(ik-1)+n)
ELSE
   CALL python_error(-1)
ENDIF
!
END SUBROUTINE get_data
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE get_data_3d(ik,xpy,ypy,zpy,npx,npy,n)
!
!  Return data from KUPLOT x,y,z arrays
!
USE kuplot_config
USE kuplot_mod
USE errlist_mod
!
IMPLICIT NONE
!
INTEGER :: ik
INTEGER :: n
INTEGER, INTENT(OUT) :: npx(1)
INTEGER, INTENT(OUT) :: npy(1)
REAL,    INTENT(OUT) :: xpy(n)
REAL,    INTENT(OUT) :: ypy(n)
REAL,    INTENT(OUT) :: zpy(n)
!
INTEGER :: ix
INTEGER :: iy
INTEGER :: i
!
IF (ik.le.(iz - 1) .and. ik.ge.1) THEN
   IF(lni(ik)) THEN
      i=1
      DO ix = 1, nx(ik)
      DO iy = 1, ny(ik)
         xpy(i)=x(offxy (ik-1) + ix) 
         ypy(i)=y(offxy (ik-1) + iy) 
         zpy(i)=z(offz(ik-1) + (ix-1) * ny(ik)+iy)
         i=i+1
      ENDDO
      ENDDO
      npx = nx(ik)
      npy = ny(ik)
   ELSE
      CALL python_error(-2)
   ENDIF
ELSE
   CALL python_error(-1)
ENDIF
!
END SUBROUTINE get_data_3d
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE set_value(val, variable)
!
!  Sets the DISCUS suite real variable 'value' to val.
!
USE do_variable_mod
USE lib_errlist_func
USE errlist_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*)     , INTENT(INOUT) :: variable
REAL                 , INTENT(IN)    :: val
!
REAL(kind=PREC_DP)    :: wert
INTEGER               :: laenge
CHARACTER(LEN=1)      :: dummy
INTEGER               :: length
INTEGER, DIMENSION(2) :: substr
!
laenge = LEN_TRIM(variable)
dummy  = " "
wert   = val
!
CALL upd_variable (variable, laenge, wert, dummy, length, substr)
!
IF (ier_num.ne.0) THEN
   CALL errlist
ENDIF
!
END SUBROUTINE set_value
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
REAL FUNCTION get_value (expression)
!
! The outer routine gets the value of a variable (think eval command)
!
USE ber_params_mod
USE lib_errlist_func
USE errlist_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: expression
INTEGER, DIMENSION(1) :: length
REAL(kind=PREC_DP), DIMENSION(1) :: values
!
length = LEN_TRIM(expression)
CALL ber_params (1, expression, length, values, 1)
get_value = values(1)
!
IF (ier_num.ne.0) THEN
   CALL errlist
ENDIF
!
END FUNCTION get_value
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
SUBROUTINE execute_macro(line)
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
CHARACTER(LEN=*), INTENT(INOUT) :: line
!f2py threadsafe
!  (releases the Python GIL while the macro runs, so that a Python thread can
!   show the output as it is written)
!
INTEGER :: length
!
length = LEN_TRIM(line)
IF(line(1:1) == '@' ) THEN
   line = line(2:length)
   length = length - 1
ENDIF
!
WRITE(output_io,'(a5,''@''a,a5)') COLOR_INFO,line(1:length),COLOR_FG_DEFAULT
CALL file_kdo(line,length)
CALL suite_loop()
!
END SUBROUTINE execute_macro
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
SUBROUTINE python_error(num)
!
!  Error when caling from Python
!
USE errlist_mod
USE lib_length
USE terminal_mod
USE prompt_mod
!
IMPLICIT NONE
!
INTEGER :: num
!
CHARACTER(LEN=45)  ERROR(-5:-1)
!
DATA ERROR (-5:-1) /                              &
   &  ' ',                                        & !-5  ! python
   &  ' ',                                        & !-4  ! python
   &  ' ',                                        & !-3  ! python
   &  'Not a 3d dataset',                         & !-2  ! python
   &  'Invalid data set selected'                 & !-1  ! python
/
!
WRITE(error_io,1500) TRIM(color_err),'ER_PYTH',error(num),num,TRIM(color_fg)
1500  FORMAT(a,' ***',a,'*** ',a45,' ***',i4,' ***',a)
!
END SUBROUTINE python_error
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE flush_output()
!
!  Write out anything the Fortran runtime still holds for the output and
!  error units. Python redirects these file descriptors to capture the
!  text; buffered output would otherwise appear later, outside the capture.
!
USE prompt_mod
USE ISO_FORTRAN_ENV, ONLY: OUTPUT_UNIT, ERROR_UNIT   ! plain form: f2py prints a debug line for "USE, INTRINSIC ::"
!
IMPLICIT NONE
!
FLUSH(OUTPUT_UNIT)
FLUSH(ERROR_UNIT)
IF(output_io /= OUTPUT_UNIT .AND. output_io /= ERROR_UNIT) FLUSH(output_io)
IF(error_io  /= OUTPUT_UNIT .AND. error_io  /= ERROR_UNIT) FLUSH(error_io)
!
END SUBROUTINE flush_output
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE dc_init()
!
!  Initialize the suite for use as a library: an error never ends the
!  process (no STOP, also not with "set error,exit" or a read error on
!  standard input); the first error is kept for dc_last_error.
!
USE errlist_mod
!
IMPLICIT NONE
!
ier_api = .TRUE.
CALL initialize_suite()
CALL dc_clear_error()
!
END SUBROUTINE dc_init
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE dc_reset()
!
!  Put DISCUS, KUPLOT, DIFFEV and REFINE back into their start-up state,
!  as the suite command "reset" does. Variables (r[], i[], user variables)
!  and the random number generator are kept.
!
USE discus_reset_all_mod
USE kuplot_reset_mod
USE diffev_reset
USE refine_reset
!
IMPLICIT NONE
!
CHARACTER(LEN=3) :: zeile
INTEGER          :: lp
!
zeile = 'all'
lp    = 3
CALL kuplot_do_reset(zeile, lp)
CALL discus_reset_all
CALL diffev_do_reset
CALL refine_do_reset
CALL dc_clear_error()
!
END SUBROUTINE dc_reset
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE dc_clear_error()
!
USE errlist_mod
!
IMPLICIT NONE
!
ier_last_num  = 0
ier_last_typ  = ER_NONE
ier_last_cls  = ' '
ier_last_text = ' '
ier_last_msg  = ' '
!
END SUBROUTINE dc_clear_error
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE dc_last_error(num, typ, cls, text, msgs)
!
!  The first error shown since dc_clear_error: its number (< 0; 0 if there
!  was none), type (ER_COMM, ER_IO, ... in errlist_mod), type name, the
!  standard message and DISCUS's extra lines (7 x 80 characters).
!
USE errlist_mod
!
IMPLICIT NONE
!
INTEGER           , INTENT(OUT) :: num
INTEGER           , INTENT(OUT) :: typ
CHARACTER(LEN=4)  , INTENT(OUT) :: cls
CHARACTER(LEN=45) , INTENT(OUT) :: text
CHARACTER(LEN=560), INTENT(OUT) :: msgs
!f2py intent(out) num, typ, cls, text, msgs
!
INTEGER :: i
!
num  = ier_last_num
typ  = ier_last_typ
cls  = ier_last_cls
text = ier_last_text
msgs = ' '
DO i = 1, 7
   msgs((i-1)*80+1:i*80) = ier_last_msg(i)
ENDDO
!
END SUBROUTINE dc_last_error
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE dc_seed(value)
!
!  Same as the command "seed value"; value 0 seeds from the clock.
!
USE random_state_mod
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: value
INTEGER, DIMENSION(1) :: iwerte
!
iwerte(1) = value
CALL ini_ran_ix(1, iwerte, 0)
!
END SUBROUTINE dc_seed
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
INTEGER FUNCTION dc_random_size()
!
USE random_state_mod
!
IMPLICIT NONE
!
dc_random_size = random_nseeds()
!
END FUNCTION dc_random_size
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE dc_get_random_state(state, n)
!
USE random_state_mod
!
IMPLICIT NONE
!
INTEGER, INTENT(IN)  :: n
INTEGER, INTENT(OUT) :: state(n)
!f2py intent(out) state
!f2py intent(in) n
INTEGER :: nseeds
!
CALL random_current(nseeds, state)
!
END SUBROUTINE dc_get_random_state
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE dc_set_random_state(state, n)
!
!  Restore a state from dc_get_random_state. Also clears the second
!  Gaussian number that gasdev keeps, as "seed" does.
!
USE random_mod
USE random_state_mod
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: n
INTEGER, INTENT(IN) :: state(n)
!
CALL alloc_random()
CALL put_seeds(n, state)
iset = 0
!
END SUBROUTINE dc_set_random_state
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE dc_run_commands(text, nchar_max, nlines)
!
!  Run command lines (separated by newlines in "text") as a macro held in
!  memory: the lines go into DISCUS's macro storage under a fixed name and
!  are run from there, as "@file" would, but without writing a file. The
!  storage is emptied when the macro ends. nchar_max returns the longest
!  line and nlines the number of lines (the caller checks them against
!  the line length limit, PREC_STRING).
!
USE class_macro_internal
USE errlist_mod
USE precision_mod
USE blanks_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(IN)  :: text
INTEGER         , INTENT(OUT) :: nchar_max
INTEGER         , INTENT(OUT) :: nlines
!f2py intent(out) nchar_max, nlines
!f2py threadsafe
!
CHARACTER(LEN=*), PARAMETER :: NAME = '/diffuse/python/commands.mac'
CHARACTER(LEN=PREC_STRING), DIMENSION(:), ALLOCATABLE :: content
CHARACTER(LEN=PREC_STRING) :: line
INTEGER :: i, start, length, n
!
!  Count the lines (the last one may lack a newline)
n = 0
DO i = 1, LEN(text)
   IF(text(i:i) == ACHAR(10)) n = n + 1
ENDDO
IF(LEN(text) > 0) THEN
   IF(text(LEN(text):LEN(text)) /= ACHAR(10)) n = n + 1
ENDIF
nlines    = n
nchar_max = 0
IF(n == 0) RETURN
ALLOCATE(content(n))
content = ' '
start = 1
n = 0
DO i = 1, LEN(text) + 1
   IF(i > LEN(text)) THEN
      IF(start > LEN(text)) EXIT
   ELSEIF(text(i:i) /= ACHAR(10)) THEN
      CYCLE
   ENDIF
   n = n + 1
   length = i - start
   nchar_max = MAX(nchar_max, length)
   line = ' '
   IF(length > 0) line = text(start:start + MIN(length, PREC_STRING) - 1)
   length = MIN(length, PREC_STRING)
   CALL rem_leading_bl(line, length)          ! as file_kdo does for files
   content(n) = line
   start = i + 1
ENDDO
IF(nchar_max > PREC_STRING) RETURN            ! caller reports the long line
CALL macro_store(NAME, nlines, content)
DEALLOCATE(content)
line = '@' // NAME
CALL execute_macro(line)
!
END SUBROUTINE dc_run_commands
!
END MODULE suite
