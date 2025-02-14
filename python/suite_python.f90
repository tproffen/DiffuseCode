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
lstandalone       = .false.      ! No standalone for DIFFEV, DISCUS, KUPLOT
!
IF( .NOT. lsetup_done ) THEN    ! ! If necessary do initial setup
   CALL run_mpi_init    ! Do the initial MPI configuration
   CALL setup_suite_start     ! Define initial parameter, array values
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
END MODULE suite
