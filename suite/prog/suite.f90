PROGRAM discus_suite 
!                                                                       
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
USE envir_mod
USE errlist_mod
use lib_errlist_func
USE gen_mpi_mod
USE prompt_mod
USE variable_mod
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
EXTERNAL :: suite_sigint
gen_mpi_myid      = 0
lstandalone       = .false.      ! No standalone for DIFFEV, DISCUS, KUPLOT
!lstandalone       = .true.      ! No standalone for DIFFEV, DISCUS, KUPLOT
!
CALL setup_suite_start           ! Do start up procedures
if(ier_num/=0) call errlist
!
CALL run_mpi_init    ! Do the initial MPI configuration for slave DIFFEV
IF(ier_num/=0) THEN
   IF(gen_mpi_myid == master) THEN   !  "DIFFEV" slave, directly go to diffev
      CALL run_mpi_finalize
      WRITE(*,*)
      WRITE(*,'(a)') ' ***MPI *** MPI System could not be initialized            ***'
      IF(ier_msg(1)/=' ') WRITE(*,'(3a)') ' ***MPI *** ',ier_msg(1)(1:46),' ***'
      IF(ier_msg(2)/=' ') WRITE(*,'(3a)') ' ***MPI *** ',ier_msg(2)(1:46),' ***'
      IF(ier_msg(3)/=' ') WRITE(*,'(3a)') ' ***MPI *** ',ier_msg(3)(1:46),' ***'
      WRITE(*,*)
      STOP
   ENDIF
ENDIF
CALL setup_suite     ! Define initial parameter, array values
CALL suite_set_sub   ! Point to specific subroutines
CALL set_signal
!CALL SIGNAL(2, -1)
!CALL SIGNAL(2, suite_sigint)
!
IF(gen_mpi_myid /= master) THEN   !  "DIFFEV" slave, directly go to diffev
   CALL program_files ()
   CALL discus_setup   (lstandalone)
   CALL kuplot_setup   (lstandalone)
   CALL diffev_setup   (lstandalone)
   suite_discus_init = .TRUE.
   suite_kuplot_init = .TRUE.
   suite_diffev_init = .TRUE.
   suite_refine_init = .TRUE.
   CALL diffev_set_sub ()
   CALL suite_set_sub_cost ()
   pname     = 'diffev'
   pname_cap = 'DIFFEV'
   prompt    = pname
   var_val(VAR_STATE)   = var_val(VAR_IS_SECTION)
   var_val(VAR_PROGRAM) = var_val(VAR_DIFFEV)
   IF(gen_mpi_active) THEN
      var_val(VAR_MPI)     = var_val(VAR_MPI_ON)
   ELSE
      var_val(VAR_MPI)     = var_val(VAR_MPI_OFF)
   ENDIF
   CALL diffev_loop    ()
ELSE
   CALL lib_f90_init_updates
   CALL program_files ()
   CALL discus_setup(lstandalone)
   CALL kuplot_setup(lstandalone)
   CALL diffev_setup(lstandalone)
   CALL refine_setup(lstandalone)
   suite_discus_init = .TRUE.
   suite_kuplot_init = .TRUE.
   suite_diffev_init = .TRUE.
   suite_refine_init = .TRUE.
   pname     = 'suite'
   pname_cap = 'SUITE'
   prompt    = pname
   hlpfile   = hlpdir(1:hlp_dir_l)//pname(1:LEN(TRIM(pname)))//'.hlp'
   hlpfile_l = LEN(TRIM(hlpfile))
   var_val(VAR_STATE)   = var_val(VAR_IS_TOP)
   var_val(VAR_PROGRAM) = var_val(VAR_SUITE)
   IF(gen_mpi_active) THEN
      var_val(VAR_MPI)     = var_val(VAR_MPI_ON)
   ELSE
      var_val(VAR_MPI)     = var_val(VAR_MPI_OFF)
   ENDIF
   IF(.NOT.gen_mpi_active) THEN
      CALL suite_set_sub_cost ()
   ENDIF
   CALL suite_loop      ! Perform the normal main loop
ENDIF
!                                                                       
CALL run_mpi_finalize
!                                                                       
END PROGRAM discus_suite
!
!*******************************************************************************
!
RECURSIVE SUBROUTINE suite_sigint
!
USE discus_exit_mod
USE diffev_do_exit_mod
USE diffev_mpi_mod
USE run_mpi_mod
!
USE appl_env_mod
USE class_macro_internal
USE doact_mod
USE errlist_mod
USE exit_mod
USE prompt_mod
USE str_comp_mod
USE sup_mod
USE terminal_mod
!
IMPLICIT NONE
CHARACTER(LEN=8)  :: dummy
INTEGER           :: lbef
INTEGER           :: prompt_keep
!
!
CALL exit_all
CALL set_signal
CALL color_set_scheme(.TRUE.,0)
!
prompt_keep = prompt_status
prompt_status = PROMPT_ON
!
WRITE(*,*) 
WRITE(*,'(a,a,a,a)') TRIM(color_err),' DISCUS SUITE closed by User Request CTRL-C ',TRIM(color_fg),CHAR(7)
WRITE(*,'(a,a,a,a)') TRIM(color_err),' Type the desired keyword and hit the ENTER key',TRIM(color_fg),CHAR(7)
WRITE(*,*) ' '
WRITE(*,*) '          : Empty ENTER is the same as continue.'
WRITE(*,*) ' continue : Continue at the interactive prompt.'
WRITE(*,*) '          : If you are inside a lengthy calculation like '
WRITE(*,*) '          : domain, fourier, mmc, pdf, powder, rmc, ...,'
WRITE(*,*) '          : it may take a while before the prompt shows up.'
WRITE(*,*) '          : No calculation, structure, current folder etc. is reliable! '
WRITE(*,*) ' save     : Save the current structure to EMERGENCY.STRU and continue. '
WRITE(*,*) ' resume   : Pick up at the interrupt, pretend there was no CTRL-C '
WRITE(*,*) '          : and let the calculation continue.'
WRITE(*,*) '          : There is no warranty that the calculation will be correct!'
WRITE(*,*) ' exit     : Stop the discus_suite'
dummy = ' '
WRITE (*, '(1X,A6,'' > '')',advance='no') 'ctrl-c'
READ(*,'(a)') dummy
lbef = LEN_TRIM(dummy)
IF(dummy==' ') THEN
   prompt_status  = PROMPT_ON
   ier_num = -14
   ier_typ = ER_COMM
   ier_ctrlc = .TRUE.   ! Flag to interrupt lengthy calculations
   RETURN
ELSEIF(str_comp (dummy, 'continue', 3, lbef, 8)) THEN
   prompt_status  = PROMPT_ON
   ier_num = -14
   ier_typ = ER_COMM
   ier_ctrlc = .TRUE.   ! Flag to interrupt lengthy calculations
   RETURN
ELSEIF(str_comp (dummy, 'save', 3, lbef, 4)) THEN
   prompt_status  = PROMPT_ON
   CALL discus_emergency_save
   CALL diffev_emergency_save
!
   ier_num = -14
   ier_typ = ER_COMM
   ier_ctrlc = .TRUE.   ! Flag to interrupt lengthy calculations
   RETURN
ELSEIF(str_comp (dummy, 'resume', 3, lbef, 6)) THEN
   prompt_status = prompt_keep
   ier_num = 0
   ier_typ = ER_NONE
   ier_msg(:) = ' '
   RETURN
ELSEIF(str_comp (dummy, 'exit', 3, lbef, 4)) THEN
   CALL diffev_emergency_mpi
   STOP
ELSE
   CALL diffev_emergency_mpi
   STOP
ENDIF
!
STOP
!
END SUBROUTINE suite_sigint
!
!*******************************************************************************
!
SUBROUTINE set_signal
!
EXTERNAL :: suite_sigint
!
!CALL SIGNAL(2,  0)
!CALL SIGNAL(2,  1)
CALL SIGNAL(2, suite_sigint)
!
END SUBROUTINE set_signal

