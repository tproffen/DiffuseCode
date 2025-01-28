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
PUBLIC get_data            ! Gets data from DISCUS
PUBLIC send_i            
PUBLIC send_r            
PUBLIC get_i            
PUBLIC get_r            
!
CONTAINS
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE get_data(xpy,ypy,n) 
!
!  Return data from DISCUS output command
!
USE kuplot_config
USE kuplot_mod
!
INTEGER :: n
REAL, INTENT(OUT) :: xpy(*)
REAL, INTENT(OUT) :: ypy(*)
!
xpy(1:n) = x(1:n)
ypy(1:n) = y(1:n)
!
END SUBROUTINE get_data
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
IF(lower>0 .and. upper<UBOUND(rpara,1) .and. lower <= upper) THEN
   rout(lower:upper) = rpara(lower:upper)
ENDIF
!
END SUBROUTINE get_r
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
END MODULE suite
