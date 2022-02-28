MODULE suite_python_support
!
!  Generic support routines for the python / discus_suite interface
!
CONTAINS 
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE suite_prae
!
!  Switch to suite
!
USE suite_setup_mod
USE suite_init_mod
USE suite_set_sub_mod
USE discus_setup_mod
USE kuplot_setup_mod
USE diffev_setup_mod
USE diffev_mpi_mod
USE run_mpi_mod
USE gen_mpi_mod
USE appl_env_mod
!
USE prompt_mod
USE envir_mod
!
IMPLICIT NONE
!
lstandalone       = .false.      ! No standalone for DIFFEV, DISCUS, KUPLOT
!
IF(lsetup_done ) THEN
   pname     = 'suite'
   pname_cap = 'SUITET'
   prompt    = pname
   oprompt   = pname
   CALL program_files ()
   lsetup_done = .TRUE.
ELSE
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
ENDIF
CALL suite_set_sub ()
!
END SUBROUTINE suite_prae
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE discus_prae
!
!  Switch from suite to discus section
!
USE suite_setup_mod
USE suite_init_mod
USE suite_set_sub_mod
USE discus_setup_mod
USE prompt_mod
USE appl_env_mod
IMPLICIT NONE
!
lstandalone       = .false.      ! No standalone for DIFFEV, DISCUS, KUPLOT
!
IF(suite_discus_init) THEN
   pname     = 'discus'
   pname_cap = 'DISCUS'
   prompt    = pname
   oprompt   = pname
   CALL program_files ()
ELSE
   CALL discus_setup   (lstandalone)
   suite_discus_init = .FALSE.
ENDIF
CALL discus_set_sub ()
CALL suite_set_sub_branch
!
END SUBROUTINE discus_prae
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE diffev_prae
!
!  Switch from suite to diffev section
!
USE suite_setup_mod
USE suite_init_mod
USE suite_set_sub_mod
USE diffev_setup_mod
USE prompt_mod
USE appl_env_mod
IMPLICIT NONE
!
lstandalone       = .false.      ! No standalone for DIFFEV, DISCUS, KUPLOT
!
IF(suite_diffev_init) THEN
   pname     = 'diffev'
   pname_cap = 'DIFFEV'
   prompt    = pname
   oprompt   = pname
   CALL program_files ()
ELSE
   CALL diffev_setup   (lstandalone)
   suite_diffev_init = .FALSE.
ENDIF
CALL diffev_set_sub ()
CALL suite_set_sub_branch
!
END SUBROUTINE diffev_prae
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE kuplot_prae
!
!  Switch from suite to discus section
!
USE suite_setup_mod
USE suite_init_mod
USE suite_set_sub_mod
USE kuplot_setup_mod
USE appl_env_mod
USE prompt_mod
IMPLICIT NONE
!
lstandalone       = .false.      ! No standalone for DIFFEV, DISCUS, KUPLOT
!
IF(suite_kuplot_init) THEN
   pname     = 'kuplot'
   pname_cap = 'KUPLOT'
   prompt    = pname
   oprompt   = pname
   CALL program_files ()
ELSE
   CALL kuplot_setup   (lstandalone)
   suite_kuplot_init = .FALSE.
ENDIF
CALL kuplot_set_sub ()
CALL suite_set_sub_branch
!
END SUBROUTINE kuplot_prae
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE back_to_suite
!
!  Go back to the suite level
!
USE suite_setup_mod
USE suite_set_sub_mod
USE prompt_mod
USE appl_env_mod
!
IMPLICIT NONE
!
pname     = 'suite'
pname_cap = 'SUITE'
prompt    = pname
oprompt   = pname
CALL suite_set_sub
CALL program_files ()
!
END SUBROUTINE back_to_suite
!
INTEGER FUNCTION log2int(log_string)
!
IMPLICIT NONE
!
LOGICAL, INTENT(IN) :: log_string
!
IF(log_string) THEN
    log2int = 1
ELSE
   log2int = 0
ENDIF
!
END FUNCTION log2int
!
INTEGER FUNCTION prop2int(property,i)
!
IMPLICIT NONE
!
INTEGER, DIMENSION(0:1), INTENT(IN) :: property
INTEGER                , INTENT(IN) :: i
!
IF(         btest(property(1),i).AND.btest(property(0),i)) THEN
   prop2int = 1
ELSEIF(.NOT.btest(property(1),i).AND.btest(property(0),i)) THEN
   prop2int = 2
ELSE
   prop2int = 0
ENDIF
!
END FUNCTION prop2int
!
!
END MODULE suite_python_support
