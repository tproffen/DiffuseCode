MODULE suite_python_support
!
!  Generic support routines for the python / discus_suite interface
!
CONTAINS 
!
SUBROUTINE discus_prae
!
!  Switch from suite to discus section
!
USE suite_setup_mod
USE suite_init_mod
USE discus_setup_mod
USE prompt_mod
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
SUBROUTINE back_to_suite
!
!  Go back to the suite level
!
USE suite_setup_mod
USE prompt_mod
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
END MODULE suite_python_support
