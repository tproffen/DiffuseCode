SUBROUTINE suite_top(from_section)
!
USE suite_loop_mod
USE suite_init_mod
USE suite_set_sub_mod
!
USE refine_setup_mod
USE kuplot_setup_mod
USE fit_set_sub_mod
!
USE appl_env_mod
USE errlist_mod
USE prompt_mod
USE variable_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(IN) :: from_section
!
INTEGER :: ier_num_copy
INTEGER :: ier_typ_copy
LOGICAL :: lend
!
ier_num_copy = 0
ier_typ_copy = 0
lend   = .false.
pname  = 'suite'
pname_cap = 'SUITE'
prompt = pname
oprompt   = pname
var_val(VAR_STATE)   = var_val(VAR_IS_TOP)
var_val(VAR_PROGRAM) = var_val(VAR_SUITE)
CALL suite_set_sub
IF(ier_num == -9 .AND. ier_typ == 1) THEN
   CALL program_files ()
   ier_num = -9
   ier_typ = ER_COMM
ELSE
   CALL program_files ()
ENDIF
!
l_to_top = .TRUE.
CALL suite_loop
l_to_top = .FALSE.
IF(ier_num/=0) THEN
   ier_num_copy = ier_num
   ier_typ_copy = ier_typ
ENDIF
!
IF(from_section=='kuplot_fit') THEN    !Return to KUPLOT
!
   IF(suite_kuplot_init) then
      pname     = 'kuplot'
      pname_cap = 'KUPLOT'
      prompt    = pname
      oprompt   = pname
      CALL program_files ()
   ELSE
      CALL kuplot_setup   (lstandalone)
      suite_kuplot_init = .TRUE.
   ENDIF
   var_val(VAR_STATE)   = var_val(VAR_IS_SECTION)
   var_val(VAR_PROGRAM) = var_val(VAR_KUPLOT)
   CALL kuplot_set_sub ()
   CALL fit_set_sub
   CALL suite_set_sub_branch
ELSEIF(from_section=='kuplot') THEN    !Return to KUPLOT
!
   IF(suite_kuplot_init) then
      pname     = 'kuplot'
      pname_cap = 'KUPLOT'
      prompt    = pname
      oprompt   = pname
      CALL program_files ()
   ELSE
      CALL kuplot_setup   (lstandalone)
      suite_kuplot_init = .TRUE.
   ENDIF
   var_val(VAR_STATE)   = var_val(VAR_IS_SECTION)
   var_val(VAR_PROGRAM) = var_val(VAR_KUPLOT)
   CALL kuplot_set_sub ()
   CALL suite_set_sub_branch
ELSEIF(from_section=='refine_fit') THEN    !Return to KUPLOT
!
   IF(suite_refine_init) then
      pname     = 'refine'
      pname_cap = 'REFINE'
      prompt    = pname
      oprompt   = pname
      CALL program_files ()
   ELSE
      CALL refine_setup   (lstandalone)
      suite_refine_init = .TRUE.
   ENDIF
   var_val(VAR_STATE)   = var_val(VAR_IS_SECTION)
   var_val(VAR_PROGRAM) = var_val(VAR_REFINE)
   CALL refine_set_sub ()
!  CALL fit_set_sub
   CALL suite_set_sub_branch
ENDIF
!
IF(ier_num_copy/=0) THEN ! hand back errors that occured while in suite_loop
   ier_num = ier_num_copy
   ier_typ = ier_typ_copy
ENDIF
!
END SUBROUTINE suite_top
