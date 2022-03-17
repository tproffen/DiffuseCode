RECURSIVE SUBROUTINE suite_branch(zeile, length, lreset, lloop)
!
!  Specific SUITE Version of a branch subroutine
!  Call a section via a branch
!
USE diffev_setup_mod
USE diffev_setup_sub_mod
USE diffev_loop_mod
!
USE discus_setup_mod
USE discus_setup_sub_mod
USE discus_loop_mod
USE discus_reset_all_mod
!
USE kuplot_setup_mod
USE kuplot_setup_sub_mod
USE kuplot_loop_mod
use kuplot_reset_mod
!
USE refine_setup_mod
USE refine_setup_sub_mod
USE refine_loop_mod
!
USE blanks_mod
USE appl_env_mod
USE errlist_mod
USE prompt_mod
USE suite_init_mod
USE suite_setup_mod
USE suite_set_sub_mod
USE charact_mod
USE lib_macro_func
USE precision_mod
USE str_comp_mod
USE variable_mod
!
IMPLICIT NONE
!
CHARACTER (LEN=*), INTENT(IN) :: zeile
INTEGER          , INTENT(IN) :: length
LOGICAL          , INTENT(IN) :: lreset
integer          , INTENT(IN) :: lloop 
!
CHARACTER(LEN=MAX(PREC_STRING,LEN(zeile)))  :: line
CHARACTER(LEN= 7   ) :: br_pname_old,br_pname_cap_old
INTEGER              :: br_prompt_status_old
INTEGER              :: br_ier_sta_old, br_state_old, br_progr_old
INTEGER              :: indxt, indxb, indxm
INTEGER              :: lbef, laenge
INTEGER              :: lcomm
LOGICAL              :: lmacro
!
if(lloop>=0) then
!
IF(str_comp(zeile, pname, LEN_TRIM(pname), length, LEN_TRIM(pname))) THEN
   ier_num = -15         ! branch to identical section
   ier_typ = ER_COMM
   ier_msg(1) = 'A branch from section '// pname(1:LEN_TRIM(pname)) // &
                ' to ' // pname(1:LEN_TRIM(pname))
   RETURN
ENDIF
!
!
! Get parameters on command line
!
lmacro = .false.
indxt  = INDEX (zeile, tab)       ! find a tabulator
IF(indxt==0) indxt = length + 1
indxb  = index (zeile, ' ') 
IF(indxb==0) indxb = length + 1
indxb  = MIN(indxb,indxt)         ! find first "white" character
lbef   = min (indxb - 1, 9) 
!
line  = ' ' 
lcomm = 0 
IF (indxb + 1.le.length) THEN     ! There is space for parameters
   line  = zeile(indxb + 1:length) 
   lcomm = length - indxb 
   CALL rem_leading_bl(line , lcomm)
ENDIF 
!
indxm = INDEX(zeile,'-macro')
IF(indxm > 0) THEN                   ! Found '-macro' qualifier
   lmacro = .true.
   indxt  = INDEX (line, tab)       ! find a tabulator
   IF(indxt==0) indxt = lcomm + 1
   indxb  = index (line, ' ') 
   IF(indxb==0) indxb = lcomm + 1
   indxb  = MIN(indxb,indxt)
   lbef   = min (indxb - 1, 9) 
   IF (indxb + 1.le.lcomm) THEN      ! There is space for a macro name
      line = line(indxb+1:lcomm)     ! copy macro name and parameter list
      CALL rem_leading_bl(line , lcomm)
      lcomm = LEN_TRIM(line)
   ENDIF
ENDIF
!
br_pname_old         = pname         ! Store old prompt/error status
br_pname_cap_old     = pname_cap
br_prompt_status_old = prompt_status
br_ier_sta_old       = ier_sta
br_state_old         = NINT(var_val(VAR_STATE))
br_progr_old         = NINT(var_val(VAR_PROGRAM))
!
IF(str_comp(zeile, 'kuplot', 2, length, 6)) THEN
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
   CALL kuplot_set_sub ()
   CALL suite_set_sub_branch ()
   IF(lreset) THEN
      line = 'all'
      laenge = 3
      call kuplot_do_reset(line, length)
   ENDIF
   var_val(VAR_PROGRAM) = var_val(VAR_KUPLOT)
   var_val(VAR_STATE)   = var_val(VAR_IS_BRANCH)
   IF(lmacro) THEN           ! Execute "command line macro"
      line = line(1:lcomm)
      CALL file_kdo(line,lcomm)
   ENDIF
   IF(lloop==0 .and. ier_num == 0) THEN     ! If no error in macro do interactive session
       CALL kuplot_loop    ()
   ENDIF
ELSEIF(str_comp(zeile, 'discus', 3, length, 6)) THEN
   IF(suite_discus_init) then
      pname     = 'discus'
      pname_cap = 'DISCUS'
      prompt    = pname
      oprompt   = pname
      CALL program_files ()
   ELSE
      CALL discus_setup   (lstandalone)
      suite_discus_init = .TRUE.
   ENDIF
   CALL discus_set_sub ()
   CALL suite_set_sub_branch ()
   IF(lreset) CALL discus_reset_all
   var_val(VAR_PROGRAM) = var_val(VAR_DISCUS)
   var_val(VAR_STATE)   = var_val(VAR_IS_BRANCH)
   IF(lmacro) THEN           ! Execute "command line macro"
      CALL file_kdo(line(1:lcomm),lcomm)
   ENDIF
   IF(lloop==0 .and. ier_num == 0) THEN     ! If no error in macro do interactive session
       CALL discus_loop    ()
   ENDIF
ELSEIF(str_comp(zeile, 'diffev', 3, length, 6)) THEN
   IF(suite_diffev_init) then
      pname     = 'diffev'
      pname_cap = 'DIFFEV'
      prompt    = pname
      oprompt   = pname
      CALL program_files ()
   ELSE
      CALL diffev_setup   (lstandalone)
      suite_diffev_init = .TRUE.
   ENDIF
   CALL diffev_set_sub ()
   CALL suite_set_sub_branch ()
   var_val(VAR_PROGRAM) = var_val(VAR_DIFFEV)
   var_val(VAR_STATE)   = var_val(VAR_IS_BRANCH)
   IF(lmacro) THEN           ! Execute "command line macro"
      CALL file_kdo(line(1:lcomm),lcomm)
   ENDIF
   IF(lloop==0 .and. ier_num == 0) THEN     ! If no error in macro do interactive session
       CALL diffev_loop    ()
   ENDIF
ELSEIF(str_comp(zeile, 'refine', 3, length, 6)) THEN
   IF(suite_diffev_init) then
      pname     = 'refine'
      pname_cap = 'REFINE'
      prompt    = pname
      oprompt   = pname
      CALL program_files ()
   ELSE
      CALL refine_setup   (lstandalone)
      suite_refine_init = .TRUE.
   ENDIF
   CALL refine_set_sub ()
   CALL suite_set_sub_branch ()
   var_val(VAR_PROGRAM) = var_val(VAR_REFINE)
   var_val(VAR_STATE)   = var_val(VAR_IS_BRANCH)
   IF(lmacro) THEN           ! Execute "command line macro"
      CALL file_kdo(line(1:lcomm),lcomm)
   ENDIF
   IF(lloop==0 .and. ier_num == 0) THEN     ! If no error in macro do interactive session
       CALL refine_loop    ()
   ENDIF
ELSE
   ier_num = -7
   ier_typ = ER_COMM
ENDIF
endif
!
if(lloop<=0) then
   pname      = br_pname_old
   pname_cap  = br_pname_cap_old
   prompt     = pname
   oprompt    = pname
   prompt_status = br_prompt_status_old
   ier_sta       = br_ier_sta_old
   var_val(VAR_PROGRAM) = br_progr_old
   var_val(VAR_STATE)   = br_state_old
!
IF(pname=='discus') THEN      ! Return to DISCUS branch
      CALL discus_set_sub ()
ELSEIF(pname=='diffev') THEN  ! Return to DIFFEV branch
      CALL diffev_set_sub ()
ELSEIF(pname=='kuplot') THEN  ! Return to KUPLOT branch
      CALL kuplot_set_sub ()
ELSEIF(pname=='refine') THEN  ! Return to REFINE branch
      CALL refine_set_sub ()
ENDIF
CALL suite_set_sub_branch ()
IF(ier_num == -9 .AND. ier_typ == 1) THEN
   CALL program_files ()
   ier_num = -9
   ier_typ = ER_COMM
ELSE
   CALL program_files ()
ENDIF
endif
!
END SUBROUTINE suite_branch
!
!*******************************************************************************
!
RECURSIVE SUBROUTINE suite_branch_io(zeile, length, lreset, lloop)
!
!  Specific SUITE Version of a branch subroutine
!  Call a section via a branch
!
USE diffev_setup_mod
USE diffev_setup_sub_mod
USE diffev_loop_mod
!
USE discus_setup_mod
USE discus_setup_sub_mod
USE discus_loop_mod
USE discus_reset_all_mod
!
USE kuplot_setup_mod
USE kuplot_setup_sub_mod
USE kuplot_loop_mod
use kuplot_reset_mod
!
USE refine_setup_mod
USE refine_setup_sub_mod
USE refine_loop_mod
!
USE blanks_mod
USE appl_env_mod
USE errlist_mod
USE prompt_mod
USE suite_init_mod
USE suite_setup_mod
USE suite_set_sub_mod
USE charact_mod
USE lib_macro_func
USE precision_mod
USE str_comp_mod
USE variable_mod
!
IMPLICIT NONE
save
!
CHARACTER (LEN=*), INTENT(IN) :: zeile
INTEGER          , INTENT(IN) :: length
LOGICAL          , INTENT(IN) :: lreset
integer          , INTENT(IN) :: lloop 
!
CHARACTER(LEN=MAX(PREC_STRING,LEN(zeile)))  :: line
CHARACTER(LEN= 7   ) :: br_pname_old,br_pname_cap_old
INTEGER              :: br_prompt_status_old
INTEGER              :: br_ier_sta_old, br_state_old, br_progr_old
INTEGER              :: indxt, indxb, indxm
INTEGER              :: lbef, laenge
INTEGER              :: lcomm
LOGICAL              :: lmacro
!
if(lloop>=0) then
!
IF(str_comp(zeile, pname, LEN_TRIM(pname), length, LEN_TRIM(pname))) THEN
   ier_num = -15         ! branch to identical section
   ier_typ = ER_COMM
   ier_msg(1) = 'A branch from section '// pname(1:LEN_TRIM(pname)) // &
                ' to ' // pname(1:LEN_TRIM(pname))
   RETURN
ENDIF
!
!
! Get parameters on command line
!
lmacro = .false.
indxt  = INDEX (zeile, tab)       ! find a tabulator
IF(indxt==0) indxt = length + 1
indxb  = index (zeile, ' ') 
IF(indxb==0) indxb = length + 1
indxb  = MIN(indxb,indxt)         ! find first "white" character
lbef   = min (indxb - 1, 9) 
!
line  = ' ' 
lcomm = 0 
IF (indxb + 1.le.length) THEN     ! There is space for parameters
   line  = zeile(indxb + 1:length) 
   lcomm = length - indxb 
   CALL rem_leading_bl(line , lcomm)
ENDIF 
!
indxm = INDEX(zeile,'-macro')
IF(indxm > 0) THEN                   ! Found '-macro' qualifier
   lmacro = .true.
   indxt  = INDEX (line, tab)       ! find a tabulator
   IF(indxt==0) indxt = lcomm + 1
   indxb  = index (line, ' ') 
   IF(indxb==0) indxb = lcomm + 1
   indxb  = MIN(indxb,indxt)
   lbef   = min (indxb - 1, 9) 
   IF (indxb + 1.le.lcomm) THEN      ! There is space for a macro name
      line = line(indxb+1:lcomm)     ! copy macro name and parameter list
      CALL rem_leading_bl(line , lcomm)
      lcomm = LEN_TRIM(line)
   ENDIF
ENDIF
!
br_pname_old         = pname         ! Store old prompt/error status
br_pname_cap_old     = pname_cap
br_prompt_status_old = prompt_status
br_ier_sta_old       = ier_sta
br_state_old         = NINT(var_val(VAR_STATE))
br_progr_old         = NINT(var_val(VAR_PROGRAM))
!
IF(str_comp(zeile, 'kuplot', 2, length, 6)) THEN
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
   CALL kuplot_set_sub ()
   CALL suite_set_sub_branch ()
   IF(lreset) THEN
      line = 'all'
      laenge = 3
      call kuplot_do_reset(line, length)
   ENDIF
   var_val(VAR_PROGRAM) = var_val(VAR_KUPLOT)
   var_val(VAR_STATE)   = var_val(VAR_IS_BRANCH)
   IF(lmacro) THEN           ! Execute "command line macro"
      line = line(1:lcomm)
      CALL file_kdo(line,lcomm)
   ENDIF
   IF(lloop==0 .and. ier_num == 0) THEN     ! If no error in macro do interactive session
       CALL kuplot_loop    ()
   ENDIF
ELSEIF(str_comp(zeile, 'discus', 3, length, 6)) THEN
   IF(suite_discus_init) then
      pname     = 'discus'
      pname_cap = 'DISCUS'
      prompt    = pname
      oprompt   = pname
      CALL program_files ()
   ELSE
      CALL discus_setup   (lstandalone)
      suite_discus_init = .TRUE.
   ENDIF
   CALL discus_set_sub ()
   CALL suite_set_sub_branch ()
   IF(lreset) CALL discus_reset_all
   var_val(VAR_PROGRAM) = var_val(VAR_DISCUS)
   var_val(VAR_STATE)   = var_val(VAR_IS_BRANCH)
   IF(lmacro) THEN           ! Execute "command line macro"
      CALL file_kdo(line(1:lcomm),lcomm)
   ENDIF
   IF(lloop==0 .and. ier_num == 0) THEN     ! If no error in macro do interactive session
       CALL discus_loop    ()
   ENDIF
ELSEIF(str_comp(zeile, 'diffev', 3, length, 6)) THEN
   IF(suite_diffev_init) then
      pname     = 'diffev'
      pname_cap = 'DIFFEV'
      prompt    = pname
      oprompt   = pname
      CALL program_files ()
   ELSE
      CALL diffev_setup   (lstandalone)
      suite_diffev_init = .TRUE.
   ENDIF
   CALL diffev_set_sub ()
   CALL suite_set_sub_branch ()
   var_val(VAR_PROGRAM) = var_val(VAR_DIFFEV)
   var_val(VAR_STATE)   = var_val(VAR_IS_BRANCH)
   IF(lmacro) THEN           ! Execute "command line macro"
      CALL file_kdo(line(1:lcomm),lcomm)
   ENDIF
   IF(lloop==0 .and. ier_num == 0) THEN     ! If no error in macro do interactive session
       CALL diffev_loop    ()
   ENDIF
ELSEIF(str_comp(zeile, 'refine', 3, length, 6)) THEN
   IF(suite_diffev_init) then
      pname     = 'refine'
      pname_cap = 'REFINE'
      prompt    = pname
      oprompt   = pname
      CALL program_files ()
   ELSE
      CALL refine_setup   (lstandalone)
      suite_refine_init = .TRUE.
   ENDIF
   CALL refine_set_sub ()
   CALL suite_set_sub_branch ()
   var_val(VAR_PROGRAM) = var_val(VAR_REFINE)
   var_val(VAR_STATE)   = var_val(VAR_IS_BRANCH)
   IF(lmacro) THEN           ! Execute "command line macro"
      CALL file_kdo(line(1:lcomm),lcomm)
   ENDIF
   IF(lloop==0 .and. ier_num == 0) THEN     ! If no error in macro do interactive session
       CALL refine_loop    ()
   ENDIF
ELSE
   ier_num = -7
   ier_typ = ER_COMM
ENDIF
endif
!
if(lloop<=0) then
   pname      = br_pname_old
   pname_cap  = br_pname_cap_old
   prompt     = pname
   oprompt    = pname
   prompt_status = br_prompt_status_old
   ier_sta       = br_ier_sta_old
   var_val(VAR_PROGRAM) = br_progr_old
   var_val(VAR_STATE)   = br_state_old
!
IF(pname=='discus') THEN      ! Return to DISCUS branch
      CALL discus_set_sub ()
ELSEIF(pname=='diffev') THEN  ! Return to DIFFEV branch
      CALL diffev_set_sub ()
ELSEIF(pname=='kuplot') THEN  ! Return to KUPLOT branch
      CALL kuplot_set_sub ()
ELSEIF(pname=='refine') THEN  ! Return to REFINE branch
      CALL refine_set_sub ()
ENDIF
CALL suite_set_sub_branch ()
IF(ier_num == -9 .AND. ier_typ == 1) THEN
   CALL program_files ()
   ier_num = -9
   ier_typ = ER_COMM
ELSE
   CALL program_files ()
ENDIF
endif
!
END SUBROUTINE suite_branch_io
