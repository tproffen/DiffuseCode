!*****7*****************************************************************
!                                                                       
SUBROUTINE suite_mache_kdo (line, lend, length) 
!+                                                                      
!     This is the main routine for command interpretation, each         
!     command is identified here and the corresponding subroutine       
!     executed. A leading @ indicates a macro.                          
!-                                                                      
!USE discus,  discus_interactive => interactive
USE diffev_setup_mod
USE diffev_loop_mod
USE discus_setup_mod
USE discus_loop_mod
USE kuplot_setup_mod
USE kuplot_loop_mod
USE suite_init_mod
USE suite_parallel_mod
USE suite_setup_mod
!
USE charact_mod 
!USE allocate_appl
!
USE appl_env_mod
USE doact_mod
USE errlist_mod 
USE class_macro_internal
USE learn_mod 
USE prompt_mod
USE envir_mod
USE variable_mod
!USE set_sub_generic_mod
!
IMPLICIT none 
!                                                                       
CHARACTER (LEN= *  ), INTENT(INOUT) :: line 
LOGICAL             , INTENT(  OUT) :: lend 
INTEGER             , INTENT(INOUT) :: length 
!
CHARACTER (LEN=1024)                  :: zeile 
CHARACTER (LEN=   9)                  :: befehl 
INTEGER                               :: indxb, indxg, lcomm, lbef, indxt 
!                                                                       
LOGICAL, EXTERNAL                     :: str_comp 
!                                                                       
lstandalone = .false.  ! Switch to slave mode for DIFFEV/DISCUS/KUPLOT
!                                                                       
CALL no_error 
!                                                                 
!-------If a commentary return immediately                        
!                                                                 
IF (line (1:1) .EQ.' ' .or. line (1:1) .eq.'#'.or. &
    line (1:1) .eq.'!' .or. length.eq.0           ) RETURN                                     
!                                                                 
!     Only the first 5 characters are significant. The command    
!     consists of the four nonblank characters                    
!                                                                 
befehl = '    ' 
indxt  = INDEX (line, tab)       ! find a tabulator
IF(indxt==0) indxt = length + 1
indxb  = index (line, ' ') 
IF(indxb==0) indxb = length + 1
indxb  = MIN(indxb,indxt)
lbef   = min (indxb - 1, 9) 
befehl = line (1:lbef) 
!                                                                 
!------ command parameters start at the first character following 
!     the blank                                                   
!                                                                 
zeile = ' ' 
lcomm = 0 
IF (indxb + 1.le.length) then 
   zeile = line (indxb + 1:length) 
   lcomm = length - indxb 
   call rem_leading_bl(zeile, lcomm)
ENDIF 
!                                                                 
!-------Suche nach einem "="                                      
!                                                                 
indxg = index (line, '=') 
IF (indxg.ne.0.and.                                              &
         .not. (str_comp (befehl, 'echo',  2, lbef, 4) ) .and.   &
         .not. (str_comp (befehl, 'syst',  2, lbef, 4) ) .and.   &
         .not. (str_comp (befehl, 'fput',  2, lbef, 4) ) .and.   &
         .not. (str_comp (befehl, 'socket',2, lbef, 5) ) .and.   &
         .not. (str_comp (befehl, 'help',  2, lbef, 4) .or.      &
          str_comp (befehl, '?   ',  2, lbef, 4) )       .AND.   &
         INDEX(line,'==') == 0                                ) THEN      
!                                                                 
!-------Zuweisung eines Funktionswertes                           
!                                                                 
   CALL do_math (line, indxg, length) 
ELSE 
!                                                                 
!     --execute a macro file                                      
!                                                                 
   IF (befehl (1:1) .eq.'@') then 
      IF (length.ge.2) then 
         CALL file_kdo (line (2:length), length - 1) 
      ELSE 
         ier_num = - 13 
         ier_typ = ER_MAC 
      ENDIF 
!                                                                 
!-------Terminate DISCUS_SUITE 'exit'                                   
!                                                                 
   ELSEIF (str_comp (befehl, 'exit', 2, lbef, 4) ) then 
      LEND = .TRUE. 
!                                                                 
!     ----Start of DISCUS_SUITE special commands                        
!                                                                 

!                                                                 
!     -- branch to DIFFEV
!
   ELSEIF ((linteractive.OR.lblock.OR.lmakro) .AND. str_comp (befehl, 'diffev', 3, lbef, 6) ) then
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
       var_val(VAR_STATE)   = var_val(VAR_IS_SECTION)
       var_val(VAR_PROGRAM) = var_val(VAR_DIFFEV)
       CALL diffev_set_sub ()
       CALL suite_set_sub_branch
       CALL diffev_loop    ()
       lend      = .false.
       pname     = 'suite'
       pname_cap = 'SUITE'
       prompt    = pname
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
!     -- branch to DISCUS
!
   ELSEIF ((linteractive.OR.lblock.OR.lmakro) .AND. str_comp (befehl, 'discus', 3, lbef, 6) ) then
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
       var_val(VAR_STATE)   = var_val(VAR_IS_SECTION)
       var_val(VAR_PROGRAM) = var_val(VAR_DISCUS)
       CALL discus_set_sub ()
       CALL suite_set_sub_branch
       CALL discus_loop    ()
       lend = .false.
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
!     -- branch to KUPLOT
!
   ELSEIF ((linteractive.OR.lblock.OR.lmakro) .AND. str_comp (befehl, 'kuplot', 3, lbef, 6) ) then
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
       CALL kuplot_loop    ()
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
!     -- Run a parallel version of discus_suite
!
   ELSEIF (str_comp (befehl, 'parallel', 3, lbef, 8) ) then
         CALL suite_do_parallel (zeile, lcomm)
!!!      CALL do_deallocate_appl (zeile, lcomm)
!                                                                 
!------   Try general commands                                    
!                                                                 
   ELSE 
      CALL kdo_all (befehl, lbef, zeile, lcomm) 
   ENDIF 
ENDIF 
!                                                                       
END SUBROUTINE suite_mache_kdo                      
