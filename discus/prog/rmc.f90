MODULE rmc_menu
!
CONTAINS
!*****7*****************************************************************
!                                                                       
      SUBROUTINE rmc 
!-                                                                      
!     This sublevel includes all commands and functions for the         
!     Reverse-Monte-Carlo simulations in DISCUS.                        
!+                                                                      
      USE discus_config_mod 
      USE discus_allocate_appl_mod
      USE crystal_mod 
      USE modify_mod
      USE molecule_mod 
      USE rmc_mod 
      USE rmc_sup_mod
      USE random_mod
!
      USE calc_expr_mod
      USE doact_mod 
      USE do_eval_mod
      USE do_wait_mod
      USE errlist_mod 
      USE get_params_mod
      USE learn_mod 
USE lib_do_operating_mod
USE lib_echo
USE lib_errlist_func
USE lib_help
USE lib_length
USE lib_macro_func
      USE class_macro_internal
      USE param_mod 
USE precision_mod
      USE prompt_mod 
USE str_comp_mod
      USE string_convert_mod
      USE sup_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER, PARAMETER :: MAXW = 2
!                                                                       
      CHARACTER(len=8) befehl 
      CHARACTER(LEN=LEN(prompt)) :: orig_prompt
      CHARACTER(40) cdummy 
      CHARACTER(LEN=PREC_STRING) :: line, zeile
      CHARACTER(LEN=PREC_STRING) :: cpara(MAXW) !  (MAXSCAT) 
      INTEGER             :: lpara(MAXW) ! (MAXSCAT)
      INTEGER lp, length 
      INTEGER indxg, ianz, lbef 
!                                                                       
!
!     maxw = MAXSCAT
!
      orig_prompt = prompt
      prompt = prompt (1:len_str (prompt) ) //'/rmc' 
!                                                                       
   10 CONTINUE 
!                                                                       
      CALL no_error 
      CALL get_cmd (line, length, befehl, lbef, zeile, lp, prompt) 
      IF (ier_num.eq.0) then 
         IF (line (1:1)  == ' '.or.line (1:1)  == '#' .or.   & 
             line == char(13) .or. line(1:1) == '!'  ) GOTO 10
!                                                                       
!------ search for "="                                                  
!                                                                       
indxg = index (line, '=') 
IF (indxg.ne.0.AND..NOT. (str_comp (befehl, 'echo',   2, lbef, 4) ) &
              .AND..NOT. (str_comp (befehl, 'system', 2, lbef, 6) )    &
              .AND..NOT. (str_comp (befehl, 'help',   2, lbef, 4) .OR. &
                          str_comp (befehl, '?   ',   2, lbef, 4) )    &
              .AND. INDEX(line,'==') == 0                            ) THEN
!
            CALL do_math (line, indxg, length) 
!                                                                       
!------ execute a macro file                                            
!                                                                       
         ELSEIF (befehl (1:1) .eq.'@') then 
            line(1:length-1) = line(2:length)
            line(length:length) = ' '
            length = length - 1
            CALL file_kdo(line, length)
!                                                                       
!------ continues a macro 'continue'                                    
!                                                                       
         ELSEIF (str_comp (befehl, 'continue', 3, lbef, 8) ) then 
            CALL macro_continue (zeile, lp) 
!                                                                       
!------ Echo a string, just for interactive check in a macro 'echo'     
!                                                                       
         ELSEIF (str_comp (befehl, 'echo', 2, lbef, 4) ) then 
            CALL echo (zeile, lp) 
!                                                                       
!------ Evaluate an expression                                          
!                                                                       
         ELSEIF (str_comp (befehl, 'evaluate', 2, lbef, 8) ) then 
            CALL do_eval (zeile, lp, .TRUE.) 
!                                                                       
!     exit 'exit'                                                       
!                                                                       
         ELSEIF (str_comp (befehl, 'exit', 3, lbef, 4) ) then 
            GOTO 9999 
!                                                                       
!     help 'help','?'                                                   
!                                                                       
      ELSEIF (str_comp (befehl, 'help', 2, lbef, 4) .or.str_comp (befehl&
     &, '?   ', 1, lbef, 4) ) then                                      
            IF (str_comp (zeile, 'errors', 2, lp, 6) ) then 
               lp = lp + 7 
               CALL do_hel ('discus '//zeile, lp) 
            ELSE 
               lp = lp + 12 
               CALL do_hel ('discus rmc '//zeile, lp) 
            ENDIF 
!                                                                       
!-------Operating System Kommandos 'syst'                               
!                                                                       
         ELSEIF (str_comp (befehl, 'system', 2, lbef, 6) ) then 
            cdummy = ' ' 
            IF (zeile.ne.' ') then 
               cdummy (1:lp) = zeile (1:lp) 
               CALL do_operating (cdummy, lp) 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
!                                                                       
!     Waiting for user input                                            
!                                                                       
         ELSEIF (str_comp (befehl, 'wait', 3, lbef, 4) ) then 
            CALL do_input (zeile, lp) 
         ELSE
!
!------- RMC specific commands
!
      IF( cr_nscat > RMC_MAXSCAT .OR. cr_ncatoms> RMC_MAXSITE) THEN
         CALL alloc_rmc ( cr_nscat, cr_ncatoms )
         IF ( ier_num < 0 ) THEN
            RETURN
         ENDIF
      ENDIF
!                                                                       
!------ reset rmc setting 'rese'                                        
!                                                                       
             IF (str_comp (befehl, 'reset', 3, lbef, 5) ) then 
            CALL rmc_reset
!                                                                       
!------ selecting/deselecting atoms                                     
!                                                                       
         ELSEIF (str_comp (befehl, 'select',   3, lbef, 6) .or. &
                 str_comp (befehl, 'deselect', 2, lbef, 8) ) then
!                                                                       
            CALL atom_select (zeile, lp, 0, RMC_MAXSCAT, rmc_allowed, &
            rmc_lsite, 0, RMC_MAXSITE,                                &
            rmc_sel_atom, .false., str_comp (  &
            befehl, 'select', 3, lbef, 6) )                               
!                                                                       
!------ selecting/deselecting of molecules                              
!                                                                       
         ELSEIF (str_comp (befehl, 'mselect', 2, lbef,   7) .or.  &
                 str_comp (befehl, 'mdeselect', 2, lbef, 8) ) then                             
!                                                                       
            CALL mole_select (zeile, lp, 0, RMC_MAXSCAT, rmc_allowed, &
            rmc_sel_atom, str_comp (  &
            befehl, 'mselect', 2, lbef, 7) )                               
!                                                                       
!------ read experimental data 'data'                                   
!                                                                       
         ELSEIF (str_comp (befehl, 'data', 3, lbef, 4) ) then 
            CALL rmc_readdata (zeile, lp) 
!                                                                       
!------ run 'run'                                                       
!                                                                       
         ELSEIF (str_comp (befehl, 'run', 2, lbef, 3) ) then 
            CALL rmc_run 
!                                                                       
!-----      save structure to file 'save'                               
!                                                                       
         ELSEIF (str_comp (befehl, 'save', 2, lbef, 4) ) then 
            CALL rmc_save (zeile, lp) 
!                                                                       
!------ set rmc parameters 'set'                                        
!                                                                       
         ELSEIF (str_comp (befehl, 'set', 3, lbef, 3) ) then 
            CALL rmc_set (zeile, lp) 
!                                                                       
!------ show rmc parameters 'show'                                      
!                                                                       
         ELSEIF (str_comp (befehl, 'show', 3, lbef, 4) ) then 
            CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
            IF (ier_num.eq.0) then 
               IF (ianz.eq.0) then 
                  CALL rmc_show ('ALL') 
               ELSEIF (ianz.eq.1) then 
                  CALL do_cap (cpara (1) ) 
                  IF (cpara (1) (1:2) .eq.'RV') then 
                     CALL rmc_rvalue 
                  ELSE 
                     CALL rmc_show (cpara (1) ) 
                  ENDIF 
               ELSE 
                  ier_typ = ER_COMM 
                  ier_num = - 6 
               ENDIF 
            ENDIF 
!                                                                       
!------ no command found                                                
!                                                                       
         ELSE 
            ier_num = - 8 
            ier_typ = ER_COMM 
         ENDIF 
         ENDIF 
      ENDIF 
!                                                                       
!------ any errors ?                                                    
!                                                                       
      IF (ier_num.ne.0) THEN 
         CALL errlist 
         IF (ier_sta.ne.ER_S_LIVE) THEN 
            IF (lmakro .OR. lmakro_error) THEN  ! Error within macro or termination errror
               IF(sprompt /= prompt ) THEN
                  ier_num = -10
                  ier_typ = ER_COMM
                  ier_msg(1) = ' Error occured in rmc menu'
                  prompt_status = PROMPT_ON 
                  prompt = orig_prompt
                  RETURN
               ELSE
                  IF(lmacro_close) THEN
                     CALL macro_close 
                     prompt_status = PROMPT_ON 
                  ENDIF 
               ENDIF 
            ENDIF 
            IF (lblock) THEN 
               ier_num = - 11 
               ier_typ = ER_COMM 
               prompt_status = PROMPT_ON 
               prompt = orig_prompt
               RETURN 
            ENDIF 
            CALL no_error 
            lmakro_error = .FALSE.
            sprompt = ' '
         ENDIF 
      ENDIF 
      GOTO 10 
!                                                                       
 9999 CONTINUE 
!
      prompt = orig_prompt
!                                                                       
      END SUBROUTINE rmc                            
!
!*****7*****************************************************************
!
SUBROUTINE rmc_reset
!
USE discus_allocate_appl_mod
USE rmc_mod
!
IMPLICIT NONE
!
CALL alloc_rmc      ( 1,  1        )
CALL alloc_rmc_data ( 1            )
CALL alloc_rmc_istl ( 1,  1, 1     )
CALL alloc_rmc_q    ( 1,  1        )
CALL alloc_rmc_planes(1, 48        )
!
rmc_allowed= .FALSE.
rmc_lsite  = .TRUE.
rmc_nplane = 0 
rmc_calc_f = .true. 
rmc_qmin   =  9999.0 
rmc_qmax   = - 9999.0 
offq       = 0   ! offq (*)
offsq      = 0   ! offsq(*)
!                                                                       
!           DO i = 1, RMC_MAX_PLANES 
!           offq (i) = 0 
!           DO j = 1, RMC_MAX_SYM 
!           offsq (i, j) = 0 
!           ENDDO 
!           ENDDO 
END SUBROUTINE rmc_reset
!
!*****7*****************************************************************
!
END MODULE rmc_menu
