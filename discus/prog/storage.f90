MODULE storage_menu_mod
!
!USE crystal_mod
!
PRIVATE
PUBLIC storage
!
CONTAINS
!
SUBROUTINE storage
!
!   Main menu to administrate internal storage of structures
!
USE discus_allocate_appl_mod
!
USE calc_expr_mod
USE doact_mod
USE do_eval_mod
USE do_wait_mod
USE errlist_mod
USE class_macro_internal
USE lib_do_operating_mod
USE lib_echo
USE lib_errlist_func
USE lib_help
USE lib_length
USE lib_macro_func
USE precision_mod
USE prompt_mod
USE str_comp_mod
USE sup_mod
!
IMPLICIT none
!
CHARACTER (LEN=8)                       :: befehl! command on input line
CHARACTER(LEN=LEN(prompt))              :: orig_prompt  ! original prompt
CHARACTER (LEN=PREC_STRING)             :: line  ! input line
CHARACTER (LEN=PREC_STRING)             :: zeile ! remainder with parameters
INTEGER                                 :: indxg ! location of "="
INTEGER                                 :: lp    ! lengtz of zeile
INTEGER laenge, lbef
LOGICAL                                 :: lalloc = .false. ! Need to allocate
LOGICAL                                 :: lend  ! condition of EOF
!
!
lend   = .FALSE.
lalloc = .FALSE.
!
!
orig_prompt = prompt
prompt = prompt (1:len_str (prompt) ) //'/storage'
!
main_loop: do
   CALL no_error
   CALL get_cmd (line, laenge, befehl, lbef, zeile, lp, prompt)
   no_err: IF (ier_num.eq.0) THEN
      no_com: IF (line /= ' '      .and. line(1:1) /= '#' .and.      &
          line /= char(13) .and. line(1:1) /= '!'        ) THEN
!                                                                       
!     ----search for "="                                                
!                                                                       
      indxg = index (line, '=') 
      is_math: IF (indxg.ne.0.AND..NOT. (str_comp (befehl, 'echo', 2, lbef, 4) ) &
                    .AND..NOT. (str_comp (befehl, 'system', 2, lbef, 6) )    &
                    .AND..NOT. (str_comp (befehl, 'help', 2, lbef, 4) .OR. &
                                str_comp (befehl, '?   ', 2, lbef, 4) )    &
                    .AND. INDEX(line,'==') == 0                            ) THEN
!                                                                       
! ------evaluate an expression and assign the value to a variabble      
!                                                                       
            CALL do_math (line, indxg, laenge)
         ELSE
!                                                                       
!------ ----execute a macro file                                        
!                                                                       
           is_com: IF (befehl (1:1) .eq.'@') THEN 
               IF (laenge.ge.2) THEN 
                  line(1:laenge-1) = line(2:laenge)
                  laenge = 1
                  CALL file_kdo(line, laenge)
               ELSE 
                  ier_num = - 13 
                  ier_typ = ER_MAC 
               ENDIF 
!                                                                       
!     ----continues a macro 'continue'                                  
!                                                                       
           ELSEIF (str_comp (befehl, 'continue', 3, lbef, 8) ) THEN 
               CALL macro_continue (zeile, lp) 
!                                                                       
!     ----Echo a string, just for interactive check in a macro 'echo'   
!                                                                       
           ELSEIF (str_comp (befehl, 'echo', 2, lbef, 4) ) THEN 
               CALL echo (zeile, lp) 
!                                                                       
!      ---Evaluate an expression, just for interactive check 'eval'     
!                                                                       
           ELSEIF (str_comp (befehl, 'evaluate', 2, lbef, 8) ) THEN 
               CALL do_eval (zeile, lp, .TRUE.) 
!                                                                       
!     ----exit 'exit'                                                   
!                                                                       
           ELSEIF (str_comp (befehl, 'exit', 2, lbef, 4) ) THEN 
               lend = .true. 
               EXIT main_loop
!                                                                       
!     ----help 'help' , '?'                                             
!                                                                       
           ELSEIF (str_comp (befehl, 'help', 1, lbef, 4) .or.  &
                    str_comp (befehl, '?   ', 1, lbef, 4) ) THEN                                      
               IF (str_comp (zeile, 'errors', 2, lp, 6) ) THEN 
                  lp = lp + 7 
                  CALL do_hel ('discus '//zeile, lp) 
               ELSE 
                  lp = lp + 15 
                  CALL do_hel ('discus storage '//zeile, lp) 
               ENDIF 
!                                                                       
!------- -Operating System Kommandos 'syst'                             
!                                                                       
           ELSEIF (str_comp (befehl, 'system', 2, lbef, 6) ) THEN 
               IF (zeile.ne.' '.and.zeile.ne.char (13) ) THEN 
                  CALL do_operating (zeile (1:lp), lp) 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
!                                                                       
!------  -----waiting for user input                                    
!                                                                       
           ELSEIF (str_comp (befehl, 'wait', 3, lbef, 4) ) THEN 
               CALL do_input (zeile, lp) 
!
!     ----Original storage command
!
!
           ELSEIF (str_comp (befehl, 'remove',3, lbef, 5)) THEN
              CALL storage_reset(zeile, lp)
!
           ELSEIF (str_comp (befehl, 'reset', 3, lbef, 5)) THEN
              CALL storage_reset(zeile, lp)
!
           ELSEIF (str_comp (befehl, 'show', 3, lbef, 4)) THEN
              CALL storage_show(zeile, lp)
!
           ELSE is_com
              ier_num = -8
              ier_typ = ER_COMM
!
           ENDIF is_com ! END IF BLOCK actual commands
        ENDIF is_math   ! END IF BLOCK math equation or specific command
     ENDIF no_com       ! END IF BLOCK no comment
   ENDIF no_err         ! END IF BLOCK no error reading input
!
   IF (ier_num.ne.0) THEN 
      CALL errlist 
      IF (ier_sta.ne.ER_S_LIVE) THEN 
         IF (lmakro .OR. lmakro_error) THEN  ! Error within macro or termination errror
            IF(sprompt /= prompt ) THEN
               ier_num = -10
               ier_typ = ER_COMM
               ier_msg(1) = ' Error occured in decorate menu'
               prompt_status = PROMPT_ON 
               prompt = orig_prompt
               RETURN
            ELSE
               IF(lmacro_close) THEN
                  CALL macro_close(-1)
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
!
!
ENDDO main_loop     ! END DO main loop of menu 
!
prompt = orig_prompt
!
END SUBROUTINE storage
!
!*******************************************************************************
!
SUBROUTINE storage_show(zeile, lp)
!
USE class_internal
USE get_params_mod
USE errlist_mod
USE precision_mod
USE str_comp_mod
USE take_param_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: zeile
INTEGER         , INTENT(INOUT) :: lp
!
INTEGER, PARAMETER :: MAXW = 2
LOGICAL, PARAMETER :: LALL = .TRUE.
LOGICAL, PARAMETER :: LSINGLE = .FALSE.
!
CHARACTER(LEN=MAX(PREC_STRING,LEN(zeile))), DIMENSION(1:MAXW) :: cpara
INTEGER            , DIMENSION(1:MAXW) :: lpara
INTEGER                                :: ianz
!
INTEGER, PARAMETER :: NOPTIONAL = 1
CHARACTER(LEN=7), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=MAX(PREC_STRING,LEN(zeile))), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent  !opt. para present
REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 0 ! Number of values to calculate 
!
!
DATA oname  / 'display' /
DATA loname /  7        /
opara  =  (/ 'short'    /)   ! Always provide fresh default values
lopara =  (/  5         /)
owerte =  (/  0.0       /)
!
cpara(:) = ' '
lpara(:) =  1
IF(zeile== ' ') THEN
   CALL store_list_node(LALL, 'all', opara(1))
ELSE
   CALL get_params (zeile, ianz, cpara, lpara, MAXW, lp)
   CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                     oname, loname, opara, lopara, lpresent, owerte)
!
   IF(ier_num ==0) THEN
      IF(ianz ==0) THEN
         CALL store_list_node(LALL, 'all', opara(1))
      ELSE
         IF(str_comp(cpara(1), 'all', 3, lpara(1) , 3) ) THEN
            CALL store_list_node(LALL, 'all', opara(1))
         ELSE
            CALL store_list_node(LSINGLE, cpara(1)(1:lpara(1)), opara(1))
         ENDIF
      ENDIF
   ENDIF
ENDIF
!
END SUBROUTINE storage_show
!
!*******************************************************************************
!
SUBROUTINE storage_reset(zeile, lp)
!
USE class_internal
USE get_params_mod
USE errlist_mod
USE precision_mod
USE str_comp_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: zeile
INTEGER         , INTENT(INOUT) :: lp
!
INTEGER, PARAMETER :: MAXW = 1
!
CHARACTER(LEN=MAX(PREC_STRING,LEN(zeile))), DIMENSION(1:MAXW) :: cpara
INTEGER            , DIMENSION(1:MAXW) :: lpara
INTEGER                                :: ianz
!
!
IF(zeile== ' ') THEN
   CALL store_remove_all(store_root)
ELSE
   CALL get_params (zeile, ianz, cpara, lpara, MAXW, lp)
   IF(ier_num ==0) THEN
      IF(str_comp(cpara(1), 'all', 3, lpara(1) , 3) ) THEN
         CALL store_remove_all(store_root)
      ELSE
         CALL store_remove_single(cpara(1)(1:lpara(1)), ier_num)
         IF(ier_num/=0) ier_typ = ER_APPL
      ENDIF
   ENDIF
ENDIF
!
END SUBROUTINE storage_reset
!
END MODULE storage_menu_mod
