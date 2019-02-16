MODULE kuplot_2dm
!
USE kuplot_2dm_mod
!
PRIVATE
PUBLIC kuplot_2dm_menu
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE kuplot_2dm_menu
!
USE kuplot_config
USE kuplot_mod
!
!
USE ber_params_mod
USE calc_expr_mod
USE class_macro_internal
USE do_eval_mod
USE do_wait_mod
USE doact_mod
USE errlist_mod
USE get_params_mod
USE learn_mod
USE prompt_mod
USE sup_mod
USE take_param_mod
!
IMPLICIT NONE
!
CHARACTER (LEN=5)                       :: befehl! command on input line
CHARACTER(LEN=LEN_TRIM(prompt))         :: orig_prompt  ! original prompt
CHARACTER (LEN=1024)                    :: line  ! input line
CHARACTER (LEN=1024)                    :: zeile ! remainder with parameters
INTEGER                                 :: indxg ! location of "="
INTEGER                                 :: lp    ! length of zeile
INTEGER                                 :: laenge
INTEGER                                 :: lbef
LOGICAL                                 :: lend
!
!
INTEGER, PARAMETER :: MAXW = 20
CHARACTER(LEN=1024), DIMENSION(MAXW) :: cpara
INTEGER            , DIMENSION(MAXW) :: lpara
REAL               , DIMENSION(MAXW) :: werte
INTEGER                              :: ianz
!
!INTEGER, EXTERNAL :: len_str
LOGICAL, EXTERNAL :: str_comp
!
!
orig_prompt = prompt
prompt = prompt (1:LEN_TRIM(prompt) ) //'/2dm'
!
!
main_loop: DO
   CALL no_error
   CALL get_cmd (line, laenge, befehl, lbef, zeile, lp, prompt)
   no_err: IF(ier_num.eq.0) THEN
      no_com: IF(line /= ' '      .AND. line(1:1) /= '#' .AND.      &
                 line /= char(13) .AND. line(1:1) /= '!'        ) THEN
!                                                                       
!     ----search for "="                                                
!                                                                       
         indxg = index (line, '=') 
         is_math: IF(indxg.ne.0                                             &
                     .AND..NOT. (str_comp (befehl, 'echo', 2, lbef, 4) )    &
                     .AND..NOT. (str_comp (befehl, 'syst', 2, lbef, 4) )    &
                     .AND..NOT. (str_comp (befehl, 'help', 2, lbef, 4) .OR. &
                                 str_comp (befehl, '?   ', 2, lbef, 4) )    &
                     .AND. INDEX(line,'==') == 0                            ) THEN
!                                                                       
! ------evaluate an expression and assign the value to a variabble      
!                                                                       
               CALL do_math (line, indxg, laenge)
         ELSE is_math                            ! is_math, al other commands
!                                                                       
!------ ----execute a macro file                                        
!                                                                       
            is_generic: IF (befehl (1:1) .eq.'@') THEN     ! macro, reset or all other commands
               IF (laenge.ge.2) THEN 
                  CALL file_kdo (line (2:laenge), laenge-1) 
               ELSE 
                  ier_num = - 13 
                  ier_typ = ER_MAC 
               ENDIF 
!                                                                       
!     ----continues a macro 'continue'                                  
!                                                                       
            ELSEIF (str_comp (befehl, 'continue', 3, lbef, 8) ) THEN  is_generic
               CALL macro_continue (zeile, lp) 
!                                                                       
!     ----Echo a string, just for interactive check in a macro 'echo'   
!                                                                       
            ELSEIF (str_comp (befehl, 'echo', 2, lbef, 4) ) THEN  is_generic
               CALL echo (zeile, lp) 
!                                                                       
!      ---Evaluate an expression, just for interactive check 'eval'     
!                                                                       
            ELSEIF (str_comp (befehl, 'eval', 2, lbef, 4) ) THEN 
               CALL do_eval (zeile, lp, .TRUE.) 
!                                                                       
!     ----exit 'exit'                                                   
!                                                                       
            ELSEIF (str_comp (befehl, 'exit', 2, lbef, 4) ) THEN  is_generic
               lend = .true. 
               EXIT main_loop
!                                                                       
!     ----help 'help' , '?'                                             
!                                                                       
            ELSEIF (str_comp (befehl, 'help', 1, lbef, 4) .OR.  &
                    str_comp (befehl, '?   ', 1, lbef, 4) ) THEN is_generic
               IF (str_comp (zeile, 'errors', 2, lp, 6) ) THEN 
                  lp = lp + 7 
                  CALL do_hel ('kuplot '//zeile, lp) 
               ELSE 
                  lp = lp + 12 
                  CALL do_hel ('kuplot 2d   '//zeile, lp) 
               ENDIF 
!                                                                       
!------- -Operating System Kommandos 'syst'                             
!                                                                       
            ELSEIF (str_comp (befehl, 'syst', 2, lbef, 4) ) THEN is_generic
                IF (zeile.ne.' '.and.zeile.ne.char (13) ) THEN
                   CALL do_operating (zeile (1:lp), lp) 
                ELSE 
                   ier_num = - 6 
                   ier_typ = ER_COMM 
                ENDIF 
!                                                                       
!------  -----waiting for user input                                    
!                                                                       
            ELSEIF (str_comp (befehl, 'wait', 3, lbef, 4) ) THEN  is_generic
                CALL do_input (zeile, lp) 
!
            ELSEIF (str_comp (befehl, 'reset', 3, lbef, 4)) THEN is_generic
               CALL do_rese(zeile,lp)
            ELSE is_generic    ! macro, reset or all other commands is_generic
!
               is_com: IF(str_comp (befehl, 'show', 3, lbef, 4) ) THEN 
!
                  WRITE(output_io,'(a)') ' Will show 2dm '
               ELSEIF(str_comp (befehl, 'load', 3, lbef, 4) ) THEN  is_com
                  CALL get_load_2dm(zeile, lp)
               ELSEIF(str_comp (befehl, 'numbers', 3, lbef, 7) ) THEN  is_com
                  CALL get_numbers(zeile, lp)
               ELSEIF(str_comp (befehl, 'xrange', 3, lbef, 6) ) THEN  is_com
                  CALL get_xrange(zeile, lp)
               ELSEIF(str_comp (befehl, 'run', 3, lbef, 3) ) THEN  is_com
                  CALL k2dm_run
               ELSE is_com
                  ier_num = -8
                  ier_typ = ER_COMM
               ENDIF is_com ! END IF BLOCK actual commands
!
            ENDIF is_generic
         ENDIF is_math
      ENDIF no_com
   ENDIF no_err
!
   IF (ier_num.ne.0) THEN 
      CALL errlist 
      IF (ier_sta.ne.ER_S_LIVE) THEN 
         IF (lmakro .OR. lmakro_error) THEN  ! Error within macro or termination errror
            IF(sprompt /= prompt ) THEN
               ier_num = -10
               ier_typ = ER_COMM
               ier_msg(1) = ' Error occured in demolec menu'
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
!
ENDDO main_loop
!
prompt = orig_prompt
!
END SUBROUTINE kuplot_2dm_menu
!
!*******************************************************************************
!
SUBROUTINE get_load_2dm(line, length)
!
USE errlist_mod
USE get_params_mod
USE take_param_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: line
INTEGER         , INTENT(INOUT) :: length
!
INTEGER, PARAMETER :: O_SKIP  = 1
INTEGER, PARAMETER :: O_COLX  = 2
INTEGER, PARAMETER :: O_COLY  = 3
INTEGER, PARAMETER :: O_COLDX = 4
INTEGER, PARAMETER :: O_COLDY = 5
INTEGER, PARAMETER :: O_SEP   = 6
!
INTEGER, PARAMETER :: NOPTIONAL = 6
CHARACTER(LEN=1024), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=1024), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent!opt. para present
REAL               , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 5 ! Number of values to calculate 
!
INTEGER, PARAMETER :: MAXW = NOPTIONAL + 20
CHARACTER(LEN=1024), DIMENSION(MAXW) :: cpara
INTEGER            , DIMENSION(MAXW) :: lpara
INTEGER                              :: ianz
!
LOGICAL, EXTERNAL :: str_comp
!
DATA oname  / 'skip', 'colx',  'coly',  'coldx', 'coldy', 'separator'  /
DATA loname /  4    ,  4    ,   4    ,   5     ,  5     ,  9           /
opara  =  (/ '25.000', '1.0000', '2.0000', '0.0000', '0.0000', ';     ' /)   ! Always provide fresh default values
lopara =  (/  6,        6,        6      ,  6      ,  6      ,  6       /)
owerte =  (/ 25.0,      1.0,      2.0    ,  0.0    ,  0.0    ,  0.0     /)
!
write(*,*) 'LINE ', line(1:len_trim(line))
CALL no_error 
CALL get_params (line, ianz, cpara, lpara, MAXW, length) 
IF (ier_num /= 0) RETURN 
write(*,*) ' cpara 1 ',cpara(1)(1:lpara(1))
CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
!
IF(str_comp (cpara(1), 'csv', 3, lpara(1), 3) ) THEN
   k2dm_type  = K2DM_CSV
   k2dm_line  = line
   k2dm_loop1 = INDEX(line, 'LOOP')
!ELSEIF(str_comp (cpara(1), 'xy', 2, lpara(1), 2) ) THEN
!   k2dm_type = K2DM_XY
ELSE
   ier_num = -9999
   ier_typ = ER_APPL
   ier_msg(1) = '2DM load requires CSV format'
   RETURN
ENDIF
CALL del_params (1, ianz, cpara, lpara, maxw)
!
IF(ianz >= 2) THEN
   k2dm_nanz = ianz
   IF(ALLOCATED(k2dm_cpara)) THEN
      DEALLOCATE(k2dm_cpara)
      DEALLOCATE(k2dm_lpara)
   ENDIF
   ALLOCATE(k2dm_cpara(1:k2dm_nanz))
   ALLOCATE(k2dm_lpara(1:k2dm_nanz))
   k2dm_cpara(1:k2dm_nanz) = cpara(1:k2dm_nanz)
   k2dm_lpara(1:k2dm_nanz) = lpara(1:k2dm_nanz)
ELSE
   ier_num = -9999
   ier_typ = ER_APPL
ENDIF 
!
END SUBROUTINE get_load_2dm
!
!*******************************************************************************
!
SUBROUTINE get_numbers(line, length)
!
USE errlist_mod
USE ber_params_mod
USE get_params_mod
USE take_param_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: line
INTEGER         , INTENT(INOUT) :: length
!
INTEGER, PARAMETER :: O_START = 1
INTEGER, PARAMETER :: O_END   = 2
INTEGER, PARAMETER :: O_STEP  = 3
!
INTEGER, PARAMETER :: NOPTIONAL = 3
CHARACTER(LEN=1024), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=1024), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent!opt. para present
REAL               , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 3 ! Number of values to calculate 
!
INTEGER, PARAMETER :: MAXW = MIN(20,NOPTIONAL)
CHARACTER(LEN=1024), DIMENSION(MAXW) :: cpara
INTEGER            , DIMENSION(MAXW) :: lpara
REAL               , DIMENSION(MAXW) :: werte
INTEGER                              :: ianz
!
LOGICAL, EXTERNAL :: str_comp
!
DATA oname  / 'start', 'end', 'step' /
DATA loname /  4     ,  3   ,  4     /
opara  =  (/ '1'     , '1'  , '1'    /)   ! Always provide fresh default values
lopara =  (/  1      ,  1   ,  1     /)
owerte =  (/  1.0    ,  1.0 ,  1.0   /)
!
CALL no_error 
CALL get_params (line, ianz, cpara, lpara, MAXW, length) 
IF (ier_num /= 0) RETURN 
CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
!
k2dm_start = owerte(O_START)
k2dm_end   = owerte(O_END)
k2dm_step  = owerte(O_STEP)
!
END SUBROUTINE get_numbers
!
!*******************************************************************************
!
SUBROUTINE get_xrange(line, length)
!
USE errlist_mod
USE ber_params_mod
USE get_params_mod
USE take_param_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: line
INTEGER         , INTENT(INOUT) :: length
!
INTEGER, PARAMETER :: O_XMIN = 1
INTEGER, PARAMETER :: O_XMAX = 2
!
INTEGER, PARAMETER :: NOPTIONAL = 2
CHARACTER(LEN=1024), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=1024), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent!opt. para present
REAL               , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 0 ! Number of values to calculate 
!
INTEGER, PARAMETER :: MAXW = MIN(20,NOPTIONAL)
CHARACTER(LEN=1024), DIMENSION(MAXW) :: cpara
INTEGER            , DIMENSION(MAXW) :: lpara
REAL               , DIMENSION(MAXW) :: werte
INTEGER                              :: ianz
!
LOGICAL, EXTERNAL :: str_comp
!
DATA oname  / 'xmin', 'xmax' /
DATA loname /  4    ,  4     /
opara  =  (/ 'xmin',   'xmax' /)   ! Always provide fresh default values
lopara =  (/  4,        4     /)
owerte =  (/  0.0,      0.0   /)
!
CALL no_error 
CALL get_params (line, ianz, cpara, lpara, MAXW, length) 
IF (ier_num /= 0) RETURN 
CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
!
IF(opara(O_XMIN)=='xmin') THEN
   k2dm_lxmin = .TRUE.
ELSE
   k2dm_lxmin = .FALSE.
   cpara(1) = opara(O_XMIN)
   lpara(1) = lopara(O_XMIN)
   ianz = 1
   CALL ber_params (ianz, cpara, lpara, werte, MAXW) 
   IF(ier_num/=0) THEN
      k2dm_xmin = werte(1)
   ENDIF
ENDIF
IF(opara(O_XMAX)=='xmax') THEN
   k2dm_lxmax = .TRUE.
ELSE
   k2dm_lxmax = .FALSE.
   cpara(1) = opara(O_XMAX)
   lpara(1) = lopara(O_XMAX)
   ianz = 1
   CALL ber_params (ianz, cpara, lpara, werte, MAXW) 
   IF(ier_num/=0) THEN
      k2dm_xmax = werte(1)
   ENDIF
ENDIF
!
END SUBROUTINE get_xrange
!
!*******************************************************************************
!
SUBROUTINE k2dm_run
!
IMPLICIT NONE
!
CHARACTER(LEN=1024) :: string
CHARACTER(LEN=10  ) :: cnumber
INTEGER :: i
INTEGER :: lp
!
DO i=k2dm_start, k2dm_end, k2dm_step
   WRITE(cnumber,'(i10)') i
   string = k2dm_line(1:k2dm_loop1-1) // cnumber // k2dm_line(k2dm_loop1+4:LEN_TRIM(k2dm_line))
   lp = LEN_TRIM(string)
write(*,*) ' STRING : ', string(1:len_trim(string)) 
   CALL do_load(string,lp, .TRUE.)
ENDDO
!
END SUBROUTINE k2dm_run
!
!*******************************************************************************
!
END MODULE kuplot_2dm
