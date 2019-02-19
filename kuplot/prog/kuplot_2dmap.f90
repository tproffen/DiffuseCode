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
USE calc_expr_mod
USE class_macro_internal
USE do_eval_mod
USE do_wait_mod
USE doact_mod
USE errlist_mod
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
               ELSEIF(str_comp (befehl, 'back', 3, lbef, 4) ) THEN  is_com
                  CALL get_back_2dm(zeile, lp)
               ELSEIF(str_comp (befehl, 'load', 3, lbef, 4) ) THEN  is_com
                  CALL get_load_2dm(zeile, lp)
               ELSEIF(str_comp (befehl, 'loop', 3, lbef, 4) ) THEN is_com
                  CALL get_numbers(zeile, lp)
               ELSEIF(str_comp (befehl, 'xrange', 3, lbef, 6) ) THEN  is_com
                  CALL get_xrange(zeile, lp)
               ELSEIF(str_comp (befehl, 'yfunc', 3, lbef, 5) ) THEN  is_com
                  CALL get_yfunc(zeile, lp)
               ELSEIF(str_comp (befehl, 'run', 3, lbef, 3) ) THEN  is_com
                  CALL k2dm_run
               ELSEIF(str_comp (befehl, 'reset', 3, lbef, 5) ) THEN  is_com
                  CALL k2dm_reset
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
               ier_msg(1) = 'Error occured in 2dm menu'
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
CALL no_error 
CALL get_params (line, ianz, cpara, lpara, MAXW, length) 
IF (ier_num /= 0) RETURN 
CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
!
IF(str_comp (cpara(1), 'csv', 3, lpara(1), 3) ) THEN
   k2dm_type  = K2DM_CSV
   k2dm_line  = line
!   k2dm_loop1 = INDEX(line, 'LOOP')
!ELSEIF(str_comp (cpara(1), 'xy', 2, lpara(1), 2) ) THEN
!   k2dm_type = K2DM_XY
ELSE
   ier_num = -6
   ier_typ = ER_FORT
   ier_msg(1) = '2DM load requires CSV format'
   RETURN
ENDIF
CALL del_params (1, ianz, cpara, lpara, maxw)
!
END SUBROUTINE get_load_2dm
!
!*******************************************************************************
!
SUBROUTINE get_back_2dm(line, length)
!
USE errlist_mod
USE build_name_mod
USE get_params_mod
USE take_param_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: line
INTEGER         , INTENT(INOUT) :: length
!
!
INTEGER, PARAMETER :: O_SKIP  = 1
INTEGER, PARAMETER :: O_COLX  = 2
INTEGER, PARAMETER :: O_COLY  = 3
INTEGER, PARAMETER :: O_COLDX = 4
INTEGER, PARAMETER :: O_COLDY = 5
INTEGER, PARAMETER :: O_SCALE = 6
INTEGER, PARAMETER :: O_SEP   = 7
!
INTEGER, PARAMETER :: NOPTIONAL = 7
CHARACTER(LEN=1024), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=1024), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent!opt. para present
REAL               , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 6 ! Number of values to calculate 
!
INTEGER, PARAMETER :: MAXW = MAX(20,NOPTIONAL)
CHARACTER(LEN=1024), DIMENSION(MAXW) :: cpara
INTEGER            , DIMENSION(MAXW) :: lpara
REAL               , DIMENSION(MAXW) :: werte
INTEGER                              :: ianz
INTEGER                              :: i
!
LOGICAL, EXTERNAL :: str_comp
!
DATA oname  / 'skip', 'colx',  'coly',  'coldx', 'coldy', 'scale', 'separator'  /
DATA loname /  4    ,  4    ,   4    ,   5     ,  5     ,  5     ,  9           /
opara  =  (/ '25.000', '1.0000', '2.0000', '0.0000', '0.0000', '1.0000', ';     ' /)   ! Always provide fresh default values
lopara =  (/  6,        6,        6      ,  6      ,  6      ,  6      ,  6       /)
owerte =  (/ 25.0,      1.0,      2.0    ,  0.0    ,  0.0    ,  1.0    ,  0.0     /)
!
CALL no_error 
CALL get_params (line, ianz, cpara, lpara, MAXW, length) 
write(*,*) ' ier ', ier_num, ier_typ
IF (ier_num /= 0) RETURN 
CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
write(*,*) ' ier ', ier_num, ier_typ
!
k2dm_scale = owerte(O_SCALE)
!
IF(str_comp (cpara(1), 'off', 3, lpara(1), 3) ) THEN
   k2dm_line_b = ' '
   RETURN
ELSEIF(str_comp (cpara(1), 'csv', 3, lpara(1), 3) ) THEN
   k2dm_type   = K2DM_CSV
   k2dm_scale  = owerte(O_SCALE)
!ELSEIF(str_comp (cpara(1), 'xy', 2, lpara(1), 2) ) THEN
!   k2dm_type = K2DM_XY
ELSE
   ier_num = -6
   ier_typ = ER_FORT
   ier_msg(1) = '2DM load requires CSV format'
   RETURN
ENDIF
!
k2dm_line_b = cpara(1)(1:lpara(1))
DO i=2,ianz
   IF(INDEX(cpara(i)(1:lpara(i)),'scale:')==0) THEN
      k2dm_line_b(LEN_TRIM(k2dm_line_b)+1:) = ','//cpara(i)(1:lpara(i))
   ENDIF
ENDDO
opts: DO i=1,NOPTIONAL
   IF(i==O_SCALE) CYCLE opts
   k2dm_line_b(LEN_TRIM(k2dm_line_b)+1:) = ','//oname(i)(1:loname(i))//':'//opara(i)(1:lopara(i))
ENDDO opts
write(*,*) ' ier ', ier_num, ier_typ
!
END SUBROUTINE get_back_2dm
!
!*******************************************************************************
!
SUBROUTINE get_numbers(line, length)
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
INTEGER, PARAMETER :: O_START = 1
INTEGER, PARAMETER :: O_END   = 2
INTEGER, PARAMETER :: O_STEP  = 3
INTEGER, PARAMETER :: O_MISS  = 4
!
INTEGER, PARAMETER :: NOPTIONAL = 4
CHARACTER(LEN=1024), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=1024), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent!opt. para present
REAL               , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 3 ! Number of values to calculate 
!
INTEGER, PARAMETER :: MAXW = MAX(20,NOPTIONAL)
CHARACTER(LEN=1024), DIMENSION(MAXW) :: cpara
INTEGER            , DIMENSION(MAXW) :: lpara
REAL               , DIMENSION(MAXW) :: werte
INTEGER                              :: ianz
!
LOGICAL, EXTERNAL :: str_comp
!
DATA oname  / 'start', 'end', 'step', 'miss' /
DATA loname /  4     ,  3   ,  4    ,  4     /
opara  =  (/ '1    ', '1    ', '1    ', 'error' /)   ! Always provide fresh default values
lopara =  (/  5     ,  5     ,  5     ,  5      /)
owerte =  (/  1.0   ,  1.0   ,  1.0   ,  0.0    /)
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
IF(opara(O_MISS) == 'error') THEN
   k2dm_miss = K2DM_ERROR
ELSEIF(opara(O_MISS) == 'ignore') THEN
   k2dm_miss = K2DM_IGNORE
ELSEIF(opara(O_MISS) == 'blank') THEN
   k2dm_miss = K2DM_BLANK
ELSE
   ier_num = -6
   ier_typ = ER_FORT
   ier_msg(1) = 'The optional parameter ''miss:'' must be'
   ier_msg(2) = '''error'', ''ignore'', or ''blank'' '
ENDIF
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
INTEGER, PARAMETER :: MAXW = MAX(20,NOPTIONAL)
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
   IF(ier_num==0) THEN
      k2dm_xmin = werte(1)
   ELSE
      RETURN
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
   IF(ier_num==0) THEN
      k2dm_xmax = werte(1)
   ELSE
      RETURN
   ENDIF
ENDIF
!
END SUBROUTINE get_xrange
!
!*******************************************************************************
!
SUBROUTINE get_yfunc (line, length)
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
INTEGER, PARAMETER :: O_LOOP = 1
INTEGER, PARAMETER :: O_MISS = 2
!
INTEGER, PARAMETER :: NOPTIONAL = 1
CHARACTER(LEN=1024), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=1024), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent!opt. para present
REAL               , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 0 ! Number of values to calculate 
!
INTEGER, PARAMETER :: MAXW = MAX(20,NOPTIONAL)
CHARACTER(LEN=1024), DIMENSION(MAXW) :: cpara
INTEGER            , DIMENSION(MAXW) :: lpara
REAL               , DIMENSION(MAXW) :: werte
INTEGER                              :: ianz
!
LOGICAL, EXTERNAL :: str_comp
!
DATA oname  / 'loop'/
DATA loname /  4    /
opara  =  (/ 'LOOP' /)   ! Always provide fresh default values
lopara =  (/  4     /)
owerte =  (/  0.0   /)
!
CALL no_error 
CALL get_params (line, ianz, cpara, lpara, MAXW, length) 
IF (ier_num /= 0) RETURN 
CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
!
IF(opara(O_LOOP)=='LOOP') THEN
   k2dm_type_yf = K2DM_YF_LOOP
   k2dm_line_yf = 'LOOP'
ELSEIF(opara(O_LOOP)=='inc') THEN
   k2dm_type_yf = K2DM_YF_INC
   k2dm_line_yf = 'inc'
ELSE
   k2dm_type_yf = K2DM_YF_FUNC
   k2dm_line_yf = opara(O_LOOP)
   continue
ENDIF
!
END SUBROUTINE get_yfunc 
!
!*******************************************************************************
!
SUBROUTINE k2dm_run
!
USE kuplot_mod
!
USE ber_params_mod
USE errlist_mod
USE get_params_mod
USE variable_mod
!
IMPLICIT NONE
!
LOGICAL, PARAMETER :: l_no_echo = .FALSE.
CHARACTER(LEN=1024) :: string
CHARACTER(LEN=10  ) :: cnumber
INTEGER :: length
INTEGER :: i, ik
INTEGER :: ikm                  ! Data set with map
INTEGER :: ikb                  ! Data set with background
INTEGER :: ix, iy
INTEGER :: iix                  ! Running index for map x
INTEGER :: lp
INTEGER :: imin, imax
LOGICAL :: lblank               ! IO error and k2dm_miss=blank
REAL    :: xxmin, xxmax
REAL    :: yy                   ! background corrected y value of a single 1D
REAL    :: DELTA
!
INTEGER, PARAMETER :: MAXW = 20
CHARACTER(LEN=1024), DIMENSION(MAXW) :: cpara
INTEGER            , DIMENSION(MAXW) :: lpara
REAL               , DIMENSION(MAXW) :: werte
INTEGER                              :: ianz
!
IF(k2dm_type_yf==K2DM_YF_FUNC) THEN                       ! yfunc is a function get params
   IF(k2dm_line_yf(1:1) == '''') k2dm_line_yf(1:1) = ' '
   i = LEN_TRIM(k2dm_line_yf)
   IF(k2dm_line_yf(i:i) == '''') k2dm_line_yf(i:i) = ' '
   length = LEN_TRIM(k2dm_line_yf)
   CALL get_params (k2dm_line_yf, ianz, cpara, lpara, MAXW, length) 
   IF(ier_num /= 0) THEN
      ier_msg(1) = 'Error getting prameters of yfunction'
      RETURN
   ENDIF 
ENDIF
!
ikm = 0
ikb = 0
!
var_val(VAR_LOOP) = k2dm_start
string = k2dm_line
lp = LEN_TRIM(string)
CALL do_load(string,lp, l_no_echo)
IF(ier_num/=0) THEN
   ier_msg(1) = 'Error loading first data set '
   RETURN
ENDIF
ik   = iz - 1                 ! This is the number of the current data set
imin = 1                      ! Default to first data point
imax = len(ik)                ! Default to length of data sets
IF(k2dm_lxmin) THEN           ! user specified xmin:xmin
   xxmin = x(offxy(ik-1)+1)
   imin  = 1
ELSE
   xxmin = k2dm_xmin
ENDIF
IF(k2dm_lxmax) THEN           ! user specified xmin:xmin
   xxmax = x(offxy(ik-1)+len(ik))
   imin  = 1
ELSE
   xxmax = k2dm_xmax
ENDIF
DELTA = ABS(x(offxy(ik-1)+2)-x(offxy(ik-1)+1))*0.01    ! Calculate a sigma as (x(2)-x(1))/10
minmax:DO i=1,len(ik)                                  ! Find number of user xmin, xmax
   IF(ABS(x(offxy(ik-1)+i)-xxmin)<DELTA) imin = i
   IF(ABS(x(offxy(ik-1)+i)-xxmax)<DELTA) THEN
      imax = i
      EXIT minmax
   ENDIF
ENDDO minmax
!
iz = iz - 1                                            ! This data set was temporary
!
! Now reserve space for actual map
WRITE(string,'(a,i10,a,i10)') 'map.data,', imax-imin+1, ',',INT((k2dm_end - k2dm_start)/k2dm_step)+1
length = LEN_TRIM(string)
CALL do_allocate(string,length, .FALSE.)               ! Allocate space for the map
IF(ier_num/=0) THEN
   ier_msg(1) = 'Error allocating space for map'
   write(*,*) ' STRING ', string(1:len_trim(string)), ' ' , length
   RETURN
ENDIF
ikm = iz - 1                                           ! Data set with map
!
IF(k2dm_line_b/=' ') THEN                              ! User requested background subtraction
   string = k2dm_line_b
   length = LEN_TRIM(string)
   CALL do_load(string,length, l_no_echo)
   IF(ier_num/=0) THEN
      ier_msg(1) = 'Error reading background file '
      RETURN
   ENDIF
   ikb = iz - 1
ELSE
   ikb = 0
ENDIF
!
! Main loop
!
iy = 0                                        ! Counter for input files
main:DO i=k2dm_start, k2dm_end, k2dm_step
   var_val(VAR_LOOP) = i
   string = k2dm_line                         ! Copy 'load' line for processing
   lp = LEN_TRIM(string)
   CALL do_load(string,lp, l_no_echo)         ! Regular kuplot/load
   IF(ier_num/=0) THEN
      IF(ier_typ==ER_IO .AND. ier_num>=-3) THEN    ! Handle input I/O errors
         IF(k2dm_miss == K2DM_ERROR) THEN          ! Stop processing
            WRITE(ier_msg(1),'(a,i5)') 'Error loading data set no ', i
            RETURN
         ELSEIF(k2dm_miss == K2DM_IGNORE) THEN     ! Skip this file, do a cycle
            CALL no_error
            CYCLE main
         ELSEIF(k2dm_miss == K2DM_BLANK) THEN      ! Flag a blank line
            CALL no_error
            lblank=.TRUE.
         ENDIF
      ELSE                                         ! All other errors
         WRITE(ier_msg(1),'(a,i5)') 'Error loading data set no', i
         RETURN
      ENDIF
   ENDIF
   ik = iz - 1                                     ! Current data set number is ik
   iy = iy + 1                                     ! Increment y-counter
   iix = 0
   DO ix=imin, imax                                    ! Copy current line into map
      iix = iix + 1
      IF(lblank) THEN                                  ! IO error and blank
         yy = 0
      ELSE
         IF(ikb/=0) THEN                                  ! Subtract background
            yy = y(offxy(ik-1) + ix) - y(offxy(ikb-1) + ix)* k2dm_scale
         ELSE                                             ! No background
            yy = y(offxy(ik-1) + ix)
         ENDIF
      ENDIF
      z(offz(ikm - 1) + (iix- 1) * ny(ikm) + iy) = yy  ! Finally place value
   ENDDO
   IF(k2dm_type_yf==K2DM_YF_LOOP) THEN
      y(offxy(ikm - 1)+iy) = i                         ! Set map's y-position
   ELSEIF(k2dm_type_yf==K2DM_YF_INC) THEN
      y(offxy(ikm - 1)+iy) = iy                        ! Set map's y-position
   ELSEIF(k2dm_type_yf==K2DM_YF_FUNC) THEN
      CALL ber_params (ianz, cpara, lpara, werte, MAXW) 
      y(offxy(ikm - 1)+iy) = werte(1)                  ! Set map's y-position
   ENDIF
!
   IF(lblank) THEN                                     ! Handle blank status
      lblank = .FALSE.
   ELSE
      iz = iz - 1                                      ! This data set was temporary
   ENDIF
ENDDO main
!
iix = 0
DO ix=imin, imax                                       ! Take x-values from lst data set
   iix = iix + 1
   x(offxy(ikm - 1) + iix) = x(offxy(ik - 1) + ix)
ENDDO
IF(ikb/=0) THEN
  iz = iz - 1                                          ! Background was additional file, reduze iz by one
ENDIF
!
!
CALL show_data(ikm)
!
END SUBROUTINE k2dm_run
!
!*******************************************************************************
!
SUBROUTINE k2dm_reset
!
IMPLICIT NONE
!
k2dm_ctype  = ' '
k2dm_line   = ' '
k2dm_line_b = ' '
k2dm_line_yf = ' '
k2dm_type      = 0
k2dm_miss      = -1
k2dm_start     = 0
k2dm_end       = 0
k2dm_step      = 0
k2dm_lxmin     = .TRUE.
k2dm_lxmax     = .TRUE.
k2dm_scale     = 1.0
k2dm_xmin      = 0.0
k2dm_xmax      = 0.0
k2dm_type_yf   = K2DM_YF_LOOP
!
END SUBROUTINE k2dm_reset
!
!*******************************************************************************
!
END MODULE kuplot_2dm
