MODULE kuplot_3dm
!
USE kuplot_3dm_mod
!
PRIVATE
PUBLIC kuplot_3dm_menu
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE kuplot_3dm_menu(string, length)
!
CHARACTER(LEN=*), INTENT(INOUT) :: stringth
INTEGER         , INTENT(INOUT) :: length
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
USE lib_errlist_func
USE lib_help
USE lib_do_operating_mod
USE lib_macro_func
USE precision_mod
USE prompt_mod
USE str_comp_mod
USE sup_mod
USE take_param_mod
!
IMPLICIT NONE
!
CHARACTER (LEN=5)                       :: befehl! command on input line
CHARACTER(LEN=LEN_TRIM(prompt))         :: orig_prompt  ! original prompt
CHARACTER (LEN=PREC_STRING)             :: line  ! input line
CHARACTER (LEN=PREC_STRING)             :: zeile ! remainder with parameters
INTEGER                                 :: indxg ! location of "="
INTEGER                                 :: lp    ! length of zeile
INTEGER                                 :: laenge
INTEGER                                 :: lbef
LOGICAL                                 :: lend
!
!
INTEGER, PARAMETER :: MAXW = 20
CHARACTER(LEN=PREC_STRING), DIMENSION(MAXW) :: cpara
INTEGER            , DIMENSION(MAXW) :: lpara
REAL               , DIMENSION(MAXW) :: werte
INTEGER                              :: ianz
!
!INTEGER, EXTERNAL :: len_str
!
CALL get_params (line, ianz, cpara, lpara, MAXW, length) 
IF (ier_num /= 0) RETURN 
CALL ber_params (ianz, cpara, lpara, werte, MAXW) 
IF (ier_num /= 0) RETURN 
k3dm_ik = NINT(werte(1))
!
orig_prompt = prompt
prompt = prompt (1:LEN_TRIM(prompt) ) //'/3dm'
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
                  line = line(2:laenge)
                  laenge = 1
                  CALL file_kdo(line, laenge)
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
               CALL k3dm_reset(zeile,lp)
            ELSE is_generic    ! macro, reset or all other commands is_generic
!
               is_com: IF(str_comp (befehl, 'show', 3, lbef, 4) ) THEN 
!
                  WRITE(output_io,'(a)') ' Will show 3dm '
               ELSEIF(str_comp (befehl, 'direction', 3, lbef, 9) ) THEN is_com
                  CALL get_direction(zeile, lp)
               ELSEIF(str_comp (befehl, 'plot', 3, lbef, 4) ) THEN  is_com
                  CALL k3dm_run(zeile, lp)
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
               ier_msg(1) = 'Error occured in 3dm menu'
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
END SUBROUTINE kuplot_3dm_menu
!
!*******************************************************************************
!
SUBROUTINE get_direction(line, length)
!
USE errlist_mod
USE get_params_mod
USE lib_errlist_func
USE precision_mod
USE str_comp_mod
USE take_param_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: line
INTEGER         , INTENT(INOUT) :: length
!
INTEGER, PARAMETER :: O_PHI   = 1
INTEGER, PARAMETER :: O_RHO   = 2
!
INTEGER, PARAMETER :: NOPTIONAL = 2
CHARACTER(LEN=   3), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=PREC_STRING), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent!opt. para present
REAL               , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 2 ! Number of values to calculate 
!
INTEGER, PARAMETER :: MAXW = MAX(20,NOPTIONAL)
CHARACTER(LEN=PREC_STRING), DIMENSION(MAXW) :: cpara
INTEGER            , DIMENSION(MAXW) :: lpara
REAL               , DIMENSION(MAXW) :: werte
INTEGER                              :: ianz
!
!
DATA oname  / 'rho', 'phi' /
DATA loname /  3   ,  3    /
opara  =  (/ '0.000', '0.000' /)   ! Always provide fresh default values
lopara =  (/  5     ,  5      /)
owerte =  (/  1.0   ,  1.0    /)
!
CALL no_error 
CALL get_params (line, ianz, cpara, lpara, MAXW, length) 
IF (ier_num /= 0) RETURN 
CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
!
k3dm_phi   = owerte(O_PHI)
k3dm_rho   = owerte(O_RHO)
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
USE lib_errlist_func
USE precision_mod
USE str_comp_mod
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
CHARACTER(LEN=   4), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=PREC_STRING), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent!opt. para present
REAL               , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 0 ! Number of values to calculate 
!
INTEGER, PARAMETER :: MAXW = MAX(20,NOPTIONAL)
CHARACTER(LEN=PREC_STRING), DIMENSION(MAXW) :: cpara
INTEGER            , DIMENSION(MAXW) :: lpara
REAL               , DIMENSION(MAXW) :: werte
INTEGER                              :: ianz
!
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
   k3dm_lxmin = .TRUE.
ELSE
   k3dm_lxmin = .FALSE.
   cpara(1) = opara(O_XMIN)
   lpara(1) = lopara(O_XMIN)
   ianz = 1
   CALL ber_params (ianz, cpara, lpara, werte, MAXW) 
   IF(ier_num==0) THEN
      k3dm_xmin = werte(1)
   ELSE
      RETURN
   ENDIF
ENDIF
IF(opara(O_XMAX)=='xmax') THEN
   k3dm_lxmax = .TRUE.
ELSE
   k3dm_lxmax = .FALSE.
   cpara(1) = opara(O_XMAX)
   lpara(1) = lopara(O_XMAX)
   ianz = 1
   CALL ber_params (ianz, cpara, lpara, werte, MAXW) 
   IF(ier_num==0) THEN
      k3dm_xmax = werte(1)
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
USE lib_errlist_func
USE precision_mod
USE str_comp_mod
USE take_param_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: line
INTEGER         , INTENT(INOUT) :: length
!
INTEGER, PARAMETER :: O_LOOP = 1
!
INTEGER, PARAMETER :: NOPTIONAL = 1
CHARACTER(LEN=   4), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=PREC_STRING), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent!opt. para present
REAL               , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 0 ! Number of values to calculate 
!
INTEGER, PARAMETER :: MAXW = MAX(20,NOPTIONAL)
CHARACTER(LEN=PREC_STRING), DIMENSION(MAXW) :: cpara
INTEGER            , DIMENSION(MAXW) :: lpara
REAL               , DIMENSION(MAXW) :: werte
INTEGER                              :: ianz
!
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
   k3dm_type_yf = K3DM_YF_LOOP
   k3dm_line_yf = 'LOOP'
ELSEIF(opara(O_LOOP)=='inc') THEN
   k3dm_type_yf = K3DM_YF_INC
   k3dm_line_yf = 'inc'
ELSE
   k3dm_type_yf = K3DM_YF_FUNC
   k3dm_line_yf = opara(O_LOOP)
   continue
ENDIF
!
END SUBROUTINE get_yfunc 
!
!*******************************************************************************
!
SUBROUTINE k3dm_run(line, length)
!
USE kuplot_mod
!
USE ber_params_mod
USE errlist_mod
USE get_params_mod
USE lib_errlist_func
USE precision_mod
USE take_param_mod
USE variable_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: line
INTEGER         , INTENT(INOUT) :: length
!
LOGICAL, PARAMETER :: l_no_echo = .FALSE.
CHARACTER(LEN=PREC_STRING) :: string
CHARACTER(LEN=10  ) :: cnumber
CHARACTER(LEN=4   ) :: bef
INTEGER :: i, ik
INTEGER :: ikm                  ! Data set with map
INTEGER :: ikb                  ! Data set with background
INTEGER :: ix, iy
INTEGER :: iix                  ! Running index for map x
INTEGER :: lp
INTEGER :: imin, imax
LOGICAL :: lblank               ! IO error and k3dm_miss=blank
REAL    :: xxmin, xxmax
REAL    :: yy                   ! background corrected y value of a single 1D
REAL    :: DELTA
!
INTEGER, PARAMETER :: O_PLOT = 1
!
INTEGER, PARAMETER :: NOPTIONAL = 1
CHARACTER(LEN=   4), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=PREC_STRING), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent!opt. para present
REAL               , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 0 ! Number of values to calculate 
!
!
INTEGER, PARAMETER :: MAXW = 20
CHARACTER(LEN=PREC_STRING), DIMENSION(MAXW) :: cpara
INTEGER            , DIMENSION(MAXW) :: lpara
REAL               , DIMENSION(MAXW) :: werte
INTEGER                              :: ianz
!
DATA oname  / 'plot'/
DATA loname /  4    /
opara  =  (/ 'auto' /)   ! Always provide fresh default values
lopara =  (/  4     /)
owerte =  (/  0.0   /)
!
CALL no_error 
CALL get_params (line, ianz, cpara, lpara, MAXW, length) 
IF (ier_num /= 0) RETURN 
CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
!
!
IF(k3dm_type_yf==K3DM_YF_FUNC) THEN                       ! yfunc is a function get params
   IF(k3dm_line_yf(1:1) == '''') k3dm_line_yf(1:1) = ' '
   i = LEN_TRIM(k3dm_line_yf)
   IF(k3dm_line_yf(i:i) == '''') k3dm_line_yf(i:i) = ' '
   length = LEN_TRIM(k3dm_line_yf)
   CALL get_params (k3dm_line_yf, ianz, cpara, lpara, MAXW, length) 
   IF(ier_num /= 0) THEN
      ier_msg(1) = 'Error getting prameters of yfunction'
      RETURN
   ENDIF 
ENDIF
!
ikm = 0
ikb = 0
!
var_val(VAR_LOOP) = k3dm_start
string = k3dm_line
lp = LEN_TRIM(string)
CALL do_load(string,lp, l_no_echo)
IF(ier_num/=0) THEN
   ier_msg(1) = 'Error loading first data set '
   RETURN
ENDIF
ik   = iz - 1                 ! This is the number of the current data set
imin = 1                      ! Default to first data point
imax = lenc(ik)                ! Default to length of data sets
IF(k3dm_lxmin) THEN           ! user specified xmin:xmin
   xxmin = x(offxy(ik-1)+1)
   imin  = 1
ELSE
   xxmin = k3dm_xmin
ENDIF
IF(k3dm_lxmax) THEN           ! user specified xmin:xmin
   xxmax = x(offxy(ik-1)+lenc(ik))
   imin  = 1
ELSE
   xxmax = k3dm_xmax
ENDIF
DELTA = ABS(x(offxy(ik-1)+2)-x(offxy(ik-1)+1))*0.01    ! Calculate a sigma as (x(2)-x(1))/10
minmax:DO i=1,lenc(ik)                                  ! Find number of user xmin, xmax
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
WRITE(string,'(a,i10,a,i10)') 'map.data,', imax-imin+1, ',',INT((k3dm_end - k3dm_start)/k3dm_step)+1
length = LEN_TRIM(string)
CALL do_allocate(string,length, .FALSE.)               ! Allocate space for the map
IF(ier_num/=0) THEN
   ier_msg(1) = 'Error allocating space for map'
   RETURN
ENDIF
ikm = iz - 1                                           ! Data set with map
!
IF(k3dm_line_b/=' ') THEN                              ! User requested background subtraction
   string = k3dm_line_b
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
main:DO i=k3dm_start, k3dm_end, k3dm_step
   var_val(VAR_LOOP) = i
   string = k3dm_line                         ! Copy 'load' line for processing
   lp = LEN_TRIM(string)
   CALL do_load(string,lp, l_no_echo)         ! Regular kuplot/load
   IF(ier_num/=0) THEN
      IF(ier_typ==ER_IO .AND. ier_num>=-3) THEN    ! Handle input I/O errors
         IF(k3dm_miss == K3DM_ERROR) THEN          ! Stop processing
            WRITE(ier_msg(1),'(a,i5)') 'Error loading data set no ', i
            RETURN
         ELSEIF(k3dm_miss == K3DM_IGNORE) THEN     ! Skip this file, do a cycle
            CALL no_error
            CYCLE main
         ELSEIF(k3dm_miss == K3DM_BLANK) THEN      ! Flag a blank line
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
            yy = y(offxy(ik-1) + ix) - y(offxy(ikb-1) + ix)* k3dm_scale
         ELSE                                             ! No background
            yy = y(offxy(ik-1) + ix)
         ENDIF
      ENDIF
      z(offz(ikm - 1) + (iix- 1) * ny(ikm) + iy) = yy  ! Finally place value
   ENDDO
   IF(k3dm_type_yf==K3DM_YF_LOOP) THEN
      y(offxy(ikm - 1)+iy) = i                         ! Set map's y-position
   ELSEIF(k3dm_type_yf==K3DM_YF_INC) THEN
      y(offxy(ikm - 1)+iy) = iy                        ! Set map's y-position
   ELSEIF(k3dm_type_yf==K3DM_YF_FUNC) THEN
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
IF(opara(O_PLOT) == 'auto') THEN
   string = ' '
   length = 0
   bef = 'aver'
   CALL para_setr(string, length, yskal_u(iwin, iframe), bef, -9999.)
   lyskal(iwin, iframe) = (yskal_u(iwin, iframe) /=  -9999.)
!
   string = ' '
   length = 0
   CALL set_skal(string, length)
!
   string = ' '
   length = 0
   CALL set_mark(string, length)
!
   string = 'fire'
   length = 4
   CALL set_cmap(string, length) 
!
   string = '1, 2'
   length = 4
   CALL para_seti(string, length, hlineart, 1, maxkurvtot, bef, 1, 4, .FALSE.)
!
   string = '1, 1,1,100,%'
   length = 12
   CALL set_hlin(string, length) 
!
   string = 'Test of title 1 '
   length = 16
   CALL para_set_title (string, length, titel (iwin, iframe, 1) ) 
!
   CALL do_plot(.FALSE.)
!
ENDIF
!
END SUBROUTINE k3dm_run
!
!*******************************************************************************
!
SUBROUTINE k3dm_reset
!
IMPLICIT NONE
!
k3dm_ctype  = ' '
k3dm_line   = ' '
k3dm_line_b = ' '
k3dm_line_yf = ' '
k3dm_type      = 0
k3dm_miss      = -1
k3dm_start     = 0
k3dm_end       = 0
k3dm_step      = 0
k3dm_lxmin     = .TRUE.
k3dm_lxmax     = .TRUE.
k3dm_scale     = 1.0
k3dm_xmin      = 0.0
k3dm_xmax      = 0.0
k3dm_type_yf   = K3DM_YF_LOOP
!
END SUBROUTINE k3dm_reset
!
!*******************************************************************************
!
END MODULE kuplot_3dm
