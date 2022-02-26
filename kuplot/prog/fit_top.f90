MODULE kuplot_fit6
!
! Revised as of DISCUS Version 6.0.0 
! Levenberg Marquadt type Least Squares fits within KUPLOT
! The procedure is adapted from the Numerical recipes.
!
!
PRIVATE
PUBLIC do_f66!, kupl_theory, write_fit
!
CONTAINS
!MODULE kuplot_fit_mod
!
!CONTAINS
!*****7*****************************************************************
!     This sublevel includes all commands and functions for the         
!     least square fits in KUPLOT.                                      
!*****7*****************************************************************
SUBROUTINE do_f66(zei, lp)
!                                                                       
!     Main fitting menu                                                 
!                                                                       
USE kuplot_theory_macro_mod
!
USE kuplot_config 
USE kuplot_fit_para
USE kuplot_mod 
!use kuplot_draw_mod
use kuplot_extrema_mod
use kuplot_para_mod
use kuplot_fit6_low_mod
use kuplot_plot_mod
!
USE ber_params_mod
USE calc_expr_mod
USE doact_mod 
USE do_eval_mod
USE do_set_mod
USE do_wait_mod
USE errlist_mod 
USE get_params_mod
USE learn_mod 
USE class_macro_internal
USE prompt_mod 

USE lib_do_operating_mod
USE lib_echo
USE lib_errlist_func
USE lib_help
USE lib_length
USE lib_macro_func
USE str_comp_mod
USE string_convert_mod
USE sup_mod
USE precision_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER, PARAMETER :: maxw = 10 
!                                                                       
CHARACTER (LEN=*), INTENT(INOUT) :: zei 
INTEGER          , INTENT(INOUT) :: lp
!
CHARACTER(LEN=PREC_STRING), DIMENSION(MAXW) ::cpara
CHARACTER(LEN=PREC_STRING) :: line, zeile
CHARACTER(LEN=40  ) :: orig_prompt 
CHARACTER(LEN=40)   ::cdummy 
CHARACTER(LEN=4)    ::befehl 
CHARACTER(LEN=1)    :: empty 
INTEGER, DIMENSION(MAXW) :: lpara
INTEGER :: ll
integer :: i
INTEGER :: ianz, indxg, lbef 
INTEGER :: maxpkt, maxzz 
REAL(KIND=PREC_DP)   , DIMENSION(MAXW) :: werte
LOGICAL, DIMENSION(3) :: flag
LOGICAL               :: sel_func 
!                                                                       
real(kind=PREC_DP) :: f, df(maxpara)
!
empty = ' '
CALL no_error 
sel_func = ftyp (1:4) /= 'NONE' 
!                                                                       
!------ check data set number and free space                            
!                                                                       
CALL get_params (zei, ianz, cpara, lpara, maxw, lp) 
IF (ier_num.ne.0) RETURN 
!                                                                       
      IF (ianz.eq.1) THEN 
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         IF (ier_num.ne.0) RETURN 
         ikfit = nint (werte (1) ) 
         IF (ikfit.lt.1.or.ikfit.gt. (iz - 1) ) THEN 
            ier_num = - 4 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
         IF ( (iz + 1) .gt.maxkurvtot) THEN 
            ier_num = - 1 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
         RETURN 
      ENDIF 
!                                                                       
      IF (ikfirst (ikfit) ) THEN 
!                                                                       
!----- -- Check if there is enough space left                           
!                                                                       
         maxpkt = maxarray - offxy (iz - 1) 
         maxzz = maxarray - offz (iz - 1) 
!                                                                       
         IF (lni (ikfit) ) THEN 
            IF ( (2.0 * nx (iz - 1) * ny (iz - 1) .gt.maxzz) .or. (max (&
            nx (iz - 1), ny (iz - 1) ) .gt.maxpkt) ) THEN               
               ier_num = - 6 
               ier_typ = ER_APPL 
               RETURN 
            ENDIF 
!                                                                       
            ikcal = iz 
            ikdif = iz + 1 
            offxy (ikcal) = offxy (iz - 1) + max (nx (iz - 1), ny (iz - &
            1) )                                                        
            offxy (ikdif) = offxy (ikcal) + max (nx (iz - 1), ny (iz -  &
            1) )                                                        
            offz (ikcal) = offz (iz - 1) + nx (iz - 1) * ny (iz - 1) 
            offz (ikdif) = offz (ikcal) + nx (iz - 1) * ny (iz - 1) 
            nx (ikcal) = nx (ikfit) 
            nx (ikdif) = nx (ikfit) 
            ny (ikcal) = ny (ikfit) 
            ny (ikdif) = ny (ikfit) 
            xmin (ikcal) = xmin (ikfit) 
            xmin (ikdif) = xmin (ikfit) 
            ymin (ikcal) = ymin (ikfit) 
            ymin (ikdif) = ymin (ikfit) 
            xmax (ikcal) = xmax (ikfit) 
            xmax (ikdif) = xmax (ikfit) 
            ymax (ikcal) = ymax (ikfit) 
            ymax (ikdif) = ymax (ikfit) 
            iz = iz + 2 
         ELSE 
            IF (2.0 * lenc(iz - 1) .gt.maxpkt) THEN 
               ier_num = - 6 
               ier_typ = ER_APPL 
               RETURN 
            ENDIF 
!                                                                       
            ikcal = iz 
            ikdif = iz + 1 
            offxy (ikcal) = offxy (iz - 1) + lenc(iz - 1) 
            offxy (ikdif) = offxy (ikcal) + lenc(iz - 1) 
            offz (ikcal) = offz (iz - 1) 
            offz (ikdif) = offz (ikcal) 
            lenc(ikcal) = lenc(ikfit) 
            lenc(ikdif) = lenc(ikfit) 
            xmin (ikcal) = xmin (ikfit) 
            xmin (ikdif) = xmin (ikfit) 
            xmax (ikcal) = xmax (ikfit) 
            xmax (ikdif) = xmax (ikfit) 
            iz = iz + 2 
         ENDIF 
!                                                                       
!----- -- Everything ok so far ..                                       
!                                                                       
         ikfirst (ikfit) = .false. 
         fit_ikcal (ikfit) = ikcal 
         fit_ikdif (ikfit) = ikdif 
!                                                                       
!----- -- Restore previous fit setting                                  
!                                                                       
      ELSE 
         ikcal = fit_ikcal (ikfit) 
         ikdif = fit_ikdif (ikfit) 
      ENDIF 
      orig_prompt = prompt
      prompt = prompt (1:len_str (prompt) ) //'/fit' 
!                                                                       
!------ here starts sublevel fit                                        
!                                                                       
   10 CONTINUE 
!                                                                       
      CALL get_cmd (line, ll, befehl, lbef, zeile, lp, prompt) 
      IF (ier_num.eq.0) THEN 
         IF (line.eq.' '.or.line (1:1) .eq.'#' .OR. line(1:1)=='!') goto 10 
!                                                                       
!------ search for "="                                                  
!                                                                       
         indxg = index (line, '=') 
      IF (indxg.ne.0                                            &
         .and..not. (str_comp (befehl, 'echo', 2, lbef, 4) )    &
         .and..not. (str_comp (befehl, 'help', 2, lbef, 4) .or. &
                     str_comp (befehl, '?   ', 2, lbef, 4) )    &
         .AND. INDEX(line,'==') == 0                          )THEN
            CALL do_math (line, indxg, ll) 
!                                                                       
!------ execute a macro file                                            
!                                                                       
         ELSEIF (befehl (1:1) .eq.'@') THEN 
            line(1:ll-1) = line(2:ll)
            ll   = ll - 1
            CALL file_kdo (line, ll) 
!                                                                       
!     continues a macro 'continue'                                      
!                                                                       
         ELSEIF (str_comp (befehl, 'continue', 3, lbef, 8) ) THEN 
            CALL macro_continue (zeile, lp) 
!                                                                       
!-------Set number of cycles 'cyc'                                      
!                                                                       
         ELSEIF (str_comp (befehl, 'cycle', 2, lbef, 5) ) THEN 
            CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
            IF (ier_num.eq.0) THEN 
               IF (ianz.eq.1) THEN 
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  IF (ier_num.eq.0) THEN 
                     ncycle = nint (werte (1) ) 
                  ENDIF 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
            ENDIF 
!                                                                       
!------ Echo a string, just for interactive check in a macro 'echo'     
!                                                                       
         ELSEIF (str_comp (befehl, 'echo', 2, lbef, 4) ) THEN 
            CALL echo (zeile, lp) 
!                                                                       
!------ Evaluate an expression                                          
!                                                                       
         ELSEIF (str_comp (befehl, 'eval', 2, lbef, 4) ) THEN 
            CALL do_eval (zeile, lp, .TRUE.) 
!                                                                       
!     exit 'exit'                                                       
!                                                                       
         ELSEIF (str_comp (befehl, 'exit', 2, lbef, 4) ) THEN 
            GOTO 9999 
!                                                                       
!     Define fit function 'func'                                        
!                                                                       
         ELSEIF (str_comp (befehl, 'func', 3, lbef, 4) ) THEN 
            lturn_off = .FALSE.
            CALL do_fit_fkt (zeile, lp) 
            IF (ier_num.eq.0) sel_func = .true. 
!                                                                       
!     help 'help','?'                                                   
!                                                                       
      ELSEIF (str_comp (befehl, 'help', 2, lbef, 4) .or.str_comp (befehl&
     &, '?   ', 1, lbef, 4) ) THEN                                      
            IF (str_comp (zeile, 'errors', 2, lp, 6) ) THEN 
               lp = lp + 7 
               CALL do_hel ('kuplot '//zeile, lp) 
            ELSE 
               lp = lp + 12 
               CALL do_hel ('kuplot fit '//zeile, lp) 
            ENDIF 
!                                                                       
!-------Save current parameters to a macro file                         
!                                                                       
         ELSEIF (str_comp (befehl, 'macro', 2, lbef, 5) ) THEN 
            CALL do_fit_macro (zeile, lp) 
!                                                                       
!-------Set parameter 'ifen' for determination of maxima                
!                                                                       
         ELSEIF (str_comp (befehl, 'mfen', 2, lbef, 4) ) THEN 
            CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
            IF (ier_num.eq.0) THEN 
               IF (ianz.eq.1) THEN 
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  IF (ier_num.eq.0) THEN 
                     fit_ifen = nint (werte (1) ) 
                  ENDIF 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
            ENDIF 
!                                                                       
!-------Toogle fit screen output 'output'                               
!                                                                       
         ELSEIF (str_comp (befehl, 'output', 2, lbef, 6) ) THEN 
            CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
            IF (ianz.eq.1) THEN 
               fstart = str_comp (cpara (1) , 'on', 2, lpara (1) , 2) 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
!                                                                       
!-------Toogle fit range 'range'                                        
!                                                                       
         ELSEIF (str_comp (befehl, 'range', 2, lbef, 5) ) THEN 
            CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
            IF (ianz.eq.1) THEN 
               frall = str_comp (cpara (1) , 'all', 2, lpara (1) , 2) 
            ELSE 
               IF (frall) THEN 
                  WRITE (output_io, 2100) 'complete range' 
               ELSE 
                  WRITE (output_io, 2100) 'plot range only' 
               ENDIF 
            ENDIF 
!                                                                       
!-------Set parameters                                                  
!                                                                       
         ELSEIF (str_comp (befehl, 'par', 2, lbef, 3) ) THEN 
            CALL do_fit_par (zeile, lp) 
!                                                                       
!-------Plot result                                                     
!                                                                       
         ELSEIF (str_comp (befehl, 'plot', 2, lbef, 4) ) THEN 
            CALL do_plot (.false.) 
!                                                                       
!-------Set lamda  'relax'                                              
!                                                                       
         ELSEIF (str_comp (befehl, 'relax', 2, lbef, 5) ) THEN 
            CALL kuplot_set_lamda(zeile, lp)
!                                                                       
!-------Set convergence criteria 'conv'
!                                                                       
         ELSEIF (str_comp (befehl, 'conv', 2, lbef, 5) ) THEN 
            CALL kuplot_set_convergence(zeile, lp)
!                                                                       
!------ Set scale                                                       
!                                                                       
         ELSEIF (str_comp (befehl, 'skal', 2, lbef, 4) ) THEN 
            CALL set_skal (zeile, lp) 
!                                                                       
!-------Run fit                                                         
!                                                                       
         ELSEIF (str_comp (befehl, 'run', 2, lbef, 3) ) THEN 
!
            lturn_off = .TRUE.
            IF (.not.sel_func) THEN 
               ier_num = - 25 
               ier_typ = ER_APPL 
            ELSE 
               IF(ex(iwin,iframe,1)==ex(iwin,iframe,2)) THEN
                  i = 0
                  CALL set_skal(empty,i)
               ENDIF
               CALL do_fit
!              IF (lni (ikfit) ) THEN 
!                 CALL do_fit_z 
!              ELSE 
!                 CALL do_fit_y 
!              ENDIF 
               CALL get_extrema 
            ENDIF 
            lturn_off = .FALSE.
!                                                                  
!-------Test a recursive call to kuplot_loop                            
!                                                                       
         ELSEIF (str_comp (befehl, 'test', 2, lbef, 4) ) THEN 
            CALL theory_macro(1.0D0,  f, df, 1)
!                                                                  
!-------Save fit results                                                
!                                                                       
         ELSEIF (str_comp (befehl, 'save', 2, lbef, 4) ) THEN 
            CALL do_fit_save 
!                                                                       
!-------Show settings                                                   
!                                                                       
         ELSEIF (str_comp (befehl, 'show', 2, lbef, 4) ) THEN 
            CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
            IF (ier_num.eq.0) THEN 
               IF (ianz.eq.0) THEN 
                  CALL do_fit_info (output_io, .true., .true., .true.) 
               ELSEIF (ianz.eq.1) THEN 
                  CALL do_cap (cpara (1) ) 
                  flag (1) = cpara (1) (1:2) .eq.'GE' 
                  flag (2) = cpara (1) (1:2) .eq.'FI' 
                  flag (3) = cpara (1) (1:2) .eq.'PA' 
                  IF (flag (1) .or.flag (2) .or.flag (3) ) THEN 
                     CALL do_fit_info (output_io, flag (1), flag (2),   &
                     flag (3) )                                         
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
            ENDIF 
!                                                                       
!-------Operating System Kommandos 'syst'                               
!                                                                       
         ELSEIF (str_comp (befehl, 'syst', 2, lbef, 4) ) THEN 
            cdummy = ' ' 
            IF (zeile.ne.' ') THEN 
               cdummy (1:lp) = zeile (1:lp) 
               CALL do_operating (cdummy, lp) 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
!                                                                       
!-------Set URF 'urf'                                                   
!                                                                       
         ELSEIF (str_comp (befehl, 'urf', 2, lbef, 3) ) THEN 
            CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
            IF (ier_num.eq.0) THEN 
               IF (ianz.eq.1) THEN 
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  IF (ier_num.eq.0) THEN 
                     urf = werte (1) 
!                     kup_fit6_lamda_s = werte(1)
                  ENDIF 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
            ENDIF 
!                                                                       
!     Waiting for user input                                            
!                                                                       
         ELSEIF (str_comp (befehl, 'wait', 3, lbef, 4) ) THEN 
            CALL do_input (zeile, lp) 
!                                                                       
!     Set weighting scheme                                              
!                                                                       
         ELSEIF (str_comp (befehl, 'wic', 3, lbef, 3) ) THEN 
            CALL do_fit_wichtung (zeile, lp) 
!                                                                       
!------ no command found                                                
!                                                                       
         ELSE 
            ier_num = - 8 
            ier_typ = ER_COMM 
         ENDIF 
      ENDIF 
!                                                                       
!------ any errors ?                                                    
!                                                                       
      IF (ier_num.ne.0) THEN 
         CALL errlist 
         IF (ier_sta.ne.ER_S_LIVE) THEN 
            IF (lmakro .OR. lmakro_error) THEN 
               IF(sprompt /= prompt ) THEN
                  ier_num = -10
                  ier_typ = ER_COMM
                  ier_msg(1) = ' Error occured in fit menu'
                  prompt_status = PROMPT_ON 
                  prompt = orig_prompt
                  RETURN
               ELSE
                  CALL macro_close 
                  prompt_status = PROMPT_ON 
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
 2100 FORMAT     (1x,'Refinement mode: ',a) 
      END SUBROUTINE do_f66                         
!
end MODULE kuplot_fit6
