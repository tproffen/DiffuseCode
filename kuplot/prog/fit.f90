MODULE kuplot_fit_const
!
! Numbers that are needed by some fit routines
!
CHARACTER(LEN=1024) :: f6_fit_func  = ' '  ! User defined fit function
INTEGER             :: f6_fit_lfunc = 0    ! Fit function string length
REAL                :: pp_origin    = 0.0  ! Origin of backgroud polynomial for Pseudovoigt
INTEGER             :: nn_backgrd   = 0    ! number of background params for Pseudovoigt
!
END MODULE kuplot_fit_const
!
!MODULE kuplot_fit_mod
!
!CONTAINS
!*****7*****************************************************************
!     This sublevel includes all commands and functions for the         
!     least square fits in KUPLOT.                                      
!*****7*****************************************************************
SUBROUTINE do_fit (zei, lp)
!                                                                       
!     Main fitting menu                                                 
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
      USE kuplot_config 
      USE kuplot_mod 
      USE string_convert_mod
      USE sup_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 10) 
!                                                                       
 CHARACTER (LEN=*), INTENT(INOUT) :: zei 
 INTEGER          , INTENT(INOUT) :: lp
!
      CHARACTER(1024) cpara (maxw) 
      CHARACTER(1024) line, zeile
      CHARACTER(LEN=40         ) :: orig_prompt 
      CHARACTER(40) cdummy 
      CHARACTER(4) befehl 
      CHARACTER(LEN=1) :: empty 
      INTEGER lpara (maxw) 
      INTEGER ll
      INTEGER ianz, indxg, lbef 
      INTEGER maxpkt, maxzz 
      REAL werte (maxw) 
      LOGICAL flag (3), sel_func 
!                                                                       
      INTEGER len_str 
      LOGICAL str_comp 
real :: f, df(maxpara)
!
      empty = ' '
      CALL no_error 
      sel_func = ftyp (1:4) .ne.'NONE' 
!                                                                       
!------ check data set number and free space                            
!                                                                       
      CALL get_params (zei, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
!                                                                       
      IF (ianz.eq.1) then 
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         IF (ier_num.ne.0) return 
         ikfit = nint (werte (1) ) 
         IF (ikfit.lt.1.or.ikfit.gt. (iz - 1) ) then 
            ier_num = - 4 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
         IF ( (iz + 1) .gt.maxkurvtot) then 
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
      IF (ikfirst (ikfit) ) then 
!                                                                       
!----- -- Check if there is enough space left                           
!                                                                       
         maxpkt = maxarray - offxy (iz - 1) 
         maxzz = maxarray - offz (iz - 1) 
!                                                                       
         IF (lni (ikfit) ) then 
            IF ( (2.0 * nx (iz - 1) * ny (iz - 1) .gt.maxzz) .or. (max (&
            nx (iz - 1), ny (iz - 1) ) .gt.maxpkt) ) then               
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
            IF (2.0 * len (iz - 1) .gt.maxpkt) then 
               ier_num = - 6 
               ier_typ = ER_APPL 
               RETURN 
            ENDIF 
!                                                                       
            ikcal = iz 
            ikdif = iz + 1 
            offxy (ikcal) = offxy (iz - 1) + len (iz - 1) 
            offxy (ikdif) = offxy (ikcal) + len (iz - 1) 
            offz (ikcal) = offz (iz - 1) 
            offz (ikdif) = offz (ikcal) 
            len (ikcal) = len (ikfit) 
            len (ikdif) = len (ikfit) 
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
      IF (ier_num.eq.0) then 
         IF (line.eq.' '.or.line (1:1) .eq.'#' .OR. line(1:1)=='!') goto 10 
!                                                                       
!------ search for "="                                                  
!                                                                       
         indxg = index (line, '=') 
      IF (indxg.ne.0                                            &
         .and..not. (str_comp (befehl, 'echo', 2, lbef, 4) )    &
         .and..not. (str_comp (befehl, 'help', 2, lbef, 4) .or. &
                     str_comp (befehl, '?   ', 2, lbef, 4) )    &
         .AND. INDEX(line,'==') == 0                          )then
            CALL do_math (line, indxg, ll) 
!                                                                       
!------ execute a macro file                                            
!                                                                       
         ELSEIF (befehl (1:1) .eq.'@') then 
            CALL file_kdo (line (2:ll), ll - 1) 
!                                                                       
!     continues a macro 'continue'                                      
!                                                                       
         ELSEIF (str_comp (befehl, 'continue', 3, lbef, 8) ) then 
            CALL macro_continue (zeile, lp) 
!                                                                       
!-------Set number of cycles 'cyc'                                      
!                                                                       
         ELSEIF (str_comp (befehl, 'cycle', 2, lbef, 5) ) then 
            CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
            IF (ier_num.eq.0) then 
               IF (ianz.eq.1) then 
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  IF (ier_num.eq.0) then 
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
         ELSEIF (str_comp (befehl, 'echo', 2, lbef, 4) ) then 
            CALL echo (zeile, lp) 
!                                                                       
!------ Evaluate an expression                                          
!                                                                       
         ELSEIF (str_comp (befehl, 'eval', 2, lbef, 4) ) then 
            CALL do_eval (zeile, lp, .TRUE.) 
!                                                                       
!     exit 'exit'                                                       
!                                                                       
         ELSEIF (str_comp (befehl, 'exit', 2, lbef, 4) ) then 
            GOTO 9999 
!                                                                       
!     Define fit function 'func'                                        
!                                                                       
         ELSEIF (str_comp (befehl, 'func', 3, lbef, 4) ) then 
            lturn_off = .FALSE.
            CALL do_fit_fkt (zeile, lp) 
            IF (ier_num.eq.0) sel_func = .true. 
!                                                                       
!     help 'help','?'                                                   
!                                                                       
      ELSEIF (str_comp (befehl, 'help', 2, lbef, 4) .or.str_comp (befehl&
     &, '?   ', 1, lbef, 4) ) then                                      
            IF (str_comp (zeile, 'errors', 2, lp, 6) ) then 
               lp = lp + 7 
               CALL do_hel ('kuplot '//zeile, lp) 
            ELSE 
               lp = lp + 12 
               CALL do_hel ('kuplot fit '//zeile, lp) 
            ENDIF 
!                                                                       
!-------Save current parameters to a macro file                         
!                                                                       
         ELSEIF (str_comp (befehl, 'macro', 2, lbef, 5) ) then 
            CALL do_fit_macro (zeile, lp) 
!                                                                       
!-------Set parameter 'ifen' for determination of maxima                
!                                                                       
         ELSEIF (str_comp (befehl, 'mfen', 2, lbef, 4) ) then 
            CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
            IF (ier_num.eq.0) then 
               IF (ianz.eq.1) then 
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  IF (ier_num.eq.0) then 
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
         ELSEIF (str_comp (befehl, 'output', 2, lbef, 6) ) then 
            CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
            IF (ianz.eq.1) then 
               fstart = str_comp (cpara (1) , 'on', 2, lpara (1) , 2) 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
!                                                                       
!-------Toogle fit range 'range'                                        
!                                                                       
         ELSEIF (str_comp (befehl, 'range', 2, lbef, 5) ) then 
            CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
            IF (ianz.eq.1) then 
               frall = str_comp (cpara (1) , 'all', 2, lpara (1) , 2) 
            ELSE 
               IF (frall) then 
                  WRITE (output_io, 2100) 'complete range' 
               ELSE 
                  WRITE (output_io, 2100) 'plot range only' 
               ENDIF 
            ENDIF 
!                                                                       
!-------Set parameters                                                  
!                                                                       
         ELSEIF (str_comp (befehl, 'par', 2, lbef, 3) ) then 
            CALL do_fit_par (zeile, lp) 
!                                                                       
!-------Plot result                                                     
!                                                                       
         ELSEIF (str_comp (befehl, 'plot', 2, lbef, 4) ) then 
            CALL do_plot (.false.) 
!                                                                       
!------ Set scale                                                       
!                                                                       
         ELSEIF (str_comp (befehl, 'skal', 2, lbef, 4) ) then 
            CALL set_skal (zeile, lp) 
!                                                                       
!-------Run fit                                                         
!                                                                       
         ELSEIF (str_comp (befehl, 'run', 2, lbef, 3) ) then 
!
            lturn_off = .TRUE.
            IF (.not.sel_func) then 
               ier_num = - 25 
               ier_typ = ER_APPL 
            ELSE 
               IF(ex(iwin,iframe,1)==ex(iwin,iframe,2)) THEN
                  CALL set_skal(empty,0)
               ENDIF
               IF (lni (ikfit) ) then 
                  CALL do_fit_z 
               ELSE 
                  CALL do_fit_y 
               ENDIF 
               CALL get_extrema 
            ENDIF 
            lturn_off = .FALSE.
!                                                                  
!-------Test a recursive call to kuplot_loop                            
!                                                                       
         ELSEIF (str_comp (befehl, 'test', 2, lbef, 4) ) then 
            CALL theory_macro(1.0,  f, df, 1)
!                                                                  
!-------Save fit results                                                
!                                                                       
         ELSEIF (str_comp (befehl, 'save', 2, lbef, 4) ) then 
            CALL do_fit_save 
!                                                                       
!-------Show settings                                                   
!                                                                       
         ELSEIF (str_comp (befehl, 'show', 2, lbef, 4) ) then 
            CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
            IF (ier_num.eq.0) then 
               IF (ianz.eq.0) then 
                  CALL do_fit_info (output_io, .true., .true., .true.) 
               ELSEIF (ianz.eq.1) then 
                  CALL do_cap (cpara (1) ) 
                  flag (1) = cpara (1) (1:2) .eq.'GE' 
                  flag (2) = cpara (1) (1:2) .eq.'FI' 
                  flag (3) = cpara (1) (1:2) .eq.'PA' 
                  IF (flag (1) .or.flag (2) .or.flag (3) ) then 
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
         ELSEIF (str_comp (befehl, 'syst', 2, lbef, 4) ) then 
            cdummy = ' ' 
            IF (zeile.ne.' ') then 
               cdummy (1:lp) = zeile (1:lp) 
               CALL do_operating (cdummy, lp) 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
!                                                                       
!-------Set URF 'urf'                                                   
!                                                                       
         ELSEIF (str_comp (befehl, 'urf', 2, lbef, 3) ) then 
            CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
            IF (ier_num.eq.0) then 
               IF (ianz.eq.1) then 
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  IF (ier_num.eq.0) then 
                     urf = werte (1) 
                  ENDIF 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
            ENDIF 
!                                                                       
!     Waiting for user input                                            
!                                                                       
         ELSEIF (str_comp (befehl, 'wait', 3, lbef, 4) ) then 
            CALL do_input (zeile, lp) 
!                                                                       
!     Set weighting scheme                                              
!                                                                       
         ELSEIF (str_comp (befehl, 'wic', 3, lbef, 3) ) then 
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
      END SUBROUTINE do_fit                         
!*****7*****************************************************************
      SUBROUTINE do_fit_fkt (zeile, lp) 
!+                                                                      
!     Set theory function                                               
!-                                                                      
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
      USE kuplot_config 
      USE kuplot_mod 
      USE string_convert_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 8) 
!                                                                       
      CHARACTER ( * ) zeile 
      CHARACTER(1024) cpara (maxw) 
      CHARACTER(80) iname 
      REAL werte (maxw) 
      INTEGER lpara (maxw), lp 
      INTEGER ianz, iianz
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
!                                                                       
      IF (ianz.ge.1) then 
         ftyp = cpara (1) (1:4) 
         CALL do_cap (ftyp) 
         IF (ftyp (1:2) .eq.'GS') then 
            IF (ianz.ne.4.and.ianz.ne.5) then 
               ier_num = - 6 
               ier_typ = ER_COMM 
               RETURN 
            ELSE 
               iname = cpara (4) (1:lpara (4) ) 
               cpara (4) = '0.0' 
               lpara (4) = 3 
            ENDIF 
         ENDIF 
         IF (ftyp(1:2)=='MA' .AND. ianz==3) THEN
            CALL del_params (1, ianz, cpara, lpara, maxw) 
            iianz = 1
            CALL ber_params (iianz, cpara, lpara, werte, maxw) 
         ELSEIF (ianz.gt.1) then 
            CALL del_params (1, ianz, cpara, lpara, maxw) 
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         ELSE 
            ianz = 0 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
      IF (ier_num.ne.0) return 
!                                                                       
!------ Polynom (only 2D)                                               
!                                                                       
      IF (ftyp (1:2) .eq.'PO') then 
         IF (.not.lni (ikfit) ) then 
            CALL setup_poly (ianz, werte, maxw) 
         ELSE 
            ier_num = - 25 
            ier_typ = ER_APPL 
         ENDIF 
!                                                                       
!------ Scale + Background Polynom (only 2D)                            
!                                                                       
      ELSEIF (ftyp (1:2) .eq.'BA') then 
         IF (.not.lni (ikfit) ) then 
            CALL setup_backpoly (ianz, werte, maxw) 
         ELSE 
            ier_num = - 25 
            ier_typ = ER_APPL 
         ENDIF 
!                                                                       
!------ Chebyshev Polynom (only 2D)                                     
!                                                                       
      ELSEIF (ftyp (1:2) .eq.'CH') then 
         IF (.not.lni (ikfit) ) then 
            CALL setup_poly_cheb (ianz, werte, maxw) 
         ELSE 
            ier_num = - 25 
            ier_typ = ER_APPL 
         ENDIF 
!                                                                       
!------ GSAS profile function(s) (only 2D)                              
!                                                                       
      ELSEIF (ftyp (1:2) .eq.'GS') then 
         IF (.not.lni (ikfit) ) then 
            CALL setup_gsas (ianz, werte, maxw, iname) 
         ELSE 
            ier_num = - 25 
            ier_typ = ER_APPL 
         ENDIF 
!                                                                       
!------ Gaussian (xy and xyz data sets)                                 
!                                                                       
      ELSEIF (ftyp (1:2) .eq.'GA') then 
         IF (.not.lni (ikfit) ) then 
            CALL setup_gauss (ianz, werte, maxw) 
         ELSE 
            CALL setup_gauss_2d (ianz, werte, maxw) 
         ENDIF 
!                                                                       
!------ Lorenzian (only 2D)                                             
!                                                                       
      ELSEIF (ftyp (1:2) .eq.'LO') then 
         IF (.not.lni (ikfit) ) then 
            CALL setup_lor (ianz, werte, maxw) 
         ELSE 
            ier_num = - 25 
            ier_typ = ER_APPL 
         ENDIF 
!                                                                       
!------ Pseudo-Voigt (only 2D)                                          
!                                                                       
      ELSEIF (ftyp (1:2) .eq.'PS') then 
         IF (.not.lni (ikfit) ) then 
            CALL setup_psvgt (ianz, werte, maxw) 
         ELSE 
            ier_num = - 25 
            ier_typ = ER_APPL 
         ENDIF 
!                                                                       
!------ User defined function (xy and xyz data sets)                    
!                                                                       
      ELSEIF (ftyp (1:2) .eq.'FX') then 
         CALL setup_user (ianz, werte, maxw, cpara, lpara) 
!                                                                       
!------ User defined macro (xy and xyz data sets)                    
!                                                                       
      ELSEIF (ftyp (1:2) .eq.'MA') then 
         CALL setup_user_macro (ianz, werte, maxw, cpara, lpara) 
!                                                                       
      ELSE 
         ier_num = - 25 
         ier_typ = ER_APPL 
      ENDIF 
!                                                                       
      END SUBROUTINE do_fit_fkt                     
!*****7*****************************************************************
      SUBROUTINE do_fit_macro (zeile, lp) 
!+                                                                      
!     fitparameter und einstellungen als macro speichern                
!-                                                                      
      USE errlist_mod 
      USE get_params_mod
      USE prompt_mod 
      USE kuplot_config 
      USE kuplot_mod 
!
      USE build_name_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 10) 
!                                                                       
      CHARACTER ( * ) zeile 
      CHARACTER(1024) cpara (maxw) 
      REAL werte (maxw) 
      INTEGER lpara (maxw), lp 
      INTEGER ianz, i 
!                                                                       
      INTEGER len_str 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
      CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
      IF (ier_num.ne.0) return 
!                                                                       
      CALL oeffne (77, cpara (1) , 'unknown') 
      IF (ier_num.ne.0) return 
!                                                                       
      WRITE (output_io, 3000) cpara (1) (1:len_str (cpara (1) ) ) 
      DO i = 1, npara 
      WRITE (77, 1000) i 
      WRITE (77, 2000) i, pinc (i), p (i) 
      ENDDO 
      CLOSE (77) 
!                                                                       
 1000 FORMAT ('# parameter number ',i3) 
 2000 FORMAT ('par ',i3,',',g12.6,',',f2.0) 
 3000 FORMAT (' ---------- > Saving current parameters to macro file ', &
     &                   a,' ..')                                       
      END SUBROUTINE do_fit_macro                   
!*****7*****************************************************************
      SUBROUTINE do_fit_pmerk 
!+                                                                      
!     speichern der fitparameter                                        
!-                                                                      
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER i 
!                                                                       
      DO i = 1, npara 
      p_merk (i) = p (i) 
      pinc_merk (i) = pinc (i) 
      ENDDO 
      END SUBROUTINE do_fit_pmerk                   
!*****7*****************************************************************
      SUBROUTINE do_fit_prueck 
!+                                                                      
!     speichern der fitparameter                                        
!-                                                                      
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER i 
!                                                                       
      DO i = 1, npara 
      p (i) = p_merk (i) 
      pinc (i) = pinc_merk (i) 
      dp (i) = 0.0 
      ENDDO 
      END SUBROUTINE do_fit_prueck                  
!*****7*****************************************************************
      SUBROUTINE do_fit_save 
!+                                                                      
!     speichern der fitergebnisse                                       
!-                                                                      
      USE errlist_mod 
      USE prompt_mod 
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 4) 
!                                                                       
      CHARACTER(80) filname 
      CHARACTER(2) cdummy 
      REAL werte (maxw) 
      INTEGER ianz 
      INTEGER len_str 
!                                                                       
      filname = fname (ikfit)(1:MIN(60,LEN_TRIM(fname(ikfit))))
      filname = filname (1:len_str (filname) ) //'.erg' 
!                                                                       
      CALL oeffne (77, filname, 'unknown') 
      IF (ier_num.ne.0) return 
!                                                                       
      CALL do_fit_info (77, .true., .true., .true.) 
      WRITE (output_io, 1000) filname (1:len_str (filname) ) 
      CLOSE (77) 
!                                                                       
      WRITE (output_io, 2000) fname (ikcal) (1:len_str (fname (ikcal) ) &
      )                                                                 
      cdummy = fform (ikcal) 
      ianz = 0 
      CALL check_form (cdummy, ianz, werte, maxw) 
      CALL do_save (ikcal, fname (ikcal), fform (ikcal), ianz, werte,   &
      maxw, .false.)                                                    
!                                                                       
      WRITE (output_io, 3000) fname (ikdif) (1:len_str (fname (ikdif) ) &
      )                                                                 
      cdummy = fform (ikdif) 
      ianz = 0 
      CALL check_form (cdummy, ianz, werte, maxw) 
      CALL do_save (ikdif, fname (ikdif), fform (ikdif), ianz, werte,   &
      maxw, .false.)                                                    
!                                                                       
 1000 FORMAT (' ---------- > Saving results in file ',a,' ..') 
 2000 FORMAT (' ---------- > Saving calculated data in file ',a,' ..') 
 3000 FORMAT (' ---------- > Saving difference data in file ',a,' ..') 
      END SUBROUTINE do_fit_save                    
!*****7*****************************************************************
      SUBROUTINE do_fit_wichtung (zeile, lp) 
!+                                                                      
!     aendern der wichtung                                              
!-                                                                      
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
      USE kuplot_config 
      USE kuplot_mod 
      USE string_convert_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 2) 
!                                                                       
      CHARACTER ( * ) zeile 
      CHARACTER(1024) cpara (maxw) 
      INTEGER lpara (maxw), lp 
      INTEGER ianz 
      REAL werte (maxw) 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
!                                                                       
      IF (ianz.eq.1) then 
         CALL do_cap (cpara (1) ) 
         IF (cpara (1) (1:3) .ne.'LOG'.and.cpara (1) (1:3)              &
         .ne.'SQR'.and.cpara (1) (1:3) .ne.'ONE'.and.cpara (1) (1:3)    &
         .ne.'LIN'.and.cpara (1) (1:3) .ne.'SQA'.and.cpara (1) (1:3)    &
         .ne.'INV'.and.cpara (1) (1:3) .ne.'BCK'.and.cpara (1) (1:3)    &
         .ne.'ISQ'.and.cpara (1) (1:3) .ne.'DAT'                        &
                   ) then                  
            ier_num = - 27 
            ier_typ = ER_APPL 
         ELSE 
            wtyp = cpara (1) (1:3) 
            wval = 0.01 
         ENDIF 
      ELSEIF (ianz.eq.2) then 
         CALL do_cap (cpara (1) ) 
         IF (cpara (1) (1:3) .eq.'BCK') then 
            wtyp = cpara (1) (1:3) 
            cpara (1) = '(0)' 
            lpara (1) = 3 
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
            IF (ier_num.eq.0) then 
               wval = werte (2) 
            ENDIF 
         ELSE 
            ier_num = - 27 
            ier_typ = ER_APPL 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
      END SUBROUTINE do_fit_wichtung                
!*****7*****************************************************************
      SUBROUTINE do_fit_par (zeile, lp) 
!+                                                                      
!     aendern der fit-parameter                                         
!-                                                                      
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
      USE prompt_mod 
      USE param_mod
      USE kuplot_config 
      USE kuplot_mod 
      USE string_convert_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 3) 
!                                                                       
      CHARACTER ( * ) zeile 
      CHARACTER(1024) cpara (maxw) 
      REAL werte (maxw) 
      INTEGER lpara (maxw), lp 
      INTEGER ianz, ip 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
!                                                                       
      IF (ianz.eq.0) then 
         CALL do_fit_info (output_io, .false., .false., .true.) 
!                                                                       
      ELSEIF (ianz.eq.1) then 
         CALL do_cap (cpara (1) ) 
         IF (cpara (1) (1:2) .eq.'SA') then 
            CALL do_fit_pmerk 
         ELSEIF (cpara (1) (1:2) .eq.'LO') then 
            CALL do_fit_prueck 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
!                                                                       
      ELSEIF (ianz.eq.2.or.ianz.eq.3) then 
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         IF (ier_num.ne.0) return 
         ip = nint (werte (1) ) 
         IF (ip.gt.0.and.ip.le.MAXPARA) then 
            pinc (ip) = werte (2) 
            IF (ianz.eq.3) p (ip) = werte (3) 
            WRITE (output_io, 1000) ip, p (ip), pinc (ip) 
         ELSE 
            ier_num = - 26 
            ier_typ = ER_APPL 
         ENDIF 
!                                                                       
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!
      IF(ier_num == 0) THEN 
         DO ip=1,MAXPARA
            kupl_para(ip) = p(ip)
         ENDDO
      ENDIF
!                                                                       
 1000 FORMAT     (' ---------- > Setting parameter p(',i2,') = ',g12.6, &
     &                   '  pinc = ',f4.1)                              
      END SUBROUTINE do_fit_par                     
!*****7*****************************************************************
      SUBROUTINE show_fit_para (idout) 
!+                                                                      
!     anzeigen der fit-parameter                                        
!-                                                                      
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER(30) fitfkt, wictyp 
      INTEGER idout, lt1, lt2, lf, lfn, lw, ipkt 
!                                                                       
      INTEGER len_str 
!                                                                       
      IF (ftyp (1:2) .eq.'PO') then 
         fitfkt = 'Polynom' 
      ELSEIF (ftyp (1:2) .eq.'BA') then 
         fitfkt = 'Background Polynom' 
      ELSEIF (ftyp (1:2) .eq.'CH') then 
         fitfkt = 'Chebyshev Polynom' 
      ELSEIF (ftyp (1:2) .eq.'FX') then 
         fitfkt = 'f = '//fit_func (1:fit_lfunc) 
      ELSEIF (ftyp (1:2) .eq.'LO') then 
         fitfkt = 'Lorenzian' 
      ELSEIF (ftyp (1:2) .eq.'PS') then 
         fitfkt = 'Pseudo-Voigt' 
      ELSEIF (ftyp (1:2) .eq.'GA') then 
         IF (lni (ikfit) ) then 
            fitfkt = 'Gaussian (2D)' 
         ELSE 
            fitfkt = 'Gaussian (1D)' 
         ENDIF 
      ELSEIF (ftyp (1:2) .eq.'MA') then 
         fitfkt = 'f = '//fit_func (1:fit_lfunc) 
      ELSE 
         fitfkt = 'not defined' 
      ENDIF 
!                                                                       
      IF (wtyp (1:3) .eq.'LOG') then 
         wictyp = 'w(i) = log(i)' 
      ELSEIF (wtyp (1:3) .eq.'SQR') then 
         wictyp = 'w(i) = sqrt(i)' 
      ELSEIF (wtyp (1:3) .eq.'ONE') then 
         wictyp = 'w(i) = 1.0' 
      ELSEIF (wtyp (1:3) .eq.'LIN') then 
         wictyp = 'w(i) = i' 
      ELSEIF (wtyp (1:3) .eq.'SQA') then 
         wictyp = 'w(i) = i**2' 
      ELSEIF (wtyp (1:3) .eq.'INV') then 
         wictyp = 'w(i) = 1.0/i' 
      ELSEIF (wtyp (1:3) .eq.'ISQ') then 
         wictyp = 'w(i) = 1.0/sqrt(i)' 
      ELSEIF (wtyp (1:3) .eq.'DAT'.and..not.lni (ikfit) ) then 
         wictyp = 'w(i) = 1/(dy(i))^2 from data ' 
      ELSEIF (wtyp (1:3) .eq.'BCK') then 
         wictyp = 'w(i) = exp(-WVAL*(Fobs-Fcalc))' 
      ELSE 
         wictyp = 'unknown' 
      ENDIF 
!                                                                       
      lt1 = max (1, len_str (titel (iwin, iframe, 1) ) ) 
      lt2 = max (1, len_str (titel (iwin, iframe, 2) ) ) 
      lf = len_str (fitfkt) 
      lfn = len_str (fname (ikfit) ) 
      lw = len_str (wictyp) 
!                                                                       
      IF (lni (ikfit) ) then 
         ipkt = nx (ikfit) * ny (ikfit) 
      ELSE 
         ipkt = len (ikfit) 
      ENDIF 
!                                                                       
      WRITE (idout, 1000) titel (iwin, iframe, 1) (1:lt1), titel (iwin, &
      iframe, 2) (1:lt2), fitfkt (1:lf), fname (ikfit) (1:lfn), ipkt,   &
      npara, urf, ncycle, fit_ifen, wictyp (1:lw)                       
!                                                                       
 1000 FORMAT (1x,'General fit parameter settings : ',/,                 &
     &        3x,'Title            : ',a,/11x,a,/,                      &
     &        3x,'Fit function     : ',a,/                              &
     &        3x,'Data file name   : ',a,/                              &
     &        3x,'# of data pts.   : ',i6,/                             &
     &        3x,'# of parameters  : ',i6,/                             &
     &        3x,'Urf              : ',f9.4,/                           &
     &        3x,'Max. cycle       : ',i6,/                             &
     &        3x,'MFEN for maxima  : ',i6,/                             &
     &        3x,'Weighting scheme : ',a,/)                             
      END SUBROUTINE show_fit_para                  
!*****7*****************************************************************
      SUBROUTINE show_fit_erg (idout) 
!+                                                                      
!     anzeigen der fit-ergebnisse                                       
!-                                                                      
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER idout, i, j 
      LOGICAL kor 
!                                                                       
      WRITE (idout, 1040) zalt, zwert, zdif, fend, r4, re 
      WRITE (idout, 1050) 
      kor = .false. 
      DO i = 2, npara 
      DO j = 1, i - 1 
      IF (abs (cl (i, j) ) .gt.0.8) then 
         WRITE (idout, 1070) i, j, cl (i, j) 
         kor = .true. 
      ENDIF 
      ENDDO 
      ENDDO 
      IF (.not.kor) write (idout, 1060) 
      WRITE (idout, * ) ' ' 
!                                                                       
 1040 FORMAT (' Information about the fit : ',/,                        &
     &        3x,'Sum n-1    : ',g12.6,12x,' Sum n   : ',g12.6/,        &
     &        3x,'Difference : ',g12.6,/,                               &
     &        3x,'Urf final  : ',g12.6,/,                               &
     &        3x,'R4 value   : ',g12.6,12x,' R exp   : ',g12.6/)        
 1050 FORMAT (' Correlations larger than 0.8 :') 
 1060 FORMAT (3x,'** none **') 
 1070 FORMAT (3x,'Between p(',i2,') - p(',i2,') : ',f6.3) 
!                                                                       
      END SUBROUTINE show_fit_erg                   
!*****7*****************************************************************
      SUBROUTINE write_fit 
!+                                                                      
!     kupl.fit schreiben fuer textframe                                 
!-                                                                      
      USE errlist_mod 
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER i, j 
      LOGICAL kor 
!                                                                       
      CALL oeffne (22, 'kupl.fit', 'unknown') 
      IF (ier_num.ne.0) return 
!                                                                       
      WRITE (22, 1000) ftyp, r4 * 100., re * 100. 
      kor = .false. 
      DO i = 2, npara 
      DO j = 1, i - 1 
      IF (abs (cl (i, j) ) .gt.0.8) then 
         WRITE (22, 1070) i, j, cl (i, j) 
         kor = .true. 
      ENDIF 
      ENDDO 
      ENDDO 
      IF (.not.kor) write (22, 1060) 
      WRITE (22, 1100) 
      DO i = 1, npara 
      IF (pinc (i) .ne.1) then 
         WRITE (22, 1200) i, p (i) 
      ELSE 
         WRITE (22, 1210) i, p (i), dp (i) 
      ENDIF 
      ENDDO 
      CLOSE (22) 
!                                                                       
 1000 FORMAT (/1x,'F i t  -  r e s u l t s',/,                          &
     &         1x,'--------------------------------------',/,           &
     &         1x,'Fit function : ',a4,/,                               &
     &         1x,'R value      : ',f5.1,' %',/                         &
     &         1x,'Rexp value   : ',f5.1,' %',/                         &
     &         1x,'--------------------------------------',/,           &
     &         1x,'Correlations > 0.8 : ',/)                            
 1060 FORMAT ( 1x,'** none **') 
 1070 FORMAT ( 1x,'betw. p(',i2,') - p(',i2,') : ',f6.3) 
 1100 FORMAT ( 1x,'--------------------------------------',/,           &
     &         1x,'Resulting parameters: ',/)                           
 1200 FORMAT ( 1x,'p(',i2,') = ',g12.6,' fixed') 
 1210 FORMAT ( 1x,'p(',i2,') = ',g12.6,' +- ',g12.6) 
!                                                                       
      CLOSE (22) 
!                                                                       
      END SUBROUTINE write_fit                      
!*****7*****************************************************************
      SUBROUTINE do_fit_info (idout, f_se, f_er, f_pa) 
!+                                                                      
!     ausgabe von fitinformationen                                      
!-                                                                      
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER idout 
      LOGICAL f_se, f_er, f_pa 
!                                                                       
      IF (f_se) call show_fit_para (idout) 
      IF (f_er) call show_fit_erg (idout) 
!                                                                       
      IF (f_pa) then 
         IF (ftyp (1:2) .eq.'PO') then 
            CALL show_poly (idout) 
         ELSEIF (ftyp (1:2) .eq.'BA') then 
            CALL show_backpoly (idout) 
         ELSEIF (ftyp (1:2) .eq.'GS') then 
            CALL show_gsas (idout) 
         ELSEIF (ftyp (1:2) .eq.'CH') then 
            CALL show_poly_cheb (idout) 
         ELSEIF (ftyp (1:2) .eq.'GA') then 
            IF (.not.lni (ikfit) ) then 
               CALL show_gauss (idout) 
            ELSE 
               CALL show_gauss_2d (idout) 
            ENDIF 
         ELSEIF (ftyp (1:2) .eq.'LO') then 
            CALL show_lor (idout) 
         ELSEIF (ftyp (1:2) .eq.'PS') then 
            CALL show_psvgt (idout) 
         ELSEIF (ftyp (1:2) .eq.'FX') then 
            CALL show_user (idout) 
         ELSEIF (ftyp (1:2) .eq.'MA') then 
            CALL show_user_macro (idout) 
         ENDIF 
      ENDIF 
!                                                                       
      END SUBROUTINE do_fit_info                    
!*****7*****************************************************************
SUBROUTINE do_fit_y 
!+                                                                      
!     der eigentliche fit fuer xy-files                                 
!-                                                                      
USE prompt_mod 
USE kuplot_config 
USE kuplot_mod 
!                                                                       
IMPLICIT none 
!                                                                       
CHARACTER(LEN=60) :: filname 
REAL              :: xx, f, df (maxpara) 
INTEGER           :: i, ii, jj, kk 
INTEGER           :: len_str 
!                                                                       
CALL wichtung (y) 
      IF (ncycle.gt.0) call fit_kupl (y) 
!                                                                       
      ii = offxy (ikfit - 1) 
      jj = offxy (ikcal - 1) 
      kk = offxy (ikdif - 1) 
!                                                                       
      DO i = 1, len (ikfit) 
      xx = x (ii + i) 
      CALL kupl_theory (xx, f, df, - i) 
      x (jj + i) = xx 
      y (jj + i) = f 
      dx (jj + i) = 0.0 
      dy (jj + i) = 0.0 
      x (kk + i) = xx 
      y (kk + i) = y (ii + i) - f 
      dx (kk + i) = 0.0 
      dy (kk + i) = 0.0 
      ENDDO 
!                                                                       
      len (ikcal) = len (ikfit) 
      len (ikdif) = len (ikfit) 
      fform (ikcal) = fform (ikfit) 
      fform (ikdif) = fform (ikfit) 
      filname = fname (ikfit)(1:MIN(60,LEN_TRIM(fname(ikfit)))) 
      fname (ikcal) = filname (1:len_str (filname) ) //'.fit' 
      fname (ikdif) = filname (1:len_str (filname) ) //'.dif' 
      CALL get_extrema 
!                                                                       
      CALL do_fit_info (output_io, .false., .false., .true.) 
      END SUBROUTINE do_fit_y                       
!*****7*****************************************************************
      SUBROUTINE do_fit_z 
!+                                                                      
!     der eigentliche fit fuer xyz-files                                
!-                                                                      
      USE prompt_mod 
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER(60) filname 
      REAL xx, f, df (maxpara) 
      INTEGER i, iii 
      INTEGER len_str 
!                                                                       
      CALL wichtung (z) 
      IF (ncycle.gt.0) call fit_kupl (z) 
!                                                                       
      DO i = 1, nx (ikfit) * ny (ikfit) 
      xx = REAL(i) 
      CALL kupl_theory (xx, f, df, - i) 
      z (offz (ikcal - 1) + i) = f 
      IF (z (offz (ikfit - 1) + i) .ne. - 9999) then 
         z (offz (ikdif - 1) + i) = z (offz (ikfit - 1) + i) - f 
      ELSE 
         z (offz (ikdif - 1) + i) = - 9999.0 
      ENDIF 
      ENDDO 
      DO iii = 1, nx (ikfit) 
      x (offxy (ikcal - 1) + iii) = x (offxy (ikfit - 1) + iii) 
      x (offxy (ikdif - 1) + iii) = x (offxy (ikfit - 1) + iii) 
      ENDDO 
      DO iii = 1, ny (ikfit) 
      y (offxy (ikcal - 1) + iii) = y (offxy (ikfit - 1) + iii) 
      y (offxy (ikdif - 1) + iii) = y (offxy (ikfit - 1) + iii) 
      ENDDO 
!                                                                       
      lni (ikcal) = .true. 
      lni (ikdif) = .true. 
      len (ikcal) = len (ikfit) 
      len (ikdif) = len (ikfit) 
      nx (ikcal) = nx (ikfit) 
      ny (ikcal) = ny (ikfit) 
      nx (ikdif) = nx (ikfit) 
      ny (ikdif) = ny (ikfit) 
      fform (ikcal) = fform (ikfit) 
      fform (ikdif) = fform (ikfit) 
!                                                                       
      filname = fname (ikfit)(1:MIN(60,LEN_TRIM(fname(ikfit)))) 
      fname (ikcal) = filname (1:len_str (filname) ) //'.fit' 
      fname (ikdif) = filname (1:len_str (filname) ) //'.dif' 
      CALL get_extrema 
!                                                                       
      CALL do_fit_info (output_io, .false., .false., .true.) 
      END SUBROUTINE do_fit_z                       
!*****7*****************************************************************
      SUBROUTINE wichtung (a) 
!+                                                                      
!     Calculation of weights. Values outside plotting range are         
!     set to zero if frall is .false.                                   
!-                                                                      
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      REAL a (maxarray) 
      INTEGER i, ii 
!                                                                       
      REAL calc_wic 
!                                                                       
      IF (lni (ikfit) ) then 
         ii = offz (ikfit - 1) 
         DO i = 1, nx (ikfit) * ny (ikfit) 
         w (ii + i) = calc_wic (a (ii + i), dy (ii + i) ) 
         ENDDO 
      ELSE 
         ii = offxy (ikfit - 1) 
         DO i = 1, len (ikfit) 
         IF (frall) then 
            w (ii + i) = calc_wic (a (ii + i), dy (ii + i) ) 
         ELSE 
            IF (x (ii + i) .lt.ex (iwin, iframe, 1) .or.x (ii + i)      &
            .gt.ex (iwin, iframe, 2) ) then                             
               w (ii + i) = 0.0 
            ELSE 
               w (ii + i) = calc_wic (a (ii + i), dy (ii + i) ) 
            ENDIF 
         ENDIF 
         ENDDO 
      ENDIF 
!                                                                       
!                                                                       
      END SUBROUTINE wichtung                       
!*****7*****************************************************************
      REAL function calc_wic (val, sig) 
!+                                                                      
!     Calculation of weights.                                           
!-                                                                      
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      REAL val, aval, sig, wic 
!                                                                       
      wic = 0.0 
      IF (val.ne. - 9999.0) then 
         aval = abs (val) 
         IF (wtyp (1:3) .eq.'LOG') then 
            IF (aval.gt.0.0) wic = log (aval) 
         ELSEIF (wtyp (1:3) .eq.'SQR') then 
            wic = sqrt (aval) 
         ELSEIF (wtyp (1:3) .eq.'ONE') then 
            wic = 1.0 
         ELSEIF (wtyp (1:3) .eq.'LIN') then 
            wic = aval 
         ELSEIF (wtyp (1:3) .eq.'SQA') then 
            wic = aval**2 
         ELSEIF (wtyp (1:3) .eq.'INV') then 
            IF (aval.ne.0.0) wic = 1.0 / aval 
         ELSEIF (wtyp (1:3) .eq.'ISQ') then 
            wic = 1.0 / sqrt (aval) 
         ELSEIF (wtyp (1:3) .eq.'DAT'.and..not.lni (ikfit) ) then 
            wic = 1.0 / (sig *sig)
         ENDIF 
      ENDIF 
!                                                                       
      calc_wic = wic 
!                                                                       
      END FUNCTION calc_wic                         
!*****7*****************************************************************
      SUBROUTINE kupl_theory (xx, f, df, iwert) 
!                                                                       
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      REAL xx, f, df (maxpara) 
      INTEGER iwert 
!                                                                       
      IF (ftyp (1:2) .eq.'PO') then 
         CALL theory_poly (xx, f, df, iwert) 
      ELSEIF (ftyp (1:2) .eq.'BA') then 
         CALL theory_backpoly (xx, f, df, iwert) 
      ELSEIF (ftyp (1:2) .eq.'CH') then 
         CALL theory_poly_cheb (xx, f, df, iwert) 
      ELSEIF (ftyp (1:2) .eq.'GS') then 
         CALL theory_gsas (xx, f, df, iwert) 
      ELSEIF (ftyp (1:2) .eq.'LO') then 
         CALL theory_lor (xx, f, df, iwert) 
      ELSEIF (ftyp (1:2) .eq.'PS') then 
         CALL theory_psvgt (xx, f, df, iwert) 
      ELSEIF (ftyp (1:2) .eq.'FX') then 
         CALL theory_user (xx, f, df, iwert) 
      ELSEIF (ftyp (1:2) .eq.'GA') then 
         IF (.not.lni (ikfit) ) then 
            CALL theory_gauss (xx, f, df, iwert) 
         ELSE 
            CALL theory_gauss_2d (xx, f, df, iwert) 
         ENDIF 
      ELSEIF (ftyp (1:2) .eq.'MA') then 
         CALL theory_macro (xx, f, df, iwert) 
      ENDIF 
!                                                                       
      END SUBROUTINE kupl_theory                    
!*****7*****************************************************************
!       User defined fit function                                       
!*****7*****************************************************************
      SUBROUTINE setup_user (ianz, werte, maxw, cpara, lpara) 
!                                                                       
      USE  berechne_mod
      USE errlist_mod 
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
!                                                                       
      CHARACTER ( * ) cpara (maxw) 
      CHARACTER(1024) cdummy 
      REAL dummy, werte (maxw) 
      INTEGER lpara (maxw) 
      INTEGER ianz, ip , length
!                                                                       
!                                                                       
      IF (ianz.eq.2) then 
         ip = nint (werte (1) ) 
         IF (ip.lt.1.or.ip.gt.maxpara) then 
            ier_num = - 31 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
!                                                                       
         npara = ip 
         fit_func = cpara (2) (1:lpara (2) ) 
         fit_lfunc = lpara (2) 
         cdummy = '('//fit_func (1:fit_lfunc) //')' 
         length = fit_lfunc + 2
         dummy = berechne (cdummy, length)
!                                                                       
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
      END SUBROUTINE setup_user                     
!*****7*****************************************************************
!       User defined fit macro                                       
!*****7*****************************************************************
      SUBROUTINE setup_user_macro (ianz, werte, maxw, cpara, lpara) 
!                                                                       
      USE  berechne_mod
      USE errlist_mod 
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER , INTENT(IN) :: maxw 
      INTEGER , INTENT(IN) :: ianz 
      REAL    , DIMENSION(MAXW), INTENT(IN) ::  werte (maxw) 
!                                                                       
      CHARACTER (LEN=*), DIMENSION(MAXW), INTENT(IN) :: cpara
      INTEGER          , DIMENSION(MAXW), INTENT(IN) :: lpara
      INTEGER :: ip 
      INTEGER :: m,i,iiw,iix
REAL, DIMENSION(maxpara) :: df
REAL                   :: f
      REAL :: xx
!                                                                       
!                                                                       
      IF (ianz.eq.2) then 
         ip = nint (werte (1) ) 
         IF (ip.lt.1.or.ip.gt.maxpara) then 
            ier_num = - 31 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
!                                                                       
         npara = ip 
         fit_func = cpara (2) (1:lpara (2) ) 
         fit_lfunc = lpara (2) 
!        cdummy = '('//fit_func (1:fit_lfunc) //')' 
!        length = fit_lfunc + 2
!        dummy = berechne (cdummy, length)
         IF (lni (ikfit) ) then 
            iiw = offz (ikfit - 1) 
            iix = offxy (ikfit - 1) 
            m = nx (ikfit) * ny (ikfit) 
         ELSE 
            iiw = offxy (ikfit - 1) 
            iix = offxy (ikfit - 1) 
            m = len (ikfit) 
         ENDIF 
         DO i = 1, 1  !!! m 
            xx = x (iix + i) 
            CALL theory_macro (xx, f, df, i) 
            IF(ier_num/=0) RETURN
         ENDDO
!                                                                       
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
      END SUBROUTINE setup_user_macro
!*****7*****************************************************************
      SUBROUTINE show_user (idout) 
!                                                                       
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER idout, i 
!                                                                       
      WRITE (idout, 1000) fit_func (1:fit_lfunc) 
      DO i = 1, npara 
      WRITE (idout, 1100) i, p (i), dp (i), pinc (i) 
      ENDDO 
      WRITE (idout, * ) ' ' 
!                                                                       
 1000 FORMAT (1x,'Fit function : F = ',a,/) 
 1100 FORMAT (3x,'p(',i2,') : ',g12.6,' +- ',g12.6,4x,'pinc : ',f2.0) 
!                                                                       
      END SUBROUTINE show_user                      
!*****7*****************************************************************
      SUBROUTINE show_user_macro (idout) 
!                                                                       
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER idout, i 
!                                                                       
      WRITE (idout, 1000) fit_func (1:fit_lfunc) 
      DO i = 1, npara 
      WRITE (idout, 1100) i, p (i), dp (i), pinc (i) 
      ENDDO 
      WRITE (idout, * ) ' ' 
!                                                                       
 1000 FORMAT (1x,'Fit function : F = ',a,/) 
 1100 FORMAT (3x,'p(',i2,') : ',g12.6,' +- ',g12.6,4x,'pinc : ',f2.0) 
!                                                                       
      END SUBROUTINE show_user_macro                      
!***********************************************************************
      SUBROUTINE theory_user (xx, f, df, i) 
!                                                                       
      USE  berechne_mod
      USE param_mod 
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER(1024) cdummy 
      REAL xx, f, df (maxpara) 
      REAL pb, h, err 
      INTEGER i, ix, iy, ind 
      INTEGER :: length
!                                                                       
      REAL dfridr 
!                                                                       
      DO ind = 1, maxpara 
      df (ind) = 0.0 
      ENDDO 
!                                                                       
!-------x und y bestimmen                                               
!                                                                       
      IF (lni (ikfit) ) then 
         ix = (abs (i) - 1) / ny (ikfit) + 1 
         iy = abs (i) - (ix - 1) * ny (ikfit) 
         rpara (0) = x (offxy (ikfit - 1) + ix) 
         rpara (1) = y (offxy (ikfit - 1) + iy) 
      ELSE 
         rpara (0) = xx 
      ENDIF 
!                                                                       
      cdummy = '('//fit_func (1:fit_lfunc) //')' 
      length = fit_lfunc + 2
      f = berechne (cdummy, length)
!                                                                       
!-------ableitungen                                                     
!                                                                       
      IF (i.gt.0) then 
         DO ind = 1, npara 
         IF (pinc (ind) .ne.0) then 
            h = 0.1 * abs (p (ind) ) 
            IF (h.eq.0) h = 0.5 
            pb = p (ind) 
            np1 = ind 
            df (ind) = dfridr (p (ind), h, err) 
            p (ind) = pb 
         ENDIF 
         ENDDO 
      ENDIF 
!                                                                       
      END SUBROUTINE theory_user                    
!*****7*****************************************************************
      REAL function func (xx) 
!                                                                       
      USE  berechne_mod
      USE param_mod 
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER(LEN=1024) :: cdummy 
      INTEGER :: length
      REAL xx
!                                                                       
      p (np1) = xx 
      cdummy = '('//fit_func (1:fit_lfunc) //')' 
      length = fit_lfunc + 2
      func = berechne (cdummy, length)
!                                                                       
      END FUNCTION func                             
!*****7*****************************************************************
!       GSAS profile functions                                          
!*****7*****************************************************************
      SUBROUTINE show_gsas (idout) 
!                                                                       
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER(40) cout (maxpara) 
      INTEGER idout, i, j, k 
!                                                                       
      DO j = 1, np2 
      cout ( (j - 1) * np3 + 1) = 'Peak position' 
      cout ( (j - 1) * np3 + 2) = 'Intensity' 
      DO i = 3, npara - 2 
      WRITE (cout ( (j - 1) * np3 + i), 2000) i - 2 
      ENDDO 
!                                                                       
      WRITE (idout, 1000) j, np1 
      DO i = 1, np3 
      k = (j - 1) * np3 + i 
      WRITE (idout, 1100) k, cout (k), p (k), dp (k), pinc (k) 
      ENDDO 
      WRITE (idout, * ) ' ' 
      ENDDO 
!                                                                       
      WRITE (idout, 1200) 
      cout (npara - 1) = 'Background const.' 
      cout (npara) = 'Background slope' 
      DO i = npara - 1, npara 
      WRITE (idout, 1100) i, cout (i), p (i), dp (i), pinc (i) 
      ENDDO 
!                                                                       
 1000 FORMAT   (1x,'GSAS profile ',i2,': Type ',i2,/) 
 1100 FORMAT   (3x,'p(',i2,') : ',a20,' : ',                            &
     &                       g12.6,' +- ',g12.6,4x,'pinc : ',f2.0)      
 1200 FORMAT     (1x,'Global parameters:',/) 
 2000 FORMAT     ('Profile Coeff. ',i2) 
!                                                                       
      END SUBROUTINE show_gsas                      
!*****7*****************************************************************
      SUBROUTINE setup_gsas (ianz, werte, maxw, iname) 
!                                                                       
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw, maxmax 
      PARAMETER (maxmax = 50) 
!                                                                       
      CHARACTER ( * ) iname 
      REAL werte (maxw) 
      INTEGER ianz 
!                                                                       
      REAL pcoff (maxpara) 
      REAL wmax (maxmax) 
      REAL stheta, dspace, xpeak, inten, delt 
      REAL x1, x2, y1, y2 
      INTEGER ixm (maxmax) 
      INTEGER ncoff, itype, ibank, ipeaks, ima, i, j 
!                                                                       
      itype = nint (werte (1) ) 
      ibank = nint (werte (2) ) 
      IF (ianz.eq.4) then 
         ipeaks = nint (werte (4) ) 
      ELSE 
         ipeaks = 1 
      ENDIF 
!                                                                       
      ifen = fit_ifen 
      CALL do_fmax_xy (ikfit, wmax, ixm, maxmax, ima) 
      IF (ima.ge.1) then 
         xpeak = x (offxy (ikfit - 1) + ixm (1) ) 
      ELSE 
         xpeak = x (offxy (ikfit - 1) + len (ikfit) / 2) 
      ENDIF 
!                                                                       
      CALL read_prof (iname, ibank, itype, pcoff, ncoff, stheta, dspace,&
      xpeak)                                                            
      CALL cnvptp1 (itype, pcoff, ncoff, p, dspace, stheta, maxpara) 
!                                                                       
      inten = 0.0 
      DO i = 1, len (ikfit) - 1 
      delt = x (offxy (ikfit - 1) + i + 1) - x (offxy (ikfit - 1)       &
      + i)                                                              
      inten = inten + y (offxy (ikfit - 1) + i) * delt 
      ENDDO 
!                                                                       
      np1 = itype 
      np2 = ipeaks 
      np3 = ncoff 
      npara = ncoff * ipeaks + 2 
!                                                                       
      DO i = 1, ipeaks 
      IF (ima.ge.i) then 
         p ( (i - 1) * np3 + 1) = x (offxy (ikfit - 1) + ixm (i) ) 
      ELSE 
         p ( (i - 1) * np3 + 1) = x (offxy (ikfit - 1) + len (ikfit)    &
         / 2)                                                           
      ENDIF 
      p ( (i - 1) * np3 + 2) = inten 
      DO j = 3, ncoff 
      p ( (i - 1) * np3 + j) = p (j) 
      ENDDO 
      ENDDO 
!                                                                       
      y1 = y (offxy (ikfit - 1) + 1) 
      y2 = y (offxy (ikfit - 1) + len (ikfit) ) 
      x1 = x (offxy (ikfit - 1) + 1) 
      x2 = x (offxy (ikfit - 1) + len (ikfit) ) 
      p (npara - 1) = y1 
      p (npara) = (y2 - y1) / (x2 - x1) 
!                                                                       
      END SUBROUTINE setup_gsas                     
!*****7*****************************************************************
      SUBROUTINE read_prof (iname, ibank, itype, pcoff, ncoff, stheta,  &
      dspace, tof)                                                      
!                                                                       
      USE debug_mod 
      USE errlist_mod 
      USE prompt_mod 
      USE kuplot_config 
      USE trig_degree_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      CHARACTER ( * ) iname 
      REAL pcoff (maxpara) 
      REAL stheta, dspace, tof 
      INTEGER ibank, itype, ncoff 
!                                                                       
      CHARACTER(80) line 
      CHARACTER(6) search 
      CHARACTER(5) key 
      REAL tmp (4) 
      REAL difc, difa, zero, tth, l2 
      REAL secondterm 
      INTEGER k, ib, itmp, ll 
!                                                                       
      INTEGER len_str 
!     REAL sind 
!                                                                       
      CALL oeffne (12, iname, 'old') 
      IF (ier_num.ne.0) return 
!                                                                       
      WRITE (key, 1000) itype 
 1000 FORMAT    ('PRCF',i1) 
!                                                                       
!------ Read instrument parameter file information                      
!                                                                       
   20 CONTINUE 
      READ (12, '(a)', end = 40) line 
      ll = len_str (line) 
      IF (ll.eq.0) goto 20 
!                                                                       
      READ (line (1:ll) , '(4x,i2,a6)', err = 20) ib, search 
      IF (ib.ne.ibank) goto 20 
!                                                                       
      IF (search.eq.' ICONS') then 
         READ (line (13:ll), *, err = 998, end = 998) difc, difa, zero 
      ELSEIF (search.eq.'BNKPAR') then 
         READ (line (13:ll), *, err = 998) l2, tth 
      ELSEIF (search (1:5) .eq.key) then 
         IF (search (6:6) .eq.' ') then 
            READ (line (13:ll), *, err = 998) itmp, ncoff 
         ELSE 
            READ (search (6:6), * ) itmp 
            READ (line (13:ll), *, err = 998) tmp 
            DO k = 1, 4 
            pcoff ( (itmp - 1) * 4 + k) = tmp (k) 
            ENDDO 
         ENDIF 
      ENDIF 
      GOTO 20 
!                                                                       
   40 CONTINUE 
      CLOSE (12) 
!                                                                       
      stheta = sind (0.5 * tth) 
      IF (difa.ne.0.0) then 
         secondterm = sqrt (4.0 * difa * (tof - zero) + difc * difc) 
         dspace = ( - difc + secondterm) / 2. / difa 
      ELSE 
         dspace = (tof - zero) / difc 
      ENDIF 
!                                                                       
      RETURN 
!                                                                       
      ier_num = - 46 
      ier_typ = ER_APPL 
      CLOSE (12) 
      RETURN 
!                                                                       
  998 CONTINUE 
      ier_num = - 47 
      ier_typ = ER_APPL 
      CLOSE (12) 
!                                                                       
      END SUBROUTINE read_prof                      
!*****7*****************************************************************
      SUBROUTINE theory_gsas (xx, f, df, i) 
!                                                                       
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER(4) htype 
      REAL deriv (maxpara) 
      REAL xx, f, df (maxpara) 
      REAL pp (maxpara) 
      REAL tth, c, dtof, delt 
      REAL xnext 
      INTEGER ind, ip, i, j, k, ptype 
!                                                                       
      DO ind = 1, maxpara 
      df (ind) = 0.0 
      deriv (ind) = 0.0 
      ENDDO 
!                                                                       
      j = abs (i) 
      IF (j.ge.len (ikfit) ) return 
!                                                                       
      tth = 0. 
      ptype = np1 
      htype = 'PNT ' 
!                                                                       
      xnext = x (offxy (ikfit - 1) + j + 1) 
      delt = xnext - xx 
      f = p (npara - 1) + p (npara) * xx 
!                                                                       
      DO ip = 1, np2 
      DO k = 1, np3 
      pp (k) = p ( (ip - 1) * np3 + k) 
      ENDDO 
      dtof = (xx - pp (1) ) / 1000. 
      CALL prpcalc (htype, tth, ptype, pp, dtof, c, deriv, MAXPARA) 
      f = f + pp (2) * c 
!                                                                       
      DO k = 1, np3 
      p ( (ip - 1) * np3 + k) = pp (k) 
      ENDDO 
!                                                                       
      IF (i.ge.1) then 
!tep   write(*,'(5g15.6)') deriv(1),deriv(2),deriv(3),deriv(4),deriv(5) 
         IF (pinc ( (ip - 1) * np3 + 1) .ne.0.0) df ( (ip - 1) * np3 +  &
         1) = pp (2) * deriv (1) / 1000.                                
         IF (pinc ( (ip - 1) * np3 + 2) .ne.0.0) df ( (ip - 1) * np3 +  &
         2) = deriv (2)                                                 
         DO k = 3, npara - 2 
         IF (pinc ( (ip - 1) * np3 + k) .ne.0.0) df ( (ip - 1) * np3 +  &
         k) = deriv (k) * pp (2)                                        
         ENDDO 
      ENDIF 
      ENDDO 
!                                                                       
      IF (i.ge.1) then 
         IF (pinc (npara - 1) .ne.0.0) df (npara - 1) = 1.0 
         IF (pinc (npara) .ne.0.0) df (npara) = xx 
      ENDIF 
!                                                                       
      END SUBROUTINE theory_gsas                    
!*****7*****************************************************************
!       Lorenzian                                                       
!*****7*****************************************************************
      SUBROUTINE show_lor (idout) 
!                                                                       
      USE wink_mod
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER idout, i, iii 
      REAL zz, sint, ds0, ds2, ds3, dsint 
!                                                                       
      WRITE (idout, 1000) np1 
      WRITE (idout, 1100) 1, p (1), dp (1), pinc (1) 
      WRITE (idout, 1200) 2, p (2), dp (2), pinc (2) 
      DO i = 1, np1 
      WRITE (idout, 1300) i 
      iii = 2 + (i - 1) * 4 
      WRITE (idout, 1400) iii + 1, p (iii + 1), dp (iii + 1), pinc (iii &
      + 1)                                                              
      WRITE (idout, 1500) iii + 2, p (iii + 2), dp (iii + 2), pinc (iii &
      + 2)                                                              
      WRITE (idout, 1600) iii + 3, p (iii + 3), dp (iii + 3), pinc (iii &
      + 3)                                                              
      WRITE (idout, 1700) iii + 4, p (iii + 4), dp (iii + 4), pinc (iii &
      + 4)                                                              
!---------integral berechnen                                            
      zz = p (iii + 3) / p (iii + 4) + p (iii + 3) * p (iii + 4) 
      sint = p (iii + 1) * REAL(zpi) * zz 
      ds0 = REAL(zpi) * zz 
      ds2 = p (iii + 1) * REAL(zpi) * (1.0 / p (iii + 4) + p (iii + 4) ) 
      ds3 = p (iii + 1) * REAL(zpi) * ( - p (iii + 3) / p (iii + 4) **2 + p ( &
      iii + 3) )                                                        
      dsint = dp (iii + 1) * ds0 + dp (iii + 3) * ds2 + dp (iii + 4)    &
      * ds3                                                             
      WRITE (idout, 1800) sint, dsint 
      ENDDO 
      WRITE (idout, * ) ' ' 
!                                                                       
 1000 FORMAT     (1x,'Fitted',i3,' Lorenzian(s) : '/) 
 1100 FORMAT     (3x,'p(',i2,') : backgr. 1 : ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1200 FORMAT     (3x,'p(',i2,') : backgr. 2 : ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1300 FORMAT     (/,1x,'Lorenzian : ',i3,/) 
 1400 FORMAT     (3x,'p(',i2,') : peak      : ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1500 FORMAT     (3x,'p(',i2,') : position  : ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1600 FORMAT     (3x,'p(',i2,') : fwhm      : ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1700 FORMAT     (3x,'p(',i2,') : asymmetry : ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1800 FORMAT     (3x,'        integral  : ',g12.6,' +- ',g12.6) 
!                                                                       
      END SUBROUTINE show_lor                       
!*****7*****************************************************************
      SUBROUTINE setup_lor (ianz, werte, maxw) 
!                                                                       
      USE errlist_mod 
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxmax 
      PARAMETER (maxmax = 50) 
!                                                                       
      INTEGER maxw 
      REAL werte (maxw) 
      REAL wmax (maxmax) 
      INTEGER ixm (maxmax) 
      INTEGER ianz, ii, jj, ima, i 
!                                                                       
      IF (ianz.eq.0) then 
         npara = 6 
         np1 = 1 
      ELSEIF (ianz.eq.1) then 
         ii = nint (werte (1) ) 
         IF (ii.gt.0.and. (2 + 4 * ii) .le.maxpara) then 
            np1 = ii 
            npara = 2 + 4 * np1 
         ELSE 
            ier_num = - 31 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
         RETURN 
      ENDIF 
!                                                                       
      ii = offxy (ikfit - 1) + 1 
      jj = offxy (ikfit - 1) + len (ikfit) 
!                                                                       
      p (1) = y (ii) 
      pinc (1) = 1.0 
      p (2) = (y (jj) - y (ii) ) / (x (jj) - x (ii) ) 
      pinc (2) = 1.0 
!                                                                       
      ifen = fit_ifen 
      CALL do_fmax_xy (ikfit, wmax, ixm, maxmax, ima) 
      IF (ima.lt.np1) then 
         ier_num = - 30 
         ier_typ = ER_APPL 
         CALL errlist 
         DO i = 1, np1 
         wmax (i) = 1.0 
         ixm (i) = 1 
         ENDDO 
      ELSE 
         DO i = ima + 1, np1 
         wmax (i) = wmax (ima) 
         ixm (i) = ixm (ima) 
         ENDDO 
      ENDIF 
      CALL no_error 
!                                                                       
      IF (ier_num.eq.0) then 
         DO i = 1, np1 
         p (2 + (i - 1) * 4 + 1) = wmax (i) 
         pinc (2 + (i - 1) * 4 + 1) = 1.0 
         p (2 + (i - 1) * 4 + 2) = x (ii + ixm (i) - 1) 
         pinc (2 + (i - 1) * 4 + 2) = 1.0 
         p (2 + (i - 1) * 4 + 3) = 0.2 * abs (x (jj) - x (ii) ) 
         pinc (2 + (i - 1) * 4 + 3) = 1.0 
         p (2 + (i - 1) * 4 + 4) = 1.0 
         pinc (2 + (i - 1) * 4 + 4) = 0.0 
         ENDDO 
!                                                                       
         DO i = 1, npara 
         dp (i) = 0.0 
         ENDDO 
      ENDIF 
!                                                                       
      END SUBROUTINE setup_lor                      
!***********************************************************************
      SUBROUTINE theory_lor (xx, f, df, i) 
!                                                                       
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      REAL xx, f, df (maxpara) 
      REAL o1, fwf, fw, xw, rn 
      INTEGER i, ind, na, nu, np, nl, nlauf 
!                                                                       
      DO ind = 1, maxpara 
      df (ind) = 0.0 
      ENDDO 
!                                                                       
      o1 = 2.0 * atan (1.0) 
      nu = 2 
      np = 4 
      nl = np1 
!-------untergrund                                                      
      f = p (1) + p (2) * xx 
      DO nlauf = 1, nl 
      na = nu + (nlauf - 1) * np + 1 
!---------halbwertsbreiten                                              
      IF (xx.le.p (na + 1) ) then 
         fwf = 1.0 / p (na + 3) 
      ELSE 
         fwf = 1.0 * p (na + 3) 
      ENDIF 
      fw = p (na + 2) * fwf 
      xw = xx - p (na + 1) 
      rn = (fw * fw + 4.0 * xw * xw) 
!---------funktionswert berechnen                                       
      f = f + p (na) * fw * fw / rn 
!---------ableitungen berechnen                                         
      IF (i.gt.0) then 
         IF (pinc (na) .ne.0) df (na) = fw * fw / rn 
         IF (pinc (na + 1) .ne.0) df (na + 1) = 8.0 * p (na) * fw * fw *&
         xw / (rn**2)                                                   
         IF (pinc (na + 2) .ne.0) df (na + 2) = 2 * p (na) * fw / rn -  &
         2 * p (na) * fw * fw * fw / (rn**2)                            
         IF (pinc (na + 3) .ne.0) then 
            IF (xx.le.p (na + 1) ) then 
               df (na + 3) = - df (na + 2) * (2.0 * (xx - p (na + 1) )  &
               **2 * o1 / (fw**3) ) * p (na + 2) / (p (na + 3) **2)     
            ELSE 
               df (na + 3) = df (na + 2) * (2.0 * (xx - p (na + 1) ) ** &
               2 * o1 / (fw**3) ) * p (na + 3)                          
            ENDIF 
         ENDIF 
      ENDIF 
      ENDDO 
!-------untergrundsableitungen                                          
      IF (i.gt.0) then 
         IF (pinc (1) .ne.0) df (1) = 1.0 
         IF (pinc (2) .ne.0) df (2) = xx 
      ENDIF 
!                                                                       
      END SUBROUTINE theory_lor                     
!*****7*****************************************************************
!     Gaussian (1D)                                                     
!*****7*****************************************************************
      SUBROUTINE show_gauss (idout) 
!                                                                       
      USE wink_mod
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      REAL o, sqpio, zz, sint, ds0, ds2, ds3, dsint 
      INTEGER idout, i, iii 
!                                                                       
      o = sqrt (4.0 * alog (2.0) ) 
      sqpio = sqrt (REAL(pi)) / o * 0.5 
!                                                                       
      WRITE (idout, 1000) np1 
      WRITE (idout, 1100) 1, p (1), dp (1), pinc (1) 
      WRITE (idout, 1200) 2, p (2), dp (2), pinc (2) 
      DO i = 1, np1 
      WRITE (idout, 1300) i 
      iii = 2 + (i - 1) * 4 
      WRITE (idout, 1400) iii + 1, p (iii + 1), dp (iii + 1), pinc (iii &
      + 1)                                                              
      WRITE (idout, 1500) iii + 2, p (iii + 2), dp (iii + 2), pinc (iii &
      + 2)                                                              
      WRITE (idout, 1600) iii + 3, p (iii + 3), dp (iii + 3), pinc (iii &
      + 3)                                                              
      WRITE (idout, 1700) iii + 4, p (iii + 4), dp (iii + 4), pinc (iii &
      + 4)                                                              
!---------integral berechnen                                            
      zz = p (iii + 3) / p (iii + 4) + p (iii + 3) * p (iii + 4) 
      sint = p (iii + 1) * sqpio * zz 
      ds0 = sqpio * zz 
      ds2 = p (iii + 1) * sqpio * (1.0 / p (iii + 4) + p (iii + 4) ) 
      ds3 = p (iii + 1) * sqpio * ( - p (iii + 3) / p (iii + 4) **2 + p &
      (iii + 3) )                                                       
      dsint = dp (iii + 1) * ds0 + dp (iii + 3) * ds2 + dp (iii + 4)    &
      * ds3                                                             
      WRITE (idout, 1800) sint, dsint 
      ENDDO 
      WRITE (idout, * ) ' ' 
!                                                                       
 1000 FORMAT     (1x,'Fitted',i3,' Gaussian(s) : '/) 
 1100 FORMAT     (3x,'p(',i2,') : backgr. 1 : ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1200 FORMAT     (3x,'p(',i2,') : backgr. 2 : ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1300 FORMAT     (/,1x,'Gaussian : ',i3,/) 
 1400 FORMAT     (3x,'p(',i2,') : peak      : ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1500 FORMAT     (3x,'p(',i2,') : position  : ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1600 FORMAT     (3x,'p(',i2,') : fwhm      : ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1700 FORMAT     (3x,'p(',i2,') : asymmetry : ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1800 FORMAT     (3x,'        integral  : ',g12.6,' +- ',g12.6) 
!                                                                       
      END SUBROUTINE show_gauss                     
!***7*******************************************************************
      SUBROUTINE setup_gauss (ianz, werte, maxw) 
!                                                                       
      USE errlist_mod 
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxmax 
      PARAMETER (maxmax = 50) 
!                                                                       
      INTEGER maxw 
      REAL werte (maxw) 
      REAL wmax (maxmax) 
      INTEGER ixm (maxmax) 
      INTEGER ianz, ii, jj, ima, i 
!                                                                       
      IF (ianz.eq.0) then 
         npara = 6 
         np1 = 1 
      ELSEIF (ianz.eq.1) then 
         ii = nint (werte (1) ) 
         IF (ii.gt.0.and. (2 + 4 * ii) .le.maxpara) then 
            np1 = ii 
            npara = 2 + 4 * np1 
         ELSE 
            ier_num = - 31 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
         RETURN 
      ENDIF 
!                                                                       
      ii = offxy (ikfit - 1) + 1 
      jj = offxy (ikfit - 1) + len (ikfit) 
!                                                                       
      p (1) = y (ii) 
      pinc (1) = 1.0 
      p (2) = (y (jj) - y (ii) ) / (x (jj) - x (ii) ) 
      pinc (2) = 1.0 
!                                                                       
      ifen = fit_ifen 
      CALL do_fmax_xy (ikfit, wmax, ixm, maxmax, ima) 
      IF (ima.lt.np1) then 
         ier_num = - 30 
         ier_typ = ER_APPL 
         CALL errlist 
         DO i = 1, np1 
         wmax (i) = 1.0 
         ixm (i) = 1 
         ENDDO 
      ELSE 
         DO i = ima + 1, np1 
         wmax (i) = wmax (ima) 
         ixm (i) = ixm (ima) 
         ENDDO 
      ENDIF 
      CALL no_error 
!                                                                       
      IF (ier_num.eq.0) then 
         DO i = 1, np1 
         p (2 + (i - 1) * 4 + 1) = wmax (i) 
         pinc (2 + (i - 1) * 4 + 1) = 1.0 
         p (2 + (i - 1) * 4 + 2) = x (ii + ixm (i) - 1) 
         pinc (2 + (i - 1) * 4 + 2) = 1.0 
         p (2 + (i - 1) * 4 + 3) = 0.2 * abs (x (jj) - x (ii) ) 
         pinc (2 + (i - 1) * 4 + 3) = 1.0 
         p (2 + (i - 1) * 4 + 4) = 1.0 
         pinc (2 + (i - 1) * 4 + 4) = 0.0 
         ENDDO 
!                                                                       
         DO i = 1, npara 
         dp (i) = 0.0 
         ENDDO 
      ENDIF 
!                                                                       
      END SUBROUTINE setup_gauss                    
!*********************************************************************  
      SUBROUTINE theory_gauss (xx, f, df, iwert) 
!                                                                       
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      REAL xx, f, df (maxpara) 
      REAL o1, fw, fwf, exx 
      INTEGER iwert, ind, nunt, npar 
      INTEGER nlauf, na 
!                                                                       
      o1 = 4.0 * alog (2.0) 
      nunt = 2 
      npar = 4 
!                                                                       
      DO ind = 1, npara 
      df (ind) = 0.0 
      ENDDO 
!-------untergrund                                                      
      f = p (1) + p (2) * xx 
      DO nlauf = 1, np1 
      na = nunt + (nlauf - 1) * npar + 1 
!---------halbwertsbreiten                                              
      IF (xx.le.p (na + 1) ) then 
         fwf = 1.0 / p (na + 3) 
      ELSE 
         fwf = 1.0 * p (na + 3) 
      ENDIF 
      fw = p (na + 2) * fwf 
!---------funktionswert berechnen                                       
      exx = exp ( - (xx - p (na + 1) ) **2 * o1 / (fw**2) ) 
      f = f + p (na) * exx 
!---------ableitungen berechnen                                         
      IF (iwert.gt.0) then 
         IF (pinc (na) .ne.0) df (na) = exx 
         IF (pinc (na + 1) .ne.0) df (na + 1) = p (na) * exx * (2.0 *   &
         o1 * (xx - p (na + 1) ) / fw**2)                               
         IF (pinc (na + 2) .ne.0) df (na + 2) = p (na) * exx * (2.0 *   &
         (xx - p (na + 1) ) **2 * o1 * fwf / (fw**3) )                  
         IF (pinc (na + 3) .ne.0) then 
            IF (xx.le.p (na + 1) ) then 
               df (na + 3) = - p (na) * exx * (2.0 * (xx - p (na + 1) ) &
               **2 * o1 / (fw**3) ) * p (na + 2) / (p (na + 3) **2)     
            ELSE 
               df (na + 3) = p (na) * exx * (2.0 * (xx - p (na + 1) ) **&
               2 * o1 / (fw**3) ) * p (na + 3)                          
            ENDIF 
         ENDIF 
      ENDIF 
      ENDDO 
!-------untergrundsableitungen                                          
      IF (iwert.gt.0) then 
         IF (pinc (1) .ne.0) df (1) = 1.0 
         IF (pinc (2) .ne.0) df (2) = xx 
      ENDIF 
!                                                                       
      END SUBROUTINE theory_gauss                   
!***7*******************************************************************
!       Pseudo-Voigt                                                    
!***7*******************************************************************
      SUBROUTINE show_psvgt (idout) 
!                                                                       
      USE wink_mod
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER idout, i, iii 
!                                                                       
      WRITE (idout, 1000) np1 
      DO i = 1, n_backgrd 
      WRITE (idout, 1100) i, i, p (i), dp (i), pinc (i) 
      ENDDO 
!     write (idout,1200) 2,p(2),dp(2),pinc(2)                           
!     write (idout,1210) 3,p(3),dp(3),pinc(3)                           
      DO i = 1, np1 
      WRITE (idout, 1300) i 
      iii = n_backgrd+ (i - 1) * 6 
      WRITE (idout, 1400) iii + 1, p (iii + 1), dp (iii + 1), pinc (iii &
      + 1)                                                              
      WRITE (idout, 1500) iii + 2, p (iii + 2), dp (iii + 2), pinc (iii &
      + 2)                                                              
      WRITE (idout, 1600) iii + 3, p (iii + 3), dp (iii + 3), pinc (iii &
      + 3)                                                              
      WRITE (idout, 1700) iii + 4, p (iii + 4), dp (iii + 4), pinc (iii &
      + 4)                                                              
      WRITE (idout, 1800) iii + 5, p (iii + 5), dp (iii + 5), pinc (iii &
      + 5)                                                              
      WRITE (idout, 1900) iii + 6, p (iii + 6), dp (iii + 6), pinc (iii &
      + 6)                                                              
      ENDDO 
      WRITE (idout, * ) ' ' 
!                                                                       
 1000 FORMAT     (1x,'Fitted',i3,' Pseudo-Voigt(s) : '/) 
 1100 FORMAT     (3x,'p(',i2,') : backgr. ',i1,' : ',g12.6,' +- ',      &
     &                   g12.6,4x,'pinc : ',f2.0)                       
!1200 FORMAT     (3x,'p(',i2,') : backgr. 2 : ',g12.6,' +- ',g12.6,     &
!    &                   4x,'pinc : ',f2.0)                             
!1210 FORMAT     (3x,'p(',i2,') : backgr. 3 : ',g12.6,' +- ',g12.6,     &
!    &                   4x,'pinc : ',f2.0)                             
 1300 FORMAT     (/,1x,'Pseudo-Voigt : ',i3,/) 
 1400 FORMAT     (3x,'p(',i2,') : eta       : ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1500 FORMAT     (3x,'p(',i2,') : int. Inten: ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1600 FORMAT     (3x,'p(',i2,') : position  : ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1700 FORMAT     (3x,'p(',i2,') : fwhm      : ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1800 FORMAT     (3x,'p(',i2,') : asymmetry1: ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1900 FORMAT     (3x,'p(',i2,') : asymmetry2: ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
!                                                                       
      END SUBROUTINE show_psvgt                     
!***7*******************************************************************
      SUBROUTINE setup_psvgt (ianz, werte, maxw) 
!                                                                       
      USE errlist_mod 
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxmax 
      PARAMETER (maxmax = 50) 
!                                                                       
      INTEGER nu 
      INTEGER maxw 
      REAL werte (maxw) 
      REAL wmax (maxmax) 
      INTEGER ixm (maxmax) 
      INTEGER ianz, ii, jj, ima, i 
!                                                                       
      IF (ianz.eq.0) then 
         npara = 8 
         np1 = 1 
      ELSEIF (ianz.le.3) then 
         ii = nint (werte (1) ) 
         IF (ii.gt.0.and. (3 + 6 * ii) .le.maxpara) then 
            np1 = ii 
            npara = 3 + 6 * np1 
         ELSE 
            ier_num = - 31 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
         p_origin = 0.0 
         n_backgrd = 2 
         IF (ianz.ge.2) then 
            p_origin = werte (2) 
         ENDIF 
         IF (ianz.eq.3) then 
            n_backgrd = nint (werte (3) ) 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
         RETURN 
      ENDIF 
      nu = n_backgrd 
!                                                                       
      ii = offxy (ikfit - 1) + 1 
      jj = offxy (ikfit - 1) + len (ikfit) 
!                                                                       
      p (1) = y (ii) 
      pinc (1) = 1.0 
      p (2) = (y (jj) - y (ii) ) / (x (jj) - x (ii) ) 
      pinc (2) = 1.0 
      DO i = 3, n_backgrd 
      p (i) = 0.0 
      pinc (i) = 1.0 
      ENDDO 
!                                                                       
      ifen = fit_ifen 
      CALL do_fmax_xy (ikfit, wmax, ixm, maxmax, ima) 
      IF (ima.lt.np1) then 
         ier_num = - 30 
         ier_typ = ER_APPL 
         CALL errlist 
         DO i = 1, np1 
         wmax (i) = 1.0 
         ixm (i) = 1 
         ENDDO 
      ELSE 
         DO i = ima + 1, np1 
         wmax (i) = wmax (ima) 
         ixm (i) = ixm (ima) 
         ENDDO 
      ENDIF 
      CALL no_error 
!                                                                       
      IF (ier_num.eq.0) then 
         DO i = 1, np1 
         p (nu + (i - 1) * 6 + 1) = 0.5 
         pinc (nu + (i - 1) * 6 + 1) = 1.0 
         p (nu + (i - 1) * 6 + 2) = wmax (i) * 0.1 * abs (x (jj)        &
         - x (ii) ) * 3.14 / 2.                                         
         pinc (nu + (i - 1) * 6 + 2) = 1.0 
         p (nu + (i - 1) * 6 + 3) = x (ii + ixm (i) - 1) 
         pinc (nu + (i - 1) * 6 + 3) = 1.0 
         p (nu + (i - 1) * 6 + 4) = 0.1 * abs (x (jj) - x (ii) ) 
         pinc (nu + (i - 1) * 6 + 4) = 1.0 
         p (nu + (i - 1) * 6 + 5) = 0.0 
         pinc (nu + (i - 1) * 6 + 5) = 0.0 
         p (nu + (i - 1) * 6 + 6) = 0.0 
         pinc (nu + (i - 1) * 6 + 6) = 0.0 
         ENDDO 
!                                                                       
         DO i = 1, npara 
         dp (i) = 0.0 
         ENDDO 
      ENDIF 
!                                                                       
      END SUBROUTINE setup_psvgt                    
!*********************************************************************  
      SUBROUTINE theory_psvgt (xx, f, df, i) 
!                                                                       
      USE kuplot_config 
      USE kuplot_mod 
      USE trig_degree_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      REAL xx, f, df (maxpara) 
      INTEGER i 
!                                                                       
      REAL fw, xw 
      REAL eta, pseudo 
      REAL asym, lore, lorn, pi 
      REAL gaus 
      REAL vln2, gpre 
      REAL zz, fa, fb 
      REAL dgdpos, dldpos 
      REAL dgdfw, dldfw 
      INTEGER j, ind, na, nu, np, nl, nlauf 
!                                                                       
!     REAL tand 
!                                                                       
      DO ind = 1, maxpara 
      df (ind) = 0.0 
      ENDDO 
!                                                                       
      pi = 4.0 * atan (1.0) 
      vln2 = 4. * alog (2.) 
      gpre = 2. * sqrt (alog (2.) / pi) 
      nu = n_backgrd 
      np = 6 
      nl = np1 
!-------untergrund                                                      
      f = 0 
      DO j = 1, nu 
      f = f + p (j) * (xx - p_origin) ** (j - 1) 
      ENDDO 
!                                                                       
      DO nlauf = 1, nl 
      na = nu + (nlauf - 1) * np + 1 
      eta = p (na) 
      xw = xx - p (na + 2) 
      fw = p (na + 3) 
!---------asymmetry                                                     
      zz = xw / fw 
      fa = 2. * zz * exp ( - zz**2) 
      fb = 2. * (2 * zz**2 - 3.) * fa 
      asym = 1.0 
      asym = asym + (p (na + 4) * fa + p (na + 5) * fb) / tand (p (na + &
      2) )                                                              
!DBG        asym = asym+(p(na+4)*fa + p(na+5)*fb)/tand(p(na+2)*0.5)     
!     --Lorentzian                                                      
      lorn = (fw * fw + 4.0 * xw * xw) 
      lore = 2. / pi * fw / lorn 
      gaus = gpre / fw * exp ( - vln2 / fw**2 * xw**2) 
      pseudo = (eta * lore+ (1 - eta) * gaus) 
!---------calculate pseudo Voigt                                        
      f = f + p (na + 1) * asym * pseudo 
!---------calculate derivatives                                         
      IF (i.gt.0) then 
!     ---- eta                                                          
         IF (pinc (na) .ne.0) then 
            df (na) = p (na + 1) * asym * (lore-gaus) 
         ENDIF 
!     ---- intensity                                                    
         IF (pinc (na + 1) .ne.0) then 
            df (na + 1) = asym * pseudo 
         ENDIF 
!     ---- position                                                     
         IF (pinc (na + 2) .ne.0) then 
            dldpos = 2. / pi * 8. * fw * xw / lorn**2 
            dgdpos = gaus * (2. * vln2 / fw**2 * xw) 
            df (na + 2) = asym * p (na + 1) * (eta * dldpos + (1. - eta)&
            * dgdpos)                                                   
         ENDIF 
!     ---- FWHM                                                         
         IF (pinc (na + 3) .ne.0) then 
            dldfw = 2. / pi * ( - 1. * fw * fw + 4. * xw * xw) / lorn** &
            2                                                           
            dgdfw = gaus * ( - 1. / fw + 2. * vln2 / fw**3 * xw * xw) 
            df (na + 3) = asym * p (na + 1) * (eta * dldfw + (1. - eta) &
            * dgdfw)                                                    
         ENDIF 
!     ---- asymmetry parameter 1                                        
         IF (pinc (na + 4) .ne.0) then 
            df (na + 4) = p (na + 1) * pseudo * fa / tand (p (na + 2) ) 
         ENDIF 
!     ---- asymmetry parameter 2                                        
         IF (pinc (na + 5) .ne.0) then 
            df (na + 5) = p (na + 1) * pseudo * fb / tand (p (na + 2) ) 
         ENDIF 
      ENDIF 
      ENDDO 
!-------derivative of background                                        
      IF (i.gt.0) then 
         DO j = 1, nu 
         IF (pinc (j) .ne.0) df (j) = (xx - p_origin) ** (j - 1) 
         ENDDO 
      ENDIF 
!                                                                       
      END SUBROUTINE theory_psvgt                   
!***7*******************************************************************
!     Gaussian (2D)                                                     
!***7*******************************************************************
      SUBROUTINE show_gauss_2d (idout) 
!                                                                       
      USE wink_mod
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      REAL o, sqpio, sqpio2 
      REAL zz_x, zz_y, sint, ds0, ds3, ds4, ds6, ds7, dsint 
      INTEGER idout, i, iii 
!                                                                       
      o = sqrt (4.0 * alog (2.0) ) 
      sqpio = sqrt (REAL(pi)) / o * 0.5 
      sqpio2 = sqpio * sqpio 
!                                                                       
      WRITE (idout, 1000) np1 
      WRITE (idout, 1100) 1, p (1), dp (1), pinc (1) 
      DO i = 1, np1 
      WRITE (idout, 1300) i 
      iii = 1 + (i - 1) * 8 
      WRITE (idout, 1400) iii + 1, p (iii + 1), dp (iii + 1), pinc (iii &
      + 1)                                                              
      WRITE (idout, 1410) iii + 2, p (iii + 2), dp (iii + 2), pinc (iii &
      + 2)                                                              
      WRITE (idout, 1420) iii + 3, p (iii + 3), dp (iii + 3), pinc (iii &
      + 3)                                                              
      WRITE (idout, 1430) iii + 4, p (iii + 4), dp (iii + 4), pinc (iii &
      + 4)                                                              
      WRITE (idout, 1440) iii + 5, p (iii + 5), dp (iii + 5), pinc (iii &
      + 5)                                                              
      WRITE (idout, 1450) iii + 6, p (iii + 6), dp (iii + 6), pinc (iii &
      + 6)                                                              
      WRITE (idout, 1460) iii + 7, p (iii + 7), dp (iii + 7), pinc (iii &
      + 7)                                                              
      WRITE (idout, 1470) iii + 8, p (iii + 8), dp (iii + 8), pinc (iii &
      + 8)                                                              
      zz_x = p (iii + 4) / p (iii + 7) + p (iii + 4) * p (iii + 7) 
      zz_y = p (iii + 5) / p (iii + 8) + p (iii + 5) * p (iii + 8) 
      sint = p (iii + 1) * sqpio2 * zz_x * zz_y 
      ds0 = sqpio2 * zz_x * zz_y 
      ds3 = p (iii + 1) * sqpio2 * zz_y * (1.0 / p (iii + 7) + p (iii + &
      7) )                                                              
      ds4 = p (iii + 1) * sqpio2 * zz_x * (1.0 / p (iii + 8) + p (iii + &
      8) )                                                              
      ds6 = p (iii + 1) * sqpio2 * zz_y * ( - p (iii + 4) / p (iii + 7) &
      **2 + p (iii + 4) )                                               
      ds7 = p (iii + 1) * sqpio2 * zz_x * ( - p (iii + 5) / p (iii + 8) &
      **2 + p (iii + 5) )                                               
      dsint = dp (iii + 1) * ds0 + dp (iii + 4) * ds3 + dp (iii + 5)    &
      * ds4 + dp (iii + 7) * ds6 + dp (iii + 8) * ds7                   
      WRITE (idout, 1800) sint, dsint 
      ENDDO 
      WRITE (idout, * ) ' ' 
!                                                                       
 1000 FORMAT     (1x,'Fitted',i3,' Gaussian(s) : '/) 
 1100 FORMAT     (3x,'p(',i2,') : backgr. 1 : ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1300 FORMAT     (/,1x,'Gaussian : ',i3,/) 
 1400 FORMAT     (3x,'p(',i2,') : peak      : ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1410 FORMAT     (3x,'p(',i2,') : position x: ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1420 FORMAT     (3x,'p(',i2,') : position y: ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1430 FORMAT     (3x,'p(',i2,') : fwhm a    : ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1440 FORMAT     (3x,'p(',i2,') : fwhm b    : ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1450 FORMAT     (3x,'p(',i2,') : angle a,x : ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1460 FORMAT     (3x,'p(',i2,') : asym. a   : ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1470 FORMAT     (3x,'p(',i2,') : asym. b   : ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1800 FORMAT     (3x,'        integral  : ',g12.6,' +- ',g12.6) 
!                                                                       
      END SUBROUTINE show_gauss_2d                  
!***7*******************************************************************
      SUBROUTINE setup_gauss_2d (ianz, werte, maxw) 
!                                                                       
      USE errlist_mod 
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxmax 
      PARAMETER (maxmax = 20) 
!                                                                       
      INTEGER maxw 
      REAL werte (maxw) 
      REAL wmax (maxmax) 
      INTEGER ixm (maxmax), iym (maxmax) 
      INTEGER ianz, ima, ii, jj, i 
!                                                                       
      IF (ianz.eq.0) then 
         npara = 1 + 8 
         np1 = 1 
      ELSEIF (ianz.eq.1) then 
         ii = nint (werte (1) ) 
         IF (ii.gt.0.and. (1 + 8 * ii) .le.maxpara) then 
            np1 = ii 
            npara = 1 + 8 * np1 
         ELSE 
            ier_num = - 31 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
         RETURN 
      ENDIF 
!                                                                       
      ii = offxy (ikfit - 1) 
      jj = offz (ikfit - 1) 
!                                                                       
!-------untergrund = wert li unten                                      
!                                                                       
      p (1) = z (jj + 1) 
      pinc (1) = 1.0 
!                                                                       
!-------nach lokalen maxima suchen                                      
!                                                                       
      ifen = fit_ifen 
      CALL do_fmax_z (ikfit, wmax, ixm, iym, maxmax, ima) 
      IF (ima.lt.np1) then 
         ier_num = - 30 
         ier_typ = ER_APPL 
         CALL errlist 
         DO i = 1, np1 
         wmax (i) = 1.0 
         ixm (i) = 1 
         iym (i) = 1 
         ENDDO 
      ELSE 
         DO i = ima + 1, np1 
         wmax (i) = wmax (ima) 
         ixm (i) = ixm (ima) 
         iym (i) = iym (ima) 
         ENDDO 
      ENDIF 
      CALL no_error 
!                                                                       
!-------parameter setzen                                                
!                                                                       
      IF (ier_num.eq.0) then 
         DO i = 1, np1 
         p (1 + (i - 1) * 8 + 1) = wmax (i) 
         pinc (1 + (i - 1) * 8 + 1) = 1.0 
         p (1 + (i - 1) * 8 + 2) = x (ii + ixm (i) ) 
         pinc (1 + (i - 1) * 8 + 2) = 1.0 
         p (1 + (i - 1) * 8 + 3) = y (ii + iym (i) ) 
         pinc (1 + (i - 1) * 8 + 3) = 1.0 
         p (1 + (i - 1) * 8 + 4) = 0.2 * abs (x (ii + nx (ikfit) )      &
         - x (ii + 1) )                                                 
         pinc (1 + (i - 1) * 8 + 4) = 1.0 
         p (1 + (i - 1) * 8 + 5) = 0.2 * abs (y (ii + ny (ikfit) )      &
         - y (ii + 1) )                                                 
         pinc (1 + (i - 1) * 8 + 5) = 1.0 
         p (1 + (i - 1) * 8 + 6) = 0.0 
         pinc (1 + (i - 1) * 8 + 6) = 0.0 
         p (1 + (i - 1) * 8 + 7) = 1.0 
         pinc (1 + (i - 1) * 8 + 7) = 0.0 
         p (1 + (i - 1) * 8 + 8) = 1.0 
         pinc (1 + (i - 1) * 8 + 8) = 0.0 
         ENDDO 
!                                                                       
         DO i = 1, npara 
         dp (i) = 0.0 
         ENDDO 
      ENDIF 
!                                                                       
      END SUBROUTINE setup_gauss_2d                 
!*********************************************************************  
      SUBROUTINE theory_gauss_2d (xx, f, df, i) 
!                                                                       
      USE wink_mod
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      REAL xx, f, o1, df (maxpara) 
      REAL rx, ry, rxs, rys, fwfx, fwx, fwfy, fwy, cosp, sinp 
      REAL exx, eyy, dfxs, dfys 
      INTEGER i, ind, nlauf, na, nu, np, ng, ix, iy 
!                                                                       
      DO ind = 1, maxpara 
      df (ind) = 0.0 
      ENDDO 
!                                                                       
      o1 = 4.0 * alog (2.0) 
      nu = 1 
      np = 8 
      ng = np1 
!-------x und y bestimmen                                               
      ix = (abs (i) - 1) / ny (ikfit) + 1 
      iy = abs (i) - (ix - 1) * ny (ikfit) 
      rx = x (offxy (ikfit - 1) + ix) 
      ry = y (offxy (ikfit - 1) + iy) 
!-------untergrund                                                      
      f = p (1) 
      DO nlauf = 1, ng 
      na = nu + (nlauf - 1) * np + 1 
!---------sinus und cosinus berechnen                                   
      cosp = cos (REAL(rad) * p (na + 5) ) 
      sinp = sin (REAL(rad) * p (na + 5) ) 
!---------transformation in hautachsensystem                            
      rxs = cosp * (rx - p (na + 1) ) + sinp * (ry - p (na + 2) ) 
      rys = - sinp * (rx - p (na + 1) ) + cosp * (ry - p (na + 2) ) 
!---------halbwertsbreiten                                              
      IF (rxs.le.0.0) then 
         fwfx = 1.0 / p (na + 6) 
      ELSE 
         fwfx = 1.0 * p (na + 6) 
      ENDIF 
      fwx = p (na + 3) * fwfx 
      IF (rys.le.0.0) then 
         fwfy = 1.0 / p (na + 7) 
      ELSE 
         fwfy = 1.0 * p (na + 7) 
      ENDIF 
      fwy = p (na + 4) * fwfy 
!---------funktionswert berechnen                                       
      exx = exp ( - rxs**2 * o1 / (fwx**2) ) 
      eyy = exp ( - rys**2 * o1 / (fwy**2) ) 
      f = f + p (na) * exx * eyy 
!---------ableitungen berechnen                                         
      IF (i.gt.0) then 
         IF (pinc (na) .ne.0) df (na) = exx * eyy 
         dfxs = p (na) * eyy * exx * ( - 2 * rxs * o1 / (fwx**2) ) 
         dfys = p (na) * eyy * exx * ( - 2 * rys * o1 / (fwy**2) ) 
         IF (pinc (na + 1) .ne.0) df (na + 1) = - cosp * dfxs + sinp *  &
         dfys                                                           
         IF (pinc (na + 2) .ne.0) df (na + 2) = - sinp * dfxs - cosp *  &
         dfys                                                           
         IF (pinc (na + 3) .ne.0) df (na + 3) = p (na) * eyy * exx *    &
         (2 * rxs**2 * o1 * fwfx / (fwx**3) )                           
         IF (pinc (na + 4) .ne.0) df (na + 4) = p (na) * eyy * exx *    &
         (2 * rys**2 * o1 * fwfy / (fwy**3) )                           
         IF (pinc (na + 5) .ne.0) df (na + 5) = ( - (rx - p (na + 1) )  &
         * sinp + (ry - p (na + 2) ) * cosp) * dfxs * REAL(rad)+ ( - (rx - p (&
         na + 1) ) * cosp - (ry - p (na + 2) ) * sinp) * dfys * REAL(rad)     
         IF (pinc (na + 6) .ne.0) then 
            IF (rxs.le.0.0) then 
               df (na + 6) = - p (na) * eyy * exx * (2 * rxs**2 * o1 /  &
               (fwx**3) ) * p (na + 3) / (p (na + 6) **2)               
            ELSE 
               df (na + 6) = p (na) * eyy * exx * (2 * rxs**2 * o1 /    &
               (fwx**3) ) * p (na + 3)                                  
            ENDIF 
         ENDIF 
         IF (pinc (na + 7) .ne.0) then 
            IF (rys.le.0.0) then 
               df (na + 7) = - p (na) * eyy * exx * (2 * rys**2 * o1 /  &
               (fwy**3) ) * p (na + 4) / (p (na + 7) **2)               
            ELSE 
               df (na + 7) = p (na) * eyy * exx * (2 * rys**2 * o1 /    &
               (fwy**3) ) * p (na + 4)                                  
            ENDIF 
         ENDIF 
      ENDIF 
      ENDDO 
!-------untergrundsableitungen                                          
      IF (i.gt.0) then 
         IF (pinc (1) .ne.0) df (1) = 1.0 
      ENDIF 
!                                                                       
      END SUBROUTINE theory_gauss_2d                
!***7*******************************************************************
!     Chebyshev polynom                                                 
!***7*******************************************************************
      SUBROUTINE show_poly_cheb (idout) 
!                                                                       
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER idout, i 
!                                                                       
      WRITE (idout, 1000) np1 
      DO i = 0, np1 
      WRITE (idout, 1100) i + 1, i, p (i + 1), dp (i + 1), pinc (i + 1) 
      ENDDO 
      WRITE (idout, * ) ' ' 
!                                                                       
 1000 FORMAT     (1x,'Fitted Chebyshev polynom of order ',i2,' : '/) 
 1100 FORMAT     (3x,'p(',i2,') : coeff. A(',i2,') : ',g12.6,           &
     &                   ' +- ',g12.6,4x,'pinc : ',f2.0)                
!                                                                       
      END SUBROUTINE show_poly_cheb                 
!***7*******************************************************************
      SUBROUTINE setup_poly_cheb (ianz, werte, maxw) 
!                                                                       
      USE errlist_mod 
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      REAL werte (maxw) 
      INTEGER ianz, ii, jj, i 
!                                                                       
      IF (ianz.eq.0) then 
         npara = 1 
         np1 = 1 
      ELSEIF (ianz.eq.1) then 
         ii = nint (werte (1) ) 
         IF (ii.ge.0.and. (1 + ii) .le.maxpara.and. (1 + ii) .le.5)     &
         then                                                           
            np1 = ii 
            npara = ii + 1 
         ELSE 
            ier_num = - 31 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
         RETURN 
      ENDIF 
!                                                                       
      ii = offxy (ikfit - 1) + 1 
      jj = offxy (ikfit - 1) + len (ikfit) 
!                                                                       
      p (1) = y (ii) 
      pinc (1) = 1.0 
      p (2) = (y (jj) - y (ii) ) / (x (jj) - x (ii) ) 
      pinc (2) = 1.0 
      DO i = 3, npara 
      p (i) = 0.0 
      pinc (i) = 1.0 
      ENDDO 
!                                                                       
      DO i = 1, npara 
      dp (i) = 0.0 
      ENDDO 
!                                                                       
      END SUBROUTINE setup_poly_cheb                
!***7*******************************************************************
      SUBROUTINE theory_poly_cheb (xx, f, df, iwert) 
!                                                                       
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      REAL xx, f, df (maxpara) 
      REAL xnew 
      INTEGER iwert, ind 
!                                                                       
      DO ind = 1, npara 
      df (ind) = 0.0 
      ENDDO 
!                                                                       
!------ Map x value on [-1,1] scale                                     
!                                                                       
      xnew = ( (xx - xmin (ikfit) ) - (xmax (ikfit) - xx) ) / (xmax (   &
      ikfit) - xmin (ikfit) )                                           
!                                                                       
      f = p (1) 
      IF (np1.ge.1) f = f + p (2) * xnew 
      IF (np1.ge.2) f = f + p (3) * (2 * xnew**2 - 1) 
      IF (np1.ge.3) f = f + p (4) * (4 * xnew**3 - 3 * xnew) 
      IF (np1.ge.4) f = f + p (5) * (8 * xnew**4 - 8 * xnew**2 + 1) 
!                                                                       
!-------Derivatives                                                     
!                                                                       
      IF (iwert.gt.0) then 
         IF (pinc (1) .ne.0.0) df (1) = 1.0 
         IF (pinc (2) .ne.0.0.and.np1.ge.1) df (2) = xnew 
         IF (pinc (3) .ne.0.0.and.np1.ge.2) df (3) = 2 * xnew**2 - 1 
         IF (pinc (4) .ne.0.0.and.np1.ge.3) df (4) = 4 * xnew**3 - 3 *  &
         xnew                                                           
         IF (pinc (5) .ne.0.0.and.np1.ge.4) df (5) = 8 * xnew**4 - 8 *  &
         xnew**2 + 1                                                    
      ENDIF 
!                                                                       
      END SUBROUTINE theory_poly_cheb               
!***7*******************************************************************
!     Polynom                                                           
!***7*******************************************************************
      SUBROUTINE show_poly (idout) 
!                                                                       
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER idout, i 
!                                                                       
      WRITE (idout, 1000) np1 
      DO i = 0, np1 
      WRITE (idout, 1100) i + 1, i, p (i + 1), dp (i + 1), pinc (i + 1) 
      ENDDO 
      WRITE (idout, * ) ' ' 
!                                                                       
 1000 FORMAT     (1x,'Fitted polynom of order ',i2,' : '/) 
 1100 FORMAT     (3x,'p(',i2,') : coeff. for x**',i2,' : ',g12.6,       &
     &                   ' +- ',g12.6,4x,'pinc : ',f2.0)                
!                                                                       
      END SUBROUTINE show_poly                      
!***7*******************************************************************
      SUBROUTINE setup_poly (ianz, werte, maxw) 
!                                                                       
      USE errlist_mod 
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      REAL werte (maxw) 
      INTEGER ianz, ii, jj, i 
!                                                                       
      IF (ianz.eq.0) then 
         npara = 1 
         np1 = 1 
      ELSEIF (ianz.eq.1) then 
         ii = nint (werte (1) ) 
         IF (ii.ge.0.and. (1 + ii) .le.maxpara) then 
            np1 = ii 
            npara = ii + 1 
         ELSE 
            ier_num = - 31 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
         RETURN 
      ENDIF 
!                                                                       
      ii = offxy (ikfit - 1) + 1 
      jj = offxy (ikfit - 1) + len (ikfit) 
!                                                                       
      p (1) = y (ii) 
      pinc (1) = 1.0 
      p (2) = (y (jj) - y (ii) ) / (x (jj) - x (ii) ) 
      pinc (2) = 1.0 
      DO i = 3, npara 
      p (i) = 0.0 
      pinc (i) = 1.0 
      ENDDO 
!                                                                       
      DO i = 1, npara 
      dp (i) = 0.0 
      ENDDO 
!                                                                       
      END SUBROUTINE setup_poly                     
!***7*******************************************************************
      SUBROUTINE theory_poly (xx, f, df, iwert) 
!                                                                       
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      REAL xx, f, df (maxpara) 
      INTEGER iwert, ind 
!                                                                       
      DO ind = 1, npara 
      df (ind) = 0.0 
      ENDDO 
!                                                                       
      f = p (1) 
      DO ind = 1, np1 
      IF (xx.ne.0) f = f + p (ind+1) * (xx**ind) 
      ENDDO 
!                                                                       
!-------Derivatives                                                     
!                                                                       
      IF (iwert.gt.0) then 
         DO ind = 0, np1 
         IF (pinc (ind+1) .ne.0) then 
            df (ind+1) = (xx**ind) 
         ENDIF 
         ENDDO 
      ENDIF 
!                                                                       
      END SUBROUTINE theory_poly                    
!***7*******************************************************************
!     Scale factor + Background Polynom                                 
!***7*******************************************************************
      SUBROUTINE show_backpoly (idout) 
!                                                                       
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER idout, i 
!                                                                       
      WRITE (idout, 1000) np1 - 2 
      WRITE (idout, 1020) ikfit2 
      WRITE (idout, 1050) p (1), dp (1), pinc (1) 
      DO i = 2, np1 
      WRITE (idout, 1100) i, i - 2, p (i), dp (i), pinc (i) 
      ENDDO 
      WRITE (idout, * ) ' ' 
!                                                                       
 1000 FORMAT     (1x,'Fitted Background polynom of order ',i2,' : '/) 
 1020 FORMAT     (1x,'Constant part from data set        ',i2) 
 1050 FORMAT     (1x,'Scale factor               : ',g12.6,             &
     &                   ' +- ',g12.6,4x,'pinc : ',f2.0)                
 1100 FORMAT     (3x,'p(',i2,') : coeff. for x**',i2,' : ',g12.6,       &
     &                   ' +- ',g12.6,4x,'pinc : ',f2.0)                
!                                                                       
      END SUBROUTINE show_backpoly                  
!***7*******************************************************************
      SUBROUTINE setup_backpoly (ianz, werte, maxw) 
!                                                                       
      USE errlist_mod 
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      REAL werte (maxw) 
      INTEGER ianz, ii, jj, i 
!                                                                       
      IF (ianz.eq.1) then 
         ikfit2 = nint (werte (1) ) 
         npara = 1 
         np1 = 1 
      ELSEIF (ianz.eq.2) then 
         ikfit2 = nint (werte (1) ) 
         ii = nint (werte (2) ) 
         IF (ii.ge.0.and. (1 + ii) .le.maxpara) then 
            np1 = ii 
            npara = ii + 1 
         ELSE 
            ier_num = - 31 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
         RETURN 
      ENDIF 
!                                                                       
      ii = offxy (ikfit - 1) + 1 
      jj = offxy (ikfit2 - 1) + 1 
!                                                                       
      IF (np1.ge.2) then 
         p (2) = y (ii) - y (jj) 
         pinc (2) = 1.0 
      ELSE 
         pinc (2) = 0.0 
      ENDIF 
!                                                                       
      ii = offxy (ikfit - 1) + len (ikfit) / 2 
      jj = offxy (ikfit2 - 1) + len (ikfit2) / 2 
      IF (y (jj) .ne.0) then 
         p (1) = (y (ii) - p (2) ) / y (jj) 
      ENDIF 
      pinc (1) = 1.0 
      DO i = 3, npara 
      p (i) = 0.0 
      pinc (i) = 1.0 
      ENDDO 
!                                                                       
      DO i = 1, npara 
      dp (i) = 0.0 
      ENDDO 
!                                                                       
      END SUBROUTINE setup_backpoly                 
!***7*******************************************************************
      SUBROUTINE theory_backpoly (xx, f, df, iwert) 
!                                                                       
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      REAL xx, f, df (maxpara) 
      INTEGER iwert, ind 
      INTEGER iii 
!                                                                       
      INTEGER idata, ipoint 
      REAL arg 
!                                                                       
      DATA iii / 0 / 
!                                                                       
      DO ind = 1, npara 
      df (ind) = 0.0 
      ENDDO 
!                                                                       
      idata = ikfit2 
      ipoint = iabs (iwert) 
      arg = (xx - xmin (idata) ) 
!                                                                       
      f = p (1) * y (offxy (idata - 1) + ipoint) 
      IF (np1.ge.2) then 
         f = f + p (2) 
      ENDIF 
      DO ind = 3, np1 
      f = f + p (ind) * (arg** (ind-2) ) 
      ENDDO 
!                                                                       
!-------Derivatives                                                     
!                                                                       
      IF (iwert.gt.0) then 
         ind = 1 
         IF (pinc (ind) .ne.0) then 
            df (ind) = y (offxy (idata - 1) + ipoint) 
         ENDIF 
         DO ind = 2, np1 
         IF (pinc (ind) .ne.0) then 
            df (ind) = (arg** (ind-2) ) 
         ENDIF 
         ENDDO 
      ENDIF 
!                                                                       
      END SUBROUTINE theory_backpoly                
!***7*******************************************************************
!     Altered subroutine FITTE for KUPLOT use                           
!***7*******************************************************************
SUBROUTINE fit_kupl (ww) 
!+                                                                      
!     eigentliche fit-routine speziell fuer kupl ..                     
!-                                                                      
USE param_mod 
USE prompt_mod 
USE kuplot_config 
USE kuplot_mod 
USE errlist_mod
!                                                                       
IMPLICIT NONE
!
REAL, DIMENSION(MAXARRAY), INTENT(IN) ::  ww !(maxarray) 
!
INTEGER :: i, j, k,l, l1, nl   ! Loop indices
INTEGER :: m, n, nf
INTEGER :: iiw, iix
REAL :: f
REAL :: g
REAL :: h
REAL :: s
REAL :: sqsum
REAL :: sum3
REAL :: sum4
REAL :: syi
REAL :: ttt
REAL :: xx
REAL :: zl
REAL, DIMENSION(MAXPARA) ::  df, dl, dz, pl
REAL, DIMENSION(MAXPARA) ::  pinc_old
!                                                                       
zl = 0.0
IF (lni (ikfit) ) THEN 
   iiw = offz (ikfit - 1) 
   iix = offxy (ikfit - 1) 
   m = nx (ikfit) * ny (ikfit) 
ELSE 
   iiw = offxy (ikfit - 1) 
   iix = offxy (ikfit - 1) 
   m = len (ikfit) 
ENDIF 
!                                                                       
DO i = 1, npara 
   pinc_old (i) = pinc (i) 
ENDDO 
!                                                                       
IF (fstart) write (output_io, 5000) 
n = npara 
j = 0 
!
cycles: DO
   zwert = 0.0 
   DO k = 1, n 
      dz (k) = 0.0 
      DO l = k, n 
         cl (l, k) = 0.0 
      ENDDO 
   ENDDO 
   points: DO i = 1, m 
      xx = x (iix + i) 
      CALL kupl_theory (xx, f, df, i) 
      s = ww (iiw + i) - f 
      IF (wtyp.eq.'BCK') THEN 
         IF (ww (iiw + i) .ne. - 9999.) THEN 
            w (iiw + i) = EXP ( - wval * s) 
         ELSE 
            w (iiw + i) = 0.0 
         ENDIF 
      ENDIF 
      zwert = zwert + w (iiw + i) * s * s 
      DO k = 1, n 
         IF (abs (df (k) ) .lt.0.0000001) df (k) = 0.0 
         ttt = w (iiw + i) * df (k) 
         dz (k) = dz (k) - s * ttt 
         DO l = k, n 
            cl (l, k) = cl (l, k) + ttt * df (l) 
         ENDDO 
      ENDDO 
   ENDDO  points
!
   f = 1.0 
   IF (urf >  0) THEN
      h = 0.0 
      DO k = 1, n 
         IF(dz(k).ne.0.0) h = h + (dz(k) * dz(k)) / (zwert * cl(k,k))
      ENDDO 
      f = 1.000001 + urf * h 
   ENDIF
!
   IF (fstart) write (output_io, 5010) j, zwert, f 
   zalt = zl 
   zdif = zwert - zl 
   fend = f 
   IF (j >  0) THEN
      IF (zwert - zl .ge. 0) EXIT cycles
   ENDIF
!
      zl = zwert 
      DO k = 1, n 
         pl (k) = p (k) 
      ENDDO 
      s = 0
      loopk: DO k = 1, n 
         loopl: DO l = 1, k 
            l1 = l - 1 
            s = cl (k, l) 
            IF (l - k .eq. 0) EXIT loopl
            IF (l1 >    0) THEN
               DO i = 1, l1 
                  s = s - cl (i, l) * cl (i, k) 
               ENDDO 
            ENDIF
            cl (l, k) = s 
         ENDDO  loopl
!
         IF (s == 0) THEN
            dl (k) = 0.0 
         ELSE
            s = s * f 
            IF (l1 >  0) THEN
               DO i = 1, l1 
                  ttt = cl (i, k) 
                  cl (i, k) = ttt * dl (i) 
                  s = s - ttt * cl (i, k) 
               ENDDO 
            ENDIF
            IF (s >  0) THEN
               dl (k) = 1.0 / s 
            ENDIF
         ENDIF
!
      ENDDO  loopk
!
      IF (j - ncycle .gt. 0) EXIT cycles
      j = j + 1 
!
      IF (n - 1 >   0) THEN
         DO l = 2, n 
            l1 = l - 1 
            DO k = 1, l1 
               dz (l) = dz (l) - cl (k, l) * dz (k) 
            ENDDO 
         ENDDO 
      ENDIF
!
      dz (n) = dz (n) * dl (n) 
      IF (n - 1 >   0) THEN
         DO nl = 2, n 
            l = n - nl + 1 
            l1 = l + 1 
            dz (l) = dz (l) * dl (l) 
            DO k = l1, n 
               dz (l) = dz (l) - cl (l, k) * dz (k) 
            ENDDO 
         ENDDO 
      ENDIF
!
      DO k = 1, n 
         pinc (k) = 0.01 * dz (k) 
         p (k) = p (k) - dz (k) 
         kupl_para(k) = p(k)            ! Update global parameter
      ENDDO 
   IF(ier_ctrlc) THEN
      ier_num = -14
      ier_typ = ER_COMM
      RETURN
   ENDIF
   IF(ier_num/=0) RETURN      ! An error occured or CTRL-C
!                                                                       
ENDDO cycles
!                                                                       
g = 1.0 
h = REAL(m - n) 
IF (h * zl.gt.0.0) g = zl / h 
DO k = 1, n 
   dl (k) = dl (k) * g 
   cl (k, k) = 1.0 
ENDDO 
!
IF (n - 1 >   0) THEN
   DO l = 2, n 
      l1 = l - 1 
      DO k = 1, l1 
         s = 0.0 
         DO i = k, l1 
            s = s - cl (i, l) * cl (i, k) 
         ENDDO 
         cl (l, k) = s 
      ENDDO 
   ENDDO 
ENDIF
!
DO k = 1, n 
   DO l = k, n 
      s = 0.0 
      IF (l - k ==   0) THEN
         DO i = l, n 
            ttt = cl (i, k) 
            cl (i, k) = ttt * dl (i) 
            s = s + ttt * cl (i, k) 
         ENDDO 
      ELSE
         DO i = l, n 
            s = s + cl (i, l) * cl (i, k) 
         ENDDO 
      ENDIF
      cl (l, k) = s 
   ENDDO 
ENDDO 
!
DO l = 1, n 
   IF (dl (l) ==   0 ) THEN
      dp (l) = 1.0 
      dp (l) = 0.0 
   ELSE
      f = SQRT (cl (l, l) ) 
      dp (l) = f 
      DO k = 1, l 
         dz (k) = cl (l, k) / f 
      ENDDO 
      DO k = l, n 
         dz (k) = cl (k, l) / f 
      ENDDO 
      DO k = 1, l 
         IF (dp (k) .eq.0.0) then 
            cl (l, k) = dz (k) 
         ELSE 
            cl (l, k) = dz (k) / dp (k) 
         ENDIF 
      ENDDO 
   ENDIF
ENDDO 
!                                                                       
DO k = 1, n 
   p (k) = pl (k) 
   kupl_para(k) = p(k)            ! Update global parameter
ENDDO 
sqsum = 0.0 
syi = 0.0 
sum3 = 0.0 
sum4 = 0.0 
nf = n 
DO i = 1, n 
   IF (df (i) .eq.0.) nf = nf - 1 
ENDDO 
DO i = 1, m 
   xx = x (iix + i) 
   CALL kupl_theory (xx, f, df, - i) 
   h = ww (iiw + i) - f 
   IF (wtyp.eq.'BCK') then 
      IF (ww (iiw + i) .ne. - 9999.) then 
         w (iiw + i) = EXP ( - wval * h) 
      ELSE 
         w (iiw + i) = 0.0 
      ENDIF 
   ENDIF 
   sum3 = sum3 + w (iiw + i) * h * h 
   sum4 = sum4 + w (iiw + i) * ww (iiw + i) * ww (iiw + i) 
ENDDO 
r4 = SQRT (sum3 / sum4) 
re = SQRT ( (m - nf) / sum4) 
res_para (0) = 1 
res_para (1) = r4 
res_para (2) = re 
!                                                                       
DO i = 1, npara 
   pinc (i) = pinc_old (i) 
ENDDO 
IF (fstart) write (output_io, * ) 
!                                                                       
 5000 FORMAT     (' Starting least square fit ... ',/) 
 5010 FORMAT     (3x,'cycle:',i4,3x,'sum: ',g12.6,3x,' urf: ',g12.6) 
!                                                                       
END SUBROUTINE fit_kupl                       
!
!***7*******************************************************************
!     Routine zur Berechnung der Ableitung aus RECIPES                  
!***7*******************************************************************
      FUNCTION dfridr (x, h, err) 
      INTEGER ntab 
      REAL dfridr, err, h, x, func, con, con2, big, safe 
      PARAMETER (con = 1.4, con2 = con * con, big = 1.e30, ntab = 10,   &
      safe = 2.)                                                        
      INTEGER i, j 
      REAL errt, fac, hh, a (ntab, ntab) 
!                                                                       
      dfridr = 0.0
!     IF (h.eq.0.) pause 'h must be nonzero in dfridr'   ! h is checked before
      hh = h 
      a (1, 1) = (func (x + hh) - func (x - hh) ) / (2.0 * hh) 
      err = big 
      DO i = 2, ntab 
      hh = hh / con 
      a (1, i) = (func (x + hh) - func (x - hh) ) / (2.0 * hh) 
      fac = con2 
      DO j = 2, i 
      a (j, i) = (a (j - 1, i) * fac - a (j - 1, i - 1) ) / (fac - 1.) 
      fac = con2 * fac 
      errt = max (abs (a (j, i) - a (j - 1, i) ), abs (a (j, i) - a (j -&
      1, i - 1) ) )                                                     
      IF (errt.le.err) then 
         err = errt 
         dfridr = a (j, i) 
      ENDIF 
      ENDDO 
      IF (abs (a (i, i) - a (i - 1, i - 1) ) .ge.SAFE * err) return 
      ENDDO 
      RETURN 
      END FUNCTION dfridr                           
!
!END MODULE kuplot_fit_mod
!
!*****FIT6******FIT6******FIT6******FIT6******FIT6******FIT6******FIT6******
!
MODULE kuplot_fit_para
!
! Variables to communicate the MRQMIN Results into KUPLOT
!
IMPLICIT NONE
!
INTEGER :: kup_fit6_MAXP
INTEGER :: kup_fit6_MAXF
INTEGER :: kup_fit6_npara
INTEGER :: kup_fit6_nfixed
INTEGER :: kup_fit6_ndata
REAL    :: kup_fit6_chisq
REAL    :: kup_fit6_conf
REAL    :: kup_fit6_lamda
REAL    :: kup_fit6_r4
REAL    :: kup_fit6_re
CHARACTER(LEN=16), DIMENSION(:),   ALLOCATABLE :: kup_fit6_params
REAL             , DIMENSION(:),   ALLOCATABLE :: kup_fit6_pp
REAL             , DIMENSION(:),   ALLOCATABLE :: kup_fit6_dpp
REAL             , DIMENSION(:,:), ALLOCATABLE :: kup_fit6_prange
REAL             , DIMENSION(:,:), ALLOCATABLE :: kup_fit6_cl
CHARACTER(LEN=16), DIMENSION(:),   ALLOCATABLE :: kup_fit6_fixed
REAL             , DIMENSION(:),   ALLOCATABLE :: kup_fit6_pf
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE kup_fit6_set(MAXP, MAXF, npara, nfixed, ndata, chisq, conf, lamda,  &
           r4, re, &
           params, pp, dpp, prange, cl, fixed, pf)
!
!  Allocate and set the variables MRW ==> KUPLOT
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: MAXP
INTEGER, INTENT(IN) :: MAXF
INTEGER, INTENT(IN) :: npara
INTEGER, INTENT(IN) :: nfixed
INTEGER, INTENT(IN) :: ndata
REAL, INTENT(IN) :: chisq
REAL, INTENT(IN) :: conf
REAL, INTENT(IN) :: lamda
REAL, INTENT(IN) :: r4
REAL, INTENT(IN) :: re
CHARACTER(LEN=*), DIMENSION(MAXP), INTENT(IN) :: params
REAL            , DIMENSION(MAXP), INTENT(IN) :: pp
REAL            , DIMENSION(MAXP), INTENT(IN) :: dpp
REAL            , DIMENSION(MAXP,2), INTENT(IN) :: prange
REAL            , DIMENSION(NPARA,NPARA), INTENT(IN) :: cl
CHARACTER(LEN=*), DIMENSION(MAXF), INTENT(IN) :: fixed
REAL            , DIMENSION(MAXF), INTENT(IN) :: pf
!                                                                       
IF(ALLOCATED(kup_fit6_params)) DEALLOCATE(kup_fit6_params)
IF(ALLOCATED(kup_fit6_pp    )) DEALLOCATE(kup_fit6_pp    )
IF(ALLOCATED(kup_fit6_dpp   )) DEALLOCATE(kup_fit6_dpp   )
IF(ALLOCATED(kup_fit6_prange)) DEALLOCATE(kup_fit6_prange)
IF(ALLOCATED(kup_fit6_cl    )) DEALLOCATE(kup_fit6_cl    )
IF(ALLOCATED(kup_fit6_fixed )) DEALLOCATE(kup_fit6_fixed )
IF(ALLOCATED(kup_fit6_pf    )) DEALLOCATE(kup_fit6_pf    )
ALLOCATE(kup_fit6_params(MAXP))
ALLOCATE(kup_fit6_pp    (MAXP))
ALLOCATE(kup_fit6_dpp   (MAXP))
ALLOCATE(kup_fit6_prange(MAXP,2))
ALLOCATE(kup_fit6_cl    (NPARA, NPARA))
ALLOCATE(kup_fit6_fixed (MAXF))
ALLOCATE(kup_fit6_pf    (MAXF))
!
kup_fit6_MAXP   = MAXP
kup_fit6_MAXF   = MAXF
kup_fit6_npara  = NPARA
kup_fit6_nfixed = nfixed
kup_fit6_ndata  = ndata
kup_fit6_chisq  = chisq
kup_fit6_conf   = conf 
kup_fit6_lamda  = lamda
kup_fit6_r4     = r4
kup_fit6_re     = re
kup_fit6_params(:) = params(:)
kup_fit6_pp    (:) = pp    (:)
kup_fit6_dpp   (:) = dpp   (:)
kup_fit6_prange(:,:) = prange(:,:)
kup_fit6_cl    (:,:) = cl    (:,:)
kup_fit6_fixed (:) = fixed (:)
kup_fit6_pf    (:) = pf    (:)
!
END SUBROUTINE kup_fit6_set
!
END MODULE kuplot_fit_para
!
!*******************************************************************************
!
MODULE kuplot_fit6_set_theory
!
! MRQCOF uses a procedure pointer to call the theory function
! Here a generic interface and the pointer variable are defined.
!
! The interface is set very broad, to reflect as much equality to the refine 
! section and to allow for easy expansion if needed.
! Most function do not really need all the parameters, mostly they are needed
! for the User macro, which incidentally is done better with REFINE...
!
! The kuplot Version 6 fit routine is still needed to be able to perform a 
! LeastSquares fit within a user macro done by REFINE.
!
INTERFACE
   SUBROUTINE kuplot_theory(MAXP, ix, iy, xx, yy, NPARA, params, par_names,     &
                            prange, l_do_deriv, data_dim, data_data,  &
                            data_weight, data_x, data_y, data_calc, kupl_last, &
                            ymod, dyda, LDERIV)
   INTEGER                                              , INTENT(IN)  :: MAXP    ! Parameter array sizes
   INTEGER                                              , INTENT(IN)  :: ix      ! Point number along x
   INTEGER                                              , INTENT(IN)  :: iy      ! Point number along y
   REAL                                                 , INTENT(IN)  :: xx      ! Point value  along x
   REAL                                                 , INTENT(IN)  :: yy      ! Point value  along y
   INTEGER                                              , INTENT(IN)  :: NPARA   ! Number of refined parameters
   REAL            , DIMENSION(MAXP )                   , INTENT(IN)  :: params  ! Parameter values
   CHARACTER(LEN=*), DIMENSION(MAXP)                    , INTENT(IN)  :: par_names    ! Parameter names
   REAL            , DIMENSION(MAXP, 2                 ), INTENT(IN)  :: prange      ! Allowed parameter range
   LOGICAL         , DIMENSION(MAXP )                   , INTENT(IN)  :: l_do_deriv  ! Parameter needs derivative
   INTEGER         , DIMENSION(2)                       , INTENT(IN)  :: data_dim     ! Data array dimensions
   REAL            , DIMENSION(data_dim(1), data_dim(2)), INTENT(IN)  :: data_data    ! Data array
   REAL            , DIMENSION(data_dim(1), data_dim(2)), INTENT(IN)  :: data_weight  ! Data sigmas
   REAL            , DIMENSION(data_dim(1))             , INTENT(IN)  :: data_x       ! Data coordinates x
   REAL            , DIMENSION(data_dim(1))             , INTENT(IN)  :: data_y       ! Data coordinates y
   REAL            , DIMENSION(data_dim(1), data_dim(2)), INTENT(OUT) :: data_calc    ! Data array
   INTEGER                                              , INTENT(IN)  :: kupl_last    ! Last KUPLOT DATA that are needed
   REAL                                                 , INTENT(OUT) :: ymod    ! Function value at (ix,iy)
   REAL            , DIMENSION(NPARA)                   , INTENT(OUT) :: dyda    ! Function derivatives at (ix,iy)
   LOGICAL                                              , INTENT(IN)  :: LDERIV  ! TRUE if derivative is needed
   END SUBROUTINE kuplot_theory
END INTERFACE
!
PROCEDURE(kuplot_theory )   , POINTER :: p_kuplot_theory  => NULL()
!
END MODULE kuplot_fit6_set_theory
!
!*******************************************************************************
!
MODULE kuplot_fit6
!
! Revised as of DISCUS Version 6.0.0 
! Levenberg Marquadt type Least Squares fits within KUPLOT
! The procedure is adapted from the Numerical recipes.
!
!
PRIVATE
PUBLIC do_f66
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
USE kuplot_config 
USE kuplot_mod 
USE string_convert_mod
USE sup_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER, PARAMETER :: maxw = 10 
!                                                                       
CHARACTER (LEN=*), INTENT(INOUT) :: zei 
INTEGER          , INTENT(INOUT) :: lp
!
      CHARACTER(1024) cpara (maxw) 
      CHARACTER(1024) line, zeile
      CHARACTER(LEN=40         ) :: orig_prompt 
      CHARACTER(40) cdummy 
      CHARACTER(4) befehl 
      CHARACTER(LEN=1) :: empty 
      INTEGER lpara (maxw) 
      INTEGER ll
      INTEGER ianz, indxg, lbef 
      INTEGER maxpkt, maxzz 
      REAL werte (maxw) 
      LOGICAL flag (3), sel_func 
!                                                                       
      INTEGER len_str 
      LOGICAL str_comp 
real :: f, df(maxpara)
!
      empty = ' '
      CALL no_error 
      sel_func = ftyp (1:4) .ne.'NONE' 
!                                                                       
!------ check data set number and free space                            
!                                                                       
      CALL get_params (zei, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
!                                                                       
      IF (ianz.eq.1) then 
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         IF (ier_num.ne.0) return 
         ikfit = nint (werte (1) ) 
         IF (ikfit.lt.1.or.ikfit.gt. (iz - 1) ) then 
            ier_num = - 4 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
         IF ( (iz + 1) .gt.maxkurvtot) then 
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
      IF (ikfirst (ikfit) ) then 
!                                                                       
!----- -- Check if there is enough space left                           
!                                                                       
         maxpkt = maxarray - offxy (iz - 1) 
         maxzz = maxarray - offz (iz - 1) 
!                                                                       
         IF (lni (ikfit) ) then 
            IF ( (2.0 * nx (iz - 1) * ny (iz - 1) .gt.maxzz) .or. (max (&
            nx (iz - 1), ny (iz - 1) ) .gt.maxpkt) ) then               
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
            IF (2.0 * len (iz - 1) .gt.maxpkt) then 
               ier_num = - 6 
               ier_typ = ER_APPL 
               RETURN 
            ENDIF 
!                                                                       
            ikcal = iz 
            ikdif = iz + 1 
            offxy (ikcal) = offxy (iz - 1) + len (iz - 1) 
            offxy (ikdif) = offxy (ikcal) + len (iz - 1) 
            offz (ikcal) = offz (iz - 1) 
            offz (ikdif) = offz (ikcal) 
            len (ikcal) = len (ikfit) 
            len (ikdif) = len (ikfit) 
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
      IF (ier_num.eq.0) then 
         IF (line.eq.' '.or.line (1:1) .eq.'#' .OR. line(1:1)=='!') goto 10 
!                                                                       
!------ search for "="                                                  
!                                                                       
         indxg = index (line, '=') 
      IF (indxg.ne.0                                            &
         .and..not. (str_comp (befehl, 'echo', 2, lbef, 4) )    &
         .and..not. (str_comp (befehl, 'help', 2, lbef, 4) .or. &
                     str_comp (befehl, '?   ', 2, lbef, 4) )    &
         .AND. INDEX(line,'==') == 0                          )then
            CALL do_math (line, indxg, ll) 
!                                                                       
!------ execute a macro file                                            
!                                                                       
         ELSEIF (befehl (1:1) .eq.'@') then 
            CALL file_kdo (line (2:ll), ll - 1) 
!                                                                       
!     continues a macro 'continue'                                      
!                                                                       
         ELSEIF (str_comp (befehl, 'continue', 3, lbef, 8) ) then 
            CALL macro_continue (zeile, lp) 
!                                                                       
!-------Set number of cycles 'cyc'                                      
!                                                                       
         ELSEIF (str_comp (befehl, 'cycle', 2, lbef, 5) ) then 
            CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
            IF (ier_num.eq.0) then 
               IF (ianz.eq.1) then 
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  IF (ier_num.eq.0) then 
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
         ELSEIF (str_comp (befehl, 'echo', 2, lbef, 4) ) then 
            CALL echo (zeile, lp) 
!                                                                       
!------ Evaluate an expression                                          
!                                                                       
         ELSEIF (str_comp (befehl, 'eval', 2, lbef, 4) ) then 
            CALL do_eval (zeile, lp, .TRUE.) 
!                                                                       
!     exit 'exit'                                                       
!                                                                       
         ELSEIF (str_comp (befehl, 'exit', 2, lbef, 4) ) then 
            GOTO 9999 
!                                                                       
!     Define fit function 'func'                                        
!                                                                       
         ELSEIF (str_comp (befehl, 'func', 3, lbef, 4) ) then 
            lturn_off = .FALSE.
            CALL do_fit_fkt (zeile, lp) 
            IF (ier_num.eq.0) sel_func = .true. 
!                                                                       
!     help 'help','?'                                                   
!                                                                       
      ELSEIF (str_comp (befehl, 'help', 2, lbef, 4) .or.str_comp (befehl&
     &, '?   ', 1, lbef, 4) ) then                                      
            IF (str_comp (zeile, 'errors', 2, lp, 6) ) then 
               lp = lp + 7 
               CALL do_hel ('kuplot '//zeile, lp) 
            ELSE 
               lp = lp + 12 
               CALL do_hel ('kuplot fit '//zeile, lp) 
            ENDIF 
!                                                                       
!-------Save current parameters to a macro file                         
!                                                                       
         ELSEIF (str_comp (befehl, 'macro', 2, lbef, 5) ) then 
            CALL do_fit_macro (zeile, lp) 
!                                                                       
!-------Set parameter 'ifen' for determination of maxima                
!                                                                       
         ELSEIF (str_comp (befehl, 'mfen', 2, lbef, 4) ) then 
            CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
            IF (ier_num.eq.0) then 
               IF (ianz.eq.1) then 
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  IF (ier_num.eq.0) then 
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
         ELSEIF (str_comp (befehl, 'output', 2, lbef, 6) ) then 
            CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
            IF (ianz.eq.1) then 
               fstart = str_comp (cpara (1) , 'on', 2, lpara (1) , 2) 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
!                                                                       
!-------Toogle fit range 'range'                                        
!                                                                       
         ELSEIF (str_comp (befehl, 'range', 2, lbef, 5) ) then 
            CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
            IF (ianz.eq.1) then 
               frall = str_comp (cpara (1) , 'all', 2, lpara (1) , 2) 
            ELSE 
               IF (frall) then 
                  WRITE (output_io, 2100) 'complete range' 
               ELSE 
                  WRITE (output_io, 2100) 'plot range only' 
               ENDIF 
            ENDIF 
!                                                                       
!-------Set parameters                                                  
!                                                                       
         ELSEIF (str_comp (befehl, 'par', 2, lbef, 3) ) then 
            CALL do_fit_par (zeile, lp) 
!                                                                       
!-------Plot result                                                     
!                                                                       
         ELSEIF (str_comp (befehl, 'plot', 2, lbef, 4) ) then 
            CALL do_plot (.false.) 
!                                                                       
!------ Set scale                                                       
!                                                                       
         ELSEIF (str_comp (befehl, 'skal', 2, lbef, 4) ) then 
            CALL set_skal (zeile, lp) 
!                                                                       
!-------Run fit                                                         
!                                                                       
         ELSEIF (str_comp (befehl, 'run', 2, lbef, 3) ) then 
!
            lturn_off = .TRUE.
            IF (.not.sel_func) then 
               ier_num = - 25 
               ier_typ = ER_APPL 
            ELSE 
               IF(ex(iwin,iframe,1)==ex(iwin,iframe,2)) THEN
                  CALL set_skal(empty,0)
               ENDIF
               IF (lni (ikfit) ) then 
                  CALL do_fit_z 
               ELSE 
                  CALL do_fit_y 
               ENDIF 
               CALL get_extrema 
            ENDIF 
            lturn_off = .FALSE.
!                                                                  
!-------Test a recursive call to kuplot_loop                            
!                                                                       
         ELSEIF (str_comp (befehl, 'test', 2, lbef, 4) ) then 
            CALL theory_macro(1.0,  f, df, 1)
!                                                                  
!-------Save fit results                                                
!                                                                       
         ELSEIF (str_comp (befehl, 'save', 2, lbef, 4) ) then 
            CALL do_fit_save 
!                                                                       
!-------Show settings                                                   
!                                                                       
         ELSEIF (str_comp (befehl, 'show', 2, lbef, 4) ) then 
            CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
            IF (ier_num.eq.0) then 
               IF (ianz.eq.0) then 
                  CALL do_fit_info (output_io, .true., .true., .true.) 
               ELSEIF (ianz.eq.1) then 
                  CALL do_cap (cpara (1) ) 
                  flag (1) = cpara (1) (1:2) .eq.'GE' 
                  flag (2) = cpara (1) (1:2) .eq.'FI' 
                  flag (3) = cpara (1) (1:2) .eq.'PA' 
                  IF (flag (1) .or.flag (2) .or.flag (3) ) then 
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
         ELSEIF (str_comp (befehl, 'syst', 2, lbef, 4) ) then 
            cdummy = ' ' 
            IF (zeile.ne.' ') then 
               cdummy (1:lp) = zeile (1:lp) 
               CALL do_operating (cdummy, lp) 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
!                                                                       
!-------Set URF 'urf'                                                   
!                                                                       
         ELSEIF (str_comp (befehl, 'urf', 2, lbef, 3) ) then 
            CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
            IF (ier_num.eq.0) then 
               IF (ianz.eq.1) then 
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  IF (ier_num.eq.0) then 
                     urf = werte (1) 
                  ENDIF 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
            ENDIF 
!                                                                       
!     Waiting for user input                                            
!                                                                       
         ELSEIF (str_comp (befehl, 'wait', 3, lbef, 4) ) then 
            CALL do_input (zeile, lp) 
!                                                                       
!     Set weighting scheme                                              
!                                                                       
         ELSEIF (str_comp (befehl, 'wic', 3, lbef, 3) ) then 
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
!*****7*****************************************************************
      SUBROUTINE do_fit_fkt (zeile, lp) 
!+                                                                      
!     Set theory function                                               
!-                                                                      
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
      USE kuplot_config 
      USE kuplot_mod 
      USE string_convert_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 8) 
!                                                                       
      CHARACTER ( * ) zeile 
      CHARACTER(1024) cpara (maxw) 
      CHARACTER(80) iname 
      REAL werte (maxw) 
      INTEGER lpara (maxw), lp 
      INTEGER ianz, iianz
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
!                                                                       
      IF (ianz.ge.1) then 
         ftyp = cpara (1) (1:4) 
         CALL do_cap (ftyp) 
         IF (ftyp (1:2) .eq.'GS') then 
            IF (ianz.ne.4.and.ianz.ne.5) then 
               ier_num = - 6 
               ier_typ = ER_COMM 
               RETURN 
            ELSE 
               iname = cpara (4) (1:lpara (4) ) 
               cpara (4) = '0.0' 
               lpara (4) = 3 
            ENDIF 
         ENDIF 
         IF (ftyp(1:2)=='MA' .AND. ianz==3) THEN
            CALL del_params (1, ianz, cpara, lpara, maxw) 
            iianz = 1
            CALL ber_params (iianz, cpara, lpara, werte, maxw) 
         ELSEIF (ianz.gt.1) then 
            CALL del_params (1, ianz, cpara, lpara, maxw) 
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         ELSE 
            ianz = 0 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
      IF (ier_num.ne.0) return 
!                                                                       
!------ Polynom (only 2D)                                               
!                                                                       
      IF (ftyp (1:2) .eq.'PO') then 
         IF (.not.lni (ikfit) ) then 
            CALL setup_poly (ianz, werte, maxw) 
         ELSE 
            ier_num = - 25 
            ier_typ = ER_APPL 
         ENDIF 
!                                                                       
!------ Scale + Background Polynom (only 2D)                            
!                                                                       
      ELSEIF (ftyp (1:2) .eq.'BA') then 
         IF (.not.lni (ikfit) ) then 
            CALL setup_backpoly (ianz, werte, maxw) 
         ELSE 
            ier_num = - 25 
            ier_typ = ER_APPL 
         ENDIF 
!                                                                       
!------ Chebyshev Polynom (only 2D)                                     
!                                                                       
      ELSEIF (ftyp (1:2) .eq.'CH') then 
         IF (.not.lni (ikfit) ) then 
            CALL setup_poly_cheb (ianz, werte, maxw) 
         ELSE 
            ier_num = - 25 
            ier_typ = ER_APPL 
         ENDIF 
!                                                                       
!------ GSAS profile function(s) (only 2D)                              
!                                                                       
      ELSEIF (ftyp (1:2) .eq.'GS') then 
         IF (.not.lni (ikfit) ) then 
            CALL setup_gsas (ianz, werte, maxw, iname) 
         ELSE 
            ier_num = - 25 
            ier_typ = ER_APPL 
         ENDIF 
!                                                                       
!------ Gaussian (xy and xyz data sets)                                 
!                                                                       
      ELSEIF (ftyp (1:2) .eq.'GA') then 
         IF (.not.lni (ikfit) ) then 
            CALL setup_gauss (ianz, werte, maxw) 
         ELSE 
            CALL setup_gauss_2d (ianz, werte, maxw) 
         ENDIF 
!                                                                       
!------ Lorenzian (only 2D)                                             
!                                                                       
      ELSEIF (ftyp (1:2) .eq.'LO') then 
         IF (.not.lni (ikfit) ) then 
            CALL setup_lor (ianz, werte, maxw) 
         ELSE 
            ier_num = - 25 
            ier_typ = ER_APPL 
         ENDIF 
!                                                                       
!------ Pseudo-Voigt (only 2D)                                          
!                                                                       
      ELSEIF (ftyp (1:2) .eq.'PS') then 
         IF (.not.lni (ikfit) ) then 
            CALL setup_psvgt (ianz, werte, maxw) 
         ELSE 
            ier_num = - 25 
            ier_typ = ER_APPL 
         ENDIF 
!                                                                       
!------ User defined function (xy and xyz data sets)                    
!                                                                       
      ELSEIF (ftyp (1:2) .eq.'FX') then 
         CALL setup_user (ianz, werte, maxw, cpara, lpara) 
!                                                                       
!------ User defined macro (xy and xyz data sets)                    
!                                                                       
      ELSEIF (ftyp (1:2) .eq.'MA') then 
         CALL setup_user_macro (ianz, werte, maxw, cpara, lpara) 
!                                                                       
      ELSE 
         ier_num = - 25 
         ier_typ = ER_APPL 
      ENDIF 
!                                                                       
      END SUBROUTINE do_fit_fkt                     
!*****7*****************************************************************
      SUBROUTINE do_fit_macro (zeile, lp) 
!+                                                                      
!     fitparameter und einstellungen als macro speichern                
!-                                                                      
      USE errlist_mod 
      USE get_params_mod
      USE prompt_mod 
      USE kuplot_config 
      USE kuplot_mod 
!
      USE build_name_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 10) 
!                                                                       
      CHARACTER ( * ) zeile 
      CHARACTER(1024) cpara (maxw) 
      REAL werte (maxw) 
      INTEGER lpara (maxw), lp 
      INTEGER ianz, i 
!                                                                       
      INTEGER len_str 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
      CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
      IF (ier_num.ne.0) return 
!                                                                       
      CALL oeffne (77, cpara (1) , 'unknown') 
      IF (ier_num.ne.0) return 
!                                                                       
      WRITE (output_io, 3000) cpara (1) (1:len_str (cpara (1) ) ) 
      DO i = 1, npara 
      WRITE (77, 1000) i 
      WRITE (77, 2000) i, pinc (i), p (i) 
      ENDDO 
      CLOSE (77) 
!                                                                       
 1000 FORMAT ('# parameter number ',i3) 
 2000 FORMAT ('par ',i3,',',g12.6,',',f2.0) 
 3000 FORMAT (' ---------- > Saving current parameters to macro file ', &
     &                   a,' ..')                                       
      END SUBROUTINE do_fit_macro                   
!*****7*****************************************************************
      SUBROUTINE do_fit_pmerk 
!+                                                                      
!     speichern der fitparameter                                        
!-                                                                      
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER i 
!                                                                       
      DO i = 1, npara 
      p_merk (i) = p (i) 
      pinc_merk (i) = pinc (i) 
      ENDDO 
      END SUBROUTINE do_fit_pmerk                   
!*****7*****************************************************************
      SUBROUTINE do_fit_prueck 
!+                                                                      
!     speichern der fitparameter                                        
!-                                                                      
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER i 
!                                                                       
      DO i = 1, npara 
      p (i) = p_merk (i) 
      pinc (i) = pinc_merk (i) 
      dp (i) = 0.0 
      ENDDO 
      END SUBROUTINE do_fit_prueck                  
!*****7*****************************************************************
      SUBROUTINE do_fit_save 
!+                                                                      
!     speichern der fitergebnisse                                       
!-                                                                      
      USE errlist_mod 
      USE prompt_mod 
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 4) 
!                                                                       
      CHARACTER(80) filname 
      CHARACTER(2) cdummy 
      REAL werte (maxw) 
      INTEGER ianz 
      INTEGER len_str 
!                                                                       
      filname = fname (ikfit)(1:MIN(60,LEN_TRIM(fname(ikfit))))
      filname = filname (1:len_str (filname) ) //'.erg' 
!                                                                       
      CALL oeffne (77, filname, 'unknown') 
      IF (ier_num.ne.0) return 
!                                                                       
      CALL do_fit_info (77, .true., .true., .true.) 
      WRITE (output_io, 1000) filname (1:len_str (filname) ) 
      CLOSE (77) 
!                                                                       
      WRITE (output_io, 2000) fname (ikcal) (1:len_str (fname (ikcal) ) &
      )                                                                 
      cdummy = fform (ikcal) 
      ianz = 0 
      CALL check_form (cdummy, ianz, werte, maxw) 
      CALL do_save (ikcal, fname (ikcal), fform (ikcal), ianz, werte,   &
      maxw, .false.)                                                    
!                                                                       
      WRITE (output_io, 3000) fname (ikdif) (1:len_str (fname (ikdif) ) &
      )                                                                 
      cdummy = fform (ikdif) 
      ianz = 0 
      CALL check_form (cdummy, ianz, werte, maxw) 
      CALL do_save (ikdif, fname (ikdif), fform (ikdif), ianz, werte,   &
      maxw, .false.)                                                    
!                                                                       
 1000 FORMAT (' ---------- > Saving results in file ',a,' ..') 
 2000 FORMAT (' ---------- > Saving calculated data in file ',a,' ..') 
 3000 FORMAT (' ---------- > Saving difference data in file ',a,' ..') 
      END SUBROUTINE do_fit_save                    
!*****7*****************************************************************
      SUBROUTINE do_fit_wichtung (zeile, lp) 
!+                                                                      
!     aendern der wichtung                                              
!-                                                                      
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
      USE kuplot_config 
      USE kuplot_mod 
      USE string_convert_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 2) 
!                                                                       
      CHARACTER ( * ) zeile 
      CHARACTER(1024) cpara (maxw) 
      INTEGER lpara (maxw), lp 
      INTEGER ianz 
      REAL werte (maxw) 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
!                                                                       
      IF (ianz.eq.1) then 
         CALL do_cap (cpara (1) ) 
         IF (cpara (1) (1:3) .ne.'LOG'.and.cpara (1) (1:3)              &
         .ne.'SQR'.and.cpara (1) (1:3) .ne.'ONE'.and.cpara (1) (1:3)    &
         .ne.'LIN'.and.cpara (1) (1:3) .ne.'SQA'.and.cpara (1) (1:3)    &
         .ne.'INV'.and.cpara (1) (1:3) .ne.'BCK'.and.cpara (1) (1:3)    &
         .ne.'ISQ'.and.cpara (1) (1:3) .ne.'DAT'                        &
                   ) then                  
            ier_num = - 27 
            ier_typ = ER_APPL 
         ELSE 
            wtyp = cpara (1) (1:3) 
            wval = 0.01 
         ENDIF 
      ELSEIF (ianz.eq.2) then 
         CALL do_cap (cpara (1) ) 
         IF (cpara (1) (1:3) .eq.'BCK') then 
            wtyp = cpara (1) (1:3) 
            cpara (1) = '(0)' 
            lpara (1) = 3 
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
            IF (ier_num.eq.0) then 
               wval = werte (2) 
            ENDIF 
         ELSE 
            ier_num = - 27 
            ier_typ = ER_APPL 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
      END SUBROUTINE do_fit_wichtung                
!*****7*****************************************************************
      SUBROUTINE do_fit_par (zeile, lp) 
!+                                                                      
!     aendern der fit-parameter                                         
!-                                                                      
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
      USE prompt_mod 
      USE param_mod
      USE kuplot_config 
      USE kuplot_mod 
      USE string_convert_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 3) 
!                                                                       
      CHARACTER ( * ) zeile 
      CHARACTER(1024) cpara (maxw) 
      REAL werte (maxw) 
      INTEGER lpara (maxw), lp 
      INTEGER ianz, ip 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
!                                                                       
      IF (ianz.eq.0) then 
         CALL do_fit_info (output_io, .false., .false., .true.) 
!                                                                       
      ELSEIF (ianz.eq.1) then 
         CALL do_cap (cpara (1) ) 
         IF (cpara (1) (1:2) .eq.'SA') then 
            CALL do_fit_pmerk 
         ELSEIF (cpara (1) (1:2) .eq.'LO') then 
            CALL do_fit_prueck 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
!                                                                       
      ELSEIF (ianz.eq.2.or.ianz.eq.3) then 
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         IF (ier_num.ne.0) return 
         ip = nint (werte (1) ) 
         IF (ip.gt.0.and.ip.le.MAXPARA) then 
            pinc (ip) = werte (2) 
            IF (ianz.eq.3) p (ip) = werte (3) 
            WRITE (output_io, 1000) ip, p (ip), pinc (ip) 
         ELSE 
            ier_num = - 26 
            ier_typ = ER_APPL 
         ENDIF 
!                                                                       
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!
      IF(ier_num == 0) THEN 
         DO ip=1,MAXPARA
            kupl_para(ip) = p(ip)
         ENDDO
      ENDIF
!                                                                       
 1000 FORMAT     (' ---------- > Setting parameter p(',i2,') = ',g12.6, &
     &                   '  pinc = ',f4.1)                              
      END SUBROUTINE do_fit_par                     
!*****7*****************************************************************
      SUBROUTINE show_fit_para (idout) 
!+                                                                      
!     anzeigen der fit-parameter                                        
!-                                                                      
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER(30) fitfkt, wictyp 
      INTEGER idout, lt1, lt2, lf, lfn, lw, ipkt 
!                                                                       
      INTEGER len_str 
!                                                                       
      IF (ftyp (1:2) .eq.'PO') then 
         fitfkt = 'Polynom' 
      ELSEIF (ftyp (1:2) .eq.'BA') then 
         fitfkt = 'Background Polynom' 
      ELSEIF (ftyp (1:2) .eq.'CH') then 
         fitfkt = 'Chebyshev Polynom' 
      ELSEIF (ftyp (1:2) .eq.'FX') then 
         fitfkt = 'f = '//fit_func (1:fit_lfunc) 
      ELSEIF (ftyp (1:2) .eq.'LO') then 
         fitfkt = 'Lorenzian' 
      ELSEIF (ftyp (1:2) .eq.'PS') then 
         fitfkt = 'Pseudo-Voigt' 
      ELSEIF (ftyp (1:2) .eq.'GA') then 
         IF (lni (ikfit) ) then 
            fitfkt = 'Gaussian (2D)' 
         ELSE 
            fitfkt = 'Gaussian (1D)' 
         ENDIF 
      ELSEIF (ftyp (1:2) .eq.'MA') then 
         fitfkt = 'f = '//fit_func (1:fit_lfunc) 
      ELSE 
         fitfkt = 'not defined' 
      ENDIF 
!                                                                       
      IF (wtyp (1:3) .eq.'LOG') then 
         wictyp = 'w(i) = log(i)' 
      ELSEIF (wtyp (1:3) .eq.'SQR') then 
         wictyp = 'w(i) = sqrt(i)' 
      ELSEIF (wtyp (1:3) .eq.'ONE') then 
         wictyp = 'w(i) = 1.0' 
      ELSEIF (wtyp (1:3) .eq.'LIN') then 
         wictyp = 'w(i) = i' 
      ELSEIF (wtyp (1:3) .eq.'SQA') then 
         wictyp = 'w(i) = i**2' 
      ELSEIF (wtyp (1:3) .eq.'INV') then 
         wictyp = 'w(i) = 1.0/i' 
      ELSEIF (wtyp (1:3) .eq.'ISQ') then 
         wictyp = 'w(i) = 1.0/sqrt(i)' 
      ELSEIF (wtyp (1:3) .eq.'DAT'.and..not.lni (ikfit) ) then 
         wictyp = 'w(i) = 1/(dy(i))^2 from data ' 
      ELSEIF (wtyp (1:3) .eq.'BCK') then 
         wictyp = 'w(i) = exp(-WVAL*(Fobs-Fcalc))' 
      ELSE 
         wictyp = 'unknown' 
      ENDIF 
!                                                                       
      lt1 = max (1, len_str (titel (iwin, iframe, 1) ) ) 
      lt2 = max (1, len_str (titel (iwin, iframe, 2) ) ) 
      lf = len_str (fitfkt) 
      lfn = len_str (fname (ikfit) ) 
      lw = len_str (wictyp) 
!                                                                       
      IF (lni (ikfit) ) then 
         ipkt = nx (ikfit) * ny (ikfit) 
      ELSE 
         ipkt = len (ikfit) 
      ENDIF 
!                                                                       
      WRITE (idout, 1000) titel (iwin, iframe, 1) (1:lt1), titel (iwin, &
      iframe, 2) (1:lt2), fitfkt (1:lf), fname (ikfit) (1:lfn), ipkt,   &
      npara, urf, ncycle, fit_ifen, wictyp (1:lw)                       
!                                                                       
 1000 FORMAT (1x,'General fit parameter settings : ',/,                 &
     &        3x,'Title            : ',a,/11x,a,/,                      &
     &        3x,'Fit function     : ',a,/                              &
     &        3x,'Data file name   : ',a,/                              &
     &        3x,'# of data pts.   : ',i6,/                             &
     &        3x,'# of parameters  : ',i6,/                             &
     &        3x,'Urf              : ',f9.4,/                           &
     &        3x,'Max. cycle       : ',i6,/                             &
     &        3x,'MFEN for maxima  : ',i6,/                             &
     &        3x,'Weighting scheme : ',a,/)                             
      END SUBROUTINE show_fit_para                  
!*****7*****************************************************************
SUBROUTINE show_fit_erg(idout)
!(MAXP, MAXF, npara, nfixed, ndata, chisq, conf, lamda,  &
!           r4, re, &
!           params, pp, dpp, prange, cl, fixed, pf)
!+                                                                      
!     Display fit results
!-                                                                      
!
USE kuplot_fit_para
USE prompt_mod
!                                                                       
IMPLICIT none 
!
INTEGER, INTENT(IN) :: idout
!                                                                       
INTEGER :: i, j 
LOGICAL :: kor 
!
IF(kup_fit6_chisq>=0.0) THEN
   WRITE(idout, 1040) kup_fit6_chisq, kup_fit6_chisq/kup_fit6_ndata,            &
      kup_fit6_conf, kup_fit6_chisq/(kup_fit6_ndata-kup_fit6_npara),            &
      kup_fit6_ndata, kup_fit6_npara, kup_fit6_lamda, kup_fit6_r4, kup_fit6_re 
   WRITE(idout, 1050) 
   kor = .false. 
   DO i = 2, kup_fit6_npara 
      DO j = 1, i - 1 
         IF(ABS(kup_fit6_cl(i, j)) >=  0.8) THEN 
            WRITE(idout, 1070) i, j, kup_fit6_cl(i, j) 
            kor = .TRUE. 
         ENDIF 
      ENDDO 
   ENDDO 
   IF (.NOT.kor) WRITE(idout, 1060) 
ENDIF
WRITE(idout, * ) ' ' 
!
!DO i=1, npara
!   IF(prange(i,1)<=prange(i,2)) THEN
!      WRITE(output_io,1110) params(i), pp(i), dpp(i), prange(i,1), prange(i,2)
!   ELSE
!      WRITE(output_io,1100) params(i), pp(i), dpp(i)
!   ENDIF
!ENDDO
!DO i=1, nfixed
!   WRITE(output_io,1200) fixed(i), pf(i)
!ENDDO
!RITE(idout, * ) ' ' 
!                                                                       
 1040 FORMAT (/,                                                        &
              ' Information about the fit : ',/,                        &
              3x,'Chi^2      : ',g12.6, 5x,' Chi^2/N : ',g12.6 /,       &
              3x,'Conf. level: ',g12.6, 1x,' Chi^2/(N-P) : ',g12.6 /,   &
              3x,'No.Data    : ',i12  , 4x,' No.Params: ',i12/,          &
              3x,'MRQ final  : ',g12.6,/,                               &
              3x,'wR value   : ',g12.6, 5x,' R exp   : ',g12.6/)        
 1050 FORMAT (' Correlations larger than 0.8 :') 
 1060 FORMAT (3x,'** none **') 
 1070 FORMAT (3x,'Between p(',i2,') - p(',i2,') : ',f6.3) 
! 1100 FORMAT (3x,a16,' : ',g15.6E3, ' +- ',g15.6,4x)
! 1110 FORMAT (3x,a16,' : ',g15.6E3, ' +- ',g15.6,4x, '[',g15.6E3,',',g15.6E3,']')
! 1200 FORMAT (3x,a16,' : ',g15.6E3)
!                                                                       
END SUBROUTINE show_fit_erg                   
!
!*******************************************************************************
!!
!      SUBROUTINE show_fit_erg (idout) 
!!+                                                                      
!!     anzeigen der fit-ergebnisse                                       
!!-                                                                      
!      USE kuplot_config 
!      USE kuplot_mod 
!!                                                                       
!      IMPLICIT none 
!!                                                                       
!      INTEGER idout, i, j 
!      LOGICAL kor 
!!                                                                       
!      WRITE (idout, 1040) zalt, zwert, zdif, fend, r4, re 
!      WRITE (idout, 1050) 
!      kor = .false. 
!      DO i = 2, npara 
!      DO j = 1, i - 1 
!      IF (abs (cl (i, j) ) .gt.0.8) then 
!         WRITE (idout, 1070) i, j, cl (i, j) 
!         kor = .true. 
!      ENDIF 
!      ENDDO 
!      ENDDO 
!      IF (.not.kor) write (idout, 1060) 
!      WRITE (idout, * ) ' ' 
!!                                                                       
! 1040 FORMAT (' Information about the fit : ',/,                        &
!     &        3x,'Sum n-1    : ',g12.6,12x,' Sum n   : ',g12.6/,        &
!     &        3x,'Difference : ',g12.6,/,                               &
!     &        3x,'Urf final  : ',g12.6,/,                               &
!     &        3x,'R4 value   : ',g12.6,12x,' R exp   : ',g12.6/)        
! 1050 FORMAT (' Correlations larger than 0.8 :') 
! 1060 FORMAT (3x,'** none **') 
! 1070 FORMAT (3x,'Between p(',i2,') - p(',i2,') : ',f6.3) 
!!                                                                       
!      END SUBROUTINE show_fit_erg                   
!*****7*****************************************************************
      SUBROUTINE write_fit 
!+                                                                      
!     kupl.fit schreiben fuer textframe                                 
!-                                                                      
      USE errlist_mod 
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER i, j 
      LOGICAL kor 
!                                                                       
      CALL oeffne (22, 'kupl.fit', 'unknown') 
      IF (ier_num.ne.0) return 
!                                                                       
      WRITE (22, 1000) ftyp, r4 * 100., re * 100. 
      kor = .false. 
      DO i = 2, npara 
      DO j = 1, i - 1 
      IF (abs (cl (i, j) ) .gt.0.8) then 
         WRITE (22, 1070) i, j, cl (i, j) 
         kor = .true. 
      ENDIF 
      ENDDO 
      ENDDO 
      IF (.not.kor) write (22, 1060) 
      WRITE (22, 1100) 
      DO i = 1, npara 
      IF (pinc (i) .ne.1) then 
         WRITE (22, 1200) i, p (i) 
      ELSE 
         WRITE (22, 1210) i, p (i), dp (i) 
      ENDIF 
      ENDDO 
      CLOSE (22) 
!                                                                       
 1000 FORMAT (/1x,'F i t  -  r e s u l t s',/,                          &
     &         1x,'--------------------------------------',/,           &
     &         1x,'Fit function : ',a4,/,                               &
     &         1x,'R value      : ',f5.1,' %',/                         &
     &         1x,'Rexp value   : ',f5.1,' %',/                         &
     &         1x,'--------------------------------------',/,           &
     &         1x,'Correlations > 0.8 : ',/)                            
 1060 FORMAT ( 1x,'** none **') 
 1070 FORMAT ( 1x,'betw. p(',i2,') - p(',i2,') : ',f6.3) 
 1100 FORMAT ( 1x,'--------------------------------------',/,           &
     &         1x,'Resulting parameters: ',/)                           
 1200 FORMAT ( 1x,'p(',i2,') = ',g12.6,' fixed') 
 1210 FORMAT ( 1x,'p(',i2,') = ',g12.6,' +- ',g12.6) 
!                                                                       
      CLOSE (22) 
!                                                                       
      END SUBROUTINE write_fit                      
!*****7*****************************************************************
      SUBROUTINE do_fit_info (idout, f_se, f_er, f_pa) 
!+                                                                      
!     ausgabe von fitinformationen                                      
!-                                                                      
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER idout 
      LOGICAL f_se, f_er, f_pa 
!                                                                       
      IF (f_se) call show_fit_para (idout) 
      IF (f_er) call show_fit_erg (idout) 
!                                                                       
      IF (f_pa) then 
         IF (ftyp (1:2) .eq.'PO') then 
            CALL show_poly (idout) 
         ELSEIF (ftyp (1:2) .eq.'BA') then 
            CALL show_backpoly (idout) 
         ELSEIF (ftyp (1:2) .eq.'GS') then 
            CALL show_gsas (idout) 
         ELSEIF (ftyp (1:2) .eq.'CH') then 
            CALL show_poly_cheb (idout) 
         ELSEIF (ftyp (1:2) .eq.'GA') then 
            IF (.not.lni (ikfit) ) then 
               CALL show_gauss (idout) 
            ELSE 
               CALL show_gauss_2d (idout) 
            ENDIF 
         ELSEIF (ftyp (1:2) .eq.'LO') then 
            CALL show_lor (idout) 
         ELSEIF (ftyp (1:2) .eq.'PS') then 
            CALL show_psvgt (idout) 
         ELSEIF (ftyp (1:2) .eq.'FX') then 
            CALL show_user (idout) 
         ELSEIF (ftyp (1:2) .eq.'MA') then 
            CALL show_user_macro (idout) 
         ENDIF 
      ENDIF 
!                                                                       
      END SUBROUTINE do_fit_info                    
!*****7*****************************************************************
!
SUBROUTINE do_fit_y 
!+                                                                      
!   Transfer routine to kuplot_mrq
!   Copies the data set into the new data, calls the fit routine
!   kuplot_mrq and copies the fit and difference into the 
!   last kuplot curves
!     der eigentliche fit fuer xy-files                                 
!-                                                                      
USE prompt_mod 
USE kuplot_config 
USE kuplot_mod 
USE kuplot_fit_para
!                                                                       
IMPLICIT none 
!                                                                       
CHARACTER(LEN=60) :: filname 
!REAL              :: xx, f, df (maxpara) 
INTEGER           :: i, ii, jj, kk 
INTEGER           :: j
!INTEGER           :: len_str 
!
INTEGER           :: MAXP     ! Maximum number of parameters 
INTEGER           :: MAXf     ! MAXIMUM number of fixed parameters 
INTEGER           :: nparams  ! number of refined parameter, NPARA is para number from KUPLOT
INTEGER           :: nfixed   ! number of fixed parameter
INTEGER           :: kupl_last! Last KUPLOT data that are needed
CHARACTER(LEN=60), DIMENSION(:), ALLOCATABLE   :: par_names 
CHARACTER(LEN=60), DIMENSION(:), ALLOCATABLE   :: fixed
INTEGER          , DIMENSION(:), ALLOCATABLE   :: par_ind      ! Index of refined parameter
INTEGER          , DIMENSION(:), ALLOCATABLE   :: fixed_ind    ! Index of fixed parameter
INTEGER          , DIMENSION(2)                :: data_dim     ! No of data points
REAL             , DIMENSION(:,:), ALLOCATABLE :: data_data    ! The data
REAL             , DIMENSION(:,:), ALLOCATABLE :: data_weight  ! The weights
REAL             , DIMENSION(:)  , ALLOCATABLE :: data_x       ! X-coordinates
REAL             , DIMENSION(:)  , ALLOCATABLE :: data_y       ! y-coordinates(for xyz data
REAL             , DIMENSION(:,:), ALLOCATABLE :: data_calc    ! The calculated data
REAL                                           :: conv_dp_sig  ! Max Dp/sigma 
REAL                                           :: conv_dchi2   ! Max DELTA Chi2
REAL                                           :: conv_chi2    ! Max Chi2
REAL                                           :: conv_conf    ! Min confidence level
LOGICAL                                        :: lconvergence ! Convergence was reached
REAL                                           :: chisq        ! Final Chi2
REAL                                           :: conf         ! Final confidence level
REAL                                           :: lamda_fin    ! Final LevenbergMarq Lamda
REAL                                           :: rval         ! Final weighted rvalue
REAL                                           :: rexp         ! R expected
REAL             , DIMENSION(:,:), ALLOCATABLE :: prange       ! Allowed parameter range
REAL             , DIMENSION(:)  , ALLOCATABLE :: par_value    ! Parameters
REAL             , DIMENSION(:)  , ALLOCATABLE :: dpp          ! Parameter uncertainty
REAL             , DIMENSION(:,:), ALLOCATABLE :: covar        ! Covariance matrix
REAL             , DIMENSION(:)  , ALLOCATABLE :: pf           ! Parameters fixed
!
INTEGER :: ifree, ifixed  ! temporary variables Number of free and fixed params
!                                                                       
CALL wichtung (y) 
!
data_dim(1) = len(ikfit)                        ! Transfer KUPLOT Dimensions
data_dim(2) = 1
ALLOCATE(data_calc(  data_dim(1), data_dim(2)))   ! Allocate CALC  at proper dimensions
ALLOCATE(data_data(  data_dim(1), data_dim(2)))   ! Allocate DATA  at proper dimensions
ALLOCATE(data_weight(data_dim(1), data_dim(2)))   ! Allocate SIGMA at proper dimensions
ALLOCATE(data_x     (data_dim(1))             )   ! Allocate x     at proper dimensions
ALLOCATE(data_y     (data_dim(2))             )   ! Allocate y     at proper dimensions
DO i=1, data_dim(1)
   data_x(i)       = x(offxy(ikfit-1) + i)
   data_data(i, 1) = y(offxy(ikfit-1) + i)
   IF(dy(offxy(ikfit-1) + i) /= 0.0) THEN
      data_weight(i,1) = dy(offxy(ikfit-1) + i)
   ELSE
      data_weight(i,1) = 1.0
   ENDIF
ENDDO

data_y(1:data_dim(2)) = 1.0                       ! Set dummy y

conv_dp_sig   = 0.005
conv_dchi2    = 0.5
conv_chi2     = 0.5
conv_conf     = 0.001
lconvergence = .FALSE.
!
! Determine and sort refined/fixed parameters
!
MAXP = npara
MAXF = 1
nfixed = 0
DO i=1, MAXP
   IF(pinc(i)==0) nfixed = nfixed + 1
ENDDO
nparams = NPARA - nfixed               ! Adjust number of free parameters
MAXP    = NPARA - nfixed               ! Adjust number of free parameters
MAXF    = MAX(1, nfixed)               ! Adjust number of fixed parameters
write(*,*) ' DO_FIT_Y', npara, MAXP, nparams, MAXF, nfixed
ALLOCATE(par_names(MAXP))
ALLOCATE(par_ind  (MAXP))
ALLOCATE(par_value(MAXP))
ALLOCATE(dpp      (MAXP))
ALLOCATE(prange   (MAXP,2))
ALLOCATE(covar    (MAXP, MAXP))
ALLOCATE(fixed    (MAXF))
ALLOCATE(fixed_ind(MAXF))
ALLOCATE(pf       (MAXF))
ifree   = 0
ifixed  = 0
DO i=1, NPARA                          ! Loop over all KUPLOT parameters
   IF(pinc(i)/=0.0) THEN               ! Refined parameter
      ifree            = ifree + 1
      par_value(ifree) = p(i)          ! Copy KUPLOT value
      WRITE(par_names(ifree),'(a,i2.2,a)') 'p[',i,']'
      par_ind(ifree)   = i             ! Set index
      prange(ifree,1)  = 1.0
      prange(ifree,2)  =-1.0
   ELSE                                ! Fixed parameter
      ifixed           = ifixed + 1
      pf(ifixed)       = p(i)          ! Copy KUPLOT value
      WRITE(fixed(ifixed),'(a,i2.2,a)') 'p[',i,']'
      fixed_ind(ifixed) = i            ! Set index
   ENDIF
ENDDO
!
dpp(:)     = 0.0
covar(:,:) = 0.0
kupl_last  = iz-1
CALL kuplot_mrq(MAXP, nparams, ncycle, kupl_last, par_names, par_ind,           &
                MAXF, nfixed, fixed, fixed_ind, pf, &
                data_dim, &
                      data_data, data_weight, data_x, data_y,              &
                      data_calc,                                           &
                      conv_dp_sig, conv_dchi2, conv_chi2, conv_conf,       &
                      lconvergence,                                        &
                      chisq, conf, lamda_fin, rval, rexp, par_value, prange, dpp, covar)
!OLD  IF (ncycle.gt.0) call fit_kupl (y) 
!
! Copy values back into KUPLOT scheme
dp(:) = 0.0
DO i=1, nparams
   ii    = par_ind(i)
   p(ii)  = par_value(i) 
   dp(ii) = dpp(i) 
   DO j=1, nparams
      jj      = par_ind(j)
      cl(ii,jj) = covar(i,j)
   ENDDO
ENDDO
!                                                                       
ii = offxy (ikfit - 1) 
jj = offxy (ikcal - 1) 
kk = offxy (ikdif - 1) 
!                                                                       
DO i=1, len(ikfit)
   x(jj+i) = data_x(i)          ! Copy x-values into result
   y(jj+i) = data_calc(i,1)     ! Copy calculated data into result
   dx(jj+i) = 0.0
   dy(jj+i) = 0.0
!
   x(kk+i)  = data_x(i)          ! Copy x-values into difference curve
   y(kk+i)  = y(ii+i) - data_calc(i, 1)  ! Calc difference
   dx(kk+i) = 0.0
   dy(kk+i) = 0.0
ENDDO
!                                                                       
len (ikcal) = len (ikfit) 
len (ikdif) = len (ikfit) 
fform (ikcal) = fform (ikfit) 
fform (ikdif) = fform (ikfit) 
filname = fname (ikfit)(1:MIN(60,LEN_TRIM(fname(ikfit)))) 
fname (ikcal) = filname (1:LEN_TRIM(filname) ) //'.fit' 
fname (ikdif) = filname (1:LEN_TRIM(filname) ) //'.dif' 
CALL get_extrema 
!

CALL kup_fit6_set(MAXP, MAXF, nparams, nfixed, data_dim(1)*data_dim(2), chisq, conf, lamda_fin,  &
                  rval, rexp, par_names, par_value, dpp, prange, covar, fixed, pf)
!                                                                       
CALL do_fit_info (output_io, .false., .false., .true.) 
!
DEALLOCATE(data_calc)
DEALLOCATE(data_data)
DEALLOCATE(data_weight)
DEALLOCATE(data_x)
DEALLOCATE(data_y)
DEALLOCATE(par_names)
DEALLOCATE(par_ind)
DEALLOCATE(par_value)
DEALLOCATE(dpp)
DEALLOCATE(prange)
DEALLOCATE(covar)
DEALLOCATE(fixed)
DEALLOCATE(fixed_ind)
DEALLOCATE(pf   )
!
END SUBROUTINE do_fit_y                       
!
!*****7*****************************************************************
      SUBROUTINE do_fit_z 
!+                                                                      
!     der eigentliche fit fuer xyz-files                                
!-                                                                      
      USE prompt_mod 
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER(60) filname 
      REAL xx, f, df (maxpara) 
      INTEGER i, iii 
      INTEGER len_str 
!                                                                       
      CALL wichtung (z) 
      IF (ncycle.gt.0) call fit_kupl (z) 
!                                                                       
      DO i = 1, nx (ikfit) * ny (ikfit) 
      xx = REAL(i) 
!      CALL kupl_theory (xx, f, df, - i) 
      z (offz (ikcal - 1) + i) = f 
      IF (z (offz (ikfit - 1) + i) .ne. - 9999) then 
         z (offz (ikdif - 1) + i) = z (offz (ikfit - 1) + i) - f 
      ELSE 
         z (offz (ikdif - 1) + i) = - 9999.0 
      ENDIF 
      ENDDO 
      DO iii = 1, nx (ikfit) 
      x (offxy (ikcal - 1) + iii) = x (offxy (ikfit - 1) + iii) 
      x (offxy (ikdif - 1) + iii) = x (offxy (ikfit - 1) + iii) 
      ENDDO 
      DO iii = 1, ny (ikfit) 
      y (offxy (ikcal - 1) + iii) = y (offxy (ikfit - 1) + iii) 
      y (offxy (ikdif - 1) + iii) = y (offxy (ikfit - 1) + iii) 
      ENDDO 
!                                                                       
      lni (ikcal) = .true. 
      lni (ikdif) = .true. 
      len (ikcal) = len (ikfit) 
      len (ikdif) = len (ikfit) 
      nx (ikcal) = nx (ikfit) 
      ny (ikcal) = ny (ikfit) 
      nx (ikdif) = nx (ikfit) 
      ny (ikdif) = ny (ikfit) 
      fform (ikcal) = fform (ikfit) 
      fform (ikdif) = fform (ikfit) 
!                                                                       
      filname = fname (ikfit)(1:MIN(60,LEN_TRIM(fname(ikfit)))) 
      fname (ikcal) = filname (1:len_str (filname) ) //'.fit' 
      fname (ikdif) = filname (1:len_str (filname) ) //'.dif' 
      CALL get_extrema 
!                                                                       
      CALL do_fit_info (output_io, .false., .false., .true.) 
      END SUBROUTINE do_fit_z                       
!*****7*****************************************************************
      SUBROUTINE wichtung (a) 
!+                                                                      
!     Calculation of weights. Values outside plotting range are         
!     set to zero if frall is .false.                                   
!-                                                                      
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      REAL a (maxarray) 
      INTEGER i, ii 
!                                                                       
      REAL calc_wic 
!                                                                       
      IF (lni (ikfit) ) then 
         ii = offz (ikfit - 1) 
         DO i = 1, nx (ikfit) * ny (ikfit) 
         w (ii + i) = calc_wic (a (ii + i), dy (ii + i) ) 
         ENDDO 
      ELSE 
         ii = offxy (ikfit - 1) 
         DO i = 1, len (ikfit) 
         IF (frall) then 
            w (ii + i) = calc_wic (a (ii + i), dy (ii + i) ) 
         ELSE 
            IF (x (ii + i) .lt.ex (iwin, iframe, 1) .or.x (ii + i)      &
            .gt.ex (iwin, iframe, 2) ) then                             
               w (ii + i) = 0.0 
            ELSE 
               w (ii + i) = calc_wic (a (ii + i), dy (ii + i) ) 
            ENDIF 
         ENDIF 
         ENDDO 
      ENDIF 
!                                                                       
!                                                                       
      END SUBROUTINE wichtung                       
!*****7*****************************************************************
      REAL function calc_wic (val, sig) 
!+                                                                      
!     Calculation of weights.                                           
!-                                                                      
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      REAL val, aval, sig, wic 
!                                                                       
      wic = 0.0 
      IF (val.ne. - 9999.0) then 
         aval = abs (val) 
         IF (wtyp (1:3) .eq.'LOG') then 
            IF (aval.gt.0.0) wic = log (aval) 
         ELSEIF (wtyp (1:3) .eq.'SQR') then 
            wic = sqrt (aval) 
         ELSEIF (wtyp (1:3) .eq.'ONE') then 
            wic = 1.0 
         ELSEIF (wtyp (1:3) .eq.'LIN') then 
            wic = aval 
         ELSEIF (wtyp (1:3) .eq.'SQA') then 
            wic = aval**2 
         ELSEIF (wtyp (1:3) .eq.'INV') then 
            IF (aval.ne.0.0) wic = 1.0 / aval 
         ELSEIF (wtyp (1:3) .eq.'ISQ') then 
            wic = 1.0 / sqrt (aval) 
         ELSEIF (wtyp (1:3) .eq.'DAT'.and..not.lni (ikfit) ) then 
            wic = 1.0 / (sig *sig)
         ENDIF 
      ENDIF 
!                                                                       
      calc_wic = wic 
!                                                                       
      END FUNCTION calc_wic                         
!*****7*****************************************************************
!      SUBROUTINE kupl_theory (xx, f, df, iwert) 
!!                                                                       
!      USE kuplot_config 
!      USE kuplot_mod 
!!                                                                       
!      IMPLICIT none 
!!                                                                       
!      REAL xx, f, df (maxpara) 
!      INTEGER iwert 
!!                                                                       
!      IF (ftyp (1:2) .eq.'PO') then 
!DONE     CALL theory_poly (xx, f, df, iwert) 
!      ELSEIF (ftyp (1:2) .eq.'BA') then 
!         CALL theory_backpoly (xx, f, df, iwert) 
!      ELSEIF (ftyp (1:2) .eq.'CH') then 
!         CALL theory_poly_cheb (xx, f, df, iwert) 
!      ELSEIF (ftyp (1:2) .eq.'GS') then 
!         CALL theory_gsas (xx, f, df, iwert) 
!      ELSEIF (ftyp (1:2) .eq.'LO') then 
!         CALL theory_lor (xx, f, df, iwert) 
!      ELSEIF (ftyp (1:2) .eq.'PS') then 
!         CALL theory_psvgt (xx, f, df, iwert) 
!      ELSEIF (ftyp (1:2) .eq.'FX') then 
!         CALL theory_user (xx, f, df, iwert) 
!      ELSEIF (ftyp (1:2) .eq.'GA') then 
!         IF (.not.lni (ikfit) ) then 
!            CALL theory_gauss (xx, f, df, iwert) 
!         ELSE 
!            CALL theory_gauss_2d (xx, f, df, iwert) 
!         ENDIF 
!      ELSEIF (ftyp (1:2) .eq.'MA') then 
!         CALL theory_macro (xx, f, df, iwert) 
!      ENDIF 
!!                                                                       
!      END SUBROUTINE kupl_theory                    
!*****7*****************************************************************
!       User defined fit function                                       
!*****7*****************************************************************
SUBROUTINE setup_user (ianz, werte, maxw, cpara, lpara) 
!                                                                       
      USE  berechne_mod
      USE errlist_mod 
      USE kuplot_config 
      USE kuplot_mod 
USE kuplot_fit6_set_theory
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
!                                                                       
      CHARACTER ( * ) cpara (maxw) 
      CHARACTER(1024) cdummy 
      REAL dummy, werte (maxw) 
      INTEGER lpara (maxw) 
      INTEGER ianz, ip , length
!                                                                       
!                                                                       
      IF (ianz.eq.2) then 
         ip = nint (werte (1) ) 
         IF (ip.lt.1.or.ip.gt.maxpara) then 
            ier_num = - 31 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
!                                                                       
         npara = ip 
         fit_func = cpara (2) (1:lpara (2) ) 
         fit_lfunc = lpara (2) 
         cdummy = '('//fit_func (1:fit_lfunc) //')' 
         length = fit_lfunc + 2
         dummy = berechne (cdummy, length)
!                                                                       
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!
p_kuplot_theory => theory_user
!                                                                       
END SUBROUTINE setup_user                     
!
!*****7*****************************************************************
!       User defined fit macro                                       
!*****7*****************************************************************
      SUBROUTINE setup_user_macro (ianz, werte, maxw, cpara, lpara) 
!                                                                       
      USE  berechne_mod
      USE errlist_mod 
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER , INTENT(IN) :: maxw 
      INTEGER , INTENT(IN) :: ianz 
      REAL    , DIMENSION(MAXW), INTENT(IN) ::  werte (maxw) 
!                                                                       
      CHARACTER (LEN=*), DIMENSION(MAXW), INTENT(IN) :: cpara
      INTEGER          , DIMENSION(MAXW), INTENT(IN) :: lpara
      INTEGER :: ip 
      INTEGER :: m,i,iiw,iix
REAL, DIMENSION(maxpara) :: df
REAL                   :: f
      REAL :: xx
!                                                                       
!                                                                       
      IF (ianz.eq.2) then 
         ip = nint (werte (1) ) 
         IF (ip.lt.1.or.ip.gt.maxpara) then 
            ier_num = - 31 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
!                                                                       
         npara = ip 
         fit_func = cpara (2) (1:lpara (2) ) 
         fit_lfunc = lpara (2) 
!        cdummy = '('//fit_func (1:fit_lfunc) //')' 
!        length = fit_lfunc + 2
!        dummy = berechne (cdummy, length)
         IF (lni (ikfit) ) then 
            iiw = offz (ikfit - 1) 
            iix = offxy (ikfit - 1) 
            m = nx (ikfit) * ny (ikfit) 
         ELSE 
            iiw = offxy (ikfit - 1) 
            iix = offxy (ikfit - 1) 
            m = len (ikfit) 
         ENDIF 
         DO i = 1, 1  !!! m 
            xx = x (iix + i) 
            CALL theory_macro (xx, f, df, i) 
            IF(ier_num/=0) RETURN
         ENDDO
!                                                                       
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
      END SUBROUTINE setup_user_macro
!*****7*****************************************************************
      SUBROUTINE show_user (idout) 
!                                                                       
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER idout, i 
!                                                                       
      WRITE (idout, 1000) fit_func (1:fit_lfunc) 
      DO i = 1, npara 
      WRITE (idout, 1100) i, p (i), dp (i), pinc (i) 
      ENDDO 
      WRITE (idout, * ) ' ' 
!                                                                       
 1000 FORMAT (1x,'Fit function : F = ',a,/) 
 1100 FORMAT (3x,'p(',i2,') : ',g12.6,' +- ',g12.6,4x,'pinc : ',f2.0) 
!                                                                       
      END SUBROUTINE show_user                      
!*****7*****************************************************************
      SUBROUTINE show_user_macro (idout) 
!                                                                       
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER idout, i 
!                                                                       
      WRITE (idout, 1000) fit_func (1:fit_lfunc) 
      DO i = 1, npara 
      WRITE (idout, 1100) i, p (i), dp (i), pinc (i) 
      ENDDO 
      WRITE (idout, * ) ' ' 
!                                                                       
 1000 FORMAT (1x,'Fit function : F = ',a,/) 
 1100 FORMAT (3x,'p(',i2,') : ',g12.6,' +- ',g12.6,4x,'pinc : ',f2.0) 
!                                                                       
      END SUBROUTINE show_user_macro                      
!***********************************************************************
!
SUBROUTINE theory_user(MAXP, ix, iy, xx, yy, NPARA, params, par_names,          &
                       prange, l_do_deriv, data_dim, &
                       data_data, data_weight, data_x, data_y, &
                       data_calc, kupl_last,      &
                       ymod, dyda, LDERIV)
!
! Calculate a user defined function
! ymod = SUM p[i]*xx^ind
!
USE kuplot_fit_const
!
USE berechne_mod
USE matrix_mod
USE param_mod 
!
IMPLICIT NONE
!
INTEGER                                              , INTENT(IN)  :: MAXP    ! Parameter array sizes
INTEGER                                              , INTENT(IN)  :: ix      ! Point number along x
INTEGER                                              , INTENT(IN)  :: iy      ! Point number along y
REAL                                                 , INTENT(IN)  :: xx      ! Point value  along x
REAL                                                 , INTENT(IN)  :: yy      ! Point value  along y
INTEGER                                              , INTENT(IN)  :: NPARA   ! Number of refined parameters
REAL            , DIMENSION(MAXP )                   , INTENT(IN)  :: params  ! Parameter values
CHARACTER(LEN=*), DIMENSION(MAXP)                    , INTENT(IN)  :: par_names    ! Parameter names
REAL            , DIMENSION(MAXP, 2                 ), INTENT(IN)  :: prange      ! Allowed parameter range
LOGICAL         , DIMENSION(MAXP )                   , INTENT(IN)  :: l_do_deriv  ! Parameter needs derivative
INTEGER         , DIMENSION(2)                       , INTENT(IN)  :: data_dim     ! Data array dimensions
REAL            , DIMENSION(data_dim(1), data_dim(2)), INTENT(IN)  :: data_data    ! Data array
REAL            , DIMENSION(data_dim(1), data_dim(2)), INTENT(IN)  :: data_weight  ! Data sigmas
REAL            , DIMENSION(data_dim(1))             , INTENT(IN)  :: data_x       ! Data coordinates x
REAL            , DIMENSION(data_dim(1))             , INTENT(IN)  :: data_y       ! Data coordinates y
REAL            , DIMENSION(data_dim(1), data_dim(2)), INTENT(OUT) :: data_calc    ! Data array
INTEGER                                              , INTENT(IN)  :: kupl_last    ! Last KUPLOT DATA that are needed
REAL                                                 , INTENT(OUT) :: ymod    ! Function value at (ix,iy)
REAL            , DIMENSION(NPARA)                   , INTENT(OUT) :: dyda    ! Function derivatives at (ix,iy)
LOGICAL                                              , INTENT(IN)  :: LDERIV  ! TRUE if derivative is needed
!     SUBROUTINE theory_user (xx, f, df, i) 
!                                                                       
!     USE kuplot_config 
!     USE kuplot_mod 
!                                                                       
!     IMPLICIT none 
!                                                                       
REAL, PARAMETER      :: SCALEF = 0.05 ! Scalefactor for parameter modification 0.01 is good?
!
CHARACTER(LEN=1024) :: cdummy 
!     REAL xx, f, df (maxpara) 
!     REAL pb, h, err 
      INTEGER :: ind 
      INTEGER :: length
!
INTEGER              :: nder          ! Numper of points for derivative
REAL                 :: delta         ! Shift to calculate derivatives
!REAL                 :: p_d           ! Shifted  parameter
REAL, DIMENSION(3)   :: dvec          ! Parameter at P, P+delta and P-delta
REAL, DIMENSION(3)   :: yvec          ! Chisquared at P, P+delta and P-delta
REAL, DIMENSION(3)   :: avec          ! Params for derivative y = a + bx + cx^2
REAL, DIMENSION(3,3) :: xmat          ! Rows are : 1, P, P^2
REAL, DIMENSION(3,3) :: imat          ! Inverse to xmat
!                                                                       
      REAL, EXTERNAL :: dfridr
!                                                                       
DO ind = 1, npara 
   dyda (ind) = 0.0 
   kupl_para(ind) = params(ind)
ENDDO 
!                                                                       
!-------x und y bestimmen                                               
!                                                                       
!     IF (lni (ikfit) ) then 
!        ix = (abs (i) - 1) / ny (ikfit) + 1 
!        iy = abs (i) - (ix - 1) * ny (ikfit) 
!        rpara (0) = x (offxy (ikfit - 1) + ix) 
!        rpara (1) = y (offxy (ikfit - 1) + iy) 
!     ELSE 
!        rpara (0) = xx 
!     ENDIF 
!
rpara(0) = data_x(ix)
rpara(1) = data_y(iy)
!                                                                       
!     cdummy = '('//fit_func (1:fit_lfunc) //')' 
cdummy = '('//f6_fit_func(1:f6_fit_lfunc) //')' 
length = f6_fit_lfunc + 2
ymod   = berechne (cdummy, length)
!                                                                       
!-------ableitungen                                                     
!                                                                       
IF(LDERIV) THEN
   DO ind = 1, npara 
      IF(l_do_deriv(ind)) THEN
!!!!
!
         nder = 1                     ! First point is at P(IND)
         dvec(2) = 0.0
         dvec(3) = 0.0
         dvec(1) = params(ind)             ! Store parameter value
         yvec(1) = ymod                    ! Store function at P(IND)
         IF(params(ind)/=0.0) THEN
            delta = params(IND)*SCALEF     ! A multiplicative variation of the parameter seems best
         ELSE
            delta = 0.010
         ENDIF
!                                     ! Test at P + DELTA
         IF(prange(ind,1)<=prange(ind,2)) THEN     ! User provided parameter range
            kupl_para(ind)     = MIN(prange(ind,2),MAX(prange(ind,1),params(ind)+delta))
         ELSE
            kupl_para(ind)     = params(ind) + delta
         ENDIF
!        p_d = p(k) + delta
         IF(params(ind)/=kupl_para(ind)) THEN           ! Parameter is not at edge of range
            nder = nder + 1
            dvec(2) = kupl_para(ind)            ! Store parameter value at P+DELTA
            yvec(2) = berechne (cdummy, length)
!           CALL refine_set_param(NPARA, par_names(k), k, p_d )  ! Set modified value
!           CALL refine_macro(MAXP, refine_mac, refine_mac_l, NPARA, kupl_last, par_names, p, &
!                             data_dim, refine_temp)
            IF(ier_num /= 0) RETURN
!           DO iiy=1, data_dim(2)
!              DO iix=1, data_dim(1)
!                 refine_derivs(iix,iiy,k) =  refine_temp(iix,iiy)
!              ENDDO
!           ENDDO
         ENDIF
!                                     ! Test at P - DELTA
         IF(prange(ind,1)<=prange(ind,2)) THEN
            kupl_para(ind)     = MIN(prange(ind,2),MAX(prange(ind,1),params(ind)-delta))
         ELSE
            kupl_para(ind)     = params(ind) - delta
         ENDIF
!        p_d = p(k) - delta
         IF(params(ind)/=kupl_para(ind)) THEN           ! Parameter is not at edge of range
            nder = nder + 1
            dvec(3) = kupl_para(ind)            ! Store parameter value at P-DELTA
            yvec(3) = berechne (cdummy, length)
!           CALL refine_set_param(NPARA, par_names(k), k, p_d )  ! Set modified value
!           CALL refine_macro(MAXP, refine_mac, refine_mac_l, NPARA, kupl_last, par_names, p, &
!                             data_dim, refine_temp)
            IF(ier_num /= 0) RETURN
!
         ENDIF
!        CALL refine_set_param(NPARA, par_names(k), k, p(k))  ! Return to original value
         kupl_para(ind) = params(ind) ! Return to original value
!
         IF(nder==3) THEN             ! Got all three points for derivative
!           DO iiy=1, data_dim(2)
!              DO iix=1, data_dim(1)
!
!              Derivative is calculated as a fit of a parabola at P, P+delta, P-delta
!                 yvec(1) = refine_calc  (iix, iiy)
!                 yvec(2) = refine_derivs(iix, iiy,k)
!                 yvec(3) = refine_temp  (iix, iiy)
                  xmat(:,1) =  1.0
                  xmat(1,2) =  dvec(1) !p(IND)
                  xmat(2,2) =  dvec(2) !p(IND) + delta
                  xmat(3,2) =  dvec(3) !p(k) - delta
                  xmat(1,3) = (dvec(1))**2 !(p(IND)        ) **2
                  xmat(2,3) = (dvec(2))**2 !(p(IND) + delta) **2
                  xmat(3,3) = (dvec(3))**2 !(p(IND) - delta) **2
                  CALL matinv3(xmat, imat)
                  avec = MATMUL(imat, yvec)
                  dyda(ind) = avec(2) + 2.*avec(3)*params(ind)
!
!                 refine_derivs(iix, iiy, k) = avec(2) + 2.*avec(3)*p(k)
!              ENDDO
!           ENDDO
         ELSEIF(nder==2) THEN
            IF(dvec(2)==0) THEN          ! P + Delta failed
               dyda(ind) = (yvec(3)-yvec(1))/(dvec(3)-dvec(1))
!              DO iiy=1, data_dim(2)
!                 DO iix=1, data_dim(1)
!                    refine_derivs(iix, iiy, k) = (refine_temp  (iix, iiy)-refine_calc  (iix, iiy))/ &
!                                                 (dvec(3)-dvec(1))
!                 ENDDO
!              ENDDO
            ELSEIF(dvec(3)==0) THEN      ! P - Delta failed
               dyda(ind) = (yvec(2)-yvec(1))/(dvec(2)-dvec(1))
!              DO iiy=1, data_dim(2)
!                 DO iix=1, data_dim(1)
!                    refine_derivs(iix, iiy, k) = (refine_derivs(iix, iiy,k)-refine_calc  (iix, iiy))/ &
!                                                 (dvec(2)-dvec(1))
!                 ENDDO
!              ENDDO
            ENDIF
!
         ELSE
            ier_num = -6
            ier_typ = ER_APPL
            ier_msg(1) = par_names(ind)
            RETURN
         ENDIF
!
!!!!
      ENDIF
   ENDDO
ENDIF
!     IF (i.gt.0) then 
!        DO ind = 1, npara 
!        IF (pinc (ind) .ne.0) then 
!           h = 0.1 * abs (p (ind) ) 
!           IF (h.eq.0) h = 0.5 
!           pb = p (ind) 
!           np1 = ind 
!           df (ind) = dfridr (p (ind), h, err) 
!           p (ind) = pb 
!        ENDIF 
!        ENDDO 
!     ENDIF 
!                                                                       
END SUBROUTINE theory_user                    
!
!*****7*****************************************************************
      REAL function func (xx) 
!                                                                       
      USE  berechne_mod
      USE param_mod 
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER(LEN=1024) :: cdummy 
      INTEGER :: length
      REAL xx
!                                                                       
      p (np1) = xx 
      cdummy = '('//fit_func (1:fit_lfunc) //')' 
      length = fit_lfunc + 2
      func = berechne (cdummy, length)
!                                                                       
      END FUNCTION func                             
!*****7*****************************************************************
!       GSAS profile functions                                          
!*****7*****************************************************************
      SUBROUTINE show_gsas (idout) 
!                                                                       
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER(40) cout (maxpara) 
      INTEGER idout, i, j, k 
!                                                                       
      DO j = 1, np2 
      cout ( (j - 1) * np3 + 1) = 'Peak position' 
      cout ( (j - 1) * np3 + 2) = 'Intensity' 
      DO i = 3, npara - 2 
      WRITE (cout ( (j - 1) * np3 + i), 2000) i - 2 
      ENDDO 
!                                                                       
      WRITE (idout, 1000) j, np1 
      DO i = 1, np3 
      k = (j - 1) * np3 + i 
      WRITE (idout, 1100) k, cout (k), p (k), dp (k), pinc (k) 
      ENDDO 
      WRITE (idout, * ) ' ' 
      ENDDO 
!                                                                       
      WRITE (idout, 1200) 
      cout (npara - 1) = 'Background const.' 
      cout (npara) = 'Background slope' 
      DO i = npara - 1, npara 
      WRITE (idout, 1100) i, cout (i), p (i), dp (i), pinc (i) 
      ENDDO 
!                                                                       
 1000 FORMAT   (1x,'GSAS profile ',i2,': Type ',i2,/) 
 1100 FORMAT   (3x,'p(',i2,') : ',a20,' : ',                            &
     &                       g12.6,' +- ',g12.6,4x,'pinc : ',f2.0)      
 1200 FORMAT     (1x,'Global parameters:',/) 
 2000 FORMAT     ('Profile Coeff. ',i2) 
!                                                                       
      END SUBROUTINE show_gsas                      
!*****7*****************************************************************
      SUBROUTINE setup_gsas (ianz, werte, maxw, iname) 
!                                                                       
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw, maxmax 
      PARAMETER (maxmax = 50) 
!                                                                       
      CHARACTER ( * ) iname 
      REAL werte (maxw) 
      INTEGER ianz 
!                                                                       
      REAL pcoff (maxpara) 
      REAL wmax (maxmax) 
      REAL stheta, dspace, xpeak, inten, delt 
      REAL x1, x2, y1, y2 
      INTEGER ixm (maxmax) 
      INTEGER ncoff, itype, ibank, ipeaks, ima, i, j 
!                                                                       
      itype = nint (werte (1) ) 
      ibank = nint (werte (2) ) 
      IF (ianz.eq.4) then 
         ipeaks = nint (werte (4) ) 
      ELSE 
         ipeaks = 1 
      ENDIF 
!                                                                       
      ifen = fit_ifen 
      CALL do_fmax_xy (ikfit, wmax, ixm, maxmax, ima) 
      IF (ima.ge.1) then 
         xpeak = x (offxy (ikfit - 1) + ixm (1) ) 
      ELSE 
         xpeak = x (offxy (ikfit - 1) + len (ikfit) / 2) 
      ENDIF 
!                                                                       
      CALL read_prof (iname, ibank, itype, pcoff, ncoff, stheta, dspace,&
      xpeak)                                                            
      CALL cnvptp1 (itype, pcoff, ncoff, p, dspace, stheta, maxpara) 
!                                                                       
      inten = 0.0 
      DO i = 1, len (ikfit) - 1 
      delt = x (offxy (ikfit - 1) + i + 1) - x (offxy (ikfit - 1)       &
      + i)                                                              
      inten = inten + y (offxy (ikfit - 1) + i) * delt 
      ENDDO 
!                                                                       
      np1 = itype 
      np2 = ipeaks 
      np3 = ncoff 
      npara = ncoff * ipeaks + 2 
!                                                                       
      DO i = 1, ipeaks 
      IF (ima.ge.i) then 
         p ( (i - 1) * np3 + 1) = x (offxy (ikfit - 1) + ixm (i) ) 
      ELSE 
         p ( (i - 1) * np3 + 1) = x (offxy (ikfit - 1) + len (ikfit)    &
         / 2)                                                           
      ENDIF 
      p ( (i - 1) * np3 + 2) = inten 
      DO j = 3, ncoff 
      p ( (i - 1) * np3 + j) = p (j) 
      ENDDO 
      ENDDO 
!                                                                       
      y1 = y (offxy (ikfit - 1) + 1) 
      y2 = y (offxy (ikfit - 1) + len (ikfit) ) 
      x1 = x (offxy (ikfit - 1) + 1) 
      x2 = x (offxy (ikfit - 1) + len (ikfit) ) 
      p (npara - 1) = y1 
      p (npara) = (y2 - y1) / (x2 - x1) 
!                                                                       
      END SUBROUTINE setup_gsas                     
!*****7*****************************************************************
      SUBROUTINE read_prof (iname, ibank, itype, pcoff, ncoff, stheta,  &
      dspace, tof)                                                      
!                                                                       
      USE debug_mod 
      USE errlist_mod 
      USE prompt_mod 
      USE kuplot_config 
      USE trig_degree_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      CHARACTER ( * ) iname 
      REAL pcoff (maxpara) 
      REAL stheta, dspace, tof 
      INTEGER ibank, itype, ncoff 
!                                                                       
      CHARACTER(80) line 
      CHARACTER(6) search 
      CHARACTER(5) key 
      REAL tmp (4) 
      REAL difc, difa, zero, tth, l2 
      REAL secondterm 
      INTEGER k, ib, itmp, ll 
!                                                                       
      INTEGER len_str 
!     REAL sind 
!                                                                       
      CALL oeffne (12, iname, 'old') 
      IF (ier_num.ne.0) return 
!                                                                       
      WRITE (key, 1000) itype 
 1000 FORMAT    ('PRCF',i1) 
!                                                                       
!------ Read instrument parameter file information                      
!                                                                       
   20 CONTINUE 
      READ (12, '(a)', end = 40) line 
      ll = len_str (line) 
      IF (ll.eq.0) goto 20 
!                                                                       
      READ (line (1:ll) , '(4x,i2,a6)', err = 20) ib, search 
      IF (ib.ne.ibank) goto 20 
!                                                                       
      IF (search.eq.' ICONS') then 
         READ (line (13:ll), *, err = 998, end = 998) difc, difa, zero 
      ELSEIF (search.eq.'BNKPAR') then 
         READ (line (13:ll), *, err = 998) l2, tth 
      ELSEIF (search (1:5) .eq.key) then 
         IF (search (6:6) .eq.' ') then 
            READ (line (13:ll), *, err = 998) itmp, ncoff 
         ELSE 
            READ (search (6:6), * ) itmp 
            READ (line (13:ll), *, err = 998) tmp 
            DO k = 1, 4 
            pcoff ( (itmp - 1) * 4 + k) = tmp (k) 
            ENDDO 
         ENDIF 
      ENDIF 
      GOTO 20 
!                                                                       
   40 CONTINUE 
      CLOSE (12) 
!                                                                       
      stheta = sind (0.5 * tth) 
      IF (difa.ne.0.0) then 
         secondterm = sqrt (4.0 * difa * (tof - zero) + difc * difc) 
         dspace = ( - difc + secondterm) / 2. / difa 
      ELSE 
         dspace = (tof - zero) / difc 
      ENDIF 
!                                                                       
      RETURN 
!                                                                       
      ier_num = - 46 
      ier_typ = ER_APPL 
      CLOSE (12) 
      RETURN 
!                                                                       
  998 CONTINUE 
      ier_num = - 47 
      ier_typ = ER_APPL 
      CLOSE (12) 
!                                                                       
      END SUBROUTINE read_prof                      
!*****7*****************************************************************
      SUBROUTINE theory_gsas (xx, f, df, i) 
!                                                                       
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER(4) htype 
      REAL deriv (maxpara) 
      REAL xx, f, df (maxpara) 
      REAL pp (maxpara) 
      REAL tth, c, dtof, delt 
      REAL xnext 
      INTEGER ind, ip, i, j, k, ptype 
!                                                                       
      DO ind = 1, maxpara 
      df (ind) = 0.0 
      deriv (ind) = 0.0 
      ENDDO 
!                                                                       
      j = abs (i) 
      IF (j.ge.len (ikfit) ) return 
!                                                                       
      tth = 0. 
      ptype = np1 
      htype = 'PNT ' 
!                                                                       
      xnext = x (offxy (ikfit - 1) + j + 1) 
      delt = xnext - xx 
      f = p (npara - 1) + p (npara) * xx 
!                                                                       
      DO ip = 1, np2 
      DO k = 1, np3 
      pp (k) = p ( (ip - 1) * np3 + k) 
      ENDDO 
      dtof = (xx - pp (1) ) / 1000. 
      CALL prpcalc (htype, tth, ptype, pp, dtof, c, deriv, MAXPARA) 
      f = f + pp (2) * c 
!                                                                       
      DO k = 1, np3 
      p ( (ip - 1) * np3 + k) = pp (k) 
      ENDDO 
!                                                                       
      IF (i.ge.1) then 
!tep   write(*,'(5g15.6)') deriv(1),deriv(2),deriv(3),deriv(4),deriv(5) 
         IF (pinc ( (ip - 1) * np3 + 1) .ne.0.0) df ( (ip - 1) * np3 +  &
         1) = pp (2) * deriv (1) / 1000.                                
         IF (pinc ( (ip - 1) * np3 + 2) .ne.0.0) df ( (ip - 1) * np3 +  &
         2) = deriv (2)                                                 
         DO k = 3, npara - 2 
         IF (pinc ( (ip - 1) * np3 + k) .ne.0.0) df ( (ip - 1) * np3 +  &
         k) = deriv (k) * pp (2)                                        
         ENDDO 
      ENDIF 
      ENDDO 
!                                                                       
      IF (i.ge.1) then 
         IF (pinc (npara - 1) .ne.0.0) df (npara - 1) = 1.0 
         IF (pinc (npara) .ne.0.0) df (npara) = xx 
      ENDIF 
!                                                                       
      END SUBROUTINE theory_gsas                    
!*****7*****************************************************************
!       Lorenzian                                                       
!*****7*****************************************************************
      SUBROUTINE show_lor (idout) 
!                                                                       
      USE wink_mod
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER idout, i, iii 
      REAL zz, sint, ds0, ds2, ds3, dsint 
!                                                                       
      WRITE (idout, 1000) np1 
      WRITE (idout, 1100) 1, p (1), dp (1), pinc (1) 
      WRITE (idout, 1200) 2, p (2), dp (2), pinc (2) 
      DO i = 1, np1 
      WRITE (idout, 1300) i 
      iii = 2 + (i - 1) * 4 
      WRITE (idout, 1400) iii + 1, p (iii + 1), dp (iii + 1), pinc (iii &
      + 1)                                                              
      WRITE (idout, 1500) iii + 2, p (iii + 2), dp (iii + 2), pinc (iii &
      + 2)                                                              
      WRITE (idout, 1600) iii + 3, p (iii + 3), dp (iii + 3), pinc (iii &
      + 3)                                                              
      WRITE (idout, 1700) iii + 4, p (iii + 4), dp (iii + 4), pinc (iii &
      + 4)                                                              
!---------integral berechnen                                            
      zz = p (iii + 3) / p (iii + 4) + p (iii + 3) * p (iii + 4) 
      sint = p (iii + 1) * REAL(zpi) * zz 
      ds0 = REAL(zpi) * zz 
      ds2 = p (iii + 1) * REAL(zpi) * (1.0 / p (iii + 4) + p (iii + 4) ) 
      ds3 = p (iii + 1) * REAL(zpi) * ( - p (iii + 3) / p (iii + 4) **2 + p ( &
      iii + 3) )                                                        
      dsint = dp (iii + 1) * ds0 + dp (iii + 3) * ds2 + dp (iii + 4)    &
      * ds3                                                             
      WRITE (idout, 1800) sint, dsint 
      ENDDO 
      WRITE (idout, * ) ' ' 
!                                                                       
 1000 FORMAT     (1x,'Fitted',i3,' Lorenzian(s) : '/) 
 1100 FORMAT     (3x,'p(',i2,') : backgr. 1 : ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1200 FORMAT     (3x,'p(',i2,') : backgr. 2 : ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1300 FORMAT     (/,1x,'Lorenzian : ',i3,/) 
 1400 FORMAT     (3x,'p(',i2,') : peak      : ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1500 FORMAT     (3x,'p(',i2,') : position  : ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1600 FORMAT     (3x,'p(',i2,') : fwhm      : ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1700 FORMAT     (3x,'p(',i2,') : asymmetry : ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1800 FORMAT     (3x,'        integral  : ',g12.6,' +- ',g12.6) 
!                                                                       
      END SUBROUTINE show_lor                       
!*****7*****************************************************************
!
SUBROUTINE setup_lor (ianz, werte, maxw) 
!                                                                       
      USE errlist_mod 
      USE kuplot_config 
      USE kuplot_mod 
USE kuplot_fit6_set_theory
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxmax 
      PARAMETER (maxmax = 50) 
!                                                                       
      INTEGER maxw 
      REAL werte (maxw) 
      REAL wmax (maxmax) 
      INTEGER ixm (maxmax) 
      INTEGER ianz, ii, jj, ima, i 
!                                                                       
      IF (ianz.eq.0) then 
         npara = 6 
         np1 = 1 
      ELSEIF (ianz.eq.1) then 
         ii = nint (werte (1) ) 
         IF (ii.gt.0.and. (2 + 4 * ii) .le.maxpara) then 
            np1 = ii 
            npara = 2 + 4 * np1 
         ELSE 
            ier_num = - 31 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
         RETURN 
      ENDIF 
!                                                                       
      ii = offxy (ikfit - 1) + 1 
      jj = offxy (ikfit - 1) + len (ikfit) 
!                                                                       
      p (1) = y (ii) 
      pinc (1) = 1.0 
      p (2) = (y (jj) - y (ii) ) / (x (jj) - x (ii) ) 
      pinc (2) = 1.0 
!                                                                       
      ifen = fit_ifen 
      CALL do_fmax_xy (ikfit, wmax, ixm, maxmax, ima) 
      IF (ima.lt.np1) then 
         ier_num = - 30 
         ier_typ = ER_APPL 
         CALL errlist 
         DO i = 1, np1 
         wmax (i) = 1.0 
         ixm (i) = 1 
         ENDDO 
      ELSE 
         DO i = ima + 1, np1 
         wmax (i) = wmax (ima) 
         ixm (i) = ixm (ima) 
         ENDDO 
      ENDIF 
      CALL no_error 
!                                                                       
      IF (ier_num.eq.0) then 
         DO i = 1, np1 
         p (2 + (i - 1) * 4 + 1) = wmax (i) 
         pinc (2 + (i - 1) * 4 + 1) = 1.0 
         p (2 + (i - 1) * 4 + 2) = x (ii + ixm (i) - 1) 
         pinc (2 + (i - 1) * 4 + 2) = 1.0 
         p (2 + (i - 1) * 4 + 3) = 0.2 * abs (x (jj) - x (ii) ) 
         pinc (2 + (i - 1) * 4 + 3) = 1.0 
         p (2 + (i - 1) * 4 + 4) = 1.0 
         pinc (2 + (i - 1) * 4 + 4) = 0.0 
         ENDDO 
!                                                                       
         DO i = 1, npara 
         dp (i) = 0.0 
         ENDDO 
      ENDIF 
!
p_kuplot_theory => theory_lor
!                                                                       
END SUBROUTINE setup_lor                      
!***********************************************************************
!
SUBROUTINE theory_lor(MAXP, ix, iy, xx, yy, NPARA, params, par_names,          &
                      prange, l_do_deriv, data_dim, &
                       data_data, data_weight, data_x, data_y, &
data_calc, kupl_last,      &
                      ymod, dyda, LDERIV)
!
! Calculate a polynomial function
! ymod = SUM p[i]*xx^ind
!
IMPLICIT NONE
!
INTEGER                                              , INTENT(IN)  :: MAXP    ! Parameter array sizes
INTEGER                                              , INTENT(IN)  :: ix      ! Point number along x
INTEGER                                              , INTENT(IN)  :: iy      ! Point number along y
REAL                                                 , INTENT(IN)  :: xx      ! Point value  along x
REAL                                                 , INTENT(IN)  :: yy      ! Point value  along y
INTEGER                                              , INTENT(IN)  :: NPARA   ! Number of refined parameters
REAL            , DIMENSION(MAXP )                   , INTENT(IN)  :: params  ! Parameter values
CHARACTER(LEN=*), DIMENSION(MAXP)                    , INTENT(IN)  :: par_names    ! Parameter names
REAL            , DIMENSION(MAXP, 2                 ), INTENT(IN)  :: prange      ! Allowed parameter range
LOGICAL         , DIMENSION(MAXP )                   , INTENT(IN)  :: l_do_deriv  ! Parameter needs derivative
INTEGER         , DIMENSION(2)                       , INTENT(IN)  :: data_dim     ! Data array dimensions
REAL            , DIMENSION(data_dim(1), data_dim(2)), INTENT(IN)  :: data_data    ! Data array
REAL            , DIMENSION(data_dim(1), data_dim(2)), INTENT(IN)  :: data_weight  ! Data sigmas
REAL            , DIMENSION(data_dim(1))             , INTENT(IN)  :: data_x       ! Data coordinates x
REAL            , DIMENSION(data_dim(1))             , INTENT(IN)  :: data_y       ! Data coordinates y
REAL            , DIMENSION(data_dim(1), data_dim(2)), INTENT(OUT) :: data_calc    ! Data array
INTEGER                                              , INTENT(IN)  :: kupl_last    ! Last KUPLOT DATA that are needed
REAL                                                 , INTENT(OUT) :: ymod    ! Function value at (ix,iy)
REAL            , DIMENSION(NPARA)                   , INTENT(OUT) :: dyda    ! Function derivatives at (ix,iy)
LOGICAL                                              , INTENT(IN)  :: LDERIV  ! TRUE if derivative is needed
!                                                                       
REAL    :: o1, fwf, fw, xw, rn 
INTEGER :: ind, na, nu, np, nl, nlauf 
!                                                                       
DO ind = 1, npara 
   dyda(ind) = 0.0 
ENDDO 
!                                                                       
o1 = 2.0 * atan (1.0) 
nu = 2 
np = 4 
nl = (npara-nu)/np  ! np1 
!-------untergrund                                                      
ymod = params (1) + params (2) * xx 
DO nlauf = 1, nl 
   na = nu + (nlauf - 1) * np + 1 
!---------halbwertsbreiten                                              
   IF (xx.le.params (na + 1) ) then 
      fwf = 1.0 / params (na + 3) 
   ELSE 
      fwf = 1.0 * params (na + 3) 
   ENDIF 
   fw = params (na + 2) * fwf 
   xw = xx - params (na + 1) 
   rn = (fw * fw + 4.0 * xw * xw) 
!---------funktionswert berechnen                                       
   ymod = ymod + params (na) * fw * fw / rn 
!---------ableitungen berechnen                                         
   IF(LDERIV) then 
      IF(l_do_deriv(na) )    dyda(na) = fw * fw / rn 
      IF(l_do_deriv(na + 1)) dyda(na + 1) = 8.0 * params(na) * fw * fw * xw / (rn**2)
      IF(l_do_deriv(na + 2)) dyda(na + 2) = 2 * params(na) * fw / rn -  &
                                            2 * params (na) * fw * fw * fw / (rn**2)
      IF(l_do_deriv(na + 3)) THEN 
         IF (xx.le.params (na + 1) ) THEN 
               dyda(na + 3) = -dyda (na + 2) * (2.0 * (xx - params(na + 1))**2  &
                              * o1 / (fw**3)) * params(na + 2) / (params(na + 3)**2)     
         ELSE 
               dyda(na + 3) = dyda(na + 2) * (2.0 * (xx - params(na + 1))**2     &
                              * o1 / (fw**3) ) * params(na + 3)                          
         ENDIF 
      ENDIF 
   ENDIF 
ENDDO 
!-------untergrundsableitungen                                          
IF (LDERIV) THEN 
   IF(l_do_deriv(1)) dyda(1) = 1.0 
   IF(l_do_deriv(2)) dyda(2) = xx 
ENDIF 
!
data_calc(ix,iy) = ymod
!                                                                       
END SUBROUTINE theory_lor                     
!
!*****7*****************************************************************
!     Gaussian (1D)                                                     
!*****7*****************************************************************
      SUBROUTINE show_gauss (idout) 
!                                                                       
      USE wink_mod
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      REAL o, sqpio, zz, sint, ds0, ds2, ds3, dsint 
      INTEGER idout, i, iii 
!                                                                       
      o = sqrt (4.0 * alog (2.0) ) 
      sqpio = sqrt (REAL(pi)) / o * 0.5 
!                                                                       
      WRITE (idout, 1000) np1 
      WRITE (idout, 1100) 1, p (1), dp (1), pinc (1) 
      WRITE (idout, 1200) 2, p (2), dp (2), pinc (2) 
      DO i = 1, np1 
      WRITE (idout, 1300) i 
      iii = 2 + (i - 1) * 4 
      WRITE (idout, 1400) iii + 1, p (iii + 1), dp (iii + 1), pinc (iii &
      + 1)                                                              
      WRITE (idout, 1500) iii + 2, p (iii + 2), dp (iii + 2), pinc (iii &
      + 2)                                                              
      WRITE (idout, 1600) iii + 3, p (iii + 3), dp (iii + 3), pinc (iii &
      + 3)                                                              
      WRITE (idout, 1700) iii + 4, p (iii + 4), dp (iii + 4), pinc (iii &
      + 4)                                                              
!---------integral berechnen                                            
      zz = p (iii + 3) / p (iii + 4) + p (iii + 3) * p (iii + 4) 
      sint = p (iii + 1) * sqpio * zz 
      ds0 = sqpio * zz 
      ds2 = p (iii + 1) * sqpio * (1.0 / p (iii + 4) + p (iii + 4) ) 
      ds3 = p (iii + 1) * sqpio * ( - p (iii + 3) / p (iii + 4) **2 + p &
      (iii + 3) )                                                       
      dsint = dp (iii + 1) * ds0 + dp (iii + 3) * ds2 + dp (iii + 4)    &
      * ds3                                                             
      WRITE (idout, 1800) sint, dsint 
      ENDDO 
      WRITE (idout, * ) ' ' 
!                                                                       
 1000 FORMAT     (1x,'Fitted',i3,' Gaussian(s) : '/) 
 1100 FORMAT     (3x,'p(',i2,') : backgr. 1 : ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1200 FORMAT     (3x,'p(',i2,') : backgr. 2 : ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1300 FORMAT     (/,1x,'Gaussian : ',i3,/) 
 1400 FORMAT     (3x,'p(',i2,') : peak      : ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1500 FORMAT     (3x,'p(',i2,') : position  : ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1600 FORMAT     (3x,'p(',i2,') : fwhm      : ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1700 FORMAT     (3x,'p(',i2,') : asymmetry : ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1800 FORMAT     (3x,'        integral  : ',g12.6,' +- ',g12.6) 
!                                                                       
      END SUBROUTINE show_gauss                     
!***7*******************************************************************
SUBROUTINE setup_gauss (ianz, werte, maxw) 
!                                                                       
      USE errlist_mod 
      USE kuplot_config 
      USE kuplot_mod 
USE kuplot_fit6_set_theory
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxmax 
      PARAMETER (maxmax = 50) 
!                                                                       
      INTEGER maxw 
      REAL werte (maxw) 
      REAL wmax (maxmax) 
      INTEGER ixm (maxmax) 
      INTEGER ianz, ii, jj, ima, i 
!                                                                       
      IF (ianz.eq.0) then 
         npara = 6 
         np1 = 1 
      ELSEIF (ianz.eq.1) then 
         ii = nint (werte (1) ) 
         IF (ii.gt.0.and. (2 + 4 * ii) .le.maxpara) then 
            np1 = ii 
            npara = 2 + 4 * np1 
         ELSE 
            ier_num = - 31 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
         RETURN 
      ENDIF 
!                                                                       
      ii = offxy (ikfit - 1) + 1 
      jj = offxy (ikfit - 1) + len (ikfit) 
!                                                                       
      p (1) = y (ii) 
      pinc (1) = 1.0 
      p (2) = (y (jj) - y (ii) ) / (x (jj) - x (ii) ) 
      pinc (2) = 1.0 
!                                                                       
      ifen = fit_ifen 
      CALL do_fmax_xy (ikfit, wmax, ixm, maxmax, ima) 
      IF (ima.lt.np1) then 
         ier_num = - 30 
         ier_typ = ER_APPL 
         CALL errlist 
         DO i = 1, np1 
         wmax (i) = 1.0 
         ixm (i) = 1 
         ENDDO 
      ELSE 
         DO i = ima + 1, np1 
         wmax (i) = wmax (ima) 
         ixm (i) = ixm (ima) 
         ENDDO 
      ENDIF 
      CALL no_error 
!                                                                       
      IF (ier_num.eq.0) then 
         DO i = 1, np1 
         p (2 + (i - 1) * 4 + 1) = wmax (i) 
         pinc (2 + (i - 1) * 4 + 1) = 1.0 
         p (2 + (i - 1) * 4 + 2) = x (ii + ixm (i) - 1) 
         pinc (2 + (i - 1) * 4 + 2) = 1.0 
         p (2 + (i - 1) * 4 + 3) = 0.2 * abs (x (jj) - x (ii) ) 
         pinc (2 + (i - 1) * 4 + 3) = 1.0 
         p (2 + (i - 1) * 4 + 4) = 1.0 
         pinc (2 + (i - 1) * 4 + 4) = 0.0 
         ENDDO 
!                                                                       
         DO i = 1, npara 
         dp (i) = 0.0 
         ENDDO 
      ENDIF 
!                                                                       
p_kuplot_theory => theory_gauss
!
      END SUBROUTINE setup_gauss                    
!*********************************************************************  
SUBROUTINE theory_gauss(MAXP, ix, iy, xx, yy, NPARA, params, par_names,         &
                       prange, l_do_deriv, data_dim,  &
                       data_data, data_weight, data_x, data_y, &
                            data_calc, kupl_last,      &
                       ymod, dyda, LDERIV)
!
! Calculate a polynomial function
! ymod = SUM p[i]*xx^ind
!
IMPLICIT NONE
!
INTEGER                                              , INTENT(IN)  :: MAXP    ! Parameter array sizes
INTEGER                                              , INTENT(IN)  :: ix      ! Point number along x
INTEGER                                              , INTENT(IN)  :: iy      ! Point number along y
REAL                                                 , INTENT(IN)  :: xx      ! Point value  along x
REAL                                                 , INTENT(IN)  :: yy      ! Point value  along y
INTEGER                                              , INTENT(IN)  :: NPARA   ! Number of refined parameters
REAL            , DIMENSION(MAXP )                   , INTENT(IN)  :: params  ! Parameter values
CHARACTER(LEN=*), DIMENSION(MAXP)                    , INTENT(IN)  :: par_names    ! Parameter names
REAL            , DIMENSION(MAXP, 2                 ), INTENT(IN)  :: prange      ! Allowed parameter range
LOGICAL         , DIMENSION(MAXP )                   , INTENT(IN)  :: l_do_deriv  ! Parameter needs derivative
INTEGER         , DIMENSION(2)                       , INTENT(IN)  :: data_dim     ! Data array dimensions
REAL            , DIMENSION(data_dim(1), data_dim(2)), INTENT(IN)  :: data_data    ! Data array
REAL            , DIMENSION(data_dim(1), data_dim(2)), INTENT(IN)  :: data_weight  ! Data sigmas
REAL            , DIMENSION(data_dim(1))             , INTENT(IN)  :: data_x       ! Data coordinates x
REAL            , DIMENSION(data_dim(1))             , INTENT(IN)  :: data_y       ! Data coordinates y
REAL            , DIMENSION(data_dim(1), data_dim(2)), INTENT(OUT) :: data_calc    ! Data array
INTEGER                                              , INTENT(IN)  :: kupl_last    ! Last KUPLOT DATA that are needed
REAL                                                 , INTENT(OUT) :: ymod    ! Function value at (ix,iy)
REAL            , DIMENSION(NPARA)                   , INTENT(OUT) :: dyda    ! Function derivatives at (ix,iy)
LOGICAL                                              , INTENT(IN)  :: LDERIV  ! TRUE if derivative is needed
!     SUBROUTINE theory_gauss (xx, f, df, iwert) 
!                                                                       
!     USE kuplot_config 
!     USE kuplot_mod 
!                                                                       
!     IMPLICIT none 
!                                                                       
!     REAL xx, f, df (maxpara) 
      REAL o1, fw, fwf, exx 
      INTEGER ind, nunt, npar 
      INTEGER nlauf, na 
INTEGER :: np1   !  Number of Gauss peaks to fit
!                                                                       
o1   = 4.0 * alog (2.0) 
nunt = 2 
npar = 4 
np1  = (npara-nunt)/npar   ! Number of Gauss Peaks
!                                                                       
DO ind = 1, npara 
   dyda(ind) = 0.0 
write(*,*) ' PARAMS ', ind, params(ind)
ENDDO 
!-------Background
ymod = params (1) + params (2) * xx 
!
DO nlauf = 1, np1 
   na = nunt + (nlauf - 1) * npar + 1 
!---------Full Width at Half Maximum
   IF (xx.le.params(na + 1) ) then 
      fwf = 1.0 / params(na + 3) 
   ELSE 
      fwf = 1.0 * params(na + 3) 
   ENDIF 
   fw = params(na + 2) * fwf 
!---------Calculate function
   exx = exp ( - (xx - params (na + 1) ) **2 * o1 / (fw**2) ) 
   ymod = ymod + params (na) * exx 
!---------Calculate derivatives
   IF(LDERIV) THEN 
      IF(l_do_deriv (na)) dyda(na) = exx 
      IF(l_do_deriv (na + 1)) dyda (na + 1) = params (na) * exx * (2.0 *   &
         o1 * (xx - params (na + 1) ) / fw**2)                               
      IF (l_do_deriv (na + 2)) dyda (na + 2) = params (na) * exx * (2.0 *   &
         (xx - params (na + 1) ) **2 * o1 * fwf / (fw**3) )                  
      IF (l_do_deriv(na + 3)) then 
         IF (xx.le.params (na + 1) ) then 
               dyda (na + 3) = - params (na) * exx * (2.0 * (xx - params (na + 1) ) &
               **2 * o1 / (fw**3) ) * params (na + 2) / (params (na + 3) **2)     
         ELSE 
               dyda (na + 3) = params (na) * exx * (2.0 * (xx - params (na + 1) ) **&
               2 * o1 / (fw**3) ) * params (na + 3)                          
         ENDIF 
      ENDIF 
   ENDIF 
ENDDO 
!-------untergrundsableitungen                                          
IF(LDERIV) THEN 
   IF(l_do_deriv(1)) dyda (1) = 1.0 
   IF(l_do_deriv(2)) dyda (2) = xx 
ENDIF 
!
data_calc(ix, iy) = ymod
!                                                                       
END SUBROUTINE theory_gauss                   
!***7*******************************************************************
!       Pseudo-Voigt                                                    
!***7*******************************************************************
      SUBROUTINE show_psvgt (idout) 
!                                                                       
      USE wink_mod
      USE kuplot_config 
      USE kuplot_mod 
USE kuplot_fit_const
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER idout, i, iii 
!                                                                       
      WRITE (idout, 1000) np1 
      DO i = 1, nn_backgrd 
      WRITE (idout, 1100) i, i, p (i), dp (i), pinc (i) 
      ENDDO 
!     write (idout,1200) 2,p(2),dp(2),pinc(2)                           
!     write (idout,1210) 3,p(3),dp(3),pinc(3)                           
      DO i = 1, np1 
      WRITE (idout, 1300) i 
      iii = nn_backgrd+ (i - 1) * 6 
      WRITE (idout, 1400) iii + 1, p (iii + 1), dp (iii + 1), pinc (iii &
      + 1)                                                              
      WRITE (idout, 1500) iii + 2, p (iii + 2), dp (iii + 2), pinc (iii &
      + 2)                                                              
      WRITE (idout, 1600) iii + 3, p (iii + 3), dp (iii + 3), pinc (iii &
      + 3)                                                              
      WRITE (idout, 1700) iii + 4, p (iii + 4), dp (iii + 4), pinc (iii &
      + 4)                                                              
      WRITE (idout, 1800) iii + 5, p (iii + 5), dp (iii + 5), pinc (iii &
      + 5)                                                              
      WRITE (idout, 1900) iii + 6, p (iii + 6), dp (iii + 6), pinc (iii &
      + 6)                                                              
      ENDDO 
      WRITE (idout, * ) ' ' 
!                                                                       
 1000 FORMAT     (1x,'Fitted',i3,' Pseudo-Voigt(s) : '/) 
 1100 FORMAT     (3x,'p(',i2,') : backgr. ',i1,' : ',g12.6,' +- ',      &
     &                   g12.6,4x,'pinc : ',f2.0)                       
!1200 FORMAT     (3x,'p(',i2,') : backgr. 2 : ',g12.6,' +- ',g12.6,     &
!    &                   4x,'pinc : ',f2.0)                             
!1210 FORMAT     (3x,'p(',i2,') : backgr. 3 : ',g12.6,' +- ',g12.6,     &
!    &                   4x,'pinc : ',f2.0)                             
 1300 FORMAT     (/,1x,'Pseudo-Voigt : ',i3,/) 
 1400 FORMAT     (3x,'p(',i2,') : eta       : ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1500 FORMAT     (3x,'p(',i2,') : int. Inten: ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1600 FORMAT     (3x,'p(',i2,') : position  : ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1700 FORMAT     (3x,'p(',i2,') : fwhm      : ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1800 FORMAT     (3x,'p(',i2,') : asymmetry1: ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1900 FORMAT     (3x,'p(',i2,') : asymmetry2: ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
!                                                                       
      END SUBROUTINE show_psvgt                     
!***7*******************************************************************
SUBROUTINE setup_psvgt (ianz, werte, maxw) 
!                                                                       
      USE errlist_mod 
      USE kuplot_config 
      USE kuplot_mod 
USE kuplot_fit6_set_theory
USE kuplot_fit_const
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxmax 
      PARAMETER (maxmax = 50) 
!                                                                       
      INTEGER nu 
      INTEGER maxw 
      REAL werte (maxw) 
      REAL wmax (maxmax) 
      INTEGER ixm (maxmax) 
      INTEGER ianz, ii, jj, ima, i 
!                                                                       
nn_backgrd = 2
pp_origin  = 0.0
!
IF(ianz.eq.0) then 
   npara = 8 
   np1 = 1 
ELSEIF(ianz.le.3) then 
   IF (ianz.eq.3) then 
      nn_backgrd = nint (werte (3) ) 
   ENDIF 
   IF (ianz.ge.2) then 
      pp_origin = werte (2) 
   ENDIF 
   ii = nint (werte (1) ) 
   IF (ii.gt.0.and. (nn_backgrd + 6 * ii) .le.maxpara) then 
      np1 = ii 
      npara = nn_backgrd + 6 * np1 
   ELSE 
      ier_num = - 31 
      ier_typ = ER_APPL 
      RETURN 
   ENDIF 
!        pp_origin = 0.0 
!        nn_backgrd = 2 
ELSE 
   ier_num = - 6 
   ier_typ = ER_COMM 
   RETURN 
ENDIF 
nu = nn_backgrd 
!                                                                       
      ii = offxy (ikfit - 1) + 1 
      jj = offxy (ikfit - 1) + len (ikfit) 
!                                                                       
      p (1) = y (ii) 
      pinc (1) = 1.0 
      p (2) = (y (jj) - y (ii) ) / (x (jj) - x (ii) ) 
      pinc (2) = 1.0 
      DO i = 3, nn_backgrd 
      p (i) = 0.0 
      pinc (i) = 1.0 
      ENDDO 
!                                                                       
      ifen = fit_ifen 
      CALL do_fmax_xy (ikfit, wmax, ixm, maxmax, ima) 
      IF (ima.lt.np1) then 
         ier_num = - 30 
         ier_typ = ER_APPL 
         CALL errlist 
         DO i = 1, np1 
         wmax (i) = 1.0 
         ixm (i) = 1 
         ENDDO 
      ELSE 
         DO i = ima + 1, np1 
         wmax (i) = wmax (ima) 
         ixm (i) = ixm (ima) 
         ENDDO 
      ENDIF 
      CALL no_error 
!                                                                       
      IF (ier_num.eq.0) then 
         DO i = 1, np1 
         p (nu + (i - 1) * 6 + 1) = 0.5 
         pinc (nu + (i - 1) * 6 + 1) = 1.0 
         p (nu + (i - 1) * 6 + 2) = wmax (i) * 0.1 * abs (x (jj)        &
         - x (ii) ) * 3.14 / 2.                                         
         pinc (nu + (i - 1) * 6 + 2) = 1.0 
         p (nu + (i - 1) * 6 + 3) = x (ii + ixm (i) - 1) 
         pinc (nu + (i - 1) * 6 + 3) = 1.0 
         p (nu + (i - 1) * 6 + 4) = 0.1 * abs (x (jj) - x (ii) ) 
         pinc (nu + (i - 1) * 6 + 4) = 1.0 
         p (nu + (i - 1) * 6 + 5) = 0.0 
         pinc (nu + (i - 1) * 6 + 5) = 0.0 
         p (nu + (i - 1) * 6 + 6) = 0.0 
         pinc (nu + (i - 1) * 6 + 6) = 0.0 
         ENDDO 
!                                                                       
         DO i = 1, npara 
         dp (i) = 0.0 
         ENDDO 
      ENDIF 
write(*,*) ' PSEUDO ', npara, nn_backgrd,p_origin
write(*,*) ' PARAMS ', p(1:npara)
p_kuplot_theory => theory_psvgt
!                                                                       
      END SUBROUTINE setup_psvgt                    
!
!*********************************************************************  
!
SUBROUTINE theory_psvgt(MAXP, ix, iy, xx, yy, NPARA, params, par_names,          &
                       prange, l_do_deriv, data_dim, &
                       data_data, data_weight, data_x, data_y, &
 data_calc, kupl_last,      &
                       ymod, dyda, LDERIV)
!
! Calculate a Pseudo Voigt Function
! ymod = eta*L + (1-eta)*G
!
USE trig_degree_mod
USE wink_mod
USE kuplot_fit_const
!
IMPLICIT NONE
!
INTEGER                                              , INTENT(IN)  :: MAXP    ! Parameter array sizes
INTEGER                                              , INTENT(IN)  :: ix      ! Point number along x
INTEGER                                              , INTENT(IN)  :: iy      ! Point number along y
REAL                                                 , INTENT(IN)  :: xx      ! Point value  along x
REAL                                                 , INTENT(IN)  :: yy      ! Point value  along y
INTEGER                                              , INTENT(IN)  :: NPARA   ! Number of refined parameters
REAL            , DIMENSION(MAXP )                   , INTENT(IN)  :: params  ! Parameter values
CHARACTER(LEN=*), DIMENSION(MAXP)                    , INTENT(IN)  :: par_names    ! Parameter names
REAL            , DIMENSION(MAXP, 2                 ), INTENT(IN)  :: prange      ! Allowed parameter range
LOGICAL         , DIMENSION(MAXP )                   , INTENT(IN)  :: l_do_deriv  ! Parameter needs derivative
INTEGER         , DIMENSION(2)                       , INTENT(IN)  :: data_dim     ! Data array dimensions
REAL            , DIMENSION(data_dim(1), data_dim(2)), INTENT(IN)  :: data_data    ! Data array
REAL            , DIMENSION(data_dim(1), data_dim(2)), INTENT(IN)  :: data_weight  ! Data sigmas
REAL            , DIMENSION(data_dim(1))             , INTENT(IN)  :: data_x       ! Data coordinates x
REAL            , DIMENSION(data_dim(1))             , INTENT(IN)  :: data_y       ! Data coordinates y
REAL            , DIMENSION(data_dim(1), data_dim(2)), INTENT(OUT) :: data_calc    ! Data array
INTEGER                                              , INTENT(IN)  :: kupl_last    ! Last KUPLOT DATA that are needed
REAL                                                 , INTENT(OUT) :: ymod    ! Function value at (ix,iy)
REAL            , DIMENSION(NPARA)                   , INTENT(OUT) :: dyda    ! Function derivatives at (ix,iy)
LOGICAL                                              , INTENT(IN)  :: LDERIV  ! TRUE if derivative is needed
!     SUBROUTINE theory_psvgt (xx, f, df, i) 
!                                                                       
!     USE kuplot_config 
!     USE kuplot_mod 
!                                                                       
!     IMPLICIT none 
!                                                                       
!     REAL xx, f, df (maxpara) 
!     INTEGER i 
!                                                                       
      REAL fw, xw 
      REAL eta, pseudo 
      REAL asym, lore, lorn
      REAL gaus 
      REAL vln2, gpre 
      REAL zz, fa, fb 
      REAL dgdpos, dldpos 
      REAL dgdfw, dldfw 
      INTEGER j, ind, na, nu, np, nl, nlauf 
!
!REAL :: p_origin = 0.0   ! Needs work!!
!INTEGER :: n_backgrd  = 2
!                                                                       
!     REAL tand 
!                                                                       
DO ind = 1, npara 
   dyda(ind) = 0.0 
ENDDO 
!                                                                       
!     pi = 4.0 * atan (1.0) 
vln2 = 4. * alog (2.) 
gpre = 2. * sqrt (alog (2.) / pi) 
nu = nn_backgrd 
np = 6 
nl = (npara-nu)/np         ! Number of peaks
!-------untergrund                                                      
ymod = 0 
DO j = 1, nu 
   ymod = ymod + params (j) * (xx - pp_origin) ** (j - 1) 
ENDDO 
!                                                                       
DO nlauf = 1, nl 
   na = nu + (nlauf - 1) * np + 1 
   eta = params (na) 
   xw = xx - params (na + 2) 
   fw = params (na + 3) 
!---------asymmetry                                                     
   zz = xw / fw 
   fa = 2. * zz * exp ( - zz**2) 
   fb = 2. * (2 * zz**2 - 3.) * fa 
   asym = 1.0 
   asym = asym + (params (na + 4) * fa + params (na + 5) * fb) / tand (params (na + 2) )
!DBG        asym = asym+(p(na+4)*fa + p(na+5)*fb)/tand(p(na+2)*0.5)     
!     --Lorentzian                                                      
   lorn = (fw * fw + 4.0 * xw * xw) 
   lore = 2. / pi * fw / lorn 
   gaus = gpre / fw * exp ( - vln2 / fw**2 * xw**2) 
   pseudo = (eta * lore+ (1 - eta) * gaus) 
!---------calculate pseudo Voigt                                        
   ymod = ymod + params(na + 1) * asym * pseudo 
!---------calculate derivatives                                         
   IF (LDERIV) then 
!     ---- eta                                                          
      IF (l_do_deriv (na)) then 
         dyda(na) = params(na + 1) * asym * (lore-gaus) 
      ENDIF 
!     ---- intensity                                                    
      IF (l_do_deriv (na + 1)) then 
         dyda(na + 1) = asym * pseudo 
      ENDIF 
!     ---- position                                                     
      IF (l_do_deriv (na + 2)) then 
         dldpos = 2. / pi * 8. * fw * xw / lorn**2 
         dgdpos = gaus * (2. * vln2 / fw**2 * xw) 
         dyda(na + 2) = asym * params(na + 1) * (eta * dldpos + (1. - eta) * dgdpos)
      ENDIF 
!     ---- FWHM                                                         
      IF (l_do_deriv (na + 3)) then 
         dldfw = 2. / pi * ( - 1. * fw * fw + 4. * xw * xw) / lorn**2
         dgdfw = gaus * ( - 1. / fw + 2. * vln2 / fw**3 * xw * xw) 
         dyda(na + 3) = asym * params(na + 1) * (eta * dldfw + (1. - eta) * dgdfw)
      ENDIF 
!     ---- asymmetry parameter 1                                        
      IF (l_do_deriv (na + 4)) then 
         dyda (na + 4) = params (na + 1) * pseudo * fa / tand (params (na + 2) ) 
      ENDIF 
!     ---- asymmetry parameter 2                                        
      IF (l_do_deriv (na + 5)) then 
         dyda (na + 5) = params (na + 1) * pseudo * fb / tand (params (na + 2) ) 
      ENDIF 
   ENDIF 
ENDDO 
!-------derivative of background                                        
IF (LDERIV) then 
   DO j = 1, nu 
      IF (l_do_deriv (j)) dyda (j) = (xx - pp_origin) ** (j - 1) 
   ENDDO 
ENDIF 
!
data_calc(ix, iy) = ymod
!                                                                       
END SUBROUTINE theory_psvgt                   
!***7*******************************************************************
!     Gaussian (2D)                                                     
!***7*******************************************************************
      SUBROUTINE show_gauss_2d (idout) 
!                                                                       
      USE wink_mod
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      REAL o, sqpio, sqpio2 
      REAL zz_x, zz_y, sint, ds0, ds3, ds4, ds6, ds7, dsint 
      INTEGER idout, i, iii 
!                                                                       
      o = sqrt (4.0 * alog (2.0) ) 
      sqpio = sqrt (REAL(pi)) / o * 0.5 
      sqpio2 = sqpio * sqpio 
!                                                                       
      WRITE (idout, 1000) np1 
      WRITE (idout, 1100) 1, p (1), dp (1), pinc (1) 
      DO i = 1, np1 
      WRITE (idout, 1300) i 
      iii = 1 + (i - 1) * 8 
      WRITE (idout, 1400) iii + 1, p (iii + 1), dp (iii + 1), pinc (iii &
      + 1)                                                              
      WRITE (idout, 1410) iii + 2, p (iii + 2), dp (iii + 2), pinc (iii &
      + 2)                                                              
      WRITE (idout, 1420) iii + 3, p (iii + 3), dp (iii + 3), pinc (iii &
      + 3)                                                              
      WRITE (idout, 1430) iii + 4, p (iii + 4), dp (iii + 4), pinc (iii &
      + 4)                                                              
      WRITE (idout, 1440) iii + 5, p (iii + 5), dp (iii + 5), pinc (iii &
      + 5)                                                              
      WRITE (idout, 1450) iii + 6, p (iii + 6), dp (iii + 6), pinc (iii &
      + 6)                                                              
      WRITE (idout, 1460) iii + 7, p (iii + 7), dp (iii + 7), pinc (iii &
      + 7)                                                              
      WRITE (idout, 1470) iii + 8, p (iii + 8), dp (iii + 8), pinc (iii &
      + 8)                                                              
      zz_x = p (iii + 4) / p (iii + 7) + p (iii + 4) * p (iii + 7) 
      zz_y = p (iii + 5) / p (iii + 8) + p (iii + 5) * p (iii + 8) 
      sint = p (iii + 1) * sqpio2 * zz_x * zz_y 
      ds0 = sqpio2 * zz_x * zz_y 
      ds3 = p (iii + 1) * sqpio2 * zz_y * (1.0 / p (iii + 7) + p (iii + &
      7) )                                                              
      ds4 = p (iii + 1) * sqpio2 * zz_x * (1.0 / p (iii + 8) + p (iii + &
      8) )                                                              
      ds6 = p (iii + 1) * sqpio2 * zz_y * ( - p (iii + 4) / p (iii + 7) &
      **2 + p (iii + 4) )                                               
      ds7 = p (iii + 1) * sqpio2 * zz_x * ( - p (iii + 5) / p (iii + 8) &
      **2 + p (iii + 5) )                                               
      dsint = dp (iii + 1) * ds0 + dp (iii + 4) * ds3 + dp (iii + 5)    &
      * ds4 + dp (iii + 7) * ds6 + dp (iii + 8) * ds7                   
      WRITE (idout, 1800) sint, dsint 
      ENDDO 
      WRITE (idout, * ) ' ' 
!                                                                       
 1000 FORMAT     (1x,'Fitted',i3,' Gaussian(s) : '/) 
 1100 FORMAT     (3x,'p(',i2,') : backgr. 1 : ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1300 FORMAT     (/,1x,'Gaussian : ',i3,/) 
 1400 FORMAT     (3x,'p(',i2,') : peak      : ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1410 FORMAT     (3x,'p(',i2,') : position x: ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1420 FORMAT     (3x,'p(',i2,') : position y: ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1430 FORMAT     (3x,'p(',i2,') : fwhm a    : ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1440 FORMAT     (3x,'p(',i2,') : fwhm b    : ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1450 FORMAT     (3x,'p(',i2,') : angle a,x : ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1460 FORMAT     (3x,'p(',i2,') : asym. a   : ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1470 FORMAT     (3x,'p(',i2,') : asym. b   : ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1800 FORMAT     (3x,'        integral  : ',g12.6,' +- ',g12.6) 
!                                                                       
      END SUBROUTINE show_gauss_2d                  
!***7*******************************************************************
      SUBROUTINE setup_gauss_2d (ianz, werte, maxw) 
!                                                                       
      USE errlist_mod 
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxmax 
      PARAMETER (maxmax = 20) 
!                                                                       
      INTEGER maxw 
      REAL werte (maxw) 
      REAL wmax (maxmax) 
      INTEGER ixm (maxmax), iym (maxmax) 
      INTEGER ianz, ima, ii, jj, i 
!                                                                       
      IF (ianz.eq.0) then 
         npara = 1 + 8 
         np1 = 1 
      ELSEIF (ianz.eq.1) then 
         ii = nint (werte (1) ) 
         IF (ii.gt.0.and. (1 + 8 * ii) .le.maxpara) then 
            np1 = ii 
            npara = 1 + 8 * np1 
         ELSE 
            ier_num = - 31 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
         RETURN 
      ENDIF 
!                                                                       
      ii = offxy (ikfit - 1) 
      jj = offz (ikfit - 1) 
!                                                                       
!-------untergrund = wert li unten                                      
!                                                                       
      p (1) = z (jj + 1) 
      pinc (1) = 1.0 
!                                                                       
!-------nach lokalen maxima suchen                                      
!                                                                       
      ifen = fit_ifen 
      CALL do_fmax_z (ikfit, wmax, ixm, iym, maxmax, ima) 
      IF (ima.lt.np1) then 
         ier_num = - 30 
         ier_typ = ER_APPL 
         CALL errlist 
         DO i = 1, np1 
         wmax (i) = 1.0 
         ixm (i) = 1 
         iym (i) = 1 
         ENDDO 
      ELSE 
         DO i = ima + 1, np1 
         wmax (i) = wmax (ima) 
         ixm (i) = ixm (ima) 
         iym (i) = iym (ima) 
         ENDDO 
      ENDIF 
      CALL no_error 
!                                                                       
!-------parameter setzen                                                
!                                                                       
      IF (ier_num.eq.0) then 
         DO i = 1, np1 
         p (1 + (i - 1) * 8 + 1) = wmax (i) 
         pinc (1 + (i - 1) * 8 + 1) = 1.0 
         p (1 + (i - 1) * 8 + 2) = x (ii + ixm (i) ) 
         pinc (1 + (i - 1) * 8 + 2) = 1.0 
         p (1 + (i - 1) * 8 + 3) = y (ii + iym (i) ) 
         pinc (1 + (i - 1) * 8 + 3) = 1.0 
         p (1 + (i - 1) * 8 + 4) = 0.2 * abs (x (ii + nx (ikfit) )      &
         - x (ii + 1) )                                                 
         pinc (1 + (i - 1) * 8 + 4) = 1.0 
         p (1 + (i - 1) * 8 + 5) = 0.2 * abs (y (ii + ny (ikfit) )      &
         - y (ii + 1) )                                                 
         pinc (1 + (i - 1) * 8 + 5) = 1.0 
         p (1 + (i - 1) * 8 + 6) = 0.0 
         pinc (1 + (i - 1) * 8 + 6) = 0.0 
         p (1 + (i - 1) * 8 + 7) = 1.0 
         pinc (1 + (i - 1) * 8 + 7) = 0.0 
         p (1 + (i - 1) * 8 + 8) = 1.0 
         pinc (1 + (i - 1) * 8 + 8) = 0.0 
         ENDDO 
!                                                                       
         DO i = 1, npara 
         dp (i) = 0.0 
         ENDDO 
      ENDIF 
!                                                                       
      END SUBROUTINE setup_gauss_2d                 
!*********************************************************************  
      SUBROUTINE theory_gauss_2d (xx, f, df, i) 
!                                                                       
      USE wink_mod
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      REAL xx, f, o1, df (maxpara) 
      REAL rx, ry, rxs, rys, fwfx, fwx, fwfy, fwy, cosp, sinp 
      REAL exx, eyy, dfxs, dfys 
      INTEGER i, ind, nlauf, na, nu, np, ng, ix, iy 
!                                                                       
      DO ind = 1, maxpara 
      df (ind) = 0.0 
      ENDDO 
!                                                                       
      o1 = 4.0 * alog (2.0) 
      nu = 1 
      np = 8 
      ng = np1 
!-------x und y bestimmen                                               
      ix = (abs (i) - 1) / ny (ikfit) + 1 
      iy = abs (i) - (ix - 1) * ny (ikfit) 
      rx = x (offxy (ikfit - 1) + ix) 
      ry = y (offxy (ikfit - 1) + iy) 
!-------untergrund                                                      
      f = p (1) 
      DO nlauf = 1, ng 
      na = nu + (nlauf - 1) * np + 1 
!---------sinus und cosinus berechnen                                   
      cosp = cos (REAL(rad) * p (na + 5) ) 
      sinp = sin (REAL(rad) * p (na + 5) ) 
!---------transformation in hautachsensystem                            
      rxs = cosp * (rx - p (na + 1) ) + sinp * (ry - p (na + 2) ) 
      rys = - sinp * (rx - p (na + 1) ) + cosp * (ry - p (na + 2) ) 
!---------halbwertsbreiten                                              
      IF (rxs.le.0.0) then 
         fwfx = 1.0 / p (na + 6) 
      ELSE 
         fwfx = 1.0 * p (na + 6) 
      ENDIF 
      fwx = p (na + 3) * fwfx 
      IF (rys.le.0.0) then 
         fwfy = 1.0 / p (na + 7) 
      ELSE 
         fwfy = 1.0 * p (na + 7) 
      ENDIF 
      fwy = p (na + 4) * fwfy 
!---------funktionswert berechnen                                       
      exx = exp ( - rxs**2 * o1 / (fwx**2) ) 
      eyy = exp ( - rys**2 * o1 / (fwy**2) ) 
      f = f + p (na) * exx * eyy 
!---------ableitungen berechnen                                         
      IF (i.gt.0) then 
         IF (pinc (na) .ne.0) df (na) = exx * eyy 
         dfxs = p (na) * eyy * exx * ( - 2 * rxs * o1 / (fwx**2) ) 
         dfys = p (na) * eyy * exx * ( - 2 * rys * o1 / (fwy**2) ) 
         IF (pinc (na + 1) .ne.0) df (na + 1) = - cosp * dfxs + sinp *  &
         dfys                                                           
         IF (pinc (na + 2) .ne.0) df (na + 2) = - sinp * dfxs - cosp *  &
         dfys                                                           
         IF (pinc (na + 3) .ne.0) df (na + 3) = p (na) * eyy * exx *    &
         (2 * rxs**2 * o1 * fwfx / (fwx**3) )                           
         IF (pinc (na + 4) .ne.0) df (na + 4) = p (na) * eyy * exx *    &
         (2 * rys**2 * o1 * fwfy / (fwy**3) )                           
         IF (pinc (na + 5) .ne.0) df (na + 5) = ( - (rx - p (na + 1) )  &
         * sinp + (ry - p (na + 2) ) * cosp) * dfxs * REAL(rad)+ ( - (rx - p (&
         na + 1) ) * cosp - (ry - p (na + 2) ) * sinp) * dfys * REAL(rad)     
         IF (pinc (na + 6) .ne.0) then 
            IF (rxs.le.0.0) then 
               df (na + 6) = - p (na) * eyy * exx * (2 * rxs**2 * o1 /  &
               (fwx**3) ) * p (na + 3) / (p (na + 6) **2)               
            ELSE 
               df (na + 6) = p (na) * eyy * exx * (2 * rxs**2 * o1 /    &
               (fwx**3) ) * p (na + 3)                                  
            ENDIF 
         ENDIF 
         IF (pinc (na + 7) .ne.0) then 
            IF (rys.le.0.0) then 
               df (na + 7) = - p (na) * eyy * exx * (2 * rys**2 * o1 /  &
               (fwy**3) ) * p (na + 4) / (p (na + 7) **2)               
            ELSE 
               df (na + 7) = p (na) * eyy * exx * (2 * rys**2 * o1 /    &
               (fwy**3) ) * p (na + 4)                                  
            ENDIF 
         ENDIF 
      ENDIF 
      ENDDO 
!-------untergrundsableitungen                                          
      IF (i.gt.0) then 
         IF (pinc (1) .ne.0) df (1) = 1.0 
      ENDIF 
!                                                                       
      END SUBROUTINE theory_gauss_2d                
!***7*******************************************************************
!     Chebyshev polynom                                                 
!***7*******************************************************************
      SUBROUTINE show_poly_cheb (idout) 
!                                                                       
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER idout, i 
!                                                                       
      WRITE (idout, 1000) np1 
      DO i = 0, np1 
      WRITE (idout, 1100) i + 1, i, p (i + 1), dp (i + 1), pinc (i + 1) 
      ENDDO 
      WRITE (idout, * ) ' ' 
!                                                                       
 1000 FORMAT     (1x,'Fitted Chebyshev polynom of order ',i2,' : '/) 
 1100 FORMAT     (3x,'p(',i2,') : coeff. A(',i2,') : ',g12.6,           &
     &                   ' +- ',g12.6,4x,'pinc : ',f2.0)                
!                                                                       
      END SUBROUTINE show_poly_cheb                 
!***7*******************************************************************
      SUBROUTINE setup_poly_cheb (ianz, werte, maxw) 
!                                                                       
      USE errlist_mod 
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      REAL werte (maxw) 
      INTEGER ianz, ii, jj, i 
!                                                                       
      IF (ianz.eq.0) then 
         npara = 1 
         np1 = 1 
      ELSEIF (ianz.eq.1) then 
         ii = nint (werte (1) ) 
         IF (ii.ge.0.and. (1 + ii) .le.maxpara.and. (1 + ii) .le.5)     &
         then                                                           
            np1 = ii 
            npara = ii + 1 
         ELSE 
            ier_num = - 31 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
         RETURN 
      ENDIF 
!                                                                       
      ii = offxy (ikfit - 1) + 1 
      jj = offxy (ikfit - 1) + len (ikfit) 
!                                                                       
      p (1) = y (ii) 
      pinc (1) = 1.0 
      p (2) = (y (jj) - y (ii) ) / (x (jj) - x (ii) ) 
      pinc (2) = 1.0 
      DO i = 3, npara 
      p (i) = 0.0 
      pinc (i) = 1.0 
      ENDDO 
!                                                                       
      DO i = 1, npara 
      dp (i) = 0.0 
      ENDDO 
!                                                                       
      END SUBROUTINE setup_poly_cheb                
!***7*******************************************************************
      SUBROUTINE theory_poly_cheb (xx, f, df, iwert) 
!                                                                       
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      REAL xx, f, df (maxpara) 
      REAL xnew 
      INTEGER iwert, ind 
!                                                                       
      DO ind = 1, npara 
      df (ind) = 0.0 
      ENDDO 
!                                                                       
!------ Map x value on [-1,1] scale                                     
!                                                                       
      xnew = ( (xx - xmin (ikfit) ) - (xmax (ikfit) - xx) ) / (xmax (   &
      ikfit) - xmin (ikfit) )                                           
!                                                                       
      f = p (1) 
      IF (np1.ge.1) f = f + p (2) * xnew 
      IF (np1.ge.2) f = f + p (3) * (2 * xnew**2 - 1) 
      IF (np1.ge.3) f = f + p (4) * (4 * xnew**3 - 3 * xnew) 
      IF (np1.ge.4) f = f + p (5) * (8 * xnew**4 - 8 * xnew**2 + 1) 
!                                                                       
!-------Derivatives                                                     
!                                                                       
      IF (iwert.gt.0) then 
         IF (pinc (1) .ne.0.0) df (1) = 1.0 
         IF (pinc (2) .ne.0.0.and.np1.ge.1) df (2) = xnew 
         IF (pinc (3) .ne.0.0.and.np1.ge.2) df (3) = 2 * xnew**2 - 1 
         IF (pinc (4) .ne.0.0.and.np1.ge.3) df (4) = 4 * xnew**3 - 3 *  &
         xnew                                                           
         IF (pinc (5) .ne.0.0.and.np1.ge.4) df (5) = 8 * xnew**4 - 8 *  &
         xnew**2 + 1                                                    
      ENDIF 
!                                                                       
      END SUBROUTINE theory_poly_cheb               
!***7*******************************************************************
!     Polynom                                                           
!***7*******************************************************************
      SUBROUTINE show_poly (idout) 
!                                                                       
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER idout, i 
!                                                                       
      WRITE (idout, 1000) np1 
      DO i = 0, np1 
      WRITE (idout, 1100) i + 1, i, p (i + 1), dp (i + 1), pinc (i + 1) 
      ENDDO 
      WRITE (idout, * ) ' ' 
!                                                                       
 1000 FORMAT     (1x,'Fitted polynom of order ',i2,' : '/) 
 1100 FORMAT     (3x,'p(',i2,') : coeff. for x**',i2,' : ',g12.6,       &
     &                   ' +- ',g12.6,4x,'pinc : ',f2.0)                
!                                                                       
      END SUBROUTINE show_poly                      
!***7*******************************************************************
SUBROUTINE setup_poly(ianz, werte, maxw) 
!                                                                       
USE errlist_mod 
USE kuplot_config 
USE kuplot_mod 
USE kuplot_fit6_set_theory
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER,                  INTENT(IN) :: ianz 
INTEGER,                  INTENT(IN) :: maxw 
REAL   , DIMENSION(maxw), INTENT(IN) :: werte  !(maxw) 
!
INTEGER ::    ii, jj, i 
!
!                                                                       
IF (ianz == 0) THEN 
   npara = 1 
   np1 = 1 
ELSEIF (ianz == 1) THEN 
   ii = nint (werte (1) ) 
   IF (ii >= 0 .AND. (1 + ii) <= maxpara) THEN 
      np1 = ii 
      npara = ii + 1 
   ELSE 
      ier_num = - 31 
      ier_typ = ER_APPL 
      RETURN 
   ENDIF 
ELSE 
   ier_num = - 6 
   ier_typ = ER_COMM 
   RETURN 
ENDIF 
!                                                                       
ii = offxy (ikfit - 1) + 1 
jj = offxy (ikfit - 1) + len (ikfit) 
!                                                                       
p (1) = y (ii) 
pinc (1) = 1.0 
p (2) = (y (jj) - y (ii) ) / (x (jj) - x (ii) ) 
pinc (2) = 1.0 
DO i = 3, npara 
   p (i) = 0.0 
   pinc (i) = 1.0 
ENDDO 
!                                                                       
DO i = 1, npara 
   dp (i) = 0.0 
ENDDO 
!
p_kuplot_theory => theory_poly
!                                                                       
END SUBROUTINE setup_poly                     
!
!*******************************************************************************
!
SUBROUTINE theory_poly(MAXP, ix, iy, xx, yy, NPARA, params, par_names,          &
                       prange, l_do_deriv, data_dim, &
                       data_data, data_weight, data_x, data_y, &
data_calc, kupl_last,      &
                       ymod, dyda, LDERIV)
!
! Calculate a polynomial function
! ymod = SUM p[i]*xx^ind
!
IMPLICIT NONE
!
INTEGER                                              , INTENT(IN)  :: MAXP    ! Parameter array sizes
INTEGER                                              , INTENT(IN)  :: ix      ! Point number along x
INTEGER                                              , INTENT(IN)  :: iy      ! Point number along y
REAL                                                 , INTENT(IN)  :: xx      ! Point value  along x
REAL                                                 , INTENT(IN)  :: yy      ! Point value  along y
INTEGER                                              , INTENT(IN)  :: NPARA   ! Number of refined parameters
REAL            , DIMENSION(MAXP )                   , INTENT(IN)  :: params  ! Parameter values
CHARACTER(LEN=*), DIMENSION(MAXP)                    , INTENT(IN)  :: par_names    ! Parameter names
REAL            , DIMENSION(MAXP, 2                 ), INTENT(IN)  :: prange      ! Allowed parameter range
LOGICAL         , DIMENSION(MAXP )                   , INTENT(IN)  :: l_do_deriv  ! Parameter needs derivative
INTEGER         , DIMENSION(2)                       , INTENT(IN)  :: data_dim     ! Data array dimensions
REAL            , DIMENSION(data_dim(1), data_dim(2)), INTENT(IN)  :: data_data    ! Data array
REAL            , DIMENSION(data_dim(1), data_dim(2)), INTENT(IN)  :: data_weight  ! Data sigmas
REAL            , DIMENSION(data_dim(1))             , INTENT(IN)  :: data_x       ! Data coordinates x
REAL            , DIMENSION(data_dim(1))             , INTENT(IN)  :: data_y       ! Data coordinates y
REAL            , DIMENSION(data_dim(1), data_dim(2)), INTENT(OUT) :: data_calc    ! Data array
INTEGER                                              , INTENT(IN)  :: kupl_last    ! Last KUPLOT DATA that are needed
REAL                                                 , INTENT(OUT) :: ymod    ! Function value at (ix,iy)
REAL            , DIMENSION(NPARA)                   , INTENT(OUT) :: dyda    ! Function derivatives at (ix,iy)
LOGICAL                                              , INTENT(IN)  :: LDERIV  ! TRUE if derivative is needed
!
INTEGER :: ind
!
ymod = params(1)
!
IF(xx /= 0.0) THEN
   DO ind = 2, npara
      ymod = ymod + params(ind) * (xx**(ind-1))
   ENDDO
ENDIF
IF(LDERIV) THEN
   DO ind = 1, npara
      IF(l_do_deriv(ind)) THEN
         dyda(ind) = xx**(ind-1)
      ENDIF
   ENDDO
ENDIF
data_calc(ix, iy) = ymod
!
END SUBROUTINE theory_poly
!!***7*******************************************************************
!      SUBROUTINE theory_poly (xx, f, df, iwert) 
!!                                                                       
!      USE kuplot_config 
!      USE kuplot_mod 
!!                                                                       
!      IMPLICIT none 
!!                                                                       
!      REAL xx, f, df (maxpara) 
!      INTEGER iwert, ind 
!!                                                                       
!      DO ind = 1, npara 
!      df (ind) = 0.0 
!      ENDDO 
!!                                                                       
!      f = p (1) 
!      DO ind = 1, np1 
!      IF (xx.ne.0) f = f + p (ind+1) * (xx**ind) 
!      ENDDO 
!!                                                                       
!!-------Derivatives                                                     
!!                                                                       
!      IF (iwert.gt.0) then 
!         DO ind = 0, np1 
!         IF (pinc (ind+1) .ne.0) then 
!            df (ind+1) = (xx**ind) 
!         ENDIF 
!         ENDDO 
!      ENDIF 
!                                                                       
!      END SUBROUTINE theory_poly                    
!***7*******************************************************************
!     Scale factor + Background Polynom                                 
!***7*******************************************************************
      SUBROUTINE show_backpoly (idout) 
!                                                                       
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER idout, i 
!                                                                       
      WRITE (idout, 1000) np1 - 2 
      WRITE (idout, 1020) ikfit2 
      WRITE (idout, 1050) p (1), dp (1), pinc (1) 
      DO i = 2, np1 
      WRITE (idout, 1100) i, i - 2, p (i), dp (i), pinc (i) 
      ENDDO 
      WRITE (idout, * ) ' ' 
!                                                                       
 1000 FORMAT     (1x,'Fitted Background polynom of order ',i2,' : '/) 
 1020 FORMAT     (1x,'Constant part from data set        ',i2) 
 1050 FORMAT     (1x,'Scale factor               : ',g12.6,             &
     &                   ' +- ',g12.6,4x,'pinc : ',f2.0)                
 1100 FORMAT     (3x,'p(',i2,') : coeff. for x**',i2,' : ',g12.6,       &
     &                   ' +- ',g12.6,4x,'pinc : ',f2.0)                
!                                                                       
      END SUBROUTINE show_backpoly                  
!***7*******************************************************************
SUBROUTINE setup_backpoly (ianz, werte, maxw) 
!                                                                       
      USE errlist_mod 
      USE kuplot_config 
      USE kuplot_mod 
USE kuplot_fit6_set_theory
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      REAL werte (maxw) 
      INTEGER ianz, ii, jj, i 
!                                                                       
      IF (ianz.eq.1) then 
         ikfit2 = nint (werte (1) ) 
         npara = 1 
         np1 = 1 
      ELSEIF (ianz.eq.2) then 
         ikfit2 = nint (werte (1) ) 
         ii = nint (werte (2) ) 
         IF (ii.ge.0.and. (1 + ii) .le.maxpara) then 
            np1 = ii 
            npara = ii + 1 
         ELSE 
            ier_num = - 31 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
         RETURN 
      ENDIF 
!                                                                       
      ii = offxy (ikfit - 1) + 1 
      jj = offxy (ikfit2 - 1) + 1 
!                                                                       
      IF (np1.ge.2) then 
         p (2) = y (ii) - y (jj) 
         pinc (2) = 1.0 
      ELSE 
         pinc (2) = 0.0 
      ENDIF 
!                                                                       
      ii = offxy (ikfit - 1) + len (ikfit) / 2 
      jj = offxy (ikfit2 - 1) + len (ikfit2) / 2 
      IF (y (jj) .ne.0) then 
         p (1) = (y (ii) - p (2) ) / y (jj) 
      ENDIF 
      pinc (1) = 1.0 
      DO i = 3, npara 
      p (i) = 0.0 
      pinc (i) = 1.0 
      ENDDO 
!                                                                       
      DO i = 1, npara 
      dp (i) = 0.0 
      ENDDO 
!
p_kuplot_theory => theory_backpoly
!                                                                       
      END SUBROUTINE setup_backpoly                 
!***7*******************************************************************
SUBROUTINE theory_backpoly(MAXP, ix, iy, xx, yy, NPARA, params, par_names,      &
                       prange, l_do_deriv, data_dim,                            &
                       data_data, data_weight, data_x, data_y,                  &
                       data_calc, kupl_last,      &
                       ymod, dyda, LDERIV)
!
! Calculate a background polynomial function
! ymod = SUM p[i]*xx^ind
!
IMPLICIT NONE
!
INTEGER                                              , INTENT(IN)  :: MAXP    ! Parameter array sizes
INTEGER                                              , INTENT(IN)  :: ix      ! Point number along x
INTEGER                                              , INTENT(IN)  :: iy      ! Point number along y
REAL                                                 , INTENT(IN)  :: xx      ! Point value  along x
REAL                                                 , INTENT(IN)  :: yy      ! Point value  along y
INTEGER                                              , INTENT(IN)  :: NPARA   ! Number of refined parameters
REAL            , DIMENSION(MAXP )                   , INTENT(IN)  :: params  ! Parameter values
CHARACTER(LEN=*), DIMENSION(MAXP)                    , INTENT(IN)  :: par_names    ! Parameter names
REAL            , DIMENSION(MAXP, 2                 ), INTENT(IN)  :: prange      ! Allowed parameter range
LOGICAL         , DIMENSION(MAXP )                   , INTENT(IN)  :: l_do_deriv  ! Parameter needs derivative
INTEGER         , DIMENSION(2)                       , INTENT(IN)  :: data_dim     ! Data array dimensions
REAL            , DIMENSION(data_dim(1), data_dim(2)), INTENT(IN)  :: data_data    ! Data array
REAL            , DIMENSION(data_dim(1), data_dim(2)), INTENT(IN)  :: data_weight  ! Data sigmas
REAL            , DIMENSION(data_dim(1))             , INTENT(IN)  :: data_x       ! Data coordinates x
REAL            , DIMENSION(data_dim(1))             , INTENT(IN)  :: data_y       ! Data coordinates y
REAL            , DIMENSION(data_dim(1), data_dim(2)), INTENT(OUT) :: data_calc    ! Data array
INTEGER                                              , INTENT(IN)  :: kupl_last    ! Last KUPLOT DATA that are needed
REAL                                                 , INTENT(OUT) :: ymod    ! Function value at (ix,iy)
REAL            , DIMENSION(NPARA)                   , INTENT(OUT) :: dyda    ! Function derivatives at (ix,iy)
LOGICAL                                              , INTENT(IN)  :: LDERIV  ! TRUE if derivative is needed
!     SUBROUTINE theory_backpoly (xx, f, df, iwert) 
!                                                                       
!     USE kuplot_config 
!     USE kuplot_mod 
!                                                                       
!     IMPLICIT none 
!                                                                       
!     REAL xx, f, df (maxpara) 
      INTEGER ind 
      INTEGER iii 
!                                                                       
!     INTEGER idata, ipoint 
      REAL arg 
!                                                                       
      DATA iii / 0 / 
!                                                                       
DO ind = 1, npara 
   dyda(ind) = 0.0 
ENDDO 
!                                                                       
!     idata = ikfit2 
!     ipoint = iabs (iwert) 
!     arg = (xx - xmin (idata) ) 
arg = (xx-data_x(1))
!                                                                       
ymod = params(1) * data_data(ix,iy)
IF (npara >= 2) THEN 
   ymod = ymod + params(2) 
ENDIF 
DO ind = 3, npara 
   ymod = ymod + params(ind) * (arg**(ind-2))
ENDDO 
!                                                                       
!-------Derivatives                                                     
!                                                                       
IF(LDERIV) THEN 
   ind = 1 
   IF (l_do_deriv(ind)) then 
      dyda(ind) = data_data(ix,iy)
   ENDIF 
   DO ind = 2, npara 
      IF (l_do_deriv(ind)) then 
         dyda (ind) = (arg** (ind-2) ) 
      ENDIF 
   ENDDO 
ENDIF 
!
data_calc(ix, iy) = ymod
!                                                                       
END SUBROUTINE theory_backpoly                
!
!*******************************************************************************
!
!
!*******************************************************************************
!
SUBROUTINE kuplot_mrq(MAXP, NPARA, ncycle, kupl_last, par_names, par_ind,  &
                      MAXF, nfixed, fixed, fixed_ind, pf, &
                      data_dim, &
                      data_data, data_weight, data_x, data_y,              &
                      data_calc,                                           &
                      conv_dp_sig, conv_dchi2, conv_chi2, conv_conf,       &
                      lconvergence,                                        &
                      chisq, conf, lamda_fin, rval, rexp, p, prange, dp, cl)
!+                                                                      
!   This routine runs the refinement cycles, interfaces with the 
!   Levenberg-Marquardt routine modified after Numerical Recipes
!-                                                                      
!
!USE refine_set_param_mod
!
USE errlist_mod
USE ber_params_mod
USE calc_expr_mod
USE gamma_mod
USE param_mod 
USE prompt_mod 
!                                                                       
IMPLICIT NONE
!
INTEGER                                              , INTENT(IN)  :: MAXP        ! Parameter array size
INTEGER                                              , INTENT(IN)  :: NPARA       ! Number of parameters
INTEGER                                              , INTENT(IN)  :: ncycle      ! maximum cycle number
INTEGER                                              , INTENT(IN)  :: kupl_last   ! Last KUPLOT DATA that are needed
CHARACTER(LEN=*), DIMENSION(MAXP)                    , INTENT(IN)  :: par_names   ! Parameter names
INTEGER         , DIMENSION(MAXP)                    , INTENT(IN)  :: par_ind     ! Index of param in KUPLOT list
INTEGER                                              , INTENT(IN)  :: MAXF        ! Fixed Parameter array size
INTEGER                                              , INTENT(IN)  :: nfixed      ! Number of fixed parameters
CHARACTER(LEN=*), DIMENSION(MAXF)                    , INTENT(IN)  :: fixed       ! Fixed Parameter names
INTEGER         , DIMENSION(MAXF)                    , INTENT(IN)  :: fixed_ind   ! Index of fixed param in KUPLOT list
REAL            , DIMENSION(MAXF)                    , INTENT(IN)  :: pf          ! Fixed Parameter array
INTEGER         , DIMENSION(2)                       , INTENT(IN)  :: data_dim    ! Data array dimensions
REAL            , DIMENSION(data_dim(1), data_dim(2)), INTENT(IN)  :: data_data   ! Data array
REAL            , DIMENSION(data_dim(1), data_dim(2)), INTENT(IN)  :: data_weight ! Data sigmas
REAL            , DIMENSION(data_dim(1))             , INTENT(IN)  :: data_x      ! Data coordinates x
REAL            , DIMENSION(data_dim(1))             , INTENT(IN)  :: data_y      ! Data coordinates y
REAL            , DIMENSION(data_dim(1), data_dim(2)), INTENT(OUT) :: data_calc   ! Calculated
REAL                                                 , INTENT(IN)  :: conv_dp_sig ! Max parameter shift
REAL                                                 , INTENT(IN)  :: conv_dchi2  ! Max Chi^2     shift
REAL                                                 , INTENT(IN)  :: conv_chi2   ! Min Chi^2     value 
REAL                                                 , INTENT(IN)  :: conv_conf   ! Min confidence level
LOGICAL                                              , INTENT(INOUT) :: lconvergence ! Convergence criteria were reached
REAL                                                 , INTENT(OUT) :: chisq       ! Chi^2 value
REAL                                                 , INTENT(OUT) :: conf        ! Confidence test from Gamma function
REAL                                                 , INTENT(OUT) :: lamda_fin   ! Final Marquardt lambda
REAL                                                 , INTENT(OUT) :: rval        ! Weighted R-value
REAL                                                 , INTENT(OUT) :: rexp        ! Expected R-value
!
REAL            , DIMENSION(MAXP)                    , INTENT(INOUT) :: p         ! Parameter array
REAL            , DIMENSION(MAXP,2)                  , INTENT(INOUT) :: prange    ! Parameter range
REAL            , DIMENSION(MAXP)                    , INTENT(INOUT) :: dp        ! Parameter sigmas
REAL            , DIMENSION(NPARA, NPARA)            , INTENT(INOUT) :: cl        ! Covariance matrix
!
INTEGER :: icyc
INTEGER :: k
INTEGER :: last_i, prev_i
REAL   , DIMENSION(NPARA, NPARA) :: alpha
REAL   , DIMENSION(NPARA       ) :: beta
REAL    :: alamda
!
INTEGER :: MAXW
INTEGER :: ianz
!
CHARACTER(LEN=1024), DIMENSION(:), ALLOCATABLE :: cpara
INTEGER            , DIMENSION(:), ALLOCATABLE :: lpara
REAL               , DIMENSION(:), ALLOCATABLE :: werte
!
REAL               , DIMENSION(0:3) :: last_chi
REAL               , DIMENSION(0:3) :: last_shift
REAL               , DIMENSION(0:3) :: last_conf 
REAL               , DIMENSION(:,:), ALLOCATABLE :: last_p 
!REAL :: shift
!
alamda     = -0.01    ! Negative lamda initializes MRQ routine
rval       = 0.0
rexp       = 0.0
cl(:,:)    = 0.0
alpha(:,:) = 0.0
beta(:)    = 0.0
chisq      = 0.0
last_chi(:)   = HUGE(0.0)
last_shift(:) = -1.0
last_conf(:)  = HUGE(0.0)
last_i        = 0
!
! Set initial parameter values
!
MAXW = NPARA
ALLOCATE(cpara(MAXW))
ALLOCATE(lpara(MAXW))
ALLOCATE(werte(MAXW))
DO k = 1, NPARA                  ! Copy parameter names to prepare calculation
   cpara(k) = par_names(k)
   lpara(k) = LEN_TRIM(par_names(k))
ENDDO
ianz = NPARA
CALL ber_params(ianz, cpara, lpara, werte, MAXW)
DO k = 1, NPARA                  ! Copy values of the parameters
   p(k) = werte(k)
ENDDO
!
ALLOCATE(last_p(0:2,NPARA))
last_p(2,:) = p(1:NPARA)
prev_i = 2
!
WRITE(output_io,'(a,10x,a,7x,a,6x,a)') 'Cyc Chi^2/(N-P)   MAX(dP/sig)   Conf',&
                                       'Lambda', 'wRvalue', 'Rexp'
!
lconvergence = .FALSE.
icyc = 0
cycles:DO
   CALL kuplot_mrqmin(MAXP, data_dim, data_data, data_weight, data_x, data_y,  &
               data_calc, p, NPARA, &
               par_names, par_ind, &
               MAXF, nfixed, fixed, fixed_ind, pf, &
               prange, kupl_last, cl, alpha, beta, chisq, alamda)
   IF(ier_num/=0) EXIT cycles
   CALL kuplot_rvalue(data_dim, data_data, data_weight, data_calc, rval, rexp, NPARA)
   last_shift(last_i) = 0.0
   DO k=1, NPARA
!     CALL kuplot_set_param(NPARA, par_names(k), k, p(k))   ! Updata parameters
      IF(ier_num/=0) EXIT cycles
      last_p(last_i,:) = p(1:NPARA)
      IF(cl(k,k)/= 0.0) THEN
         last_shift(last_i) = MAX(last_shift(last_i), ABS((p(k)-last_p(prev_i,k))/SQRT(cl(k,k))))
      ENDIF
   ENDDO
   conf = gammaq(REAL(data_dim(1)*data_dim(2)-2)*0.5, chisq*0.5)
   last_chi( last_i) = chisq/(data_dim(1)*data_dim(2)-NPARA)        ! Store Chi^2/(NDATA-NPARA)
   last_conf(last_i) = conf
   WRITE(*,'(i3,6g13.5e2)') icyc,chisq/(data_dim(1)*data_dim(2)-NPARA), last_shift(last_i), conf, alamda, rval, rexp
!
   CALL kuplot_rvalue(data_dim, data_data, data_weight, data_calc, rval, rexp, NPARA)
!  CALL kuplot_best(rval)                                              ! Write best macro
   IF(ABS(last_chi( last_i)-last_chi(prev_i))<conv_dchi2 .AND.   &
      last_shift(last_i) < conv_dp_sig                   .AND.   &
      last_conf(last_i)  > conv_conf                     .OR.    &
      last_chi( last_i)  < conv_chi2                   ) THEN
      WRITE(output_io, '(/,a)') 'Convergence reached '        
      lconvergence = .TRUE.
      EXIT cycles
   ENDIF
   IF(alamda > 0.5*HUGE(0.0)) THEN
      lconvergence = .FALSE.
      EXIT cycles
   ENDIF
   prev_i = last_i
   last_i = MOD(last_i+1, 3)
   icyc   = icyc + 1
   IF(icyc==ncycle) EXIT cycles
ENDDO cycles
!
IF(.NOT. lconvergence) THEN
   IF(icyc==ncycle) THEN
      WRITE(output_io,'(/,a)') ' Maximum cycle is reached without convergence'
   ELSEIF(alamda > 0.5*HUGE(0.0)) THEN
      WRITE(output_io,'(/,a)') ' MRQ Parameter reached infinity'
   ENDIF
ENDIF
!
! Make a final call with lamda = 0 to ensure that model is up to date
!
IF(ier_num==0) THEN
   lamda_fin = alamda
   alamda = 0.0
!
   CALL kuplot_mrqmin(MAXP, data_dim, data_data, data_weight, data_x, data_y,  &
               data_calc, p, NPARA, &
               par_names, par_ind,  &
               MAXF, nfixed, fixed, fixed_ind, pf, &
               prange, kupl_last, cl, alpha, beta, chisq, alamda)
   CALL kuplot_rvalue(data_dim, data_data, data_weight, data_calc, rval, rexp, NPARA)
!  CALL kuplot_best(rval)                                              ! Write best macro
   DO k = 1, NPARA
      dp(k) = SQRT(ABS(cl(k,k)))
   ENDDO
ENDIF
!
DEALLOCATE(cpara)
DEALLOCATE(lpara)
DEALLOCATE(werte)
DEALLOCATE(last_p)
!
END SUBROUTINE kuplot_mrq
!
!*******************************************************************************
!
SUBROUTINE kuplot_mrqmin(MAXP, data_dim, data_data, data_weight, data_x, data_y, &
               data_calc, a, NPARA, &
    par_names, par_ind, &
                      MAXF, nfixed, fixed, fixed_ind, pf, &
    prange, kupl_last, covar, alpha, beta, chisq, alamda)
!
!  Least squares routine adapted from Numerical Recipes Chapter 14
!  p 526
!
USE errlist_mod
USE gaussj_mod
!
IMPLICIT NONE
!
INTEGER                                             , INTENT(IN)    :: MAXP        ! Array size for parameters
INTEGER         , DIMENSION(2)                      , INTENT(IN)    :: data_dim    ! no of ata points along x, y
REAL            , DIMENSION(data_dim(1),data_dim(2)), INTENT(IN)    :: data_data   ! Data values
REAL            , DIMENSION(data_dim(1),data_dim(2)), INTENT(IN)    :: data_weight ! Data sigmas
REAL            , DIMENSION(data_dim(1)           ) , INTENT(IN)    :: data_x      ! Actual x-coordinates
REAL            , DIMENSION(data_dim(2)           ) , INTENT(IN)    :: data_y      ! Actual y-coordinates
REAL            , DIMENSION(data_dim(1), data_dim(2)), INTENT(OUT)  :: data_calc   ! Calculated
INTEGER                                             , INTENT(IN)    :: NPARA       ! number of parameters
REAL            , DIMENSION(MAXP, 2              )  , INTENT(IN)    :: prange      ! Allowed parameter range
INTEGER                                             , INTENT(IN)    :: kupl_last   ! Last KUPLOT DATA that are needed
REAL            , DIMENSION(MAXP)                   , INTENT(INOUT) :: a           ! Parameter values
CHARACTER(LEN=*), DIMENSION(MAXP)                   , INTENT(IN)    :: par_names   ! Parameter names
INTEGER         , DIMENSION(MAXP)                   , INTENT(IN)    :: par_ind     ! Index of param in KUPLOT list
INTEGER                                             , INTENT(IN)    :: MAXF        ! Fixed Parameter array size
INTEGER                                             , INTENT(IN)    :: nfixed      ! Number of fixed parameters
CHARACTER(LEN=*), DIMENSION(MAXF)                   , INTENT(IN)    :: fixed       ! Fixed Parameter names
INTEGER         , DIMENSION(MAXF)                   , INTENT(IN)    :: fixed_ind   ! Index of fixed param in KUPLOT list
REAL            , DIMENSION(MAXF)                   , INTENT(IN)    :: pf          ! Fixed Parameter array
REAL            , DIMENSION(NPARA, NPARA)           , INTENT(OUT)   :: covar       ! Covariance matrix
REAL            , DIMENSION(NPARA, NPARA)           , INTENT(INOUT) :: alpha       ! Temp arrays
REAL            , DIMENSION(NPARA     )             , INTENT(INOUT) :: beta        ! Temp arrays
REAL                                                , INTENT(INOUT) :: chisq       ! Chi Squared
REAL                                                , INTENT(INOUT) :: alamda      ! Levenberg parameter
!
INTEGER                         :: j     = 0
INTEGER                         :: k     = 0
REAL   , DIMENSION(MAXP)        :: atry != 0.0
REAL   , DIMENSION(NPARA)       :: da != 0.0
REAL   , DIMENSION(NPARA,NPARA) :: cl != 0.0
REAL                            :: ochisq = 0.0
LOGICAL, PARAMETER              :: LDERIV = .TRUE.
LOGICAL, PARAMETER              :: NDERIV = .FALSE.
!
ochisq = chisq
IF(alamda < 0) THEN                   ! Initialization
   alamda = 0.001
   CALL kuplot_mrqcof(MAXP, data_dim, data_data, data_weight, data_x, data_y,   &
               data_calc, a, &
               NPARA, par_names, par_ind, prange, &
                      MAXF, nfixed, fixed, fixed_ind, pf, &
               kupl_last, alpha, beta, &
               chisq, LDERIV) !, funcs)
   IF(ier_num/=0) THEN
      RETURN
   ENDIF
   ochisq = chisq
   DO j=1, NPARA
      IF(prange(j,1)<=prange(j,2)) THEN
         atry(j) = MIN(prange(j,2),MAX(prange(j,1),a(j)))
      ELSE
         atry(j) = a(j)
      ENDIF
   ENDDO
ENDIF
!
DO j=1, NPARA
   DO k=1, NPARA
      covar(j,k) = alpha(j,k)
   ENDDO
   covar(j,j) = alpha(j,j)*(1.0+alamda)       ! Adjust by Levenberg-Marquard algorithm
   da(j) = beta(j)
ENDDO
!
CALL gausj(NPARA, covar, da)
IF(ier_num/=0) THEN
   RETURN
ENDIF
!
IF(alamda==0) THEN
!
!  Last call to update the structure, no derivative is needed
!
   cl(:,:) = covar(:,:)       ! Back up of covariance
   CALL kuplot_mrqcof(MAXP, data_dim, data_data, data_weight, data_x, data_y,   &
               data_calc, a, &
               NPARA, par_names, par_ind, prange,  &
               MAXF, nfixed, fixed, fixed_ind, pf, &
               kupl_last, covar, da, chisq, NDERIV) !, funcs)
   IF(ier_num/=0) THEN
      RETURN
   ENDIF
   covar(:,:) = cl(:,:)       ! Restore    covariance
   RETURN
ENDIF
!
DO j=1,NPARA
   IF(prange(j,1)<=prange(j,2)) THEN
      atry(j) = MIN(prange(j,2),MAX(prange(j,1),a(j)+da(j)))
   ELSE
      atry(     (j)) = a(     (j)) + da(j)
   ENDIF
ENDDO
!
CALL kuplot_mrqcof(MAXP, data_dim, data_data, data_weight, data_x, data_y,   &
               data_calc, atry, &
            NPARA, par_names, par_ind, prange, &
               MAXF, nfixed, fixed, fixed_ind, pf, &
            kupl_last, covar, da, chisq, LDERIV) !, funcs)
IF(ier_num/=0) THEN
   RETURN
ENDIF
!
!
IF(chisq < ochisq) THEN            ! Success, accept solution
   alamda = 0.5*alamda
   ochisq = chisq
   DO j=1,NPARA
      DO k=1,NPARA
         alpha(j,k) = covar(j,k)
      ENDDO
      beta(j) = da(j)
      a(     (j)) = atry(     (j))
   ENDDO
ELSE                               ! Failure, reject, keep old params and 
   alamda =  2.0*alamda            ! increment alamda
   chisq = ochisq                  ! Keep old chis squared
ENDIF
!
END SUBROUTINE kuplot_mrqmin
!
!*******************************************************************************
!
SUBROUTINE kuplot_mrqcof(MAXP, data_dim, data_data, data_weight, data_x, data_y,   &
               data_calc, &
                  params, NPARA, par_names, par_ind, prange, &
               MAXF, nfixed, fixed, fixed_ind, pf, &
               kupl_last, alpha, beta, chisq, LDERIV)
!
! Modified after NumRec 14.4
! Version for 2-d data
!
USE kuplot_fit6_set_theory
USE errlist_mod
!
IMPLICIT NONE
!
INTEGER                                             , INTENT(IN)  :: MAXP        ! Maximum parameter number
INTEGER         , DIMENSION(2)                      , INTENT(IN)  :: data_dim    ! Data dimensions
REAL            , DIMENSION(data_dim(1),data_dim(2)), INTENT(IN)  :: data_data   ! Observables
REAL            , DIMENSION(data_dim(1),data_dim(2)), INTENT(IN)  :: data_weight ! Observables
REAL            , DIMENSION(data_dim(1)           ) , INTENT(IN)  :: data_x      ! x-coordinates
REAL            , DIMENSION(data_dim(2)           ) , INTENT(IN)  :: data_y      ! y-coordinates
REAL            , DIMENSION(data_dim(1),data_dim(2)), INTENT(OUT) :: data_calc   ! Calculated data
INTEGER                                             , INTENT(IN)  :: NPARA       ! Number of refine parameters
REAL            , DIMENSION(MAXP)                   , INTENT(IN)  :: params      ! Current parameter values
CHARACTER(LEN=*), DIMENSION(MAXP)                   , INTENT(IN)  :: par_names   ! Parameter names
INTEGER         , DIMENSION(MAXP)                   , INTENT(IN)  :: par_ind     ! Index of param in KUPLOT list
REAL            , DIMENSION(MAXP, 2              )  , INTENT(IN)  :: prange      ! Allowed parameter range
INTEGER                                             , INTENT(IN)  :: MAXF        ! Fixed Parameter array size
INTEGER                                             , INTENT(IN)  :: nfixed      ! Number of fixed parameters
CHARACTER(LEN=*), DIMENSION(MAXF)                   , INTENT(IN)  :: fixed       ! Fixed Parameter names
INTEGER         , DIMENSION(MAXF)                   , INTENT(IN)  :: fixed_ind   ! Index of fixed param in KUPLOT list
REAL            , DIMENSION(MAXF)                   , INTENT(IN)  :: pf          ! Fixed Parameter array
INTEGER                                             , INTENT(IN)  :: kupl_last   ! Last KUPLOT DATA that are needed
REAL            , DIMENSION(NPARA, NPARA)           , INTENT(OUT) :: alpha
REAL            , DIMENSION(NPARA)                  , INTENT(OUT) :: beta
REAL                                                , INTENT(OUT) :: chisq
LOGICAL                                             , INTENT(IN)  :: LDERIV      ! Derivatives are needed?
!
!
INTEGER                  :: i,j, jj
INTEGER                  :: k, kk
INTEGER                  :: ix, iy
!
REAL                     :: xx    = 0.0     ! current x-coordinate
REAL                     :: yy    = 0.0     ! current y-coordinate
REAL                     :: ymod  = 0.0     ! zobs(calc)
REAL, DIMENSION(1:npara+nfixed) :: dyda            ! Derivatives
REAL                     :: sig2i = 0.0     ! sigma squared
REAL                     :: dy    = 0.0     ! Difference zobs(obs) - zobs(calc)
REAL                     :: wt    = 0.0     ! Weight 
!
INTEGER :: T_MAXP
INTEGER :: t_npara
REAL             , DIMENSION(:),   ALLOCATABLE :: t_params      ! Current parameter values
CHARACTER(LEN=60), DIMENSION(:),   ALLOCATABLE :: t_par_names   ! Parameter names
REAL             , DIMENSION(:,:), ALLOCATABLE :: t_prange      ! Allowed parameter range
LOGICAL          , DIMENSION(:),   ALLOCATABLE :: t_deriv       ! TRUE if para needs derivative
!
! Unscramble Parameters into KUPLOT sequence
!
T_MAXP  = npara+nfixed
t_npara = T_MAXP
ALLOCATE(t_params   (T_MAXP))
ALLOCATE(t_par_names(T_MAXP))
ALLOCATE(t_prange   (T_MAXP,2))
ALLOCATE(t_deriv    (T_MAXP))
DO i=1,npara                     ! Copy all free parameters
  j = par_ind(i)                 ! par_ind(i) is the location in KUPLOT
  t_params(j)    = params(i)
  t_par_names(j) = par_names(i)
  t_prange(j,1)  = prange(i,1)
  t_prange(j,2)  = prange(i,2)
  t_deriv(j)     = .TRUE.
ENDDO
DO i=1,nfixed                    ! Copy all fixed parameters
  j = fixed_ind(i)               ! par_ind(i) is the location in KUPLOT
  t_params(j)    = pf   (i)
  t_par_names(j) = fixed(i)
  t_prange(j,1)  = pf    (i)
  t_prange(j,2)  = pf    (i)
  t_deriv(j)     = .FALSE.
ENDDO
!
alpha(:,:) = 0.0
beta(:)    = 0.0
chisq      = 0.0
!
loopix: do ix=1, data_dim(1)
   xx = data_x(ix)
   loopiy: do iy=1, data_dim(2)
      yy = data_y(iy)
!
      CALL p_kuplot_theory(T_MAXP, ix, iy, xx, yy, t_npara, t_params, t_par_names,   &
                         t_prange, t_deriv, &
                         data_dim, data_data, data_weight, data_x, data_y, &
                         data_calc, &
                         kupl_last, &
                         ymod, dyda, LDERIV)
      IF(ier_num/=0) THEN
         RETURN
      ENDIF
!
      sig2i =1./(data_weight(ix,iy)*data_weight(ix,iy))
      dy = data_data(ix,iy) - ymod
      DO j=1, NPARA
         jj = par_ind(j)          ! KUPLOT parameter jj is Fit parameter j
         wt = dyda(jj)*sig2i
         DO k=1, j
            kk = par_ind(k)       ! KUPLOT parameter kk is Fit parameter k
            alpha(j,k) = alpha(j,k) + wt*dyda(kk)
         ENDDO
         beta(j) = beta(j) + dy*wt
      ENDDO
      chisq = chisq + dy*dy*sig2i
!
   ENDDO loopiy
ENDDO loopix
!
DO j=2, NPARA                   ! Fill in upper diagonal
   DO k=1,j-1
      alpha(k,j) = alpha(j,k)
   ENDDO
ENDDO
!
DEALLOCATE(t_params   )
DEALLOCATE(t_par_names)
DEALLOCATE(t_prange   )
DEALLOCATE(t_deriv    )
!
END SUBROUTINE kuplot_mrqcof
!
!*******************************************************************************
!
SUBROUTINE kuplot_rvalue(data_dim, data_data, data_weight, data_calc, rval, rexp, npar)
!
IMPLICIT NONE
!
INTEGER, DIMENSION(2), INTENT(IN) :: data_dim
REAL   , DIMENSION(data_dim(1),data_dim(2)), INTENT(IN) :: data_data
REAL   , DIMENSION(data_dim(1),data_dim(2)), INTENT(IN) :: data_weight
REAL   , DIMENSION(data_dim(1),data_dim(2)), INTENT(IN) :: data_calc
REAL   , INTENT(OUT) :: rval
REAL   , INTENT(OUT) :: rexp
INTEGER, INTENT(IN)  :: npar
!
INTEGER :: iix, iiy
REAL    :: sumz
REAL    :: sumn
REAL    :: wght
!
sumz = 0.0
sumn = 0.0
!
DO iiy = 1, data_dim(2)
   DO iix = 1, data_dim(1)
      IF(data_weight(iix,iiy)/=0.0) THEN
         wght = 1./(data_weight(iix,iiy))**2
      ELSE
         wght = 1.0
      ENDIF
      sumz = sumz + wght*(data_data(iix,iiy)-data_calc(iix,iiy))**2
      sumn = sumn + wght*(data_data(iix,iiy)                   )**2
   ENDDO
ENDDO
IF(sumn/=0.0) THEN
   rval = SQRT(sumz/sumn)
   rexp = SQRT((data_dim(1)*data_dim(2)-npar)/sumn)
ELSE
   rval = -1.0
   rexp = -1.0
ENDIF
!
END SUBROUTINE kuplot_rvalue
!
!*******************************************************************************
!
END MODULE kuplot_fit6
