!*****7*****************************************************************
!     This sublevel includes all commands and functions for the         
!     least square fits in KUPLOT.                                      
!*****7*****************************************************************
      SUBROUTINE do_fit (zei, lp) 
!                                                                       
!     Main fitting menu                                                 
!                                                                       
      USE doact_mod 
      USE errlist_mod 
      USE learn_mod 
      USE macro_mod 
      USE prompt_mod 
      USE config_mod 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 10) 
!                                                                       
      CHARACTER ( * ) zei 
      CHARACTER(1024) cpara (maxw) 
      CHARACTER(1024) line, zeile 
      CHARACTER(50) prom 
      CHARACTER(40) cdummy 
      CHARACTER(4) befehl 
      INTEGER lpara (maxw) 
      INTEGER ll, lp 
      INTEGER ianz, indxg, lbef 
      INTEGER maxpkt, maxzz 
      REAL werte (maxw) 
      LOGICAL flag (3), sel_func 
!                                                                       
      INTEGER len_str 
      LOGICAL str_comp 
!                                                                       
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
!                                                                       
!------ here starts sublevel fit                                        
!                                                                       
   10 CONTINUE 
!                                                                       
      prom = prompt (1:len_str (prompt) ) //'/fit' 
      CALL get_cmd (line, ll, befehl, lbef, zeile, lp, prom) 
      IF (ier_num.eq.0) then 
         IF (line.eq.' '.or.line (1:1) .eq.'#') goto 10 
!                                                                       
!------ search for "="                                                  
!                                                                       
         indxg = index (line, '=') 
      IF (indxg.ne.0.and..not. (str_comp (befehl, 'echo', 2, lbef, 4) ) &
     &.and..not. (str_comp (befehl, 'help', 2, lbef, 4) .or.str_comp (be&
     &fehl, '?   ', 2, lbef, 4) ) ) then                                
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
            CALL do_eval (zeile, lp) 
!                                                                       
!     exit 'exit'                                                       
!                                                                       
         ELSEIF (str_comp (befehl, 'exit', 2, lbef, 4) ) then 
            GOTO 9999 
!                                                                       
!     Define fit function 'func'                                        
!                                                                       
         ELSEIF (str_comp (befehl, 'func', 3, lbef, 4) ) then 
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
            IF (.not.sel_func) then 
               ier_num = - 25 
               ier_typ = ER_APPL 
            ELSE 
               IF (lni (ikfit) ) then 
                  CALL do_fit_z 
               ELSE 
                  CALL do_fit_y 
               ENDIF 
               CALL get_extrema 
            ENDIF 
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
      IF (ier_num.ne.0) then 
         CALL errlist 
         IF (ier_sta.ne.ER_S_LIVE) then 
            IF (lmakro.and.ier_sta.ne.ER_S_LIVE) then 
               CALL macro_close 
               prompt_status = PROMPT_ON 
            ENDIF 
            IF (lblock) then 
               ier_num = - 11 
               ier_typ = ER_COMM 
               RETURN 
            ENDIF 
            CALL no_error 
         ENDIF 
      ENDIF 
      GOTO 10 
!                                                                       
 9999 CONTINUE 
!                                                                       
 2000 FORMAT     (a) 
 2100 FORMAT     (1x,'Refinement mode: ',a) 
      END SUBROUTINE do_fit                         
!*****7*****************************************************************
      SUBROUTINE do_fit_fkt (zeile, lp) 
!+                                                                      
!     Set theory function                                               
!-                                                                      
      USE errlist_mod 
      USE config_mod 
      USE kuplot_mod 
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
      INTEGER ianz, i 
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
         IF (ianz.gt.1) then 
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
      USE prompt_mod 
      USE config_mod 
      USE kuplot_mod 
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
      USE config_mod 
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
      USE config_mod 
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
      USE config_mod 
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
      filname = fname (ikfit) 
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
      USE errlist_mod 
      USE config_mod 
      USE kuplot_mod 
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
         .ne.'ISQ'.and.cpara (1) (1:3) .ne.'DAT') then                  
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
      USE errlist_mod 
      USE prompt_mod 
      USE config_mod 
      USE kuplot_mod 
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
 1000 FORMAT     (' ---------- > Setting parameter p(',i2,') = ',g12.6, &
     &                   '  pinc = ',f4.1)                              
      END SUBROUTINE do_fit_par                     
!*****7*****************************************************************
      SUBROUTINE show_fit_para (idout) 
!+                                                                      
!     anzeigen der fit-parameter                                        
!-                                                                      
      USE config_mod 
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
         wictyp = 'w(i) = 1/dy(i) from data set' 
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
      USE config_mod 
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
      USE config_mod 
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
      USE config_mod 
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
      USE config_mod 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER(60) filname 
      REAL xx, f, df (maxpara) 
      INTEGER i, ii, jj, kk 
      INTEGER len_str 
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
      filname = fname (ikfit) 
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
      USE config_mod 
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
      xx = float (i) 
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
      filname = fname (ikfit) 
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
      USE config_mod 
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
      USE config_mod 
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
            wic = 1.0 / sig 
         ENDIF 
      ENDIF 
!                                                                       
      calc_wic = wic 
!                                                                       
      END FUNCTION calc_wic                         
!*****7*****************************************************************
      SUBROUTINE kupl_theory (xx, f, df, iwert) 
!                                                                       
      USE config_mod 
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
      ENDIF 
!                                                                       
      END SUBROUTINE kupl_theory                    
!*****7*****************************************************************
!       User defined fit function                                       
!*****7*****************************************************************
      SUBROUTINE setup_user (ianz, werte, maxw, cpara, lpara) 
!                                                                       
      USE errlist_mod 
      USE config_mod 
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
      INTEGER ianz, ip 
!                                                                       
      REAL berechne 
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
         dummy = berechne (cdummy, fit_lfunc + 2) 
!                                                                       
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
      END SUBROUTINE setup_user                     
!*****7*****************************************************************
      SUBROUTINE show_user (idout) 
!                                                                       
      USE config_mod 
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
!***********************************************************************
      SUBROUTINE theory_user (xx, f, df, i) 
!                                                                       
      USE param_mod 
      USE config_mod 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER(1024) cdummy 
      REAL xx, f, df (maxpara) 
      REAL pb, h, err 
      INTEGER i, ix, iy, ind 
!                                                                       
      REAL berechne, dfridr 
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
      f = berechne (cdummy, fit_lfunc + 2) 
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
      USE param_mod 
      USE config_mod 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER(1024) cdummy 
      REAL xx, berechne 
!                                                                       
      p (np1) = xx 
      cdummy = '('//fit_func (1:fit_lfunc) //')' 
      func = berechne (cdummy, fit_lfunc + 2) 
!                                                                       
      END FUNCTION func                             
!*****7*****************************************************************
!       GSAS profile functions                                          
!*****7*****************************************************************
      SUBROUTINE show_gsas (idout) 
!                                                                       
      USE config_mod 
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
      USE config_mod 
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
      USE config_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER nbank 
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
      REAL sind 
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
  997 CONTINUE 
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
      USE config_mod 
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
      USE config_mod 
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
      sint = p (iii + 1) * zpi * zz 
      ds0 = zpi * zz 
      ds2 = p (iii + 1) * zpi * (1.0 / p (iii + 4) + p (iii + 4) ) 
      ds3 = p (iii + 1) * zpi * ( - p (iii + 3) / p (iii + 4) **2 + p ( &
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
      USE config_mod 
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
      USE config_mod 
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
      USE config_mod 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      REAL o, sqpio, zz, sint, ds0, ds2, ds3, dsint 
      INTEGER idout, i, iii 
!                                                                       
      o = sqrt (4.0 * alog (2.0) ) 
      sqpio = sqrt (pi) / o * 0.5 
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
      USE config_mod 
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
      USE config_mod 
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
      USE config_mod 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER idout, i, iii 
      REAL zz, sint, ds0, ds2, ds3, dsint 
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
 1200 FORMAT     (3x,'p(',i2,') : backgr. 2 : ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
 1210 FORMAT     (3x,'p(',i2,') : backgr. 3 : ',g12.6,' +- ',g12.6,     &
     &                   4x,'pinc : ',f2.0)                             
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
      USE config_mod 
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
      USE config_mod 
      USE kuplot_mod 
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
      REAL tand 
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
      USE config_mod 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      REAL o, sqpio, sqpio2 
      REAL zz_x, zz_y, sint, ds0, ds3, ds4, ds6, ds7, dsint 
      INTEGER idout, i, iii 
!                                                                       
      o = sqrt (4.0 * alog (2.0) ) 
      sqpio = sqrt (pi) / o * 0.5 
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
      USE config_mod 
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
      USE config_mod 
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
      cosp = cos (rad * p (na + 5) ) 
      sinp = sin (rad * p (na + 5) ) 
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
         * sinp + (ry - p (na + 2) ) * cosp) * dfxs * rad+ ( - (rx - p (&
         na + 1) ) * cosp - (ry - p (na + 2) ) * sinp) * dfys * rad     
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
      USE config_mod 
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
      USE config_mod 
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
      USE config_mod 
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
      USE config_mod 
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
      USE config_mod 
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
      USE config_mod 
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
      USE config_mod 
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
      USE config_mod 
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
      USE config_mod 
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
      USE config_mod 
      USE kuplot_mod 
!                                                                       
!
      IMPLICIT NONE
!
      INTEGER :: i, j, k,l, ll, l1, nl   ! Loop indices
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
      REAL ww (maxarray) 
      REAL df (maxpara), dl (maxpara), dz (maxpara), pl (maxpara) 
      REAL pinc_old (maxpara) 
!                                                                       
      IF (lni (ikfit) ) then 
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
   10 CONTINUE 
      zwert = 0.0 
      DO k = 1, n 
      dz (k) = 0.0 
      DO l = k, n 
      cl (l, k) = 0.0 
      ENDDO 
      ENDDO 
      DO i = 1, m 
      xx = x (iix + i) 
      CALL kupl_theory (xx, f, df, i) 
      s = ww (iiw + i) - f 
      IF (wtyp.eq.'BCK') then 
         IF (ww (iiw + i) .ne. - 9999.) then 
            w (iiw + i) = exp ( - wval * s) 
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
      ENDDO 
      f = 1.0 
      IF (urf.le.0) GOTO 60
   40 CONTINUE 
      h = 0.0 
      DO k = 1, n 
      IF (dz (k) .ne.0.0) h = h + (dz (k) * dz (k) ) / (zwert * cl (k,  &
      k) )                                                              
      ENDDO 
      f = 1.000001 + urf * h 
   60 CONTINUE 
      IF (fstart) write (output_io, 5010) j, zwert, f 
      zalt = zl 
      zdif = zwert - zl 
      fend = f 
      IF (j.le.0) GOTO 80
   70 CONTINUE 
      IF (zwert - zl .ge. 0) GOTO 300
   80 CONTINUE 
      zl = zwert 
      DO k = 1, n 
      pl (k) = p (k) 
      ENDDO 
      DO k = 1, n 
      DO l = 1, k 
      l1 = l - 1 
      s = cl (k, l) 
      IF (l - k .eq. 0) GOTO 140
  100 IF (l1 .le. 0) GOTO 130
  110 DO i = 1, l1 
      s = s - cl (i, l) * cl (i, k) 
      ENDDO 
  130 cl (l, k) = s 
      ENDDO 
  140 IF (s.ne.0) GOTO 145
  142 dl (k) = 0.0 
      GOTO 190 
  145 s = s * f 
      IF (l1.le.0) GOTO 170
  150 DO i = 1, l1 
      ttt = cl (i, k) 
      cl (i, k) = ttt * dl (i) 
      s = s - ttt * cl (i, k) 
      ENDDO 
  170 IF (s.le.0) GOTO 190
  180 dl (k) = 1.0 / s 
  190 CONTINUE 
      ENDDO 
      IF (j - ncycle .gt. 0) GOTO 300 
  200 j = j + 1 
      IF (n - 1 .le.0) GOTO 230
  210 DO l = 2, n 
      l1 = l - 1 
      DO k = 1, l1 
      dz (l) = dz (l) - cl (k, l) * dz (k) 
      ENDDO 
      ENDDO 
  230 dz (n) = dz (n) * dl (n) 
      IF (n - 1 .le.0) GOTO 260
  240 DO nl = 2, n 
      l = n - nl + 1 
      l1 = l + 1 
      dz (l) = dz (l) * dl (l) 
      DO k = l1, n 
      dz (l) = dz (l) - cl (l, k) * dz (k) 
      ENDDO 
      ENDDO 
  260 DO k = 1, n 
      pinc (k) = 0.01 * dz (k) 
      p (k) = p (k) - dz (k) 
      ENDDO 
!                                                                       
      GOTO 10 
!                                                                       
  300 g = 1.0 
      h = float (m - n) 
      IF (h * zl.gt.0.0) g = zl / h 
      DO k = 1, n 
      dl (k) = dl (k) * g 
      cl (k, k) = 1.0 
      ENDDO 
      IF (n - 1 .le.0) GOTO 350
  320 DO l = 2, n 
      l1 = l - 1 
      DO k = 1, l1 
      s = 0.0 
      DO i = k, l1 
      s = s - cl (i, l) * cl (i, k) 
      ENDDO 
      cl (l, k) = s 
      ENDDO 
      ENDDO 
  350 DO k = 1, n 
      DO l = k, n 
      s = 0.0 
      IF (l - k .ne. 0) GOTO 380
  360 DO i = l, n 
      ttt = cl (i, k) 
      cl (i, k) = ttt * dl (i) 
      s = s + ttt * cl (i, k) 
      ENDDO 
      GOTO 400 
  380 DO i = l, n 
      s = s + cl (i, l) * cl (i, k) 
      ENDDO 
  400 CONTINUE 
      cl (l, k) = s 
      ENDDO 
      ENDDO 
      DO l = 1, n 
      IF (dl (l) .ne. 0 ) GOTO 502
  501 dp (l) = 1.0 
      dp (l) = 0.0 
      GOTO 570 
  502 f = sqrt (cl (l, l) ) 
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
  570 CONTINUE 
      ENDDO 
!                                                                       
      DO k = 1, n 
      p (k) = pl (k) 
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
            w (iiw + i) = exp ( - wval * s) 
         ELSE 
            w (iiw + i) = 0.0 
         ENDIF 
      ENDIF 
      sum3 = sum3 + w (iiw + i) * h * h 
      sum4 = sum4 + w (iiw + i) * ww (iiw + i) * ww (iiw + i) 
      ENDDO 
      r4 = sqrt (sum3 / sum4) 
      re = sqrt ( (m - nf) / sum4) 
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
