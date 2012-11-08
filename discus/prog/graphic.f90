!*****7*****************************************************************
      SUBROUTINE do_niplps (linverse) 
!-                                                                      
!     This sublevel contains all routines used to write the output      
!     of the Fourier transform/Patterson to an output file in           
!     various formats.                                                  
!+                                                                      
      USE config_mod 
      USE diffuse_mod 
      USE output_mod 
      USE powder, ONLY: powder_out
      IMPLICIT none 
!                                                                       
      include'doact.inc' 
      include'errlist.inc' 
      include'learn.inc' 
      include'macro.inc' 
      include'prompt.inc' 
!                                                                       
      INTEGER maxp, maxw, iff 
      PARAMETER (maxp = 11, maxw = 1, iff = 2) 
!                                                                       
      CHARACTER(5) befehl 
      CHARACTER(50) prom 
      CHARACTER(14) cvalue (0:6) 
      CHARACTER(22) cgraphik (0:8) 
      CHARACTER(1024) infile 
      CHARACTER(1024) zeile 
      CHARACTER(1024) line, cpara (maxp) 
      INTEGER lpara (maxp) 
      INTEGER ix, iy, ianz, value, lp, length, lbef 
      INTEGER indxg 
      LOGICAL laver, lread, linverse 
      REAL xmin, ymin, xmax, ymax, werte (maxp) 
!                                                                       
      INTEGER len_str 
      LOGICAL str_comp 
!                                                                       
      DATA cgraphik / 'Standard', 'Postscript', 'Pseudo Grey Map', 'Gnup&
     &lot', 'Portable Any Map', 'Powder Pattern', 'SHELX', 'SHELXL List &
     &5', 'SHELXL List 5 real HKL' /                                    
      DATA cvalue / 'undefined     ', 'Intensity     ', 'Amplitude     '&
     &, 'Phase angle   ', 'Real Part     ', 'Imaginary Part', 'Random Ph&
     &ase  ' /                                                          
!                                                                       
      DATA value / 1 / 
      DATA laver / .false. / 
!                                                                       
      zmin = ps_low * diffumax 
      zmax = ps_high * diffumax 
   10 CONTINUE 
!                                                                       
      CALL no_error 
!                                                                       
      prom = prompt (1:len_str (prompt) ) //'/output' 
      CALL get_cmd (line, length, befehl, lbef, zeile, lp, prom) 
      IF (ier_num.eq.0) then 
         IF (line (1:1) .eq.' '.or.line (1:1) .eq.'#') goto 10 
!                                                                       
!     search for "="                                                    
!                                                                       
         indxg = index (line, '=') 
      IF (indxg.ne.0.and..not. (str_comp (befehl, 'echo', 2, lbef, 4) ) &
     &.and..not. (str_comp (befehl, 'syst', 2, lbef, 4) ) .and..not. (st&
     &r_comp (befehl, 'help', 2, lbef, 4) .or.str_comp (befehl, '?   ', &
     &2, lbef, 4) ) ) then                                              
!                                                                       
!     --evaluatean expression and assign the value to a variabble       
!                                                                       
            CALL do_math (line, indxg, length) 
         ELSE 
!                                                                       
!------ execute a macro file                                            
!                                                                       
            IF (befehl (1:1) .eq.'@') then 
               IF (length.ge.2) then 
                  CALL file_kdo (line (2:length), length - 1) 
               ELSE 
                  ier_num = - 13 
                  ier_typ = ER_MAC 
               ENDIF 
!                                                                       
!     continues a macro 'continue'                                      
!                                                                       
            ELSEIF (str_comp (befehl, 'continue', 1, lbef, 8) ) then 
               CALL macro_continue (zeile, lp) 
!                                                                       
!------ Echo a string, just for interactive check in a macro 'echo'     
!                                                                       
            ELSEIF (str_comp (befehl, 'echo', 2, lbef, 4) ) then 
               CALL echo (zeile, lp) 
!                                                                       
!     Evaluate an expression, just for interactive check 'eval'         
!                                                                       
            ELSEIF (str_comp (befehl, 'eval', 2, lbef, 4) ) then 
               CALL do_eval (zeile, lp) 
!                                                                       
!     Terminate output 'exit'                                           
!                                                                       
            ELSEIF (str_comp (befehl, 'exit', 2, lbef, 4) ) then 
               GOTO 9999 
!                                                                       
!     Determine format for output 'format'                              
!                                                                       
            ELSEIF (str_comp (befehl, 'form', 1, lbef, 4) ) then 
               CALL get_params (zeile, ianz, cpara, lpara, maxp, lp) 
               IF (ier_num.eq.0) then 
                  IF (ianz.eq.1.or.ianz.eq.2) then 
!                                                                       
!     ------Switch output type to GNUPLOT 'gnup'                        
!                                                                       
                     IF (str_comp (cpara (1) , 'gnup', 1, lpara (1) , 4)&
                     ) then                                             
                        ityp = 3 
!                                                                       
!     ------Switch output type to pgm 'pgm'                             
!                                                                       
                     ELSEIF (str_comp (cpara (1) , 'pgm ', 2, lpara (1) &
                     , 4) ) then                                        
                        ityp = 2 
!                                                                       
!     ------Switch output type to postscript 'post'                     
!                                                                       
                     ELSEIF (str_comp (cpara (1) , 'post', 3, lpara (1) &
                     , 4) ) then                                        
                        ityp = 1 
!                                                                       
!     ------Switch output type to powder pattern 'powd'                 
!                                                                       
                     ELSEIF (str_comp (cpara (1) , 'powd', 3, lpara (1) &
                     , 4) ) then                                        
                        ityp = 5 
                        IF (ianz.eq.2) then 
                           IF (str_comp (cpara (2) , 'tth', 2, lpara (2)&
                           , 3) ) then                                  
                              cpow_form = 'tth' 
                           ELSEIF (str_comp (cpara (2) , 'q', 1, lpara (&
                           2) , 1) ) then                               
      cpow_form = 'q  ' 
                           ELSEIF (str_comp (cpara (2) , 'stl', 2,      &
                           lpara (2) , 3) ) then                        
                              cpow_form = 'stl' 
                           ELSEIF (str_comp (cpara (2) , 'dst', 2,      &
                           lpara (2) , 3) ) then                        
                              cpow_form = 'dst' 
                           ELSEIF (str_comp (cpara (2) , 'lop', 2,      &
                           lpara (2) , 3) ) then                        
                              cpow_form = 'lop' 
                           ELSE 
                              ier_num = - 6 
                              ier_typ = ER_COMM 
                           ENDIF 
                        ENDIF 
!                                                                       
!     ------Switch output type to ppm 'ppm'                             
!                                                                       
                     ELSEIF (str_comp (cpara (1) , 'ppm ', 2, lpara (1) &
                     , 4) ) then                                        
                        ityp = 4 
!                                                                       
!     ------Switch output type to standard  'stan'                      
!                                                                       
                     ELSEIF (str_comp (cpara (1) , 'stan', 2, lpara (1) &
                     , 4) ) then                                        
                        ityp = 0 
!                                                                       
!     ------Switch output type to Shelx 'shel', or 'hklf4'              
!                                                                       
                     ELSEIF (str_comp (cpara (1) , 'shel', 2, lpara (1) &
                     , 4) ) then                                        
                        ityp = 6 
                     ELSEIF (str_comp (cpara (1) , 'hklf4', 2, lpara (1)&
                     , 5) ) then                                        
                        ityp = 6 
!                                                                       
!     ------Switch output type to Shelx LIST 5   'list5'                
!                                                                       
                     ELSEIF (str_comp (cpara (1) , 'list5', 2, lpara (1)&
                     , 5) ) then                                        
                        ityp = 7 
!                                                                       
!     ------Switch output type to Shelx LIST 5   'list9'                
!                                                                       
                     ELSEIF (str_comp (cpara (1) , 'list9', 2, lpara (1)&
                     , 5) ) then                                        
                        ityp = 8 
                     ELSE 
                        ier_num = - 9 
                        ier_typ = ER_APPL 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ENDIF 
!                                                                       
!     help on output 'help'                                             
!                                                                       
      ELSEIF (str_comp (befehl, 'help', 2, lbef, 4) .or.str_comp (befehl&
     &, '?   ', 1, lbef, 4) ) then                                      
               IF (str_comp (zeile, 'errors', 2, lp, 6) ) then 
                  lp = lp + 7 
                  CALL do_hel ('discus '//zeile, lp) 
               ELSE 
                  lp = lp + 14 
                  CALL do_hel ('discus output '//zeile, lp) 
               ENDIF 
!                                                                       
!     read an old output file (only for standard file type' 'inpu'      
!                                                                       
            ELSEIF (str_comp (befehl, 'inpu', 1, lbef, 4) ) then 
               CALL get_params (zeile, ianz, cpara, lpara, maxp, lp) 
               IF (ier_num.eq.0) then 
                  infile = cpara (1) 
                  lread = .true. 
                  CALL oeffne (1, infile, 'old', lread) 
                  IF (ier_num.eq.0) then 
                     READ (1, * ) out_inc (1), out_inc (2) 
                     READ (1, * ) xmin, xmax, ymin, ymax 
                     READ (1, * ) zmax 
                     zmin = zmax 
                     BACKSPACE (1) 
!                                                                       
                     DO iy = 1, out_inc (2) 
                     READ (1, * ) (dsi ( (ix - 1) * out_inc (2) + iy),  &
                     ix = 1, out_inc (1) )                              
                     DO ix = 1, out_inc (1) 
                     zmax = max (zmax, dsi ( (ix - 1) * out_inc (2)     &
                     + iy) )                                            
                     zmin = min (zmin, dsi ( (ix - 1) * out_inc (2)     &
                     + iy) )                                            
                     ENDDO 
                     ENDDO 
                     WRITE (output_io, 1015) zmin, zmax 
                     READ ( *, *, end = 20) zmin, zmax 
   20                CONTINUE 
                  ENDIF 
                  CLOSE (1) 
               ENDIF 
!                                                                       
!     define name of output file 'outf'                                 
!                                                                       
            ELSEIF (str_comp (befehl, 'outf', 1, lbef, 4) ) then 
               CALL get_params (zeile, ianz, cpara, lpara, maxp, lp) 
               IF (ier_num.eq.0) then 
                  CALL do_build_name (ianz, cpara, lpara, werte, maxp,  &
                  1)                                                    
                  IF (ier_num.eq.0) then 
                     outfile = cpara (1) 
                  ENDIF 
               ENDIF 
!                                                                       
!     write output file 'run'                                           
!                                                                       
            ELSEIF (str_comp (befehl, 'run ', 1, lbef, 4) ) then 
               CALL set_output (linverse) 
               IF (ityp.eq.0) then 
                  CALL do_output (value, laver) 
               ELSEIF (ityp.eq.1) then 
                  CALL do_post (value, laver) 
               ELSEIF (ityp.eq.2) then 
                  CALL do_pgm (value, laver) 
               ELSEIF (ityp.eq.3) then 
                  CALL do_output (value, laver) 
               ELSEIF (ityp.eq.4) then 
                  CALL do_ppm (value, laver) 
               ELSEIF (ityp.eq.5) then 
                  CALL powder_out 
               ELSEIF (ityp.eq.6) then 
                  CALL do_output (value, laver) 
               ELSEIF (ityp.eq.7) then 
                  CALL do_output (value, laver) 
               ELSEIF (ityp.eq.8) then 
                  CALL do_output (value, laver) 
               ELSE 
                  ier_num = - 9 
                  ier_typ = ER_APPL 
               ENDIF 
!                                                                       
!     Show current settings for output 'show'                           
!                                                                       
            ELSEIF (str_comp (befehl, 'show', 2, lbef, 4) ) then 
               WRITE (output_io, 3000) outfile 
               IF (ityp.lt.0.or.8.lt.ityp) then 
                  WRITE (output_io, * ) 'ityp undefiniert ', ityp 
               ELSEIF (ityp.eq.5) then 
                  WRITE (output_io, 3130) cgraphik (ityp), cpow_form 
               ELSE 
                  WRITE (output_io, 3100) cgraphik (ityp) 
                  IF (laver) then 
                     WRITE (output_io, 3110) '<'//cvalue (value) //'>' 
                  ELSE 
                     WRITE (output_io, 3110) cvalue (value) 
                  ENDIF 
               ENDIF 
               WRITE (output_io, 3060) braggmin, braggmax, diffumin,    &
               diffumax, diffuave, diffusig                             
               WRITE (output_io, 3080) 100.0 * ps_high, zmax 
               WRITE (output_io, 3090) 100.0 * ps_low, zmin 
!                                                                       
!-------Operating System Kommandos 'syst'                               
!                                                                       
            ELSEIF (str_comp (befehl, 'syst', 2, lbef, 4) ) then 
               IF (zeile.ne.' ') then 
                  CALL do_operating (zeile (1:lp), lp) 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
!                                                                       
!     Set threshold for intensity written to bitmaps 'thresh'           
!                                                                       
            ELSEIF (str_comp (befehl, 'thre', 2, lbef, 4) ) then 
               CALL get_params (zeile, ianz, cpara, lpara, maxp, lp) 
               IF (ier_num.eq.0) then 
                  IF (ianz.eq.2) then 
                     IF (str_comp (cpara (1) , 'high', 1, lpara (1) , 4)&
                     ) then                                             
                        CALL del_params (1, ianz, cpara, lpara, maxp) 
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxp)                                           
                        IF (ier_num.eq.0) then 
                           ps_high = werte (1) * 0.01 
                           zmax = diffumax * ps_high 
                        ENDIF 
                     ELSEIF (str_comp (cpara (1) , 'low', 1, lpara (1) ,&
                     3) ) then                                          
                        CALL del_params (1, ianz, cpara, lpara, maxp) 
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxp)                                           
                        IF (ier_num.eq.0) then 
                           ps_low = werte (1) * 0.01 
                           zmin = diffumax * ps_low 
                        ENDIF 
                     ELSEIF (str_comp (cpara (1) , 'sigma', 1, lpara (1)&
                     , 5) ) then                                        
                        CALL del_params (1, ianz, cpara, lpara, maxp) 
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxp)                                           
                        IF (ier_num.eq.0) then 
                           zmin = max (diffumin, diffuave-werte (1)     &
                           * diffusig)                                  
                           zmax = min (diffumax, diffuave+werte (1)     &
                           * diffusig)                                  
                           IF (diffumax.ne.0) then 
                              ps_high = zmax / diffumax 
                              ps_low = zmin / diffumax 
                           ELSE 
                              ps_high = 0.0 
                              ps_low = 0.0 
                           ENDIF 
                        ENDIF 
                     ELSEIF (str_comp (cpara (1) , 'zmax', 3, lpara (1) &
                     , 4) ) then                                        
                        CALL del_params (1, ianz, cpara, lpara, maxp) 
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxp)                                           
                        IF (ier_num.eq.0) then 
                           zmax = werte (1) 
                           IF (diffumax.ne.0) then 
                              ps_high = zmax / diffumax 
                           ELSE 
                              ps_high = 0.0 
                           ENDIF 
                        ENDIF 
                     ELSEIF (str_comp (cpara (1) , 'zmin', 3, lpara (1) &
                     , 4) ) then                                        
                        CALL del_params (1, ianz, cpara, lpara, maxp) 
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxp)                                           
                        IF (ier_num.eq.0) then 
                           zmin = werte (1) 
                           IF (diffumax.ne.0) then 
                              ps_low = zmin / diffumax 
                           ELSE 
                              ps_low = 0.0 
                           ENDIF 
                        ENDIF 
                     ELSE 
                        ier_num = - 11 
                        ier_typ = ER_APPL 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ENDIF 
!                                                                       
!     Define output value 'value'                                       
!                                                                       
            ELSEIF (str_comp (befehl, 'valu', 1, lbef, 4) ) then 
               CALL get_params (zeile, ianz, cpara, lpara, maxp, lp) 
               IF (ier_num.eq.0) then 
!------ ----Check if we want the average values <F> ?                   
                  IF (cpara (1) (1:1) .eq.'<') then 
                     ix = 2 
                     laver = .true. 
                  ELSE 
                     ix = 1 
                     laver = .false. 
                  ENDIF 
!     ----Calculate intensity 'intensity'                               
                  IF (cpara (1) (ix:ix + 1) .eq.'in') then 
                     value = 1 
!     ----Calculate amplitude 'amplitude'                               
                  ELSEIF (cpara (1) (ix:ix) .eq.'a') then 
                     value = 2 
!     ----Calculate phase 'phase'                                       
                  ELSEIF (cpara (1) (ix:ix) .eq.'p') then 
                     IF (ianz.eq.1) then 
                        value = 3 
                     ELSEIF (ianz.eq.2.and.cpara (2) (1:1) .eq.'r')     &
                     then                                               
                        value = 6 
                     ENDIF 
!     ----Calculate real part 'real'                                    
                  ELSEIF (cpara (1) (ix:ix) .eq.'r') then 
                     value = 4 
!     ----Calculate imaginary part 'imaginary'                          
                  ELSEIF (cpara (1) (ix:ix + 1) .eq.'im') then 
                     value = 5 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                     value = 0 
                  ENDIF 
!------ ----check lots and allowed output                               
                  IF (nlots.ne.1.and.value.ne.1.and..not.laver) then 
                     ier_num = - 60 
                     ier_typ = ER_APPL 
                     value = 0 
                  ENDIF 
               ENDIF 
!                                                                       
!------  -waiting for user input                                        
!                                                                       
            ELSEIF (str_comp (befehl, 'wait', 3, lbef, 4) ) then 
               CALL do_input (zeile, lp) 
!                                                                       
!------ no valid subcommand found                                       
!                                                                       
            ELSE 
               ier_num = - 8 
               ier_typ = ER_COMM 
            ENDIF 
         ENDIF 
      ENDIF 
      IF (ier_num.ne.0) then 
         CALL errlist 
         IF (ier_sta.ne.ER_S_LIVE) then 
            IF (lmakro) then 
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
 9999 CONTINUE 
!                                                                       
 1015 FORMAT ( /1x,'Z-MIN = ',G20.6,/,1x,'Z-MAX = ',G20.6,//            &
     &                     1x,'Give new values zmin, zmax    : ',$)     
 3000 FORMAT( ' Output file                  : ',a) 
 3060 FORMAT(/' Bragg minimum                : ',g12.6/                 &
     &        ' Bragg maximum                : ',g12.6/                 &
     &        ' Diffuse minimum              : ',g12.6/                 &
     &        ' Diffuse maximum              : ',g12.6/                 &
     &        ' Diffuse average intensity    : ',g12.6/                 &
     &        ' Diffuse intensity sigma      : ',g12.6)                 
 3080 FORMAT(/' Maximum value for BITMAP'/                              &
     &        ' in % of highest diffuse value'/                         &
     &        ' and absolute                 : ',2x,f9.4,2x,g12.6)      
 3090 FORMAT( ' Minimum value for BITMAP'/                              &
     &        ' in % of highest diffuse value'/                         &
     &        ' and absolute                 : ',2x,f9.4,2x,g12.6)      
 3100 FORMAT( ' Graphicsformat               : ',A) 
 3130 FORMAT( ' Graphicsformat               : ',A,A) 
 3110 FORMAT( ' Output value                 : ',A) 
      END SUBROUTINE do_niplps                      
!*****7*****************************************************************
      SUBROUTINE do_post (value, laver) 
!-                                                                      
!     Writes a POSTSCRIPT file                                          
!+                                                                      
      USE config_mod 
      USE diffuse_mod 
      USE output_mod 
      IMPLICIT none 
!                                                                       
      include'envir.inc' 
      include'errlist.inc' 
!                                                                       
      INTEGER maxcol 
      PARAMETER (maxcol = 256) 
!                                                                       
      CHARACTER(6) cout (maxqxy) 
      CHARACTER(6) cfarb (maxcol) 
      INTEGER i, ix, iy, iqqq, k, value 
      LOGICAL lread, laver 
      REAL qqq 
!                                                                       
      REAL qval 
!                                                                       
!     Check whether data are 2-dimensional                              
!                                                                       
      IF (.not. (out_inc (1) .gt.1.and.out_inc (2) .gt.1) ) then 
         ier_num = - 50 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
!                                                                       
!-------Farbtabelle einlesen                                            
!                                                                       
      lread = .true. 
      CALL oeffne (2, colorfile, 'old', lread) 
      IF (ier_num.ne.0) return 
      DO i = 1, 255 
      READ (2, 100, end = 20) cfarb (i) 
  100 FORMAT      (1x,a6) 
      ENDDO 
   20 CONTINUE 
      CLOSE (2) 
      cfarb (256) = 'ffffff' 
!                                                                       
      lread = .false. 
      CALL oeffne (2, outfile, 'unknown', lread) 
      IF (ier_num.ne.0) return 
!                                                                       
      WRITE (2, 1111) '%!PS-Adobe-2.0' 
      WRITE (2, 1111) '%%Creator: DISCUS, Version 3.0' 
      WRITE (2, 1111) '50  150 translate' 
      WRITE (2, 1111) '288 288 scale' 
      WRITE (2, 2000) nint (3.0 * out_inc (1) ), out_inc (1), out_inc ( &
      2), i, 8, out_inc (1), 0, 0, out_inc (2), 0, 0                    
!                                                                       
      DO iy = 1, out_inc (2) 
      DO ix = 1, out_inc (1) 
      k = (ix - 1) * out_inc (2) + iy 
      qqq = qval (k, value, ix, iy, laver) 
      IF (qqq.lt.zmin) then 
         qqq = zmin 
      ELSEIF (qqq.gt.zmax) then 
         qqq = zmax 
      ENDIF 
      iqqq = nint ( (maxcol - 2) * (qqq - zmin) / (zmax - zmin) )       &
      + 1                                                               
      WRITE (cout (ix), 1111) cfarb (iqqq) 
      DO i = 1, 6 
      IF (cout (ix) (i:i) .eq.' '.or.cout (ix) (i:i) .eq.' ') cout (ix) &
      (i:i) = '0'                                                       
      ENDDO 
      ENDDO 
      WRITE (2, 5000) (cout (ix), ix = 1, out_inc (1) ) 
      ENDDO 
      WRITE (2, 1111) 'showpage' 
!                                                                       
      CLOSE (2) 
!                                                                       
 1111 FORMAT (a) 
 2000 FORMAT ('/DataString ',I4,' string def'/3(I3,1X),                 &
     &        ' [ ',6(I3,1X),']'/'{'/                                   &
     &        '  currentfile DataString readhexstring pop'/             &
     &        ' }  false 3 colorimage')                                 
 2200 FORMAT (Z8) 
 5000 FORMAT (10A6) 
!                                                                       
      END SUBROUTINE do_post                        
!*****7*****************************************************************
      SUBROUTINE do_pgm (value, laver) 
!-                                                                      
!     Writes the data in a PGM format                                   
!-                                                                      
      USE config_mod 
      USE diffuse_mod 
      USE output_mod 
      IMPLICIT none 
!                                                                       
      include'errlist.inc' 
!                                                                       
      INTEGER maxcol 
      PARAMETER (maxcol = 255) 
!                                                                       
      INTEGER iqqq (maxqxy), ncol, ix, iy, k, value 
      LOGICAL lread, laver 
      REAL qqq 
!                                                                       
      REAL qval 
!                                                                       
!     Check whether data are 2-dimensional                              
!                                                                       
      IF (.not. (out_inc (1) .gt.1.and.out_inc (2) .gt.1) ) then 
         ier_num = - 50 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
!                                                                       
      lread = .false. 
      ncol = maxcol 
!                                                                       
      CALL oeffne (2, outfile, 'unknown', lread) 
      IF (ier_num.ne.0) return 
!                                                                       
      WRITE (2, 1111) 'P2' 
      WRITE (2, 2000) out_inc (1), out_inc (2), ncol 
!                                                                       
      DO iy = out_inc (2), 1, - 1 
      DO ix = 1, out_inc (1) 
      k = (ix - 1) * out_inc (2) + iy 
      qqq = qval (k, value, ix, iy, laver) 
      IF (qqq.lt.zmin) then 
         qqq = zmin 
      ELSEIF (qqq.gt.zmax) then 
         qqq = zmax 
      ENDIF 
      iqqq (ix) = nint (float (ncol - 1) * (qqq - zmin) / (zmax - zmin) &
      )                                                                 
      ENDDO 
      WRITE (2, 5000) (iqqq (ix), ix = 1, out_inc (1) ) 
      ENDDO 
!                                                                       
      CLOSE (2) 
!                                                                       
 1111 FORMAT     (a) 
 2000 FORMAT    (1x,2i4/1x,i8) 
 5000 FORMAT     (7i8) 
!                                                                       
      END SUBROUTINE do_pgm                         
!*****7*****************************************************************
      SUBROUTINE do_ppm (value, laver) 
!+                                                                      
!     Writes the data in a PGM format                                   
!-                                                                      
      USE config_mod 
      USE diffuse_mod 
      USE output_mod 
      IMPLICIT none 
!                                                                       
      include'envir.inc' 
      include'errlist.inc' 
!                                                                       
      INTEGER maxcol 
      PARAMETER (maxcol = 255) 
!                                                                       
      INTEGER iqqq (maxqxy), ncol, ix, iy, k, value 
      INTEGER icolor (maxcol, 3) 
      INTEGER i, j 
      LOGICAL lread, laver 
      REAL qqq 
!                                                                       
      REAL qval 
!                                                                       
      CHARACTER(6) cfarb (256) 
!                                                                       
!                                                                       
!     Check whether data are 2-dimensional                              
!                                                                       
      IF (.not. (out_inc (1) .gt.1.and.out_inc (2) .gt.1) ) then 
         ier_num = - 50 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
!                                                                       
!-------Farbtabelle einlesen                                            
!                                                                       
      CALL set_colorfeld (cfarb) 
      lread = .true. 
      CALL oeffne (2, colorfile, 'old', lread) 
      IF (ier_num.ne.0) return 
      DO i = 1, 255 
      READ (2, 100, end = 20) (icolor (i, j), j = 1, 3) 
  100 FORMAT      (1x,3z2) 
      ENDDO 
   20 CONTINUE 
      CLOSE (2) 
!                                                                       
      icolor (255, 1) = 255 
      icolor (255, 2) = 255 
      icolor (255, 3) = 255 
!                                                                       
      lread = .false. 
      CALL oeffne (2, outfile, 'unknown', lread) 
      IF (ier_num.ne.0) return 
!                                                                       
      ncol = maxcol 
      WRITE (2, 1111) 'P3' 
      WRITE (2, 2000) out_inc (1), out_inc (2), ncol 
!                                                                       
      DO iy = out_inc (2), 1, - 1 
      DO ix = 1, out_inc (1) 
      k = (ix - 1) * out_inc (2) + iy 
      qqq = qval (k, value, ix, iy, laver) 
      IF (qqq.lt.zmin) then 
         qqq = zmin 
      ELSEIF (qqq.gt.zmax) then 
         qqq = zmax 
      ENDIF 
      iqqq (ix) = nint (float (ncol - 1) * (qqq - zmin) / (zmax - zmin) &
      ) + 1                                                             
      ENDDO 
      WRITE (2, 5000) ( (icolor (iqqq (ix), j), j = 1, 3), ix = 1,      &
      out_inc (1) )                                                     
      ENDDO 
!                                                                       
      CLOSE (2) 
!                                                                       
 1111 FORMAT     (a) 
 2000 FORMAT    (1x,2i4/1x,i8) 
 5000 FORMAT     (15i4) 
!                                                                       
      END SUBROUTINE do_ppm                         
!*****7*****************************************************************
      SUBROUTINE set_colorfeld (cfarb) 
!+                                                                      
!     This routine sets the pseudo color color map                      
!-                                                                      
      INTEGER maxcol 
      PARAMETER (maxcol = 256) 
!                                                                       
      CHARACTER ( * ) cfarb (maxcol) 
      CHARACTER(20) ccc (256) 
      INTEGER rgb (3) 
!                                                                       
      cfarb (1) = '000000' 
!                                                                       
      DO ifarb = 2, 256 
      rh = 0.1 + float (ifarb - 1) / 283.0 
      rh = 6.0 * rh 
      i = int (rh) 
      rf = rh - float (i) 
      rp = 0.0 
      rq = 1.0 - rf 
      rt = (1.0 - (1.0 - rf) ) 
!                                                                       
      IF (rt.gt.1.0) rt = 1.0 
      IF (rp.gt.1.0) rp = 1.0 
      IF (rq.gt.1.0) rq = 1.0 
!                                                                       
      IF (i.eq.0) then 
         rgb (1) = int (0.5 + 1. * 255.0) 
         rgb (2) = int (0.5 + rt * 255.0) 
         rgb (3) = int (0.5 + rp * 255.0) 
      ELSEIF (i.eq.1) then 
         rgb (1) = int (0.5 + rq * 255.0) 
         rgb (2) = int (0.5 + 1. * 255.0) 
         rgb (3) = int (0.5 + rp * 255.0) 
      ELSEIF (i.eq.2) then 
         rgb (1) = int (0.5 + rp * 255.0) 
         rgb (2) = int (0.5 + 1. * 255.0) 
         rgb (3) = int (0.5 + rt * 255.0) 
      ELSEIF (i.eq.3) then 
         rgb (1) = int (0.5 + rp * 255.0) 
         rgb (2) = int (0.5 + rq * 255.0) 
         rgb (3) = int (0.5 + 1. * 255.0) 
      ELSEIF (i.eq.4) then 
         rgb (1) = int (0.5 + rt * 255.0) 
         rgb (2) = int (0.5 + rp * 255.0) 
         rgb (3) = int (0.5 + 1. * 255.0) 
      ELSEIF (i.eq.5) then 
         rgb (1) = int (0.5 + 1. * 255.0) 
         rgb (2) = int (0.5 + rp * 255.0) 
         rgb (3) = int (0.5 + rq * 255.0) 
      ENDIF 
!                                                                       
      WRITE (ccc (ifarb), 1000) (rgb (j), j = 3, 1, - 1) 
      DO ii = 1, 6 
      IF (ccc (ifarb) (ii:ii) .eq.' ') ccc (ifarb) (ii:ii) = '0' 
      ENDDO 
!                                                                       
      WRITE (cfarb (ifarb), 1000) (rgb (j), j = 3, 1, - 1) 
      DO ii = 1, 6 
      IF (cfarb (ifarb) (ii:ii) .eq.' ') cfarb (ifarb) (ii:ii) = '0' 
      ENDDO 
      ENDDO 
      DO ifarb = 256, 2, - 1 
      WRITE (44, 2000) ccc (ifarb) 
      ENDDO 
      CLOSE (44) 
!                                                                       
 1000 FORMAT     (3(z2)) 
 2000 FORMAT    ('#',a6,'00') 
      END SUBROUTINE set_colorfeld                  
!*****7*****************************************************************
      SUBROUTINE do_output (value, laver) 
!-                                                                      
!     Writes output in standard or GNUPLOT format                       
!+                                                                      
      USE config_mod 
      USE crystal_mod 
      USE diffuse_mod 
      USE output_mod 
      IMPLICIT none 
!                                                                       
       
      include'errlist.inc' 
      include'prompt.inc' 
!                                                                       
      INTEGER iff 
      PARAMETER (iff = 2) 
!                                                                       
      INTEGER HKLF4, LIST5, LIST9 
      PARAMETER (HKLF4 = 6, LIST5 = 7, LIST9 = 8) 
!                                                                       
      INTEGER extr_ima, i, j, k, value 
      LOGICAL lread, laver 
      REAL h (3) 
      REAL sq, qq, out_fac 
!                                                                       
      INTEGER shel_inc (2) 
      INTEGER shel_value 
      REAL shel_eck (3, 3) 
      REAL shel_vi (3, 2) 
      REAL shel_000 
      COMPLEX shel_csf 
      COMPLEX shel_acsf 
      REAL shel_dsi 
      COMPLEX shel_tcsf 
      REAL factor 
!                                                                       
      REAL qval 
!                                                                       
!     If output type is shelx, calculate qval(000) for scaling          
!                                                                       
      IF (ityp.eq.HKLF4.or.ityp.eq.LIST5) then 
         DO i = 1, 2 
         shel_inc (i) = inc (i) 
         ENDDO 
         DO i = 1, 3 
         DO j = 1, 3 
         shel_eck (i, j) = eck (i, j) 
         ENDDO 
         DO j = 1, 2 
         shel_vi (i, j) = vi (i, j) 
         ENDDO 
         ENDDO 
         shel_tcsf = csf (1) 
         shel_acsf = acsf (1) 
         shel_dsi = dsi (1) 
         inc (1) = 1 
         inc (2) = 1 
         DO i = 1, 3 
         DO j = 1, 3 
         eck (i, j) = 0.0 
         ENDDO 
         DO j = 1, 2 
         vi (i, j) = 0.0 
         ENDDO 
         ENDDO 
         IF (ityp.eq.HKLF4) then 
            value = 1 
         ELSEIF (ityp.eq.LIST5) then 
            value = 2 
         ENDIF 
         CALL four_run 
         shel_csf = csf (1) 
         shel_000 = qval (1, value, 1, 1, laver) 
         qq = qval (1, value, 1, 1, laver) / cr_icc (1) / cr_icc (2)    &
         / cr_icc (3)                                                   
         IF (ityp.eq.HKLF4) then 
            factor = max (int (log (qq) / log (10.0) ) - 3, 0) 
         ELSEIF (ityp.eq.LIST5) then 
            factor = max (int (log (qq) / log (10.0) ) - 2, 0) 
         ENDIF 
         out_fac = 10** ( - factor) 
         csf (1) = shel_tcsf 
         acsf (1) = shel_acsf 
         dsi (1) = shel_dsi 
         DO i = 1, 2 
         inc (i) = shel_inc (i) 
         ENDDO 
         DO i = 1, 3 
         DO j = 1, 3 
         eck (i, j) = shel_eck (i, j) 
         ENDDO 
         DO j = 1, 2 
         vi (i, j) = shel_vi (i, j) 
         ENDDO 
         ENDDO 
      ENDIF 
!                                                                       
      extr_ima = 6 - out_extr_abs - out_extr_ord 
!                                                                       
      lread = .false. 
      CALL oeffne (iff, outfile, 'unknown', lread) 
      IF (ier_num.eq.0) then 
         IF (out_inc (1) .gt.1.and.out_inc (2) .gt.1) then 
            IF (ityp.eq.0) then 
               WRITE (iff, * ) out_inc 
               WRITE (iff, * ) out_eck (out_extr_abs, 1), out_eck (     &
               out_extr_abs, 2), out_eck (out_extr_ord, 1), out_eck (   &
               out_extr_ord, 3)                                         
               DO j = 1, out_inc (2) 
               WRITE (iff, 4) (qval ( (i - 1) * out_inc (2) + j, value, &
               i, j, laver), i = 1, out_inc (1) )                       
               WRITE (iff, 100) 
               ENDDO 
            ELSEIF (ityp.eq.HKLF4) then 
               DO j = 1, out_inc (2) 
               DO i = 1, out_inc (1) 
               DO k = 1, 3 
               h (k) = out_eck (k, 1) + out_vi (k, 1) * float (i - 1)   &
               + out_vi (k, 2) * float (j - 1)                          
               ENDDO 
               k = (i - 1) * out_inc (2) + j 
               qq = qval (k, value, i, j, laver) / cr_icc (1) / cr_icc (&
               2) / cr_icc (3) * out_fac                                
               sq = sqrt (qq) 
               WRITE (iff, 7) int (h (1) ), int (h (2) ), int (h (3) ), &
               qq, sq                                                   
               ENDDO 
               ENDDO 
            ELSEIF (ityp.eq.LIST5) then 
               shel_value = 3 
               DO j = 1, out_inc (2) 
               DO i = 1, out_inc (1) 
               DO k = 1, 3 
               h (k) = out_eck (k, 1) + out_vi (k, 1) * float (i - 1)   &
               + out_vi (k, 2) * float (j - 1)                          
               ENDDO 
               k = (i - 1) * out_inc (2) + j 
               shel_value = 2 
               qq = qval (k, value, i, j, laver) / cr_icc (1) / cr_icc (&
               2) / cr_icc (3) * out_fac                                
               shel_value = 3 
               sq = qval (k, shel_value, i, j, laver) 
               WRITE (iff, 8) int (h (1) ), int (h (2) ), int (h (3) ), &
               qq, qq, sq                                               
               ENDDO 
               ENDDO 
            ELSEIF (ityp.eq.LIST9) then 
               shel_value = 3 
               DO j = 1, out_inc (2) 
               DO i = 1, out_inc (1) 
               DO k = 1, 3 
               h (k) = out_eck (k, 1) + out_vi (k, 1) * float (i - 1)   &
               + out_vi (k, 2) * float (j - 1)                          
               ENDDO 
               k = (i - 1) * out_inc (2) + j 
               shel_value = 2 
               qq = qval (k, value, i, j, laver) / cr_icc (1) / cr_icc (&
               2) / cr_icc (3)                                          
               shel_value = 3 
               sq = qval (k, shel_value, i, j, laver) 
               WRITE (iff, 9) h (1), h (2), h (3), qq, qq, sq 
               ENDDO 
               ENDDO 
            ELSE 
               DO j = 1, out_inc (2) 
               DO i = 1, out_inc (1) 
               DO k = 1, 3 
               h (k) = out_eck (k, 1) + out_vi (k, 1) * float (i - 1)   &
               + out_vi (k, 2) * float (j - 1)                          
               ENDDO 
               k = (i - 1) * out_inc (2) + j 
               WRITE (iff, 5) h (out_extr_abs), h (out_extr_ord),       &
               qval (k, value, i, j, laver), h (extr_ima)               
               ENDDO 
               WRITE (iff, 100) 
               ENDDO 
            ENDIF 
         ELSEIF (out_inc (1) .eq.1) then 
            i = 1 
            DO j = 1, out_inc (2) 
            DO k = 1, 3 
            h (k) = out_eck (k, 1) + out_vi (k, 1) * float (i - 1)      &
            + out_vi (k, 2) * float (j - 1)                             
            ENDDO 
            k = (i - 1) * out_inc (2) + j 
            IF (ityp.eq.HKLF4) then 
               qq = qval (k, value, i, j, laver) / cr_icc (1) / cr_icc (&
               2) / cr_icc (3) * out_fac                                
               sq = sqrt (qq) 
               WRITE (iff, 7) int (h (1) ), int (h (2) ), int (h (3) ), &
               qq, sq                                                   
            ELSEIF (ityp.eq.LIST5) then 
               shel_value = 2 
               qq = qval (k, value, i, j, laver) / cr_icc (1) / cr_icc (&
               2) / cr_icc (3) * out_fac                                
               shel_value = 3 
               sq = qval (k, shel_value, i, j, laver) / cr_icc (1)      &
               / cr_icc (2) / cr_icc (3) * out_fac                      
               WRITE (iff, 8) int (h (1) ), int (h (2) ), int (h (3) ), &
               qq, qq, sq                                               
            ELSEIF (ityp.eq.LIST9) then 
               shel_value = 2 
               qq = qval (k, value, i, j, laver) / cr_icc (1) / cr_icc (&
               2) / cr_icc (3)                                          
               shel_value = 3 
               sq = qval (k, shel_value, i, j, laver) / cr_icc (1)      &
               / cr_icc (2) / cr_icc (3)                                
               WRITE (iff, 9) h (1), h (2), h (3), qq, qq, sq 
            ELSE 
               WRITE (iff, 6) h (out_extr_ord), qval (k, value, i, j,   &
               laver)                                                   
            ENDIF 
            ENDDO 
         ELSEIF (out_inc (2) .eq.1) then 
            j = 1 
            DO i = 1, out_inc (1) 
            DO k = 1, 3 
            h (k) = out_eck (k, 1) + out_vi (k, 1) * float (i - 1)      &
            + out_vi (k, 2) * float (j - 1)                             
            ENDDO 
            k = (i - 1) * out_inc (2) + j 
            IF (ityp.eq.HKLF4) then 
               qq = qval (k, value, i, j, laver) / cr_icc (1) / cr_icc (&
               2) / cr_icc (3) * out_fac                                
               sq = sqrt (qq) 
               WRITE (iff, 7) int (h (1) ), int (h (2) ), int (h (3) ), &
               qq, sq                                                   
            ELSEIF (ityp.eq.LIST5) then 
               shel_value = 2 
               qq = qval (k, value, i, j, laver) / cr_icc (1) / cr_icc (&
               2) / cr_icc (3) * out_fac                                
               shel_value = 3 
               sq = qval (k, shel_value, i, j, laver) / cr_icc (1)      &
               / cr_icc (2) / cr_icc (3) * out_fac                      
               WRITE (iff, 8) int (h (1) ), int (h (2) ), int (h (3) ), &
               qq, qq, sq                                               
            ELSEIF (ityp.eq.LIST9) then 
               shel_value = 2 
               qq = qval (k, value, i, j, laver) / cr_icc (1) / cr_icc (&
               2) / cr_icc (3)                                          
               shel_value = 3 
               sq = qval (k, shel_value, i, j, laver) / cr_icc (1)      &
               / cr_icc (2) / cr_icc (3)                                
               WRITE (iff, 9) h (1), h (2), h (3), qq, qq, sq 
            ELSE 
               WRITE (iff, 6) h (out_extr_abs), qval (k, value, i, j,   &
               laver)                                                   
            ENDIF 
            ENDDO 
         ENDIF 
      ENDIF 
!     if(ier_num.ne.0) then                                             
!       call errlist                                                    
!     endif                                                             
      IF (ityp.eq.HKLF4.or.ityp.eq.LIST5) then 
         WRITE (output_io, 1000) out_fac 
      ENDIF 
      CLOSE (iff) 
!                                                                       
    4 FORMAT (5(1x,e11.5)) 
    5 FORMAT (4(1x,e11.5)) 
    6 FORMAT (2(1x,e11.5)) 
    7 FORMAT (3i4,2f8.2) 
    8 FORMAT (3i4,2f10.2,f7.2) 
    9 FORMAT (3(f10.6,1x),2(e11.5,1x),f7.2) 
  100 FORMAT () 
 1000 FORMAT    (' Data have been scaled by ',g17.8e3) 
!100      format(/)                                                     
!                                                                       
      END SUBROUTINE do_output                      
!*****7*****************************************************************
      REAL function qval (i, value, ix, iy, laver) 
!-                                                                      
!     transforms the real and imaginary part of the Fourier transform   
!     into the desired output format                                    
!+                                                                      
      USE config_mod 
      USE diffuse_mod 
      USE output_mod 
      IMPLICIT none 
!                                                                       
      REAL DELTA 
      PARAMETER (DELTA = 0.000001) 
!                                                                       
      include'random.inc' 
!                                                                       
      INTEGER i, value, ix, iy 
      INTEGER k 
!                                                                       
      COMPLEX f 
      REAL h (3) 
!                                                                       
      REAL atan2d 
      REAL ran1 
      LOGICAL laver 
!                                                                       
!------ Get values of F or <F>                                          
!                                                                       
      IF (laver) then 
         f = acsf (i) 
      ELSE 
         f = csf (i) 
      ENDIF 
!                                                                       
!     Calculate intensity 'intensity'                                   
!                                                                       
!------ We have to store dsi() here, because if lots are                
!------ used, csf() will only contain the values for the                
!------ last lot !!                                                     
!                                                                       
      IF (value.eq.1) then 
         IF (laver) then 
            qval = real (f * conjg (f) ) 
         ELSE 
            qval = dsi (i) 
         ENDIF 
!                                                                       
!     Calculate amplitude 'amplitude'                                   
!                                                                       
      ELSEIF (value.eq.2) then 
         qval = sqrt (real (f * conjg (f) ) ) 
!                                                                       
!     Calculate phase 'phase'                                           
!                                                                       
      ELSEIF (value.eq.3) then 
         IF (f.eq. (0, 0) ) then 
            qval = 0.0 
         ELSE 
            qval = atan2d (aimag (f), real (f) ) 
         ENDIF 
!                                                                       
!     Calculate real part 'real'                                        
!                                                                       
      ELSEIF (value.eq.4) then 
         qval = real (f) 
!                                                                       
!     Calculate imaginary part 'imaginary'                              
!                                                                       
      ELSEIF (value.eq.5) then 
         qval = aimag (f) 
!                                                                       
!     Calculate phase 'phase', random, except for integer hkl           
!                                                                       
      ELSEIF (value.eq.6) then 
         DO k = 1, 3 
         h (k) = out_eck (k, 1) + out_vi (k, 1) * float (ix - 1)        &
         + out_vi (k, 2) * float (iy - 1)                               
         ENDDO 
         IF (abs (h (1) - nint (h (1) ) ) .lt.DELTA.and.abs (h (2)      &
         - nint (h (2) ) ) .lt.DELTA.and.abs (h (3) - nint (h (3) ) )   &
         .lt.DELTA) then                                                
            IF (f.eq. (0, 0) ) then 
               qval = 0.0 
            ELSE 
               qval = atan2d (aimag (f), real (f) ) 
            ENDIF 
         ELSE 
            qval = (ran1 (idum) - 0.5) * 360. 
         ENDIF 
      ENDIF 
!                                                                       
      END FUNCTION qval                             
!*****7*****************************************************************
      SUBROUTINE set_output (linverse) 
!-                                                                      
!     Sets the proper output values for either Fourier or               
!     inverse Fourier and Patterson                                     
!+                                                                      
      USE config_mod 
      USE diffuse_mod 
      USE output_mod 
      USE patters_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      LOGICAL linverse 
!                                                                       
      INTEGER i, j 
!                                                                       
      IF (linverse) then 
         out_extr_abs = rho_extr_abs 
         out_extr_ord = rho_extr_ord 
!                                                                       
         DO i = 1, 3 
         DO j = 1, 3 
         out_eck (i, j) = rho_eck (i, j) 
         ENDDO 
         ENDDO 
!                                                                       
         DO i = 1, 2 
         DO j = 1, 3 
         out_vi (j, i) = rho_vi (j, i) 
         ENDDO 
         out_inc (i) = rho_inc (i) 
         ENDDO 
      ELSE 
         out_extr_abs = extr_abs 
         out_extr_ord = extr_ord 
!                                                                       
         DO i = 1, 3 
         DO j = 1, 3 
         out_eck (i, j) = eck (i, j) 
         ENDDO 
         ENDDO 
!                                                                       
         DO i = 1, 2 
         DO j = 1, 3 
         out_vi (j, i) = vi (j, i) 
         ENDDO 
         out_inc (i) = inc (i) 
         ENDDO 
      ENDIF 
!                                                                       
      END SUBROUTINE set_output                     
