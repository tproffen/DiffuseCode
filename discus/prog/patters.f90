MODULE patters_menu
!
CONTAINS
!*****7*****************************************************************
!                                                                       
      SUBROUTINE patterson (inverse_type) 
!-                                                                      
!     patterson sets up the menu for inverse Fourier and Patterson      
!     calculations.                                                     
!+                                                                      
      USE discus_config_mod 
      USE discus_allocate_appl_mod
      USE crystal_mod 
      USE diffuse_mod 
      USE intens_mod 
      USE inverse_mod 
      USE metric_mod
      USE output_mod 
      USE patters_mod 
!
      USE doact_mod 
      USE errlist_mod 
      USE learn_mod 
      USE class_macro_internal
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER, INTENT(INOUT)  :: inverse_type
!
      INTEGER, PARAMETER :: MIN_PARA = 21  ! A command requires at leaset these no of parameters
      INTEGER maxw 
!                                                                       
      CHARACTER(LEN=1024), DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: cpara 
      INTEGER            , DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: lpara 
      REAL               , DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: werte
!
      CHARACTER(1) extr_achs (0:3) 
      CHARACTER(1) rho_extr_achs (0:3) 
      CHARACTER(5) befehl 
      CHARACTER(LEN=LEN(prompt)) :: orig_prompt
      CHARACTER(14) cvalue (0:6) 
      CHARACTER(15) cgraphik (0:5) 
      CHARACTER(1024) zeile
      CHARACTER(1024) line 
      LOGICAL lspace 
      LOGICAL ltwo_files 
      INTEGER e_io 
      INTEGER i, j, ianz, lp, length 
      INTEGER                  :: n_qxy    = 1 ! Number of data points in direct space
      INTEGER                  :: n_nscat  = 1 ! Number of different atom types this run
      INTEGER                  :: n_natom  = 1 ! Number of atoms this run
      INTEGER indxg, lbef 
      REAL divis (2) 
      REAL rho_divis (2) 
      REAL u (3), v (3), w (3), dvi1, dvi2, dvi3, dvi4, dvi5 
!                                                                       
      INTEGER len_str 
      LOGICAL str_comp 
!     REAL do_blen, do_bang 
!                                                                       
      DATA extr_achs / ' ', 'h', 'k', 'l' / 
      DATA rho_extr_achs / ' ', 'x', 'y', 'z' / 
!                                                                       
      DATA cgraphik / 'Standard', 'Postscript', 'Pseudo Grey Map',      &
      'Gnuplot', 'SHELXL', 'HKL F4' /                                   
!                                                                       
      DATA cvalue / 'undefined     ', 'Intensity     ', 'Amplitude=Fobs'&
     &, 'Phase angle   ', 'Real Part     ', 'Imaginary Part', 'Fcalc    &
     &     ' /                                                          
!
      maxw     = MAX(MIN_PARA,MAXSCAT+1)
!     n_qxy    = 1
!     n_nscat  = 1
!     n_natom  = 1
!                                                                       
      IF (cr_natoms.gt.0.and.as_natoms.gt.0) then 
         patt_scale = float (cr_ncatoms) / float (cr_natoms) 
      ELSEIF (patt_scale.eq.0) then 
         patt_scale = 1.0 
      ENDIF 
!
      orig_prompt = prompt
      IF (inverse_type.eq.INV_INV) then 
         prompt = prompt (1:len_str (prompt) ) //'/inverse' 
      ELSEIF (inverse_type.eq.INV_DIFF) then 
         prompt = prompt (1:len_str (prompt) ) //'/diff-four' 
      ELSEIF (inverse_type.eq.INV_PATT) then 
         prompt = prompt (1:len_str (prompt) ) //'/patterson' 
      ENDIF 
!                                                                       
   10 CONTINUE 
!                                                                       
      CALL no_error 
!                                                                       
      divis (1) = float (max (1, inc (1) - 1) ) 
      divis (2) = float (max (1, inc (2) - 1) ) 
      rho_divis (1) = float (max (1, rho_inc (1) - 1) ) 
      rho_divis (2) = float (max (1, rho_inc (2) - 1) ) 
!                                                                       
      CALL get_cmd (line, length, befehl, lbef, zeile, lp, prompt) 
      IF (ier_num.eq.0) then 
         IF (line (1:1)  == ' '.or.line (1:1)  == '#' .or.   & 
             line == char(13) .or. line(1:1) == '!'  ) GOTO 10
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
!     Define the ascissa 'absc'                                         
!                                                                       
            ELSEIF (str_comp (befehl, 'absc', 1, lbef, 4) ) then 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ianz.eq.1) then 
                  IF (cpara (1) .eq.'h') then 
                     extr_abs = 1 
                  ELSEIF (cpara (1) .eq.'k') then 
                     extr_abs = 2 
                  ELSEIF (cpara (1) .eq.'l') then 
                     extr_abs = 3 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
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
!      Evaluate an expression, just for interactive check 'eval'        
!                                                                       
            ELSEIF (str_comp (befehl, 'eval', 2, lbef, 4) ) then 
               CALL do_eval (zeile, lp) 
!                                                                       
!     Terminate invers Fourier / Patterson 'exit'                       
!                                                                       
            ELSEIF (str_comp (befehl, 'exit', 3, lbef, 4) ) then 
               GOTO 9999 
!                                                                       
!     define the input file names 'file'                                
!                                                                       
            ELSEIF (str_comp (befehl, 'file', 1, lbef, 4) ) then 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ier_num.eq.0) then 
                  IF (ianz.ge.2) then 
      IF (str_comp (cpara (1) , 'a    ', 1, lpara (2) , 5) ) then 
                        CALL do_build_name (ianz, cpara, lpara, werte,  &
                        maxw, 2)                                        
                        IF (ier_num.eq.0) then 
                           i = 1 
                           rho_file (i) = cpara (2) (1:lpara(2))
                        ENDIF 
      ELSEIF (str_comp (cpara (1) , 'b    ', 1, lpara (2) , 5) ) then 
                        CALL do_build_name (ianz, cpara, lpara, werte,  &
                        maxw, 2)                                        
                        IF (ier_num.eq.0) then 
                           i = 2 
                           rho_file (i) = cpara (2) (1:lpara(2))
                        ENDIF 
                     ELSE 
                        ier_num = - 41 
                        ier_typ = ER_APPL 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ENDIF 
!                                                                       
!     Determine format for input 'format'                               
!                                                                       
            ELSEIF (str_comp (befehl, 'form', 1, lbef, 4) ) then 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ier_num.eq.0) then 
                  IF (ianz.eq.1.or.ianz.eq.2) then 
!                                                                       
!     --Switch input type to GNUPLOT 'gnup'                             
!                                                                       
                     IF (str_comp (cpara (1) , 'gnup', 1, lpara (1) , 4)&
                     ) then                                             
                        ftyp = 3 
!                                                                       
!     Switch input type to standard  'stan'                             
!                                                                       
                     ELSEIF (str_comp (cpara (1) , 'stan', 2, lpara (1) &
                     , 4) ) then                                        
                        ftyp = 0 
!                                                                       
!     Switch input type to SHELXL List Type 5 'shelxl'                  
!                                                                       
                     ELSEIF (str_comp (cpara (1) , 'shelx', 2, lpara (1)&
                     , 5) ) then                                        
                        ftyp = 4 
                        rho_type (1) = AMPLITUDE 
                        IF (ianz.eq.2) then 
                           IF (str_comp (cpara (2) , 'fobs', 2, lpara ( &
                           2) , 4) ) then                               
                              rho_type (1) = AMPLITUDE 
                           ELSEIF (str_comp (cpara (2) , 'fcalc', 2,    &
                           lpara (2) , 5) ) then                        
                              rho_type (1) = FCALC 
                           ELSE 
                              ier_num = - 52 
                              ier_typ = ER_APPL 
                           ENDIF 
                        ENDIF 
!                                                                       
!     Switch input type to SHELXS 'hklf4'                               
!                                                                       
                     ELSEIF (str_comp (cpara (1) , 'hklf4', 2, lpara (1)&
                     , 5) ) then                                        
                        ftyp = 5 
                        rho_type (1) = INTENSITY 
                        inverse_type = INV_PATT 
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
!     help 'help' , '?'                                                 
!                                                                       
      ELSEIF (str_comp (befehl, 'help', 2, lbef, 4) .or.str_comp (befehl&
     &, '?   ', 1, lbef, 4) ) then                                      
               IF (str_comp (zeile, 'errors', 2, lp, 6) ) then 
                  lp = lp + 7 
                  CALL do_hel ('discus '//zeile, lp) 
               ELSE 
                  IF (inverse_type.eq.INV_INV) then 
                     lp = lp + 15 
                     CALL do_hel ('discus inverse '//zeile, lp) 
                  ELSEIF (inverse_type.eq.INV_DIFF) then 
                     lp = lp + 17 
                     CALL do_hel ('discus diff-four '//zeile, lp) 
                  ELSEIF (inverse_type.eq.INV_PATT) then 
                     lp = lp + 17 
                     CALL do_hel ('discus patterson '//zeile, lp) 
                  ENDIF 
               ENDIF 
!                                                                       
!     define the whole layer 'laye'                                     
!                                                                       
            ELSEIF (str_comp (befehl, 'laye', 2, lbef, 4) ) then 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ier_num.eq.0) then 
                  IF (ianz.eq.11) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        DO j = 1, 3 
                        DO i = 1, 3 
                        eck (i, j) = werte ( (j - 1) * 3 + i) 
                        ENDDO 
                        ENDDO 
                        inc (1) = nint (werte (10) ) 
                        inc (2) = nint (werte (11) ) 
                        divis (1) = float (max (1, inc (1) - 1) ) 
                        divis (2) = float (max (1, inc (2) - 1) ) 
                        DO i = 1, 3 
                        vi (i, 1) = (eck (i, 2) - eck (i, 1) ) / divis (&
                        1)                                              
                        vi (i, 2) = (eck (i, 3) - eck (i, 1) ) / divis (&
                        2)                                              
                        ENDDO 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ENDIF 
!                                                                       
!     define the lower left corner 'll'                                 
!                                                                       
      ELSEIF (str_comp (befehl, 'll  ', 2, lbef, 4) ) then 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ier_num.eq.0) then 
                  IF (ianz.eq.3) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        j = 1 
                        DO i = 1, 3 
                        eck (i, j) = werte (i) 
                        ENDDO 
                        DO i = 1, 3 
                        vi (i, 1) = (eck (i, 2) - eck (i, 1) ) / divis (&
                        1)                                              
                        vi (i, 2) = (eck (i, 3) - eck (i, 1) ) / divis (&
                        2)                                              
                        ENDDO 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ENDIF 
!                                                                       
!     define the lower right corner 'lr'                                
!                                                                       
      ELSEIF (str_comp (befehl, 'lr  ', 2, lbef, 4) ) then 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ier_num.eq.0) then 
                  IF (ianz.eq.3) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        j = 2 
                        DO i = 1, 3 
                        eck (i, j) = werte (i) 
                        ENDDO 
                        DO i = 1, 3 
                        vi (i, 1) = (eck (i, 2) - eck (i, 1) ) / divis (&
                        1)                                              
                        vi (i, 2) = (eck (i, 3) - eck (i, 1) ) / divis (&
                        2)                                              
                        ENDDO 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ENDIF 
!                                                                       
!     define the number of points along the abscissa 'na'               
!                                                                       
            ELSEIF (str_comp (befehl, 'nabs', 2, lbef, 4) ) then 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ier_num.eq.0) then 
                  IF (ianz.eq.1) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        IF (werte (1) .gt.0) then 
                           inc (1) = nint (werte (1) ) 
                           divis (1) = float (max (1, inc (1) - 1) ) 
                           DO i = 1, 3 
                           vi (i, 1) = (eck (i, 2) - eck (i, 1) )       &
                           / divis (1)                                  
                           vi (i, 2) = (eck (i, 3) - eck (i, 1) )       &
                           / divis (2)                                  
                           ENDDO 
                        ELSE 
                           ier_num = - 12 
                           ier_typ = ER_APPL 
                        ENDIF 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ENDIF 
!     define the number of points along the ordinate 'no'               
!                                                                       
            ELSEIF (str_comp (befehl, 'nord', 2, lbef, 4) ) then 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ier_num.eq.0) then 
                  IF (ianz.eq.1) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        IF (werte (1) .gt.0) then 
                           inc (2) = nint (werte (1) ) 
                           divis (2) = float (max (1, inc (2) - 1) ) 
                           DO i = 1, 3 
                           vi (i, 1) = (eck (i, 2) - eck (i, 1) )       &
                           / divis (1)                                  
                           vi (i, 2) = (eck (i, 3) - eck (i, 1) )       &
                           / divis (2)                                  
                           ENDDO 
                        ELSE 
                           ier_num = - 12 
                           ier_typ = ER_APPL 
                        ENDIF 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ENDIF 
!                                                                       
!     define the ordinate  'ordi'                                       
!                                                                       
            ELSEIF (str_comp (befehl, 'ordi', 2, lbef, 4) ) then 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ianz.eq.1) then 
                  IF (cpara (1) .eq.'h') then 
                     extr_ord = 1 
                  ELSEIF (cpara (1) .eq.'k') then 
                     extr_ord = 2 
                  ELSEIF (cpara (1) .eq.'l') then 
                     extr_ord = 3 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ENDIF 
!                                                                       
!     Define the real space ascissa 'rhoabsc'                           
!                                                                       
            ELSEIF (str_comp (befehl, 'rhoab', 4, lbef, 5) ) then 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ianz.eq.1) then 
                  IF (cpara (1) .eq.'x') then 
                     rho_extr_abs = 1 
                  ELSEIF (cpara (1) .eq.'y') then 
                     rho_extr_abs = 2 
                  ELSEIF (cpara (1) .eq.'z') then 
                     rho_extr_abs = 3 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
!                                                                       
!     define the whole real space layer 'rholayer'                      
!                                                                       
            ELSEIF (str_comp (befehl, 'rhola', 5, lbef, 5) ) then 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ier_num.eq.0) then 
                  IF (ianz.eq.11) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        DO j = 1, 3 
                        DO i = 1, 3 
                        rho_eck (i, j) = werte ( (j - 1) * 3 + i) 
                        ENDDO 
                        ENDDO 
                        rho_inc (1) = nint (werte (10) ) 
                        rho_inc (2) = nint (werte (11) ) 
                        rho_divis (1) = float (max (1, rho_inc (1)      &
                        - 1) )                                          
                        rho_divis (2) = float (max (1, rho_inc (2)      &
                        - 1) )                                          
                        DO i = 1, 3 
                        rho_vi (i, 1) = (rho_eck (i, 2) - rho_eck (i, 1)&
                        ) / rho_divis (1)                               
                        rho_vi (i, 2) = (rho_eck (i, 3) - rho_eck (i, 1)&
                        ) / rho_divis (2)                               
                        ENDDO 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ENDIF 
!                                                                       
!     define the real space lower left corner 'rholl'                   
!                                                                       
            ELSEIF (str_comp (befehl, 'rholl', 5, lbef, 5) ) then 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ier_num.eq.0) then 
                  IF (ianz.eq.3) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        j = 1 
                        DO i = 1, 3 
                        rho_eck (i, j) = werte (i) 
                        ENDDO 
                        DO i = 1, 3 
                        rho_vi (i, 1) = (rho_eck (i, 2) - rho_eck (i, 1)&
                        ) / rho_divis (1)                               
                        rho_vi (i, 2) = (rho_eck (i, 3) - rho_eck (i, 1)&
                        ) / rho_divis (2)                               
                        ENDDO 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ENDIF 
!                                                                       
!     define the real space lower right corner 'rholr'                  
!                                                                       
            ELSEIF (str_comp (befehl, 'rholr', 5, lbef, 5) ) then 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ier_num.eq.0) then 
                  IF (ianz.eq.3) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        j = 2 
                        DO i = 1, 3 
                        rho_eck (i, j) = werte (i) 
                        ENDDO 
                        DO i = 1, 3 
                        rho_vi (i, 1) = (rho_eck (i, 2) - rho_eck (i, 1)&
                        ) / rho_divis (1)                               
                        rho_vi (i, 2) = (rho_eck (i, 3) - rho_eck (i, 1)&
                        ) / rho_divis (2)                               
                        ENDDO 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ENDIF 
!                                                                       
!------ define the real space number of points along the                
!       abscissa 'rhonabs'                                              
!                                                                       
            ELSEIF (str_comp (befehl, 'rhona', 5, lbef, 5) ) then 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ier_num.eq.0) then 
                  IF (ianz.eq.1) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        IF (werte (1) .gt.0) then 
                           rho_inc (1) = nint (werte (1) ) 
                           rho_divis (1) = float (max (1, rho_inc (1)   &
                           - 1) )                                       
                           DO i = 1, 3 
                           rho_vi (i, 1) = (rho_eck (i, 2) - rho_eck (i,&
                           1) ) / rho_divis (1)                         
                           rho_vi (i, 2) = (rho_eck (i, 3) - rho_eck (i,&
                           1) ) / rho_divis (2)                         
                           ENDDO 
                        ELSE 
                           ier_num = - 12 
                           ier_typ = ER_APPL 
                        ENDIF 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ENDIF 
!                                                                       
!------ define the real space number of points along the                
!       ordinate 'rhonord'                                              
!                                                                       
            ELSEIF (str_comp (befehl, 'rhono', 5, lbef, 5) ) then 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ier_num.eq.0) then 
                  IF (ianz.eq.1) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        IF (werte (1) .gt.0) then 
                           rho_inc (2) = nint (werte (1) ) 
                           rho_divis (2) = float (max (1, rho_inc (2)   &
                           - 1) )                                       
                           DO i = 1, 3 
                           rho_vi (i, 1) = (rho_eck (i, 2) - rho_eck (i,&
                           1) ) / rho_divis (1)                         
                           rho_vi (i, 2) = (rho_eck (i, 3) - rho_eck (i,&
                           1) ) / rho_divis (2)                         
                           ENDDO 
                        ELSE 
                           ier_num = - 12 
                           ier_typ = ER_APPL 
                        ENDIF 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ENDIF 
!                                                                       
!     define the real space ordinate  'rhoordi'                         
!                                                                       
            ELSEIF (str_comp (befehl, 'rhoor', 5, lbef, 5) ) then 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ianz.eq.1) then 
                  IF (cpara (1) .eq.'x') then 
                     rho_extr_ord = 1 
                  ELSEIF (cpara (1) .eq.'y') then 
                     rho_extr_ord = 2 
                  ELSEIF (cpara (1) .eq.'z') then 
                     rho_extr_ord = 3 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ENDIF 
!                                                                       
!     define the real space  upper left corner 'rhoul'                  
!                                                                       
            ELSEIF (str_comp (befehl, 'rhoul', 4, lbef, 5) ) then 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ier_num.eq.0) then 
                  IF (ianz.eq.3) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        j = 3 
                        DO i = 1, 3 
                        rho_eck (i, j) = werte (i) 
                        ENDDO 
                        DO i = 1, 3 
                        rho_vi (i, 1) = (rho_eck (i, 2) - rho_eck (i, 1)&
                        ) / rho_divis (1)                               
                        rho_vi (i, 2) = (rho_eck (i, 3) - rho_eck (i, 1)&
                        ) / rho_divis (2)                               
                        ENDDO 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ENDIF 
!                                                                       
!     start the inverse Fourier transform 'run'                         
!                                                                       
            ELSEIF (str_comp (befehl, 'run ', 1, lbef, 4) ) then 
!                                                                       
               IF (cr_v.le.0.0) then 
                  ier_num = - 35 
                  ier_typ = ER_APPL 
                  ier_msg (1) = 'A proper unit cell must be defined' 
                  ier_msg (2) = 'for this command to operate       ' 
               ELSE 
                 IF (inc(1)     * inc(2)     .gt. MAXQXY  .OR.          &
                     rho_inc(1) * rho_inc(2) .gt. MAXQXY  .OR.          &
                     cr_nscat>DIF_MAXSCAT              ) THEN
                    n_qxy   = MAX(1, n_qxy,inc(1) * inc(2),rho_inc(1)*rho_inc(2),MAXQXY)
                    n_nscat = MAX(1, n_nscat,cr_nscat,DIF_MAXSCAT)
                    n_natom = MAX(1, n_natom,cr_natoms, NMAX, DIF_MAXAT)
                    call alloc_diffuse (n_qxy,  n_nscat, n_natom )
                    IF (ier_num.ne.0) THEN
                      RETURN
                    ENDIF
                  ENDIF
                  IF (rho_inc (1) * rho_inc (2) .le.MAXQXY) then 
                     IF (rho_type (2) .eq.INTENSITY.or.rho_type (2)     &
                     .eq.AMPLITUDE.or.rho_type (2) .eq.REAL_PART) then  
                        line = ' ' 
                        line = rho_file (2) 
                        i = rho_type (2) 
                        rho_file (2) = rho_file (1) 
                        rho_type (2) = rho_type (1) 
                        rho_file (1) = trim(line)
                        rho_type (1) = i 
                     ENDIF 
                     IF (ftyp.eq.4) then 
                        ltwo_files = .false. 
                        IF (rho_type (1) .eq.AMPLITUDE.or.rho_type (1)  &
                        .eq.FCALC) then                                 
                           CALL no_error 
                        ELSE 
                           ier_num = - 51 
                           ier_typ = ER_APPL 
                        ENDIF 
                     ELSE 
                        ltwo_files = .true. 
                        IF (inverse_type.eq.INV_INV) then 
                           IF (rho_type (1) .eq.INTENSITY.and.rho_type (&
                           2) .eq.PHASE_ANG.or.rho_type (1)             &
                           .eq.AMPLITUDE.and.rho_type (2)               &
                           .eq.PHASE_ANG.or.rho_type (1)                &
                           .eq.REAL_PART.and.rho_type (2) .eq.IMAG_PART)&
                           then                                         
                              CALL no_error 
                           ELSE 
                              ier_num = - 36 
                              ier_typ = ER_APPL 
                           ENDIF 
                        ELSEIF (inverse_type.eq.INV_PATT) then 
                           IF (rho_type (1) .eq.INTENSITY.or.rho_type ( &
                           1) .eq.AMPLITUDE.or.rho_type (1)             &
                           .eq.REAL_PART.and.rho_type (2) .eq.IMAG_PART)&
                           then                                         
                              CALL no_error 
                           ELSE 
                              ier_num = - 38 
                              ier_typ = ER_APPL 
                           ENDIF 
                           IF (rho_type (1) .eq.INTENSITY) then 
                              ltwo_files = .false. 
      IF (patt_mode.eq.PATT_SHARP.or.patt_mode.eq.PATT_SUPER) then 
                                 IF (ftyp.eq.5) then 
                                    cpara (1) = 'screen' 
                                    lpara (1) = 6 
                                    cpara (2) = rho_file (1) 
                                    lpara (2) = len_str (cpara (2) ) 
                                    CALL e_create (e_io, cpara (2),     &
                                    lpara (2), .false.)                 
                                 ELSE 
                                    ier_num = - 38 
                                    ier_typ = ER_APPL 
                                 ENDIF 
                              ENDIF 
                           ENDIF 
                        ELSEIF (inverse_type.eq.INV_DIFF) then 
                           ltwo_files = .false. 
                        ENDIF 
                     ENDIF 
                     IF (ier_num.eq.0) then 
                        CALL do_patters (inverse_type, ltwo_files,      &
                        cr_vr, cr_acentric)                             
                        four_was_run = .true.
                     ENDIF 
                  ELSE 
                     ier_num = - 71 
                     ier_typ = ER_APPL 
                  ENDIF 
               ENDIF 
!                                                                       
!     Run statistics on the normalized structure factor 'stat'          
!                                                                       
            ELSEIF (str_comp (befehl, 'stat', 1, lbef, 4) ) then 
               IF (cr_v.le.0.0) then 
                  ier_num = - 35 
                  ier_typ = ER_APPL 
                  ier_msg (1) = 'A proper unit cell must be defined' 
      ier_msg (2)  = 'for this command to operate       ' 
               ELSE 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (ianz.ge.1) then 
                        IF (str_comp (cpara (1) , 'screen', 2, lpara (1)&
                        , 6) ) then                                     
                           e_io = 6 
                        ELSEIF (str_comp (cpara (1) , 'file', 2, lpara (&
                        1) , 4) ) then                                  
                           e_io = 38 
                           CALL do_build_name (ianz, cpara, lpara,      &
                           werte, maxw, 2)                              
                        ELSE 
                           ier_num = - 6 
                           ier_typ = ER_COMM 
                        ENDIF 
                        IF (ier_num.eq.0) then 
                           CALL e_create (e_io, cpara (2), lpara (2),   &
                           .true.)                                      
                        ENDIF 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ENDIF 
               ENDIF 
!                                                                       
!     Define the type of the input file 'type'                          
!                                                                       
            ELSEIF (str_comp (befehl, 'type', 1, lbef, 4) ) then 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ier_num.eq.0) then 
                  IF (ianz.eq.2) then 
      IF (str_comp (cpara (1) , 'a    ', 1, lpara (1) , 5) ) then 
                        i = 1 
      ELSEIF (str_comp (cpara (1) , 'b    ', 1, lpara (1) , 5) ) then 
                        i = 2 
                     ELSE 
                        ier_num = - 41 
                        ier_typ = ER_APPL 
                     ENDIF 
                     IF (ier_num.eq.0) then 
                        IF (str_comp (cpara (2) , 'inten', 2, lpara (2) &
                        , 5) ) then                                     
                           rho_type (i) = INTENSITY 
                        ELSEIF (str_comp (cpara (2) , 'ampli', 1, lpara &
                        (2) , 5) ) then                                 
                           rho_type (i) = AMPLITUDE 
                        ELSEIF (str_comp (cpara (2) , 'phase', 1, lpara &
                        (2) , 5) ) then                                 
                           rho_type (i) = PHASE_ANG 
                        ELSEIF (str_comp (cpara (2) , 'real ', 1, lpara &
                        (2) , 5) ) then                                 
                           rho_type (i) = REAL_PART 
                        ELSEIF (str_comp (cpara (2) , 'imagi', 2, lpara &
                        (2) , 5) ) then                                 
                           rho_type (i) = IMAG_PART 
                        ELSEIF (str_comp (cpara (2) , 'fobs', 2, lpara (&
                        2) , 4) ) then                                  
                           rho_type (i) = AMPLITUDE 
                        ELSEIF (str_comp (cpara (2) , 'fcalc', 2, lpara &
                        (2) , 5) ) then                                 
                           rho_type (i) = FCALC 
                        ELSE 
                           rho_type (i) = 0 
                           ier_num = - 42 
                           ier_typ = ER_APPL 
                        ENDIF 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ENDIF 
!                                                                       
!     define the scale factor                                           
!                                                                       
            ELSEIF (str_comp (befehl, 'scale', 2, lbef, 5) ) then 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ier_num.eq.0) then 
                  IF (ianz.eq.1) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        patt_scale = werte (1) 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ENDIF 
!                                                                       
!     set various parameters 'set'                                      
!                                                                       
            ELSEIF (str_comp (befehl, 'set', 2, lbef, 3) ) then 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ier_num.eq.0) then 
                  IF (ianz.eq.2) then 
!                                                                       
!     ------set mode for accumulation 'accu'                            
!           either initialize or add to previous patterson.             
!                                                                       
                     IF (str_comp (cpara (1) , 'accu', 1, lpara (1) , 4)&
                     ) then                                             
                        IF (str_comp (cpara (2) , 'init', 1, lpara (2) ,&
                        4) ) then                                       
                           patt_accu = PATT_INIT 
                        ELSEIF (str_comp (cpara (2) , 'add', 1, lpara ( &
                        2) , 3) ) then                                  
                           patt_accu = PATT_ADD 
                        ELSE 
                           ier_num = - 6 
                           ier_typ = ER_COMM 
                        ENDIF 
                     ELSEIF (str_comp (cpara (1) , 'excl', 1, lpara (1) &
                     , 4) ) then                                        
                        IF (str_comp (cpara (2) , 'none', 1, lpara (2) ,&
                        4) ) then                                       
                           patt_excl9999 = .false. 
                        ELSE 
                           cpara (1) = '0' 
                           lpara (1) = 1 
                           CALL ber_params (ianz, cpara, lpara, werte,  &
                           maxw)                                        
                           IF (ier_num.eq.0) then 
                              patt_excl9999 = .true. 
                              patt_excl_val = werte (2) 
                           ENDIF 
                        ENDIF 
                     ELSEIF (str_comp (cpara (1) , 'mode', 1, lpara (1) &
                     , 4) ) then                                        
                        IF (str_comp (cpara (2) , 'normal', 1, lpara (2)&
                        , 4) ) then                                     
                           patt_mode = PATT_NORMAL 
                        ELSEIF (str_comp (cpara (2) , 'sharp', 1, lpara &
                        (2) , 3) ) then                                 
                           patt_mode = PATT_SHARP 
                        ELSEIF (str_comp (cpara (2) , 'super', 1, lpara &
                        (2) , 3) ) then                                 
                           patt_mode = PATT_SUPER 
                        ELSE 
                           ier_num = - 6 
                           ier_typ = ER_COMM 
                        ENDIF 
                     ELSEIF (str_comp (cpara (1) , 'origin', 1, lpara ( &
                     1) , 4) ) then                                     
                        IF (str_comp (cpara (2) , 'normal', 1, lpara (2)&
                        , 4) ) then                                     
                           patt_origin = PATT_NORMAL 
                        ELSEIF (str_comp (cpara (2) , 'subtract', 1,    &
                        lpara (2) , 3) ) then                           
                           patt_origin = PATT_SUBTRACT 
                        ELSE 
                           ier_num = - 6 
                           ier_typ = ER_COMM 
                        ENDIF 
                     ELSEIF (str_comp (cpara (1) , 'rsym', 1, lpara (1) &
                     , 4) ) then                                        
                        IF (str_comp (cpara (2) , 'igno', 1, lpara (2) ,&
                        4) ) then                                       
                           patt_rsym = .false. 
                        ELSEIF (str_comp (cpara (2) , 'appl', 1, lpara (&
                        2) , 4) ) then                                  
                           patt_rsym = .true. 
                        ENDIF 
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
!     define the sign for inverse or direct Fourier 'sign'              
!                                                                       
            ELSEIF (str_comp (befehl, 'sign', 2, lbef, 4) ) then 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ier_num.eq.0) then 
                  IF (str_comp (cpara (1) , 'inverse', 1, lpara (1) , 7)&
                  ) then                                                
                     patt_sign = - 1 
                  ELSEIF (str_comp (cpara (1) , 'fourier', 1, lpara (1) &
                  , 7) ) then                                           
                     patt_sign = 1 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ENDIF 
!                                                                       
!     Show the current settings for the Fourier 'show'                  
!                                                                       
            ELSEIF (str_comp (befehl, 'show', 2, lbef, 4) ) then 
               IF (inverse_type.eq.INV_INV) then 
      WRITE (output_io,  * ) ' technique           : inverse Fourier' 
               ELSEIF (inverse_type.eq.INV_DIFF) then 
!23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
      WRITE (output_io,  * ) ' technique           : Difference Fourier &
     &'                                                                 
               ELSEIF (inverse_type.eq.INV_PATT) then 
                  IF (patt_mode.eq.PATT_NORMAL) then 
      WRITE (output_io,  * ) ' technique           : Patterson ' 
                  ELSEIF (patt_mode.eq.PATT_SHARP) then 
      WRITE (output_io,  * ) ' technique           : ', 'sharpened Patte&
     &rson '                                                            
                  ELSEIF (patt_mode.eq.PATT_SUPER) then 
      WRITE (output_io,  * ) ' technique           : ', 'super sharpened&
     & Patterson '                                                      
                  ENDIF 
                  IF (patt_origin.eq.PATT_SUBTRACT) then 
      WRITE (output_io,  * ) ' origin peak         : removed' 
                  ENDIF 
               ENDIF 
               IF (patt_accu.eq.PATT_INIT) then 
      WRITE (output_io,  * ) ' Accumulation mode   : Initialize' 
               ELSEIF (patt_accu.eq.PATT_ADD) then 
      WRITE (output_io,  * ) ' Accumulation mode   : Add to previous' 
               ENDIF 
               IF (patt_excl9999) then 
                  WRITE (output_io, 30605) patt_excl_val 
30605 FORMAT     ('  Exclusions          : Data points equal to ',      &
     &           g13.6e2,' will be excluded')                           
               ELSE 
      WRITE (output_io,  * ) ' Exclusions          : All data points wil&
     &l ', 'be treated as data'                                         
               ENDIF 
               WRITE (output_io, * ) 
      WRITE (output_io,  * ) ' Geometry            : reciprocal space' 
               IF (ftyp.eq.4.or.ftyp.eq.5) then 
                  IF (patt_rsym) then 
      WRITE (output_io,  * ) ' Input file format   : ', cgraphik (ftyp) &
     &, ' reciprocal space symmetry is applied'                         
                  ELSE 
      WRITE (output_io,  * ) ' Input file format   : ', cgraphik (ftyp) &
     &, ' reciprocal space symmetry is ignored'                         
                  ENDIF 
               ELSE 
      WRITE (output_io,  * ) ' Input file format   : ', cgraphik (ftyp) 
               ENDIF 
      WRITE (output_io,  * ) ' Input file A        : ', rho_file (1) 
      WRITE (output_io,  * ) ' File type  A        : ', cvalue (rho_type&
     & (1) )                                                            
      WRITE (output_io,  * ) ' Input file B        : ', rho_file (2) 
      WRITE (output_io,  * ) ' File type  B        : ', cvalue (rho_type&
     & (2) )                                                            
      WRITE (output_io,  * ) ' Scale factor        : ', patt_scale 
      WRITE (output_io,  * ) ' Sign in exponent    : ', patt_sign 
!                                                                       
!     --Output of reciprocal space geometry only for non-SHELXL files   
!                                                                       
               IF (ftyp.ne.4.and.ftyp.ne.5) then 
                  WRITE (output_io, 3060) ( (eck (i, j), i = 1, 3),     &
                  j = 1, 3)                                             
                  WRITE (output_io, 3061) inc 
                  DO i = 1, 3 
                  u (i) = vi (i, 1) 
                  v (i) = 0.0 
                  w (i) = vi (i, 2) 
                  ENDDO 
                  lspace = .false. 
                  dvi1 = do_blen (lspace, u, v) 
                  dvi2 = do_blen (lspace, w, v) 
                  IF (inc (1) .gt.1.and.inc (2) .gt.1) then 
                     dvi3 = do_bang (lspace, u, v, w) 
                  ELSE 
                     dvi3 = 0.0 
                  ENDIF 
                  IF (dvi3.gt.0) then 
                     dvi4 = dvi2 / dvi1 
                     IF (abs (u (extr_abs) ) .gt.0.0.and.abs (w (       &
                     extr_ord) ) .gt.0.0) then                          
                        dvi5 = (dvi2 / w (extr_ord) ) / (dvi1 / u (     &
                        extr_abs) )                                     
                     ELSE 
                        ier_num = - 4 
                        ier_typ = ER_FOUR 
                     ENDIF 
                  ELSE 
                     dvi4 = 0.0 
                     dvi5 = 0.0 
                  ENDIF 
                  WRITE (output_io, 3062) (vi (i, 1), i = 1, 3),        &
                  dvi1, (vi (i, 2), i = 1, 3), dvi2, dvi3, dvi4, dvi5   
                  WRITE (output_io, 3063) extr_achs (extr_abs) 
                  WRITE (output_io, 3064) extr_achs (extr_ord) 
                  IF (ier_num.ne.0) then 
                     CALL errlist 
                     CALL no_error 
                  ENDIF 
               ENDIF 
!                                                                       
               WRITE (output_io, * ) 
      WRITE (output_io,  * ) ' Geometry   : real space' 
               WRITE (output_io, * ) 
               WRITE (output_io, 3260) ( (rho_eck (i, j), i = 1, 3),    &
               j = 1, 3)                                                
               WRITE (output_io, 3061) rho_inc 
               DO i = 1, 3 
               u (i) = rho_vi (i, 1) 
               v (i) = 0.0 
               w (i) = rho_vi (i, 2) 
               ENDDO 
               lspace = .true. 
               dvi1 = do_blen (lspace, u, v) 
               dvi2 = do_blen (lspace, w, v) 
               IF (rho_inc (1) .gt.1.and.rho_inc (2) .gt.1) then 
                  dvi3 = do_bang (lspace, u, v, w) 
               ELSE 
                  dvi3 = 0.0 
               ENDIF 
               IF (dvi3.gt.0) then 
                  dvi4 = dvi2 / dvi1 
                  IF (abs (u (rho_extr_abs) ) .gt.0.0.and.abs (w (      &
                  rho_extr_ord) ) .gt.0.0) then                         
                     dvi5 = (dvi2 / w (rho_extr_ord) ) / (dvi1 / u (    &
                     rho_extr_abs) )                                    
                  ELSE 
                     ier_num = - 4 
                     ier_typ = ER_FOUR 
                  ENDIF 
               ELSE 
                  dvi4 = 0.0 
                  dvi5 = 0.0 
               ENDIF 
               WRITE (output_io, 3262) (rho_vi (i, 1), i = 1, 3),       &
               dvi1, (rho_vi (i, 2), i = 1, 3), dvi2, dvi3, dvi4, dvi5  
               WRITE (output_io, 3063) rho_extr_achs (rho_extr_abs) 
               WRITE (output_io, 3064) rho_extr_achs (rho_extr_ord) 
               IF (ier_num.ne.0) then 
                  CALL errlist 
                  CALL no_error 
               ENDIF 
!                                                                       
!     define the upper left corner 'ul'                                 
!                                                                       
      ELSEIF (str_comp (befehl, 'ul  ', 1, lbef, 4) ) then 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ier_num.eq.0) then 
                  IF (ianz.eq.3) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        j = 3 
                        DO i = 1, 3 
                        eck (i, j) = werte (i) 
                        ENDDO 
                        DO i = 1, 3 
                        vi (i, 1) = (eck (i, 2) - eck (i, 1) ) / divis (&
                        1)                                              
                        vi (i, 2) = (eck (i, 3) - eck (i, 1) ) / divis (&
                        2)                                              
                        ENDDO 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ENDIF 
!                                                                       
!     prepare Wilson statistics 'wilson'                                
!                                                                       
            ELSEIF (str_comp (befehl, 'wilson', 1, lbef, 6) ) then 
               IF (cr_v.le.0.0) then 
                  ier_num = - 35 
                  ier_typ = ER_APPL 
                  ier_msg (1) = 'A proper unit cell must be defined' 
      ier_msg (2)  = 'for this command to operate       ' 
               ELSE 
                  CALL wilson_calc 
               ENDIF 
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
!                                                                       
!     waiting for user input                                            
!                                                                       
            ELSEIF (str_comp (befehl, 'wait', 3, lbef, 4) ) then 
               CALL do_input (zeile, lp) 
!                                                                       
!     unknown command                                                   
!                                                                       
            ELSE 
               ier_num = - 8 
               ier_typ = ER_COMM 
            ENDIF 
         ENDIF 
      ENDIF 
      IF (ier_num.ne.0) THEN 
         CALL errlist 
         IF (ier_sta.ne.ER_S_LIVE) THEN 
            IF (lmakro .OR. lmakro_error) THEN  ! Error within macro or termination errror
               IF(sprompt /= prompt ) THEN
                  ier_num = -10
                  ier_typ = ER_COMM
                  ier_msg(1) = ' Error occured in '//prompt(8:LEN_TRIM(prompt))
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
 9999 CONTINUE 
!
      prompt = orig_prompt
!
 3060 FORMAT    (/' reciprocal layer    :'/                             &
     &                  ' lower left  corner  : ',3(2x,f9.4)/           &
     &                  ' lower right corner  : ',3(2x,f9.4)/           &
     &                  ' upper left  corner  : ',3(2x,f9.4))           
 3061 FORMAT    (/' number of points '/                                 &
     &                  ' horizontally        : ',2x,i4/                &
     &                  ' vertically          : ',2x,i4)                
 3062 FORMAT    (' increment vectors   :',                              &
     &             7x,'h',10x,'k',10x,'l',15x,'A**-1'/                  &
     &             ' horizontally        : ',3(2x,f9.4) ,7x,F9.4/       &
     &             ' vertically          : ',3(2x,f9.4) ,7x,F9.4/       &
     &            /' Angle               : ',2x,f9.4,' degrees'/        &
     &             ' Ratio/Aver   v/h    : ',2(2x,F9.4))                
 3063 FORMAT    (' Abszisse            : ', a) 
 3064 FORMAT    (' Ordinate            : ', a) 
 3260 FORMAT    (/' real space layer    :'/                             &
     &                  ' lower left  corner  : ',3(2x,f9.4)/           &
     &                  ' lower right corner  : ',3(2x,f9.4)/           &
     &                  ' upper left  corner  : ',3(2x,f9.4))           
 3262 FORMAT    (' increment vectors   :',                              &
     &             7x,'x',10x,'y',10x,'z',15x,'A'/                      &
     &             ' horizontally        : ',3(2x,f9.4) ,7x,F9.4/       &
     &             ' vertically          : ',3(2x,f9.4) ,7x,F9.4/       &
     &            /' Angle               : ',2x,f9.4,' degrees'/        &
     &             ' Ratio/Aver    v/h   : ',2(2x,F9.4))                
!                                                                       
      END SUBROUTINE patterson                      
!*****7*****************************************************************
      SUBROUTINE do_patters (inverse_type, ltwo_files, vr, acentric) 
!-                                                                      
!     Reads the files, if necessary calculates the values and           
!     calls the patterson function                                      
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE diffuse_mod 
      USE fourier_sup
      USE intens_mod 
      USE inverse_mod 
      USE output_mod 
      USE patters_mod 
      USE recipro_mod 
      USE quad_mod
!
      USE errlist_mod 
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER, INTENT(INOUT)  :: inverse_type
!
      LOGICAL ltwo_files 
      LOGICAL acentric 
      REAL vr 
!                                                                       
      INTEGER ifa, ifb 
      PARAMETER (ifa = 21, ifb = 22) 
!                                                                       
      INTEGER KUPL, GNU, SHELXL, HKLF4
      PARAMETER (KUPL = 0, GNU = 3, SHELXL = 4, HKLF4 = 5) 
!                                                                       
      CHARACTER(1024) line 
      INTEGER i, j, ii, k, itic 
      INTEGER extr_ima 
      INTEGER nx (2), ny (2) 
      INTEGER ihkl (3) 
      LOGICAL l_incl 
      REAL h (3) 
      REAL y (3), phase 
      REAL xmin (2), xmax (2) 
      REAL ymin (2), ymax (2) 
      REAL x1, x2, y1, y2, zz1, zz2, zz3, zz4 
      REAL z1 (10000), z2 (10000) 
      REAL e_f, dummy 
      REAL dstar2 
      COMPLEX (KIND=KIND(0.0D0)) ::a_b 
!                                                                       
!     INTEGER e_hist 
!     REAL fj2 
!     REAL quad 
!                                                                       
      zz1 = 0.0 
      zz2 = 0.0 
      zz3 = 0.0 
      zz4 = 0.0 
!                                                                       
      IF (ftyp.eq.HKLF4) then 
         CALL wilson_calc 
      ENDIF 
!                                                                       
      CALL oeffne (ifa, rho_file (1) , 'old') 
      IF (ier_num.ne.0) then 
         RETURN 
      ENDIF 
      IF (ltwo_files) then 
         CALL oeffne (ifb, rho_file (2) , 'old') 
         IF (ier_num.ne.0) then 
            RETURN 
         ENDIF 
      ENDIF 
!                                                                       
      IF (patt_accu.eq.PATT_INIT) then 
         ii = 0 
         DO i = 1, rho_inc (1) 
         DO j = 1, rho_inc (2) 
         ii = ii + 1 
         csf (ii) = cmplx (0.0D0, 0.0D0) 
         ENDDO 
         ENDDO 
      ENDIF 
!                                                                       
      IF (ftyp.eq.SHELXL.and. (.not.patt_rsym) ) then 
   10    CONTINUE 
!                                                                       
!     ----read hkl,fobs,fcalc,phase                                     
!                                                                       
         READ (ifa, 1111, end = 20, err = 900) ihkl, zz1, zz2, zz3 
 1111 FORMAT   (3I4,2F10.2,F7.2) 
         y (1) = real (ihkl (1) ) 
         y (2) = real (ihkl (2) ) 
         y (3) = real (ihkl (3) ) 
         zz4 = zz3 
!                                                                       
!     ------Calculate appropriate value, contribute to invers/Patterson 
!                                                                       
         CALL set_patt_value_s (a_b, inverse_type, rho_type, zz1, zz2,  &
         zz4)                                                           
         a_b = a_b * patt_scale * vr 
         CALL calc_patters (a_b, y) 
         GOTO 10 
   20    CONTINUE 
      ELSEIF (ftyp.eq.SHELXL.and.patt_rsym) then 
   15    CONTINUE 
!                                                                       
!     ----read hkl,fobs,fcalc,phase                                     
!                                                                       
         READ (ifa, *, end = 25, err = 900) h, zz1, zz2, zz3 
!                                                                       
!     ----Apply reciprocal space symmetry operations                    
!                                                                       
         DO k = 1, rec_n_sym 
         phase = 0.0 
         DO i = 1, 3 
         y (i) = 0 
         DO j = 1, 3 
         y (i) = y (i) + rec_sym (i, j, k) * h (j) 
         ENDDO 
         phase = phase+y (i) * rec_sym (i, 4, k) 
         ENDDO 
         zz4 = zz3 + 360.0 * phase 
!                                                                       
!     ------Calculate appropriate value, contribute to invers/Patterson 
!                                                                       
         CALL set_patt_value_s (a_b, inverse_type, rho_type, zz1, zz2,  &
         zz4)                                                           
         a_b = a_b * patt_scale * vr 
         CALL calc_patters (a_b, y) 
!                                                                       
!     ------If the space group is acentric, apply -1 operation          
!                                                                       
         IF (acentric) then 
            DO i = 1, 3 
            y (i) = - y (i) 
            ENDDO 
            zz4 = - zz4 
            CALL set_patt_value_s (a_b, inverse_type, rho_type, zz1,    &
            zz2, zz4)                                                   
            a_b = a_b * patt_scale * vr 
            CALL calc_patters (a_b, y) 
         ENDIF 
         ENDDO 
         GOTO 15 
   25    CONTINUE 
      ELSEIF (ftyp.eq.HKLF4.and.patt_rsym) then 
   30    CONTINUE 
!                                                                       
!     ----read hkl,intensity,sigma(intensity)                           
!                                                                       
         READ (ifa, *, end = 40, err = 900) h, zz1 
!                                                                       
!     ----Calculate |F|**2, |EF| or |E|**2                              
!         subtract sum(fj**2)                                           
!                                                                       
         dummy = zz1 
         IF (patt_mode.eq.PATT_NORMAL) then 
            IF (patt_origin.eq.PATT_SUBTRACT) then 
               dstar2 = quad (h, h, cr_rten) 
               dummy = zz1 * wilson_scale-fj2 (dstar2) 
            ELSE 
               dummy = zz1 * wilson_scale 
            ENDIF 
            a_b = cmplx (dummy, 0.0) 
         ELSEIF (patt_mode.eq.PATT_SHARP) then 
            dstar2 = quad (h, h, cr_rten) 
            j = e_hist (h) 
            e_f = sqrt (abs (zz1) / e_aver_f2 (j) ) 
            IF (patt_origin.eq.PATT_SUBTRACT) then 
               dstar2 = quad (h, h, cr_rten) 
               dummy = sqrt (abs (zz1) ) * e_f * sqrt (wilson_scale)    &
               - sqrt (fj2 (dstar2) )                                   
            ELSE 
               dummy = sqrt (abs (zz1) ) * e_f * sqrt (wilson_scale) 
            ENDIF 
            a_b = cmplx (dummy, 0.0) 
         ELSEIF (patt_mode.eq.PATT_SUPER) then 
            dstar2 = quad (h, h, cr_rten) 
            j = e_hist (h) 
            e_f = sqrt (abs (zz1) / e_aver_f2 (j) ) 
            IF (patt_origin.eq.PATT_SUBTRACT) then 
               dummy = e_f**2 - 1.0 
            ELSE 
               dummy = e_f**2 
            ENDIF 
            a_b = cmplx (dummy, 0.0) 
         ENDIF 
!                                                                       
!     ----Apply reciprocal space symmetry operations                    
!                                                                       
         DO k = 1, rec_n_sym 
         phase = 0.0 
         DO i = 1, 3 
         y (i) = 0 
         DO j = 1, 3 
         y (i) = y (i) + rec_sym (i, j, k) * h (j) 
         ENDDO 
         phase = phase+y (i) * rec_sym (i, 4, k) 
         ENDDO 
         zz4 = zz3 + 360.0 * phase 
!                                                                       
!     ------Calculate appropriate value, contribute to invers/Patterson 
!                                                                       
         a_b = a_b * patt_scale * vr 
         CALL calc_patters (a_b, y) 
!                                                                       
!     ------If the space group is acentric, apply -1 operation          
!                                                                       
         IF (acentric) then 
            DO i = 1, 3 
            y (i) = - y (i) 
            ENDDO 
            zz4 = - zz4 
            a_b = cmplx (dummy, 0.0) 
            a_b = a_b * patt_scale * vr 
            CALL calc_patters (a_b, y) 
         ENDIF 
         ENDDO 
         GOTO 30 
   40    CONTINUE 
      ELSEIF (ftyp.eq.HKLF4.and. (.not.patt_rsym) ) then 
   50    CONTINUE 
!                                                                       
!     ----read hkl,intensity,sigma(intensity)                           
!                                                                       
         READ (ifa, *, end = 60, err = 900) y, zz1 
!                                                                       
!     ----Calculate |F|**2, |EF| or |E|**2                              
!         subtract sum(fj**2)                                           
!                                                                       
         dummy = zz1 
         IF (patt_mode.eq.PATT_NORMAL) then 
            IF (patt_origin.eq.PATT_SUBTRACT) then 
               dstar2 = quad (y, y, cr_rten) 
               dummy = zz1 * wilson_scale-fj2 (dstar2) 
            ELSE 
               dummy = zz1 * wilson_scale 
            ENDIF 
            a_b = cmplx (dummy, 0.0) 
         ELSEIF (patt_mode.eq.PATT_SHARP) then 
            dummy = 0.0 
            dstar2 = quad (y, y, cr_rten) 
            j = e_hist (y) 
            IF (e_aver_f2 (j) .gt.0.0) then 
               e_f = sqrt (abs (zz1) ) / e_aver_f2 (j) 
               IF (patt_origin.eq.PATT_SUBTRACT) then 
                  dummy = sqrt (abs (zz1) ) * e_f * sqrt (wilson_scale) &
                  - sqrt (fj2 (dstar2) )                                
               ELSE 
                  dummy = sqrt (abs (zz1) ) * e_f * sqrt (wilson_scale) 
               ENDIF 
            ENDIF 
            a_b = cmplx (dummy, 0.0) 
         ELSEIF (patt_mode.eq.PATT_SUPER) then 
            dummy = 0.0 
            dstar2 = quad (y, y, cr_rten) 
            j = e_hist (y) 
            IF (e_aver_f2 (j) .gt.0.0) then 
               e_f = sqrt (abs (zz1) ) / e_aver_f2 (j) 
               IF (patt_origin.eq.PATT_SUBTRACT) then 
                  dummy = e_f**2 - 1.0 
               ELSE 
                  dummy = e_f**2 
               ENDIF 
            ENDIF 
            a_b = cmplx (dummy, 0.0) 
         ENDIF 
         a_b = cmplx (dummy, 0.0) * patt_scale * vr 
         CALL calc_patters (a_b, y) 
         GOTO 50 
   60    CONTINUE 
      ELSE 
         extr_ima = 6 - extr_abs - extr_ord 
!                                                                       
!     --One dimensional input files                                     
!                                                                       
         IF (inc (1) .eq.1.or.inc (2) .eq.1) then 
            ii = max (inc (1), inc (2) ) 
            itic = max (5, ii / 10) 
            DO i = 1, ii 
            IF (mod (i, itic) .eq.0) write (output_io, * ) i 
            READ (ifa, *, end = 900, err = 900) h (extr_abs), z1 (1) 
            IF (ltwo_files) then 
               READ (ifb, *, end = 900, err = 900) h (extr_abs),        &
               z2 (1)                                                   
            ELSE 
               z2 (1) = 0.0 
            ENDIF 
!           write (output_io,*) h(extr_abs),z1(1),z2(1)                 
            j = 1 
            DO k = 1, 3 
            h (k) = eck (k, 1) + vi (k, 1) * float (i - 1) + vi (k, 2)  &
            * float (j - 1)                                             
            ENDDO 
            CALL set_patt_value (a_b, inverse_type, rho_type, z1, z2, 1,&
            patt_excl9999, patt_excl_val, l_incl)                       
            IF (l_incl) then 
               a_b = a_b * patt_scale * vr 
               CALL calc_patters (a_b, h) 
            ENDIF 
            ENDDO 
!                                                                       
!     --Two dimensional Files (KUPL Format)                             
!                                                                       
         ELSE 
            IF (ftyp.eq.KUPL) then 
               ier_num = - 3 
               ier_typ = ER_IO 
               READ (ifa, *, err = 900, end = 900) nx (1), ny (1) 
               READ (ifa, *, err = 900, end = 900) xmin (1), xmax (1),  &
               ymin (1), ymax (1)                                       
               IF (ltwo_files) then 
                  READ (ifb, *, err = 900, end = 900) nx (2), ny (2) 
                  READ (ifb, *, err = 900, end = 900) xmin (2), xmax (2)&
                  , ymin (2), ymax (2)                                  
               ELSE 
                  nx (2) = nx (1) 
                  ny (2) = ny (1) 
                  xmin (2) = xmin (1) 
                  xmax (2) = xmax (1) 
                  ymin (2) = ymin (1) 
                  ymax (2) = ymax (1) 
               ENDIF 
               ier_num = 0 
               ier_typ = ER_NONE 
               IF (nx (1) .eq.nx (2) .and.ny (1) .eq.ny (2) .and.xmin ( &
               1) .eq.xmin (2) .and.ymin (1) .eq.ymin (2) .and.xmax (1) &
               .eq.xmax (2) .and.ymax (1) .eq.ymax (2) ) then           
                  itic = max (5, inc (2) / 10) 
                  DO j = 1, ny (1) 
                  IF (mod (j, itic) .eq.0) write (output_io, * ) j 
                  ier_num = - 3 
                  ier_typ = ER_IO 
                  READ (ifa, *, err = 900, end = 900) (z1 (i), i = 1,   &
                  nx (1) )                                              
                  IF (ltwo_files) then 
                     READ (ifb, *, err = 900, end = 900) (z2 (i),       &
                     i = 1, nx (1) )                                    
                  ENDIF 
                  ier_num = 0 
                  ier_typ = ER_NONE 
                  DO i = 1, inc (1) 
                  CALL set_patt_value (a_b, inverse_type, rho_type, z1, &
                  z2, i, patt_excl9999, patt_excl_val, l_incl)          
                  IF (l_incl) then 
                     DO k = 1, 3 
                     h (k) = eck (k, 1) + vi (k, 1) * float (i - 1)     &
                     + vi (k, 2) * float (j - 1)                        
                     ENDDO 
                     a_b = a_b * patt_scale * vr 
                     CALL calc_patters (a_b, h) 
                  ENDIF 
                  ENDDO 
                  ENDDO 
               ELSE 
                  ier_num = - 22 
                  ier_typ = ER_APPL 
               ENDIF 
            ELSEIF (ftyp.eq.GNU) then 
               ier_num = - 23 
               ier_typ = ER_APPL 
               itic = max (5, inc (2) / 10) 
               DO j = 1, inc (2) 
               IF (mod (j, itic) .eq.0) write (output_io, * ) j 
               DO i = 1, inc (1) 
               READ (ifa, *, end = 900, err = 900) x1, y1, z1 (1),      &
               zz1                                                      
               IF (ltwo_files) then 
                  READ (ifb, *, end = 900, err = 900) x2, y2, z2 (1),   &
                  zz2                                                   
               ENDIF 
               IF (.not.ltwo_files.or. (ltwo_files.and. (               &
               x1.eq.x2.and.y1.eq.y2.and.zz1.eq.zz2) ) ) then           
                  h (extr_abs) = x1 
                  h (extr_ord) = y1 
                  h (extr_ima) = zz1 
                  CALL set_patt_value (a_b, inverse_type, rho_type, z1, &
                  z2, 1, patt_excl9999, patt_excl_val, l_incl)          
                  IF (l_incl) then 
                     a_b = a_b * patt_scale * vr 
                     CALL calc_patters (a_b, h) 
                  ENDIF 
               ELSE 
                  ier_num = - 24 
                  ier_typ = ER_APPL 
                  GOTO 900 
               ENDIF 
               ENDDO 
               READ (ifa, 1000, end = 900, err = 900) line 
               IF (line.eq.' '.or.line.eq.char (13) ) then 
                  IF (ltwo_files) then 
                     READ (ifb, 1000, end = 900, err = 900) line 
                  ENDIF 
                  IF (line.ne.' '.and.line.ne.char (13) ) then 
                     ier_num = - 23 
                     ier_typ = ER_APPL 
                     GOTO 900 
                  ENDIF 
               ENDIF 
               ENDDO 
               ier_num = 0 
               ier_typ = ER_NONE 
            ENDIF 
         ENDIF 
      ENDIF 
!                                                                       
  900 CONTINUE 
!                                                                       
      CLOSE (ifa) 
      IF (ltwo_files) then 
         CLOSE (ifb) 
      ENDIF 
!                                                                       
 1000 FORMAT    (a) 
      END SUBROUTINE do_patters                     
!*****7*****************************************************************
      SUBROUTINE calc_patters (a_b, h) 
!-                                                                      
!     Calculates the patterson or inverse Fourier at all points         
!     in real space.                                                    
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE diffuse_mod 
      USE fourier_sup
      USE patters_mod 
      USE precision_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      REAL h (3) 
!                                                                       
      REAL(PREC_DP) xarg0, xincu, xincv 
      INTEGER iarg, iarg0, iincu, iincv, iadd 
      INTEGER i, j, ii 
      COMPLEX (KIND=KIND(0.0D0)) ::a_b 
!                                                                       
      INTEGER ISHFT, IAND 
!                                                                       
      CALL four_cexpt 
      ii = max (5, rho_inc (2) / 10) 
!                                                                       
      xarg0 = rho_eck (1, 1) * h (1) + rho_eck (2, 1) * h (2) + rho_eck &
      (3, 1) * h (3)                                                    
      xincu = rho_vi (1, 1) * h (1) + rho_vi (2, 1) * h (2) + rho_vi (3,&
      1) * h (3)                                                        
      xincv = rho_vi (1, 2) * h (1) + rho_vi (2, 2) * h (2) + rho_vi (3,&
      2) * h (3)                                                        
!                                                                       
      iarg0 = nint (64 * I2PI * (xarg0 - int (xarg0) + 1.0d0) ) 
      iincu = nint (64 * I2PI * (xincu - int (xincu) + 1.0d0) ) 
      iincv = nint (64 * I2PI * (xincv - int (xincv) + 1.0d0) ) 
!                                                                       
      iarg0 = patt_sign * iarg0 
      iarg = iarg0 
!                                                                       
      ii = 0 
!                                                                       
      DO j = 1, rho_inc (1) 
      DO i = 1, rho_inc (2) 
      iadd = ISHFT (iarg, - 6) 
      iadd = IAND (iadd, MASK) 
      ii = ii + 1 
      csf (ii) = csf (ii) + cex (iadd) * a_b 
      iarg = iarg + iincv * patt_sign 
      ENDDO 
      iarg = iarg0 + iincu * j * patt_sign 
      ENDDO 
!                                                                       
      END SUBROUTINE calc_patters                   
!*****7*******1*********************************************************
      SUBROUTINE set_patt_value (a_b, inverse_type, rho_type, z1, z2, i,&
      patt_excl9999, patt_excl_val, l_incl)                             
!-                                                                      
!     Sets the appropriate values for the inverse Fourier or            
!     Patterson.                                                        
!+                                                                      
      USE intens_mod 
      USE inverse_mod 
      IMPLICIT none 
!                                                                       
      INTEGER, INTENT(INOUT)  :: inverse_type
!                                                                       
      INTEGER rho_type (2) 
      LOGICAL patt_excl9999 
      LOGICAL l_incl 
      REAL z1 (10000), z2 (10000) 
      REAL patt_excl_val 
      COMPLEX (KIND=KIND(0.0D0)) ::a_b 
!                                                                       
      INTEGER i
!                                                                       
      REAL cosd, sind 
!                                                                       
      l_incl = .true. 
      IF (inverse_type.eq.INV_INV) then 
         IF (patt_excl9999) then 
            IF (z1 (i) .eq.patt_excl_val.or.z2 (i) .eq.patt_excl_val)   &
            then                                                        
               l_incl = .false. 
               RETURN 
            ENDIF 
         ENDIF 
         IF (rho_type (1) .eq.REAL_PART.and.rho_type (2) .eq.IMAG_PART) &
         then                                                           
            a_b = cmplx (z1 (i), z2 (i) ) 
         ELSEIF (rho_type (1) .eq.INTENSITY.and.rho_type (2)            &
         .eq.PHASE_ANG) then                                            
!                                                                       
!     ----Ignore negative intensities                                   
!                                                                       
            IF (z1 (i) .lt.0.0) then 
               l_incl = .false. 
               RETURN 
            ENDIF 
            a_b = cmplx (sqrt (z1 (i) ) * cosd (z2 (i) ), sqrt (z1 (i) )&
            * sind (z2 (i) ) )                                          
         ELSEIF (rho_type (1) .eq.AMPLITUDE.and.rho_type (2)            &
         .eq.PHASE_ANG) then                                            
            a_b = cmplx (z1 (i) * cosd (z2 (i) ), z1 (i) * sind (z2 (i) &
            ) )                                                         
         ENDIF 
      ELSEIF (inverse_type.eq.INV_PATT) then 
         IF (patt_excl9999) then 
            IF (rho_type (1) .eq.INTENSITY) then 
               IF (z1 (i) .eq.patt_excl_val) then 
                  l_incl = .false. 
                  RETURN 
               ENDIF 
               a_b = cmplx (z1 (i), 0.0) 
            ELSEIF (rho_type (1) .eq.REAL_PART.and.rho_type (2)         &
            .eq.IMAG_PART) then                                         
               IF (z1 (i) .eq.patt_excl_val.or.z2 (i) .eq.patt_excl_val)&
               then                                                     
                  l_incl = .false. 
                  RETURN 
               ENDIF 
               a_b = cmplx (z1 (i) * z1 (i) + z2 (i) * z2 (i), 0.0) 
            ELSEIF (rho_type (1) .eq.AMPLITUDE) then 
               IF (z1 (i) .eq.patt_excl_val) then 
                  l_incl = .false. 
                  RETURN 
               ENDIF 
               a_b = cmplx (z1 (i) * z1 (i), 0.0) 
            ENDIF 
         ELSE 
            IF (rho_type (1) .eq.INTENSITY) then 
               a_b = cmplx (z1 (i), 0.0) 
            ELSEIF (rho_type (1) .eq.REAL_PART.and.rho_type (2)         &
            .eq.IMAG_PART) then                                         
               a_b = cmplx (z1 (i) * z1 (i) + z2 (i) * z2 (i), 0.0) 
            ELSEIF (rho_type (1) .eq.AMPLITUDE) then 
               a_b = cmplx (z1 (i) * z1 (i), 0.0) 
            ENDIF 
         ENDIF 
      ENDIF 
!                                                                       
      END SUBROUTINE set_patt_value                 
!*****7*****************************************************************
      SUBROUTINE set_patt_value_s (a_b, inverse_type, rho_type, zz1,    &
      zz2, zz3)                                                         
!-                                                                      
!     Sets the appropriate values for the inverse Fourier or            
!     Patterson. Special Version for SHELXL list file 5 Input           
!+                                                                      
      USE intens_mod 
      USE inverse_mod 
      IMPLICIT none 
!                                                                       
      INTEGER, INTENT(INOUT)  :: inverse_type
!                                                                       
      INTEGER rho_type (2) 
      REAL zz1, zz2, zz3 
      COMPLEX (KIND=KIND(0.0D0)) ::a_b 
!                                                                       
      REAL cosd, sind 
!                                                                       
      IF (inverse_type.eq.INV_INV) then 
         IF (rho_type (1) .eq.AMPLITUDE) then 
            a_b = cmplx (zz1 * cosd (zz3), zz1 * sind (zz3) ) 
         ELSEIF (rho_type (1) .eq.FCALC) then 
            a_b = cmplx (zz2 * cosd (zz3), zz2 * sind (zz3) ) 
         ENDIF 
      ELSEIF (inverse_type.eq.INV_PATT) then 
         IF (rho_type (1) .eq.INTENSITY) then 
            a_b = cmplx (zz1 * zz1, 0.0) 
         ELSEIF (rho_type (1) .eq.AMPLITUDE) then 
            a_b = cmplx (zz1 * zz1, 0.0) 
         ELSEIF (rho_type (1) .eq.FCALC) then 
            a_b = cmplx (zz2 * zz2, 0.0) 
         ENDIF 
      ELSEIF (inverse_type.eq.INV_DIFF) then 
         a_b = cmplx ( (zz1 - zz2) * cosd (zz3), (zz1 - zz2) * sind (   &
         zz3) )                                                         
      ENDIF 
!                                                                       
      END SUBROUTINE set_patt_value_s               
!*****7*****************************************************************
      SUBROUTINE e_create (e_io, cpara, lpara, io_mode) 
!-                                                                      
!     reads the input file and creates normalized structure factors     
!                                                                       
!     Version  : 1.0                                                    
!     Date     : 27 Nov 2000                                            
!                                                                       
!     Author   : R.B. Neder (reinhard.neder@fau.de)      
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE metric_mod
      USE patters_mod 
!
      USE errlist_mod 
      USE prompt_mod 
      IMPLICIT none 
       
!                                                                       
      INTEGER MAXW 
      PARAMETER (MAXW = 20) 
      INTEGER HKLF4 
      PARAMETER (HKLF4 = 5) 
!                                                                       
      CHARACTER ( * ) cpara 
      CHARACTER(1024) outfile 
      INTEGER lpara 
      INTEGER i, n, j, lcomm 
      INTEGER e_io, e_io_hkl 
      LOGICAL io_mode 
      LOGICAL lsuccess 
      REAL h (3) 
      REAL zz1, zz2, zz3, zz4 
!                                                                       
      INTEGER e_graph (0:30) 
      INTEGER e_graph_hk0 (0:30) 
      INTEGER e_graph_h0l (0:30) 
      INTEGER e_graph_0kl (0:30) 
      INTEGER e_graph_hhl (0:30) 
      INTEGER e_graph_hMh0l (0:30) 
      INTEGER e_graph_hhM2hl (0:30) 
      INTEGER e_graph_0k0 (0:30) 
!                                                                       
      INTEGER e_num_hkl (20) 
      INTEGER e_num_hk0 (20) 
      INTEGER e_num_h0l (20) 
      INTEGER e_num_0kl (20) 
      INTEGER e_num_0k0 (20) 
      INTEGER e_num_hhl (20) 
      INTEGER e_num_hMh0l (20) 
      INTEGER e_num_hhM2hl (20) 
      INTEGER n_hkl 
      INTEGER n_hk0 
      INTEGER n_h0l 
      INTEGER n_0kl 
      INTEGER n_0k0 
      INTEGER n_hhl 
      INTEGER n_hMh0l 
      INTEGER n_hhM2hl 
      INTEGER ifa 
      PARAMETER (ifa = 21) 
!      real    e_aver_f2    (MAXW)                                      
      REAL e_aver_f2_hk0 (MAXW) 
      REAL e_aver_f2_h0l (MAXW) 
      REAL e_aver_f2_0kl (MAXW) 
      REAL e_aver_f2_0k0 (MAXW) 
      REAL e_aver_f2_hhl (MAXW) 
      REAL e_aver_f2_hMh0l (MAXW) 
      REAL e_aver_f2_hhM2hl (MAXW) 
      REAL e_f 
      REAL e_e1_hkl, e_e2_hkl, e_e3_hkl, e_e4_hkl 
      REAL e_e5_hkl, e_e6_hkl, e_e2m1_hkl 
      REAL e_e2m1_2_hkl, e_e2m1_3_hkl, e_abs_e2m1_3_hkl 
      REAL e_e1_hk0, e_e2_hk0, e_e3_hk0, e_e4_hk0 
      REAL e_e5_hk0, e_e6_hk0, e_e2m1_hk0 
      REAL e_e2m1_2_hk0, e_e2m1_3_hk0, e_abs_e2m1_3_hk0 
      REAL e_e1_h0l, e_e2_h0l, e_e3_h0l, e_e4_h0l 
      REAL e_e5_h0l, e_e6_h0l, e_e2m1_h0l 
      REAL e_e2m1_2_h0l, e_e2m1_3_h0l, e_abs_e2m1_3_h0l 
      REAL e_e1_0k0, e_e2_0k0, e_e3_0k0, e_e4_0k0 
      REAL e_e5_0k0, e_e6_0k0, e_e2m1_0k0 
      REAL e_e2m1_2_0k0, e_e2m1_3_0k0, e_abs_e2m1_3_0k0 
      REAL e_e1_0kl, e_e2_0kl, e_e3_0kl, e_e4_0kl 
      REAL e_e5_0kl, e_e6_0kl, e_e2m1_0kl 
      REAL e_e2m1_2_0kl, e_e2m1_3_0kl, e_abs_e2m1_3_0kl 
      REAL e_e1_hhl, e_e2_hhl, e_e3_hhl, e_e4_hhl 
      REAL e_e5_hhl, e_e6_hhl, e_e2m1_hhl 
      REAL e_e2m1_2_hhl, e_e2m1_3_hhl, e_abs_e2m1_3_hhl 
      REAL e_e1_hMh0l, e_e2_hMh0l, e_e3_hMh0l, e_e4_hMh0l 
      REAL e_e5_hMh0l, e_e6_hMh0l, e_e2m1_hMh0l 
      REAL e_e2m1_2_hMh0l, e_e2m1_3_hMh0l, e_abs_e2m1_3_hMh0l 
      REAL e_e1_hhM2hl, e_e2_hhM2hl, e_e3_hhM2hl, e_e4_hhM2hl 
      REAL e_e5_hhM2hl, e_e6_hhM2hl, e_e2m1_hhM2hl 
      REAL e_e2m1_2_hhM2hl, e_e2m1_3_hhM2hl, e_abs_e2m1_3_hhM2hl 
!                                                                       
      REAL i_aver_all 
      REAL i_aver_C, i_aver_B, i_aver_A, i_aver_F, i_aver_I 
      REAL i_aver_RO, i_aver_RR 
!                                                                       
      INTEGER n_all 
      INTEGER n_C, n_B, n_A, n_F, n_I 
      INTEGER n_RO, n_RR 
!                                                                       
      LOGICAL latt_P 
      LOGICAL latt_C, latt_B, latt_A, latt_F, latt_I 
      LOGICAL latt_RO, latt_RR 
!                                                                       
!     INTEGER e_hist 
!     REAL do_blen 
      REAL null (3) 
      REAL dstar 
      LOGICAL lspace 
      PARAMETER (lspace = .false.) 
!                                                                       
      DATA null / 0.0, 0.0, 0.0 / 
!                                                                       
      IF (cr_v.le.0.0) then 
         ier_num = - 35 
         ier_typ = ER_APPL 
         ier_msg (1) = 'A proper unit cell must be defined' 
      ier_msg (2)  = 'for this command to operate       ' 
         RETURN 
      ENDIF 
!                                                                       
      e_io_hkl = 45 
!                                                                       
      e_e1_hkl = 0.0 
      e_e2_hkl = 0.0 
      e_e3_hkl = 0.0 
      e_e4_hkl = 0.0 
      e_e5_hkl = 0.0 
      e_e6_hkl = 0.0 
      e_e2m1_hkl = 0.0 
      e_e2m1_2_hkl = 0.0 
      e_e2m1_3_hkl = 0.0 
      e_abs_e2m1_3_hkl = 0.0 
!                                                                       
      e_e1_hk0 = 0.0 
      e_e2_hk0 = 0.0 
      e_e3_hk0 = 0.0 
      e_e4_hk0 = 0.0 
      e_e5_hk0 = 0.0 
      e_e6_hk0 = 0.0 
      e_e2m1_hk0 = 0.0 
      e_e2m1_2_hk0 = 0.0 
      e_e2m1_3_hk0 = 0.0 
      e_abs_e2m1_3_hk0 = 0.0 
!                                                                       
      e_e1_h0l = 0.0 
      e_e2_h0l = 0.0 
      e_e3_h0l = 0.0 
      e_e4_h0l = 0.0 
      e_e5_h0l = 0.0 
      e_e6_h0l = 0.0 
      e_e2m1_h0l = 0.0 
      e_e2m1_2_h0l = 0.0 
      e_e2m1_3_h0l = 0.0 
      e_abs_e2m1_3_h0l = 0.0 
!                                                                       
      e_e1_0kl = 0.0 
      e_e2_0kl = 0.0 
      e_e3_0kl = 0.0 
      e_e4_0kl = 0.0 
      e_e5_0kl = 0.0 
      e_e6_0kl = 0.0 
      e_e2m1_0kl = 0.0 
      e_e2m1_2_0kl = 0.0 
      e_e2m1_3_0kl = 0.0 
      e_abs_e2m1_3_0kl = 0.0 
!                                                                       
      e_e1_0k0 = 0.0 
      e_e2_0k0 = 0.0 
      e_e3_0k0 = 0.0 
      e_e4_0k0 = 0.0 
      e_e5_0k0 = 0.0 
      e_e6_0k0 = 0.0 
      e_e2m1_0k0 = 0.0 
      e_e2m1_2_0k0 = 0.0 
      e_e2m1_3_0k0 = 0.0 
      e_abs_e2m1_3_0k0 = 0.0 
!                                                                       
      e_e1_hhl = 0.0 
      e_e2_hhl = 0.0 
      e_e3_hhl = 0.0 
      e_e4_hhl = 0.0 
      e_e5_hhl = 0.0 
      e_e6_hhl = 0.0 
      e_e2m1_hhl = 0.0 
      e_e2m1_2_hhl = 0.0 
      e_e2m1_3_hhl = 0.0 
      e_abs_e2m1_3_hhl = 0.0 
!                                                                       
      e_e1_hMh0l = 0.0 
      e_e2_hMh0l = 0.0 
      e_e3_hMh0l = 0.0 
      e_e4_hMh0l = 0.0 
      e_e5_hMh0l = 0.0 
      e_e6_hMh0l = 0.0 
      e_e2m1_hMh0l = 0.0 
      e_e2m1_2_hMh0l = 0.0 
      e_e2m1_3_hMh0l = 0.0 
      e_abs_e2m1_3_hMh0l = 0.0 
!                                                                       
      e_e1_hhM2hl = 0.0 
      e_e2_hhM2hl = 0.0 
      e_e3_hhM2hl = 0.0 
      e_e4_hhM2hl = 0.0 
      e_e5_hhM2hl = 0.0 
      e_e6_hhM2hl = 0.0 
      e_e2m1_hhM2hl = 0.0 
      e_e2m1_2_hhM2hl = 0.0 
      e_e2m1_3_hhM2hl = 0.0 
      e_abs_e2m1_3_hhM2hl = 0.0 
!                                                                       
      zz1 = 0.0 
      zz2 = 0.0 
      zz3 = 0.0 
      zz4 = 0.0 
      n = 0 
      DO i = 1, 20 
      e_num_hkl (i) = 0 
      e_num_hk0 (i) = 0 
      e_num_h0l (i) = 0 
      e_num_0kl (i) = 0 
      e_num_0k0 (i) = 0 
      e_num_hhl (i) = 0 
      e_num_hMh0l (i) = 0 
      e_num_hhM2hl (i) = 0 
      ENDDO 
      DO i = 1, 20 
      e_aver_f2 (i) = 0.0 
      e_aver_f2_hk0 (i) = 0.0 
      e_aver_f2_h0l (i) = 0.0 
      e_aver_f2_0kl (i) = 0.0 
      e_aver_f2_hhl (i) = 0.0 
      e_aver_f2_hMh0l (i) = 0.0 
      e_aver_f2_hhM2hl (i) = 0.0 
      e_aver_f2_0k0 (i) = 0.0 
      ENDDO 
!                                                                       
      IF (ftyp.eq.HKLF4) then 
         n_all = 0 
         n_C = 0 
         n_A = 0 
         n_B = 0 
         n_I = 0 
         n_F = 0 
         n_RO = 0 
         n_RR = 0 
         i_aver_all = 0.0 
         i_aver_C = 0.0 
         i_aver_B = 0.0 
         i_aver_A = 0.0 
         i_aver_I = 0.0 
         i_aver_F = 0.0 
         i_aver_RO = 0.0 
         i_aver_RR = 0.0 
!                                                                       
         CALL oeffne (ifa, rho_file (1) , 'old') 
         n = 0 
   50    CONTINUE 
!                                                                       
!     ----read hkl,intensity,sigma(intensity)                           
!                                                                       
         READ (ifa, *, end = 60, err = 950) h, zz1 
         n = n + 1 
         IF (.not. (h (1) .eq.0.and.h (2) .eq.0.and.h (3) .eq.0) ) then 
            n_all = n_all + 1 
            i_aver_all = i_aver_all + abs (zz1) 
            IF (mod (nint (h (1) + h (2) ), 2) .ne.0) then 
               i_aver_C = i_aver_C + abs (zz1) 
               n_C = n_C + 1 
            ENDIF 
            IF (mod (nint (h (1) + h (3) ), 2) .ne.0) then 
               i_aver_B = i_aver_B + abs (zz1) 
               n_B = n_B + 1 
            ENDIF 
            IF (mod (nint (h (2) + h (3) ), 2) .ne.0) then 
               i_aver_A = i_aver_A + abs (zz1) 
               n_A = n_A + 1 
            ENDIF 
            IF (mod (nint (h (1) + h (2) + h (3) ), 2) .ne.0) then 
               i_aver_I = i_aver_I + abs (zz1) 
               n_I = n_I + 1 
            ENDIF 
            IF (.not. (mod (nint (h (1) + h (2) ), 2) .eq.0.and.mod (   &
            nint (h (1) + h (3) ), 2) .eq.0.and.mod (nint (h (2)        &
            + h (3) ), 2) .eq.0) ) then                                 
               i_aver_F = i_aver_F + abs (zz1) 
               n_F = n_F + 1 
            ENDIF 
            IF (mod (nint ( - h (1) + h (2) + h (3) ), 3) .ne.0) then 
               i_aver_RO = i_aver_RO + abs (zz1) 
               n_RO = n_RO + 1 
            ENDIF 
            IF (mod (nint (h (1) - h (2) + h (3) ), 3) .ne.0) then 
               i_aver_RR = i_aver_RR + abs (zz1) 
               n_RR = n_RR + 1 
            ENDIF 
         ENDIF 
         GOTO 50 
!                                                                       
   60    CONTINUE 
         CLOSE (ifa) 
         latt_p = .true. 
         latt_C = .false. 
         latt_A = .false. 
         latt_B = .false. 
         latt_I = .false. 
         latt_F = .false. 
         latt_RO = .false. 
         latt_RR = .false. 
         IF (n_F.gt.0.and.i_aver_F.lt.i_aver_all * 1.e-6) then 
            latt_F = .true. 
            latt_P = .false. 
         ELSEIF (n_C.gt.0.and.i_aver_C.lt.i_aver_all * 1.e-6) then 
            latt_C = .true. 
            latt_P = .false. 
         ELSEIF (n_A.gt.0.and.i_aver_A.lt.i_aver_all * 1.e-6) then 
            latt_A = .true. 
            latt_P = .false. 
         ELSEIF (n_B.gt.0.and.i_aver_B.lt.i_aver_all * 1.e-6) then 
            latt_B = .true. 
            latt_P = .false. 
         ELSEIF (n_I.gt.0.and.i_aver_I.lt.i_aver_all * 1.e-6) then 
            latt_I = .true. 
            latt_P = .false. 
         ELSEIF (n_RO.gt.0.and.i_aver_RO.lt.i_aver_all * 1.e-6) then 
            latt_RO = .true. 
            latt_P = .false. 
         ELSEIF (n_RR.gt.0.and.i_aver_RR.lt.i_aver_all * 1.e-6) then 
            latt_RR = .true. 
            latt_P = .false. 
         ENDIF 
!                                                                       
!                                                                       
         CALL oeffne (ifa, rho_file (1) , 'old') 
         n_all = 0 
         lsuccess = .false. 
         DO i = 1, n 
!                                                                       
!     ----read hkl,intensity,sigma(intensity)                           
!                                                                       
         READ (ifa, *, end = 904, err = 950) h, zz1 
         IF (.not. (h (1) .eq.0.and.h (2) .eq.0.and.h (3) .eq.0) ) then 
            IF (latt_P.or.latt_A.and.mod (nint (h (2) + h (3) ),        &
            2) .eq.0.or.latt_B.and.mod (nint (h (1) + h (3) ), 2)       &
            .eq.0.or.latt_C.and.mod (nint (h (1) + h (2) ), 2)          &
            .eq.0.or.latt_I.and.mod (nint (h (1) + h (2) + h (3) ),     &
            2) .eq.0.or.latt_F.and.mod (nint (h (1) + h (2) ), 2)       &
            .eq.0.and.mod (nint (h (1) + h (3) ), 2) .eq.0.and.mod (    &
            nint (h (2) + h (3) ), 2) .eq.0.or.latt_RO.and.mod (nint ( -&
            h (1) + h (2) + h (3) ), 3) .eq.0.or.latt_RR.and.mod (nint (&
            h (1) - h (2) + h (3) ), 3) .eq.0) then                     
               j = e_hist (h) 
               e_aver_f2 (j) = e_aver_f2 (j) + abs (zz1) 
               e_num_hkl (j) = e_num_hkl (j) + 1 
               IF (h (1) .eq.0) then 
                  e_aver_f2_0kl (j) = e_aver_f2_0kl (j) + abs (zz1) 
                  e_num_0kl (j) = e_num_0kl (j) + 1 
               ENDIF 
               IF (h (2) .eq.0) then 
                  e_aver_f2_h0l (j) = e_aver_f2_h0l (j) + abs (zz1) 
                  e_num_h0l (j) = e_num_h0l (j) + 1 
               ENDIF 
               IF (h (3) .eq.0) then 
                  e_aver_f2_hk0 (j) = e_aver_f2_hk0 (j) + abs (zz1) 
                  e_num_hk0 (j) = e_num_hk0 (j) + 1 
               ENDIF 
               IF (h (1) .eq.h (2) .or.h (1) .eq. - h (2) ) then 
                  e_aver_f2_hhl (j) = e_aver_f2_hhl (j) + abs (zz1) 
                  e_num_hhl (j) = e_num_hhl (j) + 1 
               ENDIF 
               IF (h (1) .eq. - h (2) ) then 
                  e_aver_f2_hMh0l (j) = e_aver_f2_hMh0l (j) + abs (zz1) 
                  e_num_hMh0l (j) = e_num_hMh0l (j) + 1 
               ENDIF 
               IF (h (1) .eq.h (2) ) then 
                  e_aver_f2_hhM2hl (j) = e_aver_f2_hhM2hl (j) + abs (   &
                  zz1)                                                  
                  e_num_hhM2hl (j) = e_num_hhM2hl (j) + 1 
               ENDIF 
               IF (h (1) .eq.0.and.h (3) .eq.0) then 
                  e_aver_f2_0k0 (j) = e_aver_f2_0k0 (j) + abs (zz1) 
                  e_num_0k0 (j) = e_num_0k0 (j) + 1 
               ENDIF 
            ENDIF 
         ENDIF 
         ENDDO 
         lsuccess = .true. 
  904    CONTINUE 
         IF (.not.lsuccess) then 
            RETURN 
         ENDIF 
      ELSE 
         CLOSE (ifa) 
         RETURN 
      ENDIF 
      CLOSE (ifa) 
      DO j = 1, 20 
      IF (e_num_hkl (j) .gt.0) then 
         e_aver_f2 (j) = sqrt (e_aver_f2 (j) / e_num_hkl (j) ) 
      ENDIF 
      IF (e_num_hk0 (j) .gt.0) then 
         e_aver_f2_hk0 (j) = sqrt (e_aver_f2_hk0 (j) / e_num_hk0 (j) ) 
      ENDIF 
      IF (e_num_h0l (j) .gt.0) then 
         e_aver_f2_h0l (j) = sqrt (e_aver_f2_h0l (j) / e_num_h0l (j) ) 
      ENDIF 
      IF (e_num_0kl (j) .gt.0) then 
         e_aver_f2_0kl (j) = sqrt (e_aver_f2_0kl (j) / e_num_0kl (j) ) 
      ENDIF 
      IF (e_num_0k0 (j) .gt.0) then 
         e_aver_f2_0k0 (j) = sqrt (e_aver_f2_0k0 (j) / e_num_0k0 (j) ) 
      ENDIF 
      IF (e_num_hhl (j) .gt.0) then 
         e_aver_f2_hhl (j) = sqrt (e_aver_f2_hhl (j) / e_num_hhl (j) ) 
      ENDIF 
      IF (e_num_hMh0l (j) .gt.0) then 
         e_aver_f2_hMh0l (j) = sqrt (e_aver_f2_hMh0l (j) / e_num_hMh0l (&
         j) )                                                           
      ENDIF 
      IF (e_num_hhM2hl (j) .gt.0) then 
         e_aver_f2_hhM2hl (j) = sqrt (e_aver_f2_hhM2hl (j) /            &
         e_num_hhM2hl (j) )                                             
      ENDIF 
      ENDDO 
!                                                                       
      IF (io_mode) then 
         CALL oeffne (ifa, rho_file (1) , 'old') 
         n_hkl = 0 
         n_0kl = 0 
         n_h0l = 0 
         n_hk0 = 0 
         n_0k0 = 0 
         n_hhl = 0 
         n_hMh0l = 0 
         n_hhM2hl = 0 
!                                                                       
!     --Write results                                                   
!                                                                       
         IF (e_io.ne.6) then 
            outfile = cpara (1:lpara) //'.ehkl' 
            CALL oeffne (e_io_hkl, outfile, 'unknown') 
         ENDIF 
!                                                                       
         DO i = 0, 30 
         e_graph (i) = 0 
         e_graph_hk0 (i) = 0 
         e_graph_h0l (i) = 0 
         e_graph_0kl (i) = 0 
         e_graph_hhl (i) = 0 
         e_graph_hMh0l (i) = 0 
         e_graph_hhM2hl (i) = 0 
         e_graph_0k0 (i) = 0 
         ENDDO 
         DO i = 1, n 
         READ (ifa, *, end = 903, err = 903) h, zz1, zz2 
         j = e_hist (h) 
         IF (e_aver_f2 (j) .ne.0.0) then 
            e_f = sqrt (abs (zz1) ) / e_aver_f2 (j) 
         ELSE 
            e_f = 0.0 
         ENDIF 
         dstar = do_blen (lspace, h, null) 
!                                                                       
         IF (e_io.ne.6) then 
            WRITE (e_io_hkl, 4000) nint (h (1) ), nint (h (2) ),        &
            nint (h (3) ), e_f**2, zz2 / (e_aver_f2 (j) ) **2           
         ENDIF 
!                                                                       
         IF (.not. (h (1) .eq.0.and.h (2) .eq.0.and.h (3) .eq.0) ) then 
            IF (latt_P.or.latt_A.and.mod (nint (h (2) + h (3) ),        &
            2) .eq.0.or.latt_B.and.mod (nint (h (1) + h (3) ), 2)       &
            .eq.0.or.latt_C.and.mod (nint (h (1) + h (2) ), 2)          &
            .eq.0.or.latt_I.and.mod (nint (h (1) + h (2) + h (3) ),     &
            2) .eq.0.or.latt_F.and.mod (nint (h (1) + h (2) ), 2)       &
            .eq.0.and.mod (nint (h (1) + h (3) ), 2) .eq.0.and.mod (    &
            nint (h (2) + h (3) ), 2) .eq.0.or.latt_RO.and.mod (nint ( -&
            h (1) + h (2) + h (3) ), 3) .eq.0.or.latt_RR.and.mod (nint (&
            h (1) - h (2) + h (3) ), 3) .eq.0) then                     
!                                                                       
               IF (e_f.le.3.0) then 
                  e_graph (nint (e_f * 10) ) = e_graph (nint (e_f * 10) &
                  ) + 1                                                 
                  e_graph (0) = e_graph (0) + 1 
                  IF (h (3) .eq.0) then 
                     e_graph_hk0 (nint (e_f * 10) ) = e_graph_hk0 (nint &
                     (e_f * 10) ) + 1                                   
                     e_graph_hk0 (0) = e_graph_hk0 (0) + 1 
                  ENDIF 
                  IF (h (2) .eq.0) then 
                     e_graph_h0l (nint (e_f * 10) ) = e_graph_h0l (nint &
                     (e_f * 10) ) + 1                                   
                     e_graph_h0l (0) = e_graph_h0l (0) + 1 
                  ENDIF 
                  IF (h (1) .eq.0) then 
                     e_graph_0kl (nint (e_f * 10) ) = e_graph_0kl (nint &
                     (e_f * 10) ) + 1                                   
                     e_graph_0kl (0) = e_graph_0kl (0) + 1 
                  ENDIF 
                  IF (h (1) .eq.h (2) .or.h (1) .eq. - h (2) ) then 
                     e_graph_hhl (nint (e_f * 10) ) = e_graph_hhl (nint &
                     (e_f * 10) ) + 1                                   
                     e_graph_hhl (0) = e_graph_hhl (0) + 1 
                  ENDIF 
                  IF (h (1) .eq. - h (2) ) then 
                     e_graph_hMh0l (nint (e_f * 10) ) = e_graph_hMh0l ( &
                     nint (e_f * 10) ) + 1                              
                     e_graph_hMh0l (0) = e_graph_hMh0l (0) + 1 
                  ENDIF 
                  IF (h (1) .eq.h (2) ) then 
                     e_graph_hhM2hl (nint (e_f * 10) ) = e_graph_hhM2hl &
                     (nint (e_f * 10) ) + 1                             
                     e_graph_hhM2hl (0) = e_graph_hhM2hl (0) + 1 
                  ENDIF 
                  IF (h (1) .eq.0.and.h (3) .eq.0) then 
                     e_graph_0k0 (nint (e_f * 10) ) = e_graph_0k0 (nint &
                     (e_f * 10) ) + 1                                   
                     e_graph_0k0 (0) = e_graph_0k0 (0) + 1 
                  ENDIF 
               ENDIF 
               n_hkl = n_hkl + 1 
               CALL e_accum (e_f, e_e1_hkl, e_e2_hkl, e_e3_hkl,         &
               e_e4_hkl, e_e5_hkl, e_e6_hkl, e_e2m1_hkl, e_e2m1_2_hkl,  &
               e_e2m1_3_hkl, e_abs_e2m1_3_hkl)                          
               IF (h (3) .eq.0) then 
                  n_hk0 = n_hk0 + 1 
                  e_f = sqrt (abs (zz1) ) / e_aver_f2_hk0 (j) 
                  CALL e_accum (e_f, e_e1_hk0, e_e2_hk0, e_e3_hk0,      &
                  e_e4_hk0, e_e5_hk0, e_e6_hk0, e_e2m1_hk0,             &
                  e_e2m1_2_hk0, e_e2m1_3_hk0, e_abs_e2m1_3_hk0)         
               ENDIF 
               IF (h (2) .eq.0) then 
                  n_h0l = n_h0l + 1 
                  e_f = sqrt (abs (zz1) ) / e_aver_f2_h0l (j) 
                  CALL e_accum (e_f, e_e1_h0l, e_e2_h0l, e_e3_h0l,      &
                  e_e4_h0l, e_e5_h0l, e_e6_h0l, e_e2m1_h0l,             &
                  e_e2m1_2_h0l, e_e2m1_3_h0l, e_abs_e2m1_3_h0l)         
               ENDIF 
               IF (h (1) .eq.0) then 
                  n_0kl = n_0kl + 1 
                  e_f = sqrt (abs (zz1) ) / e_aver_f2_0kl (j) 
                  CALL e_accum (e_f, e_e1_0kl, e_e2_0kl, e_e3_0kl,      &
                  e_e4_0kl, e_e5_0kl, e_e6_0kl, e_e2m1_0kl,             &
                  e_e2m1_2_0kl, e_e2m1_3_0kl, e_abs_e2m1_3_0kl)         
               ENDIF 
               IF (h (1) .eq.h (2) .or.h (1) .eq. - h (2) ) then 
                  n_hhl = n_hhl + 1 
                  e_f = sqrt (abs (zz1) ) / e_aver_f2_hhl (j) 
                  CALL e_accum (e_f, e_e1_hhl, e_e2_hhl, e_e3_hhl,      &
                  e_e4_hhl, e_e5_hhl, e_e6_hhl, e_e2m1_hhl,             &
                  e_e2m1_2_hhl, e_e2m1_3_hhl, e_abs_e2m1_3_hhl)         
               ENDIF 
               IF (h (1) .eq. - h (2) ) then 
                  n_hMh0l = n_hMh0l + 1 
                  e_f = sqrt (abs (zz1) ) / e_aver_f2_hMh0l (j) 
                  CALL e_accum (e_f, e_e1_hMh0l, e_e2_hMh0l, e_e3_hMh0l,&
                  e_e4_hMh0l, e_e5_hMh0l, e_e6_hMh0l, e_e2m1_hMh0l,     &
                  e_e2m1_2_hMh0l, e_e2m1_3_hMh0l, e_abs_e2m1_3_hMh0l)   
               ENDIF 
               IF (h (1) .eq.h (2) ) then 
                  n_hhM2hl = n_hhM2hl + 1 
                  e_f = sqrt (abs (zz1) ) / e_aver_f2_hhM2hl (j) 
                  CALL e_accum (e_f, e_e1_hhM2hl, e_e2_hhM2hl,          &
                  e_e3_hhM2hl, e_e4_hhM2hl, e_e5_hhM2hl, e_e6_hhM2hl,   &
                  e_e2m1_hhM2hl, e_e2m1_2_hhM2hl, e_e2m1_3_hhM2hl,      &
                  e_abs_e2m1_3_hhM2hl)                                  
               ENDIF 
               IF (h (1) .eq.0.and.h (3) .eq.0) then 
                  n_0k0 = n_0k0 + 1 
                  e_f = sqrt (abs (zz1) ) / e_aver_f2_0k0 (j) 
                  CALL e_accum (e_f, e_e1_0k0, e_e2_0k0, e_e3_0k0,      &
                  e_e4_0k0, e_e5_0k0, e_e6_0k0, e_e2m1_0k0,             &
                  e_e2m1_2_0k0, e_e2m1_3_0k0, e_abs_e2m1_3_0k0)         
               ENDIF 
            ENDIF 
         ENDIF 
         ENDDO 
!                                                                       
  903    CONTINUE 
         CLOSE (ifa) 
!                                                                       
!     Write results                                                     
!                                                                       
         IF (e_io.ne.6) then 
            outfile = cpara (1:lpara) //'.statistics' 
            CALL oeffne (e_io, outfile, 'unknown') 
         ENDIF 
         lcomm = 3 
         CALL e_write (n_hkl, 'HKL', lcomm, e_e1_hkl, e_e2_hkl,         &
         e_e3_hkl, e_e4_hkl, e_e5_hkl, e_e6_hkl, e_e2m1_hkl,            &
         e_e2m1_2_hkl, e_e2m1_3_hkl, e_abs_e2m1_3_hkl, e_io, .true.)    
         IF (n_0kl.gt.0) then 
            lcomm = 3 
            CALL e_write (n_0kl, '0KL', lcomm, e_e1_0kl, e_e2_0kl,      &
            e_e3_0kl, e_e4_0kl, e_e5_0kl, e_e6_0kl, e_e2m1_0kl,         &
            e_e2m1_2_0kl, e_e2m1_3_0kl, e_abs_e2m1_3_0kl, e_io, .true.) 
         ENDIF 
         IF (n_h0l.gt.0) then 
            lcomm = 3 
            CALL e_write (n_h0l, 'H0L', lcomm, e_e1_h0l, e_e2_h0l,      &
            e_e3_h0l, e_e4_h0l, e_e5_h0l, e_e6_h0l, e_e2m1_h0l,         &
            e_e2m1_2_h0l, e_e2m1_3_h0l, e_abs_e2m1_3_h0l, e_io, .true.) 
         ENDIF 
         IF (n_hk0.gt.0) then 
            lcomm = 3 
            CALL e_write (n_hk0, 'HK0', lcomm, e_e1_hk0, e_e2_hk0,      &
            e_e3_hk0, e_e4_hk0, e_e5_hk0, e_e6_hk0, e_e2m1_hk0,         &
            e_e2m1_2_hk0, e_e2m1_3_hk0, e_abs_e2m1_3_hk0, e_io, .true.) 
         ENDIF 
         IF (n_hhl.gt.0) then 
            lcomm = 9 
            CALL e_write (n_hhl, 'HHL; H-HL', lcomm, e_e1_hhl, e_e2_hhl,&
            e_e3_hhl, e_e4_hhl, e_e5_hhl, e_e6_hhl, e_e2m1_hhl,         &
            e_e2m1_2_hhl, e_e2m1_3_hhl, e_abs_e2m1_3_hhl, e_io, .true.) 
         ENDIF 
         IF (n_hMh0l.gt.0) then 
            lcomm = 3 
            CALL e_write (n_hMh0l, 'H-HL', lcomm, e_e1_hMh0l,           &
            e_e2_hMh0l, e_e3_hMh0l, e_e4_hMh0l, e_e5_hMh0l, e_e6_hMh0l, &
            e_e2m1_hMh0l, e_e2m1_2_hMh0l, e_e2m1_3_hMh0l,               &
            e_abs_e2m1_3_hMh0l, e_io, .true.)                           
         ENDIF 
         IF (n_hhM2hl.gt.0) then 
            lcomm = 3 
            CALL e_write (n_hhM2hl, 'HHL', lcomm, e_e1_hhM2hl,          &
            e_e2_hhM2hl, e_e3_hhM2hl, e_e4_hhM2hl, e_e5_hhM2hl,         &
            e_e6_hhM2hl, e_e2m1_hhM2hl, e_e2m1_2_hhM2hl,                &
            e_e2m1_3_hhM2hl, e_abs_e2m1_3_hhM2hl, e_io, .true.)         
         ENDIF 
         IF (n_0k0.gt.0) then 
            lcomm = 3 
            CALL e_write (n_0k0, '0K0', lcomm, e_e1_0k0, e_e2_0k0,      &
            e_e3_0k0, e_e4_0k0, e_e5_0k0, e_e6_0k0, e_e2m1_0k0,         &
            e_e2m1_2_0k0, e_e2m1_3_0k0, e_abs_e2m1_3_0k0, e_io, .true.) 
         ENDIF 
         IF (e_io.ne.6) then 
            CLOSE (e_io) 
         ENDIF 
!                                                                       
!     Write short results                                               
!                                                                       
         IF (e_io.ne.6) then 
            outfile = cpara (1:lpara) //'.short' 
            CALL oeffne (e_io, outfile, 'unknown') 
         ENDIF 
         lcomm = 3 
         CALL e_write (n_hkl, 'HKL', lcomm, e_e1_hkl, e_e2_hkl,         &
         e_e3_hkl, e_e4_hkl, e_e5_hkl, e_e6_hkl, e_e2m1_hkl,            &
         e_e2m1_2_hkl, e_e2m1_3_hkl, e_abs_e2m1_3_hkl, e_io, .false.)   
         IF (n_0kl.gt.0) then 
            lcomm = 3 
            CALL e_write (n_0kl, '0KL', lcomm, e_e1_0kl, e_e2_0kl,      &
            e_e3_0kl, e_e4_0kl, e_e5_0kl, e_e6_0kl, e_e2m1_0kl,         &
            e_e2m1_2_0kl, e_e2m1_3_0kl, e_abs_e2m1_3_0kl, e_io, .false.)
         ENDIF 
         IF (n_h0l.gt.0) then 
            lcomm = 3 
            CALL e_write (n_h0l, 'H0L', lcomm, e_e1_h0l, e_e2_h0l,      &
            e_e3_h0l, e_e4_h0l, e_e5_h0l, e_e6_h0l, e_e2m1_h0l,         &
            e_e2m1_2_h0l, e_e2m1_3_h0l, e_abs_e2m1_3_h0l, e_io, .false.)
         ENDIF 
         IF (n_hk0.gt.0) then 
            lcomm = 3 
            CALL e_write (n_hk0, 'HK0', lcomm, e_e1_hk0, e_e2_hk0,      &
            e_e3_hk0, e_e4_hk0, e_e5_hk0, e_e6_hk0, e_e2m1_hk0,         &
            e_e2m1_2_hk0, e_e2m1_3_hk0, e_abs_e2m1_3_hk0, e_io, .false.)
         ENDIF 
         IF (n_hhl.gt.0) then 
            lcomm = 9 
            CALL e_write (n_hhl, 'HHL; H-HL', lcomm, e_e1_hhl, e_e2_hhl,&
            e_e3_hhl, e_e4_hhl, e_e5_hhl, e_e6_hhl, e_e2m1_hhl,         &
            e_e2m1_2_hhl, e_e2m1_3_hhl, e_abs_e2m1_3_hhl, e_io, .false.)
         ENDIF 
         IF (n_hMh0l.gt.0) then 
            lcomm = 3 
            CALL e_write (n_hMh0l, 'H-HL', lcomm, e_e1_hMh0l,           &
            e_e2_hMh0l, e_e3_hMh0l, e_e4_hMh0l, e_e5_hMh0l, e_e6_hMh0l, &
            e_e2m1_hMh0l, e_e2m1_2_hMh0l, e_e2m1_3_hMh0l,               &
            e_abs_e2m1_3_hMh0l, e_io, .false.)                          
         ENDIF 
         IF (n_hhM2hl.gt.0) then 
            lcomm = 3 
            CALL e_write (n_hhM2hl, 'HHL', lcomm, e_e1_hhM2hl,          &
            e_e2_hhM2hl, e_e3_hhM2hl, e_e4_hhM2hl, e_e5_hhM2hl,         &
            e_e6_hhM2hl, e_e2m1_hhM2hl, e_e2m1_2_hhM2hl,                &
            e_e2m1_3_hhM2hl, e_abs_e2m1_3_hhM2hl, e_io, .false.)        
         ENDIF 
         IF (n_0k0.gt.0) then 
            lcomm = 3 
            CALL e_write (n_0k0, '0K0', lcomm, e_e1_0k0, e_e2_0k0,      &
            e_e3_0k0, e_e4_0k0, e_e5_0k0, e_e6_0k0, e_e2m1_0k0,         &
            e_e2m1_2_0k0, e_e2m1_3_0k0, e_abs_e2m1_3_0k0, e_io, .false.)
         ENDIF 
         IF (e_io.ne.6) then 
            CLOSE (e_io) 
         ENDIF 
      ENDIF 
!                                                                       
      IF (io_mode) then 
         IF (e_io.ne.6) then 
            outfile = cpara (1:lpara) //'.hkl.histogram' 
            CALL oeffne (e_io, outfile, 'unknown') 
         ENDIF 
         DO i = 1, 30 
         WRITE (e_io, 1000) i * 0.1 - 0.05, float (e_graph (i) )        &
         / float (e_graph (0) ) / 0.1                                   
         ENDDO 
         IF (e_io.ne.6) then 
            CLOSE (e_io) 
         ENDIF 
         IF (e_io.ne.6) then 
            outfile = cpara (1:lpara) //'.hk0.histogram' 
            CALL oeffne (e_io, outfile, 'unknown') 
         ENDIF 
         DO i = 1, 30 
         WRITE (e_io, 1000) i * 0.1 - 0.05, float (e_graph_hk0 (i) )    &
         / float (e_graph_hk0 (0) ) / 0.1                               
         ENDDO 
         IF (e_io.ne.6) then 
            CLOSE (e_io) 
         ENDIF 
         IF (e_io.ne.6) then 
            outfile = cpara (1:lpara) //'.h0l.histogram' 
            CALL oeffne (e_io, outfile, 'unknown') 
         ENDIF 
         DO i = 1, 30 
         WRITE (e_io, 1000) i * 0.1 - 0.05, float (e_graph_h0l (i) )    &
         / float (e_graph_h0l (0) ) / 0.1                               
         ENDDO 
         IF (e_io.ne.6) then 
            CLOSE (e_io) 
         ENDIF 
         IF (e_io.ne.6) then 
            outfile = cpara (1:lpara) //'.0kl.histogram' 
            CALL oeffne (e_io, outfile, 'unknown') 
         ENDIF 
         DO i = 1, 30 
         WRITE (e_io, 1000) i * 0.1 - 0.05, float (e_graph_0kl (i) )    &
         / float (e_graph_0kl (0) ) / 0.1                               
         ENDDO 
         IF (e_io.ne.6) then 
            CLOSE (e_io) 
         ENDIF 
         IF (e_io.ne.6) then 
            outfile = cpara (1:lpara) //'.hhl.histogram' 
            CALL oeffne (e_io, outfile, 'unknown') 
         ENDIF 
         DO i = 1, 30 
         WRITE (e_io, 1000) i * 0.1 - 0.05, float (e_graph_hhl (i) )    &
         / float (e_graph_hhl (0) ) / 0.1                               
         ENDDO 
         IF (e_io.ne.6) then 
            CLOSE (e_io) 
         ENDIF 
         IF (e_io.ne.6) then 
            outfile = cpara (1:lpara) //'.hMh0l.histogram' 
            CALL oeffne (e_io, outfile, 'unknown') 
         ENDIF 
         DO i = 1, 30 
         WRITE (e_io, 1000) i * 0.1 - 0.05, float (e_graph_hMh0l (i) )  &
         / float (e_graph_hMh0l (0) ) / 0.1                             
         ENDDO 
         IF (e_io.ne.6) then 
            CLOSE (e_io) 
         ENDIF 
         IF (e_io.ne.6) then 
            outfile = cpara (1:lpara) //'.hhM2Hl.histogram' 
            CALL oeffne (e_io, outfile, 'unknown') 
         ENDIF 
         DO i = 1, 30 
         WRITE (e_io, 1000) i * 0.1 - 0.05, float (e_graph_hhM2Hl (i) ) &
         / float (e_graph_hhM2Hl (0) ) / 0.1                            
         ENDDO 
         IF (e_io.ne.6) then 
            CLOSE (e_io) 
         ENDIF 
         IF (e_io.ne.6) then 
            outfile = cpara (1:lpara) //'.0k0.histogram' 
            CALL oeffne (e_io, outfile, 'unknown') 
         ENDIF 
         DO i = 1, 30 
         WRITE (e_io, 1000) i * 0.1 - 0.05, float (e_graph_0k0 (i) )    &
         / float (e_graph_0k0 (0) ) / 0.1                               
         ENDDO 
         IF (e_io.ne.6) then 
            CLOSE (e_io) 
            CLOSE (e_io_hkl) 
         ENDIF 
      ENDIF 
!                                                                       
  950 CONTINUE 
      CLOSE (ifa) 
      CLOSE (55) 
      CLOSE (56) 
!                                                                       
 1000 FORMAT    (f8.3,2x,f14.8) 
 4000 FORMAT    (3i4,2f8.2) 
!                                                                       
      END SUBROUTINE e_create                       
!*****7*****************************************************************
      INTEGER function e_hist (h) 
!-                                                                      
!     Determines the slot for dstar                                     
!                                                                       
!     Version  : 1.0                                                    
!     Date     : 27 Nov 2000                                            
!                                                                       
!     Author   : R.B. Neder (reinhard.neder@fau.de)      
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE metric_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      REAL h (3) 
!                                                                       
      INTEGER i 
      REAL null (3) 
      REAL e_hist_limit (20) 
      REAL dstar 
      LOGICAL lspace 
      PARAMETER (lspace = .false.) 
!                                                                       
!     REAL do_blen 
!                                                                       
      DATA null / 0.0, 0.0, 0.0 / 
      DATA e_hist_limit / 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8,  &
      2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 1e6 /           
!                                                                       
      dstar = do_blen (lspace, h, null) 
!                                                                       
      DO i = 1, 20 
      IF (dstar.lt.e_hist_limit (i) ) then 
         e_hist = i 
         RETURN 
      ENDIF 
      ENDDO 
      e_hist = 20 
!                                                                       
      END FUNCTION e_hist                           
!*****7*****************************************************************
      SUBROUTINE e_accum (e_f, e_e1_hkl, e_e2_hkl, e_e3_hkl, e_e4_hkl,  &
      e_e5_hkl, e_e6_hkl, e_e2m1_hkl, e_e2m1_2_hkl, e_e2m1_3_hkl,       &
      e_abs_e2m1_3_hkl)                                                 
!-                                                                      
!     Accumulates the E statistics                                      
!                                                                       
!     Version  : 1.0                                                    
!     Date     : 27 Nov 2000                                            
!                                                                       
!     Author   : R.B. Neder (reinhard.neder@fau.de)      
!+                                                                      
      IMPLICIT none 
!                                                                       
      REAL e_f 
      REAL e_e1_hkl, e_e2_hkl, e_e3_hkl, e_e4_hkl 
      REAL e_e5_hkl, e_e6_hkl, e_e2m1_hkl 
      REAL e_e2m1_2_hkl, e_e2m1_3_hkl, e_abs_e2m1_3_hkl 
!                                                                       
      e_e1_hkl = e_e1_hkl + e_f 
      e_e2_hkl = e_e2_hkl + e_f**2 
      e_e3_hkl = e_e3_hkl + e_f**3 
      e_e4_hkl = e_e4_hkl + e_f**4 
      e_e5_hkl = e_e5_hkl + e_f**5 
      e_e6_hkl = e_e6_hkl + e_f**6 
      e_e2m1_hkl = e_e2m1_hkl + abs (e_f**2 - 1.0) 
      e_e2m1_2_hkl = e_e2m1_2_hkl + (e_f**2 - 1.0) **2 
      e_e2m1_3_hkl = e_e2m1_3_hkl + (e_f**2 - 1.0) **3 
      e_abs_e2m1_3_hkl = e_abs_e2m1_3_hkl + (abs (e_f**2 - 1.0) ) **3 
!                                                                       
      END SUBROUTINE e_accum                        
!*****7*****************************************************************
      SUBROUTINE e_write (n, e_name, lcomm, e_e1_hkl, e_e2_hkl,         &
      e_e3_hkl, e_e4_hkl, e_e5_hkl, e_e6_hkl, e_e2m1_hkl, e_e2m1_2_hkl, &
      e_e2m1_3_hkl, e_abs_e2m1_3_hkl, e_io, long)                       
!-                                                                      
!     write the E statistics to the output device                       
!                                                                       
!     Version  : 1.0                                                    
!     Date     : 27 Nov 2000                                            
!                                                                       
!     Author   : R.B. Neder (reinhard.neder@fau.de)      
!+                                                                      
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      CHARACTER(6) zeile 
      CHARACTER ( * ) e_name 
      INTEGER e_io 
      INTEGER lcomm 
      INTEGER n 
      LOGICAL long 
      REAL e_e1_hkl, e_e2_hkl, e_e3_hkl, e_e4_hkl 
      REAL e_e5_hkl, e_e6_hkl, e_e2m1_hkl 
      REAL e_e2m1_2_hkl, e_e2m1_3_hkl, e_abs_e2m1_3_hkl 
!                                                                       
      REAL cent (10), acent (10) 
!                                                                       
      DATA cent / 0.798, 1.000, 1.596, 3.000, 6.383, 15.000, 0.968,     &
      2.000, 8.000, 8.691 /                                             
      DATA acent / 0.886, 1.000, 1.329, 2.000, 3.323, 6.000, 0.736,     &
      1.000, 2.000, 2.415 /                                             
!                                                                       
!                                                                       
!     Write results                                                     
!                                                                       
      WRITE (e_io, 1000) 'Number of ', e_name (1:lcomm) , ' reflections:&
     &  ', n                                                            
      IF (long) then 
      WRITE (e_io, 1010) '<|E|**1>        : ', e_e1_hkl / n, cent (1) , &
     &acent (1)                                                         
      WRITE (e_io, 1010) '<|E|**2>        : ', e_e2_hkl / n, cent (2) , &
     &acent (2)                                                         
      WRITE (e_io, 1010) '<|E|**3>        : ', e_e3_hkl / n, cent (3) , &
     &acent (3)                                                         
      WRITE (e_io, 1010) '<|E|**4>        : ', e_e4_hkl / n, cent (4) , &
     &acent (4)                                                         
      WRITE (e_io, 1010) '<|E|**5>        : ', e_e5_hkl / n, cent (5) , &
     &acent (5)                                                         
      WRITE (e_io, 1010) '<|E|**6>        : ', e_e6_hkl / n, cent (6) , &
     &acent (6)                                                         
      WRITE (e_io, 1010) '<|E**2 - 1|>    : ', e_e2m1_hkl / n, cent (7) &
     &, acent (7)                                                       
         WRITE (e_io, 1010) '<(E**2 - 1)**2> : ', e_e2m1_2_hkl / n,     &
         cent (8) , acent (8)                                           
         WRITE (e_io, 1010) '<(E**2 - 1)**3> : ', e_e2m1_3_hkl / n,     &
         cent (9) , acent (9)                                           
         WRITE (e_io, 1010) '<|E**2 - 1|**3> : ', e_abs_e2m1_3_hkl / n, &
         cent (10) , acent (10)                                         
      ELSE 
      WRITE (e_io, 1010) '<|E**2 - 1|>    : ', e_e2m1_hkl / n, cent (7) &
     &, acent (7)                                                       
      ENDIF 
      IF (e_io.eq.6) then 
         zeile = 'return' 
         lcomm = 6 
         CALL do_input (zeile, lcomm) 
      ENDIF 
!                                                                       
 1000 FORMAT    (/,3a,i6,/,                                             &
     &           ' Criterion',11x,'observed,  centric, acentric')       
 1010 FORMAT    (a,3f10.3) 
!                                                                       
      END SUBROUTINE e_write                        
!*****7*****************************************************************
      SUBROUTINE wilson_calc 
!-                                                                      
!     Calculates a Wilson plot                                          
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE diffuse_mod 
      USE fourier_sup
      USE patters_mod 
      USE wink_mod
      USE quad_mod
!
      USE errlist_mod 
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER w_io 
!                                                                       
      INTEGER MAXW 
      PARAMETER (MAXW = 40) 
      CHARACTER(1024) outfile 
      INTEGER ifa, idot 
      INTEGER i
      INTEGER h (3) 
      INTEGER w_num (MAXW) 
      REAL hh (3) 
      REAL rint, sigma 
      REAL dstar2, stl2, fjq 
      REAL dd, dd_o 
      REAL q 
      REAL w_aver (MAXW) 
      REAL sumx, sumx2, sumy, sumy2, sumxy 
      REAL m, p, number 
!                                                                       
      INTEGER HKLF4
      PARAMETER (HKLF4 = 5) 
!                                                                       
      INTEGER len_str 
!     INTEGER e_hist 
!     REAL fj2 
!     REAL quad 
      INTEGER iimax 
      iimax = 0 
!                                                                       
      ifa = 33 
      w_io = 36 
!                                                                       
      IF (cr_nscat.eq.0) then 
         ier_num = - 98 
         ier_typ = ER_APPL 
         ier_msg (1) = ' Wilson statistics can only be calculated,' 
         ier_msg (2) = ' if atom types corresponding to the ' 
         ier_msg (3) = ' chemical formula have been defined' 
         RETURN 
      ENDIF 
      IF (lambda.eq.' '.and.rlambda.eq.0.0) then 
         ier_num = - 99 
         ier_typ = ER_APPL 
         ier_msg (1) = ' Wilson statistics can only be calculated,' 
         ier_msg (2) = ' if the radiation and wave length have' 
         ier_msg (3) = ' been defined' 
         RETURN 
      ENDIF 
      CALL dlink (ano, lambda, rlambda, renergy, l_energy, &
                  diff_radiation, diff_power) 
      DO i = 1, MAXW 
      w_aver (i) = 0.0 
      w_num (i) = 0
      ENDDO 
!                                                                       
      IF (ftyp.ne.HKLF4) then 
         ier_num = - 999 
         ier_typ = ER_APPL 
         RETURN 
      ELSE 
         CALL oeffne (ifa, rho_file (1) , 'old') 
!                                                                       
!     Open output file                                                  
!                                                                       
         IF (w_io.ne.6) then 
            idot = index (rho_file (1) , '.hkl') 
            IF (idot.gt.0) then 
               outfile = rho_file (1) (1:idot) //'wilson' 
            ELSE 
               outfile = rho_file (1) (1:len_str (rho_file (1) ) ) //   &
               '.wilson'                                                
            ENDIF 
            CALL oeffne (w_io, outfile, 'unknown') 
         ENDIF 
!                                                                       
!------ --Read HKLF4 file, accumulate the Intensities in the            
!         sin(theta)/lambda                                             
!       slots                                                           
!                                                                       
   10    CONTINUE 
         READ (ifa, *, err = 902, end = 902) h, rint, sigma 
         IF (.not. (h (1) .eq.0.and.h (2) .eq.0.and.h (3) .eq.0) ) then 
            hh (1) = h (1) 
            hh (2) = h (2) 
            hh (3) = h (3) 
            i = int (0.25 * quad (hh, hh, cr_rten) * MAXW / 2.) + 1 
            iimax = max (iimax, i) 
            IF (i.lt.MAXW + 1) then 
               w_aver (i) = w_aver (i) + rint 
               w_num (i) = w_num (i) + 1 
            ENDIF 
            GOTO 10 
!                                                                       
         ENDIF 
!                                                                       
  902    CONTINUE 
!                                                                       
         CLOSE (ifa) 
!                                                                       
!     --Normalize the accumulated values, calculate the Scale factor and
!       overall atomic displacement parameter                           
!                                                                       
         sumx = 0.0 
         sumx2 = 0.0 
         sumy = 0.0 
         sumy2 = 0.0 
         sumxy = 0.0 
         number = 0 
         dd = 2. / MAXW * 4. 
         dd_o = 1. / MAXW * 4. 
         DO i = 1, MAXW 
         IF (w_num (i) .gt.0.and.w_aver (i) .gt.0) then 
            w_aver (i) = w_aver (i) / w_num (i) 
            dstar2 = (i - 1) * dd+dd_o 
            stl2 = 0.25 * dstar2 
            fjq = fj2 (dstar2) 
            q = log (fjq / abs (w_aver (i) ) ) 
            sumx = sumx + stl2 
            sumx2 = sumx2 + stl2**2 
            sumy = sumy + q 
            sumy2 = sumy2 + q**2 
            sumxy = sumxy + stl2 * q 
            number = number + 1 
            WRITE (w_io, 1000) stl2, q 
         ENDIF 
         ENDDO 
!                                                                       
         m = number * sumx2 - sumx**2 
         p = number * sumxy - sumx * sumy 
         wilson_scale = exp (0.5 * (m * sumy - p * sumx) / number / m) 
!DBG        wilson_scale = exp(    (m*sumy - p*sumx)/number/m)          
         wilson_b = 0.5 * p / m 
!                                                                       
      ENDIF 
!                                                                       
      IF (w_io.ne.6) then 
         CLOSE (w_io) 
      ENDIF 
!                                                                       
      WRITE (output_io, * ) 
      WRITE (output_io,  * ) 'Scale: Fcalc = ', sqrt (wilson_scale) , ' &
     &* Fobs'                                                           
      WRITE (output_io, * ) 'Scale: Icalc = ', wilson_scale, ' * Iobs' 
      WRITE (output_io,  * ) 'overall  B   = ', wilson_b 
      WRITE (output_io,  * ) 'overall  U   = ', wilson_b / 8. / pi**2 
      WRITE (output_io, * ) 
!                                                                       
 1000 FORMAT    (f8.5,2x,g17.8e3) 
!                                                                       
      END SUBROUTINE wilson_calc                    
!*****7*****************************************************************
      REAL function fj2 (dstar2) 
!-                                                                      
!     Calculates the sum of all formfactors squared                     
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE diffuse_mod 
      USE fourier_sup, ONLY: form
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      REAL dstar2 
      INTEGER j 
!                                                                       
!     REAL form 
!                                                                       
      fj2 = 0.0 
      DO j = 1, cr_ncatoms 
      fj2 = fj2 + (form (cr_iscat (j), cr_scat, lxray, dstar2, diff_power) ) **2 
      ENDDO 
!                                                                       
      END FUNCTION fj2                              
END MODULE patters_menu
