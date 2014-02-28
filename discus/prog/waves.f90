MODULE waves_do_menu
!
CONTAINS
SUBROUTINE waves_menu
!-                                                                      
!     calculates the displacement of the atoms due to a plane wave      
!     travelling through the crystal                                    
!+                                                                      
      USE config_mod 
      USE allocate_appl_mod
      USE crystal_mod 
      USE modify_mod
      USE molecule_mod 
      USE show_menu
      USE waves_mod 
!
      USE doact_mod 
      USE errlist_mod 
      USE learn_mod 
      USE macro_mod 
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER, PARAMETER :: MIN_PARA = 26  ! A command requires at least these no of parameters
      INTEGER maxw 
      LOGICAL lnew, lold 
!                                                                       
      PARAMETER (lnew = .true., lold = .false.) 
!                                                                       
      CHARACTER(LEN=1024), DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: cpara
      REAL               , DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: werte
      INTEGER            , DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: lpara
!
      CHARACTER(1024) line, zeile
      CHARACTER(1024) cdummy 
      CHARACTER(50) prom 
      CHARACTER(5) befehl 
      INTEGER lp, length, lbef, ldummy 
      INTEGER indxg, ianz, is 
      INTEGER         :: nscat
      LOGICAL lend 
!                                                                       
      INTEGER len_str 
      LOGICAL str_comp 
!                                                                       
      maxw = MAX(MIN_PARA,MAXSCAT+1)
      lend = .false. 
      CALL no_error 
!
      IF( cr_nscat > WV_MAXSCAT .or. mole_num_type > WV_MAXSCAT) THEN
         nscat = max ( cr_nscat, mole_num_type)
         CALL alloc_waves ( cr_nscat )
         IF ( ier_num < 0 ) THEN
            RETURN
         ENDIF
      ENDIF
!                                                                       
      DO while (.not.lend) 
      prom = prompt (1:len_str (prompt) ) //'/wave' 
      CALL get_cmd (line, length, befehl, lbef, zeile, lp, prom) 
      IF (ier_num.eq.0) then 
         IF (line.ne.' '.and.line (1:1) .ne.'#') then 
!                                                                       
!     --- search for "="                                                
!                                                                       
            indxg = index (line, '=') 
      IF (indxg.ne.0.and..not. (str_comp (befehl, 'echo', 2, lbef, 4) ) &
     &.and..not. (str_comp (befehl, 'syst', 2, lbef, 4) ) .and..not. (st&
     &r_comp (befehl, 'help', 2, lbef, 4) .or.str_comp (befehl, '?   ', &
     &2, lbef, 4) ) ) then                                              
!                                                                       
!     ----- evaluate an expression and assign value a variabble         
!                                                                       
               CALL do_math (line, indxg, length) 
            ELSE 
!                                                                       
!------ --- execute a macro file                                        
!                                                                       
               IF (befehl (1:1) .eq.'@') then 
                  IF (length.ge.2) then 
                     CALL file_kdo (line (2:length), length - 1) 
                  ELSE 
                     ier_num = - 13 
                     ier_typ = ER_MAC 
                  ENDIF 
!                                                                       
!     --- accoustic wave 'acco'                                         
!                                                                       
               ELSEIF (str_comp (befehl, 'acco', 2, lbef, 4) ) then 
                  wv_lacoust = .true. 
!                                                                       
!     --- amplitude of wave 'ampl'                                      
!                                                                       
               ELSEIF (str_comp (befehl, 'ampl', 2, lbef, 4) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (ianz.eq.1) then 
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxw)                                           
                        IF (ier_num.eq.0) then 
                           wv_amp = werte (1) 
                        ENDIF 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ENDIF 
!                                                                       
!     --- list asymmetric unit 'asym'                                   
!                                                                       
               ELSEIF (str_comp (befehl, 'asym', 2, lbef, 4) ) then 
                  CALL show_asym 
!                                                                       
!     --- list atoms present in the crystal 'chem'                      
!                                                                       
               ELSEIF (str_comp (befehl, 'chem', 2, lbef, 4) ) then 
                  CALL show_chem 
!                                                                       
!     --- continues a macro 'continue'                                  
!                                                                       
               ELSEIF (str_comp (befehl, 'continue', 2, lbef, 8) ) then 
                  CALL macro_continue (zeile, lp) 
!                                                                       
!     --- density wave 'dens'                                           
!                                                                       
               ELSEIF (str_comp (befehl, 'dens', 3, lbef, 4) ) then 
                  wv_iwave = WV_DENS 
!                                                                       
!     --- select/deselect molecules                                     
!                                                                       
               ELSEIF (str_comp (befehl, 'msel', 2, lbef, 4)            &
               .or.str_comp (befehl, 'mdes', 2, lbef, 4) ) then         
!                                                                       
                  CALL mole_select (zeile, lp, 0, WV_MAXSCAT, wv_latom, &
                  wv_sel_atom,         &
                  str_comp (befehl, 'msel', 2, lbef, 4) )               
!                                                                       
!     --- select which atoms are copied to their image 'sele'           
!                                                                       
               ELSEIF (str_comp (befehl, 'sele', 2, lbef, 4)            &
               .or.str_comp (befehl, 'dese', 2, lbef, 4) ) then         
!                                                                       
                  CALL atom_select (zeile, lp, 0, WV_MAXSCAT, wv_latom, &
                  wv_sel_atom, lold,        &
                  str_comp (befehl, 'sele', 2, lbef, 4) )               
!                                                                       
!------ --- Echo a string, just for interactive check in a macro 'echo' 
!                                                                       
               ELSEIF (str_comp (befehl, 'echo', 2, lbef, 4) ) then 
                  CALL echo (zeile, lp) 
!                                                                       
!      -- Evaluate an expression, just for interactive check 'eval'     
!                                                                       
               ELSEIF (str_comp (befehl, 'eval', 2, lbef, 4) ) then 
                  CALL do_eval (zeile, lp) 
!                                                                       
!     --- exit 'exit'                                                   
!                                                                       
               ELSEIF (str_comp (befehl, 'exit', 2, lbef, 4) ) then 
                  lend = .true. 
!                                                                       
!     --- wave function 'func'                                          
!------ --- possible functions are 'box', 'sinusoidal', 'triangular'    
!                                                                       
               ELSEIF (str_comp (befehl, 'func', 1, lbef, 4) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (ianz.eq.1.or.ianz.eq.2) then 
                        IF (str_comp (cpara (1) , 'box', 3, lpara (1) , &
                        3) ) then                                       
                           wv_func = cpara (1) (1:lpara(1))
                           wv_ifunc = WV_BOX 
                           IF (ianz.eq.2) then 
                              CALL del_params (1, ianz, cpara, lpara,   &
                              maxw)                                     
                              CALL ber_params (ianz, cpara, lpara,      &
                              werte, maxw)                              
                              IF (werte (1) .gt.0.0.and.werte (1)       &
                              .lt.1.0) then                             
                                 wv_asym = werte (1) 
                              ELSE 
                                 ier_num = - 70 
                                 ier_typ = ER_APPL 
                              ENDIF 
                           ELSE 
                              wv_asym = 0.5 
                           ENDIF 
                        ELSEIF (str_comp (cpara (1) , 'sinus', 3, lpara &
                        (1) , 5) ) then                                 
                           wv_func = cpara (1) (1:lpara(1))
                           wv_ifunc = WV_SINUS 
                        ELSEIF (str_comp (cpara (1) , 'trian', 3, lpara &
                        (1) , 5) ) then                                 
                           wv_func = cpara (1) (1:lpara(1))
                           wv_ifunc = WV_TRIANGLE 
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
!     --- help 'help','?'                                               
!                                                                       
      ELSEIF (str_comp (befehl, 'help', 2, lbef, 4) .or.str_comp (befehl&
     &, '?   ', 1, lbef, 4) ) then                                      
                  IF (str_comp (zeile, 'errors', 2, lp, 6) ) then 
                     lp = lp + 7 
                     CALL do_hel ('discus '//zeile, lp) 
                  ELSE 
                     lp = lp + 12 
                     CALL do_hel ('discus wave '//zeile, lp) 
                  ENDIF 
!                                                                       
!     --- wave length 'leng'                                            
!                                                                       
               ELSEIF (str_comp (befehl, 'leng', 2, lbef, 4) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (ianz.eq.1) then 
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxw)                                           
                        IF (ier_num.eq.0) then 
                           wv_rlam = werte (1) 
                        ENDIF 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ENDIF 
!                                                                       
!     --- longitudinal wave 'long'                                      
!                                                                       
               ELSEIF (str_comp (befehl, 'long', 2, lbef, 4) ) then 
                  wv_iwave = WV_LONG 
!                                                                       
!     ----optical wave 'opti'                                           
!                                                                       
               ELSEIF (str_comp (befehl, 'opti', 2, lbef, 4) ) then 
                  wv_lacoust = .false. 
!                                                                       
!     --- direction of oscillation 'osci'                               
!                                                                       
               ELSEIF (str_comp (befehl, 'osci', 2, lbef, 4) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (ianz.eq.3) then 
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxw)                                           
                        IF (ier_num.eq.0) then 
                           wv_swing (1) = werte (1) 
                           wv_swing (2) = werte (2) 
                           wv_swing (3) = werte (3) 
                        ENDIF 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ENDIF 
!                                                                       
!     --- initial phase 'phas'                                          
!                                                                       
               ELSEIF (str_comp (befehl, 'phas', 3, lbef, 4) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (ianz.eq.1) then 
                        IF (str_comp (cpara (1) , 'rand', 1, lpara (1) ,&
                        4) ) then                                       
                           wv_phase_typ = WV_RAND 
                        ELSE 
                           wv_phase_typ = WV_FIX 
                           CALL ber_params (ianz, cpara, lpara, werte,  &
                           maxw)                                        
                           IF (ier_num.eq.0) then 
                              wv_phase = werte (1) 
                           ENDIF 
                        ENDIF 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ENDIF 
!                                                                       
!     --- upper limit of probability for density waves 'phig'           
!                                                                       
               ELSEIF (str_comp (befehl, 'phig', 3, lbef, 4) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (ianz.eq.1) then 
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxw)                                           
                        IF (ier_num.eq.0) then 
                           wv_phigh = min (1.0, abs (werte (1) ) ) 
                        ENDIF 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ENDIF 
!                                                                       
!     --- lower limit of probability for density waves 'plow'           
!                                                                       
               ELSEIF (str_comp (befehl, 'plow', 2, lbef, 4) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (ianz.eq.1) then 
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxw)                                           
                        IF (ier_num.eq.0) then 
                           wv_plow = min (1.0, abs (werte (1) ) ) 
                        ENDIF 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ENDIF 
!                                                                       
!     --- Select which molecules are repl. by density wave 'mrepl'      
!                                                                       
               ELSEIF (str_comp (befehl, 'mrepl', 2, lbef, 5) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF(ianz==3) THEN
                        wv_viceversa=str_comp (cpara(3),'viceversa',2,lpara(3),9)
                        cpara(3) = ' '
                        lpara(3) = 1
                        ianz     = 2
                     ELSE
                        wv_viceversa=.false.
                     ENDIF
                     IF (ianz.eq.2) then 
                        cdummy = cpara (1) 
                        ldummy = lpara (1) 
                        CALL del_params (1, ianz, cpara, lpara, maxw) 
                        CALL ber_params (ianz, cpara, lpara, werte,maxw)
                        is = nint(werte(1))
                        IF (is.ge.0.and.is.le.mole_num_type) then 
                           is = nint (werte (1) ) 
                           CALL mole_select (cdummy, ldummy, 0, WV_MAXSCAT,&
                                wv_latom, wv_sel_atom, .true.,    &
                                is, wv_repl)             
                        ELSE 
                           ier_num = - 64 
                           ier_typ = ER_APPL 
                        ENDIF 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ENDIF 
!                                                                       
!     --- Select which atoms are repl. by density wave 'repl'           
!                                                                       
               ELSEIF (str_comp (befehl, 'repl', 2, lbef, 4) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF(ianz==3) THEN
                        wv_viceversa=str_comp (cpara(3),'viceversa',2,lpara(3),9)
                        cpara(3) = ' '
                        lpara(3) = 1
                        ianz     = 2
                     ELSE
                        wv_viceversa=.false.
                     ENDIF
                     IF (ianz.eq.2) then 
                        cdummy = cpara (1) (1:lpara (1) ) 
                        ldummy = lpara (1) 
                        CALL del_params (1, ianz, cpara, lpara, maxw) 
                        CALL get_iscat (ianz, cpara, lpara, werte, maxw,&
                        lnew)                                           
                        IF (ier_num.eq.0.and.is.ne. - 1) then 
                           is = nint (werte (1) ) 
                           CALL atom_select (cdummy, ldummy, 0, MAXSCAT,&
                           wv_latom,  &
                           wv_sel_atom,        lold, .true.,            &
                           is, wv_repl)
                        ENDIF 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ENDIF 
!                                                                       
!     --- select rotation mode for molecules                            
!                                                                       
               ELSEIF (str_comp (befehl, 'rot', 2, lbef, 3) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (ianz.eq.3.or.ianz.eq.6) then 
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxw)                                           
                        IF (ier_num.eq.0) then 
                           wv_iwave = WV_ROT 
                           wv_rot_uvw (1) = werte (1) 
                           wv_rot_uvw (2) = werte (2) 
                           wv_rot_uvw (3) = werte (3) 
                           IF (ianz.eq.6) then 
                              wv_rot_orig (1) = werte (4) 
                              wv_rot_orig (2) = werte (5) 
                              wv_rot_orig (3) = werte (6) 
                           ELSE 
                              wv_rot_orig (1) = 0.0 
                              wv_rot_orig (2) = 0.0 
                              wv_rot_orig (3) = 0.0 
                           ENDIF 
                        ENDIF 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ENDIF 
!                                                                       
!     --- run waves 'run'                                               
!                                                                       
               ELSEIF (str_comp (befehl, 'run ', 2, lbef, 4) ) then 
                  CALL wave_run 
!                                                                       
!     --- set constant shift 'shif'                                     
!                                                                       
               ELSEIF (str_comp (befehl, 'shif', 3, lbef, 4) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (ianz.eq.1) then 
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxw)                                           
                        IF (ier_num.eq.0) then 
                           wv_amp0 = werte (1) 
                        ENDIF 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ENDIF 
!                                                                       
!     --- show current parameters 'show'                                
!                                                                       
               ELSEIF (str_comp (befehl, 'show', 3, lbef, 4) ) then 
                  CALL wave_show 
!                                                                       
!     --- transversal wave 'tran'                                       
!                                                                       
               ELSEIF (str_comp (befehl, 'tran', 1, lbef, 4) ) then 
                  wv_iwave = WV_TRANS 
!                                                                       
!     --- direction of wave vector 'vect'                               
!                                                                       
               ELSEIF (str_comp (befehl, 'vect', 1, lbef, 4) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (ianz.eq.3) then 
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxw)                                           
                        IF (ier_num.eq.0) then 
                           wv_wave (1) = werte (1) 
                           wv_wave (2) = werte (2) 
                           wv_wave (3) = werte (3) 
                        ENDIF 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ENDIF 
!                                                                       
!------ --- operating System Kommandos 'syst'                           
!                                                                       
               ELSEIF (str_comp (befehl, 'syst', 2, lbef, 4) ) then 
                  IF (zeile.ne.' ') then 
                     CALL do_operating (zeile (1:lp), lp) 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
!                                                                       
!------ --- waiting for user input                                      
!                                                                       
               ELSEIF (str_comp (befehl, 'wait', 3, lbef, 4) ) then 
                  CALL do_input (zeile, lp) 
!                                                                       
!------ --- End of command processing                                   
!                                                                       
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_COMM 
               ENDIF 
            ENDIF 
         ENDIF 
      ENDIF 
!                                                                       
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
      ENDDO 
!                                                                       
      END SUBROUTINE waves_menu
!*****7*********************************************************        
      SUBROUTINE wave_show 
!-                                                                      
!     Shows current settiungs in WAVES segment                          
!+                                                                      
      USE config_mod 
      USE crystal_mod 
      USE atom_name
      USE molecule_mod 
      USE waves_mod 
      USE wink_mod
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      CHARACTER(15) cwave (3) 
      CHARACTER(9) at_lis (maxscat)!, at_name 
      CHARACTER(6) cphase (2) 
      INTEGER mol_lis (maxscat) 
      INTEGER i, j, k 
!                                                                       
      DATA cphase / 'random', 'fixed ' / 
      DATA cwave / 'sinusoidal wave', 'box shaped wave', 'triangular wav&
     &e' /                                                              
!                                                                       
      WRITE (output_io, 1000) 
!                                                                       
!------ Density wave selected                                           
!                                                                       
      IF (wv_iwave.eq.WV_DENS) then 
         WRITE (output_io, 2000) 'density wave' 
         IF (wv_ifunc.eq.WV_BOX) then 
            WRITE (output_io, 2100) cwave (wv_ifunc), wv_asym 
         ELSE 
            WRITE (output_io, 2150) cwave (wv_ifunc) 
         ENDIF 
         WRITE (output_io, 2200) (wv_wave (i), i = 1, 3) 
         WRITE (output_io, 2300) wv_rlam, wv_plow, wv_phigh, wv_phase,  &
         cphase (wv_phase_typ)                                          
!                                                                       
         IF (wv_sel_atom) then 
            WRITE (output_io, 2600) 
            DO i = 0, cr_nscat 
            IF (wv_latom (i) ) then 
               at_lis (1) = at_name (i) 
               at_lis (2) = at_name (wv_repl (i) ) 
               WRITE (output_io, 2650) i, at_lis (1), wv_repl (i),      &
               at_lis (2)                                               
            ENDIF 
            ENDDO 
         ELSE 
            WRITE (output_io, 2700) 
            DO i = 1, mole_num_type 
            WRITE (output_io, 2750) i, wv_repl (i) 
            ENDDO 
         ENDIF 
!                                                                       
!------ Displacement or rotational wave                                 
!                                                                       
      ELSE 
         IF (wv_lacoust) then 
            IF (wv_iwave.eq.WV_TRANS) then 
               WRITE (output_io, 2000) 'transversal acoustic wave' 
            ELSEIF (wv_iwave.eq.WV_LONG) then 
               WRITE (output_io, 2000) 'longitudinal acoustic wave' 
            ELSEIF (wv_iwave.eq.WV_ROT) then 
               WRITE (output_io, 2000) 'rotational wave (molecules)' 
            ENDIF 
         ELSE 
            IF (wv_iwave.eq.WV_TRANS) then 
               WRITE (output_io, 2000) 'transversal optical wave' 
            ELSEIF (wv_iwave.eq.WV_LONG) then 
               WRITE (output_io, 2000) 'longitudinal optical wave' 
            ELSEIF (wv_iwave.eq.WV_ROT) then 
               WRITE (output_io, 2000) 'rotational wave (molecules)' 
            ENDIF 
         ENDIF 
         IF (wv_ifunc.eq.WV_BOX) then 
            WRITE (output_io, 2100) cwave (wv_ifunc), wv_asym 
         ELSE 
            WRITE (output_io, 2150) cwave (wv_ifunc) 
         ENDIF 
         WRITE (output_io, 2200) (wv_wave (i), i = 1, 3) 
         IF (wv_iwave.eq.WV_TRANS) then 
            WRITE (output_io, 3100) (wv_swing (i), i = 1, 3) 
         ENDIF 
         IF (wv_iwave.eq.WV_ROT) then 
            WRITE (output_io, 3120) (wv_rot_uvw (i), i = 1, 3) 
            WRITE (output_io, 3140) (wv_rot_orig (i), i = 1, 3) 
         ENDIF 
         WRITE (output_io, 3200) wv_rlam, wv_amp, wv_amp0, wv_phase,    &
         cphase (wv_phase_typ)                                          
!                                                                       
         IF (wv_sel_atom) then 
            j = 0 
            DO i = 0, cr_nscat 
            IF (wv_latom (i) ) then 
               j = j + 1 
               at_lis (j) = at_name (i) 
            ENDIF 
            ENDDO 
            WRITE (output_io, 3300) (at_lis (k), k = 1, j) 
         ELSE 
            j = 0 
            DO i = 0, mole_num_type 
            IF (wv_latom (i) ) then 
               j = j + 1 
               mol_lis (j) = i 
            ENDIF 
            ENDDO 
            WRITE (output_io, 3400) (mol_lis (k), k = 1, j) 
         ENDIF 
      ENDIF 
!                                                                       
 1000 FORMAT     (  ' Current settings for wave module : ') 
!                                                                       
 2000 FORMAT     (  '   Wave type           : ',2x,a) 
 2100 FORMAT     (  '   Wave function       : ',2x,a,' - Delta : ',f4.2) 
 2150 FORMAT     (  '   Wave function       : ',2x,a) 
 2200 FORMAT     (  '   Wave vector         : ',3(2x,f9.4)) 
 2300 FORMAT     (  '   Wave length         : ',2x,f9.4/                &
     &                     '   Low  probabilty     : ',2x,f9.4/         &
     &                     '   High probabilty     : ',2x,f9.4/         &
     &                     '   Phase at 0,0,0      : ',2x,f9.4/         &
     &                     '   Phase status        : ',2x,a)            
 2600 FORMAT     (  '   Selected atoms      : ') 
 2650 FORMAT     (  '     Atom typ ',i3,' (',a9,') to be replaced ',    &
     &                 'by atom type ',i3,' (',a9,')')                  
 2700 FORMAT     (  '   Selected molecules  : ') 
 2750 FORMAT     (  '     Molecule typ ',i3,' to be replaced ',         &
     &                 'by molecule type ',i3)                          
!                                                                       
 3100 FORMAT     (  '   Oscillation vector  : ',3(2x,f9.4)) 
 3120 FORMAT     (  '   Rotation axis       : ',3(2x,f9.4)) 
 3140 FORMAT     (  '   Rotation origin     : ',3(2x,f9.4)) 
 3200 FORMAT     (  '   Wave length         : ',2x,f9.4/                &
     &                     '   Amplitude           : ',2x,f9.4/         &
     &                     '   Constant shift      : ',2x,f9.4/         &
     &                     '   Phase at 0,0,0      : ',2x,f9.4/         &
     &                     '   Phase status        : ',2x,a   )         
 3300 FORMAT     (  '   Sel. atom types     : ',2x,50(a9,1x)) 
 3400 FORMAT     (  '   Sel. molecule types : ',2x,50(i4,1x)) 
!                                                                       
      END SUBROUTINE wave_show                      
!*****7*********************************************************        
      SUBROUTINE wave_run 
!-                                                                      
!     Main waves routine                                                
!+                                                                      
      USE config_mod 
      USE allocate_appl_mod 
      USE crystal_mod 
      USE quad_mod
      USE symm_mod 
      USE update_cr_dim_mod
      USE trafo_mod
      USE waves_mod 
      USE errlist_mod 
      USE wink_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      REAL wavec (3), swingc (3) 
      REAL wavep (3), swingp (3) 
      REAL ww, rnn 
      INTEGER i, lbef 
!                                                                       
      REAL hkl (3), uvw (3), tran (3), orig (3), angle
      INTEGER start, end, power
      LOGICAL pmult, mode, new, orig_mol, typ, sel_atom
!                                                                       
!     REAL quad 
      INTEGER len_str 
      LOGICAL str_comp 
      angle = 0.0
      start = 1
      end = 1
      power = 1
      pmult=.true.
      mode=.true.
      new=.false.
      orig_mol=.true.
      typ=.true.
      sel_atom=.true.
!                                                                       
!      REAL sinus, box, triang 
!      EXTERNAL sinus, box, triang 
!                                                                       
!------ Setup for rotational mode                                       
!                                                                       
      IF (wv_iwave.eq.WV_ROT) then 
         IF (wv_sel_atom) then 
            ier_num = - 68 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
!
!--------As sym_latom is used, we need to allocate symmetry
!
         IF( cr_nscat > SYM_MAXSCAT) THEN
            CALL alloc_symmetry ( cr_nscat )
            IF ( ier_num < 0 ) THEN
               RETURN
            ENDIF
         ENDIF
!                                                                       
         DO i = 1, 3 
            hkl (i)  = sym_hkl (i) 
            uvw (i)  = sym_uvw (i) 
            tran (i) = sym_trans (i) 
            orig (i) = sym_orig (i) 
         ENDDO 
         angle = sym_angle 
         start = sym_start 
         end   = sym_end 
         power = sym_power 
         mode  = sym_mode 
         new   = sym_new 
         orig_mol = sym_orig_mol 
         sel_atom = sym_sel_atom 
         typ      = sym_type 
         pmult    = sym_power_mult 
         DO i = 0, WV_MAXSCAT 
            wv_latom (i) = sym_latom (i) 
         ENDDO 
!                                                                       
         sym_power    = 1 
         sym_power_mult = .false. 
         sym_type     = .true. 
         sym_new      = .false. 
         sym_mode     = .false. 
         sym_orig_mol = .true. 
         sym_sel_atom = .false. 
         DO i = 1, 3 
            sym_uvw (i)   = wv_rot_uvw (i) 
            sym_orig (i)  = wv_rot_orig (i) 
            sym_trans (i) = 0.0 
         ENDDO 
!                                                                       
         CALL trans (sym_uvw, cr_gten, sym_hkl, 3) 
      ENDIF 
!                                                                       
!     normalise wavevector to length 1 A,                               
!     transform to normalised space                                     
!                                                                       
      ww = sqrt (quad (wv_wave, wv_wave, cr_gten) ) 
      DO i = 1, 3 
         wavec (i) = wv_wave (i) / ww 
      ENDDO 
      CALL trans (wavec, cr_fmat, wavep, 3) 
!                                                                       
!     For density waves set correct values of amp and amp0 depending    
!     on wave function type                                             
!                                                                       
      IF (wv_iwave.eq.WV_DENS) then 
         lbef = len_str (wv_func) 
         IF (str_comp (wv_func, 'box', 1, lbef, 3) ) then 
            wv_amp = wv_phigh - wv_plow 
            wv_amp0 = wv_plow 
         ELSEIF (str_comp (wv_func, 'sinus', 1, lbef, 5) ) then 
            wv_amp = (wv_phigh - wv_plow) * 0.5 
            wv_amp0 = (wv_phigh + wv_plow) * 0.5 
         ELSEIF (str_comp (wv_func, 'trian', 1, lbef, 5) ) then 
            wv_amp = wv_phigh - wv_plow 
            wv_amp0 = wv_plow 
         ENDIF 
      ELSE 
!                                                                       
!     normalise oscillation vector to lenght 1 A,                       
!     transform to normalised space                                     
!                                                                       
         IF (wv_iwave.eq.WV_TRANS) then 
            rnn = sqrt (quad (wv_swing, wv_swing, cr_gten) ) 
            DO i = 1, 3 
            swingc (i) = wv_swing (i) / rnn 
            ENDDO 
         ELSE 
            DO i = 1, 3 
            swingc (i) = wavec (i) 
            ENDDO 
         ENDIF 
      ENDIF 
!                                                                       
      CALL trans (swingc, cr_fmat, swingp, 3) 
!                                                                       
!------ use all of the crystal (outside)                                
!                                                                       
      IF (wv_sel_atom) then 
         IF (wv_func (1:1) .eq.'b') then 
            CALL wave_run_all (wavep, swingp, box) 
         ELSEIF (wv_func (1:1) .eq.'s') then 
            CALL wave_run_all (wavep, swingp, sinus) 
         ELSEIF (wv_func (1:1) .eq.'t') then 
            CALL wave_run_all (wavep, swingp, triang) 
         ENDIF 
      ELSE 
         IF (wv_func (1:1) .eq.'b') then 
            CALL wave_run_all_mol (wavep, swingp, box) 
         ELSEIF (wv_func (1:1) .eq.'s') then 
            CALL wave_run_all_mol (wavep, swingp, sinus) 
         ELSEIF (wv_func (1:1) .eq.'t') then 
            CALL wave_run_all_mol (wavep, swingp, triang) 
         ENDIF 
      ENDIF 
!                                                                       
      CALL update_cr_dim 
!                                                                       
!------ Restore symm settings                                           
!                                                                       
      IF (wv_iwave.eq.WV_ROT) then 
         DO i = 1, 3 
            sym_hkl (i)   = hkl (i) 
            sym_uvw (i)   = uvw (i) 
            sym_trans (i) = tran (i) 
            sym_orig (i)  = orig (i) 
         ENDDO 
         sym_angle = angle 
         sym_start = start 
         sym_end   = end 
         sym_power = power 
         sym_mode  = mode 
         sym_new   = new 
         sym_orig_mol  = orig_mol 
         sym_type      = typ 
         sym_sel_atom  = sel_atom 
         sym_power_mult = pmult 
         DO i = 0, WV_MAXSCAT 
            sym_latom (i) = wv_latom (i) 
         ENDDO 
      ENDIF 
!                                                                       
      END SUBROUTINE wave_run                       
!*****7*********************************************************        
      SUBROUTINE wave_run_all (wavep, swingp, wave_func) 
!-                                                                      
!     Displaces all atoms by the corresponding                          
!     wave function.                                                    
!+                                                                      
      USE config_mod 
      USE crystal_mod 
      USE trafo_mod
      USE waves_mod 
      USE wink_mod
      USE random_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      REAL uc (3), up (3), disp (3), disc (3) 
      REAL wavep (3), swingp (3), arg, dis 
      INTEGER i, j 
      INTEGER  :: is_target   ! dummy for target atom type
!                                                                       
      REAL ran1 
!                                                                       
      REAL wave_func 
      EXTERNAL wave_func 
!                                                                       
      IF (wv_phase_typ.eq.WV_RAND) then 
         wv_phase = ran1 (idum) * 360.0 
      ENDIF 
!                                                                       
is_density: IF (wv_iwave.eq.WV_DENS) then 
   DO i = 1, cr_natoms 
!                                                                       
!------ - Check if atom is a valid selection                            
!                                                                       
      IF (wv_latom (cr_iscat (i) ) ) then 
!                                                                       
         DO j = 1, 3 
            uc (j) = cr_pos (j, i) 
         ENDDO 
!                                                                       
         CALL trans (uc, cr_fmat, up, 3) 
         arg = up(1) * wavep(1) + up(2) * wavep(2) + up(3) * wavep(3)
         arg = arg + wv_phase * wv_rlam / 360. 
         arg = amod (arg, wv_rlam) / wv_rlam 
         dis = wave_func (wv_amp, arg, wv_amp0) 
!                                                                       
         IF (ran1 (idum) .gt.dis) then 
            cr_iscat (i) = wv_repl (cr_iscat (i) ) 
         ENDIF 
      ELSEIF(wv_viceversa) THEN       ! Is atom the target of repl command ?
         is_target = -1       
         find_target: DO j=0,cr_nscat
            IF( cr_iscat(i) == wv_repl(j)) THEN
               is_target = j
               EXIT find_target
            ENDIF
         ENDDO find_target
         IF( is_target >= 0 ) THEN     ! Atom is target, replace by source
            DO j = 1, 3 
               uc (j) = cr_pos (j, i) 
            ENDDO 
!                                                                       
            CALL trans (uc, cr_fmat, up, 3) 
            arg = up(1) * wavep(1) + up(2) * wavep(2) + up(3) * wavep(3)
            arg = arg + wv_phase * wv_rlam / 360. 
            arg = amod (arg, wv_rlam) / wv_rlam 
            dis = 1.0 - wave_func (wv_amp, arg, wv_amp0) 
!                                                                       
            IF (ran1 (idum) .gt.dis) then 
               cr_iscat (i) = is_target
            ENDIF 
         ENDIF
      ENDIF
   ENDDO
ELSE is_density
!                                                                       
      DO i = 1, cr_natoms 
!                                                                       
!------ - Check if atom is a valid selection                            
!                                                                       
      IF (wv_latom (cr_iscat (i) ) ) then 
!                                                                       
         DO j = 1, 3 
            uc (j) = cr_pos (j, i) 
         ENDDO 
!                                                                       
         CALL trans (uc, cr_fmat, up, 3) 
         arg = up (1) * wavep (1) + up (2) * wavep (2) + up (3) * wavep &
         (3)                                                            
         arg = arg + wv_phase * wv_rlam / 360. 
         arg = amod (arg, wv_rlam) / wv_rlam 
         dis = wave_func (wv_amp, arg, wv_amp0) 
!                                                                       
         IF (.not.wv_lacoust.and. (index (cr_at_lis (cr_iscat (i) ) ,&
            '-') .gt.0) ) then                                          
               dis = - dis 
         ENDIF 
!                                                                       
         DO j = 1, 3 
            disp (j) = dis * swingp (j) 
         ENDDO 
         CALL trans (disp, cr_gmat, disc, 3) 
         DO j = 1, 3 
            cr_pos (j, i) = cr_pos (j, i) + disc (j) 
         ENDDO 
      ENDIF 
   ENDDO 
ENDIF is_density
!                                                                       
      END SUBROUTINE wave_run_all                   
!*****7*********************************************************        
      SUBROUTINE wave_run_all_mol (wavep, swingp, wave_func) 
!-                                                                      
!     Displaces all molecules by the corresponding                      
!     wave function.                                                    
!+                                                                      
      USE config_mod 
      USE crystal_mod 
      USE prop_para_mod 
      USE modify_mod
      USE molecule_mod 
      USE symm_mod 
      USE symm_sup_mod
      USE trafo_mod
      USE waves_mod 
      USE wink_mod
      USE random_mod
      IMPLICIT none 
!                                                                       
      REAL uc (3), up (3), disp (3), disc (3) 
      REAL wavep (3), swingp (3), arg, dis 
      INTEGER i, j, ia, ityp, im, il 
      INTEGER  :: is_repl  ! no of a target molecule for density waves
      INTEGER  :: is_src   ! no of a source molecule for density waves
!                                                                       
      REAL ran1 
!                                                                       
      REAL wave_func 
      EXTERNAL wave_func 
!                                                                       
      IF (wv_phase_typ.eq.WV_RAND) then 
         wv_phase = ran1 (idum) * 360.0 
      ENDIF 
!                                                                       
!                                                                       
!------- ---- Density wave                                              
!                                                                       
wave_type: IF (wv_iwave.eq.WV_DENS) then 
   dens_loop:   DO i = 1, mole_num_mole 
      ia   = mole_cont (mole_off (i) + 1) 
      ityp = mole_type (i) 
!                                                                       
!------ - Check if molecule is a valid selection                        
!                                                                       
      is_source: IF (wv_latom (mole_type (i) ) ) then  ! Selected molecule to be replaced
!                                                                       
         DO j = 1, 3 
            uc (j) = cr_pos (j, ia) 
         ENDDO 
!                                                                       
         CALL trans (uc, cr_fmat, up, 3) 
         arg = up(1) * wavep(1) + up(2) * wavep(2) + up(3) * wavep(3)                                                            
         arg = arg + wv_phase * wv_rlam / 360. 
         arg = amod (arg, wv_rlam) / wv_rlam 
         dis = wave_func (wv_amp, arg, wv_amp0) 
         IF (ran1 (idum) .gt.dis) then 
            IF (wv_repl (ityp) .ne.0) then 
               is_repl = i                       ! Replace by itself if no target found
find_source:   DO is_repl=1,mole_num_mole
                  IF(mole_type(is_repl) == wv_repl(ityp)) THEN
                     EXIT find_source
                  ENDIF
               ENDDO find_source
               CALL do_swap_mole (i, is_repl, .false.) 
!!!                  CALL do_swap_mole (ityp, wv_repl (ityp), .false.) 
            ELSE 
               mole_type (i) = 0 
               DO j = 1, mole_len (i) 
                  cr_iscat (mole_cont (mole_off (i) + j) ) = 0 
                  cr_prop (mole_cont (mole_off (i) + j) ) = ibclr (     &
                  cr_prop (mole_cont (mole_off (i) + j) ), PROP_NORMAL) 
               ENDDO 
            ENDIF 
         ENDIF
      ELSEIF(wv_viceversa) THEN is_source ! Selected molecule is target?
         is_src = 0
         find_target: DO is_repl = 1, mole_num_type
            IF(ityp == wv_repl(is_repl)) THEN
               is_src = is_repl
               EXIT find_target
            ENDIF
         ENDDO find_target
         is_target: IF ( is_src /=0) THEN
!
         DO j = 1, 3 
            uc (j) = cr_pos (j, ia) 
         ENDDO 
!                                                                       
         CALL trans (uc, cr_fmat, up, 3) 
         arg = up(1) * wavep(1) + up(2) * wavep(2) + up(3) * wavep(3)                                                            
         arg = arg + wv_phase * wv_rlam / 360. 
         arg = amod (arg, wv_rlam) / wv_rlam 
         dis = 1.0 - wave_func (wv_amp, arg, wv_amp0)  ! Replace by opposite probability
         IF (ran1 (idum) .gt.dis) then 
            ityp = mole_type (i) 
            IF (is_src .ne.0) then 
find_src:   DO is_repl=1,mole_num_mole
                  IF(mole_type(is_repl) == is_src) THEN
                     EXIT find_src
                  ENDIF
               ENDDO find_src
               CALL do_swap_mole (i, is_repl, .false.) 
!!!                  CALL do_swap_mole (ityp, wv_repl (ityp), .false.) 
            ELSE 
               mole_type (i) = 0 
               DO j = 1, mole_len (i) 
                  cr_iscat (mole_cont (mole_off (i) + j) ) = 0 
                  cr_prop (mole_cont (mole_off (i) + j) ) = ibclr (     &
                  cr_prop (mole_cont (mole_off (i) + j) ), PROP_NORMAL) 
               ENDDO 
            ENDIF 
         ENDIF
         ENDIF is_target
      ENDIF is_source
   ENDDO dens_loop
ELSE wave_type
      DO i = 1, mole_num_mole 
         ia = mole_cont (mole_off (i) + 1) 
!                                                                       
!------ - Check if molecule is a valid selection                        
!                                                                       
      IF (wv_latom (mole_type (i) ) ) then 
!                                                                       
         DO j = 1, 3 
            uc (j) = cr_pos (j, ia) 
         ENDDO 
!                                                                       
         CALL trans (uc, cr_fmat, up, 3) 
         arg = up (1) * wavep (1) + up (2) * wavep (2) + up (3) * wavep &
         (3)                                                            
         arg = arg + wv_phase * wv_rlam / 360. 
         arg = amod (arg, wv_rlam) / wv_rlam 
         dis = wave_func (wv_amp, arg, wv_amp0) 
!                                                                       
!------- ---- Rotational wave                                           
!                                                                       
         IF (wv_iwave.eq.WV_ROT) then 
            sym_latom (mole_type (i) ) = .true. 
            sym_angle = dis 
            sym_start = i 
            sym_end = i 
            CALL symm_setup 
            CALL symm_mole_single 
!                                                                       
!------- ---- Displacement wave                                         
!                                                                       
         ELSE 
            IF (.not.wv_lacoust.and. (index (cr_at_lis (cr_iscat (ia) ) &
            , '-') .gt.0) ) then                                        
               dis = - dis 
            ENDIF 
!                                                                       
            DO j = 1, 3 
            disp (j) = dis * swingp (j) 
            ENDDO 
            CALL trans (disp, cr_gmat, disc, 3) 
            DO il = 1, mole_len (i) 
            im = mole_cont (mole_off (i) + il) 
            DO j = 1, 3 
            cr_pos (j, im) = cr_pos (j, im) + disc (j) 
            ENDDO 
            ENDDO 
         ENDIF 
         ENDIF 
      ENDDO 
      ENDIF  wave_type
!                                                                       
      END SUBROUTINE wave_run_all_mol               
!*****7*********************************************************        
      REAL function sinus (amp, arg, amp0) 
!-                                                                      
!     Sinusoidal wave function                                          
!+                                                                      
      USE wink_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      REAL amp, arg, amp0 
!                                                                       
      arg = arg * zpi 
      sinus = amp * cos (arg) + amp0 
!                                                                       
      END FUNCTION sinus                            
!*****7*********************************************************        
      REAL function box (amp, arg, amp0) 
!-                                                                      
!     box shaped wave function                                          
!+                                                                      
      USE config_mod 
      USE waves_mod 
      USE wink_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      REAL amp, arg, amp0 
!                                                                       
!     We have a positive argument 'arg'.                                
!                                                                       
      IF (arg.ge.0.0) then 
         IF (arg.lt. (0.5 * wv_asym) .or.arg.ge. (1.0 - 0.5 * wv_asym) )&
         then                                                           
            box = amp0 
         ELSE 
            box = amp0 + amp 
         ENDIF 
!                                                                       
!     For negative argument the exclusions are reversed                 
!------ to ensure continuity at arg=0.0                                 
!                                                                       
      ELSE 
         IF ( - arg.le. (0.5 * wv_asym) .or. - arg.gt. (1.0 - 0.5 *     &
         wv_asym) ) then                                                
            box = amp0 
         ELSE 
            box = amp0 + amp 
         ENDIF 
      ENDIF 
!                                                                       
      END FUNCTION box                              
!*****7*********************************************************        
      REAL function triang (amp, arg, amp0) 
!-                                                                      
!     triangular shaped wave function                                   
!+                                                                      
      USE wink_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      REAL amp, arg, amp0 
!                                                                       
      IF (arg.ge.0) then 
         triang = amp0 + amp * arg 
      ELSE 
         triang = amp0 + amp * (1. + arg) 
      ENDIF 
!                                                                       
      END FUNCTION triang                           
END MODULE waves_do_menu
