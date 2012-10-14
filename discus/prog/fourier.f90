      SUBROUTINE fourier 
!+                                                                      
!     This subroutine 'fourier' calculates the Fourier transform        
!     of the given crystal structure. The algorithm to speed up         
!     the explicite Fourier is based on the program 'DIFFUSE' by        
!     B.D. Butler. See also: B.D. Butler & T.R. Welberry, (1992).       
!     J. Appl. Cryst. 25, 391-399.                                      
!                                                                       
!-                                                                      
      USE config_mod 
      USE allocate_appl_mod
      USE crystal_mod 
      USE diffuse_mod 
      USE modify_mod
      USE output_mod 
      IMPLICIT none 
!                                                                       
       
      include'doact.inc' 
      include'learn.inc' 
      include'macro.inc' 
      include'prompt.inc' 
      include'errlist.inc' 
!                                                                       
      INTEGER, PARAMETER :: MIN_PARA = 21  ! A command requires at leaset these no of parameters
      INTEGER maxw 
      LOGICAL lold 
      PARAMETER (lold = .false.) 
!                                                                       
      CHARACTER (LEN=1024), DIMENSION(MAX(MIN_PARA,MAXSCAT+1))   :: cpara ! (MIN(10,MAXSCAT)) 
      INTEGER             , DIMENSION(MAX(MIN_PARA,MAXSCAT))   :: lpara ! (MIN(10,MAXSCAT))
      INTEGER             , DIMENSION(MAX(MIN_PARA,MAXSCAT))   :: jj    ! (MAXSCAT) 
      REAL                , DIMENSION(MAX(MIN_PARA,MAXSCAT))   :: werte ! (MAXSCAT)
      CHARACTER(5) befehl 
      CHARACTER(50) prom 
      CHARACTER(1024) zeile
      CHARACTER(1024) line 
      INTEGER i, j, k, ianz, lp, length 
      INTEGER indxg, lbef 
      INTEGER              :: n_qxy    ! required size in reciprocal space this run
      INTEGER              :: n_nscat  ! required no of atom types right now
      INTEGER              :: n_natoms ! required no of atoms
      LOGICAL ldim 
      REAL divis (2) 
      REAL rhkl (3) 
!                                                                       
      INTEGER len_str 
      LOGICAL str_comp 
!
      maxw     = MAX(MIN_PARA,MAXSCAT+1)
      n_qxy    = 1
      n_nscat  = 1
      n_natoms = 1
!                                                                       
   10 CONTINUE 
!                                                                       
      CALL no_error 
      divis (1) = float (max (1, inc (1) - 1) ) 
      divis (2) = float (max (1, inc (2) - 1) ) 
!                                                                       
      prom = prompt (1:len_str (prompt) ) //'/fourier' 
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
!     calculate at a single reciprocal point 'calc'                     
!                                                                       
            ELSEIF (str_comp (befehl, 'calc', 2, lbef, 4) ) then 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ier_num.eq.0) then 
                  IF (ianz.eq.3) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        rhkl (1) = werte (1) 
                        rhkl (2) = werte (2) 
                        rhkl (3) = werte (3) 
                        IF (inc (1) * inc (2) .gt. MAXQXY  .OR.          &
                            cr_natoms > DIF_MAXAT          .OR.          &
                            cr_nscat>DIF_MAXSCAT              ) THEN
                          n_qxy    = MAX(n_qxy,inc(1) * inc(2),MAXQXY)
                          n_natoms = MAX(n_natoms,cr_natoms,DIF_MAXAT)
                          n_nscat  = MAX(n_nscat,cr_nscat,DIF_MAXSCAT)
                          call alloc_diffuse (n_qxy, n_nscat, n_natoms)
                          IF (ier_num.ne.0) THEN
                            RETURN
                          ENDIF
                        ENDIF
                        CALL dlink (lxray, ano, lambda, rlambda) 
                        CALL calc_000 (rhkl) 
                     ENDIF 
                  ELSEIF (ianz.eq.0) then 
                     rhkl (1) = 0.0 
                     rhkl (2) = 0.0 
                     rhkl (3) = 0.0 
                     CALL dlink (lxray, ano, lambda, rlambda) 
                     CALL calc_000 (rhkl) 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ENDIF 
!                                                                       
!     continues a macro 'continue'                                      
!                                                                       
            ELSEIF (str_comp (befehl, 'continue', 1, lbef, 8) ) then 
               CALL macro_continue (zeile, lp) 
!                                                                       
!     define the anomalous scattering curve for an element 'delf'       
!                                                                       
            ELSEIF (str_comp (befehl, 'delf', 2, lbef, 4) ) then 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ianz.eq.3) then 
                  i = 1 
                  CALL get_iscat (i, cpara, lpara, werte, maxw, lold) 
                  IF (ier_num.eq.0) then 
                     IF (werte (1) .gt.0) then 
                        DO k = 1, i 
                        jj (k) = nint (werte (1) ) 
                        ENDDO 
                        cpara (1) = '0.0' 
                        lpara (1) = 3 
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxw)                                           
                        IF (ier_num.eq.0) then 
                           DO k = 1, i 
                           cr_delfr (jj (k) ) = werte (2) 
                           cr_delfi (jj (k) ) = werte (3) 
                           cr_delf_int (jj (k) ) = .false. 
                           ENDDO 
                        ENDIF 
                     ELSE 
                        ier_num = - 27 
                        ier_typ = ER_APPL 
                     ENDIF 
                  ENDIF 
               ELSEIF (ianz.eq.2) then 
                  IF (str_comp (cpara (2) , 'internal', 2, lpara (2) ,  &
                  8) ) then                                             
                     CALL get_iscat (i, cpara, lpara, werte, maxw, lold) 
                     IF (ier_num.eq.0) then 
                        i = nint (werte (1) ) 
                        IF (i.gt.0) then 
                           cr_delf_int (i) = .true. 
                        ELSEIF (i.eq. - 1) then 
                           DO k = 1, cr_nscat 
                           cr_delf_int (k) = .true. 
                           ENDDO 
                        ELSE 
                           ier_num = - 6 
                           ier_typ = ER_COMM 
                        ENDIF 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
!                                                                       
!     Switch dispersion on/off 'disp' second parameter 'anom' or 'off'  
!                                                                       
            ELSEIF (str_comp (befehl, 'disp', 2, lbef, 4) ) then 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ier_num.eq.0) then 
                  IF (ianz.eq.1) then 
                     IF (cpara (1) (1:1) .eq.'a') then 
                        ano = .true. 
                     ELSEIF (cpara (1) (1:1) .eq.'o') then 
                        ano = .false. 
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
!     Terminate Fourier 'exit'                                          
!                                                                       
            ELSEIF (str_comp (befehl, 'exit', 3, lbef, 4) ) then 
               GOTO 9999 
!                                                                       
!     help 'help' , '?'                                                 
!                                                                       
      ELSEIF (str_comp (befehl, 'help', 2, lbef, 4) .or.str_comp (befehl&
     &, '?   ', 1, lbef, 4) ) then                                      
               IF (str_comp (zeile, 'errors', 2, lp, 6) ) then 
                  lp = lp + 7 
                  CALL do_hel ('discus '//zeile, lp) 
               ELSE 
                  lp = lp + 12 
                  CALL do_hel ('discus four '//zeile, lp) 
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
!     define the corners 'll', 'lr', 'ul'                               
!                                                                       
      ELSEIF (str_comp (befehl, 'll  ', 2, lbef, 4) .or.str_comp (befehl&
     &, 'lr  ', 2, lbef, 4) .or.str_comp (befehl, 'ul  ', 2, lbef, 4) ) &
     &then                                                              
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ier_num.eq.0) then 
                  IF (ianz.eq.3) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
      IF (str_comp (befehl, 'll  ', 2, lbef, 4) ) j = 1 
      IF (str_comp (befehl, 'lr  ', 2, lbef, 4) ) j = 2 
      IF (str_comp (befehl, 'ul  ', 2, lbef, 4) ) j = 3 
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
!     set Fourier volume (lots)                                         
!                                                                       
            ELSEIF (str_comp (befehl, 'lots', 3, lbef, 4) ) then 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ier_num.eq.0) then 
                  IF (ianz.eq.1) then 
                     CALL do_cap (cpara (1) ) 
                     IF (cpara (1) (1:1) .eq.'O') then 
                        ilots = LOT_OFF 
                        nlots = 1 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ELSEIF (ianz.eq.6) then 
                     CALL do_cap (cpara (1) ) 
                     IF (cpara (1) (1:1) .eq.'B') then 
                        ilots = LOT_BOX 
                     ELSEIF (cpara (1) (1:1) .eq.'E') then 
                        ilots = LOT_ELI 
                     ELSE 
                        ier_num = - 2 
                        ier_typ = ER_FOUR 
                     ENDIF 
                     IF (ier_num.eq.0) then 
                        CALL del_params (1, ianz, cpara, lpara, maxw) 
                        CALL ber_params (ianz - 1, cpara, lpara, werte, &
                        maxw)                                           
                        IF (ier_num.eq.0) then 
                           ldim = .true. 
                           DO i = 1, 3 
                           ldim = ldim.and. (0..lt.nint (werte (i) )    &
                           .and.nint (werte (i) ) .le.cr_icc (i) )      
                           ENDDO 
                           IF (ldim) then 
                              ls_xyz (1) = nint (werte (1) ) 
                              ls_xyz (2) = nint (werte (2) ) 
                              ls_xyz (3) = nint (werte (3) ) 
                              nlots = nint (werte (4) ) 
                              CALL do_cap (cpara (5) ) 
                              lperiod = (cpara (5) (1:1) .eq.'Y') 
                           ELSE 
                              ier_num = - 101 
                              ier_typ = ER_APPL 
                           ENDIF 
                        ENDIF 
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
!                                                                       
!     switch to neutron diffraction 'neut'                              
!                                                                       
            ELSEIF (str_comp (befehl, 'neut', 2, lbef, 4) ) then 
               lxray = .false. 
!                                                                       
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
!     start the Fourier transform 'run'                                 
!                                                                       
            ELSEIF (str_comp (befehl, 'run ', 1, lbef, 4) ) then 
               IF (inc (1) * inc (2) .gt. MAXQXY  .OR.          &
                   cr_natoms > DIF_MAXAT          .OR.          &
                   cr_nscat>DIF_MAXSCAT              ) THEN
                 n_qxy    = MAX(n_qxy,inc(1) * inc(2),MAXQXY)
                 n_natoms = MAX(n_natoms,cr_natoms,DIF_MAXAT)
                 n_nscat  = MAX(n_nscat,cr_nscat,DIF_MAXSCAT)
                 call alloc_diffuse (n_qxy, n_nscat, n_natoms)
                 IF (ier_num.ne.0) THEN
                   RETURN
                 ENDIF
               ENDIF
               IF (inc (1) * inc (2) .le.MAXQXY) then 
                  CALL dlink (lxray, ano, lambda, rlambda) 
                  IF (four_mode.eq.INTERNAL) then 
                     IF (ier_num.eq.0) then 
                        four_log = .true. 
                        CALL four_run 
                     ENDIF 
                  ELSE 
                     four_log = .true. 
                     CALL four_external 
                  ENDIF 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_APPL 
               ENDIF 
!                                                                       
!     define the scattering curve for an element 'scat'                 
!                                                                       
            ELSEIF (str_comp (befehl, 'scat', 2, lbef, 4) ) then 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ianz.eq.10) then 
                  i = 1 
                  CALL get_iscat (i, cpara, lpara, werte, maxw, lold) 
                  IF (ier_num.eq.0) then 
                     IF (werte (1) .gt.0) then 
                        DO k = 1, i 
                        jj (k) = nint (werte (1) ) 
                        ENDDO 
                        cpara (1) = '0.0' 
                        lpara (1) = 3 
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxw)                                           
                        IF (ier_num.eq.0) then 
                           DO k = 1, i 
                           DO i = 2, 9 
                           cr_scat (i, jj (k) ) = werte (i) 
                           ENDDO 
                           cr_scat (1, jj (k) ) = werte (10) 
                           cr_scat_int (jj (k) ) = .false. 
                           ENDDO 
                        ENDIF 
                     ELSE 
                        ier_num = - 27 
                        ier_typ = ER_APPL 
                     ENDIF 
                  ENDIF 
               ELSEIF (ianz.eq.2) then 
                  i = 1 
                  CALL get_iscat (i, cpara, lpara, werte, maxw, lold) 
                  IF (ier_num.eq.0) then 
                     IF (str_comp (cpara (2) , 'internal', 2, lpara (1) &
                     , 8) ) then                                        
                        IF (werte (1) .gt.0) then 
                           k = nint (werte (1) ) 
                           cr_scat_int (k) = .true. 
                           cr_scat_equ (k) = .false. 
                        ELSEIF (werte (1) .eq. - 1) then 
                           DO k = 0, cr_nscat 
                           cr_scat_int (k) = .true. 
                           cr_scat_equ (k) = .false. 
                           ENDDO 
                        ENDIF 
                     ELSE 
                        IF (werte (1) .gt.0) then 
                           k = nint (werte (1) ) 
                           cr_scat_equ (k) = .true. 
                           CALL do_cap (cpara (2) ) 
                           cr_at_equ (k) = cpara (2) 
                        ELSEIF (werte (1) .eq. - 1) then 
                           k = nint (werte (1) ) 
                           CALL do_cap (cpara (2) ) 
                           DO k = 1, cr_nscat 
                           cr_scat_equ (k) = .true. 
                           cr_at_equ (k) = cpara (2) 
                           ENDDO 
                        ELSE 
                           ier_num = - 27 
                           ier_typ = ER_APPL 
                        ENDIF 
                     ENDIF 
                  ENDIF 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
!                                                                       
!     set desired mode of Fourier transform 'set'                       
!                                                                       
            ELSEIF (str_comp (befehl, 'set', 2, lbef, 3) ) then 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ier_num.eq.0) then 
                  IF (ianz.ge.1.and.ianz.le.2) then 
                     IF (str_comp (cpara (1) , 'aver', 1, lpara (1) , 4)&
                     ) then                                             
                        IF (ianz.eq.1) then 
                           fave = 0.0 
                        ELSE 
                           CALL del_params (1, ianz, cpara, lpara, maxw) 
                           CALL ber_params (ianz, cpara, lpara, werte,  &
                           maxw)                                        
                           IF (ier_num.eq.0) then 
                              IF (werte (1) .ge.0.0.and.werte (1)       &
                              .le.100.0) then                           
                                 fave = werte (1) * 0.01 
                              ELSE 
                                 ier_num = - 1 
                                 ier_typ = ER_FOUR 
                              ENDIF 
                           ENDIF 
                        ENDIF 
                     ELSEIF (str_comp (cpara (1) , 'extern', 1, lpara ( &
                     1) , 6) ) then                                     
                        four_mode = EXTERNAL 
                     ELSEIF (str_comp (cpara (1) , 'intern', 1, lpara ( &
                     1) , 6) ) then                                     
                        four_mode = INTERNAL 
                     ELSE 
                        ier_num = - 1 
                        ier_typ = ER_FOUR 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ENDIF 
!                                                                       
!     Show the current settings for the Fourier 'show'                  
!                                                                       
            ELSEIF (str_comp (befehl, 'show', 2, lbef, 4) ) then 
               CALL dlink (lxray, ano, lambda, rlambda) 
               CALL four_show 
!                                                                       
!     Switch usage of temperature coefficients on/off 'temp'            
!                                                                       
            ELSEIF (str_comp (befehl, 'temp', 1, lbef, 4) ) then 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ier_num.eq.0) then 
                  IF (ianz.eq.1) then 
                     IF (cpara (1) (1:2) .eq.'ig') then 
                        ldbw = .false. 
                     ELSEIF (cpara (1) (1:2) .eq.'us') then 
                        ldbw = .true. 
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
               IF (zeile.ne.' ') then 
                  CALL do_operating (zeile (1:lp), lp) 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
!                                                                       
!------ waiting for user input                                          
!                                                                       
            ELSEIF (str_comp (befehl, 'wait', 3, lbef, 4) ) then 
               CALL do_input (zeile, lp) 
!                                                                       
!     set the wave length to be used 'wvle'                             
!                                                                       
            ELSEIF (str_comp (befehl, 'wvle', 1, lbef, 4) ) then 
               CALL do_cap (zeile) 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ianz.eq.1) then 
                  IF (ichar ('A') .le.ichar (cpara (1) (1:1) )          &
                  .and.ichar (cpara (1) (1:1) ) .le.ichar ('Z') ) then  
                     lambda = cpara (1) 
                  ELSE 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     rlambda = werte (1) 
                     lambda = ' ' 
                  ENDIF 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
!                                                                       
!     switch to x-ray diffraction 'xray'                                
!                                                                       
            ELSEIF (str_comp (befehl, 'xray', 1, lbef, 4) ) then 
               lxray = .true. 
            ELSE 
               ier_num = - 8 
               ier_typ = ER_COMM 
            ENDIF 
         ENDIF 
      ENDIF 
!                                                                       
!------ end of command list                                             
!                                                                       
      IF (ier_num.ne.0) then 
         CALL errlist 
         IF (ier_sta.ne.ER_S_LIVE) then 
            IF (lmakro) then 
               CALL macro_close 
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
 2000 FORMAT    (a) 
!                                                                       
      END SUBROUTINE fourier                        
!*****7*****************************************************************
      SUBROUTINE four_run 
!+                                                                      
!     Main routine to calculate the Fourier transform for the           
!     given plane in reciprocal space. Based on the program             
!     'diffuse' written by B.D. Butler.                                 
!-                                                                      
      USE config_mod 
      USE crystal_mod 
      USE diffuse_mod 
      IMPLICIT none 
!                                                                       
      include'prompt.inc' 
       
      include'errlist.inc' 
!                                                                       
      REAL ss, seknds, dnorm 
      INTEGER lbeg (3), csize (3) 
      INTEGER iscat, nlot, ncell, i 
!                                                                       
      ier_num = 0 
!                                                                       
!------ preset some values                                              
!                                                                       
      ss = seknds (0.0) 
      CALL four_layer 
!                                                                       
!------ zero some arrays                                                
!                                                                       
      DO i = 1, num (1) * num (2) 
      csf (i) = cmplx (0.0d0, 0.0d0) 
      acsf (i) = cmplx (0.0d0, 0.0d0) 
      dsi (i) = 0.0d0 
      ENDDO 
!                                                                       
!------ preset some tables, calculate average structure                 
!                                                                       
      CALL four_cexpt 
      CALL four_stltab 
      IF (ier_num.ne.0) return 
      CALL four_formtab 
      CALL four_csize (cr_icc, csize, lperiod, ls_xyz) 
      CALL four_aver (ilots, fave) 
!                                                                       
!------ loop over crystal regions (lots)                                
!                                                                       
      DO nlot = 1, nlots 
      DO i = 1, num (1) * num (2) 
      csf (i) = cmplx (0.0d0, 0.0d0) 
      ENDDO 
!                                                                       
      CALL four_ranloc (csize, lbeg) 
      IF (four_log) then 
         IF (ilots.ne.LOT_OFF) then 
            WRITE (output_io, 2000) nlot, nlots, (lbeg (i), i = 1, 3) 
         ELSE 
            WRITE (output_io, 2010) nlot, nlots 
         ENDIF 
      ENDIF 
!                                                                       
!------ - loop over all different atom types                            
!                                                                       
      DO iscat = 1, cr_nscat 
      CALL four_getatm (iscat, ilots, lbeg, csize, ncell) 
      CALL four_strucf (iscat, .true.) 
!                                                                       
!------ --- Add this part of the structur factor to the total           
!                                                                       
      DO i = 1, num (1) * num (2) 
      csf (i) = csf (i) + tcsf (i) 
      ENDDO 
      IF (four_log) then 
         WRITE (output_io, 3000) cr_at_lis (iscat), nxat 
      ENDIF 
      ENDDO 
!                                                                       
!------ - subtract average structure factor, add intensity              
!                                                                       
      DO i = 1, num (1) * num (2) 
      csf (i) = csf (i) - acsf (i) 
      dsi (i) = dsi (i) + real (csf (i) * conjg (csf (i) ) ) 
      ENDDO 
      ENDDO 
!                                                                       
!------ if we had lots, normalise the intensities                       
!                                                                       
      IF (nlots.ne.1) then 
         dnorm = float (cr_icc (1) * cr_icc (2) * cr_icc (3) ) / float (&
         ncell * nlots)                                                 
         dnorm = dnorm**2 
         DO i = 1, num (1) * num (2) 
         dsi (i) = dnorm * dsi (i) 
         ENDDO 
      ENDIF 
!                                                                       
      CALL four_qinfo 
      ss = seknds (ss) 
      IF (four_log) then 
         WRITE (output_io, 4000) ss 
      ENDIF 
!                                                                       
 2000 FORMAT     (/,' Lot # ',I4,'/',I4,' : starting at unit cell (',   &
     &                       3(I3,1X),')')                              
 2010 FORMAT     (/,' Lot # ',I4,'/',I4,' : complete crystal') 
 3000 FORMAT (  '                   Atom typ = ',A4,7X,'(# ',I9,' )') 
 4000 FORMAT     (/,' Elapsed time    : ',G12.6,' sec') 
      END SUBROUTINE four_run                       
!*****7*****************************************************************
      SUBROUTINE four_strucf (iscat, lform) 
!+                                                                      
!     Here the complex structure factor of 'nxat' identical atoms       
!     from array 'xat' is computed.                                     
!-                                                                      
      USE config_mod 
      USE diffuse_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      REAL(8) xarg0, xincu, xincv 
      INTEGER iscat 
      INTEGER i, ii, j, k, iarg, iarg0, iincu, iincv, iadd 
      LOGICAL lform 
!                                                                       
      INTEGER IAND, ISHFT 
!                                                                       
!------ zero fourier array                                              
!                                                                       
      DO i = 1, num (1) * num (2) 
      tcsf (i) = cmplx (0.0d0, 0.0d0) 
      ENDDO 
!                                                                       
!------ Loop over all atoms in 'xat'                                    
!                                                                       
      DO k = 1, nxat 
      xarg0 = xm (1) * xat (k, 1) + xm (2) * xat (k, 2) + xm (3)        &
      * xat (k, 3)                                                      
      xincu = uin (1) * xat (k, 1) + uin (2) * xat (k, 2) + uin (3)     &
      * xat (k, 3)                                                      
      xincv = vin (1) * xat (k, 1) + vin (2) * xat (k, 2) + vin (3)     &
      * xat (k, 3)                                                      
!                                                                       
      iarg0 = nint (64 * I2PI * (xarg0 - int (xarg0) + 1.0d0) ) 
      iincu = nint (64 * I2PI * (xincu - int (xincu) + 1.0d0) ) 
      iincv = nint (64 * I2PI * (xincv - int (xincv) + 1.0d0) ) 
      iarg = iarg0 
!                                                                       
!------ - Loop over all points in Q. 'iadd' is the address of the       
!------ - complex exponent table. 'IADD' divides out the 64 and         
!------ - ISHFT acts as MOD so that the argument stays in the table     
!------ - boundaries.                                                   
!                                                                       
      ii = 0 
!                                                                       
      DO j = 1, num (1) 
      DO i = 1, num (2) 
      iadd = ISHFT (iarg, - 6) 
      iadd = IAND (iadd, MASK) 
      ii = ii + 1 
      tcsf (ii) = tcsf (ii) + cex (iadd) 
      iarg = iarg + iincv 
      ENDDO 
      iarg = iarg0 + iincu * j 
      ENDDO 
      ENDDO 
!                                                                       
!------ Now we multiply with formfactor                                 
!                                                                       
      IF (lform) then 
         DO i = 1, num (1) * num (2) 
         tcsf (i) = tcsf (i) * cfact (istl (i), iscat) 
         ENDDO 
      ENDIF 
!                                                                       
      END SUBROUTINE four_strucf                    
!*****7*****************************************************************
      SUBROUTINE four_aver (lots, ave) 
!+                                                                      
!     This routine calculates the average structure factor <F>.         
!-                                                                      
      USE config_mod 
      USE crystal_mod 
      USE diffuse_mod 
      USE modify_mod
      IMPLICIT none 
!                                                                       
      include'prompt.inc' 
       
      include'random.inc' 
!                                                                       
      REAL ran1, norm, ave 
      INTEGER isite, iatom, iscat, icell (3) 
      INTEGER lots, scell, ncell, j, ii, jj, kk 
      LOGICAL sel 
!                                                                       
      IF (ave.eq.0.0) then 
         RETURN 
      ELSE 
         IF (four_log) then 
            WRITE (output_io, 1000) 100.0 * ave 
         ENDIF 
         scell = cr_icc (1) * cr_icc (2) * cr_icc (3) 
         ncell = 0 
!                                                                       
!------ - Loop over all unit cells                                      
!                                                                       
         DO ii = 1, cr_icc (1) 
         DO jj = 1, cr_icc (2) 
         DO kk = 1, cr_icc (3) 
         icell (1) = ii 
         icell (2) = jj 
         icell (3) = kk 
!                                                                       
!------ --- get only 'ave'% of those unit cells and compute average     
!------ --- unit cell                                                   
!                                                                       
         sel = (ran1 (idum) .le.ave) 
         IF (sel) then 
            ncell = ncell + 1 
!                                                                       
!------ ----- Loop over all atom types                                  
!                                                                       
            DO iscat = 1, cr_nscat 
            nxat = 0 
!                                                                       
!------ ------- Loop over all sites within unit cell                    
!                                                                       
            DO isite = 1, cr_ncatoms 
            CALL celltoindex (icell, isite, iatom) 
            IF (cr_iscat (iatom) .eq.iscat) then 
               nxat = nxat + 1 
               DO j = 1, 3 
               xat (nxat, j) = cr_pos (j, iatom) - float (icell (j)     &
               - 1) - cr_dim0 (j, 1)                                    
               ENDDO 
            ENDIF 
            ENDDO 
            CALL four_strucf (iscat, .true.) 
            DO j = 1, num (1) * num (2) 
            acsf (j) = acsf (j) + tcsf (j) 
            ENDDO 
            ENDDO 
         ENDIF 
         ENDDO 
         ENDDO 
         ENDDO 
!                                                                       
!------ - now compute the interference function of the lot shape        
!                                                                       
         CALL four_getav (lots) 
         CALL four_strucf (0, .false.) 
         norm = 1.0 / ncell 
         DO j = 1, num (1) * num (2) 
         acsf (j) = acsf (j) * tcsf (j) * cmplx (dble (norm), 0.0d0) 
         ENDDO 
!                                                                       
!------ - write how much of the crystal we actually used                
!                                                                       
         IF (four_log) then 
            WRITE (output_io, 2000) (float (ncell) / float (scell) )    &
            * 100.0                                                     
         ENDIF 
!                                                                       
      ENDIF 
!                                                                       
 1000 FORMAT     (' Calculating <F> with ',F5.1,' % of the crystal ...') 
 2000 FORMAT     (' Used ',F5.1,' % of the crystal to calulate <F> ...') 
      END SUBROUTINE four_aver                      
!*****7*****************************************************************
      SUBROUTINE four_getatm (iscat, lots, lbeg, csize, ncell) 
!+                                                                      
!     This routine creates an atom list of atoms of type 'iscat'        
!     which are within the current lot.                                 
!-                                                                      
      USE config_mod 
      USE crystal_mod 
      USE diffuse_mod 
      USE modify_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      REAL offset (3) 
      REAL x0 (3), xtest (3) 
      INTEGER cr_end 
      INTEGER lbeg (3), csize (3), icell (3), jcell (3) 
      INTEGER ncell, iscat, lots, i, j, ir, ii, jj, kk, is, ia 
!                                                                       
      nxat = 0 
      ncell = 0 
      cr_end = cr_ncatoms * cr_icc (1) * cr_icc (2) * cr_icc (3)        &
      + 1                                                               
!                                                                       
!------ No lots required, return all atoms of type 'iscat'              
!                                                                       
      IF (lots.eq.LOT_OFF) then 
         ncell = cr_icc (1) * cr_icc (2) * cr_icc (3) 
         DO i = 1, cr_natoms 
         IF (cr_iscat (i) .eq.iscat) then 
            nxat = nxat + 1 
            DO j = 1, 3 
            xat (nxat, j) = cr_pos (j, i) 
            ENDDO 
         ENDIF 
         ENDDO 
!                                                                       
!------ Box shaped lot                                                  
!                                                                       
      ELSEIF (lots.eq.LOT_BOX) then 
         DO kk = 0, ls_xyz (3) - 1 
         icell (3) = kk + lbeg (3) 
         offset (3) = cr_dim0 (3, 1) + lbeg (3) - 1 
         IF (icell (3) .gt.cr_icc (3) ) then 
            icell (3) = icell (3) - cr_icc (3) 
            offset (3) = offset (3) - float (cr_icc (3) ) 
         ENDIF 
         DO jj = 0, ls_xyz (2) - 1 
         icell (2) = jj + lbeg (2) 
         offset (2) = cr_dim0 (2, 1) + lbeg (2) - 1 
         IF (icell (2) .gt.cr_icc (2) ) then 
            icell (2) = icell (2) - cr_icc (2) 
            offset (2) = offset (2) - float (cr_icc (2) ) 
         ENDIF 
         DO ii = 0, ls_xyz (1) - 1 
         icell (1) = ii + lbeg (1) 
         offset (1) = cr_dim0 (1, 1) + lbeg (1) - 1 
         IF (icell (1) .gt.cr_icc (1) ) then 
            icell (1) = icell (1) - cr_icc (1) 
            offset (1) = offset (1) - float (cr_icc (1) ) 
         ENDIF 
!                                                                       
         ncell = ncell + 1 
         DO is = 1, cr_ncatoms 
         CALL celltoindex (icell, is, ia) 
         IF (cr_iscat (ia) .eq.iscat) then 
            nxat = nxat + 1 
            DO j = 1, 3 
            xat (nxat, j) = cr_pos (j, ia) - offset (j) 
            ENDDO 
         ENDIF 
         ENDDO 
!                                                                       
         DO ir = cr_end, cr_natoms 
         DO j = 1, 3 
         jcell (j) = icell (j) + nint (cr_dim (j, 1) ) - 1 
         ENDDO 
         IF (int (cr_pos (1, ir) ) .eq.jcell (1) .and.int (cr_pos (2,   &
         ir) ) .eq.jcell (2) .and.int (cr_pos (3, ir) ) .eq.jcell (3)   &
         .and.cr_iscat (ir) .eq.iscat) then                             
            nxat = nxat + 1 
            DO j = 1, 3 
            xat (nxat, j) = cr_pos (j, ir) - offset (j) 
            ENDDO 
         ENDIF 
         ENDDO 
!                                                                       
         ENDDO 
         ENDDO 
         ENDDO 
!                                                                       
!------ Ellipsoid shaped lot                                            
!                                                                       
      ELSEIF (lots.eq.LOT_ELI) then 
         DO ii = 1, 3 
         x0 (ii) = float (ls_xyz (ii) ) / 2.0 
         ENDDO 
!                                                                       
         DO kk = 0, ls_xyz (3) - 1 
         icell (3) = kk + lbeg (3) 
         offset (3) = cr_dim0 (3, 1) + lbeg (3) - 1 
         xtest (3) = (float (kk) - x0 (3) + 0.5) **2 / x0 (3) **2 
         IF (icell (3) .gt.cr_icc (3) ) then 
            icell (3) = icell (3) - cr_icc (3) 
            offset (3) = offset (3) - float (cr_icc (3) ) 
         ENDIF 
         DO jj = 0, ls_xyz (2) - 1 
         icell (2) = jj + lbeg (2) 
         offset (2) = cr_dim0 (2, 1) + lbeg (2) - 1 
         xtest (2) = (float (jj) - x0 (2) + 0.5) **2 / x0 (2) **2 
         IF (icell (2) .gt.cr_icc (2) ) then 
            icell (2) = icell (2) - cr_icc (2) 
            offset (2) = offset (2) - float (cr_icc (2) ) 
         ENDIF 
         DO ii = 0, ls_xyz (1) - 1 
         icell (1) = ii + lbeg (1) 
         offset (1) = cr_dim0 (1, 1) + lbeg (1) - 1 
         xtest (1) = (float (ii) - x0 (1) + 0.5) **2 / x0 (1) **2 
         IF (icell (1) .gt.cr_icc (1) ) then 
            icell (1) = icell (1) - cr_icc (1) 
            offset (1) = offset (1) - float (cr_icc (1) ) 
         ENDIF 
!                                                                       
         IF ( (xtest (1) + xtest (2) + xtest (3) ) .le.1.0) then 
            ncell = ncell + 1 
            DO is = 1, cr_ncatoms 
            CALL celltoindex (icell, is, ia) 
            IF (cr_iscat (ia) .eq.iscat) then 
               nxat = nxat + 1 
               DO j = 1, 3 
               xat (nxat, j) = cr_pos (j, ia) - offset (j) 
               ENDDO 
            ENDIF 
            ENDDO 
!                                                                       
            DO ir = cr_end, cr_natoms 
            DO j = 1, 3 
            jcell (j) = icell (j) + nint (cr_dim (j, 1) ) - 1 
            ENDDO 
            IF (int (cr_pos (1, ir) ) .eq.jcell (1) .and.int (cr_pos (2,&
            ir) ) .eq.jcell (2) .and.int (cr_pos (3, ir) ) .eq.jcell (3)&
            .and.cr_iscat (ir) .eq.iscat) then                          
               nxat = nxat + 1 
               DO j = 1, 3 
               xat (nxat, j) = cr_pos (j, ir) - offset (j) 
               ENDDO 
            ENDIF 
            ENDDO 
         ENDIF 
!                                                                       
         ENDDO 
         ENDDO 
         ENDDO 
      ENDIF 
!                                                                       
      END SUBROUTINE four_getatm                    
!*****7*****************************************************************
      SUBROUTINE four_getav (lots) 
!+                                                                      
!     This routine computes the lattice inside the lot                  
!-                                                                      
      USE config_mod 
      USE crystal_mod 
      USE diffuse_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      REAL xtest (3), x0 (3) 
      INTEGER lots, ii, jj, kk 
!                                                                       
      nxat = 0 
!                                                                       
!------ No lots required, return complete lattice                       
!                                                                       
      IF (lots.eq.LOT_OFF) then 
         DO ii = 0, cr_icc (1) - 1 
         DO jj = 0, cr_icc (2) - 1 
         DO kk = 0, cr_icc (3) - 1 
         nxat = nxat + 1 
         xat (nxat, 1) = float (ii) + cr_dim0 (1, 1) 
         xat (nxat, 2) = float (jj) + cr_dim0 (2, 1) 
         xat (nxat, 3) = float (kk) + cr_dim0 (3, 1) 
         ENDDO 
         ENDDO 
         ENDDO 
!                                                                       
!------ Box shaped lot                                                  
!                                                                       
      ELSEIF (lots.eq.LOT_BOX) then 
         DO kk = 0, ls_xyz (3) - 1 
         DO jj = 0, ls_xyz (2) - 1 
         DO ii = 0, ls_xyz (1) - 1 
         nxat = nxat + 1 
         xat (nxat, 1) = float (ii) 
         xat (nxat, 2) = float (jj) 
         xat (nxat, 3) = float (kk) 
         ENDDO 
         ENDDO 
         ENDDO 
!                                                                       
!------ Ellipsoid shaped lot                                            
!                                                                       
      ELSEIF (lots.eq.LOT_ELI) then 
         DO ii = 1, 3 
         x0 (ii) = float (ls_xyz (ii) ) / 2.0 
         ENDDO 
!                                                                       
         DO kk = 0, ls_xyz (3) - 1 
         xtest (3) = (float (kk) - x0 (3) + 0.5) **2 / x0 (3) **2 
         DO jj = 0, ls_xyz (2) - 1 
         xtest (2) = (float (jj) - x0 (2) + 0.5) **2 / x0 (2) **2 
         DO ii = 0, ls_xyz (1) - 1 
         xtest (1) = (float (ii) - x0 (1) + 0.5) **2 / x0 (1) **2 
         IF ( (xtest (1) + xtest (2) + xtest (3) ) .le.1.0) then 
            nxat = nxat + 1 
            xat (nxat, 1) = float (ii) 
            xat (nxat, 2) = float (jj) 
            xat (nxat, 3) = float (kk) 
         ENDIF 
         ENDDO 
         ENDDO 
         ENDDO 
      ENDIF 
!                                                                       
      END SUBROUTINE four_getav                     
!*****7*****************************************************************
      SUBROUTINE four_ranloc (csize, lbeg) 
!+                                                                      
!     Returns a random cell from within the simulated crystal           
!     which is limited by 'csize'.                                      
!-                                                                      
      IMPLICIT none 
!                                                                       
      include'random.inc' 
!                                                                       
      INTEGER csize (3), lbeg (3) 
      REAL ran1 
!                                                                       
      lbeg (1) = int (ran1 (idum) * csize (1) ) + 1 
      lbeg (2) = int (ran1 (idum) * csize (2) ) + 1 
      lbeg (3) = int (ran1 (idum) * csize (3) ) + 1 
!                                                                       
      END SUBROUTINE four_ranloc                    
!*****7*****************************************************************
      SUBROUTINE four_csize (cr_icc, csize, lperiod, ls_xyz) 
!+                                                                      
!     Limits crystal size in case of no periodic boundary               
!     conditions.                                                       
!-                                                                      
      IMPLICIT none 
!                                                                       
      INTEGER cr_icc (3), csize (3), ls_xyz (3) 
      LOGICAL lperiod 
!                                                                       
      IF (lperiod) then 
         csize (1) = cr_icc (1) 
         csize (2) = cr_icc (2) 
         csize (3) = cr_icc (3) 
      ELSE 
         csize (1) = cr_icc (1) - ls_xyz (1) 
         csize (2) = cr_icc (2) - ls_xyz (2) 
         csize (3) = cr_icc (3) - ls_xyz (3) 
      ENDIF 
!                                                                       
      END SUBROUTINE four_csize                     
!*****7*****************************************************************
      SUBROUTINE four_layer 
!+                                                                      
!     This routine writes the corners of the plane to be calculated     
!     in the 'diffuse.inc' variables.                                   
!-                                                                      
      USE config_mod 
      USE diffuse_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER i 
!                                                                       
      DO i = 1, 3 
      xm (i) = eck (i, 1) 
      uin (i) = vi (i, 1) 
      vin (i) = vi (i, 2) 
      ENDDO 
!                                                                       
      DO i = 1, 2 
      num (i) = inc (i) 
      ENDDO 
!                                                                       
      END SUBROUTINE four_layer                     
!*****7*****************************************************************
      SUBROUTINE four_cexpt 
!+                                                                      
!     This routine initialises the complex exponent table and           
!     is called only at the first Fourier run.                          
!-                                                                      
      USE config_mod 
      USE diffuse_mod 
      IMPLICIT none 
!                                                                       
      include'prompt.inc' 
!                                                                       
      REAL(8) twopi, xmult, xarg 
      INTEGER i 
!                                                                       
      IF (.not.ffour) then 
         WRITE (output_io, 1000) 
!                                                                       
         twopi = 8.0d0 * datan (1.0d0) 
!                                                                       
         DO i = 0, MASK 
         xmult = (dble (i) * 1.0d0) / dble (I2PI) 
         xarg = twopi * xmult 
         cex (i) = cmplx (cos (xarg), sin (xarg) ) 
         ENDDO 
         ffour = .true. 
      ENDIF 
!                                                                       
 1000 FORMAT     (' Computing complex exponent table ...') 
      END SUBROUTINE four_cexpt                     
!*****7*****************************************************************
      SUBROUTINE four_stltab 
!+                                                                      
!     Sets up an integer array containing the corresponding integers    
!     to the formfactor table for each sin(theta)/lambda.               
!-                                                                      
      USE config_mod 
      USE crystal_mod 
      USE diffuse_mod 
      IMPLICIT none 
!                                                                       
      include'prompt.inc' 
       
      include'errlist.inc' 
!                                                                       
      REAL quad, q2, h (3) 
      INTEGER i, j, k 
!                                                                       
      IF (four_log) then 
         WRITE (output_io, 1000) 
      ENDIF 
!                                                                       
      DO j = 1, num (2) 
      DO i = 1, num (1) 
      DO k = 1, 3 
      h (k) = xm (k) + uin (k) * float (i - 1) + vin (k) * float (j - 1) 
      ENDDO 
      q2 = quad (h, h, cr_rten) / 4.0 
      k = (i - 1) * num (2) + j 
      istl (k) = nint (sqrt (q2) * (1.0 / CFINC) ) 
!                                                                       
      IF (istl (k) .gt.CFPKT) then 
         ier_num = - 3 
         ier_typ = ER_FOUR 
         RETURN 
      ENDIF 
!                                                                       
      ENDDO 
      ENDDO 
!                                                                       
 1000 FORMAT     (' Computing sin(theta)/lambda table ...') 
      END SUBROUTINE four_stltab                    
!*****7*****************************************************************
      SUBROUTINE four_formtab 
!+                                                                      
!     This routine sets up the complex formfactor lookup table          
!     for all atom types. The range in sin(theta)/lambda is             
!     0 -> 2 in steps of 0.001. These values can be changed             
!     in the 'diffuse.inc' file.                                        
!-                                                                      
      USE config_mod 
      USE crystal_mod 
      USE diffuse_mod 
      IMPLICIT none 
!                                                                       
      include'prompt.inc' 
       
!                                                                       
      REAL q2, sb, sf, sfp, sfpp 
      INTEGER iq, iscat 
!                                                                       
      REAL form 
!                                                                       
      IF (four_log) then 
         WRITE (output_io, 1000) 
      ENDIF 
!                                                                       
      DO iscat = 1, cr_nscat 
      DO iq = 0, CFPKT 
      q2 = (float (iq) * CFINC) **2 
      sf = form (iscat, cr_scat, lxray, q2) 
!                                                                       
      IF (ano) then 
         sfp = cr_delfr ( (iscat) ) 
         sfpp = cr_delfi ( (iscat) ) 
      ELSE 
         sfp = 0.0 
         sfpp = 0.0 
      ENDIF 
!                                                                       
      IF (ldbw) then 
         sb = exp ( - cr_dw ( (iscat) ) * q2) 
      ELSE 
         sb = 1.0 
      ENDIF 
!                                                                       
      cfact (iq, iscat) = cmplx (sb * (sf + sfp), sb * sfpp) 
      ENDDO 
      ENDDO 
!                                                                       
 1000 FORMAT     (' Computing formfactor lookup table ...') 
      END SUBROUTINE four_formtab                   
!*****7*****************************************************************
      SUBROUTINE four_qinfo 
!+                                                                      
!     Gives information about max/min values for diffuse and            
!     Bragg scattering.                                                 
!-                                                                      
      USE config_mod 
      USE diffuse_mod 
      IMPLICIT none 
!                                                                       
      include'prompt.inc' 
!                                                                       
      INTEGER i, j, k, nd 
      REAL h (3), dsum, dsum2 
      LOGICAL lbragg 
!                                                                       
      braggmax = - 9999.0 
      braggmin = 1E19 
      diffumax = - 9999.0 
      diffumin = 1E19 
!                                                                       
      nd = 0 
      dsum = 0.0 
      dsum2 = 0.0 
!                                                                       
      DO j = 1, num (2) 
      DO i = 1, num (1) 
      lbragg = .true. 
      DO k = 1, 3 
      h (k) = eck (k, 1) + vi (k, 1) * float (i - 1) + vi (k, 2)        &
      * float (j - 1)                                                   
      lbragg = lbragg.and.amod (h (k), 1.0) .eq.0.0 
      ENDDO 
      k = (i - 1) * num (2) + j 
      IF (lbragg) then 
         braggmax = max (braggmax, dsi (k) ) 
         braggmin = min (braggmin, dsi (k) ) 
      ELSE 
         diffumax = max (diffumax, dsi (k) ) 
         diffumin = min (diffumin, dsi (k) ) 
!                                                                       
         dsum = dsum + dsi (k) 
         dsum2 = dsum2 + dsi (k) **2 
         nd = nd+1 
      ENDIF 
      ENDDO 
      ENDDO 
!                                                                       
      IF (nd.ne.0) then 
         diffuave = dsum / float (nd) 
         diffusig = sqrt (dsum2 / float (nd) - diffuave**2) 
      ELSE 
         diffuave = 0.0 
         diffusig = 0.0 
      ENDIF 
!                                                                       
      IF (four_log) then 
         WRITE (output_io, 1000) 
         IF (braggmax.ge.0.0) then 
            WRITE (output_io, 1010) braggmin, braggmax 
         ENDIF 
         WRITE (output_io, 1020) diffumin, diffumax 
         WRITE (output_io, 1030) diffuave, diffusig 
      ENDIF 
!                                                                       
 1000 FORMAT     (/,' ') 
 1010 FORMAT     (  ' Bragg scat.     : ',G12.6,'  -> ',G12.6) 
 1020 FORMAT     (  ' Diffuse scat.   : ',G12.6,'  -> ',G12.6) 
 1030 FORMAT     (  '      Average    : ',G12.6,'  +- ',G12.6) 
      END SUBROUTINE four_qinfo                     
!*****7*****************************************************************
      SUBROUTINE four_show 
!+                                                                      
!     prints summary of current fourier settings                        
!-                                                                      
      USE config_mod 
      USE diffuse_mod 
      USE output_mod 
      IMPLICIT none 
!                                                                       
      include'prompt.inc' 
      include'errlist.inc' 
!                                                                       
      CHARACTER(7) radiation 
      CHARACTER(1) extr_achs (0:3) 
      REAL u (3), v (3), w (3) 
      REAL dvi1, dvi2, dvi3, dvi4, dvi5 
      INTEGER i, j 
      LOGICAL lspace 
!                                                                       
      REAL do_blen, do_bang 
!                                                                       
      DATA extr_achs / ' ', 'h', 'k', 'l' / 
!                                                                       
      IF (fave.eq.0.0) then 
         WRITE (output_io, 1010) 
      ELSE 
         WRITE (output_io, 1020) fave * 100.0 
      ENDIF 
      IF (four_mode.eq.INTERNAL) then 
         WRITE (output_io, 1030) 'atom form factors' 
      ELSEIF (four_mode.eq.EXTERNAL) then 
         WRITE (output_io, 1030) 'object form factors' 
      ENDIF 
!                                                                       
      IF (ilots.eq.LOT_OFF) then 
         WRITE (output_io, 1100) 
      ELSEIF (ilots.eq.LOT_BOX) then 
         WRITE (output_io, 1110) nlots 
         WRITE (output_io, 1130) (ls_xyz (i), i = 1, 3), lperiod 
      ELSEIF (ilots.eq.LOT_ELI) then 
         WRITE (output_io, 1120) nlots 
         WRITE (output_io, 1130) (ls_xyz (i), i = 1, 3), lperiod 
      ENDIF 
!                                                                       
      radiation = 'neutron' 
      IF (lxray) radiation = 'x-ray' 
      IF (lambda.eq.' ') then 
         WRITE (output_io, 1200) radiation, rlambda 
      ELSE 
         WRITE (output_io, 1210) radiation, lambda, rlambda 
      ENDIF 
!                                                                       
      IF (ldbw) then 
         WRITE (output_io, 1300) 'used' 
      ELSE 
         WRITE (output_io, 1300) 'ignored' 
      ENDIF 
!                                                                       
      IF (ano) then 
         WRITE (output_io, 1310) 'used' 
      ELSE 
         WRITE (output_io, 1310) 'ignored' 
      ENDIF 
!                                                                       
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
         IF (abs (u (extr_abs) ) .gt.0.0.and.abs (w (extr_ord) )        &
         .gt.0.0) then                                                  
            dvi5 = (dvi2 / w (extr_ord) ) / (dvi1 / u (extr_abs) ) 
         ELSE 
            ier_num = - 4 
            ier_typ = ER_FOUR 
         ENDIF 
      ELSE 
         dvi4 = 0.0 
         dvi5 = 0.0 
      ENDIF 
!                                                                       
      WRITE (output_io, 1400) ( (eck (i, j), i = 1, 3), j = 1, 3) 
      WRITE (output_io, 1410) (vi (i, 1), i = 1, 3), dvi1, (vi (i, 2),  &
      i = 1, 3), dvi2                                                   
      WRITE (output_io, 1420) (inc (i), i = 1, 2), extr_achs (extr_abs),&
      extr_achs (extr_ord)                                              
      WRITE (output_io, 1430) dvi3, dvi4, dvi5 
!                                                                       
 1010 FORMAT (  ' Fourier technique    : turbo Fourier') 
 1020 FORMAT (  ' Fourier technique    : turbo Fourier, minus <F>',     &
     &          ' (based on ',F5.1,'% of cryst.)')                      
 1030 FORMAT (  ' Fourier calculated by: ',a) 
 1100 FORMAT (  '   Fourier volume     : complete crystal') 
 1110 FORMAT (  '   Fourier volume     : ',I4,' box shaped lots') 
 1120 FORMAT (  '   Fourier volume     : ',I4,' ellipsoid shaped lots') 
 1130 FORMAT (  '   Lot size           : ',I3,' x ',I3,' x ',I3,        &
     &          ' unit cells (periodic boundaries = ',L1,')')           
 1200 FORMAT (  '   Radiation          : ',A,', wavelength = ',         &
     &          F7.4,' A')                                              
 1210 FORMAT (  '   Radiation          : ',A,', wavelength = ',A4,      &
     &          ' = ',F7.4,' A')                                        
 1300 FORMAT (  '   Temp. factors      : ',A) 
 1310 FORMAT (  '   Anomalous scat.    : ',A) 
 1400 FORMAT (/,' Reciprocal layer     : ',/                            &
     &          '   lower left  corner : ',3(2x,f9.4),/                 &
     &          '   lower right corner : ',3(2x,f9.4),/                 &
     &          '   upper left  corner : ',3(2x,f9.4))                  
 1410 FORMAT (  '   hor. increment     : ',3(2x,f9.4),2x,               &
     &          ' -> ',f9.4,' A**-1',/                                  &
     &          '   vert. increment    : ',3(2x,f9.4),2x,               &
     &          ' -> ',f9.4,' A**-1')                                   
 1420 FORMAT (  '   # of points        :  ',I4,' x',I4,'  (',           &
     &          A1,',',A1,')')                                          
 1430 FORMAT (  '   Angle              : ',2x,f9.4,' degrees',/         &
     &          '   Ratio/Aver  v/h    : ',2(2x,f9.4))                  
      END SUBROUTINE four_show                      
!*****7*****************************************************************
      REAL function form (ll, scat, lxray, h2) 
!+                                                                      
!       calculates the form factor                                      
!-                                                                      
      IMPLICIT none 
!                                                                       
      INTEGER idim, maxscat 
      PARAMETER (idim = 3, maxscat = 150) 
!                                                                       
      LOGICAL lxray 
      INTEGER ll, i 
      REAL h2, scat (9, 0:maxscat) 
!                                                                       
      form = scat (1, ll) 
      IF (lxray) then 
         DO i = 1, 4 
         form = form + scat (2 * i, ll) * exp ( - scat (2 * i + 1, ll)  &
         * h2)                                                          
         ENDDO 
      ENDIF 
      END FUNCTION form                             
!*****7*****************************************************************
      REAL function quad (h, k, rten) 
!+                                                                      
!           Calculates the scalar product of h and k.                   
!           1/d**2 = h(i)*k(j)*rten(i,j)                                
!-                                                                      
      IMPLICIT none 
!                                                                       
      INTEGER idim 
      PARAMETER (idim = 3) 
!                                                                       
      INTEGER i, j 
      REAL h (idim), k (idim), rten (idim, idim) 
!                                                                       
      quad = 0.0 
      DO i = 1, idim 
      DO j = 1, idim 
      quad = quad+h (i) * k (j) * rten (i, j) 
      ENDDO 
      ENDDO 
      END FUNCTION quad                             
!*****7*****************************************************************
      SUBROUTINE dlink (lxray, ano, lambda, rlambda) 
!-                                                                      
!     This routine reads wavelength symbols, wavelength values 
!     and atomic form factors from module "element_data_mod"
!     It replaces the old dlink and LAZY_PULVERIX
!+                                                                      
      USE config_mod 
      USE crystal_mod 
      USE element_data_mod
      IMPLICIT none 
!                                                                       
      include'charact.inc' 
      include'errlist.inc' 
!                                                                       
      LOGICAL             , INTENT(IN)   :: ano
      LOGICAL             , INTENT(IN)   :: lxray 
      CHARACTER (LEN = * ), INTENT(IN)   :: lambda 
      REAL                , INTENT(OUT)  :: rlambda
!
      CHARACTER (LEN = 4 ) :: element 
      INTEGER    :: i
      INTEGER    :: j
      REAL   , DIMENSION(1:9)   :: temp_scat  ! a1,b1,---a4,b4,c
      REAL   , DIMENSION(1:2)   :: temp_delf  ! delfr, delfi
      REAL                      :: temp_bcoh  ! b_choherent
!
!     CHARACTER(4) line 
!     INTEGER element, symwl, kodlp, i, j 
!     LOGICAL ano, lxray 
!     REAL fa (4, 8), fb (4, 8), fc (8), rlambda 
!     REAL delfr, delfi, fneu 
!     REAL wave 
!                                                                       
      ier_num = -77 
      ier_typ = ER_APPL 
!
      CALL get_wave ( lambda, rlambda, ier_num, ier_typ )
!                                                                       
      IF (ier_num.ne.0) RETURN 
!                                                                       
      any_element: IF (cr_nscat.gt.0) THEN 
         DO i = 1, cr_nscat 
         IF (cr_at_lis (i) /= 'XAXI' .AND. cr_at_lis (i) /= 'YAXI' .and. & 
             cr_at_lis (i) /= 'ZAXI') THEN                  
            IF (cr_scat_int (i) ) then 
               IF (.not.lxray) then     !  neutron scattering
                  IF (cr_scat_equ (i) ) then 
                     element =         cr_at_equ (i) 
                  ELSE 
                     element =         cr_at_lis (i) 
                  ENDIF 
!                                                                       
                  CALL symbf ( element, j)
                  IF ( j /= 0 ) THEN 
                     CALL get_scat ( j, temp_scat, temp_delf, temp_bcoh )
                     cr_scat(:,i) = 0.0
                     cr_scat(1,i) = temp_bcoh
                     cr_delfr (i) = 0.0 
                     cr_delfi (i) = 0.0 
                     ier_num = 0
                     ier_typ = ER_NONE 
                  ELSE
                     ier_num = -20 
                     ier_typ = ER_APPL 
                  ENDIF
               ELSE                    ! Xray diffraction
                  IF (cr_scat_equ (i) ) then 
                     element =         cr_at_equ (i) 
                  ELSE 
                     element =         cr_at_lis (i) 
                  ENDIF 
!
                  ier_num = -20 
                  ier_typ = ER_APPL 
!                                                                       
                  CALL symbf ( element, j)
                  IF ( j /= 0 ) THEN 
                     CALL get_scat ( j, temp_scat, temp_delf, temp_bcoh )
                     cr_scat(:,i) = 0.0
                     cr_scat(:,i) = temp_scat(:)   ! copy temp array into 1st column
                     cr_delfr (i) = 0.0   ! DEVELOPMENT !
                     cr_delfi (i) = 0.0   ! DEVELOPMENT !
!                       IF (ano.and.cr_delf_int (i) ) then 
!                          CALL anoma (symwl, element, kodlp) 
!                          cr_delfr (i) = delfr 
!                          cr_delfi (i) = delfi 
!                       ENDIF 
                     ier_num = 0
                  ELSE
                     ier_typ = ER_NONE 
                     ier_typ = ER_APPL 
                  ENDIF
               ENDIF 
!                                                                       
            ENDIF 
         ENDIF 
         ENDDO 
      ELSE any_element
         ier_num = - 21 
         ier_typ = ER_APPL 
      ENDIF any_element
!                                                                       
      END SUBROUTINE dlink                          
!*****7*****************************************************************
!     SUBROUTINE dlink_old (lxray, ano, lambda, rlambda) 
!-                                                                      
!           This routine dlinks between the DISCUS package and the      
!       routines                                                        
!           that generate the atomic form factors. These routines       
!       have been                                                       
!           adapted from the program LAZY.                              
!+                                                                      
!     USE config_mod 
!     USE crystal_mod 
!     IMPLICIT none 
!                                                                       
!     include'charact.inc' 
!      
!     include'errlist.inc' 
!                                                                       
!     COMMON / anomal / delfr, delfi 
!     COMMON / neu / fneu 
!     COMMON / scat / fa, fb, fc 
!                                                                       
!     CHARACTER ( * ) lambda 
!     CHARACTER(4) line 
!     INTEGER element, symwl, kodlp, i, j 
!     LOGICAL ano, lxray 
!     REAL fa (4, 8), fb (4, 8), fc (8), rlambda 
!     REAL delfr, delfi, fneu 
!     REAL wave 
!                                                                       
!     ier_num = 0 
!     ier_typ = ER_NONE 
!                                                                       
!     WRITE (line, 100) lambda 
!     READ (line, 100) symwl 
!     IF (.not.lxray) then 
!        kodlp = 2 
!     ELSE 
!        kodlp = 1 
!     IF (lambda.ne.'    ') rlambda = wave (symwl) 
!     ENDIF 
!                                                                       
!     IF (ier_num.ne.0) return 
!                                                                       
!     IF (cr_nscat.gt.0) then 
!        DO i = 1, cr_nscat 
!        IF (cr_at_lis (i) .ne.'XAXI'.and.cr_at_lis (i)                 &
!        .ne.'YAXI'.and.cr_at_lis (i) .ne.'ZAXI') then                  
!           IF (cr_scat_int (i) ) then 
!              IF (.not.lxray) then 
!                 IF (cr_scat_equ (i) ) then 
!                    WRITE (line, 100) cr_at_equ (i) 
!                 ELSE 
!                    WRITE (line, 100) cr_at_lis (i) 
!                 ENDIF 
!                                                                       
!     --------For neutron scattering remove the charge from an ion-name 
!                                                                       
!     line (3:4)  = '  ' 
!                 j = ichar (line (2:2) ) 
!                 IF (zero.le.j.and.j.le.nine) then 
!                    line (2:2) = ' ' 
!                 ENDIF 
!                 READ (line, 100) element 
!                 CALL anoma (symwl, element, kodlp) 
!                 cr_scat (1, i) = fneu 
!                 cr_delfr (i) = 0.0 
!                 cr_delfi (i) = 0.0 
!                 DO j = 2, 9 
!                 cr_scat (j, i) = 0.0 
!                 ENDDO 
!              ELSE 
!                 IF (cr_scat_equ (i) ) then 
!                    WRITE (line, 100) cr_at_equ (i) 
!                 ELSE 
!                    WRITE (line, 100) cr_at_lis (i) 
!                 ENDIF 
!                 READ (line, 100) element 
!             if(ano .and. .not. cr_delf_int(i) ) then                  
!                 IF (ano.and.cr_delf_int (i) ) then 
!                    CALL anoma (symwl, element, kodlp) 
!                    cr_delfr (i) = delfr 
!                    cr_delfi (i) = delfi 
!                 ENDIF 
!                 CALL fsc (element, 1) 
!                 cr_scat (1, i) = fc (1) 
!                 DO j = 1, 4 
!                 cr_scat (2 * j, i) = fa (j, 1) 
!                 cr_scat (2 * j + 1, i) = fb (j, 1) 
!                 ENDDO 
!              ENDIF 
!                                                                       
!              IF (element.eq.0) then 
!                 ier_num = - 20 
!                 ier_typ = ER_APPL 
!                 GOTO 999 
!              ENDIF 
!           ENDIF 
!        ENDIF 
!        ENDDO 
!     ELSE 
!        ier_num = - 21 
!        ier_typ = ER_APPL 
!     ENDIF 
!                                                                       
! 999 CONTINUE 
!                                                                       
! 100 FORMAT      (a4) 
! 101 FORMAT      (a2) 
!     END SUBROUTINE dlink_old                          
!*****7*****************************************************************
      SUBROUTINE calc_000 (rhkl) 
!-                                                                      
!     Calculates the value of F(rh,rk,rl)                               
!+                                                                      
      USE config_mod 
      USE crystal_mod 
      USE diffuse_mod 
      IMPLICIT none 
!                                                                       
       
      include'errlist.inc' 
!     include       'output.inc'                                        
      include'param.inc' 
!                                                                       
      INTEGER i, j 
      INTEGER shel_inc (2) 
      REAL shel_eck (3, 3) 
      REAL shel_vi (3, 2) 
      COMPLEX shel_acsf 
      REAL shel_dsi 
      COMPLEX shel_tcsf 
      REAL rhkl (3) 
!                                                                       
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
      eck (i, j) = rhkl (i) 
      ENDDO 
      DO j = 1, 2 
      vi (i, j) = 0.0 
      ENDDO 
      ENDDO 
      four_log = .false. 
      CALL four_run 
      res_para (1) = real (csf (1) ) 
      res_para (2) = aimag (csf (1) ) 
      res_para (3) = res_para (1) / cr_icc (1) / cr_icc (2) / cr_icc (3) 
      res_para (4) = res_para (2) / cr_icc (1) / cr_icc (2) / cr_icc (3) 
      res_para (0) = 4 
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
!                                                                       
      END SUBROUTINE calc_000                       
!*****7*****************************************************************
      SUBROUTINE errlist_four 
!-                                                                      
!     Displays error Messages for the error type RMC                    
!+                                                                      
      IMPLICIT none 
!                                                                       
      include'errlist.inc' 
!                                                                       
      INTEGER iu, io 
      PARAMETER (IU = - 15, IO = 0) 
!                                                                       
      CHARACTER(41) ERROR (IU:IO) 
!                                                                       
      DATA ERROR / ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '&
     &, 'Component of increment vector is zero    ', 'SIN(THETA)/LAMBDA &
     &> lookup table limits  ', 'Invalid lot shape selected             &
     &  ', 'Invalid Fourier mode selected            ', ' ' /           
                  !-15                                                  
                  !-14                                                  
                  !-13                                                  
                  !-12                                                  
                  !-11                                                  
                  !-10                                                  
                  !-9                                                   
                  !-8                                                   
                  !-7                                                   
                  !-6                                                   
                  !-5                                                   
                                                          !-4           
                                                          !-3           
                                                          !-2           
                                                          !-1           
                  !  0                                                  
!                                                                       
      CALL disp_error ('FOUR', error, iu, io) 
      END SUBROUTINE errlist_four                   
