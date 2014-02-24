MODULE transform_menu
!
CONTAINS
!
!+                                                                      
!     Generalized unit cell transformation operations:                  
!     The user provides the transformation matrix between new and       
!     old coordinate system and the program transforms the old          
!     atom coordinates into the new ones.                               
!     Individual direct and reciprocal space transformation can be      
!     calculated as well.                                               
!                                                                       
!*****7*****************************************************************
      SUBROUTINE transform 
!-                                                                      
!     Main menu for generalized transformation operations               
!+                                                                      
      USE config_mod 
      USE allocate_appl_mod
      USE crystal_mod 
      USE modify_mod
      USE show_menu
      USE transfrm_mod 
      USE trans_sup_mod
!
      USE doact_mod 
      USE errlist_mod 
      USE learn_mod 
      USE macro_mod 
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER, PARAMETER :: MIN_PARA = 20  ! A command requires at leaset these no of parameters
      INTEGER maxw 
      LOGICAL lold 
      PARAMETER (lold = .false.) 
!                                                                       
      CHARACTER(LEN=1024), DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: cpara ! (MAXSCAT) 
      INTEGER            , DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: lpara ! (MAXSCAT)
      REAL               , DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: werte ! (MAXSCAT) 
!
      CHARACTER(5) befehl 
      CHARACTER(50) prom 
      CHARACTER(1024) line, zeile
      CHARACTER(1024) tran_hfile 
      INTEGER lp, length, lbef 
      INTEGER tran_hlen 
      INTEGER indxg, ianz, i, j 
      LOGICAL lend, lchange, lscreen 
      REAL hkl (4) 
!                                                                       
      INTEGER len_str 
      LOGICAL str_comp 
!                                                                       
!                                                                       
      DATA lchange / .true. / 
      DATA lscreen / .true. / 
!                                                                       
      maxw = MAX(MIN_PARA,MAXSCAT+1)
      lend = .false. 
      CALL no_error 
!
      IF( cr_nscat > TRAN_MAXSCAT) THEN
         CALL alloc_transfrm ( cr_nscat )
         IF ( ier_num < 0 ) THEN
            RETURN
         ENDIF
      ENDIF
!                                                                       
      DO while (.not.lend) 
      prom = prompt (1:len_str (prompt) ) //'/tran' 
      CALL get_cmd (line, length, befehl, lbef, zeile, lp, prom) 
      IF (ier_num.eq.0) then 
         IF (line.ne.' '.and.line (1:1) .ne.'#') then 
!                                                                       
!     ----search for "="                                                
!                                                                       
            indxg = index (line, '=') 
      IF (indxg.ne.0.and..not. (str_comp (befehl, 'echo', 2, lbef, 4) ) &
     &.and..not. (str_comp (befehl, 'syst', 2, lbef, 4) ) .and..not. (st&
     &r_comp (befehl, 'help', 2, lbef, 4) .or.str_comp (befehl, '?   ', &
     &2, lbef, 4) ) ) then                                              
!                                                                       
!     ------evaluatean expression and assign the value to a variabble   
!                                                                       
               CALL do_math (line, indxg, length) 
            ELSE 
!                                                                       
!------ ----execute a macro file                                        
!                                                                       
               IF (befehl (1:1) .eq.'@') then 
                  IF (length.ge.2) then 
                     CALL file_kdo (line (2:length), length - 1) 
                  ELSE 
                     ier_num = - 13 
                     ier_typ = ER_MAC 
                  ENDIF 
!                                                                       
!     ----a(new) in terms of old axes                                   
!                                                                       
               ELSEIF (str_comp (befehl, 'anew', 2, lbef, 4) ) then 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0.and.ianz.eq.3) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        DO j = 1, 3 
                        tran_g (1, j) = werte (j) 
                        ENDDO 
                        tran_inp = TRAN_INP_G 
                        lchange = .true. 
                     ENDIF 
                  ENDIF 
!                                                                       
!     ----a(old) in terms of new axes                                   
!                                                                       
               ELSEIF (str_comp (befehl, 'aold', 2, lbef, 4) ) then 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0.and.ianz.eq.3) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        DO j = 1, 3 
                        tran_gi (1, j) = werte (j) 
                        ENDDO 
                        tran_inp = TRAN_INP_GI 
                        lchange = .true. 
                     ENDIF 
                  ENDIF 
!                                                                       
!     ----list asymmetric unit 'asym'                                   
!                                                                       
               ELSEIF (str_comp (befehl, 'asym', 2, lbef, 4) ) then 
                  CALL show_asym 
!                                                                       
!     ----b(new) in terms of old axes                                   
!                                                                       
               ELSEIF (str_comp (befehl, 'bnew', 2, lbef, 4) ) then 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0.and.ianz.eq.3) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        DO j = 1, 3 
                        tran_g (2, j) = werte (j) 
                        ENDDO 
                        tran_inp = TRAN_INP_G 
                        lchange = .true. 
                     ENDIF 
                  ENDIF 
!                                                                       
!     ----b(old) in terms of new axes                                   
!                                                                       
               ELSEIF (str_comp (befehl, 'bold', 2, lbef, 4) ) then 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0.and.ianz.eq.3) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        DO j = 1, 3 
                        tran_gi (2, j) = werte (j) 
                        ENDDO 
                        tran_inp = TRAN_INP_GI 
                        lchange = .true. 
                     ENDIF 
                  ENDIF 
!                                                                       
!     ----calculate a single transformation operation 'c2new'           
!                                                                       
               ELSEIF (str_comp (befehl, 'c2new', 2, lbef, 5) ) then 
                  IF (lchange) then 
                     CALL tran_setup 
                  ENDIF 
                  lchange = .false. 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (ianz.eq.3) then 
                        cpara (4) = 'd' 
                        lpara (4) = 1 
                     ENDIF 
                     ianz = 3 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        DO i = 1, 3 
                        hkl (i) = werte (i) 
                        ENDDO 
                        IF (str_comp (cpara (4) , 'd', 1, lpara (4) , 1)&
                        ) then                                          
                           hkl (4) = 1.0 
                           CALL tran_ca (hkl, tran_f, lscreen) 
                        ELSEIF (str_comp (cpara (4) , 'r', 1, lpara (4) &
                        , 1) ) then                                     
                           hkl (4) = 0.0 
                           CALL tran_ca (hkl, tran_g, lscreen) 
                        ELSE 
                           ier_num = - 6 
                           ier_typ = ER_COMM 
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
!     ----calculate a single transformation operation 'c2old'           
!                                                                       
               ELSEIF (str_comp (befehl, 'c2old', 2, lbef, 5) ) then 
                  IF (lchange) then 
                     CALL tran_setup 
                  ENDIF 
                  lchange = .false. 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (ianz.eq.3) then 
                        cpara (4) = 'd' 
                        lpara (4) = 1 
                     ENDIF 
                     ianz = 3 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        DO i = 1, 3 
                        hkl (i) = werte (i) 
                        ENDDO 
                        IF (str_comp (cpara (4) , 'd', 1, lpara (4) , 1)&
                        ) then                                          
                           hkl (4) = 1.0 
                           CALL tran_ca (hkl, tran_fi, lscreen) 
                        ELSEIF (str_comp (cpara (4) , 'r', 1, lpara (4) &
                        , 1) ) then                                     
                           hkl (4) = 0.0 
                           CALL tran_ca (hkl, tran_gi, lscreen) 
                        ELSE 
                           ier_num = - 6 
                           ier_typ = ER_COMM 
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
!     ----c(new) in terms of old axes                                   
!                                                                       
               ELSEIF (str_comp (befehl, 'cnew', 2, lbef, 4) ) then 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0.and.ianz.eq.3) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        DO j = 1, 3 
                        tran_g (3, j) = werte (j) 
                        ENDDO 
                        tran_inp = TRAN_INP_G 
                        lchange = .true. 
                     ENDIF 
                  ENDIF 
!                                                                       
!     ----c(old) in terms of new axes                                   
!                                                                       
               ELSEIF (str_comp (befehl, 'cold', 2, lbef, 4) ) then 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0.and.ianz.eq.3) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        DO j = 1, 3 
                        tran_gi (3, j) = werte (j) 
                        ENDDO 
                        tran_inp = TRAN_INP_GI 
                        lchange = .true. 
                     ENDIF 
                  ENDIF 
!                                                                       
!     ----continues a macro 'continue'                                  
!                                                                       
               ELSEIF (str_comp (befehl, 'continue', 2, lbef, 8) ) then 
                  CALL macro_continue (zeile, lp) 
!                                                                       
!     ----list atoms present in the crystal 'chem'                      
!                                                                       
               ELSEIF (str_comp (befehl, 'chem', 2, lbef, 4) ) then 
                  CALL show_chem 
!                                                                       
!     ----Deselect which atoms are included in the wave 'dese'          
!                                                                       
               ELSEIF (str_comp (befehl, 'dese', 1, lbef, 4) ) then 
                  CALL atom_select (zeile, lp, 0, TRAN_MAXSCAT, tran_latom, &
                  tran_sel_atom, lold,.false.)
!
!                 ier_num = - 6 
!                 ier_typ = ER_COMM 
!                 CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
!                 IF (ier_num.eq.0) then 
!                    CALL get_iscat (ianz, cpara, lpara, werte, maxw,   &
!                    lold)                                              
!                    IF (ier_num.eq.0) then 
!                       IF (werte (1) .eq. - 1) then 
!                          DO i = 0, cr_nscat 
!                          tran_latom (i) = .false. 
!                          ENDDO 
!                       ELSE 
!                          DO i = 1, ianz 
!                          tran_latom (nint (werte (i) ) ) = .false. 
!                          ENDDO 
!                       ENDIF 
!                    ENDIF 
!                 ENDIF 
!                                                                       
!------ ----Echo a string, just for interactive check in a macro 'echo' 
!                                                                       
               ELSEIF (str_comp (befehl, 'echo', 2, lbef, 4) ) then 
                  CALL echo (zeile, lp) 
!                                                                       
!      ---Evaluate an expression, just for interactive check 'eval'     
!                                                                       
               ELSEIF (str_comp (befehl, 'eval', 2, lbef, 4) ) then 
                  CALL do_eval (zeile, lp) 
!                                                                       
!     ----exit 'exit'                                                   
!                                                                       
               ELSEIF (str_comp (befehl, 'exit', 2, lbef, 4) ) then 
                  lend = .true. 
!                                                                       
!     ----help 'help','?'                                               
!                                                                       
      ELSEIF (str_comp (befehl, 'help', 2, lbef, 4) .or.str_comp (befehl&
     &, '?   ', 1, lbef, 4) ) then                                      
                  IF (str_comp (zeile, 'errors', 2, lp, 6) ) then 
                     lp = lp + 7 
                     CALL do_hel ('discus '//zeile, lp) 
                  ELSE 
                     lp = lp + 12 
                     CALL do_hel ('discus tran '//zeile, lp) 
                  ENDIF 
!                                                                       
!     ----transforme a list of reflections to new base 'h2new'          
!                                                                       
               ELSEIF (str_comp (befehl, 'h2new', 2, lbef, 5) ) then 
                  IF (lchange) then 
                     CALL tran_setup 
                  ENDIF 
                  lchange = .false. 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (ianz.ge.1) then 
                        CALL do_build_name (ianz, cpara, lpara, werte,  &
                        maxw, 1)                                        
                        IF (ier_num.eq.0) then 
                           tran_hfile = cpara (1) 
                           tran_hlen = lpara (1) 
                           CALL tran_hkl (tran_hfile, tran_hlen, tran_g) 
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
!     ----transforme a list of reflections to old base 'h2old'          
!                                                                       
               ELSEIF (str_comp (befehl, 'h2old', 2, lbef, 5) ) then 
                  IF (lchange) then 
                     CALL tran_setup 
                  ENDIF 
                  lchange = .false. 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (ianz.eq.1) then 
                        CALL do_build_name (ianz, cpara, lpara, werte,  &
                        maxw, 1)                                        
                        IF (ier_num.eq.0) then 
                           tran_hfile = cpara (1) 
                           tran_hlen = lpara (1) 
                           CALL tran_hkl (tran_hfile, tran_hlen,        &
                           tran_gi)                                     
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
!     ----Select range of atoms within crystal to be included 'incl'    
!                                                                       
               ELSEIF (str_comp (befehl, 'incl', 1, lbef, 4) ) then 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (ianz.eq.2) then 
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxw)                                           
                        IF (ier_num.eq.0) then 
                           tran_start = nint (werte (1) ) 
                           tran_end = nint (werte (2) ) 
                        ENDIF 
                     ELSEIF (ianz.eq.1) then 
                        IF (str_comp (cpara (1) , 'all', 1, lpara (1) , &
                        3) ) then                                       
                           tran_start = 1 
                           tran_end = - 1 
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
!------ ----Select the location of the new origin in old                
!                  coordinates 'onew'                                   
!                                                                       
               ELSEIF (str_comp (befehl, 'onew', 2, lbef, 4) ) then 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (ianz.eq.3) then 
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxw)                                           
                        IF (ier_num.eq.0) then 
                           tran_oold = .true. 
                           DO i = 1, 3 
                           tran_orig (i) = werte (i) 
                           ENDDO 
                           lchange = .true. 
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
!     ----Select the location of the old origin in new                  
!                  coordinates 'oold'                                   
!                                                                       
               ELSEIF (str_comp (befehl, 'oold', 2, lbef, 4) ) then 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (ianz.eq.3) then 
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxw)                                           
                        IF (ier_num.eq.0) then 
                           tran_oold = .false. 
                           DO i = 1, 3 
                           tran_orig (i) = werte (i) 
                           ENDDO 
                           lchange = .true. 
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
!     ----run transformation 'run'                                      
!                                                                       
               ELSEIF (str_comp (befehl, 'run ', 1, lbef, 4) ) then 
                  IF (lchange) then 
                     CALL tran_setup 
                  ENDIF 
                  lchange = .false. 
                  CALL tran_op 
!                                                                       
!     ----Select which atoms are copied to their image 'sele'           
!                                                                       
               ELSEIF (str_comp (befehl, 'sele', 3, lbef, 4) ) then 
                  CALL atom_select (zeile, lp, 0, TRAN_MAXSCAT, tran_latom, &
                  tran_sel_atom, lold,.true.)
!                 ier_num = - 6 
!                 ier_typ = ER_COMM 
!                 CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
!                 IF (ier_num.eq.0) then 
!                    CALL get_iscat (ianz, cpara, lpara, werte, maxw,   &
!                    lold)                                              
!                    IF (ier_num.eq.0) then 
!                       IF (werte (1) .eq. - 1) then 
!                          DO i = 0, cr_nscat 
!                          tran_latom (i) = .true. 
!                          ENDDO 
!                       ELSE 
!                          DO i = 1, ianz 
!                          tran_latom (nint (werte (i) ) ) = .true. 
!                          ENDDO 
!                       ENDIF 
!                    ENDIF 
!                 ENDIF 
!                                                                       
!     ----Set parameters for various commands 'set'                     
!                                                                       
               ELSEIF (str_comp (befehl, 'set', 3, lbef, 3) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (ianz.gt.0) then 
                        IF (str_comp (cpara (1) , 'deltahkl', 1, lbef,  &
                        8) ) then                                       
                           cpara (1) = '0.0' 
                           lpara (1) = 3 
                           CALL ber_params (ianz, cpara, lpara, werte,  &
                           maxw)                                        
                           IF (ier_num.eq.0) then 
                              tran_deltahkl = werte (2) 
                           ENDIF 
                        ENDIF 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ENDIF 
!                                                                       
!     ----show current parameters 'show'                                
!                                                                       
               ELSEIF (str_comp (befehl, 'show', 2, lbef, 4) ) then 
                  IF (lchange) then 
                     CALL tran_setup 
                  ENDIF 
                  lchange = .false. 
                  WRITE (output_io, 3000) 
                  WRITE (output_io, 3010) ( (tran_g (i, j), j = 1, 3),  &
                  i = 1, 3)                                             
                  WRITE (output_io, 3020) ( (tran_gi (i, j), j = 1, 3), &
                  i = 1, 3)                                             
                  WRITE (output_io, 3030) ( (tran_f (i, j), j = 1, 4),  &
                  i = 1, 3)                                             
                  WRITE (output_io, 3040) ( (tran_fi (i, j), j = 1, 4), &
                  i = 1, 3)                                             
                  WRITE (output_io, 3050) ( (tran_f (i, j), j = 1, 3),  &
                  i = 1, 3)                                             
                  WRITE (output_io, 3060) ( (tran_fi (i, j), j = 1, 3), &
                  i = 1, 3)                                             
                  WRITE (output_io, 3061) ( (tran_g (i, j), j = 1, 3),  &
                  i = 1, 3)                                             
                  WRITE (output_io, 3062) ( (tran_gi (i, j), j = 1, 3), &
                  i = 1, 3)                                             
!                                                                       
                  IF (tran_det.ne.0) then 
                     WRITE (output_io, 3100) tran_det, 1. / tran_det 
                  ELSE 
                     WRITE (output_io, 3200) 
                  ENDIF 
!                                                                       
                  WRITE (output_io, 3090) 
                  WRITE (output_io, 3091) 
                  DO i = 0, cr_nscat 
                  IF (tran_latom (i) ) then 
                     WRITE (output_io, 3092) i, cr_at_lis (i) 
                  ENDIF 
                  ENDDO 
                  IF (tran_end.eq. - 1) then 
                     WRITE (output_io, 3080) 
                  ELSE 
                     WRITE (output_io, 3081) tran_start, tran_end 
                  ENDIF 
!                                                                       
!------- -Operating System Kommandos 'syst'                             
!                                                                       
               ELSEIF (str_comp (befehl, 'syst', 2, lbef, 4) ) then 
                  IF (zeile.ne.' ') then 
                     CALL do_operating (zeile (1:lp), lp) 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
!                                                                       
!------  -----waiting for user input                                    
!                                                                       
               ELSEIF (str_comp (befehl, 'wait', 3, lbef, 4) ) then 
                  CALL do_input (zeile, lp) 
!                                                                       
!     ----x(new) in terms of old axes                                   
!                                                                       
               ELSEIF (str_comp (befehl, 'xnew', 2, lbef, 4)            &
               .or.str_comp (befehl, 'asnew', 3, lbef, 5) ) then        
                  ier_num = - 6 
                  ier_typ = ER_COMM 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0.and.ianz.eq.3) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        DO j = 1, 3 
                        tran_f (1, j) = werte (j) 
                        ENDDO 
                        tran_inp = TRAN_INP_F 
                        lchange = .true. 
                     ENDIF 
                  ENDIF 
!                                                                       
!     ----x(old) in terms of new axes                                   
!                                                                       
               ELSEIF (str_comp (befehl, 'xold', 2, lbef, 4)            &
               .or.str_comp (befehl, 'asold', 3, lbef, 5) ) then        
                  ier_num = - 6 
                  ier_typ = ER_COMM 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0.and.ianz.eq.3) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        DO j = 1, 3 
                        tran_fi (1, j) = werte (j) 
                        ENDDO 
                        tran_inp = TRAN_INP_FI 
                        lchange = .true. 
                     ENDIF 
                  ENDIF 
!                                                                       
!     ----y(new) in terms of old axes                                   
!                                                                       
               ELSEIF (str_comp (befehl, 'ynew', 2, lbef, 4)            &
               .or.str_comp (befehl, 'bsnew', 3, lbef, 5) ) then        
                  ier_num = - 6 
                  ier_typ = ER_COMM 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0.and.ianz.eq.3) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        DO j = 1, 3 
                        tran_f (2, j) = werte (j) 
                        ENDDO 
                        tran_inp = TRAN_INP_F 
                        lchange = .true. 
                     ENDIF 
                  ENDIF 
!                                                                       
!     ----y(old) in terms of new axes                                   
!                                                                       
               ELSEIF (str_comp (befehl, 'yold', 2, lbef, 4)            &
               .or.str_comp (befehl, 'bsold', 3, lbef, 5) ) then        
                  ier_num = - 6 
                  ier_typ = ER_COMM 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0.and.ianz.eq.3) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        DO j = 1, 3 
                        tran_fi (2, j) = werte (j) 
                        ENDDO 
                        tran_inp = TRAN_INP_FI 
                        lchange = .true. 
                     ENDIF 
                  ENDIF 
!                                                                       
!     ----z(new) in terms of old axes                                   
!                                                                       
               ELSEIF (str_comp (befehl, 'znew', 2, lbef, 4)            &
               .or.str_comp (befehl, 'csnew', 3, lbef, 5) ) then        
                  ier_num = - 6 
                  ier_typ = ER_COMM 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0.and.ianz.eq.3) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        DO j = 1, 3 
                        tran_f (3, j) = werte (j) 
                        ENDDO 
                        tran_inp = TRAN_INP_F 
                        lchange = .true. 
                     ENDIF 
                  ENDIF 
!                                                                       
!     ----z(old) in terms of new axes                                   
!                                                                       
               ELSEIF (str_comp (befehl, 'zold', 2, lbef, 4)            &
               .or.str_comp (befehl, 'csold', 3, lbef, 5) ) then        
                  ier_num = - 6 
                  ier_typ = ER_COMM 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0.and.ianz.eq.3) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        DO j = 1, 3 
                        tran_fi (3, j) = werte (j) 
                        ENDDO 
                        tran_inp = TRAN_INP_FI 
                        lchange = .true. 
                     ENDIF 
                  ENDIF 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_COMM 
               ENDIF 
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
      ENDDO 
!                                                                       
 3000 FORMAT    (20x,'    Unit cell transformations '/                  &
     &                  20x,' ==============================='//)       
 3010 FORMAT    (                                                       &
     &           ' ( a(new) ) = ( ',2(F9.5,','),f9.5,' )   ( a(old) )'/ &
     &           ' ( b(new) ) = ( ',2(F9.5,','),f9.5,' ) * ( b(old) )'/ &
     &           ' ( c(new) ) = ( ',2(F9.5,','),f9.5,' )   ( c(old) )'/)
 3020 FORMAT    (                                                       &
     &           ' ( a(old) ) = ( ',2(F9.5,','),f9.5,' )   ( a(new) )'/ &
     &           ' ( b(old) ) = ( ',2(F9.5,','),f9.5,' ) * ( b(new) )'/ &
     &           ' ( c(old) ) = ( ',2(F9.5,','),f9.5,' )   ( c(new) )'/)
 3030 FORMAT    (                                                       &
     &           ' ( x(new) ) = ( ',2(F9.5,','),f9.5,' )   ( x(old) )', &
     &           '   (',f9.5,')'/                                       &
     &           ' ( y(new) ) = ( ',2(F9.5,','),f9.5,' ) * ( y(old) )', &
     &           ' + (',f9.5,')'/                                       &
     &           ' ( z(new) ) = ( ',2(F9.5,','),f9.5,' )   ( z(old) )', &
     &           '   (',f9.5,')'/)                                      
 3040 FORMAT    (                                                       &
     &           ' ( x(old) ) = ( ',2(F9.5,','),f9.5,' )   ( x(new) )', &
     &           '   (',f9.5,')'/                                       &
     &           ' ( y(old) ) = ( ',2(F9.5,','),f9.5,' ) * ( y(new) )', &
     &           ' + (',f9.5,')'/                                       &
     &           ' ( z(old) ) = ( ',2(F9.5,','),f9.5,' )   ( z(new) )', &
     &           '   (',f9.5,')'/)                                      
 3050 FORMAT    (                                                       &
     &           ' ( a*(new)) = ( ',2(F9.5,','),f9.5,' )   ( a*(old))'/ &
     &           ' ( b*(new)) = ( ',2(F9.5,','),f9.5,' ) * ( b*(old))'/ &
     &           ' ( c*(new)) = ( ',2(F9.5,','),f9.5,' )   ( c*(old))'/)
 3060 FORMAT    (                                                       &
     &           ' ( a*(old)) = ( ',2(F9.5,','),f9.5,' )   ( a*(new))'/ &
     &           ' ( b*(old)) = ( ',2(F9.5,','),f9.5,' ) * ( b*(new))'/ &
     &           ' ( c*(old)) = ( ',2(F9.5,','),f9.5,' )   ( c*(new))'/)
 3061 FORMAT    (                                                       &
     &           ' ( h(new) ) = ( ',2(F9.5,','),f9.5,' )   ( h(old) )'/ &
     &           ' ( k(new) ) = ( ',2(F9.5,','),f9.5,' ) * ( k(old) )'/ &
     &           ' ( l(new) ) = ( ',2(F9.5,','),f9.5,' )   ( l(old) )'/)
 3062 FORMAT    (                                                       &
     &           ' ( h(old) ) = ( ',2(F9.5,','),f9.5,' )   ( h(new) )'/ &
     &           ' ( k(old) ) = ( ',2(F9.5,','),f9.5,' ) * ( k(new) )'/ &
     &           ' ( l(old) ) = ( ',2(F9.5,','),f9.5,' )   ( l(new) )'/)
 3100 FORMAT    (                                                       &
     &           '   V(new)   = ',G12.5E3,' * V(old)'/                     &
     &           '   V(old)   = ',G12.5E3,' * V(new)'/)                    
 3200 FORMAT    (                                                       &
     &           '   Determinant of transformation is zero'/)           
 3080 FORMAT    (' Range of atoms from to    : All atoms included') 
 3081 FORMAT    (' Range of atoms from to    : ',2(2x,i9)) 
 3090 FORMAT    (' selected atoms    :') 
 3091 FORMAT    ('                      type name') 
 3092 FORMAT    (16x,2x,i8,1x,a4) 
!                                                                       
      END SUBROUTINE transform                      
!*****7*****************************************************************
      SUBROUTINE tran_setup 
!-                                                                      
!     Performs the generalized symmetry operation                       
!     See Sands, D.E. Vectors and Tensors in Crystallography Chapt. 4.7 
!+                                                                      
      USE config_mod 
      USE tensors_mod
      USE transfrm_mod 
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER i, j 
      REAL a (3, 3) 
!                                                                       
!     initialize fourth columns and rows                                
!                                                                       
      DO i = 1, 4 
      tran_g (i, 4) = 0.0 
      tran_g (4, i) = 0.0 
      tran_gi (i, 4) = 0.0 
      tran_gi (4, i) = 0.0 
      tran_f (i, 4) = 0.0 
      tran_f (4, i) = 0.0 
      tran_fi (i, 4) = 0.0 
      tran_fi (4, i) = 0.0 
      ENDDO 
      tran_f (4, 4) = 1.0 
      tran_fi (4, 4) = 1.0 
!                                                                       
!     Matrix g was defined                                              
!                                                                       
      IF (tran_inp.eq.TRAN_INP_G) then 
         DO i = 1, 4 
         DO j = 1, 4 
         tran_gi (i, j) = tran_g (i, j) 
         ENDDO 
         ENDDO 
         CALL invmat4 (tran_gi) 
         DO i = 1, 3 
         DO j = 1, 3 
         tran_f (j, i) = tran_gi (i, j) 
         tran_fi (j, i) = tran_g (i, j) 
         ENDDO 
         ENDDO 
!                                                                       
!     Matrix gi was defined                                             
!                                                                       
      ELSEIF (tran_inp.eq.TRAN_INP_GI) then 
         DO i = 1, 4 
         DO j = 1, 4 
         tran_g (i, j) = tran_gi (i, j) 
         ENDDO 
         ENDDO 
         CALL invmat4 (tran_g) 
         DO i = 1, 3 
         DO j = 1, 3 
         tran_f (j, i) = tran_gi (i, j) 
         tran_fi (j, i) = tran_g (i, j) 
         ENDDO 
         ENDDO 
!                                                                       
!     Matrix f was defined                                              
!                                                                       
      ELSEIF (tran_inp.eq.TRAN_INP_F) then 
         DO i = 1, 4 
         DO j = 1, 4 
         tran_fi (i, j) = tran_f (i, j) 
         ENDDO 
         ENDDO 
         CALL invmat4 (tran_fi) 
         DO i = 1, 3 
         DO j = 1, 3 
         tran_g (j, i) = tran_fi (i, j) 
         tran_gi (j, i) = tran_f (i, j) 
         ENDDO 
         ENDDO 
!                                                                       
!     Matrix fi was defined                                             
!                                                                       
      ELSEIF (tran_inp.eq.TRAN_INP_FI) then 
         DO i = 1, 4 
         DO j = 1, 4 
         tran_f (i, j) = tran_fi (i, j) 
         ENDDO 
         ENDDO 
         CALL invmat4 (tran_f) 
         DO i = 1, 3 
         DO j = 1, 3 
         tran_g (j, i) = tran_fi (i, j) 
         tran_gi (j, i) = tran_f (i, j) 
         ENDDO 
         ENDDO 
      ENDIF 
!                                                                       
!     Origin was defined in terms of the old coordinates                
!                                                                       
      IF (tran_oold) then 
         j = 4 
         DO i = 1, 3 
         tran_fi (i, 4) = tran_orig (i) 
         ENDDO 
         DO i = 1, 3 
         DO j = 1, 3 
         tran_f (i, 4) = tran_f (i, 4) - tran_f (i, j) * tran_fi (j, 4) 
         ENDDO 
         ENDDO 
      ELSE 
         j = 4 
         DO i = 1, 3 
         tran_f (i, 4) = tran_orig (i) 
         ENDDO 
         DO i = 1, 3 
         DO j = 1, 3 
         tran_fi (i, 4) = tran_fi (i, 4) - tran_fi (i, j) * tran_f (j,  &
         4)                                                             
         ENDDO 
         ENDDO 
      ENDIF 
!                                                                       
      DO i = 1, 3 
      DO j = 1, 3 
      a (i, j) = tran_g (i, j) 
      ENDDO 
      ENDDO 
      tran_det = a (1, 1) * (a (2, 2) * a (3, 3) - a (2, 3) * a (3, 2) )&
      + a (2, 1) * (a (3, 2) * a (1, 3) - a (1, 2) * a (3, 3) ) + a (3, &
      1) * (a (1, 2) * a (2, 3) - a (1, 3) * a (2, 2) )                 
!                                                                       
      END SUBROUTINE tran_setup                     
!*****7*****************************************************************
      SUBROUTINE tran_op 
!-                                                                      
!     Performs the actual transformation operation.                     
!     All atoms of the structure, all symmetry elements are transformed 
!     into the new cell.                                                
!+                                                                      
      USE config_mod 
      USE crystal_mod 
      USE metric_mod
      USE spcgr_apply, ONLY: setup_lattice
      USE update_cr_dim_mod
      USE transfrm_mod 
      USE trafo_mod
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER i, j 
      INTEGER i1, i2, i3 
      INTEGER i_start, i_end 
      LOGICAL lspace 
      LOGICAL lout 
      REAL usym (4), ures (4) 
      REAL werte (5) 
      REAL u (3), v (3), w (3) 
!                                                                       
!     REAL do_bang 
!     REAL skalpro 
!                                                                       
      DATA usym / 0.0, 0.0, 0.0, 1.0 / 
      DATA werte / 0.0, 0.0, 0.0, 0.0, 0.0 / 
!                                                                       
!     Set the appropriate starting end ending number for the atoms      
!                                                                       
      i_start = tran_start 
      i_end = tran_end 
      IF (tran_end.eq. - 1) i_end = cr_natoms 
!                                                                       
!     Apply transformation operation to all atoms within selected range 
!                                                                       
      DO i = i_start, i_end 
!                                                                       
!     --Select atom if:                                                 
!       type has been selected                                          
!                                                                       
      IF (tran_latom (cr_iscat (i) ) ) then 
!                                                                       
!     ----Copy atom to temporary place                                  
!                                                                       
         DO j = 1, 3 
         usym (j) = cr_pos (j, i) 
         ENDDO 
!                                                                       
!-----      ----Apply symmetry operation                                
!                                                                       
         usym (4) = 1.0 
         CALL trans (usym, tran_f, ures, 4) 
!                                                                       
!     ----Replace original atom by its image                            
!                                                                       
         DO j = 1, 3 
         cr_pos (j, i) = ures (j) 
         ENDDO 
      ENDIF 
      ENDDO 
!                                                                       
!     If all atoms of the crystal have been selected, calculate new     
!     lattice constants and symmetry operations                         
!                                                                       
      IF (i_start.eq.1.and.i_end.eq.cr_natoms) then 
         v (1) = 0.0 
         v (2) = 0.0 
         v (3) = 0.0 
         lspace = .true. 
         DO i = 1, 3 
         DO j = 1, 3 
         usym (j) = tran_g (i, j) 
         ENDDO 
         cr_a0 (i) = sqrt (skalpro (usym, usym, cr_gten) ) 
!                                                                       
         i1 = mod (i - 1, 3) + 1 
         i2 = mod (i, 3) + 1 
         i3 = mod (i + 1, 3) + 1 
         DO j = 1, 3 
         u (j) = tran_g (i1, j) 
         w (j) = tran_g (i2, j) 
         ENDDO 
         cr_win (i3) = do_bang (lspace, u, v, w) 
         ENDDO 
         CALL tran_sym 
         cr_spcgr = 'P1' 
         cr_spcgrno = 1 
         spcgr_ianz = 0 
         spcgr_para = 1 
         lout = .true. 
         CALL setup_lattice (cr_a0, cr_ar, cr_eps, cr_gten, cr_reps,    &
         cr_rten, cr_win, cr_wrez, cr_v, cr_vr, lout, cr_gmat, cr_fmat, &
         cr_cartesian)                                                  
         CALL update_cr_dim 
      ENDIF 
!                                                                       
      END SUBROUTINE tran_op                        
!*****7*****************************************************************
      SUBROUTINE tran_sym 
!-                                                                      
!     Transforms the symmetry operations                                
!+                                                                      
      USE config_mod 
      USE crystal_mod 
      USE generate_mod 
      USE gen_add_mod 
      USE sym_add_mod 
      USE tensors_mod
      USE transfrm_mod 
      USE unitcell_mod 
!                                                                       
      USE errlist_mod 
      USE prompt_mod 
!
      IMPLICIT none 
!                                                                       
      INTEGER :: i, j, k, n, nn=0, jp, l 
      LOGICAL lequal 
!                                                                       
      REAL mat (4, 4), arr (4, 4) 
      REAL eps 
!                                                                       
!                                                                       
      DATA eps / 0.001 / 
!                                                                       
!------ The original symmetry elements are going to be stored as        
!       additional                                                      
!     generators. The original additional generators are shifted up.    
!                                                                       
      IF (gen_add_n + generspcgr (0, cr_spcgrno) .gt.gen_add_MAX) then 
         ier_typ = ER_APPL 
         ier_num = - 73 
         RETURN 
      ENDIF 
!                                                                       
!------ Shift original additional generators by generspcgr(0,cr_spcgrno)
!       entries                                                         
!                                                                       
      k = generspcgr (0, cr_spcgrno) 
      DO n = gen_add_n, 1, - 1 
      nn = n + k 
      DO i = 1, 4 
      DO j = 1, 4 
      gen_add (i, j, nn) = gen_add (i, j, n) 
      ENDDO 
      ENDDO 
      gen_add_power (nn) = gen_add_power (n) 
      ENDDO 
!                                                                       
!     copy space group generators into additional generators            
!                                                                       
      DO n = 1, generspcgr (0, cr_spcgrno) 
      IF (gen_sta.eq.GEN_SYMM) then 
         nn = generspcgr (n, cr_spcgrno) 
      ELSEIF (gen_sta.eq.GEN_CENTER) then 
         nn = generspcgr_center (n, cr_spcgrno) 
      ENDIF 
      DO i = 1, 4 
      DO j = 1, 4 
      gen_add (i, j, n) = generators (i, j, nn) 
      ENDDO 
      ENDDO 
      gen_add_power (n) = generpower (nn) 
      ENDDO 
      gen_add_n = gen_add_n + generspcgr (0, cr_spcgrno) 
!                                                                       
!     transform the additional generators                               
!                                                                       
      DO n = 1, gen_add_n 
      DO i = 1, 4 
      DO j = 1, 4 
      mat (i, j) = gen_add (i, j, n) 
      ENDDO 
      ENDDO 
      CALL matmul4 (arr, tran_f, mat) 
      CALL matmul4 (mat, arr, tran_fi) 
      DO i = 1, 4 
      DO j = 1, 4 
      gen_add (i, j, n) = mat (i, j) 
      ENDDO 
      ENDDO 
      ENDDO 
!                                                                       
!     transform the additional symmetry elements                        
!                                                                       
      DO n = 1, sym_add_n 
      DO i = 1, 4 
      DO j = 1, 4 
      mat (i, j) = sym_add (i, j, n) 
      ENDDO 
      ENDDO 
      CALL matmul4 (arr, tran_f, mat) 
      CALL matmul4 (mat, arr, tran_fi) 
      DO i = 1, 4 
      DO j = 1, 4 
      sym_add (i, j, n) = mat (i, j) 
      ENDDO 
      ENDDO 
      ENDDO 
!                                                                       
!     Check transformation of primitive translations, if it results in  
!     non-integer translation, add as additional generator.             
!                                                                       
      k = gen_add_n 
      DO n = 1, 3 
      DO i = 1, 4 
      DO j = 1, 4 
      mat (i, j) = 0.0 
      ENDDO 
      mat (i, i) = 1.0 
      ENDDO 
      mat (n, 4) = 1.0 
      CALL matmul4 (arr, tran_f, mat) 
      CALL matmul4 (mat, arr, tran_fi) 
      IF (abs (amod (mat (1, 4), 1.) ) .gt.eps.or.abs (amod (mat (2, 4),&
      1.) ) .gt.eps.or.abs (amod (mat (3, 4), 1.) ) .gt.eps) then       
         IF (gen_add_n + 1.gt.gen_add_MAX) then 
            ier_typ = ER_APPL 
            ier_num = - 73 
            RETURN 
         ENDIF 
!                                                                       
!     ----Determine power of generator                                  
!                                                                       
!                                                                       
!     ----After this loop j holds the power,which transforms the        
!         generator into a primitive translation                        
!                                                                       
!                                                                       
!     ----This loop transforms the translation into 0<=t<1              
!                                                                       
         DO i = 1, 3 
         mat (i, 4) = amod (mat (i, 4), 1.0) 
         IF (mat (i, 4) .lt.0.0) mat (i, 4) = mat (i, 4) + 1.0 
         ENDDO 
         jp = 1 
         DO while (abs (jp * mat (1, 4) - nint (jp * mat (1, 4) ) )     &
         .gt.eps.or.abs (jp * mat (2, 4) - nint (jp * mat (2, 4) ) )    &
         .gt.eps.or.abs (jp * mat (3, 4) - nint (jp * mat (3, 4) ) )    &
         .gt.eps)                                                       
         jp = jp + 1 
         ENDDO 
!                                                                       
!     ----Compare to previous generators originating from primitive     
!         translations. If identical, skip this one.                    
!                                                                       
         lequal = .true. 
         DO i = k + 1, gen_add_n 
         DO j = 1, 3 
         lequal = lequal.and. (abs (gen_add (j, 4, i) - mat (j, 4) )    &
         .lt.eps)                                                       
         ENDDO 
         ENDDO 
         IF (.not.lequal.or.k.eq.gen_add_n) then 
            gen_add_n = gen_add_n + 1 
            gen_add_power (gen_add_n) = 1 
            IF (jp.gt.1) jp = jp - 1 
            gen_add_power (gen_add_n) = jp 
            DO i = 1, 4 
            DO j = 1, 4 
            gen_add (i, j, gen_add_n) = mat (i, j) 
            ENDDO 
            ENDDO 
         ENDIF 
      ENDIF 
      ENDDO 
!                                                                       
!     Shift the generators resulting from the primitive translations    
!     up front                                                          
!                                                                       
      n = gen_add_n - k 
      DO i = 1, n 
      DO j = 1, 4 
      DO k = 1, 4 
      mat (j, k) = gen_add (j, k, gen_add_n) 
      ENDDO 
      ENDDO 
      jp = gen_add_power (gen_add_n) 
      DO l = gen_add_n, 2, - 1 
      DO j = 1, 4 
      DO k = 1, 4 
      gen_add (j, k, l) = gen_add (j, k, l - 1) 
      ENDDO 
      ENDDO 
      gen_add_power (l) = gen_add_power (l - 1) 
      ENDDO 
      DO j = 1, 4 
      DO k = 1, 4 
      gen_add (j, k, 1) = mat (j, k) 
      ENDDO 
      ENDDO 
      gen_add_power (1) = jp 
      ENDDO 
!                                                                       
      END SUBROUTINE tran_sym                       
!*****7*****************************************************************
      SUBROUTINE tran_hkl (infile, infile_l, matrix) 
!-                                                                      
!     Transforms a list of reflections, read from file.                 
!+                                                                      
      USE config_mod 
      USE transfrm_mod 
      USE trafo_mod
!                                                                       
      USE errlist_mod 
      USE param_mod 
      USE prompt_mod 
!
      IMPLICIT none 
!                                                                       
      INTEGER ird, iwr, irs 
      INTEGER i, j 
      LOGICAL lread, ltest 
!                                                                       
      CHARACTER(20) inten 
      CHARACTER ( * ) infile 
      CHARACTER(1024) outfile 
      CHARACTER(1024) restfile 
      INTEGER infile_l
      INTEGER hkl (3) 
      REAL usym (4), ures (4), utest (3) 
      REAL matrix (4, 4) 
!                                                                       
      DATA ird, iwr, irs / 7, 8, 9 / 
!                                                                       
!     Open input file, and the two output files                         
!                                                                       
      lread = .true. 
      CALL oeffne (ird, infile, 'unknown', lread) 
      IF (ier_num.ne.0) then 
         CLOSE (ird) 
         RETURN 
      ENDIF 
      lread = .false. 
      outfile = infile (1:infile_l) //'.trans' 
      CALL oeffne (iwr, outfile, 'unknown', lread) 
      IF (ier_num.ne.0) then 
         CLOSE (ird) 
         CLOSE (iwr) 
         RETURN 
      ENDIF 
      lread = .false. 
      restfile = infile (1:infile_l) //'.rest' 
      CALL oeffne (irs, restfile, 'unknown', lread) 
      IF (ier_num.ne.0) then 
         CLOSE (ird) 
         CLOSE (iwr) 
         CLOSE (irs) 
         RETURN 
      ENDIF 
!                                                                       
!-----      Loop over all reflections in input file                     
!                                                                       
 1000 CONTINUE 
      READ (ird, 2000, end = 999) hkl, inten 
      DO i = 1, 3 
      usym (i) = float (hkl (i) ) 
      ENDDO 
      usym (4) = 0.0 
!                                                                       
!-----      -- Apply transformation operation                           
!                                                                       
      CALL trans (usym, matrix, ures, 4) 
!                                                                       
!     -- Write new integer reflections to output file, non-integer      
!        reflections to restfile                                        
!                                                                       
      ltest = .true. 
      DO i = 1, 3 
      utest (i) = abs (float (nint (ures (i) ) ) - ures (i) ) 
      ltest = ltest.and. (utest (i) .lt.tran_deltahkl) 
      ENDDO 
      IF (ltest) then 
         WRITE (iwr, 2000) (nint (ures (j) ), j = 1, 3), inten 
      ELSE 
         WRITE (irs, 3000) (ures (j), j = 1, 3), inten, (nint (usym (j) &
         ), j = 1, 3)                                                   
      ENDIF 
      GOTO 1000 
  999 CONTINUE 
      CLOSE (ird) 
      CLOSE (iwr) 
      CLOSE (irs) 
!                                                                       
!                                                                       
 2000 FORMAT    (3i4,a20) 
 3000 FORMAT    (3f8.3,1x,a20,1x,3i4) 
      END SUBROUTINE tran_hkl                       
END MODULE transform_menu
