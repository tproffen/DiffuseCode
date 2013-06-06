!+                                                                      
!     Generalized symmetry operations:                                  
!                                                                       
!     Rotations around a general axis by a general angle                
!     Mirror operations with free choice of normal and glide            
!     component.                                                        
!                                                                       
!*****7*****************************************************************
      SUBROUTINE symm 
!-                                                                      
!     Main menu for generalized symmetry operations                     
!+                                                                      
      USE config_mod 
      USE allocate_appl_mod
      USE crystal_mod 
      USE modify_mod
      USE molecule_mod 
      USE structur, ONLY: update_cr_dim
      USE symm_mod 
      IMPLICIT none 
!                                                                       
       
      include'doact.inc' 
      include'errlist.inc' 
      include'learn.inc' 
      include'macro.inc' 
      include'prompt.inc' 
!                                                                       
      INTEGER, PARAMETER :: MIN_PARA = 20
      INTEGER            :: maxw 
      LOGICAL lnew, lold 
!                                                                       
      PARAMETER (lnew = .true., lold = .false.) 
!                                                                       
      CHARACTER(LEN=1024), DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: cpara
      REAL               , DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: werte
      INTEGER            , DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: lpara
!
      CHARACTER(5) befehl 
      CHARACTER(50) prom 
      CHARACTER(1024) line, zeile
      INTEGER lp, length, lbef 
      INTEGER indxg, ianz, i 
      INTEGER indxc 
      INTEGER         :: nscat
      LOGICAL lend, lspace 
      LOGICAL l_need_setup 
      LOGICAL lselect 
      REAL hkl (3) 
!                                                                       
      INTEGER len_str 
      LOGICAL str_comp 
!                                                                       
!                                                                       
      DATA l_need_setup / .true. / 
!                                                                       
      maxw = MAX(MIN_PARA,MAXSCAT+1)
      lend = .false. 
      CALL no_error 
!
      IF( cr_nscat > SYM_MAXSCAT .or. mole_num_type > SYM_MAXSCAT) THEN
         nscat = max ( cr_nscat, mole_num_type)
         CALL alloc_symmetry ( nscat )
         IF ( ier_num < 0 ) THEN
            RETURN
         ENDIF
      ENDIF
!                                                                       
      DO while (.not.lend) 
      prom = prompt (1:len_str (prompt) ) //'/symm' 
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
!     ----angle of rotation 'angle'                                     
!                                                                       
               ELSEIF (str_comp (befehl, 'angle', 2, lbef, 5) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0.and.ianz.eq.1) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        sym_angle = werte (1) 
                        l_need_setup = .true. 
                     ENDIF 
                  ENDIF 
!                                                                       
!     ----list asymmetric unit 'asym'                                   
!                                                                       
               ELSEIF (str_comp (befehl, 'asym', 2, lbef, 4) ) then 
                  CALL show_asym 
!                                                                       
!     ----calculate a single symmetry operation                         
!                                                                       
               ELSEIF (str_comp (befehl, 'calc', 2, lbef, 4) ) then 
                  IF (l_need_setup) then 
                     CALL symm_setup 
                  ENDIF 
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
                           lspace = .true. 
                        ELSEIF (str_comp (cpara (4) , 'r', 1, lpara (4) &
                        , 1) ) then                                     
                           lspace = .false. 
                        ELSE 
                           ier_num = - 6 
                           ier_typ = ER_COMM 
                        ENDIF 
                        IF (ier_num.eq.0) then 
                           IF (sym_power_mult) then 
                              CALL symm_ca_mult (hkl, lspace) 
                           ELSE 
                              CALL symm_ca_single (hkl, lspace, .true.) 
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
!     ----Work on domains 'domain'                                      
!                                                                       
               ELSEIF (str_comp (befehl, 'domain', 2, lbef, 6) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (str_comp (cpara (1) , 'select', 2, lpara (1) , &
                     6) .or.str_comp (cpara (1) , 'deselect', 2, lpara (&
                     1) , 8) ) then                                     
                        line = ' ' 
                        indxc = index (zeile, ',') 
                        line = zeile (indxc + 1:lp) 
                        lp = len_str (line) 
!                                                                       
                        lselect = str_comp (cpara (1) , 'select', 2,    &
                        lpara (1) , 6) .or.str_comp (cpara (1) ,        &
                        'deselect', 2, lpara (1) , 8)                   
!                                                                       
                        CALL mole_select (line, lp, 0, SYM_MAXSCAT,     &
                             sym_latom, sym_sel_atom, lold, lselect)
                        sym_sel_mode = SYM_RUN_DOMAIN 
                     ELSEIF (str_comp (cpara (1) , 'include', 3, lpara (&
                     1) , 7) ) then                                     
                        CALL del_params (1, ianz, cpara, lpara, maxw) 
                        IF (ianz.eq.2) then 
                           CALL ber_params (ianz, cpara, lpara, werte,  &
                           maxw)                                        
                           IF (ier_num.eq.0) then 
                              sym_sel_atom = .false. 
                              sym_start = nint (werte (1) ) 
                              sym_end = nint (werte (2) ) 
                           ENDIF 
                        ELSEIF (ianz.eq.1) then 
                           IF (str_comp (cpara (1) , 'all', 1, lpara (1)&
                           , 3) ) then                                  
                              sym_sel_atom = .false. 
                              sym_start = 1 
                              sym_end = - 1 
                           ELSE 
                              ier_num = - 6 
                              ier_typ = ER_COMM 
                           ENDIF 
                        ELSE 
                           ier_num = - 6 
                           ier_typ = ER_COMM 
                        ENDIF 
                     ELSEIF (str_comp (cpara (1) , 'atoms', 2, lpara (1)&
                     , 4) ) then                                        
                        sym_dom_mode_atom = str_comp (cpara (2) ,       &
                        'apply', 2, lpara (2) , 5)                      
                     ELSEIF (str_comp (cpara (1) , 'shape', 2, lpara (2)&
                     , 5) ) then                                        
                        sym_dom_mode_shape = str_comp (cpara (2) ,      &
                        'apply', 2, lpara (2) , 5)                      
                     ENDIF 
                  ENDIF 
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
                     CALL do_hel ('discus symm '//zeile, lp) 
                  ENDIF 
!                                                                       
!     ----Select the reciprocal space direction of the symmetry         
!         axis 'hkl'                                                    
!                                                                       
               ELSEIF (str_comp (befehl, 'hkl ', 2, lbef, 4) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0.and.ianz.eq.3) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        DO i = 1, 3 
                        sym_hkl (i) = werte (i) 
                        ENDDO 
                        CALL trans (sym_hkl, cr_rten, sym_uvw, 3) 
                        l_need_setup = .true. 
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
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (ianz.eq.2) then 
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxw)                                           
                        IF (ier_num.eq.0) then 
                           sym_sel_atom = .true. 
                           sym_start = nint (werte (1) ) 
                           sym_end = nint (werte (2) ) 
                           sym_incl = 'list' 
                        ENDIF 
                     ELSEIF (ianz.eq.1) then 
                        IF (str_comp (cpara (1) , 'all', 1, lpara (1) , &
                        3) ) then                                       
                           sym_sel_atom = .true. 
                           sym_start = 1 
                           sym_end = - 1 
                           sym_incl = 'all ' 
                        ELSEIF (str_comp (cpara (1) , 'env', 1, lpara ( &
                        1) , 3) ) then                                  
                           sym_sel_atom = .true. 
                           sym_start = 1 
                           sym_end = - 1 
                           sym_incl = 'env ' 
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
!     ----Select range of molecules within crystal to be included       
!         'mincl'                                                       
!                                                                       
               ELSEIF (str_comp (befehl, 'mincl', 3, lbef, 5)           &
               .or.str_comp (befehl, 'oincl', 3, lbef, 5) ) then        
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (ianz.eq.2) then 
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxw)                                           
                        IF (ier_num.eq.0) then 
                           sym_sel_atom = .false. 
                           sym_start = nint (werte (1) ) 
                           sym_end = nint (werte (2) ) 
                        ENDIF 
                     ELSEIF (ianz.eq.1) then 
                        IF (str_comp (cpara (1) , 'all', 1, lpara (1) , &
                        3) ) then                                       
                           sym_sel_atom = .false. 
                           sym_start = 1 
                           sym_end = - 1 
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
!     ----Select the mode of symmetry operation 'mode'                  
!                                                                       
               ELSEIF (str_comp (befehl, 'mode', 1, lbef, 4) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (ianz.eq.1.or.ianz.eq.2) then 
                        IF (str_comp (cpara (1) , 'copy', 1, lpara (1) ,&
                        4) ) then                                       
                           sym_mode = .true. 
                           l_need_setup = .true. 
                        ELSEIF (str_comp (cpara (1) , 'repl', 1, lpara (&
                        1) , 4) ) then                                  
                           sym_mode = .false. 
                           l_need_setup = .true. 
                        ELSE 
                           ier_num = - 6 
                           ier_typ = ER_COMM 
                        ENDIF 
                     ENDIF 
!                                                                       
                     IF (ier_num.eq.0.and.ianz.eq.2) then 
                        IF (str_comp (cpara (2) , 'new', 1, lpara (2) , &
                        3) ) then                                       
                           sym_new = .true. 
                           l_need_setup = .true. 
                        ELSEIF (str_comp (cpara (2) , 'old', 1, lpara ( &
                        1) , 3) ) then                                  
                           sym_new = .false. 
                           l_need_setup = .true. 
                        ELSE 
                           ier_num = - 6 
                           ier_typ = ER_COMM 
                        ENDIF 
                     ENDIF 
                  ENDIF 
!                                                                       
!     ----Select/deselect molecules                                     
!                                                                       
               ELSEIF (str_comp (befehl, 'msel', 2, lbef, 4)            &
               .or.str_comp (befehl, 'mdes', 2, lbef, 4) .or.str_comp ( &
               befehl, 'osel', 2, lbef, 4) .or.str_comp (befehl, 'odes',&
               2, lbef, 4) ) then                                       
!                                                                       
                  lselect = str_comp (befehl, 'msel', 2, lbef, 4)       &
                  .or.str_comp (befehl, 'osel', 2, lbef, 4)             
!                                                                       
                  CALL mole_select (zeile, lp, 0, SYM_MAXSCAT,         &
                             sym_latom, sym_sel_atom, lold, lselect)
!                                                                       
!     ----Select the origin of the symmetry operation  'origin'         
!                                                                       
               ELSEIF (str_comp (befehl, 'origin', 1, lbef, 6) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0.and. (ianz.eq.3.or.ianz.eq.4) ) then 
                     IF (ianz.eq.4) then 
                        sym_orig_mol = str_comp (cpara (4) , 'mol', 1,  &
                        lpara (4) , 3)                                  
                        ianz = 3 
                     ENDIF 
!                                                                       
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        DO i = 1, 3 
                        sym_orig (i) = werte (i) 
                        ENDDO 
                        l_need_setup = .true. 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
!                                                                       
!     ----Set the power of the symmetry operation  'power'              
!                                                                       
               ELSEIF (str_comp (befehl, 'power', 1, lbef, 6) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (ianz.eq.1.or.ianz.eq.2) then 
                        IF (ianz.eq.2) then 
                           IF (str_comp (cpara (2) , 'mult', 1, lpara ( &
                           2) , 4) ) then                               
                              sym_power_mult = .true. 
                           ELSEIF (str_comp (cpara (2) , 'sing', 1,     &
                           lpara (2) , 4) ) then                        
                              sym_power_mult = .false. 
                           ELSE 
                              ier_num = - 6 
                              ier_typ = ER_COMM 
                           ENDIF 
                        ENDIF 
                        ianz = 1 
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxw)                                           
                        IF (ier_num.eq.0) then 
                           sym_power = nint (werte (1) ) 
                           l_need_setup = .true. 
                        ENDIF 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ENDIF 
!                                                                       
!     ----run symmetry 'run'                                            
!                                                                       
               ELSEIF (str_comp (befehl, 'run ', 1, lbef, 4) ) then 
                  IF (l_need_setup) then 
                     CALL symm_setup 
                  ENDIF 
                  IF (sym_sel_atom) then 
!                                                                       
!-----      --------Apply symmetry operation to atoms                   
!                                                                       
                     IF (sym_power_mult) then 
                        CALL symm_op_mult 
                     ELSE 
                        CALL symm_op_single 
                     ENDIF 
                  ELSE 
                     IF (sym_sel_mode.eq.SYM_RUN_MOLECULE) then 
!                                                                       
!-----      ----------Apply symmetry operation to molecules             
!                                                                       
                        IF (sym_power_mult) then 
                           CALL symm_mole_mult 
                        ELSE 
                           CALL symm_mole_single 
                        ENDIF 
                     ELSEIF (sym_sel_mode.eq.SYM_RUN_DOMAIN) then 
!                                                                       
!-----      ----------Apply symmetry operation to molecules             
!                                                                       
                        IF (sym_power_mult) then 
                           CALL symm_domain_mult 
                        ELSE 
                           CALL symm_domain_single 
                        ENDIF 
                     ENDIF 
                  ENDIF 
                  CALL update_cr_dim 
!                                                                       
!     ----Select which atoms are copied to their image 'sele'           
!                                                                       
               ELSEIF (str_comp (befehl, 'sele', 2, lbef, 4)            &
               .or.str_comp (befehl, 'dese', 2, lbef, 4) ) then         
!                                                                       
write(*,*) ' SYMM ',SYM_MAXSCAT
                  CALL atom_select (zeile, lp, 0, SYM_MAXSCAT, sym_latom, &
                  sym_sel_atom, lold,     &
                  str_comp (befehl, 'sele', 2, lbef, 4) )               
!                                                                       
!     ----show current parameters 'show'                                
!                                                                       
               ELSEIF (str_comp (befehl, 'show', 2, lbef, 4) ) then 
                  IF (l_need_setup) then 
                     CALL symm_setup 
                  ENDIF 
                  CALL symm_show 
!                                                                       
!     ----Select translational part of the symmetry operation 'trans'   
!                                                                       
               ELSEIF (str_comp (befehl, 'trans', 2, lbef, 5) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0.and.ianz.eq.3) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        DO i = 1, 3 
                        sym_trans (i) = werte (i) 
                        ENDDO 
                        l_need_setup = .true. 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
!                                                                       
!     ----Select the type of symmetry operation                  'type' 
!                                                                       
               ELSEIF (str_comp (befehl, 'type', 2, lbef, 4) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (ianz.eq.1) then 
                        IF (str_comp (cpara (1) , 'proper', 1, lpara (1)&
                        , 6) ) then                                     
                           sym_type = .true. 
                           l_need_setup = .true. 
                        ELSEIF (str_comp (cpara (1) , 'improper', 1,    &
                        lpara (1) , 8) ) then                           
                           sym_type = .false. 
                           l_need_setup = .true. 
                        ELSE 
                           ier_num = - 6 
                           ier_typ = ER_COMM 
                        ENDIF 
                     ENDIF 
                  ENDIF 
!                                                                       
!     ----Select the direct space direction of the symmetry axis 'uvw'  
!                                                                       
               ELSEIF (str_comp (befehl, 'uvw ', 1, lbef, 4) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0.and.ianz.eq.3) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        DO i = 1, 3 
                        sym_uvw (i) = werte (i) 
                        ENDDO 
                        CALL trans (sym_uvw, cr_gten, sym_hkl, 3) 
                        l_need_setup = .true. 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
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
 9999 CONTINUE 
!                                                                       
      END SUBROUTINE symm                           
!*****7*****************************************************************
      SUBROUTINE symm_show 
!-                                                                      
!     Shows current symm settings                                       
!+                                                                      
      USE config_mod 
      USE crystal_mod 
      USE molecule_mod 
      USE symm_mod 
      IMPLICIT none 
!                                                                       
      include'prompt.inc' 
       
!                                                                       
      CHARACTER(9) at_lis (maxscat+1), at_name 
      INTEGER mol_lis (maxscat+1)
      INTEGER i, j, k 
!                                                                       
      WRITE (output_io, 3000) sym_uvw 
      WRITE (output_io, 3010) sym_hkl 
      IF (sym_orig_mol) then 
         WRITE (output_io, 3020) sym_orig, ' rel.to molecule' 
      ELSE 
         WRITE (output_io, 3020) sym_orig, ' rel. to crystal' 
      ENDIF 
      WRITE (output_io, 3030) sym_angle 
      WRITE (output_io, 3040) sym_trans 
      WRITE (output_io, 3045) sym_or_tr 
      WRITE (output_io, 3050) ( (sym_mat (i, j), j = 1, 4), i = 1, 3) 
      WRITE (output_io, 3060) ( (sym_rmat (i, j), j = 1, 3), i = 1, 3) 
      WRITE (output_io, 3070) sym_power 
!                                                                       
      IF (sym_power_mult) then 
         WRITE (output_io, 3080) 'Multiple copy of original' 
      ELSE 
         WRITE (output_io, 3080) 'Single copy of original' 
      ENDIF 
!                                                                       
      IF (sym_type) then 
         WRITE (output_io, 3090) 'Proper rotation' 
      ELSE 
         WRITE (output_io, 3090) 'Improper rotation' 
      ENDIF 
!                                                                       
      IF (sym_mode) then 
         WRITE (output_io, 3100) 'Copy atom/molecule to new position' 
      ELSE 
         WRITE (output_io, 3100) 'Move atom/molecule to new position' 
      ENDIF 
!                                                                       
      IF (sym_new.and..not.sym_sel_atom) then 
         WRITE (output_io, 3110) 'Create new molecule type' 
      ELSE 
         WRITE (output_io, 3110) 'Keep molecule type' 
      ENDIF 
!                                                                       
!------ Working with atoms ...                                          
!                                                                       
      IF (sym_sel_atom) then 
!                                                                       
         j = 0 
         DO i = 0, cr_nscat 
         IF (sym_latom (i) ) then 
            j = j + 1 
            at_lis (j) = at_name (i) 
         ENDIF 
         ENDDO 
         WRITE (output_io, 3210) (at_lis (i), i = 1, j) 
!                                                                       
         IF (sym_incl.eq.'all ') then 
            WRITE (output_io, 3220) 
         ELSEIF (sym_incl.eq.'env ') then 
            WRITE (output_io, 3225) 
         ELSE 
            WRITE (output_io, 3230) sym_start, sym_end 
         ENDIF 
!                                                                       
!------ Working with molecules                                          
!                                                                       
      ELSE 
!                                                                       
         IF (sym_orig_mol) then 
            WRITE (output_io, 3250) 'Molecule' 
         ELSE 
            WRITE (output_io, 3250) 'Crystal' 
         ENDIF 
!                                                                       
         j = 0 
         DO i = 0, mole_num_type 
         IF (sym_latom (i) ) then 
            j = j + 1 
            mol_lis (j) = i 
         ENDIF 
         ENDDO 
         WRITE (output_io, 3300) (mol_lis (k), k = 1, j) 
!                                                                       
         IF (sym_end.eq. - 1) then 
            WRITE (output_io, 3310) 
         ELSE 
            WRITE (output_io, 3320) sym_start, sym_end 
         ENDIF 
         IF (sym_sel_mode.eq.SYM_RUN_DOMAIN) then 
            IF (sym_dom_mode_atom) then 
               WRITE (output_io, 4100) 
            ELSE 
               WRITE (output_io, 4150) 
            ENDIF 
            IF (sym_dom_mode_shape) then 
               WRITE (output_io, 4200) 
            ELSE 
               WRITE (output_io, 4250) 
            ENDIF 
         ENDIF 
      ENDIF 
!                                                                       
 3000 FORMAT    ( ' Generalized Symmetry Operation'/                    &
     &                   '   Axis in direct space      : ',3(2x,f9.4))  
 3010 FORMAT    ( '   Axis in reciprocal space  : ',3(2x,f9.4)) 
 3020 FORMAT    ( '   Origin of symmetry element: ',3(2x,f9.4),a) 
 3030 FORMAT    ( '   Rotation angle            : ',  2x,f9.4 ) 
 3040 FORMAT    ( '   Translation part(true)    : ',3(2x,f9.4)) 
 3045 FORMAT    ( '   Translation part(origin)  : ',3(2x,f9.4)/) 
 3050 FORMAT    ( '   Real space matrix         : ',4(2x,f9.4)/         &
     &                2( '                             : ',4(2x,f9.4)/))
 3060 FORMAT    ( '   Reciprocal space matrix   : ',3(2x,f9.4)/         &
     &                2( '                             : ',3(2x,f9.4)/))
 3070 FORMAT    ( '   Power of symmetry element : ', (2x,i4  )) 
 3080 FORMAT    ( '   Mode of power level       : ',2x,a) 
 3090 FORMAT    ( '   Type of symmetry element  : ',2x,a) 
 3100 FORMAT    ( '   Mode of symmetry operation: ',2x,a) 
 3110 FORMAT    ( '   Molecule status           : ',2x,a) 
 3120 FORMAT    (/'   Select status             : ',2x,a) 
 3210 FORMAT    ( '   Selected atom types       : ',2x,50(a9,1x)) 
 3220 FORMAT    ( '   Range of selected atoms   :   All atoms included') 
 3225 FORMAT( '   Range of selected atoms   :   Current environment') 
 3230 FORMAT    ( '   Range of selected atoms   : ',i9,' to ',i9) 
 3250 FORMAT    ( '   Given origin relative to  : ',2x,a) 
 3300 FORMAT    ( '   Selected molecule types   : ',2x,50(i4,1x)) 
 3310 FORMAT    ( '   Range of sel. molecules   : ',                    &
     &          '  All molecules included')                             
 3320 FORMAT    ( '   Range of sel. molecules   : ',i9,' to ',i9) 
 4100 FORMAT    ( '   Status of atoms in domain :   rotated') 
 4150 FORMAT    ( '   Status of shape of domain :   invariant') 
 4200 FORMAT    ( '   Status of shape of domain :   rotated') 
 4250 FORMAT    ( '   Status of atoms in domain :   invariant') 
!                                                                       
      END SUBROUTINE symm_show                      
!*****7*****************************************************************
      SUBROUTINE symm_setup 
!-                                                                      
!     Performs the generalized symmetry operation                       
!     See Sands, D.E. Vectors and Tensors in Crystallography Chapt. 4.7 
!+                                                                      
      USE config_mod 
      USE crystal_mod 
      USE symm_mod 
      IMPLICIT none 
!                                                                       
       
      include'errlist.inc' 
!                                                                       
      INTEGER i, j, k, l 
      REAL length 
!                                                                       
      REAL uij 
      REAL ctheta, stheta 
      REAL sym_d (3), sym_r (3) 
      REAL usym (4), ures (4) 
      REAL kron (3, 3) 
      REAL a (3, 3) 
      REAL b (3, 3) 
!                                                                       
      REAL cosd, sind 
      REAL skalpro 
!                                                                       
!                                                                       
      DATA kron / 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 / 
!                                                                       
!     initialize matrix and angle                                       
!                                                                       
      DO i = 1, 4 
      DO j = 1, 4 
      sym_mat (i, j) = 0.0 
      sym_rmat (i, j) = 0.0 
      ENDDO 
      ENDDO 
      sym_mat (4, 4) = 1.0 
      sym_rmat (4, 4) = 0.0 
      IF (sym_power_mult) then 
         ctheta = cosd (sym_angle) 
         stheta = sind (sym_angle) 
      ELSE 
         ctheta = cosd (sym_angle * sym_power) 
         stheta = sind (sym_angle * sym_power) 
      ENDIF 
!                                                                       
!     Create vectors of unit length in direct and reciprocal space      
!                                                                       
      length = sqrt (skalpro (sym_uvw, sym_uvw, cr_gten) ) 
      IF (length.eq.0.0) then 
         ier_num = - 32 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
      DO j = 1, 3 
      sym_d (j) = sym_uvw (j) / length 
      ENDDO 
      length = sqrt (skalpro (sym_hkl, sym_hkl, cr_rten) ) 
      IF (length.eq.0.0) then 
         ier_num = - 32 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
      DO j = 1, 3 
      sym_r (j) = sym_hkl (j) / length 
      ENDDO 
!     write (output_io,2001) sym_d,sym_r                                
!2001      format('direct axis     :',3(2x,f10.6)/                      
!    &                  'reciprocal axis :',3(2x,f10.6))                
!                                                                       
!     calculate symmetry operation                                      
!                                                                       
      DO i = 1, 3 
      DO j = 1, 3 
      sym_mat (i, j) = 0.0 
      DO k = 1, 3 
      DO l = 1, 3 
      sym_mat (i, j) = sym_mat (i, j) + cr_rten (i, k) * cr_eps (k, l,  &
      j) * sym_d (l)                                                    
      ENDDO 
      ENDDO 
      uij = sym_d (i) * sym_r (j) 
      sym_mat (i, j) = sym_mat (i, j) * stheta + uij + (kron (i, j)     &
      - uij) * ctheta                                                   
      ENDDO 
      IF (sym_power_mult) then 
         sym_mat (i, 4) = sym_trans (i) 
      ELSE 
         sym_mat (i, 4) = sym_trans (i) * sym_power 
      ENDIF 
      ENDDO 
!                                                                       
!     In case of improper rotation, multiply by -1                      
!                                                                       
      IF (.not.sym_type.and. (                                          &
      sym_power_mult.or..not.sym_power_mult.and.mod (sym_power, 2)      &
      .ne.0) ) then                                                     
         DO i = 1, 3 
         DO j = 1, 3 
         sym_mat (i, j) = - sym_mat (i, j) 
         ENDDO 
         ENDDO 
      ENDIF 
!     write (output_io,2000) ((sym_mat(i,j),j=1,3),i=1,3)               
!2000      format(3(3(2x,f10.6)/))                                      
!                                                                       
!                                                                       
!                                                                       
!     Transform symmetry operation into reciprocal space                
!                                                                       
      DO i = 1, 3 
      DO j = 1, 3 
      a (i, j) = sym_mat (i, j) 
      ENDDO 
      ENDDO 
!                                                                       
!     do transformation q = gSg*                                        
!                                                                       
      CALL matmulx (b, a, cr_rten) 
      CALL matmulx (a, cr_gten, b) 
      DO i = 1, 3 
      DO j = 1, 3 
      sym_rmat (i, j) = a (i, j) 
      ENDDO 
      ENDDO 
!                                                                       
!     ----Calculate translational component due to origin               
!                                                                       
!                                                                       
!     ----Apply rotational part of symmetry operation to -1*origin      
!                                                                       
      DO j = 1, 3 
      usym (j) = - sym_orig (j) 
      ENDDO 
!                                                                       
!-----      ----Apply symmetry operation                                
!                                                                       
      usym (4) = 0.0 
      CALL trans (usym, sym_mat, ures, 4) 
!                                                                       
!     ----Add origin to result to obtain total translational part       
!                                                                       
      DO j = 1, 3 
      sym_or_tr (j) = ures (j) + sym_orig (j) 
      ENDDO 
!     write (output_io,2000) ((sym_rmat(i,j),j=1,3),i=1,3)              
!2000      format(3(3(2x,f10.6)/))                                      
!                                                                       
      END SUBROUTINE symm_setup                     
!*****7*****************************************************************
      SUBROUTINE symm_op_mult 
!-                                                                      
!     Performs the actual symmetry operation, multiple copy version     
!+                                                                      
      USE config_mod 
      USE crystal_mod 
      USE atom_env_mod 
      USE modify_mod
      USE symm_mod 
      IMPLICIT none 
!                                                                       
       
      include'errlist.inc' 
!                                                                       
      CHARACTER(4) name 
!                                                                       
      INTEGER i, j, k, l 
      INTEGER i_start, i_end 
!                                                                       
      REAL usym (4), ures (4) 
      REAL werte (5) 
!                                                                       
      DATA usym / 0.0, 0.0, 0.0, 1.0 / 
      DATA werte / 0.0, 0.0, 0.0, 0.0, 0.0 / 
!                                                                       
!     Set the appropriate starting end ending number for the atoms      
!                                                                       
      i_start = sym_start 
      i_end = sym_end 
      IF (sym_incl.eq.'all ') then 
         i_end = cr_natoms 
      ELSEIF (sym_incl.eq.'env ') then 
         i_end = atom_env (0) 
      ENDIF 
!                                                                       
!     Apply symmetry operation to all atoms within selected range       
!                                                                       
      DO l = i_start, i_end 
      i = l 
      IF (sym_incl.eq.'env ') i = atom_env (l) 
!                                                                       
!     --Select atom if:                                                 
!       type has been selected                                          
!                                                                       
      IF (sym_latom (cr_iscat (i) ) ) then 
!                                                                       
!     ----Subtract origin                                               
!                                                                       
         DO j = 1, 3 
         usym (j) = cr_pos (j, i) - sym_orig (j) 
         ENDDO 
!                                                                       
!-----      ----Apply symmetry operation                                
!                                                                       
         usym (4) = 1.0 
         DO k = 1, sym_power 
         CALL trans (usym, sym_mat, ures, 4) 
!                                                                       
!     ----Add origin                                                    
!                                                                       
         DO j = 1, 3 
         werte (j + 1) = ures (j) + sym_orig (j) 
         ENDDO 
!                                                                       
!     ----Insert copy of atom or replace original atom by its image     
!                                                                       
         IF (sym_mode) then 
            name = cr_at_lis (cr_iscat (i) ) 
            werte (5) = cr_dw (cr_iscat (i) ) 
            CALL do_ins_atom (name, werte, 5) 
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
         ELSE 
            DO j = 1, 3 
            cr_pos (j, i) = werte (j + 1) 
            ENDDO 
         ENDIF 
         DO j = 1, 3 
         usym (j) = ures (j) 
         ENDDO 
         ENDDO 
      ENDIF 
      ENDDO 
!                                                                       
      END SUBROUTINE symm_op_mult                   
!*****7*****************************************************************
      SUBROUTINE symm_op_single 
!-                                                                      
!     Performs the actual symmetry operation, single result version     
!+                                                                      
      USE config_mod 
      USE crystal_mod 
      USE atom_env_mod 
      USE modify_mod
      USE symm_mod 
      IMPLICIT none 
!                                                                       
       
      include'errlist.inc' 
!                                                                       
      CHARACTER(4) name 
      INTEGER i, j, l 
      INTEGER i_start, i_end 
      REAL usym (4), ures (4) 
      REAL werte (5) 
!                                                                       
      DATA usym / 0.0, 0.0, 0.0, 1.0 / 
      DATA werte / 0.0, 0.0, 0.0, 0.0, 0.0 / 
!                                                                       
!     Set the appropriate starting end ending number for the atoms      
!                                                                       
      i_start = sym_start 
      i_end = sym_end 
      IF (sym_incl.eq.'all ') then 
         i_end = cr_natoms 
      ELSEIF (sym_incl.eq.'env ') then 
         i_end = atom_env (0) 
      ENDIF 
!                                                                       
!     Apply symmetry operation to all atoms within selected range       
!                                                                       
      DO l = i_start, i_end 
      i = l 
      IF (sym_incl.eq.'env ') i = atom_env (l) 
!                                                                       
!     --Select atom if:                                                 
!       type has been selected                                          
!                                                                       
      IF (sym_latom (cr_iscat (i) ) ) then 
!                                                                       
!     ----Subtract origin                                               
!                                                                       
         DO j = 1, 3 
         usym (j) = cr_pos (j, i) - sym_orig (j) 
         ENDDO 
!                                                                       
!-----      ----Apply symmetry operation                                
!                                                                       
         usym (4) = 1.0 
         CALL trans (usym, sym_mat, ures, 4) 
!                                                                       
!     ----Add origin                                                    
!                                                                       
         DO j = 1, 3 
         werte (j + 1) = ures (j) + sym_orig (j) 
         ENDDO 
!                                                                       
!     ----Insert copy of atom or replace original atom by its image     
!                                                                       
         IF (sym_mode) then 
            name = cr_at_lis (cr_iscat (i) ) 
            werte (5) = cr_dw (cr_iscat (i) ) 
            CALL do_ins_atom (name, werte, 5) 
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
         ELSE 
            DO j = 1, 3 
            cr_pos (j, i) = werte (j + 1) 
            ENDDO 
         ENDIF 
      ENDIF 
      ENDDO 
!                                                                       
      END SUBROUTINE symm_op_single                 
!*****7*****************************************************************
      SUBROUTINE symm_mole_mult 
!-                                                                      
!     Performs the actual symmetry operation, multiple copy version     
!     Operates on molecules                                             
!+                                                                      
      USE allocate_appl_mod
      USE config_mod 
      USE crystal_mod 
      USE modify_mod
      USE molecule_mod 
      USE spcgr_apply, ONLY: mole_insert_current
      USE symm_mod 
      IMPLICIT none 
!                                                                       
       
      include'errlist.inc' 
!                                                                       
      CHARACTER(4) name 
!                                                                       
      INTEGER i, j, k, l, ii 
      INTEGER i_start, i_end 
      INTEGER imole, imole_s, imole_t 
      INTEGER  :: n_gene   ! Number of molecule generators
      INTEGER  :: n_symm   ! Number of molecule symmetry operators
      INTEGER  :: n_mole   ! Number of molecules
      INTEGER  :: n_type   ! Number of molecule types
      INTEGER  :: n_atom   ! Number of atoms in molecules
!                                                                       
      REAL usym (4), ures (4), use_orig (3) 
      REAL werte (5) 
!                                                                       
      DATA usym / 0.0, 0.0, 0.0, 1.0 / 
      DATA werte / 0.0, 0.0, 0.0, 0.0, 0.0 / 
!                                                                       
!     Set the appropriate starting end ending number for the atoms      
!                                                                       
      i_start = sym_start 
      i_end = sym_end 
      IF (sym_end.eq. - 1) i_end = mole_num_mole 
      imole_s = mole_num_mole 
!                                                                       
      IF (sym_new) then 
         n_gene = MAX( 1, MOLE_MAX_GENE)
         n_symm = MAX( 1, MOLE_MAX_SYMM)
         n_mole =         MOLE_MAX_MOLE
         n_type =         MOLE_MAX_TYPE
         n_atom =         MOLE_MAX_ATOM
         IF (mole_num_type+sym_power > MOLE_MAX_TYPE ) THEN
            n_type = MAX(mole_num_type + sym_power+5,MOLE_MAX_TYPE)
            call alloc_molecule(n_gene, n_symm, n_mole, n_type, n_atom)
         ENDIF
         imole_t = mole_num_type 
         IF (mole_num_type.lt.MOLE_MAX_TYPE) then 
            mole_num_type = mole_num_type+sym_power 
         ELSE 
            ier_num = - 66 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
      ENDIF 
!                                                                       
!     Apply symmetry operation to all molecules within selected range   
!                                                                       
      DO l = i_start, i_end 
!                                                                       
!------ - Determine origin for symmetry operation                       
!                                                                       
      IF (sym_orig_mol) then 
         DO i = 1, 3 
         ii = mole_cont (mole_off (l) + 1) 
         use_orig (i) = sym_orig (i) + cr_pos (i, ii) 
         ENDDO 
      ELSE 
         DO i = 1, 3 
         use_orig (i) = sym_orig (i) 
         ENDDO 
      ENDIF 
!                                                                       
!     --Select atom if:                                                 
!       type has been selected                                          
!                                                                       
      IF (sym_latom (mole_type (l) ) ) then 
!                                                                       
!     ----Loop over all atoms in the molecule                           
!                                                                       
         DO ii = 1, mole_len (l) 
         i = mole_cont (mole_off (l) + ii) 
!                                                                       
!     ----- Subtract origin                                             
!                                                                       
         DO j = 1, 3 
         usym (j) = cr_pos (j, i) - use_orig (j) 
         ENDDO 
!                                                                       
!-----      ----- Apply symmetry operation                              
!                                                                       
         usym (4) = 1.0 
         DO k = 1, sym_power 
         CALL trans (usym, sym_mat, ures, 4) 
!                                                                       
!     ------- Add origin                                                
!                                                                       
         DO j = 1, 3 
         werte (j + 1) = ures (j) + use_orig (j) 
         ENDDO 
!                                                                       
!     ------- Insert copy of atom or replace original atom by its image 
!                                                                       
         IF (sym_mode) then 
            name = cr_at_lis (cr_iscat (i) ) 
            werte (5) = cr_dw (cr_iscat (i) ) 
            CALL do_ins_atom (name, werte, 5) 
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
!                                                                       
!     --------- Insert atom into proper new molecule                    
!                                                                       
            imole = imole_s + (l - i_start) * sym_power + (k - 1)       &
            + 1                                                         
            CALL mole_insert_current (cr_natoms, imole) 
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
         ELSE 
            DO j = 1, 3 
            cr_pos (j, i) = werte (j + 1) 
            ENDDO 
         ENDIF 
         DO j = 1, 3 
         usym (j) = ures (j) 
         ENDDO 
         ENDDO 
         ENDDO 
!                                                                       
!-----      ----Copy all properties                                     
!                                                                       
         IF (sym_mode) then 
            CALL copy_mole_char (mole_num_mole, l) 
         ENDIF 
!                                                                       
!----- ---- Set new molecule type if requested                          
!                                                                       
         IF (sym_new) then 
            DO k = 1, sym_power 
            IF (sym_mode) then 
               mole_type (mole_num_mole-sym_power + k) = imole_t + k 
            ELSE 
               mole_type (l) = imole_t + 1 
            ENDIF 
            ENDDO 
         ENDIF 
      ENDIF 
      ENDDO 
!                                                                       
      END SUBROUTINE symm_mole_mult                 
!*****7*****************************************************************
      SUBROUTINE symm_mole_single 
!-                                                                      
!     Performs the actual symmetry operation, single result version     
!     Operates on molecules                                             
!+                                                                      
      USE allocate_appl_mod
      USE config_mod 
      USE crystal_mod 
      USE modify_mod
      USE molecule_mod 
      USE spcgr_apply, ONLY: mole_insert_current
      USE symm_mod 
      IMPLICIT none 
!                                                                       
       
      include'errlist.inc' 
!                                                                       
      CHARACTER(4) name 
      INTEGER i, j, ii, l 
      INTEGER i_start, i_end 
      INTEGER imole, imole_s, imole_t 
      INTEGER  :: n_gene   ! Number of molecule generators
      INTEGER  :: n_symm   ! Number of molecule symmetry operators
      INTEGER  :: n_mole   ! Number of molecules
      INTEGER  :: n_type   ! Number of molecule types
      INTEGER  :: n_atom   ! Number of atoms in molecules
      REAL usym (4), ures (4) 
      REAL werte (5), use_orig (3) 
!                                                                       
      DATA usym / 0.0, 0.0, 0.0, 1.0 / 
      DATA werte / 0.0, 0.0, 0.0, 0.0, 0.0 / 
!                                                                       
!     Set the appropriate starting end ending number for the molecules  
!                                                                       
      i_start = sym_start 
      i_end = sym_end 
      IF (sym_end.eq. - 1) i_end = mole_num_mole 
      imole_s = mole_num_mole 
!                                                                       
      IF (sym_new) then 
         n_gene = MAX( 1, MOLE_MAX_GENE)
         n_symm = MAX( 1, MOLE_MAX_SYMM)
         n_mole =         MOLE_MAX_MOLE
         n_type =         MOLE_MAX_TYPE
         n_atom =         MOLE_MAX_ATOM
         IF (mole_num_type+sym_power > MOLE_MAX_TYPE ) THEN
            n_type = MAX(mole_num_type + sym_power+5,MOLE_MAX_TYPE)
            call alloc_molecule(n_gene, n_symm, n_mole, n_type, n_atom)
         ENDIF
         IF (mole_num_type.lt.MOLE_MAX_TYPE) then 
            mole_num_type = mole_num_type+1 
         ELSE 
            ier_num = - 66 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
         imole_t = mole_num_type 
      ENDIF 
!                                                                       
!     Apply symmetry operation to all molecules within selected range   
!                                                                       
      DO l = i_start, i_end 
!                                                                       
!------ - Determine origin for symmetry operation                       
!                                                                       
      IF (sym_orig_mol) then 
         DO i = 1, 3 
         ii = mole_cont (mole_off (l) + 1) 
         use_orig (i) = sym_orig (i) + cr_pos (i, ii) 
         ENDDO 
      ELSE 
         DO i = 1, 3 
         use_orig (i) = sym_orig (i) 
         ENDDO 
      ENDIF 
!                                                                       
!     --Select atom if:                                                 
!       type has been selected                                          
!                                                                       
      IF (sym_latom (mole_type (l) ) ) then 
!                                                                       
!     ----Loop over all atoms in the molecule                           
!                                                                       
         DO ii = 1, mole_len (l) 
         i = mole_cont (mole_off (l) + ii) 
!                                                                       
!     ----- Subtract origin                                             
!                                                                       
         DO j = 1, 3 
         usym (j) = cr_pos (j, i) - use_orig (j) 
         ENDDO 
!                                                                       
!-----      ----- Apply symmetry operation                              
!                                                                       
         usym (4) = 1.0 
         CALL trans (usym, sym_mat, ures, 4) 
!                                                                       
!     ----- Add origin                                                  
!                                                                       
         DO j = 1, 3 
         werte (j + 1) = ures (j) + use_orig (j) 
         ENDDO 
!                                                                       
!     ----- Insert copy of atom or replace original atom by its image   
!                                                                       
         IF (sym_mode) then 
            name = cr_at_lis (cr_iscat (i) ) 
            werte (5) = cr_dw (cr_iscat (i) ) 
            CALL do_ins_atom (name, werte, 5) 
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
!                                                                       
!     ------- Insert atom into proper new molecule                      
!                                                                       
            imole = imole_s + l - i_start + 1 
            CALL mole_insert_current (cr_natoms, imole) 
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
         ELSE 
            DO j = 1, 3 
            cr_pos (j, i) = werte (j + 1) 
            ENDDO 
         ENDIF 
         ENDDO 
!                                                                       
!-----      ----Copy all properties                                     
!                                                                       
         IF (sym_mode) then 
            CALL copy_mole_char (mole_num_mole, l) 
         ENDIF 
!                                                                       
!----- ---- Set new molecule type if requested                          
!                                                                       
         IF (sym_new) then 
            IF (sym_mode) then 
               mole_type (mole_num_mole) = imole_t 
            ELSE 
               mole_type (l) = imole_t 
            ENDIF 
         ENDIF 
      ENDIF 
      ENDDO 
!                                                                       
      END SUBROUTINE symm_mole_single               
!*****7*****************************************************************
      SUBROUTINE symm_domain_mult 
!-                                                                      
!     Performs the actual symmetry operation, multiple copy version     
!     Operates on molecules                                             
!+                                                                      
      USE allocate_appl_mod
      USE config_mod 
      USE crystal_mod 
      USE modify_mod
      USE molecule_mod 
      USE spcgr_apply, ONLY: mole_insert_current
      USE symm_mod 
      IMPLICIT none 
!                                                                       
       
      include'errlist.inc' 
!                                                                       
      CHARACTER(4) name 
!                                                                       
      INTEGER i, j, k, l, ii 
      INTEGER i_start, i_end 
      INTEGER imole, imole_s, imole_t 
      INTEGER  :: n_gene   ! Number of molecule generators
      INTEGER  :: n_symm   ! Number of molecule symmetry operators
      INTEGER  :: n_mole   ! Number of molecules
      INTEGER  :: n_type   ! Number of molecule types
      INTEGER  :: n_atom   ! Number of atoms in molecules
!                                                                       
      REAL mat_atom (4, 4) 
      REAL mat_dime (4, 4) 
      REAL new_atom (4, 4) 
      REAL new_dime (4, 4) 
      REAL elements (3, 8) 
      REAL usym (4), ures (4), use_orig (3) 
      REAL werte (5) 
!                                                                       
      DATA usym / 0.0, 0.0, 0.0, 1.0 / 
      DATA werte / 0.0, 0.0, 0.0, 0.0, 0.0 / 
!                                                                       
!     Set the appropriate starting end ending number for the atoms      
!                                                                       
      i_start = sym_start 
      i_end = sym_end 
      IF (sym_end.eq. - 1) i_end = mole_num_mole 
      imole_s = mole_num_mole 
!                                                                       
      IF (sym_new) then 
         n_gene = MAX( 1, MOLE_MAX_GENE)
         n_symm = MAX( 1, MOLE_MAX_SYMM)
         n_mole =         MOLE_MAX_MOLE
         n_type =         MOLE_MAX_TYPE
         n_atom =         MOLE_MAX_ATOM
         IF (mole_num_type+sym_power > MOLE_MAX_TYPE ) THEN
            n_type = MAX(mole_num_type + sym_power+5,MOLE_MAX_TYPE)
            call alloc_molecule(n_gene, n_symm, n_mole, n_type, n_atom)
         ENDIF
         imole_t = mole_num_type 
         IF (mole_num_type.lt.MOLE_MAX_TYPE) then 
            mole_num_type = mole_num_type+sym_power 
         ELSE 
            ier_num = - 66 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
      ENDIF 
!                                                                       
!     Apply symmetry operation to all molecules within selected range   
!                                                                       
      DO l = i_start, i_end 
!                                                                       
!------ - Determine origin for symmetry operation                       
!                                                                       
      IF (sym_orig_mol) then 
         DO i = 1, 3 
         ii = mole_cont (mole_off (l) + 1) 
         use_orig (i) = sym_orig (i) + cr_pos (i, ii) 
         ENDDO 
      ELSE 
         DO i = 1, 3 
         use_orig (i) = sym_orig (i) 
         ENDDO 
      ENDIF 
!                                                                       
!     --Select atom if:                                                 
!       type has been selected                                          
!                                                                       
      IF (sym_latom (mole_type (l) ) ) then 
!                                                                       
!     ----Create the matrices from the psueodoatom positions            
!                                                                       
         DO i = 1, 8 
         DO j = 1, 3 
         elements (j, i) = cr_pos (j, mole_cont (mole_off (l) + i) ) 
         ENDDO 
         ENDDO 
         DO i = 1, 3 
         DO j = 1, 3 
         mat_atom (i, j) = cr_pos (j, mole_cont (mole_off (l) + i + 1) )&
         - cr_pos (j, mole_cont (mole_off (l) + 1) )                    
         mat_dime (i, j) = cr_pos (j, mole_cont (mole_off (l) + i + 5) )&
         - cr_pos (j, mole_cont (mole_off (l) + 5) )                    
         ENDDO 
         mat_atom (4, i) = 0.0 
         mat_dime (4, i) = 0.0 
         mat_atom (i, 4) = cr_pos (i, mole_cont (mole_off (l) + 1) ) 
         mat_dime (i, 4) = cr_pos (i, mole_cont (mole_off (l) + 5) ) 
         ENDDO 
         mat_atom (4, 4) = 1.0 
         mat_dime (4, 4) = 1.0 
!                                                                       
!     ----Loop over Power of Operation                                  
!                                                                       
         DO k = 1, sym_power 
         CALL matmul4 (new_atom, sym_mat, mat_atom) 
         CALL matmul4 (new_dime, sym_mat, mat_dime) 
         DO i = 1, 3 
         DO j = 1, 3 
         mat_atom (i, j) = new_atom (i, j) 
         mat_dime (i, j) = new_dime (i, j) 
         ENDDO 
         ENDDO 
!                                                                       
!     ----Loop over all atoms in the molecule                           
!                                                                       
         DO ii = 1, 5, 4 
         i = mole_cont (mole_off (l) + ii) 
!                                                                       
!     ----- Subtract origin                                             
!                                                                       
         DO j = 1, 3 
         usym (j) = cr_pos (j, i) - use_orig (j) 
         usym (j) = elements (j, ii) - use_orig (j) 
         ENDDO 
!                                                                       
!-----      ----- Apply symmetry operation                              
!                                                                       
         usym (4) = 1.0 
         CALL trans (usym, sym_mat, ures, 4) 
!                                                                       
!     ------- Add origin                                                
!                                                                       
         DO j = 1, 3 
         werte (j + 1) = ures (j) + use_orig (j) 
         elements (j, ii) = ures (j) + use_orig (j) 
         ENDDO 
         ENDDO 
         IF (sym_dom_mode_atom) then 
            DO i = 1, 3 
            DO j = 1, 3 
            elements (j, 1 + i) = new_atom (i, j) + elements (j, 1) 
            ENDDO 
            ENDDO 
         ELSE 
            DO j = 1, 3 
            elements (j, 1) = cr_pos (j, mole_cont (mole_off (l)        &
            + 1) )                                                      
            ENDDO 
         ENDIF 
         IF (sym_dom_mode_shape) then 
            DO i = 1, 3 
            DO j = 1, 3 
            elements (j, 5 + i) = new_dime (i, j) + elements (j, 5) 
            ENDDO 
            ENDDO 
         ELSE 
            DO j = 1, 3 
            elements (j, 5) = cr_pos (j, mole_cont (mole_off (l)        &
            + 5) )                                                      
            ENDDO 
         ENDIF 
         DO ii = 1, 8 
         i = mole_cont (mole_off (l) + ii) 
         DO j = 1, 3 
         werte (j + 1) = elements (j, ii) 
         ENDDO 
!                                                                       
!     ------- Insert copy of atom or replace original atom by its image 
!                                                                       
         IF (sym_mode) then 
            name = cr_at_lis (cr_iscat (i) ) 
            werte (5) = cr_dw (cr_iscat (i) ) 
            CALL do_ins_atom (name, werte, 5) 
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
!                                                                       
!     --------- Insert atom into proper new molecule                    
!                                                                       
            imole = imole_s + (l - i_start) * sym_power + (k - 1)       &
            + 1                                                         
            CALL mole_insert_current (cr_natoms, imole) 
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
         ELSE 
            DO j = 1, 3 
            cr_pos (j, i) = werte (j + 1) 
            ENDDO 
         ENDIF 
         DO j = 1, 3 
         usym (j) = ures (j) 
         ENDDO 
         ENDDO 
         ENDDO 
!                                                                       
!-----      ----Copy all properties                                     
!                                                                       
         IF (sym_mode) then 
            CALL copy_mole_char (mole_num_mole, l) 
         ENDIF 
!                                                                       
!----- ---- Set new molecule type if requested                          
!                                                                       
         IF (sym_new) then 
            DO k = 1, sym_power 
            IF (sym_mode) then 
               mole_type (mole_num_mole-sym_power + k) = imole_t + k 
            ELSE 
               mole_type (l) = imole_t + 1 
            ENDIF 
            ENDDO 
         ENDIF 
      ENDIF 
      ENDDO 
!                                                                       
      END SUBROUTINE symm_domain_mult               
!*****7*****************************************************************
      SUBROUTINE symm_domain_single 
!-                                                                      
!     Performs the actual symmetry operation, single result version     
!     Operates on microdomain represenatation                           
!+                                                                      
      USE allocate_appl_mod
      USE config_mod 
      USE crystal_mod 
      USE modify_mod
      USE molecule_mod 
      USE spcgr_apply, ONLY: mole_insert_current
      USE symm_mod 
      IMPLICIT none 
!                                                                       
       
      include'errlist.inc' 
!                                                                       
      CHARACTER(4) name 
      INTEGER i, j, ii, l 
      INTEGER i_start, i_end 
      INTEGER a_start, a_end 
      INTEGER imole, imole_s, imole_t 
      INTEGER  :: n_gene   ! Number of molecule generators
      INTEGER  :: n_symm   ! Number of molecule symmetry operators
      INTEGER  :: n_mole   ! Number of molecules
      INTEGER  :: n_type   ! Number of molecule types
      INTEGER  :: n_atom   ! Number of atoms in molecules
      REAL mat_atom (4, 4) 
      REAL mat_dime (4, 4) 
      REAL new_atom (4, 4) 
      REAL new_dime (4, 4) 
      REAL elements (3, 8) 
      REAL usym (4), ures (4) 
      REAL werte (5), use_orig (3) 
!                                                                       
      DATA usym / 0.0, 0.0, 0.0, 1.0 / 
      DATA werte / 0.0, 0.0, 0.0, 0.0, 0.0 / 
!                                                                       
!     Set the appropriate starting end ending number for the molecules  
!                                                                       
      i_start = sym_start 
      i_end = sym_end 
      IF (sym_end.eq. - 1) i_end = mole_num_mole 
      imole_s = mole_num_mole 
!                                                                       
      IF (sym_new) then 
         n_gene = MAX( 1, MOLE_MAX_GENE)
         n_symm = MAX( 1, MOLE_MAX_SYMM)
         n_mole =         MOLE_MAX_MOLE
         n_type =         MOLE_MAX_TYPE
         n_atom =         MOLE_MAX_ATOM
         IF (mole_num_type+sym_power > MOLE_MAX_TYPE ) THEN
            n_type = MAX(mole_num_type + sym_power+5,MOLE_MAX_TYPE)
            call alloc_molecule(n_gene, n_symm, n_mole, n_type, n_atom)
         ENDIF
         IF (mole_num_type.lt.MOLE_MAX_TYPE) then 
            mole_num_type = mole_num_type+1 
         ELSE 
            ier_num = - 66 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
         imole_t = mole_num_type 
      ENDIF 
!                                                                       
!     Apply symmetry operation to all molecules within selected range   
!                                                                       
      DO l = i_start, i_end 
!                                                                       
!------ - Determine origin for symmetry operation                       
!                                                                       
      IF (sym_orig_mol) then 
         DO i = 1, 3 
         ii = mole_cont (mole_off (l) + 1) 
         use_orig (i) = sym_orig (i) + cr_pos (i, ii) 
         ENDDO 
      ELSE 
         DO i = 1, 3 
         use_orig (i) = sym_orig (i) 
         ENDDO 
      ENDIF 
!                                                                       
!     --Select atom if:                                                 
!       type has been selected                                          
!                                                                       
      IF (sym_latom (mole_type (l) ) ) then 
!                                                                       
!     ----Create the matrices from the psueodoatom positions            
!                                                                       
         DO i = 1, 8 
         DO j = 1, 3 
         elements (j, i) = cr_pos (j, mole_cont (mole_off (l) + i) ) 
         ENDDO 
         ENDDO 
         DO i = 1, 3 
         DO j = 1, 3 
         mat_atom (i, j) = cr_pos (j, mole_cont (mole_off (l) + i + 1) )&
         - cr_pos (j, mole_cont (mole_off (l) + 1) )                    
         mat_dime (i, j) = cr_pos (j, mole_cont (mole_off (l) + i + 5) )&
         - cr_pos (j, mole_cont (mole_off (l) + 5) )                    
         ENDDO 
         mat_atom (4, i) = 0.0 
         mat_dime (4, i) = 0.0 
         mat_atom (i, 4) = cr_pos (i, mole_cont (mole_off (l) + 1) ) 
         mat_dime (i, 4) = cr_pos (i, mole_cont (mole_off (l) + 5) ) 
         ENDDO 
         mat_atom (4, 4) = 1.0 
         mat_dime (4, 4) = 1.0 
         CALL matmul4 (new_atom, sym_mat, mat_atom) 
         CALL matmul4 (new_dime, sym_mat, mat_dime) 
!                                                                       
!     ----Loop over the two origins of the microdomain                  
!                                                                       
         DO ii = 1, 5, 4 
         i = mole_cont (mole_off (l) + ii) 
!                                                                       
!     ----- Subtract origin                                             
!                                                                       
         DO j = 1, 3 
         usym (j) = cr_pos (j, i) - use_orig (j) 
         ENDDO 
!                                                                       
!-----      ----- Apply symmetry operation                              
!                                                                       
         usym (4) = 1.0 
         CALL trans (usym, sym_mat, ures, 4) 
!                                                                       
!     ----- Add origin                                                  
!                                                                       
         DO j = 1, 3 
         werte (j + 1) = ures (j) + use_orig (j) 
         elements (j, ii) = ures (j) + use_orig (j) 
         ENDDO 
         ENDDO 
                                                                        
         IF (sym_dom_mode_atom) then 
            DO i = 1, 3 
            DO j = 1, 3 
            elements (j, 1 + i) = new_atom (i, j) + elements (j, 1) 
            ENDDO 
            ENDDO 
         ELSE 
            DO j = 1, 3 
            elements (j, 1) = cr_pos (j, mole_cont (mole_off (l)        &
            + 1) )                                                      
            ENDDO 
         ENDIF 
         IF (sym_dom_mode_shape) then 
            DO i = 1, 3 
            DO j = 1, 3 
            elements (j, 5 + i) = new_dime (i, j) + elements (j, 5) 
            ENDDO 
            ENDDO 
         ELSE 
            DO j = 1, 3 
            elements (j, 5) = cr_pos (j, mole_cont (mole_off (l)        &
            + 5) )                                                      
            ENDDO 
         ENDIF 
!                                                                       
!     ----- Insert copy of atom or replace original atom by its image   
!                                                                       
         DO ii = 1, 8 
         i = mole_cont (mole_off (l) + ii) 
         DO j = 1, 3 
         werte (j + 1) = elements (j, ii) 
         ENDDO 
         IF (sym_mode) then 
            name = cr_at_lis (cr_iscat (i) ) 
            werte (5) = cr_dw (cr_iscat (i) ) 
            CALL do_ins_atom (name, werte, 5) 
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
!                                                                       
!     ------- Insert atom into proper new molecule                      
!                                                                       
            imole = imole_s + l - i_start + 1 
            CALL mole_insert_current (cr_natoms, imole) 
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
         ELSE 
            DO j = 1, 3 
            cr_pos (j, i) = werte (j + 1) 
            ENDDO 
         ENDIF 
         ENDDO 
!                                                                       
!-----      ----Copy all properties                                     
!                                                                       
         IF (sym_mode) then 
            CALL copy_mole_char (mole_num_mole, l) 
         ENDIF 
!                                                                       
!----- ---- Set new molecule type if requested                          
!                                                                       
         IF (sym_new) then 
            IF (sym_mode) then 
               mole_type (mole_num_mole) = imole_t 
            ELSE 
               mole_type (l) = imole_t 
            ENDIF 
         ENDIF 
      ENDIF 
      ENDDO 
!                                                                       
      END SUBROUTINE symm_domain_single             
!*****7*****************************************************************
      SUBROUTINE symm_ca_mult (uvw, lspace) 
!-                                                                      
!     Performs the actual symmetry operation, multiple copy version     
!     Only the input vector uvw is used in direct or reciprocal space   
!+                                                                      
      USE config_mod 
      USE symm_mod 
      IMPLICIT none 
!                                                                       
      include'prompt.inc' 
      include'errlist.inc' 
      include'param.inc' 
!                                                                       
      INTEGER j, k 
      LOGICAL lspace 
!                                                                       
      REAL uvw (3) 
      REAL usym (4), ures (4) 
      REAL werte (5) 
!                                                                       
      DATA usym / 0.0, 0.0, 0.0, 1.0 / 
      DATA werte / 0.0, 0.0, 0.0, 0.0, 0.0 / 
!                                                                       
!     real space part                                                   
!                                                                       
      IF (lspace) then 
!                                                                       
!     ----Subtract origin, if in real space                             
!                                                                       
         DO j = 1, 3 
         usym (j) = uvw (j) - sym_orig (j) 
         ENDDO 
!                                                                       
!-----      --Apply symmetry operation                                  
!                                                                       
         usym (4) = 1.0 
         DO k = 1, sym_power 
         CALL trans (usym, sym_mat, ures, 4) 
!                                                                       
!     ----Add origin and store result                                   
!                                                                       
         DO j = 1, 3 
         res_para ( (k - 1) * 3 + j) = ures (j) + sym_orig (j) 
         ENDDO 
         WRITE (output_io, 3000) (res_para ( (k - 1) * 3 + j), j = 1, 3) 
!                                                                       
!     ----Replace current vector by its image                           
!                                                                       
         DO j = 1, 3 
         usym (j) = ures (j) 
         ENDDO 
         ENDDO 
      ELSE 
!                                                                       
!     ----Reciprocal space part                                         
!                                                                       
         DO j = 1, 3 
         usym (j) = uvw (j) 
         ENDDO 
!                                                                       
!-----      --Apply symmetry operation                                  
!                                                                       
         usym (4) = 0.0 
         DO k = 1, sym_power 
         CALL trans (usym, sym_rmat, ures, 4) 
!                                                                       
!     ----Store result                                                  
!                                                                       
         DO j = 1, 3 
         res_para ( (k - 1) * 3 + j) = ures (j) 
         ENDDO 
         WRITE (output_io, 3000) (res_para ( (k - 1) * 3 + j), j = 1, 3) 
!                                                                       
!     ----Replace current vector by its image                           
!                                                                       
         DO j = 1, 3 
         usym (j) = ures (j) 
         ENDDO 
         ENDDO 
      ENDIF 
!                                                                       
      res_para (0) = sym_power * 3 
!                                                                       
 3000 FORMAT    (' Result    : ',3(2x,f9.4)) 
      END SUBROUTINE symm_ca_mult                   
!*****7*****************************************************************
      SUBROUTINE symm_ca_single (uvw, lspace, loutput) 
!-                                                                      
!     Performs the actual symmetry operation, multiple copy version     
!     Only the input vector uvw is used in direct or reciprocal space   
!+                                                                      
      USE config_mod 
      USE symm_mod 
      IMPLICIT none 
!                                                                       
      include'prompt.inc' 
      include'errlist.inc' 
      include'param.inc' 
!                                                                       
      REAL   ,DIMENSION(1:3), INTENT(IN) :: uvw
      LOGICAL,                INTENT(IN) :: lspace 
      LOGICAL,                INTENT(IN) :: loutput 
!                                                                       
      INTEGER j 
!
      REAL usym (4), ures (4) 
      REAL werte (5) 
!                                                                       
      DATA usym / 0.0, 0.0, 0.0, 1.0 / 
      DATA werte / 0.0, 0.0, 0.0, 0.0, 0.0 / 
!                                                                       
!     real space part                                                   
!                                                                       
      IF (lspace) then 
!                                                                       
!     ----Subtract origin, if in real space                             
!                                                                       
         DO j = 1, 3 
         usym (j) = uvw (j) - sym_orig (j) 
         ENDDO 
!                                                                       
!-----      --Apply symmetry operation                                  
!                                                                       
         usym (4) = 1.0 
         CALL trans (usym, sym_mat, ures, 4) 
!                                                                       
!     ----Add origin and store result                                   
!                                                                       
         DO j = 1, 3 
         res_para (j) = ures (j) + sym_orig (j) 
         ENDDO 
!                                                                       
!     ----Replace current vector by its image                           
!                                                                       
         DO j = 1, 3 
         usym (j) = ures (j) 
         ENDDO 
      ELSE 
!                                                                       
!     ----Subtract origin, if in real space                             
!                                                                       
         DO j = 1, 3 
         usym (j) = uvw (j) 
         ENDDO 
!                                                                       
!-----      --Apply symmetry operation                                  
!                                                                       
         usym (4) = 0.0 
         CALL trans (usym, sym_rmat, ures, 4) 
!                                                                       
!     ----Add origin and store result                                   
!                                                                       
         DO j = 1, 3 
         res_para (j) = ures (j) 
         ENDDO 
!                                                                       
!     ----Replace current vector by its image                           
!                                                                       
         DO j = 1, 3 
         usym (j) = ures (j) 
         ENDDO 
      ENDIF 
!                                                                       
      IF( loutput) THEN
         res_para (0) = 3 
         WRITE (output_io, 3000) (res_para (j), j = 1, 3) 
      ENDIF
!                                                                       
 3000 FORMAT    (' Result    : ',3(2x,f9.4)) 
      END SUBROUTINE symm_ca_single                 
