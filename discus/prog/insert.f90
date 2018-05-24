MODULE insert_menu
!
CONTAINS
!+                                                                      
!     Input of objects (and molecules)                                  
!                                                                       
!*****7*****************************************************************
      SUBROUTINE insert (itype) 
!-                                                                      
!     Main menu for inserting objects                                   
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE insert_mod 
      USE molecule_mod 
      USE discus_show_menu
      USE update_cr_dim_mod
!
      USE ber_params_mod
      USE build_name_mod
      USE calc_expr_mod
      USE doact_mod 
      USE do_eval_mod
      USE do_wait_mod
      USE errlist_mod 
      USE get_params_mod
      USE learn_mod 
      USE class_macro_internal
      USE prompt_mod 
      USE string_convert_mod
      USE sup_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER, INTENT(IN) :: itype 
!                                                                       
      INTEGER, PARAMETER  :: maxw =  30
!                                                                       
      CHARACTER(LEN=5)                  :: befehl 
      CHARACTER(LEN=8), DIMENSION(-1:1) :: ctype
      CHARACTER(LEN=LEN(prompt))        :: orig_prompt
      CHARACTER(LEN=1024)               :: line
      CHARACTER(LEN=1024)               :: zeile
      CHARACTER(LEN=1024), DIMENSION(1:MAXW) :: cpara (MAXW) 
      INTEGER,             DIMENSION(1:MAXW) :: lpara
      INTEGER                                :: lp, length, lbef 
      INTEGER                                :: indxg, ianz, i, j
      LOGICAL                                :: lend
      REAL   ,             DIMENSION(1:MAXW) :: werte (MAXW) 
!                                                                       
      INTEGER :: len_str 
      LOGICAL :: str_comp 
!                                                                       
      DATA ctype / 'domain  ', 'molecule', 'object  ' / 
!                                                                       
      lend = .false. 
      CALL no_error 
      orig_prompt = prompt
!                                                                       
      prompt = prompt (1:len_str (prompt) ) //'/insert'//'/'//ctype (itype)
      DO while (.not.lend) 
      CALL get_cmd (line, length, befehl, lbef, zeile, lp, prompt) 
      IF (ier_num.eq.0) then 
         IF (line /= ' '      .and. line(1:1) /= '#' .and. &
             line /= char(13) .and. line(1:1) /= '!'        ) THEN
!                                                                       
!     ----search for "="                                                
!                                                                       
indxg = index (line, '=') 
IF (indxg.ne.0.AND..NOT. (str_comp (befehl, 'echo', 2, lbef, 4) ) &
              .AND..NOT. (str_comp (befehl, 'syst', 2, lbef, 4) )    &
              .AND..NOT. (str_comp (befehl, 'help', 2, lbef, 4) .OR. &
                          str_comp (befehl, '?   ', 2, lbef, 4) )    &
              .AND. INDEX(line,'==') == 0                            ) THEN
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
!     ----Select the atomic displacementparameter of the                
!           object  'density'                                           
!                                                                       
               ELSEIF (str_comp (befehl, 'adp', 1, lbef, 3) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0.and.ianz.eq.1) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        ins_adp = werte (1) 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
!                                                                       
!     ----define the character of the object                            
!                                                                       
               ELSEIF (str_comp (befehl, 'character', 2, lbef, 9) )     &
               then                                                     
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (ianz.eq.1.or.ianz.eq.2) then 
                        IF (itype.eq. - 1) then 
                           IF (str_comp (cpara (1) , 'cube', 3, lpara ( &
                           1) , 4) ) then                               
                              ins_character = MOLE_DOM_CUBE 
                           ELSEIF (str_comp (cpara (1) , 'cylinder', 3, &
                           lpara (1) , 8) ) then                        
                              ins_character = MOLE_DOM_CYLINDER 
                           ELSEIF (str_comp (cpara (1) , 'sphere', 3,   &
                           lpara (1) , 6) ) then                        
                              ins_character = MOLE_DOM_SPHERE 
                           ELSEIF (str_comp (cpara (1) , 'fuzzy', 3,    &
                           lpara (1) , 5) ) then                        
                              ins_character = MOLE_DOM_FUZZY 
                           ELSEIF (str_comp (cpara (1) , 'atom', 3,     &
                           lpara (1) , 4) ) then                        
                              ins_character = MOLE_ATOM 
                           ELSE 
                              ier_num = - 6 
                              ier_typ = ER_COMM 
                           ENDIF 
                           ins_obj_atom = 'POSI' 
                        ELSEIF (itype.eq.0) then 
                           IF (str_comp (cpara (1) , 'atom', 3, lpara ( &
                           1) , 4) ) then                               
                              ins_character = MOLE_ATOM 
                           ELSE 
                              ier_num = - 6 
                              ier_typ = ER_COMM 
                           ENDIF 
                        ELSEIF (itype.eq.1) then 
                           IF (str_comp (cpara (1) , 'cube', 3, lpara ( &
                           1) , 4) ) then                               
                              ins_character = MOLE_CUBE 
                           ELSEIF (str_comp (cpara (1) , 'cylinder', 3, &
                           lpara (1) , 8) ) then                        
                              ins_character = MOLE_CYLINDER 
                           ELSEIF (str_comp (cpara (1) , 'sphere', 3,   &
                           lpara (1) , 6) ) then                        
                              ins_character = MOLE_SPHERE 
                           ELSEIF (str_comp (cpara (1) , 'atom', 3,     &
                           lpara (1) , 4) ) then                        
                              ins_character = MOLE_ATOM 
                           ELSE 
                              ier_num = - 6 
                              ier_typ = ER_COMM 
                           ENDIF 
                           IF (ianz.eq.1) then 
                              ins_obj_atom = 'E1- ' 
                           ELSE 
                              ins_obj_atom = cpara (2) (1:4) 
                              CALL do_cap (ins_obj_atom) 
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
!     ----help 'help','?'                                               
!                                                                       
      ELSEIF (str_comp (befehl, 'help', 2, lbef, 4) .or.str_comp (befehl&
     &, '?   ', 1, lbef, 4) ) then                                      
                  IF (str_comp (zeile, 'errors', 2, lp, 6) ) then 
                     lp = lp + 7 
                     CALL do_hel ('discus '//zeile, lp) 
                  ELSE 
                     lp = lp + 12 
                     CALL do_hel ('discus insert '//zeile, lp) 
                  ENDIF 
!                                                                       
!     ----Enter the Center of atoms within a domain 'cent'              
!                                                                       
               ELSEIF (str_comp (befehl, 'cent', 3, lbef, 4) ) then 
                  IF (itype.lt.0) then 
                     CALL get_params (zeile, ianz, cpara, lpara, maxw,  &
                     lp)                                                
                     IF (ier_num.eq.0.and.ianz.eq.3) then 
!                                                                       
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxw)                                           
                        IF (ier_num.eq.0) then 
                           DO j = 1, 3 
                           ins_cent (j) = werte (j) 
                           ENDDO 
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
!     ----Select the scattering density of the object  'density'        
!                                                                       
               ELSEIF (str_comp (befehl, 'density', 1, lbef, 7) ) then 
                  IF (itype.gt.0) then 
                     CALL get_params (zeile, ianz, cpara, lpara, maxw,  &
                     lp)                                                
                     IF (ier_num.eq.0.and.ianz.eq.1) then 
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxw)                                           
                        IF (ier_num.eq.0) then 
                           ins_density = werte (1) 
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
!     ----Select the isotropic B-value  of the object  'density'        
!                                                                       
               ELSEIF (str_comp (befehl, 'biso', 1, lbef, 4) ) then 
                  IF (itype.gt.0) then 
                     CALL get_params (zeile, ianz, cpara, lpara, maxw,  &
                     lp)                                                
                     IF (ier_num.eq.0.and.ianz.eq.1) then 
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxw)                                           
                        IF (ier_num.eq.0) then 
                           ins_biso = werte (1) 
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
!     ----Select the domain file name 'file'                            
!                                                                       
               ELSEIF (str_comp (befehl, 'file', 2, lbef, 4) ) then 
                  IF (itype.lt.0) then 
                     CALL get_params (zeile, ianz, cpara, lpara, maxw,  &
                     lp)                                                
                     IF (ier_num.eq.0.and.ianz.ge.1) then 
                        CALL do_build_name (ianz, cpara, lpara, werte,  &
                        maxw, 1)                                        
                        IF (ier_num.eq.0) then 
                           ins_file = cpara (1) (1:lpara(1))
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
!     ----Select the distance between microdomain atoms and host        
!           atoms 'fuzzy'                                               
!                                                                       
               ELSEIF (str_comp (befehl, 'fuzzy', 2, lbef, 5) ) then 
                  IF (itype.lt.0) then 
                     CALL get_params (zeile, ianz, cpara, lpara, maxw,  &
                     lp)                                                
                     IF (ier_num.eq.0.and.ianz.eq.1) then 
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxw)                                           
                        IF (ier_num.eq.0) then 
                           ins_fuzzy = werte (1) 
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
!     ----Select the origin of the object  'origin'                     
!                                                                       
               ELSEIF (str_comp (befehl, 'origin', 1, lbef, 6) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0.and.ianz.eq.3) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        DO i = 1, 3 
                        ins_origin (i) = werte (i) 
                        ENDDO 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
!                                                                       
!     ----Enter the type of the  object 'type'                          
!                                                                       
               ELSEIF (str_comp (befehl, 'type', 3, lbef, 4) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0.and.ianz.eq.1) then 
!                                                                       
                     IF (str_comp (cpara (1) , 'new', 3, lpara (1) , 3) &
                     ) then                                             
                        ins_type = INS_NEWTYPE 
                     ELSE 
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxw)                                           
                        IF (ier_num.eq.0) then 
                           IF (werte (1) .eq.INS_NEWTYPE) then 
                              ins_type = INS_NEWTYPE 
                           ELSEIF (werte (1) .gt.0) then 
                              ins_type = int (werte (1) ) 
                           ENDIF 
                        ELSE 
                           ier_num = - 6 
                           ier_typ = ER_COMM 
                        ENDIF 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
!                                                                       
!     ----Enter the X-axis of an object or domain 'xaxis'               
!                                                                       
               ELSEIF (str_comp (befehl, 'xaxis', 3, lbef, 5) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0.and.ianz.eq.3) then 
!                                                                       
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        DO j = 1, 3 
                        ins_xaxis (j) = werte (j) 
                        ENDDO 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
!                                                                       
!     ----Enter the Y-axis of an object or domain 'yaxis'               
!                                                                       
               ELSEIF (str_comp (befehl, 'yaxis', 3, lbef, 5) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0.and.ianz.eq.3) then 
!                                                                       
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        DO j = 1, 3 
                        ins_yaxis (j) = werte (j) 
                        ENDDO 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
!                                                                       
!     ----Enter the Z-axis of an object or domain 'zaxis'               
!                                                                       
               ELSEIF (str_comp (befehl, 'zaxis', 3, lbef, 5) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0.and.ianz.eq.3) then 
!                                                                       
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        DO j = 1, 3 
                        ins_zaxis (j) = werte (j) 
                        ENDDO 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
!                                                                       
!     ----Enter the X-dimension of a domain 'xdim'                      
!                                                                       
               ELSEIF (str_comp (befehl, 'xdim', 3, lbef, 4) ) then 
                  IF (itype.lt.0) then 
                     CALL get_params (zeile, ianz, cpara, lpara, maxw,  &
                     lp)                                                
                     IF (ier_num.eq.0.and.ianz.eq.3) then 
!                                                                       
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxw)                                           
                        IF (ier_num.eq.0) then 
                           DO j = 1, 3 
                           ins_xdim (j) = werte (j) 
                           ENDDO 
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
!     ----Enter the Y-dimension of a domain 'ydim'                      
!                                                                       
               ELSEIF (str_comp (befehl, 'ydim', 3, lbef, 4) ) then 
                  IF (itype.lt.0) then 
                     CALL get_params (zeile, ianz, cpara, lpara, maxw,  &
                     lp)                                                
                     IF (ier_num.eq.0.and.ianz.eq.3) then 
!                                                                       
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxw)                                           
                        IF (ier_num.eq.0) then 
                           DO j = 1, 3 
                           ins_ydim (j) = werte (j) 
                           ENDDO 
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
!     ----Enter the Z-dimension of a domain 'zdim'                      
!                                                                       
               ELSEIF (str_comp (befehl, 'zdim', 3, lbef, 4) ) then 
                  IF (itype.lt.0) then 
                     CALL get_params (zeile, ianz, cpara, lpara, maxw,  &
                     lp)                                                
                     IF (ier_num.eq.0.and.ianz.eq.3) then 
!                                                                       
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxw)                                           
                        IF (ier_num.eq.0) then 
                           DO j = 1, 3 
                           ins_zdim (j) = werte (j) 
                           ENDDO 
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
!     ----run shear 'run'                                               
!                                                                       
               ELSEIF (str_comp (befehl, 'run', 1, lbef, 3) ) then 
                  IF (ins_character.eq.MOLE_ATOM) then 
!                                                                       
!-----      --------Insert a molecule into the structure                
!                                                                       
                     CALL insert_molecules 
                  ELSEIF (ins_character.lt.0) then 
!                                                                       
!-----      --------Insert an object into the structure                 
!                                                                       
                     CALL insert_domain 
                  ELSEIF (ins_character.gt.0) then 
!                                                                       
!-----      --------Insert an object into the structure                 
!                                                                       
                     CALL insert_object 
                  ENDIF 
                  CALL update_cr_dim 
!                                                                       
!     ----show current parameters 'show'                                
!                                                                       
               ELSEIF (str_comp (befehl, 'show', 2, lbef, 4) ) then 
                  CALL insert_show () 
!                                                                       
!     ----GENERAL MENU COMMANDS                                         
!                                                                       
!                                                                       
!     ----exit 'exit'                                                   
!                                                                       
               ELSEIF (str_comp (befehl, 'exit', 2, lbef, 4) ) then 
                  lend = .true. 
!     ----list asymmetric unit 'asym'                                   
!                                                                       
               ELSEIF (str_comp (befehl, 'asym', 2, lbef, 4) ) then 
                  CALL show_asym 
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
      IF (ier_num.ne.0) THEN 
         CALL errlist 
         IF (ier_sta.ne.ER_S_LIVE) THEN 
            IF (lmakro .OR. lmakro_error) THEN  ! Error within macro or termination errror
               IF(sprompt /= prompt ) THEN
                  ier_num = -10
                  ier_typ = ER_COMM
                  ier_msg(1) = ' Error occured in insert menu'
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
      ENDDO 
!
      prompt = orig_prompt
!                                                                       
      END SUBROUTINE insert                         
!*****7*****************************************************************
      SUBROUTINE insert_show ()
!-                                                                      
!     Shows current shear settings                                      
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE insert_mod 
      USE molecule_mod 
!
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
      CHARACTER(8) c_char (0:4) 
      CHARACTER(8) ctype ( - 1:1) 
!                                                                       
      DATA c_char / 'Atom    ', 'Cube    ', 'Cylinder', 'Sphere  ', 'Fuz&
     &zy   ' /                                                          
      DATA ctype / 'Domain  ', 'Molecule', 'Object  ' / 
!                                                                       
      IF (ins_character.eq.MOLE_ATOM) then 
         WRITE (output_io, 3000) 'Inserting a molecule' 
      ELSEIF (ins_character.gt.0) then 
         WRITE (output_io, 3000) 'Inserting an object' 
         WRITE (output_io, 3010) c_char ( (ins_character) ) 
         WRITE (output_io, 3012) ins_obj_atom 
         WRITE (output_io, 3013) ins_adp 
         IF (ins_type.eq.INS_NEWTYPE) then 
            WRITE (output_io, 3015) 
         ELSE 
            WRITE (output_io, 3016) ins_type 
         ENDIF 
         WRITE (output_io, 3019) ins_density 
         WRITE (output_io, 3021) ins_biso    
         WRITE (output_io, 3020) ins_origin 
         WRITE (output_io, 3030) 'x', ins_xaxis 
         WRITE (output_io, 3030) 'y', ins_yaxis 
         WRITE (output_io, 3030) 'z', ins_zaxis 
      ELSEIF (ins_character.lt.0) then 
         WRITE (output_io, 3000) 'Inserting a domain' 
         WRITE (output_io, 3010) c_char (abs (ins_character) ) 
         WRITE (output_io, 3110) ins_file 
         WRITE (output_io, 3120) ins_fuzzy 
         WRITE (output_io, 3020) ins_origin 
         WRITE (output_io, 3030) 'x', ins_xaxis 
         WRITE (output_io, 3030) 'y', ins_yaxis 
         WRITE (output_io, 3030) 'z', ins_zaxis 
         WRITE (output_io, 3140) ins_cent 
         WRITE (output_io, 3150) 'x', ins_xdim 
         WRITE (output_io, 3150) 'y', ins_ydim 
         WRITE (output_io, 3150) 'z', ins_zdim 
      ENDIF 
!                                                                       
 3000 FORMAT    ( ' Generalized Insertion Operation'/                   &
     &                   '   Type of operation  : ',a/)                 
 3010 FORMAT    ( '   character          : ',a) 
 3012 FORMAT    ( '   atom type          : ',a4) 
 3013 FORMAT    ( '   adp                : ',2x,f9.4) 
 3015 FORMAT    ( '   type number        : New type') 
 3016 FORMAT    ( '   type number        : ',i6) 
 3019 FORMAT    ( '   Scattering density : ',2x,f9.4) 
 3021 FORMAT    ( '   Isotropic B-value  : ',2x,f9.4) 
 3020 FORMAT    ( '   origin             : ',3(2x,f9.4)/) 
 3030 FORMAT    ( '   ',a1,'-axis             : ',3(2x,f9.4)) 
 3110 FORMAT    ( '   file               : ',a) 
 3120 FORMAT    ( '   fuzzy distance     : ',2x,f9.4) 
 3140 FORMAT    (/'   center of shape    : ',3(2x,f9.4)/) 
 3150 FORMAT    ( '   ',a1,'-dimension        : ',3(2x,f9.4)) 
!                                                                       
      END SUBROUTINE insert_show                    
!*****7*****************************************************************
      SUBROUTINE insert_molecules 
!-                                                                      
!     Places an individual molecule into the postion defined in the menu
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE molecule_mod 
!
      USE errlist_mod 
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      WRITE (output_io, * ) ' Not yet implemented' 
!                                                                       
      END SUBROUTINE insert_molecules               
!*****7*****************************************************************
      SUBROUTINE insert_object 
!-                                                                      
!     Places an individual object                                       
!+                                                                      
      USE discus_config_mod 
      USE discus_allocate_appl_mod
      USE crystal_mod 
      USE insert_mod 
      USE modify_mod
      USE molecule_mod 
      USE prop_para_mod 
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      CHARACTER(LEN=4), DIMENSION(4)   :: cname 
      CHARACTER(LEN=4), DIMENSION(0:3) :: oname
      INTEGER          :: i, j, ii, k, kk 
      INTEGER          :: new_nmax
      INTEGER          :: new_nscat
      INTEGER          :: n_gene
      INTEGER          :: n_symm
      INTEGER          :: n_mole
      INTEGER          :: n_type
      INTEGER          :: n_atom
      LOGICAL          :: need_alloc = .false.
!                                                                       
      DATA cname / '    ', 'XAXI', 'YAXI', 'ZAXI' / 
      DATA oname / 'atom', 'cube', 'cyli', 'sphe' / 
!                                                                       
!     Make sure there is enough space left in the crystal               
!     1) Check maximum number of molecules                              
!     2) Check maximum number of molecule types                         
!     3) Check maximum number of atoms                                  
!                                                                       
      need_alloc = .false.
      n_gene = MAX( 1, MOLE_MAX_GENE)
      n_symm = MAX( 1, MOLE_MAX_SYMM)
      n_mole =         MOLE_MAX_MOLE
      n_type =         MOLE_MAX_TYPE
      n_atom =         MOLE_MAX_ATOM
      IF (mole_num_mole == MOLE_MAX_MOLE ) THEN
         n_mole = MOLE_MAX_MOLE + 50
         need_alloc = .true.
      ENDIF
      IF (mole_num_atom == MOLE_MAX_ATOM ) THEN
         n_atom = MOLE_MAX_ATOM + 50*8
         need_alloc = .true.
      ENDIF
      IF (ins_type.eq.INS_NEWTYPE.and.mole_num_type == MOLE_MAX_MOLE) THEN
         n_mole = MOLE_MAX_MOLE + 50
         need_alloc = .true.
      ENDIF
      IF ( need_alloc ) THEN
         call alloc_molecule(n_gene, n_symm, n_mole, n_type, n_atom)
      ENDIF
      IF (mole_num_mole.lt.MOLE_MAX_MOLE) then 
      IF (ins_type.eq.INS_NEWTYPE.and.mole_num_type.lt.MOLE_MAX_MOLE.or.&
     &ins_type.lt.MOLE_MAX_MOLE) then                                   
!
!     While developing, increment crystal if neede, but keep the check
!
            need_alloc = .false.
            new_nmax   = NMAX
            new_nscat  = MAXSCAT
            IF ( NMAX <= cr_natoms + 4 ) THEN 
               new_nmax  = MAX ( NMAX    + 1, INT (NMAX    * 1.25) )
               need_alloc = .true.
            ENDIF
            IF ( MAXSCAT <= cr_nscat+4 ) THEN 
               new_nscat = MAX ( MAXSCAT + 1, INT (MAXSCAT * 1.25) )
               need_alloc = .true.
            ENDIF
            IF( need_alloc ) THEN
              call alloc_crystal(new_nscat, new_nmax)
              IF ( ier_num /= 0) RETURN
            ENDIF
            IF (cr_natoms + 4.le.NMAX) then 
!                                                                       
!     ------ Increment the molecule number                              
!                                                                       
               i = mole_num_mole+1 
!                                                                       
!     ------ Set the object character                                   
!                                                                       
               mole_char (i) = ins_character 
!                                                                       
!     ------ Set the object type, if existing type find first           
!            molecule of this type                                      
!                                                                       
               IF (ins_type.eq.INS_NEWTYPE) then 
                  mole_num_type = mole_num_type+1 
                  ii = mole_off (mole_num_mole) + mole_len (            &
                  mole_num_mole) + 1                                    
                  mole_type (i) = mole_num_type 
                  mole_dens (i) = ins_density 
                  mole_biso (mole_type(i)) = ins_biso    
               ELSE 
                  mole_type (i) = ins_type 
                  kk = 1 
                  DO while (kk.le.mole_num_type.and.mole_type (kk)      &
                  .ne.ins_type)                                         
                  kk = kk + 1 
                  ENDDO 
                                                                        
                  IF (kk.gt.mole_num_type) then 
                     ier_num = - 64 
                     ier_typ = ER_APPL 
                     RETURN 
                  ENDIF 
                  CALL copy_mole_char (i, kk) 
               ENDIF 
!                                                                       
!     ------ Set the object origin                                      
!                                                                       
               DO j = 1, 3 
               cr_pos (j, cr_natoms + 1) = ins_origin (j) 
               ENDDO 
!                                                                       
!------ ------ Set the rotation components, add origin                  
!                                                                       
               DO j = 1, 3 
               cr_pos (j, cr_natoms + 2) = ins_xaxis (j) + ins_origin ( &
               j)                                                       
               cr_pos (j, cr_natoms + 3) = ins_yaxis (j) + ins_origin ( &
               j)                                                       
               cr_pos (j, cr_natoms + 4) = ins_zaxis (j) + ins_origin ( &
               j)                                                       
               ENDDO 
!                                                                       
               DO j = 1, 4 
               cr_prop (cr_natoms + j) = 0 
               cr_prop (cr_natoms + j) = ibset (cr_prop (cr_natoms + j),&
               PROP_NORMAL)                                             
               ENDDO 
               DO j=1,4
                  cr_surf(:,cr_natoms + j) = 0
               ENDDO
!                                                                       
!     ------ Set the atom types and names                               
!                                                                       
               IF (ins_type.eq.INS_NEWTYPE) then 
!                                                                       
!     -------- New object type                                          
!                                                                       
                  cname (1) = ins_obj_atom 
                  DO k = 1, 4 
                  j = 1 
                  IF (k.eq.1) then 
                     DO while (j.le.cr_nscat.and. (cname (k)            &
                     .ne.cr_at_lis (j) .or.ins_adp.ne.cr_dw (j) ) )     !! WORK OCC
                     j = j + 1 
                     ENDDO 
                  ELSE 
                     DO while (j.le.cr_nscat.and.cname (k)              &
                     .ne.cr_at_lis (j) )                                
                     j = j + 1 
                     ENDDO 
                  ENDIF 
                  IF (j.gt.cr_nscat) then 
                     cr_iscat (cr_natoms + k) = cr_nscat + 1 
                     IF (k.eq.1) then 
                        cr_dw (cr_nscat + 1) = ins_adp 
                     ELSE 
                        cr_dw (cr_nscat + 1) = 0.0 
                     ENDIF 
                     cr_occ(cr_nscat+1) = 1.00                 !! WORK OCC
                     cr_at_lis (cr_nscat + 1) = cname (k) 
                     cr_nscat = cr_nscat + 1 
                  ELSE 
                     cr_iscat (cr_natoms + k) = j 
                  ENDIF 
                  ENDDO 
               ELSE 
!                                                                       
!     -------- Old object type                                          
!                                                                       
                  DO j = 1, 4 
                  ii = mole_cont (mole_off (kk) + j) 
                  cr_iscat (cr_natoms + j) = cr_iscat (ii) 
                  ENDDO 
               ENDIF 
               ii = mole_off (mole_num_mole) + mole_len (mole_num_mole) 
               DO j = 1, 4 
               mole_cont (ii + j) = cr_natoms + j 
               ENDDO 
               mole_off (i) = ii 
               mole_len (i) = 4 
               cr_natoms = cr_natoms + 4 
               mole_num_mole = mole_num_mole+1 
            ELSE 
               ier_num = - 10 
               ier_typ = ER_APPL 
            ENDIF 
         ELSE 
            ier_num = - 66 
            ier_typ = ER_APPL 
         ENDIF 
      ELSE 
         ier_num = - 65 
         ier_typ = ER_APPL 
      ENDIF 
!                                                                       
      END SUBROUTINE insert_object                  
!*****7*****************************************************************
      SUBROUTINE insert_domain 
!-                                                                      
!     Places an individual domain                                       
!+                                                                      
      USE discus_config_mod 
      USE discus_allocate_appl_mod
      USE crystal_mod 
      USE insert_mod 
      USE modify_mod
      USE molecule_mod 
      USE prop_para_mod 
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      CHARACTER(LEN=4), DIMENSION(8)   :: cname
      CHARACTER(LEN=4), DIMENSION(0:3) :: oname
      INTEGER          :: i, j, ii, k, kk 
      INTEGER          :: new_nmax
      INTEGER          :: new_nscat
      INTEGER          :: n_gene
      INTEGER          :: n_symm
      INTEGER          :: n_mole
      INTEGER          :: n_type
      INTEGER          :: n_atom
      LOGICAL          :: need_alloc = .false.
!                                                                       
      DATA cname / '    ', 'XAXI', 'YAXI', 'ZAXI', &
                   'CENT', 'XDIM', 'YDIM', 'ZDIM' /
      DATA oname / 'atom', 'cube', 'cyli', 'sphe' / 
!                                                                       
!     Make sure there is enough space left in the crystal               
!     1) Check maximum number of molecules                              
!     2) Check maximum number of molecule types                         
!     3) Check maximum number of atoms                                  
!                                                                       
      need_alloc = .false.
      n_gene = MAX( 1, MOLE_MAX_GENE)
      n_symm = MAX( 1, MOLE_MAX_SYMM)
      n_mole =         MOLE_MAX_MOLE
      n_type =         MOLE_MAX_TYPE
      n_atom =         MOLE_MAX_ATOM
!
      IF (mole_num_mole == MOLE_MAX_MOLE ) THEN
         n_mole = MOLE_MAX_MOLE + 50
         need_alloc = .true.
      ENDIF
      IF (mole_num_atom >  MOLE_MAX_ATOM-8) THEN
         n_atom = MOLE_MAX_ATOM + 50*8
         need_alloc = .true.
      ENDIF
      IF (ins_type.eq.INS_NEWTYPE.and.mole_num_type == MOLE_MAX_TYPE) THEN
         n_mole = MOLE_MAX_MOLE + 50
         n_type = MOLE_MAX_TYPE + 50
         need_alloc = .true.
      ENDIF
      IF ( need_alloc ) THEN
         call alloc_molecule(n_gene, n_symm, n_mole, n_type, n_atom)
      ENDIF
      enough_mole: IF (mole_num_mole.lt.MOLE_MAX_MOLE) then 
         enough_type: IF (ins_type.eq.INS_NEWTYPE.and.mole_num_type.lt.MOLE_MAX_MOLE.or.&
             ins_type.lt.MOLE_MAX_MOLE) then
!
!     While developing, increment crystal if needed, but keep the check
!

            need_alloc = .false.
            new_nmax   = NMAX
            new_nscat  = MAXSCAT
            IF ( NMAX <= cr_natoms + 8 ) THEN 
               new_nmax  = MAX ( NMAX    + 8, INT (NMAX    * 1.25) )
               need_alloc = .true.
            ENDIF
            IF ( MAXSCAT <= cr_nscat+8 ) THEN 
               new_nscat = MAX ( MAXSCAT + 8, INT (MAXSCAT * 1.25) )
               need_alloc = .true.
            ENDIF
            IF( need_alloc ) THEN
              call alloc_crystal(new_nscat, new_nmax)
              IF ( ier_num /= 0) RETURN
            ENDIF
            enough_atom: IF (cr_natoms + 8.le.NMAX) then 
!                                                                       
!     ------ Increment the molecule number                              
!                                                                       
               i = mole_num_mole+1 
!                                                                       
!     ------ Set the domain character                                   
!                                                                       
               mole_char (i) = ins_character 
!                                                                       
!     ------ Set the domain type, if existing type find first           
!            molecule of this type                                      
!                                                                       
               IF (ins_type == INS_NEWTYPE) then 
                  mole_num_type  = mole_num_type+1 
                  ii = mole_off (mole_num_mole) + mole_len (mole_num_mole) + 1
                  mole_type (i)  = mole_num_type 
                  mole_dens (i)  = 1.0 
                  mole_biso (mole_type(i))  = 1.0 
                  mole_file (i)  = ins_file 
                  mole_fuzzy (i) = ins_fuzzy 
               ELSE 
                  mole_type (i) = ins_type 
                  kk = 1 
                  DO while (kk.le.mole_num_type.and.mole_type (kk).ne.ins_type)
                      kk = kk + 1 
                  ENDDO 
                                                                        
                  IF (kk.gt.mole_num_type) then 
                     ier_num = - 64 
                     ier_typ = ER_APPL 
                     RETURN 
                  ENDIF 
!     ------ Copy all old domain characters                             
                  CALL copy_mole_char (i, kk) 
               ENDIF 
!                                                                       
!     ------ Set the domain origin                                      
!                                                                       
               DO j = 1, 3 
                  cr_pos (j, cr_natoms + 1) = ins_origin (j) 
               ENDDO 
!                                                                       
!------ ------ Set the rotation components, add origin                  
!                                                                       
               DO j = 1, 3 
                  cr_pos (j, cr_natoms + 2) = ins_xaxis (j) + ins_origin(j)
                  cr_pos (j, cr_natoms + 3) = ins_yaxis (j) + ins_origin(j)
                  cr_pos (j, cr_natoms + 4) = ins_zaxis (j) + ins_origin(j)
               ENDDO 
!                                                                       
!     ------ Set the domain shape center                                
!                                                                       
               DO j = 1, 3 
                  cr_pos (j, cr_natoms + 5) = ins_cent (j) 
               ENDDO 
!                                                                       
!------ ------ Set the rotation components, add origin                  
!                                                                       
               DO j = 1, 3 
                  cr_pos (j, cr_natoms + 6) = ins_xdim (j) + ins_origin (j) 
                  cr_pos (j, cr_natoms + 7) = ins_ydim (j) + ins_origin (j) 
                  cr_pos (j, cr_natoms + 8) = ins_zdim (j) + ins_origin (j) 
               ENDDO 
               DO k = 1, 8 
                  cr_prop (cr_natoms + k) = 0 
                  cr_prop (cr_natoms + k) = ibset (cr_prop (cr_natoms + k),PROP_NORMAL)
               ENDDO 
               DO j=1,8
                  cr_surf(:,cr_natoms + j) = 0
               ENDDO
!                                                                       
!     ------ Set the atom types and names                               
!                                                                       
               set_atom: IF (ins_type.eq.INS_NEWTYPE) then 
!                                                                       
!     -------- New object type                                          
!                                                                       
                  cname (1) = ins_obj_atom 
                  DO k = 1, 8 
                     j = 1 
                     IF (k.eq.1) then 
                        DO while (j.le.cr_nscat.and.                  &
                                  (cname (k) .ne.cr_at_lis (j) .or.ins_adp.ne.cr_dw (j) ) )!! WORK OCC
                           j = j + 1 
                        ENDDO 
                     ELSE 
                        DO while (j.le.cr_nscat.and.cname (k) .ne.cr_at_lis (j) )
                           j = j + 1 
                        ENDDO 
                     ENDIF 
                     IF (j.gt.cr_nscat) then 
                        cr_iscat (cr_natoms + k) = cr_nscat + 1 
                        IF (k.eq.1) then 
                           cr_dw (cr_nscat + 1) = ins_adp 
                        ELSE 
                           cr_dw (cr_nscat + 1) = 0.0 
                        ENDIF 
                        cr_occ(cr_nscat + 1) = 1.00                   !! WORK OCC
                        cr_at_lis (cr_nscat + 1) = cname (k) 
                        cr_nscat = cr_nscat + 1 
                     ELSE 
                        cr_iscat (cr_natoms + k) = j 
                     ENDIF 
                  ENDDO 
               ELSE set_atom
!                                                                       
!     -------- Old object type                                          
!                                                                       
                  DO j = 1, 8 
                     ii = mole_cont (mole_off (kk) + j) 
                     cr_iscat (cr_natoms + j) = cr_iscat (ii) 
                  ENDDO 
               ENDIF set_atom
               ii = mole_off (mole_num_mole) + mole_len (mole_num_mole) 
               DO j = 1, 8 
                  mole_cont (ii + j) = cr_natoms + j 
               ENDDO 
               mole_off (i) = ii 
               mole_len (i) = 8 
               cr_natoms = cr_natoms + 8 
               mole_num_mole = mole_num_mole+1 
               mole_num_atom = mole_off (mole_num_mole) + mole_len (mole_num_mole) 
            ELSE enough_atom
               ier_num = - 10 
               ier_typ = ER_APPL 
            ENDIF enough_atom
         ELSE enough_type
            ier_num = - 66 
            ier_typ = ER_APPL 
         ENDIF enough_type
      ELSE enough_mole  
         ier_num = - 65 
         ier_typ = ER_APPL 
      ENDIF enough_mole  
!                                                                       
      END SUBROUTINE insert_domain                  
END MODULE insert_menu
