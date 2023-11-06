MODULE shear
!
IMPLICIT NONE
!
PUBLIC
!
CONTAINS
!+                                                                      
!     Generalized shear operations:                                     
!                                                                       
!     Applies a deformation tensor to the selected object.              
!     Right now, all references to the inverse shear matrix have        
!     been switched off. This would need further developping            
!                                                                       
!*****7*****************************************************************
      SUBROUTINE shear_menue
!-                                                                      
!     Main menu for generalized shear operations                        
!+                                                                      
      USE discus_config_mod 
      USE discus_allocate_appl_mod
      USE crystal_mod 
      USE metric_mod
      USE modify_mod
      USE molecule_mod 
      USE shear_mod 
      USE discus_show_menu
      USE update_cr_dim_mod
!
      USE ber_params_mod
      USE calc_expr_mod
      USE doact_mod 
      USE do_eval_mod
      USE do_wait_mod
      USE errlist_mod 
      USE get_params_mod
      USE learn_mod 
USE lib_help
USE lib_do_operating_mod
USE lib_echo
USE lib_errlist_func
USE lib_length
USE lib_macro_func
      USE class_macro_internal
      USE prompt_mod 
      USE sup_mod
USE precision_mod
USE str_comp_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER, PARAMETER :: MIN_PARA = 29  ! A command requires at leaset these no of parameters
      INTEGER maxw 
      LOGICAL lold 
!                                                                       
      PARAMETER (lold = .false.) 
!                                                                       
      CHARACTER(LEN=PREC_STRING), DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: cpara
      INTEGER            , DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: lpara
      REAL(KIND=PREC_DP) , DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: werte
!
      CHARACTER(len=9) befehl 
      CHARACTER(LEN=LEN(prompt)) :: orig_prompt
      CHARACTER(LEN=PREC_STRING) :: line, zeile
      INTEGER lp, length, lbef 
      INTEGER indxg, ianz, i, j, k 
      INTEGER indxc 
      INTEGER          :: nscat = 1
      INTEGER          :: nsite = 1
      LOGICAL lend, lspace 
      LOGICAL l_need_setup 
      LOGICAL lselect 
      REAL(kind=PREC_DP) :: hkl (3), h (3), x 
real(kind=PREC_DP), dimension(3), parameter  :: NULLV = (/ 0.0D0, 0.0D0, 0.0D0 /)
!                                                                       
!                                                                       
      DATA l_need_setup / .true. / 
!                                                                       
      maxw = MAX(MIN_PARA,MAXSCAT+1)
      lend = .false. 
      CALL no_error 
      orig_prompt = prompt
      prompt = prompt (1:len_str (prompt) ) //'/shear' 
!                                                                       
      DO while (.not.lend) 
      CALL get_cmd (line, length, befehl, lbef, zeile, lp, prompt) 
      IF (ier_num.eq.0) then 
         IF (line /= ' '      .and. line(1:1) /= '#' .and. &
             line /= char(13) .and. line(1:1) /= '!'        ) THEN
!                                                                       
!     ----search for "="                                                
!                                                                       
indxg = index (line, '=') 
IF (indxg.ne.0.AND..NOT. (str_comp (befehl, 'echo',   2, lbef, 4) ) &
              .AND..NOT. (str_comp (befehl, 'system', 2, lbef, 4) )    &
              .AND..NOT. (str_comp (befehl, 'help',   2, lbef, 4) .OR. &
                          str_comp (befehl, '?   ',   2, lbef, 4) )    &
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
                     line(1:length-1) = line(2:length)
                     line(length:length) = ' '
                     length = length - 1
                     CALL file_kdo(line, length)
                  ELSE 
                     ier_num = - 13 
                     ier_typ = ER_MAC 
                  ENDIF 
!                                                                       
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
               ELSEIF (str_comp (befehl, 'chemstry', 2, lbef, 8) ) then 
                  CALL show_chem 
!                                                                       
!------ ----Echo a string, just for interactive check in a macro 'echo' 
!                                                                       
               ELSEIF (str_comp (befehl, 'echo', 2, lbef, 4) ) then 
                  CALL echo (zeile, lp) 
!                                                                       
!      ---Evaluate an expression, just for interactive check 'eval'     
!                                                                       
               ELSEIF (str_comp (befehl, 'evaluate', 2, lbef, 8) ) then 
                  CALL do_eval (zeile, lp, .TRUE.) 
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
                     CALL do_hel ('discus shear '//zeile, lp) 
                  ENDIF 
!                                                                       
!------- -Operating System Kommandos 'syst'                             
!                                                                       
               ELSEIF (str_comp (befehl, 'system', 2, lbef, 6) ) then 
                  IF (zeile.ne.' ') then 
                     CALL do_operating (zeile (1:lp), lp) 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ELSE
!
!-------------SHEAR specific commands
!
!
                 IF( cr_nscat > SHEAR_MAXSCAT .or. mole_num_type > SHEAR_MAXSCAT) THEN
                    nscat = max ( cr_nscat, mole_num_type)
                    nsite = MAX(nsite, cr_ncatoms, SHEAR_MAXSITE)
                    CALL alloc_shear ( nscat, nsite )
                               IF ( ier_num < 0 ) THEN
                       RETURN
                    ENDIF
                 ENDIF
!                                                                       
!------  -----waiting for user input                                    
!                                                                       
                   IF (str_comp (befehl, 'wait', 3, lbef, 4) ) then 
                  CALL do_input (zeile, lp) 
!
!     ----Reset shear 'rese'                                               
!
               ELSEIF (str_comp (befehl, 'reset', 2, lbef, 5) ) then 
                  CALL shear_reset
!                                                                       
!     ----calculate a single shear operation                            
!                                                                       
               ELSEIF (str_comp (befehl, 'calc', 2, lbef, 4) ) then 
                  IF (l_need_setup) then 
                     CALL shear_setup 
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
                           CALL shear_ca_single (hkl, lspace) 
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
!     ----Enter the Eigen vectors 'eigen'                               
!                                                                       
               ELSEIF (str_comp (befehl, 'eigen', 3, lbef, 5) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0.and. (ianz.eq.4.or.ianz.eq.5) ) then 
!                                                                       
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        i = nint (werte (1) ) 
                        DO j = 1, 3 
                        shear_eigenv (j, i) = werte (j + 1) 
                        shear_eigent (i, j) = werte (j + 1) 
                        h (j) = werte (j + 1) 
                        ENDDO 
                        IF (ianz.eq.5) then 
                           shear_eigenw (i) = werte (5) 
                           lspace = .true. 
                        ELSE 
                           shear_eigenw (i) = do_blen (lspace, h, NULLV) 
                        ENDIF 
                        x = do_blen (lspace, h, NULLV) 
                        DO j = 1, 3 
                        shear_eigenv (j, i) = shear_eigenv (j, i)       &
                        / x                                             
                        shear_eigent (i, j) = shear_eigent (i, j)       &
                        / x                                             
                        ENDDO 
                        shear_input = SHEAR_EIGEN 
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
!     ----Select the reciprocal space direction of the shear axis 'hkl' 
!                                                                       
               ELSEIF (str_comp (befehl, 'hkl ', 2, lbef, 4) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0.and. (ianz.eq.3.or.ianz.eq.4) ) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        DO i = 1, 3 
                        shear_hkl (i) = werte (i) 
                        ENDDO 
!                       CALL trans (shear_hkl, cr_rten, shear_uvw, 3) 
                        shear_uvw = matmul( cr_rten, shear_hkl)
                        IF (ianz.eq.3) then 
                           lspace = .true. 
                           shear_length = do_blen (lspace, NULLV,        &
                           shear_uvw)                                   
                        ELSE 
                           shear_length = werte (4) 
                        ENDIF 
                        l_need_setup = .true. 
                        shear_input = SHEAR_PLANE 
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
               ELSEIF (str_comp (befehl, 'include', 1, lbef, 7) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (ianz.eq.2) then 
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxw)                                           
                        IF (ier_num.eq.0) then 
                           shear_sel_atom = .true. 
                           shear_start = nint (werte (1) ) 
                           shear_end = nint (werte (2) ) 
                           shear_incl = 'list' 
                        ENDIF 
                     ELSEIF (ianz.eq.1) then 
                        IF (str_comp (cpara (1) , 'all', 1, lpara (1) , &
                        3) ) then                                       
                           shear_sel_atom = .true. 
                           shear_start = 1 
                           shear_end = - 1 
                           shear_incl = 'all ' 
                        ELSEIF (str_comp (cpara (1) , 'env', 1, lpara ( &
                        1) , 3) ) then                                  
                           shear_sel_atom = .true. 
                           shear_start = 1 
                           shear_end = - 1 
                           shear_incl = 'env ' 
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
!     ----Enter the deformation matrix  'matrix'                        
!                                                                       
               ELSEIF (str_comp (befehl, 'matrix', 3, lbef, 6) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0.and. (ianz.eq.9) ) then 
!                                                                       
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        k = 0 
                        DO i = 1, 3 
                        DO j = 1, 3 
                        k = k + 1 
                        shear_mat (i, j) = werte (k) 
                        ENDDO 
                        ENDDO 
                        shear_input = SHEAR_MATRIX 
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
!     ----Select range of molecules ore generalized objects within      
!         crystal to be included 'mincl'|'oincl'                        
!                                                                       
               ELSEIF (str_comp (befehl, 'minclude', 3, lbef, 8) .or.      &
                       str_comp (befehl, 'oinclude', 3, lbef, 8) ) then        
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (ianz.eq.2) then 
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxw)                                           
                        IF (ier_num.eq.0) then 
                           shear_sel_atom = .false. 
                           shear_start = nint (werte (1) ) 
                           shear_end = nint (werte (2) ) 
                        ENDIF 
                        IF (str_comp (befehl, 'minclude', 3, lbef, 8) )    &
                        then                                            
                           shear_mode = SHEAR_MOLECULE 
                        ELSEIF (str_comp (befehl, 'oinclude', 3, lbef, 8) )&
                        then                                            
                           shear_mode = SHEAR_OBJECT 
                        ENDIF 
                     ELSEIF (ianz.eq.1) then 
                        IF (str_comp (cpara (1) , 'all', 1, lpara (1) , &
                        3) ) then                                       
                           shear_sel_atom = .false. 
                           shear_start = 1 
                           shear_end = - 1 
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
!     ----Select the mode of shear operation 'mode'                     
!                                                                       
               ELSEIF (str_comp (befehl, 'mode', 3, lbef, 4) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (ianz.eq.1.or.ianz.eq.2) then 
                        IF(str_comp(cpara (1), 'atom', 1, lpara(1), 4) ) then
                           shear_mode = SHEAR_ATOM 
                           l_need_setup = .true. 
                        ELSEIF(str_comp(cpara(1), 'molecule', 1, lpara ( 1) , 8) ) then                                  
                           shear_mode = SHEAR_MOLECULE 
                        ELSEIF(str_comp(cpara(1), 'object', 1, lpara ( 1) , 6) ) then                                  
                           shear_mode = SHEAR_OBJECT 
                           l_need_setup = .true. 
                        ELSE 
                           ier_num = - 6 
                           ier_typ = ER_COMM 
                        ENDIF 
                     ENDIF 
                  ENDIF 
!                                                                       
!     ----Work on Domains                                               
!                                                                       
               ELSEIF (str_comp (befehl, 'domain', 2, lbef, 4) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF(str_comp(cpara(1), 'select', 2, lpara(1), 6) .or. &
                        str_comp(cpara(1), 'deselect', 2, lpara(1) , 8) ) then                                     
                        line = ' ' 
                        indxc = index (zeile, ',') 
                        line = zeile (indxc + 1:lp) 
                        lp = len_str (line) 
!                                                                       
                        lselect = str_comp(cpara(1), 'select',   2, lpara(1) , 6) .or. &
                                  str_comp(cpara(1), 'deselect', 2, lpara(1) , 8)                   
                        CALL mole_select (line, lp, 0, SHEAR_MAXSCAT,   &
                        shear_latom, shear_sel_atom, lselect)
                        shear_mode = SHEAR_DOMAIN 
!                                                                       
                     ELSEIF (str_comp (cpara (1) , 'include', 3, lpara ( 1) , 7) ) then                                     
                        CALL del_params (1, ianz, cpara, lpara, maxw) 
                        IF (ianz.eq.2) then 
                           CALL ber_params (ianz, cpara, lpara, werte, maxw)                                        
                           IF (ier_num.eq.0) then 
                              shear_sel_atom = .false. 
                              shear_start = nint (werte (1) ) 
                              shear_end = nint (werte (2) ) 
                           ENDIF 
                        ELSEIF (ianz.eq.1) then 
                           IF (str_comp (cpara (1) , 'all', 1, lpara (1) , 3) ) then                                  
                              shear_sel_atom = .false. 
                              shear_start = 1 
                              shear_end = - 1 
                           ELSE 
                              ier_num = - 6 
                              ier_typ = ER_COMM 
                           ENDIF 
                        ELSE 
                           ier_num = - 6 
                           ier_typ = ER_COMM 
                        ENDIF 
                     ELSEIF (str_comp (cpara (1) , 'atoms', 2, lpara (1) , 4) ) then                                        
                        shear_dom_mode_atom = str_comp (cpara (2) , 'apply', 2, lpara (2) , 5)                      
                     ELSEIF (str_comp (cpara (1) , 'shape', 2, lpara (2) , 5) ) then                                        
                        shear_dom_mode_shape = str_comp (cpara (2) , 'apply', 2, lpara (2) , 5)                      
                     ENDIF 
                  ENDIF 
!                                                                       
!     ----Select/deselect molecules                                     
!                                                                       
               ELSEIF (str_comp (befehl, 'mselect', 2, lbef, 7) .or.       &
                       str_comp (befehl, 'mdeselect', 2, lbef, 9) .or.     &
                       str_comp (befehl, 'oselect', 2, lbef, 7) .or.       &
                       str_comp (befehl, 'odeselect', 2, lbef, 9) ) then                                       
!                                                                       
                  lselect = str_comp (befehl, 'mselect', 2, lbef, 7)       &
                        .or.str_comp (befehl, 'oselect', 2, lbef, 7)             
                  CALL mole_select (zeile, lp, 0, SHEAR_MAXSCAT,        &
                       shear_latom, shear_sel_atom, lselect)
                  IF (str_comp (befehl, 'mseliect', 3, lbef, 7)             &
                  .or.str_comp (befehl, 'mdeselect', 3, lbef, 9) ) then      
                     shear_mode = SHEAR_MOLECULE 
                  ELSEIF (str_comp (befehl, 'oselect', 3, lbef, 7)         &
                      .or.str_comp (befehl, 'odeselect', 3, lbef, 9) ) then      
                     shear_mode = SHEAR_OBJECT 
                  ENDIF 
!                                                                       
!     ----Select the origin of the shear operation  'origin'            
!                                                                       
               ELSEIF (str_comp (befehl, 'origin', 1, lbef, 6) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0.and. (ianz.eq.3.or.ianz.eq.4) ) then 
                     IF (ianz.eq.4) then 
                        shear_orig_mol = str_comp(cpara(4), 'molecule', 1, lpara(4), 8) .or. &
                                         str_comp(cpara(4), 'object',   1, lpara(4), 6)                               
                        ianz = 3 
                     ENDIF 
!                                                                       
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        DO i = 1, 3 
                        shear_orig (i) = werte (i) 
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
!     ----Enter the deformation matrix  'rmatrix'                       
!                                                                       
!         ELSEIF(str_comp(befehl,'rmatrix',3,lbef,7)) then              
!           call get_params(zeile,ianz,cpara,lpara,maxw,lp)             
!           if(ier_num.eq.0 .and. (ianz.eq. 9 )) then                   
!                                                                       
!             call ber_params(ianz,cpara,lpara,werte,maxw)              
!             if(ier_num.eq.0) then                                     
!               k = 0                                                   
!               do i=1,3                                                
!                 do j=1,3                                              
!                   k = k+1                                             
!                   shear_rmat(i,j) = werte(k)                          
!                 ENDDO                                                 
!               ENDDO                                                   
!               shear_input = SHEAR_RMATRIX                             
!               l_need_setup = .true.                                   
!             ELSE                                                      
!               ier_num = -6                                            
!               ier_typ = ER_COMM                                       
!             endif                                                     
!           ELSE                                                        
!             ier_num = -6                                              
!             ier_typ = ER_COMM                                         
!           endif                                                       
!                                                                   
!     ----run shear 'run'                                               
!                                                                       
               ELSEIF (str_comp (befehl, 'run ', 2, lbef, 4) ) then 
                  IF (l_need_setup) then 
                     CALL shear_setup 
                  ENDIF 
                  l_need_setup = .false. 
                  IF (shear_mode.eq.SHEAR_ATOM) then 
!                                                                       
!-----      --------Apply shear operation to atoms                      
!                                                                       
                     CALL shear_op_single 
                  ELSEIF (shear_mode.eq.SHEAR_MOLECULE) then 
!                                                                       
!-----      --------Apply shear operation to molecules                  
!                                                                       
                     CALL shear_mole_single 
                  ELSEIF (shear_mode.eq.SHEAR_OBJECT) then 
!                                                                       
!-----      --------Apply shear operation to objects                    
!                                                                       
                     CALL shear_obj_single 
                  ELSEIF (shear_mode.eq.SHEAR_DOMAIN) then 
!                                                                       
!-----      --------Apply shear operation to domains                    
!                                                                       
                     CALL shear_dom_single 
                  ENDIF 
                  CALL update_cr_dim 
!                                                                       
!     ----Select which atoms are copied to their image 'sele'           
!                                                                       
               ELSEIF (str_comp (befehl, 'select', 2, lbef, 6)            &
               .or.str_comp (befehl, 'deselect', 2, lbef, 7) ) then         
!                                                                       
                  CALL atom_select (zeile, lp, 0, SHEAR_MAXSCAT, shear_latom, &
                  shear_lsite, 0, SHEAR_MAXSITE,                              &
                  shear_sel_atom, &
                  lold, str_comp (befehl, 'select', 2, lbef, 6) )         
                  shear_mode = SHEAR_ATOM 
!                                                                       
!     ----show current parameters 'show'                                
!                                                                       
               ELSEIF (str_comp (befehl, 'show', 2, lbef, 4) ) then 
                  IF (l_need_setup) then 
                     CALL shear_setup 
                  ENDIF 
                  l_need_setup = .false. 
                  CALL shear_show 
!                                                                       
!     ----Select the shear vector parallel to the plane 'vector'        
!                                                                       
               ELSEIF (str_comp (befehl, 'vector', 2, lbef, 6) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0.and.ianz.eq.3) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        DO i = 1, 3 
                        shear_vector (i) = werte (i) 
                        ENDDO 
                        l_need_setup = .true. 
                        shear_input = SHEAR_PLANE 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
!                                                                       
!     ----Select the direct space direction of the shear axis 'uvw'     
!                                                                       
               ELSEIF (str_comp (befehl, 'uvw ', 1, lbef, 4) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0.and. (ianz.eq.3.or.ianz.eq.4) ) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        DO i = 1, 3 
                        shear_uvw (i) = werte (i) 
                        ENDDO 
!                       CALL trans (shear_uvw, cr_gten, shear_hkl, 3) 
                        shear_hkl = matmul( cr_gten, shear_uvw)
                        IF (ianz.eq.3) then 
                           lspace = .true. 
                           shear_length = do_blen (lspace, NULLV,        &
                           shear_uvw)                                   
                        ELSE 
                           shear_length = werte (4) 
                        ENDIF 
                        l_need_setup = .true. 
                        shear_input = SHEAR_PLANE 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_COMM 
               ENDIF 
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
                  ier_msg(1) = ' Error occured in shear menu'
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
      END SUBROUTINE shear_menue                          
!*****7*****************************************************************
      SUBROUTINE shear_show 
!-                                                                      
!     Shows current shear settings                                      
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE atom_name
      USE metric_mod
      USE molecule_mod 
      USE shear_mod 
!                                                                       
      USE prompt_mod 
!
      IMPLICIT none 
!                                                                       
      CHARACTER(9) at_lis (0:maxscat+1)! , at_name 
      INTEGER mol_lis (maxscat+1) 
      INTEGER i, j, k 
      LOGICAL lspace 
      REAL(kind=PREC_DP) :: det 
real(kind=PREC_DP) ::  u (3), v (3), w (3) 
real(kind=PREC_DP), dimension(3), parameter :: NULLV = (/ 0.0D0, 0.0D0, 0.0D0 /)
      REAL(KIND=PREC_DP) :: w12, w13, w23 
!                                                                       
      DATA lspace / .true. / 
!                                                                       
      IF (shear_input.eq.SHEAR_PLANE) then 
         WRITE (output_io, 3000) shear_uvw 
         WRITE (output_io, 3010) shear_hkl 
         WRITE (output_io, 3011) shear_length 
         WRITE (output_io, 3012) shear_vector 
      ELSEIF (shear_input.eq.SHEAR_EIGEN) then 
         WRITE (output_io, 4000) ( (shear_eigenv (j, i), j = 1, 3),     &
         shear_eigenw (i), i = 1, 3)                                    
         det = shear_eigenv (1, 1) * (shear_eigenv (2, 2) *             &
         shear_eigenv (3, 3) - shear_eigenv (2, 3) * shear_eigenv (3, 2)&
         ) + shear_eigenv (2, 1) * (shear_eigenv (3, 2) * shear_eigenv (&
         1, 3) - shear_eigenv (1, 2) * shear_eigenv (3, 3) ) +          &
         shear_eigenv (3, 1) * (shear_eigenv (1, 2) * shear_eigenv (2,  &
         3) - shear_eigenv (1, 3) * shear_eigenv (2, 2) )               
         IF (det.lt.0.0) then 
            WRITE (output_io, 4010) 'Negative determinant : ', det 
         ELSEIF (det.eq.0.0) then 
      WRITE (output_io, 4010) 'Zero     determinant : ', det 
         ELSEIF (det.gt.0.0) then 
            WRITE (output_io, 4010) 'Positive determinant : ', det 
         ENDIF 
         DO j = 1, 3 
         u (j) = shear_eigenv (j, 1) 
         v (j) = shear_eigenv (j, 2) 
         w (j) = shear_eigenv (j, 3) 
         ENDDO 
         w12 = do_bang (lspace, u, NULLV, v) 
         w13 = do_bang (lspace, u, NULLV, w) 
         w23 = do_bang (lspace, v, NULLV, w) 
         WRITE (output_io, 4020) w12, w13, w23 
      ENDIF 
      IF (shear_orig_mol) then 
         WRITE (output_io, 3020) shear_orig, ' rel.to molecule' 
      ELSE 
         WRITE (output_io, 3020) shear_orig, ' rel. to crystal' 
      ENDIF 
      WRITE (output_io, 3050) ( (shear_mat (i, j), j = 1, 4), i = 1, 3) 
!     write (output_io,3060) ((shear_rmat(i,j),j=1,3),i=1,3)            
!                                                                       
      IF (shear_mode.eq.SHEAR_ATOM) then 
         WRITE (output_io, 3100) 'Shear atoms' 
      ELSEIF (shear_mode.eq.SHEAR_MOLECULE) then 
         WRITE (output_io, 3100) 'Shear molecules' 
      ELSEIF (shear_mode.eq.SHEAR_OBJECT) then 
         WRITE (output_io, 3100) 'Shear objects' 
      ELSEIF (shear_mode.eq.SHEAR_DOMAIN) then 
         WRITE (output_io, 3100) 'Shear domains' 
      ENDIF 
!                                                                       
      IF (shear_new.and..not.shear_sel_atom) then 
         WRITE (output_io, 3110) 'Create new molecule type' 
      ELSE 
         WRITE (output_io, 3110) 'Keep molecule type' 
      ENDIF 
!                                                                       
!------ Working with atoms ...                                          
!                                                                       
      IF (shear_sel_atom) then 
!                                                                       
         j = 0 
         DO i = 0, cr_nscat 
         IF (shear_latom (i) ) then 
            j = j + 1 
            at_lis (j) = at_name (i) 
         ENDIF 
         ENDDO 
         WRITE (output_io, 3210) (at_lis (i), i = 1, j) 
!                                                                       
         IF (shear_incl.eq.'all ') then 
            WRITE (output_io, 3220) 
         ELSEIF (shear_incl.eq.'env ') then 
            WRITE (output_io, 3225) 
         ELSE 
            WRITE (output_io, 3230) shear_start, shear_end 
         ENDIF 
!                                                                       
!------ Working with molecules                                          
!                                                                       
      ELSE 
!                                                                       
         IF (shear_orig_mol) then 
            WRITE (output_io, 3250) 'Molecule' 
         ELSE 
            WRITE (output_io, 3250) 'Crystal' 
         ENDIF 
!                                                                       
         j = 0 
         DO i = 0, mole_num_type 
         IF (shear_latom (i) ) then 
            j = j + 1 
            mol_lis (j) = i 
         ENDIF 
         ENDDO 
         WRITE (output_io, 3300) (mol_lis (k), k = 1, j) 
!                                                                       
         IF (shear_end.eq. - 1) then 
            WRITE (output_io, 3310) 
         ELSE 
            WRITE (output_io, 3320) shear_start, shear_end 
         ENDIF 
      ENDIF 
!                                                                       
 3000 FORMAT    ( ' Generalized Shear Operation'/                       &
     &                   ' Axis in direct space     : ',3(2x,f11.4))    
 3010 FORMAT    ( ' Axis in reciprocal space : ',3(2x,f11.4)) 
 3011 FORMAT    ( ' Length of direct axis    : ', (2x,f11.4,' A')) 
 3012 FORMAT    ( ' Shear vector             : ',3(2x,f11.4)) 
 3020 FORMAT    ( ' Origin of shear element  : ',3(2x,f11.4),1x,a) 
 3050 FORMAT    ( ' Real space matrix        : ',4(2x,f11.4)/           &
     &                2( '                          : ',4(2x,f11.4)/))  
 3100 FORMAT    ( ' Mode of shear operation  : ',2x,a) 
 3110 FORMAT    ( ' Molecule status          : ',2x,a) 
 3210 FORMAT    ( ' Selected atom types      : ',2x,50(a9,1x)) 
 3220 FORMAT    ( ' Range of selected atoms  :   All atoms included') 
 3225 FORMAT    ( ' Range of selected atoms  :   Current environment') 
 3230 FORMAT    ( ' Range of selected atoms  : ',i9,' to ',i9) 
 3250 FORMAT    ( ' Given origin relative to : ',2x,a) 
 3300 FORMAT    ( ' Selected molecule types  : ',2x,50(i4,1x)) 
 3310 FORMAT( ' Range of sel. molecules  :   All molecules included') 
 3320 FORMAT    ( ' Range of sel. molecules  : ',i9,' to ',i9) 
 4000 FORMAT    ( ' Eigenvectors, -values    : ',3(2x,f11.4),4x,f11.4/  &
     &     2( '                          : ',3(2x,f11.4),4x,f11.4/))    
 4010 FORMAT    ( ' Determinant of Eigenvect.: ',2x,a,g15.5e2) 
 4020 FORMAT    ( ' <(1,2), <(1,3), <(2,3)   : ',3(2x,f11.4)) 
!                                                                       
      END SUBROUTINE shear_show                     
!*****7*****************************************************************
!
SUBROUTINE shear_setup 
!-                                                                      
!     Defines the deformation matrix in real and reciprocal space       
!+                                                                      
USE discus_config_mod 
USE crystal_mod 
USE metric_mod
USE shear_mod 
!
USE errlist_mod 
use matrix_mod
USE param_mod 
USE precision_mod
USE prompt_mod 
!                                                                       
IMPLICIT none 
!                                                                       
CHARACTER(LEN=PREC_STRING) ::line 
INTEGER :: i, j
INTEGER :: laenge 
LOGICAL,parameter :: lspace = .true.
!                                                                       
REAL(kind=PREC_DP) :: a (3, 3) 
REAL(kind=PREC_DP) :: fd, fn, sca 
REAL(kind=PREC_DP) :: u (3), v (3), w (3) 
REAL(kind=PREC_DP), dimension(3), parameter :: NULLV = (/ 0.0D0, 0.0D0, 0.0D0 /)
REAL(kind=PREC_DP) :: inv_eigenv (3, 3) 
!
IF (shear_input.eq.SHEAR_MATRIX) then 
!                                                                       
!     --initialize matrix and angle                                     
!                                                                       
   DO i = 1, 4 
      shear_mat (4, i) = 0.0 
      shear_mat (i, 4) = 0.0 
!         shear_rmat(i,4)= 0.0                                          
!         shear_rmat(4,i)= 0.0                                          
   ENDDO 
   shear_mat (4, 4) = 1.0 
!       shear_rmat(4,4)= 0.0                                            
!                                                                       
!     --Transform shear operation into reciprocal space                 
!                                                                       
!       do i=1,3                                                        
!         do j=1,3                                                      
!           a(i,j) = shear_mat(i,j)                                     
!         ENDDO                                                         
!       ENDDO                                                           
!c                                                                      
!     --do transformation q = gSg*                                      
!                                                                       
!       call matmulx (b,a,cr_rten)                                      
!       call matmulx (a,cr_gten,b)                                      
!       do i=1,3                                                        
!         do j=1,3                                                      
!           shear_rmat(i,j) = a(i,j)                                    
!         ENDDO                                                         
!       ENDDO                                                           
!     ELSEIF(shear_input.eq.SHEAR_RMATRIX) then                         
!c                                                                      
!     --initialize matrix and angle                                     
!                                                                       
!       do i=1,4                                                        
!         shear_mat (4,i)= 0.0                                          
!         shear_mat (i,4)= 0.0                                          
!         shear_rmat(i,4)= 0.0                                          
!         shear_rmat(4,i)= 0.0                                          
!       ENDDO                                                           
!       shear_mat (4,4)= 1.0                                            
!       shear_rmat(4,4)= 0.0                                            
!                                                                       
!     --Transform shear operation into direct space                     
!                                                                       
!       do i=1,3                                                        
!         do j=1,3                                                      
!           a(i,j) = shear_rmat(i,j)                                    
!         ENDDO                                                         
!       ENDDO                                                           
!                                                                       
!     --do transformation q = g*Sg                                      
!                                                                       
!       call matmulx (b,a,cr_gten)                                      
!       call matmulx (a,cr_rten,b)                                      
!       do i=1,3                                                        
!         do j=1,3                                                      
!           shear_mat(i,j) = a(i,j)                                     
!         ENDDO                                                         
!       ENDDO                                                           
ELSEIF (shear_input.eq.SHEAR_PLANE) then 
!                                                                       
!     --Calculate projections of unit vectors onto uvw                  
!                                                                       
   fd = do_blen (lspace, NULLV, shear_uvw) 
   laenge = 95 
   IF (fd.gt.0.0) then 
      DO j = 1, 3 
         v (j) = shear_uvw (j) / fd 
         v (j) = shear_uvw (j) / fd * shear_length 
      ENDDO 
      DO i = 1, 3 
         u (1) = 1.0D0 
         DO j = 1, 3 
            u (j) = 0.0D0 
         ENDDO 
         u (i) = 1.0D0 
         WRITE (line, 4000) v, u 
         CALL do_proj (line, laenge) 
         w (1) = res_para (1) 
         w (2) = res_para (2) 
         w (3) = res_para (3) 
         fn = do_blen (lspace, NULLV, w) 
         IF (res_para (7) .le.90) then 
            sca = fn / shear_length 
         ELSE 
            sca = - fn / shear_length 
         ENDIF 
         DO j = 1, 3 
            shear_mat (j, i) = u (j) + sca * shear_vector (j) 
         ENDDO 
      ENDDO 
   ENDIF 
ELSEIF (shear_input.eq.SHEAR_EIGEN) then 
!                                                                       
!     --Transform shear operation into direct space                     
!                                                                       
!                                                                       
!     --do transformation q = (Eigenvectors)*S*(inverse Eigenvectors)   
!                                                                       
!         DO i = 1, 3 
!         DO j = 1, 3 
!!        a (i, j) = 0.0 
!           b(i,j) = 0.0                                                
!         ENDDO 
!         a (i, i) = shear_eigenw (i) 
!!         b(i,i) = shear_eigenw(i)                                      
!        ENDDO 
      a = 0.0
      do i=1,3
         a(i,i) = shear_eigenw(i)
      enddo
!                                                                       
!        CALL invmat (inv_eigenv, shear_eigenv) 
!        CALL matmulx (b, a, inv_eigenv) 
!        CALL matmulx (a, shear_eigenv, b) 
!        DO i = 1, 3 
!        DO j = 1, 3 
!        shear_mat (i, j) = a (i, j) 
!        ENDDO 
!        ENDDO 
      shear_mat = 0.0D0
      call matinv(shear_eigenv, inv_eigenv)
      shear_mat(1:3,1:3) = matmul(shear_eigenv, matmul(a, inv_eigenv))
      shear_mat(4,4) = 1.0D0
!                                                                       
!     --Transform shear operation into reciprocal space                 
!                                                                       
!       do i=1,3                                                        
!         do j=1,3                                                      
!           a(i,j) = shear_mat(i,j)                                     
!         ENDDO                                                         
!       ENDDO                                                           
!                                                                       
!     --do transformation q = gSg*                                      
!                                                                       
!       call matmulx (b,a,cr_rten)                                      
!       call matmulx (a,cr_gten,b)                                      
!       do i=1,3                                                        
!         do j=1,3                                                      
!           shear_rmat(i,j) = a(i,j)                                    
!         ENDDO                                                         
!       ENDDO                                                           
ENDIF 
!     write (output_io,2000) ((shear_rmat(i,j),j=1,3),i=1,3)            
! 2000 FORMAT    (3(3(2x,f10.6)/)) 
 4000 FORMAT    (5(e15.8e2,','),e15.8e2) 
!                                                                       
END SUBROUTINE shear_setup                    
!
!*****7*****************************************************************
!
      SUBROUTINE shear_op_single 
!-                                                                      
!     Performs the actual shear operation, single result version        
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE atom_env_mod 
      USE shear_mod 
      USE errlist_mod 
use precision_mod
!
      IMPLICIT none 
!                                                                       
      INTEGER :: i, j, l 
      INTEGER :: i_start, i_end 
      REAL(kind=PREC_DP) :: ushear (4), ures (4) 
      REAL(kind=PREC_DP) :: werte (5) 
!                                                                       
      DATA ushear / 0.0D0, 0.0D0, 0.0D0, 1.0D0 / 
      DATA werte / 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0 / 
!                                                                       
!     Set the appropriate starting end ending number for the atoms      
!                                                                       
      i_start = shear_start 
      i_end = shear_end 
      IF (shear_incl.eq.'all ') then 
         i_end = cr_natoms 
      ELSEIF (shear_incl.eq.'env ') then 
         i_end = atom_env (0) 
      ENDIF 
!                                                                       
!     Apply shear operation to all atoms within selected range          
!                                                                       
      DO l = i_start, i_end 
      i = l 
      IF (shear_incl.eq.'env ') i = atom_env (l) 
!                                                                       
!     --Select atom if:                                                 
!       type has been selected and                                      
!            all atoms are selected                            or       
!            all microdomains have been selected               or       
!            a specific microdomain type has been selected     or       
!            microdomains are deselected and atom is outside            
!                                                                       
      IF (shear_latom (cr_iscat (i,1) ) ) then 
!                                                                       
!     ----Subtract origin                                               
!                                                                       
         DO j = 1, 3 
         ushear (j) = cr_pos (j, i) - shear_orig (j) 
         ENDDO 
!                                                                       
!-----      ----Apply shear operation                                   
!                                                                       
         ushear (4) = 1.0 
!        CALL trans (ushear, shear_mat, ures, 4) 
         ures = matmul(shear_mat, ushear)
!                                                                       
!     ----Add origin                                                    
!                                                                       
         DO j = 1, 3 
         werte (j + 1) = ures (j) + shear_orig (j) 
         ENDDO 
!                                                                       
!     ----replace original atom by its image                            
!                                                                       
         DO j = 1, 3 
         cr_pos (j, i) = werte (j + 1) 
         ENDDO 
      ENDIF 
      ENDDO 
!                                                                       
      END SUBROUTINE shear_op_single                
!*****7*****************************************************************
      SUBROUTINE shear_mole_single 
!-                                                                      
!     Performs the actual shear operation, single result version        
!     Operates on molecules                                             
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE molecule_mod 
      USE shear_mod
      USE errlist_mod 
use precision_mod
!
      IMPLICIT none 
!                                                                       
      INTEGER i, j, ii, l 
      INTEGER i_start, i_end 
      INTEGER :: imole_t = 1
      REAL(kind=PREC_DP) :: ushear (4), ures (4) 
      REAL(kind=PREC_DP) :: werte (5), use_orig (3) 
      REAL(kind=PREC_DP) :: diff (3) 
!                                                                       
      DATA ushear / 0.0D0, 0.0D0, 0.0D0, 1.0D0 / 
      DATA werte / 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0 / 
!                                                                       
!     Set the appropriate starting end ending number for the molecules  
!                                                                       
      i_start = shear_start 
      i_end = shear_end 
      IF (shear_end.eq. - 1) i_end = mole_num_mole 
!                                                                       
      IF (shear_new) then 
         IF (mole_num_type.lt.MOLE_MAX_TYPE) then 
            mole_num_type = mole_num_type+1 
         ELSE 
            ier_num = - 66 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
         imole_t = mole_num_type 
      ENDIF 
      IF (shear_orig_mol) then 
!                                                                       
!     Apply shear operation to all molecules within selected range      
!                                                                       
         DO l = i_start, i_end 
!                                                                       
!------ - Determine origin for shear operation                          
!                                                                       
         DO i = 1, 3 
         ii = mole_cont (mole_off (l) + 1) 
         use_orig (i) = shear_orig (i) + cr_pos (i, ii) 
         ENDDO 
!                                                                       
!     --Select atom if:                                                 
!       type has been selected                                          
!                                                                       
         IF (shear_latom (mole_type (l) ) .and.mole_char (l)            &
         .eq.MOLE_ATOM) then                                            
!                                                                       
!     ----Loop over all atoms in the molecule                           
!                                                                       
            DO ii = 1, mole_len (l) 
            i = mole_cont (mole_off (l) + ii) 
!                                                                       
!     ----- Subtract origin                                             
!                                                                       
            DO j = 1, 3 
            ushear (j) = cr_pos (j, i) - use_orig (j) 
            ENDDO 
!                                                                       
!-----      ----- Apply shear operation                                 
!                                                                       
            ushear (4) = 1.0 
!           CALL trans (ushear, shear_mat, ures, 4) 
            ures = matmul(shear_mat, ushear)
!                                                                       
!     ----- Add origin                                                  
!                                                                       
            DO j = 1, 3 
            werte (j) = ures (j) + use_orig (j) 
            ENDDO 
!                                                                       
!     ----- replace original atom by its image                          
!                                                                       
            DO j = 1, 3 
            cr_pos (j, i) = werte (j) 
            ENDDO 
            ENDDO 
!                                                                       
!----- ---- Set new molecule type if requested                          
!                                                                       
            IF (shear_new) then 
               mole_type (l) = imole_t 
            ENDIF 
         ENDIF 
         ENDDO 
!                                                                       
!     The origin is within the crystal. Shear the molecule position,    
!     not the internal vectors. Atom 1 of each molecule is sheared,     
!     all others are shifted by the same difference vector.             
!                                                                       
      ELSE 
!                                                                       
!     Apply shear operation to all molecules within selected range      
!                                                                       
         DO l = i_start, i_end 
!                                                                       
!------ - Determine origin for shear operation                          
!                                                                       
         DO i = 1, 3 
         use_orig (i) = shear_orig (i) 
         ENDDO 
!                                                                       
!     --Select atom if:                                                 
!       type has been selected                                          
!                                                                       
         IF (shear_latom (mole_type (l) ) .and.mole_char (l)            &
         .eq.MOLE_ATOM) then                                            
!                                                                       
!     ----Loop over all first atoms in the molecule                     
!                                                                       
            DO ii = 1, 1 
            i = mole_cont (mole_off (l) + ii) 
!                                                                       
!     ----- Subtract origin                                             
!                                                                       
            DO j = 1, 3 
            ushear (j) = cr_pos (j, i) - use_orig (j) 
            ENDDO 
!                                                                       
!-----      ----- Apply shear operation                                 
!                                                                       
            ushear (4) = 1.0 
!           CALL trans (ushear, shear_mat, ures, 4) 
            ures = matmul(shear_mat, ushear)
!                                                                       
!     ----- Add origin                                                  
!                                                                       
            DO j = 1, 3 
            werte (j) = ures (j) + use_orig (j) 
            diff (j) = werte (j) - cr_pos (j, i) 
            ENDDO 
!                                                                       
!     ----- replace original atom by its image                          
!                                                                       
            DO j = 1, 3 
            cr_pos (j, i) = werte (j) 
            ENDDO 
            ENDDO 
!                                                                       
!     ----Loop over all other atoms in the molecule                     
!                                                                       
            DO ii = 2, mole_len (l) 
            i = mole_cont (mole_off (l) + ii) 
            DO j = 1, 3 
            cr_pos (j, i) = cr_pos (j, i) + diff (j) 
            ENDDO 
            ENDDO 
!                                                                       
!----- ---- Set new molecule type if requested                          
!                                                                       
            IF (shear_new) then 
               mole_type (l) = imole_t 
            ENDIF 
         ENDIF 
         ENDDO 
      ENDIF 
!                                                                       
      END SUBROUTINE shear_mole_single              
!*****7*****************************************************************
      SUBROUTINE shear_obj_single 
!-                                                                      
!     Performs the actual shear operation, single result version        
!     Operates on objects                                               
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE molecule_mod 
      USE shear_mod 
      USE errlist_mod 
use precision_mod
!
      IMPLICIT none 
!                                                                       
      INTEGER i, j, ii, l, i0 
      INTEGER i_start, i_end 
      INTEGER :: imole_t = 1
      REAL(kind=PREC_DP) :: diff (3) 
      REAL(kind=PREC_DP) :: ushear (4), ures (4) 
      REAL(kind=PREC_DP) :: werte (5) 
!                                                                       
      DATA ushear / 0.0D0, 0.0D0, 0.0D0, 1.0D0 / 
      DATA werte / 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0 / 
!                                                                       
!     Set the appropriate starting end ending number for the molecules  
!                                                                       
      i_start = shear_start 
      i_end = shear_end 
      IF (shear_end.eq. - 1) i_end = mole_num_mole 
!                                                                       
      IF (shear_new) then 
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
      IF (shear_orig_mol) then 
!                                                                       
!     Apply shear operation to all molecules within selected range      
!                                                                       
         DO l = i_start, i_end 
!                                                                       
!------ - Determine origin for shear operation                          
!                                                                       
!                                                                       
!     --Select atom if:                                                 
!       type has been selected                                          
!                                                                       
         IF (shear_latom (mole_type (l) ) .and.mole_char (l)            &
         .ne.MOLE_ATOM) then                                            
!                                                                       
!     ----Loop over all pseudoatoms that represent the deformation      
!                                                                       
            i0 = mole_cont (mole_off (l) + 1) 
            DO ii = 2, 4 
            i = mole_cont (mole_off (l) + ii) 
!                                                                       
!     ----- Subtract origin                                             
!                                                                       
            DO j = 1, 3 
            ushear (j) = cr_pos (j, i) - cr_pos (j, i0) 
            ENDDO 
!                                                                       
!-----      ----- Apply shear operation                                 
!                                                                       
            ushear (4) = 1.0 
!           CALL trans (ushear, shear_mat, ures, 4) 
            ures = matmul(shear_mat, ushear)
!                                                                       
!     ----- Add origin                                                  
!                                                                       
            DO j = 1, 3 
            werte (j + 1) = ures (j) + cr_pos (j, i0) 
            ENDDO 
!                                                                       
!     ----- replace original atom by its image                          
!                                                                       
            DO j = 1, 3 
            cr_pos (j, i) = werte (j + 1) 
            ENDDO 
            ENDDO 
!                                                                       
!----- ---- Set new molecule type if requested                          
!                                                                       
            IF (shear_new) then 
               mole_type (l) = imole_t 
            ENDIF 
         ENDIF 
         ENDDO 
!                                                                       
!     Origin at crystal, do not shear individual object but positions   
!     of objects                                                        
!                                                                       
      ELSE 
!                                                                       
!     Apply shear operation to all molecules within selected range      
!                                                                       
         DO l = i_start, i_end 
!                                                                       
!------ - Determine origin for shear operation                          
!                                                                       
!                                                                       
!     --Select atom if:                                                 
!       type has been selected                                          
!                                                                       
         IF (shear_latom (mole_type (l) ) .and.mole_char (l)            &
         .ne.MOLE_ATOM) then                                            
!                                                                       
!     ----Loop over all pseudoatoms that represent the deformation      
!                                                                       
            i0 = mole_cont (mole_off (l) + 1) 
            DO ii = 1, 1 
            i = mole_cont (mole_off (l) + ii) 
!                                                                       
!     ----- Subtract origin                                             
!                                                                       
            DO j = 1, 3 
            ushear (j) = cr_pos (j, i) - shear_orig (j) 
            ENDDO 
!                                                                       
!-----      ----- Apply shear operation                                 
!                                                                       
            ushear (4) = 1.0 
!           CALL trans (ushear, shear_mat, ures, 4) 
            ures = matmul(shear_mat, ushear)
!                                                                       
!     ----- Add origin                                                  
!                                                                       
            DO j = 1, 3 
            werte (j) = ures (j) + shear_orig (j) 
            diff (j) = werte (j) - cr_pos (j, i) 
            ENDDO 
!                                                                       
!     ----- replace original atom by its image                          
!                                                                       
            DO j = 1, 3 
            cr_pos (j, i) = werte (j) 
            ENDDO 
            ENDDO 
!                                                                       
!     ----Apply shift to atoms 2 through 4                              
!                                                                       
            DO ii = 2, 4 
            i = mole_cont (mole_off (l) + ii) 
            DO j = 1, 3 
            cr_pos (j, i) = cr_pos (j, i) + diff (j) 
            ENDDO 
            ENDDO 
!                                                                       
!----- ---- Set new molecule type if requested                          
!                                                                       
            IF (shear_new) then 
               mole_type (l) = imole_t 
            ENDIF 
         ENDIF 
         ENDDO 
      ENDIF 
!                                                                       
      END SUBROUTINE shear_obj_single               
!*****7*****************************************************************
      SUBROUTINE shear_dom_single 
!-                                                                      
!     Performs the actual shear operation, single result version        
!     Operates on Microdomains                                          
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE molecule_mod 
      USE shear_mod 
      USE errlist_mod 
use precision_mod
!
      IMPLICIT none 
!                                                                       
      INTEGER :: i, j, ii, l, i0, m 
      INTEGER :: i_start, i_end 
      INTEGER :: imole_t=1
      REAL(kind=PREC_DP) :: diff (3) 
      REAL(kind=PREC_DP) :: ushear (4), ures (4) 
      REAL(kind=PREC_DP) :: werte (5) 
!                                                                       
      DATA ushear / 0.0D0, 0.0D0, 0.0D0, 1.0 / 
      DATA werte / 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0 / 
!                                                                       
!     Set the appropriate starting end ending number for the molecules  
!                                                                       
      i_start = shear_start 
      i_end = shear_end 
      IF (shear_end.eq. - 1) i_end = mole_num_mole 
!                                                                       
      IF (shear_new) then 
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
      IF (shear_orig_mol) then 
!                                                                       
!     Apply shear operation to all molecules within selected range      
!                                                                       
         DO l = i_start, i_end 
!                                                                       
!------ - Determine origin for shear operation                          
!                                                                       
!                                                                       
!     --Select atom if:                                                 
!       type has been selected                                          
!                                                                       
         IF (shear_latom (mole_type (l) ) .and.mole_char (l)            &
         .ne.MOLE_ATOM) then                                            
!                                                                       
!     ----Loop over all pseudoatoms that represent the deformation      
!                                                                       
            DO m = 1, 5, 4 
            i0 = mole_cont (mole_off (l) + m) 
            IF ( (m.eq.1.and.shear_dom_mode_atom) .or. (                &
            m.eq.5.and.shear_dom_mode_shape) ) then                     
               DO ii = 2, 4 
               i = mole_cont (mole_off (l) + ii + m - 1) 
!                                                                       
!     ----- Subtract origin                                             
!                                                                       
               DO j = 1, 3 
               ushear (j) = cr_pos (j, i) - cr_pos (j, i0) 
               ENDDO 
!                                                                       
!-----      ----- Apply shear operation                                 
!                                                                       
               ushear (4) = 1.0 
!              CALL trans (ushear, shear_mat, ures, 4) 
               ures = matmul(shear_mat, ushear)
!                                                                       
!     ----- Add origin                                                  
!                                                                       
               DO j = 1, 3 
               werte (j + 1) = ures (j) + cr_pos (j, i0) 
               ENDDO 
!                                                                       
!     ----- replace original atom by its image                          
!                                                                       
               DO j = 1, 3 
               cr_pos (j, i) = werte (j + 1) 
               ENDDO 
               ENDDO 
            ENDIF 
            ENDDO 
!                                                                       
!----- ---- Set new molecule type if requested                          
!                                                                       
            IF (shear_new) then 
               mole_type (l) = imole_t 
            ENDIF 
         ENDIF 
         ENDDO 
!                                                                       
!     Origin at crystal, do not shear individual object but positions   
!     of objects                                                        
!                                                                       
      ELSE 
!                                                                       
!     Apply shear operation to all molecules within selected range      
!                                                                       
         DO l = i_start, i_end 
!                                                                       
!------ - Determine origin for shear operation                          
!                                                                       
!                                                                       
!     --Select atom if:                                                 
!       type has been selected                                          
!                                                                       
         IF (shear_latom (mole_type (l) ) .and.mole_char (l)            &
         .ne.MOLE_ATOM) then                                            
!                                                                       
!     ----Loop over all pseudoatoms that represent the deformation      
!                                                                       
            DO m = 1, 5, 4 
            i0 = mole_cont (mole_off (l) + m) 
            IF ( (m.eq.1.and.shear_dom_mode_atom) .or. (                &
            m.eq.5.and.shear_dom_mode_shape) ) then                     
               DO ii = 1, 1 
               i = mole_cont (mole_off (l) + ii + m - 1) 
!                                                                       
!     ----- Subtract origin                                             
!                                                                       
               DO j = 1, 3 
               ushear (j) = cr_pos (j, i) - shear_orig (j) 
               ENDDO 
!                                                                       
!-----      ----- Apply shear operation                                 
!                                                                       
               ushear (4) = 1.0 
!              CALL trans (ushear, shear_mat, ures, 4) 
               ures = matmul(shear_mat, ushear)
!                                                                       
!     ----- Add origin                                                  
!                                                                       
               DO j = 1, 3 
               werte (j) = ures (j) + shear_orig (j) 
               diff (j) = werte (j) - cr_pos (j, i) 
               ENDDO 
!                                                                       
!     ----- replace original atom by its image                          
!                                                                       
               DO j = 1, 3 
               cr_pos (j, i) = werte (j) 
               ENDDO 
               ENDDO 
!                                                                       
!     ----Apply shift to atoms 2 through 4                              
!                                                                       
               DO ii = 2, 4 
               i = mole_cont (mole_off (l) + ii + m - 1) 
               DO j = 1, 3 
               cr_pos (j, i) = cr_pos (j, i) + diff (j) 
               ENDDO 
               ENDDO 
            ENDIF 
            ENDDO 
!                                                                       
!----- ---- Set new molecule type if requested                          
!                                                                       
            IF (shear_new) then 
               mole_type (l) = imole_t 
            ENDIF 
         ENDIF 
         ENDDO 
      ENDIF 
!                                                                       
      END SUBROUTINE shear_dom_single               
!*****7*****************************************************************
      SUBROUTINE shear_ca_single (uvw, lspace) 
!-                                                                      
!     Performs the actual shear operation, multiple copy version        
!     Only the input vector uvw is used in direct or reciprocal space   
!+                                                                      
      USE discus_config_mod 
      USE shear_mod 
!                                                                       
      USE param_mod 
      USE prompt_mod 
      USE errlist_mod 
use precision_mod
      IMPLICIT none 
!                                                                       
      INTEGER j 
      LOGICAL, intent(in) :: lspace 
!                                                                       
      REAL(kind=PREC_DP), intent(in) :: uvw (3) 
      REAL(kind=PREC_DP) :: ushear (4), ures (4) 
      REAL(kind=PREC_DP) :: werte (5) 
!                                                                       
      DATA ushear / 0.0D0, 0.0D0, 0.0D0, 1.0D0 / 
      DATA werte / 0.0D0, 0.0D0, 0.D00, 0.0D0, 0.0D0 / 
!                                                                       
!     real space part                                                   
!                                                                       
      IF (lspace) then 
!                                                                       
!     ----Subtract origin, if in real space                             
!                                                                       
         DO j = 1, 3 
         ushear (j) = uvw (j) - shear_orig (j) 
         ENDDO 
!                                                                       
!-----      --Apply shear operation                                     
!                                                                       
         ushear (4) = 1.0 
!        CALL trans (ushear, shear_mat, ures, 4) 
         ures = matmul(shear_mat, ushear)
!                                                                       
!     ----Add origin and store result                                   
!                                                                       
         DO j = 1, 3 
         res_para (j) = ures (j) + shear_orig (j) 
         ENDDO 
!                                                                       
!     ----Replace current vector by its image                           
!                                                                       
         DO j = 1, 3 
         ushear (j) = ures (j) 
         ENDDO 
      ELSE 
         CONTINUE 
!                                                                       
!     ----Subtract origin, if in real space                             
!                                                                       
!       do j=1,3                                                        
!         ushear(j) = uvw(j)                                            
!       ENDDO                                                           
!                                                                       
!-----      --Apply shear operation                                     
!                                                                       
!       ushear(4) = 0.0                                                 
!         call trans(ushear,shear_rmat,ures,4)                          
!                                                                       
!     ----Add origin and store result                                   
!                                                                       
!         do j=1,3                                                      
!           res_para(j) = ures(j)                                       
!         ENDDO                                                         
!                                                                       
!     ----Replace current vector by its image                           
!                                                                       
!         do j=1,3                                                      
!           ushear(j) = ures(j)                                         
!         ENDDO                                                         
      ENDIF 
!                                                                       
      res_para (0) = 3 
      WRITE (output_io, 3000) (res_para (j), j = 1, 3) 
!                                                                       
 3000 FORMAT    (' Result    : ',3(2x,f9.4)) 
      END SUBROUTINE shear_ca_single                
!
!*******************************************************************************
!
SUBROUTINE shear_reset
!
USE discus_allocate_appl_mod
USE shear_mod
! 
IMPLICIT NONE
!
INTEGER :: ik
!
CALL alloc_shear(1, 1)
!
IF(ALLOCATED(shear_latom))  shear_latom(:) = .TRUE.  ! (0:MAXSCAT)
IF(ALLOCATED(shear_lsite))  shear_lsite(:) = .TRUE.  ! (0:MAXSCAT)
shear_incl     = ' '
shear_sel_prop =  0
shear_start    =  1
shear_end      = -1
shear_mode     = SHEAR_OBJECT
shear_input    = SHEAR_MATRIX
shear_new      = .false.
shear_orig_mol = .false.
shear_sel_atom = .true.
shear_dom_mode_atom  = .true.
shear_dom_mode_shape = .true.
shear_hkl(:)         = (/0.,1.,0./)
shear_orig(:)        = 0.0
shear_vector(:)      = (/1.,0.,0./)
shear_length         = 0.0
shear_uvw(:)         = (/0.,1.,0./)
shear_mat(:,:)       = &
         RESHAPE((/1.,(0.,0.,0.,0.,1.,ik=1,3)/),SHAPE(shear_mat ))
shear_rmat(:,:)      = &
         RESHAPE((/1.,(0.,0.,0.,0.,1.,ik=1,3)/),SHAPE(shear_rmat))
shear_eigenv(:,:)    = &
         RESHAPE((/1.,(0.,0.,0.,1.,ik=1,2)/),SHAPE(shear_eigenv ))
shear_eigent(:,:)    = &
         RESHAPE((/1.,(0.,0.,0.,1.,ik=1,2)/),SHAPE(shear_eigent ))
shear_eigenw(  :)    = 1.0
!
END SUBROUTINE shear_reset
!
!*******************************************************************************
!
END MODULE shear
