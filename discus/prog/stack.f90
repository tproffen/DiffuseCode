MODULE stack_menu
!
CONTAINS
!
SUBROUTINE stack 
!                                                                       
!+                                                                      
!     This subroutines contain the generator for stacking faults.       
!     The menu defines the names of the files containing the different  
!     levels. The subroutine generates the sequence of levels and the   
!     list of corresponding origins. Optionally the whole crystal is    
!     created.                                                          
!-                                                                      
      USE discus_config_mod 
      USE discus_allocate_appl_mod
      USE crystal_mod 
      USE diffuse_mod 
      USE stack_mod 
      USE stack_rese_mod 
      USE spcgr_apply
!
      USE doact_mod 
      USE errlist_mod 
      USE learn_mod 
      USE class_macro_internal
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER, PARAMETER :: MIN_PARA = 99  ! A command requires at least these no of parameters
                                           ! Needs work as it should also be >= ST_MAXTYPE
      INTEGER maxw 
!                                                                       
      CHARACTER(LEN=1024), DIMENSION(MIN_PARA) :: cpara
      INTEGER            , DIMENSION(MIN_PARA) :: lpara
      REAL               , DIMENSION(MIN_PARA) :: werte
!
      CHARACTER(5) befehl 
      CHARACTER(50) prom 
      CHARACTER(1024) line, zeile
      INTEGER lp, length, lbef 
      INTEGER indxg, ianz, i, j, k 
      LOGICAL lend
      LOGICAL, SAVE :: linit    = .true. 
!
      INTEGER,SAVE  :: n_types  = 0 ! Current number of layer types
      INTEGER,SAVE  :: n_layers = 0 ! Current number of layers
      INTEGER,SAVE  :: n_qxy    = 0 ! Current number of points in reciprocal space
!                                                                       
      INTEGER len_str 
      LOGICAL str_comp 
!                                                                       
      maxw = MIN_PARA
      lend = .false. 
      CALL no_error 
!                                                                       
      IF (linit) then 
         CALL alloc_stack ( st_layer_increment, 1, 1, st_rot_status)
         n_types  = ST_MAXTYPE
         n_layers = ST_MAXLAYER
         n_qxy    = ST_MAXQXY
         linit    = .false.
      ENDIF 
!                                                                       
      DO while (.not.lend) 
      prom = prompt (1:len_str (prompt) ) //'/stack' 
      CALL get_cmd (line, length, befehl, lbef, zeile, lp, prom) 
      IF (ier_num.eq.0) then 
         IF (line /= ' '      .and. line(1:1) /= '#' .and. &
             line /= char(13) .and. line(1:1) /= '!'        ) THEN
!                                                                       
!     ----search for "="                                                
!                                                                       
            indxg = index (line, '=') 
      IF (indxg.ne.0.and..not. (str_comp (befehl, 'echo', 2, lbef, 4) ) &
     &.and..not. (str_comp (befehl, 'syst', 2, lbef, 4) ) .and..not. (st&
     &r_comp (befehl, 'help', 2, lbef, 4) .or.str_comp (befehl, '?   ', &
     &2, lbef, 4) ) ) then                                              
!                                                                       
!     ------evaluate an expression and assign the value to a variable   
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
!     ----Determine if average translation is calculated or set 'aver'  
!                                                                       
               ELSEIF (str_comp (befehl, 'aver', 1, lbef, 4) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (ianz.eq.3) then 
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxw)                                           
                        IF (ier_num.eq.0) then 
                           DO i = 1, 3 
                           st_t_aver (i) = werte (i) 
                           ENDDO 
                        ENDIF 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ENDIF 
!                                                                       
!     ----Read a single collumn of the correlation matrix 'ccolumn'     
!                                                                       
               ELSEIF (str_comp (befehl, 'ccol', 2, lbef, 4) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (st_ntypes.lt.ianz.and.ianz.le.ST_MAXTYPE+1)    &
                     then                                               
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxw)                                           
                        IF (ier_num.eq.0) then 
                           i = nint (werte (1) ) 
                           IF (0.lt.i.and.i.le.ST_MAXTYPE) then 
                              DO j = 1, ianz - 1 
                              st_corr (j, i) = werte (j + 1) 
                              ENDDO 
                           ELSE 
                              ier_num = - 13 
                              ier_typ = ER_APPL 
                           ENDIF 
                        ENDIF 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ENDIF 
!                                                                       
!     ----Read a single element of the correlation matrix 'celement'    
!                                                                       
               ELSEIF (str_comp (befehl, 'cele', 2, lbef, 4) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (ianz.eq.3) then 
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxw)                                           
                        IF (ier_num.eq.0) then 
                           i = nint (werte (1) ) 
                           j = nint (werte (2) ) 
                           IF (0.lt.i.and.i.le.ST_MAXTYPE.and.          &
                               0.lt.j.and.j.le.ST_MAXTYPE      ) then
                              st_corr (i, j) = werte (3) 
                           ELSE 
                              ier_num = - 13 
                              ier_typ = ER_APPL 
                           ENDIF 
                        ENDIF 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ENDIF 
!                                                                       
!     ----continues a macro 'continue'                                  
!                                                                       
               ELSEIF (str_comp (befehl, 'continue', 3, lbef, 8) ) then 
                  CALL macro_continue (zeile, lp) 
!                                                                       
!     ----create list of origins 'create'                               
!                                                                       
               ELSEIF (str_comp (befehl, 'crea', 3, lbef, 4) ) then 
!
!                 If necessary allocate array sizes
!
                  IF ( st_nlayer >  ST_MAXLAYER ) THEN
                     CALL alloc_stack (ST_MAXTYPE, st_nlayer, ST_MAXQXY, st_rot_status )
                  ENDIF
                  CALL do_stack_create 
!                                                                       
!     ----Read a single row of the correlation matrix 'crow'            
!                                                                       
               ELSEIF (str_comp (befehl, 'crow', 3, lbef, 4) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (st_ntypes.lt.ianz.and.ianz.le.ST_MAXTYPE+1) THEN
                        CALL ber_params (ianz, cpara, lpara, werte, maxw)
                        IF (ier_num.eq.0) then 
                           j = nint (werte (1) ) 
                           IF (0.lt.j.and.j.le.ST_MAXTYPE) then 
                              DO i = 1, ianz - 1 
                              st_corr (j, i) = werte (i + 1) 
                              ENDDO 
                           ELSE 
                              ier_num = - 13 
                              ier_typ = ER_APPL 
                           ENDIF 
                        ENDIF 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ENDIF 
!                                                                       
!     ----Define distribution type                                      
!                                                                       
               ELSEIF (str_comp (befehl, 'dist', 1, lbef, 4) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ianz.ge.1) then 
                     IF (str_comp (cpara (1) , 'matrix', 1, lpara (1) , &
                     6) ) then                                          
                        st_distr = ST_DIST_MATRIX 
                     ELSEIF (str_comp (cpara (1) , 'file', 1, lpara (1) &
                     , 4) ) then                                        
                        IF (ianz.gt.1) then 
                           CALL do_build_name (ianz, cpara, lpara,      &
                           werte, maxw, 2)                              
                           IF (ier_num.eq.0) then 
                              st_distr = ST_DIST_FILE 
                              st_infile = cpara (2) (1:lpara(2))
                              st_infile_l = lpara (2) 
                           ENDIF 
                        ELSE 
                           ier_num = - 6 
                           ier_typ = ER_COMM 
                        ENDIF 
                     ELSEIF (str_comp (cpara (1) , 'list', 1, lpara (1) &
                     , 4) ) then                                        
                        IF (ianz.gt.1) then 
                           CALL do_build_name (ianz, cpara, lpara,      &
                           werte, maxw, 2)                              
                           IF (ier_num.eq.0) then 
                              st_distr = ST_DIST_LIST 
                              st_infile = cpara (2) (1:lpara(2))
                              st_infile_l = lpara (2) 
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
!      ---Calculate the Fourier transform of the decorated              
!           s.f. 'fourier'                                              
!                                                                       
               ELSEIF (str_comp (befehl, 'four', 1, lbef, 4) ) then 
                  four_log = .true. 
                  CALL st_fourier (.false.)
                  four_was_run = .true.
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
                     CALL do_hel ('discus stack '//zeile, lp) 
                  ENDIF 
!                                                                       
!     ----read the name of a new layer type                'layer'      
!                                                                       
               ELSEIF (str_comp (befehl, 'laye', 1, lbef, 4) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     CALL do_build_name (ianz, cpara, lpara, werte,     &
                     maxw, 1)                                           
                     IF (ier_num.eq.0) then 
                        IF (ianz.eq.1) then 
!
!   If necessary increase array sizes
!
                           IF ( st_ntypes >=  ST_MAXTYPE ) THEN
                              i       = ST_MAXTYPE + st_layer_increment
                              CALL alloc_stack ( i, ST_MAXLAYER, ST_MAXQXY, st_rot_status )
                              n_types = ST_MAXTYPE
                           ENDIF
                           IF (st_ntypes.lt.ST_MAXTYPE) then 
                              st_ntypes = st_ntypes + 1 
                              st_internal(st_ntypes) = cpara(1)(1:8)=='internal'
                              st_layer (st_ntypes) = cpara (1) (1:lpara(1))
                              st_llayer (st_ntypes) = lpara (1) 
                              DO i = 1, st_nchem 
                              IF (cpara (1) .eq.st_layer_c (i) ) then 
                                 st_chem (st_ntypes) = i 
                                 GOTO 5000 
                              ENDIF 
                              ENDDO 
                              st_nchem = st_nchem + 1 
                              st_chem (st_ntypes) = st_nchem 
                              st_layer_c (st_nchem) = cpara (1) (1:lpara(1))
 5000                         CONTINUE 
                           ELSE 
                              ier_num = - 53 
                              ier_typ = ER_APPL 
                           ENDIF 
                        ELSE 
                           ier_num = - 6 
                           ier_typ = ER_COMM 
                        ENDIF 
                     ENDIF 
                  ENDIF 
!                                                                       
!     ----read the modulus for the translation             'modulus'    
!                                                                       
               ELSEIF (str_comp (befehl, 'modu', 1, lbef, 4) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (ianz.eq.6) then 
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxw)                                           
                        IF (ier_num.eq.0) then 
                           DO i = 1, 3 
                           st_mod (i, 1) = werte (i) 
                           st_mod (i, 2) = werte (i + 3) 
                           ENDDO 
                        ENDIF 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ENDIF 
!                                                                       
!     ----read the number of layers to be created in the                
!           crystal 'number'                                            
!                                                                       
               ELSEIF (str_comp (befehl, 'numb', 1, lbef, 4) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (ianz.eq.1) then 
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxw)                                           
                        IF (ier_num.eq.0) then 
                           i = nint (werte (1) ) 
!                          IF (i.lt.ST_MAXLAYER) then 
                              st_nlayer = i 
!                          ELSE 
!                             ier_num = - 54 
!                             ier_typ = ER_APPL 
!                          ENDIF 
                        ENDIF 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ENDIF 
!                                                                       
!     ----reset stacking faults 'rese'                                  
!                                                                       
               ELSEIF (str_comp (befehl, 'rese', 2, lbef, 4) ) then 
                  CALL do_stack_rese 
!                                                                       
!     ----Select values for random stacking faults 'random'             
!                                                                       
               ELSEIF (str_comp (befehl, 'random', 2, lbef, 6) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (str_comp (cpara (1) , 'prob', 1, lpara (1) , 4)&
                     ) then                                             
!                                                                       
!     ----------Random stacking fault probability                       
!                                                                       
                        IF (ianz.eq.2) then 
                           cpara (1) = '0' 
                           lpara (1) = 1 
                           CALL ber_params (ianz, cpara, lpara, werte,  &
                           maxw)                                        
                           IF (ier_num.eq.0) then 
                              st_prob = werte (2) 
                           ENDIF 
                        ENDIF 
                     ELSEIF (str_comp (cpara (1) , 'offset', 1, lpara ( &
                     1) , 6) ) then                                     
!                                                                       
!     ----------Random stacking fault offset                            
!                                                                       
                        IF (ianz.eq.4) then 
                           cpara (1) = '0' 
                           lpara (1) = 1 
                           CALL ber_params (ianz, cpara, lpara, werte,  &
                           maxw)                                        
                           IF (ier_num.eq.0) then 
                              DO i = 1, 3 
                              st_off (i) = werte (i + 1) 
                              ENDDO 
                           ENDIF 
                        ENDIF 
                     ELSEIF (str_comp (cpara (1) , 'sigma', 1, lpara (1)&
                     , 5) ) then                                        
!                                                                       
!     ----------Random stacking fault sigma for offset                  
!                                                                       
                        IF (ianz.eq.4) then 
                           cpara (1) = '0' 
                           lpara (1) = 1 
                           CALL ber_params (ianz, cpara, lpara, werte,  &
                           maxw)                                        
                           IF (ier_num.eq.0) then 
                              DO i = 1, 3 
                              st_sigma_off (i) = werte (i + 1) 
                              ENDDO 
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
!     ----Select values for rotational disorder    'rotation'           
!                                                                       
               ELSEIF (str_comp (befehl, 'rotation', 2, lbef, 8) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
!
!                    Allocate rotational stacking faults
!
                     CALL alloc_stack ( ST_MAXTYPE, ST_MAXLAYER, ST_MAXQXY, .true. )
!                    CALL alloc_stack ( st_layer_increment, 1, 1, .true.)
                     IF (str_comp (cpara (1) , 'mode', 1, lpara (1) , 4)&
                     ) then                                             
!                                                                       
!     ----------Mode of rotations "all" or only in case of              
!                  "stacking fault"                                     
!                                                                       
                        IF (ianz.eq.2) then 
                           IF (str_comp (cpara (2) , 'all', 1, lpara (2)&
                           , 3) ) then                                  
                              st_rot_mode = .true. 
                           ELSEIF (str_comp (cpara (2) , 'fault', 1,    &
                           lpara (2) , 5) ) then                        
                              st_rot_mode = .false. 
                           ELSE 
                              ier_num = - 6 
                              ier_typ = ER_COMM 
                           ENDIF 
                        ELSE 
                           ier_num = - 6 
                           ier_typ = ER_COMM 
                        ENDIF 
                     ELSEIF (str_comp (cpara (1) , 'axis', 1, lpara (1) &
                     , 4) ) then                                        
                        cpara (1) = '0' 
                        lpara (1) = 1 
!                                                                       
!     ----------Direction of in-plane vector 1                          
!                                                                       
                        IF (str_comp (cpara (2) , 'mod1', 4, lpara (2) ,&
                        4) ) then                                       
                           cpara (2) = '0' 
                           lpara (2) = 1 
!                                                                       
!     ----------Direction of in-plane vector 1                          
!                                                                       
                           IF (ianz.eq.5.or.ianz.eq.6) then 
                              IF (ianz.eq.5) then 
                                 cpara (6) = 'd' 
                                 lpara (6) = 1 
                              ENDIF 
                              ianz = 5 
                              CALL ber_params (ianz, cpara, lpara,      &
                              werte, maxw)                              
                              IF (ier_num.eq.0) then 
                                 st_rot_m1 (1) = werte (3) 
                                 st_rot_m1 (2) = werte (4) 
                                 st_rot_m1 (3) = werte (5) 
                              ENDIF 
                              IF (str_comp (cpara (6) , 'd', 1, lpara ( &
                              6) , 1) ) then                            
                                 st_rot_m1_lspace = .true. 
                              ELSEIF (str_comp (cpara (6) , 'r', 1,     &
                              lpara (6) , 1) ) then                     
                                 st_rot_m1_lspace = .false. 
                              ELSE 
                                 ier_num = - 6 
                                 ier_typ = ER_COMM 
                              ENDIF 
                           ELSE 
                              ier_num = - 6 
                              ier_typ = ER_COMM 
                           ENDIF 
                        ELSEIF (str_comp (cpara (2) , 'mod2', 4, lpara (&
                        2) , 4) ) then                                  
!                                                                       
!     ----------Direction of in-plane vector 2                          
!                                                                       
                           cpara (2) = '0' 
                           lpara (2) = 1 
                           IF (ianz.eq.5.or.ianz.eq.6) then 
                              IF (ianz.eq.5) then 
                                 cpara (6) = 'd' 
                                 lpara (6) = 1 
                              ENDIF 
                              ianz = 5 
                              CALL ber_params (ianz, cpara, lpara,      &
                              werte, maxw)                              
                              IF (ier_num.eq.0) then 
                                 st_rot_m2 (1) = werte (3) 
                                 st_rot_m2 (2) = werte (4) 
                                 st_rot_m2 (3) = werte (5) 
                              ENDIF 
                              IF (str_comp (cpara (6) , 'd', 1, lpara ( &
                              6) , 1) ) then                            
                                 st_rot_m2_lspace = .true. 
                              ELSEIF (str_comp (cpara (6) , 'r', 1,     &
                              lpara (6) , 1) ) then                     
                                 st_rot_m2_lspace = .false. 
                              ELSE 
                                 ier_num = - 6 
                                 ier_typ = ER_COMM 
                              ENDIF 
                           ELSE 
                              ier_num = - 6 
                              ier_typ = ER_COMM 
                           ENDIF 
                        ELSEIF (str_comp (cpara (2) , 'normal', 1,      &
                        lpara (2) , 6) ) then                           
!                                                                       
!     ----------Direction of in-plane vector 1                          
!                                                                       
                           cpara (2) = '0' 
                           lpara (2) = 1 
                           IF (ianz.eq.5.or.ianz.eq.6) then 
                              IF (ianz.eq.5) then 
                                 cpara (6) = 'd' 
                                 lpara (6) = 1 
                              ENDIF 
                              ianz = 5 
                              CALL ber_params (ianz, cpara, lpara,      &
                              werte, maxw)                              
                              IF (ier_num.eq.0) then 
                                 st_rot_no (1) = werte (3) 
                                 st_rot_no (2) = werte (4) 
                                 st_rot_no (3) = werte (5) 
                              ENDIF 
                              IF (str_comp (cpara (6) , 'd', 1, lpara ( &
                              6) , 1) ) then                            
                                 st_rot_no_lspace = .true. 
                              ELSEIF (str_comp (cpara (6) , 'r', 1,     &
                              lpara (6) , 1) ) then                     
                                 st_rot_no_lspace = .false. 
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
                     ELSEIF (str_comp (cpara (1) , 'sigma', 1, lpara (1)&
                     , 5) ) then                                        
!                                                                       
!     ----------Sigmas of oscillation around the three axes             
!                                                                       
                        IF (ianz.eq.3) then 
                           cpara (1) = '0' 
                           lpara (1) = 1 
                           IF (str_comp (cpara (2) , 'normal', 1, lpara &
                           (2) , 5) ) then                              
                              cpara (2) = '0' 
                              lpara (2) = 1 
                              CALL ber_params (ianz, cpara, lpara,      &
                              werte, maxw)                              
                              IF (ier_num.eq.0) then 
                                 st_rot_si_no = werte (3) 
                              ENDIF 
                           ELSEIF (str_comp (cpara (2) , 'mod1', 4,     &
                           lpara (2) , 4) ) then                        
                              cpara (2) = '0' 
                              lpara (2) = 1 
                              CALL ber_params (ianz, cpara, lpara,      &
                              werte, maxw)                              
                              IF (ier_num.eq.0) then 
                                 st_rot_si_m1 = werte (3) 
                              ENDIF 
                           ELSEIF (str_comp (cpara (2) , 'mod2', 4,     &
                           lpara (2) , 4) ) then                        
                              cpara (2) = '0' 
                              lpara (2) = 1 
                              CALL ber_params (ianz, cpara, lpara,      &
                              werte, maxw)                              
                              IF (ier_num.eq.0) then 
                                 st_rot_si_m2 = werte (3) 
                              ENDIF 
                           ELSE 
                              ier_num = - 6 
                              ier_typ = ER_COMM 
                           ENDIF 
                        ELSE 
                           ier_num = - 6 
                           ier_typ = ER_COMM 
                        ENDIF 
                     ELSEIF (str_comp (cpara (1) , 'status', 2, lpara ( &
                     1) , 6) ) then                                     
!                                                                       
!     ----------Switch status of rotational disorder "on" or "off"      
!                                                                       
                        IF (ianz.eq.2) then 
                           IF (str_comp (cpara (2) , 'on', 2, lpara (2) &
                           , 2) ) then                                  
                              st_rot_status = .true. 
                           ELSEIF (str_comp (cpara (2) , 'off', 2,      &
                           lpara (2) , 3) ) then                        
                              st_rot_status = .false. 
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
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
!                                                                       
!     ----run stacking faults 'run'                                     
!                                                                       
               ELSEIF (str_comp (befehl, 'run ', 2, lbef, 4) ) then 
                  CALL do_stack_fill 
!                                                                       
!     ----Set parameters                    'set'                       
!                                                                       
               ELSEIF (str_comp (befehl, 'set', 2, lbef, 3) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (str_comp (cpara (1) , 'aver', 1, lpara (1) , 4)&
                     ) then                                             
                        IF (ianz.eq.2) then 
                           CALL del_params (1, ianz, cpara, lpara, maxw) 
                           ianz = 1 
                           CALL ber_params (ianz, cpara, lpara, werte,  &
                           maxw)                                        
                           IF (ier_num.eq.0) then 
                              st_aver = werte (1) 
                           ENDIF 
                        ELSEIF (ianz.eq.1) then 
                           st_aver = 0.0 
                        ELSE 
                           ier_num = - 6 
                           ier_typ = ER_COMM 
                        ENDIF 
                     ELSEIF (str_comp (cpara (1) , 'modu', 1, lpara (1) &
                     , 4) ) then                                        
                        IF (ianz.eq.2) then 
                           IF (str_comp (cpara (2) , 'off', 2, lpara (2)&
                           , 3) ) then                                  
                              st_mod_sta = .false. 
                           ELSEIF (str_comp (cpara (2) , 'on', 2, lpara &
                           (2) , 2) ) then                              
                              st_mod_sta = .true. 
                           ELSE 
                              ier_num = - 6 
                              ier_typ = ER_COMM 
                           ENDIF 
                        ELSE 
                           ier_num = - 6 
                           ier_typ = ER_COMM 
                        ENDIF 
                     ELSEIF (str_comp (cpara (1) , 'tran', 1, lpara (1) &
                     , 4) ) then                                        
                        IF (ianz.eq.2) then 
                           IF (str_comp (cpara (2) , 'aver', 1, lpara ( &
                           2) , 4) ) then                               
                              st_tra_aver = .true. 
                           ELSEIF (str_comp (cpara (2) , 'fixed', 1,    &
                           lpara (2) , 5) ) then                        
                              st_tra_aver = .false. 
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
                  ENDIF 
!                                                                       
!     ----show current parameters 'show'                                
!                                                                       
               ELSEIF (str_comp (befehl, 'show', 2, lbef, 4) ) then 
                  WRITE (output_io, 3000) 
                  IF (st_aver.eq.0.0) then 
                     WRITE (output_io, 3005) 
                  ELSE 
                     WRITE (output_io, 3006) 
                  ENDIF 
                  IF (st_distr.eq.ST_DIST_MATRIX) then 
                     WRITE (output_io, 3007) 
                  ELSEIF (st_distr.eq.ST_DIST_FILE) then 
                     WRITE (output_io, 3008) st_infile (1:st_infile_l) 
                  ELSEIF (st_distr.eq.ST_DIST_LIST) then 
                     WRITE (output_io, 4008) st_infile (1:st_infile_l) 
                  ENDIF 
                  WRITE (output_io, 3009) st_ntypes 
                  WRITE (output_io, 3010) 
                  DO i = 1, st_ntypes 
                  WRITE (output_io, 3020) i, st_layer (i) (1:st_llayer (&
                  i) )                                                  
                  ENDDO 
                  WRITE (output_io, 3030) 
                  DO i = 1, st_ntypes 
                  DO j = 1, st_ntypes 
                  WRITE (output_io, 3040) i, j, (st_trans (i, j, k),    &
                  k = 1, 3), (st_sigma (i, j, k), k = 1, 3)             
                  ENDDO 
                  WRITE (output_io, * ) 
                  ENDDO 
                  IF (st_tra_aver) then 
                     WRITE (output_io, 3041) 
                  ELSE 
                     WRITE (output_io, 3042) (st_t_aver (i), i = 1, 3) 
                  ENDIF 
                  IF (st_mod_sta) then 
                     WRITE (output_io, 3043) 
                  ELSE 
                     WRITE (output_io, 3044) 
                  ENDIF 
                  WRITE (output_io, 3045) ( (st_mod (i, j), i = 1, 3),  &
                  j = 1, 2)                                             
                  WRITE (output_io, 3050) 
                  DO i = 1, st_ntypes 
                  WRITE (output_io, 3060) (st_corr (i, j), j = 1,       &
                  st_ntypes)                                            
                  ENDDO 
                  WRITE (output_io, 3070) st_nlayer 
                  WRITE (output_io, 3080) st_prob 
                  WRITE (output_io, 3090) st_off, st_sigma_off 
                  IF (st_rot_status) then 
                     WRITE (output_io, 3100) 
                  ELSE 
                     WRITE (output_io, 3101) 
                  ENDIF 
                  IF (st_rot_mode) then 
                     WRITE (output_io, 3110) 
                  ELSE 
                     WRITE (output_io, 3120) 
                  ENDIF 
                  IF (st_rot_no_lspace) then 
                     WRITE (output_io, 3140) st_rot_no, 'd',            &
                     st_rot_si_no                                       
                  ELSE 
                     WRITE (output_io, 3140) st_rot_no, 'r',            &
                     st_rot_si_no                                       
                  ENDIF 
                  IF (st_rot_m1_lspace) then 
                     WRITE (output_io, 3140) st_rot_m1, 'd',            &
                     st_rot_si_m1                                       
                  ELSE 
                     WRITE (output_io, 3140) st_rot_m1, 'r',            &
                     st_rot_si_m1                                       
                  ENDIF 
                  IF (st_rot_m2_lspace) then 
                     WRITE (output_io, 3140) st_rot_m2, 'd',            &
                     st_rot_si_m2                                       
                  ELSE 
                     WRITE (output_io, 3140) st_rot_m2, 'r',            &
                     st_rot_si_m2                                       
                  ENDIF 
!                                                                       
!     ----Select sigmas of translational vector 'sigma'                 
!                                                                       
               ELSEIF (str_comp (befehl, 'sigma', 2, lbef, 5) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (ianz.eq.5) then 
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxw)                                           
                        IF (ier_num.eq.0) then 
                           i = nint (werte (1) ) 
                           j = nint (werte (2) ) 
      IF (0.lt.i.and.i.le.ST_MAXTYPE.and.0.lt.j.and.j.le.ST_MAXTYPE) then
                              DO k = 1, 3 
                              st_sigma (i, j, k) = werte (k + 2) 
                              ENDDO 
                           ELSE 
                              ier_num = - 54 
                              ier_typ = ER_APPL 
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
!     ----Select translational vector 'trans'                           
!                                                                       
               ELSEIF (str_comp (befehl, 'trans', 1, lbef, 5) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (ianz.eq.5) then 
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxw)                                           
                        IF (ier_num.eq.0) then 
                           i = nint (werte (1) ) 
                           j = nint (werte (2) ) 
      IF (0.lt.i.and.i.le.ST_MAXTYPE.and.0.lt.j.and.j.le.ST_MAXTYPE) then
                              DO k = 1, 3 
                              st_trans (i, j, k) = werte (k + 2) 
                              ENDDO 
                           ELSE 
                              ier_num = - 54 
                              ier_typ = ER_APPL 
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
!                                                                       
!------  -----waiting for user input                                    
!                                                                       
               ELSEIF (str_comp (befehl, 'wait', 3, lbef, 4) ) then 
                  CALL do_input (zeile, lp) 
!                                                                       
!     ----Unknown command                                               
!                                                                       
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
 3000 FORMAT    (/30x,' Generalized Stacking Faults'/                   &
     &                   30x,' ==========================='/)           
 3005 FORMAT    (' Fourier Mode               : ','turbo Fourier') 
 3006 FORMAT    (' Fourier Mode               : ','turbo Fourier - <F>') 
 3007 FORMAT    (' Distribution of layer types: ','correlation matrix') 
 3008 FORMAT    (' Distribution of layer types: ','file',/,             &
     &                  '                 input file : ',a)             
 4008 FORMAT    (' Distribution of layer types: ',                      &
     &                  'file; origins taken from file',/,              &
     &                  '                 input file : ',a)             
 3009 FORMAT    (' Number of different Layers : ',  2x,i4   ) 
 3010 FORMAT    (' Names of Layers            : '           ) 
 3020 FORMAT    ('                            : ',i4,1x,a   ) 
 3030 FORMAT    (' Translations between layers:    Pairs  ',12x,        &
     &                  'Vectors /  Sigma')                             
 3040 FORMAT    ('                            : ',                      &
     &                   i4,2x,i4,1x,3(2x,f9.4)/41x,3(2x,f9.4))         
 3050 FORMAT    (' Correlation Matrix         : '           ) 
 3060 FORMAT    ('                            : ',8(2x,f9.4)) 
 3041 FORMAT    (' Average translation        : calculated') 
 3042 FORMAT    (' Average translation        : ',3(2x,f9.4)) 
 3043 FORMAT    (' Apply modulo function      : YES') 
 3044 FORMAT    (' Apply modulo function      : NO') 
 3045 FORMAT    (' Modulus of translation     : ',3(2x,f9.4)/           &
     &                  '                              ',3(2x,f9.4) )   
 3070 FORMAT    (/' Number of layers in crystal:',  2x,i4   ) 
 3080 FORMAT    (/' Prob.  of random faults    : ',  2x,f9.4 ) 
 3090 FORMAT    (' Offset of random faults    : ',3(2x,f9.4)/           &
     &                  ' Sigma  of offset, random   : ',3(2x,f9.4) )   
 3100 FORMAT    (/' Rotational disorder        : on') 
 3101 FORMAT    (/' Rotational disorder        : off') 
 3110 FORMAT    ('     affects                : ',                      &
     &                  '  all layers individually')                    
 3120 FORMAT    ('     affects                :  only stacking faults ') 
 3140 FORMAT    ('     axis in-plane 1,sigma  : ',3(2x,f9.4),2x,a1,2x,  &
     &                  f9.4 )                                          
!                                                                       
      END SUBROUTINE stack                          
!*****7*****************************************************************
      SUBROUTINE do_stack_create 
!-                                                                      
!     Creates the list of origins for the crystal with stacking faults  
!+                                                                      
      USE discus_config_mod 
      USE discus_allocate_appl_mod
      USE crystal_mod 
      USE molecule_mod 
      USE stack_mod 
      USE stack_cr_mod 
      USE structur
      USE spcgr_apply
      USE symm_sup_mod
      USE tensors_mod
      USE trafo_mod
      USE errlist_mod 
      USE random_mod
      USE param_mod 
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER i, j, k 
      INTEGER         :: nlayers
      INTEGER         :: n_mole  ! number of molecules in input file
      INTEGER         :: n_type  ! number of molecule types in input file
      INTEGER         :: n_atom  ! number of molecule atoms in input file
      LOGICAL lprev 
      REAL prob (ST_MAXTYPE) 
      REAL prob_n (ST_MAXTYPE) 
      REAL ptot_n 
      REAL u (3), v (3) 
      REAL r1, r2 
      REAL st_trans_cur (3) 
!                                                                       
      REAL gasdev 
      REAL ran1 
!                                                                       
      CALL no_error 
!                                                                       
!     reset counter for number of each layer type                       
!                                                                       
      DO i = 1, ST_MAXTYPE 
         st_number (i) = 0 
      ENDDO 
!                                                                       
!     Read layer types from file or create from probability matrix      
!                                                                       
      IF (st_distr.eq.ST_DIST_MATRIX) then 
!                                                                       
!       Set up initial probabilities                                    
!                                                                       
         CALL st_init_prob (ptot_n, prob_n) 
         IF (ptot_n.eq.0.0) then 
            ier_num = - 39 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
!                                                                       
!       The first layer is determined by random weighted choice         
!                                                                       
         lprev = .false. 
         i = 1 
         CALL stack_prob (ptot_n, prob, prob_n, lprev, i) 
         CALL stack_type (prob, i) 
         st_origin (1, 1) = 0.0 
         st_origin (2, 1) = 0.0 
         st_origin (3, 1) = 0.0 
         st_number (st_type (1) ) = 1 
      ELSEIF (st_distr.eq.ST_DIST_FILE.or.st_distr.eq.ST_DIST_LIST) THEN
!                                                                       
!     --Stacking fault origins are read from file                       
!                                                                       
         CALL test_file ( st_infile, st_natoms, st_nscat, n_mole,  &
                          n_type, n_atom, -1, .true. )
         IF(st_natoms > ST_MMAX .or. st_nscat > ST_MAX_SCAT) THEN
            CALL alloc_stack_crystal (st_nscat, st_natoms)
            IF ( ier_num /= 0 ) RETURN
         ENDIF
!
!        IF more atoms have been read than layers were allocated, 
!
         IF ( st_nlayer > ST_MAXLAYER .or.        &
              st_natoms > ST_MAXLAYER      ) THEN
            nlayers = MAX ( st_nlayer, st_natoms )
            CALL alloc_stack (ST_MAXTYPE, nlayers, MAXQXY, st_rot_status )
         ENDIF
!
!        IF more molecules have been read than were allocated
!
               IF(n_mole>MOLE_MAX_MOLE .or. n_type>MOLE_MAX_TYPE .or.   &
                  n_atom>MOLE_MAX_ATOM                          ) THEN
                  n_mole = MAX(n_mole +20 ,MOLE_MAX_MOLE)
                  n_type = MAX(n_type +10 ,MOLE_MAX_TYPE)
                  n_atom = MAX(n_atom +200,MOLE_MAX_ATOM)
                  CALL alloc_molecule(1, 1,n_mole,n_type,n_atom)
                  IF ( ier_num /= 0 ) RETURN
               ENDIF
!
         CALL stack_dist_file ()
         IF ( ier_num /= 0 ) THEN
            RETURN
         ENDIF
!        CALL stack_dist_file (cr_spcgr, cr_a0, cr_win, cr_dim, st_name,&
!        st_spcgr, st_a0, st_win, st_natoms, st_nscat, st_dw, st_at_lis,&
!        st_pos, st_iscat, st_dim, sa_natoms, sa_at_lis, sa_dw, sa_pos, &
!        sa_iscat)                                                      
         st_type (1) = st_iscat (1) 
         st_nlayer = st_natoms 
                                                                    
!DBG                                                                    
!DBG      write (output_io,*) 'st_natoms ', st_natoms                   
         IF (st_distr.eq.ST_DIST_LIST) then 
            DO i = 1, st_natoms 
            st_type (i) = st_iscat (i) 
            st_number (st_type (i) ) = st_number (st_type (i) ) + 1 
            DO j = 1, 3 
            st_origin (j, i) = st_pos (j, i) 
            ENDDO 
!DBG      write (output_io,2000) ' Origins',st_type(i)                  
!DBG     &                         ,st_origin(1,i),st_origin(2,i),      
!DBG     &                                 st_origin(3,i)               
            ENDDO 
            st_nlayer = st_natoms 
            ier_num = 0 
            RETURN 
         ENDIF 
      ENDIF 
!                                                                       
!     Determine average translation                                     
!                                                                       
      IF (st_tra_aver) then 
         DO i = 1, 3 
         st_mod (i, 3) = 0.0 
         ENDDO 
         DO i = 1, st_ntypes 
         DO j = 1, st_ntypes 
         DO k = 1, 3 
         st_mod (k, 3) = st_mod (k, 3) + st_corr (i, j) * st_trans (i,  &
         j, k)                                                          
         ENDDO 
         ENDDO 
         ENDDO 
         DO i = 1, 3 
         st_mod (i, 3) = st_mod (i, 3) / ptot_n 
         ENDDO 
         IF (st_mod (1, 3) .eq.0.0.and.st_mod (2, 3) .eq.0.0.and.st_mod &
         (3, 3) .eq.0.0) then                                           
            ier_num = - 57 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
      ELSE 
         DO i = 1, 3 
         st_mod (i, 3) = st_t_aver (i) 
         ENDDO 
      ENDIF 
!                                                                       
!     Calculate inverse modulus matrix                                  
!                                                                       
      DO i = 1, 3 
      DO j = 1, 3 
      st_inv (i, j) = st_mod (i, j) 
      ENDDO 
      ENDDO 
      CALL invmat (st_inv, st_mod) 
      IF (ier_num.eq. - 1) then 
         CALL no_error 
         ier_num = - 56 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
!                                                                       
!     --Initialise rotation disorder                                    
!                                                                       
      i = 1 
      IF (st_rot_status) then 
         st_rot_ang_no (i) = gasdev (st_rot_si_no) 
         st_rot_ang_m1 (i) = gasdev (st_rot_si_m1) 
         st_rot_ang_m2 (i) = gasdev (st_rot_si_m2) 
      ENDIF 
!                                                                       
!     Loop over all layers                                              
!                                                                       
!DBG                                                                    
      WRITE (output_io, * ) ' Start Loop' , st_nlayer
      WRITE (output_io, 2000) ' Origins', st_type (1) , st_origin (1, 1)&
      , st_origin (2, 1) , st_origin (3, 1)                             
      DO i = 2, st_nlayer 
!                                                                       
!     --Determine whether random or regular stacking fault              
!                                                                       
      IF (st_distr.eq.ST_DIST_MATRIX) then 
         IF (ran1 (idum) .lt.st_prob) then 
!                                                                       
!     ----Determine layer type from probabilities disregarding          
!           previous layer                                              
!                                                                       
            lprev = .false. 
            CALL stack_prob (ptot_n, prob, prob_n, .false., i) 
            CALL stack_type (prob, i) 
         ELSE 
!                                                                       
!     ----Determine layer type from probabilities using last layer      
!                                                                       
            lprev = .true. 
            CALL stack_prob (ptot_n, prob, prob_n, lprev, i - 1) 
            CALL stack_type (prob, i) 
         ENDIF 
      ELSEIF (st_distr.eq.ST_DIST_FILE) then 
!                                                                       
!     --Stacking fault types are read from file                         
!                                                                       
         lprev = .true. 
         st_type (i) = st_iscat (i) 
      ENDIF 
!                                                                       
!     --Determine current translation vector                            
!                                                                       
      st_trans_cur (1) = st_trans (st_type (i - 1), st_type (i),        &
      1) + gasdev (st_sigma (st_type (i - 1), st_type (i), 1) )         
      st_trans_cur (2) = st_trans (st_type (i - 1), st_type (i),        &
      2) + gasdev (st_sigma (st_type (i - 1), st_type (i), 2) )         
      st_trans_cur (3) = st_trans (st_type (i - 1), st_type (i),        &
      3) + gasdev (st_sigma (st_type (i - 1), st_type (i), 3) )         
!DBG                                                                    
!DBG      if(i.eq.2) then                                               
!DBG      write (output_io,*) ' Translation vector',i,st_type(i)        
!DBG     &    ,st_trans_cur(1),st_trans_cur(2),                         
!DBG     &    st_trans_cur(3)                                           
!DBG      endif                                                         
!                                                                       
!     --Initialise rotation disorder                                    
!                                                                       
      IF (st_rot_status) then 
!DBG                                                                    
!         IF (i.eq.2) then 
!            WRITE (output_io, * ) ' im rotation mode' 
!         ENDIF 
!                                                                       
!     ----Set new rotation angles for all layers or                     
!         only if the last layer was different than the current         
!                                                                       
         IF (st_rot_mode) then 
            st_rot_ang_no (i) = gasdev (st_rot_si_no) 
            st_rot_ang_m1 (i) = gasdev (st_rot_si_m1) 
            st_rot_ang_m2 (i) = gasdev (st_rot_si_m2) 
         ELSE 
            IF (st_type (i) .ne.st_type (i - 1) ) then 
               st_rot_ang_no (i) = gasdev (st_rot_si_no) 
               st_rot_ang_m1 (i) = gasdev (st_rot_si_m1) 
               st_rot_ang_m2 (i) = gasdev (st_rot_si_m2) 
            ELSE 
               st_rot_ang_no (i) = st_rot_ang_no (i - 1) 
               st_rot_ang_m1 (i) = st_rot_ang_m1 (i - 1) 
               st_rot_ang_m2 (i) = st_rot_ang_m2 (i - 1) 
            ENDIF 
         ENDIF 
!                                                                       
!     ----Rotate translation vector around normal                       
!                                                                       
         IF (st_rot_si_no.gt.0.0) then 
            CALL stack_rot_setup ('normal', 6, 1, i - 1, .false.) 
            CALL symm_setup 
            CALL symm_ca_single (st_trans_cur, .true., .false.) 
            DO j = 1, 3 
            st_trans_cur (j) = res_para (j) 
            ENDDO 
         ENDIF 
!                                                                       
!     ----Rotate translation vector around mod1                         
!                                                                       
         IF (st_rot_si_m1.gt.0.0) then 
            CALL stack_rot_setup ('mod1', 4, 1, i - 1, .false.) 
            CALL symm_setup 
            CALL symm_ca_single (st_trans_cur, .true., .false.) 
            DO j = 1, 3 
            st_trans_cur (j) = res_para (j) 
            ENDDO 
         ENDIF 
!                                                                       
!     ----Rotate translation vector around mod2                         
!                                                                       
         IF (st_rot_si_m2.gt.0.0) then 
            CALL stack_rot_setup ('mod2', 4, 1, i - 1, .false.) 
            CALL symm_setup 
            CALL symm_ca_single (st_trans_cur, .true., .false.) 
            DO j = 1, 3 
            st_trans_cur (j) = res_para (j) 
            ENDDO 
         ENDIF 
      ENDIF 
!                                                                       
!     --Add translation vector to last origin                           
!                                                                       
      st_origin (1, i) = st_origin (1, i - 1) + st_trans_cur (1) 
      st_origin (2, i) = st_origin (2, i - 1) + st_trans_cur (2) 
      st_origin (3, i) = st_origin (3, i - 1) + st_trans_cur (3) 
      IF (.not.lprev) then 
!DBG                                                                    
!         IF (i.eq.2) then 
!            WRITE (output_io, * ) ' im random mode' 
!         ENDIF 
!                                                                       
!     ----This is a random stacking fault, add random position          
!                                                                       
         r1 = ran1 (idum) 
         r2 = ran1 (idum) 
         DO j = 1, 3 
         st_origin (j, i) = st_origin (j, i) + r1 * st_mod (j, 1)       &
         + r2 * st_mod (j, 2) + st_off (j) + gasdev (st_sigma_off (j) ) 
         ENDDO 
      ENDIF 
!                                                                       
!     --If necessary apply modulo function                              
!                                                                       
      IF (st_mod_sta) then 
         DO j = 1, 3 
         v (j) = st_origin (j, i) - st_origin (j, 1) 
         ENDDO 
         CALL trans (v, st_inv, u, 3) 
         DO j = 1, 3 
         st_origin (j, i) = st_origin (j, i) - float (int (u (1) ) )    &
         * st_mod (j, 1) - float (int (u (2) ) ) * st_mod (j, 2)        
         ENDDO 
      ENDIF 
!                                                                       
!     --Increment number of layers of this type                         
!                                                                       
      st_number (st_type (i) ) = st_number (st_type (i) ) + 1 
      WRITE (output_io, 2000) ' Origins', st_type (i) , st_origin (1, i)&
      , st_origin (2, i) , st_origin (3, i)                             
      ENDDO 
!                                                                       
!     do i=1,st_ntypes                                                  
!       write (output_io,2020) i,st_number(i)                           
!     ENDDO                                                             
 2000 FORMAT    (a,i4,3f11.4) 
      END SUBROUTINE do_stack_create                
!*****7*****************************************************************
      SUBROUTINE st_init_prob (ptot_n, prob_n) 
!-                                                                      
!     Determines the probabilities from the correlation matrix          
!+                                                                      
      USE discus_config_mod 
      USE stack_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER i, j 
      REAL prob_n (ST_MAXTYPE) 
      REAL ptot_n 
!                                                                       
      ptot_n = 0.0 
      DO i = 1, st_ntypes 
      prob_n (i) = 0.0 
      DO j = 1, st_ntypes 
      prob_n (i) = prob_n (i) + st_corr (j, i) 
      ptot_n = ptot_n + st_corr (j, i) 
      ENDDO 
      ENDDO 
!                                                                       
      END SUBROUTINE st_init_prob                   
!*****7*****************************************************************
      SUBROUTINE stack_prob (ptot_n, prob, prob_n, lprev, l) 
!-                                                                      
!            Set up the cummulative probabilities.                      
!     if a previous layer exists, add probabilities                     
!     according to the correlation matrix, ELSE                         
!     an average probality for the                                      
!     specific layer type i.e. sum corr(i,m) /st_type                   
!     is added.                                                         
!+                                                                      
      USE discus_config_mod 
      USE stack_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER m, l 
      LOGICAL lprev 
      REAL prob (ST_MAXTYPE), ptot_n 
      REAL prob_n (ST_MAXTYPE) 
!                                                                       
      IF (lprev) then 
         DO m = 1, st_ntypes 
         prob (m) = st_corr (st_type (l), m) 
         ENDDO 
      ELSE 
         DO m = 1, st_ntypes 
         prob (m) = prob_n (m) / ptot_n 
         ENDDO 
      ENDIF 
!                                                                       
      END SUBROUTINE stack_prob                     
!*****7*****************************************************************
      SUBROUTINE stack_type (prob, i) 
!-                                                                      
!     Determine type of layer from weighted probabilities               
!+                                                                      
      USE discus_config_mod 
      USE stack_mod 
      USE random_mod
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      REAL prob (ST_MAXTYPE) 
!                                                                       
      INTEGER i, m 
      REAL ad, s, ptot 
      REAL ran1 
!                                                                       
      ptot = 0.0 
      DO m = 1, st_ntypes 
      ptot = ptot + prob (m) 
      ENDDO 
      ad = ran1 (idum) * ptot 
      s = prob (1) 
      m = 1 
!                                                                       
      DO while (s.lt.ad.and.m.lt.st_ntypes) 
      m = m + 1 
      s = s + prob (m) 
      ENDDO 
      st_type (i) = m 
!     write (output_io,*) ' st Type(i) ',i,st_type(i)                   
!                                                                       
      END SUBROUTINE stack_type                     
!*****7*****************************************************************
      SUBROUTINE do_stack_fill 
!-                                                                      
!     Creates the crystal with stacking faults                          
!+                                                                      
      USE discus_config_mod 
      USE discus_allocate_appl_mod
      USE crystal_mod 
      USE gen_add_mod 
      USE sym_add_mod 
      USE molecule_mod 
      USE read_internal_mod
      USE discus_save_mod 
      USE stack_mod  
      USE structur  
      USE symm_sup_mod
      USE spcgr_apply
      USE update_cr_dim_mod
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
      INTEGER ist 
      PARAMETER (ist = 7) 
!                                                                       
       
!                                                                       
      INTEGER          :: natoms, max_natoms
      INTEGER          :: nscats, max_nscats
      INTEGER          ::         max_n_mole
      INTEGER          ::         max_n_type
      INTEGER          ::         max_n_atom
      INTEGER i, j, k
      INTEGER iatom 
      INTEGER         :: n_mole  ! number of molecules in input file
      INTEGER         :: n_type  ! number of molecule types in input file
      INTEGER         :: n_atom  ! number of molecule atoms in input file
      LOGICAL lread, lout 
      LOGICAL           :: need_alloc = .false. 
      INTEGER, EXTERNAL :: len_str
!                                                                       
!                                                                       
!     If there are any layers in the crystal read each layer            
!                                                                       
      IF ( MAXVAL(st_number) == 0 ) then
         ier_num = -55
         ier_typ = ER_APPL
         ier_msg(1) = 'No layers have been created yet'
         ier_msg(2) = 'Use the ''create'' command prior to ''run'' '
         RETURN
      ENDIF
!
more1: IF (st_nlayer.ge.1) then 
!
         max_natoms = 0
         max_nscats = 0
         max_n_mole = 0
         max_n_type = 0
         max_n_atom = 0
         DO i = 1, st_ntypes
         enddo
         DO i = 1, st_ntypes
            IF(st_internal(i) ) THEN
               CALL testfile_internal ( st_layer (i ), natoms, nscats, n_mole, &
                                n_type, n_atom)
            ELSE
               CALL test_file ( st_layer (i ), natoms, nscats, n_mole, &
                                n_type, n_atom, -1*i, .true.)
            ENDIF
            max_natoms = MAX(max_natoms, natoms)
            max_nscats = MAX(max_nscats, nscats)
            max_n_mole = MAX(max_n_mole, n_mole)
            max_n_type = MAX(max_n_type, n_type)
            max_n_atom = MAX(max_n_atom, n_atom)
         ENDDO
         natoms = st_nlayer * max_natoms
         nscats = max_nscats
         need_alloc = .false.
            IF(natoms >= NMAX ) THEN
               natoms = MAX(NINT(natoms*1.05),natoms+10,NMAX)
               need_alloc = .true.
            ENDIF
            IF(nscats >= MAXSCAT) THEN
               nscats = MAX(NINT(nscats*1.05),nscats+ 5,MAXSCAT)
               need_alloc = .true.
            ENDIF
            IF( need_alloc) THEN
               CALL alloc_crystal (nscats, natoms)
               IF ( ier_num /= 0 ) RETURN
            ENDIF
!
!        IF more molecules have been read than were allocated
!
               IF(max_n_mole>MOLE_MAX_MOLE .or. max_n_type>MOLE_MAX_TYPE .or.   &
                  max_n_atom>MOLE_MAX_ATOM                          ) THEN
                  n_mole = MAX(max_n_mole +20 ,MOLE_MAX_MOLE)
                  n_type = MAX(max_n_type +10 ,MOLE_MAX_TYPE)
                  n_atom = MAX(max_n_atom +200,MOLE_MAX_ATOM)
                  CALL alloc_molecule(1, 1,n_mole,n_type,n_atom)
                  IF ( ier_num /= 0 ) RETURN
               ENDIF
!
!                                                                       
!     --The first layer is read by the full readstru routine to         
!        ensure correct lattice parameters                              
!                                                                       
         i = 1 
         lout = .false. 
!                                                                       
!     ----read corresponding layer                                      
!                                                                       
         CALL rese_cr 
!
         IF(st_internal(st_type(i)) ) THEN
            CALL readstru_internal (st_layer (st_type (i) ))!, &
!                 NMAX, MAXSCAT, MOLE_MAX_MOLE, &
!                 MOLE_MAX_TYPE, MOLE_MAX_ATOM )
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
         ELSE 
            CALL readstru (NMAX, MAXSCAT, st_layer (st_type (i) ), cr_name,&
            cr_spcgr, cr_a0, cr_win, cr_natoms, cr_nscat, cr_dw, cr_at_lis,&
            cr_pos, cr_mole, cr_iscat, cr_prop, cr_dim, as_natoms, as_at_lis, as_dw,&
            as_pos, as_iscat, as_prop, sav_ncell, sav_r_ncell, sav_ncatoms,&
            spcgr_ianz, spcgr_para)                                        
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
         ENDIF 
!                                                                       
         CALL setup_lattice (cr_a0, cr_ar, cr_eps, cr_gten, cr_reps,    &
         cr_rten, cr_win, cr_wrez, cr_v, cr_vr, lout, cr_gmat, cr_fmat, &
         cr_cartesian)                                                  
!                                                                       
         iatom = 1 
!                                                                       
!     --Add origin of layer                                             
!                                                                       
         DO j = iatom, cr_natoms 
         DO k = 1, 3 
         cr_pos (k, j) = cr_pos (k, j) + st_origin (k, i) 
         ENDDO 
         ENDDO 
!                                                                       
!     --Initialise rotation disorder                                    
!                                                                       
         IF (st_rot_status) then 
!                                                                       
!     ----Rotate layer around normal                                    
!                                                                       
            IF (st_rot_si_no.gt.0.0) then 
               CALL stack_rot_setup ('normal', 6, iatom, i, .true.) 
               CALL symm_setup 
               CALL symm_op_single 
            ENDIF 
!                                                                       
!     ----Rotate layer around mod1                                      
!                                                                       
            IF (st_rot_si_m1.gt.0.0) then 
               CALL stack_rot_setup ('mod1', 4, iatom, i, .true.) 
               CALL symm_setup 
               CALL symm_op_single 
            ENDIF 
!                                                                       
!     ----Rotate layer around mod2                                      
!                                                                       
            IF (st_rot_si_m2.gt.0.0) then 
               CALL stack_rot_setup ('mod2', 4, iatom, i, .true.) 
               CALL symm_setup 
               CALL symm_op_single 
            ENDIF 
         ENDIF 
!                                                                       
!     --Loop over all layers in crystal                                 
!                                                                       
layers:  DO i = 2, st_nlayer 
         iatom = cr_natoms + 1 
         lread = .true. 
internal: IF(st_internal(st_type(i)) ) THEN
!                                                                       
!     ------Read header of internal structure file                               
!                                                                       
            gen_add_n = 0 
            sym_add_n = 0 
            CALL stru_readheader_internal (st_layer(st_type(i)), MAXSCAT, cr_name,   &
            cr_spcgr, cr_at_lis, cr_nscat, cr_dw, cr_a0, cr_win,        &
            sav_ncell, sav_r_ncell, sav_ncatoms, spcgr_ianz, spcgr_para, &
            GEN_ADD_MAX, gen_add_n, gen_add_power, gen_add,                 &
            SYM_ADD_MAX, sym_add_n, sym_add_power, sym_add )
!
!                                                                       
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
!                                                                       
!     ------Now only the atoms are read from an internal file
!                                                                       
            ier_num = 0 
            ier_typ = ER_NONE 
!                                                                       
            CALL struc_read_atoms_internal (st_layer(st_type(i)),NMAX,&
                  cr_natoms, cr_pos, cr_iscat, cr_prop)
            CLOSE (ist) 
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
         ELSE internal
!                                                                       
!     --Open file, read header                                          
!                                                                       
         CALL oeffne (ist, st_layer (st_type (i) ) , 'old') 
         IF (ier_num /= 0) then 
            RETURN
         ENDIF
!                                                                       
!     ------Read header of structure file                               
!                                                                       
            gen_add_n = 0 
            sym_add_n = 0 
            CALL stru_readheader (ist, MAXSCAT, cr_name,   &
            cr_spcgr, cr_at_lis, cr_nscat, cr_dw, cr_a0, cr_win,        &
            sav_ncell, sav_r_ncell, sav_ncatoms, spcgr_ianz, spcgr_para)
!                                                                       
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
!                                                                       
!     ------Now only the atoms are read                                 
!                                                                       
            ier_num = 0 
            ier_typ = ER_NONE 
!                                                                       
            CALL struc_read_atoms (NMAX, MAXSCAT, cr_natoms, cr_nscat,  &
            cr_dw, cr_at_lis, cr_pos, cr_iscat, cr_mole, cr_prop, cr_dim,&
            as_natoms, as_at_lis, as_dw, as_pos, as_iscat, as_prop)     
            CLOSE (ist) 
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
         ENDIF internal
!                                                                       
!     ------Add origin of layer                                         
!                                                                       
            DO j = iatom, cr_natoms 
            DO k = 1, 3 
            cr_pos (k, j) = cr_pos (k, j) + st_origin (k, i) 
            ENDDO 
            ENDDO 
!                                                                       
!     ------Initialise rotation disorder                                
!                                                                       
            IF (st_rot_status) then 
!                                                                       
!     --------Rotate layer around normal                                
!                                                                       
               IF (st_rot_si_no.gt.0.0) then 
                  CALL stack_rot_setup ('normal', 6, iatom, i, .true.) 
                  CALL symm_setup 
                  CALL symm_op_single 
               ENDIF 
!                                                                       
!     --------Rotate layer around mod1                                  
!                                                                       
               IF (st_rot_si_m1.gt.0.0) then 
                  CALL stack_rot_setup ('mod1', 4, iatom, i, .true.) 
                  CALL symm_setup 
                  CALL symm_op_single 
               ENDIF 
!                                                                       
!     --------Rotate layer around mod2                                  
!                                                                       
               IF (st_rot_si_m2.gt.0.0) then 
                  CALL stack_rot_setup ('mod2', 4, iatom, i, .true.) 
                  CALL symm_setup 
                  CALL symm_op_single 
               ENDIF 
            ENDIF 
         ENDDO layers
      ENDIF more1
!                                                                       
      CLOSE (ist) 
!                                                                       
!     Update Crystal dimension                                          
!                                                                       
      CALL update_cr_dim 
!                                                                       
      END SUBROUTINE do_stack_fill                  
!*****7*****************************************************************
      SUBROUTINE stack_rot_setup (line, lbef, iatom, i, lorigin) 
!-                                                                      
!     performs the setup for the rotations around the different axes    
!+                                                                      
      USE discus_allocate_appl_mod
      USE discus_config_mod 
      USE crystal_mod 
      USE symm_mod 
      USE stack_mod 
      USE trafo_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      CHARACTER ( * ) line 
      INTEGER lbef 
      INTEGER iatom 
      INTEGER         :: new_nscat   ! Dummy for allocation 
      INTEGER i, j 
      LOGICAL lorigin 
      LOGICAL str_comp 
!
      IF ( SYM_MAXSCAT < MAXSCAT ) THEN
         new_nscat = MAXSCAT
         CALL alloc_symmetry(new_nscat)
      ENDIF
!                                                                       
      sym_power = 1 
      DO j = 1, cr_nscat 
      sym_latom (j) = .true. 
      ENDDO 
      sym_mode = .false. 
      sym_power_mult = .false. 
      sym_type = .true. 
      sym_trans (1) = 0.0 
      sym_trans (2) = 0.0 
      sym_trans (3) = 0.0 
!                                                                       
      sym_start = iatom 
      sym_end = cr_natoms 
      IF (lorigin) then 
         sym_orig (1) = st_origin (1, i) 
         sym_orig (2) = st_origin (2, i) 
         sym_orig (3) = st_origin (3, i) 
      ELSE 
         sym_orig (1) = 0.0 
         sym_orig (2) = 0.0 
         sym_orig (3) = 0.0 
      ENDIF 
!                                                                       
      IF (str_comp (line, 'normal', 1, lbef, 6) ) then 
         sym_angle = st_rot_ang_no (i) 
         IF (st_rot_no_lspace) then 
            sym_uvw (1) = st_rot_no (1) 
            sym_uvw (2) = st_rot_no (2) 
            sym_uvw (3) = st_rot_no (3) 
            CALL trans (sym_uvw, cr_gten, sym_hkl, 3) 
         ELSE 
            sym_hkl (1) = st_rot_no (1) 
            sym_hkl (2) = st_rot_no (2) 
            sym_hkl (3) = st_rot_no (3) 
            CALL trans (sym_hkl, cr_rten, sym_uvw, 3) 
         ENDIF 
      ELSEIF (str_comp (line, 'mod1', 4, lbef, 4) ) then 
         sym_angle = st_rot_ang_m1 (i) 
         IF (st_rot_m1_lspace) then 
            sym_uvw (1) = st_rot_m1 (1) 
            sym_uvw (2) = st_rot_m1 (2) 
            sym_uvw (3) = st_rot_m1 (3) 
            CALL trans (sym_uvw, cr_gten, sym_hkl, 3) 
         ELSE 
            sym_hkl (1) = st_rot_m1 (1) 
            sym_hkl (2) = st_rot_m1 (2) 
            sym_hkl (3) = st_rot_m1 (3) 
            CALL trans (sym_hkl, cr_rten, sym_uvw, 3) 
         ENDIF 
      ELSEIF (str_comp (line, 'mod2', 4, lbef, 4) ) then 
         sym_angle = st_rot_ang_m2 (i) 
         IF (st_rot_m2_lspace) then 
            sym_uvw (1) = st_rot_m2 (1) 
            sym_uvw (2) = st_rot_m2 (2) 
            sym_uvw (3) = st_rot_m2 (3) 
            CALL trans (sym_uvw, cr_gten, sym_hkl, 3) 
         ELSE 
            sym_hkl (1) = st_rot_m2 (1) 
            sym_hkl (2) = st_rot_m2 (2) 
            sym_hkl (3) = st_rot_m2 (3) 
            CALL trans (sym_hkl, cr_rten, sym_uvw, 3) 
         ENDIF 
      ENDIF 
!                                                                       
      END SUBROUTINE stack_rot_setup                
!*****7*****************************************************************
      SUBROUTINE st_fourier (calc_f2aver)
!-                                                                      
!     Calculates the Fourier transform of the stacking fault decorated  
!     by the respective layers.                                         
!+                                                                      
      USE discus_config_mod 
      USE discus_allocate_appl_mod
      USE crystal_mod 
      USE diffuse_mod 
      USE fourier_sup 
      USE four_strucf_mod
      USE molecule_mod 
      USE discus_save_mod 
      USE stack_mod 
      USE powder_mod 
      USE powder_tables_mod 
      USE structur
      USE spcgr_apply
      USE errlist_mod 
      USE prompt_mod 
      USE wink_mod 
      IMPLICIT none 
!                                                                       
      LOGICAL, INTENT(IN) :: calc_f2aver
!                                                                       
      INTEGER i, j, l 
      INTEGER iscat 
      INTEGER lbeg (3)
      INTEGER ncell 
      INTEGER         :: n_layers ! Number of layers for current layer type
      INTEGER         :: n_qxy    ! Number of data points in reciprocal space
      INTEGER         :: n_nscat  ! Number of different atom types
      INTEGER         :: n_atoms  ! Number of atoms
      INTEGER         :: n_mole  ! number of molecules in input file
      INTEGER         :: n_type  ! number of molecule types in input file
      INTEGER         :: n_atom  ! number of molecule atoms in input file
      REAL      :: xstart  ! qmin  for sin Theta / lambda calculation
      REAL      :: xdelta  ! qstep for sin Theta / lambda calculation
!
      LOGICAL lout 
!                                                                       
      REAL ss 
      REAL seknds 
!
      n_qxy   = 1
      n_nscat = 1
!                                                                       
!------ preset some values                                              
!                                                                       
      ilots = LOT_OFF 
      ncell = 0 
      lout = .false. 
      ss = seknds (0.0) 
      CALL four_layer 
!
!                 ALLOCATE the fourier part of stacking faults
!
      IF ( num (1) * num (2) .gt. ST_MAXQXY .or.            &
           num (1) * num (2) .gt.    MAXQXY      ) THEN
         i = num(1)*num(2)
         CALL alloc_stack (ST_MAXTYPE, ST_MAXLAYER, i     , st_rot_status )
      ENDIF
!                                                                       
!     Now start calculation if sufficient space                         
!                                                                       
      IF (num (1) * num (2) .gt. MAXQXY  .OR.          &
          cr_nscat>DIF_MAXSCAT           .OR.          &
          MAX(cr_natoms,st_nlayer)>DIF_MAXAT   ) THEN
        n_qxy   = MAX(n_qxy,num(1) * num(2),MAXQXY)
        n_nscat = MAX(n_nscat,cr_nscat,DIF_MAXSCAT)
        n_atoms = MAX(cr_natoms,st_nlayer,DIF_MAXSCAT)
        call alloc_diffuse (n_qxy, cr_nscat, n_atoms)
        IF (ier_num.ne.0) THEN
          RETURN
        ENDIF
      ENDIF
      IF (num (1) * num (2) .le.ST_MAXQXY) then 
!                                                                       
!------ --zero some arrays                                              
!                                                                       
         DO i = 1, num (1) * num (2) 
!           st_csf(i) = cmplx(0.0d0,0.0d0)                                
            csf (i) = cmplx (0.0, 0.0) 
!           acsf(i) = cmplx(0.0d0,0.0d0)                                
         ENDDO 
!                                                                       
!     --read first layer to ensure that the metric tensors are set      
!                                                                       
         CALL rese_cr 
         CALL test_file ( st_layer(1), n_atoms, n_nscat, n_mole, n_type,&
                          n_atom,-1, .true. )
         IF(n_atoms > NMAX .or. n_nscat > MAXSCAT) THEN
            n_atoms = MAX( n_atoms, NMAX )
            n_nscat = MAX( n_nscat, MAXSCAT)
            CALL alloc_crystal (n_nscat, n_atoms)
            IF ( ier_num /= 0 ) RETURN
         ENDIF
!
!        IF more molecules have been read than were allocated
!
               IF(n_mole>MOLE_MAX_MOLE .or. n_type>MOLE_MAX_TYPE .or.   &
                  n_atom>MOLE_MAX_ATOM                          ) THEN
                  n_mole = MAX(n_mole +20 ,MOLE_MAX_MOLE)
                  n_type = MAX(n_type +10 ,MOLE_MAX_TYPE)
                  n_atom = MAX(n_atom +200,MOLE_MAX_ATOM)
                  CALL alloc_molecule(1, 1,n_mole,n_type,n_atom)
                  IF ( ier_num /= 0 ) RETURN
               ENDIF
!
         CALL readstru (NMAX, MAXSCAT, st_layer (1), cr_name, cr_spcgr, &
         cr_a0, cr_win, cr_natoms, cr_nscat, cr_dw, cr_at_lis, cr_pos,  &
         cr_mole,                                                       &
         cr_iscat, cr_prop, cr_dim, as_natoms, as_at_lis, as_dw, as_pos,&
         as_iscat, as_prop, sav_ncell, sav_r_ncell, sav_ncatoms,        &
         spcgr_ianz, spcgr_para)                                        
         CALL setup_lattice (cr_a0, cr_ar, cr_eps, cr_gten, cr_reps,    &
         cr_rten, cr_win, cr_wrez, cr_v, cr_vr, lout, cr_gmat, cr_fmat, &
         cr_cartesian)                                                  
!                                                                       
!------ --preset some tables, calculate average structure               
!                                                                       
         CALL four_cexpt 
         CALL four_stltab 
         IF (ier_num.ne.0) return 
!                                                                       
!     --Loop over all layer chemistries                                 
!                                                                       
         DO l = 1, st_nchem 
         nxat = 0 
!                                                                       
!     ----loop over all layers                                          
!                                                                       
         DO j = 1, st_nlayer 
!                                                                       
!     ------insert stacking fault distribution as electrons into crystal
!                                                                       
         IF (st_chem (st_type (j) ) .eq.l) then 
            nxat = nxat + 1 
            xat (nxat, 1) = st_origin (1, j) 
            xat (nxat, 2) = st_origin (2, j) 
            xat (nxat, 3) = st_origin (3, j) 
         ENDIF 
         ENDDO 
         IF (nxat.gt.0) then 
!                                                                       
!     ------calculate structure factor of distribution                  
!                                                                       
            CALL four_strucf (0, .false.) 
!                                                                       
!     ------copy structure factor to temporary place                    
!                                                                       
            DO i = 1, num (1) * num (2) 
            acsf (i) = tcsf (i) 
            ENDDO 
            n_layers = nxat
!                                                                       
!     ----read corresponding layer                                      
!                                                                       
            CALL rese_cr 
            CALL test_file ( st_layer_c(l), n_atoms, n_nscat,n_mole,    &
                             n_type, n_atom, -1, .true. )
            IF(n_atoms > NMAX .or. n_nscat > MAXSCAT) THEN
               n_atoms = MAX( n_atoms, NMAX)
               n_nscat = MAX( n_nscat, MAXSCAT)
               CALL alloc_crystal (n_nscat, n_atoms)
               IF ( ier_num /= 0 ) RETURN
            ENDIF
!
!        IF more molecules have been read than were allocated
!
               IF(n_mole>MOLE_MAX_MOLE .or. n_type>MOLE_MAX_TYPE .or.   &
                  n_atom>MOLE_MAX_ATOM                          ) THEN
                  n_mole = MAX(n_mole +20 ,MOLE_MAX_MOLE)
                  n_type = MAX(n_type +10 ,MOLE_MAX_TYPE)
                  n_atom = MAX(n_atom +200,MOLE_MAX_ATOM)
                  CALL alloc_molecule(1, 1,n_mole,n_type,n_atom)
                  IF ( ier_num /= 0 ) RETURN
               ENDIF
!
            CALL readstru (NMAX, MAXSCAT, st_layer_c (l), cr_name,      &
            cr_spcgr, cr_a0, cr_win, cr_natoms, cr_nscat, cr_dw,        &
            cr_at_lis, cr_pos, cr_mole, cr_iscat, cr_prop, cr_dim, as_natoms,    &
            as_at_lis, as_dw, as_pos, as_iscat, as_prop, sav_ncell,     &
            sav_r_ncell, sav_ncatoms, spcgr_ianz, spcgr_para)           
            IF (ier_num.ne.0) then 
               ier_msg (1) = 'Reading layer file: ' 
               ier_msg (2) = trim(st_layer_c (l))
               RETURN 
            ENDIF 
            CALL setup_lattice (cr_a0, cr_ar, cr_eps, cr_gten, cr_reps, &
            cr_rten, cr_win, cr_wrez, cr_v, cr_vr, lout, cr_gmat,       &
            cr_fmat, cr_cartesian)                                      
            ier_num = 0 
            ier_typ = ER_NONE 
!
!           If necessary allocate diffuse
!
            IF ( cr_nscat>DIF_MAXSCAT           .OR.          &
                 MAX(cr_natoms,st_nlayer)>DIF_MAXAT   ) THEN
               n_qxy   = MAX(n_qxy,num(1) * num(2),MAXQXY)
               n_nscat = MAX(n_nscat,cr_nscat,DIF_MAXSCAT)
               n_atoms = MAX(cr_natoms,st_nlayer,DIF_MAXSCAT)
               call alloc_diffuse (n_qxy, cr_nscat, n_atoms)
               IF (ier_num.ne.0) THEN
                 RETURN
               ENDIF
            ENDIF
            CALL dlink (ano, lambda, rlambda, diff_radiation, &
                        diff_power) 
                                                                        
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
!                                                                       
!------ ------Load form factor tables                                   
!                                                                       
            CALL four_formtab 
!                                                                       
!------ ------zero some arrays                                          
!                                                                       
            DO i = 1, num (1) * num (2) 
            st_csf (i) = cmplx (0.0, 0.0) 
            ENDDO 
!                                                                       
!------ ------loop over all different atom types                        
!                                                                       
            DO iscat = 1, cr_nscat 
            CALL four_getatm (iscat, ilots, lbeg, ncell) 
            CALL four_strucf (iscat, .true.) 
!                                                                       
!------ --------Add this part of the structur factor to the total       
!                                                                       
            DO i = 1, num (1) * num (2) 
            st_csf (i) = st_csf (i) + tcsf (i) 
            ENDDO 
            IF (four_log) then 
               WRITE (output_io, 3000) cr_at_lis (iscat), nxat 
            ENDIF 
!
!------ -------Calculate from factor squared
!
               IF(calc_f2aver) THEN
! Really only needed for <f^2> and <f>^2 for F(Q) and S(Q)
                  xstart = pow_qmin  /zpi
                  xdelta = pow_deltaq/zpi
                  CALL powder_stltab(pow_npkt,xstart,xdelta)  
                  DO i = 1, pow_npkt         ! Always over all points in powder pattern!
                     pow_f2aver (i) = pow_f2aver (i)  + &
                                real (       cfact_pure(powder_istl(i), iscat)  * &
                                      conjg (cfact_pure(powder_istl(i), iscat)))  &
                                * n_layers * nxat
                     pow_faver2 (i) = pow_faver2 (i) +  &
                           SQRT(real (       cfact_pure(powder_istl(i), iscat)  * &
                                      conjg (cfact_pure(powder_istl(i), iscat)))) &
                                * n_layers * nxat
                  ENDDO
                  pow_nreal = pow_nreal + n_layers * nxat
                  pow_faver2(:) = pow_faver2(:)**2
                  pow_u2aver    = pow_u2aver + cr_dw(iscat) * n_layers * nxat
               ENDIF
            ENDDO 
!                                                                       
!     ------Add product of acsf und st_csf to csf                       
!                                                                       
            DO i = 1, num (1) * num (2) 
            csf (i) = csf (i) + st_csf (i) * acsf (i) 
            ENDDO 
         ENDIF 
         ENDDO 
!                                                                       
!     --Calculate average scattering and subtract                       
!                                                                       
         CALL st_fourier_aver 
         DO i = 1, num (1) * num (2) 
            csf (i) = csf (i) - acsf (i) 
         ENDDO 
      ELSE 
         ier_num = - 8 
         ier_typ = ER_APPL 
      ENDIF 
!                                                                       
!     Compute intensity                                                 
!                                                                       
      DO i = 1, num (1) * num (2) 
         dsi (i) = real (csf (i) * conjg (csf (i) ) ) 
      ENDDO 
!                                                                       
      IF (four_log) then 
         CALL four_qinfo 
         ss = seknds (ss) 
         WRITE (output_io, 4000) ss 
      ENDIF 
!                                                                       
 3000 FORMAT (  '                   Atom typ = ',A4,7X,'(# ',I9,' )') 
 4000 FORMAT      (/,' Elapsed time    : ',G12.6,' sec') 
!                                                                       
      END SUBROUTINE st_fourier                     
!*****7*****************************************************************
      SUBROUTINE st_fourier_aver 
!-                                                                      
!     Calculates the Fourier transform of the average lattice.          
!+                                                                      
      USE discus_config_mod 
      USE discus_allocate_appl_mod
      USE crystal_mod 
      USE diffuse_mod 
      USE fourier_sup 
      USE four_strucf_mod
      USE molecule_mod 
      USE discus_save_mod 
      USE stack_mod 
      USE structur
      USE spcgr_apply
      USE errlist_mod 
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER i, j, l 
      INTEGER iscat 
      INTEGER lbeg (3)
      INTEGER ncell 
      LOGICAL lout 
      INTEGER         :: n_qxy    ! Number of data points in reciprocal space
      INTEGER         :: n_nscat  ! Number of different atom types
      INTEGER         :: n_natoms ! Number of atoms 
      INTEGER         :: n_mole  ! number of molecules in input file
      INTEGER         :: n_type  ! number of molecule types in input file
      INTEGER         :: n_atom  ! number of molecule atoms in input file
!
      REAL u (3) 
!
      n_qxy   = 1
      n_nscat = 1
!                                                                       
!------ preset some values                                              
!                                                                       
      ilots = LOT_OFF 
      ncell = 0 
      lout = .false. 
!                                                                       
!     Now start calculation if sufficient space                         
!                                                                       
      IF (num (1) * num (2) .gt. MAXQXY  .OR.          &
          cr_nscat>DIF_MAXSCAT              ) THEN
        n_qxy   = MAX(n_qxy,num(1) * num(2),MAXQXY)
        n_nscat = MAX(n_nscat,cr_nscat,DIF_MAXSCAT)
        call alloc_diffuse (n_qxy, cr_nscat, cr_natoms)
        IF (ier_num.ne.0) THEN
          RETURN
        ENDIF
      ENDIF
      IF (num (1) * num (2) .le.MAXQXY) then 
!                                                                       
!------ --zero some arrays                                              
!                                                                       
         DO i = 1, num (1) * num (2) 
!         st_csf(i) = cmplx(0.0d0,0.0d0)                                
         acsf (i) = cmplx (0.0, 0.0) 
         ENDDO 
!                                                                       
!     --Calculation is only performed if average is needed              
!                                                                       
         IF (st_aver.ne.0.0) then 
!                                                                       
!     ----read first layer to ensure that the metric tensors are set    
!                                                                       
            CALL rese_cr 
         CALL test_file ( st_layer(1), n_natoms, n_nscat,n_mole, n_type,&
                          n_atom, -1, .true. )
         IF(n_natoms > NMAX .or. n_nscat > MAXSCAT) THEN
            n_natoms = MAX( n_natoms, NMAX)
            n_nscat  = MAX( n_nscat, MAXSCAT)
            CALL alloc_crystal (n_nscat, n_natoms)
            IF ( ier_num /= 0 ) RETURN
         ENDIF
!
!        IF more molecules have been read than were allocated
!
               IF(n_mole>MOLE_MAX_MOLE .or. n_type>MOLE_MAX_TYPE .or.   &
                  n_atom>MOLE_MAX_ATOM                          ) THEN
                  n_mole = MAX(n_mole +20 ,MOLE_MAX_MOLE)
                  n_type = MAX(n_type +10 ,MOLE_MAX_TYPE)
                  n_atom = MAX(n_atom +200,MOLE_MAX_ATOM)
                  CALL alloc_molecule(1, 1,n_mole,n_type,n_atom)
                  IF ( ier_num /= 0 ) RETURN
               ENDIF
!
            CALL readstru (NMAX, MAXSCAT, st_layer (1), cr_name,        &
            cr_spcgr, cr_a0, cr_win, cr_natoms, cr_nscat, cr_dw,        &
            cr_at_lis, cr_pos, cr_mole, cr_iscat, cr_prop, cr_dim, as_natoms,    &
            as_at_lis, as_dw, as_pos, as_iscat, as_prop, sav_ncell,     &
            sav_r_ncell, sav_ncatoms, spcgr_ianz, spcgr_para)           
            CALL setup_lattice (cr_a0, cr_ar, cr_eps, cr_gten, cr_reps, &
            cr_rten, cr_win, cr_wrez, cr_v, cr_vr, lout, cr_gmat,       &
            cr_fmat, cr_cartesian)                                      
!                                                                       
            IF (ier_num.ne.0) return 
!                                                                       
!     ----Calculate average translation                                 
!                                                                       
            u (1) = (st_origin (1, st_nlayer) - st_origin (1, 1) )      &
            / (st_nlayer - 1)                                           
            u (2) = (st_origin (2, st_nlayer) - st_origin (2, 1) )      &
            / (st_nlayer - 1)                                           
            u (3) = (st_origin (3, st_nlayer) - st_origin (3, 1) )      &
            / (st_nlayer - 1)                                           
!                                                                       
!         write (output_io,*) ' <t> : ', u                              
!                                                                       
!     ----Determine average displacement of each layer type off the     
!         average translation direction                                 
!                                                                       
            DO l = 1, st_ntypes 
            st_disp (1, l) = 0.0 
            st_disp (2, l) = 0.0 
            st_disp (3, l) = 0.0 
            st_ndisp (l) = 0 
            ENDDO 
!                                                                       
!     ----loop over all layers, compile displacements                   
!                                                                       
            DO j = 1, st_nlayer 
            DO i = 1, 3 
            st_disp (i, st_type (j) ) = st_disp (i, st_type (j) )       &
            + st_origin (i, j) - st_origin (i, 1) - (j - 1) * u (i)     
            ENDDO 
            st_ndisp (st_type (j) ) = st_ndisp (st_type (j) ) + 1 
            ENDDO 
!                                                                       
!     ----calculate average displacement                                
!                                                                       
            DO l = 1, st_ntypes 
            IF (st_ndisp (l) .gt.0) then 
               DO i = 1, 3 
               st_disp (i, l) = st_disp (i, l) / float (st_ndisp (l) ) 
               ENDDO 
!             write (output_io,5000) (st_disp(i,l),i=1,3)               
            ENDIF 
            ENDDO 
!                                                                       
!     ----loop over all layers                                          
!                                                                       
            nxat = 0 
            DO j = 1, st_nlayer 
!                                                                       
!     ------insert average stacking fault distribution as electrons     
!             into crystal                                              
!                                                                       
            nxat = nxat + 1 
            xat (nxat, 1) = st_origin (1, 1) + (j - 1) * u (1) 
            xat (nxat, 2) = st_origin (2, 1) + (j - 1) * u (2) 
            xat (nxat, 3) = st_origin (3, 1) + (j - 1) * u (3) 
            ENDDO 
!                                                                       
!     ----calculate structure factor of distribution                    
!                                                                       
            CALL four_strucf (0, .false.) 
!                                                                       
!     ----copy structure factor to temporary place                      
!                                                                       
            DO i = 1, num (1) * num (2) 
            acsf (i) = tcsf (i) 
            ENDDO 
!                                                                       
!------ ----zero some arrays                                            
!                                                                       
            DO i = 1, num (1) * num (2) 
            st_csf (i) = cmplx (0.0, 0.0) 
            ENDDO 
!                                                                       
!     ----Loop over all layer types                                     
!                                                                       
            DO l = 1, st_ntypes 
!                                                                       
!     ------read corresponding layer                                    
!                                                                       
            CALL rese_cr 
         CALL test_file ( st_layer(l), n_natoms, n_nscat,n_mole, n_type,&
                          n_atom, -1, .true. )
         IF(n_natoms > NMAX .or. n_nscat > MAXSCAT) THEN
            n_natoms = MAX( n_natoms, NMAX )
            n_nscat  = MAX( n_nscat,  MAXSCAT)
            CALL alloc_crystal (n_nscat, n_natoms)
            IF ( ier_num /= 0 ) RETURN
         ENDIF
!
!        IF more molecules have been read than were allocated
!
               IF(n_mole>MOLE_MAX_MOLE .or. n_type>MOLE_MAX_TYPE .or.   &
                  n_atom>MOLE_MAX_ATOM                          ) THEN
                  n_mole = MAX(n_mole +20 ,MOLE_MAX_MOLE)
                  n_type = MAX(n_type +10 ,MOLE_MAX_TYPE)
                  n_atom = MAX(n_atom +200,MOLE_MAX_ATOM)
                  CALL alloc_molecule(1, 1,n_mole,n_type,n_atom)
                  IF ( ier_num /= 0 ) RETURN
               ENDIF
!
            CALL readstru (NMAX, MAXSCAT, st_layer (l), cr_name,        &
            cr_spcgr, cr_a0, cr_win, cr_natoms, cr_nscat, cr_dw,        &
            cr_at_lis, cr_pos, cr_mole, cr_iscat, cr_prop, cr_dim, as_natoms,    &
            as_at_lis, as_dw, as_pos, as_iscat, as_prop, sav_ncell,     &
            sav_r_ncell, sav_ncatoms, spcgr_ianz, spcgr_para)           
            CALL setup_lattice (cr_a0, cr_ar, cr_eps, cr_gten, cr_reps, &
            cr_rten, cr_win, cr_wrez, cr_v, cr_vr, lout, cr_gmat,       &
            cr_fmat, cr_cartesian)                                      
!                                                                       
!     ------Add average displacement to each atom                       
!                                                                       
            DO i = 1, cr_natoms 
            DO j = 1, 3 
            cr_pos (j, i) = cr_pos (j, i) + st_disp (j, l) 
            ENDDO 
            ENDDO 
!                                                                       
            ier_num = 0 
            ier_typ = ER_NONE 
            CALL dlink (ano, lambda, rlambda, diff_radiation, &
                        diff_power) 
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
!                                                                       
!------ ------Load form factor tables                                   
!                                                                       
            CALL four_formtab 
!                                                                       
!------ ------loop over all different atom types                        
!                                                                       
            DO iscat = 1, cr_nscat 
            CALL four_getatm (iscat, ilots, lbeg, ncell) 
            CALL four_strucf (iscat, .true.) 
!                                                                       
!------ --------Add this part of the structur factor to the total       
!             Wheighted by the relative amount of layers of this type   
!                                                                       
            DO i = 1, num (1) * num (2) 
            st_csf (i) = st_csf (i) + tcsf (i) * float (st_number (l) ) &
            / float (st_nlayer)                                         
            ENDDO 
            WRITE (output_io, 3000) cr_at_lis (iscat), nxat 
            ENDDO 
            ENDDO 
!                                                                       
!     ----Save product of acsf und st_csf to acsf                       
!                                                                       
            DO i = 1, num (1) * num (2) 
            acsf (i) = st_csf (i) * acsf (i) 
            ENDDO 
         ENDIF 
      ELSE 
         ier_num = - 8 
         ier_typ = ER_APPL 
      ENDIF 
!                                                                       
!                                                                       
 3000 FORMAT     ( '                   Atom typ = ',A4,7X,'(# ',I9,' )') 
!                                                                       
      END SUBROUTINE st_fourier_aver                
!*****7*****************************************************************
      SUBROUTINE stack_dist_file ()
!     SUBROUTINE stack_dist_file (cr_spcgr, cr_a0, cr_win, cr_dim,      &
!     st_name, st_spcgr, st_a0, st_win, st_natoms, st_nscat, st_dw,     &
!     st_at_lis, st_pos, st_iscat, st_dim, sa_natoms, sa_at_lis, sa_dw, &
!     sa_pos, sa_iscat)                                                 
!-                                                                      
!     A 'crystal' file is read and the atom types are interpreted as    
!     layer types.                                                      
!     Later on the origins will be interpreted as well as               
!     origins of the stacking faults.                                   
!+                                                                      
      USE discus_config_mod 
      USE stack_mod 
      USE stack_cr_mod
      USE structur
      USE spcgr_apply
      USE errlist_mod 
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER, DIMENSION(3)   :: sav_ncell (3) 
      INTEGER                 :: sav_ncatoms 
      LOGICAL                 :: sav_r_ncell 
!                                                                       
!     Initialize counters for Atom numbers                              
!                                                                       
      st_natoms     = 0 
      st_nscat      = 0 
      sa_natoms     = 0 
      st_at_lis (0) = 'VOID' 
      st_dw (0)     = 0.0 
      st_at_lis     = ' '  ! i=1,ST_MAX_SCAT (i) = ' ' 
!                                                                       
!     Now read the pseudo microdomain structure                         
!
      CALL readstru (ST_MMAX, ST_MAX_SCAT, st_infile, st_name, st_spcgr,&
      st_a0, st_win, st_natoms, st_nscat, st_dw, st_at_lis, st_pos,     &
      st_mole,                                                          &
      st_iscat, st_prop, st_dim, sa_natoms, sa_at_lis, sa_dw, sa_pos,   &
      sa_iscat, sa_prop, sav_ncell, sav_r_ncell, sav_ncatoms,           &
      st_spcgr_ianz, st_spcgr_para)                                     
      IF (ier_num.ne.0) then 
         ier_msg (1) = 'Error occured while reading' 
         ier_msg (2) = trim(st_infile)
         RETURN 
      ENDIF 
!                                                                       
      IF (st_natoms.gt.ST_MAXLAYER) then 
         ier_num = - 5 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
!                                                                       
      END SUBROUTINE stack_dist_file                
END MODULE stack_menu
