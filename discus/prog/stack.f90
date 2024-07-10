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
use fourier_conv_mod
USE stack_mod 
USE stack_rese_mod 
USE spcgr_apply
use structur, ONLY:read_to_internal
use stack_cr_mod
!
USE ber_params_mod
USE calc_expr_mod
USE doact_mod 
USE do_eval_mod
USE do_wait_mod
USE build_name_mod
USE errlist_mod 
USE get_params_mod
use kdo_all_mod
USE learn_mod 
USE lib_do_operating_mod
USE lib_echo
USE lib_errlist_func
USE lib_help
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
INTEGER, PARAMETER :: MIN_PARA = 99  ! A command requires at least these no of parameters
                                           ! Needs work as it should also be >= ST_MAXTYPE
INTEGER :: maxw 
!                                                                       
CHARACTER(LEN=PREC_STRING), DIMENSION(MIN_PARA) :: cpara
INTEGER            , DIMENSION(MIN_PARA) :: lpara
REAL(KIND=PREC_DP) , DIMENSION(MIN_PARA) :: werte
!
CHARACTER(len=12)          :: befehl 
CHARACTER(LEN=LEN(prompt)) :: orig_prompt
CHARACTER(LEN=PREC_STRING) :: line, zeile
CHARACTER(LEN=PREC_STRING) :: strucfile    ! Dummy for input files
INTEGER :: lp, length, lbef 
INTEGER :: indxg, ianz, i, j, k 
LOGICAL :: lend
LOGICAL, SAVE :: linit    = .true. 
!
INTEGER,SAVE  :: n_types  = 0 ! Current number of layer types
INTEGER,SAVE  :: n_layers = 0 ! Current number of layers
INTEGER,dimension(3),SAVE  :: n_qxy    = 0 ! Current number of points in reciprocal space
!                                                                       
!                                                                       
maxw = MIN_PARA
lend = .false. 
CALL no_error 
!
orig_prompt = prompt
prompt = prompt (1:len_str (prompt) ) //'/stack' 
!                                                                       
loop_main: DO while (.not.lend)                 ! Main stack loop
!
   if_error: IF (ier_num.ne.0) THEN 
      CALL errlist 
      IF (ier_sta.ne.ER_S_LIVE) THEN 
         IF (lmakro .OR. lmakro_error) THEN  ! Error within macro or termination errror
            IF(sprompt /= prompt ) THEN
               ier_num = -10
               ier_typ = ER_COMM
               ier_msg(1) = ' Error occured in stack menu'
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
   ENDIF if_error
!
   CALL get_cmd(line, length, befehl, lbef, zeile, lp, prompt) 
   IF (ier_num /= 0) cycle loop_main
!
   IF(line == ' '      .or. line(1:1) == '#' .or. &
      line == char(13) .or. line(1:1) == '!'        ) cycle loop_main
!                                                                       
!     ----search for "="                                                
!                                                                       
   indxg = index (line, '=') 
   IF (indxg.ne.0                                              &
        .AND..NOT. (str_comp (befehl, 'echo',   2, lbef, 4) )    &
        .AND..NOT. (str_comp (befehl, 'system', 2, lbef, 6) )    &
        .AND..NOT. (str_comp (befehl, 'help',   2, lbef, 4) .OR. &
                    str_comp (befehl, '?   ',   2, lbef, 4) )    &
        .AND. INDEX(line,'==') == 0                            ) THEN
!                                                                       
!     ------evaluate an expression and assign the value to a variable   
!                                                                       
            CALL do_math(line, indxg, length) 
      cycle loop_main
   endif
!
!--STACK specific commands
!                                                                       
   IF(linit) then 
      CALL alloc_stack ( st_layer_increment, 1, st_rot_status)
      n_types  = ST_MAXTYPE
      n_layers = ST_MAXLAYER
      n_qxy    = ST_MAXQXY
      linit    = .false.
   ENDIF 
!                                                                       
!     ----Determine if average translation is calculated or set 'aver'  
!                                                                       
   IF (str_comp (befehl, 'average', 1, lbef, 7) ) then 
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
!     ----exit 'exit'
!
   ELSEIF (str_comp (befehl, 'exit', 2, lbef, 4) ) then
      lend = .true.
!
!     ----help 'help','?'
!
   ELSEIF(str_comp (befehl, 'help', 2, lbef, 4) .or.            &
          str_comp (befehl, '?   ', 1, lbef, 4)        ) then                                      
      IF (str_comp (zeile, 'errors', 2, lp, 6) ) then 
         lp = lp + 7 
         CALL do_hel ('discus '//zeile, lp) 
      ELSE 
         lp = lp + 12 
         CALL do_hel ('discus stack '//zeile, lp) 
      ENDIF 
   ELSEIF (str_comp (befehl, 'reset', 2, lbef, 5) ) then 
      CALL do_stack_rese 
      linit = .TRUE.
!                                                                       
!     ----Read a single collumn of the correlation matrix 'ccolumn'     
!                                                                       
   ELSEIF (str_comp (befehl, 'ccolumn', 2, lbef, 7) ) then 
      CALL get_params(zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.eq.0) then 
         IF (st_ntypes.lt.ianz.and.ianz.le.ST_MAXTYPE+1) then                                               
            CALL ber_params(ianz, cpara, lpara, werte, maxw)
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
   ELSEIF (str_comp (befehl, 'celement', 2, lbef, 8) ) then 
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.eq.0) then 
         IF (ianz.eq.3) then 
            CALL ber_params (ianz, cpara, lpara, werte, maxw)
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
!     ----create list of origins 'create'                               
!                                                                       
   ELSEIF (str_comp (befehl, 'create', 3, lbef, 6) ) then 
!
!                 If necessary allocate array sizes
!
      IF ( st_nlayer >  ST_MAXLAYER ) THEN
         CALL alloc_stack (ST_MAXTYPE, st_nlayer, st_rot_status )
      ENDIF
      st_new_form = .TRUE.
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
   ELSEIF (str_comp (befehl, 'distribution', 1, lbef, 12) ) then 
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ianz.ge.1) then 
         IF(str_comp(cpara(1), 'matrix', 1, lpara(1), 6) ) then
            st_distr = ST_DIST_MATRIX 
         ELSEIF(str_comp (cpara (1) , 'file', 1, lpara (1), 4) ) then
            IF (ianz.gt.1) then 
               CALL do_build_name (ianz, cpara, lpara, werte, maxw, 2)
               IF (ier_num.eq.0) then 
                  st_distr = ST_DIST_FILE 
                  st_infile = cpara (2) (1:lpara(2))
                  st_infile_l = lpara (2) 
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
         ELSEIF (str_comp (cpara (1) , 'list', 1, lpara (1), 4) ) then
            IF (ianz.gt.1) then 
               CALL do_build_name (ianz, cpara, lpara, werte, maxw, 2)
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
         IF (st_distr.eq.ST_DIST_FILE.or.st_distr.eq.ST_DIST_LIST) THEN
!
!           Always read distribution file into internal memory
!
            if(st_infile(1:8)/='internal') then
               strucfile = st_infile
               call read_to_internal(st_infile, 'internal_stack_list.')
               st_infile = 'internal_stack_list.' // st_infile(1:len_trim(st_infile))
            endif
         endif
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
!     ----Determine type of first layer 'first' 'random'| <number>      
!                                                                       
   ELSEIF (str_comp (befehl, 'first', 1, lbef, 5) ) then 
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.eq.0) then 
         IF (ianz.eq.1) then 
            IF(str_comp(cpara(1),'random', 3, lpara(1), 6)) THEN
               st_first = 0
            ELSE
               CALL ber_params(ianz,cpara,lpara,werte,maxw)                                           
               IF (ier_num.eq.0) then 
                  IF(nint(werte(1)) > 0) THEN
                     st_first = nint (werte (1) ) 
                  ELSE 
                     ier_num = -139 
                     ier_typ = ER_APPL 
                     ier_msg(1)='Layer type must be >= 1'
                  ENDIF 
               ENDIF 
            ENDIF 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
      ENDIF 
!                                                                       
!      ---Calculate the Fourier transform of the decorated              
!           s.f. 'fourier'                                              
!                                                                       
   ELSEIF (str_comp (befehl, 'fourier', 1, lbef, 7) ) then 
      CALL four_resolution(zeile, lp)
      if(ier_num==0) then
         four_log = .true. 
         CALL st_fourier (.false.)
         if(ier_num==0) then
            four_was_run = .true.
            call four_conv           ! Convolute diffraction pattern
         endif
      endif
!                                                                       
!     ----read the name of a new layer type                'layer'      
!                                                                       
   ELSEIF (str_comp (befehl, 'layer', 1, lbef, 5) ) then 
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.eq.0) then 
         CALL do_build_name(ianz, cpara, lpara, werte, maxw, 1)
         IF (ier_num.eq.0) then 
            IF (ianz.eq.1) then 
!
!   If necessary increase array sizes
!
               IF ( st_ntypes >=  ST_MAXTYPE ) THEN
                  i       = ST_MAXTYPE + st_layer_increment
                  CALL alloc_stack ( i, ST_MAXLAYER, st_rot_status )
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
 5000             CONTINUE 
                  i = st_ntypes
                  strucfile = st_layer(i)
                  call read_to_internal(strucfile, 'internal_stack_layer.')
                  st_layer(i) = 'internal_stack_layer.'//st_layer(i)(1:len_trim(st_layer(i)))
                  st_llayer(i) = st_llayer(i) + 21
                  st_internal(i) = .true.
!
                  i = st_nchem
                  strucfile = st_layer_c(i)
                  call read_to_internal(strucfile, 'internal_stack_layer.')
                  st_layer_c(i) = 'internal_stack_layer.'//st_layer_c(i)(1:len_trim(st_layer(i)))
!
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
   ELSEIF (str_comp (befehl, 'modulus', 1, lbef, 7) ) then 
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.eq.0) then 
         IF (ianz.eq.6) then 
            CALL ber_params (ianz, cpara, lpara, werte, maxw)
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
   ELSEIF (str_comp (befehl, 'number', 1, lbef, 6) ) then 
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.eq.0) then 
         IF (ianz.eq.1) then 
            CALL ber_params (ianz, cpara, lpara, werte, maxw)
            IF (ier_num.eq.0) then 
               i = nint (werte (1) ) 
               st_nlayer = i 
            ENDIF 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
      ENDIF 
!                                                                       
!     ----Select values for random stacking faults 'random'             
!                                                                       
   ELSEIF (str_comp (befehl, 'random', 2, lbef, 6) ) then 
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.eq.0) then 
         IF(str_comp(cpara(1), 'probability', 1, lpara(1), 11)) then
!                                                                       
!     ----------Random stacking fault probability                       
!                                                                       
            IF (ianz.eq.2) then 
               cpara (1) = '0' 
               lpara (1) = 1 
               CALL ber_params (ianz, cpara, lpara, werte, maxw)
               IF (ier_num.eq.0) then 
                  st_prob = werte (2) 
               ENDIF 
            ENDIF 
         ELSEIF (str_comp (cpara (1) , 'offset', 1, lpara ( 1) , 6) ) then
!                                                                       
!     ----------Random stacking fault offset                            
!                                                                       
            IF (ianz.eq.4) then 
               cpara (1) = '0' 
               lpara (1) = 1 
               CALL ber_params (ianz, cpara, lpara, werte, maxw)
               IF (ier_num.eq.0) then 
                  DO i = 1, 3 
                     st_off (i) = werte (i + 1) 
                  ENDDO 
               ENDIF 
            ENDIF 
         ELSEIF (str_comp (cpara (1) , 'sigma', 1, lpara (1) , 5) ) then
!                                                                       
!     ----------Random stacking fault sigma for offset                  
!                                                                       
            IF (ianz.eq.4) then 
               cpara (1) = '0' 
               lpara (1) = 1 
               CALL ber_params (ianz, cpara, lpara, werte, maxw)
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
   ELSEIF (str_comp (befehl, 'rotate', 2, lbef, 7) ) then 
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.eq.0) then 
!
!        Allocate rotational stacking faults
!
         CALL alloc_stack ( ST_MAXTYPE, ST_MAXLAYER, .true. )
!
         IF(str_comp(cpara(1), 'mode', 1, lpara(1), 4)) then
!                                                                       
!     ----------Mode of rotations "all" or only in case of              
!                  "stacking fault"                                     
!                                                                       
            IF (ianz.eq.2) then 
               IF(str_comp(cpara(2), 'all', 1, lpara(2) , 3) ) then                                  
                  st_rot_mode = .true. 
               ELSEIF (str_comp (cpara (2) , 'fault', 1, lpara (2) , 5) ) then                        
                  st_rot_mode = .false. 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
         ELSEIF (str_comp (cpara (1) , 'axis', 1, lpara (1) , 4) ) then                                        
            cpara (1) = '0' 
            lpara (1) = 1 
!                                                                       
!     ----------Direction of in-plane vector 1                          
!                                                                       
            IF (str_comp (cpara (2) , 'mod1', 4, lpara (2) , 4) ) then                                       
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
                  CALL ber_params (ianz, cpara, lpara, werte, maxw)                              
                  IF (ier_num.eq.0) then 
                     st_rot_m1 (1) = werte (3) 
                     st_rot_m1 (2) = werte (4) 
                     st_rot_m1 (3) = werte (5) 
                  ENDIF 
                  IF (str_comp (cpara (6) , 'd', 1, lpara ( 6) , 1) ) then                            
                     st_rot_m1_lspace = .true. 
                  ELSEIF (str_comp (cpara (6) , 'r', 1, lpara (6) , 1) ) then                     
                     st_rot_m1_lspace = .false. 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
            ELSEIF (str_comp (cpara (2) , 'mod2', 4, lpara ( 2) , 4) ) then                                  
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
                  CALL ber_params (ianz, cpara, lpara, werte, maxw)                              
                  IF (ier_num.eq.0) then 
                     st_rot_m2 (1) = werte (3) 
                     st_rot_m2 (2) = werte (4) 
                     st_rot_m2 (3) = werte (5) 
                  ENDIF 
                  IF (str_comp (cpara (6) , 'd', 1, lpara ( 6) , 1) ) then                            
                     st_rot_m2_lspace = .true. 
                  ELSEIF (str_comp (cpara (6) , 'r', 1, lpara (6) , 1) ) then                     
                     st_rot_m2_lspace = .false. 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
            ELSEIF (str_comp (cpara (2) , 'normal', 1, lpara (2) , 6) ) then                           
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
                  CALL ber_params (ianz, cpara, lpara, werte, maxw)                              
                  IF (ier_num.eq.0) then 
                     st_rot_no (1) = werte (3) 
                     st_rot_no (2) = werte (4) 
                     st_rot_no (3) = werte (5) 
                  ENDIF 
                  IF (str_comp (cpara (6) , 'd', 1, lpara ( 6) , 1) ) then                            
                     st_rot_no_lspace = .true. 
                  ELSEIF (str_comp (cpara (6) , 'r', 1, lpara (6) , 1) ) then                     
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
         ELSEIF (str_comp (cpara (1) , 'sigma', 1, lpara (1) , 5) ) then                                        
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
                  CALL ber_params (ianz, cpara, lpara, werte, maxw)                              
                  IF (ier_num.eq.0) then 
                     st_rot_si_no = werte (3) 
                  ENDIF 
               ELSEIF (str_comp (cpara (2) , 'mod1', 4, lpara (2) , 4) ) then                        
                  cpara (2) = '0' 
                  lpara (2) = 1 
                  CALL ber_params (ianz, cpara, lpara, werte, maxw)                              
                  IF (ier_num.eq.0) then 
                     st_rot_si_m1 = werte (3) 
                  ENDIF 
               ELSEIF (str_comp (cpara (2) , 'mod2', 4, lpara (2) , 4) ) then                        
                  cpara (2) = '0' 
                  lpara (2) = 1 
                  CALL ber_params (ianz, cpara, lpara, werte, maxw)                              
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
         ELSEIF (str_comp (cpara (1) , 'status', 2, lpara ( 1) , 6) ) then                                     
!                                                                       
!     ----------Switch status of rotational disorder "on" or "off"      
!                                                                       
            IF (ianz.eq.2) then 
               IF (str_comp (cpara (2) , 'on', 2, lpara (2) , 2) ) then                                  
                  st_rot_status = .true. 
               ELSEIF (str_comp (cpara (2) , 'off', 2, lpara (2) , 3) ) then                        
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
      call stack_set(zeile, lp)
!                                                                       
!     ----show current parameters 'show'                                
!                                                                       
   ELSEIF (str_comp (befehl, 'show', 2, lbef, 4) ) then 
      call stack_show
!                                                                       
!     ----Select sigmas of translational vector 'sigma'                 
!                                                                       
   ELSEIF (str_comp (befehl, 'sigma', 2, lbef, 5) ) then 
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.eq.0) then 
         IF (ianz.eq.5) then 
            CALL ber_params (ianz, cpara, lpara, werte, maxw)                                           
            IF (ier_num.eq.0) then 
               i = nint (werte (1) ) 
               j = nint (werte (2) ) 
               IF (0.lt.i.and.i.le.ST_MAXTYPE.and.                              &
                   0.lt.j.and.j.le.ST_MAXTYPE      ) then
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
   ELSEIF (str_comp (befehl, 'translation', 1, lbef, 11) ) then 
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.eq.0) then 
         IF (ianz.eq.5) then 
            CALL ber_params (ianz, cpara, lpara, werte, maxw)                                           
            IF (ier_num.eq.0) then 
               i = nint (werte (1) ) 
               j = nint (werte (2) ) 
               IF(0.lt.i.and.i.le.ST_MAXTYPE.and.                              &
                  0.lt.j.and.j.le.ST_MAXTYPE        ) then
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
!     ----Unknown stack command, try general commands
!                                                                       
   ELSE 
      call kdo_all(befehl, lbef, zeile, lp)
   ENDIF 
!
ENDDO loop_main
!
prompt = orig_prompt
!                                                                       
END SUBROUTINE stack                          
!
!*****7*****************************************************************
!
subroutine stack_set(zeile, lp)
!-
!   Perform all 'set' commands for stack
!+
use stack_mod
!
use ber_params_mod
use errlist_mod
use get_params_mod
use precision_mod
use str_comp_mod
use take_param_mod
!
implicit none
!
character(len=*), intent(inout) :: zeile
integer         , intent(inout) :: lp
!
integer, parameter :: MAXW = 6
!                                                                       
CHARACTER(LEN=PREC_STRING), DIMENSION(MAXW) :: cpara
INTEGER                   , DIMENSION(MAXW) :: lpara
REAL(KIND=PREC_DP)        , DIMENSION(MAXW) :: werte
integer                                     :: ianz
!
!
integer, parameter :: NOPTIONAL = 4
integer, parameter :: O_AVER    = 1
integer, parameter :: O_MOD     = 2
integer, parameter :: O_TRANS   = 3
integer, parameter :: O_ATOM    = 4
character(LEN=   7), dimension(NOPTIONAL) :: oname   !Optional parameter names
character(LEN=PREC_STRING), dimension(NOPTIONAL) :: opara   !Optional parameter strings returned
integer            , dimension(NOPTIONAL) :: loname  !Lenght opt. para name
integer            , dimension(NOPTIONAL) :: lopara  !Lenght opt. para name returned
logical            , dimension(NOPTIONAL) :: lpresent!opt. para is present
real(kind=PREC_DP) , dimension(NOPTIONAL) :: owerte   ! Calculated values
integer, parameter                        :: ncalc = 1 ! Number of values to calculate 
!
data oname  / 'aver', 'modulus',  'trans',  'atoms'   /
data loname /  2    ,  7       ,   4     ,   5       /
opara  =  (/ '0.0000', 'on    ', 'fixed ', 'off   ' /)   ! Always provide fresh default values
lopara =  (/  6,        2,        5      ,  3       /)
owerte =  (/  0.0D0,    1.0D0,    1.0D0  ,  real(ST_ATOM_OFF,KIND=PREC_DP) /)
!
!
CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
IF (ier_num /= 0) then 
      ier_num = - 6 
      ier_typ = ER_COMM 
   return
ENDIF 
!
call get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
IF (ier_num /= 0) then 
      ier_num = - 6 
      ier_typ = ER_COMM 
   return
ENDIF 
!
         IF (str_comp (cpara (1) , 'aver', 1, lpara (1) , 4)) then                                             
            IF (ianz.eq.2) then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               ianz = 1 
               CALL ber_params (ianz, cpara, lpara, werte, maxw)                                        
               IF (ier_num.eq.0) then 
                  st_aver = werte (1) 
               ENDIF 
            ELSEIF (ianz.eq.1) then 
               st_aver = 0.0 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
         ELSEIF (str_comp (cpara (1) , 'modulus', 1, lpara (1) , 7) ) then                                        
            IF (ianz.eq.2) then 
               IF (str_comp (cpara (2) , 'off', 2, lpara (2) , 3) ) then                                  
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
         ELSEIF (str_comp (cpara (1) , 'translation', 1, lpara (1) , 11) ) then                                        
            IF (ianz.eq.2) then 
               IF(str_comp(cpara(2), 'average', 1, lpara(2), 7) ) then
                  st_tra_aver = .true. 
               ELSEIF (str_comp (cpara (2) , 'fixed', 1, lpara (2) , 5) ) then                        
                  st_tra_aver = .false. 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
else     ! optional parameter style was used
   if(lpresent(O_AVER)) then
      st_aver = owerte(O_AVER)
   endif
   if(lpresent(O_MOD )) then
      if(str_comp(opara(O_MOD), 'off', 3, lopara(O_MOD), 3)) then
         st_mod_sta = .false.
      elseif(str_comp(opara(O_MOD), 'on', 2, lopara(O_MOD), 2)) then
         st_mod_sta = .true.
      else
         ier_num = - 6
         ier_typ = ER_COMM
      endif
   endif
   if(lpresent(O_TRANS )) then
      if(str_comp(opara(O_TRANS), 'average', 1, lopara(O_TRANS), 7)) then
         st_tra_aver = .true.
      elseif(str_comp(opara(O_TRANS), 'fixed', 1, lopara(O_TRANS), 5)) then
         st_tra_aver = .false.
      else
         ier_num = - 6
         ier_typ = ER_COMM
      endif
   endif
ENDIF 
st_mod_atom = nint(owerte(O_ATOM))
if(lpresent(O_ATOM)) then
   if(str_comp(opara(O_ATOM), 'off', 3, lopara(O_ATOM), 3)) then
      st_mod_atom = ST_ATOM_OFF
   elseif(str_comp(opara(O_ATOM), 'strict', 3, lopara(O_ATOM), 6)) then
      st_mod_atom = ST_ATOM_STRICT
   elseif(str_comp(opara(O_ATOM), 'on', 2, lopara(O_ATOM), 2)) then
      st_mod_atom = ST_ATOM_ON 
   else
      ier_num = - 6
      ier_typ = ER_COMM
   endif
endif
!
end subroutine stack_set
!
!*****7*****************************************************************
!
subroutine stack_show
!-
!  Show the current settings
!+
use stack_mod
!
use prompt_mod
!
implicit none
!
integer :: i, j, k
!
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
   WRITE (output_io, 3020) i, st_layer (i) (1:st_llayer(i) )
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
   WRITE(output_io, 3060) (st_corr(i, j), j = 1, st_ntypes)
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
   WRITE (output_io, 3140) st_rot_no, 'd', st_rot_si_no                                       
ELSE 
   WRITE (output_io, 3140) st_rot_no, 'r', st_rot_si_no                                       
ENDIF 
IF (st_rot_m1_lspace) then 
   WRITE (output_io, 3140) st_rot_m1, 'd', st_rot_si_m1                                       
ELSE 
   WRITE (output_io, 3140) st_rot_m1, 'r', st_rot_si_m1                                       
ENDIF 
IF (st_rot_m2_lspace) then 
   WRITE (output_io, 3140) st_rot_m2, 'd', st_rot_si_m2                                       
ELSE 
   WRITE (output_io, 3140) st_rot_m2, 'r', st_rot_si_m2                                       
ENDIF 
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
end subroutine stack_show
!
!*****7*****************************************************************
!
SUBROUTINE do_stack_create 
!-                                                                      
!     Creates the list of origins for the crystal with stacking faults  
!+                                                                      
USE discus_config_mod 
USE discus_allocate_appl_mod
USE class_internal
USE crystal_mod 
USE molecule_mod 
USE read_internal_mod
USE save_menu, ONLY:save_internal
USE stack_mod 
USE stack_cr_mod 
USE structur
USE spcgr_apply
USE symm_sup_mod
!
USE errlist_mod 
USE lib_errlist_func
use matrix_mod
USE lib_random_func
USE random_mod
USE param_mod 
USE prompt_mod 
USE precision_mod
!                                                                       
IMPLICIT none 

!                                                                       
integer, parameter :: MAXMASK = 4
real(kind=PREC_DP), parameter :: EPS = 1.0D-8
CHARACTER(LEN=23), PARAMETER :: tempfile ='internal.temporary.stru'
CHARACTER(LEN=23), PARAMETER :: stackfile='internal.stacklist.stru'
CHARACTER(LEN=25), PARAMETER :: stacksimple='internal.stacksimple.list'
character(len=PREC_STRING) :: strucfile
INTEGER i, j, k 
INTEGER         :: nlayers
integer         :: nanis   ! Number of different ADP
INTEGER         :: n_mole  ! number of molecules in input file
INTEGER         :: n_type  ! number of molecule types in input file
INTEGER         :: n_atom  ! number of molecule atoms in input file
INTEGER         :: natoms
!integer, dimension(3) :: n_cells
INTEGER         :: nscats
logical,dimension(0:MAXMASK) :: uni_mask
LOGICAL         :: need_alloc = .false.
LOGICAL         :: lprev 
REAL(kind=PREC_DP), dimension(ST_MAXTYPE) :: prob  ! (ST_MAXTYPE) 
REAL(kind=PREC_DP), dimension(ST_MAXTYPE) :: prob_n! (ST_MAXTYPE) 
REAL(kind=PREC_DP) :: ptot_n 
REAL(kind=PREC_DP) :: u (3), v (3) 
REAL(kind=PREC_DP) :: r1, r2 
REAL(kind=PREC_DP) :: st_trans_cur (3) 
!                                                                       
!                                                                       
CALL no_error 
!                                                                       
!     reset counter for number of each layer type                       
!                                                                       
DO i = 1, ST_MAXTYPE 
   st_number (i) = 0 
ENDDO 
!                                                                       
uni_mask(0)   = .true.
uni_mask(1:3) = .true.
uni_mask(4)   = .false.
!
!     Read layer types from file or create from probability matrix      
!                                                                       
if_distr: IF(st_distr.eq.ST_DIST_MATRIX) then 
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
   IF(st_first==0) THEN  ! determin first layer from random
      CALL stack_prob (ptot_n, prob, prob_n, lprev, i) 
      CALL stack_type (prob, i) 
      st_origin (1, 1) = 0.0 
      st_origin (2, 1) = 0.0 
      st_origin (3, 1) = 0.0 
      st_number (st_type (1) ) = 1 
   ELSE
      IF(st_first>0 .AND. st_first<=st_ntypes) THEN
         st_type (i) = st_first
      ELSE 
         ier_num = -139 
         ier_typ = ER_APPL 
         ier_msg(1)='The first layer type must be >= 1 and'
         ier_msg(2)='less or equal to the number of layer types'
         RETURN
      ENDIF
   ENDIF
ELSEIF (st_distr.eq.ST_DIST_FILE.or.st_distr.eq.ST_DIST_LIST) THEN if_distr
!
! Always read distribution file into internal memory
! Is done upon 'distr' command
!
!                                                                       
!     --Stacking fault origins are read from file                       
!                                                                       
   IF(st_infile(1:8)=='internal') THEN
            CALL testfile_internal (st_infile,st_natoms, st_ncatoms, &        ! Get size of internal structure
             st_nscat, nanis, n_mole, n_type, n_atom)
   ELSE
      ier_num = -182
      ier_typ = ER_APPL
      ier_msg(1) = 'Stacking faults should have been internal'
      ier_msg(2) = 'Stacking fault distribution file error 01'
      return
!           CALL test_file ( st_infile,st_natoms,st_nscat, n_mole, n_type,&
!                            n_atom, n_cells, -1, .true. )
   ENDIF
!        CALL test_file ( st_infile, st_natoms, st_nscat, n_mole,  &
!                         n_type, n_atom, -1, .true. )
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
      CALL alloc_stack (ST_MAXTYPE, nlayers, st_rot_status )
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
   st_type(1) = st_iscat(1,1) 
   st_nlayer  = st_natoms 
                                                                    
!DBG                                                                    
!DBG      write (output_io,*) 'st_natoms ', st_natoms                   
   IF (st_distr.eq.ST_DIST_LIST) then 
      DO i = 1, st_natoms 
         st_type(i) = st_iscat(1,i) 
         st_number (st_type (i) ) = st_number (st_type (i) ) + 1 
         DO j = 1, 3 
            st_origin (j, i) = st_pos (j, i) 
         ENDDO 
      ENDDO 
!
!     Copy the amount of each layer type into res
!
      DO i = 1, st_nscat
         res_para(i) = st_number(i)
      ENDDO
      res_para(0) = st_nscat
      st_nlayer = st_natoms 
      ier_num = 0 
      RETURN 
   ENDIF 
ENDIF if_distr
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
            st_mod(k, 3) = st_mod(k, 3) + st_corr(i, j) * st_trans(i, j, k)
         ENDDO 
      ENDDO 
   ENDDO 
   DO i = 1, 3 
      st_mod (i, 3) = st_mod (i, 3) / ptot_n 
   ENDDO 
   IF (st_mod(1, 3) .eq.0.0 .and.   &
       st_mod(2, 3) .eq.0.0 .and.   &
       st_mod(3, 3) .eq.0.0        ) then                                           
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
!
CALL matinv3(st_mod, st_inv) 
IF (ier_num == -45) then 
   CALL no_error 
   ier_num = - 56 
   ier_typ = ER_APPL 
   ier_msg(1) = 'Check modulus vectors and average direction'
   RETURN 
ENDIF 
!                                                                       
!     --Initialise rotation disorder                                    
!                                                                       
i = 1 
IF (st_rot_status) then 
   st_rot_ang_no (i) = gasdev(DBLE(st_rot_si_no))
   st_rot_ang_m1 (i) = gasdev(DBLE(st_rot_si_m1))
   st_rot_ang_m2 (i) = gasdev(DBLE(st_rot_si_m2))
ENDIF 
!                                                                       
!     Loop over all layers                                              
!                                                                       
!DBG                                                                    
WRITE (output_io, * ) ' Start Loop' , st_nlayer, ier_num, ier_typ
WRITE (output_io, 2000) ' Origins', st_type (1), st_origin(1, 1)&
                        , st_origin(2, 1), st_origin(3, 1)                             
loop_layer: DO i = 2, st_nlayer 
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
      st_type(i) = st_iscat(1,i) 
   ENDIF 
!                                                                       
!     --Determine current translation vector                            
!                                                                       
   st_trans_cur(1) = st_trans(st_type(i - 1), st_type (i), 1) + &
                     gasdev (DBLE(st_sigma (st_type (i - 1), st_type (i), 1)))         
   st_trans_cur(2) = st_trans(st_type(i - 1), st_type (i), 2) + &
                     gasdev (DBLE(st_sigma (st_type (i - 1), st_type (i), 2)))         
   st_trans_cur(3) = st_trans(st_type(i - 1), st_type (i), 3) + &
                     gasdev (DBLE(st_sigma (st_type (i - 1), st_type (i), 3)))         
!                                                                       
!     --Initialise rotation disorder                                    
!                                                                       
   IF (st_rot_status) then 
!                                                                       
!     ----Set new rotation angles for all layers or                     
!         only if the last layer was different than the current         
!                                                                       
      IF (st_rot_mode) then 
         st_rot_ang_no (i) = gasdev (DBLE(st_rot_si_no) )
         st_rot_ang_m1 (i) = gasdev (DBLE(st_rot_si_m1) )
         st_rot_ang_m2 (i) = gasdev (DBLE(st_rot_si_m2) )
      ELSE 
         IF (st_type (i) .ne.st_type (i - 1) ) then 
            st_rot_ang_no (i) = gasdev (DBLE(st_rot_si_no) )
            st_rot_ang_m1 (i) = gasdev (DBLE(st_rot_si_m1) )
            st_rot_ang_m2 (i) = gasdev (DBLE(st_rot_si_m2) )
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
!  Adding the modulo vector ensures that origins stay in [0,1[ range
!                                                                       
   st_origin (1, i) = st_origin (1, i - 1) + st_trans_cur (1) !+ st_mod(1,1) + st_mod(1,2) 
   st_origin (2, i) = st_origin (2, i - 1) + st_trans_cur (2) !+ st_mod(2,1) + st_mod(2,2)
   st_origin (3, i) = st_origin (3, i - 1) + st_trans_cur (3) !+ st_mod(3,1) + st_mod(3,2)
   if(st_mod_sta) then      ! Stack modulo is on, add modulo 
      st_origin (1, i) = st_origin (1, i) + st_mod(1,1) + st_mod(1,2)
      st_origin (2, i) = st_origin (2, i) + st_mod(2,1) + st_mod(2,2)
      st_origin (3, i) = st_origin (3, i) + st_mod(3,1) + st_mod(3,2)
   endif
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
         st_origin(j, i) = st_origin(j, i) + r1 * st_mod(j, 1)       &
                           + r2 * st_mod (j, 2) + st_off (j)         &
                           + gasdev (DBLE(st_sigma_off (j))) 
      ENDDO 
   ENDIF 
!                                                                       
!     --If necessary apply modulo function                              
!                                                                       
   IF (st_mod_sta) then 
      DO j = 1, 3 
         v (j) = st_origin (j, i) - st_origin (j, 1) 
      ENDDO 
!     CALL trans (v, st_inv, u, 3) 
!     EPS ensures that the origins are in the range [-EPS:1-EPS]. This helps to avoid 
!     issues with origins at 0/3 ; 3/3, wich should all be at 0.
      u = matmul(st_inv, v)
      DO j = 1, 3 
         st_origin(j, i) = st_origin(j, i) - REAL(int(u(1)+EPS)) * st_mod(j, 1) &
                                           - REAL(int(u(2)+EPS)) * st_mod(j, 2)
      ENDDO 
   ENDIF 
!                                                                       
!     --Increment number of layers of this type                         
!                                                                       
   st_number(st_type(i) ) = st_number(st_type(i) ) + 1 
   WRITE (output_io, 2000) ' Origins', st_type(i) , st_origin(1, i)             &
      , st_origin (2, i) , st_origin (3, i)
!
enddo loop_layer
!
!     Copy the amount of each layer type into res
!
DO i = 1, st_ntypes
   res_para(i) = st_number(i)
ENDDO
res_para(0) = st_ntypes
! Estimate number of layers per unit cell
! Determine growth direction as set by 'aver' command
if(    abs(st_t_aver(1))<EPS .and. abs(st_t_aver(2))<EPS) then
   i = 1
   j = 3
elseif(abs(st_t_aver(1))<EPS .and. abs(st_t_aver(3))<EPS) then
   i = 1
   j = 2
elseif(abs(st_t_aver(2))<EPS .and. abs(st_t_aver(3))<EPS) then
   i = 1
   j = 1
elseif(abs(st_t_aver(1))<EPS) then
   i = 2
   j = 1
elseif(abs(st_t_aver(2))<EPS) then
   i = 2
   j = 2
elseif(abs(st_t_aver(3))<EPS) then
   i = 2
   j = 3
else
   i = 3
   j = 0
endif
if(i==1) then            ! Single growth component
   st_ncunit =nint(st_nlayer/(abs(maxval(st_origin(j,:))-minval(st_origin(j,:)))))
elseif(i==2) then        ! Double growth direction
   if(j==1) then         ! yz
      st_ncunit =nint( &
                st_nlayer/( max(abs(maxval(st_origin(2,:))-minval(st_origin(2,:))),  &
                                abs(maxval(st_origin(3,:))-minval(st_origin(3,:))))))
   elseif(j==2) then
      st_ncunit =nint( &
                st_nlayer/( max(abs(maxval(st_origin(1,:))-minval(st_origin(1,:))),  &
                                abs(maxval(st_origin(3,:))-minval(st_origin(3,:))))))
   elseif(j==3) then
      st_ncunit =nint( &
                st_nlayer/( max(abs(maxval(st_origin(1,:))-minval(st_origin(1,:))),  &
                                abs(maxval(st_origin(2,:))-minval(st_origin(2,:))))))
   endif
else
   st_ncunit =nint( &
                st_nlayer/( max(abs(maxval(st_origin(1,:))-minval(st_origin(1,:))),  &
                                abs(maxval(st_origin(2,:))-minval(st_origin(2,:))),  &
                                abs(maxval(st_origin(3,:))-minval(st_origin(3,:))))  &
               ))
endif
!
!     Save the list of origins as internal file
!
st_nanis = 0
CALL save_internal(tempfile)        ! temporarily save current structure
!
natoms = st_nlayer
nscats = st_ntypes
need_alloc = .false.
IF(natoms > NMAX ) THEN
   need_alloc = .true.
ENDIF
IF(nscats > MAXSCAT) THEN
   need_alloc = .true.
ENDIF
IF( need_alloc) THEN
   CALL alloc_crystal_scat (nscats)
   CALL alloc_crystal_nmax (natoms)
   IF ( ier_num /= 0 ) RETURN
ENDIF
cr_natoms = st_nlayer
cr_nscat  = st_ntypes
cr_nanis  = 0
cr_dw      (:) = 0.1
cr_scat_int(:) = .FALSE.
cr_scat_equ(:) = .FALSE.
cr_delf_int(:) = .FALSE.
cr_at_equ  (:) = ' '
DO i=1,st_ntypes
   WRITE(cr_at_lis(i),3000) i
   cr_scat_equ(i) = .TRUE.
   cr_at_equ(i)   = 'H'
ENDDO
DO i=1,st_nlayer
   cr_pos  (:,i) = st_origin (:, i)
   cr_iscat(1,i) = st_type(i)
   cr_mole (  i) = 0
   cr_surf (:,i) = 0
   cr_magn (:,i) = 0
   cr_prop (  i) = 1
ENDDO
mole_num_mole = 0
mole_num_type = 0
mole_num_atom = 0
CALL save_internal(stackfile)
!
!     Save simplified list for correlation analysis
!
cr_pos(:,:) = 0.0
DO i=1,st_nlayer
   cr_pos(3,i) = REAL(i-1)
ENDDO
CALL save_internal(stacksimple)
CALL readstru_internal(MAXMASK, tempfile, uni_mask)    ! Restore current structure
CALL store_remove_single(tempfile, ier_num)  ! Cleanup internal storage
IF(ier_num/=0) THEN
   ier_typ = ER_APPL
   ier_msg(1) = 'Could not remove temporary internal storage'
   ier_msg(2) = 'in stack_create'
   ier_msg(3) = 'Please document and report'
ENDIF

if (st_distr.eq.ST_DIST_FILE.or.st_distr.eq.ST_DIST_LIST) THEN
if(st_infile(1:20)/='internal_stack_list.') then
   strucfile = st_infile
   call store_remove_single(strucfile, ier_num)
   strucfile = ' '
   strucfile(1:len_trim(st_infile)) = st_infile(21:len_trim(st_infile))
   st_infile = strucfile
endif
endif
!                                                                       
!     do i=1,st_ntypes                                                  
!       write (output_io,2020) i,st_number(i)                           
!     ENDDO                                                             
 2000 FORMAT    (a,i4,3f11.4) 
 3000 FORMAT('L',I3.3)
!
END SUBROUTINE do_stack_create                
!
!*****7*****************************************************************
!
SUBROUTINE st_init_prob (ptot_n, prob_n) 
!-                                                                      
!     Determines the probabilities from the correlation matrix          
!+                                                                      
USE discus_config_mod 
USE stack_mod 
use precision_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER ::  i, j 
REAL(kind=PREC_DP) :: prob_n (ST_MAXTYPE) 
REAL(kind=PREC_DP) :: ptot_n 
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
!
!*****7*****************************************************************
!
SUBROUTINE stack_prob (ptot_n, prob, prob_n, lprev, l) 
!-                                                                      
!     Set up the cummulative probabilities.                      
!     if a previous layer exists, add probabilities                     
!     according to the correlation matrix, ELSE                         
!     an average probality for the                                      
!     specific layer type i.e. sum corr(i,m) /st_type                   
!     is added.                                                         
!+                                                                      
USE discus_config_mod 
USE stack_mod 
!
use precision_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER :: m, l 
LOGICAL :: lprev 
REAL(kind=PREC_DP), dimension(ST_MAXTYPE) :: prob  !(ST_MAXTYPE)
REAL(kind=PREC_DP), dimension(ST_MAXTYPE) :: prob_n ! (ST_MAXTYPE) 
real(kind=PREC_DP)                        :: ptot_n 
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
!
!*****7*****************************************************************
!
SUBROUTINE stack_type (prob, i) 
!-                                                                      
!     Determine type of layer from weighted probabilities               
!+                                                                      
USE discus_config_mod 
USE stack_mod 
USE random_mod
USE errlist_mod 
USE lib_random_func
!
use precision_mod
!                                                                       
IMPLICIT none 
!                                                                       
REAL(kind=PREC_DP), dimension(ST_MAXTYPE) :: prob !(ST_MAXTYPE) 
!                                                                       
INTEGER :: i, m 
REAL(kind=PREC_DP) :: ad, s, ptot 
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
!
!*****7*****************************************************************
!
SUBROUTINE do_stack_fill 
!-                                                                      
!     Creates the crystal with stacking faults                          
!+                                                                      
USE discus_config_mod 
USE discus_allocate_appl_mod
USE chem_mod
USE crystal_mod 
USE gen_add_mod 
USE sym_add_mod 
USE molecule_mod 
USE read_internal_mod
USE discus_save_mod 
use dis_estimate_mod, ONLY : estimate_ncells
!use save_menu       , only: save_store_setting, save_default_setting, save_restore_setting, save_internal
use prep_anis_mod
USE stack_mod  
use stack_cr_mod
USE structur  
USE symm_sup_mod
USE spcgr_apply
!
USE update_cr_dim_mod
USE errlist_mod 
USE lib_length
use precision_mod
USE support_mod
!
IMPLICIT none 
!                                                                       
integer, parameter :: MAXMASK = 4
INTEGER, PARAMETER                   :: ist = 7 
real(kind=PREC_DP), parameter        :: TOL=1.0D-5
!                                                                       
character(len=PREC_STRING) :: strucfile
!                                                                       
integer :: dim_natoms  ! Internal array size
integer :: dim_ncatoms  ! Internal array size
integer :: dim_nscat  ! Internal array size
integer :: dim_nanis  ! Internal array size
integer :: dim_n_mole  ! Internal array size
integer :: dim_n_type  ! Internal array size
integer :: dim_n_atom  ! Internal array size
integer          :: natom = 0
integer          :: ncatoms = 0
integer          :: nscat = 0
integer          :: nanis = 0
integer          :: max_nanis = 0
INTEGER          :: natoms=0, max_natoms=0
INTEGER          :: nscats=0, max_nscats=0
integer          :: max_ncatoms = 0
INTEGER          ::         max_n_mole=0
INTEGER          ::         max_n_type=0
INTEGER          ::         max_n_atom=0
INTEGER         :: i, j, k, kk
INTEGER         :: iatom 
INTEGER         :: n_mole=0  ! number of molecules in input file
INTEGER         :: n_type=0  ! number of molecule types in input file
INTEGER         :: n_atom=0  ! number of molecule atoms in input file
integer, dimension(3) :: n_cells
LOGICAL         :: lread, lout 
!LOGICAL           :: need_alloc = .false. 
logical         :: lsuccess   ! Found old ADP entry
logical,dimension(0:MAXMASK) :: uni_mask
!
integer, dimension(:,:), allocatable :: adp_look  ! Lookup table for st_iscat(3,:) to get correct ADP's
!
! needed for estimate_ncells
real(kind=PREC_DP), dimension(3, 2):: pdt_dims
integer           , dimension(3)   :: pdt_ilow      ! Unit cell dimensions in periodic
integer           , dimension(3)   :: pdt_ihig      ! low and high inidce
integer                            :: pdt_ncells    ! Number of cells in periodic crystal volume
!
real(kind=PREC_DP), dimension(3,2) :: st_dims       ! Crystal dimensions to restrict if module(atom) is active
integer           , dimension(3)   :: st_idims
real(kind=PREC_DP), dimension(3)   :: st_adims
integer :: ier ! allocation error
!                                                                       
!                                                                       
!     ----read corresponding layer                                      
!                                                                       
CALL rese_cr 
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
uni_mask(0)   = .true.
uni_mask(1:3) = .true.
uni_mask(4)   = .false.
!
st_dims(:,1) = -huge(1.0D0)
st_dims(:,2) =  huge(1.0D0)
st_nscat = 0
st_natoms = 0
!
nanis = 0
max_nanis = 0
!
if(st_mod_atom/=ST_ATOM_OFF) then          ! Atom modulo function is active
   do i=1,st_ntypes                        ! Get layer dimensions
     call stru_internal_get_cr_dim(st_layer(i), cr_dim)
     st_dims(:,1) = max(st_dims(:,1), cr_dim(:,1))
     st_dims(:,2) = min(st_dims(:,2), cr_dim(:,2))
   enddo                                   ! Current maximum dimensions, used for modulo_shift
   st_adims(:) = real(nint(real(ceiling(st_dims(:,2)),kind=PREC_DP) - &
                           real(floor  (st_dims(:,1)),kind=PREC_DP)))
!
   if(st_mod_atom==ST_ATOM_ON) then        ! Keep current dimensions
      st_dims(:,1) = real(floor  (st_dims(:,1)),kind=PREC_DP)
      st_dims(:,2) = real(ceiling(st_dims(:,2)),kind=PREC_DP)
!
   elseif(st_mod_atom==ST_ATOM_STRICT) then  ! Cut crystal 1 to 2 cells smaller for stricter smoothness
      st_dims(:,1) = real(ceiling(st_dims(:,1)),kind=PREC_DP)
      st_dims(:,2) = real(floor  (st_dims(:,2)),kind=PREC_DP)
   endif
   st_idims(:) = nint(st_dims(:,2)-st_dims(:,1))
endif
!
more1: IF (st_nlayer.ge.1) then 
!
   max_natoms = 0
   max_nscats = 0
   max_n_mole = 0
   max_n_type = 0
   max_n_atom = 0
!        DO i = 1, st_ntypes
!        enddo
   DO i = 1, st_ntypes
      IF(st_internal(i) ) THEN
         CALL testfile_internal(st_layer (i ), natoms, ncatoms, nscats, nanis, &
                                n_mole, n_type, n_atom)
      ELSE
         ier_num = -182
         ier_typ = ER_APPL
         ier_msg(1) = 'Stacking faults should have been internal'
         ier_msg(2) = 'Stacking fault layer file error 02'
         return
!               CALL test_file ( st_layer (i ), natoms, nscats, n_mole, &
!                                n_type, n_atom, n_cells, -1*i, .false.)
      ENDIF
      If(ier_num/=0) RETURN
!
      max_natoms = MAX(max_natoms, natoms)
      max_nscats = MAX(max_nscats, nscats)
      max_ncatoms = MAX(max_ncatoms, nscats)
      max_nanis  = max(max_nanis, nanis)
      max_n_mole = MAX(max_n_mole, n_mole)
      max_n_type = MAX(max_n_type, n_type)
      max_n_atom = MAX(max_n_atom, n_atom)
   ENDDO
   call alloc_stack_crystal(max_nscats, max_natoms)
   ST_MMAX     = max_natoms
   ST_MAX_SCAT = max_nscats
!
   natoms = st_nlayer * max_natoms
   nscats = max_nscats
!  need_alloc = .false.
!  IF(natoms >= NMAX ) THEN
!     natoms = MAX(NINT(natoms*1.05),natoms+10,NMAX)
!     need_alloc = .true.
!  ENDIF
!  IF(nscats >= MAXSCAT) THEN
!     nscats = MAX(NINT(nscats*1.05),nscats+ 5,MAXSCAT)
!     need_alloc = .true.
!  ENDIF
!  IF( need_alloc) THEN
!     CALL alloc_crystal_scat (nscats)
!     CALL alloc_crystal_nmax (natoms)
!     IF ( ier_num /= 0 ) then
!        RETURN
!     endif
!  ENDIF
!  call alloc_crystal_nmax(st_nlayer*max_natoms)
!  call alloc_crystal_scat(max_nscats)
!   call alloc_unitcell(MAX_ncatoms)
!  call alloc_anis(MAX_nanis)
!  Allocate lookup table
   allocate(adp_look(st_ntypes, max_nanis))
   adp_look = 0
!
!        IF more molecules have been read than were allocated
!
   IF(max_n_mole>MOLE_MAX_MOLE .or. max_n_type>MOLE_MAX_TYPE .or.   &
      max_n_atom>MOLE_MAX_ATOM                          ) THEN
      n_mole = MAX(max_n_mole +20 ,MOLE_MAX_MOLE)
      n_type = MAX(max_n_type +10 ,MOLE_MAX_TYPE)
      n_atom = MAX(max_n_atom +200,MOLE_MAX_ATOM)
      CALL alloc_molecule(1, 1,n_mole,n_type,n_atom)
      IF ( ier_num /= 0 ) then
         RETURN
      endif
   ENDIF
!
!                                                                       
!     --The first layer is read by the full readstru routine to         
!        ensure correct lattice parameters                              
!                                                                       
   i = 1 
   lout = .false. 
!
   IF(st_internal(st_type(i)) ) THEN
      CALL readstru_internal(MAXMASK, st_layer(st_type(i)), uni_mask)!, &
!                 NMAX, MAXSCAT, MOLE_MAX_MOLE, &
!                 MOLE_MAX_TYPE, MOLE_MAX_ATOM )
      IF (ier_num.ne.0) then 
         RETURN 
      ENDIF 
   ELSE 
      ier_num = -182
      ier_typ = ER_APPL
      ier_msg(1) = 'Stacking faults should have been internal'
      ier_msg(2) = 'Stacking fault layer file error 03'
      return
!           CALL readstru (NMAX, MAXSCAT, st_layer (st_type (i) ), cr_name,&
!           cr_spcgr, cr_set, cr_a0, cr_win, cr_natoms, cr_nscat, cr_dw, cr_occ, cr_at_lis,&
!           cr_pos, cr_mole, cr_surf, cr_magn, cr_iscat, cr_prop, cr_dim, cr_magnetic, as_natoms, as_at_lis, as_dw,&
!           as_pos, as_iscat, as_prop, sav_ncell, sav_r_ncell, sav_ncatoms,&
!           spcgr_ianz, spcgr_para)                                        
!           IF (ier_num.ne.0) then 
!              RETURN 
!           ENDIF 
   ENDIF 
!
!                                                                       
   CALL setup_lattice (cr_a0, cr_ar, cr_eps, cr_gten, cr_reps,    &
   cr_rten, cr_win, cr_wrez, cr_v, cr_vr, lout, cr_gmat, cr_fmat, &
   cr_cartesian, cr_tran_g, cr_tran_gi, cr_tran_f, cr_tran_fi)
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
!do j=1, cr_natoms
!write(*,'(a,i4,3f9.5,3i4, 6f8.4)') cr_at_lis(cr_iscat(1,j)), j, cr_pos(:,j), cr_iscat(:,j), cr_anis_full(:,cr_iscat(3,j))
!enddo 
!write(*,*) ' cr_nanis', cr_nanis
!do j=1,cr_nanis
!write(*,'(i3,6f8.4)') j, cr_anis_full(:,j)
!write(*,'(i3,6f8.4)') 1, cr_prin(:,1,j)
!write(*,'(i3,6f8.4)') 2, cr_prin(:,2,j)
!write(*,'(i3,6f8.4)') 3, cr_prin(:,3,j)
!enddo
!  First layer allocates its own size only, increase to full need
   call alloc_crystal_nmax(st_nlayer*max_natoms)
   call alloc_crystal_scat(max_nscats)
   call alloc_unitcell(MAX_ncatoms)
   call alloc_anis(MAX_nanis)
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
   call stack_atom_modulo(iatom, st_dims, st_idims, st_adims)
!write(*,*)' MOLECULES ' , mole_num_mole, mole_num_curr, mole_num_atom, mole_num_type
!do j= 1, mole_num_mole
!  write(*,*) ' mole_len, mole_typ, mole_off ',mole_len(j), mole_type(j), mole_off(j)
!  write(*,*) ' mole_cont ', (mole_cont(mole_off(j)+k),k=1, mole_len(j))
!enddo
!
!  Read all headers of all layer types once, to get global list of atom types
!
!write(*,*) ' READING LAYER TYPES 1: ', st_ntypes, st_nscat
!read(*,*) i
   loop_st_ntypes: do i=1, st_ntypes
      call readstru_size_int(st_layer(i),  dim_natoms, dim_ncatoms, &
              dim_nscat, dim_nanis, dim_n_mole, dim_n_type, dim_n_atom, &
              st_natoms, sav_ncatoms, st_nscat, st_nanis, n_mole, n_type, n_atom)
      if(allocated(st_is_sym)) deallocate(st_is_sym)
      allocate(st_is_sym(1:sav_ncatoms))
      st_is_sym = 1
      st_nanis=max(1, st_nanis)
      call alloc_anis_generic(st_nanis, st_anis_full, st_prin, ier)
      CALL stru_readheader_internal(st_layer(i), DIM_NSCAT, DIM_NCATOMS, DIM_NANIS, st_name, &
               st_spcgr, st_spcgr_set, st_set, st_iset,                         &
               st_at_lis, st_nscat, st_dw, st_occ, st_anis, st_is_sym,          &
               st_nanis, st_anis_full, st_prin,                                 &
               st_a0, st_win,                                                   &
               sav_ncell, sav_r_ncell, sav_ncatoms, spcgr_ianz, spcgr_para,     &
               GEN_ADD_MAX, gen_add_n, gen_add_power, gen_add,                  &
               SYM_ADD_MAX, sym_add_n, sym_add_power, sym_add)
!
!write(*,*) ' LAYER TYPE ', i, st_nscat, cr_nscat, ubound(st_at_lis), MAXSCAT
      do j=1, st_nscat     ! Loop over all atom types in this layer type
         kk = -1
         loop_type_i: do k=1,cr_nscat
            if(cr_at_lis(k)==st_at_lis(j) .and. abs(cr_dw(k)-st_dw(j))<TOL) then
               kk = k
               exit loop_type_i
            endif
         enddo loop_type_i
         if(kk<0) then         ! Atom type was not found, need new type 
            if(cr_nscat==MAXSCAT) then          ! Need to allocate
               nscat = max(cr_nscat+1, MAXSCAT)
               natom = max(cr_natoms+st_natoms, NMAX)
               call alloc_crystal_scat(nscat)
               call alloc_crystal_nmax(nmax)
               MAXSCAT = nscat
               NMAX    = nmax
            endif
            kk       = cr_nscat + 1
            cr_nscat = cr_nscat + 1                  ! Increment atom types
            cr_at_lis(kk) = st_at_lis(j)             ! Set new atom type
            cr_dw    (kk) = st_dw(j )                ! Set new ADP
         endif
      enddo
!
!     Build ADP Lookup table
      do j=1, st_nanis
         call lookup_anis(.TRUE., cr_nanis, cr_anis_full, cr_prin, st_anis_full(1:6,j), st_prin(:,:,j), kk, lsuccess)
         adp_look(i,j) = kk
      enddo
!
   enddo loop_st_ntypes
!write(*,*) ' CR_NANIS', cr_nanis
!do j=1,cr_nanis
!write(*,'(i3,6f8.4)') j, cr_anis_full(:,j)
!write(*,'(i3,6f8.4)') 1, cr_prin(:,1,j)
!write(*,'(i3,6f8.4)') 2, cr_prin(:,2,j)
!write(*,'(i3,6f8.4)') 3, cr_prin(:,3,j)
!enddo
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
         call readstru_size_int(st_layer(st_type(i)),  dim_natoms, dim_ncatoms, &
              dim_nscat, dim_nanis, dim_n_mole, dim_n_type, dim_n_atom, &
              st_natoms, sav_ncatoms, st_nscat, st_nanis, n_mole, n_type, n_atom)
!        call readstru_size_int(st_layer(st_type(i)), st_natoms, sav_ncatoms, &
!                st_nscat, st_nanis, n_mole, n_type, n_atom)
         if(allocated(st_is_sym)) deallocate(st_is_sym)
         allocate(st_is_sym(1:sav_ncatoms))
         st_is_sym = 1
         st_nanis=max(1, st_nanis)
         call alloc_anis_generic(st_nanis, st_anis_full, st_prin, ier)
         CALL stru_readheader_internal (st_layer(st_type(i)), MAXSCAT, DIM_NCATOMS, DIM_NANIS, st_name, &
               st_spcgr, st_spcgr_set, st_set, st_iset,                         &
               st_at_lis, st_nscat, st_dw, st_occ, st_anis, st_is_sym,          &
               st_nanis, st_anis_full, st_prin,                                 &
               st_a0, st_win,                                                   &
               sav_ncell, sav_r_ncell, sav_ncatoms, spcgr_ianz, spcgr_para,     &
               GEN_ADD_MAX, gen_add_n, gen_add_power, gen_add,                  &
               SYM_ADD_MAX, sym_add_n, sym_add_power, sym_add)
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
         j = st_natoms    ! Use value from header for array size
         st_natoms = 0
         CALL struc_read_atoms_internal (st_layer(st_type(i)),j,&
                     st_natoms, st_pos, st_iscat, st_prop, st_surf, st_magn, st_mole)
         CLOSE (ist) 
         IF (ier_num.ne.0) then 
            RETURN 
         ENDIF 
      ELSE internal
         ier_num = -182
         ier_typ = ER_APPL
         ier_msg(1) = 'Stacking faults should have been internal'
         ier_msg(2) = 'Stacking fault layer file error 04'
         return
!INTERN!                                                                       
!INTERN!     --Open file, read header                                          
!INTERN!                                                                       
!INTERN         CALL oeffne (ist, st_layer (st_type (i) ) , 'old') 
!INTERN         IF (ier_num /= 0) then 
!INTERN            RETURN
!INTERN         ENDIF
!INTERN!                                                                       
!INTERN!     ------Read header of structure file                               
!INTERN!                                                                       
!INTERN            gen_add_n = 0 
!INTERN            sym_add_n = 0 
!INTERN            CALL stru_readheader (ist, MAXSCAT, cr_name,   &
!INTERN            cr_spcgr, cr_set, cr_at_lis, cr_nscat, cr_dw, cr_occ, cr_a0, cr_win, &
!INTERN            sav_ncell, sav_r_ncell, sav_ncatoms, spcgr_ianz, spcgr_para, &
!INTERN            AT_MAXP, at_ianz, at_param)
!INTERN!                                                                       
!INTERN            IF (ier_num.ne.0) then 
!INTERN               RETURN 
!INTERN            ENDIF 
!INTERN!                                                                       
!INTERN!     ------Now only the atoms are read                                 
!INTERN!                                                                       
!INTERN            ier_num = 0 
!INTERN            ier_typ = ER_NONE 
!INTERN!                                                                       
!INTERN            CALL struc_read_atoms (NMAX, MAXSCAT, cr_natoms, cr_nscat,  &
!INTERN            cr_dw, cr_occ, cr_at_lis, cr_pos, cr_iscat, cr_mole, cr_surf, &
!INTERN            cr_magn, cr_prop, cr_dim, cr_magnetic, &
!INTERN            as_natoms, as_at_lis, as_dw, as_occ, as_pos, as_iscat, as_prop, &
!INTERN            AT_MAXP, at_ianz, at_param)
!INTERN            CLOSE (ist) 
!INTERN            IF (ier_num.ne.0) then 
!INTERN               RETURN 
!INTERN            ENDIF 
      ENDIF internal
!                                                                       
!     -------Copy atoms from st_array
!     ------Add origin of layer                                         
!                                                                       
      DO j = 1    , st_natoms 
         kk = -1
         loop_type: do k=1,cr_nscat
            if(cr_at_lis(k)==st_at_lis(st_iscat(1,j)) .and. abs(cr_dw(k)-st_dw(st_iscat(1,j)))<TOL) then
               kk = k
               exit loop_type
            endif
         enddo loop_type
         if(kk<0) then         ! Atom type was not found, need new type 
            if(cr_nscat==MAXSCAT) then          ! Need to allocate
               nscat = max(cr_nscat+1, MAXSCAT)
               natom = max(cr_natoms+st_natoms, NMAX)
               call alloc_crystal_scat(nscat)
               call alloc_crystal_nmax(nmax)
               MAXSCAT = nscat
               NMAX    = nmax
            endif
            kk       = cr_nscat + 1
            cr_nscat = cr_nscat + 1                  ! Increment atom types
            cr_at_lis(kk) = st_at_lis(st_iscat(1,j))             ! Set new atom type
            cr_dw    (kk) = st_dw(st_iscat(1,j) )                ! Set new ADP
         endif
         cr_iscat(1,iatom-1+j) = kk 
         cr_prop (iatom-1+j) = st_prop(j)
         cr_surf (:,iatom-1+j) = 0                   ! Not on a surface
         cr_magn (:,iatom-1+j) = st_magn(:,j)
         cr_mole (  iatom-1+j) = 0 !st_mole(  j)        ! MOLECULE NEEDS WORK
         DO k = 1, 3 
            cr_pos (k, iatom-1+j) = st_pos (k, j) + st_origin (k, i) 
         ENDDO 
!        Assign correct ADP from loop up table
         cr_iscat(3,iatom-1+j) = adp_look(st_type(i), st_iscat(3,j))
      ENDDO 
      cr_natoms = cr_natoms + st_natoms
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
!write(*,*)' MOLECULES ' , mole_num_mole, mole_num_curr, mole_num_atom, mole_num_type
!do j= 1, mole_num_mole
!  write(*,*) ' mole_len, mole_typ, mole_off ',mole_len(j), mole_type(j), mole_off(j)
!  write(*,*) ' mole_cont ', (mole_cont(mole_off(j)+k),k=1, mole_len(j))
!enddo
!     Check if user wants atom modulo as well
      call stack_atom_modulo(iatom, st_dims, st_idims, st_adims)
   ENDDO layers
ENDIF more1
!                                                                       
CLOSE (ist) 
!
if(ier_num==0) then           ! No errors occured
!                                                                       
!     Update Crystal dimension                                          
!                                                                       
   CALL update_cr_dim 
   chem_purge = .TRUE.     ! The crystal is most likely NOT periodic.
                              ! In the rare circumstances that it is the user
                              ! has to turn this on explicitly
   chem_quick = .FALSE.
   chem_period(:) = .FALSE.
!
   call estimate_ncells(n_cells, pdt_dims, pdt_ilow, pdt_ihig, pdt_ncells) ! Estimate number of unit cells in the crystal
   n_cells = max(n_cells,1)
   cr_icc = n_cells
   cr_ncatoms = nint(real(cr_natoms) /real(n_cells(1)*n_cells(2)*n_cells(3)))
!
! --- Clean up internal memory
!                                                                       
   loop_clean: do i = 1, st_ntypes           ! Read all layers into internal memory
      if(st_layer(i)(1:21)=='internal_stack_layer.') then         ! An internal file, delete
         strucfile = st_layer(i)
!        call store_remove_single(strucfile, ier_num)
         if(ier_num==-113) ier_num = 0          ! File might have been cleaned up already
      endif
   enddo loop_clean
endif
!
deallocate(adp_look)
!                                                                       
END SUBROUTINE do_stack_fill                  
!
!*****7*****************************************************************
!
subroutine  stack_atom_modulo(iatom, st_dims, st_idims, st_adims)
!-
!  Perform the modulo function for atom positions
!  If an atom slides outside the intended crystal size st_dim, its removes.
!+
!
use crystal_mod
use stack_mod
!
implicit none
!
integer, intent(in) :: iatom
real(kind=PREC_DP), dimension(3,2), intent(in) :: st_dims       ! Crystal dimensions to restrict if module(atom) is active
integer           , dimension(3)  , intent(in) :: st_idims
real(kind=PREC_DP), dimension(3)  , intent(in) :: st_adims
!
integer :: j,m
integer :: ndel
!
if(st_mod_sta .and. st_mod_atom/=ST_ATOM_OFF) then
   loop_atom_mod: do j = iatom, cr_natoms      ! Loop over all atoms in current layer
      loop_dim: do m=1,3                       ! xyz dimensions
         if(st_idims(m)==1) cycle loop_dim     ! Single unit cell, skip
         if(cr_pos(m,j) < st_dims(m,1)) then   ! Atom is too low shift upward
            cr_pos(m,j) = cr_pos(m,j) + st_adims(m)
            if(cr_pos(m,j)>=st_dims(m,2)) then ! Atom is outside flag to delete
               cr_iscat(1,j) = -9
               cycle loop_atom_mod
            endif
         elseif(cr_pos(m,j) >= st_dims(m,2) .and. st_idims(m)>0) then  ! Atom is too high
            cr_pos(m,j) = cr_pos(m,j) - st_adims(m)
            if(cr_pos(m,j) <st_dims(m,1)) then ! Atom is outside flag to delete
               cr_iscat(1,j) = -9
               cycle loop_atom_mod
            endif
         endif
      enddo loop_dim
   enddo loop_atom_mod
!
   ndel = 0
   j = iatom
!
   loop_atom_rem: do                          ! Loop to remove atoms
      loop_m:do                               ! Find next atoms with iscat == -9
         if(cr_iscat(1,j+ndel)==-9) then        ! Found a candidate
            ndel = ndel + 1                   ! Need to remove another one
            if(j+ndel==cr_natoms) then        ! Last atom to delete is last in crystal
               if(cr_iscat(1,j+ndel)==-9) then  !   and needs to be removed
                  ndel = ndel + 1
                  exit loop_atom_rem
               else                           !   but is a good atom
                  exit loop_m
               endif
            endif
         else
            exit loop_m
         endif
      enddo loop_m
!
      if(ndel>0) then                         ! We need to remove atoms within the list
         cr_pos  (:,j   ) = cr_pos  (:,j+ndel)
         cr_iscat(1,j   ) = cr_iscat(1,j+ndel)
         cr_mole (  j   ) = cr_mole (  j+ndel)
         cr_surf (:,j   ) = cr_surf (:,j+ndel)
         cr_magn (:,j   ) = cr_magn (:,j+ndel)
         cr_prop (  j   ) = cr_prop (  j+ndel)
      endif
      j = j + 1                               ! Test next atom
      if(j     ==cr_natoms) then              ! We are at the crystals's end
         if(cr_iscat(1,j)==-9) ndel = ndel + 1  !   But this one needs to be removed as well
         exit loop_atom_rem
      endif
   enddo loop_atom_rem
   cr_natoms = cr_natoms - ndel               ! Update actual atom number
endif
!
end subroutine  stack_atom_modulo
!
!*****7*****************************************************************
!
      SUBROUTINE stack_rot_setup (line, lbef, iatom, i, lorigin) 
!-                                                                      
!     performs the setup for the rotations around the different axes    
!+                                                                      
USE discus_allocate_appl_mod
USE discus_config_mod 
USE crystal_mod 
USE symm_mod 
USE stack_mod 
USE str_comp_mod
!                                                                       
IMPLICIT none 
!                                                                       
CHARACTER(len= * ) :: line 
INTEGER :: lbef 
INTEGER :: iatom 
INTEGER         :: new_nscat = 1! Dummy for allocation 
INTEGER         :: new_nsite = 1! Dummy for allocation 
INTEGER :: i, j 
LOGICAL :: lorigin 
!
IF ( SYM_MAXSCAT < MAXSCAT ) THEN
   new_nscat = MAXSCAT
   new_nsite = MAX(cr_ncatoms, MAXSCAT, SYM_MAXSITE)
   CALL alloc_symmetry(new_nscat, new_nsite)
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
      sym_hkl = matmul(real(cr_gten,kind=PREC_DP), sym_uvw)
   ELSE 
      sym_hkl (1) = st_rot_no (1) 
      sym_hkl (2) = st_rot_no (2) 
      sym_hkl (3) = st_rot_no (3) 
      sym_uvw = matmul(real(cr_rten,kind=PREC_DP), sym_hkl)
   ENDIF 
ELSEIF (str_comp (line, 'mod1', 4, lbef, 4) ) then 
   sym_angle = st_rot_ang_m1 (i) 
   IF (st_rot_m1_lspace) then 
      sym_uvw (1) = st_rot_m1 (1) 
      sym_uvw (2) = st_rot_m1 (2) 
      sym_uvw (3) = st_rot_m1 (3) 
      sym_hkl = matmul(real(cr_gten,kind=PREC_DP), sym_uvw)
   ELSE 
      sym_hkl (1) = st_rot_m1 (1) 
      sym_hkl (2) = st_rot_m1 (2) 
      sym_hkl (3) = st_rot_m1 (3) 
      sym_uvw = matmul(real(cr_rten,kind=PREC_DP), sym_hkl)
   ENDIF 
ELSEIF (str_comp (line, 'mod2', 4, lbef, 4) ) then 
   sym_angle = st_rot_ang_m2 (i) 
   IF (st_rot_m2_lspace) then 
      sym_uvw (1) = st_rot_m2 (1) 
      sym_uvw (2) = st_rot_m2 (2) 
      sym_uvw (3) = st_rot_m2 (3) 
      sym_hkl = matmul(real(cr_gten,kind=PREC_DP), sym_uvw)
         ELSE 
      sym_hkl (1) = st_rot_m2 (1) 
      sym_hkl (2) = st_rot_m2 (2) 
      sym_hkl (3) = st_rot_m2 (3) 
      sym_uvw = matmul(real(cr_rten,kind=PREC_DP), sym_hkl)
   ENDIF 
ENDIF 
!                                                                       
END SUBROUTINE stack_rot_setup                
!
!*****7*****************************************************************
!
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
USE fourier_lmn_mod 
USE four_strucf_mod
USE molecule_mod 
USE discus_save_mod 
USE stack_mod 
USE phases_stack_mod
USE powder_mod 
USE powder_tables_mod 
USE read_internal_mod
USE structur
USE spcgr_apply
!
USE errlist_mod 
USE prompt_mod 
USE wink_mod 
USE support_mod
use precision_mod
!
IMPLICIT none 
!                                                                       
LOGICAL, INTENT(IN) :: calc_f2aver
!                                                                       
integer, parameter :: MAXMASK = 4
!
      INTEGER :: i, j, k, l
      INTEGER :: iscat 
      INTEGER :: lbeg (3)
      INTEGER :: ncell 
      INTEGER         :: n_layers ! Number of layers for current layer type
      INTEGER, dimension(3)         :: n_qxy    = 1 ! Number of data points in reciprocal space
      INTEGER         :: n_nscat  = 1 ! Number of different atom types
      INTEGER         :: n_atoms  = 1 ! Number of atoms
integer     :: nanis ! Number of different ADP's
      INTEGER         :: ncatoms  = 1 ! Number of atoms
      INTEGER         :: n_mole  ! number of molecules in input file
      INTEGER         :: n_type  ! number of molecule types in input file
      INTEGER         :: n_atom  ! number of molecule atoms in input file
integer, dimension(3) :: fnum   ! Number of increments (reduced by Friedel)
logical,dimension(0:MAXMASK) :: uni_mask
!integer, dimension(3) :: n_cells
logical :: four_is_new  ! The reciprocal space dimensions have changed
!
      LOGICAL :: lout 
!                                                                       
      REAL(kind=PREC_DP) :: ss 
!
!     n_qxy   = 1
!     n_nscat = 1
uni_mask(0)   = .true.
uni_mask(1:3) = .true.
uni_mask(4)   = .false.
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
!------ preset some values                                              
!                                                                       
ilots = LOT_OFF 
ncell = 0 
lout = .false. 
ss = seknds (0.0) 
CALL four_layer (four_is_new)
CALL fourier_lmn(eck,vi,inc,lmn,off_shift)
fnum = num    ! Prior to TEST_FRIEDEL WORK
!
! ALLOCATE the fourier part of stacking faults
!          To ensure compatibility with four_conv, csf, st_csf must be of 
!          the exact size and shape of num
!
if(any(num/=ubound(csf))) then
   n_qxy = num
   call alloc_diffuse_four (n_qxy )
   if(ier_num/=0) return
endif
if(any(num/=ubound(st_csf))) then
   n_qxy = num
   call alloc_stack_four (n_qxy )
   if(ier_num/=0) return
endif
if(cr_nscat/=ubound(cfact,2)) then
   call alloc_diffuse_scat(cr_nscat)
   if(ier_num/=0) return
endif
if(cr_natoms/=ubound(xat,1)) then
   call alloc_diffuse_atom(cr_natoms)
   if(ier_num/=0) return
endif
!if(ubound(csf,1)/=ubound(st_csf,1) .or. ubound(csf,1)/=ubound(st_csf,1) .or.  &
!   ubound(csf,1)/=ubound(st_csf,1) .or.                                       &
!   num(1)/=ST_MAXQXY(1) .or. num(2)/=ST_MAXQXY(2) .or. num(3)/=ST_MAXQXY(3) .or. &
!   num(1)/=   MAXQXY(1) .or. num(2)/=   MAXQXY(2) .or. num(3)/=   MAXQXY(3) .or. &
!   cr_nscat>DIF_MAXSCAT           .OR.                                        &
!          MAX(cr_natoms,st_nlayer)>DIF_MAXAT                                   ) then
!
!  n_qxy = num
!  CALL alloc_stack_four (n_qxy )
!  n_nscat = MAX(n_nscat,cr_nscat,DIF_MAXSCAT)
!   n_atoms = MAX(cr_natoms,st_nlayer,DIF_MAXAT)
!  CALL alloc_diffuse_four (n_qxy)
!  CALL alloc_diffuse_scat (n_nscat)
!   CALL alloc_diffuse_atom (n_atoms)
!ENDIF
!                                                                       
!     Now start calculation if sufficient space                         
!                                                                       
!     IF (num(1)*num(2)*num(3) .gt. MAXQXY  .OR.          &
!if(num(1)>ST_MAXQXY(1) .or. num(2)>ST_MAXQXY(2) .or. num(3)>ST_MAXQXY(3) .or. &
!          cr_nscat>DIF_MAXSCAT           .OR.          &
!          MAX(cr_natoms,st_nlayer)>DIF_MAXAT   ) THEN
!        n_qxy   = MAX(n_qxy,num,MAXQXY)
!        n_nscat = MAX(n_nscat,cr_nscat,DIF_MAXSCAT)
!        n_atoms = MAX(cr_natoms,st_nlayer,DIF_MAXAT)
!        CALL alloc_diffuse (n_qxy, cr_nscat, n_atoms)
!        IF (ier_num.ne.0) THEN
!          RETURN
!        ENDIF
!ENDIF
!     IF (num(1)*num(2)*num(3) .le.ST_MAXQXY) then 
if(num(1) >ST_MAXQXY(1) .or.  num(2) >ST_MAXQXY(2) .or.  num(3) >ST_MAXQXY(3)) then
   ier_num = - 8 
   ier_typ = ER_APPL 
   return
ENDIF 
!                                                                       
!------ --zero some arrays                                              
!                                                                       
csf(1:num(1),1:num(2), 1:num(3)) = cmplx(0.0D0, 0.0D0,KIND=KIND(0.0D0))
!
!     --read first layer to ensure that the metric tensors are set      
!                                                                       
         CALL rese_cr 
         IF(st_layer(1)(1:8)=='internal') THEN
            CALL testfile_internal (st_layer(1), n_atoms, ncatoms, &        ! Get size of internal structure
              n_nscat, nanis, n_mole, n_type, n_atom)
         ELSE
            ier_num = -182
            ier_typ = ER_APPL
            ier_msg(1) = 'Stacking faults should have been internal'
            ier_msg(2) = 'Stacking fault layer file error 05'
            return
         ENDIF
         IF(n_atoms > NMAX .or. n_nscat > MAXSCAT .or. st_nlayer > NMAX) THEN
            n_atoms = MAX( n_atoms, st_nlayer, NMAX )
            n_nscat = MAX( n_nscat, MAXSCAT)
            CALL alloc_crystal_scat (n_nscat)
            CALL alloc_crystal_nmax (n_atoms)
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
         IF(st_internal(st_type(1)) ) THEN
            CALL readstru_internal(MAXMASK, st_layer (st_type (1) ), uni_mask)!, &
         ELSE
            ier_num = -182
            ier_typ = ER_APPL
            ier_msg(1) = 'Stacking faults should have been internal'
            ier_msg(2) = 'Stacking fault layer file error 06'
            return
!INTER         CALL readstru (NMAX, MAXSCAT, st_layer (1), cr_name, cr_spcgr, cr_set, &
!INTER         cr_a0, cr_win, cr_natoms, cr_nscat, cr_dw, cr_occ, cr_at_lis, cr_pos,  &
!INTER         cr_mole, cr_surf, cr_magn,                                     &
!INTER         cr_iscat, cr_prop, cr_dim, cr_magnetic, as_natoms, as_at_lis, as_dw, as_pos,&
!INTER         as_iscat, as_prop, sav_ncell, sav_r_ncell, sav_ncatoms,        &
!INTER         spcgr_ianz, spcgr_para)                                        
         ENDIF
         CALL setup_lattice (cr_a0, cr_ar, cr_eps, cr_gten, cr_reps,    &
         cr_rten, cr_win, cr_wrez, cr_v, cr_vr, lout, cr_gmat, cr_fmat, &
         cr_cartesian,                                                  &
              cr_tran_g, cr_tran_gi, cr_tran_f, cr_tran_fi)
!                                                                       
!------ --preset some tables, calculate average structure               
!                                                                       
         CALL four_cexpt 
         CALL four_stltab 
         IF (ier_num.ne.0) return 
!                                                                       
!     --Loop over all layer chemistries                                 
!                                                                       
         layers: DO l = 1, st_nchem 
         nxat = 0 
!                                                                       
!     ----loop over all layers                                          
!                                                                       
DO j = 1, st_nlayer 
   if(st_nchem/=ubound(cfact,2)) then
      call alloc_diffuse_scat(st_nchem)
      if(ier_num/=0) return
   endif
   if(st_nlayer/=ubound(xat,1)) then
      call alloc_diffuse_atom(st_nlayer)
      if(ier_num/=0) return
   endif
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
            CALL four_strucf (0, .false., .false., .false., 1, 1, fnum) 
!                                                                       
!     ------copy structure factor to temporary place                    
!                                                                       
            !DO i = 1, num(1)*num(2)*num(3)
            !   acsf (i) = tcsf (i)                                                          ! Neder's original code
            !ENDDO 
!
            DO i = 1, num (1)
               DO j = 1, num (2)
                  DO k = 1, num (3)
                     acsf (i,j,k) = tcsf (i,j,k)                                                          ! My code
                  ENDDO
               ENDDO
            ENDDO
!
            n_layers = nxat
!                                                                       
!     ----read corresponding layer                                      
!                                                                       
            CALL rese_cr 
            IF(st_layer(l)(1:8)=='internal') THEN
               CALL testfile_internal (st_layer(l), n_atoms, ncatoms, &        ! Get size of internal structure
                 n_nscat, nanis, n_mole, n_type, n_atom)
            ELSE
               ier_num = -182
               ier_typ = ER_APPL
               ier_msg(1) = 'Stacking faults should have been internal'
               ier_msg(2) = 'Stacking fault layer file error 07'
               return
!INTER               CALL test_file ( st_layer(l), n_atoms, n_nscat, n_mole, n_type,&
!INTER                                n_atom, n_cells, -1, .false. )
            ENDIF
!            CALL test_file ( st_layer_c(l), n_atoms, n_nscat,n_mole,    &
!                             n_type, n_atom, -1, .true. )
            IF(n_atoms > NMAX .or. n_nscat > MAXSCAT) THEN
               n_atoms = MAX( n_atoms, NMAX)
               n_nscat = MAX( n_nscat, MAXSCAT)
               CALL alloc_crystal_scat (n_nscat)
               CALL alloc_crystal_nmax (n_atoms)
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
            IF(st_layer_c(l)(1:8)=='internal') THEN
               CALL readstru_internal(MAXMASK, st_layer_c(l), uni_mask)
            ELSE
               ier_num = -182
               ier_typ = ER_APPL
               ier_msg(1) = 'Stacking faults should have been internal'
               ier_msg(2) = 'Stacking fault layer file error 08'
               return
!INTER            CALL readstru (NMAX, MAXSCAT, st_layer_c (l), cr_name,      &
!INTER            cr_spcgr, cr_set, cr_a0, cr_win, cr_natoms, cr_nscat, cr_dw, cr_occ,       &
!INTER            cr_at_lis, cr_pos, cr_mole, cr_surf, cr_magn, cr_iscat, cr_prop, cr_dim, cr_magnetic, as_natoms,    &
!INTER            as_at_lis, as_dw, as_pos, as_iscat, as_prop, sav_ncell,     &
!INTER            sav_r_ncell, sav_ncatoms, spcgr_ianz, spcgr_para)           
            ENDIF
            IF (ier_num.ne.0) then 
               ier_msg (1) = 'Reading layer file: ' 
               ier_msg (2) = trim(st_layer_c (l))
               RETURN 
            ENDIF 
            CALL setup_lattice (cr_a0, cr_ar, cr_eps, cr_gten, cr_reps, &
            cr_rten, cr_win, cr_wrez, cr_v, cr_vr, lout, cr_gmat,       &
            cr_fmat, cr_cartesian,                                      &
            cr_tran_g, cr_tran_gi, cr_tran_f, cr_tran_fi)
            ier_num = 0 
            ier_typ = ER_NONE 
!
!           If necessary allocate diffuse
!
            if(cr_nscat/=ubound(cfact,2)) then
               call alloc_diffuse_scat(cr_nscat)
               if(ier_num/=0) return
            endif
            if(cr_natoms/=ubound(xat,1)) then
               call alloc_diffuse_atom(cr_natoms)
               if(ier_num/=0) return
            endif
!           IF ( cr_nscat>DIF_MAXSCAT           .OR.          &
!                MAX(cr_natoms,st_nlayer)>DIF_MAXAT   ) THEN
!              n_qxy   = MAX(n_qxy,num,MAXQXY)
!              n_nscat = MAX(n_nscat,cr_nscat,DIF_MAXSCAT)
!              n_atoms = MAX(cr_natoms,st_nlayer,DIF_MAXAT)
!              CALL alloc_diffuse_four (n_qxy)                 ! Needs to be done only once at start
!              CALL alloc_diffuse_scat (n_nscat)
!              CALL alloc_diffuse_atom (n_atoms)
!              IF (ier_num.ne.0) THEN
!                RETURN
!              ENDIF
!           ENDIF
            CALL dlink (ano, lambda, rlambda, renergy, l_energy, &
                        diff_radiation, diff_table, diff_power) 
                                                                        
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
            !DO i = 1, num (1) * num (2) *num(3)
            !   st_csf (i) = cmplx (0.0D0, 0.0D0,KIND=KIND(0.0D0))                              ! Neder's original code
            !ENDDO 
!                                                                       
            DO i = 1, num (1)
               DO j = 1, num (2)
                  DO k = 1, num (3)
                     st_csf (i,j,k) = cmplx (0.0D0, 0.0D0,KIND=KIND(0.0D0))                              ! My code, Apparently these st_csf arrays are to be treated 3D as well.
                  ENDDO
               ENDDO
            ENDDO
!
!------ ------loop over all different atom types                        
!                                                                       
         DO iscat = 1, cr_nscat 
            CALL four_getatm (iscat, ilots, lbeg, ncell) 
            CALL four_strucf (iscat, .true., .false., .false., 1, 1, fnum) 
!                                                                       
!------ --------Add this part of the structur factor to the total       
!                                                                       
            !DO i = 1, num (1) * num (2) *num(3)
            !   st_csf (i) = st_csf (i) + tcsf (i)                                                       ! Neder's original code 
            !ENDDO 
!
            DO i = 1, num (1)
               DO j = 1, num (2)
                  DO k = 1, num (3)
                     st_csf (i,j,k) = st_csf (i,j,k) + tcsf (i,j,k)                                                       ! My code
                  ENDDO
               ENDDO
            ENDDO
!
            IF (four_log) then 
              WRITE (output_io, 3000) cr_at_lis (iscat), nxat 
            ENDIF 
!
!------ -------Calculate form factor squared
!
!           IF(calc_f2aver) THEN
!              signum = 1.0D0
!              IF(REAL(cfact_pure(1, iscat))< 0.0D0) signum = -1.0D0
! Really only needed for <f^2> and <f>^2 for F(Q) and S(Q) during first calculation
!              xstart = pow_qmin  /REAL(zpi)
!              xdelta = pow_deltaq/REAL(zpi)
!              CALL powder_stltab(pow_npkt,xstart,xdelta)  
!              DO i = 1, pow_npkt         ! Always over all points in powder pattern!
!                 pow_f2aver (i) = pow_f2aver (i)  + &
!                            real (       cfact_pure(powder_istl(i), iscat)  * &
!                                  conjg (cfact_pure(powder_istl(i), iscat)))  &
!                            * n_layers * nxat
!                 pow_faver2 (i) = pow_faver2 (i) +  &
!                       SQRT(real (       cfact_pure(powder_istl(i), iscat)  * &
!                                  conjg (cfact_pure(powder_istl(i), iscat)))) &
!                            * n_layers * nxat *                               &
!                              signum
!              ENDDO
!              pow_nreal = pow_nreal + n_layers * nxat
!              pow_u2aver    = pow_u2aver + cr_dw(iscat) * n_layers * nxat
!           ENDIF
         ENDDO 
         IF(st_new_form .AND. .NOT.diff_lsingle) THEN            ! New Form factros and powder diffraction
            CALL phases_place_stack_form(n_layers, st_nlayer, st_ncunit)    ! Place all form factors of current layer type
         ENDIF
!                                                                       
!     ------Add product of acsf und st_csf to csf                       
!                                                                       
         !DO i = 1, num (1) * num (2) *num(3)
         !   csf (i) = csf (i) + st_csf (i) * acsf (i)                                                 ! Neder's original code
         !ENDDO
!
         DO i = 1, num (1)
            DO j = 1, num (2)
               DO k = 1, num (3)
                  csf (i,j,k) = csf (i,j,k) + st_csf (i,j,k) * acsf (i,j,k)                            ! My code
               ENDDO
            ENDDO
         ENDDO
!
      ENDIF 
   ENDDO  layers
!                                                                       
!     --Calculate average scattering and subtract                       
!                                                                       
   CALL st_fourier_aver (fnum)
   !DO i = 1, num (1) * num (2) *num(3)
   !   csf (i) = csf (i) - acsf (i)                                                                    ! Neders original code
   !ENDDO
!
   DO i = 1, num (1)
      DO j = 1, num (2)
         DO k = 1, num (3)
            csf (i,j,k) = csf (i,j,k) - acsf (i,j,k)                                                   ! My code
         ENDDO
      ENDDO
   ENDDO
!
st_new_form = .FALSE.                            ! Form factors have been written into phases
!                                                                       
!     Compute intensity                                                 
!                                                                       
!DO i = 1, num (1) * num (2) *num(3)
!   dsi (i) = DBLE (csf (i) * conjg (csf (i) ) )                                                      ! Neder's original code
!ENDDO 
!
!ii = 0
!DO i = 1, num (1)
!   DO j = 1, num (2)
!      DO k = 1, num (3)
!         ii = ii + 1
!         dsi (i,j,k)    = DBLE (csf (i,j,k) * conjg (csf (i,j,k) ) )                                     ! My code
!      ENDDO
!   ENDDO
!ENDDO
dsi(1:num(1),1:num(2), 1:num(3)) = DBLE(csf(1:num(1),1:num(2), 1:num(3)) * &
                                  conjg(csf(1:num(1),1:num(2), 1:num(3))))
call four_weight
!
IF (four_log) then 
   CALL four_qinfo 
   ss = seknds (ss) 
   WRITE (output_io, 4000) ss 
ENDIF 
!                                                                       
 3000 FORMAT (  '                   Atom typ = ',A4,7X,'(# ',I9,' )') 
 4000 FORMAT      (/,' Elapsed time    : ',G13.6,' sec') 
!                                                                       
      END SUBROUTINE st_fourier                     
!*****7*****************************************************************
      SUBROUTINE st_fourier_aver(fnum) 
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
      USE read_internal_mod
      USE stack_mod 
      USE structur
      USE spcgr_apply
      USE errlist_mod 
use precision_mod
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
integer, dimension(3), intent(in) :: fnum ! Number increments (reduced by Friedel)       
!                                                                       
integer, parameter :: MAXMASK = 4
      INTEGER :: i, j, k, l 
      INTEGER :: iscat 
      INTEGER :: lbeg (3)
      INTEGER :: ncell 
      LOGICAL :: lout 
      INTEGER, dimension(3)         :: n_qxy    = 1! Number of data points in reciprocal space
      INTEGER         :: n_nscat  = 1! Number of different atom types
      INTEGER         :: n_natoms = 1! Number of atoms 
      INTEGER         :: ncatoms = 1! Number of atoms 
integer     :: nanis ! Number of different ADP's
      INTEGER         :: n_mole  ! number of molecules in input file
      INTEGER         :: n_type  ! number of molecule types in input file
      INTEGER         :: n_atom  ! number of molecule atoms in input file
logical, dimension(0:MAXMASK) :: uni_mask
!integer, dimension(3) :: n_cells
!
      REAL(kind=PREC_DP) :: u (3) 
!
!     n_qxy    = 1
!     n_nscat  = 1
!     n_natoms = 1
uni_mask(0)   = .true.
uni_mask(1:3) = .true.
uni_mask(4)   = .false.
!                                                                       
!------ preset some values                                              
!                                                                       
      ilots = LOT_OFF 
      ncell = 0 
      lout = .false. 
!                                                                       
!     Now start calculation if sufficient space                         
!                                                                       
if(any(num/=ubound(csf))) then
   n_qxy = num
   CALL alloc_diffuse_four (n_qxy )
endif
if(any(num/=ubound(st_csf))) then
   n_qxy = num
   CALL alloc_stack_four (n_qxy )
endif
if(cr_nscat/=ubound(cfact,2)) then
   call alloc_diffuse_scat(cr_nscat)
   if(ier_num/=0) return
endif
if(cr_natoms/=ubound(xat,1)) then
   call alloc_diffuse_atom(cr_natoms)
   if(ier_num/=0) return
endif
!     if(num(1)>   MAXQXY(1) .or. num(2)>   MAXQXY(2) .or. num(3)>   MAXQXY(3) .or. &
!         cr_nscat>DIF_MAXSCAT           .OR.          &
!         MAX(cr_natoms,st_nlayer) > DIF_MAXAT ) THEN
!       n_qxy   = MAX(n_qxy,num,MAXQXY)
!       n_nscat = MAX(n_nscat,cr_nscat,DIF_MAXSCAT)
!       n_natoms = MAX(cr_natoms,st_nlayer,DIF_MAXAT)
!       CALL alloc_diffuse_four (n_qxy)
!       CALL alloc_diffuse_scat (n_nscat)
!       CALL alloc_diffuse_atom (n_natoms)
!       IF (ier_num.ne.0) THEN
!         RETURN
!       ENDIF
!     ENDIF
      if(num(1)<=  MAXQXY(1) .and. num(2)<=  MAXQXY(2) .and. num(3)<=  MAXQXY(3)      ) then
!                                                                       
!------ --zero some arrays                                              
!                                                                       
         !DO i = 1, num (1) * num (2)*num(3) 
!         st_csf(i) = cmplx(0.0d0,0.0d0)                                
         !acsf (i) = cmplx (0.0D0, 0.0D0,KIND=KIND(0.0D0))                                               ! Neder's original code
         !ENDDO 
!                                                                       
         DO i = 1, num (1)
            DO j = 1, num (2)
               DO k = 1, num (3)
                  acsf (i,j,k) = cmplx (0.0D0, 0.0D0,KIND=KIND(0.0D0))                                   ! My code
               ENDDO
            ENDDO
         ENDDO
!
!     --Calculation is only performed if average is needed              
!                                                                       
         IF (st_aver.ne.0.0) then 
!                                                                       
!     ----read first layer to ensure that the metric tensors are set    
!                                                                       
            CALL rese_cr 
            IF(st_layer(1)(1:8)=='internal') THEN
               CALL testfile_internal (st_layer(1), n_natoms, ncatoms, &        ! Get size of internal structure
                 n_nscat, nanis, n_mole, n_type, n_atom)
            ELSE
               ier_num = -182
               ier_typ = ER_APPL
               ier_msg(1) = 'Stacking faults should have been internal'
               ier_msg(2) = 'Stacking fault layer file error 09'
               return
!INTER               CALL test_file ( st_layer(1), n_natoms, n_nscat, n_mole, n_type,&
!INTER                                n_atom, n_cells, -1, .false. )
            ENDIF
!         CALL test_file ( st_layer(1), n_natoms, n_nscat,n_mole, n_type,&
!                          n_atom, -1, .true. )
         IF(n_natoms > NMAX .or. n_nscat > MAXSCAT) THEN
            n_natoms = MAX( n_natoms, NMAX)
            n_nscat  = MAX( n_nscat, MAXSCAT)
            CALL alloc_crystal_scat (n_nscat)
            CALL alloc_crystal_nmax (n_natoms)
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
            IF(st_layer(1)(1:8)=='internal') THEN
               CALL readstru_internal(MAXMASK, st_layer(1), uni_mask)
            ELSE
               ier_num = -182
               ier_typ = ER_APPL
               ier_msg(1) = 'Stacking faults should have been internal'
               ier_msg(2) = 'Stacking fault layer file error 10'
               return
!INTER            CALL readstru (NMAX, MAXSCAT, st_layer (1), cr_name,        &
!INTER            cr_spcgr, cr_set, cr_a0, cr_win, cr_natoms, cr_nscat, cr_dw, cr_occ,&
!INTER            cr_at_lis, cr_pos, cr_mole, cr_surf, cr_magn, cr_iscat, cr_prop, cr_dim, cr_magnetic, as_natoms,    &
!INTER            as_at_lis, as_dw, as_pos, as_iscat, as_prop, sav_ncell,     &
!INTER            sav_r_ncell, sav_ncatoms, spcgr_ianz, spcgr_para)           
            ENDIF
            CALL setup_lattice (cr_a0, cr_ar, cr_eps, cr_gten, cr_reps, &
            cr_rten, cr_win, cr_wrez, cr_v, cr_vr, lout, cr_gmat,       &
            cr_fmat, cr_cartesian,                                      &
              cr_tran_g, cr_tran_gi, cr_tran_f, cr_tran_fi)
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
               st_disp (i, l) = st_disp (i, l) / REAL(st_ndisp (l) ) 
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
            CALL four_strucf (0, .false., .false., .false., 1, 1, fnum) 
!                                                                       
!     ----copy structure factor to temporary place                      
!                                                                       
            !DO i = 1, num (1) * num (2)*num(3) 
            !acsf (i) = tcsf (i)                                                                            ! Neders original code
            !ENDDO 
!                                                                       
            DO i = 1, num (1)
               DO j = 1, num (2)
                  DO k = 1, num (3)
                     acsf (i,j,k) = tcsf (i,j,k)                                                            ! My code
                  ENDDO
               ENDDO
            ENDDO
!
!------ ----zero some arrays                                            
!                                                                       
            !DO i = 1, num (1) * num (2)*num(3) 
            !st_csf (i) = cmplx (0.0D0, 0.0D0,KIND=KIND(0.0D0))                                             ! Neders original code 
            !ENDDO 
!
            DO i = 1, num (1)
               DO j = 1, num (2)
                  DO k = 1, num (3)
                     st_csf (i,j,k) = cmplx (0.0D0, 0.0D0,KIND=KIND(0.0D0))                                 ! My code
                  ENDDO
               ENDDO
            ENDDO
!                                                                       
!     ----Loop over all layer types                                     
!                                                                          
            DO l = 1, st_ntypes                                                                             ! Loop over l
!                                                                       
!     ------read corresponding layer                                    
!                                                                       
            CALL rese_cr 
            IF(st_layer(l)(1:8)=='internal') THEN
               CALL testfile_internal (st_layer(l), n_natoms, ncatoms, &        ! Get size of internal structure
                 n_nscat, nanis, n_mole, n_type, n_atom)
            ELSE
               ier_num = -182
               ier_typ = ER_APPL
               ier_msg(1) = 'Stacking faults should have been internal'
               ier_msg(2) = 'Stacking fault layer file error 11'
               return
!INTER               CALL test_file ( st_layer(l), n_natoms, n_nscat, n_mole, n_type,&
!INTER                                n_atom, n_cells, -1, .false. )
            ENDIF
!        CALL test_file ( st_layer(l), n_natoms, n_nscat,n_mole, n_type,&
!                         n_atom, -1, .true. )
         IF(n_natoms > NMAX .or. n_nscat > MAXSCAT) THEN
            n_natoms = MAX( n_natoms, NMAX )
            n_nscat  = MAX( n_nscat,  MAXSCAT)
            CALL alloc_crystal_scat (n_nscat)
            CALL alloc_crystal_nmax (n_natoms)
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
            IF(st_layer(l)(1:8)=='internal') THEN
               CALL readstru_internal(MAXMASK, st_layer(l), uni_mask)
            ELSE
               ier_num = -182
               ier_typ = ER_APPL
               ier_msg(1) = 'Stacking faults should have been internal'
               ier_msg(2) = 'Stacking fault layer file error 12'
               return
!INTER            CALL readstru (NMAX, MAXSCAT, st_layer (l), cr_name,        &
!INTER            cr_spcgr, cr_set, cr_a0, cr_win, cr_natoms, cr_nscat, cr_dw,cr_occ, &
!INTER            cr_at_lis, cr_pos, cr_mole, cr_surf, cr_magn, cr_iscat, cr_prop, cr_dim, cr_magnetic, as_natoms,    &
!INTER            as_at_lis, as_dw, as_pos, as_iscat, as_prop, sav_ncell,     &
!INTER            sav_r_ncell, sav_ncatoms, spcgr_ianz, spcgr_para)           
            ENDIF
            CALL setup_lattice (cr_a0, cr_ar, cr_eps, cr_gten, cr_reps, &
            cr_rten, cr_win, cr_wrez, cr_v, cr_vr, lout, cr_gmat,       &
            cr_fmat, cr_cartesian,                                      &
              cr_tran_g, cr_tran_gi, cr_tran_f, cr_tran_fi)
!                                                                       
!     ------Add average displacement to each atom                       
!                                                                       
            DO i = 1, cr_natoms                                                                 ! Loop over i
            DO j = 1, 3                                                                         ! Loop over j
            cr_pos (j, i) = cr_pos (j, i) + st_disp (j, l) 
            ENDDO                                                                               ! End of loop over j
            ENDDO                                                                               ! End of loop over i
!                                                                       
            ier_num = 0 
            ier_typ = ER_NONE 
            CALL dlink (ano, lambda, rlambda, renergy, l_energy, &
                        diff_radiation, diff_table, diff_power) 
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
            DO iscat = 1, cr_nscat                                                                       ! Loop over iscat
            CALL four_getatm (iscat, ilots, lbeg, ncell) 
            CALL four_strucf (iscat, .true., .false., .false., 1, 1, fnum) 
!                                                                       
!------ --------Add this part of the structur factor to the total       
!             Wheighted by the relative amount of layers of this type   
!                                                                       
            !DO i = 1, num (1) * num (2)*num(3) 
            !st_csf (i) = st_csf (i) + tcsf (i) * DBLE (st_number (l) ) &                                ! Neder's original code
            !/ DBLE (st_nlayer)                                         
            !ENDDO 
!
            DO i = 1, num (1)
               DO j = 1, num (2)
                  DO k = 1, num (3)
                     st_csf (i,j,k) = st_csf (i,j,k) + tcsf (i,j,k) * DBLE (st_number (l) ) &           ! My code
                     / DBLE (st_nlayer)                                         
                  ENDDO
               ENDDO
            ENDDO
!
            WRITE (output_io, 3000) cr_at_lis (iscat), nxat 
            ENDDO                                                                                       ! End of loop over iscat
            ENDDO                                                                                       ! End of loop over l
!                                                                       
!     ----Save product of acsf und st_csf to acsf                       
!                                                                       
            !DO i = 1, num (1) * num (2)*num(3) 
            !acsf (i) = st_csf (i) * acsf (i)                                                            ! Neder's original code
            !ENDDO
!
            DO i = 1, num (1)
               DO j = 1, num (2)
                  DO k = 1, num (3)
                     acsf (i,j,k) = st_csf (i,j,k) * acsf (i,j,k)                                        ! My code
                  ENDDO
               ENDDO
            ENDDO
!
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
!
!*****7*****************************************************************
!
SUBROUTINE stack_dist_file ()
!-
!     A 'crystal' file is read and the atom types are interpreted as    
!     layer types.                                                      
!     Later on the origins will be interpreted as well as               
!     origins of the stacking faults.                                   
!+                                                                      
USE discus_config_mod 
use discus_allocate_appl_mod
USE read_internal_mod
USE stack_mod 
USE stack_cr_mod
USE structur
USE spcgr_apply
USE errlist_mod 
USE prompt_mod 
!                                                                       
IMPLICIT none 
!                                                                       
integer :: dim_natoms  ! Internal array size
integer :: dim_ncatoms  ! Internal array size
integer :: dim_nscat  ! Internal array size
integer :: dim_nanis  ! Internal array size
integer :: dim_n_mole  ! Internal array size
integer :: dim_n_type  ! Internal array size
integer :: dim_n_atom  ! Internal array size
INTEGER, DIMENSION(3)   :: sav_ncell (3) 
INTEGER                 :: sav_ncatoms 
integer                 :: n_mole
integer                 :: n_type
integer                 :: n_atom
LOGICAL                 :: sav_r_ncell 
integer                 :: ier ! Allocation error
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
IF(st_infile(1:8)=='internal') THEN
   call readstru_size_int(st_infile,  dim_natoms, dim_ncatoms, &
              dim_nscat, dim_nanis, dim_n_mole, dim_n_type, dim_n_atom, &
              st_natoms, sav_ncatoms, st_nscat, st_nanis, n_mole, n_type, n_atom)
!  call readstru_size_int(st_infile, st_natoms, sav_ncatoms, &
!             st_nscat, st_nanis, n_mole, n_type, n_atom)
st_natoms     = 0 
st_nscat      = 0 
sa_natoms     = 0 
   if(allocated(st_is_sym)) deallocate(st_is_sym)
   allocate(st_is_sym(1:sav_ncatoms))
   st_is_sym = 1
!         CALL readstru_internal (st_infile)
   st_nanis=max(1, st_nanis)
   call alloc_anis_generic(st_nanis, st_anis_full, st_prin, ier)
   call stru_readheader_internal (st_infile, st_MAX_SCAT, DIM_NCATOMS, DIM_NANIS, st_name,             &
            st_spcgr, st_spcgr_set, st_set, st_iset,                           &
            st_at_lis, st_nscat, st_dw, st_occ, st_anis, st_is_sym,            &
               st_nanis, st_anis_full, st_prin,                                &
               st_a0, st_win,                                                  &
            sav_ncell, sav_r_ncell, sav_ncatoms, st_spcgr_ianz, st_spcgr_para, &
            st_GEN_ADD_MAX, st_gen_add_n, st_gen_add_power, st_gen_add,        &
            st_SYM_ADD_MAX, st_sym_add_n, st_sym_add_power, st_sym_add )
   call struc_read_atoms_internal(st_infile, st_MMAX, &
            st_natoms, st_pos, st_iscat, st_prop, &
            st_surf, st_magn, st_mole )
ELSE
   ier_num = -182
   ier_typ = ER_APPL
   ier_msg(1) = 'Stacking faults should have been internal'
   ier_msg(2) = 'Stacking fault distribution file error 13'
   return
!INTER      CALL readstru (ST_MMAX, ST_MAX_SCAT, st_infile, st_name, st_spcgr, st_set, &
!INTER      st_a0, st_win, st_natoms, st_nscat, st_dw, st_occ, st_at_lis, st_pos,     &
!INTER      st_mole, st_surf, st_magn,                                        &
!INTER      st_iscat, st_prop, st_dim, st_magnetic, sa_natoms, sa_at_lis, sa_dw, sa_pos,   &
!INTER      sa_iscat, sa_prop, sav_ncell, sav_r_ncell, sav_ncatoms,           &
!INTER      st_spcgr_ianz, st_spcgr_para)                                     
ENDIF
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
!
!*******************************************************************************
!
SUBROUTINE stack_reset
!
USE discus_allocate_appl_mod
USE stack_mod
!
IMPLICIT NONE
!
CALL alloc_stack(1, 1, .TRUE.)
CALL alloc_stack_four( (/1,1,1/) )
CALL alloc_stack_crystal(1, 1)
!
ST_MAXQXY    = 1
ST_MAXLAYER  = 1
ST_MAXTYPE   = 1
!
IF(ALLOCATED(st_layer))    st_layer(:)     = ' '         ! (  ST_MAXTYPE)
IF(ALLOCATED(st_layer_c))  st_layer_c(:)   = ' '         ! (  ST_MAXTYPE)
IF(ALLOCATED(st_llayer))   st_llayer(:)    = 0           ! (  ST_MAXTYPE)
IF(ALLOCATED(st_number))   st_number(:)    = 0           ! (  ST_MAXTYPE)
IF(ALLOCATED(st_ndisp))    st_ndisp(:)     = 0           ! (  ST_MAXTYPE)
IF(ALLOCATED(st_chem))     st_chem(:)      = 0           ! (  ST_MAXTYPE)
IF(ALLOCATED(st_disp))     st_disp(:,:)    = 0           ! (3,ST_MAXTYPE)
IF(ALLOCATED(st_corr))     st_corr(:,:)    = 0.0         ! (  ST_MAXTYPE, ST_MAXTYPE)
IF(ALLOCATED(st_sigma))    st_sigma(:,:,:) = 0.0         ! (  ST_MAXTYPE, ST_MAXTYPE,3)
IF(ALLOCATED(st_trans))    st_trans(:,:,:) = 0.0         ! (  ST_MAXTYPE, ST_MAXTYPE,3)
!
IF(ALLOCATED(st_type))       st_type(:)       = 0        ! (  ST_MAXLAYER)
IF(ALLOCATED(st_internal))   st_internal(:)   = .FALSE.  ! (  ST_MAXTYPE )
IF(ALLOCATED(st_origin))     st_origin(:,:)   = 0.0      ! (3,ST_MAXLAYER)
IF(ALLOCATED(st_rot_ang_no)) st_rot_ang_no(:) = 0.0      ! (  ST_MAXLAYER)
IF(ALLOCATED(st_rot_ang_m1)) st_rot_ang_m1(:) = 0.0      ! (  ST_MAXLAYER)
IF(ALLOCATED(st_rot_ang_m2)) st_rot_ang_m2(:) = 0.0      ! (  ST_MAXLAYER)
!
!IF(ALLOCATED(st_csf))        st_csf(:) = (0.0D0, 0.0D0)          ! (ST_MAXQXY)                                   ! Neders original code
!
IF(ALLOCATED(st_csf))        st_csf(:,:,:) = (0.0D0, 0.0D0)      ! (ST_MAXQXY,ST_MAXQXY,ST_MAXQXY)?????           ! My code, discuss with Neder about the dimensionality of st_csf ...
!
!
st_infile        = ' '
!
st_distr         = ST_DIST_MATRIX
st_infile_l      = 1
st_nlayer        = 0
st_ntypes        = 0
st_nchem         = 0
st_first         = 0
st_ncunit        = 1
st_mod_sta       = .false.
st_tra_aver      = .false.
st_rot_mode      = .false.
st_rot_status    = .false.
st_rot_no_lspace = .true.
st_rot_m1_lspace = .true.
st_rot_m2_lspace = .true.
st_aver          = 0.0
st_prob          = 0.0
st_mod(:,:)      = 0.0
st_inv(:,:)      = 0.0
st_t_aver(:)     = 0.0
st_off(:)        = 0.0
st_sigma_off(:)  = 0.0
st_rot_no(:)     = 0.0
st_rot_m1(:)     = 0.0
st_rot_m2(:)     = 0.0
st_rot_si_no     = 0.0
st_rot_si_m1     = 0.0
st_rot_si_m2     = 0.0
!
END SUBROUTINE stack_reset
!
!*******************************************************************************
!
END MODULE stack_menu
