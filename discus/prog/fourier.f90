MODULE fourier_menu
!
USE errlist_mod 
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE fourier 
!+                                                                      
!     This subroutine 'fourier' calculates the Fourier transform        
!     of the given crystal structure. The algorithm to speed up         
!     the explicite Fourier is based on the program 'DIFFUSE' by        
!     B.D. Butler. See also: B.D. Butler & T.R. Welberry, (1992).       
!     J. Appl. Cryst. 25, 391-399.                                      
!                                                                       
!-                                                                      
USE discus_config_mod 
USE discus_allocate_appl_mod
USE crystal_mod 
USE chem_mod
USE diffuse_mod 
USE external_four
USE fourier_sup
USE fourier_reset_mod
USE get_iscat_mod
USE modify_mod
USE output_mod 
USE discus_show_menu
USE zone
!
USE ber_params_mod
USE build_name_mod
USE calc_expr_mod
USE doact_mod 
USE do_eval_mod
USE do_wait_mod
USE get_params_mod
USE learn_mod 
USE lib_do_operating_mod
USE lib_echo
USE lib_errlist_func
USE lib_help
USE lib_length
USE lib_macro_func
USE class_macro_internal 
USE prompt_mod 
USE do_show_mod
USE precision_mod
USE str_comp_mod
USE string_convert_mod
USE sup_mod
USE take_param_mod
!
IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER, PARAMETER :: MIN_PARA = 21  ! A command requires at least these no of parameters
      INTEGER, PARAMETER :: HKLF4 = 4
      INTEGER, PARAMETER :: CIF   = 1
      INTEGER maxw 
      LOGICAL lold 
      PARAMETER (lold = .false.) 
!                                                                       
      CHARACTER (LEN=PREC_STRING), DIMENSION(MAX(MIN_PARA,MAXSCAT+1))   :: cpara ! (MIN(10,MAXSCAT)) 
      INTEGER             , DIMENSION(MAX(MIN_PARA,MAXSCAT+1))   :: lpara ! (MIN(10,MAXSCAT))
      INTEGER             , DIMENSION(MAX(MIN_PARA,MAXSCAT+1))   :: jj    ! (MAXSCAT) 
      REAL(KIND=PREC_DP)  , DIMENSION(MAX(MIN_PARA,MAXSCAT+1))   :: werte ! (MAXSCAT)
      CHARACTER(5) befehl 
      CHARACTER(LEN=LEN(prompt)) :: orig_prompt
      CHARACTER(LEN=PREC_STRING) :: zeile
      CHARACTER(LEN=PREC_STRING) :: line 
      CHARACTER(LEN=PREC_STRING)  :: infile, calcfile
      CHARACTER(LEN=PREC_STRING)  :: symbol
      INTEGER :: i, j=1, k, ianz, lp, length , lsymbol, iianz
      INTEGER indxg, lbef 
      INTEGER              :: infile_l, outfile_l
      INTEGER              :: n_qxy    ! required size in reciprocal space this run
      INTEGER              :: n_nscat  ! required no of atom types right now
      INTEGER              :: n_natoms ! required no of atoms
      INTEGER              :: four_dim ! Dimension of Fourier that was calculated
      INTEGER              :: istyle   ! Type of hkl file at 'hkl'
      INTEGER, DIMENSION(3):: csize
      LOGICAL              :: ldim 
      LOGICAL              :: ltop = .false. ! the top left corner has been defined
      REAL   , DIMENSION(3)::  divis
      REAL   , DIMENSION(3)::  rhkl
!                                                                       
INTEGER, PARAMETER :: NOPTIONAL = 3
INTEGER, PARAMETER :: O_MODE    = 1
INTEGER, PARAMETER :: O_SYMM    = 2
INTEGER, PARAMETER :: O_TABLE   = 3
CHARACTER(LEN=   5), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=PREC_STRING), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent!opt. para is present
REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 0 ! Number of values to calculate 
!
DATA oname  / 'mode', 'symm'  , 'table'/
DATA loname /  4    ,  4      ,  5     /
opara  =  (/ '0.0000', '0.0000', 'waas  ' /)   ! Always provide fresh default values
lopara =  (/  6      ,  6      ,  4       /)
owerte =  (/  0.0    ,  0.0    ,  0.0     /)
!
!
!
      maxw     = MAX(MIN_PARA,MAXSCAT+1)
      n_qxy    = 1
      n_nscat  = 1
      n_natoms = 1
      orig_prompt = prompt
      prompt = prompt (1:len_str (prompt) ) //'/fourier' 
diff_lsingle = .TRUE.
!                                                                       
   10 CONTINUE 
!                                                                       
      CALL no_error 
      divis (1) = REAL(max (1, inc (1) - 1) ) 
      divis (2) = REAL(max (1, inc (2) - 1) ) 
      divis (3) = REAL(max (1, inc (3) - 1) ) 
!                                                                       
      CALL get_cmd (line, length, befehl, lbef, zeile, lp, prompt) 
      IF (ier_num.eq.0) then 
         IF (line (1:1)  == ' '.or.line (1:1)  == '#' .or.   & 
             line == char(13) .or. line(1:1) == '!'  ) THEN
            IF(linteractive .OR. lmakro) THEN
               GOTO 10
            ELSE
               RETURN
            ENDIF
         ENDIF
!                                                                       
!     search for "="                                                    
!                                                                       
indxg = index (line, '=') 
IF (indxg.ne.0.AND..NOT. (str_comp (befehl, 'echo', 2, lbef, 4) ) &
              .AND..NOT. (str_comp (befehl, 'syst', 2, lbef, 4) )    &
              .AND..NOT. (str_comp (befehl, 'help', 2, lbef, 4) .OR. &
                          str_comp (befehl, '?   ', 2, lbef, 4) )    &
              .AND. INDEX(line,'==') == 0                            ) THEN
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
                  line(1:length-1) = line(2:length)
                  line(length:length) = ' '
                  length = length - 1
                  CALL file_kdo(line, length)
               ELSE 
                  ier_num = - 13 
                  ier_typ = ER_MAC 
               ENDIF 
!                                                                       
!     reset all fourier settings 'reset'
!                                                                       
            ELSEIF (str_comp (befehl, 'rese', 2, lbef, 4) ) then 
               CALL fourier_reset
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
               CALL do_eval (zeile, lp, .TRUE.) 
!                                                                       
!     Terminate Fourier 'exit'                                          
!                                                                       
            ELSEIF (str_comp (befehl, 'exit', 3, lbef, 4) ) then 
               GOTO 9999 
!                                                                       
!     help 'help' , '?'                                                 
!                                                                       
            ELSEIF (str_comp (befehl, 'help', 2, lbef, 4) .OR.  &
                    str_comp (befehl, '?   ', 1, lbef, 4) ) THEN                                      
               IF (str_comp (zeile, 'errors', 2, lp, 6) ) THEN 
                  lp = lp + 7 
                  CALL do_hel ('discus '//zeile, lp) 
               ELSE 
                  lp = lp + 12 
                  CALL do_hel ('discus four '//zeile, lp) 
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
!     All Fourier commands
!
            ELSE
!                                                                       
!     Define the ascissa 'absc'                                         
!                                                                       
            IF (str_comp (befehl, 'absc', 1, lbef, 4) ) then 
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
                        IF (inc(1) * inc(2) * inc(3) .gt. MAXQXY  .OR.   &
                            cr_natoms > DIF_MAXAT                 .OR.   &
                            cr_nscat>DIF_MAXSCAT              ) THEN
                          n_qxy    = MAX(n_qxy,inc(1) * inc(2)*inc(3),MAXQXY)
                          n_natoms = MAX(n_natoms,cr_natoms,DIF_MAXAT)
                          n_nscat  = MAX(n_nscat,cr_nscat,DIF_MAXSCAT)
                          call alloc_diffuse (n_qxy, n_nscat, n_natoms)
                          IF (ier_num.ne.0) THEN
                            RETURN
                          ENDIF
                        ENDIF
                        CALL dlink (ano, lambda, rlambda, renergy, l_energy, &
                                    diff_radiation, diff_table, diff_power) 
                        IF(ier_num==0) CALL calc_000 (rhkl) 
                     ENDIF 
                  ELSEIF (ianz.eq.0) then 
                     rhkl (1) = 0.0 
                     rhkl (2) = 0.0 
                     rhkl (3) = 0.0 
                     CALL dlink (ano, lambda, rlambda, renergy, l_energy,    &
                                    diff_radiation, diff_table, diff_power) 
                     IF(ier_num==0) CALL calc_000 (rhkl) 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ENDIF 
!                                                                       
!     define the anomalous scattering curve for an element 'delf'       
!                                                                       
            ELSEIF (str_comp (befehl, 'delf', 2, lbef, 4) ) then 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ianz.eq.3) then 
                  i = 1 
                  CALL get_iscat (i, cpara, lpara, werte, maxw, lold) 
                  IF (ier_num.eq.0) then 
                     IF(werte(1) .gt. 0.0D0) then 
                        DO k = 1, i 
                        jj (k) = nint (werte (1) ) 
                        ENDDO 
                        cpara (1) = '0.0' 
                        lpara (1) = 3 
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxw)                                           
                        IF (ier_num.eq.0) then 
                           DO k = 1, i 
                           cr_delfr_u (jj (k) ) = werte (2) 
                           cr_delfi_u (jj (k) ) = werte (3) 
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
                     i = 1 
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
!     set the energy of the radiation to be used 'energy'                             
!                                                                       
            ELSEIF (str_comp (befehl, 'energy', 2, lbef, 6) ) then 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ianz.eq.1) then 
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  renergy  = werte (1) 
                  lambda   = ' ' 
                  l_energy = .true.
               ELSE 
                  ier_num = -6 
                  ier_typ = ER_COMM 
               ENDIF 
!                                                                       
!     switch to electron diffraction 'electron'                                
!                                                                       
            ELSEIF (str_comp (befehl, 'electron', 2, lbef, 8) ) then 
               lxray = .true. 
               diff_radiation = RAD_ELEC
               diff_table     = RAD_INTER
               lambda = ' '
!                                                                       
!     calculate at a SHELXL list of reciprocal points 'hkl'                     
!                                                                       
            ELSEIF (str_comp (befehl, 'hkl', 2, lbef, 3) ) then 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ier_num.eq.0) then 
                  IF (ianz >= 4) then 
                     CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1)
                     infile   = cpara(1)
                     infile_l = lpara(1)
                     CALL del_params (1, ianz, cpara, lpara, maxw) 
                     CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1)
                     calcfile   = cpara(1)
                     outfile_l = lpara(1)
                     CALL del_params (1, ianz, cpara, lpara, maxw) 
                     IF(ianz == 2 ) THEN   ! Style is given
                        IF(str_comp (cpara(ianz), 'hklf4', 5, lpara(ianz), 5)) THEN
                           istyle = HKLF4
                           ianz   = ianz - 1
                           CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                        ELSEIF(str_comp (cpara(ianz), 'cif', 3, lpara(ianz), 3)) THEN
                           istyle = CIF
                           ianz   = ianz - 1
                           CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                        ELSE
                           CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                           istyle = NINT(werte(2))
                        ENDIF
                     ELSE
                        istyle = HKLF4
                        CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     ENDIF
                     CALL calc_hkl(infile,infile_l, calcfile, outfile_l, werte(1),istyle )
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
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
                        eck(1,4) = eck(1,1)       ! This is a layer
                        eck(2,4) = eck(2,1)       ! Set verticval corner
                        eck(3,4) = eck(3,1)       ! to lower left values
                        inc (1) = nint (werte (10) ) 
                        inc (2) = nint (werte (11) ) 
                        inc (3) = 1             ! No increment along vertical axis
                        divis (1) = REAL(max (1, inc (1) - 1) ) 
                        divis (2) = REAL(max (1, inc (2) - 1) ) 
                        divis (3) =             1
                        DO i = 1, 3 
                        vi (i, 1) = (eck (i, 2) - eck (i, 1) ) / divis (1)
                        vi (i, 2) = (eck (i, 3) - eck (i, 1) ) / divis (2)
                        vi (i, 3) = (eck (i, 4) - eck (i, 1) ) / divis (3)
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
      ELSEIF (str_comp (befehl, 'll  ', 2, lbef, 4) .or. &
              str_comp (befehl, 'lr  ', 2, lbef, 4) .or. &
              str_comp (befehl, 'ul  ', 2, lbef, 4) .or. &
              str_comp (befehl, 'tl  ', 2, lbef, 4) ) then                                                              
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ier_num.eq.0) then 
                  IF (ianz.eq.3) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        IF (str_comp (befehl, 'll  ', 2, lbef, 4) ) j = 1 
                        IF (str_comp (befehl, 'lr  ', 2, lbef, 4) ) j = 2 
                        IF (str_comp (befehl, 'ul  ', 2, lbef, 4) ) j = 3 
                        IF (str_comp (befehl, 'tl  ', 2, lbef, 4) ) j = 4 
                        IF (str_comp (befehl, 'tl  ', 2, lbef, 4) ) ltop = .true.
                        DO i = 1, 3 
                           eck (i, j) = werte (i) 
                        ENDDO 
                        DO i = 1, 3 
                           vi (i, 1) = (eck (i, 2) - eck (i, 1) ) / divis ( 1)                                              
                           vi (i, 2) = (eck (i, 3) - eck (i, 1) ) / divis ( 2)                                              
                           vi (i, 3) = (eck (i, 4) - eck (i, 1) ) / divis ( 3)                                              
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
                        lot_all = str_comp(cpara(4), 'all', 3, lpara(4), 3)
                        iianz = 3
                        CALL ber_params(iianz, cpara, lpara, werte, maxw)
                        IF (ier_num.eq.0) then 
                           ldim = .true. 
                           DO i = 1, 3 
                              ldim = ldim.and. (0..lt.nint (werte (i) )    &
                                         .and. nint (werte (i) ) .le.cr_icc (i) )      
                           ENDDO 
                           IF (ldim) then 
                              ls_xyz (1) = nint (werte (1) ) 
                              ls_xyz (2) = nint (werte (2) ) 
                              ls_xyz (3) = nint (werte (3) ) 
!
                              CALL do_cap (cpara (5) )
                              lperiod = (cpara (5) (1:1) .eq.'Y')
                              CALL four_csize (cr_icc, csize, lperiod, ls_xyz)
                              iianz = 1
                              IF(lot_all) THEN
                                 WRITE(cpara(1),'(I12)')  csize(1)*csize(2)*csize(3)
                                 lpara(1) = 12
                              ELSE
                                 cpara(1) = cpara(4)
                                 lpara(1) = lpara(4)
                              ENDIF
                              CALL ber_params(iianz, cpara, lpara, werte, maxw)
                              nlots = nint(werte(1) ) 
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
                        IF (werte (1) .gt.0.0D0) then 
                           inc (1) = nint (werte (1) ) 
                           divis (1) = REAL(max (1, inc (1) - 1) ) 
                           DO i = 1, 3 
                           vi (i, 1) = (eck (i, 2) - eck (i, 1) ) / divis (1)                                  
                           vi (i, 2) = (eck (i, 3) - eck (i, 1) ) / divis (2)                                  
                           vi (i, 3) = (eck (i, 4) - eck (i, 1) ) / divis (3)                                  
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
               diff_radiation = RAD_NEUT
               diff_table     = RAD_INTER
               lambda = ' '
!                                                                       
!     define the number of points along the ordinate 'no'               
!                                                                       
            ELSEIF (str_comp (befehl, 'nord', 2, lbef, 4) ) then 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ier_num.eq.0) then 
                  IF (ianz.eq.1) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        IF (werte (1) .gt.0.0D0) then 
                           inc (2) = nint (werte (1) ) 
                           divis (2) = REAL(max (1, inc (2) - 1) ) 
                           DO i = 1, 3 
                           vi (i, 1) = (eck (i, 2) - eck (i, 1) ) / divis (1)
                           vi (i, 2) = (eck (i, 3) - eck (i, 1) ) / divis (2)
                           vi (i, 3) = (eck (i, 4) - eck (i, 1) ) / divis (3)                                  
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
!     define the number of points along the vertical 'ntop'               
!                                                                       
            ELSEIF (str_comp (befehl, 'ntop', 2, lbef, 4) ) then 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ier_num.eq.0) then 
                  IF (ianz.eq.1) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        IF (werte (1) .gt.0.0D0) then 
                           inc (3) = nint (werte (1) ) 
                           divis (3) = REAL(max (1, inc (3) - 1) ) 
                           DO i = 1, 3 
                           vi (i, 1) = (eck (i, 2) - eck (i, 1) ) / divis (1)
                           vi (i, 2) = (eck (i, 3) - eck (i, 1) ) / divis (2)
                           vi (i, 3) = (eck (i, 4) - eck (i, 1) ) / divis (3)                                  
                           ENDDO 
                           ltop = .true.
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
            ELSEIF (str_comp (befehl, 'run ', 2, lbef, 4) ) then 
               IF(ilots/=LOT_OFF .AND.  &   ! Lots will be used
                  cr_natoms /= cr_icc(1)*cr_icc(2)*cr_icc(3)*cr_ncatoms) THEN
                  ier_num = -162
                  ier_typ = ER_APPL
                  ier_msg(1) = 'Lots are on but the number of atoms in the'
                  ier_msg(2) = 'crystal differs from the value set by  '
                  ier_msg(3) = '''set crystal''. Check dimensions, purge etc.    '
               ELSEIF(ilots==LOT_OFF .AND.  &   ! Lots will not be used
                  chem_period(1) .AND. chem_period(2) .AND. chem_period(3) .AND. &
                  cr_natoms /= cr_icc(1)*cr_icc(2)*cr_icc(3)*cr_ncatoms) THEN
                  ier_num = -162
                  ier_typ = ER_APPL
                  ier_msg(1) = 'Periodic boundary conditions are on but the'
                  ier_msg(2) = 'number of atoms in the crystal differs from '
                  ier_msg(3) = 'value set by ''set crystal''. Check dimensions.'
               ELSE   
               IF(.not.ltop) THEN           ! The three-D corner was never defined, assume 2D
                  eck(1,4) = eck(1,1)       ! This is a layer
                  eck(2,4) = eck(2,1)       ! Set verticval corner
                  eck(3,4) = eck(3,1)       ! to lower left values
                  vi (1,3) = 0.00
                  vi (2,3) = 0.00
                  vi (3,3) = 0.00
                  inc(3)   = 1
                  divis(3) = 1
               ENDIF
               IF (inc(1) * inc(2) *inc(3) .gt. MAXQXY  .OR.    &
                   cr_natoms > DIF_MAXAT                .OR.    &
                   cr_nscat>DIF_MAXSCAT              ) THEN
                 n_qxy    = MAX(n_qxy,inc(1)*inc(2)*inc(3),MAXQXY)
                 n_natoms = MAX(n_natoms,cr_natoms,DIF_MAXAT)
                 n_nscat  = MAX(n_nscat,cr_nscat,DIF_MAXSCAT)
                 call alloc_diffuse (n_qxy, n_nscat, n_natoms)
                 IF (ier_num.ne.0) THEN
                   RETURN
                 ENDIF
               ENDIF
               IF (inc (1) * inc (2) * inc(3) .le.MAXQXY) then 
                  CALL dlink (ano, lambda, rlambda, renergy, l_energy, &
                              diff_radiation, diff_table, diff_power) 
                 IF (ier_num.ne.0) THEN
                   RETURN
                 ENDIF
                  CALL four_resolution(zeile, lp)
                  IF(l_zone) CALL zone_setup     ! Setup zone axis pattern
                  IF (four_mode.eq.INTERNAL) then 
                     IF (ier_num.eq.0) then 
                        four_log = .true. 
                        CALL four_run 
                     ENDIF 
                  ELSE 
                     four_log = .true. 
                     CALL four_external 
                  ENDIF 
                  IF(l_zone) CALL zone_project   ! Project zone axis pattern
                  four_was_run = .true.
                  ! Specify the fourier type that was calculated
                  four_dim = 0
                  IF(inc(1)>1) four_dim = four_dim + 1
                  IF(inc(2)>1) four_dim = four_dim + 1
                  IF(inc(3)>1) four_dim = four_dim + 1
                  IF(l_zone) THEN
                     four_last = FOUR_ZA   ! Zone axis pattern
                  ELSE
                     four_last = four_dim  ! Fourier N-Dim
                  ENDIF 
                  IF(ilots>1) four_last = -four_last  ! Lots were used
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_APPL 
               ENDIF 
               ENDIF 
!                                                                       
!     define the scattering curve for an element 'scat'                 
!                                                                       
            ELSEIF (str_comp (befehl, 'scat', 2, lbef, 4) ) then 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ianz.eq.10. .or. ianz==12) then 
                  i = 1 
                  CALL get_iscat (i, cpara, lpara, werte, maxw, lold) 
                  IF (ier_num.eq.0) then 
                     IF (werte (1) .gt.0.0D0) then 
                        DO k = 1, i 
                        jj (k) = nint (werte (1) ) 
                        ENDDO 
                        cpara (1) = '0.0' 
                        lpara (1) = 3 
                        CALL ber_params (ianz, cpara, lpara, werte,     &
                        maxw)                                           
                        IF (ier_num.eq.0) then 
                           DO k = 1, i 
                           DO i = 2, ianz-1 
                           cr_scat (i, jj (k) ) = werte (i) 
                           ENDDO 
                           cr_scat (1, jj (k) ) = werte (ianz) 
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
                        IF (werte (1) .gt.0.0D0) then 
                           k = nint (werte (1) ) 
                           cr_scat_int (k) = .true. 
                           cr_scat_equ (k) = .false. 
                        ELSEIF(NINT(werte(1)) == -1) then 
                           DO k = 0, cr_nscat 
                           cr_scat_int (k) = .true. 
                           cr_scat_equ (k) = .false. 
                           ENDDO 
                        ENDIF 
                     ELSE 
                        IF(werte(1) .gt. 0.0D0) then 
                           k = nint (werte (1) ) 
                           cr_scat_equ (k) = .true. 
                           CALL do_cap (cpara (2) ) 
                           cr_at_equ (k) = cpara (2) (1:lpara(2))
                        ELSEIF(NINT(werte(1)) == -1) then 
                           k = nint (werte (1) ) 
                           CALL do_cap (cpara (2) ) 
                           DO k = 1, cr_nscat 
                           cr_scat_equ (k) = .true. 
                           cr_at_equ (k) = cpara (2) (1:lpara(2))
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
           CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                    oname, loname, opara, lopara, lpresent, owerte)
           IF (ier_num.eq.0) then 
              four_symm = .FALSE.
              IF(lpresent(O_SYMM)) THEN       ! set mode: 
                 IF(opara(O_symm)=='apply') THEN
                    four_symm = .TRUE.
                 ELSE
                    four_symm = .FALSE.
                 ENDIF
              ENDIF
              IF(lpresent(O_MODE)) THEN       ! set mode: 
                 IF(opara(O_MODE)=='single') THEN
                    four_accum = 0
!                   IF(four_symm) CALL four_accumulate  ! Call to apply symmetry
                 ELSEIF(opara(O_MODE)=='init') THEN
                    four_accum = -1
                    CALL four_accumulate      ! Call to clear arrays
                    four_accum =  1           ! Toggle to accumulate
                 ELSEIF(opara(O_MODE)=='accumulate') then
                    IF(four_accum==0) THEN    ! First accumulate, initialize
                       four_accum = -1
                       CALL four_accumulate      ! Call to clear arrays
                    ENDIF
                    four_accum =  1           ! Toggle to accumulate
                 ELSEIF(opara(O_MODE)=='finish') THEN
                    four_accum =  2
                    CALL four_accumulate      ! Call to clear arrays
                    four_accum = 0
                 ENDIF
              ELSE                            ! No 'set mode:' parameter present
                 IF(lpresent(O_SYMM)) THEN       ! set mode: 
                    four_accum = 0
!                   IF(four_symm) CALL four_accumulate  ! Call to apply symmetry
                 ENDIF
              ENDIF
                 IF((lpresent(O_MODE) .OR. lpresent(O_SYMM)).AND. ianz==0) THEN
                    CONTINUE
                 ELSE
                 IF(ianz.ge.1.and.ianz.le.2) then 
                    IF(str_comp(cpara(1), 'aver', 1, lpara(1), 4)) THEN
                      IF(ianz.eq.1) THEN 
                         fave = 0.0 
                      ELSE 
                         CALL del_params (1, ianz, cpara, lpara, maxw) 
                         CALL ber_params(ianz, cpara, lpara, werte, maxw)
                         IF (ier_num.eq.0) then 
                             IF(werte(1) .ge.0.0D0 .AND. werte(1).le.100.0D0) THEN                           
                                fave = werte (1) * 0.01 
                             ELSE 
                                ier_num = - 1 
                                ier_typ = ER_FOUR 
                             ENDIF 
!                             IF(ianz==2) THEN
!                                fave_sca = werte(2)
!                             ELSE
!                                fave_sca = 1.0
!                             ENDIF 
                          ENDIF 
                       ENDIF 
                    ELSEIF(str_comp(cpara(1), 'extern', 1, lpara(1), 6) ) THEN
                        four_mode = EXTERNAL 
                    ELSEIF(str_comp(cpara(1), 'intern', 1, lpara(1), 6) ) THEN
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
           ENDIF 
!                                                                       
!     Show the current settings for the Fourier 'show'                  
!                                                                       
            ELSEIF (str_comp (befehl, 'show', 2, lbef, 4) ) THEN 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ier_num == 0) THEN 
                  IF(str_comp(cpara(1), 'scat', 4, lpara(1), 4)) THEN
                     CALL dlink (ano, lambda, rlambda, renergy, l_energy, &
                                 diff_radiation, diff_table, diff_power) 
                     IF (ier_num.ne.0) THEN
                       RETURN
                     ENDIF
                     CALL do_show_scat (ianz, cpara, lpara, werte, maxw)
                  ELSEIF (ianz==  1) THEN 
                     CALL do_show_generic (cpara, lpara, maxw)
                  ELSEIF(ianz==0) THEN
                     IF(.not.ltop) THEN           ! The three-D corner was never defined, assume 2D
                        eck(1,4) = eck(1,1)       ! This is a layer
                        eck(2,4) = eck(2,1)       ! Set verticval corner
                        eck(3,4) = eck(3,1)       ! to lower left values
                        vi (1,3) = 0.00
                        vi (2,3) = 0.00
                        vi (3,3) = 0.00
                        inc(3)   = 1
                        divis(3) = 1
                     ENDIF
                     CALL dlink (ano, lambda, rlambda, renergy, l_energy, &
                                 diff_radiation, diff_table, diff_power) 
                     IF (ier_num.ne.0) THEN
                       RETURN
                     ENDIF
                     IF(l_zone) CALL zone_setup     ! Setup zone axis pattern
                     CALL four_show  ( ltop )
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
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
!     define the ordinate  'top'                                       
!                                                                       
            ELSEIF (str_comp (befehl, 'top', 2, lbef, 3) ) then 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ianz.eq.1) then 
                  IF (cpara (1) .eq.'h') then 
                     extr_top = 1 
                  ELSEIF (cpara (1) .eq.'k') then 
                     extr_top = 2 
                  ELSEIF (cpara (1) .eq.'l') then 
                     extr_top = 3 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ENDIF 
!                                                                       
!     set the wave length to be used 'wvle'                             
!                                                                       
            ELSEIF (str_comp (befehl, 'wvle', 1, lbef, 4) ) then 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ianz.eq.1) then 
                  symbol  = cpara(1)
                  lsymbol = lpara(1)
                  CALL do_cap (symbol) 
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  IF(ier_num == 0) THEN
                     rlambda = werte (1) 
                     lambda = ' ' 
                     l_energy = .false.
                  ELSEIF (ICHAR ('A')  <= ICHAR (symbol    (1:1) ) .AND. &
                          ICHAR (symbol    (1:1) )  <= ICHAR ('Z') ) THEN  
                     lambda = symbol    (1:lsymbol   )
                     l_energy = .false.
                     CALL no_error
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
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
               diff_radiation = RAD_XRAY
               diff_table     = RAD_INTER
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               if(ianz>0) then
               CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                    oname, loname, opara, lopara, lpresent, owerte)
               endif
               if(opara(O_TABLE)=='waas') diff_table= RAD_WAAS
!                                                                       
!     set a zone axis pattern calculation
!                                                                       
            ELSEIF (str_comp (befehl, 'zone', 2, lbef, 4) ) then 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (str_comp (cpara(1), 'OFF', 3, lpara(1), 3) ) then 
                  l_zone = .false.
                  inc(3) =   1
                  ltop   = .true.
               ELSE
                  IF (ianz.eq.4) THEN 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     zone_uvw(1)  = werte (1) 
                     zone_uvw(2)  = werte (2) 
                     zone_uvw(3)  = werte (3) 
                     zone_res     = werte (4) 
                     lxray        = .true.        ! Switch X-ray on 
                     diff_radiation = RAD_ELEC    ! Switch electron diffraction on
                     fave         = 1.0           ! Set 100%aver
                     ilots        = LOT_BOX
                     lperiod      = .true.
                     ls_xyz(1)    = MIN(10,cr_icc(1))
                     ls_xyz(2)    = MIN(10,cr_icc(2))
                     ls_xyz(3)    = MIN(10,cr_icc(3))
                     inc(1)       = 512
                     inc(2)       = 512
                     inc(3)       =   1
                     ltop         = .true.
                     renergy      = 200.00   ! Default to 200 keV
                     l_energy     = .true.
                     l_zone       = .true.
                     lambda       = ' '
                  ELSE 
                     ier_num = -6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ENDIF 
!
!      Wrong command
!
            ELSE 
               ier_num = - 8 
               ier_typ = ER_COMM 
            ENDIF 
         ENDIF 
         ENDIF 
      ENDIF 
!                                                                       
!------ end of command list                                             
!
      IF (ier_num.ne.0) then 
         CALL errlist 
         IF (ier_sta.ne.ER_S_LIVE) then 
            IF (lmakro .OR. lmakro_error) THEN  ! Error within macro or termination errror
               IF(sprompt /= prompt ) THEN
                  ier_num = -10
                  ier_typ = ER_COMM
                  ier_msg(1) = ' Error occured in fourier menu'
                  prompt = orig_prompt
                  prompt_status = PROMPT_ON 
                  RETURN 
               ELSE
                  IF(lmacro_close) THEN
                     CALL macro_close 
                     prompt_status = PROMPT_ON 
                  ENDIF 
               ENDIF 
            ENDIF 
            IF (lblock) then 
               ier_num = - 11 
               ier_typ = ER_COMM 
               prompt = orig_prompt
               prompt_status = PROMPT_ON 
               RETURN 
            ENDIF 
            CALL no_error 
            lmakro_error = .FALSE.
            sprompt = ' '
         ENDIF 
      ENDIF 
      IF(linteractive .OR. lmakro) THEN
         GOTO 10
      ELSE
         RETURN
      ENDIF
 9999 CONTINUE 
      prompt = orig_prompt
!                                                                       
      END SUBROUTINE fourier                        
!
!*****7*****************************************************************
!
SUBROUTINE four_resolution(line, lp)
!-
!  Interpret the optional parameters on the run command to set the resolution
!+
!
USE crystal_mod
USE diffuse_mod
USE get_params_mod
USE matrix_mod
USE metric_mod
USE take_param_mod
USE precision_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: line
INTEGER         , INTENT(INOUT) :: lp
!
REAL(KIND=PREC_DP) :: EPS = 1.0D-7
INTEGER, PARAMETER :: MAXW = 4
!                                                                       
INTEGER, PARAMETER :: NOPTIONAL = 3
INTEGER, PARAMETER :: O_SIGABS  = 1                  ! Current phase number
INTEGER, PARAMETER :: O_SIGORD  = 2                  ! Weight fraction
INTEGER, PARAMETER :: O_SIGTOP  = 3                  ! Single / multiple
CHARACTER(LEN=   6), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=MAX(PREC_STRING,LEN(line))), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent!opt. para is present
REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 0 ! Number of values to calculate 
!
CHARACTER(LEN=MAX(PREC_STRING,LEN(line))), DIMENSION(MAXW) :: cpara
INTEGER            , DIMENSION(MAXW) :: lpara
REAL(KIND=PREC_DP) , DIMENSION(MAXW) :: werte
INTEGER             :: ianz
INTEGER             :: iianz
!
REAL(KIND=PREC_DP), DIMENSION(3,3) :: vimat    ! Col wise increment vectors in layer
REAL(KIND=PREC_DP), DIMENSION(3,3) :: simat    ! Col
REAL(KIND=PREC_DP), DIMENSION(3,3) :: siinv    ! simat^-1
!REAL(KIND=PREC_DP), DIMENSION(3,3) :: trmat    ! trmat = siinv x vimat
REAL(KIND=PREC_DP), DIMENSION(3  ) :: u        ! Sigma along columns of simat
REAL(KIND=PREC_DP), DIMENSION(3)   :: uu       ! Length of user si vectors
REAL(KIND=PREC_DP), DIMENSION(3  ) :: v1       ! Dummy vectors
REAL(KIND=PREC_DP), DIMENSION(3  ) :: v2       ! Dummy vectors
REAL(KIND=PREC_DP), DIMENSION(3  ) :: v3       ! Dummy vectors
!
DATA oname  / 'sigabs', 'sigord', 'sigtop' /
DATA loname /  6      ,  6      ,  6       /
!
diff_res = 0.0D0
CALL get_params (line, ianz, cpara, lpara, maxw, lp)
IF (ier_num.ne.0) RETURN
IF(ianz==0) RETURN        ! No params, use  default NULL matrix
!
CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                              oname, loname, opara, lopara, lpresent, owerte)
!
CALL four_res_optional(lpresent(O_SIGABS), 1, MAXW, opara(O_SIGABS), &
                       lopara(O_SIGABS), werte, iianz)
CALL four_res_optional(lpresent(O_SIGORD), 2, MAXW, opara(O_SIGORD), &
                       lopara(O_SIGORD), werte, iianz)
CALL four_res_optional(lpresent(O_SIGTOP), 3, MAXW, opara(O_SIGTOP), &
                       lopara(O_SIGTOP), werte, iianz)
!
! Build simat
!
!
u = diff_res(2:4,1)
uu(1) = SQRT (skalpro (u, u, cr_rten) )
IF(uu(1)< EPS) THEN
   diff_res(2:4,1) = vi(1:3,1)        ! User vector was NULL, use abscissa
   u = diff_res(2:4,1)
   uu(1) = SQRT (skalpro (u, u, cr_rten) )
ENDIF
!
u = diff_res(2:4,2)
uu(2) = SQRT (skalpro (u, u, cr_rten) )
IF(uu(2)< EPS) THEN
   diff_res(2:4,2) = vi(1:3,2)        ! User vector was NULL, use ordinate
   u = diff_res(2:4,2)
   uu(2) = SQRT (skalpro (u, u, cr_rten) )
ENDIF
!
u = diff_res(2:4,3)
uu(3) = SQRT (skalpro (u, u, cr_rten) )
IF(uu(3)< EPS) THEN
   v1 = diff_res(2:4,1)
   v2 = diff_res(2:4,2)
   CALL vekprod(v1, v2, v3, cr_reps, cr_gten)   ! Build a vecor normal to sigabs and sigord
   uu(3) = SQRT (skalpro (v3, v3, cr_rten) )
   diff_res(2:4,3) = v3
   diff_res(1  ,3) = 0.001
ENDIF
!
simat(1:3, 1) = diff_res(2:4,1) * diff_res(1,1) / uu(1)
!
simat(1:3, 2) = diff_res(2:4,2) * diff_res(1,2) / uu(2)
!
simat(1:3, 3) = diff_res(2:4,3) * diff_res(1,3) / uu(3)
!
CALL matinv3(simat, siinv) 
!
! Build vimat
!
vimat = vi
diff_tr = MATMUL(siinv, vimat)
!
END SUBROUTINE four_resolution
!
!*****7*****************************************************************
!
SUBROUTINE four_res_optional(lpresent, ientry, MAXW, opara, &
                             lopara, werte, ianz)
!
USE diffuse_mod
!
USE get_params_mod
USE take_param_mod
USE precision_mod
!
IMPLICIT NONE
!
LOGICAL                            , INTENT(IN)    :: lpresent ! Optional parameter was present 
INTEGER                            , INTENT(IN)    :: ientry   ! Entry number
INTEGER                            , INTENT(IN)    :: MAXW     ! Dimension of werte
CHARACTER(LEN=*)                   , INTENT(INOUT) :: opara    ! The string with optional values
INTEGER                            , INTENT(INOUT) :: lopara   ! length of string
REAL(KIND=PREC_DP), DIMENSION(MAXW), INTENT(OUT)   :: werte    ! Numerical values
INTEGER                            , INTENT(OUT)   :: ianz    ! Number of numerical values

ianz = 0
IF(lpresent) THEN                ! Sigmais present
   CALL get_optional_multi(MAXW, opara, lopara, werte, ianz)
   IF(ier_num==0) THEN
      IF(ianz>0) THEN
         diff_res(1,ientry) =werte(1)
         IF(ianz==4) THEN
            diff_res(2:4,ientry) = werte(2:4)
         ELSE
            diff_res(2,ientry) = 0.0D0
            diff_res(3,ientry) = 0.0D0
            diff_res(4,ientry) = 0.0D0
         ENDIF
      ELSE
         diff_res(:,ientry) = 0.00D0
      ENDIF
   ELSE
      RETURN
   ENDIF
ELSE
   diff_res(:,ientry) = 0.0D0
ENDIF
!
!
END SUBROUTINE four_res_optional
!
!*****7*****************************************************************
!
      SUBROUTINE four_show ( ltop )
!+                                                                      
!     prints summary of current fourier settings                        
!-                                                                      
      USE discus_config_mod 
      USE diffuse_mod 
      USE four_angles_mod
      USE metric_mod
      USE output_mod 
use precision_mod
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
!
      LOGICAL, INTENT(IN) :: ltop
!                                                                       
      CHARACTER(len=8) :: radiation 
      CHARACTER (LEN=8), DIMENSION(3), PARAMETER :: c_rad = (/ &
         'X-ray   ', 'neutron ', 'electron' /)
      CHARACTER(LEN=1), DIMENSION(0:3)           ::  extr_achs (0:3) 
!     REAL            , DIMENSION(3)             ::  hor
!     REAL            , DIMENSION(3)             ::  ver
!     REAL            , DIMENSION(3)             ::  top
      REAL(kind=PREC_DP)                         ::  angle_vh
      REAL(kind=PREC_DP)                         ::  ratio_vh
      REAL(kind=PREC_DP)                         ::   aver_vh
      REAL(kind=PREC_DP)                         ::  angle_ht
      REAL(kind=PREC_DP)                         ::  ratio_ht
      REAL(kind=PREC_DP)                         ::   aver_ht
      REAL(kind=PREC_DP)                         ::  angle_tv
      REAL(kind=PREC_DP)                         ::  ratio_tv
      REAL(kind=PREC_DP)                         ::   aver_tv
!     REAL            , DIMENSION(3)             ::  zero = (/0.0, 0.0, 0.0/)
      REAL(kind=PREC_DP), DIMENSION(3)           ::  length = (/0.0, 0.0, 0.0/)
      INTEGER i, j 
!     LOGICAL lspace 
!                                                                       
!     REAL do_blen, do_bang 
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
      radiation = 'neutron' 
      IF (lxray) radiation = 'x-ray' 
      radiation = c_rad(diff_radiation)
      IF (lambda.eq.' ') then 
         IF(diff_radiation==2) THEN
            IF(renergy>999.) THEN
               WRITE (output_io, 1202) radiation, rlambda , renergy*0.001
            ELSE 
               WRITE (output_io, 1201) radiation, rlambda , renergy
            ENDIF 
         ELSE 
            WRITE (output_io, 1200) radiation, rlambda , renergy
         ENDIF 
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
      IF(l_zone) THEN
         WRITE( output_io, 1500) zone_uvw(:)
         WRITE( output_io, 1510) zone_res
      ENDIF
!                                                                       
!    !DO i = 1, 3 
!    !u (i) = vi (i, 1) 
!     v (i) = 0.0 
!     w (i) = vi (i, 2) 
!     ENDDO 
!     lspace = .false. 
!     dvi1 = do_blen (lspace, u, zero) 
!     dvi2 = do_blen (lspace, w, zero) 
!     IF (inc (1) .gt.1.and.inc (2) .gt.1) then 
!        dvi3 = do_bang (lspace, u, zero, w) 
!     ELSE 
!        dvi3 = 0.0 
!     ENDIF 
!     IF (dvi3.gt.0) then 
!        dvi4 = dvi2 / dvi1 
!        IF (abs (u (extr_abs) ) .gt.0.0.and.abs (w (extr_ord) )        &
!        .gt.0.0) then                                                  
!           dvi5 = (dvi2 / w (extr_ord) ) / (dvi1 / u (extr_abs) ) 
!        ELSE 
!           ier_num = - 4 
!           ier_typ = ER_FOUR 
!        ENDIF 
!     ELSE 
!        dvi4 = 0.0 
!        dvi5 = 0.0 
!     ENDIF 
      CALL four_angles(ltop, length, angle_vh, ratio_vh, aver_vh, &
                               angle_ht, ratio_ht, aver_ht, &
                               angle_tv, ratio_tv, aver_tv)
!
!     Calculate lengths in Ang-1
!
!      lspace = .false.
!      hor(:) = vi(:,1)
!      ver(:) = vi(:,2)
!      top(:) = vi(:,3)
!      length = 0.0
!      IF(inc(1)>1) length(1) = do_blen(lspace,hor,zero)
!      IF(inc(2)>1) length(2) = do_blen(lspace,ver,zero)
!      IF(inc(3)>1) length(3) = do_blen(lspace,top,zero)
!      CALL angle(angle_vh, inc(2), inc(1), ver, hor,        &
!                 length(2), length(1), extr_ord, extr_abs,  &
!                 ratio_vh, aver_vh)
!      IF(ltop .AND. inc(3)>1) THEN
!         CALL angle(angle_ht, inc(3), inc(1), hor, top,        &
!                    length(1), length(3), extr_abs, extr_top,  &
!                    ratio_ht, aver_ht)
!         CALL angle(angle_tv, inc(3), inc(2), top, ver,        &
!                    length(3), length(2), extr_top, extr_ord,  &
!                    ratio_tv, aver_tv)
!      ENDIF
!!                                                                       
      WRITE (output_io, 1400) ( (eck (i, j), i = 1, 3), j = 1, 4) 
      WRITE (output_io, 1410) (vi (i, 1), i = 1, 3), length(1), &
                              (vi (i, 2), i = 1, 3), length(2), &
                              (vi (i, 3), i = 1, 3), length(3)
      WRITE (output_io, 1420) (inc (i), i = 1, 3), extr_achs (extr_abs),&
      extr_achs (extr_ord), extr_achs(extr_top)                          
      WRITE (output_io, 1430) 'v/h',angle_vh, ratio_vh, aver_vh
      IF(ltop .AND. inc(3)>1) THEN
         WRITE (output_io, 1430) 'h/t',angle_ht, ratio_ht, aver_ht
         WRITE (output_io, 1430) 't/v',angle_tv, ratio_tv, aver_tv
      ENDIF
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
     &          F7.4,' A == ', F8.4,'keV')                                              
 1201 FORMAT (  '   Radiation          : ',A,', wavelength = ',         &
     &          F7.4,' A == ', F8.4,'meV')                                              
 1202 FORMAT (  '   Radiation          : ',A,', wavelength = ',         &
     &          F7.4,' A == ', F8.4,' eV')                                              
 1210 FORMAT (  '   Radiation          : ',A,', wavelength = ',A4,      &
     &          ' = ',F7.4,' A')                                        
 1300 FORMAT (  '   Temp. factors      : ',A) 
 1310 FORMAT (  '   Anomalous scat.    : ',A) 
 1400 FORMAT (/,' Reciprocal layer     : ',/                            &
     &          '   lower left  corner : ',3(2x,f9.4),/                 &
     &          '   lower right corner : ',3(2x,f9.4),/                 &
     &          '   upper left  corner : ',3(2x,f9.4),/                 &
     &          '   top   left  corner : ',3(2x,f9.4))                  
 1410 FORMAT (/,'   hor. increment     : ',3(2x,f9.4),2x,               &
     &          ' -> ',f9.4,' A**-1',/                                  &
     &          '   vert. increment    : ',3(2x,f9.4),2x,               &
     &          ' -> ',f9.4,' A**-1',/                                  &
     &          '   top   increment    : ',3(2x,f9.4),2x,               &
     &          ' -> ',f9.4,' A**-1')                                   
 1420 FORMAT (  '   # of points        :  ',I5,' x',I5,' x',I5,'  (',   &
     &          A1,',',A1,',',A1,')')                                          
 1430 FORMAT (  '   Angle Ratio Aver ',a3, 3x,f9.4,' degrees',3x        &
     &                                    ,2(2x,f9.4))                  
 1500 FORMAT (/,' Zone axis pattern',/                                   &
               ,'   axis [uvw]         : ',3(2x,f9.4))  
 1510 FORMAT (  '   resolution         : ',  2x,f9.4, '  2sin(theta)/lambda')  
      END SUBROUTINE four_show 
!
END MODULE fourier_menu
