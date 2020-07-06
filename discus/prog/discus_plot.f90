MODULE discus_plot_menu
!
LOGICAL, SAVE :: pl_keep= .FALSE.
CONTAINS
!
!*****7*****************************************************************
!                                                                       
SUBROUTINE plot 
!-                                                                      
!     Write the structure properly formatted for structure display      
!     programs and KUPLOT                                               
!+                                                                      
USE discus_config_mod 
USE discus_allocate_appl_mod
USE crystal_mod 
USE metric_mod
USE modify_mod
USE molecule_mod
USE discus_plot_mod 
USE discus_show_menu
USE get_iscat_mod
USE prop_para_func
USE update_cr_dim_mod
USE trafo_mod
!
USE envir_mod
USE ber_params_mod
USE build_name_mod
USE calc_expr_mod
USE doact_mod 
USE do_wait_mod
USE errlist_mod 
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
USE string_convert_mod
USE sup_mod
USE precision_mod
USE take_param_mod
USE str_comp_mod
USE support_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER, PARAMETER :: ITMP     = 78  ! temporary unit number
INTEGER, PARAMETER :: MIN_PARA = 20  ! A command requires at leaset these no of parameters
INTEGER maxw 
LOGICAL lold 
PARAMETER (lold = .false.) 
!                                                                       
CHARACTER(LEN=PREC_STRING), DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: cpara ! (MAX(10,MAXSCAT)) 
REAL(KIND=PREC_DP) , DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: werte ! (MAX(10,MAXSCAT)) 
INTEGER            , DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: lpara ! (MAX(10,MAXSCAT))
CHARACTER(LEN=PREC_STRING) :: line, zeile
CHARACTER(LEN=PREC_STRING) :: tempfile
CHARACTER(LEN=LEN(prompt)) :: orig_prompt 
CHARACTER(5) befehl 
CHARACTER(1) cdum 
REAL :: size, rr=0.0, rg=0.0, rb=0.0
INTEGER lp, length 
INTEGER :: ianz, i, j, is, it, ic, lbef 
!      INTEGER :: npoly     
INTEGER indxg 
INTEGER         :: nscat = 1
INTEGER         :: nsite = 1
INTEGER         :: ios

LOGICAL lend, l_select 
LOGICAL :: labs = .FALSE.
LOGICAL :: lord = .FALSE.
LOGICAL :: lnor = .FALSE.
!
!                                                                       
!                                                                       
INTEGER, PARAMETER :: NOPTIONAL = 10
CHARACTER(LEN=   5), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=PREC_STRING), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent!opt. para present
REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 4 ! Number of values to calculate 
!
DATA oname  / 'dmin' , 'dmax' , 'nmin' , 'nmax' , 'face' , 'hue'  , 'color', 'plot' , 'kill'  , 'keep'/
DATA loname /  4     ,  4     ,  4     ,  4     ,  4     ,  3     ,  5     ,  4     ,  4      ,  4    /
!
! Always provide fresh default values, Repeat for each command
opara  =   (/ '0.000', '0.000', '0    ', '0    ', 'flat ', 'solid', 'auto ', 'none ', 'none ', 'no   '/)
lopara =   (/  5     ,  5     ,  1     ,  1     ,  4     ,  5     ,  4     ,  4     ,  4     ,  2     /)
owerte =   (/  0.00  ,  0.00  ,  0.    ,  0.    ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   /)
!
IF(pl_init) THEN
   labs = .FALSE.
   lord = .FALSE.
   lnor = .FALSE.
ENDIF
!
!                                                                       
      maxw = MAX(MIN_PARA,MAXSCAT+1)
      lend = .false. 
!
!     Allocate the necessary arrays
!
!     IF ( cr_nscat > PL_MAXSCAT .or. mole_num_type > PL_MAXSCAT .OR. & 
!          MAXSCAT > PL_MAXSCAT ) THEN
!        nscat = max ( cr_nscat, mole_num_type, MAXSCAT)
!        CALL alloc_plot ( nscat )
!        IF(ier_num < 0) THEN
!          RETURN
!        ENDIF
!     ENDIF
!                                                                       
      orig_prompt = prompt
      prompt = prompt (1:len_str (prompt) ) //'/plot' 
main: DO while (.not.lend)                                        !  Main loop
      CALL no_error 
      CALL get_cmd (line, length, befehl, lbef, zeile, lp, prompt) 
   IF (ier_num.eq.0) then                                         ! No error get command
      IF (line /= ' '      .and. line(1:1) /= '#' .and. &
          line /= char(13) .and. line(1:1) /= '!'        ) THEN   ! Non-blank
!                                                                       
!     ----search for "="                                                
!                                                                       
            indxg = index (line, '=') 
if_gleich:  IF (indxg /= 0.AND..NOT. (str_comp (befehl, 'echo', 2, lbef, 4) ) &
                       .AND..NOT. (str_comp (befehl, 'syst', 2, lbef, 4) ) &
                       .AND..NOT. (str_comp (befehl, 'help', 2, lbef, 4)   &
                       .OR.        str_comp (befehl, '?   ', 2, lbef, 4) ) &
                       .AND. INDEX(line,'==') == 0                        ) THEN   ! DO_MATH?
!                                                                       
!     ------evaluate an expression and assign the value to a variabble   
!                                                                       
               CALL do_math (line, indxg, length) 
         ELSE  if_gleich
!                                                                       
!------ ----execute a macro file                                        
!                                                                       
            IF (befehl (1:1) .eq.'@') then                ! macro or reset or all other commands
                  IF (length.ge.2) then 
                     line(1:length-1) = line(2:length)
                     length = 1
                     CALL file_kdo(line, length)
                  ELSE 
                     ier_num = - 13 
                     ier_typ = ER_MAC 
                  ENDIF 
!                                                                 
!     --continues a macro 'continue'                                    
!                                                                       
               ELSEIF (str_comp (befehl, 'continue', 3, lbef, 8) ) then 
                  CALL macro_continue (zeile, lp) 
!                                                                       
!------ ----list asymmetric unit 'asym'                                 
!                                                                       
               ELSEIF (str_comp (befehl, 'asym', 2, lbef, 4) ) then 
                  CALL show_asym 
!                                                                       
!------ ----list atoms present in the crystal 'chem'                    
!                                                                       
               ELSEIF (str_comp (befehl, 'chem', 2, lbef, 4) ) then 
                  CALL show_chem 
!                                                                       
!------ ----Echo a string, just for interactive check in a macro 'echo' 
!                                                                       
               ELSEIF (str_comp (befehl, 'echo', 2, lbef, 4) ) then 
                  CALL echo (zeile, lp) 
!                                                                       
!     ----exit 'exit'                                                   
!                                                                       
               ELSEIF (str_comp (befehl, 'exit', 3, lbef, 4) ) then 
                  lend = .true. 
!                                                                       
!     ----help 'help','?'                                               
!                                                                       
      ELSEIF (str_comp (befehl, 'help', 2, lbef, 4) .or.  &
              str_comp (befehl, '?   ', 1, lbef, 4) )    then
                  IF (str_comp (zeile, 'errors', 2, lp, 6) ) then 
                     lp = lp + 7 
                     CALL do_hel ('discus '//zeile, lp) 
                  ELSE 
                     lp = lp + 12 
                     CALL do_hel ('discus plot '//zeile, lp) 
                  ENDIF 
!                                                                       
!-----  ------Operating System Kommandos 'syst'                         
!                                                                       
               ELSEIF (str_comp (befehl, 'syst', 2, lbef, 4) ) then 
                  IF (zeile.ne.' ') then 
                     CALL do_operating (zeile, lp) 
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
!     ----reset plot 'reset'                                                
!
            ELSEIF (str_comp (befehl, 'reset', 3, lbef, 5) ) THEN   ! macro or reset or all other commands
               nscat = max ( cr_nscat, mole_num_type, MAXSCAT)
               CALL  plot_reset
               labs = .FALSE.
               lord = .FALSE.
               lnor = .FALSE.
            ELSE                                                   ! macro or reset or all other commands
!
!     Allocate the necessary arrays
!
               IF(cr_nscat > PL_MAXSCAT .or. mole_num_type > PL_MAXSCAT .OR. & 
                  MAXSCAT > PL_MAXSCAT ) THEN
                  nscat = MAX(cr_nscat, mole_num_type, MAXSCAT)
                  nsite = MAX(cr_ncatoms, MAXSCAT)
                  CALL alloc_plot(nscat, nsite )
                  IF(ier_num < 0) THEN
                     RETURN
                  ENDIF
                  PL_MAXSITE = UBOUND(pl_lsite,1)
               ENDIF
!                                                                       
!     ----Select the abscissa for a projection onto a to plot           
!           slice 'absc'                                                
!                                                                       
               IF (str_comp (befehl, 'absc', 2, lbef, 4) ) THEN 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0.and.ianz.eq.3) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        DO i = 1, 3 
!*****7*****************************************************************
                           pl_abs (i) = werte (i) 
                        ENDDO 
                        labs = .TRUE.
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
!                                                                       
!     ----Set the sequence of columns 'colu'                            
!                                                                       
               ELSEIF (str_comp (befehl, 'colu', 3, lbef, 4) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (ianz.eq.1) then 
                        IF (cpara (1) .eq.'xyz'.or.cpara (1).eq.'yzx'.or. &
                            cpara (1) .eq.'zxy'.or.cpara (1).eq.'xzy'.or. &
                            cpara (1) .eq.'zyx'.or.cpara (1).eq.'yxz') then
                           pl_col = cpara (1)(1:lpara(1))
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
!------ --setting the scale of the plot                                 
!                                                                       
               ELSEIF (str_comp (befehl, 'scale', 3, lbef, 5) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        IF(NINT(werte (1))  ==  -1) THEN 
                           pl_scale (0) = - 1 
                        ELSEIF (werte(1) .gt. 0.0D0) then 
                           pl_scale (0) = 0 
                           pl_scale (1) = werte (1) 
                        ELSE 
                           ier_num = - 6 
                           ier_typ = ER_COMM 
                        ENDIF 
                     ENDIF 
                  ENDIF 
!                                                                       
!------ --selecting/deselecting atoms                                   
!                                                                       
               ELSEIF (str_comp (befehl, 'sele', 3, lbef, 4) .or.       &
                        str_comp (befehl, 'dese', 2, lbef, 4) ) then
!                                                                       
                  CALL atom_select (zeile, lp, 0, PL_MAXSCAT, pl_latom, &
                  pl_lsite, 0, PL_MAXSITE, &
                  pl_sel_atom, lold,        &
                  str_comp (befehl, 'sele', 3, lbef, 4) )               
!                                                                       
!------ --selecting/deselecting of molecules                            
!                                                                       
               ELSEIF (str_comp (befehl, 'msel', 2, lbef, 4) .or.       &
                       str_comp (befehl, 'mdes', 2, lbef, 4) ) then
!                   
                  CALL mole_select (zeile, lp, 0, PL_MAXSCAT, pl_latom, &
                  pl_sel_atom, str_comp (befehl, 'msel', 2, lbef, 4) )
!                                                                       
!------ --Handle property settings 'property'                           
!                                                                       
               ELSEIF (str_comp (befehl, 'property', 4, lbef, 8) ) then 
!                                                                       
                  CALL property_select (zeile, lp, pl_sel_prop) 
!                                                                       
!     ----Select the background 'back'                                       
!                                                                       
               ELSEIF (str_comp (befehl, 'back', 2, lbef, 4) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0.and. ianz.eq.3 ) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        IF(0<=NINT(werte(1)) .AND. NINT(werte(1))<=255 .AND. &
                           0<=NINT(werte(2)) .AND. NINT(werte(2))<=255 .AND. &
                           0<=NINT(werte(3)) .AND. NINT(werte(3))<=255) THEN
                           pl_back(1) = NINT(werte(1))
                           pl_back(2) = NINT(werte(2))
                           pl_back(3) = NINT(werte(3))
                        ELSE 
                           ier_num = - 6 
                           ier_typ = ER_COMM 
                        ENDIF 
                     ENDIF 
                  ENDIF 
!                                                                       
!     ----Select the bonds 'bond'                                       
!                                                                       
               ELSEIF (str_comp (befehl, 'bonds', 2, lbef, 5) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0.and. (ianz.eq.6.or.ianz.eq.8) ) then 
!                                                                       
                     DO j = 0, MAXSCAT 
                     pl_batom_a (j) = .false. 
                     pl_batom_e (j) = .false. 
                     ENDDO 
                     zeile = cpara (1) (1:lpara (1) ) 
                     CALL atom_select (zeile, lp, 0, MAXSCAT, pl_batom_a, &
                     pl_lsite, 0, PL_MAXSITE, &
                     pl_sel_atom, lold,  .true.)
                     zeile = cpara (2) (1:lpara (2) ) 
                     CALL atom_select (zeile, lp, 0, MAXSCAT, pl_batom_e, &
                     pl_lsite, 0, PL_MAXSITE, &
                     pl_sel_atom, lold, .true.)
                     CALL del_params (2, ianz, cpara, lpara, maxw) 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        DO i = 0, MAXSCAT 
                        DO j = 0, MAXSCAT 
                        IF (pl_batom_a (i) .and.pl_batom_e (j) ) then 
                           pl_bond (i, j) = .true. 
                           pl_bond_len (1, i, j) = werte (1) 
                           pl_bond_len (2, i, j) = werte (2) 
                           pl_bond_rad (i, j) = werte (3) 
                           pl_bond_col (1, i, j) = werte (4) 
                           IF (ianz.eq.4) then 
                              pl_bond_col (2, i, j) = werte (4) 
                              pl_bond_col (3, i, j) = werte (4) 
                           ELSEIF (ianz.eq.6) then 
                              pl_bond_col (2, i, j) = werte (5) 
                              pl_bond_col (3, i, j) = werte (6) 
                           ENDIF 
                        ENDIF 
                        ENDDO 
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
!     ----Select the extend of crystal space to be plotted 'exte'       
!                                                                       
               ELSEIF (str_comp (befehl, 'exte', 3, lbef, 4) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (ianz.eq.1.and.cpara (1) (1:3) .eq.'all') then 
                        pl_ext_all (1) = .true. 
                        pl_ext_all (2) = .true. 
                        pl_ext_all (3) = .true. 
                     ELSEIF (ianz.eq.3) then 
                        cdum = cpara (1) (1:1) 
                        cpara (1) = '0.0' 
                        lpara (1) = 3 
                        CALL ber_params (ianz, cpara, lpara, werte, maxw)
                        IF (ier_num.eq.0) then 
                           IF (cdum.eq.'x') then 
                              pl_dim (1, 1) = werte (2) 
                              pl_dim (1, 2) = werte (3) 
                              pl_ext_all (1) = .false. 
                           ELSEIF (cdum.eq.'y') then 
                              pl_dim (2, 1) = werte (2) 
                              pl_dim (2, 2) = werte (3) 
                              pl_ext_all (2) = .false. 
                           ELSEIF (cdum.eq.'z') then 
                              pl_dim (3, 1) = werte (2) 
                              pl_dim (3, 2) = werte (3) 
                              pl_ext_all (3) = .false. 
                           ELSE 
                              ier_num = - 6 
                              ier_typ = ER_COMM 
                           ENDIF 
                        ENDIF 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ENDIF 
!                                                                       
!     ----Select the reciprocal space direction normal to plot          
!           slice 'hkl'                                                 
!                                                                       
               ELSEIF (str_comp (befehl, 'hkl', 2, lbef, 3) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0.and.ianz.eq.3) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        DO i = 1, 3 
                        pl_hkl (i) = werte (i) 
                        ENDDO 
                        CALL trans (pl_hkl, cr_rten, pl_uvw, 3) 
                        lnor = .TRUE.
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
!                                                                       
!     ----Define output of complete molecules/only origin               
!                                                                       
               ELSEIF (str_comp (befehl, 'mole', 2, lbef, 4) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (ianz.eq.1) then 
                        CALL do_cap (cpara (1) ) 
                        pl_mol_all = (cpara (1) (1:1) .eq.'A') 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ENDIF 
!                                                                       
!     ----Select the ordinate for a projection onto a to plot           
!           slice 'ordi'                                                
!                                                                       
               ELSEIF (str_comp (befehl, 'ordi', 2, lbef, 4) ) THEN 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0.and.ianz.eq.3) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        DO i = 1, 3 
                           pl_ord (i) = werte (i) 
                        ENDDO 
                        lord = .TRUE.
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
!                                                                       
!------ define name of output file 'outf'                               
!                                                                       
               ELSEIF (str_comp (befehl, 'outf', 1, lbef, 4) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1)
                     IF (ier_num.eq.0) then 
                        pl_out = cpara (1) (1:lpara(1))
                     ENDIF 
                  ENDIF 
!
!     ----Select the Polyhedra to be plotted by Jmol 'poly'
!
               ELSEIF (str_comp (befehl, 'poly', 3, lbef, 4) ) then 
!
      ! Always provide fresh default values, Repeat for each command
      opara (1:7) =   (/ '0.000', '0.000', '0    ', '0    ', 'flat ', 'solid', 'auto ' /)
      lopara(1:7) =   (/  5     ,  5     ,  1     ,  1     ,  4     ,  5     ,  4      /)
      owerte(1:7) =   (/  0.00  ,  0.00  ,  0.    ,  0.    ,  0.0   ,  0.0   ,  0.0    /)
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                                    oname, loname, opara, lopara, lpresent, owerte)
                  IF (ier_num.eq.0) then 
                     IF(str_comp(cpara(1), 'off', 3, lpara(1), 3)) THEN
                        pl_poly_n = 0
                     ELSE
                     line = cpara(2)
                     IF(cpara(1)(1:1)=='''' .AND. cpara(1)(lpara(1):lpara(1))=='''') THEN
                        zeile = cpara(1)(2:lpara(1)-1)
                     ELSE
                        zeile = cpara(1)(1:lpara(1))
                     ENDIF
                     CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                     CALL get_iscat (ianz, cpara, lpara, werte, maxw, lold)                               
                     IF (ier_num.eq.0) then 
                        IF(NINT(werte(1))==-1) THEN
                           pl_poly_c(:) = .TRUE.
                        ELSE
                           DO i=1,ianz
                              pl_poly_c(NINT(werte(i))) = .TRUE.
                           ENDDO
                           pl_poly_c(-1) = .FALSE.
                        ENDIF 
                        i = LEN_TRIM(line)
                        zeile = ' '
                        IF(line(1:1)=='''' .AND. line(i:i)=='''') THEN
                           zeile(1:i-2) = line(2:i-1)
                        ELSE
                           zeile = line(1:i)
                        ENDIF
                        CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                        CALL get_iscat (ianz, cpara, lpara, werte, maxw, lold)                               
                        IF (ier_num.eq.0) then 
                           IF(NINT(werte(1))==-1) THEN
                              pl_poly_o(:) = .TRUE.
                           ELSE
                              DO i=1,ianz
                                 pl_poly_o(NINT(werte(i))) = .TRUE.
                              ENDDO
                              pl_poly_o(-1) = .FALSE.
                           ENDIF 
                           IF(owerte(1)>0.0) pl_poly_dmin = owerte(1)
                           IF(owerte(2)>0.0) pl_poly_dmax = owerte(2)
                           pl_poly_nmax = NINT(owerte(4))
                           IF(NINT(owerte(3))==0) owerte(3) = pl_poly_nmax
                           pl_poly_nmin = NINT(owerte(3))
                           pl_poly_face = opara(5)=='flat'
                           pl_poly_hue  = opara(6)=='trans'
                           pl_poly_col  = opara(7)(1:MIN(LEN(pl_poly_col),lopara(7)))
                           pl_poly_n    = 1
                        ENDIF 
                     ENDIF 
                     ENDIF 
                  ENDIF 
!                                                                       
!     ----Select plotting program 'prog'                                
!                                                                       
               ELSEIF (str_comp (befehl, 'prog', 4, lbef, 4) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0.and.ianz.eq.1.or.ianz.eq.2) then 
                     CALL do_cap (cpara (1) ) 
                     IF (cpara (1) (1:1) .eq.'K'.and.ianz.eq.1) then 
                        pl_prog = 'kupl' 
                     ELSEIF (cpara (1) (1:1) .eq.'G'.and.ianz.eq.1) then
                        pl_prog = 'gnuplot' 
                     ELSEIF (cpara (1) (1:1) .eq.'A'.and.ianz.eq.1) then
                        pl_prog = 'atoms' 
                     ELSEIF (cpara (1) (1:2) .eq.'DR'.and.ianz.eq.1) then
                        pl_prog = 'drawxtl' 
                     ELSEIF (cpara (1) (1:1) .eq.'X'.and.ianz.eq.1) then
                        pl_prog = 'xbs' 
                     ELSEIF (cpara (1) (1:1) .eq.'C'.and.ianz.eq.1) then
                        pl_prog = 'cif' 
                     ELSEIF (cpara (1) (1:1) .eq.'F') then 
                        pl_prog = 'frames' 
                        IF (ianz.eq.2) then 
                           IF (str_comp (cpara (2) , 'init', 1,      &
                               lpara ( 2) , 4) ) then
                              pl_append = .false. 
                           ELSEIF (str_comp (cpara (2) , 'appe', 1,  &
                           lpara (2) , 4) ) then
                              pl_append = .true. 
                           ELSEIF (ianz.eq.1) then 
                              pl_append = .true. 
                           ELSE 
                              ier_num = - 6 
                              ier_typ = ER_COMM 
                           ENDIF 
                        ENDIF 
                     ELSEIF (cpara(1)(1:1) .eq.'D'.and.ianz.eq.1) then
                        pl_prog = 'diamond' 
                     ELSEIF (cpara(1)(1:1) .eq.'J'.and.ianz.eq.1) then
                        pl_prog = 'jmol' 
!                       Test existence of jmol / java 
                        WRITE(tempfile, '(a,i10.10)') '/tmp/which_jmol.', PID
                        IF(operating=='Linux') THEN                      
                           WRITE(line,'(a,i10.10)')  'which jmol > /tmp/which_jmol.', PID
                        ELSEIF(operating=='Linux_WSL') THEN                      
                           WRITE(line,'(a,i10.10)')  'which jmol > /tmp/which_jmol.', PID
                        ELSEIF(operating=='Windows') THEN                      
                           WRITE(line,'(a,i10.10)')  'which java > /tmp/which_jmol.', PID
                        ELSEIF(operating(1:6)=='cygwin') THEN                      
                           WRITE(line,'(a,i10.10)')  'which java > /tmp/which_jmol.', PID
                        ELSEIF(operating(1:6)=='darwin') THEN                      
                           WRITE(line,'(a,i10.10)')  'which jmol > /tmp/which_jmol.', PID
                        ENDIF
                           CALL system(line)
                           CALL oeffne( ITMP, tempfile, 'old')
                           IF(ier_NUM==0) THEN
                              READ(ITMP,'(a)',IOSTAT=ios) pl_jmol
                              CLOSE(ITMP)
                              IF(IS_IOSTAT_END(ios)) THEN
                                 ier_num = -152
                                 ier_typ = ER_APPL
                                 pl_jmol = ' '
                              ELSE
                                 IF(operating(1:6)=='cygwin' ) THEN !   &
                                    i = LEN_TRIM(pl_jmol)
                                    IF(pl_jmol(i-3:i) == 'java') THEN
                                       pl_jmol = 'cd /tmp      && java -Xmx512m -jar ./Jmol.jar ' 
                                    ELSEIF(pl_jmol(1:14) == 'which: no java') THEN
                                       ier_num = -152
                                       ier_typ = ER_APPL
                                       pl_jmol = ' '
                                    ENDIF
                                 ELSEIF(operating(1:7)=='Windows'    ) THEN
                                    i = LEN_TRIM(pl_jmol)
                                    IF(pl_jmol(i-3:i) == 'java') THEN
                                       pl_jmol = 'cd /opt/jmol && java -Xmx512m -jar ./Jmol.jar ' 
                                    ELSEIF(pl_jmol(1:14) == 'which: no java') THEN
                                       ier_num = -152
                                       ier_typ = ER_APPL
                                       pl_jmol = ' '
                                    ENDIF
                                 ENDIF
                              ENDIF
                           ENDIF
                        WRITE(line,'(a,a)') 'rm -f ', tempfile(1:LEN_TRIM(tempfile))
                        CALL system(line)
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ENDIF 
!                                                                       
!     ----run plot 'run'                                                
!                                                                       
               ELSEIF (str_comp (befehl, 'run', 3, lbef, 3) ) then 
!                 Always provide fresh default values, Repeat for each command
                  opara (8:10) =   (/ 'none ', 'none ', 'no   '/)
                  lopara(8:10) =   (/  4     ,  4     ,  2     /)
                  owerte(8:10) =   (/  0.0   ,  0.0   ,  0.0   /)
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                                    oname, loname, opara, lopara, lpresent, owerte)
                  IF (ier_num.eq.0.and.ianz.eq.0) then 
                  IF (pl_out.ne.' ') then 
                     CALL update_cr_dim 
                     DO i = 1, 3 
                     IF (pl_ext_all (i) ) then 
                        DO j = 1, 2 
                        pl_dim (i, j) = cr_dim (i, j) 
                        ENDDO 
                     ENDIF 
                     ENDDO 
                     IF (.not. (pl_dim (1, 1) .eq.pl_dim (1, 2) .and. &
                                pl_dim (2, 1) .eq.pl_dim (2, 2) .and. &
                                pl_dim(3, 1) .eq.pl_dim (3, 2) ) ) then
                        l_select = .false. 
                        DO i = 0, cr_nscat 
                        IF (pl_latom (i) ) then 
                           l_select = .true. 
                        ENDIF 
                        ENDDO 
                        IF (l_select) then 
                           CALL plot_test_aon(lnor, labs, lord)
                           IF(ier_num == 0) THEN
                              CALL do_plot 
                              pl_keep = opara(10)(1:lopara(10))=='yes'
                              IF(ier_num == 0) THEN
                                 IF(opara(8)=='inter') THEN
                                    IF(pl_prog=='jmol') THEN
                                       CALL plot_inter(opara(9)(1:lopara(9)))
                                    ELSE
                                       ier_num = -152
                                       ier_typ = ER_APPL 
                                    ENDIF
                                 ENDIF
                              ENDIF
                           ENDIF
                        ELSE 
                           ier_num = - 3 
                           ier_typ = ER_APPL 
                        ENDIF 
                     ELSE 
                        ier_num = - 4 
                        ier_typ = ER_APPL 
                     ENDIF 
                  ELSE 
                     ier_num = - 37 
                     ier_typ = ER_APPL 
                  ENDIF 
                     ENDIF 
!                                                                       
!     ----set atom representation for KUPLOT                            
!                                                                       
               ELSEIF (str_comp (befehl, 'set', 3, lbef, 3) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (ianz.eq.4.or.ianz.eq.6) then 
!                                                                       
!     ----------store 1. parameter for later analysis                   
!                                                                       
                        line = cpara (1) 
                        lp = lpara (1) 
                        cpara (1) = '0' 
                        lpara (1) = 1 
                        CALL ber_params (ianz, cpara, lpara, werte, maxw)
                        IF (ier_num.eq.0) then 
                           it = nint (werte (2) ) 
                           ic = nint (werte (3) ) 
                           IF (ianz.eq.4) then 
                              rr = werte (3) 
                              rg = werte (3) 
                              rb = werte (3) 
                           ELSEIF (ianz.eq.6) then 
                              rr = werte (3) 
                              rg = werte (4) 
                              rb = werte (5) 
                           ENDIF 
                           IF ( (ic.gt.0.and.ic.le.15.and.it.gt.0.and. &
                                 it.le.12) .or.pl_prog.eq.'xbs'.or.    &
                                 pl_prog.eq.'diamond') then
                              size = werte (ianz) 
!                                                                       
!     --------------Now analyse the atom name/type                      
!                                                                       
                              cpara (1) = line 
                              lpara (1) = lp 
                              ianz = 1 
                              CALL get_iscat (ianz, cpara, lpara, werte,&
                              maxw, lold)                               
                              IF (nint (werte (1) ) .eq. - 1) then 
                                 DO is = 1, cr_nscat 
                                 pl_typ (is) = it 
                                 pl_color (is) = ic 
                                 pl_siz (is) = size 
                                 pl_rgb (1, is) = rr 
                                 pl_rgb (2, is) = rg 
                                 pl_rgb (3, is) = rb 
                                 ENDDO 
                              ELSE 
                                 DO i = 1, ianz 
                                 is = nint (werte (i) ) 
                                 pl_typ (is) = it 
                                 pl_color (is) = ic 
                                 pl_siz (is) = size 
                                 pl_rgb (1, is) = rr 
                                 pl_rgb (2, is) = rg 
                                 pl_rgb (3, is) = rb 
                                 ENDDO 
                              ENDIF 
                           ELSE 
                              ier_num = - 59 
                              ier_typ = ER_APPL 
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
                  CALL plot_show 
!                                                                       
!     ----Thickness of the plot slice in Angstroem 'thickness'          
!                                                                       
               ELSEIF (str_comp (befehl, 'thic', 2, lbef, 4) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0.and.ianz.eq.1) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        pl_width = werte (1) 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
!                                                                       
!------ define title for frames    'title'                              
!                                                                       
               ELSEIF (str_comp (befehl, 'title', 2, lbef, 5) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1)
                     IF (ier_num.eq.0) then 
                        pl_title = cpara (1) (1:lpara(1))
                     ENDIF 
                  ENDIF 
!                                                                       
!     ----Select the direct space direction normal to plot slice 'uvw'  
!                                                                       
               ELSEIF (str_comp (befehl, 'uvw', 1, lbef, 3) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0.and.ianz.eq.3) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        DO i = 1, 3 
                        pl_uvw (i) = werte (i) 
                        ENDDO 
                        CALL trans (pl_uvw, cr_gten, pl_hkl, 3) 
                        lnor = .TRUE.
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
!                                                                       
!     ----Select a point in direct space within the plot slice 'vec'    
!                                                                       
               ELSEIF (str_comp (befehl, 'vect', 1, lbef, 4) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0.and.ianz.eq.3) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        DO i = 1, 3 
                        pl_vec (i) = werte (i) 
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
!     ----set typ of output                                             
!                                                                       
               ELSEIF (str_comp (befehl, 'type', 2, lbef, 4) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     IF (ianz.eq.1) then 
                        CALL do_cap (cpara (1) ) 
                        pl_dens = (cpara (1) (1:1) .eq.'P') 
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ENDIF 
!                                                                       
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_COMM 
               ENDIF                                            ! All commands
            ENDIF                                               ! macro or reset or all other commands
         ENDIF if_gleich                                        ! DO Math ?
      ENDIF                                                     ! non_blank
   ENDIF                                                        ! No error get command
!                                                                       
   IF (ier_num.ne.0) then 
      CALL errlist 
      IF (ier_sta.ne.ER_S_LIVE) then 
         IF (lmakro .OR. lmakro_error) THEN  ! Error within macro or termination errror
            IF(sprompt /= prompt ) THEN
               ier_num = -10
               ier_typ = ER_COMM
               ier_msg(1) = ' Error occured in plot menu'
               prompt_status = PROMPT_ON 
               prompt = orig_prompt
               RETURN
            ELSE
               CALL macro_close 
               prompt_status = PROMPT_ON 
            ENDIF 
         ENDIF 
         IF (lblock) then 
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
ENDDO  main                 ! Main loop
!                                                                       
prompt = orig_prompt
!
END SUBROUTINE plot                           
!
!*****7*****************************************************************
      SUBROUTINE plot_show 
!                                                                       
      USE discus_config_mod 
      USE crystal_mod 
      USE atom_name
      USE metric_mod
      USE modify_mod
      USE molecule_mod 
      USE discus_plot_mod 
      USE prop_char_mod
!      USE prop_para_func
!      USE prop_para_mod 
!
      USE errlist_mod 
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      CHARACTER(9) at_name_i 
      CHARACTER(32) c_property 
      INTEGER i, j, length 
      LOGICAL lspace 
      REAL null (3), w_na, w_no, w_ao, aver 
      REAL u (3), w (3), angl, b_a, b_o 
!                                                                       
!     REAL do_bang, do_blen 
!                                                                       
      lspace = .true. 
      null (1) = 0.0 
      null (2) = 0.0 
      null (3) = 0.0 
!                                                                       
      WRITE (output_io, 3000) pl_prog 
      WRITE (output_io, 3005) pl_out 
      IF (pl_dens) then 
         WRITE (output_io, 3007) 'projection in single unit cell' 
      ELSE 
         WRITE (output_io, 3007) 'crystal' 
      ENDIF 
      WRITE (output_io, 3010) ( (pl_dim (i, j), j = 1, 2), i = 1, 3) 
      WRITE (output_io, 3020) pl_col 
      WRITE (output_io, 3021) pl_hkl 
      WRITE (output_io, 3022) pl_uvw 
!                                                                       
      IF (do_blen (lspace, null, pl_uvw) .ne.0.0) then 
!                                                                       
!     --A slice is to be plottet normal to a vector                     
!                                                                       
         WRITE (output_io, 3025) pl_abs 
         WRITE (output_io, 3026) pl_ord 
         w_na = do_bang (lspace, pl_uvw, null, pl_abs) 
         w_no = do_bang (lspace, pl_uvw, null, pl_ord) 
         w_ao = do_bang (lspace, pl_abs, null, pl_ord) 
         WRITE (output_io, 3027) w_na, w_no, w_ao 
!                                                                       
!     --Make sure that none of the abscissa or ordinate are zero        
!                                                                       
         b_o = do_blen (lspace, null, pl_ord) 
         IF (b_o.eq.0) then 
            ier_msg (1) = 'length of ordinate is zero' 
            ier_num = - 32 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
         b_a = do_blen (lspace, null, pl_abs) 
         IF (b_a.eq.0) then 
            ier_msg (1) = 'length of abscissa is zero' 
            ier_num = - 32 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
!                                                                       
!     --Determine correct ratio for the plot. Depends on pl_col         
!                                                                       
         IF (pl_col (2:2) .eq.'x') then 
            u (1) = pl_abs (1) 
            u (2) = pl_abs (2) 
            u (3) = pl_abs (3) 
         ELSEIF (pl_col (2:2) .eq.'y') then 
            u (1) = pl_ord (1) 
            u (2) = pl_ord (2) 
            u (3) = pl_ord (3) 
         ELSEIF (pl_col (2:2) .eq.'z') then 
            u (1) = pl_uvw (1) 
            u (2) = pl_uvw (2) 
            u (3) = pl_uvw (3) 
         ENDIF 
         b_o = do_blen (lspace, null, u) 
         IF (pl_col (1:1) .eq.'x') then 
            w (1) = pl_abs (1) 
            w (2) = pl_abs (2) 
            w (3) = pl_abs (3) 
         ELSEIF (pl_col (1:1) .eq.'y') then 
            w (1) = pl_ord (1) 
            w (2) = pl_ord (2) 
            w (3) = pl_ord (3) 
         ELSEIF (pl_col (1:1) .eq.'z') then 
            w (1) = pl_uvw (1) 
            w (2) = pl_uvw (2) 
            w (3) = pl_uvw (3) 
         ENDIF 
         b_a = do_blen (lspace, null, w) 
         aver = b_o / b_a 
         angl = do_bang (lspace, u, null, w) 
         WRITE (output_io, 3028) aver 
         WRITE (output_io, 3029) angl 
      ELSE 
!                                                                       
!     --Projection onto xy, yz or zx plane                              
!                                                                       
         IF (pl_col (2:2) .eq.'x') then 
            u (1) = 1.0 
            u (2) = 0.0 
            u (3) = 0.0 
         ELSEIF (pl_col (2:2) .eq.'y') then 
            u (1) = 0.0 
            u (2) = 1.0 
            u (3) = 0.0 
         ELSEIF (pl_col (2:2) .eq.'z') then 
            u (1) = 0.0 
            u (2) = 0.0 
            u (3) = 1.0 
         ENDIF 
         b_o = do_blen (lspace, null, u) 
         IF (pl_col (1:1) .eq.'x') then 
            w (1) = 1.0 
            w (2) = 0.0 
            w (3) = 0.0 
         ELSEIF (pl_col (1:1) .eq.'y') then 
            w (1) = 0.0 
            w (2) = 1.0 
            w (3) = 0.0 
         ELSEIF (pl_col (1:1) .eq.'z') then 
            w (1) = 0.0 
            w (2) = 0.0 
            w (3) = 1.0 
         ENDIF 
         b_a = do_blen (lspace, null, w) 
         aver = b_o / b_a 
         angl = do_bang (lspace, u, null, w) 
         WRITE (output_io, 3028) aver 
         WRITE (output_io, 3029) angl 
      ENDIF 
      WRITE (output_io, 3023) pl_vec 
      WRITE (output_io, 3024) pl_width 
      CALL char_prop_2 (c_property, pl_sel_prop (1), pl_sel_prop (0),   &
      length)                                                           
      WRITE (output_io, 3131) c_property (1:length) 
!                                                                       
      IF (pl_sel_atom) then 
         WRITE (output_io, 3050) 
         DO i = 0, cr_nscat 
         IF (pl_latom (i) ) then 
            at_name_i = at_name (i) 
            WRITE (output_io, 3060) i, at_name_i, pl_typ (i),   &
                                    pl_color (i), pl_siz (i)
         ENDIF 
         ENDDO 
      ELSE 
         IF (pl_mol_all) then 
            WRITE (output_io, 3065) 'All atoms in molecule' 
         ELSE 
            WRITE (output_io, 3065) 'Origin of molecule only' 
         ENDIF 
         WRITE (output_io, 3070) 
         DO i = 1, mole_num_type 
         IF (pl_latom (i) ) then 
            WRITE (output_io, 3080) i, pl_typ (i), pl_color (i), &
            pl_siz (i)                                                  
         ENDIF 
         ENDDO 
      ENDIF 
!                                                                       
 3000 FORMAT    ( ' Plotting Program        : ',a) 
 3005 FORMAT    ( ' Output file             : ',a) 
 3007 FORMAT    ( ' Output type             : ',a) 
 3010 FORMAT    (/' Extend of Crystal       : '/                 &
                  ' Xmin, Xmax              : ',2(g12.3E3)/      &
                  ' Ymin, Ymax              : ',2(g12.3E3)/      &
                  ' Zmin, Zmax              : ',2(g12.3E3))      
 3020 FORMAT    ( ' Sequence of coordinates : ',a) 
 3021 FORMAT    (/' HKL of Normal           : ',3(f8.3)) 
 3022 FORMAT    ( ' UVW of Normal           : ',3(f8.3)) 
 3025 FORMAT    (/' Abscissa of slice       : ',3(f8.3)) 
 3026 FORMAT    ( ' Ordinate of slice       : ',3(f8.3)) 
 3027 FORMAT    ( ' <(N,A), <(N,O), <(A,O)  : ',3(f8.3)) 
 3028 FORMAT    ( ' Aver v/h for KUPLOT     : ', (f8.3)) 
 3029 FORMAT    ( ' Angle                   : ', (f8.3),' degrees') 
 3023 FORMAT    (/' Point in plot slice     :  ',3(f8.3)) 
 3024 FORMAT    ( ' Half thickness of slice : ',  g12.3E3,' Angstroem' ) 
 3131 FORMAT    (/' Atom properties         : ','NMDOEI'/               &
                  '      absent=- ignored=. : ',a)               
 3050 FORMAT    ( ' Selected atoms          : ',                        &
                  '  #  name       mtyp  color  size ')           
 3065 FORMAT    ( ' Plotting status         : ',a) 
 3060 FORMAT    (27x,i3,2x,a9,3x,i2,4x,i2,2x,f6.2) 
 3070 FORMAT    ( ' Selected molecules      : ',                        &
                  '  #   mtyp  color  size ')                    
 3080 FORMAT    (27x,i3,4x,i2,4x,i2,2x,f6.2) 
      END SUBROUTINE plot_show                      
!*****7*****************************************************************
      SUBROUTINE do_plot 
!                                                                       
      USE discus_config_mod 
      USE crystal_mod 
      USE discus_plot_mod 
      USE discus_plot_init_mod
      USE discus_plot_export_mod
      USE errlist_mod 
USE support_mod
      IMPLICIT none 
!                                                                       
      INTEGER iff 
      PARAMETER (iff = 2) 
!                                                                       
       
!                                                                       
      LOGICAL lread 
      LOGICAL lexist 
!                                                                       
      IF (pl_prog.ne.'frames'.or..not.pl_append) then 
         lread = .false. 
         CALL oeffne (iff, pl_out, 'unknown') 
      ELSEIF (pl_prog.eq.'frames'.and.pl_append) then 
         INQUIRE (file = pl_out, exist = lexist) 
         IF (lexist) then 
            CALL oeffne_append (iff, pl_out, 'old')
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
!DBG          open (unit=iff,file=pl_out,status='old',access='append',  
!DBG     &               err=999)                                       
            ier_num = 0 
         ELSE 
            CALL oeffne (iff, pl_out, 'new') 
         ENDIF 
      ENDIF 
      IF (ier_num.ne.0) return 
!                                                                       
!------ Output program KUPLOT or GNUPLOT                                
!                                                                       
      IF (pl_prog.eq.'kupl'.or.pl_prog.eq.'gnuplot') then 
!                                                                       
!------ - Atom or molecule mode                                         
!                                                                       
         IF (pl_sel_atom) then 
            CALL plot_kuplot (iff, (pl_prog.eq.'kupl') ) 
         ELSE 
            CALL plot_kuplot_mol (iff, (pl_prog.eq.'kupl') ) 
         ENDIF 
!                                                                       
!------ Output ATOMS                                                    
!                                                                       
      ELSEIF (pl_prog.eq.'atoms') then 
!                                                                       
!------ - Atom or molecule mode                                         
!                                                                       
         IF (pl_sel_atom) then 
            CALL plot_atoms (iff) 
         ELSE 
            ier_num = - 76 
            ier_typ = ER_APPL 
         ENDIF 
!                                                                       
!------ Output DRAWxtl                                                  
!                                                                       
      ELSEIF (pl_prog.eq.'drawxtl') then 
!                                                                       
!------ - Atom or molecule mode                                         
!                                                                       
         IF (pl_sel_atom) then 
            CALL plot_drawxtl (iff) 
         ELSE 
            ier_num = - 76 
            ier_typ = ER_APPL 
         ENDIF 
!                                                                       
!------ Output CIF (subset)                                             
!                                                                       
      ELSEIF (pl_prog.eq.'cif'.or.pl_prog.eq.'jmol') then 
!                                                                       
!------ - Atom or molecule mode                                         
!                                                                       
         IF (pl_sel_atom) then 
            CALL plot_cif (iff, .TRUE., 'P1') 
         ELSE 
            ier_num = - 76 
            ier_typ = ER_APPL 
         ENDIF 
!                                                                       
!------ Output XBS                                                      
!                                                                       
      ELSEIF (pl_prog.eq.'xbs') then 
!                                                                       
!------ - Atom or molecule mode                                         
!                                                                       
         IF (pl_sel_atom) then 
            CALL plot_xbs (iff) 
         ELSE 
            ier_num = - 76 
            ier_typ = ER_APPL 
         ENDIF 
!                                                                       
!------ Output XBS Movie frames                                         
!                                                                       
      ELSEIF (pl_prog.eq.'frames') then 
!                                                                       
!------ - Atom or molecule mode                                         
!                                                                       
         IF (pl_sel_atom) then 
            CALL plot_frames (iff) 
         ELSE 
            ier_num = - 76 
            ier_typ = ER_APPL 
         ENDIF 
!                                                                       
!------ Output DIAMOND                                                  
!                                                                       
      ELSEIF (pl_prog.eq.'diamond') then 
!                                                                       
!------ - Atom or molecule mode                                         
!                                                                       
         IF (pl_sel_atom) then 
            CALL plot_diamond (iff) 
         ELSE 
            ier_num = - 76 
            ier_typ = ER_APPL 
         ENDIF 
      ENDIF 
!                                                                       
      CLOSE (iff) 
      RETURN 
!                                                                       
      CLOSE (iff) 
      ier_num = - 2 
      ier_typ = ER_IO 
!                                                                       
      END SUBROUTINE do_plot                        
!
!*****7*****************************************************************
!
SUBROUTINE plot_inter(kill)
!   
USE crystal_mod
USE discus_plot_mod
USE discus_plot_init_mod
USE metric_mod
USE symm_mod
USE symm_sup_mod
USE tensors_mod
USE trans_sup_mod
USE trafo_mod
!
USE spcgr_apply
!
USE envir_mod
USE errlist_mod
USE param_mod
USE precision_mod
USE prompt_mod
USE trig_degree_mod
USE matrix_mod
USE support_mod
!
IMPLICIT NONE
!
REAL, PARAMETER :: azero = 1.0
REAL, PARAMETER :: eps   = 1.0E-5
CHARACTER(LEN=*), INTENT(IN) :: kill
!
INTEGER, PARAMETER  :: ITMP     = 78  ! temporary unit number
LOGICAL, PARAMETER  :: lscreen = .FALSE.
CHARACTER(LEN=PREC_STRING) :: tempfile
CHARACTER(LEN=PREC_STRING) :: line
CHARACTER(LEN=PREC_STRING) :: line_move
INTEGER             :: i,j,k
INTEGER             :: length
LOGICAL             :: lsecond  = .false.
LOGICAL             :: lunit    = .false.
INTEGER, SAVE       :: jmol_no  = 0   ! Jmol instance number
REAL, DIMENSION(4)  :: hkl      !Normal in reciprocal space for augmented matrix
REAL, DIMENSION(3)  :: n_hkl    !Normal in reciprocal space
REAL, DIMENSION(3)  :: n_uvw    !Normal in direct space
REAL, DIMENSION(3)  :: v        !temporary vector
REAL, DIMENSION(3)  :: u        !temporary vector
REAL(KIND=PREC_DP), DIMENSION(3)  :: axis     !axis for moveto command
REAL, DIMENSION(3)  :: p_a      !projected abscissa
REAL, DIMENSION(3)  :: p_o      !projected abscissa
REAL(KIND=PREC_DP), DIMENSION(3,3):: roti     !Rotation matrix
REAL(KIND=PREC_DP), DIMENSION(3,3):: rotf     !Rotation matrix
REAL(KIND=PREC_DP), DIMENSION(3,3):: test     !Test for unit matrix
REAL(KIND=PREC_DP) :: det, trace
INTEGER :: ier
REAL                :: rr
REAL(KIND=PREC_DP)                :: beta
!
IF(pl_prog=='jmol') THEN
   IF(kill=='yes') THEN
      CALL jmol_kill(.TRUE., .FALSE.)
   ENDIF
!
! Determine moveto command
!
   CALL plot_ini_jmol  (azero,                        &
        pl_tran_g, pl_tran_gi, pl_tran_f, pl_tran_fi, &
        cr_gten, cr_rten, cr_eps)
!
   hkl(1:3) = pl_hkl(1:3)                      ! Take HKL of normal as this is c*, if default
   hkl(4)   = 0.0
   CALL tran_ca (hkl, pl_tran_fi, lscreen)     ! TRAN_CA works on 4x4 matrix
   n_hkl(1:3) = hkl(1:3)                       ! TRANS works on 3x3 matrix
   CALL trans (n_hkl, cr_rten, n_uvw, 3)       ! uvw is the normal in caresian space
   rr = SQRT(n_uvw(1)**2 + n_uvw(2)**2 + n_uvw(3)**2)
   n_uvw(1) = n_uvw(1) / rr
   n_uvw(2) = n_uvw(2) / rr
   n_uvw(3) = n_uvw(3) / rr
!
!!!! THE FOLLOWING DOES NOT WORK RIGHT NOW !!! PATCHED BELOW
   WRITE(line,'(6(f12.6,'',''),a)') pl_abs(1:3), n_uvw(1:3), 'dddd'
   length = 82
   CALL do_proj (line, length)                ! Project abscissa into plane normal to 'normal'
   p_a(1:3) = res_para(4:6)
   rr = SQRT(p_a(1)**2 + p_a(2)**2 + p_a(3)**2)
   p_a(1) = p_a(1) / rr
   p_a(2) = p_a(2) / rr
   p_a(3) = p_a(3) / rr
!
!  Do vector product of normal and projected abscissa to get an ordinate
!
   p_o(1) = n_uvw(2)*p_a(3) - n_uvw(3)*p_a(2)
   p_o(2) = n_uvw(3)*p_a(1) - n_uvw(1)*p_a(3)
   p_o(3) = n_uvw(1)*p_a(2) - n_uvw(2)*p_a(1)
!
! write(*,*) ' normal ', n_uvw(:), SQRT(n_uvw(1)**2 + n_uvw(2)**2 + n_uvw(3)**2)
! write(*,*) ' abscis ', p_a  (:), SQRT(p_a  (1)**2 + p_a  (2)**2 + p_a  (3)**2)
! write(*,*) ' ordina ', p_o  (:), SQRT(p_o  (1)**2 + p_o  (2)**2 + p_o  (3)**2)
! write(*,*) ' w_na   ', acosd(n_uvw(1)*p_a(1) + n_uvw(2)*p_a(2) + n_uvw(3)*p_a(3) )
! write(*,*) ' w_no   ', acosd(n_uvw(1)*p_o(1) + n_uvw(2)*p_o(2) + n_uvw(3)*p_o(3) )
! write(*,*) ' w_ao   ', acosd(p_a  (1)*p_o(1) + p_a  (2)*p_o(2) + p_a  (3)*p_o(3) )
   roti(:,1) =  p_a  (:)
   roti(:,2) =  p_o  (:)
   roti(:,3) = n_uvw(:)
   CALL matinv3(roti,rotf)
!
!  The second operation will be the tilt of the intended normal 
   axis(1)   = 0.0
   axis(2)   = 0.0
   axis(3)   = 1.0
   test(:,:) = rotf(:,:)
   test(1,1) = ABS(test(1,1)) - 1.0
   test(2,2) = ABS(test(2,2)) - 1.0
   test(3,3) = ABS(test(3,3)) - 1.0
   lunit = .TRUE.
   DO i=1,3
      DO j=1,3
         lunit = lunit .AND. ABS(test(i,j)) < eps
      ENDDO
   ENDDO
   IF(lunit) THEN                            ! We have a unit matrix
      line_move = 'moveto 0.0  {0 0 1}  0.0'
   ELSE                                      ! Regular matrix determne axis amd angle
      CALL mat_axis(rotf, axis, beta, det, trace, ier)
!      CALL get_detail_axis(rotf, w_char, axis)
      rr = SQRT(axis(1)**2 + axis(2)**2 + axis(3)**2)
      axis(1) = axis(1) / rr
      axis(2) = axis(2) / rr
      axis(3) = axis(3) / rr
!
! write(*,*) ' MATRIX              ', rotf(1,:)
! write(*,*) ' MATRIX              ', rotf(2,:)
! write(*,*) ' MATRIX              ', rotf(3,:)
! write(*,*) ' AXIS                ', axis
      v(1) = rotf(1,1)*n_uvw(1) + rotf(1,2)*n_uvw(2) + rotf(1,3)*n_uvw(3)
      v(2) = rotf(2,1)*n_uvw(1) + rotf(2,2)*n_uvw(2) + rotf(2,3)*n_uvw(3)
      v(3) = rotf(3,1)*n_uvw(1) + rotf(3,2)*n_uvw(2) + rotf(3,3)*n_uvw(3)
! write(*,*) ' TRANSFORMED normal  ', v
      u(1) = rotf(1,1)*p_a  (1) + rotf(1,2)*p_a  (2) + rotf(1,3)*p_a  (3)
      u(2) = rotf(2,1)*p_a  (1) + rotf(2,2)*p_a  (2) + rotf(2,3)*p_a  (3)
      u(3) = rotf(3,1)*p_a  (1) + rotf(3,2)*p_a  (2) + rotf(3,3)*p_a  (3)
! write(*,*) ' TRANSFORMED abscissa', u
      u(1) = rotf(1,1)*p_o  (1) + rotf(1,2)*p_o  (2) + rotf(1,3)*p_o  (3)
      u(2) = rotf(2,1)*p_o  (1) + rotf(2,2)*p_o  (2) + rotf(2,3)*p_o  (3)
      u(3) = rotf(3,1)*p_o  (1) + rotf(3,2)*p_o  (2) + rotf(3,3)*p_o  (3)
! write(*,*) ' TRANSFORMED ordinate', u
!     We need negative angle as view is into -C in JMOL !!
      beta = -ACOSD(0.5*(rotf(1,1) + rotf(2,2) + rotf(3,3)-1.)) ! trace is cos+cos+1
! write(*,*) ' Determinant, beta   ', det3(rotf), beta
      IF(ABS(beta)<EPS) THEN
         line_move = 'moveto 0.0  {0 0 1}  0.0'
      ELSE
         WRITE(line_move,'(a,3(1x,f6.3),a,f8.3)') 'moveto 0.0  {', axis(:),'}  ', beta
      ENDIF
   ENDIF
!!!!!!!!!!!PATCH FOR RIGHT NOW!!!!!!!!!!!!
! Just straighten up the normal
   IF(ABS(n_uvw(3)-1.0)<EPS) THEN
      line_move = 'moveto 0.0  {0 0 1}  0.0'
   ELSE
      beta = ACOSD(n_uvw(3))
      axis(1) =  n_uvw(2)
      axis(2) = -n_uvw(1)
      axis(3) =  0.0
      WRITE(line_move,'(a,3(1x,f6.3),a,f8.3)') 'moveto 0.0  {', axis(:),'}  ', beta
   ENDIF
! write(*,*) 'LINEMOVE ', line_move(1:LEN_TRIM(line_move))
!! END !!!!PATCH FOR RIGHT NOW!!!!!!!!!!!!
   CALL plot_ini_trans (azero,                        &
        pl_tran_g, pl_tran_gi, pl_tran_f, pl_tran_fi, &
        cr_gten, cr_rten, cr_eps)
!
   jmol_no = jmol_no + 1
   IF(operating=='Linux') THEN
      WRITE(tempfile,'(a,i10.10,a,i10.10,a)')  '/tmp/jmol.',PID,'.',jmol_no,'.mol'
   ELSEIF(operating=='Linux_WSL') THEN
      WRITE(tempfile,'(a,i10.10,a,i10.10,a)')  '/tmp/jmol.',PID,'.',jmol_no,'.mol'
   ELSEIF(operating(1:6)=='darwin') THEN
      WRITE(tempfile,'(a,i10.10,a,i10.10,a)')  '/tmp/jmol.',PID,'.',jmol_no,'.mol'
   ELSEIF(operating(1:6)=='cygwin') THEN
      WRITE(tempfile,'(a,i10.10,a,i10.10,a)')  '/tmp/jmol.',PID,'.',jmol_no,'.mol'
   ELSEIF(operating(1:7)=='Windows') THEN
      WRITE(tempfile,'(a,i10.10,a,i10.10,a)')  '/tmp/jmol.',PID,'.',jmol_no,'.mol'
   ENDIF
   CALL oeffne( ITMP, tempfile, 'unknown')
   IF(operating=='Linux'.OR.operating(1:6)=='darwin' .OR. &
      operating=='Linux_WSL'                             ) THEN
!
!     Write load command, enclose file name in apostrophes to cover 
!     blanks in filenemes
      i=LEN_TRIM(current_dir)
      IF(current_dir(i:i) /= '/') current_dir(i+1:i+1)='/'
      WRITE(ITMP,'(3a)') 'load ''',                  &
         current_dir(1:LEN_TRIM(current_dir))//      &
         pl_out(1:LEN_TRIM(pl_out)), ''' {1 1 1}'
   ELSEIF(operating(1:6)=='cygwin') THEN
      IF(current_dir(1:5)=='/home') THEN
      WRITE(ITMP,'(3a)') 'load ''',                   &
         operat_top(1:LEN_TRIM(operat_top))//        &
         current_dir(1:LEN_TRIM(current_dir))//      &
         pl_out(1:LEN_TRIM(pl_out)), ''' {1 1 1}'
      ELSEIF(current_dir(1:10)=='/cygdrive/') THEN
      WRITE(ITMP,'(3a)') 'load ''',                   &
         'C:\'//                                      &
         current_dir(13:LEN_TRIM(current_dir))//      &
         pl_out(1:LEN_TRIM(pl_out)), ''' {1 1 1}'
      ELSE
         WRITE(*,*) 'WRONG DRIVE ', current_dir(1:len_trim(current_dir))
      ENDIF
   ELSEIF(operating(1:7)=='Windows') THEN
      WRITE(ITMP,'(3a)') 'load ''',                  &
         current_dir(11:11) // ':' //                &
         current_dir(12:LEN_TRIM(current_dir)) //    &
         pl_out(1:LEN_TRIM(pl_out)), ''' {1 1 1}'
   ENDIF
!  WRITE(ITMP,'(a)') 'moveto 0.0  {0 0 1}  0.0'
   WRITE(ITMP,'(a)') line_move(1:len_TRIM(line_move))
   WRITE(ITMP,'(a)') 'boundbox off'
   WRITE(ITMP,'(a)') 'unitcell off'
   WRITE(ITMP,'(a)') 'axes 0.03'
   WRITE(ITMP,'(a)') 'hide none'
   WRITE(ITMP,'(a)') 'select * '
   WRITE(ITMP,'(a)') 'center selected'
   IF(ANY(pl_bond(0:cr_nscat, 0:cr_nscat))) THEN   ! Bonds were specified
   WRITE(ITMP,'(a)') 'connect (all) (all) DELETE'
   DO i=0, cr_nscat
      DO j=0, cr_nscat
         IF(pl_bond(i,j)) THEN
            WRITE(ITMP,'(a, f5.2, 1x, f5.2, 1x, a,a,a, a,a,a, a, f5.2, a)') &
            'connect ', pl_bond_len(1,i,j), pl_bond_len(2,i,j), &
            '(_',cr_at_lis(i)(1:len_trim(cr_at_lis(i))),')',    &
            '(_',cr_at_lis(j)(1:len_trim(cr_at_lis(j))),')',    &
            ' SINGLE radius ',pl_bond_rad(i,j), ' ModifyOrCreate'
         ENDIF
      ENDDO
   ENDDO
   ENDIF
   IF(pl_poly_n>0) THEN
!  DO i=1, pl_poly_n
      WRITE(line,'(a)') 'polyhedra '
      k = LEN_TRIM(line)
      IF(pl_poly_nmin>0 .AND. pl_poly_nmax>0) WRITE(line(k+1:k+3),'(i2,a)') pl_poly_nmin,' '
      k = LEN_TRIM(line)
      IF(pl_poly_nmax>0) WRITE(line(k+1:k+3),'(1x,i2  )') pl_poly_nmax
      k = LEN_TRIM(line)
      IF(pl_poly_dmin>0.0 .AND. pl_poly_dmax>0.0) WRITE(line(k+1:k+8),'(f8.2)') pl_poly_dmin
      k = LEN_TRIM(line)
      IF(pl_poly_dmax>0.0) WRITE(line(k+1:k+8),'(f8.2)') pl_poly_dmax
      k = LEN_TRIM(line)
      WRITE(line(k+1:k+8),'(a)') ' BONDS ('
      k = LEN_TRIM(line)
      IF(pl_poly_c(-1)) THEN
         WRITE(line(k+1:k+7),'(a)') '*) to ('
         k = LEN_TRIM(line)
      ELSE
         lsecond = .FALSE.
         DO j=0,cr_nscat
            IF(pl_poly_c(j)) THEN
               IF(lsecond) THEN
                  WRITE(line(k+1:k+8),'(a,a)') ' OR ',cr_at_lis(j)
               ELSE
                  WRITE(line(k+1:k+4),'(a)') cr_at_lis(j)
                  lsecond = .TRUE.
               ENDIF
               k = LEN_TRIM(line)
            ENDIF
         ENDDO
         WRITE(line(k+1:k+6),'(a)') ') to ('
         k = LEN_TRIM(line)
      ENDIF
      IF(pl_poly_o(-1)) THEN
         WRITE(line(k+1:k+2),'(a)') '*)'
         k = LEN_TRIM(line)
      ELSE
         lsecond = .FALSE.
         DO j=0,cr_nscat
            IF(pl_poly_o(j)) THEN
               IF(lsecond) THEN
                  WRITE(line(k+1:k+8),'(a,a)') ' OR ',cr_at_lis(j)
               ELSE
                  WRITE(line(k+1:k+4),'(a)') cr_at_lis(j)
                  lsecond = .TRUE.
               ENDIF
               k = LEN_TRIM(line)
            ENDIF
         ENDDO
         WRITE(line(k+1:k+1),'(a)') ')'
         k = LEN_TRIM(line)
      ENDIF
      IF(.NOT.pl_poly_face) WRITE(line(k+1:k+10),'(a)') ' COLLAPSED'
      WRITE(ITMP,'(a)') line(1:LEN_TRIM(line))
      IF(pl_poly_hue) THEN
         line = 'color polyhedra translucent '
         IF(pl_poly_col/='auto') line = line(1:len_trim(line)) // ' ' // pl_poly_col
         WRITE(ITMP,'(a)') line(1:LEN_TRIM(line))
      ELSE
         IF(pl_poly_col/='auto') THEN
            line = 'color polyhedra ' // pl_poly_col
            WRITE(ITMP,'(a)') line(1:LEN_TRIM(line))
         ENDIF
      ENDIF
!  ENDDO
   ENDIF
   WRITE(ITMP,'(a)') 'spacefill 25%'
   WRITE(ITMP,'(a)') 'zoom 125'
   WRITE(ITMP,'(a,i3,a,i3,a,i3,a)') 'background [',pl_back(1),',' , pl_back(2),',', pl_back(3),']'
   CLOSE(ITMP)
   IF(operating=='Linux') THEN
      WRITE(line,'(a,a,a,a)') pl_jmol(1:LEN_TRIM(pl_jmol)), ' -s ',&
            tempfile(1:LEN_TRIM(tempfile)), ' > /dev/null &'
   ELSEIF(operating=='Linux_WSL') THEN
      WRITE(line,'(a,a,a,a)') pl_jmol(1:LEN_TRIM(pl_jmol)), ' -s ',&
            tempfile(1:LEN_TRIM(tempfile)), ' > /dev/null &'
   ELSEIF(operating(1:6)=='darwin') THEN
      WRITE(line,'(a,a,a,a)') pl_jmol(1:LEN_TRIM(pl_jmol)), ' -s ',&
            tempfile(1:LEN_TRIM(tempfile)), ' > /dev/null &'
   ELSEIF(operating(1:6)=='cygwin') THEN
      WRITE(line,'(a,a,i10.10,a,i10.10,a,a,a)') pl_jmol(1:LEN_TRIM(pl_jmol)), &
            ' -s jmol.', PID, '.',jmol_no, '.mol  > /dev/null &'
   ELSEIF(operating(1:7)=='Windows') THEN
      WRITE(line,'(a,a,i10.10,a,i10.10,a,a,a)') pl_jmol(1:LEN_TRIM(pl_jmol)), &
            ' -s ../../tmp/jmol.', PID, '.',jmol_no, '.mol  > /dev/null &'
!
   ENDIF
   WRITE(output_io,'(a)') ' JMOL may take a moment to show up'
   CALL system(line)
ENDIF
END SUBROUTINE plot_inter
!
!*****7*****************************************************************
!
SUBROUTINE jmol_kill(lout, lfinal)
!-
! Kill jmol processes that were started by this discus_suite 
!
USE envir_mod
USE errlist_mod
USE prompt_mod
USE support_mod
!
IMPLICIT NONE
!
LOGICAL, INTENT(IN) :: lout   ! Flag info on screen for Linux
LOGICAL, INTENT(IN) :: lfinal ! Flag info on screen for all systems
!
INTEGER, PARAMETER  :: ITMP     = 78  ! temporary unit number
!
CHARACTER(LEN=PREC_STRING) :: kill_file
CHARACTER(LEN=PREC_STRING) :: line
!
INTEGER             :: jmol_pid   ! jmol Process identifier
INTEGER             :: jppid      ! jmol Parent Process identifier
INTEGER             :: ios
!
LOGICAL :: lpresent
LOGICAL :: did_kill
!
WRITE(kill_file, '(a,I10.10,a)') '/tmp/jmol.',PID,'.pid'
INQUIRE(FILE=kill_file, EXIST=lpresent)
IF(lpresent) THEN
   WRITE(line, '(a,a)') 'rm -f ', kill_file(1:LEN_TRIM(kill_file))
   CALL system(line)
ENDIF
!
IF(lfinal) THEN
   IF(pl_keep) RETURN
ENDIF
!
did_kill = .FALSE.
IF(operating=='Linux') THEN
   WRITE(line,'(a,i10,a,a)') 'ps j | grep ',PID,' | grep java |  grep -v grep | awk ''{print $2, $3}'' >> ', &
         kill_file(1:LEN_TRIM(kill_file))
   CALL system(line)
ELSEIF(operating=='Linux_WSL') THEN
   WRITE(line,'(a,i10,a,a)') 'ps j | grep ',PID,' | grep java |  grep -v grep | awk ''{print $2, $3}'' >> ', &
         kill_file(1:LEN_TRIM(kill_file))
   CALL system(line)
ELSEIF(operating(1:6)=='darwin') THEN
   WRITE(line,'(a,i10,a,a)') 'ps j | grep ',PID,' | grep java |  grep -v grep | awk ''{print $2, $4}'' >> ', &
         kill_file(1:LEN_TRIM(kill_file))
   CALL system(line)

ELSEIF(operating(1:6)=='cygwin' .OR. operating(1:7)=='Windows') THEN
!  WRITE(line,'(a,I5.5,a,a)') 'ps aux | grep ''jmol -s /tmp/jmol.',PID, &
!     '.mol'' | grep -v grep | awk ''{print $1, $3}''  '
!  CALL system(line)
!  WRITE(line,'(a,I5.5,a,a)') 'ps aux | grep ''jmol -s /tmp/jmol.',PID, &
!     '.mol'' | grep -v grep | awk ''{print $1, $3}'' > ', kill_file(1:LEN_TRIM(kill_file))
!  CALL system(line)
!  line = 'ps aux | grep java |                      grep -v grep | awk ''{print $1, $3}'' '
!  CALL system(line)
   line = 'ps aux | grep java |                      grep -v grep | awk ''{print $1, $3}'' >> ' //  &
         kill_file(1:LEN_TRIM(kill_file))
   CALL system(line)
ENDIF
!
CALL oeffne( ITMP, kill_file, 'old')
IF(ier_num==0) THEN
   READ(ITMP,*,IOSTAT=ios) jmol_pid, jppid
   DO WHILE (.NOT.IS_IOSTAT_END(ios)) 
      IF(jppid==PID .OR. jppid==PPID) THEN                    ! Current discus_suite has started jmol
         WRITE(line,'(a,i12,a)') 'kill -9 ',jmol_pid, ' > /dev/null'
         CALL system(line)
         did_kill = .TRUE.
      ENDIF
      READ(ITMP,*,IOSTAT=ios) jmol_pid, jppid
   ENDDO
   CLOSE(ITMP)
!!!RBN   WRITE(line, '(a,a)') 'rm -f ', kill_file(1:LEN_TRIM(kill_file))
!   CALL system(line)
   IF(operating=='Linux' .AND. lout) THEN
      IF(did_kill) WRITE(output_io,*) ' Processes were killed, hit ENTER to retrieve prompt'
   ENDIF
ENDIF
!
WRITE(line, '(a,a)') 'rm -f ', kill_file(1:LEN_TRIM(kill_file))
CALL system(line)
!
IF(did_kill .AND.lfinal) THEN
   WRITE(output_io,*) ' Closed JMOL windows, ignore ''Killed'' messages '
   WRITE(output_io,*) ' DISCUS_SUITE will terminate properly '
ENDIF
!
END SUBROUTINE jmol_kill
!
!*****7*****************************************************************
!
SUBROUTINE plot_test_aon(lnor, labs, lord)
!
USE crystal_mod
USE metric_mod
USE discus_plot_mod 
USE trafo_mod
!
USE errlist_mod
!
IMPLICIT NONE
!
LOGICAL, INTENT(IN) :: lnor
LOGICAL, INTENT(IN) :: labs
LOGICAL, INTENT(IN) :: lord
!
LOGICAL           , PARAMETER :: lspace = .TRUE.
REAL              , PARAMETER :: EPS  = 1.E-5!A small value
REAL, DIMENSION(3), PARAMETER :: NULL = (/0.00, 0.00, 0.00/)
!
REAL               :: w_na, w_no, w_ao    ! pairwise angles
REAL, DIMENSION(3) :: v                   ! a vector
!  Tests if normal, abscissa nad ordinate are parallel, try to correct
w_na = do_bang (lspace, pl_uvw, NULL, pl_abs) 
w_no = do_bang (lspace, pl_uvw, NULL, pl_ord) 
w_ao = do_bang (lspace, pl_abs, NULL, pl_ord) 
IF(ABS(w_na)<EPS)  THEN     ! Normal and Abscissa are parallel
  IF(lnor) THEN            ! User specified Normal
     IF(labs) THEN         ! User specified Abscissa, Flag error
        ier_num = -80 
        ier_typ = ER_APPL
        RETURN
     ELSE                  ! No abscissa give, lets calculate one
        IF(ABS(w_no) > EPS) THEN ! Normal and Ordinate are not parallel
           CALL vekprod (pl_ord, pl_uvw, pl_abs, cr_eps, cr_rten)
           w_na = 90.0
           w_ao = 90.0
        ELSE                     ! Normal and Ordinate are parallel
           IF(labs) THEN         ! User specified Ordinate, Flag error
              ier_num = -153
              ier_typ = ER_APPL
              RETURN
           ELSE                  ! No abscissa give, lets calculate one
              v(1) = pl_uvw(1) + 0.1
              v(2) = pl_uvw(2) + 0.01
              v(3) = pl_uvw(3) + 0.03
              CALL vekprod (v     , pl_uvw, pl_abs, cr_eps, cr_rten)
              CALL vekprod (pl_uvw, pl_abs, pl_ord, cr_eps, cr_rten)
              w_na = 90.0
              w_no = 90.0
              w_ao = 90.0
           ENDIF
        ENDIF
     ENDIF
  ELSE                           ! User did not give a normal
     IF(labs) THEN               ! User specified Abscissa
        IF(lord) THEN            ! User specified Ordinate
           CALL vekprod (pl_abs, pl_ord, pl_uvw, cr_eps, cr_rten)
           w_na = 90.0
           w_no = 90.0
        ELSE                        ! User did not give Ordinate
           v(1) = pl_abs(1) + 0.1
           v(2) = pl_abs(2) + 0.01
           v(3) = pl_abs(3) + 0.03
           CALL vekprod (v     , pl_abs, pl_uvw, cr_eps, cr_rten)
           CALL vekprod (pl_uvw, pl_abs, pl_ord, cr_eps, cr_rten)
           w_na = 90.0
           w_no = 90.0
           w_ao = 90.0
        ENDIF
     ELSE                        ! User did not give abscissa
        IF(lord) THEN            ! User specified Ordinate
           v(1) = pl_ord(1) + 0.1
           v(2) = pl_ord(2) + 0.01
           v(3) = pl_ord(3) + 0.03
           CALL vekprod (pl_ord, v     , pl_abs, cr_eps, cr_rten)
           CALL vekprod (pl_abs, pl_ord, pl_uvw, cr_eps, cr_rten)
           w_na = 90.0
           w_no = 90.0
           w_ao = 90.0
        ELSE                        ! User did not give ordinate
           pl_hkl(1) = 0.0          ! Make JMOL standard
           pl_hkl(2) = 0.0
           pl_hkl(3) = 1.0
           pl_abs(1) = 1.0
           pl_abs(2) = 0.0
           pl_abs(3) = 0.0
           CALL trans (pl_hkl, cr_rten, pl_uvw, 3) 
           CALL vekprod (pl_uvw, pl_abs, pl_ord, cr_eps, cr_rten)
           w_na = 90.0
           w_no = 90.0
           w_ao = 90.0
        ENDIF
     ENDIF
  ENDIF
ENDIF
!   ! at this point w_na is not ZERO, check Ordinate and correct
IF(ABS(w_no)<EPS .OR. ABS(w_ao)<EPS)  THEN     ! Ordinate is parallel to Normal or Abscissa 
   IF(lord) THEN            ! User specified Ordinate
      ier_num = -154
      ier_typ = ER_APPL
      RETURN
   ELSE                     ! User did not give ordinate
      CALL vekprod (pl_uvw, pl_abs, pl_ord, cr_eps, cr_rten)
      w_no = 90.0
      w_ao = 90.0
   ENDIF
ENDIF
END SUBROUTINE plot_test_aon
!
!*******************************************************************************
!
SUBROUTINE plot_reset
!
! reset all variables needed for the plotting of structures
!-
!
USE discus_allocate_appl_mod
USE discus_plot_mod
!
IMPLICIT NONE
!
INTEGER :: ik
!
CALL alloc_plot(1, 1)
!
pl_init  = .TRUE.
pl_jmol  = ' '
pl_out   = 'plot.cif'
pl_title = 'crystal structure'
pl_prog  = 'cif'
pl_vrml  = 'u'
pl_col   = 'xyz'
pl_width     = 1.e12
pl_dim       = reshape((/(-1e14,ik=1,3),(1e12,ik=1,3)/),shape(pl_dim)) ! (3,2)
pl_hkl       = (/0.0,0.0,1.0/)
pl_uvw       = (/0.0,0.0,1.0/)
pl_vec       = 0.0
pl_abs       = (/1.0,0.0,0.0/)
pl_ord       = (/0.0,1.0,0.0/)
pl_siz(:)    = 1.0      ! 
pl_rgb(:,:)  = 0.0       ! 
pl_rgb(1,:)  = 1.0       ! 
pl_back      = (/ 240, 240, 240 /)   ! plot background
pl_bond_len(:,:,:) = 0.0  ! 
pl_bond_rad(  :,:) = 0.2  !
pl_bond_col(:,:,:) = 0.0  !
pl_vrml_scaling = 0.05
pl_scale        = (/-1.0, 1.0/)
pl_sel_prop     = (/0,0/)
pl_typ(:)       = 3
pl_color(:)     = 0
pl_latom(:)     = .FALSE. 
pl_lsite(:)     = .TRUE.  
pl_batom_a(:)   = .FALSE.   ! 
pl_batom_e(:)   = .FALSE.   ! 
pl_poly_n       = 0  ! Number of polyhedra definitions
pl_poly_c(:) = .FALSE.    ! 
pl_poly_o(:) = .FALSE.    ! 
pl_poly_dmin = 0.0     ! Minimum distance to neighbor for polhedra
pl_poly_dmax = 0.0     ! Maximum distance to neighbor for polhedra
pl_poly_nmin = 0       ! Minimum neighbors for polhedra
pl_poly_nmax = 192     ! Maximum neighbors for polhedra
pl_poly_face = .TRUE.  ! Face style flat/collapsed
pl_poly_hue  = .FALSE. ! Face style solid / transparent
pl_poly_col  = 'auto'  ! Face color
pl_dens      = .FALSE.
pl_sel_atom  = .TRUE.
pl_mol_all   = .TRUE.
pl_bond(:,:) = .FALSE.      ! (0:MAXSCAT,0:MAXSCAT)
pl_append    = .FALSE.
pl_ext_all   = .TRUE.
!
END SUBROUTINE plot_reset
!
END MODULE discus_plot_menu
