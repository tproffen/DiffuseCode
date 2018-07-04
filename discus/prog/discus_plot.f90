MODULE discus_plot_menu
!
CONTAINS
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
      USE class_macro_internal
      USE prompt_mod 
      USE string_convert_mod
      USE sup_mod
      USE take_param_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER, PARAMETER :: ITMP     = 78  ! temporary unit number
      INTEGER, PARAMETER :: MIN_PARA = 20  ! A command requires at leaset these no of parameters
      INTEGER maxw 
      LOGICAL lold 
      PARAMETER (lold = .false.) 
!                                                                       
      CHARACTER(LEN=1024), DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: cpara ! (MAX(10,MAXSCAT)) 
      REAL               , DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: werte ! (MAX(10,MAXSCAT)) 
      INTEGER            , DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: lpara ! (MAX(10,MAXSCAT))
      CHARACTER(1024) line, zeile
      CHARACTER(LEN=1024) :: tempfile
      CHARACTER(LEN=LEN(prompt)) :: orig_prompt 
      CHARACTER(5) befehl 
      CHARACTER(1) cdum 
      REAL :: size, rr=0.0, rg=0.0, rb=0.0
      INTEGER lp, length 
      INTEGER :: ianz, i, j, is, it, ic, lbef 
!      INTEGER :: npoly     
      INTEGER indxg 
      INTEGER         :: nscat
      INTEGER         :: ios

      LOGICAL lend, l_select 
      LOGICAL :: labs = .FALSE.
      LOGICAL :: lord = .FALSE.
      LOGICAL :: lnor = .FALSE.
!
!                                                                       
      INTEGER len_str 
      LOGICAL str_comp 
!                                                                       
      INTEGER, PARAMETER :: NOPTIONAL = 9
      CHARACTER(LEN=1024), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
      CHARACTER(LEN=1024), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
      INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
      INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
      REAL               , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
      INTEGER, PARAMETER                        :: ncalc = 4 ! Number of values to calculate 
!
      DATA oname  / 'dmin' , 'dmax' , 'nmin' , 'nmax' , 'face' , 'hue'  , 'color', 'plot' , 'kill'  /
      DATA loname /  4     ,  4     ,  4     ,  4     ,  4     ,  3     ,  4     ,  4     ,  4      /
!
      ! Always provide fresh default values, Repeat for each command
      opara  =   (/ '0.000', '0.000', '0    ', '0    ', 'flat ', 'solid', 'auto ', 'none ', 'none ' /)
      lopara =   (/  5     ,  5     ,  1     ,  1     ,  4     ,  5     ,  4     ,  4     ,  4      /)
      owerte =   (/  0.00  ,  0.00  ,  0.    ,  0.    ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0    /)
!
!
!                                                                       
      maxw = MAX(MIN_PARA,MAXSCAT+1)
      lend = .false. 
!
!     Allocate the necessary arrays
!
      IF ( cr_nscat > PL_MAXSCAT .or. mole_num_type > PL_MAXSCAT .OR. & 
           MAXSCAT > PL_MAXSCAT ) THEN
         nscat = max ( cr_nscat, mole_num_type, MAXSCAT)
         CALL alloc_plot ( nscat )
         IF(ier_num < 0) THEN
           RETURN
         ENDIF
      ENDIF
!                                                                       
      orig_prompt = prompt
      prompt = prompt (1:len_str (prompt) ) //'/plot' 
      DO while (.not.lend) 
      CALL no_error 
      CALL get_cmd (line, length, befehl, lbef, zeile, lp, prompt) 
      IF (ier_num.eq.0) then 
         IF (line /= ' '      .and. line(1:1) /= '#' .and. &
             line /= char(13) .and. line(1:1) /= '!'        ) THEN
!                                                                       
!     ----search for "="                                                
!                                                                       
            indxg = index (line, '=') 
      IF (indxg /= 0.AND..NOT. (str_comp (befehl, 'echo', 2, lbef, 4) ) &
                    .AND..NOT. (str_comp (befehl, 'syst', 2, lbef, 4) ) &
                    .AND..NOT. (str_comp (befehl, 'help', 2, lbef, 4)   &
                    .OR.        str_comp (befehl, '?   ', 2, lbef, 4) ) &
                    .AND. INDEX(line,'==') == 0                        ) THEN
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
!     ----Select the abscissa for a projection onto a to plot           
!           slice 'absc'                                                
!                                                                       
               ELSEIF (str_comp (befehl, 'absc', 2, lbef, 4) ) THEN 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0.and.ianz.eq.3) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        DO i = 1, 3 
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
!     --continues a macro 'continue'                                    
!                                                                       
               ELSEIF (str_comp (befehl, 'continue', 3, lbef, 8) ) then 
                  CALL macro_continue (zeile, lp) 
!                                                                       
!------ --setting the scale of the plot                                 
!                                                                       
               ELSEIF (str_comp (befehl, 'scale', 3, lbef, 5) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        IF (werte (1) .eq. - 1) then 
                           pl_scale (0) = - 1 
                        ELSEIF (werte (1) .gt.0.0) then 
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
                     pl_sel_atom, lold,  .true.)
                     zeile = cpara (2) (1:lpara (2) ) 
                     CALL atom_select (zeile, lp, 0, MAXSCAT, pl_batom_e, &
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
                                    oname, loname, opara, lopara, owerte)
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
                        IF(line(1:1)=='''' .AND. line(i:i)=='''') THEN
                           zeile = line(2:i-1)
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
                        IF(operating=='Linux') THEN                      
                           tempfile = '/tmp/which_jmol'
                           line = 'which jmol > /tmp/which_jmol'
                        ELSEIF(operating=='Windows') THEN                      
                           tempfile = '/tmp/which_jmol'
                           line = 'which java > /tmp/which_jmol'
                        ELSEIF(operating(1:6)=='cygwin') THEN                      
                           tempfile = '/tmp/which_jmol'
                           line = 'which java > /tmp/which_jmol'
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
                                 IF(operating(1:6)=='cygwin' .OR.       &
                                    operating(1:7)=='Windows'    ) THEN
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
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ENDIF 
!                                                                       
!     ----run plot 'run'                                                
!                                                                       
               ELSEIF (str_comp (befehl, 'run', 1, lbef, 3) ) then 
!                 Always provide fresh default values, Repeat for each command
                  opara (8:9) =   (/ 'none ', 'none ' /)
                  lopara(8:9) =   (/  4     ,  4      /)
                  owerte(8:9) =   (/  0.0   ,  0.0    /)
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                                    oname, loname, opara, lopara, owerte)
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
!------  -----waiting for user input                                    
!                                                                       
               ELSEIF (str_comp (befehl, 'wait', 3, lbef, 4) ) then 
                  CALL do_input (zeile, lp) 
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
      ENDDO 
!                                                                       
      prompt = orig_prompt
      END SUBROUTINE plot                           
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
      USE errlist_mod 
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
            CALL oeffne_append (iff, pl_out, 'old', lread) 
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
            CALL plot_cif (iff) 
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
!*****7*****************************************************************
      SUBROUTINE plot_atoms (iff) 
!-                                                                      
!     Writes the selected atoms for plotting with ATOMS                 
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE metric_mod
      USE modify_func_mod
      USE discus_plot_mod 
      USE tensors_mod
      USE errlist_mod 
      USE wink_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      REAL d, dist, fac 
      REAL v (3) 
      INTEGER i, j, iff 
      LOGICAL lno_slice, latom 
!                                                                       
      INTEGER len_str 
!                                                                       
      WRITE (iff, 500) cr_name (1:len_str (cr_name) ) 
      WRITE (iff, 510) (cr_a0 (i) * cr_icc (i), i = 1, 3),  &
                       (cr_win (i), i = 1, 3)
      WRITE (iff, 520) 'P1' 
      WRITE (iff, 530) 
!                                                                       
      latom = .false. 
      fac = 1.0 / (8.0 * REAL(pi)**2) 
!                                                                       
      IF (pl_hkl(1) .eq.0.and.pl_hkl(2) .eq.0.and.pl_hkl(3) .eq.0) then
         lno_slice = .true. 
         d = 1.0 
      ELSE 
         lno_slice = .false. 
         d = sqrt (skalpro (pl_uvw, pl_uvw, cr_gten) ) 
         DO j = 1, 3 
         pl_mat (j, 1) = pl_abs (j) 
         pl_mat (j, 2) = pl_ord (j) 
         pl_mat (j, 3) = pl_uvw (j) 
         ENDDO 
         CALL invmat (pl_inv, pl_mat) 
         IF (ier_num.eq. - 1) then 
            RETURN 
         ENDIF 
      ENDIF 
!                                                                       
      DO i = 1, cr_natoms 
!                                                                       
!     --Select atom if:                                                 
!       type has been selected and                                      
!            all atoms are selected                            or       
!                                                                       
      IF (check_select_status (i, pl_latom (cr_iscat (i) ), cr_prop (i),   &
      pl_sel_prop) ) then                                               
!                                                                       
!     --Check dimensions of plotting space                              
!                                                                       
         IF (pl_dim (1, 1) .le.cr_pos (1, i) .and.     &
             cr_pos (1, i) .le.pl_dim (1, 2) .and.     &
             pl_dim (2, 1) .le.cr_pos (2, i) .and.     &
             cr_pos (2, i) .le.pl_dim (2, 2) .and.     &
             pl_dim (3, 1) .le.cr_pos (3, i) .and.     &
             cr_pos (3, i) .le.pl_dim (3, 2) ) then
!                                                                       
!     ------Determine distance to point in plot slice                   
!                                                                       
            DO j = 1, 3 
            v (j) = cr_pos (j, i) - pl_vec (j) 
            ENDDO 
            dist = abs (skalpro (pl_uvw, v, cr_gten) ) / d 
!                                                                       
!     ------Write atom if inside slice or if slice has been deselected  
!                                                                       
            IF (dist.le.pl_width.or.lno_slice) then 
!                                                                       
!     ------write atom position or projection onto abscissa/ordinate    
!                                                                       
               DO j = 1, 3 
               v (j) = cr_pos (j, i) 
               ENDDO 
!                                                                       
!     ------Write atom position in desired sequence of coordinates      
!                                                                       
               latom = .true. 
               WRITE (iff, 1000) cr_at_lis (cr_iscat (i) ),          &
                 cr_iscat (i), (v (j) / cr_icc (j), j = 1, 3),       &
                 fac * cr_dw (cr_iscat (i) ), 0.0, 0.0, 0.0, 0.0, 0.0                            
            ENDIF 
         ENDIF 
      ENDIF 
      ENDDO 
!                                                                       
      IF (.not.latom) then 
         ier_num = - 58 
         ier_typ = ER_APPL 
      ENDIF 
!                                                                       
  500 FORMAT    ('TITL ',a) 
  510 FORMAT    ('CELL ',6(f10.6,1x)) 
  520 FORMAT    ('SPGP ',a) 
  530 FORMAT    ('FIELDS LAB TYP COO TFU') 
 1000 FORMAT     (a4,1x,i4,3x,3(f11.6,1x),/,6(f8.6,1x)) 
      END SUBROUTINE plot_atoms                     
!*****7*****************************************************************
      SUBROUTINE plot_drawxtl (iff) 
!-                                                                      
!     Writes the selected atoms for plotting with DRAWxtl               
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE metric_mod
      USE modify_func_mod
      USE discus_plot_mod 
      USE tensors_mod
      USE errlist_mod 
      USE wink_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      REAL d, dist, fac 
      REAL v (3) 
      INTEGER i, j, iff 
      LOGICAL lno_slice, latom 
!                                                                       
      INTEGER len_str 
      CHARACTER(12) povcolor (15) 
      DATA povcolor / 'Scarlet', 'HuntersGreen', 'MediumBlue',          &
      'Magenta', 'Yellow', 'Black', 'IndianRed', 'DarkGreen',           &
      'MidnightBlue', 'Maroon', 'Gold', 'Gray20', 'Cyan', 'SkyBlue',    &
      'White' /                                                         
!                                                                       
      WRITE (iff, 501) cr_name (1:len_str (cr_name) ) 
      WRITE (iff, 511) (cr_a0 (i) * cr_icc (i), i = 1, 3),              &
                       (cr_win (i), i = 1, 3)
      WRITE (iff, 521) 'P 1' 
      WRITE (iff, 531) 
!                                                                       
      latom = .false. 
      fac = 1.0 / (8.0 * REAL(pi)**2) 
!                                                                       
      IF (pl_hkl(1) .eq.0.and.pl_hkl(2) .eq.0.and.pl_hkl(3) .eq.0) then
         lno_slice = .true. 
         d = 1.0 
      ELSE 
         lno_slice = .false. 
         d = sqrt (skalpro (pl_uvw, pl_uvw, cr_gten) ) 
         DO j = 1, 3 
         pl_mat (j, 1) = pl_abs (j) 
         pl_mat (j, 2) = pl_ord (j) 
         pl_mat (j, 3) = pl_uvw (j) 
         ENDDO 
         CALL invmat (pl_inv, pl_mat) 
         IF (ier_num.eq. - 1) then 
            RETURN 
         ENDIF 
      ENDIF 
!                                                                       
      DO i = 1, cr_natoms 
!                                                                       
!     --Select atom if:                                                 
!       type has been selected and                                      
!            all atoms are selected                            or       
!                                                                       
      IF (check_select_status (i, pl_latom (cr_iscat (i) ), cr_prop (i),   &
      pl_sel_prop) ) then                                               
!                                                                       
!     --Check dimensions of plotting space                              
!                                                                       
         IF (pl_dim (1, 1) .le.cr_pos (1, i) .and.      &
             cr_pos (1, i) .le.pl_dim (1, 2) .and.      &
             pl_dim (2, 1) .le.cr_pos (2, i) .and.      &
             cr_pos (2, i) .le.pl_dim (2, 2) .and.      &
             pl_dim (3, 1) .le.cr_pos (3, i) .and.      &
             cr_pos (3, i) .le.pl_dim (3, 2) ) then  
!                                                                       
!     ------Determine distance to point in plot slice                   
!                                                                       
            DO j = 1, 3 
            v (j) = cr_pos (j, i) - pl_vec (j) 
            ENDDO 
            dist = abs (skalpro (pl_uvw, v, cr_gten) ) / d 
!                                                                       
!     ------Write atom if inside slice or if slice has been deselected  
!                                                                       
            IF (dist.le.pl_width.or.lno_slice) then 
!                                                                       
!     ------write atom position or projection onto abscissa/ordinate    
!                                                                       
               DO j = 1, 3 
               v (j) = cr_pos (j, i) 
               ENDDO 
!                                                                       
!     ------Write atom position in desired sequence of coordinates      
!                                                                       
               latom = .true. 
               WRITE (iff, 1001) cr_at_lis (cr_iscat (i) ),   &
                                 cr_iscat(i),(v(j)/cr_icc (j), j = 1, 3)
               WRITE (iff, 1002) cr_at_lis (cr_iscat (i) ),      &
                     cr_iscat (i),fac * cr_dw (cr_iscat (i) ),   &
                     fac * cr_dw (cr_iscat (i) ),                &
                     fac * cr_dw (cr_iscat (i) ), 0.0, 0.0, 0.0, &
                     povcolor ( pl_color (i) )
!     &                          pl_rgb(1,i),pl_rgb(2,i),pl_rgb(3,i);   
            ENDIF 
         ENDIF 
      ENDIF 
      ENDDO 
!                                                                       
      IF (.not.latom) then 
         ier_num = - 58 
         ier_typ = ER_APPL 
      ENDIF 
!                                                                       
      WRITE (iff, 541) 
                                                                        
  501 FORMAT    ('title ',a) 
  511 FORMAT    ('cell ',6(f10.6,1x)) 
  521 FORMAT    ('spgp ',a) 
  531 FORMAT    ('ellipsoids 50') 
  541 FORMAT    ('end') 
 1001 FORMAT     ('atom ',a4,1x,i4,3x,3(f11.6,1x)) 
!1002      format ('uij ',a4,1x,i4,3x,6(f8.6,1x),3f7.3)                 
 1002 FORMAT     ('uij ',a4,1x,i4,3x,6(f8.6,1x),a) 
      END SUBROUTINE plot_drawxtl                   
!*****7*****************************************************************
      SUBROUTINE plot_cif (iff) 
!-                                                                      
!     Writes the selected atoms as cif file (one unit cell)             
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE metric_mod
      USE modify_func_mod
      USE discus_plot_mod 
      USE tensors_mod
      USE errlist_mod 
      USE wink_mod
      IMPLICIT none 
!                                                                       
      REAL   , PARAMETER   :: eightpisq = 8.*3.14159265**2
!                                                                       
       
!                                                                       
      REAL d, dist 
      REAL v (3) 
      INTEGER, DIMENSION(3)  :: scalef
      INTEGER i, j, iff 
      LOGICAL lno_slice, latom 
!                                                                       
!                                                                       
!     scalef(1) = MAX(cr_icc(1), INT((cr_dim(1,2)-cr_dim(1,1)))+2)
!     scalef(2) = MAX(cr_icc(2), INT((cr_dim(2,2)-cr_dim(2,1)))+2)
!     scalef(3) = MAX(cr_icc(3), INT((cr_dim(3,2)-cr_dim(3,1)))+2)
      scalef(1) =                INT((cr_dim(1,2)-cr_dim(1,1)))+2
      scalef(2) =                INT((cr_dim(2,2)-cr_dim(2,1)))+2
      scalef(3) =                INT((cr_dim(3,2)-cr_dim(3,1)))+2
      WRITE (iff, 500) 
      WRITE (iff, 510) (cr_a0 (i) * scalef (i), i = 1, 3),  &
                       (cr_win (i), i = 1, 3)
!                                                                       
      latom = .false. 
!                                                                       
      IF (pl_hkl(1).eq.0.and.pl_hkl(2).eq.0.and.pl_hkl(3).eq.0) then
         lno_slice = .true. 
         d = 1.0 
      ELSE 
         lno_slice = .false. 
         d = sqrt (skalpro (pl_uvw, pl_uvw, cr_gten) ) 
         DO j = 1, 3 
         pl_mat (j, 1) = pl_abs (j) 
         pl_mat (j, 2) = pl_ord (j) 
         pl_mat (j, 3) = pl_uvw (j) 
         ENDDO 
         CALL invmat (pl_inv, pl_mat) 
         IF (ier_num.eq. - 1) then 
            RETURN 
         ENDIF 
      ENDIF 
!                                                                       
      DO i = 1, cr_natoms 
!                                                                       
!     --Select atom if:                                                 
!       type has been selected and                                      
!            all atoms are selected                            or       
!                                                                       
      IF (check_select_status (i, pl_latom (cr_iscat (i) ), cr_prop (i),   &
      pl_sel_prop) ) then                                               
!                                                                       
!     --Check dimensions of plotting space                              
!                                                                       
         IF (pl_dim (1, 1) .le.cr_pos (1, i) .and.        &
             cr_pos (1, i) .le.pl_dim (1, 2) .and.        &
             pl_dim (2, 1) .le.cr_pos (2, i) .and.        &
             cr_pos (2, i) .le.pl_dim (2, 2) .and.        &
             pl_dim (3, 1) .le.cr_pos (3, i) .and.        &
             cr_pos (3, i) .le.pl_dim (3, 2) ) then
!                                                                       
!     ------Determine distance to point in plot slice                   
!                                                                       
            DO j = 1, 3 
            v (j) = cr_pos (j, i) - pl_vec (j) 
            ENDDO 
            dist = abs (skalpro (pl_uvw, v, cr_gten) ) / d 
!                                                                       
!     ------Write atom if inside slice or if slice has been deselected  
!                                                                       
            IF (dist.le.pl_width.or.lno_slice) then 
!                                                                       
!     ------write atom position or projection onto abscissa/ordinate    
!                                                                       
               DO j = 1, 3 
               v (j) = cr_pos (j, i) + 0.01
               ENDDO 
!                                                                       
!     ------Write atom position in desired sequence of coordinates      
!                                                                       
               latom = .true. 
               WRITE (iff, 1000) cr_at_lis (cr_iscat (i) ),             &
                   ( (v (j) - cr_dim (j, 1) ) / scalef (j), j = 1, 3),  &
                    cr_dw ( cr_iscat (i) )/eightpisq
            ENDIF 
         ENDIF 
      ENDIF 
      ENDDO 
!                                                                       
      IF (.not.latom) then 
         ier_num = - 58 
         ier_typ = ER_APPL 
      ENDIF 
!                                                                       
  500 FORMAT ('data_00001',/,                                           &
     &        '_audit_creation_method   ''DISCUS''',/)                  
  510 FORMAT ('_cell_length_a     ',f15.4,/                             &
     &        '_cell_length_b     ',f15.4,/                             &
     &        '_cell_length_c     ',f15.4,/                             &
     &        '_cell_angle_alpha  ',f10.4,/                             &
     &        '_cell_angle_beta   ',f10.4,/                             &
     &        '_cell_angle_gamma  ',f10.4,//                            &
     &        '_symmetry_space_group_name_H-M   ''P1''',//              &
     &        'loop_',/                                                 &
     &        '_atom_site_label',/                                      &
     &        '_atom_site_fract_x',/                                    &
     &        '_atom_site_fract_y',/                                    &
     &        '_atom_site_fract_z',/                                    &
     &        '_atom_site_u_iso_or_equiv')                              
 1000 FORMAT (a4,3x,3(f11.6,1x),4x,f8.6) 
      END SUBROUTINE plot_cif                       
!*****7*****************************************************************
      SUBROUTINE plot_kuplot (iff, lkupl) 
!-                                                                      
!     Writes the selected atoms in a format suitable for plotting       
!     using KUPLOT or GNUPLOT                                           
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE metric_mod
      USE modify_func_mod
      USE discus_plot_mod 
      USE tensors_mod
      USE trafo_mod
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
      INTEGER idim 
      PARAMETER (idim = 3) 
!                                                                       
       
!                                                                       
      REAL d, dist, ps 
      REAL v (3), u (3) 
      INTEGER i, j, iff, pt, pc 
      LOGICAL lno_slice, latom, lkupl 
!                                                                       
!                                                                       
      latom = .false. 
      IF (pl_hkl(1).eq.0.and.pl_hkl(2).eq.0.and.pl_hkl(3).eq.0) then
         lno_slice = .true. 
         d = 1.0 
      ELSE 
         lno_slice = .false. 
         d = sqrt (skalpro (pl_uvw, pl_uvw, cr_gten) ) 
         DO j = 1, 3 
         pl_mat (j, 1) = pl_abs (j) 
         pl_mat (j, 2) = pl_ord (j) 
         pl_mat (j, 3) = pl_uvw (j) 
         ENDDO 
         CALL invmat (pl_inv, pl_mat) 
         IF (ier_num.eq. - 1) then 
            RETURN 
         ENDIF 
      ENDIF 
!                                                                       
      DO i = 1, cr_natoms 
!                                                                       
!     --Select atom if:                                                 
!       type has been selected and                                      
!            all atoms are selected                            or       
!                                                                       
      IF (check_select_status (i, pl_latom (cr_iscat (i) ), cr_prop (i),   &
      pl_sel_prop) ) then                                               
!                                                                       
!     --Check dimensions of plotting space                              
!                                                                       
         IF (pl_dim (1, 1) .le.cr_pos (1, i) .and.          &
             cr_pos (1, i) .le.pl_dim (1, 2) .and.          &
             pl_dim (2, 1) .le.cr_pos (2, i) .and.          &
             cr_pos (2, i) .le.pl_dim (2, 2) .and.          &
             pl_dim (3, 1) .le.cr_pos (3, i) .and.          &
             cr_pos (3, i) .le.pl_dim (3, 2) ) then
!                                                                       
!     ------Determine distance to point in plot slice                   
!                                                                       
            DO j = 1, 3 
            v (j) = cr_pos (j, i) - pl_vec (j) 
            ENDDO 
            dist = abs (skalpro (pl_uvw, v, cr_gten) ) / d 
!                                                                       
!     ------Write atom if inside slice or if slice has been deselected  
!                                                                       
            IF (dist.le.pl_width.or.lno_slice) then 
!                                                                       
!     ------write atom position or projection onto abscissa/ordinate    
!                                                                       
               IF (lno_slice) then 
                  DO j = 1, 3 
                  v (j) = cr_pos (j, i) 
                  ENDDO 
               ELSE 
                  DO j = 1, 3 
                  u (j) = cr_pos (j, i) 
                  ENDDO 
                  CALL trans (u, pl_inv, v, idim) 
               ENDIF 
!                                                                       
!     ------Write atom position in desired sequence of coordinates      
!                                                                       
               latom = .true. 
               pt = pl_typ (cr_iscat (i) ) 
               pc = pl_color (cr_iscat (i) ) 
               ps = pl_siz (cr_iscat (i) ) 
!                                                                       
               IF (pl_col.eq.'xyz') then 
                  CALL write_atom (lkupl, iff, i, v (1), v (2), v (3),  &
                  pt, pc, ps)                                           
               ELSEIF (pl_col.eq.'yzx') then 
                  CALL write_atom (lkupl, iff, i, v (2), v (3), v (1),  &
                  pt, pc, ps)                                           
               ELSEIF (pl_col.eq.'zxy') then 
                  CALL write_atom (lkupl, iff, i, v (3), v (1), v (2),  &
                  pt, pc, ps)                                           
               ELSEIF (pl_col.eq.'zyx') then 
                  CALL write_atom (lkupl, iff, i, v (3), v (2), v (1),  &
                  pt, pc, ps)                                           
               ELSEIF (pl_col.eq.'yxz') then 
                  CALL write_atom (lkupl, iff, i, v (2), v (1), v (3),  &
                  pt, pc, ps)                                           
               ELSEIF (pl_col.eq.'xzy') then 
                  CALL write_atom (lkupl, iff, i, v (1), v (3), v (2),  &
                  pt, pc, ps)                                           
               ENDIF 
            ENDIF 
         ENDIF 
      ENDIF 
      ENDDO 
!                                                                       
      IF (.not.latom) then 
         ier_num = - 58 
         ier_typ = ER_APPL 
      ENDIF 
!                                                                       
      END SUBROUTINE plot_kuplot                    
!*****7*****************************************************************
      SUBROUTINE plot_kuplot_mol (iff, lkupl) 
!-                                                                      
!     Writes the selected molecules in a format suitable                
!     for plotting using KUPLOT or GNUPLOT                              
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE metric_mod
      USE molecule_mod 
      USE discus_plot_mod 
      USE tensors_mod
      USE trafo_mod
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
      INTEGER idim 
      PARAMETER (idim = 3) 
!                                                                       
       
!                                                                       
      REAL d, dist, ps 
      REAL v (3), u (3) 
      INTEGER i, j, k, i0, iff, pt, pc 
      INTEGER i_start, i_end, imol 
      LOGICAL lno_slice, latom, lkupl 
!                                                                       
!     REAL skalpro 
!                                                                       
      latom = .false. 
      IF (pl_hkl(1).eq.0.and.pl_hkl(2).eq.0.and.pl_hkl(3).eq.0) then
         lno_slice = .true. 
         d = 1.0 
      ELSE 
         lno_slice = .false. 
         d = sqrt (skalpro (pl_uvw, pl_uvw, cr_gten) ) 
         DO j = 1, 3 
         pl_mat (j, 1) = pl_abs (j) 
         pl_mat (j, 2) = pl_ord (j) 
         pl_mat (j, 3) = pl_uvw (j) 
         ENDDO 
         CALL invmat (pl_inv, pl_mat) 
         IF (ier_num.eq. - 1) then 
            RETURN 
         ENDIF 
      ENDIF 
!                                                                       
      DO i = 1, mole_num_mole 
!                                                                       
!     --Select molecule if:                                             
!       type has been selected                                          
!                                                                       
      IF (pl_latom (mole_type (i) ) ) then 
!                                                                       
!     --Check dimensions of plotting space                              
!                                                                       
         i0 = mole_cont (mole_off (i) + 1) 
         IF (pl_dim (1, 1) .le.cr_pos (1, i0) .and.             &
             cr_pos (1, i0) .le.pl_dim (1, 2) .and.             &
             pl_dim (2, 1) .le.cr_pos (2, i0) .and.             &
             cr_pos (2, i0) .le.pl_dim (2, 2) .and.             &
             pl_dim (3, 1) .le.cr_pos (3, i0) .and.             &
             cr_pos (3, i0) .le.pl_dim (3, 2) ) then
!                                                                       
!     ------Determine distance to point in plot slice                   
!------ ------This is done only for the molecules origin                
!                                                                       
            DO j = 1, 3 
            v (j) = cr_pos (j, i0) - pl_vec (j) 
            ENDDO 
            dist = abs (skalpro (pl_uvw, v, cr_gten) ) / d 
!                                                                       
!     ------Write atom if inside slice or if slice has been deselected  
!                                                                       
            IF (dist.le.pl_width.or.lno_slice) then 
!                                                                       
!     ------Write atom position in desired sequence of coordinates      
!                                                                       
               latom = .true. 
               pt = pl_typ (mole_type (i) ) 
               pc = pl_color (mole_type (i) ) 
               ps = pl_siz (mole_type (i) ) 
!                                                                       
               IF (pl_mol_all) then 
                  i_start = 1 
                  i_end = mole_len (i) 
               ELSE 
                  i_start = 1 
                  i_end = 1 
               ENDIF 
!                                                                       
               DO imol = i_start, i_end 
               j = mole_cont (mole_off (i) + imol) 
!                                                                       
!     ------write atom position or projection onto abscissa/ordinate    
!                                                                       
               IF (lno_slice) then 
                  DO k = 1, 3 
                  v (k) = cr_pos (k, j) 
                  ENDDO 
               ELSE 
                  DO k = 1, 3 
                  u (k) = cr_pos (k, j) 
                  ENDDO 
                  CALL trans (u, pl_inv, v, idim) 
               ENDIF 
               IF (pl_col.eq.'xyz') then 
                  CALL write_atom (lkupl, iff, j, v (1), v (2), v (3),  &
                  pt, pc, ps)                                           
               ELSEIF (pl_col.eq.'yzx') then 
                  CALL write_atom (lkupl, iff, j, v (2), v (3), v (1),  &
                  pt, pc, ps)                                           
               ELSEIF (pl_col.eq.'zxy') then 
                  CALL write_atom (lkupl, iff, j, v (3), v (1), v (2),  &
                  pt, pc, ps)                                           
               ELSEIF (pl_col.eq.'zyx') then 
                  CALL write_atom (lkupl, iff, j, v (3), v (2), v (1),  &
                  pt, pc, ps)                                           
               ELSEIF (pl_col.eq.'yxz') then 
                  CALL write_atom (lkupl, iff, j, v (2), v (1), v (3),  &
                  pt, pc, ps)                                           
               ELSEIF (pl_col.eq.'xzy') then 
                  CALL write_atom (lkupl, iff, j, v (1), v (3), v (2),  &
                  pt, pc, ps)                                           
               ENDIF 
               ENDDO 
            ENDIF 
         ENDIF 
      ENDIF 
      ENDDO 
!                                                                       
      IF (.not.latom) then 
         ier_num = - 58 
         ier_typ = ER_APPL 
      ENDIF 
!                                                                       
      END SUBROUTINE plot_kuplot_mol                
!*****7*****************************************************************
      SUBROUTINE write_atom (lkupl, iff, iatom, x, y, z, pt, pc, ps) 
!                                                                       
      USE discus_config_mod 
      USE crystal_mod 
      USE celltoindex_mod
!     USE modify_mod
      USE discus_plot_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      REAL x, y, z, ps 
      INTEGER cr_end 
      INTEGER i, pt, pc, iff, iatom, isite, icell (3) 
      LOGICAL lkupl 
!                                                                       
      cr_end = cr_ncatoms*cr_icc(1) * cr_icc(2)*cr_icc(3) + 1
!                                                                       
      IF (pl_dens) then 
         IF (iatom.lt.cr_end) then 
            CALL indextocell (iatom, icell, isite) 
         ELSE 
            DO i = 1, 3 
            icell (i) = int (cr_pos (i, iatom) - cr_dim0 (i, 1) ) + 1
            ENDDO 
         ENDIF 
         x = x - float (icell (1) - 1) - cr_dim0 (1, 1) 
         y = y - float (icell (2) - 1) - cr_dim0 (2, 1) 
         z = z - float (icell (3) - 1) - cr_dim0 (3, 1) 
      ENDIF 
!                                                                       
      IF (lkupl) then 
         WRITE (iff, 1000) x, y, z, pt, pc, ps 
      ELSE 
         WRITE (iff, 2000) x, y, z 
      ENDIF 
!                                                                       
 1000 FORMAT  (3(2x,f12.6),2(2x,i2),2x,f6.2) 
 2000 FORMAT  (3(2x,f12.6)) 
      END SUBROUTINE write_atom                     
!*****7*****************************************************************
      SUBROUTINE plot_xbs (iff) 
!-                                                                      
!     Writes the selected atoms in a format suitable for plotting       
!     using xbs                                                         
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE metric_mod
      USE modify_func_mod
      USE discus_plot_mod 
      USE discus_plot_init_mod
      USE tensors_mod
      USE trans_sup_mod
      USE errlist_mod 
      USE param_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      CHARACTER(1024) zeile 
      INTEGER laenge 
      INTEGER i, j, k, iff 
      LOGICAL latom, lspace 
      LOGICAL lscreen 
      REAL uvw (4) 
      REAL absz (3) 
      REAL xmin (3), xmax (3), null (3) 
      REAL xx 
!                                                                       
!                                                                       
      DATA lscreen / .false. / 
      DATA lspace / .true. / 
      DATA null / 0.0, 0.0, 0.0 / 
!                                                                       
      CALL plot_ini_trans (1.0,                                &
                 pl_tran_g, pl_tran_gi, pl_tran_f, pl_tran_fi, &
                 cr_gten, cr_rten, cr_eps)

!                                                                       
      latom = .false. 
      uvw (4) = 1.0 
      DO j = 1, 3 
      xmin (j) = 0.0 
      xmax (j) = 0.0 
      ENDDO 
      xx = 0.0 
!                                                                       
      DO i = 1, cr_natoms 
!                                                                       
!     --Select atom if:                                                 
!       type has been selected and                                      
!            all atoms are selected                            or       
!                                                                       
      IF (check_select_status (i, pl_latom (cr_iscat (i) ), cr_prop (i),   &
      pl_sel_prop) ) then                                               
!                                                                       
!     --Check dimensions of plotting space                              
!                                                                       
         IF (pl_dim (1, 1) .le.cr_pos (1, i) .and.      &
             cr_pos (1, i) .le.pl_dim (1, 2) .and.      &
             pl_dim (2, 1) .le.cr_pos (2, i) .and.      &
             cr_pos (2, i) .le.pl_dim (2, 2) .and.      &
             pl_dim (3, 1) .le.cr_pos (3, i) .and.      &
             cr_pos (3, i) .le.pl_dim (3, 2) ) then

!                                                                       
!     ------write atom position                                         
!                                                                       
            latom = .true. 
            DO j = 1, 3 
            uvw (j) = cr_pos (j, i) 
            ENDDO 
            CALL tran_ca (uvw, pl_tran_f, lscreen) 
            DO j = 1, 3 
            xmin (j) = min (xmin (j), uvw (j) ) 
            xmax (j) = max (xmax (j), uvw (j) ) 
            ENDDO 
            WRITE (iff, 2100) cr_at_lis(cr_iscat(i)),(uvw(j), j = 1, 3)
         ENDIF 
      ENDIF 
      ENDDO 
!                                                                       
!     determine scale                                                   
!                                                                       
      IF (pl_scale (0) .eq. - 1.0) then 
         xx = max (abs (xmax (1) - xmin (1) ), abs (xmax (2) - xmin (1))&
                 , abs (xmax (3) - xmin (1) ) )                                
         pl_scale (1) = 200. / xx 
      ENDIF 
!                                                                       
!     write Atom specifications                                         
!                                                                       
      WRITE (iff, * ) 
      DO i = 0, cr_nscat 
      IF (pl_latom (i) ) then 
         WRITE (iff, 2200) cr_at_lis (i), pl_siz (i),                   &
                          (pl_rgb (j, i), j = 1, 3)
      ENDIF 
      ENDDO 
!                                                                       
!     Check that normal and asbzissa are not parallel                   
!                                                                       
      xx = do_bang (lspace, pl_uvw, null, pl_abs) 
      IF (xx.eq.0.0) then 
         ier_num = - 80 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
!                                                                       
      WRITE (zeile, 1000) pl_abs, pl_uvw 
      laenge = 82 
      CALL do_proj (zeile, laenge) 
      absz (1) = res_para (4) 
      absz (2) = res_para (5) 
      absz (3) = res_para (6) 
!                                                                       
      xx = do_blen (lspace, pl_uvw, null) 
      IF (xx.eq.0.0) then 
         ier_num = - 32 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
      DO j = 1, 3 
      uvw (j) = pl_uvw (j) / xx 
      ENDDO 
      CALL tran_ca (uvw, pl_tran_f, lscreen) 
      DO j = 1, 3 
      pl_tran_gi (3, j) = uvw (j) 
      ENDDO 
!                                                                       
      xx = do_blen (lspace, absz, null) 
      IF (xx.eq.0.0) then 
         ier_num = - 32 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
      DO j = 1, 3 
      uvw (j) = absz (j) / xx 
      ENDDO 
      CALL tran_ca (uvw, pl_tran_f, lscreen) 
      DO j = 1, 3 
      pl_tran_gi (1, j) = uvw (j) 
      ENDDO 
      pl_tran_gi (2, 1) = pl_tran_gi (3, 2) * pl_tran_gi (1, 3) -       &
      pl_tran_gi (3, 3) * pl_tran_gi (1, 2)                             
      pl_tran_gi (2, 2) = pl_tran_gi (3, 3) * pl_tran_gi (1, 1) -       &
      pl_tran_gi (3, 1) * pl_tran_gi (1, 3)                             
      pl_tran_gi (2, 3) = pl_tran_gi (3, 1) * pl_tran_gi (1, 2) -       &
      pl_tran_gi (3, 2) * pl_tran_gi (1, 1)                             
!                                                                       
      DO i = 1, 4 
      DO j = 1, 4 
      pl_tran_g (i, j) = pl_tran_gi (i, j) 
      ENDDO 
      ENDDO 
!                                                                       
      CALL invmat4 (pl_tran_g) 
!                                                                       
!     write Bond types                                                  
!                                                                       
      WRITE (iff, * ) 
      DO i = 1, cr_nscat 
      DO j = 1, cr_nscat 
      IF (pl_bond (i, j) ) then 
         WRITE (iff, 2300) cr_at_lis (i), cr_at_lis (j),                &
               pl_bond_len (1, i, j), pl_bond_len (2, i, j),            &
               pl_bond_rad (i, j), (pl_bond_col (k, i, j), k = 1, 3)
      ENDIF 
      ENDDO 
      ENDDO 
!                                                                       
!     write trailer, for right now a standard                           
!                                                                       
      WRITE (iff, * ) 
      WRITE (iff, 2400) ( (pl_tran_g (i, j), i = 1, 3), j = 1, 3) 
      WRITE (iff, 2500) 12.0 
      WRITE (iff, 2600) 1.0 
      WRITE (iff, 2700) pl_scale (1) 
      WRITE (iff, 2800) 1.0 
      WRITE (iff, 2900) 1.0 
      WRITE (iff, 3000) 0.00, 0.00 
      WRITE (iff, 3100) 1, 0, 1, 0, 0, 1, 0, 0, 0 
!                                                                       
      IF (.not.latom) then 
         ier_num = - 58 
         ier_typ = ER_APPL 
      ENDIF 
!                                                                       
 1000 FORMAT    (6(f12.6,','),'dddd') 
 2100 FORMAT    ('atom      ',a4,3f11.3) 
 2200 FORMAT    ('spec      ',a4,f10.3,3f7.2) 
 2300 FORMAT    ('bonds     ',a4,4x,a4,6f9.3) 
 2400 FORMAT    ('tmat',9f7.3) 
 2500 FORMAT    ('dist   ',f7.3) 
 2600 FORMAT    ('inc    ',f7.3) 
 2700 FORMAT    ('scale  ',f7.3) 
 2800 FORMAT    ('rfac   ',f7.3) 
 2900 FORMAT    ('bfac   ',f7.3) 
 3000 FORMAT    ('pos    ',2f7.3) 
 3100 FORMAT    ('switches',9i2) 
!                                                                       
      END SUBROUTINE plot_xbs                       
!*****7*****************************************************************
      SUBROUTINE plot_frames (iff) 
!-                                                                      
!     Writes the selected atoms in a format suitable for plotting       
!     using xbs                                                         
!     this subroutine writes the frames for a movie                     
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
!     USE metric_mod
      USE modify_func_mod
      USE discus_plot_mod 
      USE discus_plot_init_mod
      USE trans_sup_mod
      USE errlist_mod 
      USE param_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER i, j, iff 
      LOGICAL latom, lspace 
      LOGICAL lscreen 
      REAL uvw (4) 
      REAL xmin (3), xmax (3), null (3) 
      REAL xx 
!                                                                       
      DATA lscreen / .false. / 
      DATA lspace / .true. / 
      DATA null / 0.0, 0.0, 0.0 / 
!                                                                       
      CALL plot_ini_trans (1.0,                                &
                 pl_tran_g, pl_tran_gi, pl_tran_f, pl_tran_fi, &
                 cr_gten, cr_rten, cr_eps)
!                                                                       
      latom = .false. 
      uvw (4) = 1.0 
      DO j = 1, 3 
      xmin (j) = 0.0 
      xmax (j) = 0.0 
      ENDDO 
      xx = 0.0 
!                                                                       
      WRITE (iff, 1000) pl_title 
!                                                                       
      DO i = 1, cr_natoms 
!                                                                       
!     --Select atom if:                                                 
!       type has been selected and                                      
!            all atoms are selected                            or       
!                                                                       
      IF (check_select_status (i, pl_latom (cr_iscat (i) ), cr_prop (i),   &
      pl_sel_prop) ) then                                               
!                                                                       
!     --Check dimensions of plotting space                              
!                                                                       
         IF (pl_dim (1, 1) .le.cr_pos (1, i) .and.            &
             cr_pos (1, i) .le.pl_dim (1, 2) .and.            &
             pl_dim (2, 1) .le.cr_pos (2, i) .and.            &
             cr_pos (2, i) .le.pl_dim (2, 2) .and.            &
             pl_dim (3, 1) .le.cr_pos (3, i) .and.            &
             cr_pos (3, i) .le.pl_dim (3, 2) ) then 
!                                                                       
!     ------write atom position                                         
!                                                                       
            latom = .true. 
            DO j = 1, 3 
            uvw (j) = cr_pos (j, i) 
            ENDDO 
            CALL tran_ca (uvw, pl_tran_f, lscreen) 
            DO j = 1, 3 
            xmin (j) = min (xmin (j), uvw (j) ) 
            xmax (j) = max (xmax (j), uvw (j) ) 
            ENDDO 
            WRITE (iff, 2100) (uvw (j), j = 1, 3) 
         ENDIF 
      ENDIF 
      ENDDO 
      WRITE (iff, * ) 
!                                                                       
      IF (.not.latom) then 
         ier_num = - 58 
         ier_typ = ER_APPL 
      ENDIF 
!                                                                       
 1000 FORMAT    ('frame ',a) 
 2100 FORMAT    (3f11.3) 
!                                                                       
      END SUBROUTINE plot_frames                    
!*****7*****************************************************************
      SUBROUTINE plot_diamond (iff) 
!-                                                                      
!     Writes the selected atoms in a format suitable for plotting       
!     using diamond                                                     
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE modify_func_mod
      USE discus_plot_mod 
      USE discus_plot_init_mod
      USE trans_sup_mod
      USE blanks_mod
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      CHARACTER(1024) zeile 
      INTEGER laenge 
      INTEGER i, j, iff 
      INTEGER :: natoms
      LOGICAL latom, lspace 
      LOGICAL lscreen 
      REAL uvw (4) 
!                                                                       
!                                                                       
      DATA lscreen / .false. / 
      DATA lspace / .true. / 
!                                                                       
      CALL plot_ini_trans (1.0,                                &
                 pl_tran_g, pl_tran_gi, pl_tran_f, pl_tran_fi, &
                 cr_gten, cr_rten, cr_eps)
!                                                                       
      latom = .false. 
!
      natoms = 0
      DO i = 1, cr_natoms 
         IF (check_select_status(i, pl_latom(cr_iscat(i)), cr_prop (i),pl_sel_prop) ) THEN
            natoms = natoms + 1
         ENDIF
      ENDDO
!                                                                       
!     Write number of atoms                                             
!                                                                       
      WRITE (zeile, 1000) natoms 
      laenge = 20 
      CALL rem_bl (zeile, laenge) 
      WRITE (iff, 1100) zeile (1:laenge) 
!                                                                       
!     Write title                                                       
!                                                                       
      WRITE (iff, 1100) cr_name 
!                                                                       
      DO i = 1, cr_natoms 
!                                                                       
!     --Select atom if:                                                 
!       type has been selected and                                      
!            all atoms are selected                            or       
!                                                                       
      IF (check_select_status (i, pl_latom (cr_iscat (i) ), cr_prop (i),   &
      pl_sel_prop) ) then                                               
!                                                                       
!     --Check dimensions of plotting space                              
!                                                                       
         IF (pl_dim (1, 1) .le.cr_pos (1, i) .and.       &
             cr_pos (1, i) .le.pl_dim (1, 2) .and.       &
             pl_dim (2, 1) .le.cr_pos (2, i) .and.       &
             cr_pos (2, i) .le.pl_dim (2, 2) .and.       &
             pl_dim (3, 1) .le.cr_pos (3, i) .and.       &
             cr_pos (3, i) .le.pl_dim (3, 2) ) then
!                                                                       
!     ------write atom position                                         
!                                                                       
            latom = .true. 
            DO j = 1, 3 
            uvw (j) = cr_pos (j, i) 
            ENDDO 
            CALL tran_ca (uvw, pl_tran_f, lscreen) 
            WRITE (iff, 2100) cr_at_lis (cr_iscat (i) ),  &
                              (uvw (j), j = 1, 3)
         ENDIF 
      ENDIF 
      ENDDO 
!                                                                       
      IF (.not.latom) then 
         ier_num = - 58 
         ier_typ = ER_APPL 
      ENDIF 
!                                                                       
 1000 FORMAT    (i20) 
 1100 FORMAT    (a) 
 2100 FORMAT    (a4,3f11.3) 
!                                                                       
      END SUBROUTINE plot_diamond                   
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
USE prompt_mod
USE trig_degree_mod
use matrix_mod
!
IMPLICIT NONE
!
REAL, PARAMETER :: azero = 1.0
REAL, PARAMETER :: eps   = 1.0E-5
CHARACTER(LEN=*), INTENT(IN) :: kill
!
INTEGER, PARAMETER  :: ITMP     = 78  ! temporary unit number
LOGICAL, PARAMETER  :: lscreen = .FALSE.
CHARACTER(LEN=1024) :: tempfile
CHARACTER(LEN=1024) :: line
CHARACTER(LEN=1024) :: line_move
INTEGER             :: i,j,k
INTEGER             :: length
INTEGER             :: jmol_pid
INTEGER             :: ios
LOGICAL             :: did_kill = .false.
LOGICAL             :: lpresent = .false.
LOGICAL             :: lsecond  = .false.
LOGICAL             :: lunit    = .false.
REAL, DIMENSION(4)  :: hkl      !Normal in reciprocal space for augmented matrix
REAL, DIMENSION(3)  :: n_hkl    !Normal in reciprocal space
REAL, DIMENSION(3)  :: n_uvw    !Normal in direct space
REAL, DIMENSION(3)  :: v        !temporary vector
REAL, DIMENSION(3)  :: u        !temporary vector
REAL, DIMENSION(3)  :: axis     !axis for moveto command
REAL, DIMENSION(3)  :: p_a      !projected abscissa
REAL, DIMENSION(3)  :: p_o      !projected abscissa
REAL, DIMENSION(3,3):: roti     !Rotation matrix
REAL, DIMENSION(3,3):: rotf     !Rotation matrix
REAL, DIMENSION(3,3):: test     !Test for unit matrix
REAL :: det, trace
INTEGER :: ier
REAL                :: rr
REAL                :: beta
!
IF(pl_prog=='jmol') THEN
   IF(kill=='yes') THEN
      INQUIRE(FILE='/tmp/jmol.pid', EXIST=lpresent)
      IF(lpresent) THEN
         line = 'rm -f /tmp/jmol.pid'
         CALL system(line)
      ENDIF
      did_kill = .FALSE.
      IF(operating=='Linux') THEN
         WRITE(line,'(a,i8,a)') 'ps j | grep ',PID,' | grep java |  grep -v grep | awk ''{print $2}'' >> /tmp/jmol.pid'
         CALL system(line)
!
      ELSEIF(operating(1:6)=='cygwin' .OR. operating(1:7)=='Windows') THEN
         line = 'ps aux | grep ''jmol -s /tmp/jmol.mol'' | grep -v grep | awk ''{print $1}'' > /tmp/jmol.pid'
         CALL system(line)
         line = 'ps aux | grep java |                      grep -v grep | awk ''{print $1}'' >> /tmp/jmol.pid'
         CALL system(line)
      ENDIF
      tempfile = '/tmp/jmol.pid'
      CALL oeffne( ITMP, tempfile, 'old')
      IF(ier_num==0) THEN
         READ(ITMP,*,IOSTAT=ios) jmol_pid
         DO WHILE (.NOT.IS_IOSTAT_END(ios)) 
            WRITE(line,'(a,i12,a)') 'kill -9 ',jmol_pid, ' > /dev/null'
            CALL system(line)
            did_kill = .TRUE.
            READ(ITMP,*,IOSTAT=ios) jmol_pid
         ENDDO
         CLOSE(ITMP)
         line = 'rm -f /tmp/jmol.pid'
         CALL system(line)
         IF(operating=='Linux') THEN
            IF(did_kill) WRITE(output_io,*) ' Processes were killed, hit ENTER to retrieve prompt'
         ENDIF
      ENDIF
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
   IF(operating=='Linux') THEN
      tempfile = '/tmp/jmol.mol'
   ELSEIF(operating(1:6)=='cygwin' .OR. operating(1:7)=='Windows') THEN
      tempfile = '/opt/jmol/jmol.mol'
   ENDIF
   CALL oeffne( ITMP, tempfile, 'unknown')
   IF(operating=='Linux') THEN
!
!     Write load command, enclose file name in apostrophes to cover 
!     blanks in filenemes
      WRITE(ITMP,'(3a)') 'load ''',                  &
         current_dir(1:LEN_TRIM(current_dir))//      &
         pl_out(1:LEN_TRIM(pl_out)), ''' {1 1 1}'
   ELSEIF(operating(1:6)=='cygwin') THEN
      WRITE(ITMP,'(3a)') 'load ''',                   &
         operat_top(1:LEN_TRIM(operat_top))//        &
         current_dir(1:LEN_TRIM(current_dir))//      &
         pl_out(1:LEN_TRIM(pl_out)), ' {1 1 1}'
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
   ELSEIF(operating(1:6)=='cygwin') THEN
      WRITE(line,'(a,a,a,a,a)') pl_jmol(1:LEN_TRIM(pl_jmol)), ' -s jmol.mol > /dev/null &'
   ELSEIF(operating(1:7)=='Windows') THEN
      WRITE(line,'(a,a,a,a,a)') pl_jmol(1:LEN_TRIM(pl_jmol)), ' -s jmol.mol > /dev/null &'
!
   ENDIF
   WRITE(output_io,'(a)') ' JMOL may take a moment to show up'
   CALL system(line)
ENDIF
END SUBROUTINE plot_inter
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
END MODULE discus_plot_menu
