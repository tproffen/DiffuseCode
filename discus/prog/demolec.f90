MODULE demolec
!
USE demolec_mod
PRIVATE
PUBLIC :: demolecularize
PUBLIC :: demol_reset
!
CONTAINS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE demolecularize
!
! Menu to remove a molecule status
!
USE calc_expr_mod
USE class_macro_internal
USE crystal_mod
USE discus_allocate_appl_mod
USE molecule_mod
USE prop_para_func
!
USE do_eval_mod
USE do_wait_mod
USE doact_mod
USE errlist_mod
USE learn_mod
USE lib_do_operating_mod
USE lib_echo
USE lib_errlist_func
USE lib_help
USE lib_length
USE lib_macro_func
USE precision_mod
USE prompt_mod
USE sup_mod
USE str_comp_mod
!
IMPLICIT NONE
!
CHARACTER (LEN=8)                       :: befehl! command on input line
CHARACTER(LEN=LEN(prompt))              :: orig_prompt  ! original prompt
CHARACTER (LEN=PREC_STRING)                    :: line  ! input line
CHARACTER (LEN=PREC_STRING)                    :: zeile ! remainder with parameters
INTEGER                                 :: indxg ! location of "="
INTEGER                                 :: lp    ! length of zeile
INTEGER                                 :: laenge
INTEGER                                 :: lbef
INTEGER                                 :: nmoletype ! temp variable for allocation
INTEGER                                 :: natomtype ! temp variable for allocation
!
LOGICAL                                 :: sel   ! condition of select or deselect
LOGICAL                                 :: lend  ! condition of EOF
!
!
!
orig_prompt = prompt
prompt = prompt (1:len_str (prompt) ) //'/demol'
!
nmoletype = MAX(DEM_MAX_MOLETYPE, mole_num_type)
natomtype = MAX(DEM_MAX_ATOMTYPE, MAXSCAT, cr_nscat)
IF(nmoletype>DEM_MAX_MOLETYPE .OR.        &
   natomtype>DEM_MAX_ATOMTYPE       ) THEN
   CALL alloc_demol(natomtype, nmoletype)
ENDIF
!
main_loop: DO
   CALL no_error
   CALL get_cmd (line, laenge, befehl, lbef, zeile, lp, prompt)
   no_err: IF(ier_num.eq.0) THEN
      no_com: IF(line /= ' '      .AND. line(1:1) /= '#' .AND.      &
                 line /= char(13) .AND. line(1:1) /= '!'        ) THEN
!                                                                       
!     ----search for "="                                                
!                                                                       
         indxg = index (line, '=') 
         is_math: IF(indxg.ne.0                                             &
                     .AND..NOT. (str_comp (befehl, 'echo',   2, lbef, 4) )    &
                     .AND..NOT. (str_comp (befehl, 'system', 2, lbef, 6) )    &
                     .AND..NOT. (str_comp (befehl, 'help',   2, lbef, 4) .OR. &
                                 str_comp (befehl, '?   ',   2, lbef, 4) )    &
                     .AND. INDEX(line,'==') == 0                            ) THEN
!                                                                       
! ------evaluate an expression and assign the value to a variabble      
!                                                                       
               CALL do_math (line, indxg, laenge)
         ELSE is_math                            ! is_math, al other commands
!                                                                       
!------ ----execute a macro file                                        
!                                                                       
            is_generic: IF (befehl (1:1) .eq.'@') THEN     ! macro, reset or all other commands
               IF (laenge.ge.2) THEN 
                  line(1:laenge-1) = line(2:laenge)
                  line(laenge:laenge) = ' '
                  laenge = laenge - 1
                  CALL file_kdo(line, laenge)
               ELSE 
                  ier_num = - 13 
                  ier_typ = ER_MAC 
               ENDIF 
!                                                                       
!     ----continues a macro 'continue'                                  
!                                                                       
            ELSEIF (str_comp (befehl, 'continue', 3, lbef, 8) ) THEN  is_generic
               CALL macro_continue (zeile, lp) 
!                                                                       
!     ----Echo a string, just for interactive check in a macro 'echo'   
!                                                                       
            ELSEIF (str_comp (befehl, 'echo', 2, lbef, 4) ) THEN  is_generic
               CALL echo (zeile, lp) 
!                                                                       
!      ---Evaluate an expression, just for interactive check 'eval'     
!                                                                       
            ELSEIF (str_comp (befehl, 'evaluate', 2, lbef, 8) ) THEN 
               CALL do_eval (zeile, lp, .TRUE.) 
!                                                                       
!     ----exit 'exit'                                                   
!                                                                       
            ELSEIF (str_comp (befehl, 'exit', 2, lbef, 4) ) THEN  is_generic
               lend = .true. 
               EXIT main_loop
!                                                                       
!     ----help 'help' , '?'                                             
!                                                                       
            ELSEIF (str_comp (befehl, 'help', 1, lbef, 4) .OR.  &
                    str_comp (befehl, '?   ', 1, lbef, 4) ) THEN is_generic
               IF (str_comp (zeile, 'errors', 2, lp, 6) ) THEN 
                  lp = lp + 7 
                  CALL do_hel ('discus '//zeile, lp) 
               ELSE 
                  lp = lp + 12 
                  CALL do_hel ('discus demo '//zeile, lp) 
               ENDIF 
!                                                                       
!------- -Operating System Kommandos 'syst'                             
!                                                                       
            ELSEIF (str_comp (befehl, 'system', 2, lbef, 6) ) THEN is_generic
                IF (zeile.ne.' '.and.zeile.ne.char (13) ) THEN
                   CALL do_operating (zeile (1:lp), lp) 
                ELSE 
                   ier_num = - 6 
                   ier_typ = ER_COMM 
                ENDIF 
!                                                                       
!------  -----waiting for user input                                    
!                                                                       
            ELSEIF (str_comp (befehl, 'wait', 3, lbef, 4) ) THEN  is_generic
                CALL do_input (zeile, lp) 
!
            ELSEIF (str_comp (befehl, 'reset', 3, lbef, 5)) THEN is_generic
               CALL demol_reset
            ELSE is_generic    ! macro, reset or all other commands is_generic
!
               is_com: IF(str_comp (befehl, 'include', 3, lbef, 7) ) THEN 
                  CALL demol_incl(zeile, lp)
!                                                                       
!------ --Handle property settings 'property'                           
!                                                                       
               ELSEIF (str_comp (befehl, 'property', 4, lbef, 8) ) THEN is_com
!                                                                       
                  CALL property_select (zeile, lp, dem_sel_prop)

               ELSEIF(str_comp (befehl, 'run', 3, lbef, 3) ) THEN is_com
                  CALL demol_run
               ELSEIF(str_comp (befehl, 'select',   3, lbef, 6) .OR.  &
                      str_comp (befehl, 'deselect', 3, lbef, 8)     ) THEN is_com
!
!--Select
!
                  sel = str_comp (befehl, 'select', 3, lbef, 6)
                  CALL demol_msel(zeile, lp, sel)
               ELSEIF(str_comp (befehl, 'show', 3, lbef, 4) ) THEN is_com
!
!--show
!
                  CALL demol_show
!
               ELSE is_com
                  ier_num = -8
                  ier_typ = ER_COMM
!
               ENDIF is_com ! END IF BLOCK actual commands
            ENDIF is_generic ! END IF BLOCK generic commands
         ENDIF is_math      ! END IF BLOCK math equation or specific command
      ENDIF no_com          ! END IF BLOCK no comment
   ENDIF no_err             ! END IF BLOCK no error reading input
!
   IF (ier_num.ne.0) THEN 
      CALL errlist 
      IF (ier_sta.ne.ER_S_LIVE) THEN 
         IF (lmakro .OR. lmakro_error) THEN  ! Error within macro or termination errror
            IF(sprompt /= prompt ) THEN
               ier_num = -10
               ier_typ = ER_COMM
               ier_msg(1) = ' Error occured in demolec menu'
               prompt_status = PROMPT_ON 
               prompt = orig_prompt
               RETURN
            ELSE
               IF(lmacro_close) THEN
                  CALL macro_close (-1)
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
!
!
ENDDO main_loop     ! END DO main loop of menu 
!
prompt = orig_prompt
!
END SUBROUTINE demolecularize
!
!*******************************************************************************
!
SUBROUTINE demol_incl(zeile, lp)
!
! Select molecule ranges, or  atom ranges
!
USE crystal_mod
USE molecule_mod
!
USE errlist_mod
USE ber_params_mod
USE get_params_mod
USE precision_mod
USE take_param_mod
USE str_comp_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: zeile
INTEGER, INTENT(INOUT)          :: lp
!
INTEGER, PARAMETER                   :: MAXW = 4
INTEGER, PARAMETER                   :: MIN_PARA = 2
CHARACTER(LEN=MAX(PREC_STRING,LEN(zeile))), DIMENSION(MAXW) :: cpara
INTEGER            , DIMENSION(MAXW) :: lpara
CHARACTER(LEN=PREC_STRING), DIMENSION(MAX(MIN_PARA,MAXSCAT+1, MOLE_MAX_TYPE)) :: ccpara
INTEGER            , DIMENSION(MAX(MIN_PARA,MAXSCAT+1, MOLE_MAX_TYPE)) :: llpara
REAL(KIND=PREC_DP) , DIMENSION(MAX(MIN_PARA,MAXSCAT+1, MOLE_MAX_TYPE)) :: wwerte
INTEGER                              :: ianz
INTEGER                              :: iianz
!
CHARACTER(LEN=MAX(PREC_STRING,LEN(zeile)))                  :: string
INTEGER                              :: j, k
INTEGER                              :: llp
INTEGER                              :: MAXWW
!
INTEGER, PARAMETER :: NOPTIONAL = 2
INTEGER, PARAMETER :: O_MOLRANG = 1
INTEGER, PARAMETER :: O_ATOMRANG= 2
CHARACTER(LEN=   9), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=MAX(PREC_STRING,LEN(zeile))), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent!opt. para is present
REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 0 ! Number of values to calculate 
!
!
DATA oname  / 'molerange', 'atomrange'  /
DATA loname /  9         ,  9           /
opara  =  (/ '0.0000', '0.0000' /)   ! Always provide fresh default values
lopara =  (/  6,        6       /)
owerte =  (/  0.0,      0.0     /)
!
MAXWW = UBOUND(ccpara,1)
!
CALL get_params(zeile, ianz, cpara, lpara, MAXW, lp)
!
CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
IF(lpresent(O_MOLRANG)) THEN                ! User provided molecule types
   IF(str_comp (opara(O_MOLRANG), 'all', 3, lopara(O_MOLRANG), 3) ) THEN
      dem_molerange(1) = 1
      dem_molerange(2) = -1
   ELSEIF(str_comp (opara(O_MOLRANG), 'none', 3, lopara(O_MOLRANG), 4) ) THEN
      dem_molerange(1) = 0
      dem_molerange(2) = 0 
   ELSE
      IF(opara(O_MOLRANG)(1:1)=='(') opara(O_MOLRANG)(1:1)= ' '
      IF(opara(O_MOLRANG)(lopara(O_MOLRANG):lopara(O_MOLRANG))==')') &
         opara(O_MOLRANG)(lopara(O_MOLRANG):lopara(O_MOLRANG))= ' '
      string = opara(O_MOLRANG)(1:lopara(O_MOLRANG))
      llp = lopara(O_MOLRANG)
      CALL get_params(string, iianz, ccpara, llpara, MAXW, llp)
      IF(ier_num/=0) RETURN
      IF(iianz == 2 .AND. ccpara(2)=='last') ccpara(2) = '  -1'
      wwerte(:) = 0
      CALL ber_params(iianz, ccpara, llpara, wwerte, MAXWW)
      IF(ier_num/=0) RETURN
      IF(iianz == 1) THEN
         j = NINT(wwerte(1))
         IF(j<1 .OR. j>mole_num_mole) THEN
            ier_num = -63
            ier_typ = ER_APPL
            RETURN
         ENDIF
         dem_molerange(1:2) = j
      ELSEIF(iianz == 2) THEN
         j = NINT(wwerte(1))
         k = NINT(wwerte(2))
         IF(j<1 .OR. j>mole_num_mole .OR.            &
            (k/=-1 .AND.(k<1 .OR. k>mole_num_mole .OR. k<j))) THEN
            ier_num = -63
            ier_typ = ER_APPL
            RETURN
         ENDIF
         dem_molerange(1) = j
         dem_molerange(2) = k
      ENDIF
   ENDIF
ELSE
   dem_molerange(1) =  1
   dem_molerange(2) = -1
ENDIF
!
IF(lpresent(O_ATOMRANG)) THEN                ! User provided molecule types
   IF(str_comp (opara(O_ATOMRANG), 'all', 3, lopara(O_ATOMRANG), 3) ) THEN
      dem_atomrange(1) =  1
      dem_atomrange(2) = -1
   ELSEIF(str_comp (opara(O_ATOMRANG), 'none', 3, lopara(O_ATOMRANG), 4) ) THEN
      dem_atomrange(1) =  0
      dem_atomrange(2) =  0
   ELSE
      IF(opara(O_ATOMRANG)(1:1)=='(') opara(O_ATOMRANG)(1:1)= ' '
      IF(opara(O_ATOMRANG)(lopara(O_ATOMRANG):lopara(O_ATOMRANG))==')') &
         opara(O_ATOMRANG)(lopara(O_ATOMRANG):lopara(O_ATOMRANG))= ' '
      string = opara(O_ATOMRANG)(1:lopara(O_ATOMRANG))
      llp = lopara(O_ATOMRANG)
      CALL get_params(string, iianz, ccpara, llpara, MAXW, llp)
      IF(ier_num/=0) RETURN
      IF(iianz == 2 .AND. ccpara(2)=='last') ccpara(2) = '  -1'
      wwerte(:) = 0
      CALL ber_params(iianz, ccpara, llpara, wwerte, MAXWW)
      IF(ier_num/=0) RETURN
      IF(iianz == 1) THEN
         j = NINT(wwerte(1))
         IF(j<1 .OR. j>mole_num_mole) THEN
            ier_num = -19
            ier_typ = ER_APPL
            RETURN
         ENDIF
         dem_atomrange(1:2) = j
      ELSEIF(iianz == 2) THEN
         j = NINT(wwerte(1))
         k = NINT(wwerte(2))
         IF(j<1 .OR. j>cr_natoms .OR.            &
            (k/=-1 .AND. (k<1 .OR. k>cr_natoms .OR. k<j))) THEN
            ier_num = -19
            ier_typ = ER_APPL
            RETURN
         ENDIF
         dem_atomrange(1) = j
         dem_atomrange(2) = k
      ENDIF
   ENDIF
ELSE
   dem_atomrange(1) =  0
   dem_atomrange(2) =  0
ENDIF
!
END SUBROUTINE demol_incl
!
!*******************************************************************************
!
SUBROUTINE demol_msel(zeile, lp, sel)
!
! Select molecule types, or  atom types
!
USE crystal_mod
USE molecule_mod
!
USE errlist_mod
USE ber_params_mod
USE get_params_mod
USE precision_mod
USE take_param_mod
USE str_comp_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: zeile
INTEGER, INTENT(INOUT)          :: lp
LOGICAL, INTENT(IN)             :: sel    ! =true for select
!
INTEGER, PARAMETER                   :: MAXW = 4
INTEGER, PARAMETER                   :: MIN_PARA = 2
CHARACTER(LEN=MAX(PREC_STRING,LEN(zeile))), DIMENSION(MAXW) :: cpara
INTEGER            , DIMENSION(MAXW) :: lpara
CHARACTER(LEN=MAX(PREC_STRING,LEN(zeile))), DIMENSION(MAX(MIN_PARA,MAXSCAT+1, MOLE_MAX_TYPE)) :: ccpara
INTEGER            , DIMENSION(MAX(MIN_PARA,MAXSCAT+1, MOLE_MAX_TYPE)) :: llpara
REAL(KIND=PREC_DP) , DIMENSION(MAX(MIN_PARA,MAXSCAT+1, MOLE_MAX_TYPE)) :: wwerte
INTEGER                              :: ianz
INTEGER                              :: iianz
!
CHARACTER(LEN=PREC_STRING)                  :: string
INTEGER                              :: i, j
INTEGER                              :: llp
INTEGER                              :: MAXWW
!
INTEGER, PARAMETER :: NOPTIONAL = 2
INTEGER, PARAMETER :: O_MOLTYPE = 1
INTEGER, PARAMETER :: O_ATOMTYPE= 2
CHARACTER(LEN=   8), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=MAX(PREC_STRING,LEN(zeile))), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent!opt. para is present
REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 0 ! Number of values to calculate 
!
!
DATA oname  / 'moletype', 'atomtype'  /
DATA loname /  8        ,  8          /
opara  =  (/ '0.0000', '0.0000' /)   ! Always provide fresh default values
lopara =  (/  6,        6       /)
owerte =  (/  0.0,      0.0     /)
!
MAXWW = UBOUND(ccpara,1)
!
CALL get_params(zeile, ianz, cpara, lpara, MAXW, lp)
!
CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
IF(lpresent(O_MOLTYPE)) THEN                ! User provided molecule types
   IF(str_comp (opara(O_MOLTYPE), 'all', 3, lopara(O_MOLTYPE), 3) ) THEN
      dem_lmoletype(:) = sel
   ELSEIF(str_comp (opara(O_MOLTYPE), 'none', 3, lopara(O_MOLTYPE), 4) ) THEN
      dem_lmoletype(:) = .NOT.sel
   ELSE
      IF(opara(O_MOLTYPE)(1:1)=='(') opara(O_MOLTYPE)(1:1)= ' '
      IF(opara(O_MOLTYPE)(lopara(O_MOLTYPE):lopara(O_MOLTYPE))==')') &
         opara(O_MOLTYPE)(lopara(O_MOLTYPE):lopara(O_MOLTYPE))= ' '
      string = opara(O_MOLTYPE)(1:lopara(O_MOLTYPE))
      llp = lopara(O_MOLTYPE)
      CALL get_params(string, iianz, ccpara, llpara, MAXW, llp)
      IF(ier_num/=0) RETURN
      CALL ber_params(iianz, ccpara, llpara, wwerte, MAXWW)
      IF(ier_num/=0) RETURN
      DO i=1, iianz
         j = NINT(wwerte(i))
         IF(j<1 .OR. j>mole_num_type) THEN
            ier_num = -64
            ier_typ = ER_APPL
            RETURN
         ENDIF
         dem_lmoletype(j) = sel
      ENDDO
   ENDIF
ELSE
   dem_lmoletype(:) = .NOT.sel
ENDIF
!
IF(lpresent(O_ATOMTYPE)) THEN                ! User provided molecule types
   IF(str_comp (opara(O_ATOMTYPE), 'all', 3, lopara(O_ATOMTYPE), 3) ) THEN
      dem_latomtype(:) = sel
   ELSEIF(str_comp (opara(O_ATOMTYPE), 'none', 3, lopara(O_ATOMTYPE), 4) ) THEN
      dem_latomtype(:) = .NOT.sel
   ELSE
      IF(opara(O_ATOMTYPE)(1:1)=='(') opara(O_ATOMTYPE)(1:1)= ' '
      IF(opara(O_ATOMTYPE)(lopara(O_ATOMTYPE):lopara(O_ATOMTYPE))==')') &
         opara(O_ATOMTYPE)(lopara(O_ATOMTYPE):lopara(O_ATOMTYPE))= ' '
      string = opara(O_ATOMTYPE)(1:lopara(O_ATOMTYPE))
      llp = lopara(O_ATOMTYPE)
      CALL get_params(string, iianz, ccpara, llpara, MAXW, llp)
      IF(ier_num/=0) RETURN
      CALL ber_params(iianz, ccpara, llpara, wwerte, MAXWW)
      IF(ier_num/=0) RETURN
      DO i=1, iianz
        j = NINT(wwerte(i))
         IF(j<0 .OR. j>cr_nscat) THEN
            ier_num = -122
            ier_typ = ER_APPL
            RETURN
         ENDIF
         dem_latomtype(j) = sel
      ENDDO
   ENDIF
ELSE
   dem_latomtype(:) = .NOT.sel
ENDIF
!
END SUBROUTINE demol_msel
!
!*******************************************************************************
!
SUBROUTINE demol_show
!
USE crystal_mod
USE molecule_mod
USE prop_char_mod
USE prompt_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=32)                      :: c_property
INTEGER, DIMENSION(1:DEM_MAX_MOLETYPE) :: moletype_list
INTEGER, DIMENSION(1:DEM_MAX_ATOMTYPE) :: atomtype_list
INTEGER :: j,k
INTEGER :: i
INTEGER :: length
!
j = 0
DO i=1, mole_num_type
  IF(dem_lmoletype(i)) THEN
     j = j+1
     moletype_list(j) = i
  ENDIF
ENDDO
WRITE(output_io, 3300) (moletype_list(k),k=1,j)
!
IF(dem_molerange(2)==-1) THEN
   IF(dem_molerange(1)==1) THEN
      WRITE(output_io, 3400)
   ELSE
      WRITE(output_io, 3500) dem_molerange(1)
   ENDIF
ELSEIF(dem_molerange(1)==0  .AND.dem_molerange(2)==0) THEN
  WRITE(output_io, 3450) 
ELSE
  WRITE(output_io, 3550) dem_molerange(:)
ENDIF
!
j = 0
DO i=0, cr_nscat
  IF(dem_latomtype(i)) THEN
     j = j+1
     atomtype_list(j) = i
  ENDIF
ENDDO
WRITE(output_io, 4300) (atomtype_list(k),k=1,j)
!
IF(dem_atomrange(1)==0 .AND.dem_atomrange(2)==0) THEN
   WRITE(output_io, 4401) 
ELSEIF(dem_atomrange(2)==-1) THEN
   IF(dem_atomrange(1)==1) THEN
      WRITE(output_io, 4400)
   ELSE
      WRITE(output_io, 4500) dem_atomrange(1)
   ENDIF
ELSE
  WRITE(output_io, 4550) dem_atomrange(:)
ENDIF
CALL char_prop_2 (c_property, dem_sel_prop (1), dem_sel_prop (0),   &
                  length)
WRITE (output_io, 5000) c_property (1:length)
!
!
3300 FORMAT( '   Selected molecule types   : ',2x,50(i4,1x))
3400 FORMAT( '   Range of incl. molecules  : ',                    &
     &          '  All molecules included')
3450 FORMAT( '   Range of incl. molecules  : ',                    &
     &          '  No  molecules included')
3500 FORMAT( '   Range of incl. molecules  : ',i9,' to all remaining')
3550 FORMAT( '   Range of incl. molecules  : ',i9,' to ',i9)
4300 FORMAT( '   Must include atom types   : ',2x,50(i4,1x))
4401 FORMAT( '   Range of incl. atoms      : ',                    &
     &          '  No range needs to be included')
4400 FORMAT( '   Atoms must be in range    : ',                    &
     &          '  All atoms must be included')
4500 FORMAT( '   Atoms must be in range    : ',i9,' to all remaining')
4550 FORMAT( '   Atoms must be in range    : ',i9,' to ',i9)
!
5000 FORMAT(/' Atom properties         : ','NMDOEI'/               &
                  '      absent=- ignored=. : ',a)
!
END SUBROUTINE demol_show
!
!*******************************************************************************
!
!
SUBROUTINE demol_run
!
! Performs the demolecularization
!
USE crystal_mod
USE modify_func_mod
USE molecule_mod
USE prop_para_mod
!
USE errlist_mod
!
IMPLICIT NONE
!
INTEGER :: i, j, k, l, mmol
INTEGER :: istart, iend
INTEGER :: nat                  ! >0 if there is a true atom in a molecule
LOGICAL :: isdel = .FALSE.
LOGICAL :: empty = .FALSE.
LOGICAL :: latomtypes = .TRUE.
!
isdel = .FALSE.
istart = dem_molerange(1)                      ! Cope include range into working variables
IF(dem_molerange(2)==-1) THEN
   iend = mole_num_mole
ELSE
   iend   = dem_molerange(2)
ENDIF
!
IF(istart>=1    .AND. istart <= mole_num_mole .AND. &      ! Test include range
   istart<=iend .AND. iend   <= mole_num_mole       ) THEN
   mainloop:DO i=istart, iend                              ! Loop over all included molecules
      moletypes:IF(dem_lmoletype(mole_type(i))) THEN       ! Moletype is selected
         IF(dem_atomrange(1)>0 .OR.                                &
            (dem_atomrange(2)>0 .AND. dem_atomrange(2)<cr_natoms)) THEN  !Limited atom range
            DO j=1, mole_len(i)                            ! Loop over all atoms in molecule
               k = mole_cont(mole_off(i)+j)                ! Absolute atom number
               IF(dem_atomrange(2)>0) THEN                 ! Upper limit is restricted
                  IF(k<dem_atomrange(1) .OR. dem_atomrange(2)<k) EXIT moletypes  ! Atom range outside , skip
               ELSE
                  IF(k<dem_atomrange(1))                         EXIT moletypes  ! Atom range outside , skip
               ENDIF
            ENDDO
         ENDIF
!
! Test if required atom types are present
!
         latomtypes=.TRUE.                                 ! Assume required atom types are there
         atomtypes:DO l=0,cr_nscat
            IF(dem_latomtype(l)) THEN                      ! This atom type is required
               DO j=1, mole_len(i)                               ! Loop over all atoms in molecule
                  k = mole_cont(mole_off(i)+j)                ! Absolute atom number
                  IF(cr_iscat(1,mole_cont(mole_off(i)+j))==l) CYCLE atomtypes  ! correct atom type
               ENDDO
               latomtypes=.FALSE.                          ! End of molecule without atom type
               EXIT atomtypes                              ! no need for further search
            ENDIF
         ENDDO atomtypes
!
! Test property masks
!
         DO j=1, mole_len(i)                               ! Loop over all atoms in molecule
            k = mole_cont(mole_off(i)+j)                   ! Absolute atom number
            IF(.NOT.check_select_status(k, .TRUE., cr_prop(k), dem_sel_prop) ) EXIT moletypes  ! Atom does not fulfil properties
         ENDDO
!
         passed: IF(latomtypes) THEN                       ! Passed atom type test
!
! All tests were passed
!
            DO j=1, mole_len(i)                            ! Loop over all atoms in molecule
               k = mole_cont(mole_off(i)+j)                ! Absolute atom number
               cr_mole(k) = 0
               cr_prop (k ) = IBCLR (cr_prop (k ), PROP_MOLECULE)   ! UNFLAG THIS ATOM AS Molecule  atom
               cr_prop (k ) = IBCLR (cr_prop (k ), PROP_LIGAND  )   ! UNFLAG THIS ATOM AS Ligand Atom
               mole_cont(mole_off(i)+j) = 0                ! Clear molecule content
               isdel = .TRUE.
            ENDDO
         ENDIF passed
      ENDIF moletypes
   ENDDO mainloop
ELSE
   ier_num = -63
   ier_typ = ER_APPL
   ier_msg(1) = 'Molecule include range is outside range'
ENDIF
!
IF(isdel) THEN                     ! A molecule was deleted
   i = 1
   clean: DO 
      IF(i > mole_num_mole) EXIT clean 
      empty = .TRUE.
      isempty: DO j=1, mole_len(i)                       ! Loop over all atoms in molecule
         IF(mole_cont(mole_off(i)+j) > 0) THEN           ! Absolute atom number
            empty = .FALSE.
            EXIT isempty
         ENDIF
      ENDDO isempty
      IF(empty) THEN                                     ! Molecule is empty
         nat = mole_len(i)
         DO j=mole_off(i)+1, mole_num_atom - mole_len(i)                             ! Shift atoms down
            mole_cont(j) = mole_cont(j+mole_len(i))
            k = mole_cont(j)
            IF(k>0) cr_mole(k) = cr_mole(k) - 1
         ENDDO
         DO mmol = i, mole_num_mole-1                    ! Shift molecule properties down
            mole_len (mmol) = mole_len (mmol + 1)
            mole_off (mmol) = mole_off (mmol + 1) - nat
            mole_type(mmol) = mole_type(mmol + 1)
            mole_file(mmol) = mole_file(mmol + 1)
            mole_char(mmol) = mole_char(mmol + 1)
         ENDDO
         mole_num_mole = mole_num_mole - 1
         mole_num_atom = mole_num_atom - nat
      ELSE
         i = i + 1
      ENDIF
   ENDDO clean
ENDIF
!
END SUBROUTINE demol_run
!
!*******************************************************************************
!
SUBROUTINE demol_reset
!
USE demolec_mod
!
IMPLICIT NONE
!
dem_molerange(1) = 0
dem_molerange(2) = 0
dem_atomrange(1) =  1
dem_atomrange(2) = -1
IF(ALLOCATED(dem_lmoletype)) dem_lmoletype(:) = .FALSE.
IF(ALLOCATED(dem_latomtype)) dem_latomtype(:) = .FALSE.
dem_sel_prop(:)  = 0
!
END SUBROUTINE demol_reset
!
!*******************************************************************************
!
END MODULE demolec
