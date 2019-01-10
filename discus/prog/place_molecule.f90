MODULE mole_surf_mod
!
!USE crystal_mod
!
PRIVATE
PUBLIC do_place_molecule, deco_reset
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE do_place_molecule
!
!   Main menu and procedures to place molecules onto a surface
!
   USE deco_mod
   USE discus_allocate_appl_mod
   USE modify_mod
   USE modify_func_mod
USE prop_para_func
USE prop_para_mod
!
   USE calc_expr_mod
   USE doact_mod
   USE do_eval_mod
USE do_wait_mod
USE errlist_mod
   USE get_params_mod
   USE learn_mod
   USE class_macro_internal
   USE prompt_mod
      USE sup_mod
!
   IMPLICIT none
!
   CHARACTER (LEN=5)                       :: befehl! command on input line
   CHARACTER(LEN=LEN(prompt))              :: orig_prompt  ! original prompt
   CHARACTER (LEN=1024)                    :: line  ! input line
   CHARACTER (LEN=1024)                    :: zeile ! remainder with parameters
   INTEGER                                 :: indxg ! location of "="
   INTEGER                                 :: lp    ! lengtz of zeile
   INTEGER laenge, lbef
   INTEGER                                 :: n_num
   LOGICAL                                 :: lalloc = .false. ! Need to allocate
   LOGICAL                                 :: ladd = .true.  ! condition add command is fine
   LOGICAL                                 :: lend  ! condition of EOF
   LOGICAL                                 :: lold =.true. ! condition of EOF
!
   INTEGER, EXTERNAL :: len_str
   LOGICAL, EXTERNAL :: str_comp
!
   lend   = .FALSE.
   lalloc = .FALSE.
!
!  IF INITIALIZATION is needed, set property flags to external surface
!
   IF (dc_init) THEN
      dc_sel_prop(0) = bit_set ( dc_sel_prop, 0, PROP_SURFACE_EXT, .true.)
      dc_sel_prop(1) = bit_set ( dc_sel_prop, 1, PROP_SURFACE_EXT, .true.)
      dc_use_conn(:) = 1 ! while developing!!!
!      dc_n_molecules = 0 ! while developing!!!
      dc_init = .false.
   ENDIF
!
   orig_prompt = prompt
   prompt = prompt (1:len_str (prompt) ) //'/deco'
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
                     .AND..NOT. (str_comp (befehl, 'echo', 2, lbef, 4) )    &
                     .AND..NOT. (str_comp (befehl, 'syst', 2, lbef, 4) )    &
                     .AND..NOT. (str_comp (befehl, 'help', 2, lbef, 4) .OR. &
                                 str_comp (befehl, '?   ', 2, lbef, 4) )    &
                     .AND. INDEX(line,'==') == 0                            ) THEN
!                                                                       
! ------evaluate an expression and assign the value to a variabble      
!                                                                       
               CALL do_math (line, indxg, laenge)
         ELSE                                    ! is_math, al other commands
!                                                                       
!------ ----execute a macro file                                        
!                                                                       
              IF (befehl (1:1) .eq.'@') THEN     ! macro, reset or all other commands
                  IF (laenge.ge.2) THEN 
                     CALL file_kdo (line (2:laenge), laenge-1) 
                  ELSE 
                     ier_num = - 13 
                     ier_typ = ER_MAC 
                  ENDIF 
!                                                                       
!     ----continues a macro 'continue'                                  
!                                                                       
              ELSEIF (str_comp (befehl, 'continue', 3, lbef, 8) ) THEN 
                  CALL macro_continue (zeile, lp) 
!                                                                       
!     ----Echo a string, just for interactive check in a macro 'echo'   
!                                                                       
              ELSEIF (str_comp (befehl, 'echo', 2, lbef, 4) ) THEN 
                  CALL echo (zeile, lp) 
!                                                                       
!      ---Evaluate an expression, just for interactive check 'eval'     
!                                                                       
              ELSEIF (str_comp (befehl, 'eval', 2, lbef, 4) ) THEN 
                  CALL do_eval (zeile, lp, .TRUE.) 
!                                                                       
!     ----exit 'exit'                                                   
!                                                                       
              ELSEIF (str_comp (befehl, 'exit', 2, lbef, 4) ) THEN 
                  lend = .true. 
                  EXIT main_loop
!                                                                       
!     ----help 'help' , '?'                                             
!                                                                       
              ELSEIF (str_comp (befehl, 'help', 1, lbef, 4) .or.  &
                       str_comp (befehl, '?   ', 1, lbef, 4) ) THEN                                      
                  IF (str_comp (zeile, 'errors', 2, lp, 6) ) THEN 
                     lp = lp + 7 
                     CALL do_hel ('discus '//zeile, lp) 
                  ELSE 
                     lp = lp + 12 
                     CALL do_hel ('discus deco '//zeile, lp) 
                  ENDIF 
!                                                                       
!------- -Operating System Kommandos 'syst'                             
!                                                                       
              ELSEIF (str_comp (befehl, 'syst', 2, lbef, 4) ) THEN 
                  IF (zeile.ne.' '.and.zeile.ne.char (13) ) THEN 
                     CALL do_operating (zeile (1:lp), lp) 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
!                                                                       
!------  -----waiting for user input                                    
!                                                                       
              ELSEIF (str_comp (befehl, 'wait', 3, lbef, 4) ) THEN 
                  CALL do_input (zeile, lp) 
!
              ELSEIF (str_comp (befehl, 'reset', 3, lbef, 4)) THEN
                  CALL deco_reset
                  if(ier_num == 0) ladd = .true.   ! allow new add command
              ELSE    ! macro, reset or all other commands
!
!---------------- All other commands
                  IF(MAXSCAT > UBOUND(dc_latom,1)) THEN
                     lalloc = .TRUE.
                  ENDIF
                  IF(lalloc) THEN
                     CALL alloc_deco(MAXSCAT, DCC_MAXNUM, DCC_MAXANCH, DCC_MAXHKL, DCC_MAXNEW, DCC_MAXMSCAT)
                  ENDIF
!
!----------------------------------------------------------------------------------
!     ----Original decorate commands                                      
!----------------------------------------------------------------------------------
!
!
              is_com: IF (str_comp (befehl, 'add', 3, lbef, 3)) THEN
                  IF(ladd) THEN
                     IF(dcc_num == DCC_MAXNUM) THEN
                        n_num = DCC_MAXNUM + 5
                        CALL alloc_deco(MAXSCAT, n_num, DCC_MAXANCH, DCC_MAXHKL, DCC_MAXNEW, DCC_MAXMSCAT)
                     ENDIF
                        IF(ier_num == 0 ) THEN
                           dcc_num = dcc_num + 1
                        ELSE
                           RETURN
                        ENDIF
                     dcc_axis(0,  dcc_num) = -1     ! Set axis default as auto
                     CALL deco_add (zeile, lp)
!                    IF(ier_num == 0) ladd = .false.   ! exclude new add command
                  ELSE
                     ier_num = -128
                     ier_typ = ER_APPL
                     ier_msg(1) = 'Currently only one decoration definition'
                     ier_msg(2) = 'can be used. Run the calculation or reset.'
                  ENDIF
!                                                                       
!------ --Handle property settings 'property'                           
!                                                                       
              ELSEIF (str_comp (befehl, 'property', 4, lbef, 8) ) THEN
!                                                                       
                  CALL property_select (zeile, lp, dc_sel_prop)
!
              ELSEIF (str_comp (befehl, 'run', 3, lbef, 3)) THEN
                 IF(dcc_num > 0 ) THEN
                     CALL deco_run
                     if(ier_num == 0) ladd = .true.   ! allow new add command
                 ELSE
                     ier_num = -129
                     ier_typ = ER_APPL
                 ENDIF
!
              ELSEIF (str_comp (befehl, 'sel', 3, lbef, 3) .or.  &
                      str_comp (befehl, 'del', 3, lbef, 3)     ) THEN
                  CALL atom_select (zeile, lp, 0,  DC_MAXSCAT,  dc_latom, &
                  dc_sel_atom , lold  ,   &
                  str_comp (befehl, 'sel', 2, lbef, 3) )
!
              ELSEIF (str_comp (befehl, 'set', 3, lbef, 3)) THEN
!                IF(dc_temp_type /= DC_NONE) THEN
                 IF(dcc_num > 0 ) THEN
                    CALL deco_set (zeile, lp)
                 ELSE
                    ier_num = -129
                    ier_typ = ER_APPL
                    ier_msg(1) = 'Before setting details, a decoration'
                    ier_msg(2) = 'must have been defined via => add'
                 ENDIF
!
!
              ELSEIF (str_comp (befehl, 'show', 3, lbef, 4)) THEN
                  CALL deco_show
                  WRITE(output_io,*)
!
              ELSE is_com
                 ier_num = -8
                 ier_typ = ER_COMM
!
              ENDIF is_com ! END IF BLOCK actual commands
              ENDIF
           ENDIF is_math   ! END IF BLOCK math equation or specific command
        ENDIF no_com       ! END IF BLOCK no comment
      ENDIF no_err         ! END IF BLOCK no error reading input
!
      IF (ier_num.ne.0) THEN 
         CALL errlist 
         IF (ier_sta.ne.ER_S_LIVE) THEN 
            IF (lmakro .OR. lmakro_error) THEN  ! Error within macro or termination errror
               IF(sprompt /= prompt ) THEN
                  ier_num = -10
                  ier_typ = ER_COMM
                  ier_msg(1) = ' Error occured in decorate menu'
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
!
!
   ENDDO main_loop     ! END DO main loop of menu 
!
   prompt = orig_prompt
!
   END SUBROUTINE do_place_molecule
!
!######################################################################
!
   SUBROUTINE deco_set (zeile,lp)
!
!  Set the definitions for the molecule decorator
!
   USE atom_env_mod
   USE chem_mod
   USE deco_mod
   USE discus_allocate_appl_mod
   USE get_iscat_mod
   USE modify_mod
   USE point_grp
!
   USE ber_params_mod
   USE get_params_mod
   USE take_param_mod
!
   IMPLICIT NONE
!
!
   CHARACTER (LEN=*), INTENT(INOUT) :: zeile
   INTEGER          , INTENT(INOUT) :: lp
!
   INTEGER, PARAMETER   :: MAXW = 20
   INTEGER, PARAMETER :: NOPTIONAL = 2
   CHARACTER (LEN=1024), DIMENSION(1:MAXW) :: cpara
   CHARACTER (LEN=1024), DIMENSION(1:2)    :: ccpara
   INTEGER             , DIMENSION(1:MAXW) :: lpara
   INTEGER             , DIMENSION(1:2)    :: llpara
!
   INTEGER, PARAMETER :: O_ANGLE  = 1
   INTEGER, PARAMETER :: O_ANCHOR = 2
   CHARACTER(LEN=   6), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
   CHARACTER(LEN=1024), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
   INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
   INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
   REAL               , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
   INTEGER, PARAMETER                        :: ncalc = 2 ! Number of values to calculate

!
   INTEGER             , DIMENSION(:,:), ALLOCATABLE :: temp_hkl
!  INTEGER             , DIMENSION(:,:), ALLOCATABLE :: temp_hkl_new
!  INTEGER             , DIMENSION(1:4)    :: hkl
   INTEGER             , DIMENSION(1:4)    :: hkl_new
   INTEGER              :: temp_nhkl   ! current number of hkl for restictions
   INTEGER              :: temp_num    ! current decoration number
   INTEGER              :: n_anch      ! current number of anchor atom types
   INTEGER              :: n_hkl       ! current number of hkl for restictions
   INTEGER              :: n_surfnew   ! current number of new surface type atoms
   INTEGER              :: i,j         ! Dummy variables
   INTEGER              :: ianz, janz, kanz, success  !, ncon
   LOGICAL              :: lnew !, lnew_mole = .false.
   LOGICAL              :: l_form
   REAL   , DIMENSION(1:MAXW) :: werte
   REAL   , DIMENSION(1:2   ) :: wwerte
!
   LOGICAL str_comp
!
   DATA oname  / 'angle ', 'anchor' /
   DATA loname /  6      ,  6       /
   opara  =  (/ '170.00' , '1.0000' /)    ! Always provide fresh default values
   lopara =  (/  6       ,  6       /)
   owerte =  (/  170.00  ,  1.0000  /)
!
   CALL get_params(zeile, ianz, cpara, lpara, maxw, lp)
   IF ( ier_num /= 0 ) RETURN              ! Error reading parameters
   IF ( ianz <  3 ) THEN                   ! All commands need three parameters
      ier_num = -6
      ier_typ = ER_COMM
      ier_msg(1) = 'All deco set commands need at least three '
      ier_msg(2) = 'parameters'
      RETURN
   ENDIF
!
!  Get optional parameters
!
!  opara  =  (/ '170.00'/)    ! Always provide fresh default values
!  lopara =  (/  6      /)
!  owerte =  (/  170.00 /)
   CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                     oname, loname, opara, lopara, owerte)
   IF(ier_num /= 0) RETURN
   success       = 1
   temp_num      = 0
   locate:DO j=1,dcc_num
      IF(cpara(1)(1:lpara(1)) == dcc_name(j)(1:dcc_lname(j))) THEN
         success = 0
         temp_num = j
         EXIT locate
      ENDIF
   ENDDO locate
   IF(success /= 0) THEN
      ier_num = -145
      ier_typ = ER_APPL
      RETURN
   ENDIF
!
   lnew          = .false.
   success       = 1
   IF ( str_comp(cpara(2),'axis',4,lpara(2),4) ) THEN
      IF ( ianz == 4 ) THEN
         CALL del_params (2, ianz, cpara, lpara, maxw)   ! delete first 2 params
         CALL ber_params (ianz, cpara, lpara, werte, maxw)
         IF(ier_num ==0) THEN
            dcc_axis(0,dcc_num) = 2
            dcc_axis(1,dcc_num) = NINT(werte(1))
            dcc_axis(2,dcc_num) = NINT(werte(2))
         ENDIF
      ELSEIF ( ianz == 3 ) THEN
         CALL del_params (2, ianz, cpara, lpara, maxw)   ! delete first 2 params
         IF( str_comp(cpara(1),'auto',4,lpara(1),4)) THEN
            dcc_axis(0,  dcc_num) = -1
            dcc_axis(1:2,dcc_num) =  0
         ELSE
            ier_num = -6
            ier_typ = ER_COMM
            ier_msg(1) = 'set axis command needs parameters:'
            ier_msg(2) = '''auto'' or two atom numbers '
         ENDIF
      ELSE
         ier_num = -6
         ier_typ = ER_COMM
         ier_msg(1) = 'set axis command needs four parameters'
      ENDIF
   ELSEIF ( str_comp(cpara(2),'bond',4,lpara(2),4) ) THEN
      IF ( ianz >= 5 ) THEN
         CALL del_params (2, ianz, cpara, lpara, maxw)   ! delete first 2 params
         IF(ier_num /= 0) RETURN
         ccpara(1) = cpara(ianz-1)
         ccpara(2) = cpara(ianz  )
         llpara(1) = lpara(ianz-1)
         llpara(2) = lpara(ianz  )
         janz = 2
         CALL ber_params (janz, ccpara, llpara, wwerte, 2)
         IF(ier_num /= 0) RETURN
!
         IF(dcc_neig(0,temp_num) == 0) THEN
            j = 1
         ELSEIF(dcc_neig(0,temp_num) == 1) THEN
            j = 2
         ELSE
            RETURN
         ENDIF
         dcc_neig(0,  temp_num) = j
         dcc_surf(0,0,temp_num) = j
         dcc_neig(j,temp_num) = NINT(wwerte(1))
         dcc_dist(j,temp_num) =      wwerte(2)
         cpara(ianz-1) = ' '
         cpara(ianz  ) = ' '
         lpara(ianz-1) = 1
         lpara(ianz  ) = 1
         janz = ianz - 2                                 ! ignore last two parameters
         kanz = ianz - 2                                 ! ignore last two parameters
         CALL get_iscat (janz, cpara, lpara, werte, maxw, .false.)
         IF(ier_num /= 0) RETURN
!        IF(dcc_type(temp_num) == DC_NORMAL .AND. kanz /= 1) THEN
!           ier_num = -143
!           ier_typ = ER_APPL
!           ier_msg(1) = 'For the NORMAL decoration we need '
!           ier_msg(2) = 'exactly one surface atom type '
!        ELSEIF(dcc_type(temp_num) == DC_DOUBLE .AND. kanz /= 1) THEN
!           ier_num = -143
!           ier_typ = ER_APPL
!           ier_msg(1) = 'For the DOUBLE decoration we need '
!           ier_msg(2) = 'exactly one surface atom type '
!           ier_msg(3) = 'and two bond statements'
!        ELSEIF(dcc_type(temp_num) == DC_BRIDGE .AND. kanz /= 1) THEN
!           ier_num = -143
!           ier_typ = ER_APPL
!           ier_msg(1) = 'For the BRIDGE decoration we need '
!           ier_msg(2) = 'exactly one surface atom types'
!           ier_msg(3) = 'and two bond statements'
!        ELSEIF(dcc_type(temp_num) == DC_ACCEPTOR .AND. kanz /= 1) THEN
!           ier_num = -143
!           ier_typ = ER_APPL
!           ier_msg(1) = 'For the ACCEPTOR decoration we need '
!           ier_msg(2) = 'exactly one surface atom type '
!        ELSEIF(dcc_type(temp_num) == DC_DONOR .AND. kanz /= 1) THEN
!           ier_num = -143
!           ier_typ = ER_APPL
!           ier_msg(1) = 'For the DONOR decoration we need '
!           ier_msg(2) = 'exactly one surface atom type '
!        ELSE ! SUCCESS
            IF(janz > UBOUND(dcc_surf,2)) THEN
               n_anch = janz+4
               CALL alloc_deco(MAXSCAT, DCC_MAXNUM, n_anch     , DCC_MAXHKL, DCC_MAXNEW, DCC_MAXMSCAT)
               IF(ier_num /= 0) RETURN
            ENDIF
            DO i=1,janz
               dcc_surf(j,i,temp_num) = NINT(werte(i))  ! Surface atom type
            ENDDO
            dcc_surf(j,0,temp_num) = janz               ! Number of surface atom types
            dcc_maxsurf = MAX(dcc_maxsurf, dcc_surf(j,0,temp_num)) ! Keep maximum number
            dcc_nanch(j, temp_num) = NINT(owerte(O_ANCHOR))
            IF(dcc_type(temp_num) == DC_MULTIPLE) THEN     ! Special treatment of 'multi'
               IF(j==1) THEN                            ! First 'set bond' command
                  IF(NINT(owerte(O_ANCHOR))==1) THEN    ! User did not give anchor:3
                     IF(kanz==3) THEN                   ! number of surface atom types == 3
                        dcc_nanch(j, temp_num) = 3
                     ELSE
                        ier_num = -143
                        ier_typ = ER_APPL
                        ier_msg(1) = 'For the MULTI decoration we need '
                        ier_msg(2) = 'exactly three surface atom types'
                        ier_msg(3) = 'or the ''anchor:3'' parameter'
                     ENDIF
                  ELSEIF(NINT(owerte(O_ANCHOR))==2) THEN ! number of surface atom types == 2
                     ier_num = -143
                     ier_typ = ER_APPL
                     ier_msg(1) = 'For the MULTI decoration we need '
                     ier_msg(2) = 'exactly three surface atom types'
                     ier_msg(3) = 'or the ''anchor:3'' parameter'
                  ELSEIF(NINT(owerte(O_ANCHOR))> 3) THEN ! number of surface atom types >  3
                     ier_num = -143
                     ier_typ = ER_APPL
                     ier_msg(1) = 'For the MULTI decoration we need '
                     ier_msg(2) = 'exactly three surface atom types'
                     ier_msg(3) = 'or the ''anchor:3'' parameter'
                  ENDIF
               ELSE                                      ! Second 'set bond' command
                  IF(NINT(owerte(O_ANCHOR))==1) THEN    ! User did not give anchor:3
                     dcc_nanch(j, temp_num) = 1
                  ELSE
                     ier_num = -143
                     ier_typ = ER_APPL
                     ier_msg(1) = 'For the MULTI second bond we need '
                     ier_msg(2) = 'exactly one surface atom type'
                     ier_msg(3) = 'or the ''anchor:1'' parameter'
                  ENDIF
               ENDIF
            ENDIF
!
!  Interpret optional values, as default values are provided, we can take it blindly
!
            dcc_angle(temp_num) = owerte(O_ANGLE)
!
!        ENDIF
      ELSE
         ier_num = -6
         ier_typ = ER_COMM
         ier_msg(1) = 'set bond command needs >= five parameters'
      ENDIF
   ELSEIF ( str_comp(cpara(2),'ligand',4,lpara(2),6) ) THEN
      IF ( ianz == 4 ) THEN
         dcc_file (temp_num) = cpara(3)
         dcc_lfile(temp_num) = lpara(3)
         CALL del_params (3, ianz, cpara, lpara, maxw)   ! delete first 3 params
         IF(ier_num /= 0) RETURN
         CALL ber_params (ianz, cpara, lpara, werte, maxw)
         IF(ier_num /= 0) RETURN
         dcc_dens (temp_num) = werte(1)
!
      ELSE
         ier_num = -6
         ier_typ = ER_COMM
         ier_msg(1) = 'set ligand command needs four parameters'
      ENDIF
   ELSEIF ( str_comp(cpara(2),'hkl',3,lpara(2),3) .OR.    &
            str_comp(cpara(2),'form',4,lpara(2),4) ) THEN
      dcc_lform(temp_num) = str_comp(cpara(2),'form',4,lpara(2),4)
      l_form = str_comp(cpara(2),'form',4,lpara(2),4)  ! Obsolete
      IF ( ianz == 3 ) THEN
         IF( str_comp(cpara(3),'none',4,lpara(2),4) ) THEN
            dcc_lrestrict(temp_num) = .FALSE.
            dcc_hkl(:,:,temp_num)   = 0
         ELSE
            ier_num = -6
            ier_typ = ER_COMM
            ier_msg(1) = 'set hkl command needs four parameters'
            ier_msg(2) = 'or must be                           '
            ier_msg(2) = 'set hkl, none                        '
         ENDIF
      ELSEIF ( ianz == 5 ) THEN
         CALL del_params (2, ianz, cpara, lpara, maxw)   ! delete first 2 params
         CALL ber_params (ianz, cpara, lpara, werte, maxw)
         IF( ier_num==0) THEN
            hkl_new(:)   = 0
            hkl_new(1:3) = NINT(werte(1:3))
            temp_nhkl    = 1
         ENDIF
      ELSE
         ier_num = -6
         ier_typ = ER_COMM
         ier_msg(1) = 'set hkl command needs 3 or 5 parameters'
      ENDIF
      IF(ier_num==0) then
         IF(ALLOCATED(temp_hkl)) THEN
            DEALLOCATE(temp_hkl)
         ENDIF
         ALLOCATE(temp_hkl(3,48))
         temp_nhkl = 1
         IF(l_form) THEN    ! If necessary expand form to symmetrically equivalent hkl
            CALL point_init(hkl_new, temp_hkl, temp_nhkl)
            CALL point_set(hkl_new, temp_hkl, temp_nhkl)
         ELSE
            temp_hkl(:,1) = hkl_new(1:3)
         ENDIF
         IF(dcc_hkl(1,0,temp_num)+temp_nhkl  > UBOUND(dcc_hkl,2)) THEN
            n_hkl = dcc_hkl(j,0,temp_num)+temp_nhkl
            CALL alloc_deco(MAXSCAT, DCC_MAXNUM, DCC_MAXANCH, n_hkl    , DCC_MAXNEW, DCC_MAXMSCAT)
            IF(ier_num /= 0) RETURN
         ENDIF
         j = dcc_hkl(j,0,temp_num)
         DO i=1,temp_nhkl
            dcc_hkl(:, j+i, temp_num)    = temp_hkl(:,i)
         ENDDO
         dcc_hkl(1,0,temp_num)  = dcc_hkl(1,0,temp_num) + temp_nhkl
         dcc_lrestrict(temp_num) = .TRUE.
      ENDIF
   ELSEIF ( str_comp(cpara(2),'surface',4,lpara(2),7) ) THEN
      CALL del_params (2, ianz, cpara, lpara, maxw)   ! delete first 2 params
      IF(ianz>UBOUND(dc_temp_surfnew,1)) THEN
         ier_num = -6999
         ier_typ = ER_COMM
         ier_msg(1) = 'The number of surface atoms is limited to'
         WRITE(ier_msg(2),'(i2)') UBOUND(dc_temp_surfnew,1)
         RETURN
      ENDIF
      CALL ber_params (ianz, cpara, lpara, werte, maxw)
      IF(ier_num/=0) RETURN
      IF(dcc_surfnew(0,temp_num)+ianz>UBOUND(dcc_surfnew,1)) THEN
         n_surfnew = dcc_surfnew(0,temp_num)+ianz
         CALL alloc_deco(MAXSCAT, DCC_MAXNUM, DCC_MAXANCH, DCC_MAXHKL, n_surfnew, DCC_MAXMSCAT)
         IF(ier_num /= 0) RETURN
      ENDIF
      i = dcc_surfnew(0,temp_num)
      dcc_surfnew(i+1:i+ianz,temp_num) = NINT(werte(1:ianz))
      dcc_surfnew(0,temp_num) = dcc_surfnew(0,temp_num) + ianz
   ELSEIF ( str_comp(cpara(2),'tilt',4,lpara(2),4) ) THEN
      CALL del_params (2, ianz, cpara, lpara, maxw)   ! delete first 2 params
      IF( str_comp(cpara(1),'angle',4,lpara(1),5) ) THEN
         CALL del_params (1, ianz, cpara, lpara, maxw)   ! delete first 2 params
         IF(ianz==1) THEN
            CALL ber_params (ianz, cpara, lpara, werte, maxw)
            IF(ier_num/=0) RETURN
            dcc_tilt(temp_num) = werte(1)
         ELSE
            ier_num = -6
            ier_typ = ER_COMM
         ENDIF
      ELSEIF( str_comp(cpara(1),'plane',4,lpara(1),5) ) THEN
         CALL del_params (1, ianz, cpara, lpara, maxw)   ! delete first 2 params
         IF(ianz==1) THEN
            IF(str_comp(cpara(1), 'auto', 4, lpara(1), 4)) THEN
               dcc_tilt_hkl (1:3,temp_num) = -1.0
               dcc_tilt_atom(1:4,temp_num) = -1
               dcc_tilt_is_auto(temp_num) = .TRUE.
            ELSE
               ier_num = -6
               ier_typ = ER_COMM
            ENDIF
         ELSEIF(ianz==3) THEN
            CALL ber_params (ianz, cpara, lpara, werte, maxw)
            IF(ier_num/=0) RETURN
            dcc_tilt_hkl (1:3,temp_num) = werte(1:3)
            dcc_tilt_atom(1:4,temp_num) = -1
            dcc_tilt_is_auto(temp_num) = .FALSE.
         ELSE
            ier_num = -6
            ier_typ = ER_COMM
         ENDIF
      ELSEIF( str_comp(cpara(1),'atoms',4,lpara(1),5) ) THEN
         CALL del_params (1, ianz, cpara, lpara, maxw)   ! delete first 2 params
         IF(ianz==1) THEN
            IF(str_comp(cpara(1), 'auto', 4, lpara(1), 4)) THEN
               dcc_tilt_hkl (1:3,temp_num) = -1.0
               dcc_tilt_atom(1:4,temp_num) = -1
               dcc_tilt_is_auto(temp_num) = .TRUE.
            ELSE
               ier_num = -6
               ier_typ = ER_COMM
            ENDIF
         ELSEIF(ianz==4) THEN
            CALL ber_params (ianz, cpara, lpara, werte, maxw)
            IF(ier_num/=0) RETURN
            dcc_tilt_hkl (1:3,temp_num) = -1.0
            dcc_tilt_atom(1:4,temp_num) = NINT(werte(1:4))
            dcc_tilt_is_auto(temp_num) = .FALSE.
         ELSE
            ier_num = -6
            ier_typ = ER_COMM
         ENDIF
      ENDIF
   ELSE
      ier_num = -6
      ier_typ = ER_COMM
   ENDIF
!
   IF(success == -1) THEN
      ier_num = -128
      ier_typ = ER_APPL
   ENDIF
!
   END SUBROUTINE deco_set
!
!######################################################################
!
   SUBROUTINE deco_add (zeile,lp)
!
!  Add a   definitions for the molecule decorator
!
   USE atom_env_mod
   USE chem_mod
USE deco_mod
   USE modify_mod
USE errlist_mod
   USE get_params_mod
!
   IMPLICIT NONE
!
!
   CHARACTER (LEN=*), INTENT(INOUT) :: zeile
   INTEGER          , INTENT(INOUT) :: lp
!
   INTEGER, PARAMETER   :: MAXW = 20
   CHARACTER (LEN=1024) :: cpara(1:MAXW)
   INTEGER              :: lpara(1:MAXW)
   INTEGER              :: ianz, success
!  LOGICAL              :: lnew
!
   LOGICAL str_comp
!
   CALL get_params(zeile, ianz, cpara, lpara, maxw, lp)
   IF ( ier_num /= 0 ) RETURN              ! Error reading parameters
!
   IF ( ianz == 2 ) THEN
      dcc_name(dcc_num)  = cpara(1)(1:lpara(1))
      dcc_lname(dcc_num) = lpara(1)
      success       = 1
      IF ( str_comp(cpara(2),'normal',3,lpara(2),6) ) THEN
         dcc_type(dcc_num) = DC_NORMAL
      ELSEIF ( str_comp(cpara(2),'bridge',3,lpara(2),6) ) THEN
         dcc_type(dcc_num) = DC_BRIDGE
      ELSEIF ( str_comp(cpara(2),'chelate',3,lpara(2),7) ) THEN
         dcc_type(dcc_num) = DC_CHELATE
      ELSEIF ( str_comp(cpara(2),'double',3,lpara(2),6) ) THEN
         dcc_type(dcc_num) = DC_DOUBLE
      ELSEIF ( str_comp(cpara(2),'multiple',3,lpara(2),8) ) THEN
         dcc_type(dcc_num)      = DC_MULTIPLE
      ELSEIF ( str_comp(cpara(2),'acceptor',3,lpara(2),8) ) THEN
         dcc_type(dcc_num) = DC_ACCEPTOR
      ELSEIF ( str_comp(cpara(2),'donor',3,lpara(2),5) ) THEN
         dcc_type(dcc_num) = DC_DONOR
      ELSE
         ier_num = -6
         ier_typ = ER_COMM
         ier_msg(1) = 'Allowed decoration types are:'
         ier_msg(2) = 'normal, bridge, double, multi, acceptor'
         ier_msg(3) = 'donor , chelate' 
!        NULLIFY(dc_def_temp) ! 
         RETURN
      ENDIF
      dcc_file   (    dcc_num) = ' '
      dcc_lfile  (    dcc_num) = 0
      dcc_surf   (:,:,dcc_num) = 0
      dcc_neig   (:,  dcc_num) = 0
      dcc_secnd  (    dcc_num) = 0
      dcc_axis   (:,  dcc_num) = 0
      dcc_axis   (0,  dcc_num) = -1
      dcc_lform  (    dcc_num) = .FALSE.
      dcc_hkl    (:,:,dcc_num) = 0
      dcc_surfnew(:,  dcc_num) = 0
      dcc_dens   (    dcc_num) = 0.0
      dcc_dist   (:,  dcc_num) = 0.0
   ELSE
      ier_num = -6
      ier_typ = ER_COMM
   ENDIF
!
   END SUBROUTINE deco_add
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   SUBROUTINE deco_run
!
!  Performs the actual decoration
!
USE crystal_mod
   USE chem_menu
   USE conn_mod
   USE conn_def_mod
   USE chem_aver_mod
   USE chem_mod
USE deco_mod
   USE discus_allocate_appl_mod
   USE discus_plot_init_mod
   USE discus_save_mod
   USE domain_menu
   USE domain_mod
   USE micro_mod
   USE mmc_menu
   USE mc_mod
   USE mmc_mod
   USE modify_mod
   USE modify_func_mod
   USE molecule_mod
USE prop_para_func
   USE prop_para_mod
   USE read_internal_mod
   USE structur , ONLY: rese_cr
   USE save_menu, ONLY: save_internal, save_store_setting, save_restore_setting, save_default_setting, save_struc, save_show
!USE surface_func_mod
   USE discus_init_mod, ONLY: mmc_init
!
   USE param_mod
   USE random_mod
USE prompt_mod
!
   IMPLICIT none
!
   INTEGER, PARAMETER         :: MAXW = 200  ! not ideal, should be dynamic ....
   REAL   , PARAMETER         :: EPS = 0.00001 ! Accuracy to locate anchor
!
   CHARACTER(LEN=1024) :: line      ! a string
   CHARACTER(LEN=1024) :: mole_name ! molecule file name
   CHARACTER(LEN= 200) :: corefile  ! original structure file name
!  CHARACTER(LEN= 200) :: corelist  ! Place holder for core position at 0,0,0
   CHARACTER(LEN= 200) :: shellfile  ! original structure file name
   CHARACTER(LEN=1024), DIMENSION(MAXW) :: cpara      ! a string
   INTEGER            , DIMENSION(MAXW) :: lpara      ! a string
   REAL               , DIMENSION(MAXW) :: werte      ! a string
!
   INTEGER   :: ier_num_deco, ier_typ_deco
   CHARACTER(LEN=LEN(ier_msg)), DIMENSION(1:3) :: ier_msg_deco

   INTEGER   :: istatus          ! status
   INTEGER   :: i,j, k, i1  ! dummy index
   INTEGER   :: ia   ! dummy index
   INTEGER   :: ianz ! dummy number of parameters
   INTEGER   :: length ! dummy length
   INTEGER   :: mole_length ! length of molecule file name
   INTEGER   :: istart,iend ! dummy index
   INTEGER   :: is, js ! scattering number of surface atom
!  INTEGER   :: idef ! connectivity definition number
   INTEGER   :: jdef ! connectivity definition number
   INTEGER   :: ncon ! number of connections defined
!   INTEGER, DIMENSION(:), ALLOCATABLE   :: c_list
   INTEGER, DIMENSION(:), ALLOCATABLE   :: temp_iatom  ! temporary storage for host atom number
   INTEGER, DIMENSION(:), ALLOCATABLE   :: temp_iscat  ! temporary storage for anchor types
   REAL   , DIMENSION(:,:), ALLOCATABLE :: temp_pos    ! temporary storage for anchor positions
   INTEGER, DIMENSION(:), ALLOCATABLE   :: temp_ident  ! temporary storage for anchor types
   INTEGER, DIMENSION(:,:), ALLOCATABLE :: anch_id  ! Lookup which anchor belongs to which env and surface 
INTEGER  :: dc_temp_natoms   ! Number of atoms in the ligand molecule
   INTEGER   :: natoms                 ! number of atoms in connectivity
   INTEGER   :: natoms_prior           ! number of atoms in shell prior to decoration
   INTEGER   :: n_scat                 ! Dummy maximum scattering types
   INTEGER   :: n_corr                 ! Dummy maximum correlation number
   INTEGER   :: nscat_old              ! maximum scattering types prior to modifying anchors
   INTEGER   :: nscat_tmp              ! maximum scattering types prior to current anchors
   INTEGER   :: n_anch                 ! Number of different anchor types
   INTEGER   :: n_repl, i_repl         ! counter for anchor atoms
   INTEGER   :: m_type_old             ! Molecule types in original crystal
   INTEGER   :: j_surf                 ! Surface atom type replaced by an achor
   INTEGER   :: temp_secnd             ! Second neighbor in the molecule
   LOGICAL   :: l_correct              ! A dummy for logical comparisons
   LOGICAL   :: temp_lrestrict         ! Local copy of dcc_lrestict
   REAL      :: temp_angle             ! Local copy of dcc_angle
INTEGER :: nanch
INTEGER, DIMENSION(:,:), ALLOCATABLE :: anchor
INTEGER, DIMENSION(:  ), ALLOCATABLE :: anchor_num
   REAL   , DIMENSION(:), ALLOCATABLE :: temp_prob ! Relative probabilities for a definition
   REAL                    :: r_anch   ! relative amount of anchor atoms
   REAL                    :: temp_grand ! sum of all definition probabilities
!  REAL                    :: temp_choose ! Actual random number for definition choice
   REAL                    :: density  ! Average ligand density
   REAL                    :: prob     ! Probability to replace by non_anchor dummy
   REAL   , DIMENSION(1:3) :: xyz
   REAL   , DIMENSION(3)   :: host_a0
   REAL   , DIMENSION(3)   :: host_win
   REAL   , DIMENSION(4,4) :: host_tran_fi
!
   REAL ran1
!
   IF(MAXVAL(cr_surf(0,:)) == 0 .AND. MINVAL(cr_surf(0,:)) == 0) THEN
      ier_num = -130
      ier_typ = ER_APPL
      ier_msg(1) = 'The structure needs to be cut in surface menu'
      ier_msg(2) = 'Use the '' boundary '' command prior to DECO'
      ier_msg(3) = 'Check use of ''set external'' in ''property'' '
      RETURN
   ENDIF
   IF(cr_natoms == 0 ) THEN
      ier_num = -27
      ier_typ = ER_CHEM
      RETURN
   ENDIF
   IF(ALLOCATED(anch_id)) THEN
      DEALLOCATE(anch_id)
   ENDIF
   ALLOCATE(anch_id(1:dcc_num*dcc_maxsurf,1:2))   ! make lookup for anchors
   anch_id(:,:) = 0
!  ALLOCATE(anch_id(1:dc_temp_id*dc_temp_maxsurf,1:2))   ! make lookup for anchors
   temp_grand = 1.0
!
!  Section for automatic distribution of all surface anchors
!
   m_type_old = mole_num_type          ! Remember molecule type in original crystal
!
   CALL save_store_setting             ! Backup user "save" setting
   CALL save_default_setting           ! Default to full saving
   line       = 'ignore, all'          ! Ignore all properties
   length     = 11
   CALL property_select(line, length, sav_sel_prop)
!
   corefile   = 'internal.decorate'             ! internal user files always start with 'internal'
   shellfile  = 'internal.decoshell'   
   CALL save_internal(corefile)        !     thus this file name is unique
   line       = 'present, external'    ! Force atom to be close to a surface
   length     = 17
   CALL property_select(line, length, sav_sel_prop)
   line       = 'absent, outside'      ! Force atom to be inside
   length     = 15
   CALL property_select(line, length, sav_sel_prop)
      sav_w_scat  = .FALSE.
      sav_w_adp   = .FALSE.
      sav_w_occ   = .FALSE.
      sav_r_ncell = .TRUE.
      sav_w_ncell = .TRUE.
      sav_w_gene  = .FALSE.
      sav_w_symm  = .FALSE.
      sav_w_mole  = .FALSE.
      sav_w_obje  = .FALSE.
      sav_w_doma  = .FALSE.
      sav_w_prop  = .TRUE.
!
   CALL save_internal(shellfile)
! RBN DECO NEEDS ERROR CHECK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Make single atom "structure" for domain list file 
! Should be obsolete
   CALL rese_cr
!  cr_natoms    = 1
!  cr_ncatoms   = 1
!  cr_nscat     = 1
!  cr_pos(:,1)  = 0.0
!  cr_iscat(1)  = 1
!  cr_at_lis(1) = 'CORE'
!  corelist     = 'internal.core.list'
!  line       = 'ignore, all'          ! Ignore all properties
!  length     = 11
!  CALL property_select(line, length, sav_sel_prop)
!  CALL save_internal(corelist)        ! Save the core list
!
!     Load the molecules into temporary structures to reduce disk I/O
!
   host_a0 (:) = cr_a0(:)
   host_win(:) = cr_win(:)
   host_tran_fi(:,:) = cr_tran_fi(:,:)
   CALL deco_get_molecules(host_a0,host_win, host_tran_fi)
   IF(ier_num /=0) THEN
      i = ier_num                         ! Keep error status
      j = ier_typ
      CALL save_restore_setting
      CALL no_error
      CALL readstru_internal( corefile)   ! Read  core file
      ier_num = i                         ! Restore error status
      ier_typ = j
      DEALLOCATE(anch_id)
      RETURN
   ENDIF
!
!     Determine average density
!
   IF(dcc_num > 0) ALLOCATE(temp_prob(0:dcc_num), STAT = istatus)
   temp_prob(:) = 0.0
   density      = 0.0
   DO i=1,dcc_num
      density = density + dcc_dens(i)
      temp_prob(i) = temp_prob(i-1) + dcc_dens(i)
      temp_grand   = temp_grand     + dcc_dens(i)
   ENDDO
   density = density/FLOAT(dcc_num)
   DO i=1,dcc_num
      temp_prob(i) = temp_prob(i)/temp_grand
   ENDDO
!
!
CALL rese_cr
CALL no_error
!
CALL readstru_internal(shellfile)   ! Read shell file
!
cr_icc(:) = 1                       ! The shell is not periodic, 
cr_ncatoms = cr_natoms              ! place all atoms into one "unit" cell
nscat_old = cr_nscat                ! Save old scattering curve number
chem_period(:) = .FALSE.
chem_quick     = .FALSE.
IF(cr_natoms > 0) THEN              ! The Shell does consist of atoms
!  IF(ASSOCIATED(dc_def_head)) THEN
!
!     Now sort those surface atoms that are anchors to the ligand
!
      n_anch = 0                                           !Initialize number of anchors
      n_repl = 0
   name_anchors: DO dc_temp_id=1, dcc_num                  ! Loop over all definitions to name anchors
     dc_temp_dens = dcc_dens(  dc_temp_id)                 ! Get density for surface coverage
     dc_temp_dist = dcc_dist(1,dc_temp_id)
     j = dcc_surf(1,0,dc_temp_id)
     dc_temp_surf(0:j) = dcc_surf(1,0:j,dc_temp_id)
         CALL chem_elem(.false.)                           ! get composition
         r_anch = 0.0
         DO j = 1, dc_temp_surf(0)
            r_anch = r_anch + res_para(dc_temp_surf(j)+1)  ! Fractional composition of the anchoring atoms
         ENDDO
         prob   = MAX(0.0,MIN(1.0,DC_AREA*dc_temp_dens/r_anch)) ! replacement probability
!        Replace anchors by a new atom type
         nscat_tmp = cr_nscat
         IF(cr_nscat+dc_temp_surf(0) > MAXSCAT) THEN                   ! Number of scattering types increased
            n_scat = MAX(cr_nscat+MAX(5,dc_temp_surf(0)), MAXSCAT)     ! Allow for several anchor types
            natoms = MAX(cr_natoms, NMAX)
            CALL alloc_crystal(n_scat,natoms)
         ENDIF
         DO k = 1, dc_temp_surf(0)
            n_anch = n_anch + 1
            WRITE(cr_at_lis(nscat_tmp + k),'(a2,i2.2)') 'AN', n_anch  ! Fixed name for anchor, numbered
            cr_dw    (nscat_tmp + k) = cr_dw(dc_temp_surf(1))
            anch_id(n_anch,1) = dc_temp_id          ! This anchor belongs to environment dc_temp_id
            anch_id(n_anch,2) = dc_temp_surf(k)     ! This anchor connects at atom type  dc_temp_surf(k)
         ENDDO
         cr_nscat = nscat_tmp + dc_temp_surf(0)
!        As surface atom number is bound to be small try several times until we get sufficient anchors
         i = 0
         i_repl = 0
!
         replace: DO i=1, cr_natoms 
            DO ia = 1, cr_natoms                        ! Loop to replace
               l_correct = .FALSE.
               find_surf: DO j=1,dc_temp_surf(0)
                  IF(cr_iscat(ia) == dc_temp_surf(j)) THEN
                     l_correct = .TRUE.
                     j_surf = j
                     EXIT find_surf
                  ENDIF
               ENDDO find_surf
               IF(l_correct                    ) THEN   ! Got a surface atom of correct type
               IF(ran1(idum) < prob) THEN               ! Randomly pick a fraction
                     cr_iscat(ia) = nscat_tmp + j_surf  ! Replace by standard Anchor type
                     cr_prop (ia) = IBSET (cr_prop (ia), PROP_DECO_ANCHOR)    ! FLAG THIS ATOM AS SURFACE ANCHOR
                     n_repl       = n_repl  + 1         ! Increment replaced atoms, grand total
                     i_repl       = i_repl  + 1         ! Increment replaced atoms, local count
                     IF(i_repl == NINT(cr_natoms*r_anch*prob) .OR. &
                        i_repl == cr_natoms                       ) EXIT replace   ! got enough anchors
                  ENDIF
               ENDIF
            ENDDO
         ENDDO replace
      ENDDO name_anchors
!
IF(n_repl==0) THEN
   CALL save_restore_setting
   CALL no_error
   CALL readstru_internal( corefile)   ! Read  core file
   ier_num = -131
   ier_typ = ER_APPL
   ier_msg(1) = 'Is the surface very small, just a few atoms?'
   ier_msg(2) = 'Is the coverage too small? '
   ier_msg(3) = 'Check the set ligand command A1'
   DEALLOCATE(anch_id)
   RETURN
ENDIF
!
   line       = 'ignore, all'          ! Ignore all properties
   length     = 11
   CALL property_select(line, length, sav_sel_prop)
   CALL save_internal('internal_anchors')        !     thus this file name is unique
!
      IF(n_repl > 0                       ) THEN    ! Need at least one anchor
         ALLOCATE(anchor_num(1:cr_nscat-nscat_old))
         anchor_num(:) = 0
         DO k=1, cr_natoms
            j= cr_iscat(k)-nscat_old
            IF(j>0) anchor_num(j) = anchor_num(j)+1 
         ENDDO
         j=MAXVAL(anchor_num)
         DEALLOCATE(anchor_num)
      IF(n_repl > 2 .AND. n_repl<cr_natoms-2 .AND. j>2) THEN  ! Need at least two anchors for sorting
!
!        Sorting requires separate loops to avoid exchange of anchor atom types
!
         line   = ' '
         length = 1
         CALL conn_do_set (code_res,line, length) ! Reset connectivities
!
!        As RMC_SELECT always uses connectivity no 1 for the selection mode
!        I first create a connectivity of the pair: surface-anchor
!        For MMC a second set of connectivities is created that includes
!        all surface types and all anchors
! 
!        PAIR Correlations for RMC select
         chem_period(:) = .FALSE.
         chem_quick     = .FALSE.
         make_conn: DO k= 1, n_anch
!
!           Prepare the pairwise connectivities, the command will be:
!           '5, 5, Zn, 0.5, 18, deco_0001'  with 5 the anchor type number and Zn the 
!                                           actual anchor type(s)
!
            WRITE(line,'(i4,a1,i4,a1,i4,a22)') &
                  nscat_old+k,     ',', nscat_old+k, ',',   &
                  dc_temp_surf(k), ', 0.5, 18.0, deco_0001'
            WRITE(line,'(i4,a1,i4,a1,i4,a22)') &
                  nscat_old+k,     ',', nscat_old+k, ',',   &
                  anch_id(k,2)   , ', 0.5, 18.0, deco_0001'
            length = 36
            CALL conn_do_set (code_add,line, length)  ! Add connectivity around "AN01"
!
!           Add connectivities around the actual anchors
!
            WRITE(line,'(i4,a1,i4,a1,i4,a22)') &
                  anch_id(k,2)   , ',', nscat_old+k, ',',   &
                  anch_id(k,2)   , ', 0.5, 18.0, deco_0001'
            length = 36
            CALL conn_do_set (code_add,line, length)  ! Add connectivity around surface
         ENDDO make_conn
!
!        All Correlations for MMC energy
!        make_conn_all: DO k= 1, dc_temp_surf(0)
         make_conn_all: DO k= 1, n_anch
            WRITE(line,'(i4,a1)') nscat_old+k,     ','
            DO j=1,n_anch
               i1 = 6 + (j-1)*10
               WRITE(line(i1+0:i1+ 4),'(i4,a1)') nscat_old+j,     ','
               WRITE(line(i1+5:i1+ 9),'(i4,a1)') anch_id(j,2)   , ','
            ENDDO
            i1 = 6 + n_anch         *10
            WRITE(line(i1+0:i1+22),'(a22)'     ) '  0.5, 18.0, deco_0002'
            length = 6 + n_anch         *10+22
            CALL conn_do_set (code_add,line, length)  ! Add connectivity around surface
!
            WRITE(line,'(i4,a1)') anch_id(k,2)   , ','
            DO j=1,n_anch
               i1 = 6 + (j-1)*10
               WRITE(line(i1+0:i1+ 4),'(i4,a1)') nscat_old+j,     ','
               WRITE(line(i1+5:i1+ 9),'(i4,a1)') anch_id(j,2)   , ','
            ENDDO
            i1 = 6 + dc_temp_surf(0)*10
            WRITE(line(i1+0:i1+22),'(a22)'     ) '  0.5, 18.0, deco_0002'
            length = 6 + n_anch         *10+22
            CALL conn_do_set (code_add,line, length)  ! Add connectivity around surface
         ENDDO make_conn_all
!
         CALL create_connectivity(.TRUE.,0,0, ' ')            ! Create actual connecitivity list
!
!        SORT ATOMS WITH MMC REPULSIVE
!
         n_corr = MAX(CHEM_MAX_COR,MMC_MAX_CORR)
         n_scat = MAX(MAXSCAT, MMC_MAX_SCAT)
         CALL alloc_mmc ( n_corr, MC_N_ENERGY, n_scat )       ! Basic mmc allocation
         CALL mmc_init
         ianz = 1
         cpara(1) = 'rese'                                    ! Prepare "set con, reset"
         lpara(1) = 4
         CALL chem_set_con (ianz, cpara, lpara, werte, maxw)  ! Reset conn list
         cpara(1) = 'rese'                                    ! Prepare "set neig, reset"
         CALL chem_set_neig(ianz, cpara, lpara, werte, maxw)  ! Reset neig list
!
         DO k= 1, n_anch
            ianz = 3
            WRITE(cpara(1),'(i4)') k                             ! Prepare "set con, 1, AN01, deco_0002"
            lpara(1) = 4
            WRITE(cpara(2),'(i4)')  nscat_old+k
            lpara(2) = 4
            cpara(3) = 'deco_0002'
            lpara(3) = 9
            CALL chem_set_con (ianz, cpara, lpara, werte, maxw)  ! Set conn list 1 'set con,1, AN01, deco_0002'
         ENDDO
         DO k= 1, n_anch         
            ianz = 3
            WRITE(cpara(1),'(i4)') n_anch          + k           ! Prepare "set con, 2, <sur>,deco_0002"
            lpara(1) = 4
            WRITE(cpara(2),'(i4)')  anch_id(k,2)
            lpara(2) = 4
            cpara(3) = 'deco_0002'
            lpara(3) = 9
            CALL chem_set_con (ianz, cpara, lpara, werte, maxw)  ! Set conn list 2 'set con,2, <sur>, deco_0002'
         ENDDO
         ianz = 1 + n_anch         *2
         cpara(1) = 'con'                                     ! Prepare "set neig, con, 1, 2"
         lpara(1) = 3
         DO k=1, 2*n_anch
            WRITE(cpara(1+k),'(i4)') k
            lpara(1+k) = 4
         ENDDO
         CALL chem_set_neig(ianz, cpara, lpara, werte, maxw)  ! Set neig 'set neigh, con, 1, 2'
         ianz = 2
         cpara(1) = 'swc'                                     ! Prepare "set mode, swchem, 1.0, all"
         lpara(1) = 3
         cpara(2) = 'con'
         lpara(2) = 1
         werte(1) = 1.000
         cpara(3) = ' '
         lpara(3) = 1
         CALL  mmc_set_mode(ianz, cpara, lpara, werte, maxw)  ! Set mode 'set mode,swchem, 1.0, all'
!
         IF(1 > MMC_REP_CORR .or.  1 > CHEM_MAX_COR  .or. &   ! Allocate Repulsive
            MAXSCAT > MMC_REP_SCAT                         ) THEN
            n_corr = MAX(n_corr, CHEM_MAX_COR, MMC_REP_CORR)
            n_scat = MAX(MAXSCAT,MMC_REP_SCAT)
            CALL alloc_mmc_rep (n_corr, n_scat)
!              IF(ier_num /= 0) THEN                             ! Does not seem to fail :-)
!                 RETURN
!              ENDIF
         ENDIF
!                                                                ! define repulsive energy
         DO k= 1, n_anch         
            mmc_allowed(anch_id(k,2)   ) = .true.                ! 'set allowed actual anchor type 
            mmc_allowed(nscat_old+k    ) = .true.                !  and current anchor type "AN01"
            is = nscat_old+k
            DO j= 1, n_anch
               js = nscat_old+j
            CALL mmc_set_disp (1, MC_REPULSIVE, is, js, 100.0, 15.0)
            CALL mmc_set_rep  (1, is, js, 15.,16000., 0.5, 1.)
            ENDDO
         ENDDO
         mmc_cor_energy (1, MC_REPULSIVE) = .true.
         mmc_cor_energy (0, MC_REPULSIVE) = .true.
!
         mo_cyc  =  10*cr_natoms                              ! Define cycles
         mo_feed =   1*cr_natoms                              ! Define feedback
         mo_kt   =   2.5                                      ! Define Temperature
!
         CALL mmc_run_multi(.FALSE. )                          ! Run actual sorting
      ENDIF ! (n_repl > 2 ) THEN                           ! Need at least two anchors for sorting
   line       = 'ignore, all'          ! Ignore all properties
   length     = 11
   CALL property_select(line, length, sav_sel_prop)
   CALL save_internal('internal_sorted')        !     thus this file name is unique
!        Change Anchors back to their original names, keep property flag
!
!        Copy anchor types and positions into temporary array
!
         ALLOCATE(temp_iatom(    1:n_repl))
         ALLOCATE(temp_iscat(    1:n_repl))
         ALLOCATE(temp_pos  (1:3,1:n_repl))
         ALLOCATE(temp_ident(    1:n_repl))
         temp_iatom(:) = 0
         temp_iscat(:) = 0
         temp_pos(:,:) = 0
         temp_ident(:) = 0
         j = 0
         DO ia = 1, cr_natoms
            DO k=1, n_anch
               IF(cr_iscat(ia)==nscat_old + k     ) THEN
                  cr_iscat(ia) = anch_id(k,2)
                  j = j + 1
                  temp_iscat(j) = anch_id(k,2)
                  temp_pos(:,j) = cr_pos(:, ia)
                  temp_ident(j) = anch_id(k,1)
               ENDIF
            ENDDO
         ENDDO
!
         cr_nscat = nscat_old                                 ! Set number of atom types back
         DO k=1, n_anch
            cr_at_lis(nscat_old+k) = ' '                      ! Have atom type disappear
         ENDDO
!
!
CALL rese_cr
CALL save_restore_setting
CALL no_error
CALL readstru_internal( corefile)   ! Read  original core file
IF(cr_natoms == 0 ) THEN
   ier_num = -27
   ier_typ = ER_CHEM
   ier_msg(1) = 'Central core is empty upon read in DECO ?'
   ier_msg(2) = 'Check the property settings. '
   ier_msg(3) = ' '
   DEALLOCATE(anch_id)
   RETURN
ENDIF
!
! Set properties to a well defined setting
!
   line       = 'ignore, all'          ! Ignore all properties
   length     = 11
   CALL property_select(line, length, sav_sel_prop)
   line       = 'present, external'    ! Force atom to be close to a surface
   length     = 17
   CALL property_select(line, length, sav_sel_prop)
   line       = 'absent, outside'      ! Force atom to be inside
   length     = 15
   CALL property_select(line, length, sav_sel_prop)
i=0
loop_anchor: DO j=1,n_repl
   find_atom: DO ia=1,cr_natoms
      IF(ABS(temp_pos(1,j)-cr_pos(1,ia))<EPS .AND. &
         ABS(temp_pos(2,j)-cr_pos(2,ia))<EPS .AND. &
         ABS(temp_pos(3,j)-cr_pos(3,ia))<EPS      ) THEN  ! Found correct position
         cr_iscat(ia) = temp_iscat(j)  ! Set to correct chemistry
         cr_prop (ia) = IBSET (cr_prop (ia), PROP_DECO_ANCHOR)    ! FLAG THIS ATOM AS SURFACE ANCHOR
         temp_iatom(j) = ia                                  ! Store atom number for this anchor
         i=i+1
      ENDIF
   ENDDO find_atom
ENDDO loop_anchor
!
IF(MAXVAL(temp_iatom)==0) THEN
   CALL save_restore_setting
   CALL no_error
   CALL readstru_internal( corefile)   ! Read  original core file
   ier_num = -131
   ier_typ = ER_APPL
   ier_msg(1) = 'Is the surface very small, just a few atoms?'
   ier_msg(2) = 'Is the coverage too small? '
   ier_msg(3) = 'Check the set ligand command A2'
   DEALLOCATE(anch_id)
   RETURN
ENDIF
!
nanch = UBOUND(dcc_surf,2)
ALLOCATE(anchor(1:2,0:UBOUND(dcc_surf,2)))
anchor(:,:) = 0
!
         istart = 1
         iend   = cr_natoms
         natoms_prior = cr_natoms
         ier_num = 0
         main_loop: DO ia=1, n_repl                  ! Loop over the anchors
              is     = temp_iscat(ia)                ! get scattering type central
              xyz(:) = temp_pos(:,ia)                ! Get atom position
!
!             Go through all definitions
              jdef = 0
      defs: DO dc_temp_id=1, dcc_num
                 is_sel: IF(dc_temp_id==temp_ident(ia)) THEN
                 jdef = jdef + 1
        dc_temp_dist = dcc_dist(1,dc_temp_id)
        j = dcc_surf(1,0,dc_temp_id)
        dc_temp_surf(0:j) = dcc_surf(1,0:j,dc_temp_id)
        temp_secnd        = dcc_secnd(dc_temp_id)
                    l_correct = .FALSE.
                    DO i=1,dc_temp_surf(0)
                       l_correct = l_correct .OR. temp_iscat(ia) == dc_temp_surf(i)
                    ENDDO
                    IF(l_correct) THEN                               ! Matching atom in current definition
           dc_temp_type = dcc_type(dc_temp_id)
           dc_temp_axis(0:2) = dcc_axis(0:2, dc_temp_id)
           dc_temp_neig      = dcc_neig(1,dc_temp_id)
           dc_temp_natoms    = dcc_natoms(dc_temp_id)
              temp_lrestrict = dcc_lrestrict(dc_temp_id)
           dc_temp_surfnew(:) = 0
           j = dcc_surfnew(0, dc_temp_id)
           dc_temp_surfnew(1:j) = dcc_surfnew(1:j, dc_temp_id)
           WRITE(mole_name,1000) dc_temp_id, dcc_file(dc_temp_id)(1:dcc_lfile(dc_temp_id))
1000 FORMAT('internal_',I4.4,'_',a)
           mole_length = 14+dcc_lfile(dc_temp_id)
           ncon = dcc_surf(0,0, dc_temp_id)
           anchor(:,:) =dcc_surf(1:2,:,dc_temp_id)
           temp_angle = dcc_angle(dc_temp_id)
                       SELECT CASE(dc_temp_type)
                          CASE ( DC_NORMAL )                 ! Molecule in normal position
                             IF(ncon == 1) THEN
                                CALL deco_place_normal(dc_temp_id, temp_iatom(ia), &
                                  dc_temp_axis, dc_temp_surfnew,             &
                                  m_type_old, mole_name, &
                                  dc_temp_natoms, &
                                  UBOUND(dcc_atom_name,1), &
                                  dcc_atom_name(:,dc_temp_id),    &
                                  dcc_adp      (:,dc_temp_id),    &
                                  dcc_biso     (  dc_temp_id),    &
                                  dcc_clin     (  dc_temp_id),    &
                                  dcc_cqua     (  dc_temp_id),    &
                                  dc_temp_neig, dc_temp_dist,      &
                                  istart, iend, temp_lrestrict,    &
                                  dcc_hkl(1,0,dc_temp_id), &
                                  dcc_hkl(1:3,1:dcc_hkl(1,0,dc_temp_id),dc_temp_id), &
                                  dcc_tilt(dc_temp_id), dcc_tilt_hkl(1:3,dc_temp_id),&
                                  dcc_tilt_atom(1:4,dc_temp_id),                     &
                                  dcc_tilt_is_atom(dc_temp_id),                      &
                                  dcc_tilt_is_auto(dc_temp_id)      )
                             ELSE
                                ier_num = -158
                                ier_msg(1) = 'The normal connection requires one bond'
                                EXIT main_loop
                             ENDIF
                          CASE ( DC_BRIDGE )                 ! Molecule in bridge position
                             IF(ncon == 2) THEN
                             CALL deco_place_bridge(dc_temp_id, temp_iatom(ia), &
                                  dc_temp_axis, dc_temp_surfnew,             &
                                  m_type_old, mole_name, &
                                  dc_temp_natoms, &
                                  UBOUND(dcc_atom_name,1), &
                                  dcc_atom_name(:,dc_temp_id),    &
                                  dcc_adp      (:,dc_temp_id),    &
                                  dcc_biso     (  dc_temp_id),    &
                                  dcc_clin     (  dc_temp_id),    &
                                  dcc_cqua     (  dc_temp_id),    &
                                  nanch, anchor,                  &
!                                 dc_temp_neig, dc_temp_dist, &
                                  dcc_neig(1:2,dc_temp_id), dcc_dist(1:2,dc_temp_id), ncon, &
                                  istart, iend, temp_lrestrict,    &
                                  dcc_hkl(1,0,dc_temp_id), &
                                  dcc_hkl(1:3,1:dcc_hkl(1,0,dc_temp_id),dc_temp_id), &
                                  dcc_tilt(dc_temp_id), dcc_tilt_hkl(1:3,dc_temp_id),&
                                  dcc_tilt_atom(1:4,dc_temp_id),                     &
                                  dcc_tilt_is_atom(dc_temp_id) ,                     &
                                  dcc_tilt_is_auto(dc_temp_id)      )
                             ELSE
                                ier_num = -158
                                ier_msg(1) = 'The bridge connection requires two bonds'
                                EXIT main_loop
                             ENDIF
                          CASE ( DC_DOUBLE   )               ! Molecule in double   connection position
                             IF(ncon >  1) THEN
                             CALL deco_place_double(dc_temp_id, temp_iatom(ia), &
                                  dc_temp_axis, dc_temp_surfnew,             &
                                  m_type_old, mole_name, &
                                  dc_temp_natoms, &
                                  UBOUND(dcc_atom_name,1), &
                                  dcc_atom_name(:,dc_temp_id),    &
                                  dcc_adp      (:,dc_temp_id),    &
                                  dcc_biso     (  dc_temp_id),    &
                                  dcc_clin     (  dc_temp_id),    &
                                  dcc_cqua     (  dc_temp_id),    &
                                  nanch, anchor, & ! UBOUND(dcc_surf,2),             &
                                  dcc_neig(1:2,dc_temp_id), dcc_dist(1:2,dc_temp_id), ncon, &
                                  istart, iend, temp_lrestrict,    &
                                  dcc_hkl(1,0,dc_temp_id), &
                                  dcc_hkl(1:3,1:dcc_hkl(1,0,dc_temp_id),dc_temp_id), &
                                  dcc_tilt(dc_temp_id), dcc_tilt_hkl(1:3,dc_temp_id),&
                                  dcc_tilt_atom(1:4,dc_temp_id),                     &
                                  dcc_tilt_is_atom(dc_temp_id) ,                     &
                                  dcc_tilt_is_auto(dc_temp_id)      )
CYCLE main_loop
                             ELSE
                                ier_num = -158
                                ier_msg(1) = 'The double connection requires two bonds'
                                EXIT main_loop
                             ENDIF
                          CASE ( DC_MULTIPLE )               ! Molecule in multiple connection position
                             IF(ncon >  1) THEN
                             CALL deco_place_multi(dc_temp_id, temp_iatom(ia), &
                                  dc_temp_axis, dc_temp_surfnew,             &
                                  m_type_old, mole_name, dc_temp_natoms,    &
                                  UBOUND(dcc_atom_name,1), &
                                  dcc_atom_name(:,dc_temp_id),    &
                                  dcc_adp      (:,dc_temp_id),    &
                                  dcc_biso     (  dc_temp_id),    &
                                  dcc_clin     (  dc_temp_id),    &
                                  dcc_cqua     (  dc_temp_id),    &
                                  nanch, anchor,                  &
                                  dcc_nanch(1:2, dc_temp_id),     &
                                  dcc_neig(1:2,dc_temp_id), dcc_dist(1:2,dc_temp_id), ncon, &
                                  istart, iend, temp_lrestrict,    &
                                  dcc_hkl(1,0,dc_temp_id), &
                                  dcc_hkl(1:3,1:dcc_hkl(1,0,dc_temp_id),dc_temp_id), &
                                  dcc_tilt(dc_temp_id), dcc_tilt_hkl(1:3,dc_temp_id),&
                                  dcc_tilt_atom(1:4,dc_temp_id),                     &
                                  dcc_tilt_is_atom(dc_temp_id) ,                     &
                                  dcc_tilt_is_auto(dc_temp_id)      )
CYCLE main_loop
                             ELSE
                                ier_num = -158
                                ier_msg(1) = 'The multiple connection requires > one bond'
                                EXIT main_loop
                             ENDIF
                          CASE ( DC_ACCEPTOR )                 ! Molecule in normal position
                             IF(ncon == 1) THEN
                                CALL deco_place_acceptor(dc_temp_id, temp_iatom(ia), &
                                  dc_temp_axis, dc_temp_surfnew,                &
                                  m_type_old, mole_name, dc_temp_natoms, &
                                  UBOUND(dcc_atom_name,1), &
                                  dcc_atom_name(:,dc_temp_id),    &
                                  dcc_adp      (:,dc_temp_id),    &
                                  dcc_biso     (  dc_temp_id),    &
                                  dcc_clin     (  dc_temp_id),    &
                                  dcc_cqua     (  dc_temp_id),    &
                                  dc_temp_neig, dc_temp_dist, temp_secnd, &
                                  istart, iend, temp_lrestrict,    &
                                  dcc_hkl(1,0,dc_temp_id), &
                                  dcc_hkl(1:3,1:dcc_hkl(1,0,dc_temp_id),dc_temp_id), &
                                  temp_angle)
                             ELSE
                                ier_num = -158
                                ier_msg(1) = 'The acceptor connection requires one bond'
                                EXIT main_loop
                             ENDIF
                          CASE ( DC_DONOR )                 ! Molecule in normal position
                             IF(ncon == 1) THEN
                                CALL deco_place_donor(dc_temp_id, temp_iatom(ia), &
                                  dc_temp_axis, dc_temp_surfnew,                &
                                  m_type_old, mole_name, dc_temp_natoms, &
                                  UBOUND(dcc_atom_name,1), &
                                  dcc_atom_name(:,dc_temp_id),    &
                                  dcc_adp      (:,dc_temp_id),    &
                                  dcc_biso     (  dc_temp_id),    &
                                  dcc_clin     (  dc_temp_id),    &
                                  dcc_cqua     (  dc_temp_id),    &
                                  dc_temp_neig, dc_temp_dist,      &
                                  istart, iend, temp_lrestrict,    &
                                  dcc_hkl(1,0,dc_temp_id), &
                                  dcc_hkl(1:3,1:dcc_hkl(1,0,dc_temp_id),dc_temp_id), &
                                  temp_angle)
                             ELSE
                                ier_num = -158
                                ier_msg(1) = 'The donor connection requires one bond'
                                EXIT main_loop
                             ENDIF
                          CASE ( DC_CHELATE  )               ! Molecule in chelate  connection position
                             IF(ncon == 2) THEN
                             CALL deco_place_chelate(dc_temp_id, temp_iatom(ia), &
                                  dc_temp_axis, dc_temp_surfnew,             &
                                  m_type_old, mole_name, &
                                  dc_temp_natoms, &
                                  UBOUND(dcc_atom_name,1), &
                                  dcc_atom_name(:,dc_temp_id),    &
                                  dcc_adp      (:,dc_temp_id),    &
                                  dcc_biso     (  dc_temp_id),    &
                                  dcc_clin     (  dc_temp_id),    &
                                  dcc_cqua     (  dc_temp_id),    &
                                  nanch, anchor, & ! UBOUND(dcc_surf,2),             &
                                  dcc_neig(1:2,dc_temp_id), dcc_dist(1:2,dc_temp_id), ncon, &
                                  istart, iend, temp_lrestrict,    &
                                  dcc_hkl(1,0,dc_temp_id), &
                                  dcc_hkl(1:3,1:dcc_hkl(1,0,dc_temp_id),dc_temp_id), &
                                  dcc_tilt(dc_temp_id), dcc_tilt_hkl(1:3,dc_temp_id),&
                                  dcc_tilt_atom(1:4,dc_temp_id),                     &
                                  dcc_tilt_is_atom(dc_temp_id) ,                     &
                                  dcc_tilt_is_auto(dc_temp_id)      )
                                  IF(ier_num/=0) EXIT main_loop
                                  CYCLE main_loop
                             ELSE
                                ier_num = -158
                                ier_msg(1) = 'The chelate connection requires two bond'
                                EXIT main_loop
                             ENDIF
                       END SELECT
                    ENDIF
                    CYCLE main_loop
               ENDIF is_sel ! END IF BLOCK is selected
            ENDDO defs
         ENDDO main_loop   ! END DO main loop over all atoms
!
         IF(ier_num==-158) THEN
           CALL save_restore_setting
           ier_num = 0
           ier_typ = 0
           CALL readstru_internal(corefile)   ! Read core file
           ier_num = -158
           ier_typ = ER_APPL
         ELSEIF(cr_natoms == natoms_prior) THEN    ! No atoms added Failure
           CALL save_restore_setting
           CALL no_error
           CALL readstru_internal(corefile)   ! Read core file
           ier_num = -131
           ier_typ = ER_APPL
           ier_msg(1) = 'Is the surface very small, just a few atoms?'
           ier_msg(2) = 'Is the coverage too small? '
           ier_msg(3) = 'Check the set ligand command A3'
         ENDIF
!
         IF(ier_num == 0)  THEN       ! Success in main_loop
!
            DO ia=1, n_repl
               cr_prop (temp_iatom(ia)) = IBCLR (cr_prop (temp_iatom(ia)), PROP_DECO_ANCHOR)  ! UNFLAG THIS ATOM AS SURFACE ANCHOR
            ENDDO
            CALL do_purge
!        ELSE     ! Error in main_loop
!          CALL save_restore_setting
!          CALL no_error
!          CALL readstru_internal(corefile)   ! Read core file
!          ier_num = -131
!          ier_typ = ER_APPL
!          ier_msg(1) = 'Is the surface very small, just a few atoms?'
!          ier_msg(2) = 'Is the coverage too small? '
!          ier_msg(3) = 'Check the set ligand command'
         ENDIF
      ELSE     ! n_repl > 0   !! No anchor atoms found
        CALL save_restore_setting
        CALL no_error
        CALL readstru_internal( corefile)   ! Read  core file
        ier_num = -131
        ier_typ = ER_APPL
        ier_msg(1) = 'Is the surface very small, just a few atoms?'
        ier_msg(2) = 'Is the coverage too small? '
        ier_msg(3) = 'Check the set ligand command A4'
      ENDIF    ! n_repl > 0   !! No anchor atoms found
   ELSE     ! SHELL has atoms
     CALL rese_cr
     CALL save_restore_setting
     CALL no_error
     CALL readstru_internal( corefile)   ! Read  core file
     ier_num = -130
     ier_typ = ER_APPL
     ier_msg(1) = 'Possible reasons: no boundary was used to cut'
     ier_msg(2) = 'Distance to external surface is too large'
     ier_msg(3) = 'Check settings and command in surface menu'
   ENDIF
!
   IF(ier_num /=0) THEN
      ier_num_deco = ier_num
      ier_typ_deco = ier_typ
      ier_msg_deco = ier_msg
   ELSE
      ier_num_deco =  0
      ier_typ_deco =  0
      ier_msg_deco = ' '
   ENDIF
!  IF(ier_num == 0 ) THEN
      CALL save_restore_setting     ! Restore user save settigns
!  ENDIF
!
!  Clean up temporary arrays
!
   IF(ALLOCATED(temp_prob )) DEALLOCATE(temp_prob, STAT=istatus)
   IF(ALLOCATED(anch_id   )) DEALLOCATE(anch_id)
   IF(ALLOCATED(temp_ident)) DEALLOCATE(temp_ident)
   IF(ALLOCATED(temp_iscat)) DEALLOCATE(temp_iscat)
   IF(ALLOCATED(temp_iatom)) DEALLOCATE(temp_iatom)
   IF(ALLOCATED(temp_pos  )) DEALLOCATE(temp_pos  )
IF(ALLOCATED(anchor    )) DEALLOCATE(anchor)
   line   = ' '
   length = 1
   CALL conn_do_set (code_res,line, length) ! Reset connectivities
!
! Clean up internal files
!
CALL store_remove_single(corefile, ier_num)
CALL store_remove_single(shellfile, ier_num)
CALL store_remove_single('internal_anchors', ier_num)
CALL store_remove_single('internal_sorted', ier_num)
rdefs: DO dc_temp_id=1, dcc_num
   WRITE(mole_name,1000) dc_temp_id, dcc_file(dc_temp_id)(1:dcc_lfile(dc_temp_id))
   CALL store_remove_single(mole_name, ier_num)
ENDDO rdefs
IF(ier_num/=0) THEN
   ier_typ = ER_APPL
   ier_msg(1) = 'Could not remove temporary internal storage'
   ier_msg(2) = 'in deco_run '
   ier_msg(3) = 'Please document and report'
ENDIF

!
!  Restore DECO ERROR SETTINGS 
!
   IF(ier_num_deco/=0) THEN
      ier_num = ier_num_deco
      ier_typ = ier_typ_deco
      ier_msg = ier_msg_deco
   ELSE
      chem_purge = .TRUE.     ! The crystal is most likely NOT periodic.
                              ! In the rare circumstances that is is the user
                              ! has to turn this on explicitly
      chem_quick = .FALSE.
      chem_period(:) = .FALSE.
   ENDIF
!
   END SUBROUTINE deco_run
!
!######################################################################
!
SUBROUTINE deco_get_molecules(host_a0, host_win, host_tran_fi)
!
!  Loops over all the environments that are defined. 
!  For each the corresponding molecule is loaded and the first
!  "bond" connectivity is read out to ensure that the "first neighbor"
!  is set properly. it shifts the molecule so that the first 
!  neighbor is at (0,0,0)
!  The molecules are stored in the internal storage
!
USE crystal_mod
USE deco_mod
USE discus_allocate_appl_mod
USE metric_mod
!USE modify_mod
USE molecule_mod
USE prop_para_func
USE spcgr_apply
USE save_menu, ONLY: save_internal
USE discus_save_mod
USE structur, ONLY: do_readstru, rese_cr
USE trafo_mod
!
use read_internal_mod
!
IMPLICIT none
!
REAL, DIMENSION(1:3), INTENT(IN) :: host_a0
REAL, DIMENSION(1:3), INTENT(IN) :: host_win
REAL, DIMENSION(4,4), INTENT(IN) :: host_tran_fi
CHARACTER (LEN=1024)                :: strufile
   CHARACTER (LEN=1024)                :: mole_name
CHARACTER (LEN=1024)                :: savefile
CHARACTER (LEN=1024)                :: line
!
   INTEGER  :: mole_length      ! Molecule name length
   INTEGER  :: i,j              ! dummy index
   INTEGER  :: length           ! dummy index
!  INTEGER  :: istatus          ! status
!  INTEGER  :: natoms           ! number of atoms in file
!  INTEGER  :: ntypes           ! number of atom types in file
!  INTEGER  :: n_mole           ! number of molecules
!  INTEGER  :: n_type           ! number of molecule types
!  INTEGER  :: n_atom           ! number of atoms in molecules
!  INTEGER  :: init   = -1      ! Initialize atom names/types
!  INTEGER  :: itype            ! atom scattering type
!  INTEGER, DIMENSION(0:3)  :: isurface            ! atom surface type
!  REAL, DIMENSION(3) :: posit  ! atom position
   REAL, DIMENSION(4) :: posit4 ! atom position
   REAL, DIMENSION(4) :: uvw4   ! atom position
!  REAL, DIMENSION(4,4) :: rd_tran_f  ! transformation matrix to cartesian
!  INTEGER  :: secnd            ! Index of second neighbor
!  INTEGER  :: iprop            ! atom property 
!  LOGICAL  :: lcell  = .true.  ! Treat atoms with equal name and B as one type
   LOGICAL  :: success= .true.  ! Found matching molecule file name
   INTEGER, DIMENSION(1:2) :: temp_axis
!  INTEGER, DIMENSION(0:4) :: temp_surf
   INTEGER                 :: temp_neig
   REAL                    :: temp_dist
!
INTEGER               :: n_mscat     ! temporary number of molecule scattering types
INTEGER               :: lll         ! Dummy index
LOGICAL, PARAMETER    :: lout = .FALSE.
REAL                  :: dist
REAL                  :: shortest
REAL                  :: longest
REAL, DIMENSION(3)    :: uvw_out
REAL, DIMENSION(3)    :: w
!
   strufile = ' '
!
main: DO i=1, dcc_num
!
   mole_name   = dcc_file(i)
   mole_length = dcc_lfile(i)
   temp_axis(1) = dcc_axis(1,i)
   temp_axis(2) = dcc_axis(2,i)
   temp_neig    = dcc_neig(1,i)
   temp_dist    = dcc_dist(1,i)
   strufile     = mole_name
   success = .FALSE.
   CALL rese_cr
   CALL do_readstru(strufile)
   IF(ier_num /= 0) RETURN
   CALL setup_lattice (cr_a0, cr_ar, cr_eps, cr_gten, cr_reps,    &
         cr_rten, cr_win, cr_wrez, cr_v, cr_vr, lout, cr_gmat, cr_fmat, &
         cr_cartesian,                                                  &
              cr_tran_g, cr_tran_gi, cr_tran_f, cr_tran_fi)
   IF(ier_num /= 0) RETURN
   IF(cr_natoms < temp_axis(1) .OR. cr_natoms < temp_axis(2)) THEN
      ier_num = -144
      ier_typ = ER_APPL
      ier_msg(1) = ' The atom numbers listed on the set axis'
      ier_msg(2) = ' command are outside the number of atoms'
      WRITE(ier_msg(3),'(a,a)') ' for molecule ', &
         strufile(1:MIN(LEN(ier_msg),LEN_TRIM(strufile)))
      RETURN
   ENDIF
   IF(cr_natoms < temp_neig   ) THEN
      ier_num = -144
      ier_typ = ER_APPL
      ier_msg(1) = ' The atom numbers listed on the set bond'
      ier_msg(2) = ' command are outside the number of atoms'
      WRITE(ier_msg(3),'(a,a)') ' for molecule ', &
              strufile(1:MIN(LEN(ier_msg),LEN_TRIM(strufile)))
      RETURN
   ENDIF
!
   dcc_natoms(i) = cr_natoms         ! Store number of atoms in this molecule
   IF(cr_nscat > DCC_MAXMSCAT .OR. cr_nscat > UBOUND(dcc_atom_name,1)) THEN
      n_mscat = cr_nscat + 5
      CALL alloc_deco(MAXSCAT, DCC_MAXNUM, DCC_MAXANCH, DCC_MAXHKL, DCC_MAXNEW, n_mscat)
      IF(ier_num /= 0) RETURN
   ENDIF
   dcc_atom_name(0:cr_nscat,i) = cr_at_lis(0:cr_nscat)
   dcc_adp      (0:cr_nscat,i) = cr_dw    (0:cr_nscat)
   dcc_mole_type(i) = mole_type(1)
   dcc_biso     (i) = mole_biso(mole_type(1))    ! Here we use dcc_biso for each molecule
   dcc_clin     (i) = mole_clin(mole_type(1))    ! Here we use dcc_clin for each molecule
   dcc_cqua     (i) = mole_cqua(mole_type(1))    ! Here we use dcc_cqua for each molecule
!
   uvw_out(:) = cr_pos(:, temp_neig) ! save position of molecule atom bonded to surface
   DO j=1, cr_natoms ! Shift the first neighbor to 0,0,0
      cr_pos(1,j) = cr_pos(1,j) - uvw_out(1)
      cr_pos(2,j) = cr_pos(2,j) - uvw_out(2)
      cr_pos(3,j) = cr_pos(3,j) - uvw_out(3)
   ENDDO
!
!  Determine closest neighbor to first atom
   dcc_secnd(i) = 0
   shortest     = 1.0E8
   longest      = 0.0
   lll          = 0
   DO j=1, cr_natoms
      w(:) = cr_pos(:,j)
      dist = skalpro (w, w, cr_gten) 
      IF(dist>0.0 .AND. dist<shortest) THEN  ! Found new closest neighbor
         dcc_secnd(i) = j                           ! store as dcc_secnd
         shortest = dist
      ENDIF
!
!     FIND longest distance for axis = AUTO mode
!
      IF(j/=dcc_neig(1,i) .AND. j/=dcc_neig(2,i)) THEN  ! Ignore neigbor atoms
         IF(dist>0.0 .AND. dist>longest ) THEN  ! Found new longest neighbor
            lll = j                             ! store as lll
            longest = dist
         ENDIF
      ENDIF
   ENDDO
   IF(dcc_axis(0,i) == -1) THEN ! We have axis == auto mode
      IF(lll > 1) THEN          ! A proper second neighbor exists 
         dcc_axis(1,i) = 1
         dcc_axis(2,i) = lll
         dcc_axis(0,i) = 2      ! set to normal mode
      ELSE
         dcc_axis(:,i) = 0      ! Flag to axis = none mode
      ENDIF
   ENDIF
!
! Transform molecule geometry into host geometry
   DO j=1, cr_natoms
      posit4(1:3) = cr_pos(1:3,j)
      posit4(4)   = 1.0
      CALL trans(posit4,cr_tran_f ,uvw4, 4)       ! Transform mol to cartesian
      CALL trans(uvw4  ,host_tran_fi,posit4, 4)   ! Transform cartesian to host_crystal
      cr_pos(1:3,j) = posit4(1:3)
   ENDDO
!
   cr_a0(:)  = host_a0
   cr_win(:) = host_win
      CALL setup_lattice (cr_a0, cr_ar, cr_eps, cr_gten, cr_reps,    &
      cr_rten, cr_win, cr_wrez, cr_v, cr_vr, lout, cr_gmat, cr_fmat, &
      cr_cartesian,                                                  &
           cr_tran_g, cr_tran_gi, cr_tran_f, cr_tran_fi)
!
! save as file internal_nnnn_originalname
   line       = 'ignore, all'          ! Ignore all properties
   length     = 11
   CALL property_select(line, length, sav_sel_prop)
   WRITE(savefile,1000) i,dcc_file(i)(1:dcc_lfile(i))
!                                                                       
   CALL save_internal(savefile)        !     thus this file name is unique
!
ENDDO main
1000 FORMAT('internal_',i4.4,'_',a)
!
!
   END SUBROUTINE deco_get_molecules
!
!######################################################################
!
   SUBROUTINE deco_show
!
   USE atom_name
   USE deco_mod
   USE prompt_mod
   IMPLICIT none
!
   CHARACTER(LEN=1024) :: line
   INTEGER :: i, j, k
!
!
   DO i=1, dcc_num
      WRITE(output_io,*)
      WRITE(output_io,1100) i, dcc_name(i)(1:dcc_lname(i))
      WRITE(output_io,1200) dcc_file(i)(1:dcc_lfile(i))
      IF(dcc_axis(0,i) == 0) THEN
         WRITE(output_io,1300) dcc_type(i),dcc_ctype(dcc_type(i))
      ELSEIF(dcc_axis(0,i) == -1) THEN
         WRITE(output_io,1305) dcc_type(i),dcc_ctype(dcc_type(i))
      ELSEIF(dcc_axis(2,i) == -1) THEN
         WRITE(output_io,1310) dcc_type(i),dcc_ctype(dcc_type(i)), dcc_axis(1,i)
      ELSE
         WRITE(output_io,1320) dcc_type(i),dcc_ctype(dcc_type(i)), dcc_axis(1,i),dcc_axis(2,i)
      ENDIF
      WRITE(output_io,1323) dcc_type(i),dcc_ctype(dcc_type(i)), dcc_tilt(i)
      IF(dcc_tilt_is_auto(i)) THEN
         WRITE(output_io,1325) dcc_type(i),dcc_ctype(dcc_type(i))
      ELSE
         IF(dcc_tilt_is_atom(i)) THEN
            WRITE(output_io,1322) dcc_type(i),dcc_ctype(dcc_type(i)), dcc_tilt_atom(1:4,i)
         ELSE
            WRITE(output_io,1321) dcc_type(i),dcc_ctype(dcc_type(i)), dcc_tilt_hkl(1:3,i)
         ENDIF
      ENDIF

      IF(.NOT.dcc_lrestrict(i)) THEN
            WRITE(output_io, 1330)
      ELSE
         IF(dcc_lform(i)) THEN
            WRITE(output_io, 1340)
         ELSE
            WRITE(output_io, 1345)
         ENDIF
         DO j=1, dcc_hkl(1,0,i)
            WRITE(output_io, 1350) dcc_hkl(:,j,i)
         ENDDO
      ENDIF
      WRITE(output_io, 1400 ) dcc_dens(i)
      WRITE(output_io, 1500 ) dcc_surf(0,0,i)
      DO k=1, dcc_surf(0,0,i)
         line = ' '
         DO j=1,dcc_surf(k,0,i)
            line((j-1)*9+1:(j-1)*9+9) = at_name (dcc_surf(k,j,i)) // ' '
         ENDDO
         WRITE(output_io, 1600) line(1: dcc_surf(k,0,i)*10), &
            dcc_neig(k,i), dcc_dist(k,i)
      ENDDO
   ENDDO
!
1100 FORMAT(' Definition        :',i4,' ',4a)
1200 FORMAT('   Molecule file   :     ',a)
1300 FORMAT('   Connection type :     ',i4,' ',a8, ' Ligand axis : NONE')
1305 FORMAT('   Connection type :     ',i4,' ',a8, ' Ligand axis : AUTO')
1310 FORMAT('   Connection type :     ',i4,' ',a8, ' Ligand axis : ',i4,' to last')
1320 FORMAT('   Connection type :     ',i4,' ',a8, ' Ligand axis : ',i4,' to ',i4 )
1321 FORMAT('   Connection type :     ',i4,' ',a8, ' Tilt plane  : ',3(f8.3,2x)  )
1322 FORMAT('   Connection type :     ',i4,' ',a8, ' Tilt atoms  : ',4(i4  ,2x)  )
1323 FORMAT('   Connection type :     ',i4,' ',a8, ' Ligand tilt : ',F8.3        )
1325 FORMAT('   Connection type :     ',i4,' ',a8, ' Tilt plane  : AUTO'         )
1330 FORMAT('   Connection open to all surfaces')
1340 FORMAT('   Connection restricted to following forms')
1345 FORMAT('   Connection restricted to following surfaces')
1350 FORMAT('              H K L:     ',2(i4,2x),i4)
1400 FORMAT('   Surface density :     ',f6.3,' Ligands/A^2')
1500 FORMAT('   No. of Connect. :     ',i4)
1600 FORMAT('   Type , lig, dist: ',a,1x,i4,1x,f8.5)
!
   END SUBROUTINE deco_show
!
!######################################################################
!
SUBROUTINE deco_reset
!
!  Set all definitions back to system default
!
USE deco_mod
use discus_allocate_appl_mod
!
IMPLICIT none
!
!  INTEGER :: istatus
CALL alloc_deco(MAXSCAT, 4, 3, 3, 2, 3)
!
dcc_num = 0
IF(ALLOCATED(dcc_name))         dcc_name      (    :) = ' '
IF(ALLOCATED(dcc_lname))        dcc_lname     (    :) = 0
IF(ALLOCATED(dcc_file))         dcc_file      (    :) = ' '
IF(ALLOCATED(dcc_lfile))        dcc_lfile     (    :) = 0
IF(ALLOCATED(dcc_natoms))       dcc_natoms    (    :) = 0
IF(ALLOCATED(dcc_atom_name))    dcc_atom_name (:,  :) = ' '
IF(ALLOCATED(dcc_adp))          dcc_adp       (:,  :) = 0
IF(ALLOCATED(dcc_biso))         dcc_biso      (    :) = 0
IF(ALLOCATED(dcc_clin))         dcc_clin      (    :) = 0
IF(ALLOCATED(dcc_cqua))         dcc_cqua      (    :) = 0
IF(ALLOCATED(dcc_mole_type))    dcc_mole_type (    :) = 0
IF(ALLOCATED(dcc_type))         dcc_type      (    :) = 0
IF(ALLOCATED(dcc_surf))         dcc_surf      (:,:,:) = 0
IF(ALLOCATED(dcc_neig))         dcc_neig      (:,  :) = 0
IF(ALLOCATED(dcc_secnd))        dcc_secnd     (    :) = 0
IF(ALLOCATED(dcc_axis))         dcc_axis      (:,  :) = 0
IF(ALLOCATED(dcc_axis))         dcc_axis      (0,  :) = -1
IF(ALLOCATED(dcc_lrestrict))    dcc_lrestrict (    :) = .FALSE.
IF(ALLOCATED(dcc_lform))        dcc_lform     (    :) = .FALSE.
IF(ALLOCATED(dcc_hkl))          dcc_hkl       (:,:,:) = 0
IF(ALLOCATED(dcc_surfnew))      dcc_surfnew   (:,  :) = 0
IF(ALLOCATED(dcc_dens))         dcc_dens      (    :) = 0.0
IF(ALLOCATED(dcc_dist))         dcc_dist      (:,  :) = 0.0
IF(ALLOCATED(dcc_angle))        dcc_angle     (    :) = 170.0  !Bond angle for Hydrogen bonds
IF(ALLOCATED(dcc_tilt))         dcc_tilt(:)         = 0.0            ! Tilt angle for ligand off axis
IF(ALLOCATED(dcc_tilt_hkl))     dcc_tilt_hkl(:,:)   = 0.0        ! Normal to molecule plane
IF(ALLOCATED(dcc_tilt_atom))    dcc_tilt_atom(:,:)  = 0       ! Atoms that form molecule plane
IF(ALLOCATED(dcc_tilt_is_atom)) dcc_tilt_is_atom(:) = .TRUE.    ! Tilt angle for ligand off axis
IF(ALLOCATED(dcc_tilt_is_auto)) dcc_tilt_is_auto(:) = .TRUE.    ! Tilt angle for ligand off axis
!
dc_init        = .true.    ! We need to initialize
dc_latom(:)    = .false.
!
END SUBROUTINE deco_reset
!
!******************************************************************************
!
   SUBROUTINE deco_place_normal(temp_id, ia, &
                            mole_axis, mole_surfnew,           &
                            m_type_old, mole_name, &
                            mole_natoms, &
                            mole_nscat, mole_atom_name, &
                            mole_dw, r_m_biso, r_m_clin, r_m_cqua,         &
                            neig, dist, istart, iend, &
                            lrestrict, nhkl, rhkl   , &
                            tilt, tilt_hkl, tilt_atom, tilt_is_atom, &
                            tilt_is_auto)
!
USE crystal_mod
   USE chem_mod
!  USE deco_mod
   USE metric_mod
   USE modify_mod
   USE molecule_mod
   USE molecule_func_mod
   USE point_grp
   USE prop_para_mod
   USE read_internal_mod
   USE surface_func_mod
   USE surface_mod
   USE symm_menu
   USE symm_mod
   USE symm_sup_mod
   USE trafo_mod
!
   USE param_mod
!
   IMPLICIT NONE
!
   INTEGER,                 INTENT(IN) ::    temp_id      ! The definition to be used
   INTEGER,                 INTENT(IN) :: ia              ! Surface atom number
   INTEGER,                 INTENT(IN) :: m_type_old      ! molecule types previous to all placements
   INTEGER, DIMENSION(0:2), INTENT(IN) :: mole_axis       ! Atoms that define molecule axis
   INTEGER, DIMENSION(1:20),INTENT(IN) :: mole_surfnew    ! Atoms that will be flagged as surface atoms
   CHARACTER (LEN=1024),    INTENT(IN) :: mole_name       ! Molecule file name
   INTEGER,                 INTENT(IN) :: mole_natoms     ! Number of atoms in molecule
   INTEGER,                 INTENT(IN) :: mole_nscat      ! Number of atoms in molecule
   CHARACTER (LEN=4   ), DIMENSION(0:mole_nscat),   INTENT(IN) :: mole_atom_name ! Atom names in the molecule
   REAL                , DIMENSION(0:mole_nscat),   INTENT(IN) :: mole_dw        ! ADPs       in the molecule
REAL,                    INTENT(IN) :: r_m_biso        ! Molecular Biso
REAL,                    INTENT(IN) :: r_m_clin        ! Molecular linear correction
REAL,                    INTENT(IN) :: r_m_cqua        ! Molecular quadratic correction
   INTEGER,                 INTENT(IN) :: neig            ! Connected to this neighbor in mole
   REAL   ,                 INTENT(IN) :: dist            ! distance to ligand molecule
   INTEGER,                 INTENT(IN) :: istart          ! First atom for surface determination
   INTEGER,                 INTENT(IN) :: iend            ! Last  atom for surface determination
   LOGICAL,                 INTENT(IN) :: lrestrict       ! Restriction to a surface type T/F
   INTEGER,                 INTENT(IN) :: nhkl            ! Number of faces for the restriction
INTEGER, DIMENSION(3,nhkl), INTENT(IN) :: rhkl            ! actual faces for the restriction
REAL                      , INTENT(IN) :: tilt            ! Molecule tilt angle
REAL,    DIMENSION(3)     , INTENT(IN) :: tilt_hkl        ! Molecule tilt plane by this normal
INTEGER, DIMENSION(4)     , INTENT(IN) :: tilt_atom       ! Molecule tilt plane defined by these atoms
LOGICAL                   , INTENT(IN) :: tilt_is_atom    ! Plane defined by atoms
LOGICAL                   , INTENT(IN) :: tilt_is_auto    ! Plane defined by atoms
!
   REAL   , PARAMETER      :: EPS = 1.0E-6
!
   CHARACTER (LEN=1024)    :: line
   INTEGER                 ::    j, im, laenge ! Dummy variables
   INTEGER                 :: nold            ! atom number previous to current molecule
   INTEGER                 :: m_type_new      ! new molecule types 
   INTEGER                 :: surf_char       ! Surface character, plane, edge, corner, ...
   INTEGER, DIMENSION(3,6) :: surface_normal  ! Set of local normals (:,1) is main normal
   INTEGER, DIMENSION(3)   :: surf_kante      ! Edge vector if not a plane
   INTEGER, DIMENSION(6)   :: surf_weight     ! Best normal has heighest weight
   INTEGER                              :: test_nhkl
   INTEGER, DIMENSION(:,:), ALLOCATABLE :: test_hkl   ! Temporary array, needed as point_test expects an allocatable
!
   REAL   , DIMENSION(1:3) :: surf_normal     ! Normal to work with
   REAL   , DIMENSION(3)   :: vnull           ! null vector
   REAL   , DIMENSION(3)   :: origin          ! Symmetry origin
   REAL                    :: normal_l        ! local normal

   REAL   , DIMENSION(3)   :: posit           ! Atom coordinates for read internal
!  REAL                    :: dw1             ! Atom ADP for read internal
   INTEGER                 :: itype           ! Atom type for read internal
   INTEGER                 :: iprop           ! Atom property for read internal
   INTEGER, DIMENSION(0:3) :: isurface        ! Atom surface for read internal
   REAL   , DIMENSION(1:3) :: axis_ligand                 ! Initial molecule orientation
!  INTEGER                 :: i_m_type
   INTEGER, DIMENSION(4)   :: hkl
   INTEGER                 :: in_mole,in_moleatom
   INTEGER                 :: n_atoms_orig
!
   vnull(:) = 0.00
!
!  Determine surface character, if growth is restricted check if we're at proper surface
!
   hkl(4) = 0
   CALL surface_character(ia, istart, iend, surf_char, surface_normal, surf_kante, surf_weight, .TRUE.)
IF(surf_char == 0) RETURN
   surf_normal(1:3) = FLOAT(surface_normal(:,1))
   hkl(1:3)         =       surface_normal(:,1)
   IF(lrestrict) THEN
      test_nhkl = nhkl
      ALLOCATE(test_hkl(3,test_nhkl))
      test_hkl(1:3,1:test_nhkl) =    rhkl(1:3,1:test_nhkl)
      IF(.NOT.point_test(hkl, test_hkl, test_nhkl, .TRUE.) ) THEN
         RETURN
      ENDIF
      DEALLOCATE(test_hkl)
   ENDIF
!
   n_atoms_orig = cr_natoms
!
!  Accept as surfaces any outside facing surface and all internal "flat" surfaces
   IF(surf_char /=0 .AND. surf_char > -SURF_EDGE) THEN    ! Surface atoms only
!
      nold = cr_natoms
      IF(mole_axis(0)==2) THEN
      im = mole_axis(2)                                ! Make mole axis from spcified atoms
      CALL struc_read_one_atom_internal(mole_name, im, posit, itype, iprop, isurface,in_mole,in_moleatom)
      axis_ligand(1) = posit(1)
      axis_ligand(2) = posit(2)
      axis_ligand(3) = posit(3)
      im = mole_axis(1)
!
      CALL struc_read_one_atom_internal(mole_name, im, posit, itype, iprop, isurface,in_mole,in_moleatom)
      axis_ligand(1) = axis_ligand(1) - posit(1)       ! this is mole axis
      axis_ligand(2) = axis_ligand(2) - posit(2)
      axis_ligand(3) = axis_ligand(3) - posit(3)
      ENDIF
!     Insert molecule atoms into crystal at correct origin in initial orientation
      normal_l = sqrt (skalpro (surf_normal, surf_normal, cr_gten))
      CALL struc_read_one_atom_internal(mole_name, neig, posit, itype, iprop, isurface,in_mole,in_moleatom)
!
      origin(1)  = cr_pos(1,ia) + surf_normal(1)/normal_l*dist - posit(1)  ! Origin is shifted
      origin(2)  = cr_pos(2,ia) + surf_normal(2)/normal_l*dist - posit(2)  ! by dist away from 
      origin(3)  = cr_pos(3,ia) + surf_normal(3)/normal_l*dist - posit(3)  ! surface atom
      sym_latom(:) = .false.                           ! Initially deselect all atomtypes
!
      atoms: DO im=1,mole_natoms                       ! Load all atoms from the molecule
         CALL struc_read_one_atom_internal(mole_name, im, posit, itype, iprop, isurface,in_mole,in_moleatom)
         posit(:) = posit(:) + origin(:)
         WRITE(line, 1000) mole_atom_name(itype), posit, mole_dw(itype)
         laenge = 60
         CALL do_ins(line, laenge)                     ! Insert into crystal
         cr_prop (cr_natoms) = ibset (cr_prop (cr_natoms), PROP_LIGAND)
         cr_prop (cr_natoms) = ibset (cr_prop (cr_natoms), PROP_SURFACE_EXT)
         cr_surf(:,cr_natoms) = 0
         CALL check_symm
         sym_latom(cr_iscat(cr_natoms)) = .true.       ! Select atom type for rotation
      ENDDO atoms
!
      IF(mole_axis(0)==2) THEN    ! Rotate upright, if two atoms are given
        CALL rotate_directly(neig, n_atoms_orig, mole_axis, surf_normal)
!! define rotation operation
!      sym_angle      = do_bang(lspace, surf_normal, vnull, axis_ligand)
!      IF(ABS(sym_angle) > EPS ) THEN                   ! Rotate if not zero degrees
!         sym_orig(:)    = origin(:)                    ! Define origin
!         sym_trans(:)   = 0.0                          ! No translation needed
!         sym_sel_atom   = .true.                       ! Select atoms
!         sym_new        = .false.                      ! No new types
!         sym_power      =  1                           ! Just need one operation
!         sym_type       = .true.                       ! Proper rotation
!         sym_mode       = .false.                      ! Move atom to new position
!         sym_orig_mol   = .false.                      ! Origin at crystal
!         sym_power_mult =.false.                       ! No multiple copies
!         sym_sel_atom   = .true.                       ! Select atoms not molecules
!         sym_start      =  cr_natoms - mole_natoms + 1 ! set range of atoms numbers
!         sym_end        =  cr_natoms
!         IF(ABS(sym_angle-180.) < EPS ) THEN           ! Ligand and surface normal are antiparallel
!            WRITE(line,1100) axis_ligand(3)+0.1,axis_ligand(2)+0.01,axis_ligand(1)+0.001, surf_normal
!         ELSE
!            WRITE(line,1100) axis_ligand, surf_normal
!         ENDIF
!         laenge = 81
!         CALL vprod(line, laenge)                      ! Make rotation axis
!         sym_uvw(:) = res_para(1:3)
!         CALL trans (sym_uvw, cr_gten, sym_hkl, 3)     ! Make reciprocal space axis
!         CALL symm_setup                               ! Symmetry setup defines matrix
!         CALL symm_show                                ! Show only in debug
!         CALL symm_op_single       ianz, werte                    ! Perform the operation
!      ENDIF  ! Rotate if not zero degrees
      ENDIF
      CALL deco_tilt(origin, tilt, tilt_hkl, tilt_atom, tilt_is_atom, &
                     tilt_is_auto,                                    &
                     surf_normal, mole_natoms, 0, 0)
      IF(ABS(dist) < EPS ) THEN                        ! Remove surface atom
         cr_iscat(ia) = 0
         cr_prop (ia) = ibclr (cr_prop (ia), PROP_NORMAL)
      ENDIF
!
      m_type_new = m_type_old  + temp_id 
      CALL molecularize_numbers(nold+1,cr_natoms, m_type_new, r_m_biso, r_m_clin, r_m_cqua)
      flagsurf: DO j=1,20
               IF(mole_surfnew(j)>0) THEN
            im = nold + mole_surfnew(j)
            cr_prop  (im) = IBSET (cr_prop (im), PROP_SURFACE_EXT)
            cr_surf(:,im) = cr_surf(:, ia)          ! Copy surface vector from anchor
         ELSE
            EXIT flagsurf
         ENDIF
      ENDDO flagsurf
   ENDIF
!
   chem_period(:) = .false.                         ! We inserted atoms, turn off periodic boundaries
   chem_quick     = .false.                         ! turn of quick search
   cr_prop (ia) = IBCLR (cr_prop (ia), PROP_DECO_ANCHOR)  ! UNFLAG THIS ATOM AS SURFACE ANCHOR
   cr_prop (ia) = IBCLR (cr_prop (ia), PROP_SURFACE_EXT)   ! Anchor is no longer at a surface
   1000 FORMAT(a4,4(2x,',',F12.6))
!
   END SUBROUTINE deco_place_normal
!
!*******************************************************************************
!
SUBROUTINE deco_place_bridge(temp_id, ia,                                       &
                         mole_axis, mole_surfnew,                               &
                         m_type_old, mole_name,                                 &
                         mole_natoms,                                           &
                         mole_nscat, mole_atom_name,                            &
                         mole_dw, r_m_biso, r_m_clin, r_m_cqua,                 &
                         nanch, anchor,                                         &
                         neig, dist, ncon, istart, iend, lrestrict, nhkl, rhkl, &
                         tilt, tilt_hkl, tilt_atom, tilt_is_atom,               &
                         tilt_is_auto)
!
USE crystal_mod
USE atom_env_mod
USE chem_mod
USE do_find_mod
USE metric_mod
USE modify_mod
USE molecule_func_mod
USE point_grp
USE prop_para_mod
USE read_internal_mod
USE surface_func_mod
USE symm_menu
USE symm_mod
USE symm_sup_mod
USE trafo_mod
!
USE param_mod
!
IMPLICIT NONE
!
INTEGER,                    INTENT(IN) :: temp_id         ! The definition to be used
INTEGER,                    INTENT(IN) :: ia              ! Surface atom number
INTEGER,                    INTENT(IN) :: m_type_old      ! molecule types previous to all placements
INTEGER, DIMENSION(0:2),    INTENT(IN) :: mole_axis       ! Atoms that define molecule axis
INTEGER, DIMENSION(1:20),   INTENT(IN) :: mole_surfnew    ! Atoms that will be flagged as surface atoms
CHARACTER (LEN=1024),       INTENT(IN) :: mole_name       ! Molecule file name
INTEGER,                    INTENT(IN) :: mole_natoms     ! Number of atoms in molecule
INTEGER,                    INTENT(IN) :: mole_nscat      ! Number of atoms in molecule
CHARACTER (LEN=4   ), DIMENSION(0:mole_nscat),   INTENT(IN) :: mole_atom_name ! Atom names in the molecule
REAL                , DIMENSION(0:mole_nscat),   INTENT(IN) :: mole_dw        ! ADPs       in the molecule
REAL,                       INTENT(IN) :: r_m_biso        ! Molecular Biso
REAL,                       INTENT(IN) :: r_m_clin        ! Molecular linear correction
REAL,                       INTENT(IN) :: r_m_cqua        ! Molecular quadratic correction
INTEGER,                    INTENT(IN) :: nanch           ! Connected to this neighbor in mole
INTEGER, DIMENSION(1:2,0:nanch), INTENT(IN) :: anchor          ! Surface atom type
INTEGER, DIMENSION(1:2),    INTENT(IN) :: neig            ! Connected to this neighbor in mole
REAL   , DIMENSION(1:2),    INTENT(IN) :: dist            ! distance to ligand molecule
INTEGER,                    INTENT(IN) :: ncon            ! Number of defined bonds
INTEGER,                    INTENT(IN) :: istart          ! First atom for surface determination
INTEGER,                    INTENT(IN) :: iend            ! Last  atom for surface determination
LOGICAL,                    INTENT(IN) :: lrestrict       ! Restriction to a surface type T/F
INTEGER,                    INTENT(IN) :: nhkl            ! Number of faces for the restriction
INTEGER, DIMENSION(3,nhkl), INTENT(IN) :: rhkl            ! actual faces for the restriction
REAL                      , INTENT(IN) :: tilt            ! Molecule tilt angle
REAL,    DIMENSION(3)     , INTENT(IN) :: tilt_hkl        ! Molecule tilt plane by this normal
INTEGER, DIMENSION(4)     , INTENT(IN) :: tilt_atom       ! Molecule tilt plane defined by these atoms
LOGICAL                   , INTENT(IN) :: tilt_is_atom    ! Plane defined by atoms
LOGICAL                   , INTENT(IN) :: tilt_is_auto    ! Plane defined by atoms
!
INTEGER                                 :: surf_char      ! Surface character, plane, edge, corner, ...
INTEGER, DIMENSION(3,6)                 :: surface_normal ! Set of local normals (:,1) is main normal
INTEGER, DIMENSION(3)                   :: surf_kante     ! Edge vector if not a plane
INTEGER, DIMENSION(6)                   :: surf_weight    ! Best normal has heighest weight
INTEGER, DIMENSION(:), ALLOCATABLE      :: all_surface         ! Surface atom type
!
REAL   , DIMENSION(1:3)                 :: surf_normal    ! Normal to work with
!
INTEGER, PARAMETER                      :: MINPARA = 2
INTEGER                                 :: MAXW = MINPARA
REAL     , DIMENSION(1:MAX(MINPARA,nanch)) :: werte
!
CHARACTER (LEN=1024)                    :: line, zeile
INTEGER                                 :: ianz
INTEGER                                 :: i,j, im, laenge
INTEGER                                 :: l
INTEGER                                 :: iprop
INTEGER                                 :: itype
INTEGER  , DIMENSION(0:3)               :: isurface
INTEGER                                 :: nold   ! atom number previous to current molecule
INTEGER  , DIMENSION(1:4)               :: hkl
!
INTEGER                                 :: test_nhkl
INTEGER  , DIMENSION(:,:), ALLOCATABLE  :: test_hkl
!
LOGICAL  , DIMENSION(1:3)               :: fp
LOGICAL                                 :: fq
REAL                                    :: c_ang_ia, c_ang_nei         ! COS/SIN of Angles in ia and nei
REAL                                    :: angle, a_test               ! Dummy angles
REAL                                    :: dist_m
REAL     , DIMENSION(1:3)               :: axis_ligand                 ! Initial molecule orientation
REAL                                    :: rmin, radius, normal_l, b_l, b_n, b_l_min
REAL     , DIMENSION(1:3)               :: x, bridge, tangent, origin, posit
REAL     , DIMENSION(1:3)               :: vnull
REAL     , DIMENSION(1:3)               :: pos1, pos2
INTEGER             :: m_type_new   ! new molecule types 
INTEGER             :: in_mole,in_moleatom
INTEGER             :: n_atoms_orig
!
maxw     = MAX(MINPARA,nanch)
!
vnull(:) = 0.00
nold = cr_natoms                                   ! Remember old atom number
!
!  Determine surface character, if growth is restricted check if we're at proper surface
!
hkl(4) = 0
CALL surface_character(ia, istart, iend, surf_char, surface_normal, surf_kante, surf_weight, .TRUE.)
surf_normal(1:3) = FLOAT(surface_normal(:,1))
IF(surf_char == 0) RETURN
hkl(1:3) = surface_normal(:,1)
IF(lrestrict) THEN
   test_nhkl =    nhkl
   ALLOCATE(test_hkl(3,test_nhkl))
   test_hkl(1:3,1:test_nhkl) =    rhkl(1:3,1:test_nhkl)
   IF(.NOT.point_test(hkl, test_hkl, test_nhkl, .TRUE.) ) THEN
      RETURN
   ENDIF
   DEALLOCATE(test_hkl)
ENDIF
!
n_atoms_orig = cr_natoms                                  ! Store original atom numbers
!
ALLOCATE(all_surface(1:ncon))
!
all_surface (1) = ia
!
!  FIND the other surface partners involved in the bonds.
!
x(1)     = cr_pos(1,ia)
x(2)     = cr_pos(2,ia)
x(3)     = cr_pos(3,ia)
!
angle   = 400.0
b_l_min = 1.E12
!
search: DO l=2,ncon
   rmin     = 0.1
   radius   = dist(1) + dist(2) - 0.1
   ianz = anchor(2,0)
   werte(1:ianz) = anchor(2,1:ianz) 
   fp (:)   = chem_period (:)
      fq    = chem_quick
   CALL do_find_env (ianz, werte, maxw, x, rmin, radius, fq, fp)  ! Find all neighbors
   j = 0
   IF(atom_env(0) >= 1 ) THEN                                     ! We need at least one neighbor
      pos1(:) = cr_pos(:,ia)                                       ! Temporary pos of atom 1=ia
      angle   = 400.0
      b_l_min = 1.E12
!     j = 0
      check_prop: DO i=1,atom_env(0)                               ! Check properties 
         IF(IBITS(cr_prop(atom_env(i)),PROP_SURFACE_EXT,1).eq.1 .and.        &  ! real Atom is near surface
            IBITS(cr_prop(atom_env(i)),PROP_OUTSIDE    ,1).eq.0       ) THEN    ! real Atom is near surface
            pos2(:) = cr_pos(:,atom_env(i))                       ! temporarily store atom 2
            bridge(:) = pos1(:) - pos2(:)
            b_l = sqrt (skalpro (bridge, bridge, cr_gten))        ! Calculate bridge length
            IF(b_l<=b_l_min) THEN
               IF(b_l<b_l_min) THEN
                  b_l_min = b_l
                  angle   = 400.0
               ENDIF
               a_test = do_bang(.TRUE., bridge,vnull, surf_normal)
               IF(ABS(a_test-90.0) < angle) THEN
                  j = i                                           ! Will use this neighbor
                  angle = ABS(a_test - 90.)
               ENDIF
            ENDIF
         ENDIF
      ENDDO check_prop
      IF(j==0) THEN                                      ! No suitable neighbor, quietly leave
         GOTO 9999
      ENDIF
      all_surface (l) = atom_env(j)
   ELSE
      GOTO 9999
   ENDIF  ! 
ENDDO search
!
!
pos2(:)   = cr_pos(:,all_surface(2))
bridge(1) = (cr_pos(1,ia)-cr_pos(1,all_surface(2)))   ! Calculate vector along bridge
bridge(2) = (cr_pos(2,ia)-cr_pos(2,all_surface(2)))
bridge(3) = (cr_pos(3,ia)-cr_pos(3,all_surface(2)))
b_l       = sqrt (skalpro (bridge, bridge, cr_gten))     ! Calculate bridge length
c_ang_ia  = (b_l**2 + dist(1)**2 -dist(2)**2) / (2.*b_l*dist(1))   ! COS(Angle in atom ia)
IF(ABS(c_ang_ia) > 1.00) GOTO 9999
!
c_ang_nei = (b_l**2 + dist(2)**2 -dist(1)**2) / (2.*b_l*dist(2))   ! COS(Angle in atom ia)
IF(ABS(c_ang_nei) > 1.00) GOTO 9999
!
dist_m    = (dist(1) * c_ang_ia) / b_l            ! relative length of 'midpoint' location
x(1) = cr_pos(1,ia) + (cr_pos(1,all_surface(2))-cr_pos(1,ia))*dist_m ! Calculate midpoint
x(2) = cr_pos(2,ia) + (cr_pos(2,all_surface(2))-cr_pos(2,ia))*dist_m
x(3) = cr_pos(3,ia) + (cr_pos(3,all_surface(2))-cr_pos(3,ia))*dist_m
!
IF(angle>0.0) THEN                                 ! Bridge is not normal to surface_ normal
   WRITE(line, 1100) bridge, surf_normal
   laenge = LEN_TRIM(line)
   CALL vprod(line, laenge)
   tangent(1:3) = res_para(1:3)
   WRITE(line, 1100) tangent, bridge
   laenge = LEN_TRIM(line)
   CALL vprod(line, laenge)
   surf_normal(:) = res_para(1:3)
ENDIF
!
IF(mole_axis(0)==2) THEN
   im   = mole_axis(2)
   CALL struc_read_one_atom_internal(mole_name, im, posit, itype, iprop, isurface,in_mole,in_moleatom)
   axis_ligand(1) = posit(1)
   axis_ligand(2) = posit(2)
   axis_ligand(3) = posit(3)
   im = mole_axis(1)
CALL struc_read_one_atom_internal(mole_name, im, posit, itype, iprop, isurface,in_mole,in_moleatom)
   axis_ligand(1) = axis_ligand(1) - posit(1)
   axis_ligand(2) = axis_ligand(2) - posit(2)
   axis_ligand(3) = axis_ligand(3) - posit(3)
ENDIF
!     Insert molecule atoms into crystal at correct origin in initial orientation
normal_l = sqrt (skalpro (surf_normal, surf_normal, cr_gten))
b_n      = dist(1) * sin(acos(c_ang_ia))                 ! Calculate distance along normal
!
!     Get position of "Neig" and subtract to ensure correct origin
!     in read_crystal this atom is always shifted to (0,0,0)
CALL struc_read_one_atom_internal(mole_name, neig(1), posit, itype, iprop, isurface,in_mole,in_moleatom)
origin(1)  = x(1) + surf_normal(1)/normal_l*b_n - posit(1)    ! Calculate ligand origin
origin(2)  = x(2) + surf_normal(2)/normal_l*b_n - posit(2)
origin(3)  = x(3) + surf_normal(3)/normal_l*b_n - posit(3)
sym_latom(:) = .false.                        ! Initially deselect all atomtypes
DO im=1,mole_natoms                           ! Insert all atoms
   CALL struc_read_one_atom_internal(mole_name, im, posit, itype, iprop, isurface,in_mole,in_moleatom)
   posit(:) = posit(:) + origin(:)
   WRITE(line, 1000) mole_atom_name(itype), posit, mole_dw(itype)
   laenge = 60
   zeile = line
   CALL do_ins(line, laenge)
   cr_prop (cr_natoms) = ibset (cr_prop (cr_natoms), PROP_LIGAND)
   cr_prop (cr_natoms) = ibset (cr_prop (cr_natoms), PROP_SURFACE_EXT)
   cr_surf(:,cr_natoms) = 0
   CALL check_symm
   sym_latom(cr_iscat(cr_natoms)) = .true.    ! Select atopm type for rotation
ENDDO
IF(cr_natoms==nold) GOTO 9999
!
IF(mole_axis(0)==2) THEN    ! Rotate upright, if two atoms are given
  CALL rotate_directly(neig(1), n_atoms_orig, mole_axis, surf_normal)
ENDIF
!
!     Tilt molecule by user request
!
CALL deco_tilt(origin, tilt, tilt_hkl, tilt_atom, tilt_is_atom, &
               tilt_is_auto,                                    &
               surf_normal, mole_natoms, 0, 0)
m_type_new = m_type_old  + temp_id
!
CALL molecularize_numbers(nold+1,cr_natoms, m_type_new, r_m_biso, r_m_clin, r_m_cqua)
flagsurf: DO j=1,20
   IF(mole_surfnew(j)>0) THEN
      im = nold + mole_surfnew(j)
      cr_prop  (im) = IBSET (cr_prop (im), PROP_SURFACE_EXT)
      cr_surf(:,im) = cr_surf(:, ia)          ! Copy surface vector from anchor
   ELSE
      EXIT flagsurf
   ENDIF
ENDDO flagsurf
!
9999 CONTINUE                                    ! Jump here from errors to ensure dealloc
IF(nold<cr_natoms) THEN                          ! We did insert a molecule
    chem_period(:) = .false.                         ! We inserted atoms, turn off periodic boundaries
    chem_quick     = .false.                         ! turn of quick search
    cr_prop (ia) = IBCLR (cr_prop (ia), PROP_DECO_ANCHOR)   ! UNFLAG THIS ATOM AS SURFACE ANCHOR
    cr_prop (ia) = IBCLR (cr_prop (ia), PROP_SURFACE_EXT)   ! Anchor is no longer at a surface
    j  = all_surface(2)
    cr_prop (j ) = IBCLR (cr_prop (j ), PROP_DECO_ANCHOR)   ! UNFLAG THIS ATOM AS SURFACE ANCHOR
    cr_prop (j ) = IBCLR (cr_prop (j ), PROP_SURFACE_EXT)   ! Anchor is no longer at a surface
ENDIF
DEALLOCATE(all_surface)
!
1000 FORMAT(a4,4(2x,',',F12.6))
1100 FORMAT(6(G15.6E3,', '),'ddd')
!
END SUBROUTINE deco_place_bridge
!
!*******************************************************************************
!
   SUBROUTINE deco_place_double(temp_id, ia, mole_axis, &
                            mole_surfnew,           &
                            m_type_old, mole_name, &
                            mole_natoms,           &
                            mole_nscat, mole_atom_name, &
                            mole_dw, r_m_biso, r_m_clin, r_m_cqua,         &
                            nanch,anchor, neig, dist,ncon, &
                            istart, iend, lrestrict, nhkl, rhkl,   &
                            tilt, tilt_hkl, tilt_atom, tilt_is_atom, &
                            tilt_is_auto)
!
!  The molecule is bound by two of its atoms to two different surfae atoms.
!  The molecule is placed such that:
!     the vector between the two binding molecule atoms is parallel to the
!     vector of the two surface atoms.
!     The resulting trapezoid is parallel to the surface normal
!     The remainder of the molecule is rotated around the vector between
!     the two molecule atoms such that the molecule axis is as parallel as
!     possible to the surface normal.
!
USE crystal_mod
   USE atom_env_mod
   USE chem_mod
USE do_find_mod
   USE metric_mod
   USE modify_mod
   USE molecule_func_mod
   USE point_grp
USE prop_para_mod
   USE read_internal_mod
USE surface_func_mod
   USE symm_menu
   USE symm_mod
   USE symm_sup_mod
   USE trafo_mod
!
   USE param_mod
   USE wink_mod
!
   IMPLICIT NONE
!
   INTEGER,                 INTENT(IN) ::    temp_id      ! The definition to be used
   INTEGER,                 INTENT(IN) :: ia              ! Surface atom number
   INTEGER, DIMENSION(0:2), INTENT(IN) :: mole_axis       ! Atoms that define molecule axis
   INTEGER, DIMENSION(1:20),INTENT(IN) :: mole_surfnew    ! Atoms that will be flagged as surface atoms
   INTEGER,                 INTENT(IN) :: m_type_old      ! molecule types previous to all placements
   CHARACTER (LEN=1024),    INTENT(IN) :: mole_name       ! Molecule file name
   INTEGER,                 INTENT(IN) :: mole_natoms     ! Molecule file name length
   INTEGER,                 INTENT(IN) :: mole_nscat      ! Number of atoms in molecule
   CHARACTER (LEN=4   ), DIMENSION(0:mole_nscat),   INTENT(IN) :: mole_atom_name ! Atom names in the molecule
   REAL                , DIMENSION(0:mole_nscat),   INTENT(IN) :: mole_dw        ! ADPs       in the molecule
REAL,                    INTENT(IN) :: r_m_biso        ! Molecular Biso
REAL,                    INTENT(IN) :: r_m_clin        ! Molecular linear correction
REAL,                    INTENT(IN) :: r_m_cqua        ! Molecular quadratic correction
   INTEGER,                 INTENT(IN) :: nanch           ! Number of anchor types
   INTEGER, DIMENSION(1:2,0:nanch), INTENT(IN) :: anchor          ! Surface atom type
   INTEGER, DIMENSION(1:2), INTENT(IN) :: neig            ! Connected to this neighbor in mole
   REAL   , DIMENSION(1:2), INTENT(IN) :: dist            ! distance to ligand molecule
   INTEGER,                 INTENT(IN) :: ncon            ! Number of defined bonds
   INTEGER,                 INTENT(IN) :: istart          ! First atom for surface determination
   INTEGER,                 INTENT(IN) :: iend            ! Last  atom for surface determination
   LOGICAL,                 INTENT(IN) :: lrestrict       ! Restriction to a surface type T/F
   INTEGER,                 INTENT(IN) :: nhkl            ! Number of faces for the restriction
   INTEGER, DIMENSION(3,nhkl), INTENT(IN) :: rhkl            ! actual faces for the restriction
REAL                      , INTENT(IN) :: tilt            ! Molecule tilt angle
REAL,    DIMENSION(3)     , INTENT(IN) :: tilt_hkl        ! Molecule tilt plane by this normal
INTEGER, DIMENSION(4)     , INTENT(IN) :: tilt_atom       ! Molecule tilt plane defined by these atoms
LOGICAL                   , INTENT(IN) :: tilt_is_atom    ! Plane defined by atoms
LOGICAL                   , INTENT(IN) :: tilt_is_auto    ! Plane defined by atoms
!
   INTEGER, PARAMETER                      :: MINPARA = 2
   INTEGER                                 :: MAXW = MINPARA
   REAL                , DIMENSION(1:MAX(MINPARA,nanch)) :: werte
!
   CHARACTER (LEN=1024)                    :: line
   INTEGER, DIMENSION(:), ALLOCATABLE  :: all_surface         ! Surface atom type
   INTEGER                                 :: ianz
   INTEGER                                 :: i,j, l,im, laenge
   INTEGER                                 :: iprop
   INTEGER                                 :: itype
   INTEGER, DIMENSION(0:3) :: isurface       ! Atom surface
   INTEGER                                 :: n_atoms_orig   ! Number of atoms prior to insertion
   INTEGER                                 :: n1,n2          ! number of mol neighbours after insertion
   INTEGER                                 :: success        ! Everything went fine
   INTEGER  , DIMENSION(1:4)               :: hkl
   LOGICAL  , DIMENSION(1:3)               :: fp
   LOGICAL                                 :: fq
   LOGICAL, PARAMETER :: lspace = .true.
   REAL                                    :: rmin, radius, b_l, t_l
   REAL                                    :: arg            ! argument for acos
   REAL     , DIMENSION(1:3)               :: x, bridge, tangent, origin, posit, v, w
   REAL     , DIMENSION(1:3)               :: shift, v1, v2
   REAL     , DIMENSION(1:3)               :: vnull
   INTEGER                 :: surf_char      ! Surface character, plane, edge, corner, ...
   INTEGER, DIMENSION(3,6) :: surface_normal ! Set of local normals (:,1) is main normal
   INTEGER, DIMENSION(3)   :: surf_kante     ! Edge vector if not a plane
   INTEGER, DIMENSION(6)   :: surf_weight    ! Best normal has heighest weight
!
   REAL   , DIMENSION(1:3) :: surf_normal    ! Normal to work with
!
   INTEGER                              :: test_nhkl
   INTEGER, DIMENSION(:,:), ALLOCATABLE :: test_hkl
!
   INTEGER             :: nold         ! atom number previous to current molecule
   INTEGER             :: m_type_new   ! new molecule types 
   INTEGER             :: in_mole,in_moleatom
!
maxw     = MAX(MINPARA,nanch)
vnull(:) = 0.00
success = -1
n1 = 1
n2 = 1
b_l = 0.0
nold = cr_natoms                           ! Remember original atom number
!
!  Determine surface character, if growth is restricted check if we're at proper surface
!
hkl(4) = 0
CALL surface_character(ia, istart, iend, surf_char, surface_normal, surf_kante, surf_weight, .TRUE.)
IF(surf_char == 0) RETURN
surf_normal(1:3) = FLOAT(surface_normal(:,1))
hkl(1:3) = surface_normal(:,1)
IF(lrestrict) THEN
   test_nhkl =    nhkl
   ALLOCATE(test_hkl(3,test_nhkl))
   test_hkl(1:3,1:test_nhkl) =    rhkl(1:3,1:test_nhkl)
   IF(.NOT.point_test(hkl, test_hkl, test_nhkl, .TRUE.) ) THEN
      RETURN
   ENDIF
   DEALLOCATE(test_hkl)
ENDIF
!
!  Load the molecule into the crystal structure
!
n_atoms_orig = cr_natoms                         ! Number of atoms prior to insertion
nold = cr_natoms                                 ! Remember original atom number
CALL struc_read_one_atom_internal(mole_name, neig(1), posit, itype, iprop, isurface,in_mole,in_moleatom)
origin(:) = 0.0 -posit(:)                     ! initially place at 0,0,0
sym_latom(:) = .false.                        ! Initially deselect all atomtypes
insert: DO im=1,mole_natoms                   ! Insert all atoms
      CALL struc_read_one_atom_internal(mole_name, im, posit, itype, iprop, isurface,in_mole,in_moleatom)
      posit(:) = posit(:) + origin(:)
      WRITE(line, 1000) mole_atom_name(itype), posit, mole_dw(itype)
      laenge = 60
      CALL do_ins(line, laenge)
      cr_prop (cr_natoms) = ibset (cr_prop (cr_natoms), PROP_LIGAND)
      cr_prop (cr_natoms) = ibset (cr_prop (cr_natoms), PROP_SURFACE_EXT)
      cr_surf(:,cr_natoms) = 0
      CALL check_symm
      sym_latom(cr_iscat(cr_natoms)) = .true.    ! Select atopm type for rotation
ENDDO insert
!
chem_period(:) = .false.                         ! We inserted atoms, turn off periodic boundaries
chem_quick     = .false.                         ! turn of quick search
!
ALLOCATE(all_surface(1:ncon))
!
all_surface (1) = ia
!
!  FIND the other surface partners involved in the bonds.
!
x(1)     = cr_pos(1,ia)
x(2)     = cr_pos(2,ia)
x(3)     = cr_pos(3,ia)
search: DO l=2,ncon
   n2 = n_atoms_orig +     neig(l)   
   n1 = n_atoms_orig +     neig(1)
   bridge(1) = cr_pos(1,n2) - cr_pos(1,n1)
   bridge(2) = cr_pos(2,n2) - cr_pos(2,n1)
   bridge(3) = cr_pos(3,n2) - cr_pos(3,n1)
   b_l      = sqrt(skalpro(bridge, bridge, cr_gten))
   rmin     = MAX( 0.1, b_l - dist(1) - dist(2) )          ! Minimum distance between surface atoms
   radius   = b_l + dist(1) + dist(2)                      ! Maximum distance between surface atoms
   ianz = anchor(2,0)
   werte(1:ianz) = anchor(2,1:ianz) 
   fp (:)   = chem_period (:)
      fq    = chem_quick
   CALL do_find_env (ianz, werte, maxw, x, rmin, radius, fq, fp)  ! Find all neighbors
   IF(atom_env(0) >= 1 ) THEN                                     ! We need at least one neighbor
     j = 0
     check_prop: DO i=1,atom_env(0)                               ! Check properties 
        IF(IBITS(cr_prop(atom_env(i)),PROP_SURFACE_EXT,1).eq.1 .and.        &  ! real Atom is near surface
           IBITS(cr_prop(atom_env(i)),PROP_OUTSIDE    ,1).eq.0       ) THEN    ! real Atom is near surface
            j = i                                                 ! Will use this neighbor
            EXIT check_prop                                       ! Found first good neighbor
         ENDIF
      ENDDO check_prop
      IF(j==0) THEN                                      ! No suitable neighbor, quietly leave
         GOTO 9999
      ENDIF
      all_surface (l) = atom_env(j)
   ELSE
      GOTO 9999
   ENDIF  ! 
ENDDO search
!  Determine rotation axis for surface vector
tangent(:) = cr_pos(:,all_surface (2)) - cr_pos(:,ia)  ! Vector between surface atoms
t_l        = sqrt(skalpro(tangent, tangent, cr_gten))  ! Distance between surface atoms
WRITE(line,1100) tangent, surf_normal                  ! Calculate rotation axis
laenge = LEN_TRIM(line)
CALL vprod(line, laenge)
sym_uvw(:) = res_para(1:3)
!  Calculate angle in first trapezoid corner
sym_angle  = ACOS(-(dist(2)**2-dist(1)**2-(t_l-b_l)**2)/ &
                   (2.*dist(1)*(t_l-b_l))) /REAL(rad)
!
sym_orig(:)    = 0.0                       ! Define origin at 0,0,0
sym_trans(:)   = 0.0                       ! No translation needed
sym_sel_atom   = .true.                    ! Select atoms
sym_new        = .false.                   ! No new types
sym_power      =  1                        ! Just need one operation
sym_type       = .true.                    ! Proper rotation
sym_mode       = .false.                   ! Move atom to new position
sym_orig_mol   = .false.                   ! Origin at crystal
sym_power_mult =.false.                    ! No multiple copies
sym_sel_atom   = .true.                    ! Select atoms not molecules
CALL trans (sym_uvw, cr_gten, sym_hkl, 3)  ! Make reciprocal space axis
CALL symm_setup
v(:) =  tangent(:)*dist(1)/t_l     ! Scale vector to bond length
CALL symm_ca_single (v, .true., .false.)
origin(:) = cr_pos(:,all_surface(1)) + res_para(1:3)  ! Origin of the molecule = surface 1 + result 
shift (:) = origin(:) - cr_pos(:,n1)       ! All molecule atoms need to be shifted by this vector
DO i=n_atoms_orig+1,cr_natoms
   cr_pos(:,i) = cr_pos(:,i) + shift(:)
ENDDO
!
!  Calculate angle in second trapezoid corner
arg        = (-(dist(1)**2-dist(2)**2-(t_l-b_l)**2)/ &
               (2.*dist(2)*(t_l-b_l)))
IF(ABS(arg) > 1.00) THEN
   GOTO 9999                               ! No solution is found skip
ENDIF
sym_angle  = ACOS(arg)/REAL(rad)
!
sym_trans(:)   = 0.0                       ! No translation needed
sym_orig(:)    = 0.0                       ! Define origin in 0,0,0
sym_uvw(:)     = -sym_uvw(:)               ! invert axis
v(:) =  -tangent(:)*dist(2)/t_l    ! invert and scale vector to bond length
CALL trans (sym_uvw, cr_gten, sym_hkl, 3)  ! Make reciprocal space axis
CALL symm_setup
CALL symm_ca_single (v, .true., .false.)   ! rotate negative surface vector
w(:) = cr_pos(:,all_surface(2)) + res_para(1:3)  ! Target position for 2nd molecule atom
!
sym_orig(:)    = cr_pos(:,n1)              ! Define origin in 1st attached molecule atom
v1(:) = cr_pos(:,n2) - cr_pos(:,n1)        ! Current vector from 1st to 2nd molecule atom
v2(:) = w(:)         - cr_pos(:,n1)        ! Vector from 1st to target  2nd molecule atom
WRITE(line,1100) v1, v2                    ! Rotation axis will be v1 x v2
laenge = LEN_TRIM(line)
CALL vprod(line, laenge)
sym_uvw(:) = res_para(1:3)
IF(res_para(1)**2+res_para(2)**2+res_para(3)**2 >  1e-5) THEN !Non-zero axis
sym_trans(:)   = 0.0                          ! No translation needed
sym_angle  = do_bang(lspace, v1, vnull, v2)   ! Calculate rotation angle = < (v1,v2)
sym_start  =  n_atoms_orig + 1                ! set range of atoms numbers
sym_end    =  cr_natoms
CALL trans (sym_uvw, cr_gten, sym_hkl, 3)     ! Make reciprocal space axis
CALL symm_setup
CALL symm_op_single                           ! Perform the operation
ELSE
   GOTO 9999
ENDIF
!
IF(mole_axis(0)==2) THEN    ! Rotate upright, if two atoms are given
   CALL rotate_projected(n1, n2, n_atoms_orig, mole_axis, surf_normal)
ENDIF
cr_prop (all_surface(2)) = IBCLR (cr_prop (all_surface(2)), PROP_DECO_ANCHOR)  ! UNFLAG THIS ATOM AS SURFACE ANCHOR
!
!     Tilt molecule by user request
!
origin(1:3) = cr_pos(1:3, n1)
CALL deco_tilt(origin, tilt, tilt_hkl, tilt_atom, tilt_is_atom, &
               tilt_is_auto,                                    &
               surf_normal, mole_natoms, n1, n2)
!
m_type_new = m_type_old  + temp_id
CALL molecularize_numbers(nold+1,cr_natoms, m_type_new, r_m_biso, r_m_clin, r_m_cqua)
flagsurf: DO j=1,20
   IF(mole_surfnew(j)>0) THEN
      im = nold + mole_surfnew(j)
      cr_prop  (im) = IBSET (cr_prop (im), PROP_SURFACE_EXT)
      cr_surf(:,im) = cr_surf(:, ia)          ! Copy surface vector from anchor
   ELSE
      EXIT flagsurf
   ENDIF
ENDDO flagsurf
success = 0 
!
9999 CONTINUE                                             ! Jump here from errors to ensure dealloc
IF(success /=0) THEN                          ! An error occurred, reset crystal
   cr_natoms = n_atoms_orig
ENDIF
IF(nold<cr_natoms) THEN                          ! We did insert a molecule
   chem_period(:) = .false.                         ! We inserted atoms, turn off periodic boundaries
   chem_quick     = .false.                         ! turn of quick search
   cr_prop (ia) = IBCLR (cr_prop (ia), PROP_DECO_ANCHOR)   ! UNFLAG THIS ATOM AS SURFACE ANCHOR
   cr_prop (ia) = IBCLR (cr_prop (ia), PROP_SURFACE_EXT)   ! Anchor is no longer at a surface
   j  = all_surface(2)
   cr_prop (j ) = IBCLR (cr_prop (j ), PROP_DECO_ANCHOR)   ! UNFLAG THIS ATOM AS SURFACE ANCHOR
   cr_prop (j ) = IBCLR (cr_prop (j ), PROP_SURFACE_EXT)   ! Anchor is no longer at a surface
ENDIF
DEALLOCATE(all_surface)
!
1000 FORMAT(a4,4(2x,',',F12.6))
1100 FORMAT(6(G15.6E3,', '),'ddd')
!
END SUBROUTINE deco_place_double
!
!*******************************************************************************
!
SUBROUTINE deco_place_chelate(temp_id, ia, &
                            mole_axis, mole_surfnew,           &
                            m_type_old, mole_name, &
                            mole_natoms, &
                            mole_nscat, mole_atom_name, &
                            mole_dw, r_m_biso, r_m_clin, r_m_cqua,         &
                            nanch, anchor, &
                            neig, dist, ncon,istart, iend, lrestrict, nhkl, rhkl, &
                            tilt, tilt_hkl, tilt_atom, tilt_is_atom, &
                            tilt_is_auto)
!
! Place the ligand in a "chelate" bond. 
! The two connecting ligand atoms are placed along a line normal to the local surface normal.
! The molecule axis is rotated around these two atoms into a plane that contains the 
! local surface normal.
!
USE crystal_mod
USE atom_env_mod
USE chem_mod
USE metric_mod
USE modify_mod
USE molecule_func_mod
USE point_grp
USE prop_para_mod
USE read_internal_mod
USE surface_func_mod
USE symm_menu
USE symm_mod
USE symm_sup_mod
USE trafo_mod
!
USE param_mod
USE trig_degree_mod
!
IMPLICIT NONE
!
INTEGER,                 INTENT(IN) :: temp_id         ! The definition to be used
INTEGER,                 INTENT(IN) :: ia              ! Surface atom number
INTEGER,                 INTENT(IN) :: m_type_old      ! molecule types previous to all placements
INTEGER, DIMENSION(0:2), INTENT(IN) :: mole_axis       ! Atoms that define molecule axis
INTEGER, DIMENSION(1:20),INTENT(IN) :: mole_surfnew    ! Atoms that will be flagged as surface atoms
CHARACTER (LEN=1024),    INTENT(IN) :: mole_name       ! Molecule file name
INTEGER,                 INTENT(IN) :: mole_natoms     ! Number of atoms in molecule
INTEGER,                 INTENT(IN) :: mole_nscat      ! Number of atoms in molecule
CHARACTER (LEN=4   ), DIMENSION(0:mole_nscat),   INTENT(IN) :: mole_atom_name ! Atom names in the molecule
REAL                , DIMENSION(0:mole_nscat),   INTENT(IN) :: mole_dw        ! ADPs       in the molecule
REAL,                    INTENT(IN) :: r_m_biso        ! Molecular Biso
REAL,                    INTENT(IN) :: r_m_clin        ! Molecular linear correction
REAL,                    INTENT(IN) :: r_m_cqua        ! Molecular quadratic correction
INTEGER,                 INTENT(IN) :: nanch           ! Connected to this neighbor in mole
INTEGER, DIMENSION(1:2,0:nanch), INTENT(IN) :: anchor          ! Surface atom type
INTEGER, DIMENSION(1:2), INTENT(IN) :: neig            ! Connected to this neighbor in mole
REAL   , DIMENSION(1:2), INTENT(IN) :: dist            ! distance to ligand molecule
INTEGER,                 INTENT(IN) :: ncon            ! Number of defined bonds
INTEGER,                 INTENT(IN) :: istart          ! First atom for surface determination
INTEGER,                 INTENT(IN) :: iend            ! Last  atom for surface determination
LOGICAL,                 INTENT(IN) :: lrestrict       ! Restriction to a surface type T/F
INTEGER,                 INTENT(IN) :: nhkl            ! Number of faces for the restriction
INTEGER, DIMENSION(3,nhkl), INTENT(IN) :: rhkl            ! actual faces for the restriction
REAL                      , INTENT(IN) :: tilt            ! Molecule tilt angle
REAL,    DIMENSION(3)     , INTENT(IN) :: tilt_hkl        ! Molecule tilt plane by this normal
INTEGER, DIMENSION(4)     , INTENT(IN) :: tilt_atom       ! Molecule tilt plane defined by these atoms
LOGICAL                   , INTENT(IN) :: tilt_is_atom    ! Plane defined by atoms
LOGICAL                   , INTENT(IN) :: tilt_is_auto    ! Plane defined by atoms
!
INTEGER                 :: surf_char      ! Surface character, plane, edge, corner, ...
INTEGER, DIMENSION(3,6) :: surface_normal ! Set of local normals (:,1) is main normal
INTEGER, DIMENSION(3)   :: surf_kante     ! Edge vector if not a plane
INTEGER, DIMENSION(6)   :: surf_weight    ! Best normal has heighest weight
!
REAL   , DIMENSION(1:3) :: surf_normal    ! Normal to work with
!
REAL, PARAMETER :: EPS = 1.0E-6
!
CHARACTER (LEN=1024)                    :: line
INTEGER                                 :: j, im, laenge
INTEGER                                 :: iprop
INTEGER                                 :: itype
INTEGER, DIMENSION(0:3)                 :: isurface
INTEGER                                 :: nold   ! atom number previous to current molecule
INTEGER  , DIMENSION(1:4)               :: hkl
INTEGER                                 :: n1, n2  ! Ligand atoms
!
INTEGER                              :: test_nhkl
INTEGER, DIMENSION(:,:), ALLOCATABLE :: test_hkl
!
LOGICAL, PARAMETER :: lspace = .true.
REAL                    :: aa, bb, cc, arg  ! Triangle sides for cosine theorem
REAL                                    :: normal_l
REAL     , DIMENSION(1:3)               :: x, origin, posit
REAL     , DIMENSION(1:3)               :: vnull
INTEGER             :: m_type_new   ! new molecule types 
INTEGER             :: in_mole,in_moleatom
!
vnull(:) = 0.00
nold     = cr_natoms
!
!  Determine surface character, if growth is restricted check if we're at proper surface
!
hkl(4) = 0
CALL surface_character(ia, istart, iend, surf_char, surface_normal, surf_kante, surf_weight, .TRUE.)
IF(surf_char == 0) RETURN
surf_normal(1:3) = FLOAT(surface_normal(:,1))
hkl(1:3) = surface_normal(:,1)
IF(lrestrict) THEN
   test_nhkl =    nhkl
   ALLOCATE(test_hkl(3,test_nhkl))
   test_hkl(1:3,1:test_nhkl) =    rhkl(1:3,1:test_nhkl)
      IF(.NOT.point_test(hkl, test_hkl, test_nhkl, .TRUE.) ) THEN
         RETURN
      ENDIF
      DEALLOCATE(test_hkl)
ENDIF
!
!  Insert molecule atoms into crystal at initial origin along normal
normal_l = sqrt (skalpro (surf_normal, surf_normal, cr_gten))
origin(1)  = cr_pos(1,ia) + surf_normal(1)/normal_l*dist(1)
origin(2)  = cr_pos(2,ia) + surf_normal(2)/normal_l*dist(1)
origin(3)  = cr_pos(3,ia) + surf_normal(3)/normal_l*dist(1)
sym_latom(:) = .false.                        ! Initially deselect all atomtypes
DO im=1,mole_natoms                           ! Insert all atoms
   CALL struc_read_one_atom_internal(mole_name, im, posit, itype, iprop, isurface,in_mole,in_moleatom)
   posit(:) = posit(:) + origin(:)
   WRITE(line, 1000) mole_atom_name(itype), posit, mole_dw(itype)
   laenge = 60
   CALL do_ins(line, laenge)
   cr_prop (cr_natoms) = ibset (cr_prop (cr_natoms), PROP_LIGAND)
   cr_prop (cr_natoms) = ibset (cr_prop (cr_natoms), PROP_SURFACE_EXT)
   cr_surf(:,cr_natoms) = 0
   CALL check_symm
   sym_latom(cr_iscat(cr_natoms)) = .true.    ! Select atom type for rotation
ENDDO
n1 = nold + neig(1)                           ! Atom number for 1st
n2 = nold + neig(2)                           ! and 2nd bonded ligand atoms
!
! 1.st step rotate ligand molecule in neig(1) to achieve the intended 
! distance dist(2) between surface atom ia and neig(2)
!
!  Find the vector from ligand atoms neig(1) and neig(2), the ligand_chelate axis
x(1) = cr_pos(1,n2) - cr_pos(1,n1)
x(2) = cr_pos(2,n2) - cr_pos(2,n1)
x(3) = cr_pos(3,n2) - cr_pos(3,n1)
cc = sqrt (skalpro (x, x, cr_gten))           ! Distance between chelate atoms
aa = dist(2)                                  ! Intended distance surface to ligand(2)
bb = dist(1)                                  ! Intended distance surface to ligand(1)
arg = (aa*aa - bb*bb -cc*cc)/(-2.*bb*cc)      ! Arg of acos
IF(arg > 1) THEN                              ! No solution exists
   cr_natoms = nold
   RETURN
ENDIF
!
! Define rotation operation
! Determine rotation axis by vector product ligand_chelate x surface_normal
! Angle is 180 - current angle (chelate to Normal) - intended angle
sym_angle = -(180. - do_bang(lspace, surf_normal, vnull, x) - acosd(arg))
IF(ABS(sym_angle) > EPS ) THEN                ! Rotate if not zero degrees
   sym_orig(:)    = cr_pos(:,n1)              ! Define origin in neig(1)
   sym_trans(:)   = 0.0                       ! No translation needed
   sym_sel_atom   = .true.                    ! Select atoms
   sym_new        = .false.                   ! No new types
   sym_power      =  1                        ! Just need one operation
   sym_type       = .true.                    ! Proper rotation
   sym_mode       = .false.                   ! Move atom to new position
   sym_orig_mol   = .false.                   ! Origin at crystal
   sym_power_mult =.false.                    ! No multiple copies
   sym_sel_atom   = .true.                    ! Select atoms not molecules
   sym_start      =  cr_natoms - mole_natoms + 1 ! set range of atoms numbers
   sym_end        =  cr_natoms
   IF(ABS(sym_angle-180.) < EPS ) THEN        ! Ligand and surface normal are antiparallel
      WRITE(line,1100) x(1)+0.1,x(2)+0.1,x(3)+0.1, surf_normal
   ELSE
      WRITE(line,1100) x,surf_normal
   ENDIF
   laenge = LEN_TRIM(line)
   CALL vprod(line, laenge)
   sym_uvw(:) = res_para(1:3)
   CALL trans (sym_uvw, cr_gten, sym_hkl, 3)
   CALL symm_setup
!  CALL symm_show
   CALL symm_op_single
ENDIF
!
! 2nd step: Rotate ligand molecule to place chelate axis normal to the surface
! Determine rotation axis by vector product ligand_chelate x surface_normal
! define rotation operation
!  Find the vector from ligand atoms neig(1) and neig(2), the ligand_chelate axis
!  Needs new calculation, as molecule was rotated in step 1
x(1) = cr_pos(1,n2) - cr_pos(1,n1)
x(2) = cr_pos(2,n2) - cr_pos(2,n1)
x(3) = cr_pos(3,n2) - cr_pos(3,n1)
sym_angle = do_bang(lspace, surf_normal, vnull, x) - 90.0
IF(ABS(sym_angle) > EPS ) THEN                ! Rotate if not zero degrees
   sym_orig(:)    = cr_pos(:,ia)              ! Define origin in surface atom
   sym_trans(:)   = 0.0                       ! No translation needed
   sym_sel_atom   = .true.                    ! Select atoms
   sym_new        = .false.                   ! No new types
   sym_power      =  1                        ! Just need one operation
   sym_type       = .true.                    ! Proper rotation
   sym_mode       = .false.                   ! Move atom to new position
   sym_orig_mol   = .false.                   ! Origin at crystal
   sym_power_mult =.false.                    ! No multiple copies
   sym_sel_atom   = .true.                    ! Select atoms not molecules
   sym_start      =  cr_natoms - mole_natoms + 1 ! set range of atoms numbers
   sym_end        =  cr_natoms
   IF(ABS(sym_angle-180.) < EPS ) THEN        ! Ligand and surface normal are antiparallel
      WRITE(line,1100) x(1)+0.1,x(2)+0.1,x(3)+0.1, surf_normal
   ELSE
      WRITE(line,1100) x,surf_normal
   ENDIF
   laenge = LEN_TRIM(line)
   CALL vprod(line, laenge)
   sym_uvw(:) = res_para(1:3)
   CALL trans (sym_uvw, cr_gten, sym_hkl, 3)
   CALL symm_setup
!  CALL symm_show
   CALL symm_op_single
ENDIF
!
IF(mole_axis(0)==2) THEN    ! Rotate upright, if two atoms are given
   CALL rotate_projected(n1, n2, nold, mole_axis, surf_normal)
ENDIF
!
!     Tilt molecule by user request
!
origin(1:3) = cr_pos(1:3, n1)
CALL deco_tilt(origin, tilt, tilt_hkl, tilt_atom, tilt_is_atom, &
               tilt_is_auto,                                    &
               surf_normal, mole_natoms, n1, n2)
!
m_type_new = m_type_old  + temp_id
!
CALL molecularize_numbers(nold+1,cr_natoms, m_type_new, r_m_biso, r_m_clin, r_m_cqua)
!
flagsurf: DO j=1,20
   IF(mole_surfnew(j)>0) THEN
      im = nold + mole_surfnew(j)
      cr_prop  (im) = IBSET (cr_prop (im), PROP_SURFACE_EXT)
      cr_surf(:,im) = cr_surf(:, ia)          ! Copy surface vector from anchor
   ELSE
      EXIT flagsurf
   ENDIF
ENDDO flagsurf
!
IF(nold<cr_natoms) THEN                          ! We did insert a molecule
   chem_period(:) = .false.                         ! We inserted atoms, turn off periodic boundaries
   chem_quick     = .false.                         ! turn of quick search
   cr_prop (ia) = IBCLR (cr_prop (ia), PROP_DECO_ANCHOR)  ! UNFLAG THIS ATOM AS SURFACE ANCHOR
   cr_prop (ia) = IBCLR (cr_prop (ia), PROP_SURFACE_EXT)   ! Anchor is no longer at a surface
ENDIF
!
1000 FORMAT(a4,4(2x,',',F12.6))
1100 FORMAT(6(G15.6E3,', '),'ddd')
!
END SUBROUTINE deco_place_chelate
!
!*****7*****************************************************************
!
   SUBROUTINE deco_place_multi(temp_id, ia, &
                            mole_axis, mole_surfnew,           &
                            m_type_old, mole_name, mole_natoms, &
                            mole_nscat, mole_atom_name, &
                            mole_dw, r_m_biso, r_m_clin, r_m_cqua,         &
                            nanch, anchor, &
                            nsites,                                               &
                            neig, dist,ncon, istart, iend, lrestrict, nhkl, rhkl, &
                            tilt, tilt_hkl, tilt_atom, tilt_is_atom, &
                            tilt_is_auto)
!
!  Places a molecule that has multiple bonds to the surface.
!  The first bond should be the one that carries multiple connections to
!  the surface. 
!  The first molecule atom is in a uniquely specified position.
!  If further bonds are specified, the molecule is rotated around this
!  first anchor point to fulfill as best as possible the further conditions.
!
USE crystal_mod
   USE atom_env_mod
   USE chem_mod
USE do_find_mod
   USE metric_mod
   USE modify_mod
   USE molecule_func_mod
   USE point_grp
USE prop_para_mod
USE read_internal_mod
USE surface_func_mod
   USE symm_menu
   USE symm_mod
   USE symm_sup_mod
   USE trafo_mod
!
   USE param_mod
   USE wink_mod
!
   IMPLICIT NONE
!
   INTEGER,                 INTENT(IN) :: temp_id         ! The definition to be used
   INTEGER,                 INTENT(IN) :: ia              ! Surface atom number
   INTEGER, DIMENSION(0:2), INTENT(IN) :: mole_axis       ! Atoms that define molecule axis
   INTEGER, DIMENSION(1:20),INTENT(IN) :: mole_surfnew    ! Atoms that will be flagged as surface atoms
   INTEGER,                 INTENT(IN) :: m_type_old      ! molecule types previous to all placements
   CHARACTER (LEN=1024),    INTENT(IN) :: mole_name       ! Molecule file name
   INTEGER,                 INTENT(IN) :: mole_natoms     ! Number of atoms in the molecule
   INTEGER,                 INTENT(IN) :: mole_nscat      ! Number of atoms in molecule
   CHARACTER (LEN=4   ), DIMENSION(0:mole_nscat),   INTENT(IN) :: mole_atom_name ! Atom names in the molecule
   REAL                , DIMENSION(0:mole_nscat),   INTENT(IN) :: mole_dw        ! ADPs       in the molecule
REAL,                    INTENT(IN) :: r_m_biso        ! Molecular Biso
REAL,                    INTENT(IN) :: r_m_clin        ! Molecular linear correction
REAL,                    INTENT(IN) :: r_m_cqua        ! Molecular quadratic correction
   INTEGER,                 INTENT(IN) :: nanch           ! Number of anchor types
   INTEGER, DIMENSION(1:2,0:nanch), INTENT(IN) :: anchor          ! Surface atom type
!  INTEGER, DIMENSION(0:4), INTENT(IN) :: surf            ! Surface atom type
   INTEGER, DIMENSION(1:2), INTENT(IN) :: nsites          ! Number of surface anchor positions per bond
   INTEGER, DIMENSION(1:2), INTENT(IN) :: neig            ! Connected to this neighbor in mole
   REAL   , DIMENSION(1:2), INTENT(IN) :: dist            ! distance to ligand molecule
   INTEGER,                 INTENT(IN) :: ncon            ! Number of defined bonds
   INTEGER,                 INTENT(IN) :: istart          ! First atom for surface determination
   INTEGER,                 INTENT(IN) :: iend            ! Last  atom for surface determination
   LOGICAL,                 INTENT(IN) :: lrestrict       ! Restriction to a surface type T/F
   INTEGER,                 INTENT(IN) :: nhkl            ! Number of faces for the restriction
   INTEGER, DIMENSION(3,nhkl), INTENT(IN) :: rhkl            ! actual faces for the restriction
REAL                      , INTENT(IN) :: tilt            ! Molecule tilt angle
REAL,    DIMENSION(3)     , INTENT(IN) :: tilt_hkl        ! Molecule tilt plane by this normal
INTEGER, DIMENSION(4)     , INTENT(IN) :: tilt_atom       ! Molecule tilt plane defined by these atoms
LOGICAL                   , INTENT(IN) :: tilt_is_atom    ! Plane defined by atoms
LOGICAL                   , INTENT(IN) :: tilt_is_auto    ! Plane defined by atoms
!
   INTEGER, PARAMETER                      :: MINPARA = 2
   INTEGER                                 :: MAXW = MINPARA
   REAL                , DIMENSION(1:MAX(MINPARA,nanch  )) :: werte
!
   INTEGER                              :: test_nhkl
   INTEGER, DIMENSION(:,:), ALLOCATABLE :: test_hkl
   CHARACTER (LEN=1024)                    :: line
   INTEGER   :: surf_char ! Surface character, plane, edge, corner, ...
   INTEGER, DIMENSION(0:4)             :: surface         ! Surface atom type
   INTEGER                                 :: ianz
   INTEGER                                 :: i,j, l,im, laenge
   INTEGER                                 :: iprop
   INTEGER                                 :: itype
   INTEGER, DIMENSION(0:3)                 :: isurface       ! Atom surface
   INTEGER                                 :: n_atoms_orig   ! Number of atoms prior to insertion
   INTEGER                                 :: n1,n2          ! number of mol neighbours after insertion
   INTEGER                                 :: a1,a2          ! number of mol axis atoms after rotations
   INTEGER                                 :: success        ! Everything went fine
   INTEGER  , DIMENSION(1:2)               :: is_good        ! Atom numbers of 2nd and 3rd surface sites
   INTEGER  , DIMENSION(1:4)               :: hkl
   LOGICAL  , DIMENSION(1:3)               :: fp
   LOGICAL                                 :: fq
   LOGICAL, PARAMETER :: lspace = .true.
   REAL                                    :: rmin, radius, b_l
   REAL                                    :: arg
   REAL                                    :: alpha, beta, v_l
   REAL     , DIMENSION(1:3)               :: x, bridge, base, origin, posit, v, w, u
   REAL     , DIMENSION(1:3)               :: shift, v1, v2, v3
   REAL     , DIMENSION(1:3)               :: surf_normal
   REAL     , DIMENSION(1:3)               :: vnull
   INTEGER, DIMENSION(3,6) :: surface_normal ! Set of local normals (:,1) is main normal
   INTEGER, DIMENSION(3)   :: surf_kante     ! Edge vector if not a plane
   INTEGER, DIMENSION(6)   :: surf_weight    ! Best normal has heighest weight
   INTEGER             :: nold         ! atom number previous to current molecule
   INTEGER             :: m_type_new   ! new molecule types 
   INTEGER             :: in_mole,in_moleatom
!  INTEGER             :: i_m_mole
!   INTEGER             :: i_m_type
!  INTEGER             :: i_m_char
!  CHARACTER (LEN=200) :: c_m_file
!  REAL                :: r_m_fuzzy
!  REAL                :: r_m_dens
!
!
maxw     = MAX(MINPARA,nanch)
vnull(:) = 0.00
success = -1
n2 = 1
nold = cr_natoms                           ! Remember original atom number
!
!  Determine surface character, if growth is restricted check if we're at proper surface
!
hkl(4) = 0
CALL surface_character(ia, istart, iend, surf_char, surface_normal, surf_kante, surf_weight, .TRUE.)
IF(surf_char == 0) RETURN
surf_normal(1:3) = FLOAT(surface_normal(:,1))
hkl(1:3) = surface_normal(:,1)
IF(lrestrict) THEN
   test_nhkl =    nhkl
   ALLOCATE(test_hkl(3,test_nhkl))
   test_hkl(1:3,1:test_nhkl) =    rhkl(1:3,1:test_nhkl)
   IF(.NOT.point_test(hkl, test_hkl, test_nhkl, .TRUE.) ) THEN
      RETURN
   ENDIF
   DEALLOCATE(test_hkl)
ENDIF
!
!  Load the molecule into the crystal structure
!
n_atoms_orig = cr_natoms                      ! Number of atoms prior to insertion
nold         = cr_natoms                      ! Remember original atom number
origin(:)    = 0.0                            ! initially place at 0,0,0
sym_latom(:) = .false.                        ! Initially deselect all atomtypes
!
insert: DO im=1,mole_natoms                   ! Insert all atoms
   CALL struc_read_one_atom_internal(mole_name, im, posit, itype, iprop, isurface,in_mole,in_moleatom)
   posit(:) = posit(:) + origin(:)
   WRITE(line, 1000) mole_atom_name(itype), posit, mole_dw(itype)
   laenge = 60
   CALL do_ins(line, laenge)
   cr_prop (cr_natoms) = ibset (cr_prop (cr_natoms), PROP_LIGAND)
   cr_prop (cr_natoms) = ibset (cr_prop (cr_natoms), PROP_SURFACE_EXT)
   cr_surf(:,cr_natoms) = 0
   CALL check_symm
   sym_latom(cr_iscat(cr_natoms)) = .true.    ! Select atom type for rotation
ENDDO insert
!
chem_period(:) = .false.                         ! We inserted atoms, turn off periodic boundaries
chem_quick     = .false.                         ! turn of quick search
!
!  Find the other surface atoms involved in this bond
j = anchor(1,0)
surface(1:j) = anchor(1,1:j)
surface(0) = j
n1 = n_atoms_orig +     neig(1)  !Absolute number for 1st neighbor in molecule
CALL deco_find_anchor(nsites(1), surface(0), surface, dist(1), ia,  &
                      surf_normal, posit, is_good, base, success)
IF(success/=0) THEN ! DID not find a suitable anchor, flag error
   GOTO 9999
ENDIF
!
!   Move molecule to anchor position
!
shift (:) = posit(:) - cr_pos(:,n1)
DO i=n_atoms_orig+1,cr_natoms
   cr_pos(:,i) = cr_pos(:,i) + shift(:)
ENDDO
success = 1                                     ! Start with error flag
!
IF(ncon == 2) THEN                            ! We have the second connection
   l = 2
   surface = 0
   n2 = n_atoms_orig +     neig(2)   
   bridge(1) = cr_pos(1,n2) - cr_pos(1,n1)
   bridge(2) = cr_pos(2,n2) - cr_pos(2,n1)
   bridge(3) = cr_pos(3,n2) - cr_pos(3,n1)
   b_l      = sqrt(skalpro(bridge, bridge, cr_gten))
   rmin     = MAX( 0.0, b_l - dist(1) - dist(2) )          ! Minimum distance between surface atoms
   radius   = b_l + dist(1) + dist(2)                      ! Maximum distance between surface atoms
   ianz     = 1
   j            = anchor(2,0)
   surface(1:j) = anchor(2,1:j)
   surface(0) = j
   werte(1:j) = surface(1:j)
   ianz = j
   x(:)  = cr_pos(:,ia)                               ! Search around 1.st surface atom
   fp(:) = chem_period (:)
   fq    = chem_quick
   CALL do_find_env (ianz, werte, maxw, x, rmin, radius, fq, fp)  ! Find all neighbors
   IF(atom_env(0) >= 1 ) THEN                                     ! We need at least one neighbor
     j = 0
     check_prop: DO i=1,atom_env(0)                               ! Check properties 
        IF(IBITS(cr_prop(atom_env(i)),PROP_SURFACE_EXT,1).eq.1 .and.        &  ! real Atom is near surface
           IBITS(cr_prop(atom_env(i)),PROP_OUTSIDE    ,1).eq.0       ) THEN    ! real Atom is near surface
            j = i                                                 ! Will use this neighbor
            EXIT check_prop                                       ! Found first good neighbor
         ENDIF
      ENDDO check_prop
      IF(j==0) THEN                                      ! No suitable neighbor, quietly leave
         GOTO 9999
      ENDIF
   ELSE
      GOTO 9999
   ENDIF  ! 
ENDIF
!
bridge(1) = cr_pos(1,n2) - cr_pos(1,n1)    ! Current vector Mole 1 to mole 2
bridge(2) = cr_pos(2,n2) - cr_pos(2,n1)
bridge(3) = cr_pos(3,n2) - cr_pos(3,n1)
b_l      = sqrt(skalpro(bridge, bridge, cr_gten))
u(:) = base(:)      - cr_pos(:,ia)         ! Vector from 1st surface to point below 1st mole
v(:) = cr_pos(:,n1) - cr_pos(:,ia)         ! Vector from 1st surface to             1st mole
v_l      = sqrt(skalpro(v, v, cr_gten))    ! Bond length 1st surface to 1st mole
WRITE(line,1100) u, v                      ! Do vector product 
laenge = LEN_TRIM(line)
CALL vprod(line, laenge)
sym_uvw(:) =  res_para(1:3)
arg = (dist(1)**2 + dist(2)**2 - b_l**2)/(2.*dist(1)*dist(2))
IF(ABS(arg)> 1.0) GOTO 9999                ! No suitable solution
sym_angle  = acos( (dist(1)**2 + dist(2)**2 - b_l**2)/(2.*dist(1)*dist(2)))/REAL(rad)
v(:) = v(:) *dist(2)/v_l                        ! Scale vector 1st surface to 1st mole to distance2
!
sym_orig(:)    = 0.0                       ! Define origin at 0,0,0
sym_trans(:)   = 0.0                       ! No translation needed
sym_sel_atom   = .true.                    ! Select atoms
sym_new        = .false.                   ! No new types
sym_power      =  1                        ! Just need one operation
sym_type       = .true.                    ! Proper rotation
sym_mode       = .false.                   ! Move atom to new position
sym_orig_mol   = .false.                   ! Origin at crystal
sym_power_mult =.false.                    ! No multiple copies
sym_sel_atom   = .true.                    ! Select atoms not molecules
CALL trans (sym_uvw, cr_gten, sym_hkl, 3)  ! Make reciprocal space axis
CALL symm_setup
CALL symm_ca_single (v, .true., .false.)
IF(ier_num /= 0) GOTO 9999                 ! Rotation is erroneous, |Axis| = 0 or similar 
posit(:) = cr_pos(:,ia) + res_para(1:3)    ! Add rotated vector to 1st surface
!
! next step rotate molecule for 2nd mole to fall onto target posit
u(:) = posit(:) - cr_pos(:,n1)             ! Vector from 1st mole to target
WRITE(line,1100) bridge, u                 ! Do vector product (1st to 2nd mole) x (1st mole to target)
laenge = LEN_TRIM(line)
CALL vprod(line, laenge)
sym_uvw(:) =  res_para(1:3)
sym_angle  = do_bang(lspace, bridge, vnull, u) ! Angle (1st to 2nd mole) and (1st mole to target)
sym_orig(:) = cr_pos(:,n1)                 ! Set origin in 1st mole
sym_start  =  n_atoms_orig + 1             ! set range of atoms numbers
sym_end    =  cr_natoms
CALL trans (sym_uvw, cr_gten, sym_hkl, 3)  ! Make reciprocal space axis
CALL symm_setup
CALL symm_op_single                        ! Perform the operation
IF(ier_num /= 0) GOTO 9999                 ! Rotation is erroneous, |Axis| = 0 or similar 
!
IF(mole_axis(0)==2) THEN    ! Rotate upright, if two atoms are given
!   Rotate molecule up to straighten molecule axis out
   sym_uvw(:) = cr_pos(:,n2) - cr_pos(:,n1)       ! Rotation axis
   a1 = n_atoms_orig + mole_axis(1)  ! dc_temp_axis(1)
   a2 = n_atoms_orig + mole_axis(2)  ! dc_temp_axis(2)
   v1(:)      = cr_pos(:,a2) - cr_pos(:,a1)      ! Current molecule axis
   WRITE(line,1200) v1, sym_uvw                  ! First project molecule axis into 
   laenge = LEN_TRIM(line)                       !   plane normal to the 
   CALL do_proj(line, laenge)                    !   vector between connected molecule atoms
   v3(:) = res_para(4:6)                         ! This is the projection
   WRITE(line,1100) sym_uvw, surf_normal         ! Find normal to plane defined by
   laenge = LEN_TRIM(line)                       !   vector between connected molecule atoms
   CALL vprod(line, laenge)                      !   and surface normal
   w(:) = res_para(1:3)                          ! Need to project (projected) mol axis into plane normal to w
   WRITE(line,1200) v3, w                        ! Prepare projection
   laenge = LEN_TRIM(line)
   CALL do_proj(line, laenge)                    ! Project axis into plane 
   v2(:) = res_para(4:6)                         ! This is the projection
   alpha      = do_bang(lspace, surf_normal, vnull, v2)   ! Calculate angle normal and projection
   WRITE(line,1100) v3,v2                        ! Do vector product (mol_axis) x (projection)
   laenge = LEN_TRIM(line)
   CALL vprod(line, laenge)
   u(:) =  res_para(1:3)
   beta = do_bang(lspace, sym_uvw, vnull, u)     ! Calculate angle (rot-axis) to vector product 
   IF(beta < 90) THEN                            ! Need to invert rotation axis
      IF(alpha < 90) THEN
         sym_angle  = do_bang(lspace, v3, vnull, v2)   ! Calculate rotation angle = < (v1,v2)
      ELSE
         sym_angle  =-180.+do_bang(lspace, v3, vnull, v2)   ! Calculate rotation angle = < (v1,v2)
      ENDIF
   ELSE
      sym_uvw(:) = -sym_uvw(:)
      IF(alpha < 90) THEN
         sym_angle  = do_bang(lspace, v3, vnull, v2)   ! Calculate rotation angle = < (v1,v2)
      ELSE
         sym_angle =-180.+do_bang(lspace, v3, vnull, v2)
      ENDIF
   ENDIF
   sym_orig(:) = cr_pos(:,n1)                    ! Origin in 1st bonded atom
   CALL trans (sym_uvw, cr_gten, sym_hkl, 3)     ! Make reciprocal space axis
   CALL symm_setup
   CALL symm_op_single                           ! Perform the operation
ENDIF
!
!     Tilt molecule by user request
!
origin(1:3) = cr_pos(1:3, n1)
CALL deco_tilt(origin, tilt, tilt_hkl, tilt_atom, tilt_is_atom, &
               tilt_is_auto,                                    &
               surf_normal, mole_natoms, n1, n2)
!
cr_prop (ia) = IBCLR (cr_prop (ia), PROP_DECO_ANCHOR)  ! UNFLAG THIS ATOM AS SURFACE ANCHOR
cr_prop (n1) = IBCLR (cr_prop (n1), PROP_DECO_ANCHOR)  ! UNFLAG THIS ATOM AS SURFACE ANCHOR
cr_prop (n2) = IBCLR (cr_prop (n2), PROP_DECO_ANCHOR)  ! UNFLAG THIS ATOM AS SURFACE ANCHOR
cr_prop (ia) = IBCLR (cr_prop (ia), PROP_SURFACE_EXT)   ! Anchor is no longer at a surface
cr_prop (n1) = IBCLR (cr_prop (n1), PROP_SURFACE_EXT)   ! Anchor is no longer at a surface
cr_prop (n2) = IBCLR (cr_prop (n2), PROP_SURFACE_EXT)   ! Anchor is no longer at a surface
DO j=1, 2
  n1 = is_good(j)
  cr_prop (n1) = IBCLR (cr_prop (n1), PROP_DECO_ANCHOR)  ! UNFLAG THIS ATOM AS SURFACE ANCHOR
  cr_prop (n1) = IBCLR (cr_prop (n1), PROP_SURFACE_EXT)   ! Anchor is no longer at a surface
ENDDO
flagsurf: DO j=1,20
   IF(mole_surfnew(j)>0) THEN
      im = nold + mole_surfnew(j)
      cr_prop  (im) = IBSET (cr_prop (im), PROP_SURFACE_EXT)
      cr_surf(:,im) = cr_surf(:, ia)          ! Copy surface vector from anchor
   ELSE
      EXIT flagsurf
   ENDIF
ENDDO flagsurf
!
m_type_new = m_type_old  + temp_id
CALL molecularize_numbers(nold+1,cr_natoms, m_type_new, r_m_biso, r_m_clin, r_m_cqua)
success = 0                                   ! Clear error flag
!
9999 CONTINUE                                 ! Jump here from errors to ensure dealloc
   IF(success /=0) THEN                       ! An error occurred, reset crystal
      cr_natoms = n_atoms_orig
   ENDIF
!
1000 FORMAT(a4,4(2x,',',F12.6))
1100 FORMAT(6(G15.6E3,', '),'ddd')
1200 FORMAT(6(G15.6E3,', '),'dddd')
!
   END SUBROUTINE deco_place_multi
!
!******************************************************************************
!
   SUBROUTINE deco_place_acceptor(temp_id, ia, &
                            mole_axis, mole_surfnew,          &
                            m_type_old, mole_name, mole_natoms, &
                            mole_nscat, mole_atom_name, &
                            mole_dw, r_m_biso, r_m_clin, r_m_cqua,         &
                            neig, dist, temp_secnd, istart, iend, lrestrict, nhkl, rhkl, dha_angle)
!
!  Place molecules via a hydrogen bond onto the surface acceptor atom
!
USE crystal_mod
   USE chem_mod
   USE metric_mod
   USE modify_mod
   USE molecule_func_mod
   USE point_grp
USE prop_para_mod
   USE read_internal_mod
USE surface_func_mod
   USE surface_mod
   USE symm_menu
   USE symm_mod
   USE symm_sup_mod
   USE trafo_mod
!
   USE param_mod
!
   IMPLICIT NONE
!
   INTEGER,                 INTENT(IN) :: temp_id         ! The definition to be used
   INTEGER,                 INTENT(IN) :: ia              ! Surface atom number
   INTEGER,                 INTENT(IN) :: m_type_old      ! molecule types previous to all placements
   INTEGER, DIMENSION(0:2), INTENT(IN) :: mole_axis       ! Atoms that define molecule axis
   INTEGER, DIMENSION(1:20),INTENT(IN) :: mole_surfnew    ! Atoms that will be flagged as surface atoms
   CHARACTER (LEN=1024),    INTENT(IN) :: mole_name       ! Molecule file name
   INTEGER,                 INTENT(IN) :: mole_natoms     ! Number of atoms in the molecule
   INTEGER,                 INTENT(IN) :: mole_nscat      ! Number of atoms in molecule
   CHARACTER (LEN=4   ), DIMENSION(0:mole_nscat),   INTENT(IN) :: mole_atom_name ! Atom names in the molecule
   REAL                , DIMENSION(0:mole_nscat),   INTENT(IN) :: mole_dw        ! ADPs       in the molecule
REAL,                    INTENT(IN) :: r_m_biso        ! Molecular Biso
REAL,                    INTENT(IN) :: r_m_clin        ! Molecular linear correction
REAL,                    INTENT(IN) :: r_m_cqua        ! Molecular quadratic correction
   INTEGER,                 INTENT(IN) :: neig            ! Connected to this neighbor in mole
   REAL   ,                 INTENT(IN) :: dist            ! distance to ligand molecule
   INTEGER,                 INTENT(IN) :: temp_secnd      ! Second neighbor atom in molecule
   INTEGER,                 INTENT(IN) :: istart          ! First atom for surface determination
   INTEGER,                 INTENT(IN) :: iend            ! Last  atom for surface determination
   LOGICAL,                 INTENT(IN) :: lrestrict       ! Restriction to a surface type T/F
   INTEGER,                 INTENT(IN) :: nhkl            ! Number of faces for the restriction
   INTEGER, DIMENSION(3,nhkl), INTENT(IN) :: rhkl         ! actual faces for the restriction
   REAL   ,                 INTENT(IN) :: dha_angle       ! hydrogen-Bond angle in Hydrogen atom
!
!  REAL, PARAMETER         :: DIST_A_H     = 1.920   ! Average Acceptor Hydrogon distance
!  REAL, PARAMETER         :: SIGMA_A_H    = 0.001   ! Sigma for Acceptor Hydrogon distance
!  REAL, PARAMETER         :: ANGLE_A_H_D  = 170.0   ! Average Angle in Hydrogen bond
   REAL, PARAMETER         :: SIGMA_A_H_D  =   0.0001! Sigma for Angle in Hydrogen bond
   REAL, DIMENSION(3), PARAMETER :: VNULL = (/ 0.0, 0.0, 0.0 /) 
   LOGICAL, PARAMETER      :: lspace=.TRUE.
   CHARACTER (LEN=1024)    :: line
   INTEGER                 ::    j, im, laenge  ! Dummy index
   INTEGER                 :: n1, n2         ! Atoms that define molecule axis
   INTEGER                 :: itype, iprop   ! Atom types, properties
   INTEGER, DIMENSION(0:3) :: isurface       ! Atom surface
   INTEGER                 :: nold           ! atom number previous to current molecule
   INTEGER                 :: surf_char      ! Surface character, plane, edge, corner, ...
   INTEGER, DIMENSION(3,6) :: surface_normal ! Set of local normals (:,1) is main normal
   INTEGER, DIMENSION(3)   :: surf_kante     ! Edge vector if not a plane
   INTEGER, DIMENSION(6)   :: surf_weight    ! Best normal has heighest weight
INTEGER                              :: test_nhkl
   INTEGER, DIMENSION(:,:), ALLOCATABLE :: test_hkl
   REAL                    :: normal_l       ! length of normal vector
   REAL                    :: hbond          ! actual hydrogen bond A..H
   REAL                    :: angle          ! Temporary angle
   REAL                    :: solution_1, solution_2 ! Temporary angles
   REAL   , DIMENSION(1:3) :: surf_normal    ! Normal to work with
   REAL   , DIMENSION(3)   :: posit          ! Temporary atom position
   REAL   , DIMENSION(3)   :: origin         ! Temporary origin for symmetry operations
   REAL   , DIMENSION(3)   :: u,v,w,up,wp,p, prot  ! Temporary vectors

   INTEGER, DIMENSION(4) :: hkl
   INTEGER             :: m_type_new   ! new molecule types 
   INTEGER             :: in_mole,in_moleatom
!
   REAL :: gaslim, ran1
!
!  Determine surface character, if growth is restricted check if we're at proper surface
!
   hkl(4) = 0
   CALL surface_character(ia, istart, iend, surf_char, surface_normal, surf_kante, surf_weight, .TRUE.)
IF(surf_char == 0) RETURN
!
   surf_normal(1:3) = FLOAT(surface_normal(:,1))
   hkl(1:3)         =       surface_normal(:,1)
IF(lrestrict) THEN
   test_nhkl =    nhkl
   ALLOCATE(test_hkl(3,test_nhkl))
   test_hkl(1:3,1:test_nhkl) =    rhkl(1:3,1:test_nhkl)
      IF(.NOT.point_test(hkl, test_hkl, test_nhkl, .TRUE.) ) THEN
         RETURN
      ENDIF
      DEALLOCATE(test_hkl)
ENDIF
!
!
!  Accept as surfaces any outside facing surface and all internal "flat" surfaces
IF(surf_char /=0 .AND. surf_char > -SURF_EDGE) THEN    ! Surface atoms only
   nold = cr_natoms
!
!  Insert molecule atoms into crystal at correct origin in initial orientation
   normal_l = SQRT (skalpro (surf_normal, surf_normal, cr_gten))
   hbond = dist
!  hbond = DIST_A_H + gaslim(SIGMA_A_H, 2.0)     ! Make distributed distances
   origin(1)  = cr_pos(1,ia) + surf_normal(1)/normal_l*hbond  ! Origin is shifted
   origin(2)  = cr_pos(2,ia) + surf_normal(2)/normal_l*hbond  ! by hbond away from
   origin(3)  = cr_pos(3,ia) + surf_normal(3)/normal_l*hbond  ! surface atom
   sym_latom(:) = .false.                        ! Initially deselect all atomtypes
!
   atoms: DO im=1,mole_natoms                    ! Load all atoms from the molecule
      CALL struc_read_one_atom_internal(mole_name, im, posit, itype, iprop, isurface,in_mole,in_moleatom)
      posit(:) = posit(:) + origin(:)
      WRITE(line, 1000) mole_atom_name(itype), posit, mole_dw(itype)
      laenge = 60
      CALL do_ins(line, laenge)                  ! Insert into crystal
      cr_prop  (cr_natoms) = IBSET (cr_prop (cr_natoms), PROP_LIGAND)
      cr_prop  (cr_natoms) = IBCLR (cr_prop (cr_natoms), PROP_SURFACE_EXT)
      cr_surf(:,cr_natoms) = 0
      CALL check_symm
      sym_latom(cr_iscat(cr_natoms)) = .true.    ! Select atom type for rotation
   ENDDO atoms
!
!           Rotate ligand to achieve Hydrogen bond angle: ANGLE_A_H_D
!
            u(:) = cr_pos(:, nold+neig) - cr_pos(:,ia)          ! Vector Hydrogen to surface atom
            v(:) = cr_pos(:, nold+neig) - cr_pos(:,nold+temp_secnd)  ! Vector Hydrogen to Neighbor in ligand
            WRITE(line,1100) u,v                          ! Do vector product (mol_axis) x (projection)
            laenge = LEN_TRIM(line)
            CALL vprod(line, laenge)
            sym_uvw(:)     =  res_para(1:3)
!           sym_angle      = ANGLE_A_H_D - do_bang(lspace, u, VNULL, v) &
            sym_angle      = dha_angle   - do_bang(lspace, u, VNULL, v) &
                                         + gaslim(SIGMA_A_H_D, 2.0)! 
            sym_orig(:)    = cr_pos(:,nold+neig)          ! Rotate in Hydrogen
            sym_trans(:)   = 0.0                          ! No translation needed
            sym_sel_atom   = .true.                       ! Select atoms
            sym_new        = .false.                      ! No new types
            sym_power      =  1                           ! Just need one operation
            sym_type       = .true.                       ! Proper rotation
            sym_mode       = .false.                      ! Move atom to new position
            sym_orig_mol   = .false.                      ! Origin at crystal
            sym_power_mult =.false.                       ! No multiple copies
            sym_sel_atom   = .true.                       ! Select atoms not molecules
            sym_start      =  cr_natoms - mole_natoms + 1 ! set range of atoms numbers
            sym_end        =  cr_natoms
            sym_uvw(:) = res_para(1:3)
            CALL trans (sym_uvw, cr_gten, sym_hkl, 3)     ! Make reciprocal space axis
            CALL symm_setup                               ! Symmetry setup defines matrix
!           CALL symm_show                                ! Show only in debug
            CALL symm_op_single                           ! Perform the operation
!
!           Rotate ligand around Acceptor-Hydrogon bond by random degree
!
            sym_uvw(:)     = u(:)
            sym_angle      = 360.0*ran1(0)
            CALL trans (sym_uvw, cr_gten, sym_hkl, 3)     ! Make reciprocal space axis
            CALL symm_setup                               ! Symmetry setup defines matrix
!           CALL symm_show                                ! Show only in debug
            CALL symm_op_single                           ! Perform the operation
!           Group molecule
            m_type_new = m_type_old  + temp_id 
!
!  Rotate molecule around the H-O vector to align molecule axis along the surface normal
!
IF(mole_axis(0)==2) THEN    ! Rotate upright, if two atoms are given
!   Rotate molecule up to straighten molecule axis out
   n1 = nold + mole_axis(1)                      ! First atom on molecule axis
   n2 = nold + mole_axis(2)                      ! Second atom on molecule axis
   v(:) = cr_pos(:, nold+neig) - cr_pos(:,nold+temp_secnd)  ! Vector Hydrogen to Neighbor in ligand
   w(:) = cr_pos(:, nold+neig) - cr_pos(:,ia)               ! Vector anchor to Hydrogen
   u(:) = cr_pos(:,n2) - cr_pos(:,n1)                       ! Current molecule axis
   WRITE(line,1200) w, v                         ! Prepare projection normal onto H=>O vector
   laenge = LEN_TRIM(line)
   CALL do_proj(line, laenge)                    ! Project normal into plane normal to H=> O vector
   wp(:) =  res_para(4:6)                        ! Normal projected into plane normal to H=>O vector
   WRITE(line,1200) u, v                         ! Prepare projection axis   onto H=>O vector
   laenge = LEN_TRIM(line)
   CALL do_proj(line, laenge)                    ! Project axis into plane normal to H=> O vector
   up(:) =  res_para(4:6)                        ! Axis projected into plane normal to H=>O vector
   angle = do_bang(lspace, wp, VNULL, up)        ! angle between the two projections
!
!  Two rotations will move the second axis atom into the plane defined by surface normal and
!  H==>O vector, a rotation by alpha and (180-alpha) We have to test which of these creates the 
!  smaller angle between the normal and the resulting axis.
   sym_uvw(:) = cr_pos(:, nold+neig) - cr_pos(:,nold+temp_secnd)  ! Rotation axis: Vector Hydrogen to Neighbor in ligand
   sym_angle     = angle
   sym_orig(:)    = cr_pos(:,nold+temp_secnd)    ! Rotate in Oxygen
   sym_trans(:)   = 0.0                          ! No translation needed
   sym_sel_atom   = .true.                       ! Select atoms
   sym_new        = .false.                      ! No new types
   sym_power      =  1                           ! Just need one operation
   sym_type       = .true.                       ! Proper rotation
   sym_mode       = .false.                      ! Move atom to new position
   sym_orig_mol   = .false.                      ! Origin at crystal
   sym_power_mult = .false.                      ! No multiple copies
   sym_sel_atom   = .true.                       ! Select atoms not molecules
   sym_start      =  cr_natoms - mole_natoms + 1 ! set range of atoms numbers
   sym_end        =  cr_natoms
   CALL trans (sym_uvw, cr_gten, sym_hkl, 3)     ! Make reciprocal space axis
   CALL symm_setup                               ! Symmetry setup defines matrix
!
   p(:) = cr_pos(:,n2)                           ! absolute position of axis atom number two
   CALL symm_ca_single(p,lspace,.FALSE.)
   prot(:) = res_para(1:3) - cr_pos(:,n1)        ! Resulting rotated axis solution 1
   solution_1 = do_bang(lspace, w, VNULL, prot)
!
   sym_angle = 180.-angle
   CALL symm_setup                               ! Symmetry setup defines matrix
   p(:) = cr_pos(:,n2)                           ! absolute position of axis atom number two
   CALL symm_ca_single(p,lspace,.FALSE.)
   prot(:) = res_para(1:3) - cr_pos(:,n1)        ! Resulting rotated axis solution 2
   solution_2 = do_bang(lspace, w, VNULL, prot)
   IF( solution_1 <= solution_2) THEN
      sym_angle = angle
   ELSE
      sym_angle = 180.0 - angle
   ENDIF
   CALL symm_setup                               ! Symmetry setup defines matrix
   CALL symm_op_single                           ! Perform the operation on the whole molecule
   
ENDIF
!
   CALL molecularize_numbers(nold+1,cr_natoms, m_type_new, r_m_biso, r_m_clin, r_m_cqua)
   flagsurf: DO j=1,20
      IF(mole_surfnew(j)>0) THEN
         im = nold + mole_surfnew(j)
         cr_prop  (im) = IBSET (cr_prop (im), PROP_SURFACE_EXT)
         cr_surf(:,im) = cr_surf(:, ia)          ! Copy surface vector from anchor
      ELSE
         EXIT flagsurf
      ENDIF
   ENDDO flagsurf
ENDIF
chem_period(:) = .false.                         ! We inserted atoms, turn off periodic boundaries
chem_quick     = .false.                         ! turn of quick search
cr_prop (ia) = IBCLR (cr_prop (ia), PROP_DECO_ANCHOR)  ! UNFLAG THIS ATOM AS SURFACE ANCHOR
cr_prop (ia) = IBCLR (cr_prop (ia), PROP_SURFACE_EXT)  ! Anchor is no longer at a surface
!
1000 FORMAT(a4,4(2x,',',F12.6))
1100 FORMAT(6(G15.6E3,', '),'ddd')
1200 FORMAT(6(G15.6E3,', '),'dddd')
!
END SUBROUTINE deco_place_acceptor
!
!******************************************************************************
!
SUBROUTINE deco_place_donor(temp_id, ia, &
                            mole_axis, mole_surfnew,          &
                            m_type_old, mole_name, mole_natoms, &
                            mole_nscat, mole_atom_name, &
                            mole_dw, r_m_biso, r_m_clin, r_m_cqua,         &
                            neig, dist, istart, iend, lrestrict, nhkl, rhkl, dha_angle)
!
!  Place molecules via a hydrogn bond onto the surface donor atom
!
USE crystal_mod
USE atom_env_mod
USE chem_mod
USE do_find_mod
USE metric_mod
USE modify_mod
USE molecule_func_mod
USE point_grp
USE prop_para_mod
USE read_internal_mod
USE surface_func_mod
USE surface_mod
USE symm_menu
USE symm_mod
USE symm_sup_mod
USE trafo_mod
!
use molecule_mod
!
USE param_mod
USE random_mod
!
IMPLICIT NONE
!
INTEGER,                 INTENT(IN) :: temp_id         ! The definition to be used
INTEGER,                 INTENT(IN) :: ia              ! Surface atom number
INTEGER,                 INTENT(IN) :: m_type_old      ! molecule types previous to all placements
INTEGER, DIMENSION(0:2), INTENT(IN) :: mole_axis       ! Atoms that define molecule axis
INTEGER, DIMENSION(1:20),INTENT(IN) :: mole_surfnew    ! Atoms that will be flagged as surface atoms
CHARACTER (LEN=1024),    INTENT(IN) :: mole_name       ! Molecule file name
INTEGER,                 INTENT(IN) :: mole_natoms     ! Number of atoms in the molecule
   INTEGER,                 INTENT(IN) :: mole_nscat      ! Number of atoms in molecule
   CHARACTER (LEN=4   ), DIMENSION(0:mole_nscat),   INTENT(IN) :: mole_atom_name ! Atom names in the molecule
   REAL                , DIMENSION(0:mole_nscat),   INTENT(IN) :: mole_dw        ! ADPs       in the molecule
REAL,                    INTENT(IN) :: r_m_biso        ! Molecular Biso
REAL,                    INTENT(IN) :: r_m_clin        ! Molecular linear correction
REAL,                    INTENT(IN) :: r_m_cqua        ! Molecular quadratic correction
INTEGER,                 INTENT(IN) :: neig            ! Connected to this neighbor in mole
REAL   ,                 INTENT(IN) :: dist            ! distance to ligand molecule
INTEGER,                 INTENT(IN) :: istart          ! First atom for surface determination
INTEGER,                 INTENT(IN) :: iend            ! Last  atom for surface determination
LOGICAL,                 INTENT(IN) :: lrestrict       ! Restriction to a surface type T/F
INTEGER,                 INTENT(IN) :: nhkl            ! Number of faces for the restriction
INTEGER, DIMENSION(3,nhkl), INTENT(IN) :: rhkl         ! actual faces for the restriction
REAL   ,                 INTENT(IN) :: dha_angle       ! hydrogen-Bond angle in Hydrogen atom
!
INTEGER, PARAMETER      :: MAXW = 2
!  REAL, PARAMETER         :: DIST_A_H     = 1.920   ! Average Acceptor Hydrogon distance
!  REAL, PARAMETER         :: SIGMA_A_H    = 0.001   ! Sigma for Acceptor Hydrogon distance
REAL, PARAMETER         :: EPS = 1.0E-7
!REAL, PARAMETER         :: ANGLE_A_H_D  = 170.0   ! Average Angle in Hydrogen bond
!   REAL, PARAMETER         :: SIGMA_A_H_D  =   0.0001! Sigma for Angle in Hydrogen bond
   REAL, DIMENSION(3), PARAMETER :: VNULL = (/ 0.0, 0.0, 0.0 /) 
   LOGICAL, PARAMETER      :: lspace=.TRUE.
   CHARACTER (LEN=1024)    :: line
   INTEGER                 :: ianz
   INTEGER                 ::    j, im, laenge  ! Dummy index
   INTEGER                 :: n1, n2         ! Atom that define teh molecule axis
   INTEGER                 :: itype, iprop   ! Atom types, properties
   INTEGER, DIMENSION(0:3) :: isurface       ! Atom surface
   INTEGER                 :: nold           ! atom number previous to current molecule
   INTEGER                 :: surf_char      ! Surface character, plane, edge, corner, ...
   INTEGER                 :: success        ! Failure or success ?
   INTEGER, DIMENSION(3,6) :: surface_normal ! Set of local normals (:,1) is main normal
   INTEGER, DIMENSION(3)   :: surf_kante     ! Edge vector if not a plane
   INTEGER, DIMENSION(6)   :: surf_weight    ! Best normal has heighest weight
   INTEGER                              :: test_nhkl
   INTEGER, DIMENSION(:,:), ALLOCATABLE :: test_hkl
   LOGICAL  , DIMENSION(3) :: fp
   LOGICAL                 :: fq
   REAL                    :: normal_l       ! length of normal vector
   REAL                    :: hbond          ! actual hydrogen bond A..H
   REAL                    :: rmin, rmax
   REAL                    :: angle          !Temporary angle
   REAL                    :: d2             !Temporary distance
   REAL   , DIMENSION(1:3) :: surf_normal    ! Normal to work with
   REAL   , DIMENSION(3)   :: posit          ! Temporary atom position
   REAL   , DIMENSION(3)   :: origin         ! Temporary origin for symmetry operations
   REAL   , DIMENSION(3)   :: u,v,x, w       ! Temporary vectors
   REAL   , DIMENSION(1:MAXW) :: werte

   INTEGER, DIMENSION(4) :: hkl
   INTEGER             :: m_type_new   ! new molecule types 
   INTEGER             :: in_mole,in_moleatom
!
   REAL :: ran1
!
success = -1
fp(:) = .FALSE.
fq    = .FALSE.
!
!  Determine surface character, if growth is restricted check if we're at proper surface
!
hkl(4) = 0
CALL surface_character(ia, istart, iend, surf_char, surface_normal, surf_kante, surf_weight, .TRUE.)
IF(surf_char == 0) RETURN
!
surf_normal(1:3) = FLOAT(surface_normal(:,1))
hkl(1:3)         =       surface_normal(:,1)
IF(lrestrict) THEN
   test_nhkl =    nhkl
   ALLOCATE(test_hkl(3,test_nhkl))
   test_hkl(1:3,1:test_nhkl) =    rhkl(1:3,1:test_nhkl)
   IF(.NOT.point_test(hkl, test_hkl, test_nhkl, .TRUE.) ) THEN
      RETURN
   ENDIF
   DEALLOCATE(test_hkl)
ENDIF
!
!
nold = cr_natoms
!  Accept as surfaces any outside facing surface and all internal "flat" surfaces
IF(surf_char /=0 .AND. surf_char > -SURF_EDGE) THEN    ! Surface atoms only
!
!  Insert molecule atoms into crystal at correct origin in initial orientation
   normal_l = SQRT (skalpro (surf_normal, surf_normal, cr_gten))
!  hbond = DIST_A_H + gaslim(SIGMA_A_H, 2.0)     ! Make distributed distances
   hbond = dist
   origin(1)  = cr_pos(1,ia) + surf_normal(1)/normal_l*hbond  ! Origin is shifted
   origin(2)  = cr_pos(2,ia) + surf_normal(2)/normal_l*hbond  ! by hbond away from
   origin(3)  = cr_pos(3,ia) + surf_normal(3)/normal_l*hbond  ! surface Hydrogen atom
   sym_latom(:) = .false.                        ! Initially deselect all atomtypes
!
   atoms: DO im=1,mole_natoms                    ! Load all atoms from the molecule
      CALL struc_read_one_atom_internal(mole_name, im, posit, itype, iprop, isurface,in_mole,in_moleatom)
      posit(:) = posit(:) + origin(:)
      WRITE(line, 1000) mole_atom_name(itype), posit, mole_dw(itype)
      laenge = 60
      CALL do_ins(line, laenge)                  ! Insert into crystal
      cr_prop  (cr_natoms) = IBSET (cr_prop (cr_natoms), PROP_LIGAND)
      cr_prop  (cr_natoms) = IBCLR (cr_prop (cr_natoms), PROP_SURFACE_EXT)
      cr_surf(:,cr_natoms) = 0
      CALL check_symm
      sym_latom(cr_iscat(cr_natoms)) = .true.    ! Select atom type for rotation
   ENDDO atoms
!
!  Rotate ligand to achieve Hydrogen bond angle: ANGLE_A_H_D
!
!  Find the donor neighbor to the Hydrogen
   x(:) = cr_pos(:,ia)      ! Copy Hydrogen atom position for formal reasons
   rmin = 0.1
   rmax = 1.2               ! Small limit to find donor atom only
   werte(:) = -1
   CALL do_find_env (ianz, werte, MAXW, x, rmin, rmax, fq, fp)
!
!  Did we find a proper covalent neighbor at less than 1.2 A?
!
   IF(atom_env(0)==0) GOTO 9999                  !Failure, leave
   u(:) = cr_pos(:, atom_env((1))) - cr_pos(:,ia)     ! Vector Hydrogen to Donor
   v(:) = cr_pos(:, nold+neig) - cr_pos(:,ia)          ! Vector Hydrogen to Neighbor in ligand
   WRITE(line,1100) u,v                          ! Do vector product (mol_axis) x (projection)
   laenge = LEN_TRIM(line)
   CALL vprod(line, laenge)
   d2 = res_para(1)**2 + res_para(2)**2 + res_para(3)**2
   angle = do_bang(lspace, u, VNULL, v) ! + gaslim(SIGMA_A_H_D, 2.0)
   IF(d2 < EPS    ) THEN   !D==>H and H==>A are parallel, seek alternative solution
      ier_num = 0
      ier_typ = 0
      v(1) = v(1) + ran1(idum)
      v(2) = v(2) - ran1(idum)
      WRITE(line,1100) u,v                          ! Do vector product (mol_axis) x (projection)
      laenge = LEN_TRIM(line)
      CALL vprod(line, laenge)
      d2 = res_para(1)**2 + res_para(2)**2 + res_para(3)**2
      IF(d2 < EPS    ) THEN   !D==>H and H==>A are still parallel, silently give up, should not happen ?
         GOTO 9999
      ENDIF
   ENDIF
   sym_uvw(:)     =  res_para(1:3)
!  sym_angle      = ANGLE_A_H_D - angle 
   sym_angle      = dha_angle   - angle 
   sym_orig(:)    = cr_pos(:,ia)                 ! Rotate in Hydrogen
   sym_trans(:)   = 0.0                          ! No translation needed
   sym_sel_atom   = .true.                       ! Select atoms
   sym_new        = .false.                      ! No new types
   sym_power      =  1                           ! Just need one operation
   sym_type       = .true.                       ! Proper rotation
   sym_mode       = .false.                      ! Move atom to new position
   sym_orig_mol   = .false.                      ! Origin at crystal
   sym_power_mult =.false.                       ! No multiple copies
   sym_sel_atom   = .true.                       ! Select atoms not molecules
   sym_start      =  cr_natoms - mole_natoms + 1 ! set range of atoms numbers
   sym_end        =  cr_natoms
   sym_uvw(:) = res_para(1:3)
   CALL trans (sym_uvw, cr_gten, sym_hkl, 3)     ! Make reciprocal space axis
   CALL symm_setup                               ! Symmetry setup defines matrix
!  CALL symm_show                                ! Show only in debug
   CALL symm_op_single                           ! Perform the operation
!
!  Rotate ligand around Donor-Hydrogen bond by random degree
!
   sym_uvw(:)     = u(:)
   sym_orig(:)    = cr_pos(:,nold+neig)          ! Rotate in Acceptor
   sym_angle      = 360.0*ran1(0)
   CALL trans (sym_uvw, cr_gten, sym_hkl, 3)     ! Make reciprocal space axis
   CALL symm_setup                               ! Symmetry setup defines matrix
!  CALL symm_show                                ! Show only in debug
   CALL symm_op_single                           ! Perform the operation
!
!  Rotate molecule around the H-O vector to align molecule axis along the surface normal
!
IF(mole_axis(0)==2) THEN    ! Rotate upright, if two atoms are given
!   Rotate molecule up to straighten molecule axis out
   n1 = nold + mole_axis(1)                      ! First atom on molecule axis
   n2 = nold + mole_axis(2)                      ! Second atom on molecule axis
!  v(:) = cr_pos(:, nold+neig) - cr_pos(:,ia)    ! Vector Hydrogen to Neighbor in ligand
   w(:) = surf_normal(:)                         ! Vector anchor to Hydrogen
   u(:) = cr_pos(:,n2) - cr_pos(:,n1)            ! Current molecule axis
   angle = do_bang(lspace,u, VNULL, w)
   IF(angle /= 0.0) THEN
      WRITE(line,1100) u,w                          ! Do vector product (mol_axis) x (projection)
      laenge = LEN_TRIM(line)
      CALL vprod(line, laenge)
      sym_uvw(:) = res_para(1:3)
   sym_angle     = angle
   sym_orig(:)    = cr_pos(:,nold+neig)          ! Rotate in Oxygen
   sym_trans(:)   = 0.0                          ! No translation needed
   sym_sel_atom   = .true.                       ! Select atoms
   sym_new        = .false.                      ! No new types
   sym_power      =  1                           ! Just need one operation
   sym_type       = .true.                       ! Proper rotation
   sym_mode       = .false.                      ! Move atom to new position
   sym_orig_mol   = .false.                      ! Origin at crystal
   sym_power_mult = .false.                      ! No multiple copies
   sym_sel_atom   = .true.                       ! Select atoms not molecules
   sym_start      =  cr_natoms - mole_natoms + 1 ! set range of atoms numbers
   sym_end        =  cr_natoms
   CALL trans (sym_uvw, cr_gten, sym_hkl, 3)     ! Make reciprocal space axis
   CALL symm_setup                               ! Symmetry setup defines matrix
!
   CALL symm_op_single                           ! Perform the operation on the whole molecule
   
   ENDIF
ENDIF
!
!  Group molecule
   m_type_new = m_type_old  + temp_id
   CALL molecularize_numbers(nold+1,cr_natoms, m_type_new, r_m_biso, r_m_clin, r_m_cqua)
   flagsurf: DO j=1,20
      IF(mole_surfnew(j)>0) THEN
         im = nold + mole_surfnew(j)
         cr_prop  (im) = IBSET (cr_prop (im), PROP_SURFACE_EXT)
         cr_surf(:,im) = cr_surf(:, ia)          ! Copy surface vector from anchor
      ELSE
         EXIT flagsurf
      ENDIF
   ENDDO flagsurf
!
ENDIF
chem_period(:) = .false.                         ! We inserted atoms, turn off periodic boundaries
chem_quick     = .false.                         ! turn of quick search
cr_prop (ia) = IBCLR (cr_prop (ia), PROP_DECO_ANCHOR)  ! UNFLAG THIS ATOM AS SURFACE ANCHOR
cr_prop (ia) = IBCLR (cr_prop (ia), PROP_SURFACE_EXT)  ! Anchor is no longer at a surface
success = 0
!
9999 CONTINUE                                 ! Jump here from errors to ensure dealloc
   IF(success /=0) THEN                       ! An error occurred, reset crystal
      cr_natoms = nold
   ENDIF
!
!
1000 FORMAT(a4,4(2x,',',F12.6))
1100 FORMAT(6(G15.6E3,', '),'ddd')
!
END SUBROUTINE deco_place_donor
!
!*****7*****************************************************************
!
   SUBROUTINE deco_find_anchor(MAXAT,MAXTYPE, surface, distance, ia, &
                               normal, posit, is_good, &
                               base, ierror)
!-                                                                      
!  Find a common point around MAXAT atom types in surface
!  ia is the first surface atom number
!
USE crystal_mod
   USE atom_env_mod
   USE chem_mod
USE do_find_mod
   USE metric_mod
   USE modify_mod
USE prop_para_mod
!
USE errlist_mod
   USE param_mod
   USE wink_mod
!
   IMPLICIT NONE
   INTEGER                    , INTENT(IN)  :: MAXAT
   INTEGER                    , INTENT(IN)  :: MAXTYPE
   INTEGER, DIMENSION(0:MAXTYPE), INTENT(IN)  :: surface
   REAL                       , INTENT(IN)  :: distance
   INTEGER                    , INTENT(IN)  :: ia
   REAL   , DIMENSION(1:3)    , INTENT(IN)  :: normal
   REAL   , DIMENSION(1:3)    , INTENT(OUT) :: posit
   REAL   , DIMENSION(1:3)    , INTENT(OUT) :: base
   INTEGER, DIMENSION(1:2)    , INTENT(OUT) :: is_good
   INTEGER                    , INTENT(INOUT) :: ierror
!
   LOGICAL, PARAMETER         :: lspace = .true. 
!
   CHARACTER (LEN=1024)                    :: line
   INTEGER :: ianz, i, j, k, l, kgood, lgood
   INTEGER :: good1, good2, good3
   INTEGER :: laenge
   INTEGER, DIMENSION(0:6,2:MAXAT) :: neig
   LOGICAL, DIMENSION(1:3) :: fp
   LOGICAL                 :: fq
   REAL                    :: rmin, radius
   REAL   , DIMENSION(1:MAXTYPE) :: werte
   REAL   , DIMENSION(1:3)    :: x, u,v,w, e1,e2,e3
   REAL                    :: u_l, v_l, w_l    ! length of vectors in triangle
   REAL                    :: av, sig, av_min, sig_min ! average length  and sigma
   REAL                    :: tx,ty, tz        ! Cartesion coordinates of target position
   REAL                    :: g2x, g2y         ! Cartesion coordinates of atom 3
   REAL                    :: arg
   REAL     , DIMENSION(1:3)               :: vnull
!
   vnull(:) = 0.00
!
   good2     = ia
   good3     = ia
   ierror    = 0
   neig(:,:) = 0
   rmin     = 0.0                  ! Minimum distance between surface atoms
   radius   = 2.0 * distance       ! Maximum distance between surface atoms
   ianz     = 1
   ianz     = MAXTYPE
   x(:)     = cr_pos(:, ia)        ! Seach around atom ia
   fp (:)   = chem_period (:)
   fq       = chem_quick
   find: DO l=2, MAXAT
      werte(1) = surface(l)           ! Find this atom type
      werte(1:MAXTYPE) = surface(1:MAXTYPE)           ! Find this atom type
      CALL do_find_env (ianz, werte, MAXTYPE, x, rmin, radius, fq, fp)  ! Find all neighbors
      IF(atom_env(0) >= 1 ) THEN                                     ! We need at least one neighbor
        j = 0
        check_prop: DO i=1,atom_env(0)                               ! Check properties 
           IF(IBITS(cr_prop(atom_env(i)),PROP_SURFACE_EXT,1).eq.1 .and.        &  ! real Atom is near surface
              IBITS(cr_prop(atom_env(i)),PROP_OUTSIDE    ,1).eq.0       ) THEN    ! real Atom is near surface
               j = j +1                                              ! Will use this neighbor
               neig(j,l) = atom_env(i)
               neig(0,l) = j
               IF(j==6) EXIT check_prop                              ! Found first six good neighbors
            ENDIF
         ENDDO check_prop
         IF(j==0) THEN                                      ! No suitable neighbor, quietly leave
            ier_num = -1118
            ier_typ = ER_APPL
            RETURN
         ENDIF
      ELSE
         ier_num = -1118
         ier_typ = ER_APPL
         RETURN
      ENDIF  ! 
   ENDDO find
!
!
   IF(maxat==3) THEN                        ! Find a suitable triangle
       av_min = 1.E10
      sig_min = 1.E10
      lgood   = 0
      kgood   = 0
      first: DO l=1,neig(0,2)
         u(:) = cr_pos(:,neig(l,2)) - x(:)         ! Difference vector central to neighbor
         u_l  = SQRT(skalpro(u,u,cr_gten))         ! Calculate length ia to neighbor
         second: DO k = 1,neig(0,3)                         ! Loop over all secon neighbors
            IF(neig(l,2) /= neig(k,3) ) THEN          ! Exclude identical neighbors
               v(:) = cr_pos(:,neig(k,3)) - x(:)         ! Difference vector central to neighbor
               v_l  = SQRT(skalpro(v,v,cr_gten))         ! Calculate length ia to neighbor
               w(:) = cr_pos(:,neig(k,3)) - cr_pos(:,neig(l,2)) ! Difference vector neigh to neighbor
               w_l  = SQRT(skalpro(w,w,cr_gten))         ! Calculate length neigh to neighbor
               av   = (u_l+v_l+w_l)/3.
               arg  =     ((u_l-av)**2+(v_l-av)**2+(w_l-av)**2)
               IF(arg<0) THEN                     ! Negative arg, no solution return with error status
                  ierror = -1 
                  RETURN
               ENDIF
               sig  = SQRT(arg                                )/3
               IF(sig < sig_min) THEN
                  lgood = l
                  kgood = k
                  sig_min = sig
                   av_min = av
               ENDIF
            ENDIF
         ENDDO second
      ENDDO first
   good1 = ia                             ! Atom is is 0,0,0 corner in cartesian space
   good2 = neig(lgood,2)                  ! Atom is along x-axis in cartesian space
   good3 = neig(kgood,3)                  ! Atoms definex x-y plane in cartesian space
   u (:) = cr_pos(:,good2) - cr_pos(:,good1)  ! Cartesian x-axis
   u_l  = SQRT(skalpro(u ,u ,cr_gten))
   IF(u_l<=0) THEN                        ! Negative arg, no solution return with error status
      ierror = -2 
      RETURN
   ENDIF
   e1(:) = u (:) / u_l                    ! Normalize to 1 angstroem
   v (:) = cr_pos(:,good3) - cr_pos(:,good1)  ! Temporary vector
   WRITE(line,1100) e1,v                         ! Do vector product (e1) x (atom good 3)
   laenge = LEN_TRIM(line)
   CALL vprod(line, laenge)
   e3(:) =  res_para(1:3)                 ! Result is cartesian z-axis
   v_l  = SQRT(skalpro(e3,e3,cr_gten))
   IF(v_l<=0) THEN                        ! Negative arg, no solution return with error status
      ierror = -2 
      RETURN
   ENDIF
   IF(do_bang(lspace, e3, vnull, normal) <= 90) THEN
      e3(:) =  e3(:) / v_l                ! Normalize to 1 angstroem
   is_good(1) = good2
   ELSE
      e3(:) = -e3(:) / v_l                ! invert and Normalize to 1 angstroem
   ENDIF
   WRITE(line,1100) e3,e1                 ! Do vector product (e3) x (e1 )
   laenge = LEN_TRIM(line)
   CALL vprod(line, laenge)
   e2(:) =  res_para(1:3)                 ! Result is cartesian y-axis
   v_l  = SQRT(skalpro(e2,e2,cr_gten))    ! cartesian x-coordinate of atom 2
   IF(v_l<=0) THEN                        ! Negative arg, no solution return with error status
      ierror = -2 
      RETURN
   ENDIF
   e2(:) = e2(:) / v_l                    ! Normalize to 1 angstroem
   g2x =     (skalpro(v,e1,cr_gten))      ! cartesian x-coordinate of atom 3
   g2y =     (skalpro(v,e2,cr_gten))      ! cartesian y-coordinate of atom 3
!  Calculate target coordinates from trilateration
   tx = 0.5 * u_l
   ty = 0.5 * (g2x**2+g2y**2)/g2y - g2x/g2y*tx
   arg =     (distance**2-tx**2 - ty**2)
   IF(arg<=0) THEN                     ! Negative arg, no solution return with error status
      ierror = -1 
      RETURN
   ENDIF
   tz = sqrt(arg)
   posit(:) = cr_pos(:,good1) + tx*e1(:) + ty*e2(:) + tz*e3(:)
   
   base(:) = (cr_pos(:,neig(lgood,2))+ cr_pos(:,neig(kgood,3)))*0.5
   ENDIF
   is_good(1) = good2
   is_good(2) = good3
!
1100 FORMAT(6(G15.6E3,', '),'ddd')
!
   END SUBROUTINE deco_find_anchor
!
!###############################################################################
!
SUBROUTINE deco_tilt(origin, tilt, tilt_hkl, tilt_atom, tilt_is_atom, &
                     tilt_is_auto,                                    &
                     surf_normal, mole_natoms, n1, n2)
!
USE crystal_mod
USE fit_mod
USE symm_menu
USE symm_mod
USE symm_sup_mod
USE metric_mod
USE trafo_mod
!
USE param_mod
!
IMPLICIT NONE
!
REAL,    DIMENSION(3), INTENT(IN) :: origin          ! Molecule origin
REAL                 , INTENT(IN) :: tilt            ! Molecule tilt angle
REAL,    DIMENSION(3), INTENT(IN) :: tilt_hkl        ! Molecule tilt plane by this normal normal
INTEGER, DIMENSION(4), INTENT(IN) :: tilt_atom       ! Molecule tilt plane defined by these atoms
LOGICAL              , INTENT(IN) :: tilt_is_atom    ! Plane defined by atoms
LOGICAL              , INTENT(IN) :: tilt_is_auto    ! Plane defined automatically
REAL,    DIMENSION(3), INTENT(IN) :: surf_normal     ! Local surface normal
INTEGER              , INTENT(IN) :: mole_natoms     ! number of atoms in molecule
INTEGER              , INTENT(IN) :: n1              ! number of atom that defines rotation axis
INTEGER              , INTENT(IN) :: n2              ! number of atom that defines rotation axis
!
CHARACTER(LEN=1024)   :: line
INTEGER               :: nold      ! Original number of atoms
INTEGER               :: laenge
INTEGER               :: i
INTEGER, DIMENSION(:), ALLOCATABLE :: list
REAL   , DIMENSION(3) :: u, v, hkl ! Dummy vectors
REAL                  :: dist
!
IF(ABS(tilt)>0.001) THEN
nold = cr_natoms - mole_natoms
!
   IF(n1==0 .AND. n2==0) THEN
      IF(tilt_is_auto) THEN                      ! Determine molecular plane automatically
         ALLOCATE(list(1:mole_natoms))
         DO i=1,mole_natoms
            list(i) = nold + i
         ENDDO
         CALL dis_fit_plane(mole_natoms, list, hkl, dist)
         DEALLOCATE(list)
      ELSE
      IF(tilt_is_atom) THEN
         u(:) = cr_pos(:,nold+tilt_atom(2)) - cr_pos(:,nold+tilt_atom(1))
         v(:) = cr_pos(:,nold+tilt_atom(4)) - cr_pos(:,nold+tilt_atom(3))
         WRITE(line,1100) u, v                  ! Do vector product (e3) x (e1 )
         WRITE(*   ,1100) u, v                  ! Do vector product (e3) x (e1 )
         laenge = LEN_TRIM(line)
         CALL vprod(line, laenge)
         hkl(:) =  res_para(1:3)                 ! Result is cartesian y-axis
      ELSE
         hkl(:) = tilt_hkl(:)
      ENDIF
      ENDIF
   ELSE
      hkl(:) = cr_pos(:,n2) - cr_pos(:,n1)
   ENDIF
!
   sym_uvw(:)     = hkl(1:3)
   sym_angle      = tilt
   sym_orig(:)    = origin(:)                    ! Rotate in origin
   sym_trans(:)   = 0.0                          ! No translation needed
   sym_sel_atom   = .true.                       ! Select atoms
   sym_new        = .false.                      ! No new types
   sym_power      =  1                           ! Just need one operation
   sym_type       = .true.                       ! Proper rotation
   sym_mode       = .false.                      ! Move atom to new position
   sym_orig_mol   = .false.                      ! Origin at crystal
   sym_power_mult = .false.                       ! No multiple copies
   sym_sel_atom   = .true.                       ! Select atoms not molecules
   sym_start      =  cr_natoms - mole_natoms + 1 ! set range of atoms numbers
   sym_end        =  cr_natoms
   CALL trans (sym_uvw, cr_gten, sym_hkl, 3)     ! Make reciprocal space axis
   CALL symm_setup                               ! Symmetry setup defines matrix
!  CALL symm_show                                ! Show only in debug
   CALL symm_op_single                           ! Perform the operation
ENDIF
!
1100 FORMAT(6(G15.6E3,', '),'ddd')
!
END SUBROUTINE deco_tilt
!
!###############################################################################
!
   SUBROUTINE check_symm
!
USE crystal_mod
   USE discus_allocate_appl_mod
   USE molecule_mod
   USE symm_mod
!
   IMPLICIT NONE
   INTEGER   :: nscat     ! dummy for allocation
!
   IF( cr_nscat > SYM_MAXSCAT .or. mole_num_type > SYM_MAXSCAT) THEN
      nscat = max ( cr_nscat+10, mole_num_type)
      CALL alloc_symmetry ( nscat )
      IF ( ier_num < 0 ) THEN
         RETURN
      ENDIF
   ENDIF
   END SUBROUTINE check_symm
!
!###############################################################################
!
SUBROUTINE rotate_projected(n1, n2, n_atoms_orig, mole_axis, surf_normal)
!
USE crystal_mod
USE metric_mod
USE symm_menu
USE symm_mod
USE symm_sup_mod
USE trafo_mod
USE param_mod
!
IMPLICIT NONE
!
INTEGER                , INTENT(IN) :: n1 ! Neighbor 1 defines rotation axis
INTEGER                , INTENT(IN) :: n2 ! Neighbor 2 defines rotation axis
INTEGER                , INTENT(IN) :: n_atoms_orig ! Number of atoms prior to deco
INTEGER, DIMENSION(0:2), INTENT(IN) :: mole_axis    ! Axis to straighten up
REAL   , DIMENSION(1:3), INTENT(IN) :: surf_normal  ! Surface normal at anchor atom
!
LOGICAL, PARAMETER  :: lspace = .TRUE.
CHARACTER(LEN=1024) :: line
INTEGER             :: a1, a2          ! absolute atom numbers for rotation axix
INTEGER             :: laenge
REAL                :: alpha, beta
REAL, DIMENSION(3)  :: v1, v2, v3, u, w, vnull   ! Dummy vectors
!
vnull(:) = 0.0
!
!   Rotate molecule up to straighten molecule axis out
sym_uvw(:) = cr_pos(:,n2) - cr_pos(:,n1)      ! Rotation axis
a1 = n_atoms_orig + mole_axis(1)              ! Absolute number for axis atom 1
a2 = n_atoms_orig + mole_axis(2)              ! Absolute number for axis atom 2
v1(:)      = cr_pos(:,a2) - cr_pos(:,a1)      ! Current molecule axis
WRITE(line,1200) v1, sym_uvw                  ! First project molecule axis into 
laenge = LEN_TRIM(line)                       !   plane normal to the 
CALL do_proj(line, laenge)                    !   vector between connected molecule atoms
v3(:) = res_para(4:6)                         ! This is the projection
WRITE(line,1100) sym_uvw, surf_normal         ! Find normal to plane defined by
laenge = LEN_TRIM(line)                       !   vector between connected molecule atoms
CALL vprod(line, laenge)                      !   and surface normal
w(:) = res_para(1:3)                          ! Need to project (projected) mol axis into plane normal to w
WRITE(line,1200) v3, w                        ! Prepare projection
laenge = LEN_TRIM(line)
CALL do_proj(line, laenge)                    ! Project axis into plane 
v2(:) = res_para(4:6)                         ! This is the projection
alpha      = do_bang(lspace, surf_normal, vnull, v2)   ! Calculate angle normal and projection
WRITE(line,1100) v3,v2                        ! Do vector product (mol_axis) x (projection)
laenge = LEN_TRIM(line)
CALL vprod(line, laenge)
u(:) =  res_para(1:3)
beta = do_bang(lspace, sym_uvw, vnull, u)     ! Calculate angle (rot-axis) to vector product 
IF(beta < 90) THEN                            ! Need to invert rotation axis
   IF(alpha < 90) THEN
      sym_angle  = do_bang(lspace, v3, vnull, v2)   ! Calculate rotation angle = < (v1,v2)
   ELSE
      sym_angle  =-180.+do_bang(lspace, v3, vnull, v2)   ! Calculate rotation angle = < (v1,v2)
   ENDIF
ELSE
   sym_uvw(:) = -sym_uvw(:)
   IF(alpha < 90) THEN
      sym_angle  = do_bang(lspace, v3, vnull, v2)   ! Calculate rotation angle = < (v1,v2)
   ELSE
      sym_angle =-180.+do_bang(lspace, v3, vnull, v2)
   ENDIF
ENDIF
sym_orig(:) = cr_pos(:,n1)                 ! Origin in 1st bonded atom
sym_trans(:)   = 0.0                       ! No translation needed
sym_sel_atom   = .true.                    ! Select atoms
sym_new        = .false.                   ! No new types
sym_power      =  1                        ! Just need one operation
sym_type       = .true.                    ! Proper rotation
sym_mode       = .false.                   ! Move atom to new position
sym_orig_mol   = .false.                   ! Origin at crystal
sym_power_mult =.false.                    ! No multiple copies
sym_sel_atom   = .true.                    ! Select atoms not molecules
sym_start      =  n_atoms_orig + 1         ! set range of atoms numbers
sym_end        =  cr_natoms
!
CALL trans (sym_uvw, cr_gten, sym_hkl, 3)     ! Make reciprocal space axis
CALL symm_setup
CALL symm_op_single                           ! Perform the operation
!
1100 FORMAT(6(G15.6E3,', '),'ddd')
1200 FORMAT(6(G15.6E3,', '),'dddd')
!
END SUBROUTINE rotate_projected
!
!###############################################################################
!
SUBROUTINE rotate_directly(neig, n_atoms_orig, mole_axis, surf_normal)
!
USE crystal_mod
USE metric_mod
USE symm_menu
USE symm_mod
USE symm_sup_mod
USE trafo_mod
USE param_mod
use errlist_mod
!
IMPLICIT NONE
!
INTEGER                , INTENT(IN) :: neig ! Neighbor 1 defines rotation axis
INTEGER                , INTENT(IN) :: n_atoms_orig ! Number of atoms prior to deco
INTEGER, DIMENSION(0:2), INTENT(IN) :: mole_axis    ! Axis to straighten up
REAL   , DIMENSION(1:3), INTENT(IN) :: surf_normal  ! Surface normal at anchor atom
!
REAL   , PARAMETER  :: EPS = 1.0E-6
LOGICAL, PARAMETER  :: lspace = .TRUE.
CHARACTER(LEN=1024) :: line
INTEGER             :: n1              ! absolute atom numbers for neighbor
INTEGER             :: a1, a2          ! absolute atom numbers for rotation axix
INTEGER             :: laenge
REAL, DIMENSION(3)  :: axis_ligand, vnull   ! Dummy vectors
!
vnull(:) = 0.0
!
! define rotation operation
n1 = n_atoms_orig + neig                      ! Absolute number for neighboring atom
a1 = n_atoms_orig + mole_axis(1)              ! Absolute number for axis atom 1
a2 = n_atoms_orig + mole_axis(2)              ! Absolute number for axis atom 2
axis_ligand(:)      = cr_pos(:,a2) - cr_pos(:,a1)      ! Current molecule axis
      sym_angle      = do_bang(lspace, surf_normal, vnull, axis_ligand)
      IF(ABS(sym_angle) > EPS ) THEN                   ! Rotate if not zero degrees
         sym_orig(:)    = cr_pos(:,n1)                 ! Rotate in 1st neighbor atom
         sym_trans(:)   = 0.0                          ! No translation needed
         sym_sel_atom   = .true.                       ! Select atoms
         sym_new        = .false.                      ! No new types
         sym_power      =  1                           ! Just need one operation
         sym_type       = .true.                       ! Proper rotation
         sym_mode       = .false.                      ! Move atom to new position
         sym_orig_mol   = .false.                      ! Origin at crystal
         sym_power_mult =.false.                       ! No multiple copies
         sym_sel_atom   = .true.                       ! Select atoms not molecules
         sym_start      =  n_atoms_orig + 1            ! set range of atoms numbers
         sym_end        =  cr_natoms
         IF(ABS(sym_angle-180.) < EPS ) THEN           ! Ligand and surface normal are antiparallel
            WRITE(line,1100) axis_ligand(3)+0.1,axis_ligand(2)+0.01,axis_ligand(1)+0.001, surf_normal
         ELSE
            WRITE(line,1100) axis_ligand, surf_normal
         ENDIF
         laenge = LEN_TRIM(line)
         CALL vprod(line, laenge)                      ! Make rotation axis
         sym_uvw(:) = res_para(1:3)
         CALL trans (sym_uvw, cr_gten, sym_hkl, 3)     ! Make reciprocal space axis
         CALL symm_setup                               ! Symmetry setup defines matrix
!        CALL symm_show                                ! Show only in debug
         CALL symm_op_single                           ! Perform the operation
ENDIF
!
1100 FORMAT(6(G15.6E3,', '),'ddd')
!
END SUBROUTINE rotate_directly
!
!*******************************************************************************
!
END MODULE mole_surf_mod
