MODULE mole_surf_mod
!
USE crystal_mod
USE deco_mod
USE prop_para_mod
USE errlist_mod
!
PRIVATE
PUBLIC do_place_molecule
!
CONTAINS
   SUBROUTINE do_place_molecule
!
!   Main menu and procedures to place molecules onto a surface
!
   USE discus_allocate_appl_mod
   USE modify_mod
   USE modify_func_mod
!
   USE doact_mod
   USE learn_mod
   USE class_macro_internal
   USE prompt_mod
   IMPLICIT none
!
!
   CHARACTER (LEN=5)                       :: befehl! command on input line
   CHARACTER (LEN=50)                      :: prom  ! Menu prompt
   CHARACTER (LEN=1024)                    :: line  ! input line
   CHARACTER (LEN=1024)                    :: zeile ! remainder with parameters
   INTEGER                                 :: indxg ! location of "="
   INTEGER                                 :: lp    ! lengtz of zeile
   INTEGER laenge, lbef
   LOGICAL                                 :: ladd = .true.  ! condition add command is fine
   LOGICAL                                 :: lend  ! condition of EOF
   LOGICAL                                 :: lold =.true. ! condition of EOF
!
   INTEGER, EXTERNAL :: len_str
   LOGICAL, EXTERNAL :: str_comp
!
   lend = .false.
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
   IF(MAXSCAT > UBOUND(dc_latom,1)) THEN
      CALL alloc_deco(MAXSCAT)
   ENDIF
!
   main_loop: do
      CALL no_error
      prom = prompt (1:len_str (prompt) ) //'/deco'
      CALL get_cmd (line, laenge, befehl, lbef, zeile, lp, prom)
      no_err: IF (ier_num.eq.0) THEN
         no_com: IF (line /= ' '      .and. line(1:1) /= '#' .and.      &
             line /= char(13) .and. line(1:1) /= '!'        ) THEN
!                                                                       
!     ----search for "="                                                
!                                                                       
            indxg = index (line, '=')
            is_math: IF (indxg.ne.0.and.                            &
               .not. (str_comp (befehl, 'echo', 2, lbef, 4) ) .and. &
               .not. (str_comp (befehl, 'syst', 2, lbef, 4) ) .and. &
               .not. (str_comp (befehl, 'help', 2, lbef, 4)   .or.  &
                      str_comp (befehl, '?   ', 2, lbef, 4) ) ) THEN
!                                                                       
! ------evaluate an expression and assign the value to a variabble      
!                                                                       
               CALL do_math (line, indxg, laenge)
            ELSE
!                                                                       
!------ ----execute a macro file                                        
!                                                                       
              is_com: IF (befehl (1:1) .eq.'@') THEN 
                  IF (laenge.ge.2) THEN 
                     CALL file_kdo (line (2:laenge), laenge-1) 
                  ELSE 
                     ier_num = - 13 
                     ier_typ = ER_MAC 
                  ENDIF 
!                                                                       
!     ----continues a macro 'continue'                                  
!                                                                       
              ELSEIF (str_comp (befehl, 'continue', 5, lbef, 8) ) THEN 
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
                  CALL do_eval (zeile, lp) 
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
!     ----Original decorate commands                                      
!
!
              ELSEIF (str_comp (befehl, 'add', 3, lbef, 3)) THEN
                  IF(ladd) THEN
                     CALL deco_add (zeile, lp)
                     IF(ier_num == 0) ladd = .false.   ! exclude new add command
                  ELSE
                     ier_num = -128
                     ier_typ = ER_COMM
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
              ELSEIF (str_comp (befehl, 'reset', 3, lbef, 4)) THEN
                  CALL deco_reset
                  if(ier_num == 0) ladd = .true.   ! allow new add command
!
              ELSEIF (str_comp (befehl, 'run', 3, lbef, 3)) THEN
                  IF(ASSOCIATED(dc_def_temp)) THEN         ! We need at least one definition
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
                 IF(dc_temp_type /= DC_NONE) THEN
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
           ENDIF is_math   ! END IF BLOCK math equation or specific command
        ENDIF no_com       ! END IF BLOCK no comment
      ENDIF no_err         ! END IF BLOCK no error reading input
!
      IF (ier_num.ne.0) then 
         CALL errlist 
         IF (ier_sta.ne.ER_S_LIVE) then 
            IF (lmakro) then 
               CALL macro_close 
               prompt_status = PROMPT_ON 
            ENDIF 
            IF (lblock) then 
               RETURN 
            ENDIF 
            CALL no_error 
         ENDIF 
      ENDIF 
!
   ENDDO main_loop     ! END DO main loop of menu 
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
   USE modify_mod
!
   IMPLICIT NONE
!
!
   CHARACTER (LEN=*), INTENT(INOUT) :: zeile
   INTEGER          , INTENT(INOUT) :: lp
!
   INTEGER, PARAMETER   :: MAXW = 20
   CHARACTER (LEN=1024), DIMENSION(1:MAXW) :: cpara
   CHARACTER (LEN=1024), DIMENSION(1:2)    :: ccpara
   INTEGER             , DIMENSION(1:MAXW) :: lpara
   INTEGER             , DIMENSION(1:2)    :: llpara
   INTEGER              :: ianz, janz, success,i, ncon
   LOGICAL              :: lnew, lnew_mole = .false.
   REAL   , DIMENSION(1:MAXW) :: werte
   REAL   , DIMENSION(1:2   ) :: wwerte
!
   LOGICAL str_comp
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
   dc_temp_name  = cpara(1)
   dc_temp_lname = lpara(1)
   lnew          = .false.
   success       = 1
   IF ( str_comp(cpara(2),'axis',4,lpara(2),4) ) THEN
      IF ( ianz == 4 ) THEN
         CALL del_params (2, ianz, cpara, lpara, maxw)   ! delete first 2 params
         CALL ber_params (ianz, cpara, lpara, werte, maxw)
         IF(ier_num ==0) THEN
            dc_temp_axis(1) = NINT(werte(1))
            dc_temp_axis(2) = NINT(werte(2))
            CALL dc_find_def(dc_def_head,dc_def_temp, dc_temp_lname, dc_temp_name,dc_temp_id,lnew,success)
            IF(success==0) CALL dc_set_axis(dc_def_temp, dc_temp_axis)
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
         dc_temp_neig  = NINT(wwerte(1))
         dc_temp_dist  =      wwerte(2)
         cpara(ianz-1) = ' '
         cpara(ianz  ) = ' '
         lpara(ianz-1) = 1
         lpara(ianz  ) = 1
         janz = ianz - 2                                 ! ignore last two parameters
         CALL get_iscat (janz, cpara, lpara, werte, maxw, .false.)
         IF(ier_num /= 0) RETURN
         IF(dc_temp_type == DC_NORMAL .AND. janz /= 1) THEN
            ier_num = -132
            ier_typ = ER_APPL
            ier_msg(1) = 'For the NORMAL decoration we need one '
            ier_msg(2) = 'exactly surface atom type '
         ELSEIF(dc_temp_type == DC_BRIDGE .AND. janz /= 2) THEN
            ier_num = -132
            ier_typ = ER_APPL
            ier_msg(1) = 'For the BRIDGE decoration we need '
            ier_msg(2) = 'exactly two surface atom types'
         ELSEIF(dc_temp_type == DC_DOUBLE .AND. janz <  3) THEN
            ier_num = -132
            ier_typ = ER_APPL
            ier_msg(1) = 'For the MULTI decoration we need at least '
            ier_msg(2) = 'three surface atom types '
         ELSE ! SUCCESS
            DO i=1,janz
               dc_temp_surf(i)  = NINT(werte(i))  ! Surface atom type
            ENDDO
            dc_temp_surf(0) = janz
            dc_temp_id      = 0
            dc_def_temp => dc_def_head
!
            CALL dc_find_def(dc_def_head,dc_def_temp, dc_temp_lname, dc_temp_name,dc_temp_id,lnew,success)
            IF(success==0) THEN
               CALL dc_set_con(dc_def_temp%dc_def_con, dc_temp_surf, dc_temp_neig, dc_temp_dist)
               IF(ier_num == 0) THEN
                  CALL dc_inc_ncon(dc_def_temp, ncon)
!               DO i=1, dc_temp_surf(0)
!                  dc_latom(dc_temp_surf(i)) = .true.              ! Select this atom type
!               ENDDO
          
                  IF(ncon ==1 ) THEN
                     dc_latom(dc_temp_surf(1)) = .true.              ! Select first atom type
                  ENDIF
               ENDIF
            ELSE
               ier_num = -128
               ier_typ = ER_APPL
               ier_msg(1) = 'Check the definition type. '
               ier_msg(2) = 'Did you use the add command first ?'
            ENDIF
         ENDIF
      ELSE
         ier_num = -6
         ier_typ = ER_COMM
         ier_msg(1) = 'set bond command needs >= five parameters'
      ENDIF
   ELSEIF ( str_comp(cpara(2),'ligand',4,lpara(2),6) ) THEN
      IF ( ianz == 4 ) THEN
         dc_temp_file  = cpara(3)
         dc_temp_lfile = lpara(3)
         CALL del_params (3, ianz, cpara, lpara, maxw)   ! delete first 3 params
         IF(ier_num /= 0) RETURN
         CALL ber_params (ianz, cpara, lpara, werte, maxw)
         IF(ier_num /= 0) RETURN
         dc_temp_dens = werte(1)
         dc_temp_id    = 0
         dc_def_temp => dc_def_head
!
         CALL dc_find_def(dc_def_head,dc_def_temp, dc_temp_lname, dc_temp_name,dc_temp_id,lnew,success)
         IF(success==0) THEN
            CALL dc_set_file(dc_def_temp, dc_temp_lfile, dc_temp_file)
            CALL dc_set_dens(dc_def_temp, dc_temp_dens)
            lnew_mole = .true.
            search: DO i=1, dc_n_molecules
               IF(dc_input(i)(1:LEN_TRIM(dc_input(i))) == dc_temp_file(1:dc_temp_lfile)) THEN
                  lnew_mole = .false.
                  EXIT search
               ENDIF
            ENDDO search
            IF(lnew_mole) THEN
               dc_n_molecules = dc_n_molecules + 1
               dc_input(dc_n_molecules) = dc_temp_file(1:dc_temp_lfile)
            ENDIF
         ELSE
            ier_num = -128
            ier_typ = ER_APPL
            ier_msg(1) = 'Check the definition type. '
            ier_msg(2) = 'Did you use the add command first ?'
         ENDIF
      ELSE
         ier_num = -6
         ier_typ = ER_COMM
         ier_msg(1) = 'set ligand command needs four parameters'
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
   USE modify_mod
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
   LOGICAL              :: lnew
!
   LOGICAL str_comp
!
   CALL get_params(zeile, ianz, cpara, lpara, maxw, lp)
   IF ( ier_num /= 0 ) RETURN              ! Error reading parameters
!
   IF ( ianz == 2 ) THEN
      dc_temp_name  = cpara(1)
      dc_temp_lname = lpara(1)
      dc_temp_id    = 0
      success       = 1
      IF ( str_comp(cpara(2),'normal',4,lpara(2),6) ) THEN
         dc_temp_type = DC_NORMAL
      ELSEIF ( str_comp(cpara(2),'bridge',4,lpara(2),6) ) THEN
         dc_temp_type = DC_BRIDGE
      ELSEIF ( str_comp(cpara(2),'double',4,lpara(2),6) ) THEN
         dc_temp_type = DC_DOUBLE
      ELSEIF ( str_comp(cpara(2),'multiple',4,lpara(2),8) ) THEN
         dc_temp_type = DC_MULTIPLE
      ELSE
         ier_num = -6
         ier_typ = ER_COMM
         ier_msg(1) = 'Allowed decoration types are:'
         ier_msg(2) = 'normal, bridge, double, multi'
         NULLIFY(dc_def_temp) ! => dc_def_head
         RETURN
      ENDIF
      lnew          = .true.
      NULLIFY(dc_def_temp) ! => dc_def_head
      CALL dc_find_def(dc_def_head,dc_def_temp, dc_temp_lname, dc_temp_name,dc_temp_id,lnew,success)
      IF(success==0) THEN  ! set the type
         CALL dc_set_type(dc_def_temp, dc_temp_type)
      ELSE
         ier_num = -128
         ier_typ = ER_APPL
      ENDIF
   ELSE
      ier_num = -6
      ier_typ = ER_COMM
   ENDIF
!
   END SUBROUTINE deco_add
!
!
!######################################################################
!
   SUBROUTINE find_surface_character(ia, surf_char, surf_normal)
!
   USE atom_env_mod
   USE chem_mod
   USE metric_mod
   USE modify_mod
   USE tensors_mod
   USE prompt_mod
!
   IMPLICIT NONE
!
   INTEGER,INTENT(IN)  :: ia     ! dummy index
!  INTEGER,INTENT(IN)  :: istart ! dummy index
!  INTEGER,INTENT(IN)  :: iend   ! dummy index
   INTEGER,INTENT(OUT) :: surf_char   ! Surface character
   REAL   ,DIMENSION(1:3),INTENT(OUT) :: surf_normal   ! Surface normal

   INTEGER, PARAMETER :: MAXW = 1
   LOGICAL, PARAMETER :: LNEW = .false.
   REAL   , PARAMETER :: EPS  = 1E-7
!
   CHARACTER  (LEN=4), DIMENSION(1:MAXW) :: cpara
   INTEGER  ,          DIMENSION(1:MAXW) :: lpara
   REAL     ,          DIMENSION(1:MAXW) :: werte
!  INTEGER   :: ia ! dummy index
   INTEGER   :: ib ! dummy index
   INTEGER                 :: ianz  ! number of parameters in cpara
   INTEGER                 :: istat
   LOGICAL  , DIMENSION(3) :: fp
   LOGICAL                 :: fq
   REAL     , DIMENSION(3) :: x
   REAL      :: rmin
   REAL      :: radius
   REAL                      :: sigma
   REAL     , DIMENSION(3,3) :: direct
   REAL     , DIMENSION(3,3) :: recipr
   REAL     , DIMENSION(3  ) :: vector
   REAL     , DIMENSION(3  ) :: finalv
   REAL     , DIMENSION(:), ALLOCATABLE :: deviat  ! distances to plane
!
   surf_char = SURF_NONE
!
   fp (1) = chem_period (1)
   fp (2) = chem_period (2)
   fp (3) = chem_period (3)
   fq = chem_quick
   ALLOCATE(deviat(0:MAX_ATOM_ENV), STAT = istat)
!
!   atomloop: DO ia=istart,iend
!
     IF(IBITS(cr_prop(ia),PROP_SURFACE_EXT,1).eq.1 .and.        &  ! real Atom is near surface
        IBITS(cr_prop(ia),PROP_OUTSIDE    ,1).eq.0       ) THEN    ! real Atom is near surface
        cpara(1) = 'all'
        lpara(1) = 3
        ianz     = 1
        x(1)     = cr_pos(1,ia)
        x(2)     = cr_pos(2,ia)
        x(3)     = cr_pos(3,ia)
        rmin     = 0.0
        radius   = 6.5
        CALL get_iscat (ianz, cpara, lpara, werte, maxw, lnew)
        ianz = 1
        werte(1) = -1                                              ! Find all atom types
        CALL do_find_env (ianz, werte, maxw, x, rmin,&
                          radius, fq, fp)
        IF(atom_env(0) > 3 ) THEN
        direct = 0.0
        vector = 0.0
        DO ib = 1, atom_env(0) 
           IF(IBITS(cr_prop(atom_env(ib)),PROP_SURFACE_EXT,1).eq.1 ) THEN          ! real Atom is near surface
           direct(1,1) = direct(1,1) + atom_pos(1,ib)**2             ! Accumulate x^2 etc
           direct(2,2) = direct(2,2) + atom_pos(2,ib)**2
           direct(3,3) = direct(3,3) + atom_pos(3,ib)**2
           direct(1,2) = direct(1,2) + atom_pos(1,ib)*atom_pos(2,ib)
           direct(1,3) = direct(1,3) + atom_pos(1,ib)*atom_pos(3,ib)
           direct(2,3) = direct(2,3) + atom_pos(2,ib)*atom_pos(3,ib)
           vector(1  ) = vector(1  ) + atom_pos(1,ib)
           vector(2  ) = vector(2  ) + atom_pos(2,ib)
           vector(3  ) = vector(3  ) + atom_pos(3,ib)
           ENDIF
        ENDDO
        direct(2,1) = direct(1,2) ! copy into transposed elements
        direct(3,1) = direct(1,3)
        direct(3,2) = direct(2,3)
!
        CALL invmat( recipr, direct )  ! calculate inverse matrix
        finalv      = MATMUL(recipr, vector) ! finalv contains the hkl of the plane
        sigma = 0.0
        DO ib = 1, atom_env(0) 
     IF(IBITS(cr_prop(atom_env(ib)),PROP_SURFACE_EXT,1).eq.1 ) THEN          ! real Atom is near surface
              deviat(ib) = (atom_pos(1,ib)*finalv(1)+      &  ! accumulate the distances
                            atom_pos(2,ib)*finalv(2)+      &
                            atom_pos(3,ib)*finalv(3)-1)
              sigma      = sigma + deviat(ib)**2              ! calculate the sigma of the distance distribution
     ENDIF
        ENDDO
        sigma = SQRT(sigma)
        surf_normal(:) = finalv(:)
        surf_char = SURF_PLANE
        IF(sigma <    0.05) THEN
!           cr_iscat(ia) = cr_iscat(ia) + 4 !cr_nscat
!          CYCLE atomloop
        ENDIF
        ENDIF
!
!       Atom is not in plane, try a row
!
     ENDIF
!
!  ENDDO atomloop
!
   DEALLOCATE(deviat, STAT = istat)
!
   IF(skalpro(surf_normal,surf_normal, cr_gten) < EPS) THEN
      surf_normal(:) = cr_pos(:,ia)
      surf_char = SURF_ATOM
   ENDIF
!
   END SUBROUTINE find_surface_character
!
!######################################################################
!
   SUBROUTINE deco_run
!
!  Performs the actual decoration
!
   USE chem_menu
   USE conn_mod
   USE conn_def_mod
   USE chem_aver_mod
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
   USE prop_para_mod
   USE read_internal_mod
   USE structur , ONLY: rese_cr
   USE save_menu, ONLY: save_internal, save_store_setting, save_restore_setting, save_default_setting
!
   USE param_mod
   USE random_mod
!
   IMPLICIT none
!
   INTEGER, PARAMETER         :: MAXW = 200  ! not ideal, should be dynamic ....
!
   CHARACTER(LEN=1024) :: line      ! a string
   CHARACTER(LEN=1024) :: mole_name ! molecule file name
   CHARACTER(LEN= 200) :: corefile  ! original structure file name
   CHARACTER(LEN= 200) :: corelist  ! Place holder for core position at 0,0,0
   CHARACTER(LEN= 200) :: shellfile  ! original structure file name
   CHARACTER(LEN=1024), DIMENSION(MAXW) :: cpara      ! a string
   INTEGER            , DIMENSION(MAXW) :: lpara      ! a string
   REAL               , DIMENSION(MAXW) :: werte      ! a string
   INTEGER   :: istatus          ! status
   INTEGER   :: i  ! dummy index
   INTEGER   :: ia ! dummy index
   INTEGER   :: ianz ! dummy number of parameters
   INTEGER   :: length ! dummy length
   INTEGER   :: mole_length ! length of molecule file name
   INTEGER   :: istart,iend ! dummy index
   INTEGER   :: is ! scattering number of surface atom
   INTEGER   :: idef ! connectivity definition number
   INTEGER   :: ncon ! number of connections defined
   INTEGER, DIMENSION(:), ALLOCATABLE :: c_list
   INTEGER, DIMENSION(:,:), ALLOCATABLE :: c_offs ! Result of connectivity search
   INTEGER   :: natoms                 ! number of atoms in connectivity
   INTEGER   :: n_scat                 ! Dummy maximum scattering types
   INTEGER   :: n_corr                 ! Dummy maximum correlation number
   INTEGER   :: nscat_old              ! maximum scattering types prior to modifying anchors
   INTEGER   :: n_repl                 ! counter for anchor aatoms
   REAL                    :: r_anch   ! relative amount of anchor atoms
   REAL                    :: prob     ! Probability to replace by non_anchor dummy
   REAL   , DIMENSION(1:3) :: xyz
!
   REAL ran1
!
!  Section for automatic distribution of all surface anchors
!
   CALL save_store_setting             ! Backup user "save" setting
   CALL save_default_setting           ! Default to full saving
   corefile   = 'internal.decorate'             ! internal user files always start with 'internal'
   CALL save_internal(corefile)        !     thus this file name is unique
   shellfile  = 'internal.decoshell'   
   line       = 'ignore, all'          ! Ignore all properties
   length     = 11
   CALL property_select(line, length, sav_sel_prop)
   line       = 'present, external'    ! Force atom to be close to a surface
   length     = 17
   CALL property_select(line, length, sav_sel_prop)
   line       = 'absent, outside'      ! Force atom to be inside
   length     = 15
write(*,*) ' ALL     ATOMS,             ', cr_natoms
   CALL property_select(line, length, sav_sel_prop)
   CALL save_internal(shellfile)
!
!  Make single atom "structure" for domain list file 
   CALL rese_cr
   cr_natoms    = 1
   cr_ncatoms   = 1
   cr_nscat     = 1
   cr_pos(:,1)  = 0.0
   cr_iscat(1)  = 1
   cr_at_lis(1) = 'CORE'
   corelist     = 'internal.core.list'
   CALL save_internal(corelist)        ! Save the core list
!
!     Load the molecules into temporary structures to reduce disk I/O
!
   CALL deco_get_molecules
!
   CALL rese_cr
!
   CALL readstru_internal(shellfile)   ! Read shell file
   IF(cr_natoms > 0) THEN              ! The Shell does consist of atoms
!
!     Now sort those surface atoms that are anchors to the ligand
!
      dc_def_temp => dc_def_head
      CALL dc_get_dens(dc_def_temp, dc_temp_dens)    ! Get density for surface coverage 
      CALL dc_get_con(dc_con_temp, dc_temp_surf, dc_temp_neig, dc_temp_dist)
      CALL chem_elem(.false.)             ! get composition
      r_anch = res_para(dc_temp_surf(1)+1)           ! Fractional composition of the anchoring atoms
      prob   = MAX(0.0,MIN(1.0,DC_AREA*dc_temp_dens/r_anch)) ! replacement probability
!     Replace anchors by a new atom type
      nscat_old = cr_nscat
      IF(cr_nscat == MAXSCAT) THEN                   ! Number of scattering types increased
         n_scat = MAX(cr_nscat+5,MAXSCAT)
         natoms = MAX(cr_natoms, NMAX)
         CALL alloc_crystal(n_scat,natoms)
      ENDIF
      cr_at_lis(nscat_old + 1) = 'AN01'              ! Fixed name for anchor
      cr_dw    (nscat_old + 1) = cr_dw(dc_temp_surf(1))
      cr_nscat = nscat_old + 1
      n_repl = 0
!     As surface atom number is bound to be small try several tims unit we get sufficient anchors
      i = 0
      replace: DO i=1, cr_natoms 
         DO ia = 1, cr_natoms                        ! Loop to replace
            IF(cr_iscat(ia)==dc_temp_surf(1)) THEN   ! Got a surface atom of correct type
            IF(ran1(idum) < prob) THEN            ! Randomly pick a fraction
                  cr_iscat(ia) = nscat_old + 1       
                  cr_prop (ia) = IBSET (cr_prop (ia), PROP_DECO_ANCHOR)  ! FLAG THIS ATOM AS SURFACE ANCHOR
                  n_repl       = n_repl  + 1         ! Increment replaced atoms
                  IF(n_repl == NINT(cr_natoms*r_anch*prob)) EXIT replace   ! got enough anchors
               ENDIF
            ENDIF
         ENDDO
      ENDDO replace
      IF(n_repl > 0 ) THEN                        ! Need at least one anchor
!
!        Prepare a connectivity 
!
         WRITE(line,1000), 'AN01','AN01', cr_at_lis(dc_temp_surf(1)), 'deco_0001'
1000 FORMAT(3(a4,','),' 0.5, 18.0, ',a9)
         length = 36
         CALL conn_do_set (code_add,line, length)          ! Add connectivity around "AN01"
         WRITE(line,1000),cr_at_lis(dc_temp_surf(1)),'AN01', cr_at_lis(dc_temp_surf(1)), 'deco_0002'
         CALL conn_do_set (code_add,line, length)          ! Add connectivity around  other surface atoms
!        CALL conn_show
         CALL create_connectivity                             ! Create actual connecitivity list
!
!        SORT ATOMS WITH MMC REPULSIVE
!
         n_corr = MAX(CHEM_MAX_COR,MMC_MAX_CORR)
         n_scat = MAX(MAXSCAT, MMC_MAX_SCAT)
         CALL alloc_mmc ( n_corr, MC_N_ENERGY, n_scat )       ! Basic mmc allocation
         ianz = 1
         cpara(1) = 'rese'                                    ! Prepare "set con, reset"
         lpara(1) = 4
         CALL chem_set_con (ianz, cpara, lpara, werte, maxw)  ! Reset conn list
         cpara(1) = 'rese'                                    ! Prepare "set neig, reset"
         CALL chem_set_neig(ianz, cpara, lpara, werte, maxw)  ! Reset neig list
         ianz = 3
         cpara(1) = '1'                                       ! Prepare "set con, 1, AN01, deco_0001"
         lpara(1) = 1
         cpara(2) = 'AN01'
         lpara(2) = 4
         cpara(3) = 'deco_0001'
         lpara(3) = 9
         CALL chem_set_con (ianz, cpara, lpara, werte, maxw)  ! Set conn list 1 'set con,1, AN01, deco_0001'
         ianz = 3
         cpara(1) = '2'                                       ! Prepare "set con, 2, <sur>, deco_0002"
         lpara(1) = 1
         cpara(2) = cr_at_lis(dc_temp_surf(1))
         lpara(2) = 4
         cpara(3) = 'deco_0002'
         lpara(3) = 9
         CALL chem_set_con (ianz, cpara, lpara, werte, maxw)  ! Set conn list 2 'set con,2, <surf>, deco_0002'
         ianz = 3
         cpara(1) = 'con'                                     ! Prepare "set neig, con, 1, 2"
         lpara(1) = 3
         cpara(2) = '1'
         lpara(2) = 1
         cpara(3) = '2'
         lpara(3) = 1
         CALL chem_set_neig(ianz, cpara, lpara, werte, maxw)  ! Set neig 'set neigh,1, 2'
         ianz = 2
         cpara(1) = 'swc'                                     ! Prepare "set mode, swchem, 1.0, all"
         lpara(1) = 3
         cpara(2) = 'all'
         lpara(2) = 1
         werte(1) = 1.000
         cpara(3) = ' '
         lpara(3) = 1
         CALL  mmc_set_mode(ianz, cpara, lpara, werte, maxw)  ! Set mode 'set mode,swchem, 1.0, all'
!
         mmc_allowed(dc_temp_surf(1)) = .true.                ! 'set allowed AN01   and
         mmc_allowed(nscat_old+1    ) = .true.                !  actual anchor type
!
         IF(1 > MMC_REP_CORR .or.  1 > CHEM_MAX_COR  .or. &   ! Allocate Repulsive
            MAXSCAT > MMC_REP_SCAT                         ) THEN
            n_corr = MAX(n_corr, CHEM_MAX_COR, MMC_REP_CORR)
            n_scat = MAX(MAXSCAT,MMC_REP_SCAT)
            CALL alloc_mmc_rep (n_corr, n_scat)
!        IF(ier_num /= 0) THEN                             ! Does not seem to fail :-)
!           RETURN
!        ENDIF
         ENDIF
!                                                          ! define repulsive energy
         is = nscat_old+1
         CALL mmc_set_disp (1, MC_REPULSIVE, is, is, 100.0, 15.0)
         CALL mmc_set_rep  (1, is, is, 15.,16000., 0.5, 1.)
         mmc_cor_energy (1, MC_REPULSIVE) = .true.
         mmc_cor_energy (0, MC_REPULSIVE) = .true.
!
         mo_cyc  = 100*cr_natoms                              ! Define cycles
         mo_feed =   5*cr_natoms                              ! Define feedback
         mo_kt   =   2.5                                      ! Define Temperature
!
!        CALL mmc_show
         CALL mmc_run_multi                                   ! Run actual sorting

!        Change Anchors back to their original names, keep property flag
!
         DO ia = 1, cr_natoms
            IF(cr_at_lis(cr_iscat(ia))=='AN01') THEN
               cr_iscat(ia) = dc_temp_surf(1)
            ENDIF
         ENDDO
         cr_nscat = nscat_old                                 ! Set number of atom types back
         cr_at_lis(nscat_old+1) = ' '                         ! Have atom type disappear
!
!  Transform atom coordinates into cartesian space to ease computations
!
!   CALL plot_ini_trans (1.0)
!   CALL trans_atoms_tocart(uvw_out)
!
         istart = 1
         iend   = cr_natoms
         ier_num = 0
         main_loop: DO ia=istart,iend
           is_sel: IF(check_select_status (dc_latom (cr_iscat (ia) ), cr_prop (ia),   &
                                        dc_sel_prop)              ) THEN
              is     = cr_iscat(ia)                  ! get scattering type central
              xyz(:) = cr_pos(:,ia)                  ! Get atom position
              idef   = dc_use_conn(is)               ! use this connectivity list for surface
              CALL get_connectivity_list (ia, is, idef, maxw, c_list, c_offs, natoms )
!             Go through all definitions
              dc_def_temp => dc_def_head        
              defs: DO WHILE(ASSOCIATED(dc_def_temp))
                 dc_con_temp => dc_def_temp%dc_def_con
                 cons: DO WHILE(ASSOCIATED(dc_con_temp))           ! A connectivity exists
                    CALL dc_get_con(dc_con_temp, dc_temp_surf, dc_temp_neig, dc_temp_dist)
                    IF(cr_iscat(ia) == dc_temp_surf(1).AND.        &
                       BTEST(cr_prop(ia),PROP_DECO_ANCHOR) ) THEN    ! Matching atom in current definition
                       CALL dc_get_type(dc_def_temp, dc_temp_type)
                       CALL dc_get_axis(dc_def_temp, dc_temp_axis)
                       CALL dc_get_mole_name(dc_def_temp, mole_name, mole_length)
                       CALL dc_get_ncon(dc_def_temp, ncon)
                       SELECT CASE(dc_temp_type)
                          CASE ( DC_NORMAL )                 ! Molecule in normal position
                             IF(ncon == 1) THEN
                                CALL deco_place_normal(dc_def_temp, ia, is, xyz, &
                                  dc_temp_axis, mole_name, mole_length,       &
                                  dc_temp_surf, dc_temp_neig, dc_temp_dist)
                             ELSE
                                ier_num = -1118
                                ier_msg(1) = 'The bridge connection requires one bond'
                                EXIT main_loop
                             ENDIF
                          CASE ( DC_BRIDGE )                 ! Molecule in bridge position
                             IF(ncon == 1) THEN
                             CALL deco_place_bridge(dc_def_temp, ia, is, xyz, &
                                  dc_temp_axis, mole_name, mole_length,       &
                                  dc_temp_surf, dc_temp_neig, dc_temp_dist)
                             ELSE
                                ier_num = -1118
                                ier_msg(1) = 'The bridge connection requires one bond'
                                EXIT main_loop
                             ENDIF
                          CASE ( DC_DOUBLE   )               ! Molecule in double   connection position
                             IF(ncon >  1) THEN
                             CALL deco_place_double(dc_def_temp, ia, is, xyz, &
                                  dc_temp_axis, mole_name, mole_length,       &
                                  dc_temp_surf, dc_temp_neig, dc_temp_dist, ncon)
                                EXIT cons
                             ELSE
                                ier_num = -1118
                                ier_msg(1) = 'The mult   connection requires > one bond'
                                EXIT main_loop
                             ENDIF
                          CASE ( DC_MULTIPLE )               ! Molecule in multiple connection position
                             IF(ncon >  1) THEN
                             CALL deco_place_multi(dc_def_temp, ia, is, xyz, &
                                  dc_temp_axis, mole_name, mole_length,       &
                                  dc_temp_surf, dc_temp_neig, dc_temp_dist, ncon)
                                EXIT cons
                             ELSE
                                ier_num = -1118
                                ier_msg(1) = 'The mult   connection requires > one bond'
                                EXIT main_loop
                             ENDIF
                       END SELECT
                    ENDIF
                    dc_con_temp => dc_con_temp%next
                 ENDDO cons
                 dc_def_temp => dc_def_temp%next
              ENDDO defs
           ENDIF is_sel ! END IF BLOCK is selected
         ENDDO main_loop   ! END DO main loop over all atoms
!
         IF(ier_num == 0)  THEN       ! Success in main_loop
!
!           Use domain to insert core back into the structure
!
            CALL alloc_domain ( clu_increment )
            MK_MAX_SCAT = MAX(MK_MAX_SCAT, MAXSCAT)
            MK_MAX_ATOM = MAX(MK_MAX_ATOM, NMAX)
            CALL alloc_micro  ( MK_MAX_SCAT , MK_MAX_ATOM)
!
            clu_infile          = corelist
            clu_infile_internal = .true.
            clu_mode            = CLU_IN_PSEUDO
            clu_remove_mode     = CLU_REMOVE_STRICT
            clu_remove_end      = iend     ! we do not need to check the ligand
            clu_remove_dist     = 0.50
            clu_content(1)      = corefile
            clu_name(1)         = 'CORE'
            clu_character(1)    = CLU_CHAR_FUZZY
            clu_fuzzy(1)        = 0.5
            clu_orient(1,:,:)   = 0.0
            clu_orient(1,1,1)   = 1.0
            clu_orient(1,2,2)   = 1.0
            clu_orient(1,3,3)   = 1.0
            clu_shape (1,:,:)   = 0.0
            clu_shape (1,1,1)   = 1.0
            clu_shape (1,2,2)   = 1.0
            clu_shape (1,3,3)   = 1.0
            clu_sigma (1,:  )   = 0.0
            clu_index           = 1
            clu_number          = 1
!
            CALL micro_filereading
            DO ia=istart, iend
               cr_prop (ia) = IBCLR (cr_prop (ia), PROP_DECO_ANCHOR)  ! FLAG THIS ATOM AS SURFACE ANCHOR
            ENDDO
            CALL do_purge
         ELSE     ! Error in main_loop
           CALL readstru_internal(shellfile)   ! Read shell file
         ENDIF
      ELSE     ! n_repl > 0   !! No anchor atoms found
        CALL readstru_internal( corefile)   ! Read  core file
        ier_num = -131
        ier_typ = ER_APPL
        ier_msg(1) = 'Is the surface very small, just a few atoms?'
        ier_msg(2) = 'Is the coverage too small? '
        ier_msg(3) = 'Check the set ligand command'
      ENDIF
   ELSE     ! SHELL has atoms
     CALL readstru_internal( corefile)   ! Read  core file
     ier_num = -130
     ier_typ = ER_APPL
     ier_msg(1) = 'Possible reasons: no boundary was used to cut'
     ier_msg(2) = 'Distance to external surface is too large'
     ier_msg(3) = 'Check settings and command in surface menu'
   ENDIF
!
   IF(ier_num == 0 ) THEN
      CALL save_restore_setting     ! Restore user save settigns
   ENDIF
!
!  Clean up temporary arrays
!
   DO i=1, dc_n_molecules
      CALL dc_molecules(i)%finalize_atoms()
   ENDDO
   DEALLOCATE(dc_molecules)
   DEALLOCATE(m_ntypes, STAT=istatus)
   DEALLOCATE(m_length, STAT=istatus)
   DEALLOCATE(c_list  , STAT=istatus)
!
!  Transform atom coordinates back into crystal space 
!
!   CALL trans_atoms_fromcart()
!
   END SUBROUTINE deco_run
!
!######################################################################
   SUBROUTINE deco_get_molecules
!
   USE structur, ONLY: test_file
   USE trafo_mod
   IMPLICIT none
!
   CHARACTER (LEN=1024)                :: strufile
!
   INTEGER  :: i,j,k            ! dummy index
   INTEGER  :: istatus          ! status
   INTEGER  :: natoms           ! number of atoms in file
   INTEGER  :: ntypes           ! number of atom types in file
   INTEGER  :: n_mole           ! number of molecules
   INTEGER  :: n_type           ! number of molecule types
   INTEGER  :: n_atom           ! number of atoms in molecules
   INTEGER  :: init   = -1      ! Initialize atom names/types
   INTEGER  :: itype            ! atom scattering type
   REAL, DIMENSION(3) :: posit  ! atom position
   REAL, DIMENSION(3) :: uvw    ! atom position
   REAL, DIMENSION(4,4) :: rd_tran_f  ! transformation matrix to cartesian
   INTEGER  :: iprop            ! atom property 
   LOGICAL  :: lcell  = .true.  ! Treat atoms with equal name and B as one type
!
   ALLOCATE(m_name  (dc_n_molecules), STAT=istatus)
   ALLOCATE(m_lname (dc_n_molecules), STAT=istatus)
   ALLOCATE(m_ntypes(dc_n_molecules), STAT=istatus)
   ALLOCATE(m_length(dc_n_molecules), STAT=istatus)
!
   ALLOCATE(dc_molecules(dc_n_molecules), STAT = istatus)
!
   m_name(:) = ' '
   m_lname(:) = 0
   m_ntypes(:) = 0
   m_length(:) = 0
!
   DO i=1,dc_n_molecules        ! load all molecules
      strufile = dc_input(i)
      CALL test_file(strufile, natoms, ntypes, n_mole, n_type, &
                     n_atom, init, lcell)
      CALL dc_molecules(i)%alloc_arrays(natoms, ntypes, n_mole, n_atom)
      CALL read_crystal ( dc_molecules(i), strufile )
      CALL dc_molecules(i)%get_cryst_tran_f(rd_tran_f)  ! Get transformation matrix to cartesian
!
      DO j=1, natoms
         CALL dc_molecules(i)%get_cryst_atom ( j, itype, posit, iprop)
         CALL trans(posit,rd_tran_f ,uvw, 4)       ! Transform mol to cartesian
         CALL trans(uvw  ,cr_tran_fi,posit, 4)     ! Transfrom cartesian to crystal
         CALL dc_molecules(i)%set_cryst_atom ( j, itype, posit, iprop)
      ENDDO
      m_length(i) = natoms
      m_ntypes(i) = ntypes
      m_lname(i)  = LEN_TRIM(strufile)
      m_name(i)   = strufile(1:m_lname(i))
   ENDDO
!
   END SUBROUTINE deco_get_molecules
!
!######################################################################
!
   SUBROUTINE deco_show
!
   IMPLICIT none
!
   dc_def_temp => dc_def_head
   CALL dc_show_def(dc_def_temp, ier_num)
!
   END SUBROUTINE deco_show
!
!######################################################################
!
   SUBROUTINE deco_reset
!
!  Set all definitions back to system default
!
   IMPLICIT none
   INTEGER :: istatus
!
   dc_init        = .true.    ! We need to initialize
   dc_n_molecules = 0         ! There are no molecules
   dc_latom(:)    = .false.
   dc_use_conn(:) = 0
   CALL dc_reset_def ( dc_def_head)
   NULLIFY(dc_def_head)
   NULLIFY(dc_def_temp)
   DEALLOCATE(m_name      , STAT=istatus)
   DEALLOCATE(m_lname     , STAT=istatus)
   DEALLOCATE(m_ntypes    , STAT=istatus)
   DEALLOCATE(m_ntypes    , STAT=istatus)
   DEALLOCATE(m_length    , STAT=istatus)
   DEALLOCATE(dc_molecules,STAT = istatus)
   dc_temp_type = DC_NONE
!
   END SUBROUTINE deco_reset
!
!******************************************************************************
!
   SUBROUTINE read_crystal ( this, infile)
!
!  Read a crystal structure from file
!  This procedure interfaces to the old "readtru" in "structur.f90"
!
   USE discus_plot_init_mod
   USE spcgr_apply
   USE inter_readstru
   USE structur, ONLY: readstru
   USE trans_cart_mod
!
   IMPLICIT none
!
   LOGICAL, PARAMETER :: lout = .true.
   TYPE (cl_cryst)                  :: this   ! The current crystal
   CHARACTER (LEN=*   ), INTENT(IN)  :: infile ! Disk file
   INTEGER                           :: inum   ! dummy index
   REAL   , DIMENSION(3)             :: posit  ! dummy position vector
   REAL   , DIMENSION(3)             :: uvw_out ! dummy position vector
   INTEGER                           :: istat  ! status variable
   integer :: i,j
   
!
   rd_strucfile = infile
   rd_NMAX      = this%get_natoms() ! cr_natoms
   rd_MAXSCAT   = this%get_nscat()  ! cr_nscat
!
!  Allocate temporary arrays
!
   ALLOCATE ( rd_cr_dw    (0:rd_MAXSCAT)    , STAT = istat )
   ALLOCATE ( rd_cr_at_lis(0:rd_MAXSCAT)    , STAT = istat )
   ALLOCATE ( rd_cr_as_lis(0:rd_MAXSCAT)    , STAT = istat )
   ALLOCATE ( rd_cr_at_equ(0:rd_MAXSCAT)    , STAT = istat )
   ALLOCATE ( rd_cr_scat_equ(0:rd_MAXSCAT)    , STAT = istat )
   ALLOCATE ( rd_cr_scat_int(0:rd_MAXSCAT)    , STAT = istat )
   ALLOCATE ( rd_cr_delf_int(0:rd_MAXSCAT)    , STAT = istat )
   ALLOCATE ( rd_sav_latom(0:rd_MAXSCAT)    , STAT = istat )
   ALLOCATE ( rd_cr_delfi (0:rd_MAXSCAT)    , STAT = istat )
   ALLOCATE ( rd_cr_delfr (0:rd_MAXSCAT)    , STAT = istat )
   ALLOCATE ( rd_cr_scat(11,0:rd_MAXSCAT)    , STAT = istat )
   ALLOCATE ( rd_cr_pos   (1:3,1:rd_NMAX)   , STAT = istat )
   ALLOCATE ( rd_cr_iscat (1:rd_NMAX)       , STAT = istat )
   ALLOCATE ( rd_cr_prop  (1:rd_NMAX)       , STAT = istat )
   ALLOCATE ( rd_cr_mole  (1:rd_NMAX)       , STAT = istat )
   ALLOCATE ( rd_as_at_lis(0:rd_MAXSCAT)    , STAT = istat )
   ALLOCATE ( rd_as_dw    (0:rd_MAXSCAT)    , STAT = istat )
   ALLOCATE ( rd_as_pos   (1:3,1:rd_MAXSCAT), STAT = istat )
   ALLOCATE ( rd_as_iscat (1:rd_MAXSCAT)    , STAT = istat )
   ALLOCATE ( rd_as_prop  (1:rd_MAXSCAT)    , STAT = istat )
!
   rd_cr_dw    (:) = 0.0
   rd_cr_at_lis(:) = ' '
   rd_cr_as_lis(:) = ' '
   rd_cr_at_equ(:) = ' '
   rd_cr_scat_equ(:) = .false.
   rd_cr_scat_int(:) = .true.
   rd_cr_delf_int(:) = .true.
   rd_sav_latom  (:) = .true.
   rd_cr_delfi (:) = 0.0
   rd_cr_delfr (:) = 0.0
   rd_cr_delfr (:) = 0.0
   rd_cr_pos (:,:) = 0.0
   rd_cr_iscat (:) = 0
   rd_cr_prop  (:) = 0
   rd_cr_mole  (:) = 0
   rd_as_at_lis(:) = ' '
   rd_as_dw    (:) = 0.0
   rd_as_pos (:,:) = 0.0
   rd_as_iscat (:) = 0
   rd_as_prop  (:) = 0
   rd_cr_natoms    = 0
   rd_cr_nscat     = 0
!
   CALL readstru (rd_NMAX, rd_MAXSCAT, rd_strucfile, rd_cr_name,        &
               rd_cr_spcgr, rd_cr_a0, rd_cr_win, rd_cr_natoms, rd_cr_nscat, rd_cr_dw,     &
               rd_cr_at_lis, rd_cr_pos, rd_cr_mole, rd_cr_iscat, rd_cr_prop, rd_cr_dim, rd_as_natoms, &
               rd_as_at_lis, rd_as_dw, rd_as_pos, rd_as_iscat, rd_as_prop, rd_sav_ncell,  &
               rd_sav_r_ncell, rd_sav_ncatoms, rd_spcgr_ianz, rd_spcgr_para)
!
   CALL lattice (rd_cr_a0, rd_cr_ar, rd_cr_eps, rd_cr_gten, rd_cr_reps, rd_cr_rten,  &
                 rd_cr_win, rd_cr_wrez, rd_cr_vol, rd_cr_vr, lout,                   &
                 rd_tran_g, rd_tran_gi, rd_tran_f, rd_tran_fi)
!
   rd_n_latom = rd_MAXSCAT
   CALL this%set_crystal_save_flags (rd_sav_scat, & 
            rd_sav_adp, rd_sav_gene, rd_sav_symm,                     &
            rd_sav_w_ncell, rd_sav_obje, rd_sav_doma, rd_sav_mole, rd_sav_prop, &
            rd_sav_sel_prop,rd_n_latom,rd_sav_latom)
   CALL this%set_crystal_from_local   ( rd_strucfile, &
                                         rd_NMAX, rd_MAXSCAT, rd_cr_name,      &
            rd_cr_natoms, rd_cr_ncatoms, rd_cr_n_REAL_atoms, rd_cr_spcgrno, rd_cr_syst, &
            rd_cr_spcgr, rd_cr_at_lis, rd_cr_at_equ, rd_cr_as_lis,                      &
            rd_cr_nscat, rd_cr_dw, rd_cr_a0, rd_cr_win,                                 &
            rd_cr_ar, rd_cr_wrez, rd_cr_vol, rd_cr_vr, rd_cr_dim, rd_cr_dim0, rd_cr_icc,  &
            rd_sav_ncell, rd_sav_r_ncell, rd_sav_ncatoms, rd_spcgr_ianz, rd_spcgr_para, &
            rd_tran_g, rd_tran_gi, rd_tran_f, rd_tran_fi, rd_cr_gmat, rd_cr_fmat, &
            rd_cr_gten, rd_cr_rten, rd_cr_eps, rd_cr_reps,                              &
            rd_cr_acentric, rd_cr_newtype, rd_cr_cartesian, rd_cr_sel_prop,             &
            rd_cr_scat, rd_cr_delfi , rd_cr_delfr, rd_cr_delf_int,                    &
            rd_cr_scat_int, rd_cr_scat_equ,                                             &
            rd_cr_pos, rd_cr_iscat, rd_cr_prop                                          &
            )
!
   DEALLOCATE ( rd_cr_dw    , STAT = istat )
   DEALLOCATE ( rd_cr_at_lis, STAT = istat )
   DEALLOCATE ( rd_cr_as_lis, STAT = istat )
   DEALLOCATE ( rd_cr_at_equ, STAT = istat )
   DEALLOCATE ( rd_cr_scat_equ, STAT = istat )
   DEALLOCATE ( rd_cr_scat_int, STAT = istat )
   DEALLOCATE ( rd_cr_delf_int, STAT = istat )
   DEALLOCATE ( rd_sav_latom, STAT = istat )
   DEALLOCATE ( rd_cr_delfi , STAT = istat )
   DEALLOCATE ( rd_cr_delfr , STAT = istat )
   DEALLOCATE ( rd_cr_scat  , STAT = istat )
   DEALLOCATE ( rd_cr_pos   , STAT = istat )
   DEALLOCATE ( rd_cr_iscat , STAT = istat )
   DEALLOCATE ( rd_cr_prop  , STAT = istat )
   DEALLOCATE ( rd_cr_mole  , STAT = istat )
   DEALLOCATE ( rd_as_at_lis, STAT = istat )
   DEALLOCATE ( rd_as_dw    , STAT = istat )
   DEALLOCATE ( rd_as_pos   , STAT = istat )
   DEALLOCATE ( rd_as_iscat , STAT = istat )
   DEALLOCATE ( rd_as_prop  , STAT = istat )
!
   END SUBROUTINE read_crystal
!
!******************************************************************************
!
   SUBROUTINE deco_place_normal(dc_def_temp, ia, ia_scat, xyz, &
                            mole_axis, mole_name, mole_length, &
                            surf, neig, dist)
!
   USE chem_mod
   USE metric_mod
   USE modify_mod
   USE prop_para_mod
   USE symm_menu
   USE symm_mod
   USE symm_sup_mod
   USE trafo_mod
!
   USE param_mod
!
   IMPLICIT NONE
!
   TYPE (dc_def), POINTER              :: dc_def_temp     ! The definition to be used
   INTEGER,                 INTENT(IN) :: ia              ! Surface atom number
   INTEGER,                 INTENT(IN) :: ia_scat         ! Scattering type of surface atom
   REAL   , DIMENSION(1:3), INTENT(IN) :: xyz             ! Position of surface atom
   INTEGER, DIMENSION(1:2), INTENT(IN) :: mole_axis       ! Atoms that define molecule axis
   CHARACTER (LEN=1024),    INTENT(IN) :: mole_name       ! Molecule file name
   INTEGER,                 INTENT(IN) :: mole_length     ! Molecule file name length
   INTEGER, DIMENSION(0:4), INTENT(IN) :: surf            ! Surface atom type
   INTEGER,                 INTENT(IN) :: neig            ! Connected to this neighbor in mole
   REAL   ,                 INTENT(IN) :: dist            ! distance to ligand molecule
!
   REAL, PARAMETER :: EPS = 1.0E-6
!
   CHARACTER (LEN=4) :: atom_name
   CHARACTER (LEN=1024) :: line
   INTEGER :: i, im, laenge
   INTEGER :: itype
   INTEGER   :: surf_char ! Surface character, plane, edge, corner, ...
   REAL   , DIMENSION(1:3) :: surf_normal
   REAL, DIMENSION(3) :: posit
   REAL, DIMENSION(3) :: vnull
   REAL, DIMENSION(3) :: origin
   REAL               :: normal_l
   REAL               :: dw1
   INTEGER :: iprop
   LOGICAL, PARAMETER :: lspace = .true.
   REAL   , DIMENSION(1:3) :: axis_ligand                 ! Initial molecule orientation
!
   vnull(:) = 0.00
!
   CALL find_surface_character(ia,surf_char, surf_normal)
!  IF(surf_char == SURF_PLANE ) THEN                      ! Ignore other than planar surfaces
!
      moles: DO i=1, dc_n_molecules                       ! Loop over all loaded molecules
         IF(mole_name(1:mole_length) == m_name(i)(1:m_lname(i))) THEN
            im = mole_axis(2)                             ! Make mole axis from spcified atoms
            CALL dc_molecules(i)%get_cryst_atom(im, itype, posit, iprop)
            axis_ligand(1) = posit(1)
            axis_ligand(2) = posit(2)
            axis_ligand(3) = posit(3)
            im = mole_axis(1)
            CALL dc_molecules(i)%get_cryst_atom(im, itype, posit, iprop)
            axis_ligand(1) = axis_ligand(1) - posit(1)   ! this is mole axis
            axis_ligand(2) = axis_ligand(2) - posit(2)
            axis_ligand(3) = axis_ligand(3) - posit(3)
!           Insert molecule atoms into crystal at correct origin in initial orientation
            normal_l = sqrt (skalpro (surf_normal, surf_normal, cr_gten))
            origin(1)  = cr_pos(1,ia) + surf_normal(1)/normal_l*dist  ! Origin is shifted
            origin(2)  = cr_pos(2,ia) + surf_normal(2)/normal_l*dist  ! by dist away from 
            origin(3)  = cr_pos(3,ia) + surf_normal(3)/normal_l*dist  ! surface atom
            sym_latom(:) = .false.                        ! Initially deselect all atomtypes
            atoms: DO im=1,m_length(i)                    ! Load all atoms from the molecule
               CALL dc_molecules(i)%get_cryst_atom(im, itype, posit, iprop)
               CALL dc_molecules(i)%get_cryst_scat(im, itype, atom_name, dw1)
               posit(:) = posit(:) + origin(:)
               WRITE(line, 1000) atom_name, posit, dw1
               laenge = 60
               CALL do_ins(line, laenge)                  ! Insert into crystal
               CALL check_symm
               sym_latom(cr_iscat(cr_natoms)) = .true.    ! Select atopm type for rotation
            ENDDO atoms
! define rotation operation
            sym_angle      = do_bang(lspace, surf_normal, vnull, axis_ligand)
            IF(ABS(sym_angle) > EPS ) THEN                ! Rotate if not zero degrees
               sym_orig(:)    = origin(:)                    ! Define origin
               sym_trans(:)   = 0.0                          ! No translation needed
               sym_sel_atom   = .true.                       ! Select atoms
               sym_new        = .false.                      ! No new types
               sym_power      =  1                           ! Just need one operation
               sym_type       = .true.                       ! Proper rotation
               sym_mode       = .false.                      ! Move atom to new position
               sym_orig_mol   = .false.                      ! Origin at crystal
               sym_power_mult =.false.                       ! No multiple copies
               sym_sel_atom   = .true.                       ! Select atoms not molecules
               sym_start      =  cr_natoms - m_length(i) + 1 ! set range of atoms numbers
               sym_end        =  cr_natoms
               IF(ABS(sym_angle-180.) < EPS ) THEN           ! Ligand and surface normal are antiparallel
                  WRITE(line,1100) axis_ligand(3)+0.1,axis_ligand(2)+0.01,axis_ligand(1)+0.001, surf_normal
               ELSE
                  WRITE(line,1100) axis_ligand, surf_normal
               ENDIF
               laenge = 81
               CALL vprod(line, laenge)                      ! Make rotation axis
               sym_uvw(:) = res_para(1:3)
               CALL trans (sym_uvw, cr_gten, sym_hkl, 3)     ! Make reciprocal space axis
               CALL symm_setup                               ! Symmetry setup defines matrix
!              CALL symm_show                                ! Show only in debug
               CALL symm_op_single                           ! Perform the operation
            ENDIF
            IF(ABS(dist) < EPS ) THEN                      ! Remove surface atom
               cr_iscat(ia) = 0
               cr_prop (ia) = ibclr (cr_prop (ia), PROP_NORMAL)
            ENDIF
         ENDIF
      ENDDO moles
!   ENDIF
!
   chem_period(:) = .false.                         ! We inserted atoms, turn off periodic boundaries
   chem_quick     = .false.                         ! turn of quick search
   1000 FORMAT(a4,4(2x,',',F12.6))
   1100 FORMAT(6(F12.6,', '),'ddd')
   END SUBROUTINE deco_place_normal
!
!*******************************************************************************
!
   SUBROUTINE deco_place_bridge(dc_def_temp, ia, ia_scat, xyz, &
                            mole_axis, mole_name, mole_length, &
                            surf, neig, dist)
!
   USE atom_env_mod
   USE chem_mod
   USE metric_mod
   USE modify_mod
   USE symm_menu
   USE symm_mod
   USE symm_sup_mod
   USE trafo_mod
!
   USE param_mod
!
   IMPLICIT NONE
!
   TYPE (dc_def), POINTER              :: dc_def_temp     ! The definition to be used
   INTEGER,                 INTENT(IN) :: ia              ! Surface atom number
   INTEGER,                 INTENT(IN) :: ia_scat         ! Scattering type of surface atom
   REAL   , DIMENSION(1:3), INTENT(IN) :: xyz             ! Position of surface atom
   INTEGER, DIMENSION(1:2), INTENT(IN) :: mole_axis       ! Atoms that define molecule axis
   CHARACTER (LEN=1024),    INTENT(IN) :: mole_name       ! Molecule file name
   INTEGER,                 INTENT(IN) :: mole_length     ! Molecule file name length
   INTEGER, DIMENSION(0:4), INTENT(IN) :: surf            ! Surface atom type
   INTEGER,                 INTENT(IN) :: neig            ! Connected to this neighbor in mole
   REAL   ,                 INTENT(IN) :: dist            ! distance to ligand molecule
!
   REAL, PARAMETER :: EPS = 1.0E-6
   INTEGER, PARAMETER                      :: MAXW = 2
   REAL                , DIMENSION(1:MAXW) :: werte
!
   CHARACTER (LEN=4) :: atom_name
   CHARACTER (LEN=1024)                    :: line
   INTEGER                                 :: ianz
   INTEGER                                 :: i,j, im, laenge
   INTEGER                                 :: iprop
   INTEGER                                 :: itype
   LOGICAL  , DIMENSION(1:3)               :: fp
   LOGICAL                                 :: fq
   LOGICAL, PARAMETER :: lspace = .true.
   REAL   , DIMENSION(1:3) :: axis_ligand                 ! Initial molecule orientation
   REAL                                    :: rmin, radius, normal_l, dw1, b_l, b_n
   REAL     , DIMENSION(1:3)               :: x, bridge, tangent, origin, posit
   REAL     , DIMENSION(1:3)               :: surf_normal
   REAL     , DIMENSION(1:3)               :: vnull
!
   vnull(:) = 0.00
!
!  Find the other partners involved in the bridge 
   x(1)     = cr_pos(1,ia)
   x(2)     = cr_pos(2,ia)
   x(3)     = cr_pos(3,ia)
   rmin     = 0.1
   radius   = dist*2.0
   ianz     = 1
   werte(1) = surf(2)
   fp (:)   = chem_period (:)
   fq       = chem_quick
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
         RETURN
      ENDIF
      bridge(1) = (cr_pos(1,ia)-cr_pos(1,atom_env(j)))   ! Calculate vector along bridge
      bridge(2) = (cr_pos(2,ia)-cr_pos(2,atom_env(j)))
      bridge(3) = (cr_pos(3,ia)-cr_pos(3,atom_env(j)))
      x(1) = (cr_pos(1,ia)+cr_pos(1,atom_env(j)))*0.5    ! Calculate midpoint
      x(2) = (cr_pos(2,ia)+cr_pos(2,atom_env(j)))*0.5
      x(3) = (cr_pos(3,ia)+cr_pos(3,atom_env(j)))*0.5
      WRITE(line,1100) x, bridge                         ! Calculate vector parallel to surface
      laenge = 81
      CALL vprod(line, laenge)
      tangent(:) = res_para(1:3)
      WRITE(line,1100) bridge, tangent                   ! Calculate surface normal
      laenge = 81
      CALL vprod(line, laenge)
      surf_normal(:) = res_para(1:3)
      moles: DO i=1, dc_n_molecules
         IF(mole_name(1:mole_length) == m_name(i)(1:m_lname(i))) THEN
            im = mole_axis(2)
            CALL dc_molecules(i)%get_cryst_atom(im, itype, posit, iprop)
            axis_ligand(1) = posit(1)
            axis_ligand(2) = posit(2)
            axis_ligand(3) = posit(3)
            im = mole_axis(1)
            CALL dc_molecules(i)%get_cryst_atom(im, itype, posit, iprop)
            axis_ligand(1) = axis_ligand(1) - posit(1)
            axis_ligand(2) = axis_ligand(2) - posit(2)
            axis_ligand(3) = axis_ligand(3) - posit(3)
!           Insert molecule atoms into crystal at correct origin in initial orientation
            normal_l = sqrt (skalpro (surf_normal, surf_normal, cr_gten))
            b_l = sqrt (skalpro (bridge, bridge, cr_gten))     ! Calculate bridge length
            b_n = sqrt(dist**2-b_l**2/4.)                      ! Calculate distance along normal
            origin(1)  = x(1) + surf_normal(1)/normal_l*b_n    ! Calculate ligand origin
            origin(2)  = x(2) + surf_normal(2)/normal_l*b_n
            origin(3)  = x(3) + surf_normal(3)/normal_l*b_n
            sym_latom(:) = .false.                        ! Initially deselect all atomtypes
            DO im=1,m_length(i)                           ! Insert all atoms
               CALL dc_molecules(i)%get_cryst_atom(im, itype, posit, iprop)
               CALL dc_molecules(i)%get_cryst_scat(im, itype, atom_name, dw1)
               posit(:) = posit(:) + origin(:)
               WRITE(line, 1000) atom_name, posit, dw1
               laenge = 60
               CALL do_ins(line, laenge)
               CALL check_symm
               sym_latom(cr_iscat(cr_natoms)) = .true.    ! Select atopm type for rotation
            ENDDO
! define rotation operation
            sym_angle      = do_bang(lspace, surf_normal, vnull, axis_ligand)
            IF(ABS(sym_angle) > EPS ) THEN                ! Rotate if not zero degrees
               sym_orig(:)    = origin(:)                 ! Define origin
               sym_trans(:)   = 0.0                       ! No translation needed
               sym_sel_atom   = .true.                    ! Select atoms
               sym_new        = .false.                   ! No new types
               sym_power      =  1                        ! Just need one operation
               sym_type       = .true.                    ! Proper rotation
               sym_mode       = .false.                   ! Move atom to new position
               sym_orig_mol   = .false.                   ! Origin at crystal
               sym_power_mult =.false.                    ! No multiple copies
               sym_sel_atom   = .true.                    ! Select atoms not molecules
               sym_start      =  cr_natoms - m_length(i) + 1 ! set range of atoms numbers
               sym_end        =  cr_natoms
               IF(ABS(sym_angle-180.) < EPS ) THEN        ! Ligand and surface normal are antiparallel
                  WRITE(line,1100) axis_ligand(3)+0.1,axis_ligand(2)+0.01,axis_ligand(1)+0.001, surf_normal
               ELSE
                  WRITE(line,1100) axis_ligand, surf_normal
               ENDIF
               laenge = 81
               CALL vprod(line, laenge)                   ! Make rotation axis
               sym_uvw(:) = res_para(1:3)
               CALL trans (sym_uvw, cr_gten, sym_hkl, 3)
               CALL symm_setup
!              CALL symm_show
               CALL symm_op_single
            ENDIF
         ENDIF
      ENDDO moles
   ENDIF
!
   chem_period(:) = .false.                         ! We inserted atoms, turn off periodic boundaries
   chem_quick     = .false.                         ! turn of quick search
!
1000 FORMAT(a4,4(2x,',',F12.6))
1100 FORMAT(6(F12.6,', '),'ddd')
!
   END SUBROUTINE deco_place_bridge
!
!*******************************************************************************
!
   SUBROUTINE deco_place_double(dc_def_temp, ia, ia_scat, xyz, &
                            mole_axis, mole_name, mole_length, &
                            surf, neig, dist,ncon)
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
   USE atom_env_mod
   USE chem_mod
   USE metric_mod
   USE modify_mod
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
   TYPE (dc_def), POINTER              :: dc_def_temp     ! The definition to be used
   INTEGER,                 INTENT(IN) :: ia              ! Surface atom number
   INTEGER,                 INTENT(IN) :: ia_scat         ! Scattering type of surface atom
   REAL   , DIMENSION(1:3), INTENT(IN) :: xyz             ! Position of surface atom
   INTEGER, DIMENSION(1:2), INTENT(IN) :: mole_axis       ! Atoms that define molecule axis
   CHARACTER (LEN=1024),    INTENT(IN) :: mole_name       ! Molecule file name
   INTEGER,                 INTENT(IN) :: mole_length     ! Molecule file name length
   INTEGER, DIMENSION(0:4), INTENT(IN) :: surf            ! Surface atom type
   INTEGER,                 INTENT(IN) :: neig            ! Connected to this neighbor in mole
   REAL   ,                 INTENT(IN) :: dist            ! distance to ligand molecule
   INTEGER,                 INTENT(IN) :: ncon            ! Number of defined bonds
!
   TYPE (dc_con), POINTER              :: dc_con_temp     ! A connectivity definition to be used
   REAL, PARAMETER :: EPS = 1.0E-6
   INTEGER, PARAMETER                      :: MAXW = 2
   REAL                , DIMENSION(1:MAXW) :: werte
!
   CHARACTER (LEN=4) :: atom_name
   CHARACTER (LEN=1024)                    :: line
   INTEGER   :: surf_char ! Surface character, plane, edge, corner, ...
   INTEGER, DIMENSION(0:4)             :: surface         ! Surface atom type
   INTEGER                             :: neighbor        ! Connected to this neighbor in mole
   REAL                                :: distance        ! distance to ligand molecule
   INTEGER, DIMENSION(:), ALLOCATABLE  :: all_surface         ! Surface atom type
   INTEGER, DIMENSION(:), ALLOCATABLE  :: all_neighbor        ! Connected to this neighbor in mole
   REAL   , DIMENSION(:), ALLOCATABLE  :: all_distance        ! distance to ligand molecule
   INTEGER                                 :: ianz
   INTEGER                                 :: i,j, l,im, laenge
   INTEGER                                 :: iprop
   INTEGER                                 :: itype
   INTEGER                                 :: n_atoms_orig   ! Number of atoms prior to insertion
   INTEGER                                 :: n1,n2          ! number of mol neighbours after insertion
   INTEGER                                 :: a1,a2          ! number of mol axis atoms after rotations
   INTEGER                                 :: success        ! Everything went fine
   LOGICAL  , DIMENSION(1:3)               :: fp
   LOGICAL                                 :: fq
   LOGICAL, PARAMETER :: lspace = .true.
   REAL                                    :: rmin, radius, dw1, b_l, t_l
   REAL                                    :: alpha, beta
   REAL     , DIMENSION(1:3)               :: x, bridge, tangent, origin, posit, v, w, u
   REAL     , DIMENSION(1:3)               :: shift, v1, v2, v3
   REAL     , DIMENSION(1:3)               :: surf_normal
   REAL     , DIMENSION(1:3)               :: vnull
!
   vnull(:) = 0.00
   success = -1
!
!  Load the molecule into the crystal structure
!
   n_atoms_orig = cr_natoms                         ! Number of atoms prior to insertion
   moles: DO i=1, dc_n_molecules
      IF(mole_name(1:mole_length) == m_name(i)(1:m_lname(i))) THEN
      im = mole_axis(2)
      CALL dc_molecules(i)%get_cryst_atom(im, itype, posit, iprop)
      origin(:) = 0.0                               ! initially place at 0,0,0
      sym_latom(:) = .false.                        ! Initially deselect all atomtypes
      insert: DO im=1,m_length(i)                   ! Insert all atoms
         CALL dc_molecules(i)%get_cryst_atom(im, itype, posit, iprop)
         CALL dc_molecules(i)%get_cryst_scat(im, itype, atom_name, dw1)
         posit(:) = posit(:) + origin(:)
         WRITE(line, 1000) atom_name, posit, dw1
         laenge = 60
         CALL do_ins(line, laenge)
         CALL check_symm
         sym_latom(cr_iscat(cr_natoms)) = .true.    ! Select atopm type for rotation
      ENDDO insert
      ENDIF
   ENDDO moles
!
   chem_period(:) = .false.                         ! We inserted atoms, turn off periodic boundaries
   chem_quick     = .false.                         ! turn of quick search
!
   ALLOCATE(all_surface(1:ncon))
   ALLOCATE(all_neighbor(1:ncon))
   ALLOCATE(all_distance(1:ncon))
!
   all_surface (1) = ia
   all_neighbor(1) = neig
   all_distance(1) = dist
!
!  FIND the other surface partners involved in the bonds.
!
   x(1)     = cr_pos(1,ia)
   x(2)     = cr_pos(2,ia)
   x(3)     = cr_pos(3,ia)
   dc_con_temp => dc_def_temp%dc_def_con%next   ! Point to second connectivity
   search: DO l=2,ncon
      IF(.NOT. ASSOCIATED(dc_con_temp) ) GOTO 9999
      CALL dc_get_con(dc_con_temp, surface, neighbor, distance)
      n2 = n_atoms_orig +     neighbor   
      n1 = n_atoms_orig + all_neighbor(1)
      bridge(1) = cr_pos(1,n2) - cr_pos(1,n1)
      bridge(2) = cr_pos(2,n2) - cr_pos(2,n1)
      bridge(3) = cr_pos(3,n2) - cr_pos(3,n1)
      b_l      = sqrt(skalpro(bridge, bridge, cr_gten))
      rmin     = MAX( 0.1, b_l - dist - distance )          ! Minimum distance between surface atoms
      radius   = b_l + dist + distance                      ! Maximum distance between surface atoms
      ianz     = 1
      werte(1) = surface(1)
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
         all_neighbor(l) = neighbor
         all_distance(l) = distance
      ELSE
         GOTO 9999
      ENDIF  ! 
   ENDDO search
!  Find surface character and local normal
   CALL find_surface_character(ia, surf_char, surf_normal)
!  Determine rotation axis for surface vector
   tangent(:) = cr_pos(:,all_surface (2)) - cr_pos(:,ia)  ! Vector between surface atoms
   t_l        = sqrt(skalpro(tangent, tangent, cr_gten))  ! Distance between surface atoms
   WRITE(*   ,1100) tangent, surf_normal                  ! Calculate rotation axis
   WRITE(line,1100) tangent, surf_normal                  ! Calculate rotation axis
   laenge = 81
   CALL vprod(line, laenge)
   sym_uvw(:) = res_para(1:3)
!  Calculate angle in first trapezoid corner
   sym_angle  = ACOS(-(all_distance(2)**2-all_distance(1)**2-(t_l-b_l)**2)/ &
                      (2.*all_distance(1)*(t_l-b_l))) /rad
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
   v(:) =  tangent(:)*all_distance(1)/t_l     ! Scale vector to bond length
   CALL symm_ca_single (v, .true., .false.)
   origin(:) = cr_pos(:,all_surface(1)) + res_para(1:3)  ! Origin of the molecule = surface 1 + result 
   shift (:) = origin(:) - cr_pos(:,n1)       ! All molecule atoms need to be shifted by this vector
   DO i=n_atoms_orig+1,cr_natoms
      cr_pos(:,i) = cr_pos(:,i) + shift(:)
   ENDDO
!
!  Calculate angle in second trapezoid corner
   sym_angle  = ACOS(-(all_distance(1)**2-all_distance(2)**2-(t_l-b_l)**2)/ &
                      (2.*all_distance(2)*(t_l-b_l))) /rad
!
   sym_trans(:)   = 0.0                       ! No translation needed
   sym_orig(:)    = 0.0                       ! Define origin in 0,0,0
   sym_uvw(:)     = -sym_uvw(:)               ! invert axis
   v(:) =  -tangent(:)*all_distance(2)/t_l    ! invert and scale vector to bond length
   CALL trans (sym_uvw, cr_gten, sym_hkl, 3)  ! Make reciprocal space axis
   CALL symm_setup
   CALL symm_ca_single (v, .true., .false.)   ! rotate negative surface vector
   w(:) = cr_pos(:,all_surface(2)) + res_para(1:3)  ! Target position for 2nd molecule atom
!
   sym_orig(:)    = cr_pos(:,n1)              ! Define origin in 1st attached molecule atom
   v1(:) = cr_pos(:,n2) - cr_pos(:,n1)        ! Current vector from 1st to 2nd molecule atom
   v2(:) = w(:)         - cr_pos(:,n1)        ! Vector from 1st to target  2nd molecule atom
   WRITE(line,1100) v1, v2                    ! Rotation axis will be v1 x v2
   laenge = 81
   CALL vprod(line, laenge)
   sym_uvw(:) = res_para(1:3)
   sym_trans(:)   = 0.0                       ! No translation needed
   sym_angle  = do_bang(lspace, v1, vnull, v2)  ! Calculate rotation angle = < (v1,v2)
   sym_start  =  n_atoms_orig + 1             ! set range of atoms numbers
   sym_end    =  cr_natoms
   CALL trans (sym_uvw, cr_gten, sym_hkl, 3)  ! Make reciprocal space axis
   CALL symm_setup
   CALL symm_op_single                           ! Perform the operation
!
!   Rotate molecule up to straighten molecule axis out
   sym_uvw(:) = cr_pos(:,n2) - cr_pos(:,n1)       ! Rotation axis
   CALL dc_get_axis(dc_def_temp, dc_temp_axis)
   a1 = n_atoms_orig + dc_temp_axis(1)
   a2 = n_atoms_orig + dc_temp_axis(2)
   v1(:)      = cr_pos(:,a2) - cr_pos(:,a1)      ! Current molecule axis
   WRITE(line,1200) v1, sym_uvw                  ! First project molecule axis into 
   laenge = 82                                   !   plane normal to the 
   CALL do_proj(line, laenge)                    !   vector between connected molecule atoms
   v3(:) = res_para(4:6)                         ! This is the projection
   WRITE(line,1100) sym_uvw, surf_normal         ! Find normal to plane defined by
   laenge = 81                                   !   vector between connected molecule atoms
   CALL vprod(line, laenge)                      !   and surface normal
   w(:) = res_para(1:3)                          ! Need to project (projected) mol axis into plane normal to w
   WRITE(line,1200) v3, w                        ! Prepare projection
   laenge = 82
   CALL do_proj(line, laenge)                    ! Project axis into plane 
   v2(:) = res_para(4:6)                         ! This is the projection
   alpha      = do_bang(lspace, surf_normal, vnull, v2)   ! Calculate angle normal and projection
   WRITE(line,1100) v3,v2                        ! Do vector product (mol_axis) x (projection)
   laenge = 81
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
   success = 0 

!
9999 CONTINUE                                             ! Jump here from errors to ensure dealloc
   DEALLOCATE(all_surface)
   DEALLOCATE(all_neighbor)
   DEALLOCATE(all_distance)
   IF(success /=0) THEN                          ! An error occurred, reset crystal
      cr_natoms = n_atoms_orig
   ENDIF
!
1000 FORMAT(a4,4(2x,',',F12.6))
1100 FORMAT(6(F12.6,', '),'ddd')
1200 FORMAT(6(F12.6,', '),'dddd')
!
   END SUBROUTINE deco_place_double
!
!*****7*****************************************************************
!
   SUBROUTINE deco_place_multi(dc_def_temp, ia, ia_scat, xyz, &
                            mole_axis, mole_name, mole_length, &
                            surf, neig, dist,ncon)
!
!  Places a molecule that has multiple bonds to the surface.
!  The first bond should be the one that carries multiple connections to
!  the surface. 
!  The first molecule atom is in a uniquely specified position.
!  If further bonds are specified, the molecule is rotated around this
!  first anchor point to fulfill as best as possible the further conditions.
!
   USE atom_env_mod
   USE chem_mod
   USE metric_mod
   USE modify_mod
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
   TYPE (dc_def), POINTER              :: dc_def_temp     ! The definition to be used
   INTEGER,                 INTENT(IN) :: ia              ! Surface atom number
   INTEGER,                 INTENT(IN) :: ia_scat         ! Scattering type of surface atom
   REAL   , DIMENSION(1:3), INTENT(IN) :: xyz             ! Position of surface atom
   INTEGER, DIMENSION(1:2), INTENT(IN) :: mole_axis       ! Atoms that define molecule axis
   CHARACTER (LEN=1024),    INTENT(IN) :: mole_name       ! Molecule file name
   INTEGER,                 INTENT(IN) :: mole_length     ! Molecule file name length
   INTEGER, DIMENSION(0:4), INTENT(IN) :: surf            ! Surface atom type
   INTEGER,                 INTENT(IN) :: neig            ! Connected to this neighbor in mole
   REAL   ,                 INTENT(IN) :: dist            ! distance to ligand molecule
   INTEGER,                 INTENT(IN) :: ncon            ! Number of defined bonds
!
   TYPE (dc_con), POINTER              :: dc_con_temp     ! A connectivity definition to be used
   REAL, PARAMETER :: EPS = 1.0E-6
   INTEGER, PARAMETER                      :: MAXW = 2
   REAL                , DIMENSION(1:MAXW) :: werte
!
   CHARACTER (LEN=4) :: atom_name
   CHARACTER (LEN=1024)                    :: line
   INTEGER   :: surf_char ! Surface character, plane, edge, corner, ...
   INTEGER, DIMENSION(0:4)             :: surface         ! Surface atom type
   INTEGER                             :: neighbor        ! Connected to this neighbor in mole
   REAL                                :: distance        ! distance to ligand molecule
   INTEGER, DIMENSION(:), ALLOCATABLE  :: all_surface         ! Surface atom type
   INTEGER, DIMENSION(:), ALLOCATABLE  :: all_neighbor        ! Connected to this neighbor in mole
   REAL   , DIMENSION(:), ALLOCATABLE  :: all_distance        ! distance to ligand molecule
   INTEGER                                 :: ianz
   INTEGER                                 :: i,j, l,im, laenge
   INTEGER                                 :: iprop
   INTEGER                                 :: itype
   INTEGER                                 :: n_atoms_orig   ! Number of atoms prior to insertion
   INTEGER                                 :: n1,n2          ! number of mol neighbours after insertion
   INTEGER                                 :: a1,a2          ! number of mol axis atoms after rotations
   INTEGER                                 :: success        ! Everything went fine
   LOGICAL  , DIMENSION(1:3)               :: fp
   LOGICAL                                 :: fq
   LOGICAL, PARAMETER :: lspace = .true.
   REAL                                    :: rmin, radius, dw1, b_l
   REAL                                    :: alpha, beta, d1,d2, v_l
   REAL     , DIMENSION(1:3)               :: x, bridge, base, origin, posit, v, w, u
   REAL     , DIMENSION(1:3)               :: shift, v1, v2, v3
   REAL     , DIMENSION(1:3)               :: surf_normal
   REAL     , DIMENSION(1:3)               :: vnull
!
!
   vnull(:) = 0.00
   success = -1
!
!  Load the molecule into the crystal structure
!
   n_atoms_orig = cr_natoms                         ! Number of atoms prior to insertion
   moles: DO i=1, dc_n_molecules
      IF(mole_name(1:mole_length) == m_name(i)(1:m_lname(i))) THEN
      im = mole_axis(2)
      CALL dc_molecules(i)%get_cryst_atom(im, itype, posit, iprop)
      origin(:) = 0.0                               ! initially place at 0,0,0
      sym_latom(:) = .false.                        ! Initially deselect all atomtypes
      insert: DO im=1,m_length(i)                   ! Insert all atoms
         CALL dc_molecules(i)%get_cryst_atom(im, itype, posit, iprop)
         CALL dc_molecules(i)%get_cryst_scat(im, itype, atom_name, dw1)
         posit(:) = posit(:) + origin(:)
         WRITE(line, 1000) atom_name, posit, dw1
         laenge = 60
         CALL do_ins(line, laenge)
         CALL check_symm
         sym_latom(cr_iscat(cr_natoms)) = .true.    ! Select atom type for rotation
      ENDDO insert
      ENDIF
   ENDDO moles
!
   chem_period(:) = .false.                         ! We inserted atoms, turn off periodic boundaries
   chem_quick     = .false.                         ! turn of quick search
!
   ALLOCATE(all_surface(1:ncon))
   ALLOCATE(all_neighbor(1:ncon))
   ALLOCATE(all_distance(1:ncon))
!
   all_surface (1) = ia
   all_neighbor(1) = neig
   all_distance(1) = dist
!
!  Find surface character and local normal
   CALL find_surface_character(ia, surf_char, surf_normal)
!  Find the other surface atoms involved in this bond
   dc_con_temp => dc_def_temp%dc_def_con
   CALL dc_get_con(dc_con_temp, surface, neighbor, distance)
   n1 = n_atoms_orig +     neighbor   
   CALL deco_find_anchor(surface(0), surface, distance, ia, surf_normal, posit, base)
!
!   Move molecule to anchor position
!
   shift (:) = posit(:) - cr_pos(:,n1)
   DO i=n_atoms_orig+1,cr_natoms
      cr_pos(:,i) = cr_pos(:,i) + shift(:)
   ENDDO
   success = 0
!
   IF(ncon == 2) THEN                            ! We have the second connection
      dc_con_temp => dc_def_temp%dc_def_con%next   ! Point to second connectivity
      l = 2
!  search: DO l=2,ncon
      IF(.NOT. ASSOCIATED(dc_con_temp) ) GOTO 9999
      surface = 0
      neighbor= 0
      distance= 0.0
      CALL dc_get_con(dc_con_temp, surface, neighbor, distance)
      n2 = n_atoms_orig +     neighbor   
      bridge(1) = cr_pos(1,n2) - cr_pos(1,n1)
      bridge(2) = cr_pos(2,n2) - cr_pos(2,n1)
      bridge(3) = cr_pos(3,n2) - cr_pos(3,n1)
      b_l      = sqrt(skalpro(bridge, bridge, cr_gten))
      rmin     = MAX( 0.0, b_l - dist - distance )          ! Minimum distance between surface atoms
      radius   = b_l + dist + distance                      ! Maximum distance between surface atoms
      ianz     = 1
      werte(1) = surface(1)
      x(:)     = cr_pos(:,ia)                               ! Search around 1.st surface atom
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
         all_neighbor(l) = neighbor
         all_distance(l) = distance
      ELSE
         GOTO 9999
      ENDIF  ! 
!   ENDDO search
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
   laenge = 81
   CALL vprod(line, laenge)
   sym_uvw(:) =  res_para(1:3)
   d1 = all_distance(1)
   d2 = all_distance(2)
   sym_angle  = acos( (d1**2 + d2**2 - b_l**2)/(2.*d1*d2))/rad
   v(:) = v(:) *d2/v_l                        ! Scale vector 1st surface to 1st mole to distance2
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
   posit(:) = cr_pos(:,ia) + res_para(1:3)    ! Add rotated vector to 1st surface
!
! next step rotate molecule for 2nd mole to fall onto target posit
   u(:) = posit(:) - cr_pos(:,n1)             ! Vector from 1st mole to target
   WRITE(line,1100) bridge, u                 ! Do vector product (1st to 2nd mole) x (1st mole to target)
   laenge = 81
   CALL vprod(line, laenge)
   sym_uvw(:) =  res_para(1:3)
   sym_angle  = do_bang(lspace, bridge, vnull, u) ! Angle (1st to 2nd mole) and (1st mole to target)
   sym_orig(:) = cr_pos(:,n1)                 ! Set origin in 1st mole
   sym_start  =  n_atoms_orig + 1             ! set range of atoms numbers
   sym_end    =  cr_natoms
   CALL trans (sym_uvw, cr_gten, sym_hkl, 3)  ! Make reciprocal space axis
   CALL symm_setup
   CALL symm_op_single                        ! Perform the operation
!
!   Rotate molecule up to straighten molecule axis out
   sym_uvw(:) = cr_pos(:,n2) - cr_pos(:,n1)       ! Rotation axis
   CALL dc_get_axis(dc_def_temp, dc_temp_axis)
   a1 = n_atoms_orig + dc_temp_axis(1)
   a2 = n_atoms_orig + dc_temp_axis(2)
   v1(:)      = cr_pos(:,a2) - cr_pos(:,a1)      ! Current molecule axis
   WRITE(line,1200) v1, sym_uvw                  ! First project molecule axis into 
   laenge = 82                                   !   plane normal to the 
   CALL do_proj(line, laenge)                    !   vector between connected molecule atoms
   v3(:) = res_para(4:6)                         ! This is the projection
   WRITE(line,1100) sym_uvw, surf_normal         ! Find normal to plane defined by
   laenge = 81                                   !   vector between connected molecule atoms
   CALL vprod(line, laenge)                      !   and surface normal
   w(:) = res_para(1:3)                          ! Need to project (projected) mol axis into plane normal to w
   WRITE(line,1200) v3, w                        ! Prepare projection
   laenge = 82
   CALL do_proj(line, laenge)                    ! Project axis into plane 
   v2(:) = res_para(4:6)                         ! This is the projection
   alpha      = do_bang(lspace, surf_normal, vnull, v2)   ! Calculate angle normal and projection
   WRITE(line,1100) v3,v2                        ! Do vector product (mol_axis) x (projection)
   laenge = 81
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
!
9999 CONTINUE                                 ! Jump here from errors to ensure dealloc
   DEALLOCATE(all_surface)
   DEALLOCATE(all_neighbor)
   DEALLOCATE(all_distance)
   IF(success /=0) THEN                       ! An error occurred, reset crystal
      cr_natoms = n_atoms_orig
   ENDIF
!
1000 FORMAT(a4,4(2x,',',F12.6))
1100 FORMAT(6(F12.6,', '),'ddd')
1200 FORMAT(6(F12.6,', '),'dddd')
!
   END SUBROUTINE deco_place_multi
!
!*****7*****************************************************************
   SUBROUTINE deco_find_anchor(MAXAT,surface, distance, ia, normal, posit, base)
!-                                                                      
!  Find a common point around MAXAT atom types in surface
!  ia is the first surface atom number
!
   USE atom_env_mod
   USE chem_mod
   USE metric_mod
   USE modify_mod
!
   USE param_mod
   USE wink_mod
!
   IMPLICIT NONE
   INTEGER                    , INTENT(IN)  :: MAXAT
   INTEGER, DIMENSION(0:MAXAT), INTENT(IN)  :: surface
   REAL                       , INTENT(IN)  :: distance
   INTEGER                    , INTENT(IN)  :: ia
   REAL   , DIMENSION(1:3)    , INTENT(IN)  :: normal
   REAL   , DIMENSION(1:3)    , INTENT(OUT) :: posit
   REAL   , DIMENSION(1:3)    , INTENT(OUT) :: base
!
   INTEGER, PARAMETER         :: MAXW = 3
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
   REAL   , DIMENSION(1:MAXW) :: werte
   REAL   , DIMENSION(1:3)    :: x, u,v,w, e1,e2,e3
   REAL                    :: u_l, v_l, w_l    ! length of vectors in triangle
   REAL                    :: av, sig, av_min, sig_min ! average length  and sigma
   REAL                    :: tx,ty, tz        ! Cartesion coordinates of target position
   REAL                    :: g2x, g2y         ! Cartesion coordinates of atom 3
   REAL     , DIMENSION(1:3)               :: vnull
!
   vnull(:) = 0.00
!
   neig(:,:) = 0
   rmin     = 0.0                  ! Minimum distance between surface atoms
   radius   = 2.0 * distance       ! Maximum distance between surface atoms
   ianz     = 1
   x(:)     = cr_pos(:, ia)        ! Seach around atom ia
   fp (:)   = chem_period (:)
   fq       = chem_quick
   find: DO l=2, MAXAT
      werte(1) = surface(l)           ! Find this atom type
      CALL do_find_env (ianz, werte, maxw, x, rmin, radius, fq, fp)  ! Find all neighbors
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
               sig  = SQRT((u_l-av)**2+(v_l-av)**2+(w_l-av)**2)/3
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
   e1(:) = u (:) / u_l                    ! Normalize to 1 angstroem
   v (:) = cr_pos(:,good3) - cr_pos(:,good1)  ! Temporary vector
   WRITE(line,1100) e1,v                         ! Do vector product (e1) x (atom good 3)
   laenge = 81
   CALL vprod(line, laenge)
   e3(:) =  res_para(1:3)                 ! Result is cartesian z-axis
   v_l  = SQRT(skalpro(e3,e3,cr_gten))
   IF(do_bang(lspace, e3, vnull, normal) <= 90) THEN
      e3(:) =  e3(:) / v_l                ! Normalize to 1 angstroem
   ELSE
      e3(:) = -e3(:) / v_l                ! invert and Normalize to 1 angstroem
   ENDIF
   WRITE(line,1100) e3,e1                 ! Do vector product (e3) x (e1 )
   laenge = 81
   CALL vprod(line, laenge)
   e2(:) =  res_para(1:3)                 ! Result is cartesian y-axis
   v_l  = SQRT(skalpro(e2,e2,cr_gten))    ! cartesian x-coordinate of atom 2
   e2(:) = e2(:) / v_l                    ! Normalize to 1 angstroem
   g2x =     (skalpro(v,e1,cr_gten))      ! cartesian x-coordinate of atom 3
   g2y =     (skalpro(v,e2,cr_gten))      ! cartesian y-coordinate of atom 3
!  Calculate target coordinates from trilateration
   tx = 0.5 * u_l
   ty = 0.5 * (g2x**2+g2y**2)/g2y - g2x/g2y*tx
   tz = sqrt(distance**2-tx**2 - ty**2)
   posit(:) = cr_pos(:,good1) + tx*e1(:) + ty*e2(:) + tz*e3(:)
   
   base(:) = (cr_pos(:,neig(lgood,2))+ cr_pos(:,neig(kgood,3)))*0.5
   ENDIF
!
1100 FORMAT(6(F12.6,', '),'ddd')
!
   END SUBROUTINE deco_find_anchor
!
!*****7*****************************************************************
!
      SUBROUTINE trans_atoms_tocart (uvw_out)
!-                                                                      
!     transforms atom coordinates into a cartesian space                
!     Warning, only the fractional coordinates are transformed,         
!     the unit cell and space group information is not touched.         
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE discus_plot_mod 
      USE trans_sup_mod
      IMPLICIT none 
!                                                                       
      REAL ,DIMENSION(1:3), INTENT(OUT) :: uvw_out !(3)
!
      INTEGER              ::  i
      LOGICAL, PARAMETER   :: lscreen = .false. 
      REAL, DIMENSION(1:4) :: uvw
      REAL             :: xmin
      REAL             :: xmax
      REAL             :: ymin
      REAL             :: ymax
      REAL             :: zmin
      REAL             :: zmax
!                                                                       
      xmin = 0.0
      xmax = 0.0
      ymin = 0.0
      ymax = 0.0
      zmin = 0.0
      zmax = 0.0
      uvw(4) = 1.0
!         
      DO i = 1, cr_natoms 
         uvw (1) = cr_pos (1, i) 
         uvw (2) = cr_pos (2, i) 
         uvw (3) = cr_pos (3, i) 
         CALL tran_ca (uvw, pl_tran_f, lscreen) 
         cr_pos (1, i) = uvw (1) 
         cr_pos (2, i) = uvw (2) 
         cr_pos (3, i) = uvw (3) 
         xmin = MIN(xmin,uvw(1))
         xmax = MAX(xmax,uvw(1))
         ymin = MIN(ymin,uvw(2))
         ymax = MAX(ymax,uvw(2))
         zmin = MIN(zmin,uvw(3))
         zmax = MAX(zmax,uvw(3))
      ENDDO
      uvw_out (1) = ABS(xmax-xmin)
      uvw_out (2) = ABS(ymax-ymin)
      uvw_out (3) = ABS(zmax-zmin) 
!                                                                       
      END SUBROUTINE trans_atoms_tocart      
!*****7*****************************************************************
      SUBROUTINE trans_atoms_fromcart 
!-                                                                      
!     transforms atom coordinates from a cartesian space back           
!     to the original coordinates                                       
!     Warning, only the fractional coordinates are transformed,         
!     the unit cell and space group information is not touched.         
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE discus_plot_mod 
      USE trans_sup_mod
      IMPLICIT none 
!                                                                       
      INTEGER              :: i 
      LOGICAL, PARAMETER   :: lscreen = .false.
!                                                                       
      REAL, DIMENSION(1:4) ::  uvw !(4) 
!                                                                       
!                                                                       
      uvw(4) = 1.0
      DO i = 1, cr_natoms 
         uvw (1) = cr_pos (1, i) 
         uvw (2) = cr_pos (2, i) 
         uvw (3) = cr_pos (3, i) 
         CALL tran_ca (uvw, pl_tran_fi, lscreen) 
         cr_pos (1, i) = uvw (1) 
         cr_pos (2, i) = uvw (2) 
         cr_pos (3, i) = uvw (3) 
      ENDDO 
   END SUBROUTINE trans_atoms_fromcart
!
   SUBROUTINE check_symm
!
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
END MODULE mole_surf_mod
