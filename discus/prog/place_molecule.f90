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
!  CHARACTER (LEN=*)   :: line
   CHARACTER (LEN=5)                       :: befehl! command on input line
   CHARACTER (LEN=50)                      :: prom  ! Menu prompt
   CHARACTER (LEN=1024)                    :: line  ! input line
   CHARACTER (LEN=1024)                    :: zeile ! remainder with parameters
   INTEGER                                 :: indxg ! location of "="
   INTEGER                                 :: lp    ! lengtz of zeile
   INTEGER laenge, lbef
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
      dc_use_conn    = 1 ! while developing!!!
      dc_n_molecules = 1 ! while developing!!!
      dc_init = .false.
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
                  CALL deco_add (zeile, lp)
!                                                                       
!------ --Handle property settings 'property'                           
!                                                                       
              ELSEIF (str_comp (befehl, 'property', 4, lbef, 8) ) then
!                                                                       
                  CALL property_select (zeile, lp, dc_sel_prop)
!
              ELSEIF (str_comp (befehl, 'reset', 3, lbef, 4)) THEN
                  CALL deco_reset
!
              ELSEIF (str_comp (befehl, 'run', 3, lbef, 3)) THEN
                  CALL deco_run
!
              ELSEIF (str_comp (befehl, 'sel', 3, lbef, 3) .or.  &
                      str_comp (befehl, 'del', 3, lbef, 3)     ) THEN
                  CALL atom_select (zeile, lp, 0,  DC_MAXSCAT,  dc_latom, &
                  dc_sel_atom , lold  ,   &
                  str_comp (befehl, 'sel', 2, lbef, 3) )
!
              ELSEIF (str_comp (befehl, 'set', 3, lbef, 3)) THEN
                  CALL deco_set (zeile, lp)
!
!
              ELSEIF (str_comp (befehl, 'show', 3, lbef, 4)) THEN
                  CALL deco_show
                  WRITE(output_io,*)
!
!                  CALL find_surface_character
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
   INTEGER              :: ianz, janz, success,i
   LOGICAL              :: lnew
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
         dc_temp_axis(1) = werte(1)
         dc_temp_axis(2) = werte(2)
         CALL dc_find_def(dc_def_head,dc_def_temp, dc_temp_lname, dc_temp_name,dc_temp_id,lnew,success)
         IF(success==0) CALL dc_set_axis(dc_def_temp, dc_temp_axis)
      ELSE
         ier_num = -6
         ier_typ = ER_COMM
      ENDIF
   ELSEIF ( str_comp(cpara(2),'bond',4,lpara(2),4) ) THEN
      IF ( ianz >= 5 ) THEN
         CALL del_params (2, ianz, cpara, lpara, maxw)   ! delete first 2 params
         ccpara(1) = cpara(ianz-1)
         ccpara(2) = cpara(ianz  )
         llpara(1) = lpara(ianz-1)
         llpara(2) = lpara(ianz  )
         janz = 2
         CALL ber_params (janz, ccpara, llpara, wwerte, 2)
         dc_temp_neig  = NINT(wwerte(1))
         dc_temp_dist  =      wwerte(2)
         cpara(ianz-1) = ' '
         cpara(ianz  ) = ' '
         lpara(ianz-1) = 1
         lpara(ianz  ) = 1
         janz = ianz - 2                                 ! ignore last two parameters
         CALL get_iscat (janz, cpara, lpara, werte, maxw, .false.)
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
!            DO i=1, dc_temp_surf(0)
!               dc_latom(dc_temp_surf(i)) = .true.              ! Select this atom type
!            ENDDO
               dc_latom(dc_temp_surf(1)) = .true.              ! Select first atom type
         ENDIF
      ELSE
         ier_num = -6
         ier_typ = ER_COMM
      ENDIF
   ELSEIF ( str_comp(cpara(2),'ligand',4,lpara(2),6) ) THEN
      IF ( ianz == 3 ) THEN
         dc_temp_file  = cpara(3)
         dc_temp_lfile = lpara(3)
         dc_temp_id    = 0
         dc_def_temp => dc_def_head
!
         CALL dc_find_def(dc_def_head,dc_def_temp, dc_temp_lname, dc_temp_name,dc_temp_id,lnew,success)
         IF(success==0) CALL dc_set_file(dc_def_temp, dc_temp_lfile, dc_temp_file)
      ELSE
         ier_num = -6
         ier_typ = ER_COMM
      ENDIF
!   ELSEIF ( str_comp(cpara(2),'connection',4,lpara(2),10) ) THEN
!      IF ( ianz == 6 ) THEN
!         dc_temp_name  = cpara(2)
!         dc_temp_lname = lpara(2)
!         dc_temp_id    = 0
!         dc_def_temp => dc_def_head
!         CALL dc_find_def(dc_def_head,dc_def_temp, dc_temp_lname, dc_temp_name,dc_temp_id,lnew,success)
!!         CALL dc_set_connection()
!      ELSE
!         ier_num = -6
!         ier_typ = ER_COMM
!      ENDIF
!  ELSEIF ( str_comp(cpara(1),'type',4,lpara(1),4) ) THEN
!     IF ( ianz == 2 ) THEN
!        IF ( str_comp(cpara(3),'normal',4,lpara(1),6) ) THEN
!           dc_temp_type = DC_NORMAL
!        ELSE
!           ier_num = -6
!           ier_typ = ER_COMM
!           RETURN
!        ENDIF
!        dc_temp_name  = cpara(2)
!        dc_temp_lname = lpara(2)
!        dc_temp_id    = 0
!        NULLIFY(dc_def_temp)
!        lnew          = .false.
!        CALL dc_find_def(dc_def_head,dc_def_temp, dc_temp_lname, dc_temp_name,dc_temp_id,lnew,ier_typ)
!         CALL dc_set_type(dc_def_temp, dc_temp_type)
!     ELSE
!        ier_num = -6
!        ier_typ = ER_COMM
!     ENDIF
   ELSE
      ier_num = -6
      ier_typ = ER_COMM
   ENDIF
!
   IF(success == -1) THEN
      ier_num = -1116
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
      IF ( str_comp(cpara(2),'normal',4,lpara(1),6) ) THEN
         dc_temp_type = DC_NORMAL
      ELSEIF ( str_comp(cpara(2),'bridge',4,lpara(1),6) ) THEN
         dc_temp_type = DC_BRIDGE
      ELSE
         ier_num = -6
         ier_typ = ER_COMM
         RETURN
      ENDIF
      ALLOCATE(dc_def_temp)
      NULLIFY(dc_def_temp) ! => dc_def_head
      lnew          = .true.
      NULLIFY(dc_def_temp) ! => dc_def_head
      CALL dc_find_def(dc_def_head,dc_def_temp, dc_temp_lname, dc_temp_name,dc_temp_id,lnew,success)
      IF(success==0) THEN  ! set the type
         CALL dc_set_type(dc_def_temp, dc_temp_type)
      ELSE
         ier_num = -1117
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
   SUBROUTINE find_surface_character(istart, iend, surf_char, surf_normal)
!
   USE atom_env_mod
   USE chem_mod
   USE modify_mod
   USE tensors_mod
   USE prompt_mod
!
   IMPLICIT NONE
!
   INTEGER,INTENT(IN)  :: istart ! dummy index
   INTEGER,INTENT(IN)  :: iend   ! dummy index
   INTEGER,INTENT(OUT) :: surf_char   ! Surface character
   REAL   ,DIMENSION(1:3),INTENT(OUT) :: surf_normal   ! Surface normal

   INTEGER, PARAMETER :: MAXW = 1
   LOGICAL, PARAMETER :: LNEW = .false.
!
   CHARACTER  (LEN=4), DIMENSION(1:MAXW) :: cpara
   INTEGER  ,          DIMENSION(1:MAXW) :: lpara
   REAL     ,          DIMENSION(1:MAXW) :: werte
   INTEGER   :: ia ! dummy index
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
   atomloop: DO ia=istart,iend
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
        call do_find_env (ianz, werte, maxw, x, rmin,&
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
        write(output_io,*) ' row 1 ',direct(1,:), ' V ',vector(1)
        write(output_io,*) ' row 2 ',direct(2,:), ' V ',vector(2)
        write(output_io,*) ' row 3 ',direct(3,:), ' V ',vector(3)
        call invmat( recipr, direct )  ! calculate inverse matrix
        write(output_io,*) ' row 1 ',recipr(1,:)
        write(output_io,*) ' row 2 ',recipr(2,:)
        write(output_io,*) ' row 3 ',recipr(3,:)
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
        write(*,*) ' Normal ',finalv,' sigma ',sigma
        surf_normal(:) = finalv(:)
        surf_char = SURF_PLANE
        IF(sigma <    0.05) THEN
!           cr_iscat(ia) = cr_iscat(ia) + 4 !cr_nscat
           CYCLE atomloop
        ENDIF
        ENDIF
!
!       Atom is not in plane, try a row
!
     ENDIF
!
   ENDDO atomloop
!
   DEALLOCATE(deviat, STAT = istat)
!
   END SUBROUTINE find_surface_character
!
!######################################################################
!
   SUBROUTINE deco_run
!
!  Performs the actual decoration
!
   USE conn_mod
   USE modify_func_mod
!
   IMPLICIT none
!
   INTEGER, PARAMETER         :: MAXW = 2000  ! not ideal, should be dynamic ....
!
   CHARACTER(LEN=1024) :: mole_name ! molecule file name
   INTEGER   :: istatus          ! status
   INTEGER   :: i  ! dummy index
   INTEGER   :: ia ! dummy index
   INTEGER   :: mole_length ! length of molecule file name
   INTEGER   :: istart,iend ! dummy index
   INTEGER   :: is ! scattering number of surface atom
   INTEGER   :: idef ! connectivity definition number
   INTEGER, DIMENSION(:), ALLOCATABLE :: c_list
   INTEGER, DIMENSION(:,:), ALLOCATABLE :: c_offs ! Result of connectivity search
   INTEGER   :: natoms                 ! number of atoms in connectivity
   REAL   , DIMENSION(1:3) :: xyz
!
character(len=4) :: atom_name
integer ::itype
real, dimension(3) :: posit
real :: dw1
integer iprop
character(len=1) :: dummy
!
!  Load the molecules into temporary structures to reduce disk I/O
!
write(*,*) ' LOADING THE MOLECULES'
   CALL deco_get_molecules
write(*,*) ' Content of molecules'
DO i=1, dc_n_molecules
   DO ia=1, dc_molecules(i)%get_natoms()
   CALL dc_molecules(i)%get_cryst_atom(ia, itype, posit, iprop)
   CALL dc_molecules(i)%get_cryst_scat(ia, itype, atom_name, dw1)
   write(*,*) atom_name, itype, posit,iprop, dw1
   ENDDO
ENDDO
write(*,*) 'Weiter mit a return'
read(*,'(a)') dummy

!
   istart = 1
   iend   = cr_natoms
write(*,*) ' STARTING MAIN LOOP', istart, iend
   main_loop: DO ia=istart,iend
     is_sel: IF(check_select_status (dc_latom (cr_iscat (ia) ), cr_prop (ia),   &
                                     dc_sel_prop)              ) THEN
write(*,*)
write(*,*) ' Central : ',ia, cr_iscat(ia),cr_pos(:,ia), cr_prop(ia)
write(*,*) 'Weiter mit  return'
read(*,'(a)') dummy
        is     = cr_iscat(ia)                  ! get scattering type central
        xyz(:) = cr_pos(:,ia)                  ! Get atom position
        idef   = dc_use_conn(is)               ! use this connectivity list for surface
        CALL get_connectivity_list (ia, is, idef, maxw, c_list, c_offs, natoms )
write(*,*) 'CONNECT length : ', natoms, ' > ',c_list(1:natoms)
!       Go through all definitions
        dc_def_temp => dc_def_head        
        DO WHILE(ASSOCIATED(dc_def_temp))
write(*,*) ' Found a definition '
           dc_con_temp => dc_def_temp%dc_def_con
           DO WHILE(ASSOCIATED(dc_con_temp))           ! A connectivity exists
              CALL dc_get_con(dc_con_temp, dc_temp_surf, dc_temp_neig, dc_temp_dist)
write(*,*) ' Zentral, neigh, dist ', dc_temp_surf(1), dc_temp_neig, dc_temp_dist
              IF(cr_iscat(ia) == dc_temp_surf(1)) THEN    ! Matching atom in current definition
                 CALL dc_get_type(dc_def_temp, dc_temp_type)
                 CALL dc_get_axis(dc_def_temp, dc_temp_axis)
                 CALL dc_get_mole_name(dc_def_temp, mole_name, mole_length)
write(*,*) ' MATCHING PAIR TYPE', dc_temp_type
                 SELECT CASE(dc_temp_type)
                    CASE ( DC_NORMAL )                 ! Molecule in normal position
                       CALL deco_place_normal(dc_def_temp, ia, is, xyz, &
                            dc_temp_axis, mole_name, mole_length,       &
                            dc_temp_surf, dc_temp_neig, dc_temp_dist)
                    CASE ( DC_BRIDGE )                 ! Molecule in normal position
                       CALL deco_place_bridge(dc_def_temp, ia, is, xyz, &
                            dc_temp_axis, mole_name, mole_length,       &
                            dc_temp_surf, dc_temp_neig, dc_temp_dist)
                 END SELECT
              ENDIF
              dc_con_temp => dc_def_temp%dc_def_con%next
           ENDDO
           dc_def_temp => dc_def_temp%next
        ENDDO
     ENDIF is_sel ! END IF BLOCK is selected
   ENDDO main_loop   ! END DO main loop over all atoms
write(*,*) 'CLEANING UP'
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
   END SUBROUTINE deco_run
!
!######################################################################
   SUBROUTINE deco_get_molecules
!
   USE structur, ONLY: test_file
   IMPLICIT none
!
   CHARACTER (LEN=1024)                :: strufile
!
   INTEGER  :: i                ! dummy index
   INTEGER  :: istatus          ! status
   INTEGER  :: natoms           ! number of atoms in file
   INTEGER  :: ntypes           ! number of atom types in file
   INTEGER  :: n_mole           ! number of molecules
   INTEGER  :: n_type           ! number of molecule types
   INTEGER  :: n_atom           ! number of atoms in molecules
   INTEGER  :: init   = -1      ! Initialize atom names/types
   LOGICAL  :: lcell  = .true.  ! Treat atoms with equal name and B as one type
!
   ALLOCATE(m_name  (dc_n_molecules), STAT=istatus)
   ALLOCATE(m_lname (dc_n_molecules), STAT=istatus)
   ALLOCATE(m_ntypes(dc_n_molecules), STAT=istatus)
   ALLOCATE(m_ntypes(dc_n_molecules), STAT=istatus)
   ALLOCATE(m_length(dc_n_molecules), STAT=istatus)
!
   ALLOCATE(dc_molecules(dc_n_molecules), STAT = istatus)
!
   DO i=1,dc_n_molecules        ! load all molecules
      strufile = dc_input(i)
      CALL test_file(strufile, natoms, ntypes, n_mole, n_type, &
                     n_atom, init, lcell)
      CALL dc_molecules(i)%alloc_arrays(natoms, ntypes, n_mole, n_atom)
      CALL read_crystal ( dc_molecules(i), strufile )
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
!
   dc_init        = .true.    ! We need to initialize
   dc_n_molecules = 0         ! There are no molecules
   CALL dc_reset_def ( dc_def_head)
!
   END SUBROUTINE deco_reset
!
!******************************************************************************
   SUBROUTINE read_crystal ( this, infile)
!
!  Read a crystal structure from file
!  This procedure interfaces to the old "reastru" in "structur.f90"
!
   USE inter_readstru
   USE structur, ONLY: readstru
!
   IMPLICIT none
!
   TYPE (cl_cryst)                  :: this   ! The current crystal
   CHARACTER (LEN=*   ), INTENT(IN)  :: infile ! Disk file
   INTEGER                           :: inum   ! dummy index
   REAL   , DIMENSION(3)             :: posit  ! dummy position vector
   INTEGER                           :: istat  ! status variable
!
   rd_strucfile = infile
   rd_NMAX      = this%get_natoms() ! cr_natoms
   rd_MAXSCAT   = this%get_nscat()  ! cr_nscat
!
!  Allocate temporary arrays
!
   ALLOCATE ( rd_cr_dw    (0:rd_MAXSCAT)    , STAT = istat )
   ALLOCATE ( rd_cr_at_lis(0:rd_MAXSCAT)    , STAT = istat )
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
   rd_cr_dw   = 0.0
   rd_cr_at_lis = ' '
   rd_cr_pos    = 0.0
   rd_cr_iscat  = 0
   rd_cr_prop   = 0
   rd_as_at_lis = ' '
   rd_as_dw     = 0.0
   rd_as_pos    = 0.0
   rd_as_iscat  = 0
   rd_as_prop   = 0
!
   CALL readstru (rd_NMAX, rd_MAXSCAT, rd_strucfile, rd_cr_name,        &
               rd_cr_spcgr, rd_cr_a0, rd_cr_win, rd_cr_natoms, rd_cr_nscat, rd_cr_dw,     &
               rd_cr_at_lis, rd_cr_pos, rd_cr_mole, rd_cr_iscat, rd_cr_prop, rd_cr_dim, rd_as_natoms, &
               rd_as_at_lis, rd_as_dw, rd_as_pos, rd_as_iscat, rd_as_prop, rd_sav_ncell,  &
               rd_sav_r_ncell, rd_sav_ncatoms, rd_spcgr_ianz, rd_spcgr_para)
!
!
   DO inum=1, rd_NMAX
     posit = rd_cr_pos(:,inum)
!     CALL this%atoms(inum)%set_atom ( rd_cr_iscat(inum), posit, rd_cr_prop(inum) )
     CALL this%set_cryst_atom ( inum, rd_cr_iscat(inum), posit , rd_cr_prop(inum) )
   ENDDO
   CALL this%set_cryst_at_lis( rd_MAXSCAT, rd_cr_nscat, rd_cr_at_lis)
   CALL this%set_cryst_dw    ( rd_MAXSCAT, rd_cr_nscat, rd_cr_dw)
!
   DEALLOCATE ( rd_cr_dw    , STAT = istat )
   DEALLOCATE ( rd_cr_at_lis, STAT = istat )
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
   SUBROUTINE deco_place_normal(dc_def_temp, ia, ia_scat, xyz, &
                            mole_axis, mole_name, mole_length, &
                            surf, neig, dist)
!
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
!
   CHARACTER (LEN=4) :: atom_name
   CHARACTER (LEN=1024) :: line
   INTEGER :: i, im, laenge
   INTEGER :: itype
   INTEGER   :: surf_char ! connectivity definition number
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
   CALL find_surface_character(ia,ia, surf_char, surf_normal)
   IF(surf_char == SURF_PLANE ) THEN                      ! Ignore other than planar surfaces
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
!           CALL symm_show                                ! Show only in debug
            CALL symm_op_single                           ! Perform the operation
            ENDIF
         ENDIF
      ENDDO moles
   ENDIF
   1000 FORMAT(a4,4(2x,',',F12.6))
   1100 FORMAT(6(F12.6,', '),'ddd')
   END SUBROUTINE deco_place_normal
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
   CHARACTER (LEN=1024), DIMENSION(1:MAXW) :: cpara
   INTEGER             , DIMENSION(1:MAXW) :: lpara
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
!  Find the second partner involved in the bridge 
   x(1)     = cr_pos(1,ia)
   x(2)     = cr_pos(2,ia)
   x(3)     = cr_pos(3,ia)
   rmin     = 0.1
   radius   = dist*2.0
   ianz     = 1
   werte(1) = surf(2)
   fp (:)   = chem_period (:)
   fq       = chem_quick
   call do_find_env (ianz, werte, maxw, x, rmin, radius, fq, fp)  ! Find all neighbors
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
               CALL symm_show
               CALL symm_op_single
            ENDIF
         ENDIF
      ENDDO moles
   ENDIF
!
1000 FORMAT(a4,4(2x,',',F12.6))
1100 FORMAT(6(F12.6,', '),'ddd')
!
   END SUBROUTINE deco_place_bridge
END MODULE mole_surf_mod
