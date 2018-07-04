MODULE conn_mod
!
USE conn_def_mod
USE conn_sup_mod
USE conn_type_mod
USE crystal_mod
!
USE errlist_mod
USE random_mod
!
IMPLICIT none
!
PRIVATE
PUBLIC  conn_menu              !  Main menu to interact with user
!PUBLIC  get_connectivity_list  !  Read out the actual list of atoms around a central atom
!PUBLIC  get_connectivity_identity ! Identify a connectivity definition
!PUBLIC  get_connectivity_numbers  ! Read out number of connectivities for atom type
PUBLIC  do_show_connectivity   !  Show the current definitions
PRIVATE allocate_conn_list     !  Allocate memory for the connectivity
PRIVATE deallocate_conn        !  Free memory
PUBLIC  create_connectivity    !  Create the actual list of neighbors around each atom
PUBLIC  conn_do_set            !  Set parameters for the connectivity definitions
PUBLIC  conn_show              !  Main show routine
PRIVATE conn_test              !  While developing, a routine to test functionality
PRIVATE bond_switch_para       !  Define  a bond witching operation
PRIVATE do_bond_switch         !  Perform a bond witching operation
PRIVATE get_connect_pointed    !  Read out connectivity list return the pointer to list
PRIVATE do_exchange            !  Helper to exchange atoms between connectivities
PUBLIC  conn_update            !  Update the connectivity for an atom
!
INTEGER, PARAMETER              :: MAX_ATOM=10
!
CONTAINS
!
!  Module procedures to allocate, deallocate,
!
   SUBROUTINE allocate_conn_list(MAX_ATOM)
!
!  Simply allocates the connectivity list
!
   IMPLICIT  none
!
   INTEGER, INTENT(in)  :: MAX_ATOM
!
   INTEGER              :: i
!
   IF(.not. ALLOCATED(at_conn)) THEN
      ALLOCATE (at_conn(1:MAX_ATOM))                 ! allocate the array
      DO i=1,MAX_ATOM                                ! initialise the array
        at_conn(i)%number       = 0
        NULLIFY (at_conn(i)%liste)                   ! Initialy no NEIGHBORHOOD
      END DO
   ENDIF
   END SUBROUTINE allocate_conn_list
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   SUBROUTINE deallocate_conn(MAX_ATOM)
!
!  Properly deallocates the connectivity list
!
   IMPLICIT  none
!
   INTEGER, INTENT(in)  :: MAX_ATOM
!
   INTEGER              :: i
!
!
! Do proper removal of connectivity list
!
   IF(ALLOCATED(at_conn)) THEN
      cln_nghb : DO i=1,MAX_ATOM                         ! Loop over all atoms
        liste: IF ( ASSOCIATED(at_conn(i)%liste) ) THEN  ! If neighborhood was created
          hood_temp => at_conn(i)%liste                  ! point to the first neighborhood
          hood_head => at_conn(i)%liste
          hood: DO WHILE ( ASSOCIATED(hood_temp) )       ! While there are further neighborhood
             temp => hood_temp%nachbar                   ! point to the first neighbor within current neigborhood
             head => hood_temp%nachbar
             DO WHILE ( ASSOCIATED(temp) )               ! While there are further neighbors
               temp => temp%next                         ! Point to next neighbor
!              IF(ALLOCATED(head)) THEN
               IF(ASSOCIATED(head)) THEN
                  DEALLOCATE ( head )                       ! deallocate previous neighbor
               ENDIF
               head => temp                              ! also point to next neighbor
             END DO
             hood_temp => hood_temp%next_neighborhood    ! Point to next neighborhood
!            IF(ALLOCATED(hood_head)) THEN
             IF(ASSOCIATED(hood_head)) THEN
                DEALLOCATE ( hood_head )                 ! deallocate previous neighborhood
             ENDIF
             hood_head => hood_temp
          END DO hood                                    ! Loop over neighborhoods
        END IF liste                                     ! If atom has neigborhoods
      END DO cln_nghb                                    ! Loop over atoms
!
!
      DEALLOCATE (at_conn )                     ! Deallocate the initial array
!
      conn_status = .false.
   ENDIF
!
   END SUBROUTINE deallocate_conn
!
   SUBROUTINE create_connectivity
!
!  Performs a loop over all atoms and creates the individual connectivities
!
   USE chem_mod
   USE crystal_mod
   USE atom_env_mod
!  USE modify_mod
   USE do_find_mod
!
   IMPLICIT NONE
!
!  INTEGER, INTENT(IN)  :: i
!
   INTEGER, PARAMETER  :: MIN_PARA = 1
   INTEGER             :: maxw
!
   REAL   , DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: werte ! Array for neighbors
!
   INTEGER              :: j,i
   INTEGER              :: is  ! dummies for scattering types
   INTEGER              :: ianz
   INTEGER              :: n_neig      ! Actual number of neighboring atoms
   INTEGER,DIMENSION(:), ALLOCATABLE :: valid_neig  ! Valid neigbors under the scope
   LOGICAL, DIMENSION(3):: fp    ! periodic boundary conditions
   LOGICAL              :: fq    ! quick search algorithm
   REAL                 :: rmin        ! Minimum bond length
   REAL                 :: rmax        ! Maximum bond length
   REAL   , DIMENSION(3)     :: x      ! Atom position
!
   maxw = MAX(MIN_PARA, MAXSCAT+1)
!
   CALL deallocate_conn(conn_nmax)                     ! Deallocate old connectivity
   conn_nmax = cr_natoms                               ! Remember current atom number
   CALL allocate_conn_list(conn_nmax)                  ! Allocate connectivity
!
   fp (1) = chem_period (1)
   fp (2) = chem_period (2)
   fp (3) = chem_period (3)
   fq     = chem_quick
!
   atome: DO i = 1,cr_natoms                                ! Check all atoms in the structure    
      x(1) = cr_pos(1,i)
      x(2) = cr_pos(2,i)
      x(3) = cr_pos(3,i)
      ianz = 1
      is   = cr_iscat(i)                                    ! Keep atom type
      allowed: IF ( ASSOCIATED(def_main(is)%def_liste )) THEN  ! def.s exist
         def_temp => def_main(is)%def_liste
         neighs: DO
            werte(1:def_temp%valid_no) =      &
               def_temp%valid_types(1:def_temp%valid_no)    ! Copy valid atom types
            ianz  = def_temp%valid_no                       ! Copy no. of valid atom types
            rmin     = def_temp%def_rmin                    ! Copy distance limits
            rmax     = def_temp%def_rmax
            CALL do_find_env (ianz, werte, maxw, x, rmin,rmax, fq, fp)
            at_conn(i)%number = i                           ! Just set atom no
!
            IF ( atom_env(0) > 0) THEN                      ! The atom has neighbors
!
!              FIND valid neighbors
               ALLOCATE(valid_neig(1:atom_env(0)))
               n_neig = 0
               list: DO j = 1, atom_env(0)                  ! Add all (intended) neighbors to list
                  ! Add current atom if scope matches the intended conditions
                  IF( def_temp%mole_scope==MOLE_SCOPE_IGN                   .OR. &
                     (def_temp%mole_scope==MOLE_SCOPE_WITH  .AND.                &
                       cr_mole(i)/=0 .AND. cr_mole(i)==cr_mole(atom_env(j))).OR. &
                     (def_temp%mole_scope==MOLE_SCOPE_OUT   .AND.                &
                       cr_mole(i)/=0 .AND. cr_mole(i)/=cr_mole(atom_env(j)))     &
                    ) THEN
                    IF(i==atom_env(j)) THEN
                       IF(NINT(atom_pos(1,j)-cr_pos(1,atom_env(j)))==0 .AND.  &
                          NINT(atom_pos(2,j)-cr_pos(2,atom_env(j)))==0 .AND.  &
                          NINT(atom_pos(3,j)-cr_pos(3,atom_env(j)))==0       ) THEN
                          CYCLE list         ! got identical atom without offset
                       ENDIF
                    ENDIF
                    n_neig = n_neig + 1
                    valid_neig(n_neig) = atom_env(j)
                  ENDIF
               ENDDO list
!
!              Properly set the pointer hood_temp
!
               IF ( ASSOCIATED(at_conn(i)%liste)) THEN      ! A previous NEIGHBORHOOD exists
                  ALLOCATE(hood_temp%next_neighborhood)     ! Create one NEIGHBORHOOD
                  hood_temp => hood_temp%next_neighborhood  ! Point to the new NEIGHBORHOOD
               ELSE
                  ALLOCATE (at_conn(i)%liste)               ! Create one NEIGHBORHOOD
                  hood_temp => at_conn(i)%liste             ! Point to current NEIGHBORHOOD
               ENDIF
!
!           Now we set parameters of current NEIGHBORHOOD
!
!              IF(def_temp%intend_no == -1) THEN
!                 n_neig = atom_env(0)
!              ELSE
!                 n_neig = MIN(atom_env(0),def_temp%intend_no)
!              ENDIF
               IF(def_temp%intend_no /= -1) THEN
                  n_neig = MIN(n_neig     ,def_temp%intend_no)
               ENDIF
               NULLIFY (hood_temp%next_neighborhood)        ! No further NEIGHBORHOODs
               hood_temp%central_number = i                 ! Just set central atom no.
               hood_temp%central_type   = cr_iscat(i)
               hood_temp%neigh_type     = def_temp%valid_id   ! Set definition type number
               hood_temp%conn_name      = def_temp%def_name   ! Set name from definition type
               hood_temp%conn_name_l    = def_temp%def_name_l ! Set name length from definition type
               hood_temp%mmc_sel        = def_temp%mmc_sel    ! Copy mmc selection mode
               hood_temp%mmc_ene        = def_temp%mmc_ene    ! Copy mmc energy    mode
               hood_temp%natoms         = n_neig              ! Set number of neighbors
               NULLIFY (hood_temp%nachbar)                  ! Initially there are no NEIGHBORS
!
               ALLOCATE (hood_temp%nachbar)                 ! create the first NEIGHBOR slot
               j = 1
               tail => hood_temp%nachbar                    ! tail points to the first NEIGHBOR
!              tail%atom_number = atom_env(j)               ! I store the atom_no of the neighbor
               tail%atom_number = valid_neig(j)             ! I store the atom_no of the neighbor
               tail%offset(1)   = NINT(atom_pos(1,j)-cr_pos(1,atom_env(j)))
               tail%offset(2)   = NINT(atom_pos(2,j)-cr_pos(2,atom_env(j)))
               tail%offset(3)   = NINT(atom_pos(3,j)-cr_pos(3,atom_env(j)))
               NULLIFY (tail%next)                          ! No further neighbors
!
               DO j = 2, n_neig                             ! Add all (intended) neighbors to list
                     ALLOCATE (tail%next)                      ! create a further NEIGHBOR
                     tail => tail%next                         ! reassign tail to new end of list
                     tail%atom_number = valid_neig(j)          ! I store the atom_no of the neighbor
                     tail%offset(1)   = NINT(atom_pos(1,j)-cr_pos(1,atom_env(j)))
                     tail%offset(2)   = NINT(atom_pos(2,j)-cr_pos(2,atom_env(j)))
                     tail%offset(3)   = NINT(atom_pos(3,j)-cr_pos(3,atom_env(j)))
                     NULLIFY (tail%next)                       ! No further neighbors
               ENDDO
               DEALLOCATE(valid_neig)
            ENDIF
            IF ( .NOT. ASSOCIATED(def_temp%def_next)) THEN  ! No more def.s
               CYCLE atome
            ENDIF
            def_temp => def_temp%def_next
         ENDDO neighs                                       ! Loop over def.s
      ENDIF allowed                                         ! Atom has def.s
   ENDDO atome                                              ! Loop over all atoms in structure
!
   conn_status = .true.
!
   END SUBROUTINE create_connectivity
!
   SUBROUTINE recreate_connectivity(itype, ino, c_name)
!
!  Performs a loop over all atoms and creates the individual connectivities
!
   USE chem_mod
   USE crystal_mod
   USE do_find_mod
   USE atom_env_mod
!  USE modify_mod
!
   IMPLICIT NONE
!
   INTEGER           , INTENT(IN)     :: itype   ! Atom type
   INTEGER           , INTENT(INOUT)  :: ino     ! Connectivity def. no.
   CHARACTER(LEN=256), INTENT(INOUT)  :: c_name  ! Connectivity name
!
!  INTEGER, INTENT(IN)  :: i
!
   INTEGER, PARAMETER  :: MIN_PARA = 1
   INTEGER             :: maxw
!
   REAL   , DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: werte ! Array for neighbors
!
   INTEGER              :: j,i
   INTEGER              :: is  ! dummies for scattering types
   INTEGER              :: ianz
   INTEGER              :: n_neig      ! Actual number of neighboring atoms
   LOGICAL, DIMENSION(3):: fp    ! periodic boundary conditions
   LOGICAL              :: fq    ! quick search algorithm
   LOGICAL              :: found_def   ! Found correct neiborhood def to replace
   LOGICAL              :: found       ! Found correct neiborhood to replace
   REAL                 :: rmin        ! Minimum bond length
   REAL                 :: rmax        ! Maximum bond length
   REAL   , DIMENSION(3)     :: x      ! Atom position
!
   maxw = MAX(MIN_PARA, MAXSCAT+1)
   found = .FALSE.
!
!  CALL deallocate_conn(conn_nmax)                     ! Deallocate old connectivity
   conn_nmax = cr_natoms                               ! Remember current atom number
   IF ( .NOT.ALLOCATED(at_conn)) THEN                 ! No previous NEIGHBORHOOD exists
      CALL allocate_conn_list(conn_nmax)               ! Allocate connectivity
   ENDIF
!
   fp (1) = chem_period (1)
   fp (2) = chem_period (2)
   fp (3) = chem_period (3)
   fq     = chem_quick
!
   atome: DO i = 1,cr_natoms                                ! Check all atoms in the structure    
      is   = cr_iscat(i)                                    ! Keep atom type
      IF(is/=itype) CYCLE atome
      x(1) = cr_pos(1,i)
      x(2) = cr_pos(2,i)
      x(3) = cr_pos(3,i)
      ianz = 1
      allowed: IF ( ASSOCIATED(def_main(is)%def_liste )) THEN  ! def.s exist
         def_temp => def_main(is)%def_liste
         found_def = .FALSE.
         neighs: DO
            IF(.NOT.(ino==def_temp%valid_id .OR. c_name==def_temp%def_name)) THEN
               def_temp => def_temp%def_next                ! not the right definition,
               CYCLE neighs                                 ! go to next definition
            ENDIF
            found_def = .TRUE.
            werte(1:def_temp%valid_no) =      &
               def_temp%valid_types(1:def_temp%valid_no)    ! Copy valid atom types
            ianz  = def_temp%valid_no                       ! Copy no. of valid atom types
            rmin     = def_temp%def_rmin                    ! Copy distance limits
            rmax     = def_temp%def_rmax
            CALL do_find_env (ianz, werte, maxw, x, rmin,rmax, fq, fp)
            at_conn(i)%number = i                           ! Just set atom no
!
            IF ( atom_env(0) > 0) THEN                      ! The atom has neighbors
!
!              Properly set the pointer hood_temp
!
               found = .FALSE.
               IF ( ASSOCIATED(at_conn(i)%liste)) THEN      ! A previous NEIGHBORHOOD exists
                  hood_temp => at_conn(i)%liste             ! Point to current NEIGHBORHOOD
find_hood:        DO WHILE(ASSOCIATED(hood_temp))
                     IF(ino==hood_temp%neigh_type .OR. &
                        c_name==hood_temp%conn_name   ) THEN ! Found the correct neighborhood
                        found = .TRUE.
                        EXIT find_hood
                     ENDIF
                     hood_temp => hood_temp%next_neighborhood  ! Point to the next NEIGHBORHOOD
                  ENDDO find_hood
                  IF(found) THEN                            ! We are at correct neighborhood
                     CONTINUE                               ! Deallocate old neighbors
                     temp => hood_temp%nachbar              ! point to the first neighbor within current neigborhood
                     head => hood_temp%nachbar
                     DO WHILE ( ASSOCIATED(temp) )          ! While there are further neighbors
                       temp => temp%next                    ! Point to next neighbor
                       IF(ASSOCIATED(head)) THEN
                          DEALLOCATE ( head )               ! deallocate previous neighbor
                       ENDIF
                       head => temp                         ! also point to next neighbor
                     END DO
                  ELSE                                      ! No neighborood found thus:
                     ALLOCATE(hood_temp%next_neighborhood)     ! Create one NEIGHBORHOOD
                     hood_temp => hood_temp%next_neighborhood  ! Point to the new NEIGHBORHOOD
                  ENDIF
               ELSE
                  ALLOCATE (at_conn(i)%liste)               ! Create one NEIGHBORHOOD
                  hood_temp => at_conn(i)%liste             ! Point to current NEIGHBORHOOD
                  NULLIFY (hood_temp%next_neighborhood)     ! No further NEIGHBORHOODs
               ENDIF
!
!           Now we set parameters of current NEIGHBORHOOD
!
               IF(def_temp%intend_no == -1) THEN
                  n_neig = atom_env(0)
               ELSE
                  n_neig = MIN(atom_env(0),def_temp%intend_no)
               ENDIF
               hood_temp%central_number = i                 ! Just set central atom no.
               hood_temp%central_type   = cr_iscat(i)
               hood_temp%neigh_type     = def_temp%valid_id   ! Set definition type number
               hood_temp%conn_name      = def_temp%def_name   ! Set name from definition type
               hood_temp%conn_name_l    = def_temp%def_name_l ! Set name length from definition type
               hood_temp%mmc_sel        = def_temp%mmc_sel    ! Copy mmc selection mode
               hood_temp%mmc_ene        = def_temp%mmc_ene    ! Copy mmc energy    mode
               hood_temp%natoms         = n_neig              ! Set number of neighbors
               NULLIFY (hood_temp%nachbar)                  ! Initially there are no NEIGHBORS
!
               ALLOCATE (hood_temp%nachbar)                 ! create the first NEIGHBOR slot
               j = 1
               tail => hood_temp%nachbar                    ! tail points to the first NEIGHBOR
               tail%atom_number = atom_env(j)               ! I store the atom_no of the neighbor
               tail%offset(1)   = NINT(atom_pos(1,j)-cr_pos(1,atom_env(j)))
               tail%offset(2)   = NINT(atom_pos(2,j)-cr_pos(2,atom_env(j)))
               tail%offset(3)   = NINT(atom_pos(3,j)-cr_pos(3,atom_env(j)))
               NULLIFY (tail%next)                          ! No further neighbors
!
               DO j = 2, n_neig                             ! Add all (intended) neighbors to list
                  ALLOCATE (tail%next)                      ! create a further NEIGHBOR
                  tail => tail%next                         ! reassign tail to new end of list
                  tail%atom_number = atom_env(j)            ! I store the atom_no of the neighbor
                  tail%offset(1)   = NINT(atom_pos(1,j)-cr_pos(1,atom_env(j)))
                  tail%offset(2)   = NINT(atom_pos(2,j)-cr_pos(2,atom_env(j)))
                  tail%offset(3)   = NINT(atom_pos(3,j)-cr_pos(3,atom_env(j)))
                  NULLIFY (tail%next)                       ! No further neighbors
               ENDDO
            ENDIF
            IF(found) CYCLE atome
            IF ( .NOT. ASSOCIATED(def_temp%def_next)) THEN  ! No more def.s
               CYCLE atome
            ENDIF
            def_temp => def_temp%def_next
         ENDDO neighs                                       ! Loop over def.s
            IF(.NOT.found_def) THEN                             ! found no correct definition
               ier_num = -109
               ier_typ = ER_APPL
               ier_msg(1) = 'None of the connectivity definitions for the'
               ier_msg(2) = 'atom type to be renewed matches the number or name'
               RETURN
            ENDIF
      ENDIF allowed                                         ! Atom has def.s
   ENDDO atome                                              ! Loop over all atoms in structure
!
   conn_status = .true.
!
   END SUBROUTINE recreate_connectivity
!
!
   SUBROUTINE conn_do_set ( code, zeile, length)
!-                                                                      
!     Set the parameters for the connectivity
!+                                                                      
      USE discus_allocate_appl_mod 
      USE discus_config_mod 
      USE crystal_mod 
      USE get_iscat_mod
      USE modify_mod
      USE variable_test
      USE berechne_mod
      USE ber_params_mod
      USE errlist_mod
      USE get_params_mod
      USE take_param_mod
!
      IMPLICIT none
!
!
      INTEGER          , INTENT(IN   )  :: code 
      CHARACTER (LEN=*), INTENT(INOUT)  :: zeile
      INTEGER          , INTENT(INOUT)  :: length 
!
      INTEGER, PARAMETER  :: MIN_PARA = 5
      INTEGER             :: maxw
      INTEGER, PARAMETER  :: maxw2 = 2
      INTEGER, PARAMETER  :: NOPTIONAL = 2
!
      CHARACTER(LEN=1024), DIMENSION(MAX(MIN_PARA,MAXSCAT+5)) :: cpara
      INTEGER            , DIMENSION(MAX(MIN_PARA,MAXSCAT+5)) :: lpara
      REAL               , DIMENSION(MAX(MIN_PARA,MAXSCAT+5)) :: werte
!
      CHARACTER(LEN=1024), DIMENSION(2)  :: ccpara
      INTEGER            , DIMENSION(2)  :: llpara
      REAL               , DIMENSION(2)  :: wwerte
!
      INTEGER             :: ianz         ! number of command line parameters
      INTEGER             :: iianz        ! dummy number
      INTEGER             :: is1          ! first atom type
      INTEGER             :: is2          ! second atom type
      INTEGER             :: temp_id      ! temporary definition ID
      INTEGER             :: temp_number  ! temporary Number of neighbors
      INTEGER             :: work_id      ! ID of the definition to change/delete
      CHARACTER(LEN=256)  :: work_name    ! Name of the definition to change/delete
      INTEGER             :: work_name_l  ! Length of name for the definition to change/delete
      LOGICAL             :: lnew         ! require atom type to exist
      LOGICAL             :: l_exist      ! TRUE if conn. name is a variable name
      LOGICAL             :: l_type       ! Unused DUMMY argument
      INTEGER             :: is_no        ! Unused DUMMY argument
      INTEGER             :: all_status   ! Allocation status
      INTEGER             :: mole_scope   ! Scope is limited to a molecule
      REAL                :: rmin         ! minimum bond distance
      REAL                :: rmax         ! maximum bond distance
!
      CHARACTER(LEN=   9), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
      CHARACTER(LEN=1024), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
      INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
      INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
      REAL               , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
      INTEGER, PARAMETER                        :: ncalc = 1 ! Number of values to calculate 
!
      LOGICAL :: str_comp
!
      DATA oname  / 'first', 'molescope' /
      DATA loname /  5     ,  9          /
      opara  =  (/ '-1     ', 'ignore ' /)   ! Always provide fresh default values
      lopara =  (/  7       ,  7        /)
      owerte =  (/  -1.0    ,  0.0      /)
!
!                                                                       
      rmin = 0.0
      rmax = 0.5
      maxw = MAX(MIN_PARA, MAXSCAT+5)     ! (MAXSCAT+ void) + 4 Parameters  
      temp_number = NINT(owerte(1))       ! Take default from optional parameter list
!                                                                       
!     Check definitions array
!
      IF ( .NOT. ALLOCATED(def_main)) THEN
         ALLOCATE (def_main(0:MAXSCAT), stat = all_status) 
         DO is1 = 0, MAXSCAT
            NULLIFY(def_main(is1)%def_liste)
            def_main(is1)%def_id     = is1
            def_main(is1)%def_number = 0
         ENDDO
!        conn_max_def = MAXSCAT
         IF ( all_status /= 0 ) THEN
            ier_num = -3
            ier_typ = ER_COMM
            ier_msg(1) = ' Error allocating the definitions'
            WRITE(ier_msg(2),'(a,i6)') 'Error no:',all_status
            RETURN
         ENDIF
      ENDIF
!
!     Get parameters from input line
!
      CALL get_params (zeile, ianz, cpara, lpara, maxw, length) 
      IF (ier_num.ne.0) return 
!
!     Sort optional parameters
!
      opara  =  (/ '-1     ', 'ignore ' /)   ! Always provide fresh default values
      lopara =  (/  7       ,  7        /)
      owerte =  (/  -1.0    ,  0.0      /)
      CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                        oname, loname, opara, lopara, owerte)
      IF(     str_comp(opara(2), 'ignore',  2, lopara(2), 6) ) THEN
         mole_scope = MOLE_SCOPE_IGN
      ELSEIF( str_comp(opara(2), 'within',  2, lopara(2), 6) ) THEN
         mole_scope = MOLE_SCOPE_WITH
      ELSEIF( str_comp(opara(2), 'outside', 2, lopara(2), 7) ) THEN
         mole_scope = MOLE_SCOPE_OUT
      ELSE
         ier_num = -13
         ier_typ = ER_COMM
         ier_msg(1) = 'Offending parameter'
         ier_msg(2) = oname(2)(1:MIN(43,loname(2)))
         ier_msg(3) = opara(2)(1:MIN(43,lopara(2)))
         RETURN
      ENDIF

!
!     Set mmc behaviour
!
!     IF( str_comp (cpara(1), 'mmc',    2, lpara(1), 3) ) THEN
!        IF( str_comp (cpara(2), 'select', 2, lpara(2), 6) ) THEN
!           IF( str_comp (cpara(3), 'on', 2, lpara(3), 2) ) THEN
!              conn_mmc_sel =  STATUS_ON
!           ELSEIF( str_comp (cpara(3), 'off', 2, lpara(3), 3) ) THEN
!              conn_mmc_sel =  STATUS_OFF
!           ELSEIF( str_comp (cpara(3), 'ignore', 2, lpara(3), 6) ) THEN
!              conn_mmc_sel =  STATUS_IGN
!           ELSE
!              ier_num = -6
!              ier_typ = ER_FORT
!           ENDIF
!           RETURN
!        ELSEIF( str_comp (cpara(2), 'energy', 2, lpara(2), 6) ) THEN
!           IF( str_comp (cpara(3), 'on', 2, lpara(3), 2) ) THEN
!              conn_mmc_ene = STATUS_ON
!           ELSEIF( str_comp (cpara(3), 'off', 2, lpara(3), 3) ) THEN
!              conn_mmc_ene = STATUS_OFF
!           ELSEIF( str_comp (cpara(3), 'ignore', 2, lpara(3), 6) ) THEN
!              conn_mmc_ene = STATUS_IGN
!           ELSE
!              ier_num = -6
!              ier_typ = ER_FORT
!           ENDIF
!           RETURN
!        ENDIF
!     ENDIF
!
!     Remove all definitions 
!
      reset: IF ( code == code_res ) THEN        ! remove all definitions
         exist_def: IF ( ALLOCATED(def_main)) THEN    ! Are there any definitions
            DO is1 = 0, MAXSCAT
               is_range: IF (is1 <= UBOUND(def_main,1) ) THEN 
               is_reset: IF ( ASSOCIATED(def_main(is1)%def_liste) ) THEN  ! A list of definitions exists
                  def_head => def_main(is1)%def_liste
                  def_temp => def_main(is1)%def_liste
                  work_id = 1
                  DO WHILE ( ASSOCIATED(def_temp%def_next) )        ! A further definition exists
                     IF ( work_id == 1 ) THEN                       ! This is the first def.
                        def_main(is1)%def_liste => def_temp%def_next  ! properly point in def_main
                        DEALLOCATE(def_temp, stat=all_status)       ! Remove this node
                        def_temp => def_main(is1)%def_liste         ! Point to next node
                     ELSE                                           ! Not the first definition
                        def_head%def_next   => def_temp%def_next    ! Point previous next to further
                        DEALLOCATE(def_temp, stat=all_status)       ! Remove this node
                        def_temp => def_head%def_next               ! Point to next node
                        work_id = work_id + 1                       ! Increment node number
                     ENDIF
                  ENDDO
                  DEALLOCATE(def_main(is1)%def_liste, stat=all_status) !Remove the whole list for this scat type
               ENDIF is_reset
               ENDIF is_range
            ENDDO
            DEALLOCATE(def_main, stat=all_status)                   ! Deallocate the main structure
         ENDIF exist_def
         RETURN                           ! All definitions have been removed
      ENDIF reset
!
!     All other codes
!
      IF ((code==code_del .AND. ianz/=2) .OR. (code/=code_del .AND.ianz < 4)) THEN
         ier_num = -6         ! Wrong number of input parameters
         ier_typ = ER_COMM
         return      ! At least four parameters
      ENDIF
!                                                                       
      iianz = 1
      lnew  = .false.
!
!     Get the first atom type
!
      CALL get_iscat (iianz, cpara, lpara, werte, maxw, lnew)
      is1 = NINT(werte(1))
      CALL del_params (1, ianz, cpara, lpara, maxw) 
!
      IF ( code /= code_add ) then
         iianz = 1
         CALL ber_params (iianz, cpara, lpara, werte, maxw)
         IF(ier_num == 0) THEN
            work_id     = NINT(werte(1))
            work_name   = ' '
            work_name_l = 1
         ELSE
            work_id     = -2
            work_name   = cpara(iianz)(1:lpara(iianz))
            work_name_l = lpara(iianz)
            call no_error
         ENDIF
         CALL variable_exist (work_name, work_name_l,0, l_exist, l_type, is_no)
         IF(l_exist) THEN
            ier_num = -120
            ier_typ = ER_APPL
            RETURN
         ENDIF
         CALL discus_validate_var_spec(work_name, work_name_l)
         IF( ier_num == -25 ) THEN
            ier_num = -120
            ier_typ = ER_APPL
            RETURN
         ENDIF
         CALL del_params (1, ianz, cpara, lpara, maxw) 
      ELSE
         work_id     = -1
         work_name   = cpara(ianz)(1:lpara(ianz))
         work_name_l = lpara(ianz)
         ianz        = ianz - 1
         IF(cpara(ianz)(1:6)=='first_') THEN
            cpara(ianz) = cpara(ianz)(7:lpara(ianz))
            temp_number = NINT(berechne(cpara(ianz), lpara(ianz)))
            ianz        = ianz - 1
         ELSE
            temp_number = NINT(owerte(1))  ! Take default from optional parameter list
         ENDIF
      ENDIF
!
      IF ( code /= code_del ) THEN
!
!     Get minimum and maximum bond length, last two parameters
!
         ccpara(1) = cpara(ianz-1)
         llpara(1) = lpara(ianz-1)
         ccpara(2) = cpara(ianz  )
         llpara(2) = lpara(ianz  )
         iianz     = 2
         CALL ber_params (iianz, ccpara, llpara, wwerte, maxw2)
         rmin = wwerte(1)
         rmax = wwerte(2)
         ianz = ianz - 2
!
!     get scattering types of neighbors
!
         CALL get_iscat (ianz, cpara, lpara, werte, maxw, lnew)
      ENDIF
!
      is_there: IF ( ASSOCIATED(def_main(is1)%def_liste) ) THEN  ! A list of definitions exists
         is_work: IF ( work_id /= -1 ) THEN                      ! Work on an existing definition
            IF ( work_id > def_main(is1)%def_number ) THEN       ! Definition does not exist
               ier_num = -109
               ier_typ = ER_APPL
               RETURN
            ENDIF
            def_head => def_main(is1)%def_liste
            def_temp => def_main(is1)%def_liste
            search: DO                                           ! search for working definition
               IF ( .NOT. ASSOCIATED(def_temp)) THEN             ! target is not associated ERROR
                  ier_num = -109
                  ier_typ = ER_APPL
                  RETURN
               ENDIF
               IF ( work_id   == def_temp%valid_id  .OR. &
                    work_name == def_temp%def_name       ) THEN  ! Found working definition
                  work_id   = def_temp%valid_id                  ! Make sure ID matches
                  work_name = def_temp%def_name                  ! Make sure name matches
                  EXIT search
               ENDIF
               def_head => def_temp
               def_temp => def_temp%def_next
            ENDDO search                                         ! End search for working definition
            IF ( code == code_set ) THEN                         ! Replace current entries
               DEALLOCATE( def_temp%valid_types, stat = all_status) ! delete old arry
               IF ( all_status /= 0 ) THEN
                  ier_num = -3
                  ier_typ = ER_COMM
                  ier_msg(1) = ' Error allocating the definitions'
                  WRITE(ier_msg(2),'(a,i6)') 'Error no:',all_status
                  RETURN
               ENDIF
               ALLOCATE ( def_temp%valid_types(1:ianz), stat = all_status) ! alloc array of neighbor types
               IF ( all_status /= 0 ) THEN
                  ier_num = -3
                  ier_typ = ER_COMM
                  ier_msg(1) = ' Error allocating the definitions'
                  WRITE(ier_msg(2),'(a,i6)') 'Error no:',all_status
                  RETURN
               ENDIF
               DO is2 = 1, ianz                                  ! Set all neighbor types
                  def_temp%valid_types(is2) = NINT(werte(is2))
               ENDDO
               def_temp%valid_no   = ianz                        ! Set number of neighb or types
               def_temp%intend_no  = temp_number                 ! Set intended number of neighbor atoms
               def_temp%mole_scope = mole_scope                  ! Set scope
               def_temp%def_rmin   = rmin                        ! Set bond length limits
               def_temp%def_rmax   = rmax                        ! Set bond length limits
            ELSEIF ( code == code_del ) THEN                     ! Remove this definition
               IF ( ASSOCIATED(def_temp%def_next) ) THEN         ! A further definition exists
                  IF ( work_id == 1 ) THEN                       ! This is the first def.
                     def_main(is1)%def_liste => def_temp%def_next  ! properly point in def_main
                     DEALLOCATE(def_temp)                        ! Remove this node
                     def_temp => def_main(is1)%def_liste         ! Point to next node
                  ELSE                                           ! Not the first definition
                     def_head%def_next   => def_temp%def_next    ! Point previous next to further
                     DEALLOCATE(def_temp)                        ! Remove this node
                     def_temp => def_head%def_next               ! Point to next node
                  ENDIF
                  def_temp%valid_id = def_temp%valid_id -1       ! Decrement the ID of current node
                  def_temp => def_temp%def_next                  ! Point to next node
                  search3: DO                                    ! Decrement ID's of all further nodes
                     IF ( ASSOCIATED(def_temp)) THEN
                        def_temp%valid_id = def_temp%valid_id -1
                     ELSE
                        EXIT search3
                     ENDIF
                     def_temp => def_temp%def_next
                  ENDDO search3
               ELSE                                              ! No further definition exists
                  IF ( work_id == 1 ) THEN                       ! This is the first def.
                     DEALLOCATE(def_main(is1)%def_liste)         ! Remove the whole list for this scat type
                  ELSE                                           ! Not the first definition
                     NULLIFY(def_head%def_next)                  ! Nullify previous pointer to next
                     DEALLOCATE(def_temp)                        ! Remove current definition
                  ENDIF
               ENDIF
            ENDIF
         ELSE is_work                                            ! Add a new definition
            def_head => def_main(is1)%def_liste                  ! Point to first def.
            def_temp => def_main(is1)%def_liste                  ! Point to first def.
            temp_id  = 0                                         ! In case we have no def.s at all
            search2: DO                                          ! search for working definition
               IF ( .NOT. ASSOCIATED(def_temp)) THEN             ! TRUE if at end of definitions
                  EXIT search2
               ENDIF
               def_head => def_temp                              ! Increment point to next def.
               def_temp => def_temp%def_next                     ! Increment point to next def.
               temp_id  = def_head%valid_id                      ! Store previous ID number
            ENDDO search2                                        ! End search for working definition
            ALLOCATE(def_temp, stat=all_status)                  ! Allocate a next node
            IF ( all_status /= 0 ) THEN
               ier_num = -3
               ier_typ = ER_COMM
               ier_msg(1) = ' Error allocating the definitions'
               WRITE(ier_msg(2),'(a,i6)') 'Error no:',all_status
               RETURN
            ENDIF
            ALLOCATE ( def_temp%valid_types(1:ianz), stat = all_status) ! alloc array of neighbor types
            IF ( all_status /= 0 ) THEN
               ier_num = -3
               ier_typ = ER_COMM
               ier_msg(1) = ' Error allocating the definitions'
               WRITE(ier_msg(2),'(a,i6)') 'Error no:',all_status
               RETURN
            ENDIF
            def_head%def_next => def_temp                     ! Point previous def. to current
            DO is2 = 1, ianz                                  ! Set all neighbor types
               def_temp%valid_types(is2) = NINT(werte(is2))
            ENDDO
            def_temp%valid_id   = temp_id + 1                 ! Set number of neighb or types
            def_temp%def_name   = work_name                   ! Set definition name
            def_temp%def_name_l = work_name_l                 ! Set definition name length
            def_temp%mmc_sel    = conn_mmc_sel                ! Set mmc selection type
            def_temp%mmc_ene    = conn_mmc_ene                ! Set mmc energy    type
            def_temp%valid_no   = ianz                        ! Set number of neighbor types
            def_temp%intend_no  = temp_number                 ! Set intended number of neighbor atoms
            def_temp%mole_scope = mole_scope                  ! Set scope
            def_temp%def_rmin   = rmin                        ! Set bond length limits
            def_temp%def_rmax   = rmax                        ! Set bond length limits
            NULLIFY(def_temp%def_next)                        ! No further definition
            def_main(is1)%def_number = def_main(is1)%def_number + 1
         ENDIF is_work
      ELSE  is_there                                          ! No list exists yet, add first node
         IF ( code == code_add ) THEN                         ! Replace current entries
            ALLOCATE(def_main(is1)%def_liste, stat = all_status) ! allocate a list of defs.
            def_main(is1)%def_id     = is1                    ! Set identifier
            def_main(is1)%def_number = 0                      ! No def.s exist yet.
            IF ( all_status /= 0 ) THEN
               ier_num = -3
               ier_typ = ER_COMM
               ier_msg(1) = ' Error allocating the definitions'
               WRITE(ier_msg(2),'(a,i6)') 'Error no:',all_status
               RETURN
            ENDIF
            def_temp => def_main(is1)%def_liste
            ALLOCATE ( def_temp%valid_types(1:ianz), stat = all_status) ! alloc array of neighbor types
            IF ( all_status /= 0 ) THEN
               ier_num = -3
               ier_typ = ER_COMM
               ier_msg(1) = ' Error allocating the definitions'
               WRITE(ier_msg(2),'(a,i6)') 'Error no:',all_status
               RETURN
            ENDIF
            DO is2 = 1, ianz                                  ! Set all neighbor types
               def_temp%valid_types(is2) = NINT(werte(is2))
            ENDDO
            def_temp%valid_id   = 1                           ! Set number of neighb or types
            def_temp%def_name   = work_name                   ! Set definition name
            def_temp%def_name_l = work_name_l                 ! Set definition name length
            def_temp%mmc_sel    = conn_mmc_sel                ! Set mmc selection type
            def_temp%mmc_ene    = conn_mmc_ene                ! Set mmc energy    type
            def_temp%valid_no   = ianz                        ! Set number of neighb or types
            def_temp%intend_no  = temp_number                 ! Set intended number of neighbor atoms
            def_temp%mole_scope = mole_scope                  ! Set scope
            def_temp%def_rmin   = rmin                        ! Set bond length limits
            def_temp%def_rmax   = rmax                        ! Set bond length limits
            NULLIFY(def_temp%def_next)                        ! No further definition
            def_main(is1)%def_number = def_main(is1)%def_number + 1
         ENDIF
      ENDIF is_there
!
!
   END SUBROUTINE conn_do_set
!
   SUBROUTINE conn_menu
!-                                                                      
!     Main menu for connectivity related operations                          
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE get_iscat_mod
      USE modify_mod
!
      USE ber_params_mod
      USE calc_expr_mod
      USE doact_mod 
      USE do_eval_mod
      USE do_wait_mod
      USE learn_mod 
      USE class_macro_internal
      USE get_params_mod
      USE prompt_mod 
      USE sup_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER, PARAMETER :: MIN_PARA = 4
      INTEGER maxw 
!                                                                       
      CHARACTER(5) befehl 
!     CHARACTER(50) prom 
      CHARACTER(LEN=LEN(PROMPT)) :: orig_prompt 
      CHARACTER(1024) line, zeile
      CHARACTER(LEN=1024), DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: cpara ! (MAXSCAT) 
      INTEGER            , DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) :: lpara ! (MAXSCAT)
      CHARACTER (LEN=256)  :: c_name   ! Connectivity name
      INTEGER              :: c_name_l ! connectivity name length
      INTEGER              :: ino      ! connectivity no
      INTEGER              :: iatom    ! atoms no for show
      INTEGER              :: itype    ! atomR type for recreate
      INTEGER lp, length, lbef 
      INTEGER indxg, ianz, iianz
      LOGICAL              :: long     ! make long output
      LOGICAL              :: lnew     ! Do not make new atom type
      LOGICAL lend
      REAL               , DIMENSION(MAX(MIN_PARA,MAXSCAT+1)) ::  werte ! (MAXSCAT) 
!                                                                       
      INTEGER len_str 
      LOGICAL str_comp 
!                                                                       
      maxw = MAX(MIN_PARA,MAXSCAT+1)
      lend = .false. 
      CALL no_error 
      orig_prompt = prompt
      prompt = prompt (1:len_str (prompt) ) //'/conn' 
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
            IF (indxg.ne.0.AND.                                      &
                .NOT. (str_comp (befehl, 'echo', 2, lbef, 4) ) .AND. &
                .NOT. (str_comp (befehl, 'syst', 2, lbef, 4) ) .AND. &
                .NOT. (str_comp (befehl, 'help', 2, lbef, 4) .OR.    &
                       str_comp (befehl, '?   ', 2, lbef, 4) ) .AND. &
                INDEX(line, '==') == 0                               &
               ) then                                              
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
!     ----continues a macro 'continue'                                  
!                                                                       
               ELSEIF (str_comp (befehl, 'continue', 2, lbef, 8) ) then 
                  CALL macro_continue (zeile, lp) 
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
!     ----help 'help','?'                                               
!                                                                       
               ELSEIF (str_comp (befehl, 'help', 2, lbef, 4) .or.  &
                       str_comp (befehl, '?   ', 1, lbef, 4) ) then                                      
                  IF (str_comp (zeile, 'errors', 2, lp, 6) ) then 
                     lp = lp + 7 
                     CALL do_hel ('discus '//zeile, lp) 
                  ELSE 
                     lp = lp + 12 
                     CALL do_hel ('discus conn '//zeile, lp) 
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
!------  -----waiting for user input                                    
!                                                                       
               ELSEIF (str_comp (befehl, 'wait', 3, lbef, 4) ) then 
                  CALL do_input (zeile, lp) 
!                                                                       
!     ----create a connectivity list 'create'                                     
!                                                                       
               ELSEIF (str_comp (befehl, 'create', 2, lbef, 6) ) then 
                  CALL create_connectivity
!                                                                       
!     ----delete a connectivity list 'delete'                                     
!                                                                       
               ELSEIF (str_comp (befehl, 'delete', 2, lbef, 6) ) then 
                  CALL deallocate_conn (conn_nmax)
!                                                                       
!     ----add a new connectivity definition       'add'
!                                                                       
               ELSEIF (str_comp (befehl, 'add', 2, lbef, 3) ) then 
                  CALL conn_do_set (code_add,zeile, lp) 
!                                                                       
!     ----remove an old    connectivity definition 'remove'                 
!                                                                       
               ELSEIF (str_comp (befehl, 'remove', 3, lbef, 6) ) then 
                  CALL conn_do_set (code_del,zeile, lp) 
!                                                                       
!     ----recreate a connectivity list 'recreate'                                     
!                                                                       
               ELSEIF (str_comp (befehl, 'recreate', 3, lbef, 8) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, length) 
                  IF(ier_num==0) THEN
                     itype  = 0
                     ino    = 0
                     c_name = ' '
                     IF(ianz==2) THEN
                        iianz = 1
                        lnew  = .false.
                        CALL get_iscat (iianz, cpara, lpara, werte, maxw, lnew)
                        IF(ier_num==0) THEN
                           itype = NINT(werte(1))
                           CALL del_params (1, ianz, cpara, lpara, maxw)
                           c_name = cpara(1)(1:lpara(1))
                           CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                           IF(ier_num==0) THEN
                              ino = NINT(werte(1))
                              c_name = ' '
                           ELSEif(ier_num==-1 .AND. ier_typ==ER_FORT) THEN
                              CALL no_error   ! assume 
                           ENDIF
                           CALL recreate_connectivity(itype, ino, c_name)
                        ENDIF
                     ELSE
                        ier_num = -6
                        ier_typ = ER_COMM
                     ENDIF
                  ENDIF
!                                                                       
!     ----reset to no      connectivity definition 'reset'                 
!                                                                       
               ELSEIF (str_comp (befehl, 'reset', 3, lbef, 5) ) then 
                  CALL conn_do_set (code_res,zeile, lp) 
                  CALL deallocate_conn (conn_nmax)
!                                                                       
!     ----overwrite  a new connectivity definition 'set'                 
!                                                                       
               ELSEIF (str_comp (befehl, 'set', 2, lbef, 3) ) then 
                  CALL conn_do_set (code_set,zeile, lp) 
!                                                                       
!     ----show current parameters 'show'                                
!                                                                       
               ELSEIF (str_comp (befehl, 'show', 2, lbef, 4) ) then 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, length) 
                  IF ( ianz==0) THEN
                     CALL conn_show 
!                  CALL conn_test  For debug only
                  ELSE
                     IF (str_comp (cpara(1), 'connect', 3, lpara(1), 7) ) then
                        CALL del_params (1, ianz, cpara, lpara, maxw)
                        IF(str_comp (cpara(ianz), 'long',3, lpara(ianz), 4)) THEN
                           long = .true.
                           ianz = ianz - 1
                        ELSE
                           long = .false.
                        ENDIF
                        iianz = 1
                        CALL ber_params (iianz, cpara, lpara, werte, maxw) 
                        iatom = NINT(werte(1))
                        CALL del_params (1, ianz, cpara, lpara, maxw)
                        CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                        IF(ier_num/=0) THEN
                           c_name_l = MIN(256,lpara(1))
                           c_name   = cpara(1)(1:c_name_l)
                           ino      = 0
                           CALL no_error
                        ELSE                                               ! Success set to value
                           ino = nint (werte (1) ) 
                           c_name   = ' '
                           c_name_l = 1
                        ENDIF
                        CALL get_connectivity_identity( cr_iscat(iatom), ino, c_name, c_name_l)
                        CALL do_show_connectivity ( iatom, ino, c_name, long)
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ENDIF 
!                                                                       
!     ----Perform bond switching 'switch'                 
!                                                                       
               ELSEIF (str_comp (befehl, 'switch', 2, lbef, 6) ) then 
                  CALL bond_switch_para(zeile,lp)
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_COMM 
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
                  ier_msg(1) = ' Error occured in connectivity menu'
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
      prompt = orig_prompt
!                                                                       
   END SUBROUTINE conn_menu
!
   SUBROUTINE conn_show
!-                                                                      
!     Show connectivity definitions
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE atom_name
      USE prompt_mod
      IMPLICIT none 
!
      INTEGER   :: is
      INTEGER   :: i
!
!     IF(conn_mmc_sel == STATUS_ON) THEN
!        WRITE(output_io, 3000) 'selection', 'on'
!     ELSEIF(conn_mmc_sel == STATUS_OFF) THEN
!        WRITE(output_io, 3000) 'selection', 'off'
!     ELSE
!        WRITE(output_io, 3000) 'selection', 'ignore'
!     ENDIF
!     IF(conn_mmc_ene == STATUS_ON) THEN
!        WRITE(output_io, 3000) 'energy', 'on'
!     ELSEIF(conn_mmc_ene == STATUS_OFF) THEN
!        WRITE(output_io, 3000) 'energy', 'off'
!     ELSE
!        WRITE(output_io, 3000) 'energy', 'ignore'
!     ENDIF
      exist_def: IF ( ALLOCATED(def_main)) THEN    ! Are there any definitions
        scats: DO is=0,maxscat                     ! Loop over all atom types
           IF ( .NOT. ASSOCIATED(def_main(is)%def_liste)) THEN  ! This type has no def.s
              CYCLE scats
           ENDIF
           def_temp => def_main(is)%def_liste
           WRITE(output_io, 1000) at_name(is)
           DO
              IF ( .NOT. ASSOCIATED(def_temp)) THEN  ! This type has no more def.s
                 CYCLE scats
              ENDIF
              WRITE(output_io, 2000) def_temp%valid_id, &
                  def_temp%def_name(1:def_temp%def_name_l), def_temp%valid_no,   &
                  (at_name  (def_temp%valid_types(i)),i=1,def_temp%valid_no)
!                 (cr_at_lis(def_temp%valid_types(i)),                                &
!                            def_temp%valid_types(i) ,i=1,def_temp%valid_no)
              WRITE(output_io, 2300) def_temp%intend_no,def_temp%def_rmin,&
                                                        def_temp%def_rmax
!             IF(def_temp%mmc_sel==STATUS_ON) THEN
!                WRITE(output_io, 4000) 'selection', 'on'
!             ELSEIF(def_temp%mmc_sel==STATUS_OFF) THEN
!                WRITE(output_io, 4000) 'selection', 'off'
!             ELSE
!                WRITE(output_io, 4000) 'selection', 'ignore'
!             ENDIF
!             IF(def_temp%mmc_ene==STATUS_ON) THEN
!                WRITE(output_io, 4000), 'energy', 'on'
!             ELSEIF(def_temp%mmc_ene==STATUS_OFF) THEN
!                WRITE(output_io, 4000), 'energy', 'off'
!             ELSE
!                WRITE(output_io, 4000), 'energy', 'ignore'
!             ENDIF
              def_temp => def_temp%def_next
           ENDDO
        ENDDO scats
      ELSE exist_def                               ! No def.s exist
         WRITE(output_io, 7000) 
      ENDIF exist_def
!
1000  FORMAT(' Central atom type       : ',a9)
2000  FORMAT('     Def.no; Name; No.of.types; Types : ',i4,1x, a,1x,i4,1x, &
             ': ',20(a9:,',',2x))
!            20(a4,'(',i4,')',2x))
2300  FORMAT('     Max neig, Bond length range', 4x,i8,2x,f8.4, 2x, f8.4)
7000  FORMAT(' No connectivity definitions set')
!
   END SUBROUTINE conn_show
!
!
   SUBROUTINE conn_test
!-                                                                      
!     Mainly used while developing code, tests the connectivity
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      IMPLICIT none 
!
      INTEGER    :: i
!
!
      IF ( ALLOCATED(at_conn) ) THEN
         atoms:DO i=1,cr_natoms
        IF ( ASSOCIATED(at_conn(i)%liste) ) THEN      ! If neighborhood was created
          hood_temp => at_conn(i)%liste               ! point to the first neighborhood
!         hood_head => at_conn(i)%liste
          DO WHILE ( ASSOCIATED(hood_temp) )          ! While there are further neighborhood
             temp => hood_temp%nachbar                    ! point to the first neighbor within current neigborhood
             head => hood_temp%nachbar
             head => at_conn(i)%liste%nachbar             ! head points to the first NEIGHBOR
             DO WHILE ( ASSOCIATED(temp) )               ! While there are further neighbors
               temp => temp%next                         ! Point to next neighbor
               head => temp                              ! also point to next neighbor
             END DO
             hood_temp => hood_temp%next_neighborhood     ! Point to next neighborhood
          END DO
        END IF
         ENDDO atoms
      ENDIF
!
   END SUBROUTINE conn_test
!
!
   SUBROUTINE do_show_connectivity ( iatom, idef, c_name, long )
!-                                                                      
!     Shows the connectivity no. idef around atom no iatom
!+                                                                      
      USE crystal_mod
      USE atom_name
!     USE modify_mod
      USE prop_char_mod
      USE metric_mod
      USE param_mod 
      USE prompt_mod 
      USE lib_f90_allocate_mod
      IMPLICIT none 
!                                                                       
!
      INTEGER          , INTENT(IN)    :: iatom
      INTEGER          , INTENT(INOUT) :: idef
      CHARACTER(LEN=*) , INTENT(IN)    :: c_name
      LOGICAL          , INTENT(IN)    :: long
!
      CHARACTER (LEN=9)          :: at_name_d
      CHARACTER (LEN=32)         :: c_property
      INTEGER                    :: is1
      INTEGER, DIMENSION(:), ALLOCATABLE :: c_list
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: c_offs
      INTEGER                    :: natoms
!                                                                       
      INTEGER                    :: i, j
      INTEGER                    :: length
      INTEGER                    :: n_res
      REAL   , DIMENSION(3)      :: u, v
      REAL                       :: distance
!
      is1 = cr_iscat(iatom)
      u   = cr_pos(:,iatom)
      CALL get_connectivity_list (iatom, is1, idef, c_list, c_offs, natoms )
!
      WRITE(output_io,1000) iatom, at_name(is1), idef,c_name(1:LEN_TRIM(c_name)), natoms
      IF ( natoms > 0 ) THEN
        DO j=1,(natoms-1)/6+1
           WRITE(output_io, 1100) (c_list(i),i=(j-1)*6+1, MIN((j-1)*6 + 6,natoms))
        ENDDO
      ELSE
         WRITE(output_io, 1200)
      ENDIF
      IF( natoms > MAXPAR_RES) THEN
        n_res = MAX(natoms,MAXPAR_RES,CHEM_MAX_NEIG)
        CALL alloc_param(n_res)
        MAXPAR_RES = n_res
      ENDIF
      IF ( natoms>0 .AND. natoms <= MAXPAR_RES ) THEN
         res_para(0) = FLOAT(natoms)
         res_para(1:natoms) = FLOAT(c_list(1:natoms))
      ENDIF
!
      IF(long) THEN
         DO j= 1,natoms
            i = c_list(j)
            at_name_d = at_name(cr_iscat(i))
            CALL char_prop_1 (c_property, cr_prop (i), length) 
            WRITE (output_io, 3010) at_name_d, cr_pos (1, i), cr_pos (2,&
            i), cr_pos (3, i), cr_dw (cr_iscat (i) ), c_property (1:    &   !! WORK OCC
            length)
            v(:) = cr_pos(:,i) + c_offs(:,j)
            distance = do_blen(.TRUE., u,v)
            WRITE (output_io, 3020) c_offs(:,j), distance
         ENDDO
      ENDIF
!
      IF(ALLOCATED(c_list)) THEN
         DEALLOCATE(c_list)
      ENDIF
!
1000  FORMAT( ' Connectivity for atom No. ',I6,' of type ', &
             a9,' No : ',i4,1x,a,/,                         &
              '      Neighbor number:', i3)
1100  FORMAT( '      Neighbors are        ',20(i6:,2x))
1200  FORMAT( '      Atom has no neighbours')
3010  FORMAT(1x,a9,3(2x,f12.6),4x,f10.6,2x,a) 
3020  FORMAT( 3x  ,3(8x,i6   ),9x,f12.6    ) 
!
      END SUBROUTINE do_show_connectivity 
!
!
   SUBROUTINE bond_switch_para(zeile,length)
!
   USE discus_config_mod 
   USE crystal_mod 
   USE chem_aver_mod
   USE ber_params_mod
   USE errlist_mod 
   USE get_params_mod
   USE param_mod 
!
   IMPLICIT NONE
!
   CHARACTER (LEN=*), INTENT(INOUT) :: zeile
   INTEGER          , INTENT(INOUT) :: length
!
   INTEGER, PARAMETER                   :: MAXW = 3
   CHARACTER(LEN=1024), DIMENSION(MAXW) :: cpara
   INTEGER            , DIMENSION(MAXW) :: lpara
   REAL               , DIMENSION(MAXW) :: werte
   INTEGER                                 :: ianz , iianz
!
   INTEGER                            :: iatom   ! Selected atom number
   INTEGER                            :: natom   ! Dummy atom number
   INTEGER                            :: i       ! Dummy atom number
   INTEGER                            :: iscat   ! Selected atom type
   INTEGER                            :: success ! Selected atom type
   INTEGER                            :: ino     ! Connectivity def. no.
   CHARACTER(LEN=256)                 :: c_name  ! Connectivity name
!
   REAL       :: ran1
!
   iscat  = 0
   ino    = 0
   c_name = ' '
   CALL get_params (zeile, ianz, cpara, lpara, maxw, length) 
   IF ( ianz==0) THEN
!
!     No parameters, choose any atom, use connectivity number one
!
      iatom  = INT(ran1(idum)*cr_natoms) + 1
      iscat  = cr_iscat(iatom)
      ino    = 1
      c_name = ' '
   ELSE
      iianz = 1
      CALL ber_params (iianz, cpara, lpara, werte, maxw) 
      IF(ier_num==0) THEN
         IF(ianz==2) THEN   ! Atom type is specified 
            iscat = NINT(werte(1))             ! use this ytom type
            res_para(:) = 0
            CALL chem_elem(.FALSE.)            ! Check if any atoms are present
            IF(iscat<0) THEN
               ier_num = -27
               ier_typ = ER_APPL
               RETURN
            ENDIF
            IF(res_para(iscat+1)>0.0) THEN     ! Do we have atoms iscat?
               natom = INT(ran1(idum)*cr_natoms*res_para(iscat+1))+1 ! random choice
               iatom = 0
search:        DO i=1, cr_natoms               ! loop until we find the natoms atom
                  IF(cr_iscat(i)==iscat) THEN  !      of type iscat
                     iatom = iatom + 1
                     IF(iatom==natom) THEN
                        iatom = i              ! Found, this is the abolute atom number
                        EXIT search
                     ENDIF
                  ENDIF
               ENDDO search
            ELSE
               ier_num = -27
               ier_typ = ER_APPL
               RETURN
            ENDIF
         ELSEIF(IANZ==3) THEN
            iatom = NINT(werte(1))
            IF(iatom<1 .OR.iatom>cr_natoms) THEN
               ier_num = -19
               ier_typ = ER_APPL
               ier_msg(1) = 'atom number must be > 0 and < n[1]'
               RETURN
            ENDIF
            CALL del_params (1, ianz, cpara, lpara, maxw)
            CALL ber_params (iianz, cpara, lpara, werte, maxw) 
            IF(ier_num/=0) THEN
               ier_num = -6
               ier_typ = ER_FORT
               RETURN
            ENDIF
            iscat = NINT(werte(1))
            IF(cr_iscat(iatom)/= iscat) THEN
               ier_num = -6
               ier_typ = ER_FORT
               ier_msg(1) = 'Atom type does not match selected atom'
               RETURN
            ENDIF
         ELSE
            ier_num = -6
            ier_typ = ER_COMM
            ier_msg(1) = 'Switch command needs 2 or 3 parameters'
            RETURN
         ENDIF
         CALL del_params (1, ianz, cpara, lpara, maxw)
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         IF(ier_num/=0) THEN
            c_name   = cpara(1)(1:lpara(1))
            ino      = 0
            CALL no_error
         ELSE                                               ! Success set to value
            ino = nint (werte (1) ) 
            IF(ino < 0 ) THEN
               ier_num = -6
               ier_typ = ER_COMM
               ier_msg(1) = 'Connectivity number needs to be positive'
               RETURN
            ENDIF
            c_name   = ' '
         ENDIF
      ENDIF
   ENDIF
   IF(ier_num == 0) THEN
      IF(iatom > 0 .AND. (ino>0.OR.c_name/=' ')) THEN
         CALL do_bond_switch (iatom, ino, c_name, success) 
         IF(success == 0) THEN
            CALL no_error
         ELSEIF(success == -1) THEN
            ier_num = -134
            ier_typ = ER_APPL
         ELSEIF(success == -2) THEN
            ier_num = -135
            ier_typ = ER_APPL
         ELSEIF(success == -3) THEN
            ier_num = -136
            ier_typ = ER_APPL
         ENDIF
      ENDIF
   ENDIF
!
   END SUBROUTINE bond_switch_para
!
!
   SUBROUTINE do_bond_switch (jatom, ino, c_name, success)
!-                                                                      
!     Get the list of neighbors for central atom jatom of type is1
!+                                                                      
   USE discus_config_mod 
   USE crystal_mod 
   USE metric_mod
!
   IMPLICIT none 
!
   INTEGER           , INTENT(IN)     :: jatom   ! central atom number
   INTEGER           , INTENT(INOUT)  :: ino     ! Connectivity def. no.
   CHARACTER(LEN=256), INTENT(INOUT)  :: c_name  ! Connectivity name
   INTEGER           , INTENT(OUT)    :: success
!
   INTEGER, DIMENSION(:),   ALLOCATABLE :: c_list  ! List of all neighbors 
   INTEGER, DIMENSION(:,:), ALLOCATABLE :: c_offs  ! Offsets from periodic boundary
   INTEGER, DIMENSION(:),   ALLOCATABLE :: s_list  ! List of all neighbors 
   INTEGER, DIMENSION(:,:), ALLOCATABLE :: s_offs  ! Offsets from periodic boundary
   INTEGER, DIMENSION(:),   ALLOCATABLE :: j_list  ! List of all neighbors 
   INTEGER, DIMENSION(:,:), ALLOCATABLE :: j_offs  ! Offsets from periodic boundary
   INTEGER, DIMENSION(:),   ALLOCATABLE :: k_list  ! List of all neighbors 
   INTEGER, DIMENSION(:,:), ALLOCATABLE :: k_offs  ! Offsets from periodic boundary
   REAL   , DIMENSION(:,:), ALLOCATABLE :: test_par  ! Test for parallel vectors to the neighbors
   INTEGER, DIMENSION(2)                :: test_max  ! Location ov the maximum scalar product
   INTEGER                              :: c_natoms  ! number of atoms in connectivity list
   INTEGER                              :: s_natoms  ! number of atoms in connectivity list
   INTEGER                              :: j_natoms  ! number of atoms in connectivity list
   INTEGER                              :: k_natoms  ! number of atoms in connectivity list
!
   TYPE (NEIGHBORHOOD), POINTER       :: hood_j
   TYPE (NEIGHBORHOOD), POINTER       :: hood_k
   INTEGER    :: c_neig, s_neig, c_ex, s_ex
   INTEGER, DIMENSION(3) :: t_offs, in_offs
   INTEGER, DIMENSION(3) :: c2ex_old_offs, s2ex_old_offs
   INTEGER, DIMENSION(3) :: c2ex_new_offs, s2ex_new_offs
   REAL   , DIMENSION(3) :: u,v, un, vn
   REAL                  :: distance, unn, vnn
   INTEGER               :: t_ex, in_ex, in_ref
   INTEGER    :: katom, j_ex, k_ex
   INTEGER    :: i,k
   INTEGER :: j_in_j, k_in_k
   LOGICAL :: j_ex_in_s_list
   LOGICAL :: k_ex_in_c_list
   REAL       :: ran1
!
   NULLIFY(hood_j)
   NULLIFY(hood_k)
   c_natoms = 0
!
!
   IF ( ALLOCATED(at_conn) ) THEN
      CALL get_connect_pointed(hood_central, jatom, ino, c_name, c_list, c_offs, c_natoms, success)
      IF(success/=0) RETURN
! We should have found a central atom and its connectivity list.
! Randomly pick a neighbor atom
      c_neig = INT(ran1(idum)*c_natoms) + 1
      katom  = c_list(c_neig) 
      CALL get_connect_pointed(hood_second, katom, ino, c_name, s_list, s_offs, s_natoms, success)
      IF(success/=0) RETURN
!
      u(:) = cr_pos(:,jatom)
      v(:) = cr_pos(:,katom) + c_offs(:,c_neig)
      distance = do_blen(.TRUE., u,v)
      IF(distance > 2.0) THEN       !   DEVELOPMENT !! MUST BE A USER DEFINED VALUE  !!
!do i=1, 100
!write(*,'(a,2i4,f8.3)') ' ATOMS far apart ', jatom, katom, distance
!ENDDO
         RETURN
      ENDIF
!
      s_neig = 0                                          ! Central not yet found
search_c: DO i=1, s_natoms
         IF(s_list(i) == jatom) THEN
            s_neig = i
            EXIT search_c
         ENDIF
      ENDDO search_c
      IF(s_neig == 0) THEN
!        write(output_io,*) ' CENTRAL atom is NOT a neighbor to second', jatom, c_list(c_neig)
!        ier_num = -6
!        ier_typ = ER_FORT
         DEALLOCATE(s_list)
         DEALLOCATE(c_list)
         DEALLOCATE(s_offs)
         DEALLOCATE(c_offs)
         RETURN
      ENDIF
!
! Test for parallel
!
!write(*,'(a,i5,3f6.2,2x, 3f6.2,i5)') ' jatom, pos', jatom, cr_pos(:,jatom), u, c_neig
!write(*,'(a,i5,3f6.2,2x, 3f6.2,i5)') ' katom, pos', katom, cr_pos(:,katom), v, s_neig
!i=1
!c_name = 'c_first'
!call do_show_connectivity ( jatom, i, c_name, .true. )
!call do_show_connectivity ( katom, i, c_name, .true. )
      ALLOCATE(test_par(1:c_natoms-1,1:s_natoms-1))
      test_par(:,:) = 0.0
      DO i=1, c_natoms-1
         c_ex  = MOD(i+c_neig-1, c_natoms)+1
         un(:) = cr_pos(:,c_list(c_ex)) + c_offs(:,c_ex)  -  u(:)  ! Vector jatom ==> neig
         unn   = sqrt (skalpro (un, un, cr_gten) )
         DO k=1, s_natoms-1
            s_ex  = MOD(k+s_neig-1, s_natoms)+1
            vn(:) = v(:) - c_offs(:, c_neig) - cr_pos(:,s_list(s_ex)) - s_offs(:,s_ex) ! Vector neig ==> katom
            vnn   = sqrt (skalpro (vn, vn, cr_gten) )
            test_par(i,k) =  skalpro (un, vn, cr_gten)/(unn*vnn)
!write(*,'(4i4,3f6.2,2x,3f6.2)') c_list(c_ex), s_list(s_ex), c_ex, s_ex, &
!   cr_pos(:,c_list(c_ex)) + c_offs(:,c_ex),                             &
!   cr_pos(:,s_list(s_ex)) + s_offs(:,s_ex)
!write(*,'(16x,3f6.2,2x, 3f6.2,2x, f6.2)') un, vn, test_par(i,k)
         ENDDO
      ENDDO
!
!     Choose the pair that is as parallel as possible
!
      test_max = MAXLOC(test_par)
      i= test_max(1)
      k= test_max(2)
      c_ex = MOD(i+c_neig-1, c_natoms ) + 1
      s_ex = MOD(k+s_neig-1, s_natoms ) + 1
      DEALLOCATE(test_par)
!write(*,'(a,2i5,4i4,2i5)') ' jatom, katom', jatom, katom, test_max, c_ex, s_ex, c_list(c_ex), s_list(s_ex)
!
!     c_ex = MOD(INT(ran1(idum)*(c_natoms-1)) + (c_neig), c_natoms ) + 1
!     s_ex = MOD(INT(ran1(idum)*(s_natoms-1)) + (s_neig), s_natoms ) + 1
      j_ex = c_list(c_ex)
      k_ex = s_list(s_ex)
!
!     TEST for distance original to new partner and second to new partner
!
      vn(:) = cr_pos(:, k_ex) + c_offs(:,c_neig) + s_offs(:,s_ex)
      distance = do_blen(.TRUE., u, vn)
      IF(distance .GT. 3.00) THEN       !   DEVELOPMENT !! MUST BE A USER DEFINED VALUE  !! 
!write(*,'(a,3f6.2,2x,3f6.2, 2x, f6.2)') ' Neighbor to A far apart rejected',u,vn, distance
         RETURN
      ENDIF
      un(:) = cr_pos(:, j_ex) + c_offs(:,c_ex)  !+ s_offs(:,s_neig) !+ c_offs(:,c_ex)
      distance = do_blen(.TRUE., v, un)
      IF(distance .GT. 3.00) THEN       !   DEVELOPMENT !! MUST BE A USER DEFINED VALUE  !! 
!write(*,'(a,3f6.2,2x,3f6.2, 2x, f6.2)') ' Neighbor to B far apart rejected',v,un, distance
         RETURN
      ENDIF
!return
!
      j_ex_in_s_list = .FALSE.
      DO i=1, s_natoms
         IF(s_list(i) == j_ex) THEN
            j_ex_in_s_list = .TRUE.
         ENDIF
      ENDDO
      k_ex_in_c_list = .FALSE.
      DO i=1, c_natoms
         IF(c_list(i) == k_ex) THEN
            k_ex_in_c_list = .TRUE.
         ENDIF
      ENDDO
      if(j_ex_in_s_list .OR. k_ex_in_c_list) THEN
!        write(output_io,*) 'EXCHANGE WOULD PRODUCE DOUBLE PARTNER '
!ier_num = -6
!ier_typ = ER_FORT
         DEALLOCATE(s_list)
         DEALLOCATE(c_list)
         DEALLOCATE(s_offs)
         DEALLOCATE(c_offs)
         RETURN
      ENDIF
      if(j_ex == k_ex) then
!write(*,*) 'Trying to exchange identical neighbors'
!ier_num = -6
!ier_typ = ER_FORT
         DEALLOCATE(s_list)
         DEALLOCATE(c_list)
         DEALLOCATE(s_offs)
         DEALLOCATE(c_offs)
         RETURN
      ENDIF

      IF(.NOT.(j_ex_in_s_list .OR. k_ex_in_c_list .OR. j_ex == k_ex )) THEN ! SUCCESS
!
!     Get neighbors for the atoms that will be exchanged
         CALL get_connect_pointed(hood_j, j_ex, ino, c_name, j_list, j_offs, j_natoms, success)
         IF(success/=0) RETURN
         CALL get_connect_pointed(hood_k, k_ex, ino, c_name, k_list, k_offs, k_natoms, success)
         IF(success/=0) RETURN
         j_in_j = 0                                          ! Central not yet found
search_k: DO i=1, j_natoms
           IF(j_list(i) == jatom) THEN
              j_in_j = i
              EXIT search_k
           ENDIF
         ENDDO search_k
         k_in_k = 0                                          ! Central not yet found
search_j: DO i=1, k_natoms
            IF(k_list(i) == katom) THEN
               k_in_k = i
               EXIT search_j
            ENDIF
         ENDDO search_j
!
! now do the exchange
!
         c2ex_old_offs(:) = c_offs(:,c_ex)   ! offset central to original neig
         s2ex_old_offs(:) = s_offs(:,s_ex)   ! offset second  to original neig
         c2ex_new_offs(:) = c_offs(:,c_neig) + s2ex_old_offs(:) ! offs central to new
         s2ex_new_offs(:) = s_offs(:,s_neig) + c2ex_old_offs(:) ! offs second  to new
!     Modify second atom
         in_ref     = s_list(s_ex)    ! Find this partner
         in_ex      = c_list(c_ex)    ! store this atom number instead of old
         in_offs(:) = s2ex_new_offs(:)  ! store this offset instead of old
         CALL do_exchange(hood_second%nachbar, in_ref, in_ex, in_offs, t_ex, t_offs)
!     Modify central atom
         in_ref     = c_list(c_ex)    ! Find this partner
         in_ex      = t_ex            ! store this atom number instead of old
         in_offs(:) = c2ex_new_offs(:)   ! store this offset instead of old
!     Modify new neig to second atom
         CALL do_exchange(hood_central%nachbar, in_ref, in_ex, in_offs, t_ex, t_offs)
         in_ref     = jatom           ! Find this partner
         in_ex      = katom           ! store this atom number instead of old
         in_offs(:) = -s2ex_new_offs(:)! store this offset instead of old
         CALL do_exchange(hood_j%nachbar, in_ref, in_ex, in_offs, t_ex, t_offs)
!     Modify new neig to central atom
         in_ref     = katom           ! Find this partner
         in_ex      = jatom           ! store this atom number instead of old
         in_offs(:) = -c2ex_new_offs(:)! store this offset instead of old
         CALL do_exchange(hood_k%nachbar, in_ref, in_ex, in_offs, t_ex, t_offs)
!
!
      ENDIF
   ENDIF
   DEALLOCATE(s_list)
   DEALLOCATE(c_list)
   DEALLOCATE(s_offs)
   DEALLOCATE(c_offs)
   NULLIFY(hood_j)
   NULLIFY(hood_k)
!
   END SUBROUTINE do_bond_switch
!
   SUBROUTINE get_connect_pointed(hood_p, jatom, ino, c_name, c_list, c_offs, c_natoms, success)
!
   IMPLICIT NONE
!
!
   TYPE (NEIGHBORHOOD)    , POINTER                  :: hood_p
   INTEGER                             , INTENT(IN)  :: jatom   ! central atom number
   INTEGER                             , INTENT(INOUT)  :: ino     ! Connectivity def. no.
   CHARACTER(LEN=256)                  , INTENT(INOUT)  :: c_name  ! Connectivity name
   INTEGER, DIMENSION(:),   ALLOCATABLE, INTENT(OUT) :: c_list  ! List of all neighbors 
   INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: c_offs  ! Offsets from periodic boundary
   INTEGER                             , INTENT(OUT) :: c_natoms  ! number of neigbor atoms
   INTEGER                             , INTENT(out) :: success   ! sucess flag
!
   TYPE (NEIGHBORS), POINTER  :: p_atoms
   INTEGER :: i, k
!
   i = jatom
   success = -1                                  ! No connectivity list exists
   IF(ASSOCIATED(at_conn(i)%liste)) THEN
      hood_p => at_conn(i)%liste                 ! point to the first neighborhood
      success = -2                               ! No neighborhood exists
search1:  DO WHILE ( ASSOCIATED(hood_p) )        ! While there are further neighborhood
         p_atoms => hood_p%nachbar               ! point to the first neighbor within current neigborhood
         success = -2                            ! No connectivity matches
         IF ( hood_p%neigh_type == ino  .OR. &
            hood_p%conn_name  == c_name   ) THEN ! This is the right neighborhood
            ino    = hood_p%neigh_type           ! Return actual number and name
            c_name = hood_p%conn_name
            c_natoms = hood_p%natoms             ! Store number of neigboring atoms
            IF(c_natoms==0) THEN
               success = -3                      ! List is empty
               RETURN
            ENDIF
            IF(ALLOCATED(c_list)) THEN           ! Just in case, if the list exists already
               DEALLOCATE(c_list)
               ALLOCATE(c_list(0:c_natoms))
            ELSE
               ALLOCATE(c_list(0:c_natoms))
            ENDIF
            IF(ALLOCATED(c_offs)) THEN
               DEALLOCATE(c_offs)
               ALLOCATE(c_offs(1:3,0:c_natoms))
            ELSE
               ALLOCATE(c_offs(1:3,0:c_natoms))
            ENDIF
            c_list = 0                           ! clear connectivity list
            k = 0
            DO WHILE ( ASSOCIATED(p_atoms) )     ! While there are further neighbors
               k         = k+ 1
               c_list(k) = p_atoms%atom_number   ! Add atoms
               c_offs(1,k) = p_atoms%offset(1)   ! and relative offsets
               c_offs(2,k) = p_atoms%offset(2)
               c_offs(3,k) = p_atoms%offset(3)
               p_atoms => p_atoms%next           ! Point to next neighbor
            END DO
            success =  0                         ! Finally, success
            EXIT search1                         ! End of connectivity list
         ENDIF
         hood_p => hood_p%next_neighborhood      ! Point to next neighborhood
      END DO search1                             ! end of search central atom
   ENDIF
!
   END SUBROUTINE get_connect_pointed
!
   SUBROUTINE  do_exchange(hood_start, in_ref, in_ex, in_offs, t_ex, t_offs)
!
   IMPLICIT NONE
   TYPE (NEIGHBORS), POINTER                    :: hood_start
   INTEGER                        , INTENT(IN ) :: in_ref  ! Reference atom
   INTEGER                        , INTENT(IN ) :: in_ex   ! input exchange partner
   INTEGER, DIMENSION(3         ) , INTENT(IN ) :: in_offs ! Input offsets from periodic boundary
   INTEGER                        , INTENT(OUT) :: t_ex    ! Output echange partner
   INTEGER, DIMENSION(3         ) , INTENT(OUT) :: t_offs  ! Output offsets from periodic boundary
!
   TYPE (NEIGHBORS), POINTER  :: p_atoms
!
   p_atoms => hood_start               ! point to the first neighbor within current neigborhood
search_s_ex: DO WHILE(ASSOCIATED(p_atoms))
      IF(p_atoms%atom_number == in_ref      ) THEN ! Found second exchange partner
         t_ex                = p_atoms%atom_number
         t_offs(:)           = p_atoms%offset(:)
         p_atoms%atom_number = in_ex
         p_atoms%offset(:)   = in_offs(:)
         EXIT search_s_ex
      ENDIF
      p_atoms => p_atoms%next           ! Point to next neighborhood
   ENDDO search_s_ex
   NULLIFY(p_atoms)
   END SUBROUTINE  do_exchange
!
SUBROUTINE conn_update(isel, shift)
!
IMPLICIT NONE
!
INTEGER              , INTENT(IN) :: isel
REAL   , DIMENSION(3), INTENT(IN) :: shift
!
CHARACTER(LEN=256)                   :: c_name  ! Connectivity name
INTEGER, DIMENSION(:),   ALLOCATABLE :: c_list  ! List of all neighbors 
INTEGER, DIMENSION(:,:), ALLOCATABLE :: c_offs  ! Offsets from periodic boundary
CHARACTER(LEN=256)                   :: j_name  ! Connectivity name
INTEGER, DIMENSION(:),   ALLOCATABLE :: j_list  ! List of all neighbors 
INTEGER, DIMENSION(:,:), ALLOCATABLE :: j_offs  ! Offsets from periodic boundary
INTEGER :: c_natoms
INTEGER :: j_natoms
INTEGER :: ino
INTEGER :: i, success
INTEGER :: iatom
INTEGER :: is
TYPE (NEIGHBORHOOD), POINTER         :: hood_c
TYPE (NEIGHBORHOOD), POINTER         :: hood_j
   TYPE (NEIGHBORS), POINTER         :: p_atoms
!
NULLIFY(hood_c)
NULLIFY(hood_j)
!
ino = 0
IF(.NOT.ALLOCATED(at_conn)) RETURN
IF(UBOUND(at_conn,1)==0) RETURN
IF(ASSOCIATED(at_conn(isel)%liste)) THEN     ! A connectivity list has been created
   is = cr_iscat(isel) 
   IF(ASSOCIATED(def_main(is)%def_liste)) THEN  ! This type has a definition
      def_temp => def_main(is)%def_liste
search_defs:      DO WHILE (ASSOCIATED(def_temp))           ! There are definitions to follow
         c_name = def_temp%def_name(1:def_temp%def_name_l)
         ino = 0
         CALL get_connect_pointed(hood_c, isel, ino, c_name, c_list, c_offs, c_natoms, success)
         IF(success/=0) RETURN
         DO i=1,c_natoms                          ! Update all offsets
            c_offs(:,i) = c_offs(:,i) + INT(shift(:))
         ENDDO
         p_atoms => hood_c%nachbar
         i = 1
         DO WHILE(ASSOCIATED(p_atoms))            ! Place into structure
            IF(p_atoms%atom_number == c_list(i)) THEN
               p_atoms%offset(:) = c_offs(:,i)
            ENDIF
            i = i + 1
            p_atoms => p_atoms%next
         ENDDO
         def_temp => def_temp%def_next
      ENDDO search_defs
   ENDIF
!DBG_PERIOD j= 1
!DBG_PERIOD call do_show_connectivity ( isel, j, c_name, .TRUE. )
!
!  Loop over all atoms, and find out if atom isel is a neighbor to any, if so update
!  Might have to be replaced by a connectivity list that indicates for a given atom
!  which other atoms have this listed as neighbor....
   DO iatom = 1, cr_natoms
      IF(ASSOCIATED(at_conn(iatom)%liste)) THEN     ! A connectivity list has been created
         is = cr_iscat(iatom) 
         IF(ASSOCIATED(def_main(is)%def_liste)) THEN  ! This type has a definition
            def_temp => def_main(is)%def_liste
search_def2:DO WHILE (ASSOCIATED(def_temp))           ! There are definitions to follow
               j_name = def_temp%def_name(1:def_temp%def_name_l)
               ino = 0
               CALL get_connect_pointed(hood_j, iatom, ino, j_name, j_list, j_offs, j_natoms, success)
               IF(success/=0) RETURN
               p_atoms => hood_j%nachbar
search_neig:   DO WHILE(ASSOCIATED(p_atoms))            ! Place into structure
                  IF(p_atoms%atom_number == isel     ) THEN
                     p_atoms%offset(:) = p_atoms%offset(:) - INT(shift(:))
                     EXIT search_neig
                  ENDIF
                  p_atoms => p_atoms%next
               ENDDO search_neig
               def_temp => def_temp%def_next
            ENDDO search_def2
         ENDIF
      ENDIF
   ENDDO
ENDIF
IF(ALLOCATED(c_list)) DEALLOCATE(c_list)
IF(ALLOCATED(c_offs)) DEALLOCATE(c_offs)
IF(ALLOCATED(j_list)) DEALLOCATE(j_list)
IF(ALLOCATED(j_offs)) DEALLOCATE(j_offs)
NULLIFY(hood_c)
NULLIFY(hood_j)
!
!
END SUBROUTINE conn_update
!
END MODULE conn_mod
