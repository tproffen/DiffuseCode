MODULE conn_mod
!
USE conn_def_mod
USE crystal_mod
!
USE errlist_mod
!
IMPLICIT none
!
PRIVATE
PUBLIC  conn_menu              !  Main menu to interact with user
PUBLIC  get_connectivity_list  !  Read out the actual list of atoms around a central atom
PUBLIC  get_connectivity_identity ! Identify a connectivity definition
PUBLIC  do_show_connectivity   !  Show the current definitions
PRIVATE allocate_conn_list     !  Allocate memory for the connectivity
PRIVATE deallocate_conn        !  Free memory
PRIVATE create_connectivity    !  Create the actual list of neighbors around each atom
PRIVATE conn_do_set            !  Set parameters for the connectivity definitions
PRIVATE conn_show              !  Main show routine
PRIVATE conn_test              !  While developing, a routine to test functionality
!
INTEGER, PARAMETER              :: MAX_ATOM=10
!
! TYPE NEIGHBORS  is a linear chain to be filled with the actual neighboring atoms
!                 currently only the neighbor number is stored no further info, as
!                 I do not want to update this list all the time.
TYPE :: NEIGHBORS
   INTEGER                      :: atom_number
   TYPE (NEIGHBORS), POINTER    :: next
END TYPE
!
! TYPE NEIGHBORHOOD is a linear chain of the possible neighborhoods.
!                   Each node branches off to one side with the actual neighboring atoms
TYPE :: NEIGHBORHOOD
   INTEGER                      :: central_number     ! absolute number of central atom
   INTEGER                      :: central_type       ! central atom is of this type
   INTEGER                      :: neigh_type         ! this neighbors belongs to this definition
   CHARACTER (LEN=256)          :: conn_name          ! Connectivity name
   INTEGER                      :: conn_name_l        ! Connectivity name length
   REAL                         :: distance_min       ! minimum distance to neighbors
   REAL                         :: distance_max       ! maximum distance to neighbors
   TYPE (NEIGHBORS), POINTER    :: nachbar            ! The actual list of neighboring atoms
   TYPE (NEIGHBORHOOD), POINTER :: next_neighborhood  ! A next neighborhood
END TYPE
!
! TYPE MAIN_LIST is a structure that contains info on the central atom, as
!                well as a pointer to the NEIGHBORHOOD
TYPE :: MAIN_LIST
   INTEGER                      :: number
   TYPE (NEIGHBORHOOD), POINTER :: liste
END TYPE
!
! In order to have FAST access to the neighborhood of any atom, an
! allocatable array is defined. Entry at_conn(i) gives access to the 
! neighborhood of atom i
TYPE (main_list), DIMENSION(:), ALLOCATABLE :: at_conn
!
! (temporary) pointers of TYPE NEIGHBORS. This allows to move along the 
! neighbors in an individual neighborhood.
TYPE (NEIGHBORS), POINTER       :: head, tail, temp
!
! (temporary) pointers of TYPE NEIGHBORHOOD. This allows to move along the 
! neighborhoods of an individual atom.
TYPE (NEIGHBORHOOD), POINTER       :: hood_head
TYPE (NEIGHBORHOOD), POINTER       :: hood_temp
!
LOGICAL                            :: conn_status = .false.
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
   USE modify_mod
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
   INTEGER              :: j,i,jj
   INTEGER              :: is,js ! dummies for scattering types
   INTEGER              :: ianz
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
               NULLIFY (hood_temp%next_neighborhood)        ! No further NEIGHBORHOODs
               hood_temp%central_number = i                 ! Just set central atom no.
               hood_temp%central_type   = cr_iscat(i)
               hood_temp%neigh_type     = def_temp%valid_id   ! Set definition type number
               hood_temp%conn_name      = def_temp%def_name   ! Set name from definition type
               hood_temp%conn_name_l    = def_temp%def_name_l ! Set name length from definition type
               NULLIFY (hood_temp%nachbar)                  ! Initially there are no NEIGHBORS
!
               ALLOCATE (hood_temp%nachbar)                 ! create the first NEIGHBOR slot
               j = atom_env(1)
               tail => hood_temp%nachbar                    ! tail points to the first NEIGHBOR
               tail%atom_number = j                         ! I store the atom_no of the neighbor
               NULLIFY (tail%next)                          ! No further neighbors
!
               DO j = 2, atom_env(0)                        ! Add all neighbors to list
                  ALLOCATE (tail%next)                      ! create a further NEIGHBOR
                  tail => tail%next                         ! reassign tail to new end of list
                  tail%atom_number = atom_env(j)            ! I store the atom_no of the neighbor
                  NULLIFY (tail%next)                       ! No further neighbors
               ENDDO
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
   SUBROUTINE conn_do_set ( code, zeile, length)
!-                                                                      
!     Set the parameters for the connectivity
!+                                                                      
      USE allocate_appl_mod 
      USE config_mod 
      USE crystal_mod 
      USE modify_mod
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
!
      CHARACTER(LEN=1024), DIMENSION(MAX(MIN_PARA,MAXSCAT+4)) :: cpara
      INTEGER            , DIMENSION(MAX(MIN_PARA,MAXSCAT+4)) :: lpara
      REAL               , DIMENSION(MAX(MIN_PARA,MAXSCAT+4)) :: werte
!
      CHARACTER(LEN=1024), DIMENSION(2)  :: ccpara
      INTEGER            , DIMENSION(2)  :: llpara
      REAL               , DIMENSION(2)  :: wwerte
!
      INTEGER             :: ianz         ! number of command line parameters
      INTEGER             :: iianz        ! dummy number
      INTEGER             :: max_scat     ! dummy number
      INTEGER             :: i,j          ! dummy number
      INTEGER             :: is1          ! first atom type
      INTEGER             :: is2          ! second atom type
      INTEGER             :: temp_id      ! temporary definition ID
      INTEGER             :: work_id      ! ID of the definition to change/delete
      CHARACTER(LEN=256)  :: work_name    ! Name of the definition to change/delete
      INTEGER             :: work_name_l  ! Length of name for the definition to change/delete
      LOGICAL             :: lnew         ! require atom type to exist
      INTEGER             :: all_status   ! Allocation status
      REAL                :: rmin         ! minimum bond distance
      REAL                :: rmax         ! maximum bond distance
!                                                                       
      LOGICAL :: str_comp 
      REAL    :: berechne 
!
      maxw = MAX(MIN_PARA, MAXSCAT+4)
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
         CALL del_params (1, ianz, cpara, lpara, maxw) 
      ELSE
         work_id     = -1
         work_name   = cpara(ianz)(1:lpara(ianz))
         work_name_l = lpara(ianz)
         ianz        = ianz - 1
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
               def_temp%valid_no = ianz                          ! Set number of neighb or types
               def_temp%def_rmin = rmin                          ! Set bond length limits
               def_temp%def_rmax = rmax                          ! Set bond length limits
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
            def_temp%valid_no   = ianz                        ! Set number of neighb or types
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
            def_temp%valid_no   = ianz                        ! Set number of neighb or types
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
      USE config_mod 
      USE crystal_mod 
!
      USE doact_mod 
      USE learn_mod 
      USE macro_mod 
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER maxw 
      LOGICAL lnew, lold 
!                                                                       
      PARAMETER (lnew = .true., lold = .false.) 
!                                                                       
      CHARACTER(5) befehl 
      CHARACTER(50) prom 
      CHARACTER(1024) line, zeile, cpara (MAXSCAT) 
      INTEGER lpara (MAXSCAT), lp, length, lbef 
      INTEGER indxg, ianz, i 
      INTEGER indxc 
      LOGICAL lend, lspace 
      LOGICAL lselect 
      REAL werte (MAXSCAT) 
!                                                                       
      INTEGER len_str 
      LOGICAL str_comp 
!                                                                       
      maxw = MAXSCAT
      lend = .false. 
      CALL no_error 
!                                                                       
      DO while (.not.lend) 
      prom = prompt (1:len_str (prompt) ) //'/conn' 
      CALL get_cmd (line, length, befehl, lbef, zeile, lp, prom) 
      IF (ier_num.eq.0) then 
         IF (line.ne.' '.and.line (1:1) .ne.'#') then 
!                                                                       
!     ----search for "="                                                
!                                                                       
            indxg = index (line, '=') 
            IF (indxg.ne.0.and.                                      &
                .not. (str_comp (befehl, 'echo', 2, lbef, 4) ) .and. &
                .not. (str_comp (befehl, 'syst', 2, lbef, 4) ) .and. &
                .not. (str_comp (befehl, 'help', 2, lbef, 4) .or.    &
                       str_comp (befehl, '?   ', 2, lbef, 4) ) ) then                                              
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
                     IF (str_comp (cpara(1), 'connect', 4, lpara(1), 7) ) then
                        CALL del_params (1, ianz, cpara, lpara, maxw)
                        CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                        CALL do_show_connectivity ( NINT(werte(1)), NINT(werte(2)))
                     ELSE 
                        ier_num = - 6 
                        ier_typ = ER_COMM 
                     ENDIF 
                  ENDIF 
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
 9999 CONTINUE 
!
   END SUBROUTINE conn_menu
!
   SUBROUTINE conn_show
!-                                                                      
!     Show connectivity definitions
!+                                                                      
      USE config_mod 
      USE crystal_mod 
      USE atom_name
      USE prompt_mod
      IMPLICIT none 
!
      INTEGER   :: is
      INTEGER   :: i
      CHARACTER(LEN=9)           :: at_name_i
!     CHARACTER(LEN=9), EXTERNAL :: at_name
!
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
              WRITE(output_io, 2300) def_temp%def_rmin,def_temp%def_rmax
              def_temp => def_temp%def_next
           ENDDO
        ENDDO scats
      ELSE exist_def                               ! No def.s exist
         WRITE(output_io, 7000) 
      ENDIF exist_def
!
1000  FORMAT(' Central atom type       : ',a9)
2000  FORMAT('     Def.no; Name; No.of.neigh; Types : ',i4,1x, a,1x,i4,1x, &
             ': ',20(a9:,',',2x))
!            20(a4,'(',i4,')',2x))
2300  FORMAT('     Bond length range',23x     ,f8.4, 2x, f8.4)
7000  FORMAT(' No connectivity definitions set')
!
   END SUBROUTINE conn_show
!
!
   SUBROUTINE conn_test
!-                                                                      
!     Mainly used while developing code, tests the connectivity
!+                                                                      
      USE config_mod 
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
   SUBROUTINE get_connectivity_list (jatom, is1, ino, maxw, c_list, natoms )
!-                                                                      
!     Get the list of neighbors for central atom jatom of type is1
!+                                                                      
      USE config_mod 
      USE crystal_mod 
      IMPLICIT none 
!
      INTEGER, INTENT(IN)  :: jatom   ! central atom number
      INTEGER, INTENT(IN)  :: is1     ! central atom type
      INTEGER, INTENT(IN)  :: ino     ! Connectivity def. no.
      CHARACTER(LEN=256)   :: c_name  ! Connectivity name
      INTEGER              :: c_name_l! Connectivity name length
      INTEGER, INTENT(IN)  :: maxw    ! Size of array c_list 
      INTEGER, DIMENSION(1:maxw), INTENT(OUT) :: c_list    ! Size of array c_list 
      INTEGER, INTENT(OUT) :: natoms  ! number of atoms in connectivity list
!
      INTEGER    :: i
!
      natoms = 0
      c_list = 0   ! clear connectivity list
!
!
      IF ( ALLOCATED(at_conn) ) THEN
         i = jatom 
        IF ( ASSOCIATED(at_conn(i)%liste) ) THEN      ! If neighborhood was created
          hood_temp => at_conn(i)%liste               ! point to the first neighborhood
!         hood_head => at_conn(i)%liste
          DO WHILE ( ASSOCIATED(hood_temp) )          ! While there are further neighborhood
             temp => hood_temp%nachbar                    ! point to the first neighbor within current neigborhood
             IF ( hood_temp%neigh_type == ino  .OR. &
                  hood_temp%conn_name  == c_name   ) THEN   ! This is the right neighborhood
                DO WHILE ( ASSOCIATED(temp) )               ! While there are further neighbors
                    natoms         = natoms + 1
                    c_list(natoms) = temp%atom_number
                  temp => temp%next                         ! Point to next neighbor
                END DO
                RETURN                                      ! End of connectivity list
             ENDIF
             hood_temp => hood_temp%next_neighborhood     ! Point to next neighborhood
          END DO
        END IF
      ENDIF
!
   END SUBROUTINE get_connectivity_list
!
!
   SUBROUTINE get_connectivity_identity (is1, work_id, work_name, work_name_l)
!-                                                                      
!     Get the identity of a connectivity from central atom and number or name
!+                                                                      
      USE config_mod 
      USE crystal_mod 
      IMPLICIT none 
!
!
      INTEGER,            INTENT(IN)      :: is1        ! central atom type
      INTEGER,            INTENT(INOUT)   :: work_id    ! Connectivity def. no.
      CHARACTER(LEN=256), INTENT(INOUT)   :: work_name  ! Connectivity name
      INTEGER           , INTENT(INOUT)   :: work_name_l! Connectivity name length
!
      INTEGER    :: i
!
!
      IF ( ALLOCATED(def_main) ) THEN
         is_there: IF ( ASSOCIATED(def_main(is1)%def_liste) ) THEN  ! A list of definitions exists
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
                  work_id     = def_temp%valid_id                ! Make sure ID matches
                  work_name   = def_temp%def_name                ! Make sure name matches
                  work_name_l = def_temp%def_name_l              ! Make sure name matches
                  EXIT search
               ENDIF
               def_head => def_temp
               def_temp => def_temp%def_next
            ENDDO search
         ELSE
            ier_num = -109
            ier_typ = ER_APPL
         ENDIF is_there
      ELSE
         ier_num = -110
         ier_typ = ER_APPL
      ENDIF
!
   END SUBROUTINE get_connectivity_identity
!
!
   SUBROUTINE do_show_connectivity ( iatom, idef )
!-                                                                      
!     Shows the connectivity no. idef around atom no iatom
!+                                                                      
      USE param_mod 
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
!
      INTEGER, INTENT(in)        :: iatom
      INTEGER, INTENT(in)        :: idef
!
      INTEGER, PARAMETER         :: maxw = 2000
      INTEGER                    :: is1
      INTEGER, DIMENSION(1:maxw) :: c_list
      INTEGER                    :: natoms
!                                                                       
      INTEGER                    :: i
!
      is1 = cr_iscat(iatom)
      CALL get_connectivity_list (iatom, is1, idef, maxw, c_list, natoms )
!
      WRITE(output_io,1000) iatom, is1
      IF ( natoms > 0 ) THEN
        WRITE(output_io, 1100) (c_list(i),i=1,natoms)
      ELSE
         WRITE(output_io, 1200)
      ENDIF
      IF ( natoms <= MAXPAR_RES ) THEN
         res_para(0) = FLOAT(natoms)
         res_para(1:natoms) = FLOAT(c_list(1:natoms))
      ENDIF
!
1000  FORMAT( ' Connectivity for atom No. ',I6,' of type ',i4)
1100  FORMAT( '      Neighbors are        ',20(i6,2x))
1200  FORMAT( '      Atom has no neighbours')
!
      END SUBROUTINE do_show_connectivity 
END MODULE conn_mod
