MODULE read_internal_mod
!
USE class_internal
USE errlist_mod
!
IMPLICIT NONE
!
!  readstru_size_int  ! Read header, get sizes
!  readstru_internal  ! Read internal structure
!  readcell_internal  ! Read an internal asymmetric unit
!  stru_readheader_internal     ! Read complete internal header
!  struc_read_atoms_internal    ! Read all atoms from internal file
!  struc_read_one_atom_internal ! Read one atom from internal file
!  stru_internal_molecules      ! Read molecules from internal file
!  testfile_internal  ! Test internal file for number of atoms, scattering curves
!
CONTAINS
!*******************************************************************************
   SUBROUTINE readstru_size_int(strucfile, natoms, & 
              nscat, n_mole, n_type, n_atom)
!
!  Reads the header section of a structure from an internal cystal. 
!
   USE discus_allocate_appl_mod
   IMPLICIT NONE
!
!
   CHARACTER (LEN=*), INTENT(IN) :: strucfile
!
   INTEGER       , INTENT(INOUT) :: natoms
   INTEGER       , INTENT(INOUT) :: nscat
   INTEGER       , INTENT(INOUT) :: n_mole
   INTEGER       , INTENT(INOUT) :: n_type
   INTEGER       , INTENT(INOUT) :: n_atom
!
integer ier
!
   NULLIFY(read_from)
   NULLIFY(read_parent)
   IF(.NOT.ASSOCIATED(store_root)) THEN
      ier_num = -113
      ier_typ = ER_APPL
      RETURN
   ENDIF
   CALL store_find_node(store_root, read_from, strucfile, read_temp, read_parent, ier ) ! Find the proper node
   IF ( ier /= 0 .OR. .NOT.ASSOCIATED(read_temp)) THEN
      ier_num = -113
      ier_typ = ER_APPL
      RETURN
   ENDIF
!
   natoms = read_temp%crystal%get_natoms()
   nscat  = read_temp%crystal%get_nscat ()
   n_mole = read_temp%crystal%get_n_mole()
   n_type = read_temp%crystal%get_n_type()
   n_atom = read_temp%crystal%get_n_atom()
!
   END SUBROUTINE readstru_size_int
!*******************************************************************************
   SUBROUTINE readstru_internal(strucfile) !,NMAX, MAXSCAT, MOLE_MAX_MOLE, &
!              mole_max_TYPE, MOLE_MAX_ATOM)
!
!  Reads a structure from an internal cystal. The old crystal is overwritten
!
   USE discus_allocate_appl_mod
   USE chem_mod
   USE crystal_mod
   USE prop_para_mod
   USE molecule_mod
!   USE class_internal
   IMPLICIT NONE
!
!
   CHARACTER (LEN=*), INTENT(IN) :: strucfile
!  INTEGER          , INTENT(IN) :: NMAX
!  INTEGER          , INTENT(IN) :: MAXSCAT
!  INTEGER          , INTENT(IN) :: MOLE_MAX_MOLE
!  INTEGER          , INTENT(IN) :: MOLE_MAX_TYPE
!  INTEGER          , INTENT(IN) :: MOLE_MAX_ATOM
!
   INTEGER                       :: istatus
   INTEGER                       :: natoms
   INTEGER                       :: nscat
   INTEGER                       :: n_mole
   INTEGER                       :: n_type
   INTEGER                       :: n_atom
   INTEGER                       :: i,j
   INTEGER                       :: iatom
!
   CALL readstru_size_int(strucfile, natoms, & 
                   nscat, n_mole, n_type, n_atom)
   IF ( ier_num /= 0) THEN                        ! Could not find the internal storage file
      RETURN
   ENDIF
!
!  Allocate sufficient space
   IF(natoms > NMAX .or. nscat > MAXSCAT ) THEN
      natoms = MAX(natoms, NMAX)
      nscat  = MAX(nscat , MAXSCAT)
      CALL alloc_crystal (nscat, natoms)
      IF ( ier_num /= 0 ) THEN
         ier_num = -114
         ier_typ = ER_APPL
         ier_msg(1) = 'Standard crystal could not be allocated'
         ier_msg(2) = 'In readstru_internal'
         RETURN
      ENDIF
   ENDIF
   IF ( n_mole > MOLE_MAX_MOLE .or. n_type > MOLE_MAX_TYPE .or. &
        n_atom > MOLE_MAX_MOLE                                ) THEN  ! If molecules were present in internal file
      n_mole = MAX(n_mole, MOLE_MAX_MOLE)
      n_type = MAX(n_type, MOLE_MAX_TYPE)
      n_atom = MAX(n_atom, MOLE_MAX_ATOM)
      CALL alloc_molecule(1, 1,n_mole,n_type,n_atom)
      IF ( ier_num /= 0 ) THEN
         ier_num = -114
         ier_typ = ER_APPL
         ier_msg(1) = 'Standard crystal could not be allocated'
         ier_msg(2) = 'In readstru_internal'
         RETURN
      ENDIF
   ENDIF
!
!  Now copy from crystal to standard
!
   CALL read_temp%crystal%get_header_from_crystal()
   CALL read_temp%crystal%get_atoms_from_crystal()
   CALL read_temp%crystal%get_molecules_from_crystal(mole_max_mole,       &
              mole_max_type, mole_max_atom, mole_num_mole, mole_num_type, &
              mole_num_atom, mole_len, mole_off, mole_type, mole_char,    &
              mole_file, mole_dens, mole_biso, mole_clin, mole_cqua,      &
              mole_fuzzy, mole_cont)
   DO i = 1, mole_num_mole       ! set molecule number for each atom
      DO j = 1, mole_len (i)
         iatom          = mole_cont (mole_off (i) + j)
         cr_prop(iatom) = ibset(cr_prop(iatom),PROP_MOLECULE)
         cr_mole(iatom) = i
      ENDDO
   ENDDO
!
   chem_purge = .FALSE.                          ! No purge, period boundary is OK
!
   END SUBROUTINE readstru_internal
!*******************************************************************************
   SUBROUTINE readcell_internal ( strucfile )
!
!  Reads a unit cell from an internal cystal. The old crystal is overwritten
!
   USE discus_allocate_appl_mod
!   USE class_internal
   USE chem_mod
   USE crystal_mod
   USE molecule_mod
   USE spcgr_apply, ONLY: get_symmetry_matrices, firstcell, symmetry
   USE wyckoff_mod
   IMPLICIT NONE
!
!
   CHARACTER (LEN=*), INTENT(IN) :: strucfile
!
   INTEGER                       :: i,j,k        ! Dummy
   INTEGER                       :: ia         ! Dummy; atoms in internal crystal
   INTEGER                       :: istatus
   INTEGER                       :: natoms
   INTEGER                       :: nscat
   INTEGER                       :: n_mole
   INTEGER                       :: n_type
   INTEGER                       :: n_atom
   INTEGER                       :: new_mole
   INTEGER                       :: new_nmax
   INTEGER                       :: new_nscat
   INTEGER, DIMENSION(3)         :: rd_icc     ! Crystal size on 'cell' command line
   INTEGER                       :: itype      ! type of current atom
   REAL   , DIMENSION(3)         :: posit      ! position of current atom
   INTEGER, DIMENSION(0:3)       :: isurface   ! surface  of current atom
   INTEGER                       :: iprop      ! property of current atom
   REAL   , DIMENSION(5)         :: werte      ! temporary array
!
   CHARACTER (LEN=4)             :: at_name    ! temporary atom name
   REAL                          :: dw1        ! temporary DW factor
   REAL                          :: occ1       ! temporary occupancy factor
!
   INTEGER                       :: i_mole     ! Atoom is part of molecule i_mole
   INTEGER                       :: i_type     ! and this molecule is of type i_type
   INTEGER                       :: i_char     ! and this molecule is of type i_type
   CHARACTER (LEN=200)           :: c_file     ! and this molecule is of type i_type
   REAL                          :: r_fuzzy    ! and this molecule is of type i_type
   REAL                          :: r_dens     ! and this molecule is of type i_type
   REAL                          :: r_biso     ! and this molecule is of type i_type
   REAL                          :: r_clin     ! and this molecule is of type i_type
   REAL                          :: r_cqua     ! and this molecule is of type i_type
!

   INTEGER                              :: temp_num_mole
   INTEGER                              :: temp_num_type
   INTEGER                              :: temp_num_atom
   INTEGER, DIMENSION(:  ), ALLOCATABLE :: temp_len
   INTEGER, DIMENSION(:  ), ALLOCATABLE :: temp_off
   INTEGER, DIMENSION(:  ), ALLOCATABLE :: temp_type
   CHARACTER (LEN=200), DIMENSION(:  ), ALLOCATABLE :: temp_file
   INTEGER, DIMENSION(:  ), ALLOCATABLE :: temp_char
   REAL   , DIMENSION(:  ), ALLOCATABLE :: temp_dens
   REAL   , DIMENSION(:  ), ALLOCATABLE :: temp_biso
   REAL   , DIMENSION(:  ), ALLOCATABLE :: temp_clin
   REAL   , DIMENSION(:  ), ALLOCATABLE :: temp_cqua
   REAL   , DIMENSION(:  ), ALLOCATABLE :: temp_fuzz
   INTEGER, DIMENSION(:  ), ALLOCATABLE :: temp_cont
   INTEGER, DIMENSION(:  ), ALLOCATABLE :: temp_look
   INTEGER, DIMENSION(1,2)              :: iin_mole
!
   LOGICAL                       :: need_alloc ! we need to allocate something
   LOGICAL                       :: new_type   ! Each atom is a new type
!
   rd_icc = cr_icc               ! Save crystal dimensions
!
   new_type = cr_newtype         ! Was defined via the 'cell' or 'lcell' command
!  ALLOCATE(read_temp, STAT = istatus )            ! Allocate a temporary storage
!  IF ( istatus /= 0) THEN
!     ier_num = -114
!     ier_typ = ER_APPL
!     RETURN
!  ENDIF
!   ALLOCATE(read_temp, STAT = istatus )           ! Allocate a temporary storage
!write(*,*) ' STATUS of ALLOC ', istatus
   CALL readstru_size_int(strucfile, natoms, &     ! Get the sizes of the internal crystal
                   nscat, n_mole, n_type, n_atom)
   IF ( ier_num /= 0) THEN
      ier_typ = ER_APPL
      ier_msg(1) = 'Could not get size of stored crystal'
      RETURN
   ENDIF
!  Get header
   CALL read_temp%crystal%get_header_from_crystal() ! Read the header
   cr_icc = rd_icc               ! Restore crystal dimensions
!
   CALL get_symmetry_matrices                       ! Setup symmetry
!
!  Allocate enough space for one unit cell
   need_alloc = .false.
   IF ( NMAX < spc_n*natoms ) THEN
      new_nmax = spc_n*natoms + 1
      need_alloc = .true.
   ELSE
      new_nmax = NMAX
   ENDIF
   IF ( MAXSCAT < nscat                    ) THEN
      new_nscat =     nscat
      need_alloc = .true.
   ELSE
      new_nscat= MAXSCAT
   ENDIF
   IF( new_type .and. natoms > new_nscat ) THEN     ! Each atom is a new scattering type make enough space
      new_nscat = MAX(new_nscat, natoms)
      need_alloc = .true.
   ENDIF
   IF ( need_alloc ) THEN
      call alloc_crystal(new_nscat, new_nmax)
      IF ( ier_num /= 0 ) THEN
         ier_msg(1) = 'Could not allocate space for actual crystal'
         RETURN
      ENDIF
   ENDIF
!  Allocate space for molecules
found: IF ( n_mole > 0 ) THEN      ! FOUND MOLECULES
   need_alloc = .false.
   IF ( MOLE_MAX_MOLE < n_mole*spc_n ) THEN
      new_mole = n_mole*spc_n
      need_alloc = .true.
   ELSE
      new_mole = MOLE_MAX_MOLE
   ENDIF
   IF ( n_type > MOLE_MAX_TYPE ) THEN
      need_alloc = .true.
   ENDIF
!
   IF ( MOLE_MAX_ATOM < n_atom*spc_n ) THEN
      new_nmax = n_atom*spc_n
      need_alloc = .true.
   ELSE
      new_nmax = MOLE_MAX_ATOM
   ENDIF
   IF ( need_alloc ) THEN
      call alloc_molecule(1,1,new_mole, n_type, new_nmax)
      IF ( ier_num /= 0 ) THEN
         ier_msg(1) = 'Could not allocate space for molecules'
         WRITE(ier_msg(2), '(i15,1x,i15)')  new_mole, n_type
         WRITE(ier_msg(3), '(i15       )')  new_nmax
         RETURN
      ENDIF
   ENDIF
!
!  Allocate temporary space for molecule info
! 
      ALLOCATE(temp_len (0:n_mole))
      ALLOCATE(temp_off (0:n_mole))
      ALLOCATE(temp_type(0:n_mole))
      ALLOCATE(temp_file(0:n_mole))
      ALLOCATE(temp_char(0:n_mole))
      ALLOCATE(temp_dens(0:n_mole))
      ALLOCATE(temp_biso(0:n_type))
      ALLOCATE(temp_clin(0:n_type))
      ALLOCATE(temp_cqua(0:n_type))
      ALLOCATE(temp_fuzz(0:n_mole))
      ALLOCATE(temp_cont(0:natoms))
      ALLOCATE(temp_look(0:natoms))
      temp_len (0:n_mole) = 0
      temp_off (0:n_mole) = 0
      temp_type(0:n_mole) = 0
      temp_file(0:n_mole) = ' '
      temp_char(0:n_mole) = 0
      temp_dens(0:n_mole) = 0.0
      temp_biso(0:n_type) = 0.0
      temp_clin(0:n_type) = 0.0
      temp_cqua(0:n_type) = 0.0
      temp_fuzz(0:n_mole) = 0.0
      temp_cont(0:natoms) = 0
      temp_look(0:natoms) = 0
!
!
!  Now copy from internal crystal to local variables
!
      CALL read_temp%crystal%get_molecules_from_crystal(n_mole,       &
              n_type, natoms, temp_num_mole, temp_num_type, &
              temp_num_atom, temp_len, temp_off, temp_type, temp_char,    &
              temp_file, temp_dens, temp_biso, temp_clin, temp_cqua,    &
              temp_fuzz, temp_cont)
!
!  Build lookup table for original molecules
!
   DO i=1, temp_num_mole
      DO j=1, temp_len(i)
         ia = temp_cont(temp_off(i)+j)
         temp_look(ia) = i            !atom(iatom) is in molecule i
      ENDDO
   ENDDO
   ENDIF found
!
!  Main loop over all atoms in the asymmetric unit   
!
   mole_l_on = .false.
   cr_natoms = 0
   cr_nscat  = 0
   mole_len  = 0
   mole_off  = 0
   mole_type = 0
   mole_cont = 0
   i_mole    = 0
   mole_num_atom = 0
   mole_num_mole = 0
!
main: do ia = 1, natoms
      werte    = 0.0
      cr_natoms = cr_natoms + 1
      CALL read_temp%crystal%get_cryst_atom ( ia, itype, posit, iprop, isurface, iin_mole)
      CALL read_temp%crystal%get_cryst_scat ( ia, itype, at_name , dw1, occ1  )
mole_exist: if(n_mole > 0) THEN
      CALL read_temp%crystal%get_cryst_mole ( ia, i_mole, i_type,  &
                 i_char, c_file, r_fuzzy, r_dens, r_biso, r_clin, r_cqua)
in_mole: IF ( temp_look(ia) > 0 ) THEN            ! This atom belongs to a molecule
            IF ( .not. mole_l_on .OR. temp_look(ia)> mole_num_curr ) THEN  ! Right now we are not in a molecule
               mole_l_on    = .true.              ! Turn molecule on
               mole_l_first = .true.              ! This is the first atom in the molecule
               mole_gene_n  = 0                   ! No molecule generators
               mole_symm_n  = 0                   ! no molecule symmetry
               i_mole = mole_num_mole + 1         ! Need to work on a new (=next) molecule
               mole_off  (i_mole) = mole_off(mole_num_mole)+mole_len(mole_num_mole) ! Set offset
               mole_num_mole      = MAX(mole_num_mole, i_mole) ! Adjust no of molecules
               mole_num_type      = MAX(mole_num_type, i_type) ! Adjust no of molecule types
               mole_type (i_mole) = i_type        ! Set current molecule type
               mole_char (i_mole) = i_char        ! Set current molecule character
               mole_file (i_mole) = c_file        ! Set current molecule file
               mole_fuzzy(i_mole) = r_fuzzy       ! Set current molecule Fuzzy distance
               mole_dens (i_mole) = r_dens        ! Set current molecule density
               mole_biso (mole_type(i_mole)) = r_biso ! Set current molecule b-value
               mole_clin (mole_type(i_mole)) = r_clin ! Set current molecule b-value
               mole_cqua (mole_type(i_mole)) = r_cqua ! Set current molecule b-value
               mole_num_curr      = i_mole        ! Set molecule no we are working on
            ENDIF
      ELSE  in_mole
         mole_l_on = .false.                      ! Next atom is not inside a molecule
      ENDIF in_mole
      ENDIF mole_exist
scat_dw: IF ( new_type ) THEN                     ! force each atom to be a new type 
         cr_nscat = cr_nscat + 1
         cr_at_lis(cr_nscat) = at_name
         cr_dw    (cr_nscat) = dw1
         cr_occ   (cr_nscat) = occ1
         itype               = cr_nscat
      ELSE scat_dw                                ! check for previous atom types
         itype = -1
do_scat_dw: DO k = 1,cr_nscat
            IF ( at_name == cr_at_lis(k) .AND.  &
                 dw1     == cr_dw    (k) .AND.  &
                 occ1    == cr_occ   (k)   ) THEN
               itype = k
               EXIT do_scat_dw
            ENDIF
         ENDDO do_scat_dw
         IF ( itype == -1 ) THEN
            cr_nscat = cr_nscat + 1               ! Atom type does not exist, create a new one
            cr_at_lis(cr_nscat) = at_name
            cr_dw    (cr_nscat) = dw1
            cr_occ   (cr_nscat) = occ1
            itype               = cr_nscat
         ENDIF
      ENDIF scat_dw                               ! End check for atom type
!
      werte(1)  = posit(1)                        ! Just to be compatible to firstcell...
      werte(2)  = posit(2)
      werte(3)  = posit(3)
      IF ( .not. mole_l_on ) THEN
         CALL firstcell ( werte, 5)
      ENDIF
      DO j=1,3
         cr_pos(j,cr_natoms) = werte(j)           ! strore in actual crystal
      ENDDO
      cr_iscat(cr_natoms) = itype                 ! set the atom type
!     cr_mole (cr_natoms) = i_mole                ! set the molecule number
      cr_mole (cr_natoms) = 0                     ! set the molecule number
      cr_surf(:,cr_natoms) = isurface               ! set the property flag
      cr_prop (cr_natoms) = iprop                 ! set the property flag
      CALL symmetry
   ENDDO main
!
   CALL no_error 
!                                                                       
!  move first unit cell into lower left corner of crystal          
!                                                                       
   DO i = 1, cr_natoms 
      DO j = 1, 3 
         cr_pos (j, i) = cr_pos (j, i) - int ( (cr_icc (j) ) / 2) 
         cr_dim (j, 1) = amin1 (cr_dim (j, 1), cr_pos (j, i) ) 
         cr_dim (j, 2) = amax1 (cr_dim (j, 2), cr_pos (j, i) ) 
      ENDDO 
   ENDDO 
!
!  DEALLOCATE(read_temp, STAT = istatus )         ! Deallocate a temporary storage
IF(n_mole > 0) THEN
!
!  Deallocate temporary space for molecule info
! 
   DEALLOCATE(temp_len )
   DEALLOCATE(temp_off )
   DEALLOCATE(temp_type)
   DEALLOCATE(temp_file)
   DEALLOCATE(temp_char)
   DEALLOCATE(temp_dens)
   DEALLOCATE(temp_biso)
   DEALLOCATE(temp_clin)
   DEALLOCATE(temp_cqua)
   DEALLOCATE(temp_fuzz)
   DEALLOCATE(temp_cont)
   DEALLOCATE(temp_look)
ENDIF
!
   chem_purge = .FALSE.                          ! No purge, period boundary is OK
!
   END SUBROUTINE readcell_internal
!*******************************************************************************
   SUBROUTINE stru_readheader_internal (rd_strucfile, rd_MAXSCAT, rd_cr_name,   &
            rd_cr_spcgr, rd_cr_at_lis, rd_cr_nscat, rd_cr_dw, rd_cr_occ, rd_cr_a0, rd_cr_win, &
            rd_sav_ncell, rd_sav_r_ncell, rd_sav_ncatoms, rd_spcgr_ianz, rd_spcgr_para,   &
            rd_GEN_ADD_MAX, rd_gen_add_n, rd_gen_add_power, rd_gen_add,                 &
            rd_SYM_ADD_MAX, rd_sym_add_n, rd_sym_add_power, rd_sym_add )
!
!  Copies the header info from the internal crystal 'rd_strucfile' into local variables
!
   IMPLICIT NONE
!
   CHARACTER (LEN=  * )                         , INTENT(IN   ) :: rd_strucfile 
   INTEGER                                      , INTENT(INOUT) :: rd_MAXSCAT 
!
   CHARACTER (LEN=  80)                         , INTENT(INOUT) :: rd_cr_name 
   CHARACTER (LEN=  16)                         , INTENT(INOUT) :: rd_cr_spcgr 
   REAL                , DIMENSION(3)           , INTENT(INOUT) :: rd_cr_a0
   REAL                , DIMENSION(3)           , INTENT(INOUT) :: rd_cr_win
   INTEGER                                      , INTENT(INOUT) :: rd_cr_nscat 
   REAL                , DIMENSION(0:rd_MAXSCAT), INTENT(INOUT) :: rd_cr_dw     ! (0:MAXSCAT) 
   REAL                , DIMENSION(0:rd_MAXSCAT), INTENT(INOUT) :: rd_cr_occ    ! (0:MAXSCAT) 
   CHARACTER (LEN=   4), DIMENSION(0:rd_MAXSCAT), INTENT(INOUT) :: rd_cr_at_lis ! (0:MAXSCAT) 
   INTEGER             , DIMENSION(3)           , INTENT(INOUT) :: rd_sav_ncell ! (3) 
   LOGICAL                                      , INTENT(INOUT) :: rd_sav_r_ncell 
   INTEGER                                      , INTENT(INOUT) :: rd_sav_ncatoms 
   INTEGER                                      , INTENT(INOUT) :: rd_spcgr_ianz 
   INTEGER                                      , INTENT(INOUT) :: rd_spcgr_para 
!
   INTEGER             ::  rd_GEN_ADD_MAX
   INTEGER             ::  rd_gen_add_n
   INTEGER             ::  rd_gen_add_power(rd_GEN_ADD_MAX)
   REAL                ::  rd_gen_add(4,4,0:rd_GEN_ADD_MAX)
!
   INTEGER             ::  rd_SYM_ADD_MAX
   INTEGER             ::  rd_sym_add_n
   INTEGER             ::  rd_sym_add_power(rd_SYM_ADD_MAX)
   REAL                ::  rd_sym_add(4,4,0:rd_SYM_ADD_MAX)
!
   INTEGER             :: istatus
!
   NULLIFY(read_from)
   NULLIFY(read_parent)
   CALL store_find_node(store_root, read_from, rd_strucfile, read_temp, read_parent, ier_num ) ! Find the proper node
   IF ( ier_num /= 0) THEN
      ier_typ = ER_APPL
      RETURN
   ENDIF
!
   CALL read_temp%crystal%get_header_to_local (rd_MAXSCAT, rd_cr_name,      &
            rd_cr_spcgr, rd_cr_at_lis, rd_cr_nscat, rd_cr_dw, rd_cr_occ, rd_cr_a0, rd_cr_win, &
            rd_sav_ncell, rd_sav_r_ncell, rd_sav_ncatoms, rd_spcgr_ianz, rd_spcgr_para, &
            rd_GEN_ADD_MAX, rd_gen_add_n, rd_gen_add_power, rd_gen_add,                 &
            rd_SYM_ADD_MAX, rd_sym_add_n, rd_sym_add_power, rd_sym_add )
!
   END SUBROUTINE stru_readheader_internal
!*******************************************************************************
   SUBROUTINE struc_read_atoms_internal(strucfile, RD_NMAX, &
              rd_cr_natoms, rd_cr_pos, rd_cr_iscat, rd_cr_prop, &
              rd_cr_surf, rd_cr_mole )
!
!  This subroutine adds all atoms from the internal storage to the local
!  crystal rd_cr_pos. The number of atoms cr_natoms is incremented accordingly
!
   USE discus_allocate_appl_mod
   IMPLICIT NONE
!
   CHARACTER (LEN=  * )                         , INTENT(IN   ) :: strucfile 
   INTEGER                                      , INTENT(IN   ) :: rd_NMAX 
   INTEGER                                      , INTENT(INOUT) :: rd_cr_natoms
   REAL                , DIMENSION(3,1:RD_NMAX) , INTENT(INOUT) :: rd_cr_pos
   INTEGER             , DIMENSION(  1:RD_NMAX) , INTENT(INOUT) :: rd_cr_iscat
   INTEGER             , DIMENSION(0:3,1:RD_NMAX),INTENT(INOUT) :: rd_cr_surf
   INTEGER             , DIMENSION(  1:RD_NMAX) , INTENT(INOUT) :: rd_cr_prop
   INTEGER             , DIMENSION(  1:RD_NMAX) , INTENT(INOUT) :: rd_cr_mole
!
   INTEGER                       :: i
   INTEGER                       :: istatus
   INTEGER                       :: natoms
   INTEGER                       :: nscat
   INTEGER                       :: n_mole
   INTEGER                       :: n_type
   INTEGER                       :: n_atom
   INTEGER                       :: itype
   REAL   , DIMENSION(3)         :: posit
   INTEGER, DIMENSION(0:3)       :: isurface
   INTEGER, DIMENSION(1:2)       :: iin_mole
   INTEGER                       :: iprop
!
   NULLIFY(read_from)
   NULLIFY(read_parent)
   CALL store_find_node(store_root, read_from, strucfile, read_temp, read_parent, ier_num ) ! Find the proper node
   IF ( ier_num /= 0) THEN
      ier_typ = ER_APPL
      RETURN
   ENDIF
   CALL readstru_size_int(strucfile, natoms, &           ! Get the number of atoms
              nscat, n_mole, n_type, n_atom)
   DO i=1,natoms
      CALL read_temp%crystal%get_cryst_atom(i, itype, posit, iprop, isurface, iin_mole)
      rd_cr_natoms = rd_cr_natoms + 1
      rd_cr_pos(1,rd_cr_natoms) = posit(1)
      rd_cr_pos(2,rd_cr_natoms) = posit(2)
      rd_cr_pos(3,rd_cr_natoms) = posit(3)
      rd_cr_iscat(rd_cr_natoms) = itype
      rd_cr_prop (rd_cr_natoms) = iprop 
      rd_cr_mole (rd_cr_natoms) = iin_mole(1) 
   ENDDO
!
   END SUBROUTINE struc_read_atoms_internal
!*******************************************************************************
   SUBROUTINE struc_read_one_atom_internal(strucfile, iatom,  &
              rd_cr_pos, rd_cr_iscat, rd_cr_prop, rd_cr_surf, &
              rd_cr_mole, rd_cr_moleatom )
!
!  This subroutine reads just atom "iatom" from the internal storage.
!  The atom is copied into internal storage
!
   IMPLICIT NONE
!
   CHARACTER (LEN=  * )               , INTENT(IN   ) :: strucfile 
   INTEGER                            , INTENT(IN   ) :: iatom
   REAL                , DIMENSION(3) , INTENT(INOUT) :: rd_cr_pos
   INTEGER                            , INTENT(INOUT) :: rd_cr_iscat
   INTEGER                            , INTENT(INOUT) :: rd_cr_prop
   INTEGER             , DIMENSION(0:3),INTENT(INOUT) :: rd_cr_surf
   INTEGER                             ,INTENT(INOUT) :: rd_cr_mole
   INTEGER                             ,INTENT(INOUT) :: rd_cr_moleatom
!
   INTEGER                       :: i
   INTEGER                       :: istatus
   INTEGER                       :: natoms
   INTEGER                       :: nscat
   INTEGER                       :: n_mole
   INTEGER                       :: n_type
   INTEGER                       :: n_atom
   INTEGER, DIMENSION(1:2)       :: iin_mole
!
   NULLIFY(read_from)
   NULLIFY(read_parent)
   CALL store_find_node(store_root, read_from, strucfile, read_temp, read_parent, ier_num ) ! Find the proper node
   IF ( ier_num /= 0) THEN
      ier_typ = ER_APPL
      RETURN
   ENDIF
   CALL readstru_size_int(strucfile, natoms, &           ! Get the number of atoms
              nscat, n_mole, n_type, n_atom)
   IF ( iatom <= natoms ) THEN
      i = iatom
      CALL read_temp%crystal%get_cryst_atom(i, rd_cr_iscat, rd_cr_pos, &
           rd_cr_prop, rd_cr_surf, iin_mole)
      rd_cr_mole     = iin_mole(1)
      rd_cr_moleatom = iin_mole(2)
   ELSE
      ier_num = -105
      ier_typ = ER_APPL
   ENDIF
!
   END SUBROUTINE struc_read_one_atom_internal
!*******************************************************************************
   SUBROUTINE stru_internal_molecules(strucfile, MOLE_MAX_MOLE,           &
              MOLE_MAX_TYPE, MOLE_MAX_ATOM, mole_num_mole, mole_num_type, &
              mole_num_atom, mole_len, mole_off, mole_type, mole_char,    &
              mole_file, mole_dens, mole_biso, mole_clin, mole_cqua,      &
              mole_fuzzy, mole_cont)
!
!  Read all molecules from internal storage into local variables.
!
!
!  Reads a structure from an internal cystal. The old crystal is overwritten
!
   USE discus_allocate_appl_mod
   IMPLICIT NONE
!
!
   CHARACTER (LEN=*), INTENT(IN) :: strucfile
   INTEGER,                                       INTENT(IN ) :: mole_max_mole
   INTEGER,                                       INTENT(IN ) :: mole_max_type
   INTEGER,                                       INTENT(IN ) :: mole_max_atom
   INTEGER,                                       INTENT(OUT) :: mole_num_mole
   INTEGER,                                       INTENT(OUT) :: mole_num_type
   INTEGER,                                       INTENT(OUT) :: mole_num_atom
   INTEGER,           DIMENSION(0:MOLE_MAX_MOLE), INTENT(OUT) :: mole_len
   INTEGER,           DIMENSION(0:MOLE_MAX_MOLE), INTENT(OUT) :: mole_off
   INTEGER,           DIMENSION(0:MOLE_MAX_MOLE), INTENT(OUT) :: mole_type
   INTEGER,           DIMENSION(0:MOLE_MAX_MOLE), INTENT(OUT) :: mole_char
   CHARACTER (LEN=*), DIMENSION(0:MOLE_MAX_MOLE), INTENT(OUT) :: mole_file
   REAL   ,           DIMENSION(0:MOLE_MAX_MOLE), INTENT(OUT) :: mole_dens
   REAL   ,           DIMENSION(0:MOLE_MAX_TYPE), INTENT(OUT) :: mole_biso
   REAL   ,           DIMENSION(0:MOLE_MAX_TYPE), INTENT(OUT) :: mole_clin
   REAL   ,           DIMENSION(0:MOLE_MAX_TYPE), INTENT(OUT) :: mole_cqua
   REAL   ,           DIMENSION(0:MOLE_MAX_MOLE), INTENT(OUT) :: mole_fuzzy
   INTEGER,           DIMENSION(0:MOLE_MAX_ATOM), INTENT(OUT) :: mole_cont
!
   INTEGER                       :: natoms
   INTEGER                       :: nscat
   INTEGER                       :: istatus
!
!  ALLOCATE(read_temp, STAT = istatus )                  ! Allocate a temporary storage
!  IF ( istatus /= 0) THEN
!     ier_num = -114
!     ier_typ = ER_APPL
!     RETURN
!  ENDIF
   CALL readstru_size_int(strucfile, natoms, & 
                   nscat, mole_num_mole, mole_num_type, mole_num_atom)
   IF ( ier_num /= 0) THEN                        ! Could not find the internal storage file
      RETURN
   ENDIF
!
!
!  Now copy from crystal to local variables
!
   CALL read_temp%crystal%get_molecules_from_crystal(MOLE_MAX_MOLE,       &
              MOLE_MAX_TYPE, MOLE_MAX_ATOM, mole_num_mole, mole_num_type, &
              mole_num_atom, mole_len, mole_off, mole_type, mole_char,    &
              mole_file, mole_dens, mole_biso, mole_clin, mole_cqua,      &
              mole_fuzzy, mole_cont)
!
!  DEALLOCATE(read_temp, STAT = istatus )        ! Deallocate a temporary storage
!
   END SUBROUTINE stru_internal_molecules
!*******************************************************************************
   SUBROUTINE testfile_internal(strucfile , natoms, & 
              nscat, n_mole, n_type, n_atom)
!
!  Reads the header section of an internal file and determines 
!  the required sizes
!
   IMPLICIT NONE
!
!
   CHARACTER (LEN=*), INTENT(IN) :: strucfile
!
   INTEGER                       :: istatus
   INTEGER, INTENT(INOUT)        :: natoms
   INTEGER, INTENT(INOUT)        :: nscat
   INTEGER, INTENT(INOUT)        :: n_mole
   INTEGER, INTENT(INOUT)        :: n_type
   INTEGER, INTENT(INOUT)        :: n_atom
!
!  ALLOCATE(read_temp, STAT = istatus )                  ! Allocate a temporary storage
!  IF ( istatus /= 0) THEN
!     ier_num = -114
!     ier_typ = ER_APPL
!     RETURN
!  ENDIF
   CALL readstru_size_int(strucfile, natoms, & 
                   nscat, n_mole, n_type, n_atom)
!
!  DEALLOCATE(read_temp, STAT = istatus )        ! Deallocate a temporary storage
!
   END SUBROUTINE testfile_internal
!*******************************************************************************

END MODULE read_internal_mod
