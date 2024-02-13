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
!
!*******************************************************************************
!
SUBROUTINE readstru_size_int(strucfile, dim_natoms, dim_ncatoms,     &
           dim_nscat, dim_nanis, dim_n_mole, dim_n_type, dim_n_atom, &
           natoms, ncatoms, nscat, nanis, n_mole, n_type, n_atom)
!
!  Read the array sizes and number of atoms, scattering types etc from internal file
!
USE discus_allocate_appl_mod
IMPLICIT NONE
!
!
CHARACTER (LEN=*), INTENT(IN) :: strucfile
!
INTEGER       , INTENT(  OUT) :: dim_natoms     ! Total number of atoms 
INTEGER       , intent(  out) :: dim_ncatoms    ! Number of atoms per uni cell
INTEGER       , INTENT(  OUT) :: dim_nscat      ! Number of chemical species
INTEGER       , INTENT(  OUT) :: dim_nanis      ! Number of ADPs
INTEGER       , INTENT(  OUT) :: dim_n_mole     ! Number of molecules
INTEGER       , INTENT(  OUT) :: dim_n_type     ! Number of molecule types
INTEGER       , INTENT(  OUT) :: dim_n_atom     ! Total number of atoms in molecules
INTEGER       , INTENT(INOUT) :: natoms     ! Total number of atoms 
INTEGER       , intent(  out) :: ncatoms    ! Number of atoms per uni cell
INTEGER       , INTENT(INOUT) :: nscat      ! Number of chemical species
INTEGER       , INTENT(INOUT) :: nanis      ! Number of ADPs
INTEGER       , INTENT(INOUT) :: n_mole     ! Number of molecules
INTEGER       , INTENT(INOUT) :: n_type     ! Number of molecule types
INTEGER       , INTENT(INOUT) :: n_atom     ! Total number of atoms in molecules
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
call read_temp%crystal%get_dimensions(dim_natoms, dim_nscat, dim_ncatoms, dim_nanis, &
                                dim_n_mole, dim_n_type, dim_n_atom)
!
natoms  = read_temp%crystal%get_natoms()
ncatoms = read_temp%crystal%get_ncatoms()
nscat   = read_temp%crystal%get_nscat ()
nanis   = read_temp%crystal%get_nanis ()
n_mole  = read_temp%crystal%get_n_mole()
n_type  = read_temp%crystal%get_n_type()
n_atom  = read_temp%crystal%get_n_atom()
!
END SUBROUTINE readstru_size_int
!
!*******************************************************************************
!
SUBROUTINE readstru_internal(MAXMASK, strucfile, uni_mask) !,NMAX, MAXSCAT, MOLE_MAX_MOLE, &
!              mole_max_TYPE, MOLE_MAX_ATOM)
!
!  Reads a structure from an internal cystal. The old crystal is overwritten
!
USE discus_allocate_appl_mod
USE chem_mod
USE crystal_mod
USE prop_para_mod
USE molecule_mod
use reduce_atoms_mod
use reduce_atoms_mod
!   USE class_internal
!
IMPLICIT NONE
!
integer          , intent(in) :: MAXMASK
CHARACTER (LEN=*), INTENT(IN) :: strucfile
logical,dimension(0:MAXMASK), intent(in) :: uni_mask   ! Unique atom type mask
!
!
integer                       :: dim_natoms   ! Number of atoms in the structure
integer                       :: dim_ncatoms  ! Number of atoms per unit call in the structure
integer                       :: dim_nscat    ! Number of different atom types
integer                       :: dim_nanis    ! Number of ADPs
integer                       :: dim_n_mole
integer                       :: dim_n_type
integer                       :: dim_n_atom
INTEGER                       :: natoms   ! Number of atoms in the structure
INTEGER                       :: ncatoms  ! Number of atoms per unit call in the structure
INTEGER                       :: nscat    ! Number of different atom types
INTEGER                       :: nanis    ! Number of ADPs
INTEGER                       :: n_mole
INTEGER                       :: n_type
INTEGER                       :: n_atom
INTEGER                       :: i,j
INTEGER                       :: iatom
!
CALL readstru_size_int(strucfile, dim_natoms, dim_ncatoms, dim_nscat, dim_nanis, &
                                  dim_n_mole, dim_n_type, dim_n_atom,            &
                                  natoms, ncatoms, nscat, nanis,                 &
                                  n_mole, n_type, n_atom)
IF ( ier_num /= 0) THEN                        ! Could not find the internal storage file
   RETURN
ENDIF
!write(*,*) ' READ ', dim_natoms, dim_nscat, dim_ncatoms, dim_nanis,&
!dim_n_mole, dim_n_type, dim_n_atom
!write(*,*) ' Read ', natoms, nscat, ncatoms, nanis,&
!n_mole, n_type, n_atom
!write(*,*) ' read ', NMAX, MAXSCAT
call alloc_crystal_nmax(dim_natoms)
call alloc_crystal_scat(dim_nscat)
call alloc_unitcell(dim_ncatoms)
call alloc_anis(dim_nanis)
call alloc_molecule(1, 1,dim_n_mole, dim_n_type, dim_n_atom)
cr_natoms = natoms
cr_ncatoms= ncatoms
cr_nscat  = nscat
cr_nanis  = nanis
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
!
chem_purge = .FALSE.                          ! No purge, period boundary is OK
!
!enddo
if(uni_mask(0)) then            ! User requested reduction
   call reduce_atoms(MAXMASK, uni_mask)
endif
!
END SUBROUTINE readstru_internal
!
!*******************************************************************************
!
SUBROUTINE readcell_internal (MAXMASK, strucfile, uni_mask)
!-
!  Reads a unit cell from an internal cystal. The old crystal is overwritten
!+
!
USE discus_allocate_appl_mod
use atom_line_mod
USE chem_mod
USE cryst_class
USE crystal_mod
USE molecule_mod
USE spcgr_apply, ONLY: get_symmetry_matrices, firstcell, symmetry
USE wyckoff_mod
USE lib_errlist_func
USE precision_mod
!
IMPLICIT NONE
!
integer          , intent(in) :: MAXMASK
CHARACTER (LEN=*), INTENT(IN) :: strucfile
logical, dimension(0:MAXMASK) :: uni_mask
!
INTEGER                       :: i,j,k        ! Dummy
INTEGER                       :: ia         ! Dummy; atoms in internal crystal
integer                       :: dim_natoms   ! Number of atoms in the structure
integer                       :: dim_ncatoms  ! Number of atoms per unit call in the structure
integer                       :: dim_nscat    ! Number of different atom types
integer                       :: dim_nanis    ! Number of ADPs
integer                       :: dim_n_mole
integer                       :: dim_n_type
integer                       :: dim_n_atom
INTEGER                       :: natoms
INTEGER                       :: ncatoms
INTEGER                       :: nscat
INTEGER                       :: nanis
INTEGER                       :: n_mole
INTEGER                       :: n_type
INTEGER                       :: n_atom
INTEGER                       :: new_mole
INTEGER                       :: new_nmax
INTEGER                       :: new_nscat
INTEGER, DIMENSION(3)         :: rd_icc     ! Crystal size on 'cell' command line
INTEGER           , dimension(3)         :: itype      ! type of current atom
REAL(kind=PREC_DP), dimension(3)         :: posit      ! position of current atom
INTEGER, DIMENSION(0:3)       :: isurface   ! surface  of current atom
REAL(kind=PREC_DP), DIMENSION(0:3)       :: magn_mom   ! Magnetic moment
INTEGER                       :: iprop      ! property of current atom
REAL(KIND=PREC_DP), DIMENSION(5)         :: werte      ! temporary array
!
CHARACTER (LEN=4)             :: at_name    ! temporary atom name
CHARACTER (LEN=4)             :: nw_name    ! temporary atom name
REAL(kind=PREC_DP)                          :: dw1        ! temporary DW factor
REAL(kind=PREC_DP)                          :: occ1       ! temporary occupancy factor
!
INTEGER                       :: i_mole     ! Atoom is part of molecule i_mole
INTEGER                       :: i_type     ! and this molecule is of type i_type
INTEGER                       :: i_char     ! and this molecule is of type i_type
CHARACTER (LEN=200)           :: c_file     ! and this molecule is of type i_type
REAL(kind=PREC_DP)                          :: r_fuzzy    ! and this molecule is of type i_type
REAL(kind=PREC_DP)                          :: r_dens     ! and this molecule is of type i_type
REAL(kind=PREC_DP)                          :: r_biso     ! and this molecule is of type i_type
REAL(kind=PREC_DP)                          :: r_clin     ! and this molecule is of type i_type
REAL(kind=PREC_DP)                          :: r_cqua     ! and this molecule is of type i_type
!

INTEGER                              :: temp_num_mole
INTEGER                              :: temp_num_type
INTEGER                              :: temp_num_atom
INTEGER, DIMENSION(:  ), ALLOCATABLE :: temp_len
INTEGER, DIMENSION(:  ), ALLOCATABLE :: temp_off
INTEGER, DIMENSION(:  ), ALLOCATABLE :: temp_type
CHARACTER (LEN=200), DIMENSION(:  ), ALLOCATABLE :: temp_file
INTEGER, DIMENSION(:  ), ALLOCATABLE :: temp_char
REAL(kind=PREC_DP)   , DIMENSION(:  ), ALLOCATABLE :: temp_dens
REAL(kind=PREC_DP)   , DIMENSION(:  ), ALLOCATABLE :: temp_biso
REAL(kind=PREC_DP)   , DIMENSION(:  ), ALLOCATABLE :: temp_clin
REAL(kind=PREC_DP)   , DIMENSION(:  ), ALLOCATABLE :: temp_cqua
REAL(kind=PREC_DP)   , DIMENSION(:  ), ALLOCATABLE :: temp_fuzz
INTEGER, DIMENSION(:  ), ALLOCATABLE :: temp_cont
INTEGER, DIMENSION(:  ), ALLOCATABLE :: temp_look
INTEGER, DIMENSION(:  ), ALLOCATABLE :: temp_upd
INTEGER, DIMENSION(1,2)              :: iin_mole
!
LOGICAL                       :: need_alloc ! we need to allocate something
LOGICAL                       :: new_type   ! Each atom is a new type
!
!
rd_icc = cr_icc               ! Save crystal dimensions
!
new_type = cr_newtype         ! Was defined via the 'cell' or 'lcell' command
!
CALL readstru_size_int(strucfile, dim_natoms, dim_ncatoms, dim_nscat, dim_nanis, &
                                  dim_n_mole, dim_n_type, dim_n_atom,            &
                                  natoms, ncatoms, nscat, nanis,                 &
                                  n_mole, n_type, n_atom)
!CALL readstru_size_int(strucfile, natoms, ncatoms, &     ! Get the sizes of the internal crystal
!                   nscat, nanis, n_mole, n_type, n_atom)
IF ( ier_num /= 0) THEN
   ier_typ = ER_APPL
   ier_msg(1) = 'Could not get size of stored crystal'
   RETURN
ENDIF
!
!  Get number of symmetry operations for propper allocations
spc_n = read_temp%crystal%get_spc_n()
!
dim_n_mole, dim_n_type, dim_n_atom
n_mole, n_type, n_atom
!
!  Allocate enough space for one unit cell
new_nmax  = spc_n*natoms + 1
new_nscat =     nscat
call alloc_crystal_scat(new_nscat)
call alloc_crystal_nmax(new_nmax)
call alloc_unitcell(ncatoms)
call alloc_anis(nanis)

!need_alloc = .false.
!IF ( NMAX < spc_n*natoms ) THEN
!   new_nmax = spc_n*natoms + 1
!   need_alloc = .true.
!ELSE
!   new_nmax = NMAX
!ENDIF
!IF ( MAXSCAT < nscat                    ) THEN
!   new_nscat =     nscat
!   need_alloc = .true.
!ELSE
!   new_nscat= MAXSCAT
!ENDIF
!IF( new_type .and. natoms > new_nscat ) THEN     ! Each atom is a new scattering type make enough space
!   new_nscat = MAX(new_nscat, natoms)
!   need_alloc = .true.
!ENDIF
!IF ( need_alloc ) THEN
!   call alloc_crystal_scat(new_nscat)
!   call alloc_crystal_nmax(new_nmax)
!   IF ( ier_num /= 0 ) THEN
!      ier_msg(1) = 'Could not allocate space for actual crystal'
!      RETURN
!   ENDIF
!ENDIF
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
   ALLOCATE(temp_upd (0:n_mole))
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
   temp_upd (0:n_mole) = 0
!
!  Get header
   CALL read_temp%crystal%get_header_from_crystal() ! Read the header
   cr_icc = rd_icc               ! Restore crystal dimensions
!
   CALL get_symmetry_matrices                       ! Setup symmetry
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
      temp_upd(i) = i                 ! Initially each molecule is at the corect number
      DO j=1, temp_len(i)
         ia = temp_cont(temp_off(i)+j)
         temp_look(ia) = i            !atom(iatom) is in molecule i
!write(*,*) ' UPD, LOOK ', temp_upd(i), temp_look(ia), i, j
      ENDDO
   ENDDO
ELSE
!
!  Get header
   CALL read_temp%crystal%get_header_from_crystal() ! Read the header
   cr_icc = rd_icc               ! Restore crystal dimensions
!
   CALL get_symmetry_matrices                       ! Setup symmetry
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
mole_num_curr = 0                                  ! Start with no molecules
main: do ia = 1, natoms
   werte    = 0.0
   cr_natoms = cr_natoms + 1
   CALL read_temp%crystal%get_cryst_atom ( ia, itype, posit, iprop, isurface, magn_mom, iin_mole)
   CALL read_temp%crystal%get_cryst_scat ( ia, itype, at_name , dw1, occ1  )
   mole_exist: if(n_mole > 0) THEN
   CALL read_temp%crystal%get_cryst_mole ( ia, i_mole, i_type,  &
                 i_char, c_file, r_fuzzy, r_dens, r_biso, r_clin, r_cqua)
      in_mole: IF ( temp_look(ia) /= 0 ) THEN            ! This atom belongs to a molecule
!        IF ( .not. mole_l_on .OR. temp_upd(temp_look(ia))> mole_num_curr ) THEN  ! Right now we are not in a molecule
         IF ( .not. mole_l_on .OR. temp_upd(temp_look(ia))> mole_num_mole ) THEN  ! Right now we are not in a molecule
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
      nw_name = at_name
      if(uni_mask(0)) then       ! Use unique mask
               itype(1) = atom_get_type(MAXSCAT, 0, cr_nscat, MAXMASK,   &
                                 cr_at_lis, cr_dw, cr_occ,        &
                                 nw_name, dw1, occ1, uni_mask)
      else
         itype(1) = -1                  ! No unique mask, ==> New atom type
         itype(2:3)          = 1
      endif
      IF ( itype(1) == -1 ) THEN
      cr_nscat = cr_nscat + 1
      if(cr_nscat>ubound(cr_at_lis,1)) then
         new_nscat =  cr_nscat
         call alloc_crystal_scat(new_nscat)
      endif
      cr_at_lis(cr_nscat) = nw_name
      cr_dw    (cr_nscat) = dw1
      cr_occ   (cr_nscat) = occ1
      itype(1)            = cr_nscat
      itype(2:3)          = 1
      endif
   ELSE scat_dw                                ! check for previous atom types
      nw_name = at_name
      if(uni_mask(0)) then       ! Use unique mask
               itype(1) = atom_get_type(MAXSCAT, 0, cr_nscat, MAXMASK,   &
                                 cr_at_lis, cr_dw, cr_occ,        &
                                 nw_name, dw1, occ1, uni_mask)
      else
         itype(1) = -1                  ! No unique mask, ==> New atom type
         itype(2:3)          = 1
      endif
!     itype = -1
!     do_scat_dw: DO k = 1,cr_nscat
!        IF( at_name == cr_at_lis(k) .AND.  &
!            dw1     == cr_dw    (k) .AND.  &
!            occ1    == cr_occ   (k)   ) THEN
!           itype = k
!           EXIT do_scat_dw
!        ENDIF
!     ENDDO do_scat_dw
      IF ( itype(1) == -1 ) THEN
         cr_nscat = cr_nscat + 1               ! Atom type does not exist, create a new one
         cr_at_lis(cr_nscat) = at_name
         cr_dw    (cr_nscat) = dw1
         cr_occ   (cr_nscat) = occ1
         itype(1)               = cr_nscat
         itype(2:3)          = 1
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
   cr_iscat(:,cr_natoms) = itype(:)                 ! set the atom type
   if(mole_l_on) then
      mole_num_curr = temp_look(ia)
   endif
!     cr_mole (cr_natoms) = i_mole                ! set the molecule number
!  else
      cr_mole (cr_natoms) = 0                     ! set the molecule number
!  endif
   cr_surf(:,cr_natoms) = isurface               ! set the property flag
   cr_magn(:,cr_natoms) = magn_mom               ! set magnetic vector
   cr_prop (cr_natoms) = iprop                 ! set the property flag
   j = mole_num_mole                           ! temporarily store number of molecules
   CALL symmetry
   if(mole_num_mole>j) then                    ! New molecules were generated update lookup
      do k=(temp_look(ia))+1, temp_num_mole
         temp_upd(k) = temp_upd(k) + mole_num_mole - j
      enddo
   endif
ENDDO main
!
!
CALL no_error 
!                                                                       
!  move first unit cell into lower left corner of crystal          
!                                                                       
DO i = 1, cr_natoms 
   DO j = 1, 3 
      cr_pos (j, i) = cr_pos (j, i) - int ( (cr_icc (j) ) / 2) 
      cr_dim (j, 1) = min (cr_dim (j, 1), cr_pos (j, i) ) 
      cr_dim (j, 2) = max (cr_dim (j, 2), cr_pos (j, i) ) 
   ENDDO 
ENDDO 
!
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
   DEALLOCATE(temp_upd )
ENDIF
!
chem_purge = .FALSE.                          ! No purge, period boundary is OK
!write(*,*) 'ICELL ', dim_natoms, dim_nscat, dim_ncatoms, dim_nanis,&
!dim_n_mole, dim_n_type, dim_n_atom
!write(*,*) 'ICell ', natoms, nscat, ncatoms, nanis,&
!n_mole, n_type, n_atom
!write(*,*) 'Icell ', NMAX, MAXSCAT
!
END SUBROUTINE readcell_internal
!
!*******************************************************************************
   SUBROUTINE stru_readheader_internal (rd_strucfile, DIM_MAXSCAT, DIM_NCATOMS, &
            DIM_NANIS, &
            rd_cr_name,   &
            rd_cr_spcgr, rd_cr_spcgr_set, rd_cr_set, rd_cr_iset,                &
            rd_cr_at_lis, rd_cr_nscat, rd_cr_dw, rd_cr_occ, rd_cr_anis, rd_cr_is_sym,      &
            rd_nanis, rd_cr_anis_full, rd_cr_prin, &
            rd_cr_a0, rd_cr_win, &
            rd_sav_ncell, rd_sav_r_ncell, rd_sav_ncatoms, rd_spcgr_ianz, rd_spcgr_para,   &
            rd_GEN_ADD_MAX, rd_gen_add_n, rd_gen_add_power, rd_gen_add,                 &
            rd_SYM_ADD_MAX, rd_sym_add_n, rd_sym_add_power, rd_sym_add )
!
!  Copies the header info from the internal crystal 'rd_strucfile' into local variables
!
use precision_mod
!
   IMPLICIT NONE
!
   CHARACTER (LEN=  * )                         , INTENT(IN   ) :: rd_strucfile 
   INTEGER                                      , INTENT(INOUT) :: DIM_MAXSCAT 
   INTEGER                                      , INTENT(INOUT) :: DIM_NCATOMS 
   INTEGER                                      , INTENT(INOUT) :: DIM_NANIS 
   INTEGER                                      , INTENT(INOUT) :: rd_sav_ncatoms 
!
   CHARACTER (LEN=  80)                         , INTENT(INOUT) :: rd_cr_name 
   CHARACTER (LEN=  16)                         , INTENT(INOUT) :: rd_cr_spcgr 
   CHARACTER (LEN=  16)                         , INTENT(INOUT) :: rd_cr_spcgr_set
   CHARACTER (LEN=   3)                         , INTENT(INOUT) :: rd_cr_set
   INTEGER                                      , INTENT(INOUT) :: rd_cr_iset
   REAL(kind=PREC_DP)  , DIMENSION(3)           , INTENT(INOUT) :: rd_cr_a0
   REAL(kind=PREC_DP)  , DIMENSION(3)           , INTENT(INOUT) :: rd_cr_win
   INTEGER                                      , INTENT(INOUT) :: rd_cr_nscat 
   REAL(kind=prec_DP)  , DIMENSION(0:DIM_MAXSCAT), INTENT(INOUT) :: rd_cr_dw     ! (0:MAXSCAT) 
   REAL(kind=prec_DP), DIMENSION(6,0:DIM_MAXSCAT), intent(inout) :: rd_cr_anis   ! (0:MAXSCAT) 
   integer             ,dimension(DIM_NCATOMS)   ,intent(inout) :: rd_cr_is_sym
   integer                                      , intent(in)    :: rd_nanis
   real(kind=prec_DP), dimension(6,1:DIM_NANIS)     , intent(inout) :: rd_cr_anis_full   ! (0:MAXSCAT) 
   real(kind=prec_DP), dimension(4,3,1:DIM_NANIS)   , intent(inout) :: rd_cr_prin        ! (0:MAXSCAT) 
   REAL(kind=PREC_DP)  , DIMENSION(0:DIM_MAXSCAT), INTENT(INOUT) :: rd_cr_occ    ! (0:MAXSCAT) 
   CHARACTER (LEN=   4), DIMENSION(0:DIM_MAXSCAT), INTENT(INOUT) :: rd_cr_at_lis ! (0:MAXSCAT) 
   INTEGER             , DIMENSION(3)           , INTENT(INOUT) :: rd_sav_ncell ! (3) 
   LOGICAL                                      , INTENT(INOUT) :: rd_sav_r_ncell 
   INTEGER                                      , INTENT(INOUT) :: rd_spcgr_ianz 
   INTEGER                                      , INTENT(INOUT) :: rd_spcgr_para 
!
   INTEGER             ::  rd_GEN_ADD_MAX
   INTEGER             ::  rd_gen_add_n
   INTEGER             ::  rd_gen_add_power(rd_GEN_ADD_MAX)
   REAL(kind=PREC_DP)  ::  rd_gen_add(4,4,0:rd_GEN_ADD_MAX)
!
   INTEGER             ::  rd_SYM_ADD_MAX
   INTEGER             ::  rd_sym_add_n
   INTEGER             ::  rd_sym_add_power(rd_SYM_ADD_MAX)
   REAL(kind=PREC_DP)  ::  rd_sym_add(4,4,0:rd_SYM_ADD_MAX)
!
   NULLIFY(read_from)
   NULLIFY(read_parent)
   CALL store_find_node(store_root, read_from, rd_strucfile, read_temp, read_parent, ier_num ) ! Find the proper node
   IF ( ier_num /= 0) THEN
      ier_typ = ER_APPL
      RETURN
   ENDIF
!
CALL read_temp%crystal%get_header_to_local (DIM_MAXSCAT, DIM_NCATOMS, DIM_NANIS, &
   rd_cr_name,      &
            rd_cr_spcgr, rd_cr_spcgr_set, rd_cr_set, rd_cr_iset,            &
            rd_cr_at_lis, rd_cr_nscat, rd_cr_dw, rd_cr_anis, rd_cr_occ, &
            rd_cr_is_sym, rd_nanis, rd_cr_anis_full, rd_cr_prin, rd_cr_a0, rd_cr_win, &
            rd_sav_ncell, rd_sav_r_ncell, rd_sav_ncatoms, rd_spcgr_ianz, rd_spcgr_para, &
            rd_GEN_ADD_MAX, rd_gen_add_n, rd_gen_add_power, rd_gen_add,                 &
            rd_SYM_ADD_MAX, rd_sym_add_n, rd_sym_add_power, rd_sym_add )
!
   END SUBROUTINE stru_readheader_internal
!
!*******************************************************************************
!
subroutine stru_get_anis(rd_strucfile, DIM_NANIS, rd_nanis, rd_cr_anis_full, rd_cr_prin)
!-
!  Get the files related to anisotropic parameters
!+
!
use errlist_mod
use precision_mod
!
implicit none
!
character(len=*)                              , intent(in)    :: rd_strucfile 
integer                                       , intent(in)    :: DIM_NANIS
integer                                       , intent(out)   :: rd_nanis
real(kind=prec_DP), dimension(6,1:DIM_NANIS)  , intent(inout) :: rd_cr_anis_full   ! (0:MAXSCAT) 
real(kind=prec_DP), dimension(4,3,1:DIM_NANIS), intent(inout) :: rd_cr_prin        ! (0:MAXSCAT) 
!
integer :: lanis   ! Local querry for nanis
!
nullify(read_from)
nullify(read_parent)
call store_find_node(store_root, read_from, rd_strucfile, read_temp, read_parent, ier_num ) ! Find the proper node
if(ier_num /= 0) THEN
   ier_typ = ER_APPL
   return
endif
!
lanis = read_temp%crystal%get_nanis()
if(lanis>rd_nanis) then
   ier_msg(1) = 'More ADPs in internal storage than expected'
   ier_num = -6
   ier_typ = ER_COMM
   return
endif
call read_temp%crystal%get_anis(DIM_NANIS, rd_nanis, rd_cr_anis_full, rd_cr_prin)
!
end subroutine stru_get_anis
!
!*******************************************************************************
!
   SUBROUTINE struc_read_atoms_internal(strucfile, RD_NMAX, &
              rd_cr_natoms, rd_cr_pos, rd_cr_iscat, rd_cr_prop, &
              rd_cr_surf, rd_cr_magn, rd_cr_mole )
!
!  This subroutine adds all atoms from the internal storage to the local
!  crystal rd_cr_pos. The number of atoms cr_natoms is incremented accordingly
!
   USE discus_allocate_appl_mod
!
use precision_mod
!
   IMPLICIT NONE
!
   CHARACTER (LEN=  * )                         , INTENT(IN   ) :: strucfile 
   INTEGER                                      , INTENT(IN   ) :: rd_NMAX 
   INTEGER                                      , INTENT(INOUT) :: rd_cr_natoms
   REAL(kind=PREC_DP)  , DIMENSION(3,1:RD_NMAX) , INTENT(INOUT) :: rd_cr_pos
   INTEGER             , DIMENSION(3,1:RD_NMAX) , INTENT(INOUT) :: rd_cr_iscat
   INTEGER             , DIMENSION(0:3,1:RD_NMAX),INTENT(INOUT) :: rd_cr_surf
   REAL(kind=PREC_DP)  , DIMENSION(0:3,1:RD_NMAX),INTENT(INOUT) :: rd_cr_magn
   INTEGER             , DIMENSION(  1:RD_NMAX) , INTENT(INOUT) :: rd_cr_prop
   INTEGER             , DIMENSION(  1:RD_NMAX) , INTENT(INOUT) :: rd_cr_mole
!
integer                       :: dim_natoms   ! Number of atoms in the structure
integer                       :: dim_ncatoms  ! Number of atoms per unit call in the structure
integer                       :: dim_nscat    ! Number of different atom types
integer                       :: dim_nanis    ! Number of ADPs
integer                       :: dim_n_mole
integer                       :: dim_n_type
integer                       :: dim_n_atom
   INTEGER                       :: i
   INTEGER                       :: natoms
   INTEGER                       :: ncatoms
   INTEGER                       :: nscat
   INTEGER                       :: nanis
   INTEGER                       :: n_mole
   INTEGER                       :: n_type
   INTEGER                       :: n_atom
   INTEGER           , dimension(3) :: itype
   REAL(kind=PREC_DP), DIMENSION(3)         :: posit
   INTEGER, DIMENSION(0:3)       :: isurface
   REAL(kind=PREC_DP)   , DIMENSION(0:3)       :: magn_mom
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
   CALL readstru_size_int(strucfile, dim_natoms, dim_ncatoms, dim_nscat, dim_nanis, &
                                  dim_n_mole, dim_n_type, dim_n_atom,            &
                                  natoms, ncatoms, nscat, nanis,                 &
                                  n_mole, n_type, n_atom)
!  CALL readstru_size_int(strucfile, natoms, ncatoms, &           ! Get the number of atoms
!             nscat, nanis, n_mole, n_type, n_atom)
   DO i=1,natoms
      CALL read_temp%crystal%get_cryst_atom(i, itype, posit, iprop, isurface, magn_mom, iin_mole)
      rd_cr_natoms = rd_cr_natoms + 1
      rd_cr_pos(1,rd_cr_natoms) = posit(1)
      rd_cr_pos(2,rd_cr_natoms) = posit(2)
      rd_cr_pos(3,rd_cr_natoms) = posit(3)
      rd_cr_iscat(:,rd_cr_natoms) = itype(:)
      rd_cr_prop (rd_cr_natoms) = iprop 
      rd_cr_magn(:,rd_cr_natoms) = magn_mom
      rd_cr_mole (rd_cr_natoms) = iin_mole(1) 
   ENDDO
!
   END SUBROUTINE struc_read_atoms_internal
!*******************************************************************************
   SUBROUTINE struc_read_one_atom_internal(strucfile, iatom,  &
              rd_cr_pos, rd_cr_iscat, rd_cr_prop, rd_cr_surf, &
              rd_cr_magn, rd_cr_mole, rd_cr_moleatom )
!
!  This subroutine reads just atom "iatom" from the internal storage.
!  The atom is copied into internal storage
!
use precision_mod
!
   IMPLICIT NONE
!
   CHARACTER (LEN=  * )               , INTENT(IN   ) :: strucfile 
   INTEGER                            , INTENT(IN   ) :: iatom
   REAL(kind=PREC_DP)  , DIMENSION(3) , INTENT(INOUT) :: rd_cr_pos
   INTEGER             , dimension(3) , INTENT(INOUT) :: rd_cr_iscat
   INTEGER                            , INTENT(INOUT) :: rd_cr_prop
   INTEGER             , DIMENSION(0:3),INTENT(INOUT) :: rd_cr_surf
   REAL(kind=PREC_DP)  , DIMENSION(0:3),INTENT(INOUT) :: rd_cr_magn
   INTEGER                             ,INTENT(INOUT) :: rd_cr_mole
   INTEGER                             ,INTENT(INOUT) :: rd_cr_moleatom
!
   INTEGER                       :: i
integer                       :: dim_natoms   ! Number of atoms in the structure
integer                       :: dim_ncatoms  ! Number of atoms per unit call in the structure
integer                       :: dim_nscat    ! Number of different atom types
integer                       :: dim_nanis    ! Number of ADPs
integer                       :: dim_n_mole
integer                       :: dim_n_type
integer                       :: dim_n_atom
   INTEGER                       :: natoms
   INTEGER                       :: ncatoms
   INTEGER                       :: nscat
   INTEGER                       :: nanis
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
   CALL readstru_size_int(strucfile, dim_natoms, dim_ncatoms, dim_nscat, dim_nanis, &
                                  dim_n_mole, dim_n_type, dim_n_atom,            &
                                  natoms, ncatoms, nscat, nanis,                 &
                                  n_mole, n_type, n_atom)
!  CALL readstru_size_int(strucfile, natoms, ncatoms, &           ! Get the number of atoms
!             nscat, nanis, n_mole, n_type, n_atom)
   IF ( iatom <= natoms ) THEN
      i = iatom
      CALL read_temp%crystal%get_cryst_atom(i, rd_cr_iscat, rd_cr_pos, &
           rd_cr_prop, rd_cr_surf, rd_cr_magn, iin_mole)
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
!  USE discus_allocate_appl_mod
!
use precision_mod
!
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
   REAL(kind=PREC_DP)   ,           DIMENSION(0:MOLE_MAX_MOLE), INTENT(OUT) :: mole_dens
   REAL(kind=PREC_DP)   ,           DIMENSION(0:MOLE_MAX_TYPE), INTENT(OUT) :: mole_biso
   REAL(kind=PREC_DP)   ,           DIMENSION(0:MOLE_MAX_TYPE), INTENT(OUT) :: mole_clin
   REAL(kind=PREC_DP)   ,           DIMENSION(0:MOLE_MAX_TYPE), INTENT(OUT) :: mole_cqua
   REAL(kind=PREC_DP)   ,           DIMENSION(0:MOLE_MAX_MOLE), INTENT(OUT) :: mole_fuzzy
   INTEGER,           DIMENSION(0:MOLE_MAX_ATOM), INTENT(OUT) :: mole_cont
!
integer                       :: dim_natoms   ! Number of atoms in the structure
integer                       :: dim_ncatoms  ! Number of atoms per unit call in the structure
integer                       :: dim_nscat    ! Number of different atom types
integer                       :: dim_nanis    ! Number of ADPs
integer                       :: dim_n_mole
integer                       :: dim_n_type
integer                       :: dim_n_atom
   INTEGER                       :: natoms
   INTEGER                       :: ncatoms
   INTEGER                       :: nscat
   INTEGER                       :: nanis
!
!  ALLOCATE(read_temp, STAT = istatus )                  ! Allocate a temporary storage
!  IF ( istatus /= 0) THEN
!     ier_num = -114
!     ier_typ = ER_APPL
!     RETURN
!  ENDIF
   CALL readstru_size_int(strucfile, dim_natoms, dim_ncatoms, dim_nscat, dim_nanis, &
                                  dim_n_mole, dim_n_type, dim_n_atom,            &
                                  natoms, ncatoms, nscat, nanis,                 &
                                  mole_num_mole, mole_num_type, mole_num_atom)
!  CALL readstru_size_int(strucfile, natoms, ncatoms, & 
!                  nscat, nanis, mole_num_mole, mole_num_type, mole_num_atom)
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
   SUBROUTINE testfile_internal(strucfile , natoms, ncatoms, & 
              nscat, nanis, n_mole, n_type, n_atom)
!
!  Reads the header section of an internal file and determines 
!  the required sizes
!
   IMPLICIT NONE
!
!
   CHARACTER (LEN=*), INTENT(IN) :: strucfile
!
integer                       :: dim_natoms   ! Number of atoms in the structure
integer                       :: dim_ncatoms  ! Number of atoms per unit call in the structure
integer                       :: dim_nscat    ! Number of different atom types
integer                       :: dim_nanis    ! Number of ADPs
integer                       :: dim_n_mole
integer                       :: dim_n_type
integer                       :: dim_n_atom
   INTEGER, INTENT(INOUT)        :: natoms
   INTEGER, INTENT(INOUT)        :: ncatoms
   INTEGER, INTENT(INOUT)        :: nscat
   INTEGER, intent(inout)        :: nanis
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
   CALL readstru_size_int(strucfile, dim_natoms, dim_ncatoms, dim_nscat, dim_nanis, &
                                  dim_n_mole, dim_n_type, dim_n_atom,            &
                                  natoms, ncatoms, nscat, nanis,                 &
                                  n_mole, n_type, n_atom)
!  CALL readstru_size_int(strucfile, natoms, ncatoms, & 
!                  nscat, nanis, n_mole, n_type, n_atom)
!
!  DEALLOCATE(read_temp, STAT = istatus )        ! Deallocate a temporary storage
!
   END SUBROUTINE testfile_internal
!
!*******************************************************************************
!
subroutine stru_internal_get_cr_dim(strucfile, cr_dim)
!-
!  Read the crystal dimensions of the internal structure strucfile
!+
!
use precision_mod
!
character(len=*), intent(in) :: strucfile
real(kind=PREC_DP), dimension(3,2), intent(out) :: cr_dim
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
call read_temp%crystal%get_cr_dim_from_crystal(cr_dim)
nullify(read_from)
nullify(read_parent)
nullify(read_temp)
!
end subroutine stru_internal_get_cr_dim
!*******************************************************************************

END MODULE read_internal_mod
