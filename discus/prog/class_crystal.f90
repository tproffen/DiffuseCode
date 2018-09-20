MODULE cryst_class
!
!  Defines a type "cl_cryst".
!  This very incomplete example just has the natoms i.e. 
!  the number of atoms and an array of type "cl_atom"
!
!  Just a few procedures are bound to this class:
!      alloc_atoms      allocate the atom array
!      get_natoms       Returns the number of atoms in this crystal
!      set_cryst_atom   Sets all values for an atom
!      get_cryst_atom   Gets all values for an atom
!
!  A working class will need the full suite of methods....
!
USE atom_class          ! this way we know all interfaces to the "cl_atom" class
USE errlist_mod
!
IMPLICIT none 
!
TYPE :: cl_cryst        ! Define a type "cl_cryst"
   PRIVATE              ! hide local variables from outside; access through procedures only
!
   CHARACTER (LEN=200)                              ::  cr_file     = 'crystal.stru' ! Crystal file name
   CHARACTER (LEN=80)                               ::  cr_name     = 'crystal' ! Crystal name
   CHARACTER (LEN=16)                               ::  cr_spcgr    = 'P1'      ! Space group symbol
   INTEGER                                          ::  cr_spcgrno  = 1 ! Space group number
   INTEGER                                          ::  cr_syst     = 1 ! Crystal system
   LOGICAL                                          ::  cr_acentric = .true.
   INTEGER                                          ::  cr_ncatoms  = 1 ! Atoms in unit cell
   INTEGER           , DIMENSION(3)                 ::  cr_icc      = 1 ! Number of unit cells
!
   INTEGER                                          ::  spcgr_ianz  = 0 ! Number of parameters on spcg line
   INTEGER                                          ::  spcgr_para  = 0 ! 
   INTEGER           , DIMENSION(0:1)               ::  cr_sel_prop = 0 ! Global selection properties
!
   REAL              , DIMENSION(3)                 ::  cr_a0   =  1.0  ! unit cell lengths
   REAL              , DIMENSION(3)                 ::  cr_win  = 90.0  ! unit cell angles
   REAL              , DIMENSION(3,2)               ::  cr_dim  = RESHAPE((/0,0,0, 1,1,1/),(/3,2/))
   REAL              , DIMENSION(3,2)               ::  cr_dim0 = RESHAPE((/0,0,0, 1,1,1/),(/3,2/))
   REAL                                             ::  cr_v    = 1.    ! Unit cell volume
   REAL              , DIMENSION(3,3)               ::  cr_gten = RESHAPE((/1,0,0, 0,1,0, 0,0,1/),(/3,3/))
   REAL              , DIMENSION(3,3,3)             ::  cr_eps          ! Epsilon tensor
   REAL              , DIMENSION(3,3)               ::  cr_gmat = RESHAPE((/1,0,0, 0,1,0, 0,0,1/),(/3,3/))
   REAL              , DIMENSION(3)                 ::  cr_ar   = 1.0   ! reciprocal unit cell lengths
   REAL              , DIMENSION(3)                 ::  cr_wrez = 90.00 ! reciprocal unit cell angles
   REAL                                             ::  cr_vr   = 1.    ! ceiprocal unit cell volume
   REAL              , DIMENSION(3,3)               ::  cr_rten = RESHAPE((/1,0,0, 0,1,0, 0,0,1/),(/3,3/))
   REAL              , DIMENSION(3,3,3)             ::  cr_reps         ! reciprocel epsilon tensor
   REAL              , DIMENSION(3,3)               ::  cr_fmat = RESHAPE((/1,0,0, 0,1,0, 0,0,1/),(/3,3/))
   REAL              , DIMENSION(4,4)               ::  cr_tran_g       ! trans to cartesian
   REAL              , DIMENSION(4,4)               ::  cr_tran_gi      ! trans to cartesian
   REAL              , DIMENSION(4,4)               ::  cr_tran_f       ! trans to cartesian
   REAL              , DIMENSION(4,4)               ::  cr_tran_fi      ! trans to cartesian
!
   INTEGER                                          ::  cr_gen_add_n         = 0
   INTEGER                                          ::  cr_gen_add_power(14) = 1
   REAL              , DIMENSION(:,:,:),ALLOCATABLE ::  cr_gen_add   ! Additional Generator matrices (4,4,0:14)
!
   INTEGER                                          ::  cr_sym_add_n         = 0
   INTEGER                                          ::  cr_sym_add_power(14) = 1
   REAL              , DIMENSION(:,:,:),ALLOCATABLE ::  cr_sym_add   ! Additional Symmetry matrices (4,4,0:14)
!
   REAL              , DIMENSION(:,:), ALLOCATABLE  ::  cr_scat   ! Scattering curves (9,0:MAXSCAT)
   REAL              , DIMENSION(  :), ALLOCATABLE  ::  cr_delfr  ! anomalous scattering real part (  0:MAXSCAT)
   REAL              , DIMENSION(  :), ALLOCATABLE  ::  cr_delfi  ! anomalous scattering imag part (  0:MAXSCAT)
!
   LOGICAL           , DIMENSION(  :), ALLOCATABLE  ::  cr_scat_int  ! (  0:MAXSCAT)
   LOGICAL           , DIMENSION(  :), ALLOCATABLE  ::  cr_scat_equ  ! (  0:MAXSCAT)
   LOGICAL           , DIMENSION(  :), ALLOCATABLE  ::  cr_delf_int  ! (  0:MAXSCAT)
   CHARACTER (LEN=4 ), DIMENSION(  :), ALLOCATABLE  ::  cr_at_equ    ! (  0:MAXSCAT)
   LOGICAL           , DIMENSION(  :), ALLOCATABLE  ::  cr_sav_atom  ! (  0:MAXSCAT)
   LOGICAL                                          ::  cr_newtype = .true.
   LOGICAL                                          ::  cr_cartesian = .false.
   LOGICAL                                          ::  cr_sav_scat  = .false.
   LOGICAL                                          ::  cr_sav_adp   = .false.
   LOGICAL                                          ::  cr_sav_occ   = .false.
   LOGICAL                                          ::  cr_sav_gene  = .true.
   LOGICAL                                          ::  cr_sav_symm  = .true.
   LOGICAL                                          ::  cr_sav_ncell = .false.
   LOGICAL                                          ::  cr_sav_obje  = .true.
   LOGICAL                                          ::  cr_sav_doma  = .true.
   LOGICAL                                          ::  cr_sav_mole  = .true.
   LOGICAL                                          ::  cr_sav_prop  = .true.
   INTEGER           , DIMENSION(0:1)               ::  cr_sav_sel_prop  = (/0,0/)
!
   INTEGER                                          ::  cr_natoms  !the number of atoms
   INTEGER                                          ::  cr_nscat   !the number of atom types
   INTEGER                                          ::  cr_n_REAL_atoms = 0

   REAL              , DIMENSION(  :), ALLOCATABLE  ::  cr_dw      ! (  0:MAXSCAT)
   REAL              , DIMENSION(  :), ALLOCATABLE  ::  cr_occ     ! (  0:MAXSCAT)
   CHARACTER (LEN=4 ), DIMENSION(  :), ALLOCATABLE  ::  cr_at_lis  ! (  0:MAXSCAT)
   TYPE(cl_atom)     , DIMENSION(:)  , POINTER      ::  atoms            !the actual atoms
!
   INTEGER                                          ::  cr_mole_max_mole
   INTEGER                                          ::  cr_mole_max_type
   INTEGER                                          ::  cr_mole_max_atom
   INTEGER                                          ::  cr_num_mole
   INTEGER                                          ::  cr_num_type
   INTEGER                                          ::  cr_num_atom
   INTEGER, DIMENSION(:), ALLOCATABLE               ::  cr_mole_len     ! (  0:MOLE_MAX_MOLE)
   INTEGER, DIMENSION(:), ALLOCATABLE               ::  cr_mole_off     ! (  0:MOLE_MAX_MOLE)
   INTEGER, DIMENSION(:), ALLOCATABLE               ::  cr_mole_type    ! (  0:MOLE_MAX_MOLE)
   INTEGER, DIMENSION(:), ALLOCATABLE               ::  cr_mole_char    ! (  0:MOLE_MAX_MOLE)
   CHARACTER(LEN=200), DIMENSION(:), ALLOCATABLE    ::  cr_mole_file    ! (  0:MOLE_MAX_MOLE)
   INTEGER, DIMENSION(:), ALLOCATABLE               ::  cr_mole_cont    ! (  0:MOLE_MAX_ATOM)
   REAL   , DIMENSION(:), ALLOCATABLE               ::  cr_mole_dens    ! (  0:MOLE_MAX_MOLE)
   REAL   , DIMENSION(:), ALLOCATABLE               ::  cr_mole_biso    ! (  0:MOLE_MAX_TYPE)
   REAL   , DIMENSION(:), ALLOCATABLE               ::  cr_mole_fuzzy   ! (  0:MOLE_MAX_MOLE)
!
   LOGICAL                                          ::  latom = .false.  !have atoms been allocated?
!
CONTAINS
!
   PROCEDURE, PUBLIC, PASS :: alloc_arrays         ! allocate a set of atoms
   PROCEDURE, PUBLIC, PASS :: get_natoms           ! Return number of atoms in this crystal
   PROCEDURE, PUBLIC, PASS :: get_nscat            ! Return number of atom types in this crystal
   PROCEDURE, PUBLIC, PASS :: get_n_mole           ! Return number of molecules types in this crystal
   PROCEDURE, PUBLIC, PASS :: get_n_type           ! Return number of molecule types in this crystal
   PROCEDURE, PUBLIC, PASS :: get_n_atom           ! Return number of molecule atoms in this crystal
   PROCEDURE, PUBLIC, PASS :: set_cryst_atom       ! Set iscat, posit, property
   PROCEDURE, PUBLIC, PASS :: set_cryst_at_lis     ! Set cr_at_lis
   PROCEDURE, PUBLIC, PASS :: set_cryst_dw         ! Set cr_dw
   PROCEDURE, PUBLIC, PASS :: get_cryst_tran_f     ! Get transformation matrix
   PROCEDURE, PUBLIC, PASS :: get_cryst_atom       ! Get iscat, posit, property
   PROCEDURE, PUBLIC, PASS :: get_cryst_mole       ! Get molecule number in which inum is
   PROCEDURE, PUBLIC, PASS :: get_cryst_scat       ! Get iscat, atom_name and DW
!   PROCEDURE, PUBLIC, PASS :: read_crystal         ! Read from disk file
   PROCEDURE, PUBLIC, PASS :: set_crystal_from_standard ! Copy from DISCUS standard into a crystal type
   PROCEDURE, PUBLIC, PASS :: set_crystal_from_local    ! Copy from local copy      into a crystal type
   PROCEDURE, PUBLIC, PASS :: set_crystal_save_flags    ! Set the "save" flags in the crystal type
   PROCEDURE, PUBLIC, PASS :: get_header_from_crystal   ! Copy from a crystal type into DISCUS standard
   PROCEDURE, PUBLIC, PASS :: get_header_to_local       ! Copy from a crystal type into a local copy
   PROCEDURE, PUBLIC, PASS :: get_atoms_from_crystal    ! Copy from a crystal type into DISCUS standard
   PROCEDURE, PUBLIC, PASS :: get_atoms_to_local        ! Copy from a crystal type into a local copy
   PROCEDURE, PUBLIC, PASS :: get_molecules_from_crystal ! Copy from a crystal type into DISCUS standard
   PROCEDURE, PUBLIC, PASS :: is_alloc_atom        ! Has atom array been allocated ?
   PROCEDURE, PUBLIC       :: finalize_atoms
!   FINAL :: finalize_atoms                        ! NOT YET IN GFORTRAN!
END TYPE cl_cryst
!
CONTAINS
!
!******************************************************************************
   SUBROUTINE alloc_arrays   ( this, natoms, nscat, n_mole, n_type, n_atom)
!
!  Allocate the arrays for "this" crystal 
!  Initialize all variables
!
   IMPLICIT none
!
   CLASS (cl_cryst)                 :: this        ! Work on "this" crystal
   INTEGER,              INTENT(IN) :: natoms      ! Number of atoms for this crystal
   INTEGER,              INTENT(IN) :: nscat       ! Number of atom types for this crystal
   INTEGER,              INTENT(IN) :: n_mole      ! Number of molecules  for this crystal
   INTEGER,              INTENT(IN) :: n_type      ! Number of molecule types  for this crystal
   INTEGER,              INTENT(IN) :: n_atom      ! Number of atoms in molecules for this crystal
!
   INTEGER  :: istatus
!
   IF(this%latom) THEN                             ! Crystal was allocated 
      DEALLOCATE ( this%cr_gen_add , STAT=istatus ) ! Allocate equivalent atom names
      DEALLOCATE ( this%cr_sym_add , STAT=istatus ) ! Allocate equivalent atom names
      DEALLOCATE ( this%cr_scat    , STAT=istatus ) ! Allocate equivalent atom names
      DEALLOCATE ( this%cr_delfr   , STAT=istatus ) ! Allocate equivalent atom names
      DEALLOCATE ( this%cr_delfi   , STAT=istatus ) ! Allocate equivalent atom names
      DEALLOCATE ( this%cr_scat_int, STAT=istatus ) ! Allocate equivalent atom names
      DEALLOCATE ( this%cr_scat_equ, STAT=istatus ) ! Allocate equivalent atom names
      DEALLOCATE ( this%cr_delf_int, STAT=istatus ) ! Allocate equivalent atom names
      DEALLOCATE ( this%cr_at_equ  , STAT=istatus ) ! Always deallocate
      DEALLOCATE ( this%cr_sav_atom, STAT=istatus ) ! Always deallocate
      DEALLOCATE ( this%cr_at_lis  , STAT=istatus ) ! Always deallocate
      DEALLOCATE ( this%cr_dw      , STAT=istatus ) ! Always deallocate
      DEALLOCATE ( this%cr_occ     , STAT=istatus ) ! Always deallocate
      DEALLOCATE ( this%atoms      , STAT=istatus ) ! Always deallocate
      DEALLOCATE ( this%cr_mole_len  , STAT=istatus ) ! Always deallocate molecules
      DEALLOCATE ( this%cr_mole_off  , STAT=istatus ) ! Always deallocate molecules
      DEALLOCATE ( this%cr_mole_type , STAT=istatus ) ! Always deallocate molecules
      DEALLOCATE ( this%cr_mole_char , STAT=istatus ) ! Always deallocate molecules
      DEALLOCATE ( this%cr_mole_file , STAT=istatus ) ! Always deallocate molecules
      DEALLOCATE ( this%cr_mole_cont , STAT=istatus ) ! Always deallocate molecules
      DEALLOCATE ( this%cr_mole_dens , STAT=istatus ) ! Always deallocate molecules
      DEALLOCATE ( this%cr_mole_biso , STAT=istatus ) ! Always deallocate molecules
      DEALLOCATE ( this%cr_mole_fuzzy, STAT=istatus ) ! Always deallocate molecules
   ENDIF
   ALLOCATE ( this%cr_gen_add (4,4,0:14 ), STAT=istatus )  ! Allocate equivalent atom names
   ALLOCATE ( this%cr_sym_add (4,4,0:14 ), STAT=istatus )  ! Allocate equivalent atom names
   ALLOCATE ( this%cr_scat    (11,0:nscat), STAT=istatus ) ! Allocate equivalent atom names
   ALLOCATE ( this%cr_delfr   (  0:nscat), STAT=istatus ) ! Allocate equivalent atom names
   ALLOCATE ( this%cr_delfi   (  0:nscat), STAT=istatus ) ! Allocate equivalent atom names
   ALLOCATE ( this%cr_scat_int(  0:nscat), STAT=istatus ) ! Allocate equivalent atom names
   ALLOCATE ( this%cr_scat_equ(  0:nscat), STAT=istatus ) ! Allocate equivalent atom names
   ALLOCATE ( this%cr_delf_int(  0:nscat), STAT=istatus ) ! Allocate equivalent atom names
   ALLOCATE ( this%cr_at_equ  (  0:nscat), STAT=istatus ) ! Allocate equivalent atom names
   ALLOCATE ( this%cr_sav_atom(  0:nscat), STAT=istatus ) ! Allocate equivalent atom names
   ALLOCATE ( this%cr_at_lis  (  0:nscat), STAT=istatus ) ! Allocate atom names
   ALLOCATE ( this%cr_dw      (  0:nscat), STAT=istatus ) ! Allocate Debye Waller terms
   ALLOCATE ( this%cr_occ     (  0:nscat), STAT=istatus ) ! Allocate Debye Waller terms
   ALLOCATE ( this%atoms      (natoms   ), STAT=istatus ) ! Allocate list of atoms
   ALLOCATE ( this%cr_mole_len  (0:n_mole), STAT=istatus ) ! Allocate molecules
   ALLOCATE ( this%cr_mole_off  (0:n_mole), STAT=istatus ) ! Allocate molecules
   ALLOCATE ( this%cr_mole_type (0:n_mole), STAT=istatus ) ! Allocate molecules
   ALLOCATE ( this%cr_mole_char (0:n_mole), STAT=istatus ) ! Allocate molecules
   ALLOCATE ( this%cr_mole_file (0:n_mole), STAT=istatus ) ! Allocate molecules
   ALLOCATE ( this%cr_mole_cont (0:n_atom), STAT=istatus ) ! Allocate molecules
   ALLOCATE ( this%cr_mole_dens (0:n_mole), STAT=istatus ) ! Allocate molecules
   ALLOCATE ( this%cr_mole_biso (0:n_type), STAT=istatus ) ! Allocate molecules
   ALLOCATE ( this%cr_mole_fuzzy(0:n_mole), STAT=istatus ) ! Allocate molecules
!
   this%cr_gen_add (4,4,0:14 ) = 0.0
   this%cr_sym_add (4,4,0:14 ) = 0.0
   this%cr_scat    (11,0:nscat) = 0.0
   this%cr_delfr   (  0:nscat) = 0.0
   this%cr_delfi   (  0:nscat) = 0.0
   this%cr_scat_int(  0:nscat) = .TRUE.
   this%cr_scat_equ(  0:nscat) = .FALSE.
   this%cr_delf_int(  0:nscat) = .FALSE.
   this%cr_at_equ  (  0:nscat) = ' '
   this%cr_sav_atom(  0:nscat) = .TRUE.
   this%cr_at_lis  (  0:nscat) = ' '
   this%cr_dw      (  0:nscat) = 0.0
   this%cr_occ     (  0:nscat) = 0.0
!  this%atoms      (natoms   )
   this%cr_mole_len  (0:n_mole) = 0
   this%cr_mole_off  (0:n_mole) = 0
   this%cr_mole_type (0:n_mole) = 0
   this%cr_mole_char (0:n_mole) = 0
   this%cr_mole_file (0:n_mole) = ' '
   this%cr_mole_cont (0:n_atom) = 0
   this%cr_mole_dens (0:n_mole) = 0.0
   this%cr_mole_biso (0:n_type) = 0.0
   this%cr_mole_fuzzy(0:n_mole) = 0.0
!
   this%cr_natoms = natoms                            ! Store number of atoms
   this%cr_nscat  = nscat                             ! Store number of atom types
   this%latom  = .true.                            ! Flag that crystal is allocated
!
   END SUBROUTINE alloc_arrays
!******************************************************************************
   INTEGER FUNCTION get_natoms (this )
!
!  Return the number of atoms in "this" crystal
!
   IMPLICIT none
!
   CLASS (cl_cryst) :: this
!
   get_natoms = this%cr_natoms
!
   END FUNCTION get_natoms
!******************************************************************************
   INTEGER FUNCTION get_nscat (this )
!
!  Return the number of atoms in "this" crystal
!
   IMPLICIT none
!
   CLASS (cl_cryst) :: this
!
   get_nscat = this%cr_nscat
!
   END FUNCTION get_nscat
!******************************************************************************
   INTEGER FUNCTION get_n_mole (this )
!
!  Return the number of molecules in "this" crystal
!
   IMPLICIT none
!
   CLASS (cl_cryst) :: this
!
   get_n_mole = this%cr_num_mole
!
   END FUNCTION get_n_mole
!******************************************************************************
   INTEGER FUNCTION get_n_type (this )
!
!  Return the number of molecule types in "this" crystal
!
   IMPLICIT none
!
   CLASS (cl_cryst) :: this
!
   get_n_type = this%cr_num_type
!
   END FUNCTION get_n_type
!******************************************************************************
   INTEGER FUNCTION get_n_atom (this )
!
!  Return the number of molecule atoms in "this" crystal
!
   IMPLICIT none
!
   CLASS (cl_cryst) :: this
!
   get_n_atom = this%cr_num_atom
!
   END FUNCTION get_n_atom
!******************************************************************************
   SUBROUTINE set_cryst_atom ( this, inum, itype, posit, iprop, isurface, iin_mole)
!
!  Fully define atom No. <inum> for "this" crystal
!
   IMPLICIT none
!
   CLASS (cl_cryst)                 :: this
   INTEGER,              INTENT(IN) :: inum
   INTEGER,              INTENT(IN) :: itype
   REAL   , DIMENSION(3),INTENT(IN) :: posit
   INTEGER,              INTENT(IN) :: iprop
   INTEGER, DIMENSION(0:3),INTENT(IN) :: isurface
   INTEGER, DIMENSION(1:2),INTENT(IN) :: iin_mole
!
   CALL this%atoms(inum)%set_atom ( itype, posit, iprop, isurface, iin_mole )
!
   END SUBROUTINE set_cryst_atom 
!******************************************************************************
   SUBROUTINE set_cryst_at_lis ( this, MAXSCAT, inum, cr_at_lis)
!
!  Set the atom name array for "this" crystal
!
   IMPLICIT none
!
   CLASS (cl_cryst)                 :: this
   INTEGER,              INTENT(IN) :: MAXSCAT
   INTEGER,              INTENT(IN) :: inum
   CHARACTER (LEN=4),DIMENSION(0:MAXSCAT), INTENT(IN) :: cr_at_lis
!
   INTEGER :: i
!
   DO i=0, inum
      this%cr_at_lis(i) = cr_at_lis(i)
   ENDDO
!
   END SUBROUTINE set_cryst_at_lis 
!******************************************************************************
   SUBROUTINE set_cryst_dw ( this, MAXSCAT, inum, cr_dw)
!
!  Set the atom name array for "this" crystal
!
   IMPLICIT none
!
   CLASS (cl_cryst)                 :: this
   INTEGER,              INTENT(IN) :: MAXSCAT
   INTEGER,              INTENT(IN) :: inum
   REAL   ,DIMENSION(0:MAXSCAT), INTENT(IN) :: cr_dw
!
   INTEGER :: i
!
   DO i=0, inum
      this%cr_dw(i) = cr_dw(i)
   ENDDO
!
   END SUBROUTINE set_cryst_dw 
!******************************************************************************
   SUBROUTINE set_cryst_occ( this, MAXSCAT, inum, cr_occ)
!
!  Set the atom name array for "this" crystal
!
   IMPLICIT none
!
   CLASS (cl_cryst)                 :: this
   INTEGER,              INTENT(IN) :: MAXSCAT
   INTEGER,              INTENT(IN) :: inum
   REAL   ,DIMENSION(0:MAXSCAT), INTENT(IN) :: cr_occ
!
   INTEGER :: i
!
   DO i=0, inum
      this%cr_occ(i) = cr_occ(i)
   ENDDO
!
   END SUBROUTINE set_cryst_occ 
!******************************************************************************
   SUBROUTINE get_cryst_tran_f ( this, cr_tran_f)
!
!  Get the transformation matrix to cartesian space
!
   IMPLICIT NONE
!
!
   CLASS (cl_cryst)                  :: this
   REAL, DIMENSION(4,4), INTENT(OUT) :: cr_tran_f
!
   cr_tran_f(:,:) = this%cr_tran_f(:,:)
!
   END SUBROUTINE get_cryst_tran_f
!******************************************************************************
   SUBROUTINE get_cryst_scat ( this, inum, itype, at_name, dw1, occ1)
!
!  Get the atom name and Debye Waller factor for atom Nr. inum
!
   IMPLICIT none
!
   CLASS (cl_cryst)                  :: this
   INTEGER,              INTENT(IN)  :: inum
   INTEGER,              INTENT(OUT) :: itype
   CHARACTER (LEN=4),    INTENT(OUT) :: at_name
   REAL   ,              INTENT(OUT) :: dw1
   REAL   ,              INTENT(OUT) :: occ1   !!WORK OCC
!
   REAL   , DIMENSION(3)             :: posit
   INTEGER                           :: iprop
   INTEGER, DIMENSION(0:3)           :: isurface
   INTEGER, DIMENSION(1:2)           :: iin_mole
!
   CALL this%atoms(inum)%get_atom ( itype, posit, iprop, isurface, iin_mole )
   at_name = this%cr_at_lis(itype)
   dw1     = this%cr_dw    (itype)
   occ1    = this%cr_occ   (itype)   !!WORK OCC
!
   END SUBROUTINE get_cryst_scat 
!******************************************************************************
   SUBROUTINE get_cryst_atom ( this, inum, itype, posit, iprop, isurface, iin_mole)
!
!  Return everything about atom <inum> in "this" crystal
!
   IMPLICIT none
!
   CLASS (cl_cryst)                  :: this
   INTEGER,              INTENT(IN ) :: inum
   INTEGER,              INTENT(OUT) :: itype
   REAL   , DIMENSION(3),INTENT(OUT) :: posit
   INTEGER,              INTENT(OUT) :: iprop
   INTEGER, DIMENSION(0:3),INTENT(OUT) :: isurface
   INTEGER, DIMENSION(1:2),INTENT(OUT) :: iin_mole
!
   CALL this%atoms(inum)%get_atom ( itype, posit, iprop, isurface, iin_mole )
!
   END SUBROUTINE get_cryst_atom 
!******************************************************************************
   SUBROUTINE get_cryst_mole ( this, inum, i_m_mole, i_m_type, i_m_char, &
                               c_m_file, r_m_fuzzy, r_m_dens, r_m_biso)
!
!  Get the atom name and Debye Waller factor for atom Nr. inum
!
   IMPLICIT none
!
   CLASS (cl_cryst)                  :: this
   INTEGER,              INTENT(IN)  :: inum
   INTEGER,              INTENT(OUT) :: i_m_mole
   INTEGER,              INTENT(OUT) :: i_m_type
   INTEGER,              INTENT(OUT) :: i_m_char
   CHARACTER (LEN=200),  INTENT(OUT) :: c_m_file
   REAL   ,              INTENT(OUT) :: r_m_fuzzy
   REAL   ,              INTENT(OUT) :: r_m_dens
   REAL   ,              INTENT(OUT) :: r_m_biso
!
   INTEGER                           :: i,j  ! Dummy
!
   i_m_mole = 0
   i_m_type = 0
   DO i=1, this%cr_num_mole
      DO j=1,this%cr_mole_len(i)
         IF ( this%cr_mole_cont(this%cr_mole_off(i)+j) == inum ) THEN
            i_m_mole = i
            i_m_type = this%cr_mole_type(i)
            i_m_char = this%cr_mole_char(i)
            c_m_file = this%cr_mole_file(i)
            r_m_fuzzy= this%cr_mole_fuzzy(i)
            r_m_dens = this%cr_mole_dens(i)
            r_m_biso = this%cr_mole_biso(this%cr_mole_type(i))
            RETURN
         ENDIF
      ENDDO
   ENDDO
!
   END SUBROUTINE get_cryst_mole 
!******************************************************************************
   SUBROUTINE set_crystal_from_standard   ( this, strucfile)
!
!  Set all values for the crystal to those of the main DISCUS crystal
!
   USE crystal_mod
   USE molecule_mod
   USE discus_save_mod
   USE gen_add_mod
   USE sym_add_mod
   USE modify_func_mod

   IMPLICIT none
!
   CLASS (cl_cryst)                 :: this        ! Work on "this" crystal
   CHARACTER (LEN=*), INTENT(IN)    :: strucfile
!
   INTEGER, DIMENSION(:), ALLOCATABLE :: iscat_table
   INTEGER               :: istatus
   INTEGER               :: inum
   INTEGER               :: itype
   REAL   , DIMENSION(3) :: posit
   INTEGER, DIMENSION(0:3) :: isurface
   INTEGER, DIMENSION(1:2) :: iin_mole
   INTEGER               :: iprop
   INTEGER               :: i,j
   INTEGER               :: ia
   
!
   this%cr_file         = strucfile
   this%cr_name         = cr_name
   this%cr_spcgr        = cr_spcgr
   this%cr_spcgrno      = cr_spcgrno
   this%cr_syst         = cr_syst
   this%cr_acentric     = cr_acentric
   IF(this%cr_sav_ncell) THEN
      this%cr_ncatoms      = cr_ncatoms
      this%cr_icc          = cr_icc
   ENDIF
!
   this%spcgr_ianz      = spcgr_ianz
   this%spcgr_para      = spcgr_para
   this%cr_sel_prop     = cr_sel_prop
!
   this%cr_a0           = cr_a0
   this%cr_win          = cr_win
   this%cr_dim          = cr_dim
   this%cr_dim0         = cr_dim0
   this%cr_v            = cr_v
   this%cr_gten         = cr_gten
   this%cr_eps          = cr_eps
   this%cr_gmat         = cr_gmat
   this%cr_ar           = cr_ar
   this%cr_wrez         = cr_wrez
   this%cr_vr           = cr_vr
   this%cr_rten         = cr_rten
   this%cr_reps         = cr_reps
   this%cr_fmat         = cr_fmat
   this%cr_tran_g       = cr_tran_g
   this%cr_tran_gi      = cr_tran_gi
   this%cr_tran_f       = cr_tran_f
   this%cr_tran_fi      = cr_tran_fi
!
   IF(this%cr_sav_gene) THEN
      this%cr_gen_add_n    = gen_add_n
      this%cr_gen_add_power= gen_add_power
      this%cr_gen_add      = gen_add
   ELSE
      this%cr_gen_add_n    = 0
      this%cr_gen_add_power= 1
      this%cr_gen_add      = 0.0
   ENDIF
   IF(this%cr_sav_symm) THEN
      this%cr_sym_add_n    = sym_add_n
      this%cr_sym_add_power= sym_add_power
      this%cr_sym_add      = sym_add
   ELSE
      this%cr_sym_add_n    = 0
      this%cr_sym_add_power= 1
      this%cr_sym_add      = 0.0
   ENDIF
!
!  Save scattering curves and names if "WRITE" flag was set
!
   ALLOCATE(iscat_table(0:MAXSCAT), STAT = istatus)
   iscat_table = 0
   IF(this%cr_sav_scat) THEN
      ia = -1
      DO i=0,cr_nscat
!        IF ( this%cr_sav_atom(i)) THEN
            ia = ia + 1
            iscat_table(i) = ia
            DO j=1,11
               this%cr_scat(j,ia)     = cr_scat(j,ia)
            ENDDO
!
            this%cr_delfr(ia)     = cr_delfr(i)
            this%cr_delfi(ia)     = cr_delfi(i)
!
            this%cr_scat_int(ia)  = cr_scat_int(i)
            this%cr_scat_equ(ia)  = cr_scat_equ(i)
            this%cr_delf_int(ia)  = cr_delf_int(i)
            this%cr_at_equ(ia)    = cr_at_equ(i)
            this%cr_at_lis(ia)    = cr_at_lis(i)
!        ENDIF
      ENDDO
   ELSE
      this%cr_scat         = 0.0
      this%cr_delfr        = 0.0
      this%cr_delfi        = 0.0
!
      this%cr_scat_int     = .true.
      this%cr_scat_equ     = .false.
      this%cr_delf_int     = .true.
      this%cr_at_equ       = ' '
      this%cr_at_lis(0)    = cr_at_lis(0)   ! always save void name
!     Always save atom names
      ia = 0
      DO i=1,cr_nscat
         IF ( this%cr_sav_atom(i)) THEN
            ia = ia + 1
            iscat_table(i) = ia
            this%cr_at_lis(ia)    = cr_at_lis(i)
         ENDIF
      ENDDO
   ENDIF
!  Always save ADP values 
!
   ia =  0
   this%cr_dw(0)       = cr_dw(0)           ! always save void ADP
   IF(this%cr_sav_adp) THEN
      DO i=1,cr_nscat
!        IF ( this%cr_sav_atom(i)) THEN
            ia = ia + 1
            iscat_table(i) = ia
            this%cr_dw(ia)       = cr_dw(ia)
!        ENDIF
      ENDDO
   ELSE
      DO i=1,cr_nscat
         IF ( this%cr_sav_atom(i)) THEN
            ia = ia + 1
            iscat_table(i) = ia
            this%cr_dw(ia)       = cr_dw(ia)
         ENDIF
      ENDDO
   ENDIF
!  save Occupancy values 
!
   ia =  0
   this%cr_occ(0)      = cr_occ(0)          ! always save void OCC
   IF(this%cr_sav_occ) THEN
      DO i=1,cr_nscat
!        IF ( this%cr_sav_atom(i)) THEN
            ia = ia + 1
            iscat_table(i) = ia
            this%cr_occ(ia)      = cr_occ(ia)
!        ENDIF
      ENDDO
   ELSE
      DO i=1,cr_nscat
         IF ( this%cr_sav_atom(i)) THEN
            ia = ia + 1
            iscat_table(i) = ia
            this%cr_occ(ia)      = 1.00    ! Use default instead
         ENDIF
      ENDDO
   ENDIF
   this%cr_nscat = ia
!
   this%cr_newtype      = cr_newtype
   this%cr_cartesian    = cr_cartesian
!
   this%cr_natoms       = cr_natoms
   this%cr_n_REAL_atoms = cr_n_REAL_atoms
!
! Save molecules if selected
!
   IF(this%cr_sav_mole .or. this%cr_sav_doma .or. this%cr_sav_obje) THEN
      this%cr_mole_max_mole = mole_max_mole
      this%cr_mole_max_type = mole_max_type
      this%cr_mole_max_atom = mole_max_atom
      this%cr_num_mole      = mole_num_mole
      this%cr_num_type      = mole_num_type
      this%cr_num_atom      = mole_num_atom
      FORALL ( i=0:mole_num_type)
         this%cr_mole_biso (i) = mole_biso (i)
      END FORALL
      FORALL ( i=0:mole_num_mole)
         this%cr_mole_len  (i) = mole_len  (i)
         this%cr_mole_off  (i) = mole_off  (i)
         this%cr_mole_type (i) = mole_type (i)
         this%cr_mole_char (i) = mole_char (i)
         this%cr_mole_file (i) = mole_file (i)
         this%cr_mole_dens (i) = mole_dens (i)
         this%cr_mole_fuzzy(i) = mole_fuzzy(i)
      END FORALL
      FORALL ( i=0:mole_num_atom)
         this%cr_mole_cont (i) = mole_cont (i)
      END FORALL
   ELSE
      this%cr_mole_max_mole = 0
      this%cr_mole_max_type = 0
      this%cr_mole_max_atom = 0
      this%cr_num_mole      = 0
      this%cr_num_type      = 0
      this%cr_num_atom      = 0
      this%cr_mole_len      = 0
      this%cr_mole_off      = 0
      this%cr_mole_type     = 0
      this%cr_mole_char     = 0
      this%cr_mole_file     = ' '
      this%cr_mole_dens     = 0.0
      this%cr_mole_biso     = 0.0
      this%cr_mole_fuzzy    = 0.0
      this%cr_mole_cont     = 0
   ENDIF
!
!  Save atoms, if selected
!
   ia = 0
   DO inum=1,cr_natoms
!     IF(this%cr_sav_atom(cr_iscat(inum))) THEN
      IF(check_select_status(inum,this%cr_sav_atom(cr_iscat(inum)),   &
                                  cr_prop (inum),                &
                             this%cr_sav_sel_prop) ) THEN
         ia = ia + 1
         itype = iscat_table(cr_iscat(inum))
         posit = cr_pos(:,inum)
         isurface(:) = cr_surf(:, inum)
         IF(this%cr_sav_prop) THEN
            iprop = cr_prop (inum)
         ELSE
            iprop = 1
         ENDIF
         IF(this%cr_sav_mole .OR. this%cr_sav_doma .OR. this%cr_sav_obje) THEN
            IF(cr_mole(inum)/=0) THEN
               iin_mole(1) = cr_mole(inum)
               check_mole: DO j = 1, mole_len (cr_mole(inum))
                  IF(mole_cont (mole_off(cr_mole(inum))+j) == inum) THEN
                     iin_mole(2) = j
                     this%cr_mole_cont (mole_off(cr_mole(inum))+j) = ia
                     EXIT check_mole
                  ENDIF
               ENDDO check_mole
            ELSE
               iin_mole(:) = 0
            ENDIF
         ENDIF
         CALL this%atoms(ia)%set_atom ( itype, posit, iprop, isurface, iin_mole )
      ELSE
         IF(this%cr_sav_mole .OR. this%cr_sav_doma .OR. this%cr_sav_obje) THEN
            IF(cr_mole(inum)/=0) THEN
               ier_num = -157
               ier_typ = ER_APPL
               WRITE(ier_msg(1), '(a,i7)') 'Atom number ',inum
               ier_msg(2) = 'The atom is deseclected but part of a molecule'
               ier_msg(3) = 'Remove this atom type and purge prior to save '
            ENDIF
         ENDIF
      ENDIF
   ENDDO
   this%cr_natoms = MIN(cr_natoms,ia)
   DEALLOCATE(iscat_table, STAT=istatus)

!
   END SUBROUTINE set_crystal_from_standard
!******************************************************************************
   SUBROUTINE set_crystal_from_local   ( this, strucfile, &
            rd_NMAX, rd_MAXSCAT, rd_n_mole, rd_n_mole_type, rd_n_atom, rd_cr_name,       &
            rd_cr_natoms, rd_cr_ncatoms, rd_cr_n_REAL_atoms, rd_cr_spcgrno, rd_cr_syst, &
            rd_cr_spcgr, rd_cr_at_lis, rd_cr_at_equ, rd_cr_as_lis,                      &
            rd_cr_nscat, rd_cr_dw, rd_cr_occ, rd_cr_a0, rd_cr_win,                      &
            rd_cr_ar, rd_cr_wrez, rd_cr_v, rd_cr_vr, rd_cr_dim, rd_cr_dim0, rd_cr_icc,  &
            rd_sav_ncell, rd_sav_r_ncell, rd_sav_ncatoms, rd_spcgr_ianz, rd_spcgr_para, &
            rd_cr_tran_g, rd_cr_tran_gi, rd_cr_tran_f, rd_cr_tran_fi, rd_cr_gmat, rd_cr_fmat, &
            rd_cr_gten, rd_cr_rten, rd_cr_eps, rd_cr_reps,                              &
            rd_cr_acentric, rd_cr_newtype, rd_cr_cartesian, rd_cr_sel_prop,             &
            rd_cr_scat, rd_cr_delfi , rd_cr_delfr, rd_cr_delf_int,                    &
            rd_cr_scat_int, rd_cr_scat_equ,                                             &
            rd_cr_pos, rd_cr_iscat, rd_cr_prop, rd_cr_surf                              &
            )
!
!
!  Set all values for the crystal to those of the local copy
!
!  USE crystal_mod
   USE molecule_mod
   USE discus_save_mod
   USE gen_add_mod
   USE sym_add_mod
   USE modify_func_mod

   IMPLICIT none
!
   CLASS (cl_cryst)                 :: this        ! Work on "this" crystal
   CHARACTER (LEN=*), INTENT(IN)    :: strucfile
!
   INTEGER                                      , INTENT(IN) :: rd_NMAX
   INTEGER                                      , INTENT(IN) :: rd_MAXSCAT 
!
   CHARACTER (LEN=  80)                         , INTENT(IN) :: rd_cr_name 
   INTEGER                                      , INTENT(IN) :: rd_cr_natoms 
   INTEGER                                      , INTENT(IN) :: rd_cr_ncatoms 
   INTEGER                                      , INTENT(IN) :: rd_n_mole      ! Molecule number in this crystal
   INTEGER                                      , INTENT(IN) :: rd_n_mole_type ! Molecule number in this crystal
   INTEGER                                      , INTENT(IN) :: rd_n_atom      ! Atoms in molecule in this crystal 
   INTEGER                                      , INTENT(IN) :: rd_cr_n_REAL_atoms 
   INTEGER                                      , INTENT(IN) :: rd_cr_spcgrno 
   CHARACTER (LEN=  16)                         , INTENT(IN) :: rd_cr_spcgr 
   INTEGER                                      , INTENT(IN) :: rd_cr_syst 
   REAL                , DIMENSION(3)           , INTENT(IN) :: rd_cr_a0
   REAL                , DIMENSION(3)           , INTENT(IN) :: rd_cr_win
   REAL                , DIMENSION(3)           , INTENT(IN) :: rd_cr_ar
   REAL                , DIMENSION(3)           , INTENT(IN) :: rd_cr_wrez
   REAL                                         , INTENT(IN) :: rd_cr_v
   REAL                                         , INTENT(IN) :: rd_cr_vr
   REAL                , DIMENSION(3,2)         , INTENT(IN) :: rd_cr_dim
   REAL                , DIMENSION(3,2)         , INTENT(IN) :: rd_cr_dim0
   INTEGER             , DIMENSION(3)           , INTENT(IN) :: rd_cr_icc
   INTEGER                                      , INTENT(IN) :: rd_cr_nscat 
   REAL                , DIMENSION(0:rd_MAXSCAT), INTENT(IN) :: rd_cr_dw     ! (0:MAXSCAT) 
   REAL                , DIMENSION(0:rd_MAXSCAT), INTENT(IN) :: rd_cr_occ    ! (0:MAXSCAT)   !! WORK OCC
   CHARACTER (LEN=   4), DIMENSION(0:rd_MAXSCAT), INTENT(IN) :: rd_cr_at_lis ! (0:MAXSCAT) 
   CHARACTER (LEN=   4), DIMENSION(0:rd_MAXSCAT), INTENT(IN) :: rd_cr_at_equ ! (0:MAXSCAT) 
   CHARACTER (LEN=   4), DIMENSION(0:rd_MAXSCAT), INTENT(IN) :: rd_cr_as_lis ! (0:MAXSCAT) 
   INTEGER             , DIMENSION(3)           , INTENT(IN) :: rd_sav_ncell ! (3) 
   LOGICAL                                      , INTENT(IN) :: rd_sav_r_ncell 
   INTEGER                                      , INTENT(IN) :: rd_sav_ncatoms 
   INTEGER                                      , INTENT(IN) :: rd_spcgr_ianz 
   INTEGER                                      , INTENT(IN) :: rd_spcgr_para 
!  INTEGER                                      , INTENT(IN) :: rd_GEN_ADD_MAX
!  INTEGER                                      , INTENT(IN) :: rd_GEN_ADD_n
!  INTEGER           , DIMENSION(rd_GEN_ADD_MAX), INTENT(IN) :: rd_gen_add_power
!  REAL        , DIMENSION(4,4,0:rd_GEN_ADD_MAX), INTENT(IN) :: rd_gen_add
!  INTEGER                                      , INTENT(IN) :: rd_SYM_ADD_MAX
!  INTEGER                                      , INTENT(IN) :: rd_SYM_ADD_n
!  INTEGER           , DIMENSION(rd_SYM_ADD_MAX), INTENT(IN) :: rd_sym_add_power
!  REAL        , DIMENSION(4,4,0:rd_SYM_ADD_MAX), INTENT(IN) :: rd_sym_add
   REAL                , DIMENSION(3,3)         , INTENT(IN) :: rd_cr_gten
   REAL                , DIMENSION(3,3)         , INTENT(IN) :: rd_cr_rten
   REAL                , DIMENSION(3,3,3)       , INTENT(IN) :: rd_cr_eps
   REAL                , DIMENSION(3,3,3)       , INTENT(IN) :: rd_cr_reps
   REAL                , DIMENSION(3,3)         , INTENT(IN) :: rd_cr_gmat
   REAL                , DIMENSION(3,3)         , INTENT(IN) :: rd_cr_fmat
   REAL                , DIMENSION(4,4)         , INTENT(IN) :: rd_cr_tran_g
   REAL                , DIMENSION(4,4)         , INTENT(IN) :: rd_cr_tran_gi
   REAL                , DIMENSION(4,4)         , INTENT(IN) :: rd_cr_tran_f
   REAL                , DIMENSION(4,4)         , INTENT(IN) :: rd_cr_tran_fi
   LOGICAL                                      , INTENT(IN) :: rd_cr_acentric
   LOGICAL                                      , INTENT(IN) :: rd_cr_newtype
   LOGICAL                                      , INTENT(IN) :: rd_cr_cartesian
   INTEGER             , DIMENSION(0:1)         , INTENT(IN) :: rd_cr_sel_prop
   REAL              ,DIMENSION(11,0:rd_MAXSCAT), INTENT(IN) :: rd_cr_scat   ! (11,0:MAXSCAT)
   REAL                , DIMENSION(0:rd_MAXSCAT), INTENT(IN) :: rd_cr_delfr  ! (  0:MAXSCAT)
   REAL                , DIMENSION(0:rd_MAXSCAT), INTENT(IN) :: rd_cr_delfi  ! (  0:MAXSCAT)
   LOGICAL             , DIMENSION(0:rd_MAXSCAT), INTENT(IN) :: rd_cr_scat_int  ! (  0:MAXSCAT)
   LOGICAL             , DIMENSION(0:rd_MAXSCAT), INTENT(IN) :: rd_cr_scat_equ  ! (  0:MAXSCAT)
   LOGICAL             , DIMENSION(0:rd_MAXSCAT), INTENT(IN) :: rd_cr_delf_int  ! (  0:MAXSCAT)

   REAL                , DIMENSION(3,rd_NMAX)   , INTENT(IN) :: rd_cr_pos
   INTEGER             , DIMENSION(  rd_NMAX)   , INTENT(IN) :: rd_cr_iscat
   INTEGER             , DIMENSION(  rd_NMAX)   , INTENT(IN) :: rd_cr_prop
   INTEGER             , DIMENSION(0:3,rd_NMAX) , INTENT(IN) :: rd_cr_surf
!
   INTEGER, DIMENSION(:), ALLOCATABLE :: iscat_table
   INTEGER               :: istatus
   INTEGER               :: inum
   INTEGER               :: itype
   REAL   , DIMENSION(3) :: posit
   INTEGER, DIMENSION(0:3) :: isurface
   INTEGER, DIMENSION(1:2) :: iin_mole
   INTEGER               :: iprop
   INTEGER               :: i,j
   INTEGER               :: ia
   
!
   this%cr_file         = strucfile
   this%cr_name         = rd_cr_name
   this%cr_spcgr        = rd_cr_spcgr
   this%cr_spcgrno      = rd_cr_spcgrno
   this%cr_syst         = rd_cr_syst
   this%cr_acentric     = rd_cr_acentric
   IF(this%cr_sav_ncell) THEN
      this%cr_ncatoms      = rd_cr_ncatoms
      this%cr_icc          = rd_cr_icc
   ENDIF
!
   this%spcgr_ianz      = rd_spcgr_ianz
   this%spcgr_para      = rd_spcgr_para
   this%cr_sel_prop     = rd_cr_sel_prop
!
   this%cr_a0           = rd_cr_a0
   this%cr_win          = rd_cr_win
   this%cr_dim          = rd_cr_dim
   this%cr_dim0         = rd_cr_dim0
   this%cr_v            = rd_cr_v
   this%cr_gten         = rd_cr_gten
   this%cr_eps          = rd_cr_eps
   this%cr_gmat         = rd_cr_gmat
   this%cr_ar           = rd_cr_ar
   this%cr_wrez         = rd_cr_wrez
   this%cr_vr           = rd_cr_vr
   this%cr_rten         = rd_cr_rten
   this%cr_reps         = rd_cr_reps
   this%cr_fmat         = rd_cr_fmat
   this%cr_tran_g       = rd_cr_tran_g
   this%cr_tran_gi      = rd_cr_tran_gi
   this%cr_tran_f       = rd_cr_tran_f
   this%cr_tran_fi      = rd_cr_tran_fi
!
   IF(this%cr_sav_gene) THEN
      this%cr_gen_add_n    = gen_add_n
      this%cr_gen_add_power= gen_add_power
      this%cr_gen_add      = gen_add
   ELSE
      this%cr_gen_add_n    = 0
      this%cr_gen_add_power= 1
      this%cr_gen_add      = 0.0
   ENDIF
   IF(this%cr_sav_symm) THEN
      this%cr_sym_add_n    = sym_add_n
      this%cr_sym_add_power= sym_add_power
      this%cr_sym_add      = sym_add
   ELSE
      this%cr_sym_add_n    = 0
      this%cr_sym_add_power= 1
      this%cr_sym_add      = 0.0
   ENDIF
!
!  Save scattering curves and names if "WRITE" flag was set
!
   ALLOCATE(iscat_table(0:rd_MAXSCAT), STAT = istatus)
   iscat_table(:) = 0
   IF(this%cr_sav_scat) THEN
      ia = -1
      DO i=0,rd_cr_nscat
!        IF ( this%cr_sav_atom(i)) THEN
            ia = ia + 1
            iscat_table(i) = ia
            DO j=1,11
               this%cr_scat(j,ia)     = rd_cr_scat(j,i)
            ENDDO
!
            this%cr_delfr(ia)     = rd_cr_delfr(i)
            this%cr_delfi(ia)     = rd_cr_delfi(i)
!
            this%cr_scat_int(ia)  = rd_cr_scat_int(i)
            this%cr_scat_equ(ia)  = rd_cr_scat_equ(i)
            this%cr_delf_int(ia)  = rd_cr_delf_int(i)
            this%cr_at_equ(ia)    = rd_cr_at_equ(i)
            this%cr_at_lis(ia)    = rd_cr_at_lis(i)
!        ENDIF
      ENDDO
   ELSE
      this%cr_scat         = 0.0
      this%cr_delfr        = 0.0
      this%cr_delfi        = 0.0
!
      this%cr_scat_int     = .true.
      this%cr_scat_equ     = .false.
      this%cr_delf_int     = .true.
      this%cr_at_equ       = ' '
      this%cr_at_lis(0)    = rd_cr_at_lis(0)   ! always save void name
!     Always save atom names
      ia = 0
      DO i=1,rd_cr_nscat
         IF ( this%cr_sav_atom(i)) THEN
            ia = ia + 1
            iscat_table(i) = ia
            this%cr_at_lis(ia)    = rd_cr_at_lis(i)
         ENDIF
      ENDDO
   ENDIF
!  Always save ADP values 
!
   ia =  0
   this%cr_dw(0)       = rd_cr_dw(0)           ! always save void ADP
   IF(this%cr_sav_adp) THEN
      DO i=1,rd_cr_nscat
!        IF ( this%cr_sav_atom(i)) THEN
            ia = ia + 1
            iscat_table(i) = ia
            this%cr_dw(ia)       = rd_cr_dw(ia)
!        ENDIF
      ENDDO
   ELSE
      DO i=1,rd_cr_nscat
         IF ( this%cr_sav_atom(i)) THEN
            ia = ia + 1
            iscat_table(i) = ia
            this%cr_dw(ia)       = rd_cr_dw(ia)
         ENDIF
      ENDDO
   ENDIF
   this%cr_nscat = ia
!
!  save occupancy values 
!
   ia =  0
   this%cr_occ(0)       = rd_cr_occ(0)           ! always save void ADP
   IF(this%cr_sav_occ) THEN
      DO i=1,rd_cr_nscat
!        IF ( this%cr_sav_atom(i)) THEN
            ia = ia + 1
            iscat_table(i) = ia
            this%cr_occ(ia)       = rd_cr_occ(ia)
!        ENDIF
      ENDDO
   ELSE
      DO i=1,rd_cr_nscat
         IF ( this%cr_sav_atom(i)) THEN
            ia = ia + 1
            iscat_table(i) = ia
            this%cr_occ(ia)       = 1.00
         ENDIF
      ENDDO
   ENDIF
   this%cr_nscat = ia
!
   this%cr_newtype      = rd_cr_newtype
   this%cr_cartesian    = rd_cr_cartesian
!
   this%cr_natoms       = rd_cr_natoms
   this%cr_n_REAL_atoms = rd_cr_n_REAL_atoms
!
!  Save atoms, if selected
!
   ia = 0
   DO inum=1,rd_cr_natoms
!     IF(this%cr_sav_atom(cr_iscat(inum))) THEN
      IF(check_select_status(inum,this%cr_sav_atom(rd_cr_iscat(inum)),   &
                                  rd_cr_prop (inum),                &
                             this%cr_sav_sel_prop) ) THEN
         ia = ia + 1
         itype = iscat_table(rd_cr_iscat(inum))
         posit = rd_cr_pos(:,inum)
         isurface(:) = rd_cr_surf(:,inum)
         IF(this%cr_sav_prop) THEN
            iprop = rd_cr_prop (inum)
         ELSE
            iprop = 1
         ENDIF
         CALL this%atoms(ia)%set_atom ( itype, posit, iprop, isurface, iin_mole )
      ENDIF
   ENDDO
   this%cr_natoms = MIN(rd_cr_natoms,ia)
   DEALLOCATE(iscat_table, STAT=istatus)

!
   IF(this%cr_sav_mole .or. this%cr_sav_doma .or. this%cr_sav_obje .AND. rd_n_mole > 0) THEN
      this%cr_mole_max_mole = rd_n_mole
      this%cr_mole_max_type = rd_n_mole_type
      this%cr_mole_max_atom = mole_max_atom
      this%cr_num_mole      = rd_n_mole
      this%cr_num_type      = rd_n_mole_type
      this%cr_num_atom      = rd_n_atom
      FORALL ( i=0:rd_n_mole_type)
         this%cr_mole_biso (i) = mole_biso (i)
      END FORALL
      FORALL ( i=0:rd_n_mole)
         this%cr_mole_len  (i) = mole_len  (i)
         this%cr_mole_off  (i) = mole_off  (i)
         this%cr_mole_type (i) = mole_type (i)
         this%cr_mole_char (i) = mole_char (i)
         this%cr_mole_file (i) = mole_file (i)
         this%cr_mole_dens (i) = mole_dens (i)
         this%cr_mole_fuzzy(i) = mole_fuzzy(i)
      END FORALL
      FORALL ( i=0:rd_n_atom)
         this%cr_mole_cont (i) = mole_cont (i)
      END FORALL
   ELSE
      this%cr_mole_max_mole = 0
      this%cr_mole_max_type = 0
      this%cr_mole_max_atom = 0
      this%cr_num_mole      = 0
      this%cr_num_type      = 0
      this%cr_num_atom      = 0
      this%cr_mole_len      = 0
      this%cr_mole_off      = 0
      this%cr_mole_type     = 0
      this%cr_mole_char     = 0
      this%cr_mole_file     = ' '
      this%cr_mole_dens     = 0.0
      this%cr_mole_biso     = 0.0
      this%cr_mole_fuzzy    = 0.0
      this%cr_mole_cont     = 0
   ENDIF
!
   END SUBROUTINE set_crystal_from_local
!******************************************************************************
   SUBROUTINE set_crystal_save_flags ( this,  sav_scat, sav_adp, sav_occ, &
                                       sav_gene, sav_symm, &
                                       sav_w_ncell, sav_obje, sav_doma, sav_mole, &
                                       sav_prop,sav_sel_prop,n_latom,sav_latom)
!
   USE crystal_mod
   IMPLICIT none
!
   CLASS (cl_cryst)                 :: this        ! Work on "this" crystal
!
!   INTEGER, INTENT(IN)              ::   n_nscat
   LOGICAL, INTENT(IN)              ::  sav_scat
   LOGICAL, INTENT(IN)              ::  sav_adp
   LOGICAL, INTENT(IN)              ::  sav_occ
   LOGICAL, INTENT(IN)              ::  sav_gene
   LOGICAL, INTENT(IN)              ::  sav_symm
   LOGICAL, INTENT(IN)              ::  sav_w_ncell
   LOGICAL, INTENT(IN)              ::  sav_obje
   LOGICAL, INTENT(IN)              ::  sav_doma
   LOGICAL, INTENT(IN)              ::  sav_mole
   LOGICAL, INTENT(IN)              ::  sav_prop
   INTEGER, DIMENSION(0:1),INTENT(IN) ::  sav_sel_prop
   INTEGER, INTENT(IN)              ::  n_latom
   LOGICAL, DIMENSION(0:n_latom),INTENT(IN) ::  sav_latom
!
   INTEGER                          :: i
!
   this%cr_sav_scat  = sav_scat
   this%cr_sav_adp   = sav_adp
   this%cr_sav_occ   = sav_occ
   this%cr_sav_gene  = sav_gene
   this%cr_sav_symm  = sav_symm
   this%cr_sav_ncell = sav_w_ncell
   this%cr_sav_obje  = sav_obje
   this%cr_sav_doma  = sav_doma
   this%cr_sav_mole  = sav_mole
   this%cr_sav_prop  = sav_prop
   this%cr_sav_sel_prop(:)  = sav_sel_prop(:)
   DO i = 0, n_latom
      this%cr_sav_atom(i) = sav_latom(i)
   ENDDO
!
   END SUBROUTINE set_crystal_save_flags
!******************************************************************************
   SUBROUTINE get_header_from_crystal   ( this)
!
!  Get all values from "this" crystal and place into the main DISCUS crystal!
!  The standard crystal is completely overwritten!
!
   USE crystal_mod
   USE molecule_mod
   USE discus_save_mod    ! should be a read_mod !?!?
   USE gen_add_mod
   USE sym_add_mod
   IMPLICIT none
!
   CLASS (cl_cryst)                 :: this        ! Work on "this" crystal
!
   INTEGER               :: i,j
   INTEGER               :: ia
!
   cr_name         = this%cr_name
   cr_spcgr        = this%cr_spcgr
   cr_spcgrno      = this%cr_spcgrno
   cr_syst         = this%cr_syst
   cr_acentric     = this%cr_acentric
!
   spcgr_ianz      = this%spcgr_ianz
   spcgr_para      = this%spcgr_para
   cr_sel_prop     = this%cr_sel_prop
!
   cr_a0           = this%cr_a0
   cr_win          = this%cr_win
   cr_dim          = this%cr_dim
   cr_dim0         = this%cr_dim0
   cr_v            = this%cr_v
   cr_gten         = this%cr_gten
   cr_eps          = this%cr_eps
   cr_gmat         = this%cr_gmat
   cr_ar           = this%cr_ar
   cr_wrez         = this%cr_wrez
   cr_vr           = this%cr_vr
   cr_rten         = this%cr_rten
   cr_reps         = this%cr_reps
   cr_fmat         = this%cr_fmat
   cr_tran_g       = this%cr_tran_g
   cr_tran_gi      = this%cr_tran_gi
   cr_tran_f       = this%cr_tran_f
   cr_tran_fi      = this%cr_tran_fi
!
   cr_natoms       = this%cr_natoms
!
   IF ( sav_w_ncell ) THEN
      cr_ncatoms   = this%cr_ncatoms
      cr_icc       = this%cr_icc
   ELSE
!                                                                       
!     Define initial crystal size in number of unit cells         
!                                                                       
      DO ia = 1, 3
         cr_icc(ia) = MAX (1,INT(cr_dim0(ia,2) - cr_dim0 (ia,1) ) )
      ENDDO
!                                                                       
!     Define (average) number of atoms per unit cell              
!                                                                       
      cr_ncatoms = cr_natoms / (cr_icc (1) * cr_icc (2) *   &
                                cr_icc (3) )
   ENDIF
!
   IF ( sav_w_gene ) THEN
      gen_add_n    = this%cr_gen_add_n
      gen_add_power= this%cr_gen_add_power
      gen_add      = this%cr_gen_add
   ELSE
      gen_add_n    = 0
      gen_add_power= 1
      gen_add      = 0.0
   ENDIF
   IF ( sav_w_symm ) THEN
      sym_add_n    = this%cr_sym_add_n
      sym_add_power= this%cr_sym_add_power
      sym_add      = this%cr_sym_add
   ELSE
      sym_add_n    = 0
      sym_add_power= 1
      sym_add      = 0.0
   ENDIF
!
   FORALL (j=1:11,i=0:this%cr_nscat)
      cr_scat(j,i)     = this%cr_scat(j,i)
   END FORALL
!
   FORALL (i=0:this%cr_nscat)
      cr_delfr(i)     = this%cr_delfr(i)
      cr_delfi(i)     = this%cr_delfi(i)
!
      cr_scat_int(i)  = this%cr_scat_int(i)
      cr_scat_equ(i)  = this%cr_scat_equ(i)
      cr_delf_int(i)  = this%cr_delf_int(i)
      cr_at_equ(i)    = this%cr_at_equ(i)
      cr_at_lis(i)    = this%cr_at_lis(i)
      cr_dw(i)        = this%cr_dw(i)
      cr_occ(i)       = this%cr_occ(i)
   END FORALL
!
   cr_newtype      = this%cr_newtype
   cr_cartesian    = this%cr_cartesian
!
   cr_nscat        = this%cr_nscat
   cr_n_REAL_atoms = this%cr_n_REAL_atoms
!
   END SUBROUTINE get_header_from_crystal
!******************************************************************************
   SUBROUTINE get_header_to_local (this, rd_MAXSCAT, rd_cr_name,      &
            rd_cr_spcgr, rd_cr_at_lis, rd_cr_nscat, rd_cr_dw, rd_cr_occ, rd_cr_a0, rd_cr_win,      &
            rd_sav_ncell, rd_sav_r_ncell, rd_sav_ncatoms, rd_spcgr_ianz, rd_spcgr_para, &
            rd_GEN_ADD_MAX, rd_gen_add_n, rd_gen_add_power, rd_gen_add,                 &
            rd_SYM_ADD_MAX, rd_sym_add_n, rd_sym_add_power, rd_sym_add )
!
   IMPLICIT NONE
!
   INTEGER                                      , INTENT(INOUT) :: rd_MAXSCAT 
!
   CHARACTER (LEN=  80)                         , INTENT(INOUT) :: rd_cr_name 
   CHARACTER (LEN=  16)                         , INTENT(INOUT) :: rd_cr_spcgr 
   REAL                , DIMENSION(3)           , INTENT(INOUT) :: rd_cr_a0
   REAL                , DIMENSION(3)           , INTENT(INOUT) :: rd_cr_win
   INTEGER                                      , INTENT(INOUT) :: rd_cr_nscat 
   REAL                , DIMENSION(0:rd_MAXSCAT), INTENT(INOUT) :: rd_cr_dw     ! (0:MAXSCAT) 
   REAL                , DIMENSION(0:rd_MAXSCAT), INTENT(INOUT) :: rd_cr_occ    ! (0:MAXSCAT)   !! WORK OCC
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
!  Get header variables from an internal crystal to local variables
!
   CLASS (cl_cryst)                 :: this        ! Work on "this" crystal
!
   INTEGER                          :: i
!
   rd_cr_name         = this%cr_name
   rd_cr_spcgr        = this%cr_spcgr
!
   rd_spcgr_ianz      = this%spcgr_ianz
   rd_spcgr_para      = this%spcgr_para
!
   rd_cr_a0           = this%cr_a0
   rd_cr_win          = this%cr_win
!
   rd_sav_ncell       = this%cr_icc
   rd_sav_ncatoms     = this%cr_ncatoms
   rd_sav_r_ncell     = .true.
!
!
   FORALL (i=0:this%cr_nscat)
      rd_cr_at_lis(i)    = this%cr_at_lis(i)
      rd_cr_dw(i)        = this%cr_dw(i)
      rd_cr_occ(i)       = this%cr_occ(i)  !! WORK OCC
   END FORALL
!
   rd_cr_nscat        = MAX(rd_cr_nscat,this%cr_nscat)
!
   rd_gen_add_n    = this%cr_gen_add_n
   rd_gen_add_power= this%cr_gen_add_power
   rd_gen_add      = this%cr_gen_add
   rd_sym_add_n    = this%cr_sym_add_n
   rd_sym_add_power= this%cr_sym_add_power
   rd_sym_add      = this%cr_sym_add
!
   END SUBROUTINE get_header_to_local
!******************************************************************************
   SUBROUTINE get_atoms_from_crystal   ( this)
!
!  Get all values from "this" crystal and place into the main DISCUS crystal!
!  The standard crystal is completely overwritten!
!
   USE crystal_mod
   USE molecule_mod
   USE discus_save_mod    ! should be a read_mod !?!?
   USE gen_add_mod
   USE sym_add_mod
   IMPLICIT none
!
   CLASS (cl_cryst)                 :: this        ! Work on "this" crystal
!
   INTEGER               :: inum
   INTEGER               :: itype
   REAL   , DIMENSION(3) :: posit
   INTEGER, DIMENSION(0:3) :: isurface
   INTEGER, DIMENSION(1:2) :: iin_mole
   INTEGER               :: iprop
   INTEGER               :: ia
!
   ia = 0
   DO inum=1,this%cr_natoms
      CALL this%atoms(inum)%get_atom ( itype, posit, iprop, isurface, iin_mole )
      ia = ia + 1
      cr_iscat(ia) = itype
      cr_pos(:,ia) = posit
      cr_surf(:,ia) = isurface(:)
      cr_prop (ia) = iprop
      cr_mole (ia) = iin_mole(1)
   ENDDO
   cr_natoms = ia
!
   END SUBROUTINE get_atoms_from_crystal
!******************************************************************************
   SUBROUTINE get_atoms_to_local( this, RD_NMAX, rd_cr_natoms,  &
            rd_cr_pos, rd_cr_iscat, rd_cr_prop, rd_cr_surf, rd_cr_mole )
!
!  Get atom positions from internal crystal and copy them into
!  local variables.
!
   IMPLICIT none
!
   INTEGER,                          INTENT(IN)     :: RD_NMAX
   INTEGER,                          INTENT(INOUT)  :: rd_cr_natoms
   INTEGER, DIMENSION(1:RD_NMAX),    INTENT(INOUT)  :: rd_cr_iscat
   INTEGER, DIMENSION(1:RD_NMAX),    INTENT(INOUT)  :: rd_cr_prop
   REAL   , DIMENSION(1:3,1:RD_NMAX),INTENT(INOUT)  :: rd_cr_pos
   INTEGER, DIMENSION(0:3,1:RD_NMAX),INTENT(INOUT)  :: rd_cr_surf
   INTEGER, DIMENSION(    1:RD_NMAX),INTENT(INOUT)  :: rd_cr_mole
!                                                                       

   CLASS (cl_cryst)                 :: this        ! Work on "this" crystal
!
   INTEGER               :: inum
   INTEGER               :: itype
   REAL   , DIMENSION(3) :: posit
   INTEGER, DIMENSION(0:3) :: isurface
   INTEGER, DIMENSION(1:2) :: iin_mole
   INTEGER               :: iprop
   INTEGER               :: ia
!
   ia = 0
   DO inum=1,this%cr_natoms
      CALL this%atoms(inum)%get_atom ( itype, posit, iprop, isurface, iin_mole )
      IF ( this%cr_sav_atom(itype)) THEN
         ia = ia + 1
         rd_cr_iscat(ia) = itype
         rd_cr_pos(:,ia) = posit
         rd_cr_prop (ia) = iprop
         rd_cr_surf(:,ia)= isurface(:)
         rd_cr_mole(  ia)= iin_mole(1)
      ENDIF
   ENDDO
   rd_cr_natoms = ia
!
   END SUBROUTINE get_atoms_to_local
!******************************************************************************
   SUBROUTINE get_molecules_from_crystal   ( this, mole_max_mole,         &
              mole_max_type, mole_max_atom, mole_num_mole, mole_num_type, &
              mole_num_atom, mole_len, mole_off, mole_type, mole_char,    &
              mole_file, mole_dens, mole_biso, mole_fuzzy, mole_cont)
!
!  Get all values from "this" crystal and place into the main DISCUS crystal!
!  The standard crystal is completely overwritten!
!
   USE discus_save_mod    ! should be a read_mod !?!?
   IMPLICIT none
!
   CLASS (cl_cryst)                 :: this        ! Work on "this" crystal
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
   REAL   ,           DIMENSION(0:MOLE_MAX_MOLE), INTENT(OUT) :: mole_fuzzy
   INTEGER,           DIMENSION(0:MOLE_MAX_ATOM), INTENT(OUT) :: mole_cont
!
   INTEGER               :: ia
!
   mole_num_mole = this%cr_num_mole
   mole_num_type = this%cr_num_type
   mole_num_atom = this%cr_num_atom
   IF ( sav_w_mole .or. sav_w_obje .or. sav_w_doma ) THEN
      FORALL ( ia=0:mole_num_type)
         mole_biso (ia) = this%cr_mole_biso (ia)
      END FORALL
      FORALL ( ia=0:mole_num_mole)
         mole_len  (ia) = this%cr_mole_len  (ia)
         mole_off  (ia) = this%cr_mole_off  (ia)
         mole_type (ia) = this%cr_mole_type (ia)
         mole_char (ia) = this%cr_mole_char (ia)
         mole_file (ia) = this%cr_mole_file (ia)
         mole_dens (ia) = this%cr_mole_dens (ia)
         mole_fuzzy(ia) = this%cr_mole_fuzzy(ia)
      END FORALL
      FORALL ( ia=0:mole_num_atom)
         mole_cont (ia) = this%cr_mole_cont (ia)
      END FORALL
   ENDIF
!
   END SUBROUTINE get_molecules_from_crystal
!******************************************************************************
   LOGICAL FUNCTION is_alloc_atom (this)
!
   IMPLICIT NONE
!
   CLASS (cl_cryst)                 :: this        ! Work on "this" crystal
!
   is_alloc_atom = this%latom
!
   END FUNCTION is_alloc_atom 
!******************************************************************************
   SUBROUTINE finalize_atoms   ( this)
!
!  Finalize the class cl_cryst
!  Deallocate the array of atoms for "this" crystal 
!
   IMPLICIT none
!
   CLASS (cl_cryst)                 :: this        ! Work on "this" crystal
!
   INTEGER  :: istatus
!
   IF(this%latom) THEN
      DEALLOCATE ( this%cr_gen_add , STAT=istatus ) ! Allocate equivalent atom names
      DEALLOCATE ( this%cr_sym_add , STAT=istatus ) ! Allocate equivalent atom names
      DEALLOCATE ( this%cr_scat    , STAT=istatus ) ! Allocate equivalent atom names
      DEALLOCATE ( this%cr_delfr   , STAT=istatus ) ! Allocate equivalent atom names
      DEALLOCATE ( this%cr_delfi   , STAT=istatus ) ! Allocate equivalent atom names
      DEALLOCATE ( this%cr_scat_int, STAT=istatus ) ! Allocate equivalent atom names
      DEALLOCATE ( this%cr_scat_equ, STAT=istatus ) ! Allocate equivalent atom names
      DEALLOCATE ( this%cr_delf_int, STAT=istatus ) ! Allocate equivalent atom names
      DEALLOCATE ( this%cr_at_equ  , STAT=istatus ) ! Always deallocate
      DEALLOCATE ( this%cr_sav_atom, STAT=istatus ) ! Always deallocate
      DEALLOCATE ( this%cr_at_lis  , STAT=istatus ) ! Always deallocate
      DEALLOCATE ( this%cr_dw      , STAT=istatus ) ! Always deallocate
      DEALLOCATE ( this%cr_occ     , STAT=istatus ) ! Always deallocate
      DEALLOCATE ( this%atoms      , STAT=istatus ) ! Always deallocate
      DEALLOCATE ( this%cr_mole_len  , STAT=istatus ) ! Always deallocate molecules
      DEALLOCATE ( this%cr_mole_off  , STAT=istatus ) ! Always deallocate molecules
      DEALLOCATE ( this%cr_mole_type , STAT=istatus ) ! Always deallocate molecules
      DEALLOCATE ( this%cr_mole_char , STAT=istatus ) ! Always deallocate molecules
      DEALLOCATE ( this%cr_mole_file , STAT=istatus ) ! Always deallocate molecules
      DEALLOCATE ( this%cr_mole_cont , STAT=istatus ) ! Always deallocate molecules
      DEALLOCATE ( this%cr_mole_dens , STAT=istatus ) ! Always deallocate molecules
      DEALLOCATE ( this%cr_mole_biso , STAT=istatus ) ! Always deallocate molecules
      DEALLOCATE ( this%cr_mole_fuzzy, STAT=istatus ) ! Always deallocate molecules
!
!     IF(ALLOCATED(this%atoms)) DEALLOCATE ( this%atoms)
   ENDIF
   this%cr_natoms = 0
   this%latom  = .false.
!
   END SUBROUTINE finalize_atoms
!******************************************************************************
END MODULE cryst_class
