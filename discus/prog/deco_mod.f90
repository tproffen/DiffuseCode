MODULE deco_mod
!+
!
!     variables needed for the decoration of structures
!-
USE discus_config_mod
USE cryst_class
USE dc_def_class
!
SAVE
!
INTEGER                                   :: DC_MAXSCAT = 1
!
!  Temporary arrays
!
TYPE(cl_cryst), DIMENSION(:), POINTER     :: dc_molecules  ! temporary molecules 
INTEGER, DIMENSION(:), ALLOCATABLE        :: m_ntypes ! number of atom tyes in molecules
INTEGER, DIMENSION(:), ALLOCATABLE        :: m_length ! number of atoms in molecules
!
!  INTPUT definitions
!
INTEGER                                   :: dc_temp_id = 0  ! Temporary number for definition structure
CHARACTER (LEN=1024)                      :: dc_temp_name    ! Temporary name for definition structure
INTEGER                                   :: dc_temp_lname   ! Temporary name length for definition structure
CHARACTER (LEN=1024)                      :: dc_temp_file    ! Temporary filename for input
INTEGER                                   :: dc_temp_lfile   ! Temporary filename length for input
!
!TYPE dc_def
!   INTEGER              :: dc_def_lname
!   CHARACTER(LEN=1024)  :: dc_def_name
!   INTEGER              :: dc_def_lfile
!   CHARACTER(LEN=1024)  :: dc_def_file
!   TYPE (dc_def), POINTER :: next
!END TYPE dc_def
TYPE (dc_def), POINTER                    :: dc_def_head => NULL()
TYPE (dc_def), POINTER                    :: dc_def_tail => NULL()
TYPE (dc_def), POINTER                    :: dc_def_temp => NULL()
INTEGER                                   :: dc_def_number = 0  ! number of definitions
!CHARACTER(LEN= 1024), DIMENSION(:), ALLOCATABLE  :: dc_input = 'molecule.stru'
CHARACTER(LEN= 1024), DIMENSION(1:200)  :: dc_input = 'molecule.stru'
!
INTEGER, DIMENSION(0:1)                   :: dc_sel_prop  = (/0,0/)

INTEGER, DIMENSION(0:200)                 :: dc_use_conn      ! Use this conn to check surface atom
INTEGER                                   :: dc_n_molecules = 0  ! number of molecule types
LOGICAL                                   :: dc_init      = .true. ! do we need to initialize
!LOGICAL, DIMENSION(:), ALLOCATABLE        :: dc_latom     ! (0:MAXSCAT)
LOGICAL, DIMENSION(0:200    )             :: dc_latom     ! (0:MAXSCAT)
LOGICAL                                   :: dc_sel_atom  = .true.
LOGICAL                                   :: dc_mol_all   = .true.
!
! The molecules
!
TYPE :: is_atom
   INTEGER                :: is_iscat
   REAL   , DIMENSION(3)  :: is_pos
   INTEGER                :: is_prop
END TYPE is_atom
!
!
END MODULE deco_mod
