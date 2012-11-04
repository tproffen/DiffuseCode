MODULE deco_mod
!+
!
!     variables needed for the decoration of structures
!-
USE config_mod
USE cryst_class
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
!
!CHARACTER(LEN= 1024), DIMENSION(:), ALLOCATABLE  :: dc_input = 'molecule.stru'
CHARACTER(LEN= 1024), DIMENSION(1:200)  :: dc_input = 'molecule.stru'
!
INTEGER, DIMENSION(0:1)                   :: dc_sel_prop  = (/0,0/)

INTEGER, DIMENSION(0:200)                 :: dc_use_conn      ! Use this conn to check surface atom
INTEGER                                   :: dc_n_molecules   ! number of molecule types
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
