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
INTEGER, PARAMETER                        :: DC_NONE     = 0
INTEGER, PARAMETER                        :: DC_NORMAL   = 1
INTEGER, PARAMETER                        :: DC_BRIDGE   = 2
INTEGER, PARAMETER                        :: DC_DOUBLE   = 3
INTEGER, PARAMETER                        :: DC_MULTIPLE = 4
!
INTEGER, PARAMETER                        :: SURF_NONE   = 0
INTEGER, PARAMETER                        :: SURF_PLANE  = 1
!
!  Temporary arrays
!
TYPE(cl_cryst), DIMENSION(:), POINTER     :: dc_molecules  ! temporary molecules 
CHARACTER (LEN=1024), DIMENSION(:), ALLOCATABLE :: m_name   ! Name of molecule file
INTEGER, DIMENSION(:), ALLOCATABLE        :: m_lname  ! Length of molecule file length
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
INTEGER                                   :: dc_temp_type    ! Temporary connection type for input
INTEGER, DIMENSION(0:4)                   :: dc_temp_surf    ! Temporary surface atom type for input
INTEGER                                   :: dc_temp_neig    ! Temporary neighbor atom number for input
INTEGER, DIMENSION(1:2)                   :: dc_temp_axis    ! Temporary axis for ligand      for input
REAL                                      :: dc_temp_dist    ! Temporary neighbor distance    for input
!
TYPE (dc_def), POINTER                    :: dc_def_head => NULL()
!TYPE (dc_def), POINTER                    :: dc_def_tail => NULL()
TYPE (dc_def), POINTER                    :: dc_def_temp => NULL()
TYPE (dc_con), POINTER                    :: dc_con_temp => NULL()
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
