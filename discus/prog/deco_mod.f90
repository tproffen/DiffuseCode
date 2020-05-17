MODULE deco_mod
!+
!
!     variables needed for the decoration of structures
!-
USE precision_mod
USE discus_config_mod
!
INTEGER, PARAMETER        :: DC_MAXMODE  = 7  ! We have six decoration modes
!
SAVE
!
INTEGER                                   :: DC_MAXSCAT = 1
INTEGER                                   :: DC_MAXSITE = 1
!
INTEGER, PARAMETER                        :: DC_NONE     = 0
INTEGER, PARAMETER                        :: DC_NORMAL   = 1
INTEGER, PARAMETER                        :: DC_BRIDGE   = 2
INTEGER, PARAMETER                        :: DC_DOUBLE   = 3
INTEGER, PARAMETER                        :: DC_MULTIPLE = 4
INTEGER, PARAMETER                        :: DC_ACCEPTOR = 5
INTEGER, PARAMETER                        :: DC_DONOR    = 6
INTEGER, PARAMETER                        :: DC_CHELATE  = DC_MAXMODE
!
INTEGER, PARAMETER                        :: SURF_ATOM   =-1
INTEGER, PARAMETER                        :: SURF_NONE   = 0
INTEGER, PARAMETER                        :: SURF_PLANE  = 1
!
REAL   , PARAMETER                        :: DC_AREA     = 11.00  ! Estimates area per surface atom in [A^2]
!
CHARACTER(LEN=8),DIMENSION(0:DC_MAXMODE)  :: dcc_ctype
!
!  All the variables for the surface decoration definitions
!
INTEGER                                            :: dcc_MAXNUM          ! Maximum Number of definitions
INTEGER                                            :: dcc_MAXANCH         ! Maximum Number of anchors
INTEGER                                            :: dcc_MAXHKL          ! Maximum Number of faces
INTEGER                                            :: dcc_MAXNEW          ! Maximum Number of new surface atoms
INTEGER                                            :: dcc_MAXMSCAT        ! Maximum Number of molecule atom types
INTEGER                                            :: dcc_num = 0         ! Number of definitions
INTEGER                                            :: dcc_nhkl = 0        ! Number of faces
INTEGER                                            :: dcc_nnew = 0        ! Number of new surface atoms
INTEGER                                            :: dcc_maxsurf = 0     ! Max no of surface sites
CHARACTER(LEN=PREC_STRING), DIMENSION(    :), ALLOCATABLE :: dcc_name            ! Definition names
CHARACTER(LEN=PREC_STRING), DIMENSION(    :), ALLOCATABLE :: dcc_file            ! Content file name
INTEGER            , DIMENSION(    :), ALLOCATABLE :: dcc_natoms          ! Content number of atoms
CHARACTER(LEN=4)   , DIMENSION(:,  :), ALLOCATABLE :: dcc_atom_name       ! Content atom names
REAL               , DIMENSION(:,  :), ALLOCATABLE :: dcc_adp             ! Content ADP's
REAL               , DIMENSION(    :), ALLOCATABLE :: dcc_biso            ! Content molecule ADP's
REAL               , DIMENSION(    :), ALLOCATABLE :: dcc_clin            ! Content molecule ADP's
REAL               , DIMENSION(    :), ALLOCATABLE :: dcc_cqua            ! Content molecule ADP's
INTEGER            , DIMENSION(    :), ALLOCATABLE :: dcc_mole_type       ! Content molecule type
INTEGER            , DIMENSION(    :), ALLOCATABLE :: dcc_type            ! definitions type
INTEGER            , DIMENSION(    :), ALLOCATABLE :: dcc_lname           ! Length definitions name
INTEGER            , DIMENSION(    :), ALLOCATABLE :: dcc_lfile           ! Length content file name
INTEGER            , DIMENSION(:  ,:), ALLOCATABLE :: dcc_nanch           ! Number of anchor atoms
INTEGER            , DIMENSION(:,:,:), ALLOCATABLE :: dcc_surf            ! Surface anchors
INTEGER            , DIMENSION(:,  :), ALLOCATABLE :: dcc_neig            ! Molecule neighbor
INTEGER            , DIMENSION(    :), ALLOCATABLE :: dcc_secnd           ! Molecule neighbor closest to dcc_neig
INTEGER            , DIMENSION(:,  :), ALLOCATABLE :: dcc_axis            ! Molecule axis
LOGICAL            , DIMENSION(    :), ALLOCATABLE :: dcc_spread          ! Spread anchors evenly, or keep fixed
LOGICAL            , DIMENSION(    :), ALLOCATABLE :: dcc_lrestrict       ! Restiction by form or hkl
LOGICAL            , DIMENSION(    :), ALLOCATABLE :: dcc_lform           ! Restiction by form or hkl
INTEGER            , DIMENSION(:,:,:), ALLOCATABLE :: dcc_hkl             ! Surface restriction
INTEGER            , DIMENSION(:,  :), ALLOCATABLE :: dcc_surfnew         ! Molecule atoms as new surface
REAL               , DIMENSION(    :), ALLOCATABLE :: dcc_dens            ! Molecule density per A^2
REAL               , DIMENSION(:,  :), ALLOCATABLE :: dcc_dist            ! Bond length in       A
REAL               , DIMENSION(    :), ALLOCATABLE :: dcc_angle           ! Bond angle for Hydrogen bonds
REAL               , DIMENSION(    :), ALLOCATABLE :: dcc_tilt            ! Tilt angle for ligand off axis
REAL               , DIMENSION(:,  :), ALLOCATABLE :: dcc_tilt_hkl        ! Normal to molecule plane
INTEGER            , DIMENSION(:,  :), ALLOCATABLE :: dcc_tilt_atom       ! Atoms that form molecule plane
LOGICAL            , DIMENSION(    :), ALLOCATABLE :: dcc_tilt_is_atom    ! Tilt angle for ligand off axis
LOGICAL            , DIMENSION(    :), ALLOCATABLE :: dcc_tilt_is_auto    ! Tilt angle for ligand off axis
!
!  INTPUT definitions
!
INTEGER                                   :: dc_temp_id = 0  ! Temporary number for definition structure
CHARACTER (LEN=PREC_STRING)               :: dc_temp_name    ! Temporary name for definition structure
INTEGER                                   :: dc_temp_lname   ! Temporary name length for definition structure
CHARACTER (LEN=PREC_STRING)               :: dc_temp_file    ! Temporary filename for input
INTEGER                                   :: dc_temp_lfile   ! Temporary filename length for input
INTEGER                                   :: dc_temp_type = DC_NONE   ! Temporary connection type for input
INTEGER, DIMENSION(0:40)                  :: dc_temp_surf = 0! Temporary surface atom type for input
INTEGER                                   :: dc_temp_maxsurf ! Maximum number of surface atom typers
INTEGER                                   :: dc_temp_neig    ! Temporary neighbor atom number for input
INTEGER, DIMENSION(0:2)                   :: dc_temp_axis    ! Temporary axis for ligand      for input
INTEGER, DIMENSION(1:40)                  :: dc_temp_surfnew=0 ! Temporary list of new surface atoms
REAL                                      :: dc_temp_dist    ! Temporary neighbor distance    for input
REAL                                      :: dc_temp_dens    ! Temporary ligand density       for input
LOGICAL                                   :: dc_temp_restrict = .FALSE. ! Restriction yes / no
LOGICAL                                   :: dc_temp_l_form   = .FALSE. ! Single hkl or form
INTEGER                                   :: dc_temp_n_hkl    = 0       ! Number of hkls that make special form
INTEGER,DIMENSION(:,:), ALLOCATABLE       :: dc_temp_hkl     ! Surface restrictions
!
INTEGER                                   :: dc_def_number = 0  ! number of definitions
CHARACTER(LEN= PREC_STRING), DIMENSION(1:200)  :: dc_input = 'molecule.stru'
!
INTEGER, DIMENSION(0:1)                   :: dc_sel_prop  = (/0,0/)

INTEGER, DIMENSION(0:200)                 :: dc_use_conn      ! Use this conn to check surface atom
INTEGER                                   :: dc_n_molecules = 0  ! number of molecule types
LOGICAL                                   :: dc_init      = .true. ! do we need to initialize
LOGICAL, DIMENSION(:), ALLOCATABLE        :: dc_latom     ! (0:MAXSCAT)
LOGICAL, DIMENSION(:), ALLOCATABLE        :: dc_lsite     ! (0:MAXSCAT)
LOGICAL                                   :: dc_sel_atom  = .true.
LOGICAL                                   :: dc_mol_all   = .true.
!
!
DATA dcc_ctype /'NONE    ', 'NORMAL  ', 'BRIDGE  ', 'DOUBLE  ', 'MULTIPLE', &
                'ACCEPTOR', 'DONOR   ', 'CHELATE ' /
!
END MODULE deco_mod
