MODULE shear_mod
!+
!     variables needed for the generalized symmetry operations
!-
!
USE discus_config_mod
!
SAVE
!
INTEGER, PRIVATE    :: ik
INTEGER             ::  SHEAR_MAXSCAT = 1
INTEGER             ::  SHEAR_MAXSITE = 1
!
LOGICAL,          DIMENSION(:), ALLOCATABLE  ::  shear_latom  ! (0:MAXSCAT)
LOGICAL,          DIMENSION(:), ALLOCATABLE  ::  shear_lsite  ! (0:MAXSCAT)
!
!
INTEGER, PARAMETER  ::    SHEAR_ATOM      =  0
INTEGER, PARAMETER  ::    SHEAR_MOLECULE  =  1
INTEGER, PARAMETER  ::    SHEAR_OBJECT    =  2
INTEGER, PARAMETER  ::    SHEAR_DOMAIN    =  3
!
INTEGER, PARAMETER  ::    SHEAR_MATRIX    =  0
INTEGER, PARAMETER  ::    SHEAR_RMATRIX   =  1
INTEGER, PARAMETER  ::    SHEAR_PLANE     =  2
INTEGER, PARAMETER  ::    SHEAR_EIGEN     =  3
!
CHARACTER(LEN=4)        :: shear_incl
INTEGER, DIMENSION(0:1) :: shear_sel_prop =  0
INTEGER                 :: shear_start    =  1
INTEGER                 :: shear_end      = -1
INTEGER                 :: shear_mode     = SHEAR_OBJECT
INTEGER                 :: shear_input    = SHEAR_MATRIX
LOGICAL                 :: shear_new      = .false.
LOGICAL                 :: shear_orig_mol = .false.
LOGICAL                 :: shear_sel_atom = .true.
LOGICAL                 :: shear_dom_mode_atom  = .true.
LOGICAL                 :: shear_dom_mode_shape = .true.
REAL   , DIMENSION(3)   :: shear_hkl    = (/0.,1.,0./)
REAL   , DIMENSION(3)   :: shear_orig   = 0.0
REAL   , DIMENSION(3)   :: shear_vector = (/1.,0.,0./)
REAL                    :: shear_length = 0.0
REAL   , DIMENSION(3)   :: shear_uvw    = (/0.,1.,0./)
REAL   , DIMENSION(4,4) :: shear_mat    = &
         RESHAPE((/1.,(0.,0.,0.,0.,1.,ik=1,3)/),SHAPE(shear_mat ))
REAL   , DIMENSION(4,4) :: shear_rmat   = &
         RESHAPE((/1.,(0.,0.,0.,0.,1.,ik=1,3)/),SHAPE(shear_rmat))
REAL   , DIMENSION(3,3) :: shear_eigenv = &
         RESHAPE((/1.,(0.,0.,0.,1.,ik=1,2)/),SHAPE(shear_eigenv ))
REAL   , DIMENSION(3,3) :: shear_eigent = &
         RESHAPE((/1.,(0.,0.,0.,1.,ik=1,2)/),SHAPE(shear_eigent ))
REAL   , DIMENSION(3)   :: shear_eigenw = 1.0
!
INTEGER                 :: shear_size_of  = 0 ! Bytes allocated for shear
!
END MODULE shear_mod
