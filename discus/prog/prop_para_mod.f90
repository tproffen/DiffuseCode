MODULE prop_para_mod
!+
!
!     Parameter definitions for property settings
!-
SAVE
!
INTEGER, PARAMETER ::  MINPROP           =  1
INTEGER, PARAMETER ::  MAXPROP           =  8
!
INTEGER, PARAMETER ::  PROP_NORMAL       =  0
INTEGER, PARAMETER ::  PROP_MOLECULE     =  1
INTEGER, PARAMETER ::  PROP_DOMAIN       =  2
INTEGER, PARAMETER ::  PROP_OUTSIDE      =  3
INTEGER, PARAMETER ::  PROP_SURFACE_EXT  =  4
INTEGER, PARAMETER ::  PROP_SURFACE_INT  =  5
INTEGER, PARAMETER ::  PROP_LIGAND       =  6
INTEGER, PARAMETER ::  PROP_TEMP         =  7
!
INTEGER, PARAMETER ::  PROP_DECO_ANCHOR  =  8
!
INTEGER, PARAMETER ::  PROP_IGNORE       =  MAXPROP+1
!
CHARACTER(LEN=8)   ::  c_prop_letter     = 'NMDOEILT'
CHARACTER(LEN=8)   ::  c_prop_small      = 'nmdoeilt'
!
INTEGER            ::  prop_user_no      = 0
!
TYPE :: prop_templ
   INTEGER :: act
   INTEGER :: at_type
   INTEGER :: conn_no
   CHARACTER(LEN=256) :: conn_name
   INTEGER :: n_min
   INTEGER :: n_max
   INTEGER :: e_min
   INTEGER :: e_max
END TYPE prop_templ
!
TYPE(prop_templ), DIMENSION(:), ALLOCATABLE :: prop_user
!
END MODULE prop_para_mod
