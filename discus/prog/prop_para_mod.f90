MODULE prop_para_mod
!+
!
!     Parameter definitions for property settings
!-
SAVE
!
INTEGER, PARAMETER ::  MINPROP           =  1
INTEGER, PARAMETER ::  MAXPROP           =  6
!
INTEGER, PARAMETER ::  PROP_NORMAL       =  0
INTEGER, PARAMETER ::  PROP_MOLECULE     =  1
INTEGER, PARAMETER ::  PROP_DOMAIN       =  2
INTEGER, PARAMETER ::  PROP_OUTSIDE      =  3
INTEGER, PARAMETER ::  PROP_SURFACE_EXT  =  4
INTEGER, PARAMETER ::  PROP_SURFACE_INT  =  5
!
INTEGER, PARAMETER ::  PROP_IGNORE       =  MAXPROP+1
!
CHARACTER(LEN=8)   ::  c_prop_letter     = 'NMDOEI  '
!
!     COMMON /property_params/ c_prop_letter
!
END MODULE prop_para_mod
