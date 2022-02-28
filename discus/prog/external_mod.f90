MODULE external_mod
!+
!
!     variables needed for the externally computed structure factors
!-
!
USE molecule_mod
use precision_mod
!
SAVE
!
INTEGER, PARAMETER  ::  EXTE_HLINES  =  10
!
CHARACTER(LEN=200), DIMENSION(:), ALLOCATABLE  ::  exte_names ! (MOLE_MAX_TYPE)
CHARACTER(LEN=200)                             ::  exte_filename
INTEGER, DIMENSION(:), ALLOCATABLE ::  exte_length! (MOLE_MAX_TYPE)
!
INTEGER, DIMENSION(:), ALLOCATABLE ::  exte_type ! (MOLE_MAX_TYPE)
!
INTEGER                 ::  exte_version
INTEGER                 ::  exte_hdrblks
CHARACTER(LEN=144)      ::  exte_title
INTEGER                 ::  exte_nrows
INTEGER                 ::  exte_ncols
INTEGER                 ::  exte_layer
INTEGER                 ::  exte_npixelb
INTEGER                 ::  exte_wordord
INTEGER                 ::  exte_longord
INTEGER, DIMENSION(3)   ::  exte_orig ! (3)
REAL(kind=PREC_DP)                    ::  exte_scale
!
INTEGER, DIMENSION(3)   ::  exte_iii  ! (3)
REAL(kind=PREC_DP)   , DIMENSION(4,4) ::  exte_mat  ! (4,4)
REAL(kind=PREC_DP)   , DIMENSION(4,4) ::  exte_rmat ! (4,4)
REAL(kind=PREC_DP)   , DIMENSION(4,4) ::  exte_rot  ! (4,4)
REAL(kind=PREC_DP)   , DIMENSION(4,4) ::  exte_rrot ! (4,4)
!
END MODULE external_mod
