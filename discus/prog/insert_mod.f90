MODULE insert_mod
!-
!     Variables needed to insert objects
!+
!
use precision_mod
!
SAVE
!
INTEGER, PARAMETER               ::  INS_NEWTYPE  =  -1
!
CHARACTER(LEN=4  )               ::  ins_obj_atom  = 'VOID'
CHARACTER(LEN=200)               ::  ins_file      = ' '
INTEGER                          ::  ins_CHARACTER = 0
INTEGER                          ::  ins_type      = 0
REAL(kind=PREC_DP), DIMENSION(3) ::  ins_origin    = (/0.0, 0.0, 0.0/) ! (3)
REAL(kind=PREC_DP), DIMENSION(3) ::  ins_cent      = (/0.0, 0.0, 0.0/) ! (3)
REAL(kind=PREC_DP)               ::  ins_density   = 0.0
REAL(kind=PREC_DP)               ::  ins_biso      = 0.0
REAL(kind=PREC_DP)               ::  ins_clin      = 0.0
REAL(kind=PREC_DP)               ::  ins_cqua      = 0.0
REAL(kind=PREC_DP)               ::  ins_fuzzy     = 0.0
REAL(kind=PREC_DP)               ::  ins_adp       = 0.0
REAL(kind=PREC_DP), DIMENSION(3) ::  ins_xaxis     = (/1.0, 0.0, 0.0/) ! (3)
REAL(kind=PREC_DP), DIMENSION(3) ::  ins_yaxis     = (/0.0, 1.0, 0.0/) ! (3)
REAL(kind=PREC_DP), DIMENSION(3) ::  ins_zaxis     = (/0.0, 0.0, 1.0/) ! (3)
REAL(kind=PREC_DP), DIMENSION(3) ::  ins_xdim      = (/1.0, 0.0, 0.0/) ! (3)
REAL(kind=PREC_DP), DIMENSION(3) ::  ins_ydim      = (/0.0, 1.0, 0.0/) ! (3)
REAL(kind=PREC_DP), DIMENSION(3) ::  ins_zdim      = (/0.0, 0.0, 1.0/) ! (3)
!
END MODULE insert_mod
