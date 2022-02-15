MODULE patters_mod
!+
!
!     This file contains variables for inverse Fourier and
!     Patterson input
!-
use precision_mod
SAVE
!
INTEGER, PARAMETER  ::  PATT_INIT      =  0
INTEGER, PARAMETER  ::  PATT_ADD       =  1
INTEGER, PARAMETER  ::  PATT_NORMAL    =  0
INTEGER, PARAMETER  ::  PATT_SHARP     =  1
INTEGER, PARAMETER  ::  PATT_SUPER     =  2
INTEGER, PARAMETER  ::  PATT_SUBTRACT  =  1
!
CHARACTER(LEN=200), DIMENSION(2)  ::   rho_file = (/' ',' '/)
INTEGER                 ::  patt_accu     = 0
INTEGER, DIMENSION(3)   ::  rho_inc       = (/121,121, 1/)
INTEGER, DIMENSION(2)   ::  rho_type      = (/4,5/)
INTEGER                 ::  patt_sign     = -1
INTEGER                 ::  patt_mode     = PATT_NORMAL
INTEGER                 ::  patt_origin   = 0
INTEGER                 ::  ftyp          = 0
LOGICAL                 ::  patt_excl9999 = .false.
LOGICAL                 ::  patt_rsym     = .false.
REAL(kind=PREC_DP)   , DIMENSION(3,3) ::  rho_eck       = reshape((/0,0,0, 1,0,0, 0,1,0/),shape(rho_eck))
REAL(kind=PREC_DP)   , DIMENSION(3,2) ::  rho_vi        = reshape((/0.02,0.0,0.0, 0.0, 0.02,0.0/),shape(rho_vi))
REAL(kind=PREC_DP)                    ::  patt_scale    = 1.0
REAL(kind=PREC_DP)                    ::  patt_excl_val = -9999.0
REAL(kind=PREC_DP)   , DIMENSION(20)  ::  e_aver_f2     = 0.0
REAL(kind=PREC_DP)                    ::  wilson_scale  = 1.0
REAL(kind=PREC_DP)                    ::  wilson_b      = 0.0
!
END MODULE patters_mod
