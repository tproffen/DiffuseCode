MODULE mole_env_mod
!+
!
!     include file for molecular environment variables
!-
USE atom_env_mod
use precision_mod
!
SAVE
!
INTEGER, PARAMETER ::  MAX_MOLE_ENV  =  MAX_ATOM_ENV
!
INTEGER           , DIMENSION(0:MAX_MOLE_ENV) ::  mole_env ! (0:MAX_MOLE_ENV)
REAL(kind=PREC_DP), DIMENSION(3,MAX_MOLE_ENV) ::  mole_pos ! (3,MAX_MOLE_ENV)
!
END MODULE mole_env_mod
