MODULE atom_env_mod
!+
!
!     include file for atomic environment variables
!-
USE config_mod
!
SAVE
!
INTEGER, PARAMETER  ::  MAX_ATOM_ENV_D  =  MAX_ATOM_ENV+1
!
INTEGER             ::  atom_env(0:MAX_ATOM_ENV) = 0
REAL                ::  atom_pos(3,MAX_ATOM_ENV) = 0.0
!
!     COMMON  /at_env/ atom_env,atom_pos
!
END MODULE atom_env_mod
