MODULE atom_env_mod
!+
!
!     include file for atomic environment variables
!-
SAVE
!
INTEGER, PARAMETER  ::  MAX_ATOM_ENV    =  6000
INTEGER, PARAMETER  ::  MAX_ATOM_ENV_D  =  MAX_ATOM_ENV+1
!
INTEGER             ::  atom_env(  0:MAX_ATOM_ENV) = 0
REAL                ::  atom_pos(3,0:MAX_ATOM_ENV) = 0.0
REAL                ::  atom_dis(  0:MAX_ATOM_ENV) = 0.0
!
END MODULE atom_env_mod
