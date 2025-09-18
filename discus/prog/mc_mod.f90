MODULE mc_mod
!+
!
!     Variables for MONTE-CARLO level
!-
!
USE discus_config_mod
!
USE precision_mod
!
SAVE
!
INTEGER, PARAMETER  ::  MC_NONE     = 0
INTEGER, PARAMETER  ::  MC_OCC      = 1
INTEGER, PARAMETER  ::  MC_DISP     = 2
INTEGER, PARAMETER  ::  MC_SPRING   = 3
INTEGER, PARAMETER  ::  MC_ANGLE    = 4
INTEGER, PARAMETER  ::  MC_VECTOR   = 5
INTEGER, PARAMETER  ::  MC_BLEN     = 6
INTEGER, PARAMETER  ::  MC_LENNARD  = 7
INTEGER, PARAMETER  ::  MC_BUCKING  = 8
INTEGER, PARAMETER  ::  MC_REPULSIVE= 9
INTEGER, PARAMETER  ::  MC_COORDNUM = 10
INTEGER, PARAMETER  ::  MC_UNI      = 11
INTEGER, PARAMETER  ::  MC_GROUP    = 12
INTEGER, PARAMETER  ::  MC_PREF     = 13
INTEGER, PARAMETER  ::  MC_VALUE    = 14
!
INTEGER, PARAMETER  ::  MC_N_ENERGY        = 14
!
LOGICAL, PARAMETER  :: MMC_IS_ATOM = .TRUE.
LOGICAL, PARAMETER  :: MMC_IS_MOLE = .FALSE.
!
CHARACTER(LEN=200)  ::  mo_atom(3)
!
INTEGER             ::  mo_energy
INTEGER             ::  mo_mode,mo_local
INTEGER(KIND=PREC_INT_LARGE)  ::  mo_cyc,mo_feed
!
REAL(kind=PREC_DP), DIMENSION(:)    , ALLOCATABLE                ::  mo_ach_corr    ! (200)
REAL(kind=PREC_DP), DIMENSION(:,:  ), ALLOCATABLE ::  mo_maxmove     ! (4,0:DEF_MAXSCAT)
REAL(kind=PREC_DP), DIMENSION(:,:  ), ALLOCATABLE ::  mo_maxmove_mole! (4,0:DEF_MAXSCAT)
REAL(kind=PREC_DP), DIMENSION(:,:  ), ALLOCATABLE ::  mo_maxrota_mole! (4,0:DEF_MAXSCAT)
REAL(kind=PREC_DP)                ::  mo_kt
LOGICAL             ::  mo_sel_atom
!
END MODULE mc_mod
