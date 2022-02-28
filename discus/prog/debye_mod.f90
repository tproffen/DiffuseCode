MODULE debye_mod
!+
!
!     Contains all variables for Debye formula
!-
!
use precision_mod
SAVE
!
INTEGER                               ::  MAXLOOK   =1 ! = MAXSCAT*(MAXSCAT+1)/2
INTEGER                               ::  MAXHIST   =1 ! 
INTEGER                               ::  MAXDSCAT  =1 ! 
INTEGER                               ::  MAXDQXY   =1 ! 
INTEGER                               ::  DEB_MAXMASK   =1 ! 
!
INTEGER                                 ::  nlook
!
REAL(kind=PREC_DP), DIMENSION(:),   ALLOCATABLE    ::  rsf          ! (1:MAXQXY)
REAL(kind=PREC_DP), DIMENSION(:),   ALLOCATABLE    ::  sinetab      ! (0:2**16-1)
REAL(kind=PREC_DP)                      ::  pow_del_hist = 0.0100000000 ! = 0.001 can be changed by user
INTEGER                                 ::  deb_size_of  = 0 ! Bytes allocated for DEBYE
LOGICAL                                 ::  deb_conv     = .FALSE.  ! convolute ADP's
!
END MODULE debye_mod
