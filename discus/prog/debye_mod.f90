MODULE debye_mod
!+
!
!     Contains all variables for Debye formula
!-
!
SAVE
!
INTEGER                               ::  MAXLOOK   =1 ! = MAXSCAT*(MAXSCAT+1)/2
INTEGER                               ::  MAXHIST   =1 ! 
INTEGER                               ::  MAXDSCAT  =1 ! 
INTEGER                               ::  MAXDQXY   =1 ! 
INTEGER                               ::  DEB_MAXMASK   =1 ! 
!
INTEGER, DIMENSION(:,:), ALLOCATABLE  ::  histogram    ! (MAXHIST,MAXLOOK)
INTEGER                               ::  nlook
INTEGER, DIMENSION(:,:), ALLOCATABLE  ::  look         ! (MAXSCAT,MAXSCAT)
!
REAL   , DIMENSION(:),   ALLOCATABLE  ::  rsf          ! (1:MAXQXY)
REAL   , DIMENSION(:,:), ALLOCATABLE  ::  partial      ! (32001,MAXLOOK)
!!!   REAL                                  ::  sinetab(0:MASK) !!! MASK is in diffuse.inc!
REAL   , DIMENSION(:),   ALLOCATABLE  ::  sinetab      ! (0:2**16-1)
REAL                                  ::  pow_del_hist = 0.0100000000 ! = 0.001 can be changed by user
INTEGER                               ::  deb_size_of  = 0 ! Bytes allocated for DEBYE
!
END MODULE debye_mod
