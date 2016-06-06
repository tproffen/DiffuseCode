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
INTEGER                                 ::  nlook
!
REAL (KIND(0.0D0)), DIMENSION(:),   ALLOCATABLE    ::  rsf          ! (1:MAXQXY)
!!!   REAL                                  ::  sinetab(0:MASK) !!! MASK is in diffuse.inc!
REAL (KIND(0.0D0)), DIMENSION(:),   ALLOCATABLE    ::  sinetab      ! (0:2**16-1)
REAL                                    ::  pow_del_hist = 0.0100000000 ! = 0.001 can be changed by user
INTEGER                                 ::  deb_size_of  = 0 ! Bytes allocated for DEBYE
!
END MODULE debye_mod
