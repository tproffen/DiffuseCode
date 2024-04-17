MODULE discus_config_mod
!
SAVE
!#######################################################################
!
!
! DISCUS Configuration file
!#######################################################################
! (1) Crystal dimensions
!
!     MAXSCAT : Maximum number of different atomtypes
!     NMAX    : Maximum number of atoms
!
      INTEGER              ::  MAXSCAT               = 1
      INTEGER              ::  NMAX                  = 1
!
!#######################################################################
! (2) Fourier transform
!
!     MAXQXY  : Maximum number of points in Q 
!     CFPKT   : Number of points in SIN(THETA/LAMBDA) lookup table
!     CFINC   : Increment for SIN(THETA/LAMBDA) table. The lookup
!               table ranges from 0 to (CFPKT+1)*CFINC
!
!
INTEGER, dimension(3) ::   MAXQXY
      INTEGER, PARAMETER  ::   CFPKT  = 9999
      REAL(KIND=KIND(0.0D0)), PARAMETER  ::   CFINC  =    0.001D0
!
!#######################################################################
!
END MODULE discus_config_mod
