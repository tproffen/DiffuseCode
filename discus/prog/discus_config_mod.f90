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
!     MAXAT_CELL : Maximum number of atoms per unit cell
!
      INTEGER,  PARAMETER  ::  DEF_NMAX              = 1100000
      INTEGER,  PARAMETER  ::  DEF_MAXSCAT           =     100
      INTEGER,  PARAMETER  ::  DEF_MAXAT_CELL        =     250
!
      INTEGER        MAXAT_CELL
      INTEGER              ::  MAXSCAT               = 1
      INTEGER              ::  NMAX                  = 1
!
!     PARAMETER      (NMAX       = 1000000)
!     PARAMETER      (MAXSCAT    =  25)
      PARAMETER      (MAXAT_CELL = 250)
!
!#######################################################################
! (1b) Molecules
!
!     Currently the array sizes for molecules are set in the
!     file 'molecule_mod.f90'.
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
      INTEGER,  PARAMETER  ::  DEF_MAXQXY            = 301*301
      INTEGER,  PARAMETER  ::  DEF_CFPKT             =    4999
      REAL   ,  PARAMETER  ::  DEF_CFINC             =       0.001
!
      INTEGER             ::   MAXQXY
      INTEGER, PARAMETER  ::   CFPKT  = 4999
      REAL   (KIND=KIND(0.0D0)), PARAMETER  ::   CFINC  =    0.001D0
!      REAL   , PARAMETER  ::   CFINC  =    0.0001
!
!#######################################################################
! (3) Reverse-Monte-Carlo
!
!     RMC_MAX_Q      : Maximum number of input data points. 
!                      : This is the value for all points of all planes,
!                      : the size of a single plane must be <= MAXQXY !
!
!     RMC_MAX_SQ     : Size of Fourier-Array (RMC_MAX_Q * SYM)
!
!     RMC_MAX_PLANES : Maximum number of input data planes
!     RMC_MAX_SYM    : Maximum number of symmetrically equivalent planes
!     RMC_MAX_ATOM   : Maximum number of atoms moved in one RMC move
!     RMC_MAX_LOTS   : Maximum number of 'lots' allowed
!
!     INTEGER            RMC_MAX_Q ,RMC_MAX_PLANES,RMC_MAX_SYM
      INTEGER            RMC_MAX_ATOM
!     INTEGER            RMC_MAX_ATOM,RMC_MAX_LOTS,RMC_MAX_SQ
!
!     PARAMETER      (RMC_MAX_PLANES = 1)
!     PARAMETER      (RMC_MAX_SYM    = 1)
!     PARAMETER      (RMC_MAX_Q      = 100000 )
!     PARAMETER      (RMC_MAX_SQ     = RMC_MAX_SYM*RMC_MAX_Q) ! RMC_MAX_SYM*RMC_MAX_Q)
      PARAMETER      (RMC_MAX_ATOM   = 20)
!     PARAMETER      (RMC_MAX_LOTS   = 100)
!
!#######################################################################
! (4) CHEM and MC 
!
!     Note that the array sizes defined here are used in the 
!     CHEM sublevel as well as the Monte-Carlo (MC, AMC) levels.
!
!     CHEM_MAX_VEC  : Maximum number of neighbour vector definitions
!     CHEM_MAX_RAN  : Maximum number of neighbour range definitions
!     CHEM_MAX_ANG  : Maximum number of neighbour angle definitions
!     CHEM_MAX_ENV  : Maximum number of neighbours in environment def.
!     CHEM_MAX_COR  : Maximum number of correlations
!     CHEM_MAX_BIN  : Maximum number of points for histograms
!     CHEM_MAX_ATOM : Maximum number of different atoms
!     CHEM_MAX_NEIG : Maximum number of neighbouring atoms/molecules
!     CHEM_MAX_CENT : Maximum number of central atoms
!
!       INTEGER         CHEM_MAX_VEC
!       INTEGER         CHEM_MAX_ANG
        INTEGER         CHEM_MAX_COR
!       INTEGER         CHEM_MAX_RAN
!       INTEGER         CHEM_MAX_ENV
        INTEGER         CHEM_MAX_BIN,CHEM_MAX_ATOM
      INTEGER            CHEM_MAX_NEIG,MAX_ATOM_ENV
      INTEGER            CHEM_MAX_CENT
!
!       PARAMETER       (CHEM_MAX_VEC  = 200 )
!       PARAMETER       (CHEM_MAX_RAN  = 200 )
!       PARAMETER       (CHEM_MAX_ANG  = 200 )
!       PARAMETER       (CHEM_MAX_ENV  = 5000)
        PARAMETER       (CHEM_MAX_COR  = 200 )
        PARAMETER       (CHEM_MAX_BIN  = 2001)
        PARAMETER       (CHEM_MAX_ATOM = 20  )
!     Warning!      MAXPAR_RES  in "param_mod.f90" must always be of 
!                                 identical size to
!                 CHEM_MAX_NEIG in "config_mod.f90"
        PARAMETER       (CHEM_MAX_NEIG = 6000 )
       PARAMETER        (MAX_ATOM_ENV  = CHEM_MAX_NEIG)
      PARAMETER      (CHEM_MAX_CENT = 200 )
!
!#######################################################################
! (5) Microdomains
!
!     MMAX         : Maximum number of atoms within a microdomain
!     MK_MAX_SCAT  : Maximum number of atom types within a microdomain
!
!      INTEGER            MMAX,MK_MAX_SCAT
!
!      PARAMETER      (MMAX         = 250   )
!      PARAMETER      (MK_MAX_SCAT  = DEF_MAXSCAT)
!
!#######################################################################
! (6) Stacking faults
!
!     ST_MAXQXY    : Maximum number of points of Fourier (<= MAXQXY !)
!     ST_MAXLAYER  : Maximum number of layers
!     ST_MAXTYPE   : Maximum number of layer types
!
!     INTEGER, PARAMETER  ::  ST_MAXQXY    = DEF_MAXQXY
!     INTEGER, PARAMETER  ::  ST_MAXLAYER  = 5000
!     INTEGER, PARAMETER  ::  ST_MAXTYPE   = 6
!     INTEGER, PARAMETER  ::  ST_MMAX      = 5000
!     INTEGER, PARAMETER  ::  ST_MAX_SCAT  = DEF_MAXSCAT
!
!#######################################################################
! (7) Pair distribution function
!
!     MAXDAT       : Maximum number of points in PDF
!     MAXBND       : Maximum extend of periodic boundary array
! 
!     INTEGER            MAXDAT,MAXBND
!
!     PARAMETER      (MAXDAT = 15000)
!     PARAMETER      (MAXBND =  1000)
!
!#######################################################################
! (8) Powder module
!
!     POW_MAXPKT   : Maximum number of points in powder pattern
!     MAXHIST      : Maximum number of points in histogram used by debye
!
!     INTEGER        MAXHIST
!     INTEGER        POW_MAXPKT
!
!     PARAMETER      (POW_MAXPKT=1000  )
!     PARAMETER      (MAXHIST   =3201  )
!     PARAMETER      (POW_MAXPKT=100000)
!     PARAMETER      (MAXHIST   =320001)
!

! create kind for portable double precision real number
!      integer, parameter:: dp=kind(0.d0)  ! double precision

END MODULE discus_config_mod
