MODULE param_mod
!+
!     Include file for free variables
!
!     Warning!      MAXPAR_RES  in "param.inc" must always be of 
!                                 identical size to
!                 CHEM_MAX_NEIG in DISCUS/"config.inc"
!     This warning is obsolete.
!-
   IMPLICIT NONE
   PUBLIC
   SAVE
!
   INTEGER, PARAMETER :: MAXPAR     =  500
!  INTEGER, PARAMETER :: MAXPAR_RES = 12000
   INTEGER            :: MAXPAR_RES = 6000
!
   INTEGER, DIMENSION(0:MAXPAR)     :: inpara   = 0   ! (0:MAXPAR)
   REAL   , DIMENSION(0:MAXPAR)     :: rpara    = 0.0 ! (0:MAXPAR)
!  REAL   , DIMENSION(0:MAXPAR_RES) :: res_para = 0.0 ! (0:MAXPAR_RES)
   REAL   , DIMENSION(:),ALLOCATABLE:: res_para       ! (0:MAXPAR_RES)
   REAL   , DIMENSION(1:2)          :: rvalues  = 0.0 ! (1:2)
   LOGICAL                          :: rvalue_yes = .false.
!
END MODULE param_mod
