MODULE param_mod
!+
!     Include file for free variables
!
!     Warning!      MAXPAR_RES  in "param.inc" must always be of 
!                                 identical size to
!                 CHEM_MAX_NEIG in DISCUS/"config.inc"
!     This warning is obsolete.
!-
USE precision_mod
!
IMPLICIT NONE
PUBLIC
SAVE
!
INTEGER, PARAMETER :: MAXPAR     =  500
!  INTEGER, PARAMETER :: MAXPAR_RES = 12000
INTEGER            :: MAXPAR_RES = 6000
INTEGER            :: MAXPAR_REF =    1
!
INTEGER           , DIMENSION(0:MAXPAR)     :: inpara   = 0   ! (0:MAXPAR)
REAL(kind=PREC_DP), DIMENSION(0:MAXPAR)     :: rpara    = 0.0 ! (0:MAXPAR)
!  REAL   , DIMENSION(0:MAXPAR_RES) :: res_para = 0.0 ! (0:MAXPAR_RES)
REAL(KIND=PREC_DP), DIMENSION(:),ALLOCATABLE:: res_para       ! (0:MAXPAR_RES)
   REAL   , DIMENSION(:),ALLOCATABLE:: ref_para       ! Defined by DIFFEV
   REAL   , DIMENSION(0:MAXPAR)     :: kupl_para = 0.0 ! (0:MAXPAR)
   REAL   , DIMENSION(0:MAXPAR)     :: kupl_deriv= 0.0 ! (0:MAXPAR)
   INTEGER                          :: nrvalues = 0
   REAL   , DIMENSION(1:2, 0:15)    :: rvalues  = 0.0 ! (1:2)
   LOGICAL                          :: rvalue_yes = .false.
!
END MODULE param_mod
