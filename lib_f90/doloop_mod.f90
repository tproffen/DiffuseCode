MODULE doloop_mod
!+
!     This file contains Variable definitions for do loop's, if's
!-
USE precision_mod
!
IMPLICIT NONE
!
PUBLIC
SAVE
!
INTEGER, PARAMETER           :: MAXLEV = 20
INTEGER, PARAMETER           :: MAXCOM = 2000
!
INTEGER, DIMENSION(0:MAXLEV) :: nloop    = 0       ! (0:MAXLEV)
INTEGER, DIMENSION(0:MAXLEV) :: iloop    = 0       ! (0:MAXLEV)
INTEGER, DIMENSION(1:1     ) :: do_kpara = 0       ! (1)
!
LOGICAL, DIMENSION(0:MAXLEV) :: ldostart = .false. ! (0:MAXLEV)
!
REAL(KIND=PREC_DP), DIMENSION(0:MAXLEV) :: ghigh    = 0.0     ! (0:MAXLEV)
REAL(KIND=PREC_DP), DIMENSION(0:MAXLEV) :: glow     = 0.0     ! (0:MAXLEV)
REAL(KIND=PREC_DP), DIMENSION(0:MAXLEV) :: ginc     = 1.0     ! (0:MAXLEV)
!
END MODULE doloop_mod
