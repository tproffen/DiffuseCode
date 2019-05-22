MODULE sine_table_mod
!
! Lookup tables for SINE and COSINE functions 
! Needed for a (reasonabley) fast explicit FOURIER
! The argument is in radians
!
USE precision_mod
!
!
INTEGER(KIND=PREC_INT_LARGE) , PARAMETER  :: ST_I2PI     = 2**16
INTEGER(KIND=PREC_INT_LARGE) , PARAMETER  :: ST_MASK     = ST_I2PI-1

REAL(KIND=PREC_DP), DIMENSION(:), ALLOCATABLE :: sine
REAL(KIND=PREC_DP), DIMENSION(:), ALLOCATABLE :: cosine
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE set_sine
!
USE wink_mod
!
IMPLICIT NONE
!
INTEGER :: i
!
IF(.NOT.ALLOCATED(sine)) THEN
!
   ALLOCATE(sine(0:ST_MASK))
!
   DO i=0, ST_MASK
     sine(i) = SIN(REAL(i*zpi/ST_I2PI,KIND(1.0D0)))
   ENDDO
ENDIF
!
END SUBROUTINE set_sine
!
!*******************************************************************************
!
SUBROUTINE set_cosine
!
USE wink_mod
!
IMPLICIT NONE
!
INTEGER :: i
!
IF(.NOT.ALLOCATED(cosine)) THEN
!
   ALLOCATE(cosine(0:ST_MASK))
!
   DO i=0, ST_MASK
     cosine(i) = COS(REAL(i*zpi/ST_I2PI,KIND(1.0D0)))
   ENDDO
ENDIF
!
END SUBROUTINE set_cosine
!
!*******************************************************************************
!
END MODULE sine_table_mod
