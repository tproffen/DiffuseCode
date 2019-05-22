MODULE discus_fft_mod
!
! Routines to calculate a fourier for 3D PDF
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE do_fft_2d_cos(npkt1, npkt2, zwrt, out_eck, out_vi, out_inc)
!
USE precision_mod
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: npkt1
INTEGER, INTENT(IN) :: npkt2
REAL(KIND=PREC_SP), DIMENSION(npkt1,npkt2), INTENT(IN) :: zwrt
REAL(KIND=PREC_SP), DIMENSION(3, 4)       , INTENT(IN) :: out_eck
REAL(KIND=PREC_SP), DIMENSION(3, 3)       , INTENT(IN) :: out_vi
INTEGER           , DIMENSION(3)          , INTENT(IN) :: out_inc
!
!
END SUBROUTINE do_fft_2d_cos
!
!*******************************************************************************
!
END MODULE discus_fft_mod
