MODULE discus_fft_mod
!
! Routines to calculate a fourier for 3D PDF
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE do_fft_2d_cos(npkt1, npkt2, zwrt, out_eck, out_vi, out_inc, &
                         nnew1, nnew2, znew, pdf3d_eck, pdf3d_vi, pdf3d_inc)
!
USE precision_mod
USE sine_table_mod
USE wink_mod
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: npkt1
INTEGER, INTENT(IN) :: npkt2
REAL(KIND=PREC_SP), DIMENSION(npkt1,npkt2), INTENT(IN) :: zwrt
REAL(KIND=PREC_SP), DIMENSION(3, 4)       , INTENT(IN) :: out_eck
REAL(KIND=PREC_SP), DIMENSION(3, 3)       , INTENT(IN) :: out_vi
INTEGER           , DIMENSION(3)          , INTENT(IN) :: out_inc
INTEGER                                   , INTENT(IN) :: nnew1
INTEGER                                   , INTENT(IN) :: nnew2
REAL(KIND=PREC_SP), DIMENSION(nnew1,nnew2), INTENT(OUT) :: znew
REAL(KIND=PREC_SP), DIMENSION(3, 4)       , INTENT(IN) :: pdf3d_eck
REAL(KIND=PREC_SP), DIMENSION(3, 3)       , INTENT(IN) :: pdf3d_vi
INTEGER           , DIMENSION(3)          , INTENT(IN) :: pdf3d_inc
!
!REAL(KIND=PREC_SP), DIMENSION(3,npkt1) :: hwrt
!REAL(KIND=PREC_SP), DIMENSION(3,npkt2) :: kwrt
!
INTEGER :: i,j, ih,ik
REAL(KIND=PREC_DP) :: xx
REAL(KIND=PREC_DP) :: yy
REAL(KIND=PREC_DP) :: hh,kk
REAL(KIND=PREC_DP) :: arg
!
CALL set_cosine
!
DO j=0,nnew2-1
   DO i=0,nnew1-1
      xx = pdf3d_eck(1,1) + (i   )*pdf3d_vi(1,1) + (i   )*pdf3d_vi(1,2)
      yy = pdf3d_eck(2,1) + (j   )*pdf3d_vi(2,1) + (j   )*pdf3d_vi(2,2)
      DO ik=1,NPKT2
         DO ih=1,NPKT1
            hh = out_eck(1,1) + (ih-1)*out_vi(1,1) + (ik-1)*out_vi(1,2)
            kk = out_eck(2,1) + (ih-1)*out_vi(2,1) + (ik-1)*out_vi(2,2)
            arg = ABS(hh*xx + kk*yy)
            znew(i+1,j+1) = znew(i+1,j+1) +                              &
                zwrt(ih,ik)*cosine(MOD(INT( arg*REAL(ST_MASK,KIND=PREC_DP), KIND=PREC_INT_LARGE), ST_MASK))
         ENDDO
      ENDDO
   ENDDO
ENDDO
!
END SUBROUTINE do_fft_2d_cos
!
!*******************************************************************************
!
END MODULE discus_fft_mod
