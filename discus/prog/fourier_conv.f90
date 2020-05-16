MODULE fourier_conv_mod
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE four_conv        ! Convolute original diffraction pattern
!-
!   Convolutes the original single crystal pattern with a "resolution" function
!+
USE diffuse_mod
!
IMPLICIT NONE
REAL(KIND=PREC_DP), PARAMETER :: EPS = 1.0E-6
INTEGER :: isdim
INTEGER, DIMENSION(3) :: dsort
!
IF(diff_res(1,1)<EPS .AND. diff_res(1,2)<EPS .AND. diff_res(1,3)<EPS) RETURN
!
!  Determine sequence of array dimensions 
!
dsort(1)      = MAXLOC(num, 1)
num(dsort(1)) = -num(dsort(1))
dsort(2)      = MAXLOC(num, 1)
num(dsort(2)) = -num(dsort(2))
dsort(3)      = MAXLOC(num, 1)
num(dsort(3)) = -num(dsort(3))
num = -num
!
!  Determine dimensions that we need 
!
isdim = 3
IF(num(1)==1) isdim = isdim - 1
IF(num(2)==1) isdim = isdim - 1
IF(num(3)==1) isdim = isdim - 1
!
IF(isdim==3) THEN
   CALL four_conv_3D(dsort)
ELSEIF(isdim==2) THEN
   CALL four_conv_2D(dsort)
ELSEIF(isdim==1) THEN
   CALL four_conv_1D(dsort)
ENDIF
!
END SUBROUTINE four_conv        ! Convolute original diffraction pattern
!
!*******************************************************************************
!
SUBROUTINE four_conv_1D(dsort)     ! Convolute original diffraction pattern
!
USE diffuse_mod
!
USE map_1dtofield
USE precision_mod
USE singleton
USE wink_mod
!
IMPLICIT NONE
!
INTEGER, DIMENSION(3), INTENT(IN) :: dsort
!
COMPLEX(KIND=PREC_DP) , DIMENSION(:), ALLOCATABLE  :: profile  ! the profile function
COMPLEX(KIND=PREC_DP) , DIMENSION(:), ALLOCATABLE  :: pattern  ! the diffraction pattern
COMPLEX(KIND=PREC_DP) , DIMENSION(:), ALLOCATABLE  :: temp     ! A temporary pattern
REAL(KIND=PREC_DP) :: dreal
REAL(KIND=PREC_DP) :: dimag
REAL(KIND=PREC_DP), DIMENSION(3  ) :: vector   ! Position = vimat x vector; vector [i,j,k]
REAL(KIND=PREC_DP), DIMENSION(3  ) :: posit    ! Position = vimat x vector; vector [i,j,k]
REAL(KIND=PREC_DP), DIMENSION(3  ) :: sigma    ! Sigma along columns of simat
!
INTEGER :: i, i1
INTEGER :: ipos
!
ALLOCATE(profile(num(1) ))        ! Allocate profile in regular dimensions
ALLOCATE(pattern(num(dsort(1)) )) ! Allocated pattern in sorted dimensions
!
!  Build profile function
!
sigma(1) = diff_res(1, 1)
sigma(2) = diff_res(1, 2)
sigma(3) = diff_res(1, 3)
!
ipos = 1
dimag = 0.0D0
vector(3) = 0.0
vector(2) = 0.0
!
DO i=1, num(1)
   i1 = i-ipos
   IF(i1>num(1)/2) i1 = i1 - num(1)
   vector(1) = REAL(i1, KIND=PREC_DP)
   posit = MATMUL(diff_tr, vector)
   dreal = 1.0D0/SQRT(     ZPI)/sigma(1)*EXP(-0.50D0*(posit(1) )**2)
   profile(i) = COMPLEX(dreal, dimag)
ENDDO
!
profile   = fft(profile  ) / SQRT(REAL(num(1)))    ! FFT profile
!
IF(ilots.eq.LOT_OFF) THEN
   CALL maptofftfd(num, dsort, csf, pattern)          ! Use complex structure factor
   pattern   = fft(pattern)   / SQRT(REAL(num(1)))    ! FFT pattern
   temp      = profile  *pattern                      ! Multiply the Fouriers
   temp      = fft(temp)      / SQRT(REAL(num(1)))    ! FFT multiplied pattern
   CALL mapfftfdtoline(num, dsort, csf, temp)
ENDIF
!
CALL maptofftfd(num, dsort, dsi, pattern)          ! Use intensities
pattern   = fft(pattern)   / SQRT(REAL(num(1)))    ! FFT pattern
temp      = profile  *pattern                      ! Multiply the Fouriers
temp      = fft(temp)      / SQRT(REAL(num(1)))    ! FFT multiplied pattern
CALL mapfftfdtoline(num, dsort, dsi, temp)
!
DEALLOCATE(pattern)
DEALLOCATE(temp)
DEALLOCATE(profile)
!
END SUBROUTINE four_conv_1D     ! Convolute original diffraction pattern
!
!*******************************************************************************
!
SUBROUTINE four_conv_2D(dsort)     ! Convolute original diffraction pattern
!
USE diffuse_mod
!
USE map_1dtofield
USE precision_mod
USE singleton
USE wink_mod
!
IMPLICIT NONE
!
INTEGER, DIMENSION(3), INTENT(IN) :: dsort
!
COMPLEX(KIND=PREC_DP) , DIMENSION(:,:), ALLOCATABLE  :: profile  ! the profile function
COMPLEX(KIND=PREC_DP) , DIMENSION(:,:), ALLOCATABLE  :: profile_t! the profile function
COMPLEX(KIND=PREC_DP) , DIMENSION(:,:), ALLOCATABLE  :: pattern  ! the diffraction pattern
COMPLEX(KIND=PREC_DP) , DIMENSION(:,:), ALLOCATABLE  :: temp     ! A temporary pattern
REAL(KIND=PREC_DP) :: dreal
REAL(KIND=PREC_DP) :: dimag
REAL(KIND=PREC_DP), DIMENSION(3  ) :: vector   ! Position = vimat x vector; vector [i,j,k]
REAL(KIND=PREC_DP), DIMENSION(3  ) :: posit    ! Position = vimat x vector; vector [i,j,k]
REAL(KIND=PREC_DP), DIMENSION(3  ) :: sigma    ! Sigma along columns of simat
!
INTEGER :: i,j, i1, j1
INTEGER :: ipos
!
ALLOCATE(profile(num(1), num(2) ))               ! Allocate profile in regular dimensions
ALLOCATE(pattern(num(dsort(1)), num(dsort(2)) )) ! Allocated pattern in sorted dimensions
!
!  Build profile function
!
sigma(1) = diff_res(1, 1)
sigma(2) = diff_res(1, 2)
sigma(3) = diff_res(1, 3)
!
ipos = 1
dimag = 0.0D0
vector(3) = 0.0
!
DO j=1, num(2)
   j1 = j-ipos
   IF(j1>num(2)/2) j1 = j1 - num(2)
   vector(2) = REAL(j1, KIND=PREC_DP)
   DO i=1, num(1)
      i1 = i-ipos
      IF(i1>num(1)/2) i1 = i1 - num(1)
      vector(1) = REAL(i1, KIND=PREC_DP)
      posit = MATMUL(diff_tr, vector)
      dreal = 1.0D0/SQRT(     ZPI)/sigma(1)*EXP(-0.50D0*(posit(1) )**2) * &
              1.0D0/SQRT(     ZPI)/sigma(2)*EXP(-0.50D0*(posit(2) )**2)
      profile(i,j) = COMPLEX(dreal, dimag)
   ENDDO
ENDDO
!
ALLOCATE(profile_t(num(dsort(1)), num(dsort(2))))     ! Allocate array for FFT
IF(dsort(1)==1) THEN   ! 1st coordinate is largest
   profile_t = profile
ELSE                   ! 2nd coordinate is largest, transpose profile
   profile_t = TRANSPOSE(profile)
ENDIF
!
profile_t = fft(profile_t) / SQRT(REAL(num(1)*num(2)))    ! FFT profile
!
IF(ilots.eq.LOT_OFF) THEN
   CALL maptofftfd(num, dsort, csf, pattern)                 ! Use complex structure factor
   pattern   = fft(pattern)   / SQRT(REAL(num(1)*num(2)))    ! FFT pattern
   temp      = profile_t*pattern                             ! Multiply the Fouriers
   temp      = fft(temp)      / SQRT(REAL(num(1)*num(2)))    ! FFT multiplied pattern
   CALL mapfftfdtoline(num, dsort, csf, temp)
ENDIF
!
CALL maptofftfd(num, dsort, dsi, pattern)                 ! Use intensities
pattern   = fft(pattern)   / SQRT(REAL(num(1)*num(2)))    ! FFT pattern
temp      = profile_t*pattern                             ! Multiply the Fouriers
temp      = fft(temp)      / SQRT(REAL(num(1)*num(2)))    ! FFT multiplied pattern
CALL mapfftfdtoline(num, dsort, dsi, temp)                ! Resore convoluted intensities
!
DEALLOCATE(pattern)
DEALLOCATE(temp)
DEALLOCATE(profile_t)
DEALLOCATE(profile)
!
END SUBROUTINE four_conv_2D     ! Convolute original diffraction pattern
!
!*******************************************************************************
!
SUBROUTINE four_conv_3D(dsort)     ! Convolute original diffraction pattern
!
USE diffuse_mod
!
USE map_1dtofield
USE precision_mod
USE singleton
USE wink_mod
!
IMPLICIT NONE
!
INTEGER, DIMENSION(3), INTENT(IN) :: dsort
!
COMPLEX(KIND=PREC_DP) , DIMENSION(:,:,:), ALLOCATABLE  :: profile  ! the profile function
COMPLEX(KIND=PREC_DP) , DIMENSION(:,:,:), ALLOCATABLE  :: profile_t! the profile function
COMPLEX(KIND=PREC_DP) , DIMENSION(:,:,:), ALLOCATABLE  :: pattern  ! the diffraction pattern
COMPLEX(KIND=PREC_DP) , DIMENSION(:,:,:), ALLOCATABLE  :: temp     ! A temporary pattern
REAL(KIND=PREC_DP) :: dreal
REAL(KIND=PREC_DP) :: dimag
REAL(KIND=PREC_DP), DIMENSION(3  ) :: vector   ! Position = vimat x vector; vector [i,j,k]
REAL(KIND=PREC_DP), DIMENSION(3  ) :: posit    ! Position = vimat x vector; vector [i,j,k]
REAL(KIND=PREC_DP), DIMENSION(3  ) :: sigma    ! Sigma along columns of simat
!
INTEGER :: i,j, k, i1, j1, k1
INTEGER, DIMENSION(3) :: ientry
INTEGER :: ipos
!
ALLOCATE(profile(num(1), num(2), num(3) ))                     ! Allocate profile in regular dimensions
ALLOCATE(pattern(num(dsort(1)), num(dsort(2)), num(dsort(3)))) ! Allocated pattern in sorted dimensions
!
!  Build profile function
!
sigma(1) = diff_res(1, 1)
sigma(2) = diff_res(1, 2)
sigma(3) = diff_res(1, 3)
!
ipos = 1
dimag = 0.0D0
!
DO k=1, num(3)
   k1 = k-ipos
   IF(k1>num(3)/2) k1 = k1 - num(3)
   vector(3) = REAL(k1, KIND=PREC_DP)
   DO j=1, num(2)
      j1 = j-ipos
      IF(j1>num(2)/2) j1 = j1 - num(2)
      vector(2) = REAL(j1, KIND=PREC_DP)
      DO i=1, num(1)
         i1 = i-ipos
         IF(i1>num(1)/2) i1 = i1 - num(1)
         vector(1) = REAL(i1, KIND=PREC_DP)
         posit = MATMUL(diff_tr, vector)
         dreal = 1.0D0/SQRT(     ZPI)/sigma(1)*EXP(-0.50D0*(posit(1) )**2) * &
                 1.0D0/SQRT(     ZPI)/sigma(2)*EXP(-0.50D0*(posit(2) )**2) * &
                 1.0D0/SQRT(     ZPI)/sigma(3)*EXP(-0.50D0*(posit(3) )**2)
         profile(i,j,k) = COMPLEX(dreal, dimag)
      ENDDO
   ENDDO
ENDDO
!
ALLOCATE(profile_t(num(dsort(1)), num(dsort(2)), num(dsort(3))))    ! Allocate array for FFT
!
DO k=1, num(3)
   ientry(dsort(3)) = k
   DO j=1, num(2)
      ientry(dsort(2)) = j
      DO i=1, num(1)
         ientry(dsort(1)) = i
         profile_t(ientry(1), ientry(2), ientry(3)) = profile(i,j,k)
      ENDDO
   ENDDO
ENDDO
!
profile_t = fft(profile_t) / SQRT(REAL(num(1)*num(2)*num(3)))    ! FFT profile
!
IF(ilots.eq.LOT_OFF) THEN
   CALL maptofftfd(num, dsort, csf, pattern)              ! Use complex structure factor
   pattern   = fft(pattern)   / SQRT(REAL(num(1)*num(2)*num(3)))    ! FFT pattern
   temp      = profile_t*pattern                                    ! Multiply the Fouriers
   temp      = fft(temp)      / SQRT(REAL(num(1)*num(2)*num(3)))    ! FFT multiplied pattern
   CALL mapfftfdtoline(num, dsort, csf, temp)
ENDIF
!
CALL maptofftfd(num, dsort, dsi, pattern)              ! Use intensities
pattern   = fft(pattern)   / SQRT(REAL(num(1)*num(2)*num(3)))    ! FFT pattern
temp      = profile_t*pattern                                    ! Multiply the Fouriers
temp      = fft(temp)      / SQRT(REAL(num(1)*num(2)*num(3)))    ! FFT multiplied pattern
CALL mapfftfdtoline(num, dsort, dsi, temp)
!
DEALLOCATE(pattern)
DEALLOCATE(temp)
DEALLOCATE(profile_t)
DEALLOCATE(profile)
!
END SUBROUTINE four_conv_3D     ! Convolute original diffraction pattern
!
!*******************************************************************************
!
!
!*******************************************************************************
!
END MODULE fourier_conv_mod
