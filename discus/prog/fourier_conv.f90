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
use lib_f90_fftw3
USE map_1dtofield
USE precision_mod
USE wink_mod
!
use iso_c_binding
!
IMPLICIT NONE
!
INTEGER, DIMENSION(3), INTENT(IN) :: dsort
!
COMPLEX(KIND=C_DOUBLE_COMPLEX) , DIMENSION(:), ALLOCATABLE  :: profile  ! the profile function
COMPLEX(KIND=C_DOUBLE_COMPLEX) , DIMENSION(:), ALLOCATABLE  ::  in_pattern  ! the diffraction pattern
COMPLEX(KIND=C_DOUBLE_COMPLEX) , DIMENSION(:), ALLOCATABLE  :: out_pattern  ! the diffraction pattern
COMPLEX(KIND=C_DOUBLE_COMPLEX) , DIMENSION(:), ALLOCATABLE  :: temp     ! A temporary pattern
REAL(KIND=PREC_DP) :: dreal
REAL(KIND=PREC_DP) :: dimag
REAL(KIND=PREC_DP), DIMENSION(3  ) :: vector   ! Position = vimat x vector; vector [i,j,k]
REAL(KIND=PREC_DP), DIMENSION(3  ) :: posit    ! Position = vimat x vector; vector [i,j,k]
REAL(KIND=PREC_DP), DIMENSION(3  ) :: sigma    ! Sigma along columns of simat
!
INTEGER :: i, i1
INTEGER :: ipos
!
type(c_ptr) :: plan    ! FFWT3 plan
!
ALLOCATE(profile(num(1) ))        ! Allocate profile in regular dimensions
ALLOCATE(temp   (num(1) ))        ! Allocate profile in regular dimensions
ALLOCATE( in_pattern(num(dsort(1)) )) ! Allocated pattern in sorted dimensions
ALLOCATE(out_pattern(num(dsort(1)) )) ! Allocated pattern in sorted dimensions
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
   profile(i) = CMPLX(dreal, dimag)
ENDDO
!
plan = fftw_plan_dft_1d(num(1)        , in_pattern, out_pattern, FFTW_FORWARD, FFTW_ESTIMATE)
!
!profile   = fft(profile  ) / SQRT(REAL(num(1)))    ! FFT profile
call   fftw_execute_dft(plan, profile, temp)          ! FFT profile
temp = temp / SQRT(REAL(num(1)))                      ! Normalize Fourier of profile
!
!
IF(ilots.eq.LOT_OFF) THEN
   CALL maptofftfd(num, dsort, csf, in_pattern)       ! Use complex structure factor
!   pattern   = fft(pattern)   / SQRT(REAL(num(1)))    ! FFT pattern
   call   fftw_execute_dft(plan, in_pattern, out_pattern)
   out_pattern = temp     * out_pattern               ! Multiply the Fouriers
!  temp      = fft(temp)      / SQRT(REAL(num(1)))    ! FFT multiplied pattern
   call   fftw_execute_dft(plan, out_pattern,  in_pattern)
   in_pattern =  in_pattern/(real(num(1),kind=PREC_DP)) ! Scale for H 
   CALL mapfftfdtoline(num, dsort, csf, in_pattern)
ENDIF
!
CALL maptofftfd(num, dsort, dsi, in_pattern)       ! Use intensities
call   fftw_execute_dft(plan, in_pattern, out_pattern)
out_pattern = temp  * out_pattern                 ! Multiply the Fouriers
call   fftw_execute_dft(plan, out_pattern,  in_pattern)
in_pattern =  in_pattern/(real(num(1),kind=PREC_DP)) ! Scale for H 
CALL mapfftfdtoline(num, dsort, dsi, in_pattern)
!
call   fftw_destroy_plan(plan)
!
DEALLOCATE( in_pattern)
DEALLOCATE(out_pattern)
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
use lib_f90_fftw3
USE map_1dtofield
USE precision_mod
USE wink_mod
!
use iso_c_binding
!
IMPLICIT NONE
!
INTEGER, DIMENSION(3), INTENT(IN) :: dsort
!
COMPLEX(KIND=C_DOUBLE_COMPLEX) , DIMENSION(:,:), ALLOCATABLE  :: profile  ! the profile function
COMPLEX(KIND=C_DOUBLE_COMPLEX) , DIMENSION(:,:), ALLOCATABLE  :: profile_t! the profile function
COMPLEX(KIND=C_DOUBLE_COMPLEX) , DIMENSION(:,:), ALLOCATABLE  :: in_pattern  ! the diffraction pattern
COMPLEX(KIND=C_DOUBLE_COMPLEX) , DIMENSION(:,:), ALLOCATABLE  :: out_pattern  ! the diffraction pattern
COMPLEX(KIND=C_DOUBLE_COMPLEX) , DIMENSION(:,:), ALLOCATABLE  :: temp     ! A temporary pattern
REAL(KIND=PREC_DP) :: dreal
REAL(KIND=PREC_DP) :: dimag
REAL(KIND=PREC_DP), DIMENSION(3  ) :: vector   ! Position = vimat x vector; vector [i,j,k]
REAL(KIND=PREC_DP), DIMENSION(3  ) :: posit    ! Position = vimat x vector; vector [i,j,k]
REAL(KIND=PREC_DP), DIMENSION(3  ) :: sigma    ! Sigma along columns of simat
!
INTEGER :: i,j, i1, j1
INTEGER :: ipos
!
type(c_ptr) :: plan    ! FFWT3 plan
type(c_ptr) :: plan_p  ! FFWT3 plan for PROFILE
!
ALLOCATE(profile(num(1), num(2) ))               ! Allocate profile in regular dimensions
ALLOCATE( in_pattern(num(dsort(1)), num(dsort(2)) )) ! Allocated pattern in sorted dimensions
ALLOCATE(out_pattern(num(dsort(1)), num(dsort(2)) )) ! Allocated pattern in sorted dimensions
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
      profile(i,j) = CMPLX(dreal, dimag)
   ENDDO
ENDDO
!
ALLOCATE(profile_t(num(dsort(1)), num(dsort(2))))     ! Allocate array for FFT
ALLOCATE(temp     (num(dsort(1)), num(dsort(2))))     ! Allocate array for FFT
IF(dsort(1)==1) THEN   ! 1st coordinate is largest
   profile_t = profile
ELSE                   ! 2nd coordinate is largest, transpose profile
   profile_t = TRANSPOSE(profile)
ENDIF
!
plan_p = fftw_plan_dft_2d(num(dsort(2)), num(dsort(1)), profile_t, temp, FFTW_FORWARD, FFTW_ESTIMATE)  ! Plan for PROFILE
plan   = fftw_plan_dft_2d(num(dsort(2)), num(dsort(1)), in_pattern, out_pattern, FFTW_FORWARD, FFTW_ESTIMATE)  ! Plan for diffraction
!profile_t = fft(profile_t) / SQRT(REAL(num(1)*num(2)))    ! FFT profile
call   fftw_execute_dft(plan_p, profile_t, temp)
temp = temp/sqrt(real(num(dsort(1))*num(dsort(2)),kind=PREC_DP))
!
!
IF(ilots.eq.LOT_OFF) THEN
   CALL maptofftfd(num, dsort, csf, in_pattern)              ! Use complex structure factor
!   pattern   = fft(pattern)   / SQRT(REAL(num(1)*num(2)))    ! FFT pattern
   call   fftw_execute_dft(plan, in_pattern, out_pattern)
   out_pattern      = temp*out_pattern                       ! Multiply the Fouriers
!   temp      = fft(temp)      / SQRT(REAL(num(1)*num(2)))   ! FFT multiplied pattern
   call   fftw_execute_dft(plan, out_pattern, in_pattern)    ! FFT multiplied pattern
   in_pattern = in_pattern / (REAL(num(1)*num(2)))
   CALL mapfftfdtoline(num, dsort, csf, in_pattern)          ! Map back into complex structure factor
ENDIF
!
CALL maptofftfd(num, dsort, dsi, in_pattern)              ! Use intensities
!pattern   = fft(pattern)   / SQRT(REAL(num(1)*num(2)))    ! FFT pattern
!temp      = profile_t*pattern                             ! Multiply the Fouriers
!temp      = fft(temp)      / SQRT(REAL(num(1)*num(2)))    ! FFT multiplied pattern
call   fftw_execute_dft(plan, in_pattern, out_pattern)    ! FFT Diffraction pattern
out_pattern      = temp*out_pattern                       ! Multiply the Fouriers
call   fftw_execute_dft(plan, out_pattern, in_pattern)    ! FFT multiplied pattern
in_pattern = in_pattern / (REAL(num(1)*num(2)))
CALL mapfftfdtoline(num, dsort, dsi, in_pattern)          ! Restore convoluted intensities
!
!
CALL maptofftfd(num, dsort, acsf, in_pattern)             ! Use average complex structure factor
!pattern   = fft(pattern)   / SQRT(REAL(num(1)*num(2)))    ! FFT pattern
!temp      = profile_t*pattern                             ! Multiply the Fouriers
!temp      = fft(temp)      / SQRT(REAL(num(1)*num(2)))    ! FFT multiplied pattern
call   fftw_execute_dft(plan, in_pattern, out_pattern)    ! FFT Diffraction pattern
out_pattern      = temp*out_pattern                       ! Multiply the Fouriers
call   fftw_execute_dft(plan, out_pattern, in_pattern)    ! FFT multiplied pattern
in_pattern = in_pattern / (REAL(num(1)*num(2)))
CALL mapfftfdtoline(num, dsort, acsf, in_pattern)         ! Restore convoluted average complex structure factor
!
DEALLOCATE(in_pattern)
DEALLOCATE(out_pattern)
DEALLOCATE(temp)
DEALLOCATE(profile_t)
DEALLOCATE(profile)
call   fftw_destroy_plan(plan  )
call   fftw_destroy_plan(plan_p)
!
END SUBROUTINE four_conv_2D     ! Convolute original diffraction pattern
!
!*******************************************************************************
!
SUBROUTINE four_conv_3D(dsort)     ! Convolute original diffraction pattern
!
USE diffuse_mod
!
use lib_f90_fftw3
USE map_1dtofield
USE precision_mod
!USE singleton
USE wink_mod
!
use iso_c_binding
!
IMPLICIT NONE
!
INTEGER, DIMENSION(3), INTENT(IN) :: dsort
!
COMPLEX(KIND=C_DOUBLE_COMPLEX) , DIMENSION(:,:,:), ALLOCATABLE  :: profile  ! the profile function
COMPLEX(KIND=C_DOUBLE_COMPLEX) , DIMENSION(:,:,:), ALLOCATABLE  :: profile_t! the profile function
COMPLEX(KIND=C_DOUBLE_COMPLEX) , DIMENSION(:,:,:), ALLOCATABLE  ::  in_pattern  ! the diffraction pattern
COMPLEX(KIND=C_DOUBLE_COMPLEX) , DIMENSION(:,:,:), ALLOCATABLE  :: out_pattern  ! the diffraction pattern
COMPLEX(KIND=C_DOUBLE_COMPLEX) , DIMENSION(:,:,:), ALLOCATABLE  :: temp     ! A temporary pattern
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
type(c_ptr) :: plan    ! FFWT3 plan
type(c_ptr) :: plan_p  ! FFWT3 plan for PROFILE
!
ALLOCATE(profile(num(1), num(2), num(3) ))                     ! Allocate profile in regular dimensions
ALLOCATE(in_pattern(num(dsort(1)), num(dsort(2)), num(dsort(3)))) ! Allocated pattern in sorted dimensions
ALLOCATE(out_pattern(num(dsort(1)), num(dsort(2)), num(dsort(3)))) ! Allocated pattern in sorted dimensions
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
         profile(i,j,k) = CMPLX(dreal, dimag)
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
!profile_t = fft(profile_t) / SQRT(REAL(num(1)*num(2)*num(3)))    ! FFT profile
plan_p = fftw_plan_dft_3d(num(dsort(3)), num(dsort(2)), num(dsort(1)), profile_t, temp, FFTW_FORWARD, FFTW_ESTIMATE)  ! Plan for PROFILE
plan   = fftw_plan_dft_3d(num(dsort(3)), num(dsort(2)), num(dsort(1)), in_pattern, out_pattern, FFTW_FORWARD, FFTW_ESTIMATE)  ! Plan for diffraction
temp = temp/sqrt(real(num(dsort(1))*num(dsort(2))*num(dsort(3)),kind=PREC_DP))
!
IF(ilots.eq.LOT_OFF) THEN
   CALL maptofftfd(num, dsort, csf, in_pattern)              ! Use complex structure factor
!  pattern   = fft(pattern)   / SQRT(REAL(num(1)*num(2)*num(3)))    ! FFT pattern
!  temp      = profile_t*pattern                                    ! Multiply the Fouriers
!  temp      = fft(temp)      / SQRT(REAL(num(1)*num(2)*num(3)))    ! FFT multiplied pattern
   call   fftw_execute_dft(plan, in_pattern, out_pattern)
   out_pattern      = temp*out_pattern                       ! Multiply the Fouriers
   call   fftw_execute_dft(plan, out_pattern, in_pattern)    ! FFT multiplied pattern
   in_pattern = in_pattern / (REAL(num(1)*num(2)*num(3)))
   CALL mapfftfdtoline(num, dsort, csf, in_pattern)
ENDIF
!
CALL maptofftfd(num, dsort, dsi, in_pattern)              ! Use intensities
!pattern   = fft(pattern)   / SQRT(REAL(num(1)*num(2)*num(3)))    ! FFT pattern
!temp      = profile_t*pattern                                    ! Multiply the Fouriers
!temp      = fft(temp)      / SQRT(REAL(num(1)*num(2)*num(3)))    ! FFT multiplied pattern
call   fftw_execute_dft(plan, in_pattern, out_pattern)    ! FFT pattern
out_pattern      = temp*out_pattern                       ! Multiply the Fouriers
call   fftw_execute_dft(plan, out_pattern, in_pattern)    ! FFT multiplied pattern
in_pattern = in_pattern / (REAL(num(1)*num(2)*num(3)))
CALL mapfftfdtoline(num, dsort, dsi, in_pattern)
!
DEALLOCATE( in_pattern)
DEALLOCATE(out_pattern)
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
