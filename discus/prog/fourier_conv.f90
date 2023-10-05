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
real(kind=PREC_DP) :: weight, weight2
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
   dreal = 1.0D0/SQRT(     ZPI)/sigma(1)*EXP(-0.50D0*(posit(1)          )**2)
   profile(i) = CMPLX(dreal, dimag)
ENDDO
!
plan = fftw_plan_dft_1d(num(1)        , in_pattern, out_pattern, FFTW_FORWARD, FFTW_ESTIMATE)
!
call   fftw_execute_dft(plan, profile, temp)          ! FFT profile
temp = temp /                 (REAL(num(1),kind=PREC_DP))       ! Normalize Fourier of profile
!
IF(ilots.eq.LOT_OFF) THEN
   CALL maptofftfd(num, dsort, csf, in_pattern)       ! Use complex structure factor
   call   fftw_execute_dft(plan, in_pattern, out_pattern)
   out_pattern = temp     * out_pattern               ! Multiply the Fouriers
   call   fftw_execute_dft(plan, out_pattern,  in_pattern)
   in_pattern =  in_pattern/(real(num(1),kind=PREC_DP)) ! Scale for H 
   CALL mapfftfdtoline(num, dsort, csf, in_pattern)
ENDIF
!
weight = sum(dsi)
CALL maptofftfd(num, dsort, dsi, in_pattern)       ! Use intensities
call   fftw_execute_dft(plan, in_pattern, out_pattern)
out_pattern = temp  * out_pattern                 ! Multiply the Fouriers
call   fftw_execute_dft(plan, out_pattern,  in_pattern)
in_pattern =  in_pattern/(real(num(1),kind=PREC_DP)) ! Scale for H 
CALL mapfftfdtoline(num, dsort, dsi, in_pattern)
!
! Scale resulting pattern to maintain integral value
!
weight2 = sum(dsi)
csf = csf*weight/weight2
dsi = dsi*weight/weight2
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
real(kind=PREC_DP) :: weight, weight2
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
!
weight = sum(dsi)
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
! Scale resulting pattern to maintain integral value
!
weight2 = sum(dsi)
acsf=acsf*weight/weight2
csf = csf*weight/weight2
dsi = dsi*weight/weight2
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
real(kind=PREC_DP) :: weight, weight2
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
ALLOCATE(temp     (num(dsort(1)), num(dsort(2)), num(dsort(3))))    ! Allocate array for FFT
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
call   fftw_execute_dft(plan_p, profile_t, temp)
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
weight = sum(dsi)
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
! Scale resulting pattern to maintain integral value
!
weight2 = sum(dsi)
acsf=acsf*weight/weight2
csf = csf*weight/weight2
dsi = dsi*weight/weight2
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
SUBROUTINE four_resolution(line, lp)
!-
!  Interpret the optional parameters on the run command to set the resolution
!+
!
USE crystal_mod
USE diffuse_mod
USE get_params_mod
USE matrix_mod
USE metric_mod
use lib_metric_mod
USE take_param_mod
USE precision_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: line
INTEGER         , INTENT(INOUT) :: lp
!
REAL(KIND=PREC_DP) :: EPS = 1.0D-7
INTEGER, PARAMETER :: MAXW = 4
!                                                                       
INTEGER, PARAMETER :: NOPTIONAL = 3
INTEGER, PARAMETER :: O_SIGABS  = 1                  ! Current phase number
INTEGER, PARAMETER :: O_SIGORD  = 2                  ! Weight fraction
INTEGER, PARAMETER :: O_SIGTOP  = 3                  ! Single / multiple
CHARACTER(LEN=   6), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=MAX(PREC_STRING,LEN(line))), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent!opt. para is present
REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 0 ! Number of values to calculate 
!
CHARACTER(LEN=MAX(PREC_STRING,LEN(line))), DIMENSION(MAXW) :: cpara
INTEGER            , DIMENSION(MAXW) :: lpara
REAL(KIND=PREC_DP) , DIMENSION(MAXW) :: werte
INTEGER             :: ianz
INTEGER             :: iianz
integer             :: isdim
integer             :: j
!
real(kind=PREC_DP)                 :: r        ! A real number
REAL(KIND=PREC_DP), DIMENSION(3,3) :: vimat    ! Col wise increment vectors in layer
REAL(KIND=PREC_DP), DIMENSION(3,3) :: simat    ! Col
REAL(KIND=PREC_DP), DIMENSION(3,3) :: siinv    ! simat^-1
!REAL(KIND=PREC_DP), DIMENSION(3,3) :: trmat    ! trmat = siinv x vimat
REAL(KIND=PREC_DP), DIMENSION(3  ) :: u        ! Sigma along columns of simat
REAL(KIND=PREC_DP), DIMENSION(3)   :: uu       ! Length of user si vectors
REAL(KIND=PREC_DP), DIMENSION(3  ) :: v1       ! Dummy vectors
REAL(KIND=PREC_DP), DIMENSION(3  ) :: v2       ! Dummy vectors
REAL(KIND=PREC_DP), DIMENSION(3  ) :: v3       ! Dummy vectors
!
DATA oname  / 'sigabs', 'sigord', 'sigtop' /
DATA loname /  6      ,  6      ,  6       /
!
diff_res = 0.0D0
CALL get_params (line, ianz, cpara, lpara, maxw, lp)
IF (ier_num.ne.0) RETURN
IF(ianz==0) RETURN        ! No params, use  default NULL matrix
!
CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
      oname, loname, opara, lopara, lpresent, owerte)
!
!  Determine dimensions that we need 
!
isdim = 3
IF(inc(1)==1) isdim = isdim - 1
IF(inc(2)==1) isdim = isdim - 1
IF(inc(3)==1) isdim = isdim - 1
!
CALL four_res_optional(lpresent(O_SIGABS), 1, MAXW, opara(O_SIGABS), &
       lopara(O_SIGABS), werte, iianz)
if(isdim>1) then
   CALL four_res_optional(lpresent(O_SIGORD), 2, MAXW, opara(O_SIGORD), &
          lopara(O_SIGORD), werte, iianz)
endif
if(isdim==3) then
   CALL four_res_optional(lpresent(O_SIGTOP), 3, MAXW, opara(O_SIGTOP), &
          lopara(O_SIGTOP), werte, iianz)
endif
!
! Build simat
!
!
u = diff_res(2:4,1)
uu(1) = SQRT (skalpro (u, u, cr_rten) )
IF(uu(1)< EPS) THEN
diff_res(2:4,1) = vi(1:3,1)        ! User vector was NULL, use abscissa
u = diff_res(2:4,1)
uu(1) = SQRT (skalpro (u, u, cr_rten) )
ENDIF
!
if(isdim==1) then                  ! 1D diffraction pattern ; generate ord and top
   diff_res(2:4,2) = 0.0D0
   diff_res(2:4,3) = 0.0D0
   diff_res(1  ,2) = 0.0001D0
   diff_res(1  ,3) = 0.0001D0
   v1 = diff_res(2:4,1)
   j = minloc(abs(v1), 1)
   uu(1) = SQRT (skalpro (v1, v1, cr_rten) )
   if(v1(j)<1.0D-8) then
      diff_res(j+1  ,2) = 1.0D0
      uu(2) = 1.0D0
      v1(j) = 1.0D9
      j = minloc(abs(v1), 1)
      if(v1(j)<1.0D-8) then
         diff_res(j+1  ,3) = 1.0D0
         uu(3) = 1.0D0
      else
         v1 = diff_res(2:4,1)
         v2 = diff_res(2:4,2)
         call vekprod(v1, v2, v3, cr_reps, cr_gten)   ! Build a vector normal to sigabs and sigord
         r = lib_blen(cr_rten, v3)
         diff_res(2:4,3) = v3/r
         u = diff_res(2:4,3)
         uu(3) = SQRT (skalpro (u, u, cr_rten) )
      endif
   else     ! No zeros in abscissa
!      
      v1 = diff_res(2:4,1)
      call random_number(r)
      v2(1) = 2.0D0*r - 1.0D0
      call random_number(r)
      v2(2) = 2.0D0*r - 1.0D0
      call random_number(r)
      v2(3) = 2.0D0*r - 1.0D0
      call vekprod(v1, v2, v3, cr_reps, cr_gten)   ! Build a vector normal to sigabs and sigord
      r = lib_blen(cr_rten, v3)
!write(*,*) ' VECTORS ', v1
!write(*,*) ' VECTORS ', v2
!write(*,*) ' VECTORS ', v3, '    ', r
!write(*,*) ' ANGLE 1, 3 ', lib_bang(cr_rten, v1, v3)
      v3(1) = v3(1)/r
      v3(2) = v3(2)/r
      v3(3) = v3(3)/r
!write(*,*) ' VECTORS ', v3
!write(*,*) ' ANGLE 1, 3 ', lib_bang(cr_rten, v1, v3)
      diff_res(2:4,2) = v3
      v2 = diff_res(2:4,2)
      call vekprod(v1, v2, v3, cr_reps, cr_gten)   ! Build a vector normal to sigabs and sigord
      r = lib_blen(cr_rten, v3)
      diff_res(2:4,3) = v3/r
!     uu(1) = 1.0D0
      uu(2) = 1.0D0
      uu(3) = 1.0D0
      diff_res(1  ,2) = 0.0001D0
      diff_res(1  ,3) = 0.0001D0
   endif
!write(*,*) ' DIFF_RES ', diff_res(:,1)
!write(*,*) ' DIFF_RES ', diff_res(:,2)
!write(*,*) ' DIFF_RES ', diff_res(:,3)
!
else
!
   u = diff_res(2:4,2)
   uu(2) = SQRT (skalpro (u, u, cr_rten) )
   IF(uu(2)< EPS) THEN
      diff_res(2:4,2) = vi(1:3,2)        ! User vector was NULL, use ordinate
      u = diff_res(2:4,2)
      uu(2) = SQRT (skalpro (u, u, cr_rten) )
      diff_res(1  ,2) = 0.001
   ENDIF
!
   u = diff_res(2:4,3)
   uu(3) = SQRT (skalpro (u, u, cr_rten) )
   IF(uu(3)< EPS) THEN
      v1 = diff_res(2:4,1)
      v2 = diff_res(2:4,2)
      CALL vekprod(v1, v2, v3, cr_reps, cr_gten)   ! Build a vecor normal to sigabs and sigord
      uu(3) = SQRT (skalpro (v3, v3, cr_rten) )
      diff_res(2:4,3) = v3
      diff_res(1  ,3) = 0.001
   ENDIF
endif
!
simat(1:3, 1) = diff_res(2:4,1) * diff_res(1,1) / uu(1)
!
simat(1:3, 2) = diff_res(2:4,2) * diff_res(1,2) / uu(2)
!
simat(1:3, 3) = diff_res(2:4,3) * diff_res(1,3) / uu(3)
!write(*,*) ' SIMAT ', simat(:,1)
!write(*,*) ' SIMAT ', simat(:,2)
!write(*,*) ' SIMAT ', simat(:,3)
!
CALL matinv3(simat, siinv) 
!
! Build vimat
!
if(isdim==1) then
   vimat(:,1) = vi(:,1)
   vimat(:,2) = diff_res(2:4,2)
   vimat(:,3) = diff_res(2:4,3)
else
   vimat = vi
endif
!write(*,*) ' viMAT ', vimat(:,1)
!write(*,*) ' viMAT ', vimat(:,2)
!write(*,*) ' viMAT ', vimat(:,3)
diff_tr = MATMUL(siinv, vimat)
!write(*,*) 'diff_tr', diff_tr(:,1)
!write(*,*) 'diff_tr', diff_tr(:,2)
!write(*,*) 'diff_tr', diff_tr(:,3)
!
END SUBROUTINE four_resolution
!
!*****7*****************************************************************
!
SUBROUTINE four_res_optional(lpresent, ientry, MAXW, opara, &
           lopara, werte, ianz)
!
USE diffuse_mod
!
use errlist_mod
USE get_params_mod
USE take_param_mod
USE precision_mod
!
IMPLICIT NONE
!
LOGICAL                            , INTENT(IN)    :: lpresent ! Optional parameter was present 
INTEGER                            , INTENT(IN)    :: ientry   ! Entry number
INTEGER                            , INTENT(IN)    :: MAXW     ! Dimension of werte
CHARACTER(LEN=*)                   , INTENT(INOUT) :: opara    ! The string with optional values
INTEGER                            , INTENT(INOUT) :: lopara   ! length of string
REAL(KIND=PREC_DP), DIMENSION(MAXW), INTENT(OUT)   :: werte    ! Numerical values
INTEGER                            , INTENT(OUT)   :: ianz    ! Number of numerical values

ianz = 0
IF(lpresent) THEN                ! Sigmais present
CALL get_optional_multi(MAXW, opara, lopara, werte, ianz)
IF(ier_num==0) THEN
IF(ianz>0) THEN
 diff_res(1,ientry) =werte(1)
 IF(ianz==4) THEN
    diff_res(2:4,ientry) = werte(2:4)
 ELSE
    diff_res(2,ientry) = 0.0D0
    diff_res(3,ientry) = 0.0D0
    diff_res(4,ientry) = 0.0D0
 ENDIF
ELSE
 diff_res(:,ientry) = 0.00D0
ENDIF
ELSE
RETURN
ENDIF
ELSE
diff_res(:,ientry) = 0.0D0
ENDIF
!
!
END SUBROUTINE four_res_optional
!
!*******************************************************************************
!
END MODULE fourier_conv_mod
