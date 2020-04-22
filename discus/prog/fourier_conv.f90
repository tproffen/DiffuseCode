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
USE errlist_mod
USE fftpack_mod
USE precision_mod
USE wink_mod
!
IMPLICIT NONE
!
REAL(KIND=PREC_DP), PARAMETER :: EPS = 1.0E-6
COMPLEX(KIND=KIND(0.0E0)) , DIMENSION(:,:), ALLOCATABLE  :: profile  ! the profile function
COMPLEX(KIND=KIND(0.0E0)) , DIMENSION(:,:), ALLOCATABLE  :: pattern  ! the diffraction pattern
COMPLEX(KIND=KIND(0.0E0)) , DIMENSION(:,:), ALLOCATABLE  :: temp     ! temporary pattern
REAL(KIND=PREC_SP), DIMENSION(:), ALLOCATABLE :: work   ! temporary array
REAL(KIND=PREC_SP), DIMENSION(:), ALLOCATABLE :: wsave  ! temporary array
!
INTEGER :: dimen
INTEGER :: i, j
INTEGER :: i1, j1
INTEGER :: n1, n2
INTEGER :: ii
INTEGER :: ipos
INTEGER                 :: lenwrk   ! Length of work array
INTEGER                 :: lensav   ! Length of wsave array
REAL(KIND=PREC_DP) :: rdimen
REAL(KIND=PREC_DP) :: dreal
REAL(KIND=PREC_DP) :: dimag
REAL(KIND=PREC_DP), DIMENSION(3  ) :: vector   ! Position = vimat x vector; vector [i,j,k]
REAL(KIND=PREC_DP), DIMENSION(3  ) :: posit    ! Position = vimat x vector; vector [i,j,k]
REAL(KIND=PREC_DP), DIMENSION(3  ) :: sigma    ! Sigma along columns of simat
!
IF(diff_res(1,1)<EPS .AND. diff_res(1,2)<EPS .AND. diff_res(1,3)<EPS) RETURN
dimen = MAX(num(1), num(2), num(3))
rdimen = REAL(dimen,KIND=4)
n1    = dimen
n2    = dimen
lenwrk = 2*n1*n2
lensav = 2*(n1*n2) + INT(LOG(REAL(n1))/LOG( 2.0E+00 ))       &
                   + INT(LOG(REAL(n2))/LOG( 2.0E+00 )) + 8.
ALLOCATE(profile(dimen, dimen))    
ALLOCATE(pattern(dimen, dimen))    
ALLOCATE(temp   (dimen, dimen))    
ALLOCATE(work(1:lenwrk))
ALLOCATE(wsave(1:lensav))
!
sigma(1) = diff_res(1, 1)
sigma(2) = diff_res(1, 2)
sigma(3) = diff_res(1, 3)
write(*,*) ' sigma ', sigma
!
ipos = 1
dimag = 0.0D0
vector(3) = 0.0
!
DO j=1, dimen
   j1 = j-ipos
   IF(j1>dimen/2) j1 = j1 - dimen
   vector(2) = REAL(j1, KIND=PREC_DP)
   DO i=1, dimen
      i1 = i-ipos
      IF(i1>dimen/2) i1 = i1 - dimen
      vector(1) = REAL(i1, KIND=PREC_DP)
      posit = MATMUL(diff_tr, vector)
      dreal = 1.0D0/SQRT(     ZPI)/sigma(1)*EXP(-0.50D0*(posit(1)         )**2) * &
              1.0D0/SQRT(     ZPI)/sigma(2)*EXP(-0.50D0*(posit(2)         )**2)
      profile(i,j) = COMPLEX(dreal, dimag)
   ENDDO
ENDDO
CALL output2(dimen, profile, 'INPUT/profile')
!
! Set profile
!
pattern = COMPLEX(0.0D0, 0.0D0)
ii = 0
DO j = 1, num (1) 
   j1 = MOD(j-1+dimen/2, dimen) + 1
   DO i = 1, num (2)
      i1 = MOD(i-1+dimen/2, dimen) + 1
      ii = ii + 1
      pattern(j1,i1) = csf(ii)
   ENDDO
ENDDO
CALL output2(dimen, pattern, 'INPUT/pattern')
!
CALL cfft2i (n1, n2, wsave, lensav, ier_num ) ! Initialize wsave
CALL cfft2f ( dimen, n1, n2, profile, wsave, lensav, work, lenwrk, ier_num )
CALL cfft2f ( dimen, n1, n2, pattern, wsave, lensav, work, lenwrk, ier_num )
CALL output2(dimen, profile, 'INPUT/profile_four')
CALL output2(dimen, pattern, 'INPUT/pattern_four')
!
temp = profile*pattern          ! Do pointwise multiplication of transforms
CALL output2(dimen, pattern, 'INPUT/temp_four')
CALL cfft2b ( dimen, n1, n2, temp   , wsave, lensav, work, lenwrk, ier_num )
CALL output2(dimen, pattern, 'INPUT/temp')
!
ii = 0
csf = COMPLEX(0.0D0, 0.0D0)
DO j = 1, num (1)
   j1 = MOD(j-1+dimen/2, dimen) + 1
   DO i = 1, num (2)
      i1 = MOD(i-1+dimen/2, dimen) + 1
      ii = ii + 1
      csf(ii) = temp(j1,i1)
   ENDDO
ENDDO
!
DEALLOCATE(profile)
DEALLOCATE(pattern)
DEALLOCATE(temp)
DEALLOCATE(work)
DEALLOCATE(wsave)
!
END SUBROUTINE four_conv        ! Convolute original diffraction pattern
!
!*******************************************************************************
!
SUBROUTINE output2(dimen, datum, base)
!
IMPLICIT NONE
!

INTEGER(KIND=4)                            , INTENT(IN) :: dimen
COMPLEX(KIND=4), DIMENSION(1:dimen,1:dimen), INTENT(IN) :: datum
CHARACTER(LEN=*)                           , INTENT(IN) :: base
!
INTEGER, PARAMETER:: IWR = 7
!
CHARACTER(len=1024) :: ofile
INTEGER (KIND=4) :: i, j, i1,j1
COMPLEX(KIND=4), DIMENSION(1:dimen) :: temp
!
ofile = base(1:LEN_TRIM(base)) // '.nipl.real'
OPEN(UNIT=IWR, FILE=ofile, STATUS='unknown')
WRITE(IWR, '(i4,2x,i4)') dimen, dimen
WRITE(IWR, '(a)') ' 0.00, 1.00,  0.00, 1.00'
DO j=1,dimen
  j1 = MOD(j-1+dimen/2, dimen) + 1
  DO i=1, dimen
    i1 = MOD(i-1+dimen/2, dimen) + 1
    temp(i) = datum(i1, j1)
  ENDDO
  WRITE(IWR, '(8(2x,F16.7))') (REAL(temp(i)), i=1, dimen)
!  WRITE(IWR, '(8(2x,F16.7))') (datum(i,j)*CONJG(datum(i,j)), j=1, dimen)
ENDDO
CLOSE(IWR)
!
ofile = base(1:LEN_TRIM(base)) // '.nipl.imag'
OPEN(UNIT=IWR, FILE=ofile, STATUS='unknown')
WRITE(IWR, '(i4,2x,i4)') dimen, dimen
WRITE(IWR, '(a)') ' 0.00, 1.00,  0.00, 1.00'
DO j=1,dimen
  j1 = MOD(j-1+dimen/2, dimen) + 1
  DO i=1, dimen
    i1 = MOD(i-1+dimen/2, dimen) + 1
    temp(i) = datum(i1, j1)
  ENDDO
  WRITE(IWR, '(8(2x,F16.7))') (IMAG(temp(i)), i=1, dimen)
!  WRITE(IWR, '(8(2x,F16.7))') (datum(i,j)*CONJG(datum(i,j)), j=1, dimen)
ENDDO
CLOSE(IWR)
!
ofile = base(1:LEN_TRIM(base)) // '.nipl.inte'
OPEN(UNIT=IWR, FILE=ofile, STATUS='unknown')
WRITE(IWR, '(i4,2x,i4)') dimen, dimen
WRITE(IWR, '(a)') ' 0.00, 1.00,  0.00, 1.00'
DO j=1,dimen
  j1 = MOD(j-1+dimen/2, dimen) + 1
  DO i=1, dimen
    i1 = MOD(i-1+dimen/2, dimen) + 1
    temp(i) = datum(i1, j1)
  ENDDO
  WRITE(IWR, '(8(2x,F16.7))') (REAL(temp(i))**2+IMAG(temp(i))**2, i=1, dimen)
!  WRITE(IWR, '(8(2x,F16.7))') (datum(i,j)*CONJG(datum(i,j)), j=1, dimen)
ENDDO
CLOSE(IWR)
!
ofile = base(1:LEN_TRIM(base)) // '.nipl.ampl'
OPEN(UNIT=IWR, FILE=ofile, STATUS='unknown')
WRITE(IWR, '(i4,2x,i4)') dimen, dimen
WRITE(IWR, '(a)') ' 0.00, 1.00,  0.00, 1.00'
DO j=1,dimen
  j1 = MOD(j-1+dimen/2, dimen) + 1
  DO i=1, dimen
    i1 = MOD(i-1+dimen/2, dimen) + 1
    temp(i) = datum(i1, j1)
  ENDDO
  WRITE(IWR, '(8(2x,F16.7))') (SQRT(REAL(temp(i))**2+IMAG(temp(i))**2), i=1, dimen)
!  WRITE(IWR, '(8(2x,F16.7))') (datum(i,j)*CONJG(datum(i,j)), j=1, dimen)
ENDDO
CLOSE(IWR)
!
ofile = base(1:LEN_TRIM(base)) // '.nipl.phase'
OPEN(UNIT=IWR, FILE=ofile, STATUS='unknown')
WRITE(IWR, '(i4,2x,i4)') dimen, dimen
WRITE(IWR, '(a)') ' 0.00, 1.00,  0.00, 1.00'
DO j=1,dimen
  j1 = MOD(j-1+dimen/2, dimen) + 1
  DO i=1, dimen
    i1 = MOD(i-1+dimen/2, dimen) + 1
    temp(i) = datum(i1, j1)
  ENDDO
  WRITE(IWR, '(8(2x,F16.7))') (ATAN2(IMAG(temp(i)),REAL(temp(i)))+3.14159, i=1, dimen)
!  WRITE(IWR, '(8(2x,F16.7))') (datum(i,j)*CONJG(datum(i,j)), j=1, dimen)
ENDDO
CLOSE(IWR)
!
END SUBROUTINE output2
!
!*******************************************************************************
!
END MODULE fourier_conv_mod
