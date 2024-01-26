MODULE fourier_lmn_mod
!
CONTAINS
SUBROUTINE fourier_lmn(eck_in,vi_in,inc,lmn,off_in)
!
!   Calculate integer "phase" values for the lower left corner
!   in reciprocal space, which are used in four_strucf.
!   This calculation ensures that the the phases at hkl and -hkl
!   are identical.
!
!   If a rod or plane does not intersect hkl=0, offset phases 
!   are calculated.
!
USE errlist_mod
use matrix_mod
USE lib_random_func
use precision_mod
USE random_mod
!
IMPLICIT NONE
!
REAL(kind=PREC_DP), DIMENSION(1:3, 1:4), INTENT(IN)  ::  eck_in  ! Corners in reciprocal space
REAL(kind=PREC_DP), DIMENSION(1:3, 1:3), INTENT(IN)  ::  vi_in   ! increment vectors
INTEGER           , DIMENSION(1:3)     , INTENT(IN)  ::  inc   ! Number data points
INTEGER           , DIMENSION(1:6)     , INTENT(OUT) ::  lmn   ! eck(:,1) is at lmn*vi
REAL(kind=PREC_DP), DIMENSION(1:3,1:3) , INTENT(OUT) :: off_in   ! Additional phase shift vector
!
REAL(kind=PREC_QP), DIMENSION(1:3, 1:4)              ::  eck   ! Corners in reciprocal space
REAL(kind=PREC_QP), DIMENSION(1:3, 1:3)              ::  vi    ! increment vectors
REAL(kind=PREC_QP), DIMENSION(1:3,1:3)               :: off    ! Additional phase shift vector
!INTEGER, DIMENSION(1:3)     ::  opq
INTEGER                     :: i,j, j1, j2
INTEGER                     :: dimen   ! dimension spanned by vi vectors (1, 2, 3)
INTEGER, DIMENSION(1:3)     :: direc   ! index of directions with inc > 1
LOGICAL, DIMENSION(1:3)     :: nooff   ! offset vectors
REAL(kind=PREC_QP)          :: dummy
REAL(kind=PREC_QP)   , DIMENSION(1:3,1:3) :: mat_a
REAL(kind=PREC_QP)   , DIMENSION(1:3)     :: vec_r
!REAL(kind=PREC_QP)   , DIMENSION(1:3)     :: vec_s
REAL(kind=PREC_QP)   , DIMENSION(1:3)     :: vec_t
REAL(kind=PREC_QP)   , DIMENSION(1:3,1:3) :: mat_i
!character(len=1) :: cdummy
real(kind=PREC_QP) :: big=10000000.0_PREC_QP
!
!write(*,*) big
eck = real(real(nint(eck_in*big), kind=PREC_QP)/big, kind=PREC_QP)
vi  = real(real(nint( vi_in*big), kind=PREC_QP)/big, kind=PREC_QP)
!write(*,*) ' LMN CALCULATION '
!write(*,*) ' vi in1 ',vi_in(:,1) , ' | ', 1.0D0
!write(*,*) ' vi in2 ',vi_in(:,2) 
!write(*,*) ' vi in3 ',vi_in(:,3) 
!write(*,*) ' vi     ',real(nint(vi_in(:,1)*1.0D8)/1.0D8 , kind=PREC_DP)
!write(*,*) ' vi     ',nint(vi_in(:,2)*1.0D8)
!write(*,*) ' vi     ',nint(vi_in(:,3)*1.0D8) 
!write(*,*) ' vi     ',vi   (:,1) 
!write(*,*) ' vi     ',vi   (:,2) 
!write(*,*) ' vi     ',vi   (:,3) 
!write(*,*) ' vi     ',real(vi   (:,1) , kind=PREC_DP)
!write(*,*) ' vi     ',real(vi   (:,2) , kind=PREC_DP)
!write(*,*) ' vi     ',real(vi   (:,3) , kind=PREC_DP)
!write(*,*) ' ECK_IN ', eck_in(:,1)
!write(*,*) ' ECK    ', eck   (:,1)
!
mat_a(:,:) = 0.0d0
mat_i(:,:) = 0.0d0
dimen      = 0
nooff(:)   = .false.
off  (:,:) = 0.0D0
direc(:)   = 0
DO i=1,3
   IF(inc(i) > 1) THEN    ! reciprocal space has non_zero size along this increment
      dummy = SQRT(ABS(vi(1,i)**2 + vi(2,i)**2 +vi(3,i)**2))
      IF(dummy > 1.0E-7) THEN    ! Increment vector is non zero
         dimen = dimen + 1   ! Increment total dimension
         direc(dimen) = i    ! Keep track which increments are "thick"
         nooff(i) = .true.   ! No offset needed along this direction
         DO j=1,3
            mat_a(j,i) = vi(j,i)  ! Store increment vector
         ENDDO
      ELSE
         ier_num = -132
         ier_typ = ER_APPL
         WRITE(ier_msg(1),'(a,i2,a,i8)') 'No. of points at axis',i,' is ',inc(1)
         WRITE(ier_msg(2),'(a,3f8.4)') 'Zero increment vec:',vi(:,i)
         WRITE(ier_msg(3),'(a      )') 'Check corners ll,lr,ul,tl and na,no,nt'
         RETURN
      ENDIF
   ENDIF
ENDDO
IF(dimen==0) THEN       ! Single spot in reciprocal space
   off(1,1) = 1.0D0
   off(2,2) = 1.0D0
   off(3,3) = 1.0D0
   mat_a = off
ELSEIF(dimen==1) THEN   ! 1D line in reciprocal space
   j  = direc(1)
   j1 = MOD(j  ,3) + 1 ! Cyclically increment j by one
   j2 = MOD(j+1,3) + 1 ! Cyclically increment j by two
   off(1, j2) = eck(2,1) *vi(3,j) - eck(3,1) *vi(2,j)  ! Vector product of indices
   off(2, j2) = eck(3,1) *vi(1,j) - eck(1,1) *vi(3,j)  ! Vector product of indices
   off(3, j2) = eck(1,1) *vi(2,j) - eck(2,1) *vi(1,j)  ! Vector product of indices
   dummy = SQRT(ABS(off(1,j2)**2+off(2,j2)**2+off(3,j2)**2))
   IF(dummy > 1.E-7) THEN           ! Nonzero offset
      off(:,j2) = off(:,j2)*0.001/dummy
      off(1, j1) = off(2,j2)*vi(3,j) - off(3,j2)*vi(2,j) 
      off(2, j1) = off(3,j2)*vi(1,j) - off(1,j2)*vi(3,j) 
      off(3, j1) = off(1,j2)*vi(2,j) - off(2,j2)*vi(1,j) 
!
      dummy = SQRT(ABS(off(1,j1)**2+off(2,j1)**2+off(3,j1)**2))
      IF(dummy > 1.E-7) THEN           ! Nonzero offset
         off(:,j1) = off(:,j1)*0.001/dummy                     ! scale to 0.001
      ELSE                             ! Should never occur ?
         off(1,j1) = ran1(idum)
         off(2,j1) = ran1(idum)
         off(3,j1) = ran1(idum)
      ENDIF
   ELSE                          ! Zero offset == line through origin
      off(1,j1) = ran1(idum)
      off(2,j1) = ran1(idum)
      off(3,j1) = ran1(idum)
      off(1,j2) = ran1(idum)
      off(2,j2) = ran1(idum)
      off(3,j2) = ran1(idum)
   ENDIF
   mat_a(:,j1) = off(:,j1)
   mat_a(:,j2) = off(:,j2)
ELSEIF(dimen==2) THEN  ! 2D layer in reciprocal space
   j  = direc(1)
   j1 = direc(2)
   j2 = 3 - MOD(j+j1,3)  ! last index 1,2=>3 1,3=>2 2,3=>1
   off(1, j2) = vi(2,j)*vi(3,j1) - vi(3,j)*vi(2,j1)
   off(2, j2) = vi(3,j)*vi(1,j1) - vi(1,j)*vi(3,j1)
   off(3, j2) = vi(1,j)*vi(2,j1) - vi(2,j)*vi(1,j1)
   dummy = SQRT(off(1,j2)**2+off(2,j2)**2+off(3,j2)**2)
   off(:,j2) = off(:,j2)*0.001D0/dummy
   mat_a(:,j2) = off(:,j2)
ENDIF
!write(*,*) ' NUM , dimension ', inc, dimen
!write(*,*) ' DIREC, offsets  ', direc(:), nooff(:), j, j1, j2
!do j=1,3
!   write(*,1000) (mat_a(j,i),i=1,3),'   ',(mat_i(j,i),i=1,3), vi(j,1), vi(j,2),vi(j,3), &
!   (off(j,i),i=1,3)
!enddo
!
call matinv_q(mat_a, mat_i)
!write(*,*) ' vi     ',vi   (:,1) 
!write(*,*) ' vi     ',vi   (:,2) 
!write(*,*) ' vi     ',vi   (:,3) 
!write(*,*) ' mat  1 ',mat_a(1,:) 
!write(*,*) ' mat  2 ',mat_a(2,:) 
!write(*,*) ' mat  3 ',mat_a(3,:) 
!write(*,*) ' mat_i1 ',mat_i(1,:) 
!write(*,*) ' mat_i2 ',mat_i(2,:) 
!write(*,*) ' mat_i3 ',mat_i(3,:) 
!do j=1,3
!   write(*,1000) (mat_a(j,i),i=1,3),'   ',(mat_i(j,i),i=1,3), vi(j,1), vi(j,2),vi(j,3), &
!   (off(j,i),i=1,3)
!enddo
!1000 FORMAT(3(f9.5,1x),a,3(f9.5,1x),2(2x, 3(f9.5,1x)))
!
lmn(:) = 0
!vec_r = 0.0D0
!vec_r(1) = mat_i(1,1)*eck(1,1) + mat_i(1,2)*eck(2,1) + mat_i(1,3)*eck(3,1)
!vec_r(2) = mat_i(2,1)*eck(1,1) + mat_i(3,2)*eck(2,1) + mat_i(2,3)*eck(3,1)
!vec_r(3) = mat_i(3,1)*eck(1,1) + mat_i(2,2)*eck(2,1) + mat_i(3,3)*eck(3,1)
vec_r = matmul(mat_i, eck(:,1))
!write(*,*) ' vec_r ', vec_r
if(abs(vec_r(1)-real(nint(vec_r(1)),kind=PREC_QP))<1.0D-9) vec_r(1) = real(nint(vec_r(1)),kind=PREC_QP)
if(abs(vec_r(2)-real(nint(vec_r(2)),kind=PREC_QP))<1.0D-9) vec_r(2) = real(nint(vec_r(2)),kind=PREC_QP)
if(abs(vec_r(3)-real(nint(vec_r(3)),kind=PREC_QP))<1.0D-9) vec_r(3) = real(nint(vec_r(3)),kind=PREC_QP)
!write(*,*) ' vec_r ', vec_r
!lmn(1:3) = NINT(vec_r(1:3))
!
lmn(1:3) =  int(vec_r(1:3))        ! lmn(1:3) is integer part only, frac goes into off
vec_t = matmul(mat_a, vec_r - int(vec_r))   ! Fractional part of vi's
!write(*,*) ' vec_t ', vec_t
!
DO j=1,3
   IF(.NOT.nooff(j)) THEN     ! We need an offset vector
      lmn(j+3) = lmn(j)
      lmn(j  ) = 0
   else                       ! No offset, just regular phase shift;
!                             !  eck(:,1) is at non-integer multiple of vi's
      off(j,j) = off(j,j) + vec_t(j)
      lmn(3+j) = 1
      if(abs(off(j,j))<1.0D-9) then
         off(j,j) = 0.0D0
         lmn(3+j) = 0
      endif
   ENDIF
ENDDO
!write(*,*) ' OFF ', off(1,1), nooff(1)
!write(*,*) ' OFF ', off(2,2), nooff(2)
!write(*,*) ' OFF ', off(3,3), nooff(3)
!write(*,*)
!write(*,*) ' eck(j,1) =   vi(j,1) *  lambda +   vi(j,2) *      my +   vi(j,3) *      ny'
!do j=1,3
!   write(*,2000) eck(j,1), vi(j,1), lmn(1), vi(j,2), lmn(2), vi(j,3), lmn(3)
!enddo
!do j=1,3
!   write(*,3000)          off(j,1), lmn(4),off(j,2), lmn(5),off(j,3), lmn(6)
!enddo
!2000 FORMAT('(',f9.5,')=(',f9.5,')*',i8,' +(',f9.5,')*',i8,' +(',f9.5,')*',i8)
!3000 FORMAT(12x,       '(',f9.5,')*',i8,' +(',f9.5,')*',i8,' +(',f9.5,')*',i8)
!read(*,'(a)') cdummy
off_in = real(off, kind=PREC_DP)

END SUBROUTINE fourier_lmn
END MODULE fourier_lmn_mod
