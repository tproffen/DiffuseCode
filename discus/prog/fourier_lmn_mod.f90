MODULE fourier_lmn_mod
!
CONTAINS
SUBROUTINE fourier_lmn(eck,vi,inc,lmn,off)
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
USE tensors_mod
USE random_mod
!
IMPLICIT NONE
!
REAL   , DIMENSION(1:3, 1:4), INTENT(IN)  ::  eck
REAL   , DIMENSION(1:3, 1:3), INTENT(IN)  ::  vi
INTEGER, DIMENSION(1:3)     , INTENT(IN)  ::  inc
INTEGER, DIMENSION(1:6)     , INTENT(OUT) ::  lmn
REAL   , DIMENSION(1:3,1:3) , INTENT(OUT) :: off
!
INTEGER                     :: i,j, j1, j2
INTEGER                     :: dimen   ! dimension spanned by vi vectors (1, 2, 3)
INTEGER, DIMENSION(1:3)     :: direc   ! index of directions with inc > 1
LOGICAL, DIMENSION(1:3)     :: nooff   ! offset vectors
REAL                        :: dummy
REAL   , DIMENSION(1:3,1:3) :: mat_a
REAL   , DIMENSION(1:3)     :: vec_r
REAL   , DIMENSION(1:3,1:3) :: mat_i
REAL     :: ran1
!
mat_a(:,:) = 0.0
mat_i(:,:) = 0.0
dimen      = 0
nooff(:)   = .false.
off  (:,:) = 0.0
direc(:)   = 0
DO i=1,3
   IF(inc(i) > 1) THEN    ! reciprocal space has non_zero size along this increment
      dummy = SQRT(ABS(vi(1,i)**2 + vi(2,i)**2 +vi(3,i)**2))
      IF(dummy > 1.0E-7) THEN    ! Increment vector is non zero
         dimen = dimen + 1   ! Increment total dimension
         direc(dimen) = i    ! Keep track which increments are "thick"
         nooff(i) = .true.   ! No offest needed along this direction
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
   off(1,1) = 1.0
   off(2,2) = 1.0
   off(3,3) = 1.0
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
   off(:,j2) = off(:,j2)*0.001/dummy
   mat_a(:,j2) = off(:,j2)
ENDIF
!write(*,*) ' NUM , dimension ', inc, dimen
!write(*,*) ' DIREC, offsets  ', direc(:), nooff(:), j, j1, j2
!do j=1,3
!   write(*,1000) (mat_a(j,i),i=1,3),'   ',(mat_i(j,i),i=1,3), vi(j,1), vi(j,2),vi(j,3), &
!   (off(j,i),i=1,3)
!enddo
CALL invmat(mat_i, mat_a)
!do j=1,3
!   write(*,1000) (mat_a(j,i),i=1,3),'   ',(mat_i(j,i),i=1,3), vi(j,1), vi(j,2),vi(j,3), &
!   (off(j,i),i=1,3)
!enddo
!1000 FORMAT(3(f9.5,1x),a,3(f9.2,1x),2(2x, 3(f9.5,1x)))
lmn(:) = 0
vec_r(:) = 0.0
!write(*,*)
!write(*,*) ' eck(j,1) =   vi(j,1) *  lambda +   vi(j,2) *      my +   vi(j,3) *      ny'
DO j=1,3
   DO i=1,3
      vec_r(j) = vec_r(j) + mat_i(j,i)*eck(i,1)
   ENDDO
ENDDO
lmn(1:3) = NINT(vec_r(1:3))
DO j=1,3
   IF(.NOT.nooff(j)) THEN     ! We need an offset vector
      lmn(j+3) = lmn(j)
      lmn(j  ) = 0
   ENDIF
ENDDO
!do j=1,3
!   write(*,2000) eck(j,1), vi(j,1), lmn(1), vi(j,2), lmn(2), vi(j,3), lmn(3)
!enddo
!do j=1,3
!   write(*,3000)          off(j,1), lmn(4),off(j,2), lmn(5),off(j,3), lmn(6)
!enddo
!2000 FORMAT('(',f9.5,')=(',f9.5,')*',i8,' +(',f9.5,')*',i8,' +(',f9.5,')*',i8)
!3000 FORMAT(12x,       '(',f9.5,')*',i8,' +(',f9.5,')*',i8,' +(',f9.5,')*',i8)
!
END SUBROUTINE fourier_lmn
END MODULE fourier_lmn_mod
