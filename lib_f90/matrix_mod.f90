MODULE matrix_mod
!
USE errlist_mod
IMPLICIT NONE
!
!interface deter
!  module procedure det2, det3, det4
!end interface deter
!!
!interface matinv
!  module procedure matinv2, matinv3, matinv4
!end interface matinv
!
CONTAINS
!
!*******************************************************************************
!
PURE FUNCTION determinant(a) RESULT(det)
!
USE precision_mod
REAL(KIND=PREC_DP), DIMENSION(:,:), INTENT(IN) :: a
!
REAL(KIND=PREC_DP)                             :: det
!
det = 0.0D0
if(    ubound(a,1)==2) then
  det = det2(a)
elseif(ubound(a,1)==3) then
  det = det3(a)
elseif(ubound(a,1)==4) then
  det = det4(a)
endif
!
END FUNCTION determinant
!
!*******************************************************************************
!
PURE FUNCTION det2(a) RESULT(det)
!
USE precision_mod
REAL(KIND=PREC_DP), DIMENSION(2,2), INTENT(IN) :: a
!
REAL(KIND=PREC_DP)                             :: det
!
det = (A(1,1)*A(2,2) - A(1,2)*A(2,1))
!
END FUNCTION det2
!
!*******************************************************************************
!
!PURE FUNCTION det3(a) RESULT(det)
PURE FUNCTION det3(a) RESULT(det)
!
USE precision_mod
REAL(KIND=PREC_DP), DIMENSION(3,3), INTENT(IN) :: a
!
REAL(KIND=PREC_DP)                             :: det
!
det  =   (A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2) &
        - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1) &
        + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))
!
END FUNCTION det3
!
!*******************************************************************************
!
PURE FUNCTION det4(a) RESULT(det)
!
USE precision_mod
REAL(KIND=PREC_DP), DIMENSION(4,4), INTENT(IN) :: a
!
REAL(KIND=PREC_DP)                             :: det
!
det = &
    (A(1,1)*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))&
   - A(1,2)*(A(2,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))&
   + A(1,3)*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))&
   - A(1,4)*(A(2,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(2,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(2,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1))))
!
END FUNCTION det4
!
!*******************************************************************************
!
subroutine matinv(mat, imat)
!-
! Performs matrix inversion mat => imat
! Generic interface to the 2x2, 3x2, 4x4 versions
!+
!
use precision_mod
!
implicit none
real(kind=PREC_DP), dimension(:,:), intent(in) :: mat      !! Matrix
real(kind=PREC_DP), dimension(:,:), intent(out):: imat     !! Inverse matrix
!
imat = 0.0D0
!
if(    ubound(mat,1)==2) then
  call  matinv2(mat, imat)
elseif(ubound(mat,1)==3) then
  call  matinv3(mat, imat)
elseif(ubound(mat,1)==4) then
  call  matinv4(mat, imat)
endif
!
end subroutine matinv
!
!*******************************************************************************
!
SUBROUTINE matinv2(A, B)
    !! Performs a direct calculation of the inverse of a 2×2 matrix.
USE precision_mod
REAL(KIND=PREC_DP), DIMENSION(2,2), INTENT(IN) :: A        !! Matrix
REAL(KIND=PREC_DP), DIMENSION(2,2), INTENT(OUT):: B        !! Inverse matrix
!
REAL(KIND=PREC_DP)                        :: det, detinv
!
! Calculate the inverse determinant of the matrix
det = det2(a)
!
IF(det/=0.0D0) THEN
   detinv = 1./det

    ! Calculate the inverse of the matrix
    B(1,1) = +detinv * A(2,2)
    B(2,1) = -detinv * A(2,1)
    B(1,2) = -detinv * A(1,2)
    B(2,2) = +detinv * A(1,1)
ELSE
   ier_num = -45
   ier_typ = ER_FORT
ENDIF
!
END SUBROUTINE matinv2
!
!*******************************************************************************
!
SUBROUTINE matinv3(A, B)
    !! Performs a direct calculation of the inverse of a 3×3 matrix.
USE precision_mod
REAL(KIND=PREC_DP), DIMENSION(3,3), INTENT(IN) :: A        !! Matrix
REAL(KIND=PREC_DP), DIMENSION(3,3), INTENT(OUT):: B        !! Inverse matrix
!
REAL(KIND=PREC_DP)                             :: det, detinv
!
! Calculate the inverse determinant of the matrix
det = det3(a)
!
IF(det/=0.0D0) THEN
   detinv = 1./det
!
! Calculate the inverse of the matrix
   B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
   B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
   B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
   B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
   B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
   B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
   B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
   B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
   B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
ELSE
   ier_num = -45
   ier_typ = ER_FORT
ENDIF
!
END SUBROUTINE matinv3
!
!*******************************************************************************
!
SUBROUTINE matinv4(A, B)
    !! Performs a direct calculation of the inverse of a 4×4 matrix.
USE precision_mod
REAL(KIND=PREC_DP), DIMENSION(4,4), INTENT(IN) :: A        !! Matrix
REAL(KIND=PREC_DP), DIMENSION(4,4), INTENT(OUT):: B        !! Inverse matrix
!
REAL(KIND=PREC_DP)                             :: det, detinv
!
! Calculate the inverse determinant of the matrix
det = det4(a)
!
IF(det/=0.0D0) THEN
   detinv = 1./det
    ! Calculate the inverse of the matrix
   B(1,1) = detinv*(A(2,2)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(2,3)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(2,4)*(A(3,2)*A(4,3)-A(3,3)*A(4,2)))
   B(2,1) = detinv*(A(2,1)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(2,3)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(2,4)*(A(3,3)*A(4,1)-A(3,1)*A(4,3)))
   B(3,1) = detinv*(A(2,1)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(2,2)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(2,4)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
   B(4,1) = detinv*(A(2,1)*(A(3,3)*A(4,2)-A(3,2)*A(4,3))+A(2,2)*(A(3,1)*A(4,3)-A(3,3)*A(4,1))+A(2,3)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
   B(1,2) = detinv*(A(1,2)*(A(3,4)*A(4,3)-A(3,3)*A(4,4))+A(1,3)*(A(3,2)*A(4,4)-A(3,4)*A(4,2))+A(1,4)*(A(3,3)*A(4,2)-A(3,2)*A(4,3)))
   B(2,2) = detinv*(A(1,1)*(A(3,3)*A(4,4)-A(3,4)*A(4,3))+A(1,3)*(A(3,4)*A(4,1)-A(3,1)*A(4,4))+A(1,4)*(A(3,1)*A(4,3)-A(3,3)*A(4,1)))
   B(3,2) = detinv*(A(1,1)*(A(3,4)*A(4,2)-A(3,2)*A(4,4))+A(1,2)*(A(3,1)*A(4,4)-A(3,4)*A(4,1))+A(1,4)*(A(3,2)*A(4,1)-A(3,1)*A(4,2)))
   B(4,2) = detinv*(A(1,1)*(A(3,2)*A(4,3)-A(3,3)*A(4,2))+A(1,2)*(A(3,3)*A(4,1)-A(3,1)*A(4,3))+A(1,3)*(A(3,1)*A(4,2)-A(3,2)*A(4,1)))
   B(1,3) = detinv*(A(1,2)*(A(2,3)*A(4,4)-A(2,4)*A(4,3))+A(1,3)*(A(2,4)*A(4,2)-A(2,2)*A(4,4))+A(1,4)*(A(2,2)*A(4,3)-A(2,3)*A(4,2)))
   B(2,3) = detinv*(A(1,1)*(A(2,4)*A(4,3)-A(2,3)*A(4,4))+A(1,3)*(A(2,1)*A(4,4)-A(2,4)*A(4,1))+A(1,4)*(A(2,3)*A(4,1)-A(2,1)*A(4,3)))
   B(3,3) = detinv*(A(1,1)*(A(2,2)*A(4,4)-A(2,4)*A(4,2))+A(1,2)*(A(2,4)*A(4,1)-A(2,1)*A(4,4))+A(1,4)*(A(2,1)*A(4,2)-A(2,2)*A(4,1)))
   B(4,3) = detinv*(A(1,1)*(A(2,3)*A(4,2)-A(2,2)*A(4,3))+A(1,2)*(A(2,1)*A(4,3)-A(2,3)*A(4,1))+A(1,3)*(A(2,2)*A(4,1)-A(2,1)*A(4,2)))
   B(1,4) = detinv*(A(1,2)*(A(2,4)*A(3,3)-A(2,3)*A(3,4))+A(1,3)*(A(2,2)*A(3,4)-A(2,4)*A(3,2))+A(1,4)*(A(2,3)*A(3,2)-A(2,2)*A(3,3)))
   B(2,4) = detinv*(A(1,1)*(A(2,3)*A(3,4)-A(2,4)*A(3,3))+A(1,3)*(A(2,4)*A(3,1)-A(2,1)*A(3,4))+A(1,4)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)))
   B(3,4) = detinv*(A(1,1)*(A(2,4)*A(3,2)-A(2,2)*A(3,4))+A(1,2)*(A(2,1)*A(3,4)-A(2,4)*A(3,1))+A(1,4)*(A(2,2)*A(3,1)-A(2,1)*A(3,2)))
   B(4,4) = detinv*(A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))+A(1,2)*(A(2,3)*A(3,1)-A(2,1)*A(3,3))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1)))
ELSE
   ier_num = -45
   ier_typ = ER_FORT
ENDIF
!
END SUBROUTINE matinv4
!
SUBROUTINE matinvn(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
USE precision_mod
!
IMPLICIT NONE
!
INTEGER                           , INTENT(IN)    :: n
REAL(KIND=PREC_DP), DIMENSION(n,n), INTENT(INOUT) :: a
REAL(KIND=PREC_DP), DIMENSION(n,n), INTENT(OUT)   :: c
!double precision a(n,n), c(n,n)
!double precision L(n,n), U(n,n), b(n), d(n), x(n)
!double precision coeff
REAL(KIND=PREC_DP), DIMENSION(n,n) :: L 
REAL(KIND=PREC_DP), DIMENSION(n,n) :: U 
REAL(KIND=PREC_DP), DIMENSION(n)   :: b 
REAL(KIND=PREC_DP), DIMENSION(n)   :: d 
REAL(KIND=PREC_DP), DIMENSION(n)   :: x 
REAL(KIND=PREC_DP)                 :: coeff 
INTEGER :: i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do
end SUBROUTINE matinvn
!
!*******************************************************************************
!
SUBROUTINE mat_axis(mat, axis, alpha, det, trace, ier)
!
!  Determine the rotation angle and the rotation axis for a 3x3 matrix
!  RETURNS: 
!  axis(3):    Rotation axis
!  alpha  :    Rotation angle
!  det    :    Determinant
!  trace  :    Trace
!  ier    :    0:success
!             -1:det==0
!             -2:det/=+-1
!             -3:trace outside [-3:3]
!             +1: Rotation angle is zero
!
USE precision_mod
USE trig_degree_mod
IMPLICIT NONE
!
REAL(KIND=PREC_DP), DIMENSION(3,3), INTENT(IN)  :: mat
REAL(KIND=PREC_DP), DIMENSION(3)  , INTENT(OUT) :: axis
REAL(KIND=PREC_DP)                , INTENT(OUT) :: alpha
REAL(KIND=PREC_DP)                , INTENT(OUT) :: det
REAL(KIND=PREC_DP)                , INTENT(OUT) :: trace
INTEGER             , INTENT(OUT) :: ier
!
REAL(KIND=PREC_DP), PARAMETER :: EPS = 1.D-5
!
INTEGER              :: i1,i2,i3
REAL(KIND=PREC_DP), DIMENSION(3,3) :: work
REAL(KIND=PREC_DP)                 :: ca, sa
REAL(KIND=PREC_DP)                 :: n1,n2,n3
REAL(KIND=PREC_DP)                 :: m1,m2,m3
!
axis(:) = 0.0D0
alpha   = 0.0D0
det     = det3(mat)
trace   = mat(1,1) + mat(2,2) + mat(3,3)
ca      = 0.5D0*(trace-1)
IF(ABS(det)<EPS) THEN
   ier = -1
   RETURN
ENDIF
IF(ABS(det)-1>EPS) THEN
   ier = -2
   RETURN
ENDIF
IF(ABS(trace)>3.0D0) THEN
   ier = -3
   RETURN
ENDIF
work  = mat/det                   ! Normalize the matrix
alpha = ACOSD(ca)
sa    = SIND(alpha)
!
IF(ABS(alpha)<EPS) THEN
   ier =  1
   RETURN
ENDIF
!
IF(         ( (work(1,1)-ca)/(1.0D0-ca) )<0.0D0) THEN
   n1 = 0.0D0
ELSE
   n1    = SQRT( (work(1,1)-ca)/(1.0D0-ca) )
ENDIF
IF(         ( (work(2,2)-ca)/(1.0D0-ca) )<0.0D0) THEN
   n2 = 0.0
ELSE
   n2    = SQRT( (work(2,2)-ca)/(1.0D0-ca) )
ENDIF
IF(         ( (work(3,3)-ca)/(1.0D0-ca) )<0.0D0) THEN
   n3 = 0.0
ELSE
   n3    = SQRT( (work(3,3)-ca)/(1.0D0-ca) )
ENDIF
!
loop1: DO i1=-1,1,2
   m1 = n1*REAL(i1)
   loop2: DO i2=-1,1,2
      m2 = n2*REAL(i2)
      loop3: DO i3=-1,1,2
         m3 = n3*REAL(i3)
         IF(ABS(m1*m2*(1.-ca)-m3*sa - work(1,2))<EPS .AND. &
            ABS(m1*m3*(1.-ca)-m2*sa - work(1,3))<EPS .AND. &
            ABS(m2*m3*(1.-ca)-m1*sa - work(2,3))<EPS ) THEN 
            EXIT loop1
         ENDIF
      ENDDO loop3
   ENDDO loop2
ENDDO loop1
axis(1) = m1
axis(2) = m2
axis(3) = m3
!
END SUBROUTINE mat_axis
END MODULE matrix_mod
