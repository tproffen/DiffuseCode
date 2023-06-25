MODULE math_sup
!-
! some useful general math stuff
!
!+
!
! gcd          ! Greatest common divisor
! cardano      ! Solve a cubic equation: 0 = C*x^3 + Q*x^2 + L*x + const ; only real roots
! eigen_value  ! Determine eigenvalues and Eigenvectors of a symmetric 3x3 matrix
!
INTERFACE gcd
   MODULE PROCEDURE gcd_two, gcd_three
END INTERFACE gcd
!
interface math_solve_poly
   module procedure math_cardano, math_quadratic, math_linear
end interface math_solve_poly
!
private
public gcd
public math_solve_poly
public eigen_value
public test_eigen_value
!
CONTAINS
!
!*******************************************************************************
!
INTEGER FUNCTION gcd_two(i1, i2)
!
! Determine the greatest common divisor
!
IMPLICIT NONE
INTEGER, INTENT(IN)  :: i1,i2
!
INTEGER :: j1,j2, k
!
k = 1
IF(IABS(i1) > IABS(i2)) THEN
   j1 = IABS(i1)
   j2 = IABS(i2)
   IF(j2==0) THEN
      gcd_two = j1
      RETURN
   ENDIF
ELSEIF(IABS(i1) < IABS(i2)) THEN
   j1 = IABS(i2)
   j2 = IABS(i1)
   IF(j2==0) THEN
      gcd_two = j1
      RETURN
   ENDIF
ELSE
   j1 = IABS(i1)
   j2 = IABS(i2)
   IF(j2==0) THEN
      gcd_two = 1
      RETURN
   ENDIF
ENDIF
!
search: DO
   k = MOD(j1,j2)
   IF(k == 0) EXIT search
   j1 = j2
   j2 = k
ENDDO search
!
gcd_two = j2
!
END FUNCTION gcd_two
!
!*******************************************************************************
!
INTEGER FUNCTION gcd_three(i1,i2,i3)
!
! Determine the greatest common divisor. Version for three integers
!
IMPLICIT NONE
INTEGER, INTENT(IN) :: i1
INTEGER, INTENT(IN) :: i2
INTEGER, INTENT(IN) :: i3
INTEGER             :: divisor
!
divisor   = gcd_two(i1, i2)
IF(divisor>1) THEN
   gcd_three = gcd_two(divisor, i3)
ELSE
   gcd_three = 1
   IF(i1==0 .AND. I2==0 .AND. IABS(I3)>0) THEN
      gcd_three = IABS(i3)
   ENDIF
ENDIF
!
END FUNCTION gcd_three
!
!*******************************************************************************
!
subroutine math_cardano(a, b, c, d, n, root)
!
! Determine the real roots of a cubic equation
! A x^3 + B x^2 + C x + D = 0
! with A, B, C, D real and 
! A /= 0
!
!+
!
use precision_mod
use wink_mod
!
implicit none
!
real(kind=PREC_DP)              , intent(in)  :: a     ! A x^3
real(kind=PREC_DP)              , intent(in)  :: b     ! B x^2
real(kind=PREC_DP)              , intent(in)  :: c     ! C x
real(kind=PREC_DP)              , intent(in)  :: d     ! D
integer                         , intent(out) :: n     ! number of distinct roots
real(kind=PREC_DP), dimension(3), intent(out) :: root
!
real(kind=prec_DP), parameter :: eps = 1.0D-5
real(kind=prec_DP) :: p
real(kind=prec_DP) :: q
real(kind=prec_DP) :: delta
real(kind=prec_DP) :: u
real(kind=prec_DP) :: v
real(kind=prec_DP) :: t1
real(kind=prec_DP) :: t2
real(kind=prec_DP) :: arg
!
!write(*,*) ' a,b,c,d', a,b,c,d
if(abs(a)<eps) then         ! Not a cubic equations, but a quadratic
   call math_quadratic(b, c, d, n, root)
   return
endif
p     = ( 9.0D0*a*c -3.0D0*b**2)/(9.0D0*a**2)
q     = (-9.0D0*a*b*c + 27.0D0*a**2*d + 2.0D0*b**3)/(27.0D0*a**3)
delta = (18.0D0*a*b*c*d - 4.0D0*a*c**3 - 27.0D0*a**2*d**2 + b**2*c**2 - 4.0D0*b**3*d)/(108.D0*a**4)
!write(*,*) ' p, q, d', p,q, delta
!
if(abs(delta)< EPS .and. abs(p)<EPS) then            ! Case A Triple solution
!write(*,*) ' CASE A '
   root(1) = -b/3.0D0
   root(2) = -b/3.0D0
   root(3) = -b/3.0D0
   n       = 1
elseif(abs(delta)< EPS .and. abs(p)>=EPS) then       ! Case B a single and a double solution
!write(*,*) ' CASE B '
   root(1) = (b**3 - 4.0D0*a*b*c + 9.0D0*a**2*d) / (3.0D0*a**2*c - a*b**2)
   root(2) = (b*c  - 9.0D0*a*d) / (6.0D0*a*c - 2.0D0*b**2)
   root(3) = (b*c  - 9.0D0*a*d) / (6.0D0*a*c - 2.0D0*b**2)
   n       = 2
elseif(delta<0.0D0) then                             ! Case C a single solution
!write(*,*) ' CASE C '
   t1   = -q/2. + sqrt(-delta)
   t2   = -q/2. - sqrt(-delta)
   u    =  ((abs(t1))**(1./3.))
   v    =  ((abs(t2))**(1./3.))
   if(t1<0) then
     u     = -u
   endif
   if(t2<0) then
     v     = -v
   endif
   root(1) = -b/a/3. + u + v
   root(2) = root(1)
   root(3) = root(1)
   n       = 1
elseif(delta>0.0D0) then                            ! Case D Three solutions
   arg =               (-0.5D0*q*sqrt(-27.D0/p**3))
!write(*,*) ' ARG 1 ', arg
   if(arg<-1.0D0) arg = -1.0D0
   if(arg> 1.0D0) arg =  1.0D0
   root(1) = -b/a/3.0D0 + sqrt(-4.D0/3.D0*p)*cos(1.D0/3.D0*acos(arg)          )
   root(2) = -b/a/3.0D0 - sqrt(-4.D0/3.D0*p)*cos(1.D0/3.D0*acos(arg)+ pi/3.0D0)
   root(3) = -b/a/3.0D0 - sqrt(-4.D0/3.D0*p)*cos(1.D0/3.D0*acos(arg)- pi/3.0D0)
!write(*,*) ' CASE D ', root
   n       = 3
endif
!
end subroutine math_cardano
!
!*******************************************************************************
!
subroutine math_quadratic(a, b, c, n, root)
!
! Determine the real roots of a quadratic equation
! A x^2 + B x   + C       = 0
! with A, B, C    real and 
! A /= 0
!
!+
!
use precision_mod
!
implicit none
!
real(kind=PREC_DP)              , intent(in)  :: a     ! A x^2
real(kind=PREC_DP)              , intent(in)  :: b     ! B x
real(kind=PREC_DP)              , intent(in)  :: c     ! C  
integer                         , intent(out) :: n     ! number of distinct roots
real(kind=PREC_DP), dimension(3), intent(out) :: root
!
real(kind=prec_DP), parameter :: eps = 1.0D-8
real(kind=prec_DP) :: p
real(kind=prec_DP) :: q
real(kind=prec_DP) :: arg
!
if(abs(a)<eps) then         ! Not a quadratic equation, but a linear
   call math_linear(b, c, n, root)
   return
endif
p = b/a
q = c/a
arg = 0.25*p*p - q
!
if(arg>=0.0D0) then      ! Two solutions
   root(1) = -0.5D0*p + sqrt(arg)
   root(2) = -0.5D0*p - sqrt(arg)
   root(3) = -0.5D0*p + sqrt(arg)
   n       = 2
elseif(arg < 0.0D0) then ! No real solutions
   n       = 0
else                     ! A single solution
   root(1) = -0.5D0*p + sqrt(arg)
   root(2) = -0.5D0*p + sqrt(arg)
   root(3) = -0.5D0*p + sqrt(arg)
   n       = 1
endif
!
end subroutine math_quadratic
!
!*******************************************************************************
!
subroutine math_linear(a, b, n, root)
!-
! Determine the real roots of a linear equation
! A x   + B               = 0
! with A, B       real and 
! A /= 0
!
!+
!
use precision_mod
!
implicit none
!
real(kind=PREC_DP)              , intent(in)  :: a     ! A x
real(kind=PREC_DP)              , intent(in)  :: b     ! B 
integer                         , intent(out) :: n     ! number of distinct roots
real(kind=PREC_DP), dimension(3), intent(out) :: root
!
real(kind=prec_DP), parameter :: eps = 1.0D-8
!
if(abs(a) < EPS) then    ! Constant, no solution
   n = 0
else
   root(1) = -b/a
   root(2) = -b/a
   root(3) = -b/a
   n       = 1
endif
!
end subroutine math_linear
!
!*******************************************************************************
!
!-
!  Calculate the Eigenvalues of a symmetric 3x3 matrix within a cartesian space
!  Eigenvalues based on
!  Based on O.K. Smith "Eigenvalue of a Symmetric 3 X 3 Matrix", Comm of the ACM, (1961), 4, 168
!  3m = Tr(a)
!  2q = det(a - mI)            (I = Unit matrix)
!  6p = SUM( (a-mI)i_ij^2)     ! SUm of squared elements of (a -mI)
!  Eigenvalues are 
!  l_1 = m + 2 sqrt(p)   cos(phi)
!  l_2 = m -   sqrt(p) ( cos(phi) + sqrt(3) sind(phi))
!  l_3 = m -   sqrt(p) ( cos(phi) - sqrt(3) sind(phi))
!  with
!  phi = 1/3 atan( sqrt(p^3 -q^2) / q)  ; phi : [0, pi]
!
!  Eigenvectors based on math.stackexchange.com/questions/3760388/general-formula-for-eigenvectors-of-a-3x3-matrix
! "General formula for Eigenvectors of a 3x3 matrix"
!       ( a b c )
!   A = ( d e f )
!       ( g h i )
!  test_1 = ( a-lambda, b       , c       )
!  test_2 = ( d       , e-lambda, f       )
!  test_3 = ( g       , h       , i-lambda)
!  Vector product of any two test_i gives eigenvectors
!  Additions RBN: length of any test_i must be /= null
!                 Special cases for 2 equal eigenvalues
!                 Special cases for 3 equal eigenvalues
!+
!
!*******************************************************************************
!
subroutine eigen_value(a, eigen_val, eigen_vec, gten, neigen)
!
use errlist_mod
use lib_metric_mod
use matrix_mod
use precision_mod
use prompt_mod
use wink_mod
!
implicit none
!
real(kind=PREC_DP), dimension(3,3), intent(in)  :: a          ! The matrix
real(kind=PREC_DP), dimension(3)  , intent(out) :: eigen_val  ! Eigenvalues
real(kind=PREC_DP), dimension(3,3), intent(out) :: eigen_vec  ! Eigenvectors
real(kind=PREC_DP), dimension(3,3), intent(in)  :: gten
integer                           , intent(out) :: neigen     ! Number distinct Eigenvalues
!integer                           , intent(out) :: ier_num    ! Error message 0 == success
!
real(kind=PREC_DP), parameter :: TOL = 1.0D-10
real(kind=PREC_DP), parameter :: EQL = 1.0D-5     ! Akzept values as equal if difference is smaller
integer :: i
integer :: j
integer :: k
real(kind=PREC_DP), dimension(3,3)  :: b       ! a -  m*I
real(kind=PREC_DP), dimension(3,3)  :: t       ! Trialvectors
real(kind=PREC_DP), dimension(3,3)  :: one_mat ! Trialvectors
real(kind=PREC_DP), dimension(3  )  :: u       ! Vectors, whose vecotor product will be eigenvector
real(kind=PREC_DP), dimension(3  )  :: v       ! Vectors, whose vecotor product will be eigenvector
real(kind=PREC_DP), dimension(3  )  :: t_len   ! Length of trial vectors
real(kind=PREC_DP), dimension(3  )  :: eigen   ! Eigenvalues prior to sorting
real(kind=PREC_DP), dimension(3  )  :: ang     ! Angles of eigenvector 1 to base vectors

real(kind=PREC_DP) :: phi                      ! Angle in 
real(kind=PREC_DP) :: p                        ! 6p = SUM(b(i,j)^2)
real(kind=PREC_DP) :: q                        ! 2q = DET(a - mI)
real(kind=PREC_DP) :: m                        !  m = Tr(a)
real(kind=PREC_DP) :: det                      ! a determinant
real(kind=PREC_DP) :: length2                  ! vector length ^2
real(kind=PREC_DP) :: length                   ! a vector length
real(kind=PREC_DP) :: cphi                     ! cos(phi) length
real(kind=PREC_DP) :: sphi                     ! sin(phi) length
!
ier_num = 0
eigen_val = 0.0D0
eigen_vec = 0.0D0
neigen    = 0
one_mat = 0.0D0
one_mat(1,1) = 1.0D0
one_mat(2,2) = 1.0D0
one_mat(3,3) = 1.0D0
!
det = det3(a)
if(abs(det)<TOL) then        ! Determinant(A) is zero
   ier_num = -4 
   ier_typ = ER_FORT
   ier_msg(1) = 'Determinant of A is zero'
   return
endif
!
m = (a(1,1) + a(2,2) + a(3,3))/3.0D0            ! m = Trace(A) / 3
!
b = a - m*one_mat                               ! A - mI
!
q = det3(b)/2.0D0                               ! q = determinant(1) / 2
p = ( b(1,1)**2 + b(1,2)**2 + b(1,3)**2 +    &  ! p = SUM( (A - mI)_ij^2 )
      b(2,1)**2 + b(2,2)**2 + b(2,3)**2 +    &
      b(3,1)**2 + b(3,2)**2 + b(3,3)**2 )/6.0D0
!
if((p**3 -q**2)<0.0D0) then                     ! Negative root
!  write(output_io,*) ' p^3 - q^2 is negative', p**3 -q**2
   ier_num = -5
   ier_typ = ER_FORT
   ier_msg(1) = 'Eigenvalue is complex '
   return
endif
!
if(abs(q)<TOL       ) then                      ! q == Null phi = PI/6
   phi = pi/2.0D0 / 3.0D0
   cphi = sqrt(3.0D0)*0.5D0
   sphi = 0.50D0
else
   phi = datan(sqrt(p**3 -q**2)/q)/3.
   if(phi<0) phi=phi + pi
   cphi = cos(phi)
   sphi = sin(phi)
endif
!write(*,*)  ' MPQ     ', m, p, q, abs(q)<TOL
!write(*,*)  ' phi     ', phi, p**3 -q**2
!write(*,*)  ' c s phi ', cphi, sphi, cos(phi), sin(phi)
!
eigen(1) = m + 2.0D0*sqrt(p) * dcos(phi)
eigen(2) = m - sqrt(p)*(cphi + sqrt(3.0D0)*sphi)
eigen(3) = m - sqrt(p)*(cphi - sqrt(3.0D0)*sphi)
eigen_val = eigen
!
!  Determine if : all equal, two equal or all different
!
if(abs(eigen(1)-eigen(2))<EQL .and. abs(eigen(1)-eigen(3))<EQL .and.   &
   abs(eigen(2)-eigen(3))<EQL) then           ! All eigenvalues are equal
   eigen_val = eigen
   neigen = 1
elseif(abs(eigen(1)-eigen(2))<EQL .or. abs(eigen(1)-eigen(3))<EQL .or. &
       abs(eigen(2)-eigen(3))<EQL ) then      ! one unique Eigenvalue
   neigen = 2
   if(abs(eigen(1)-eigen(2))<EQL) then        ! sort unique to first
      eigen_val(1) = eigen(3)
      eigen_val(2) = eigen(1)
      eigen_val(3) = eigen(2)
   elseif(abs(eigen(1)-eigen(3))<EQL) then
      eigen_val(1) = eigen(2)
      eigen_val(2) = eigen(1)
      eigen_val(3) = eigen(3)
   elseif(abs(eigen(2)-eigen(3))<EQL) then
      eigen_val(1) = eigen(1)
      eigen_val(2) = eigen(2)
      eigen_val(3) = eigen(3)
   endif
else                                          ! All eigenvalues differ
   neigen = 3
   do i=1, 3                                  ! Sort, largest first
      j = maxloc(eigen, 1)
      eigen_val(i) = eigen(j)
      eigen(j) = -huge(1.0D0)
   enddo
endif
!
!  Calculate eigenvectors
!  Based on stackexchange answer for three distinct eigenvalues
!  Modified to take into account that test vector length can be zero
!  If two eigenvalues are equal, the first eigenvector is OK, as
!  I take the vector product of the two non-zero test vectors
!
do i=1, 3
   t(1,1) = a(1,1) - eigen_val(i)
   t(2,1) = a(1,2)
   t(3,1) = a(1,3) 
!
   t(1,2) = a(2,1)
   t(2,2) = a(2,2) - eigen_val(i)
   t(3,2) = a(2,3) 
!
   t(1,3) = a(3,1)
   t(2,3) = a(3,2) 
   t(3,3) = a(3,3) - eigen_val(i)
   t_len(1) = sqrt(t(1,1)**2 + t(2,1)**2 + t(3,1)**2)
   t_len(2) = sqrt(t(1,2)**2 + t(2,2)**2 + t(3,2)**2)
   t_len(3) = sqrt(t(1,3)**2 + t(2,3)**2 + t(3,3)**2)
!  Determine the two longest vectors
   j = maxloc(t_len,1)
   v = t(:,j)
   t_len(j) = -1.0d0
   k = maxloc(t_len,1)
   u = t(:,k)
!  Eigenvector is vector product of test vectors
   eigen_vec(1, i) = u(2)*v(3) - u(3)*v(2)
   eigen_vec(2, i) = u(3)*v(1) - u(1)*v(3)
   eigen_vec(3, i) = u(1)*v(2) - u(2)*v(1)
   length = sqrt(eigen_vec(1,i)**2 + eigen_vec(2,i)**2 + eigen_vec(3,i)**2)
   if(length>0.0D0) eigen_vec(:,i) = eigen_vec(:,i)/length  ! Normalize to one
enddo
!
!  Special cases for equal eigenvalues
!
if(neigen==1) then                ! All eigenvalues are equal, use base vectors
   eigen_vec(1,1) = 1.00D0
   eigen_vec(2,2) = 1.00D0
   eigen_vec(3,3) = 1.00D0
elseif(neigen==2) then            ! Two eigenvalues are equal
   do i = 1, 3                    ! Determine deviation of (angle (base, eigenvector)) from 90°
      ang(i) = abs(lib_bang(gten, eigen_vec(:,i), one_mat(:,i)) - 90.0D0)
   enddo
   j = minloc(ang,1)              ! This base vector is closest to 90° off 1st eigenvector
   u = eigen_vec(:,1)             ! Vector product of 1st Eigenvector and this base vector
   v = one_mat(:,j)               ! Will create 2nd Eigenvector
   eigen_vec(1, 2) = u(2)*v(3) - u(3)*v(2)
   eigen_vec(2, 2) = u(3)*v(1) - u(1)*v(3)
   eigen_vec(3, 2) = u(1)*v(2) - u(2)*v(1)
!
   v = eigen_vec(:,2)             ! Third eigenvector is vprod(1st x 2nd )
   eigen_vec(1, 3) = u(2)*v(3) - u(3)*v(2)
   eigen_vec(2, 3) = u(3)*v(1) - u(1)*v(3)
   eigen_vec(3, 3) = u(1)*v(2) - u(2)*v(1)
endif
!
if(det3(eigen_vec)<0.0D0) eigen_vec(:,3) = -eigen_vec(:,3)  ! Ensure righthandedness
!
end subroutine eigen_value
!
!*******************************************************************************
!
subroutine test_eigen_value(a, eigen_val, eigen_vec, gten)
!
use errlist_mod
use lib_metric_mod
use matrix_mod
use precision_mod
use prompt_mod
use wink_mod
!
implicit none
!
real(kind=PREC_DP), dimension(3,3), intent(in)  :: a          ! The matrix
real(kind=PREC_DP), dimension(3)  , intent(out) :: eigen_val  ! Eigenvalues
real(kind=PREC_DP), dimension(3,3), intent(out) :: eigen_vec  ! Eigenvectors
real(kind=PREC_DP), dimension(3,3), intent(in)  :: gten
!integer                           , intent(out) :: ier_num    ! Error message 0 == success
!
integer :: i
real(kind=PREC_DP), dimension(3,3)  :: t       ! Trialvectors
!
do i= 1, 3
   t(:,i) = matmul(a, eigen_vec(:,i))
enddo
do i=1, 3
  write(output_io,'(a,i2,3f12.6)') ' Lambda i    ', i, eigen_val(i), &
                              lib_blen(gten, t(:,i)) / lib_blen(gten, eigen_vec(:,i)) , &
                              lib_bang(gten, t(:,i),                  eigen_vec(:,i))
enddo
write(output_io,*) ' Det(eigenvectors ) ', det3(eigen_vec)
write(output_io,'(a,2x, f12.6)') ' angle 1, 2 ', lib_bang(gten,eigen_vec(:,1), eigen_vec(:,2))
write(output_io,'(a,2x, f12.6)') ' angle 1, 3 ', lib_bang(gten,eigen_vec(:,1), eigen_vec(:,3))
write(output_io,'(a,2x, f12.6)') ' angle 2, 3 ', lib_bang(gten,eigen_vec(:,2), eigen_vec(:,3))
!
end subroutine test_eigen_value
!
!*******************************************************************************
!
END MODULE math_sup
