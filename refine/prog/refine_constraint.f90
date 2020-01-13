MODULE refine_constraint_mod
!
!  Rotines related to parameter constraints
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE refine_constrain_auto
!-
! Apply and set automatic constraints for special parameter names
! These are:
! P_lat...
! P_biso...
! P_dia...
!+
USE refine_params_mod
!
USE errlist_mod
USE prompt_mod
!
IMPLICIT NONE
!
INTEGER :: i    ! Dummy loop variable
!
DO i=1,refine_par_n
   IF(refine_params(i)(1:5)=='P_lat'  .OR.                                      &
      refine_params(i)(1:6)=='P_biso' .OR.                                      &
      refine_params(i)(1:5)=='P_dia'                                            &
                                           ) THEN
      IF(refine_range(i,1)<refine_range(i,2)) THEN  ! User did set a range
         refine_range(i,1) = MAX(0.0, refine_range(i,1))
      ELSE
         refine_range(i,1) = 0.0
         refine_range(i,2) = HUGE(0.0)
      ENDIF
   ENDIF
   IF(refine_range(i,1)<refine_range(i,2)) THEN  ! USER/AUTO constrain
      IF(refine_p(i)<refine_range(i,1) .OR.    &
         refine_range(i,2)<refine_p(i)) THEN
         WRITE(output_io, 1000) refine_params(i), refine_p(i), &
                                refine_range(i,1), refine_range(i,2)
         ier_num = -11
         ier_typ = ER_APPL
      ENDIF
   ENDIF
1000 FORMAT(' **** Para: ',A16, g15.6E3,4x, '[',g15.6E3,',',g15.6E3,']')
ENDDO
!
END SUBROUTINE refine_constrain_auto
!
!*******************************************************************************
!
SUBROUTINE refine_constrain_fwhm(u_in,v_in,w_in, du,dv,dw, EPS, ier_num)
!-
! Constrain DELTA(u,v,w) such that FWHM is > EPS for Q>0.0
! for (uvw) + DELTA(uvw)
!+
IMPLICIT NONE
!
REAL   , INTENT(IN)    :: u_in
REAL   , INTENT(IN)    :: v_in
REAL   , INTENT(IN)    :: w_in
REAL   , INTENT(INOUT) :: du
REAL   , INTENT(INOUT) :: dv
REAL   , INTENT(INOUT) :: dw
REAL   , INTENT(IN)    :: EPS
INTEGER, INTENT(OUT)   :: ier_num
!
INTEGER :: i  ! dummy index
REAL :: u   ! working u  = u_in + du
REAL :: v   ! working v  = v_in + dv
REAL :: w   ! working w  = w_in + dw
!
u = u_in
v = v_in
w = w_in
!
ier_num = -10
!
IF(refine_constrain_test_fwhm(u,v,w, EPS)) THEN   ! Initial values are OK
   u = u_in + du
   v = v_in + dv
   w = w_in + dw
   IF(refine_constrain_test_fwhm(u,v,w, EPS)) THEN ! Modified values are OK
      ier_num = 0
   ELSE
      loop: DO i=99,0, -1
         u = u_in + REAL(i)*0.01*du
         v = v_in + REAL(i)*0.01*dv
         w = w_in + REAL(i)*0.01*dw
         IF(refine_constrain_test_fwhm(u,v,w, EPS)) THEN ! Modified values are OK
            ier_num = 0
            EXIT loop
         ENDIF
      ENDDO loop
   ENDIF
ENDIF
!
END SUBROUTINE refine_constrain_fwhm
!
!*******************************************************************************
!
LOGICAL FUNCTION refine_constrain_test_fwhm(u,v,w_in, EPS)
!
! Tests if FHWM is > EPS for Q > 0.0
! u*Q**2 + v*Q + w - EPS > 0.0
!+
IMPLICIT NONE
!
REAL, INTENT(IN)    :: u
REAL, INTENT(IN)    :: v
REAL, INTENT(IN)    :: w_in
REAL, INTENT(IN)    :: EPS
!
REAL :: w   ! working w  = w_in - EPS
REAL :: q1  ! Lower zero point
REAL :: q2  ! Upper zero point
!
w = w_in - EPS
refine_constrain_test_fwhm = .FALSE.
q1 = -1.000E4
q2 = -1.000E4
!
IF(u>0.0 .AND. v>0.0 .AND. w>0.0) THEN     ! All positive, OK
   refine_constrain_test_fwhm = .TRUE.
ELSE                                       ! At least one is negative or zero
   IF(u/=0.0) THEN                         ! Non-zero u, square test
      IF(0.25*(v/u)**2 - w/u<0.0) THEN     ! No zero points
         refine_constrain_test_fwhm = .TRUE.
      ELSE
         q1 = -0.5*v/u - SQRT(0.25*(v/u)**2 - w/u)
         q2 = -0.5*v/u + SQRT(0.25*(v/u)**2 - w/u)
         IF(q2 < 0.0 ) THEN                ! Zero point at negative Q, tolerate
            refine_constrain_test_fwhm = .TRUE.
         ELSEIF(q1 > 200.0 ) THEN          ! Zero point at huge     Q, tolerate
            refine_constrain_test_fwhm = .TRUE.
         ELSE
            refine_constrain_test_fwhm = .FALSE.
         ENDIF
      ENDIF
   ELSE                                    ! U==0 linear test
      IF(v/=0.0) THEN                      ! Non-zero v
         q1 = -w/v
         IF(v>0.0 .AND. q1 < 0.0 ) THEN    ! positive slope Zero point at negative Q, tolerate
            refine_constrain_test_fwhm = .TRUE.
         ELSEIF(v<0.0 .AND. q1 > 200.0 ) THEN  ! negative slope Zero point at huge     Q, tolerate
            refine_constrain_test_fwhm = .TRUE.
         ELSE
            refine_constrain_test_fwhm = .FALSE.
         ENDIF
      ELSE
         IF(w>0.0) THEN
            refine_constrain_test_fwhm = .TRUE.
         ELSE
            refine_constrain_test_fwhm = .FALSE.
         ENDIF
      ENDIF
   ENDIF
ENDIF
write(*,*) ' zero points ', q1, q2
!
END FUNCTION refine_constrain_test_fwhm
!
subroutine refine_constrain_temp(line, length)
!
USE get_params_mod
USE ber_params_mod
use precision_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: line
INTEGER         , INTENT(INOUT) :: length
!
INTEGER, PARAMETER :: MAXW = 3
CHARACTER (LEN=1024), DIMENSION(MAXW) :: cpara   = ' '
INTEGER             , DIMENSION(MAXW) :: lpara = 0
REAL(KIND=PREC_DP)  , DIMENSION(MAXW) :: werte = 0.0
INTEGER :: ianz
logical :: ans = .false.
!
CALL get_params(line, ianz, cpara, lpara, MAXW, length)
CALL ber_params (ianz, cpara, lpara, werte, MAXW)
!
ans = refine_constrain_test_fwhm(real(werte(1)),real(werte(2)), real(werte(3)), 1.0E-5)
write(*,'('' TEST == '', L1)') ans
!
end subroutine refine_constrain_temp
!
!*******************************************************************************
!
END MODULE refine_constraint_mod
