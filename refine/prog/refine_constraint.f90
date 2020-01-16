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
REAL    :: r_eta, r_u, r_v,r_w
REAL, PARAMETER :: EPS=1.0E-5
!
refine_fwhm(0) = .FALSE.
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
   ELSEIF(refine_params(i)=='P_eta') THEN
      IF(refine_range(i,1)<refine_range(i,2)) THEN  ! User did set a range
         refine_range(i,1) = MAX(0.0, refine_range(i,1))
         refine_range(i,2) = MIN(1.0, refine_range(i,2))
      ELSE
         refine_range(i,1) = 0.0
         refine_range(i,2) = 1.0
      ENDIF
      refine_fwhm(0) = .TRUE.
   ELSEIF(refine_params(i)=='P_u' .OR. refine_params(i)=='P_v' .OR.             &
          refine_params(i)=='P_w'                                   ) THEN
      refine_fwhm(0) = .TRUE.
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
! Collect P_eta, P-u, P_v and P_w from refined/fixed parameters
DO i=1,refine_fix_n
   IF(refine_fixed (i)=='P_eta') THEN
      r_eta = refine_f(i)
      refine_fwhm(0) = .TRUE.
   ENDIF
   IF(refine_fixed (i)=='P_u'  ) THEN
      r_u = refine_f(i)
      refine_fwhm(0) = .TRUE.
      refine_fwhm(1) = .FALSE.
      refine_fwhm_ind(1) = i
   ENDIF
   IF(refine_fixed (i)=='P_v'  ) THEN
      r_v = refine_f(i)
      refine_fwhm(0) = .TRUE.
      refine_fwhm(2) = .FALSE.
      refine_fwhm_ind(2) = i
   ENDIF
   IF(refine_fixed (i)=='P_w'  ) THEN
      r_w = refine_f(i)
      refine_fwhm(0) = .TRUE.
      refine_fwhm(3) = .FALSE.
      refine_fwhm_ind(3) = i
   ENDIF
ENDDO
DO i=1,refine_par_n
   IF(refine_params(i)=='P_eta') THEN
      r_eta = refine_p(i)
      refine_fwhm(0) = .TRUE.
   ENDIF
   IF(refine_params(i)=='P_u'  ) THEN
      r_u = refine_p(i)
      refine_fwhm(0) = .TRUE.
      refine_fwhm(1) = .TRUE.
      refine_fwhm_ind(1) = i
   ENDIF
   IF(refine_params(i)=='P_v'  ) THEN
      r_v = refine_p(i)
      refine_fwhm(2) = .TRUE.
      refine_fwhm_ind(2) = i
   ENDIF
   IF(refine_params(i)=='P_w'  ) THEN
      r_w = refine_p(i)
      refine_fwhm(0) = .TRUE.
      refine_fwhm(3) = .TRUE.
      refine_fwhm_ind(3) = i
   ENDIF
ENDDO
      
IF(refine_fwhm(0)) THEN        ! Eta and Pu,v,w are used
   IF(.NOT.refine_constrain_test_fwhm(r_u,r_v,r_w, EPS)) THEN
!     WRITE(output_io, 2000) 'P_eta',r_eta
      WRITE(output_io, 2100) 'u,v,w', r_u, r_v, r_w
      ier_num = -10
      ier_typ = ER_APPL
   ENDIF
ENDIF
!2000 FORMAT(' **** Para: ',A16, g15.6E3,4x, '[',g15.6E3,',',g15.6E3,']')
2100 FORMAT(' **** Para: ',A16, g15.6E3,4x, g15.6E3,4x ,g15.6E3)
!
END SUBROUTINE refine_constrain_auto
!
!*******************************************************************************
!
SUBROUTINE refine_constrain_fwhm_in(MAXP, NPARA, a, da, ier_num, ier_typ)
!-
!  Copy the MRQMIN Parameter into the working u,v,w, values and adapt if needed
!+
USE refine_params_mod
IMPLICIT NONE
!
INTEGER                 , INTENT(IN)    :: MAXP        ! MAX number of parameters
INTEGER                 , INTENT(IN)    :: NPARA       ! number of parameters
REAL   , DIMENSION(MAXP), INTENT(IN)    :: a           ! Parameter values
REAL   , DIMENSION(MAXP), INTENT(INOUT) :: da          ! Parameter shifts
INTEGER                 , INTENT(OUT)   :: ier_num
INTEGER                 , INTENT(OUT)   :: ier_typ
!
REAL, DIMENSION(3) :: uvw_in
REAL, DIMENSION(3) :: uvw_del
REAL, PARAMETER :: EPS=1.0E-5
INTEGER :: j
!
IF(refine_fwhm(0)) THEN                                  ! Constrain FWHM
   DO j=1,3
      IF(refine_fwhm(j) .AND. refine_fwhm_ind(j)<=NPARA ) THEN  ! P_u is refined
         uvw_in (j) = a (refine_fwhm_ind(j))                    ! Copy entry from MRQMIN
         uvw_del(j) = da(refine_fwhm_ind(j))                    ! Intended shift P_u 
      ELSE                                                      ! P_u is fixed or absent
         IF(refine_fwhm_ind(j)<=refine_fix_n) THEN              ! P_u is fixed
            uvw_in(j)  = refine_f(refine_fwhm_ind(j))           ! Copy value from fixed list
            uvw_del(j) = 0.0                                    ! P_u may not be modified
         ELSE
            uvw_in(j)  = 0.0                                    ! P_u is absent default to 0.0
            uvw_del(j) = 0.0                                    ! P_u may not be modified
         ENDIF
      ENDIF
   ENDDO
   CALL refine_constrain_fwhm(uvw_in, uvw_del, EPS, ier_num, ier_typ)
   DO j=1,3
      IF(refine_fwhm(j) .AND. refine_fwhm_ind(j)<=NPARA ) THEN  ! P_u is refined
         da(refine_fwhm_ind(j)) = uvw_del(j)                    ! use adjusted value
      ENDIF
   ENDDO
ENDIF
!
END SUBROUTINE refine_constrain_fwhm_in
!
!*******************************************************************************
!
SUBROUTINE refine_constrain_fwhm(uvw_in, uvw_del, EPS, ier_num, ier_typ)
!-
! Constrain DELTA(u,v,w) such that FWHM is > EPS for Q>0.0
! for (uvw) + DELTA(uvw)
!+
IMPLICIT NONE
!
REAL   , DIMENSION(3), INTENT(IN)    :: uvw_in
REAL   , DIMENSION(3), INTENT(INOUT) :: uvw_del
REAL   , INTENT(IN)    :: EPS
INTEGER, INTENT(OUT)   :: ier_num
INTEGER, INTENT(OUT)   :: ier_typ
!
INTEGER :: i  ! dummy index
REAL :: u   ! working u  = u_in + du
REAL :: v   ! working v  = v_in + dv
REAL :: w   ! working w  = w_in + dw
REAL :: fac
!
u = uvw_in(1)
v = uvw_in(2)
w = uvw_in(3)
fac = 1.0
!
ier_num = -10
ier_typ = 6
!
IF(refine_constrain_test_fwhm(u,v,w, EPS)) THEN   ! Initial values are OK
   u = uvw_in(1) + uvw_del(1)
   v = uvw_in(2) + uvw_del(2)
   w = uvw_in(3) + uvw_del(3)
   fac = 0.0
   IF(refine_constrain_test_fwhm(u,v,w, EPS)) THEN ! Modified values are OK
      ier_num = 0
      ier_typ = 0
      fac = 1.0
   ELSE
      loop: DO i=99,0, -1
         fac = REAL(i)*0.01
         u = uvw_in(1) + fac*uvw_del(1)
         v = uvw_in(2) + fac*uvw_del(2)
         w = uvw_in(3) + fac*uvw_del(3)
         IF(refine_constrain_test_fwhm(u,v,w, EPS)) THEN ! Modified values are OK
            ier_num = 0
            ier_typ = 0
            EXIT loop
         ENDIF
      ENDDO loop
   ENDIF
!else
!  write(*,*) ' UVW NEGATIVE !'
ENDIF
IF(ier_num==0 .AND. fac<1.0) THEN
!  write(*,*) ' UVW restricted ', fac
   uvw_del(1) = fac*uvw_del(1)
   uvw_del(2) = fac*uvw_del(2)
   uvw_del(3) = fac*uvw_del(3)
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
!write(*,*) ' zero points ', q1, q2
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
