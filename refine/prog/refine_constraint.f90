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
use refine_data_mod
!
USE errlist_mod
USE prompt_mod
use precision_mod
!
IMPLICIT NONE
!
INTEGER :: i    ! Dummy loop variable
REAL(kind=PREC_DP)    :: r_eta, r_eta_l, r_eta_q
REAL(kind=PREC_DP)    :: r_u, r_v, r_w
REAL(kind=PREC_DP), PARAMETER :: EPS=1.0E-5
!
refine_fwhm(0) = .FALSE.
refine_eta (0) = .FALSE.
refine_eta_ind  = 0
refine_fwhm_ind = 0
r_eta   = 0.0D0
r_eta_l = 0.0D0
r_eta_q = 0.0D0
r_u     = 0.0D0
r_v     = 0.0D0
r_w     = 0.0D0
DO i=1,refine_par_n
   IF(refine_params(i)(1:5)=='P_lat'  .OR.                                      &
      refine_params(i)(1:6)=='P_biso' .OR.                                      &
      refine_params(i)(1:5)=='P_dia'                                            &
                                           ) THEN
      IF(refine_range(i,1)<refine_range(i,2)) THEN  ! User did set a range
         refine_range(i,1) = MAX(0.0D0, refine_range(i,1))
      ELSE
         refine_range(i,1) = 0.0
         refine_range(i,2) = HUGE(0.0)
      ENDIF
   ELSEIF(refine_params(i)=='P_eta') THEN
      IF(refine_range(i,1)<refine_range(i,2)) THEN  ! User did set a range
         refine_range(i,1) = MAX(0.0D0, refine_range(i,1))
         refine_range(i,2) = MIN(1.0D0, refine_range(i,2))
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
      refine_eta(0) = .TRUE.
      refine_eta(3) = .FALSE.
      refine_eta_ind(3) = i
   ENDIF
   IF(refine_fixed (i)=='P_eta_l') THEN
      r_eta_l = refine_f(i)
      refine_eta(0) = .TRUE.
      refine_eta(2) = .FALSE.
      refine_eta_ind(2) = i
   ENDIF
   IF(refine_fixed (i)=='P_eta_q') THEN
      r_eta_q = refine_f(i)
      refine_eta(0) = .TRUE.
      refine_eta(1) = .FALSE.
      refine_eta_ind(1) = i
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
      r_eta   = refine_p(i)
      refine_eta(0) = .TRUE.
      refine_eta(3) = .TRUE.
      refine_eta_ind(3) = i
   ENDIF
   IF(refine_params(i)=='P_eta_l') THEN
      r_eta   = refine_p(i)
      refine_eta(0) = .TRUE.
      refine_eta(2) = .TRUE.
      refine_eta_ind(2) = i
   ENDIF
   IF(refine_params(i)=='P_eta_q') THEN
      r_eta   = refine_p(i)
      refine_eta(0) = .TRUE.
      refine_eta(1) = .TRUE.
      refine_eta_ind(1) = i
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
   IF(.NOT.refine_constrain_test_fwhm(r_u,r_v,r_w, EPS, ref_x(1), ref_x(ref_dim(1)))) THEN
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
SUBROUTINE refine_constrain_fwhm_in(MAXP, NPARA, a, da, qmax, qmin, ier_num, ier_typ)
!-
!  Copy the MRQMIN Parameter into the working u,v,w, values and adapt if needed
!+
USE refine_params_mod
use precision_mod
IMPLICIT NONE
!
INTEGER                            , INTENT(IN)    :: MAXP        ! MAX number of parameters
INTEGER                            , INTENT(IN)    :: NPARA       ! number of parameters
REAL(kind=PREC_DP), DIMENSION(MAXP), INTENT(IN)    :: a           ! Parameter values
REAL(kind=PREC_DP), DIMENSION(MAXP), INTENT(INOUT) :: da          ! Parameter shifts
real(kind=PREC_DP)                 , intent(in)    :: qmin        ! Qmin in experimental data
real(kind=PREC_DP)                 , intent(in)    :: qmax        ! Qmax in experimental data
INTEGER                            , INTENT(OUT)   :: ier_num
INTEGER                            , INTENT(OUT)   :: ier_typ
!
REAL(kind=PREC_DP), DIMENSION(3) :: uvw_in
REAL(kind=PREC_DP), DIMENSION(3) :: uvw_del
REAL(kind=PREC_DP), PARAMETER :: EPS=1.0E-5
INTEGER :: j
!
IF(refine_fwhm(0)) THEN                                  ! Constrain FWHM
   DO j=1,3
      IF(refine_fwhm(j) .AND. refine_fwhm_ind(j)<=NPARA ) THEN  ! P_u is refined
         uvw_in (j) = a (refine_fwhm_ind(j))                    ! Copy entry from MRQMIN
         uvw_del(j) = da(refine_fwhm_ind(j))                    ! Intended shift P_u 
      ELSE                                                      ! P_u is fixed or absent
         IF(refine_fwhm_ind(j)> 0             .and.   &         ! P_u is fixed
            refine_fwhm_ind(j)<=refine_fix_n) THEN              ! P_u is fixed
            uvw_in(j)  = refine_f(refine_fwhm_ind(j))           ! Copy value from fixed list
            uvw_del(j) = 0.0                                    ! P_u may not be modified
         ELSE
            uvw_in(j)  = 0.0                                    ! P_u is absent default to 0.0
            uvw_del(j) = 0.0                                    ! P_u may not be modified
         ENDIF
      ENDIF
   ENDDO
   CALL refine_constrain_fwhm(uvw_in, uvw_del, EPS, qmin, qmax, ier_num, ier_typ)
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
SUBROUTINE refine_constrain_fwhm(uvw_in, uvw_del, EPS, qmin, qmax, ier_num, ier_typ)
!-
! Constrain DELTA(u,v,w) such that FWHM is > EPS for Q>0.0
! for (uvw) + DELTA(uvw)
!+
use precision_mod
IMPLICIT NONE
!
REAL(kind=PREC_DP), DIMENSION(3), INTENT(IN)    :: uvw_in
REAL(kind=PREC_DP), DIMENSION(3), INTENT(INOUT) :: uvw_del
REAL(kind=PREC_DP)              , iNTENT(IN)    :: EPS
real(kind=PREC_DP)              , intent(in)    :: qmin        ! Qmin in experimental data
real(kind=PREC_DP)              , intent(in)    :: qmax        ! Qmax in experimental data
INTEGER                         , INTENT(OUT)   :: ier_num
INTEGER                         , INTENT(OUT)   :: ier_typ
!
INTEGER :: i  ! dummy index
REAL(kind=PREC_DP) :: u   ! working u  = u_in + du
REAL(kind=PREC_DP) :: v   ! working v  = v_in + dv
REAL(kind=PREC_DP) :: w   ! working w  = w_in + dw
REAL(kind=PREC_DP) :: fac
!
u = uvw_in(1)
v = uvw_in(2)
w = uvw_in(3)
fac = 1.0
!
ier_num = -10
ier_typ = 6
!
IF(refine_constrain_test_fwhm(u,v,w, EPS, qmin, qmax)) THEN   ! Initial values are OK
   u = uvw_in(1) + uvw_del(1)
   v = uvw_in(2) + uvw_del(2)
   w = uvw_in(3) + uvw_del(3)
   fac = 0.0
   IF(refine_constrain_test_fwhm(u,v,w, EPS, qmin, qmax)) THEN ! Modified values are OK
      ier_num = 0
      ier_typ = 0
      fac = 1.0
   ELSE
      loop: DO i=99,0, -1
         fac = REAL(i, kind=PREC_DP)*0.01d0
         u = uvw_in(1) + fac*uvw_del(1)
         v = uvw_in(2) + fac*uvw_del(2)
         w = uvw_in(3) + fac*uvw_del(3)
         IF(refine_constrain_test_fwhm(u,v,w, EPS, qmin, qmax)) THEN ! Modified values are OK
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
LOGICAL FUNCTION refine_constrain_test_fwhm(u,v,w_in, EPS, qmin, qmax)
!
! Tests if FHWM is > EPS for Q > 0.0
! u*Q**2 + v*Q + w - EPS > 0.0
!+
use precision_mod
IMPLICIT NONE
!
REAL(kind=PREC_DP), INTENT(IN)    :: u
REAL(kind=PREC_DP), INTENT(IN)    :: v
REAL(kind=PREC_DP), INTENT(IN)    :: w_in
REAL(kind=PREC_DP), INTENT(IN)    :: EPS
real(kind=PREC_DP), intent(in)    :: qmin        ! Qmin in experimental data
real(kind=PREC_DP), intent(in)    :: qmax        ! Qmax in experimental data
!
REAL(kind=PREC_DP) :: w   ! working w  = w_in - EPS
REAL(kind=PREC_DP) :: q1  ! Lower zero point
REAL(kind=PREC_DP) :: q2  ! Upper zero point
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
         IF(q2 < qmin) THEN                ! Zero point q < Qmin     , tolerate
            refine_constrain_test_fwhm = .TRUE.
         ELSEIF(q1 > qmax  ) THEN          ! Zero point q > Qmax     , tolerate
            refine_constrain_test_fwhm = .TRUE.
         ELSE
            refine_constrain_test_fwhm = .FALSE.
         ENDIF
      ENDIF
   ELSE                                    ! U==0 linear test
      IF(v/=0.0) THEN                      ! Non-zero v
         q1 = -w/v
         IF(v>0.0 .AND. q1 < qmin) THEN    ! positive slope, Zero point q < Qmin , tolerate
            refine_constrain_test_fwhm = .TRUE.
         ELSEIF(v<0.0 .AND. q1 > qmax  ) THEN  ! negative slope Zero point q > Qmax , tolerate
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
!*******************************************************************************
!
SUBROUTINE refine_constrain_eta_in(MAXP, NPARA, a, da, qmin, qmax, ier_num, ier_typ)
!-
!  Copy the MRQMIN Parameter into the working eta  , values and adapt if needed
!+
USE refine_params_mod
use precision_mod
IMPLICIT NONE
!
INTEGER                            , INTENT(IN)    :: MAXP        ! MAX number of parameters
INTEGER                            , INTENT(IN)    :: NPARA       ! number of parameters
REAL(kind=PREC_DP), DIMENSION(MAXP), INTENT(IN)    :: a           ! Parameter values
REAL(kind=PREC_DP), DIMENSION(MAXP), INTENT(INOUT) :: da          ! Parameter shifts
real(kind=PREC_DP)                 , intent(in)    :: qmin        ! Qmin in experimental data
real(kind=PREC_DP)                 , intent(in)    :: qmax        ! Qmax in experimental data
INTEGER                            , INTENT(OUT)   :: ier_num
INTEGER                            , INTENT(OUT)   :: ier_typ
!
REAL(kind=PREC_DP), DIMENSION(3) :: eta_in
REAL(kind=PREC_DP), DIMENSION(3) :: eta_del
REAL(kind=PREC_DP), PARAMETER :: EPS=1.0E-5
INTEGER :: j
!
IF(refine_eta(0)) THEN                                          ! Constrain FWHM
   DO j=1,3
      IF(refine_eta(j) .AND. refine_eta_ind(j)<=NPARA ) THEN    ! P_eta is refined
         eta_in (j) = a (refine_eta_ind(j))                     ! Copy entry from MRQMIN
         eta_del(j) = da(refine_eta_ind(j))                     ! Intended shift P_eta 
      ELSE                                                      ! P_eta is fixed or absent
         IF(refine_eta_ind(j)> 0            .and. &             ! P_eta is fixed
            refine_eta_ind(j)<=refine_fix_n) THEN               ! P_eta is fixed
            eta_in(j)  = refine_f(refine_eta_ind(j))            ! Copy value from fixed list
            eta_del(j) = 0.0                                    ! P_eta may not be modified
         ELSE
            eta_in(j)  = 0.0                                    ! P_eta is absent default to 0.0
            eta_del(j) = 0.0                                    ! P_eta may not be modified
         ENDIF
      ENDIF
   ENDDO
   CALL refine_constrain_eta(eta_in, eta_del, EPS, qmin, qmax, ier_num, ier_typ)
   DO j=1,3
      IF(refine_eta(j) .AND. refine_eta_ind(j)<=NPARA ) THEN  ! P_u is refined
         da(refine_eta_ind(j)) = eta_del(j)                    ! use adjusted value
      ENDIF
   ENDDO
ENDIF
!
END SUBROUTINE refine_constrain_eta_in
!
!*******************************************************************************
!
SUBROUTINE refine_constrain_eta(eta_in, eta_del, EPS, qmin, qmax, ier_num, ier_typ)
!-
! Constrain DELTA(u,v,w) such that FWHM is > EPS for Q>0.0
! for (eta) + DELTA(eta)
!+
use precision_mod
IMPLICIT NONE
!
REAL(kind=PREC_DP), DIMENSION(3), INTENT(IN)    :: eta_in
REAL(kind=PREC_DP), DIMENSION(3), INTENT(INOUT) :: eta_del
REAL(kind=PREC_DP)              , INTENT(IN)    :: EPS
real(kind=PREC_DP)              , intent(in)    :: qmin        ! Qmin in experimental data
real(kind=PREC_DP)              , intent(in)    :: qmax        ! Qmax in experimental data
INTEGER                         , INTENT(OUT)   :: ier_num
INTEGER                         , INTENT(OUT)   :: ier_typ
!
INTEGER :: i  ! dummy index
REAL(kind=PREC_DP) :: u   ! working u  = u_in + du
REAL(kind=PREC_DP) :: v   ! working v  = v_in + dv
REAL(kind=PREC_DP) :: w   ! working w  = w_in + dw
REAL(kind=PREC_DP) :: fac
!
u = eta_in(1)
v = eta_in(2)
w = eta_in(3)
fac = 1.0
!
ier_num = -10
ier_typ = 6
!
IF(refine_constrain_test_eta(u,v,w, EPS, qmin, qmax)) THEN   ! Initial values are OK
   u = eta_in(1) + eta_del(1)
   v = eta_in(2) + eta_del(2)
   w = eta_in(3) + eta_del(3)
   fac = 0.0
   IF(refine_constrain_test_eta(u,v,w, EPS, qmin, qmax)) THEN ! Modified values are OK
      ier_num = 0
      ier_typ = 0
      fac = 1.0
   ELSE
      loop: DO i=99,0, -1
         fac = REAL(i, kind=PREC_DP)*0.01d0
         u = eta_in(1) + fac*eta_del(1)
         v = eta_in(2) + fac*eta_del(2)
         w = eta_in(3) + fac*eta_del(3)
         IF(refine_constrain_test_eta(u,v,w, EPS, qmin, qmax)) THEN ! Modified values are OK
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
   eta_del(1) = fac*eta_del(1)
   eta_del(2) = fac*eta_del(2)
   eta_del(3) = fac*eta_del(3)
ENDIF
!
END SUBROUTINE refine_constrain_eta
!
!*******************************************************************************
!
LOGICAL FUNCTION refine_constrain_test_eta(u,v,w_in, EPS, qmin, qmax)
!
! Tests if ETA is > EPS for Q > 0.0 and < 1.0
! u*Q**2 + v*Q + w - EPS > 0.0
!+
use precision_mod
IMPLICIT NONE
!
REAL(kind=PREC_DP), INTENT(IN)    :: u
REAL(kind=PREC_DP), INTENT(IN)    :: v
REAL(kind=PREC_DP), INTENT(IN)    :: w_in
REAL(kind=PREC_DP), INTENT(IN)    :: EPS
real(kind=PREC_DP), intent(in)    :: qmin        ! Qmin in experimental data
real(kind=PREC_DP), intent(in)    :: qmax        ! Qmax in experimental data
!
REAL(kind=PREC_DP) :: w   ! working w  = w_in - EPS
REAL(kind=PREC_DP) :: q01  ! Lower zero point
REAL(kind=PREC_DP) :: q02  ! Upper zero point
REAL(kind=PREC_DP) :: q11  ! Lower one  point
REAL(kind=PREC_DP) :: q12  ! Upper one  point
REAL(kind=PREC_DP) :: qext ! location of extremum
REAL(kind=PREC_DP) :: emax ! extremum value
!
w = w_in - EPS
refine_constrain_test_eta = .FALSE.
q01 = -1.000E4
q02 = -1.000E4
!
IF(u==0.0 .AND. v==0.0) then               ! Constant eta
   if(w>=0.0 .and. w<=1) THEN              ! All positive, OK
      refine_constrain_test_eta = .TRUE.
   else                                    ! Constant value outside [0,1]
      refine_constrain_test_eta = .false.
   endif
ELSE                                       ! At least one is negative or zero
   IF(u/=0.0) THEN                         ! Non-zero u, square test
      qext = -(v/u)*005D0                  ! Extremum
      emax = u*qext**2 + v*qext + w        ! Extremum
      if(u<0.0D0) then                     ! Negative parabola
         if(emax < 0.0D0) then             ! Negative everywhere!
            refine_constrain_test_eta = .FALSE.
         elseif(emax < 1.0D0) then         ! Less than 1.0 everywhere!
            q01 = -0.5*v/u - SQRT(0.25*(v/u)**2 - w/u)          ! Left zero intercept
            q02 = -0.5*v/u + SQRT(0.25*(v/u)**2 - w/u)          ! Right zero intercept
            IF(q01 < qmin .and. q02> qmax) THEN    ! Zero points outside Q, tolerate
               refine_constrain_test_eta = .TRUE.
            ELSE
               refine_constrain_test_eta = .FALSE.
            ENDIF
         else                              ! General solution
            q01 = -0.5*v/u - SQRT(0.25*(v/u)**2 -  w/u)         ! Left zero intercept
            q02 = -0.5*v/u + SQRT(0.25*(v/u)**2 -  w/u)         ! Right zero intercept
            q11 = -0.5*v/u - SQRT(0.25*(v/u)**2 - (w-1.0D0)/u)  ! Left ONE intercept
            q12 = -0.5*v/u + SQRT(0.25*(v/u)**2 - (w-1.0D0)/u)  ! Right ONE intercept
            if(q01<qmin .and. q11>qmax ) then         ! Left intercepts include Qrange; OK
               refine_constrain_test_eta = .true.
            elseif(q12<qmin .and. q02>qmax ) then     ! Right intercepts include Qrange; OK
               refine_constrain_test_eta = .true.
            else
               refine_constrain_test_eta = .FALSE.
            endif
         endif
      else                                 ! Positive parabola
         if(emax > 1.0D0) then             ! Larger than 1.0 everywhere!
            refine_constrain_test_eta = .FALSE.
         elseif(emax > 0.0D0) then         ! Larger than 0.0 everywhere!
            q11 = -0.5*v/u - SQRT(0.25*(v/u)**2 - (w-1.0d0)/u)   ! Left ONE intercept
            q12 = -0.5*v/u + SQRT(0.25*(v/u)**2 - (w-1.0d0)/u)   ! Right ONE intercep
            IF(q11 < qmin .and. q12>qmax   ) THEN   ! Zero point outsiede Q range; tolerate
               refine_constrain_test_eta = .TRUE.
            ELSE
               refine_constrain_test_eta = .FALSE.
            ENDIF
         else                              ! General solution
            q01 = -0.5*v/u - SQRT(0.25*(v/u)**2 -  w/u)         ! Left zero intercept
            q02 = -0.5*v/u + SQRT(0.25*(v/u)**2 -  w/u)         ! Right zero intercept
            q11 = -0.5*v/u - SQRT(0.25*(v/u)**2 - (w-1.0D0)/u)  ! Left ONE intercept
            q12 = -0.5*v/u + SQRT(0.25*(v/u)**2 - (w-1.0D0)/u)  ! Right ONE intercept
            if(q11<qmin .and. q01>qmax) then         ! Left intercepts include Qrange; OK
               refine_constrain_test_eta = .true.
            elseif(q02<qmin .and. q12>qmax) then     ! Right intercepts include Qrange; OK
               refine_constrain_test_eta = .true.
            else
               refine_constrain_test_eta = .FALSE.
            endif
         endif
      endif
   ELSE                                    ! U==0 linear test
      IF(v/=0.0) THEN                      ! Non-zero v
         q01 = (0.0D0-w)/v                 ! Q at which eta = 0
         q11 = (1.0D0-w)/v                 ! Q at which eta = 1
         IF(v>0.0) then                    ! positive slope 
            if(q01<qmin  .and. q11> qmax   ) then  ! Intercepts outside range, OK
               refine_constrain_test_eta = .TRUE.
            else
               refine_constrain_test_eta = .FALSE.
            endif
         elseif(v<0.0) then                ! negative slope 
            if(q11<qmin  .and. q01> qmax   ) then  ! Intercepts outside range, OK
               refine_constrain_test_eta = .TRUE.
            else
               refine_constrain_test_eta = .FALSE.
            endif
         ELSE
            refine_constrain_test_eta = .FALSE.
         ENDIF
      ELSE
         IF(w>=0.0 .and. w <= 1.0) THEN
            refine_constrain_test_eta = .TRUE.
         ELSE
            refine_constrain_test_eta = .FALSE.
         ENDIF
      ENDIF
   ENDIF
ENDIF
!write(*,*) ' zero points ', q1, q2
!
END FUNCTION refine_constrain_test_eta
!
!*******************************************************************************
!
subroutine refine_constrain_temp(line, length, luvw)
!
USE get_params_mod
USE ber_params_mod
use precision_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: line
INTEGER         , INTENT(INOUT) :: length
logical         , INTENT(IN   ) :: luvw
!
INTEGER, PARAMETER :: MAXW = 5
CHARACTER (LEN=MAX(PREC_STRING,LEN(line))), DIMENSION(MAXW) :: cpara   != ' '
INTEGER             , DIMENSION(MAXW) :: lpara = 0
REAL(KIND=PREC_DP)  , DIMENSION(MAXW) :: werte = 0.0
INTEGER :: ianz
logical :: ans = .false.
!
CALL get_params(line, ianz, cpara, lpara, MAXW, length)
CALL ber_params (ianz, cpara, lpara, werte, MAXW)
!
if(luvw) then
ans = refine_constrain_test_fwhm(werte(1),werte(2), werte(3), 1.0D-5, werte(4), werte(5))
else
ans = refine_constrain_test_eta(werte(1),werte(2), werte(3), 1.0D-5, werte(4), werte(5))
endif

!write(*,'('' TEST == '', L1)') ans
!
end subroutine refine_constrain_temp
!
!*******************************************************************************
!
END MODULE refine_constraint_mod
