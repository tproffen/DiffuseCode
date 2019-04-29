MODULE gaussj_mod
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE gausj(NP, a, b)
!
! Gauss-Jordan elimination with full pivoting
! NumRec 2.1 p 28-29
! Modified to work with 1-D arrays b
!
USE errlist_mod
!
IMPLICIT NONE
!
INTEGER                  , INTENT(IN)    :: NP  ! Dimensions for matrix a
REAL   , DIMENSION(NP,NP), INTENT(INOUT) :: a   ! 
REAL   , DIMENSION(NP   ), INTENT(INOUT) :: b   ! 
!
!
INTEGER, DIMENSION(NP  ) :: ipiv
INTEGER, DIMENSION(NP  ) :: indxr
INTEGER, DIMENSION(NP  ) :: indxc
INTEGER                  :: j, i, k, l, ll
REAL                     :: big, dum
INTEGER                  :: irow, icol
REAL                     :: pivinv
!
DO j= 1, NP
   ipiv(j) = 0
ENDDO
irow = 0
icol = 0
!
main: do i=1,NP
   big = 0.0
   DO j= 1, NP
      IF(ipiv(j) /= 1) THEN
         DO k= 1, NP
            IF(ipiv(k)==0) THEN
               IF(ABS(a(j,k)) >= big) THEN
                  BIG = ABS(a(j,k))
                  irow = j
                  icol = k
               ENDIF
            ELSEIF(ipiv(k) > 1.) THEN
               ier_num = -1
               ier_typ = ER_MATH
               ier_msg(1) = 'Gauss-Jordan elimination failed'
               ier_msg(2) = 'Are all derivatives zero?'
               ier_msg(3) = 'Check user macro with different values'
               RETURN
            ENDIF
         ENDDO
      ENDIF
   ENDDO
!
   ipiv(icol) = ipiv(icol) + 1
!
   IF(irow /= icol) THEN
      DO l=1,NP
         dum = a(irow,l)
         a(irow,l) = a(icol, l)
         a(icol,l) = dum
      ENDDO
      dum = b(irow  )
      b(irow  ) = b(icol   )
      b(icol  ) = dum
   ENDIF
   indxr(i) = irow
   indxc(i) = icol
   IF(a(icol, icol)== 0.0) THEN
     ier_num = -1
     ier_typ = ER_MATH
     ier_msg(1) = 'Gauss-Jordan elimination failed'
     ier_msg(2) = 'Are all derivatives zero?'
     ier_msg(3) = 'Check user macro with different values'
     RETURN
   ENDIF
   pivinv = 1./a(icol, icol)
   a(icol, icol) = 1.
   DO l=1, NP
      a(icol, l) = a(icol, l)*pivinv
   ENDDO
!
   b(icol   ) = b(icol   )*pivinv
!
   DO ll=1, NP
      IF(ll /= icol) THEN
         dum = a(ll, icol)
         a(ll, icol) = 0.0
         DO l=1, NP
            a(ll,l) = a(ll,l) - a(icol,l)*dum
         ENDDO
!
         b(ll  ) = b(ll  ) - b(icol  )*dum
      ENDIF
   ENDDO
!
ENDDO main
!
DO l=NP,1, -1
   IF(indxr(l) /= indxc(l)) THEN
      DO k=1, NP
         dum = a(k,indxr(l))
         a(k,indxr(l)) = a(k, indxc(l))
         a(k,indxc(l)) = dum
      ENDDO
   ENDIF
ENDDO
!
END SUBROUTINE gausj
!
!*******************************************************************************
!
END MODULE gaussj_mod
