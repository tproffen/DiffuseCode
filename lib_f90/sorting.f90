MODULE sorting_mod

!CONTAINS
!
interface indexx
  module procedure indexx_sp, indexx_dp
end interface indexx
!
contains
!
!*****7**************************************************************** 
!
   SUBROUTINE indexx_sp (n, arr, indx) 
!                                                                       
use precision_mod
!
   IMPLICIT none 
!
!                                                                       
   INTEGER, INTENT(IN   ) :: n
   INTEGER, INTENT(INOUT) :: indx (n)
   REAL(kind=PREC_SP)   , INTENT(IN   ) :: arr (n) 
!
   INTEGER, PARAMETER     :: M      =  7
   INTEGER, PARAMETER     :: NSTACK = 250 
!
   INTEGER i, indxt, ir, itemp, j, jstack, k, l, istack (NSTACK) 
   INTEGER                :: ii
   REAL(kind=PREC_DP):: a 
!
loop11:  DO j = 1, n 
            indx (j) = j 
         END DO loop11
         jstack = 0 
         l = 1 
         ir = n 
mainlp: DO
mainif:  IF (ir - l < m ) THEN 
loop13:     DO j = l + 1, ir 
               indxt = indx (j) 
               a = arr (indxt) 
               ii = 0
loop12:        DO i = j - 1, 1, - 1 
                  IF (arr (indx (i) ) <=  a) THEN
                     II = i
                     EXIT loop12
                  ENDIF
                  indx (i + 1) = indx (i) 
               END DO  loop12
!              i = 0 
!      2       CONTINUE
               i = ii
               indx (i + 1) = indxt 
            END DO  loop13
            IF (jstack == 0) EXIT mainlp   ! We're done !! 
            ir = istack (jstack) 
            l = istack (jstack - 1) 
            jstack = jstack - 2 
         ELSE mainif
            k = (l + ir) / 2 
            itemp = indx (k) 
            indx (k) = indx (l + 1) 
            indx (l + 1) = itemp 
            IF (arr (indx (l + 1) )  >  arr (indx (ir) ) ) THEN 
               itemp = indx (l + 1) 
               indx (l + 1) = indx (ir) 
               indx (ir) = itemp 
            ENDIF 
            IF (arr (indx (l) )  >  arr (indx (ir) ) ) THEN 
               itemp = indx (l) 
               indx (l) = indx (ir) 
               indx (ir) = itemp 
            ENDIF 
            IF (arr (indx (l + 1) )  >  arr (indx (l) ) ) THEN 
               itemp = indx (l + 1) 
               indx (l + 1) = indx (l) 
               indx (l) = itemp 
            ENDIF 
            i = l + 1 
            j = ir 
            indxt = indx (l) 
            a = arr (indxt) 
loop3:      DO
               i = i + 1 
               IF (arr (indx (i) )  <  a) CYCLE loop3 
loop4:         DO
                  j = j - 1 
                  IF (arr (indx (j) )  <= a) EXIT loop4 
               ENDDO loop4
               IF (j <  i) EXIT loop3 
               itemp = indx (i) 
               indx (i) = indx (j) 
               indx (j) = itemp 
            ENDDO loop3
            indx (l) = indx (j) 
            indx (j) = indxt 
            jstack = jstack + 2 
            IF (jstack >  NSTACK) THEN
!              ier_num = -15
!              ier_typ = ER_APPL
!              ier_msg(1) = 'NSTACK too small in indexx' 
               RETURN
            ENDIF
            IF (ir - i + 1 >= j - l) THEN 
               istack (jstack) = ir 
               istack (jstack - 1) = i 
               ir = j - 1 
            ELSE 
               istack (jstack) = j - 1 
               istack (jstack - 1) = l 
               l = i 
            ENDIF 
         ENDIF mainif
      ENDDO mainlp
   END SUBROUTINE indexx_sp
!
!*****7**************************************************************** 
!
   SUBROUTINE indexx_dp (n, arr, indx) 
!                                                                       
use precision_mod
!
   IMPLICIT none 
!
!                                                                       
   INTEGER, INTENT(IN   ) :: n
   INTEGER, INTENT(INOUT) :: indx (n)
   REAL(kind=PREC_DP)   , INTENT(IN   ) :: arr (n) 
!
   INTEGER, PARAMETER     :: M      =  7
   INTEGER, PARAMETER     :: NSTACK = 250 
!
   INTEGER i, indxt, ir, itemp, j, jstack, k, l, istack (NSTACK) 
   INTEGER                :: ii
   REAL(kind=PREC_DP):: a 
!
loop11:  DO j = 1, n 
            indx (j) = j 
         END DO loop11
         jstack = 0 
         l = 1 
         ir = n 
mainlp: DO
mainif:  IF (ir - l < m ) THEN 
loop13:     DO j = l + 1, ir 
               indxt = indx (j) 
               a = arr (indxt) 
               ii = 0
loop12:        DO i = j - 1, 1, - 1 
                  IF (arr (indx (i) ) <=  a) THEN
                     II = i
                     EXIT loop12
                  ENDIF
                  indx (i + 1) = indx (i) 
               END DO  loop12
!              i = 0 
!      2       CONTINUE
               i = ii
               indx (i + 1) = indxt 
            END DO  loop13
            IF (jstack == 0) EXIT mainlp   ! We're done !! 
            ir = istack (jstack) 
            l = istack (jstack - 1) 
            jstack = jstack - 2 
         ELSE mainif
            k = (l + ir) / 2 
            itemp = indx (k) 
            indx (k) = indx (l + 1) 
            indx (l + 1) = itemp 
            IF (arr (indx (l + 1) )  >  arr (indx (ir) ) ) THEN 
               itemp = indx (l + 1) 
               indx (l + 1) = indx (ir) 
               indx (ir) = itemp 
            ENDIF 
            IF (arr (indx (l) )  >  arr (indx (ir) ) ) THEN 
               itemp = indx (l) 
               indx (l) = indx (ir) 
               indx (ir) = itemp 
            ENDIF 
            IF (arr (indx (l + 1) )  >  arr (indx (l) ) ) THEN 
               itemp = indx (l + 1) 
               indx (l + 1) = indx (l) 
               indx (l) = itemp 
            ENDIF 
            i = l + 1 
            j = ir 
            indxt = indx (l) 
            a = arr (indxt) 
loop3:      DO
               i = i + 1 
               IF (arr (indx (i) )  <  a) CYCLE loop3 
loop4:         DO
                  j = j - 1 
                  IF (arr (indx (j) )  <= a) EXIT loop4 
               ENDDO loop4
               IF (j <  i) EXIT loop3 
               itemp = indx (i) 
               indx (i) = indx (j) 
               indx (j) = itemp 
            ENDDO loop3
            indx (l) = indx (j) 
            indx (j) = indxt 
            jstack = jstack + 2 
            IF (jstack >  NSTACK) THEN
!              ier_num = -15
!              ier_typ = ER_APPL
!              ier_msg(1) = 'NSTACK too small in indexx' 
               RETURN
            ENDIF
            IF (ir - i + 1 >= j - l) THEN 
               istack (jstack) = ir 
               istack (jstack - 1) = i 
               ir = j - 1 
            ELSE 
               istack (jstack) = j - 1 
               istack (jstack - 1) = l 
               l = i 
            ENDIF 
         ENDIF mainif
      ENDDO mainlp
   END SUBROUTINE indexx_dp
END MODULE sorting_mod
