MODULE math_sup
!
INTERFACE gcd
   MODULE PROCEDURE gcd_two, gcd_three
END INTERFACE gcd
!
CONTAINS
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
gcd_three = gcd_two(divisor, i3)
!
END FUNCTION gcd_three
END MODULE math_sup
