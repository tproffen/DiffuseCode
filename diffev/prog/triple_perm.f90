!*****7**************************************************************** 
!                                                                       
MODULE triple_perm
!
PRIVATE
PUBLIC  :: do_triple_perm
!
CONTAINS
!
   SUBROUTINE do_triple_perm (j, j1, j2, j3, n) 
!                                                                       
   IMPLICIT none 
!                                                                       
      include'random.inc' 
!                                                                       
   INTEGER, INTENT(IN   ) :: j
   INTEGER, INTENT(INOUT) :: j1
   INTEGER, INTENT(INOUT) :: j2
   INTEGER, INTENT(INOUT) :: j3
   INTEGER, INTENT(IN   ) :: n 
!
   INTEGER              :: ii 
   REAL                 :: r 
!                                                                       
   REAL                 :: ran1 
!                                                                       
   j1 = mod (j + int (ran1 (idum) * (n - 1) ), n) + 1 
   DO while (j1.eq.j) 
      j1 = mod (j1, n) + 1 
   ENDDO 
   j2 = mod (j + int (ran1 (idum) * (n - 1) ), n) + 1 
   DO while (j2.eq.j.or.j2.eq.j1) 
      j2 = mod (j2, n) + 1 
   ENDDO 
   j3 = mod (j + int (ran1 (idum) * (n - 1) ), n) + 1 
   DO while (j3.eq.j.or.j3.eq.j1.or.j3.eq.j2) 
      j3 = mod (j3, n) + 1 
   ENDDO 
!                                                                       
   END SUBROUTINE do_triple_perm                    
!
END MODULE triple_perm
