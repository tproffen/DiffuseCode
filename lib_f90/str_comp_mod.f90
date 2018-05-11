!*****7***********************************************************      
LOGICAL FUNCTION str_comp (a, b, j, la, lb) 
!-                                                                      
!     compares the first non blank characters of the two strings        
!     for equality. At least j characters must be identical.            
!+                                                                      
IMPLICIT none 
!                                                                       
CHARACTER(LEN=*), INTENT(IN) :: a, b 
INTEGER         , INTENT(IN) :: j, la, lb
!
INTEGER i, ia, ib 
!                                                                       
IF (la == 0 .OR. lb == 0) THEN 
   str_comp = .false. 
ELSE 
   ia = MIN (INDEX (a, ' ') , la) 
   ib = MIN (INDEX (b, ' ') , lb) 
   IF (ia == 0) THEN 
      ia = la 
   ENDIF 
   IF (ib == 0) THEN 
      ib = lb 
   ENDIF 
   i = MIN (ia, ib) 
   i = MIN (i, la) 
   i = MIN (i, lb) 
   IF (i <  j) THEN 
      str_comp = .false. 
   ELSE 
      str_comp = a (1:i) .eq.b (1:i) 
   ENDIF 
ENDIF 
!                                                                       
END FUNCTION str_comp                         
!*****7***********************************************************      
