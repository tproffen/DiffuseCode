MODULE str_comp_mod
!
CONTAINS
!
!*****7***********************************************************      
LOGICAL FUNCTION str_comp (a, b, j, la, lb) 
!-                                                                      
!     compares the first non blank characters of the two strings        
!     for equality. At least j characters must be identical.            
!+                                                                      
!
use precision_mod
IMPLICIT none 
!                                                                       
CHARACTER(LEN=*), INTENT(IN) :: a, b 
INTEGER         , INTENT(IN) :: j, la, lb
!
!character(len=PREC_STRING) :: aa
!character(len=PREC_STRING) :: bb
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
   i = max (ia, ib) 
   i = MIN (i, la) 
   i = MIN (i, lb) 
!
!   aa = ' '
!   bb = ' '
!   aa = a(1:len_trim(a))
!   bb = b(1:len_trim(b))
!   i  = max(len_trim(aa), len_trim(bb))
!if(bb(1:3)=='var') then
!write(*,*) ' AA >',aa(1:len_trim(aa)),'<>', i, ' >>',aa(1:i),'<<'
!write(*,*) ' BB >',bb(1:len_trim(bb)),'<>', j, ' >>',bb(1:i),'<<'
!write(*,*) ' ==  ', aa(1:i) .eq.bb(1:i)
!endif
   IF (i <  j) THEN 
      str_comp = .false. 
   ELSE 
      str_comp = a (1:i) .eq.b (1:i) 
!     str_comp = aa(1:i) .eq.bb(1:i) 
   ENDIF 
ENDIF 
!                                                                       
END FUNCTION str_comp                         
!*****7***********************************************************      
END MODULE str_comp_mod
