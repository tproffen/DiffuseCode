MODULE string_convert_mod
!
CONTAINS
!
!*****7*****************************************************************
!
SUBROUTINE do_cap (str) 
!+                                                                      
!       Converts string 'str' to upper letters                          
!-                                                                      
IMPLICIT none 
!                                                                       
CHARACTER (LEN=* ), INTENT(INOUT) :: str 
INTEGER  :: i 
INTEGER  :: len 
!                                                                       
DO i = 1, LEN (str) 
   IF (IACHAR (str (i:i) ) .ge.IACHAR ('a') .and. &
       IACHAR (str (i:i) ) .le.IACHAR ('z')       ) THEN                                            
      str(i:i) = achar(IACHAR(str(i:i)) - IACHAR('a') + IACHAR('A'))
   ENDIF 
ENDDO 
!
END SUBROUTINE do_cap                         
!
!*****7*****************************************************************
!
SUBROUTINE do_low (str) 
!+                                                                      
!       Converts string 'str' to lower case letters                          
!-                                                                      
IMPLICIT none 
!                                                                       
CHARACTER (LEN=* ), INTENT(INOUT) :: str 
INTEGER  :: i 
INTEGER  :: len 
!                                                                       
DO i = 1, LEN (str) 
   IF (IACHAR (str (i:i) ) .ge.IACHAR ('A') .and. &
       IACHAR (str (i:i) ) .le.IACHAR ('Z')       ) THEN                                            
      str(i:i) = achar(IACHAR(str(i:i)) - IACHAR('A') + IACHAR('a'))
   ENDIF 
ENDDO 
!
END SUBROUTINE do_low
!
!*****7***********************************************************      
!
END MODULE string_convert_mod
