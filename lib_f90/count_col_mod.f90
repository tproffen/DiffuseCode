MODULE count_col_mod
!
CONTAINS
!
!*****7***************************************************************  
SUBROUTINE count_col (zeile, ianz) 
!+                                                                      
!     This subroutine counts the number of columns of string 'zeile'    
!-                                                                      
USE prompt_mod 
!
IMPLICIT NONE
!                                                                       
CHARACTER (LEN=* ), INTENT(IN)    ::  zeile 
INTEGER           , INTENT(INOUT) :: ianz
!
INTEGER :: i 
LOGICAL :: ein 
!                                                                       
INTEGER :: len_str 
!                                                                       
ianz = 0 
ein = .false. 
!                                                                       
DO i = 1, len_str (zeile) 
   IF (zeile (i:i) /=  ' ') THEN 
      ein = .true. 
   ELSEIF (zeile (i:i) ==  ' '.AND.ein) THEN 
      ein = .false. 
      ianz = ianz + 1 
   ENDIF 
ENDDO 
!                                                                       
IF (ein) ianz = ianz + 1 
!
END SUBROUTINE count_col                      
!
!*****7**************************************************************** 
!
END MODULE count_col_mod
