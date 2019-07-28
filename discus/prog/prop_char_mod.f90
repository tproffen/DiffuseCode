MODULE prop_char_mod
!
CONTAINS
!
!*****7*****************************************************************
!
SUBROUTINE char_prop_1 (c_property, property, length) 
!-                                                                      
!     sets letters for true property bits                               
!     property = 1; ==> Letter, i.e. property is present                
!     property = 0; ==> "-"     i.e. property is absent                 
!+                                                                      
USE prop_para_mod 
!
IMPLICIT none 
!                                                                       
CHARACTER(LEN=*), INTENT(OUT) :: c_property 
INTEGER         , INTENT(IN)  :: property 
INTEGER         , INTENT(OUT) :: length 
!                                                                       
INTEGER :: i 
!                                                                       
DO i = 0, MAXPROP - 1
   IF (BTEST (property, i) ) THEN 
      c_property (i + 1:i + 1) = c_prop_letter (i + 1:i + 1) 
   ELSEIF (.NOT.BTEST (property, i) ) THEN 
      c_property (i + 1:i + 1) = '-' 
   ENDIF 
ENDDO 
length = MAXPROP 
!
END SUBROUTINE char_prop_1                    
!
!*****7*****************************************************************
!
SUBROUTINE char_prop_2(c_property, property, mask, length) 
!-                                                                      
!     sets letters for true property bits                               
!     property ; mask = 1; 1 ==> Letter, i.e. property is present       
!     property ; mask = 0; 1 ==> "-" i.e. property must be absent       
!     property ; mask = 0; 0 ==> "." i.e. property is ignored           
!     property ; mask = 1; 0 ==> "." i.e. property is ignored           
!      This last situation should never occur...                        
!+                                                                      
USE prop_para_mod 
!                                                                       
IMPLICIT none 
!                                                                       
CHARACTER(LEN=*), INTENT(OUT) :: c_property 
INTEGER         , INTENT(IN)  :: property 
INTEGER         , INTENT(IN)  :: mask 
INTEGER         , INTENT(OUT) :: length 
!                                                                       
INTEGER :: i 
!                                                                       
DO i = 0, MAXPROP - 1
   IF (BTEST (property, i) .AND.BTEST (mask, i) ) THEN 
      c_property (i + 1:i + 1) = c_prop_letter (i + 1:i + 1) 
   ELSEIF (.NOT.BTEST (property, i) .AND.BTEST (mask, i) ) THEN 
      c_property (i + 1:i + 1) = '-' 
   ELSEIF (.NOT.BTEST (property, i) .AND..NOT.BTEST (mask, i) ) THEN 
      c_property (i + 1:i + 1) = '.' 
   ELSEIF (BTEST (property, i) .AND..NOT.BTEST (mask, i) ) THEN 
      c_property (i + 1:i + 1) = '.' 
   ENDIF 
ENDDO 
length = MAXPROP 
!
END SUBROUTINE char_prop_2                    
!
!*****7*****************************************************************
!
END MODULE prop_char_mod
