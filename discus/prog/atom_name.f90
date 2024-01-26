MODULE atom_name
!
CONTAINS
!
!*****7***************************************************************  
!
CHARACTER(len=9) function at_name (iscat) 
!+                                                                      
!     This function builds the atom name as XX(iscat) to be             
!     able to distinguish between different atom types with             
!     the same name.                                                    
!-                                                                      
USE discus_config_mod 
USE crystal_mod 
USE lib_length
!                                                                       
IMPLICIT none 
!                                                                       
integer, intent(in) :: iscat
CHARACTER(len=5) :: istr 
INTEGER   :: il, is 
!                                                                       
!                                                                       
if(iscat==-1) then
  at_name ='(*)'
else
   IF (iscat.ge.100) then 
      WRITE (istr, 1000) iscat 
   ELSEIF (iscat.ge.10) then 
      WRITE (istr, 1100) iscat 
   ELSE 
      WRITE (istr, 1200) iscat 
   ENDIF 
!                                                                       
   il = len_str (cr_at_lis (iscat) ) 
   is = len_str (istr) 
!                                                                       
   at_name = cr_at_lis (iscat) (1:il) //istr (1:is) 
endif
!                                                                       
 1000 FORMAT     ('(',I3,')') 
 1100 FORMAT     ('(',I2,')') 
 1200 FORMAT     ('(',I1,')') 
END FUNCTION at_name                          
!
!*****7***************************************************************  
!
CHARACTER(len=17) function at_name_long(iscat) 
!+                                                                      
!     This function builds the atom name as XX(iscat) to be             
!     able to distinguish between different atom types with             
!     the same name.                                                    
!-                                                                      
USE discus_config_mod 
USE crystal_mod 
!
USE lib_length
use precision_mod
use blanks_mod
!                                                                       
IMPLICIT none 
!                                                                       
integer, dimension(3), intent(in) :: iscat
CHARACTER(len=PREC_STRING) :: istr 
INTEGER   :: i, il, is 
!                                                                       
!                                                                       
if(iscat(1)==-1) then
  at_name_long ='(*)'
else
   write(istr, '(a1,3(i5,a1))') '(', iscat(1), ',', iscat(2), ',', iscat(3), ')'
   is = 19
   call rem_bl(istr, is)
!                                                                       
   il = len_str(cr_at_lis(iscat(1)) ) 
!                                                                       
   at_name_long = cr_at_lis(iscat(1)) (1:il) //istr(1:is) 
endif
!                                                                       
END FUNCTION at_name_long
!
!*****7***************************************************************  
!
END MODULE atom_name
