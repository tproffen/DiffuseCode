MODULE exit_mod
!
CONTAINS
      SUBROUTINE do_exit 
!+                                                                      
!           Clean exit from the program DISCUS ;-)                      
!-                                                                      
      USE prompt_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      CALL exit_all 
      END SUBROUTINE do_exit                        
END MODULE exit_mod
