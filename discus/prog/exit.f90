MODULE discus_exit_mod
!
CONTAINS
      SUBROUTINE discus_do_exit 
!+                                                                      
!           Clean exit from the program DISCUS ;-)                      
!-                                                                      
      USE prompt_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      CALL exit_all 
      END SUBROUTINE discus_do_exit                        
END MODULE discus_exit_mod
