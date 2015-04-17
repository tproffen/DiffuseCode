MODULE stack_rese_mod
!
CONTAINS
!*****7*****************************************************************
      SUBROUTINE do_stack_rese 
!-                                                                      
!     Resets the stacking fault setup                                   
!+                                                                      
      USE discus_config_mod 
      USE stack_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      st_nlayer = 0 
      st_ntypes = 0 
      st_nchem = 0 
!                                                                       
      END SUBROUTINE do_stack_rese                  
END MODULE stack_rese_mod
