PROGRAM discus 
!                                                                       
USE setup_mod
USE discus_loop_mod
!
IMPLICIT none 
!
!
CALL setup
CALL discus_set_sub
CALL discus_loop
!                                                                       
END PROGRAM discus                            
