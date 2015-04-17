PROGRAM discus 
!                                                                       
USE discus_setup_mod
USE discus_loop_mod
!
IMPLICIT none 
!
LOGICAL, PARAMETER :: standalone = .true.
!
CALL discus_setup(standalone)
CALL discus_set_sub
CALL discus_loop
!                                                                       
END PROGRAM discus                            
