MODULE do_exit_mod
!
CONTAINS
!
SUBROUTINE do_exit 
!                                                                       
USE allocate_appl
USE prompt_mod 
!+                                                                      
!           Clean exit from the program DIFFEV ;-)                      
!-                                                                      
IMPLICIT none 
!                                                                       
!                                                                       
IF (output_io.ne.OUTPUT_SCREEN) then 
   CLOSE (output_io) 
ENDIF 
!
CALL do_deallocate_appl ( 'all',3)
!
END SUBROUTINE do_exit                        
END MODULE do_exit_mod
