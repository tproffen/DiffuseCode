MODULE do_exit_mod
!
CONTAINS
!
SUBROUTINE do_exit 
!                                                                       
!USE allocate_appl
USE prompt_mod 
!+                                                                      
!           Clean exit from the program DIFFEV ;-)                      
!-                                                                      
IMPLICIT none 
!                                                                       
!
!CALL do_deallocate_appl ( 'all',3)
CALL exit_all
!                                                                       
IF (output_status.ne.OUTPUT_SCREEN) then 
   CLOSE (output_io) 
ENDIF 
!
END SUBROUTINE do_exit                        
END MODULE do_exit_mod
