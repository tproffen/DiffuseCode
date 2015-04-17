MODULE diffev_do_exit_mod
!
CONTAINS
!
SUBROUTINE diffev_do_exit 
!                                                                       
USE diffev_allocate_appl
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
CALL diffev_do_deallocate_appl ( 'all',3)
!
END SUBROUTINE diffev_do_exit                        
END MODULE diffev_do_exit_mod
