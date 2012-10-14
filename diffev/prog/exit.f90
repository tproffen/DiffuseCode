SUBROUTINE do_exit 
!                                                                       
USE allocate_appl
!+                                                                      
!           Clean exit from the program DIFFEV ;-)                      
!-                                                                      
IMPLICIT none 
!                                                                       
      include'prompt.inc' 
!                                                                       
IF (output_io.ne.OUTPUT_SCREEN) then 
   CLOSE (output_io) 
ENDIF 
!
CALL alloc_appl ( 'all',3,'deallocate')
!
END SUBROUTINE do_exit                        
