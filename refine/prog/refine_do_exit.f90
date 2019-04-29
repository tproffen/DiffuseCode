MODULE refine_do_exit_mod
!
CONTAINS
!
SUBROUTINE refine_do_exit 
!                                                                       
USE refine_allocate_appl
!
USE exit_mod
USE prompt_mod 
!+                                                                      
!           Clean exit from the program REFINE ;-)                      
!-                                                                      
IMPLICIT none 
INTEGER :: length
!                                                                       
!                                                                       
length = 3
CALL refine_do_deallocate_appl ( 'all',length)
CALL exit_all
!
IF (output_io.ne.OUTPUT_SCREEN) THEN 
   CLOSE (output_io) 
ENDIF 
!
END SUBROUTINE refine_do_exit                        
!
SUBROUTINE refine_emergency_save
!
!     Write the GENERATION , PARAMETER SUMMARY files
!     Currently left blank intentionally, as these
!     Files are updated regularly at points where it
!     makes sense to update them
!
IMPLICIT NONE
!
!WRITE(*,*) ' SAVING STRUCTURE TO EMERGENCY.STRU '
!
!
END SUBROUTINE refine_emergency_save
!
SUBROUTINE refine_emergency_mpi
!
!     Closes down MPI in case of emergency. 
!     Currently just a finalize, to be developed further
!
USE diffev_mpi_mod
IMPLICIT NONE
!
CALL run_mpi_finalize
!
!
END SUBROUTINE refine_emergency_mpi
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
END MODULE refine_do_exit_mod
