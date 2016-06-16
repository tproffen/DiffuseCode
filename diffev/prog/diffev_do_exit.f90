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
!
SUBROUTINE diffev_emergency_save
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
END SUBROUTINE diffev_emergency_save
!
SUBROUTINE diffev_emergency_mpi
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
END SUBROUTINE diffev_emergency_mpi
END MODULE diffev_do_exit_mod
