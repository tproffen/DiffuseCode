MODULE do_exit_mod
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE do_exit 
!-
! Generic exit from discsu_suite
!+
!USE allocate_appl
USE discus_plot_menu, ONLY:jmol_kill
USE exit_mod
USE prompt_mod 
!
use hdf5_state    ! HDF5 was initialized yes/no?
!+                                                                      
!           Clean exit from the program DIFFEV ;-)                      
!-                                                                      
IMPLICIT none 
!
logical :: did_kill
!                                                                       
!  Terminate any Jmol processes started by this discus_suite
!
CALL jmol_kill(.FALSE., .TRUE., did_kill)
!
!CALL do_deallocate_appl ( 'all',3)
CALL exit_all
!
if(hdf5_get_state()) call hdf5_final_close    ! HDF5 was opened, perform final closure
!                                                                       
IF (output_status.ne.OUTPUT_SCREEN) then 
   CLOSE (output_io) 
ENDIF 
!
END SUBROUTINE do_exit                        
!
!*******************************************************************************
!
END MODULE do_exit_mod
