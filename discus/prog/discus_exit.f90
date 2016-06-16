MODULE discus_exit_mod
!
CONTAINS
      SUBROUTINE discus_do_exit 
!+                                                                      
!           Clean exit from the program DISCUS ;-)                      
!-                                                                      
      USE prompt_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      CALL exit_all 
      END SUBROUTINE discus_do_exit                        
!
      SUBROUTINE discus_emergency_save
!
!     Write the structure to a file EMERGENCY.STRU
!
      USE discus_config_mod
      USE crystal_mod
      USE save_menu
!
      CHARACTER (LEN=14) :: strucfile
!
      WRITE(*,*) ' SAVING STRUCTURE TO EMERGENCY.STRU '
!
      CALL save_default_setting        ! Set default save flags
      strucfile = 'EMERGENCY.STRU'
      CALL save_keyword(strucfile)
!
      END SUBROUTINE discus_emergency_save
END MODULE discus_exit_mod
