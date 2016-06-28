!*****7**************************************************************** 
      SUBROUTINE kuplot_do_exit 
!                                                                       
!       Things to do when KUPLOT exits                                  
!                                                                       
      USE errlist_mod 
      USE prompt_mod 
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
!------ call system wide exit routine                                   
!                                                                       
      CALL exit_all 
!                                                                       
!------ close PGPLOT devices                                            
!                                                                       
      CALL PGEND 
!                                                                       
      END SUBROUTINE kuplot_do_exit                        
!
      SUBROUTINE kuplot_sigint
!
!     Handle KUPLOT specific part of a CTRL-C interrupt.
!     CALLED within standalone KUPLOT only
!     This subroutine calls all KUPLOT specific emergency handlers
!     which can also be called from the SUITE
!
      IMPLICIT NONE
      CHARACTER(LEN=1) :: dummy
      WRITE(*,*)
      WRITE(*,*) ' EMERGENCY Shutdown with USER CTRL-C Interrupt'
!
      CALL exit_all
!
      WRITE(*,*) 
      WRITE(*,*) ' KUPLOT closed by User Request CTRL-C '
      WRITE(*,*) ' For final close down hit ENTER key'
      READ(*,'(a)') dummy
      STOP        ! Terminate program
!
      END SUBROUTINE kuplot_sigint
