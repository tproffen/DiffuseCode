MODULE exit_mod
!
CONTAINS
!*****7*****************************************************************
      SUBROUTINE exit_all 
!                                                                       
!     This part contains all the program independent commands.          
!                                                                       
      USE errlist_mod 
      USE prompt_mod 
      USE operating_mod
      USE terminal_mod
      IMPLICIT none 
!                                                                       
!
      INTEGER socket_send
!                                                                       
      IF (ier_num.ne.0) THEN 
         CALL errlist 
      ENDIF 
!                                                                       
!------ Close output file                                               
!                                                                       
      IF (output_status.ne.OUTPUT_SCREEN) THEN 
         CLOSE (output_io) 
      ENDIF 
!                                                                       
!------ Close sockets                                                   
!                                                                       
      IF (lremote) THEN 
         ier_num =  socket_send (s_remote, 'bye', 3) 
         CALL socket_close (s_remote) 
      ENDIF 
!                                                                       
      IF (lsocket) THEN 
         CALL socket_close (s_conid) 
         CALL socket_close (s_sock) 
      ENDIF 
!
!       For library procedure style set setup_done to false
!
      lsetup_done = .false.
!
      CALL operating_exit   ! Exit for operating system
      write(*,'(a,a)') COLOR_BG_DEFAULT, COLOR_FG_DEFAULT
!                                                                       
      END SUBROUTINE exit_all                       
!
END MODULE exit_mod
