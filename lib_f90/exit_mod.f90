MODULE exit_mod
!
CONTAINS
!
!*****7*****************************************************************
!
SUBROUTINE exit_all 
!                                                                       
!     This part contains all the program independent commands.          
!                                                                       
USE envir_mod
USE errlist_mod 
USE lib_errlist_func
USE operating_mod
USE precision_mod
USE prompt_mod 
USE terminal_mod
!                                                                       
IMPLICIT none 
!
CHARACTER(LEN=PREC_STRING) :: tempfile   ! jmol script files to remove
CHARACTER(LEN=PREC_STRING) :: line       ! temporary string
INTEGER             :: socket_send
LOGICAL             :: lpresent
!                                                                       
IF (ier_num.ne.0) THEN 
   CALL errlist 
ENDIF 
!
!------ Delete JMOL scripts
!
WRITE(tempfile, '(a,i5.5,a)') '/tmp/jmol.',PID,'.00001.mol'
INQUIRE(FILE=tempfile,EXIST=lpresent)
IF(lpresent) THEN
   line = 'rm -f ' // tempfile(1:len_trim(tempfile)-9) // '*.mol'
   CALL system(line)
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
