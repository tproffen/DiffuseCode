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
USE appl_env_mod
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
!INTEGER             :: socket_send
LOGICAL             :: lpresent
INTEGER           :: old_version
INTEGER           :: new_version
CHARACTER(LEN=10) :: cversion
INTEGER           :: since_update

!                                                                       
IF (ier_num.ne.0) THEN 
   CALL errlist 
ENDIF 
CALL lib_f90_init_updates
CALL lib_f90_test_updates(old_version, new_version, cversion, since_update)
IF(new_version > old_version ) THEN
   WRITE(output_io,*)
   WRITE(output_io,1000) cversion
   WRITE(output_io,1100) 
   WRITE(output_io,1200) 
   WRITE(output_io,1300) 
   WRITE(output_io,*)
ENDIF
IF(operating == OS_LINUX_WSL) THEN
   IF(since_update>7) THEN
      CALL lib_f90_update_ubuntu
   ENDIF
ENDIF
!
!------ Delete JMOL scripts
!
WRITE(tempfile, '(a,i10.10,a,i10.10,a)') '/tmp/jmol.',PID,'.',1,'.mol'
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
!IF (lremote) THEN 
!   ier_num =  socket_send (s_remote, 'bye', 3) 
!   CALL socket_close (s_remote) 
!ENDIF 
!                                                                       
!IF (lsocket) THEN 
!   CALL socket_close (s_conid) 
!   CALL socket_close (s_sock) 
!ENDIF 
!
!       For library procedure style set setup_done to false
!
lsetup_done = .false.
!
CALL operating_exit   ! Exit for operating system
write(*,'(a,a)') COLOR_BG_DEFAULT, COLOR_FG_DEFAULT
!
1000 FORMAT(' New DISCUS version available at GIThub : ',a)
1100 FORMAT(' To update start DISCUS again and use     '  )
1200 FORMAT(' the command:         ''update''          '  )
1300 FORMAT(' or run the installation script           '  )
!                                                                       
END SUBROUTINE exit_all                       
!
END MODULE exit_mod
