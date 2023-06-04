MODULE exit_mod
!
private
!
public exit_all
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
use blanks_mod
USE envir_mod
USE errlist_mod 
USE lib_errlist_func
use lib_config
USE operating_mod
USE precision_mod
USE prompt_mod 
use str_comp_mod
use string_convert_mod
USE terminal_mod
!                                                                       
IMPLICIT none 
!
CHARACTER(LEN=PREC_STRING) :: tempfile   ! jmol script files to remove
CHARACTER(LEN=PREC_STRING) :: line       ! temporary string
character(len=8)           :: frm
character(len=PREC_STRING) :: message
integer          :: exit_msg
integer          :: ier_cmd
!INTEGER             :: socket_send
integer          :: ios
integer          :: ll
integer          :: i
LOGICAL             :: lpresent
INTEGER           :: old_version
INTEGER           :: new_version
CHARACTER(LEN=10) :: cversion
INTEGER           :: since_update

!                                                                       
IF (ier_num.ne.0) THEN 
   CALL errlist 
ENDIF 
!CALL lib_f90_init_updates
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
   if(generic_get_interval()>0) then
      IF(since_update>generic_get_interval()) THEN
         line = 'NO'
         write(*,*)
         write(*,'(a,i4,a)' ) ' Last Ubuntu update >',generic_get_interval(), ' days ago, '
         write(*,'(a)' ) ' Operating system update is strongly recommended'
         write(*,'(a)', advance='no' ) ' Update Ubuntu now ? Yes/No/Never/New Interval in days : '
         read(*, '(a)', iostat=ios) line
         call do_cap(line)
         ll = len_trim(line)
         call rem_leading_bl(line, ll)
         if(is_iostat_end(ios) .or. line(1:3)=='YES' .or. &
            str_comp(line, 'Y', 1, ll, 1)                ) then
            CALL lib_f90_update_ubuntu
         elseif(line(1:2)=='NO' .or. str_comp(line, 'NO', 2, ll, 2)) then
            write(*,*)
            write(*,'(A)') ' Ubuntu update deferred for right now'
            write(*,*)
         elseif(line(1:2)=='NEVER' .or. str_comp(line, 'NEVER', 2, ll, 5)) then
            write(*,*)
            write(*,'(A)') ' DISCUS will not update Ubuntu any longer' 
            write(*,'(A)') ' Regular updates are strongly recommended' 
            write(*,*)
            call generic_set_interval(-1)
            call generic_write_config(discus_dir, discus_dir_l)
         else
            write(frm,'(a,i1,a)') '(i',len_trim(line),')'
            read(line,fmt=frm,iostat=ios) i
            if(ios==0) then
               if(i>0) then
                  call generic_set_interval(i)
                  call generic_write_config(discus_dir, discus_dir_l)
                  write(*,*)
                  write(*,'(A)') ' Ubuntu update deferred for right now'
                  write(*,'(a, i4, a)') ' DISCUS will ask again in ',i, ' days'
                  write(*,*)
               else
                  write(*,*)
                  write(*,'(A)') ' DISCUS will not update Ubuntu any longer' 
                  write(*,'(A)') ' Regular updates are strongly recommended' 
                  write(*,*)
                  call generic_set_interval(-1)
                  call generic_write_config(discus_dir, discus_dir_l)
               endif
            endif
         endif
      ENDIF
   endif
ENDIF
!
!------ Delete JMOL scripts
!
WRITE(tempfile, '(a,i10.10,a,i10.10,a)') '/tmp/jmol.',PID,'.',1,'.mol'
INQUIRE(FILE=tempfile,EXIST=lpresent)
IF(lpresent) THEN
   line = 'rm -f ' // tempfile(1:len_trim(tempfile)-9) // '*.mol'
   call execute_command_line(line(1:LEN_TRIM(line)), wait=.false., &
        CMDSTAT=ier_cmd, CMDMSG=message, EXITSTAT=exit_msg)
ENDIF
!
WRITE(tempfile, '(a,i10.10,a,i10.10,a)') '/tmp/jmol.',PID,'.',1,'.pid'
INQUIRE(FILE=tempfile,EXIST=lpresent)
!write(*,*) ' JMOL ', tempfile(1:len_trim(tempfile)), ' ', lpresent
IF(lpresent) THEN
   line = 'rm -f ' // tempfile(1:len_trim(tempfile)-9) // '*.pid'
   call execute_command_line(line(1:LEN_TRIM(line)), wait=.false., &
        CMDSTAT=ier_cmd, CMDMSG=message, EXITSTAT=exit_msg)
ENDIF
!
WRITE(tempfile, '(a,i10.10,a)') '/tmp/jmol.',PID,'.pid'
INQUIRE(FILE=tempfile,EXIST=lpresent)
!write(*,*) ' JMOL ', tempfile(1:len_trim(tempfile)), ' ', lpresent
IF(lpresent) THEN
   line = 'rm -f ' // tempfile(1:len_trim(tempfile))
   call execute_command_line(line(1:LEN_TRIM(line)), wait=.false., &
        CMDSTAT=ier_cmd, CMDMSG=message, EXITSTAT=exit_msg)
ENDIF
!
!-- Delete temporary DISCUS*PID files
!
WRITE(tempfile, '(a,i10.10)') '/tmp/DISCUS.UFILE.',PID
INQUIRE(FILE=tempfile,EXIST=lpresent)
!write(*,*) ' UFIL ', tempfile(1:len_trim(tempfile)), ' ', lpresent
IF(lpresent) THEN
   line = 'rm -f ' // tempfile(1:len_trim(tempfile))
   call execute_command_line(line(1:LEN_TRIM(line)), wait=.false., &
        CMDSTAT=ier_cmd, CMDMSG=message, EXITSTAT=exit_msg)
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
!write(*,'(a,a)') COLOR_BG_DEFAULT, COLOR_FG_DEFAULT
!
1000 FORMAT(' New DISCUS version available at GIThub : ',a)
1100 FORMAT(' To update start DISCUS again and use     '  )
1200 FORMAT(' the command:         ''update''          '  )
1300 FORMAT(' or run the installation script           '  )
!                                                                       
END SUBROUTINE exit_all                       
!
END MODULE exit_mod
