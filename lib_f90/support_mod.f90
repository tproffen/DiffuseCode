MODULE support_mod
!
! Contains basic support routines 
!
interface seknds
  module procedure seknds_sp, seknds_dp
end interface seknds
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE remove_comment (line, ll) 
!                                                                       
!     removes trailing in line comments                                 
!                                                                       
USE charact_mod
!                                                                       
IMPLICIT none 
!                                                                       
CHARACTER (LEN=*), INTENT(INOUT) :: line 
INTEGER          , INTENT(INOUT) ::ll 
!                                                                       
INTEGER :: i 
LOGICAL :: quote 
LOGICAL :: search 
!                                                                       
search = .true. 
DO i = 1, ll 
   quote = line (i:i) .eq.'"'.or.line (i:i) .eq.'''' 
   IF (quote) THEN 
      search = .not.search 
   ENDIF 
   IF (search) THEN 
      IF (line (i:i) .eq.'#'.or.line (i:i) .eq.'!') THEN 
         line (i:ll) = ' ' 
         ll = i - 1 
      ENDIF 
   ENDIF 
ENDDO 
ll = LEN_TRIM(line)
IF(ll>0) THEN
   DO WHILE(line(ll:ll)==TAB)      !  Remove trailing TABs
      line(ll:ll) = ' '
      ll = LEN_TRIM(line)
   ENDDO
ENDIF
!                                                                       
END SUBROUTINE remove_comment                 
!
!*******************************************************************************
!
SUBROUTINE datum () 
!                                                                       
USE times_mod
!
IMPLICIT none 
!                                                                       
CHARACTER(LEN=8) :: date 
CHARACTER(LEN=10):: time 
CHARACTER(LEN=5) :: zone 
INTEGER, DIMENSION(8) :: values
!                                                                       
CALL DATE_AND_TIME (date, time, zone, values) 
int_date (1) = values (1) 
int_date (2) = values (2) 
int_date (3) = values (3) 
int_time (1) = values (5) 
int_time (2) = values (6) 
int_time (3) = values (7) 
millisec = values (8) 
midnight = values (5) * 3600 * 1000 + values (6) * 60 * 1000 +    &
values (7) * 1000 + values (8)                                    
!                                                                       
CALL fdate (f_date) 
!                                                                       
END SUBROUTINE datum                          
!
!*****7**************************************************************** 
!
SUBROUTINE datum_intrinsic () 
!                                                                       
USE times_mod
!                                                                       
IMPLICIT none 
!                                                                       
CHARACTER(LEN=8)  :: date 
CHARACTER(LEN=10) :: time 
CHARACTER(LEN=5)  :: zone 
INTEGER, DIMENSION(8) :: values
!                                                                       
CALL DATE_AND_TIME (date, time, zone, values) 
int_date (1) = values (1) 
int_date (2) = values (2) 
int_date (3) = values (3) 
int_time (1) = values (5) 
int_time (2) = values (6) 
int_time (3) = values (7) 
millisec = values (8) 
midnight = values (5) * 3600 * 1000 + values (6) * 60 * 1000 +    &
values (7) * 1000 + values (8)                                    
!                                                                       
f_date = ' ' 
f_date (1:8) = date 
f_date (9:18) = time 
f_date (19:24)  = '    ' 
!                                                                       
END SUBROUTINE datum_intrinsic                
!
!*****7**************************************************************** 
!
SUBROUTINE holecwd (cwd, dummy) 
!-
! Get current directory
!  generic version, calls compiler specific version
!                                                                       
USE sys_compiler
!
IMPLICIT none 
!                                                                       
CHARACTER(LEN=*), INTENT(OUT) :: cwd 
INTEGER         , INTENT(OUT) :: dummy 
!                                                                       
CALL sys_holecwd (cwd, dummy )
!                                                                       
END SUBROUTINE holecwd                        
!
!*****7**************************************************************** 
!
SUBROUTINE holeenv (request, answer) 
!                                                                       
USE lib_length
!
IMPLICIT none 
!                                                                       
CHARACTER(LEN=*), INTENT(IN)  :: request 
CHARACTER(LEN=*), INTENT(OUT) :: answer 
INTEGER length 
!                                                                       
answer = ' ' 
length = len_str(request)
!
CALL GET_ENVIRONMENT_VARIABLE (request(1:length), answer) 
!                                                                       
END SUBROUTINE holeenv
!
!*****7*****************************************************************
!
SUBROUTINE file_info (ifile)
!+                                                                      
!  Gets and stores file modification time
!  Generic modification, calls compiler specifc version
!+
!                                                                       
USE sys_compiler
!
IMPLICIT none
!                                                                       
INTEGER, INTENT(IN) :: ifile
!
CALL sys_file_info (ifile)
!
END SUBROUTINE file_info
!
!*****7*****************************************************************
!
SUBROUTINE file_info_disk (filename, fmodt) 
!+                                                                      
!     Gets and stored sile modification time                            
!-                                                                      
USE errlist_mod 
USE times_mod
USE sys_compiler
IMPLICIT none 
!                                                                       
!                                                                       
CHARACTER(LEN=*) , INTENT(IN)  :: filename
CHARACTER(LEN=24), INTENT(OUT) :: fmodt
!
INTEGER buff (13) 
CHARACTER :: ctime
!
fmodt = ' '
CALL stat (filename, buff, ier_num) 
IF (ier_num.ne.0) RETURN 
fmodt = ctime (buff (10) ) 
!                                                                       
END SUBROUTINE file_info_disk
!
!*****7*****************************************************************
!
SUBROUTINE do_getargs (iarg, arg, narg) 
!+                                                                      
!     Get command line parameters                                       
!-                                                                      
IMPLICIT none 
!                                                                       
INTEGER, INTENT(IN ) :: narg
INTEGER, INTENT(OUT) :: iarg
CHARACTER(LEN=*), DIMENSION(narg), INTENT(OUT) :: arg
!
INTEGER :: i 
!                                                                       
iarg = COMMAND_ARGUMENT_COUNT()
iarg = min (iarg, narg) 
!                                                                       
IF (iarg.gt.0) then 
   DO i = 1, iarg 
   CALL GET_COMMAND_ARGUMENT (i, arg (i) ) 
   ENDDO 
ENDIF 
!                                                                       
END SUBROUTINE do_getargs                     
!
!*****7*****************************************************************
!
SUBROUTINE do_operating_comm (command) 
!+                                                                      
!     Executes operating system command                                 
!-                                                                      
USE errlist_mod 
USE param_mod 
USE precision_mod
USE prompt_mod
IMPLICIT none 
!                                                                       
!                                                                       
CHARACTER (LEN=*), INTENT(IN) :: command 
CHARACTER(LEN=PREC_STRING) :: message
INTEGER :: exit_msg
integer :: ier_cmd
!
INTEGER length
!                                                                       
!     CALL system (command(1:len_str(command)), ier_num) 
length = LEN_TRIM(command)
call execute_command_line (command(1:length), wait=.true., &
     CMDSTAT=ier_cmd, CMDMSG=message, EXITSTAT=exit_msg) 
IF (ier_num.eq.0) then 
   ier_typ = ER_NONE 
ELSE 
   WRITE ( output_io, 2000) ier_num 
   WRITE ( output_io, 2010) command(1:MIN(38,length))
   res_para (0) = - 1 
   res_para (1) = - 5 
   res_para (2) = ER_COMM 
   res_para (3) = ier_num 
   ier_num = - 5 
   ier_typ = ER_COMM 
ENDIF 
!                                                                       
 2000 FORMAT    (' ****SYST**** Operating System/Shell Error Number:',  &
     &             i5,  ' ****')                                        
 2010 FORMAT    (' ****SYST**** ',a38,' ****')
!
END SUBROUTINE do_operating_comm              
!
!*****7*****************************************************************
!
SUBROUTINE do_chdir (dir, ld, echo) 
!+                                                                      
!     Changes working directory ..                                      
!-                                                                      
USE build_name_mod
USE errlist_mod 
USE envir_mod 
USE lib_length
USE precision_mod
USE get_params_mod
USE prompt_mod 
USE string_convert_mod
USE sys_compiler
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER, PARAMETER :: MAXW = 20 
CHARACTER(LEN=*), INTENT(INOUT) :: dir 
INTEGER         , INTENT(INOUT) :: ld
LOGICAL         , INTENT(IN)    :: echo 
!
CHARACTER(LEN=PREC_STRING) :: cwd                    ! to wotk with locally
CHARACTER(LEN=PREC_STRING) :: cwd_echo               ! For display only
CHARACTER(LEN=PREC_STRING), DIMENSION(MAXW) :: cpara
CHARACTER(LEN=PREC_STRING) :: string
CHARACTER(LEN=   1)        :: drive
CHARACTER(LEN=   6)        :: mntdrive
INTEGER                    :: i
INTEGER                    :: lld
INTEGER, DIMENSION(MAXW)   :: lpara
INTEGER                    :: ianz 
INTEGER                    :: dummy 
REAL(KIND=PREC_DP), DIMENSION(MAXW) :: werte
!                                                                       
!                                                                       
IF (dir.eq.' ') THEN 
   CALL sys_holecwd(cwd, dummy)
   ld  = LEN_TRIM(cwd) 
   dir = cwd                                    ! Update user directory
   cwd_echo = cwd                               ! Make copy for display
   IF(operating==OS_LINUX_WSL) THEN
      lld = LEN_TRIM(cwd_echo)
      DO i=1, lld
         IF(cwd(i:i)=='/') cwd(i:i)='\'
      ENDDO
      IF(cwd(lld:lld) /= '\') THEN
         cwd(lld+1:lld+1) = '\'
         lld = lld + 1
      ENDIF
      IF(cwd(1:5) == '\mnt\') THEN
         drive = cwd(6:6)
         CALL do_cap(drive)
         cwd_echo = drive // ':' // cwd(7:lld)
      ENDIF
   ENDIF
   IF (echo) then 
      WRITE ( *, 1000) cwd_echo(1:LEN_TRIM(cwd_echo)) 
   ENDIF 
ELSE 
   ld = -ld
   cwd = dir
   CALL get_params (cwd, ianz, cpara, lpara, maxw, ld) 
   CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
   IF(operating==OS_LINUX_WSL) THEN
      string = cpara(1)
      IF(string(2:2)==':') THEN
         drive = string(1:1)
         CALL do_low(drive)
         cpara(1) = '/mnt/' // drive // string(3:LEN_TRIM(string))
         mntdrive = '/mnt/' // drive
         IF(drive /= 'c') THEN  !     For all other drives but C:
            CALL mount_drive(drive)
         ENDIF
      ENDIF
      lpara(1) = LEN_TRIM(cpara(1))
      DO i=1,lpara(1)
         IF(cpara(1)(i:i)=='\') cpara(1)(i:i) = '/'
      ENDDO
   ENDIF
   cwd = cpara(1)(1:lpara (1) ) 
   ld  = lpara(1) 
   CALL sys_chdir(cwd(1:ld), ier_num) 
   IF (ier_num.eq.0) then 
      ier_typ = ER_NONE 
      CALL sys_holecwd(cwd, dummy)
      ld = LEN_TRIM(cwd) 
      IF(cwd(ld:ld) /= '/' ) THEN
         cwd = cwd(1:ld) // '/'
         ld  = ld + 1
      ENDIF
      cwd_echo = cwd
      IF(operating==OS_LINUX_WSL) THEN
         lld = LEN_TRIM(cwd_echo)
         DO i=1, lld
            IF(cwd_echo(i:i)=='/') cwd_echo(i:i)='\'
         ENDDO
         IF(cwd_echo(1:5) == '\mnt\') THEN
            drive = cwd_echo(6:6)
            CALL do_cap(drive)
            cwd_echo = drive // ':' // cwd_echo(7:lld)
         ENDIF
      ELSE
         cwd_echo = cwd
      ENDIF
      IF(echo) WRITE (output_io, 1000) cwd_echo(1:LEN_TRIM(cwd_echo))
   ELSE 
      WRITE ( *, 2000) ier_num 
      ier_num = - 5 
      ier_typ = ER_COMM 
   ENDIF 
ENDIF 
!
current_dir   = cwd(1:LEN(current_dir))             ! Update the global variable
current_dir_l = MIN(ld, LEN(current_dir))
dir = cwd
ld  = LEN_TRIM(dir)
!                                                                       
 1000 FORMAT    (' ------ > Working directory : ',a) 
 2000 FORMAT    (' ****SYST**** Unable to change directory - Error :',  &
     &             i5,  ' ****')                                        
END SUBROUTINE do_chdir                       
!
!*****7*****************************************************************
!
SUBROUTINE display_dir(cwd)
!
USE envir_mod
USE string_convert_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: cwd
!
CHARACTER(LEN=   1)        :: drive
INTEGER :: i, ld
!
IF(operating==OS_LINUX_WSL) THEN
   ld = LEN_TRIM(cwd)
   DO i=1, ld
      IF(cwd(i:i)=='/') cwd(i:i)='\'
   ENDDO
   IF(cwd(ld:ld) /= '\') THEN
      cwd(ld+1:ld+1) = '\'
      ld = ld + 1
   ENDIF
   IF(cwd(1:5) == '\mnt\') THEN
      drive = cwd(6:6)
      CALL do_cap(drive)
      cwd = drive // ':' // cwd(7:ld)
      ld = LEN_TRIM(cwd)
   ENDIF
ENDIF
!
END SUBROUTINE display_dir
!
!*****7*****************************************************************
!
SUBROUTINE do_cwd (dir, ld) 
!+                                                                      
!     Gets the current working directory and stores the result in       
!     dir                                                               
!-                                                                      
USE errlist_mod 
USE lib_length
USE precision_mod
USE sys_compiler
!
IMPLICIT none 
!                                                                       
CHARACTER(LEN=*), INTENT(OUT) :: dir 
INTEGER         , INTENT(OUT) :: ld
!
CHARACTER(LEN=PREC_STRING) :: cwd 
INTEGER ::  dummy 
!                                                                       
CALL sys_holecwd( cwd, dummy)
ld = len_str (cwd)
dir = cwd (1:ld) 
!                                                                       
END SUBROUTINE do_cwd                         
!
!*****7*****************************************************************
!
SUBROUTINE do_del_file (name, wait_stat) 
!+                                                                      
!     Deletes a file                                                    
!-                                                                      
USE blanks_mod
USE errlist_mod 
USE lib_length
USE precision_mod
IMPLICIT none 
!                                                                       
!                                                                       
CHARACTER (LEN=*), INTENT(INOUT) :: name 
logical          , intent(in), optional :: wait_stat   ! User provided true/false
!
CHARACTER(LEN=PREC_STRING) :: command 
CHARACTER(LEN=PREC_STRING) :: message
INTEGER :: laenge
integer :: exit_msg 
integer :: ier_cmd
logical :: wait_local
!                                                                       
!
wait_local = .true.             ! Default to wait
if(present(wait_stat)) then     ! use optional user parameter
   wait_local = wait_stat
endif
laenge = len_str (name) 
CALL rem_bl (name, laenge) 
command = ' ' 
!                                                                       
!     sun,hp                                                            
!                                                                       
command (1:6) = 'rm -f ' 
command (7:7 + laenge) = name (1:laenge) 
call execute_command_line (command(1:7+laenge), wait=wait_local,  &
     CMDSTAT=ier_cmd, CMDMSG=message, exitstat=exit_msg) 
IF (ier_num.eq.0) then 
   ier_typ = ER_NONE 
ELSE 
   WRITE ( *, 1900) name
   WRITE ( *, 2000) ier_num 
   ier_num = - 5 
   ier_typ = ER_COMM 
ENDIF 
!                                                                       
 1900 FORMAT('Could not delete file: ', a)
 2000 FORMAT(' ****SYST****Operating System/Shell Error Number:',i5,'****')                                         
!
END SUBROUTINE do_del_file                    
!
!*****7*****************************************************************
!
SUBROUTINE do_rename_file (nameold, namenew) 
!+                                                                      
!     Renames the file <nameold> to <namenew>                           
!-                                                                      
USE blanks_mod
USE errlist_mod 
USE lib_length
USE precision_mod
!                                                                       
IMPLICIT none 
!                                                                       
CHARACTER (LEN=*), INTENT(INOUT) :: nameold
CHARACTER (LEN=*), INTENT(OUT)   :: namenew 
!
CHARACTER(LEN=80)          :: command 
CHARACTER(LEN=PREC_STRING) :: message
INTEGER :: lo, ln
integer :: exit_msg
integer :: ier_cmd
!                                                                       
lo = len_str (nameold) 
CALL rem_bl (nameold, lo) 
ln = len_str (namenew) 
CALL rem_bl (namenew, ln) 
command = ' ' 
!                                                                       
!     sun,hp                                                            
!                                                                       
command (1:3) = 'mv ' 
command (4:3 + lo) = nameold (1:lo) 
command (4 + lo:4 + lo) = ' ' 
command (5 + lo:4 + lo + ln) = namenew (1:ln) 
!                                                                       
call execute_command_line (command(1:4+lo+ln), wait=.true., &
     CMDSTAT=ier_cmd, CMDMSG=message, exitstat=exit_msg) 
IF (ier_num.eq.0) then 
   ier_typ = ER_NONE 
ELSE 
   WRITE ( *, 2000) ier_num 
   ier_num = - 5 
   ier_typ = ER_COMM 
ENDIF 
!                                                                       
 2000 FORMAT('****SYST****Operating System/Shell Error Number:',i5,'****')                                         
!
END SUBROUTINE do_rename_file                 
!
!*****7*****************************************************************
!
REAL(kind=PREC_SP) FUNCTION seknds_sp (s) 
!-                                                                      
!     returns the elapsed user time since the last call or start (s=0)  
!     im seconds                                                        
! Single precision version
!                                                                       
use precision_mod
IMPLICIT none 
!                                                                       
REAL(kind=PREC_SP), INTENT(IN) :: s 
REAL(kind=PREC_DP) :: current
!                                                                       
CALL CPU_TIME(current) 
seknds_sp = current - s
!
END FUNCTION seknds_sp                           
!
!*****7*****************************************************************
!
REAL(kind=PREC_DP) FUNCTION seknds_dp (s) 
!-                                                                      
!     returns the elapsed user time since the last call or start (s=0)  
!     im seconds                                                        
! double precision version
!                                                                       
use precision_mod
IMPLICIT none 
!                                                                       
REAL(kind=PREC_DP), INTENT(IN) :: s 
REAL(kind=PREC_DP) :: current
!                                                                       
CALL CPU_TIME(current) 
seknds_dp = current - s
!
END FUNCTION seknds_dp                           
!
!*****7*****************************************************************
!
SUBROUTINE open_def (idef) 
!                                                                       
USE envir_mod 
USE errlist_mod 
USE precision_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER, INTENT(IN) :: idef
!
CHARACTER(LEN=PREC_STRING) :: dffile 
CHARACTER(LEN=3), PARAMETER :: accessm = 'old'
INTEGER :: l1, l2 
LOGICAL lread 
!                                                                       
DATA dffile / ' ' / 
!                                                                       
!     try in present directory                                          
!                                                                       
dffile = deffile 
lread = .true. 
CALL oeffne (idef, dffile, accessm)
IF (ier_num.ne.0) then 
!                                                                       
!     not found, try in home directory                                  
!                                                                       
   dffile (1:home_dir_l) = home_dir 
   l1 = home_dir_l + 1 
   l2 = home_dir_l + 1 + deffile_l + 1 
   dffile (l1:l2) = '/'//deffile 
   CALL oeffne (idef, dffile, accessm)
   IF (ier_num.ne.0) then 
!                                                                       
!     not found, try in home/bin directory                              
!                                                                       
      dffile (1:home_dir_l) = home_dir 
      l1 = home_dir_l + 1 
      l2 = l1 + 4 
      dffile (l1:l2) = '/bin/' 
      l1 = l2 + 1 
      l2 = l1 + deffile_l 
      dffile (l1:l2) = deffile 
      CALL oeffne (idef, dffile, accessm)
   ENDIF 
ENDIF 
!
END SUBROUTINE open_def                       
!
!*****7***************************************************************  
!
SUBROUTINE oeffne (inum, in_datei, stat) 
!-                                                                      
!     opens a file 'datei' with status 'stat' and unit 'inum'           
!     if lread = .true. the file is opened readonly                     
!     The readonly parameter is VMS specific, not implemented           
!     on SUN FORTRAN                                                    
!+                                                                      
USE envir_mod 
USE errlist_mod 
USE lib_length
USE precision_mod
USE times_mod
!                                                                       
IMPLICIT none 
!                                                                       
CHARACTER (LEN=* ), INTENT(IN) :: in_datei
CHARACTER (LEN=* ), INTENT(IN)    :: stat 
INTEGER           , INTENT(IN)    :: inum
!
CHARACTER(LEN=PREC_STRING) :: datei
CHARACTER (LEN=PREC_STRING)  :: line 
CHARACTER (LEN=PREC_STRING)  :: message
INTEGER :: ios 
INTEGER :: l_datei 
LOGICAL :: lda
!                                                                       
message = ' '
!                                                                       
datei = in_datei
ios     = 0
l_datei = len_str (datei) 
IF (l_datei.gt.0) then 
   IF (datei (1:1) .eq.'~') then 
      line = home_dir (1:home_dir_l) //datei (2:l_datei) 
      datei = line 
   ENDIF 
   ier_num = - 2 
   ier_typ = ER_IO 
   IF (stat.eq.'unknown') then 
!
      OPEN (UNIT=inum, FILE=datei , STATUS=stat, IOSTAT=ios, IOMSG=message)
      IF(ios==0) THEN
         ier_num = 0 
         ier_typ = ER_NONE 
      ELSE
         ier_num = -2
         ier_typ = er_io
         ier_msg(3) = message(1:80)
      ENDIF
   ELSE 
      INQUIRE (file = datei, exist = lda) 
      IF (stat.eq.'old') then 
         IF (lda) then 
            OPEN (inum, FILE = datei, STATUS = stat,  &
            IOSTAT = ios,IOMSG=message)
            IF(ios==0) THEN
               ier_num = 0 
               ier_typ = ER_NONE 
               CALL file_info (inum) 
            ELSE
               ier_num = -2
               ier_typ = er_io
               ier_msg(3) = message(1:80)
            ENDIF
         ELSEIF (.not.lda) then 
            ier_num = - 1 
            ier_typ = ER_IO 
         ENDIF 
      ELSEIF (stat.eq.'new') then 
         IF (.not.lda) then 
            OPEN (inum, file = datei, status = stat,  &
            iostat = ios,IOMSG=message)
            IF(ios==0) THEN
               ier_num = 0 
               ier_typ = ER_NONE 
            ELSE
               ier_num = -2
               ier_typ = er_io
               ier_msg(3) = message(1:80)
            ENDIF
         ELSEIF (lda) then 
            ier_num = - 4 
            ier_typ = ER_IO 
         ENDIF 
      ENDIF 
   ENDIF 
ELSE 
   ier_num = - 14 
   ier_typ = ER_IO 
ENDIF 
!
END SUBROUTINE oeffne                         
!
!*****7***************************************************************  
!
SUBROUTINE oeffne_append (inum, in_datei, stat) 
!-                                                                      
!     opens a file 'datei' with status 'stat' and unit 'inum'           
!     in append mode.                                                   
!     if lread = .true. the file is opened readonly                     
!     The readonly parameter is VMS specific, not implemented           
!     on SUN FORTRAN                                                    
!+                                                                      
USE envir_mod 
USE errlist_mod 
USE lib_length
USE precision_mod
!                                                                       
IMPLICIT none 
!                                                                       
CHARACTER (LEN=*), INTENT(IN) :: in_datei
CHARACTER (LEN=*), INTENT(IN) :: stat 
INTEGER          , INTENT(IN) :: inum 
CHARACTER(LEN=PREC_STRING) :: line 
CHARACTER(LEN=PREC_STRING) :: message 
CHARACTER(LEN=PREC_STRING) :: datei 
INTEGER :: ios 
INTEGER :: l_datei 
LOGICAL :: lda
!                                                                       
!                                                                       
datei = in_datei
l_datei = len_str (datei) 
IF (l_datei.gt.0) then 
   IF (datei (1:1) .eq.'~') then 
      line = home_dir (1:home_dir_l) //datei (2:l_datei) 
      datei = line 
   ENDIF 
   ier_num = - 2 
   ier_typ = ER_IO 
   IF (stat.eq.'unknown') then 
      OPEN (inum, file = datei, status = stat, position =         &
      'append', iostat = ios, IOMSG = message)
            IF(ios==0) then
         ier_num = 0 
         ier_typ = ER_NONE 
      ELSE
         ier_num = -2
         ier_typ = ER_IO
         ier_msg(3) = message(1:80)
      ENDIF
   ELSE 
      INQUIRE (file = datei, exist = lda) 
      IF (stat.eq.'old') then 
         IF (lda) then 
            OPEN (inum, file = datei, status = stat, position =   &
            'append', iostat = ios, IOMSG = message)
            IF(ios==0) then
               ier_num = 0 
               ier_typ = ER_NONE 
            ELSE
               ier_num = -2
               ier_typ = ER_IO
               ier_msg(3) = message(1:80)
            ENDIF
         ELSEIF (.not.lda) then 
            ier_num = - 1 
            ier_typ = ER_IO 
         ENDIF 
      ELSEIF (stat.eq.'new') then 
         IF (.not.lda) then 
            OPEN (inum, file = datei, status = stat, position =   &
            'append', iostat = ios, IOMSG = message)
            IF(ios==0) then
               ier_num = 0 
               ier_typ = ER_NONE 
            ELSE
               ier_num = -2
               ier_typ = ER_IO
               ier_msg(3) = message(1:80)
            ENDIF
         ELSEIF (lda) then 
            ier_num = - 4 
            ier_typ = ER_IO 
         ENDIF 
      ENDIF 
   ENDIF 
ELSE 
   ier_num = - 14 
   ier_typ = ER_IO 
ENDIF 
!
END SUBROUTINE oeffne_append                  
!
!*****7***************************************************************  
!
SUBROUTINE do_sleep (seconds) 
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER,INTENT(IN) :: seconds 
!                                                                       
CALL sleep (seconds) 
!                                                                       
END SUBROUTINE do_sleep                       
!
!*****7***************************************************************  
!
SUBROUTINE do_fexist (zeile, lp, lout) 
!                                                                       
USE build_name_mod
USE errlist_mod 
USE get_params_mod
USE param_mod 
USE precision_mod
USE prompt_mod 
USE sys_compiler
IMPLICIT none 
!                                                                       
!                                                                       
INTEGER, PARAMETER :: maxw = 10 
!                                                                       
CHARACTER (LEN=*), INTENT(INOUT) :: zeile 
INTEGER , INTENT(INOUT)::lp 
LOGICAL , INTENT(IN)   ::lout
!
CHARACTER(LEN=PREC_STRING), DIMENSION(MAXW) :: cpara
REAL(KIND=PREC_DP), DIMENSION(MAXW) :: werte
INTEGER, DIMENSION(MAXW) :: lpara
INTEGER :: ianz 
LOGICAL :: lexist 
!                                                                       
CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
IF (ier_num.ne.0) return 
!                                                                       
IF (ianz.ge.1) then 
   CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
   INQUIRE (file = cpara (1) (1:lpara (1) ), exist = lexist) 
   IF (lexist) then 
      res_para (0) = 1 
      res_para (1) = 1 
      IF(lout) WRITE (output_io, 1000) cpara (1) (1:lpara (1) ) 
   ELSE 
      CALL sys_inquire_directory(MAXW, cpara, lpara,lexist)
!     INQUIRE (directory = cpara (1) (1:lpara (1) ), exist = lexist)
      IF (lexist) then
         res_para (0) = 1
         res_para (1) = 1
         WRITE (output_io, 1000) cpara (1) (1:lpara (1) )
      ELSE 
         res_para (0) = 1 
         res_para (1) = 0 
         IF(lout) WRITE (output_io, 1100) cpara (1) (1:lpara (1) ) 
      ENDIF 
   ENDIF 
ELSE 
   ier_num = - 6 
   ier_typ = ER_COMM 
ENDIF 
!                                                                       
 1000 FORMAT     (' ------ > File ',a,' exists ...') 
 1100 FORMAT     (' ------ > File ',a,' does NOT exist ...') 
END SUBROUTINE do_fexist                      
!
!*****7***********************************************************      
!
LOGICAL FUNCTION is_nan(x)
!
use precision_mod
!
REAL(kind=PREC_DP), INTENT(IN) :: x
!
is_nan = isnan(x)
!
END FUNCTION is_nan
!
!*****7***********************************************************      
!
INTEGER FUNCTION lib_f90_getpid()
!
! Interface to the compiler dependent GETPID routines
!
lib_f90_getpid = getpid()
!
END FUNCTION lib_f90_getpid
!
!*****7***********************************************************      
!
SUBROUTINE mount_drive(indrive)
!
! Mount a removable drive at WSL
!
USE blanks_mod
USE envir_mod
USE errlist_mod
USE prompt_mod
USE precision_mod
USE string_convert_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=1), INTENT(IN) :: indrive
!
INTEGER, PARAMETER :: IRD = 69
CHARACTER(LEN=1) :: drive
CHARACTER(LEN=6) :: mntdrive
CHARACTER(LEN=PREC_STRING) :: string
CHARACTER(LEN=PREC_STRING) :: message
integer :: exit_msg
integer :: ier_cmd
INTEGER :: length
INTEGER :: ios
LOGICAL :: lda
LOGICAL :: lmounted
!
IF(operating == OS_LINUX_WSL) THEN    ! Only for WSL
   drive = indrive
!
   length = LEN_TRIM(drive)
   CALL rem_bl(drive, length)
   CALL do_low(drive)
   mntdrive = '/mnt/' // drive(1:1)
!
   INQUIRE(FILE=mntdrive,EXIST=lda)
   IF(.NOT. lda) THEN     ! Need to build a '/mnt/X' file
      WRITE(output_io,*) 'Access to drives requires Ubuntu password'
      string='sudo mkdir -p '// mntdrive
      call execute_command_line (string(1:LEN_TRIM(string)), wait=.true., &
                              CMDSTAT=ier_cmd, CMDMSG=message, exitstat=exit_msg) 
      CALL do_cap(drive)
      string = 'sudo mount -t drvfs ' // drive // ': ' // mntdrive
      call execute_command_line (string(1:LEN_TRIM(string)), wait=.true., &
                                 CMDSTAT=ier_cmd, CMDMSG=message, exitstat=exit_msg) 
   ELSE
      CALL do_cap(drive)
      OPEN(UNIT=IRD, FILE='/proc/mounts', STATUS='OLD', ACTION='READ')
      lmounted = .FALSE.
      find: DO
         READ(IRD,'(a)', IOSTAT=IOS) string
         IF(IS_IOSTAT_END(ios)) EXIT find
         length = LEN_TRIM(string)
         CALL rem_leading_bl(string,length)
         IF(string(1:2) == drive//':') THEN
            lmounted=.TRUE.
            EXIT find
         ENDIF
      ENDDO find
      CLOSE(UNIT=IRD)
      IF(.NOT.lmounted) THEN
         WRITE(output_io,*) 'Access to drives may require Ubuntu password'
         string = 'sudo mount -t drvfs ' // drive // ': ' // mntdrive
         call execute_command_line (string(1:LEN_TRIM(string)), wait=.true., &
                                    CMDSTAT=ier_cmd, CMDMSG=message, exitstat=exit_msg) 
      ENDIF
   ENDIF
ENDIF
!
END SUBROUTINE mount_drive
!
!*****7***********************************************************      
!
SUBROUTINE umount_drive(indrive)
!
! un mount a WSL removable drive
!
USE blanks_mod
USE envir_mod
USE errlist_mod
USE precision_mod
USE prompt_mod
USE string_convert_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=1), INTENT(IN) :: indrive
!
CHARACTER(LEN=1) :: drive
CHARACTER(LEN=6) :: mntdrive
CHARACTER(LEN=PREC_STRING) :: string
CHARACTER(LEN=PREC_STRING) :: message
integer :: exit_msg
integer :: ier_cmd
INTEGER :: length
LOGICAL :: lda
!
IF(operating == OS_LINUX_WSL) THEN    ! Only for WSL
   drive = indrive
!
   length = LEN_TRIM(drive)
   CALL rem_bl(drive, length)
   CALL do_low(drive)
   mntdrive = '/mnt/' // drive(1:1)
!
   INQUIRE(FILE=mntdrive,EXIST=lda)
!
   IF(lda) THEN       ! Drive letter exists, otherwise silently leave
      WRITE(output_io,*) 'Access to drives may require Ubuntu password'
      string = 'sudo umount ' // mntdrive
      call execute_command_line (string(1:LEN_TRIM(string)), wait=.true., &
           CMDSTAT=ier_cmd, CMDMSG=message, exitstat=exit_msg) 
   ENDIF
ENDIF
!
END SUBROUTINE umount_drive
!
!*****7***********************************************************      
!
subroutine set_rw(infile, length)
!-
!  set user rw  privillege to  a file
!+
!
use precision_mod
!
implicit none
!
character(len=*), intent(in) :: infile
integer         , intent(in) :: length
!
character(len=PREC_STRING) :: string
character(len=PREC_STRING) :: message
integer :: exit_msg
integer :: ier_cmd
!
string = 'chmod u+rw ' // infile(1:length)
call execute_command_line (string(1:LEN_TRIM(string)), wait=.true., &
           CMDSTAT=ier_cmd, CMDMSG=message, exitstat=exit_msg)
!
end subroutine set_rw
!
!*****7**************************************************************** 
!
END MODULE support_mod
