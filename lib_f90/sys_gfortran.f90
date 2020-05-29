MODULE sys_compiler
!
CONTAINS
!                                                                       
!     This file contains subroutines for:                               
!     Compiler specific routines GNU gfortran version                   
!
!*****7**************************************************************** 
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
      IMPLICIT none 
!                                                                       
!                                                                       
      CHARACTER(8) date 
      CHARACTER(10) time 
      CHARACTER(5) zone 
      INTEGER values (8) 
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
!*****7**************************************************************** 
      SUBROUTINE holecwd (cwd, dummy) 
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER(LEN=*) cwd 
      INTEGER dummy 
      INTEGER  :: getcwd
!                                                                       
      dummy = getcwd (cwd )
!                                                                       
      END SUBROUTINE holecwd                        
!*****7**************************************************************** 
      SUBROUTINE holeenv (request, answer) 
!                                                                       
USE lib_length
      IMPLICIT none 
!                                                                       
      CHARACTER(LEN=*) request 
      CHARACTER(LEN=*) answer 
      INTEGER length 
!                                                                       
      answer = ' ' 
      length = len_str(request)
!     CALL getenv (request(1:length), answer) 
      CALL GET_ENVIRONMENT_VARIABLE (request(1:length), answer) 
!                                                                       
      END SUBROUTINE holeenv                        
!*****7*****************************************************************
      SUBROUTINE file_info (ifile) 
!+                                                                      
!     Gets and stored sile modification time                            
!-                                                                      
      USE errlist_mod 
      USE times_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER ifile 
      INTEGER buff (13) 
      CHARACTER :: ctime
!
      CALL fstat (ifile, buff, ier_num) 
      IF (ier_num.ne.0) return 
      f_modt = ctime (buff (10) ) 
!                                                                       
      END SUBROUTINE file_info                      
!*****7*****************************************************************
      SUBROUTINE file_info_disk (filename, fmodt) 
!+                                                                      
!     Gets and stored sile modification time                            
!-                                                                      
      USE errlist_mod 
      USE times_mod
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
!*****7*****************************************************************
      SUBROUTINE do_getargs (iarg, arg, narg) 
!+                                                                      
!     Get command line parameters                                       
!-                                                                      
      IMPLICIT none 
!                                                                       
      INTEGER narg, iarg, i 
      CHARACTER (LEN=*) arg (narg) 
!                                                                       
!     iarg = iargc () 
      iarg = COMMAND_ARGUMENT_COUNT()
      iarg = min (iarg, narg) 
!                                                                       
      IF (iarg.gt.0) then 
         DO i = 1, iarg 
!        CALL getarg (i, arg (i) ) 
         CALL GET_COMMAND_ARGUMENT (i, arg (i) ) 
         ENDDO 
      ENDIF 
!                                                                       
      END SUBROUTINE do_getargs                     
!*****7*****************************************************************
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
!
      INTEGER length
!                                                                       
!     CALL system (command(1:len_str(command)), ier_num) 
      length = LEN_TRIM(command)
      CALL EXECUTE_COMMAND_LINE (command(1:length), CMDSTAT=ier_num, CMDMSG=message, EXITSTAT=exit_msg) 
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
      END SUBROUTINE do_operating_comm              
!*****7*****************************************************************
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
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER MAXW 
      PARAMETER (MAXW = 20) 
      CHARACTER (LEN=*) dir 
      CHARACTER(LEN=PREC_STRING) :: cwd 
      CHARACTER(LEN=PREC_STRING) :: cpara (MAXW) 
      CHARACTER(LEN=PREC_STRING) :: string
      CHARACTER(LEN=   1) :: drive
      INTEGER :: i
      INTEGER lpara (MAXW) 
      INTEGER ianz 
      INTEGER ld, dummy 
      LOGICAL echo 
      REAL(KIND=PREC_DP) :: werte (MAXW) 
!                                                                       
      INTEGER  :: getcwd
!                                                                       
      IF (dir.eq.' ') then 
         dummy = getcwd (cwd )
         ld = len_str (cwd) 
         dir = cwd 
         IF (echo) then 
            WRITE ( *, 1000) cwd (1:ld) 
         ENDIF 
      ELSE 
         ld = -ld
         CALL get_params (dir, ianz, cpara, lpara, maxw, ld) 
         CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
         IF(operating==OS_LINUX_WSL) THEN
            string = cpara(1)
            IF(string(2:2)==':') THEN
               drive = string(1:1)
               CALL do_low(drive)
               cpara(1) = '/mnt/' // drive // string(3:LEN_TRIM(string))
            ENDIF
            lpara(1) = LEN_TRIM(cpara(1))
            DO i=1,lpara(1)
               IF(cpara(1)(i:i)=='\') cpara(1)(i:i) = '/'
            ENDDO
         ENDIF
         dir = cpara (1) (1:lpara (1) ) 
         ld = lpara (1) 
         CALL chdir (dir (1:ld), ier_num) 
         IF (ier_num.eq.0) then 
            ier_typ = ER_NONE 
            dummy = getcwd (cwd )
            ld = len_str (cwd) 
            IF(cwd(ld:ld) /= '/' ) THEN
               cwd = cwd(1:ld) // '/'
               ld  = ld + 1
            ENDIF
            IF(echo) WRITE (output_io, 1000) cwd (1:ld)
            current_dir   = cwd(1:LEN(current_dir))
            current_dir_l = MIN(ld, LEN(current_dir))
         ELSE 
            WRITE ( *, 2000) ier_num 
            ier_num = - 5 
            ier_typ = ER_COMM 
         ENDIF 
      ENDIF 
!                                                                       
 1000 FORMAT    (' ------ > Working directory : ',a) 
 2000 FORMAT    (' ****SYST**** Unable to change directory - Error :',  &
     &             i5,  ' ****')                                        
      END SUBROUTINE do_chdir                       
!*****7*****************************************************************
      SUBROUTINE do_cwd (dir, ld) 
!+                                                                      
!     Gets the current working directory and stores the result in       
!     dir                                                               
!-                                                                      
      USE errlist_mod 
USE lib_length
USE precision_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      CHARACTER (LEN=*) dir 
      CHARACTER(LEN=PREC_STRING) :: cwd 
      INTEGER ld, dummy 
!                                                                       
      INTEGER  :: getcwd
!                                                                       
      dummy = getcwd (cwd )
      ld = len_str (cwd)
      dir = cwd (1:ld) 
!                                                                       
      END SUBROUTINE do_cwd                         
!*****7*****************************************************************
      SUBROUTINE do_del_file (name) 
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
      CHARACTER (LEN=*) name 
      CHARACTER(LEN=PREC_STRING) :: command 
      CHARACTER(LEN=PREC_STRING) :: message
      INTEGER laenge 
!                                                                       
      laenge = len_str (name) 
      CALL rem_bl (name, laenge) 
      command = ' ' 
!                                                                       
!     sun,hp                                                            
!                                                                       
      command (1:6) = 'rm -f ' 
      command (7:7 + laenge) = name (1:laenge) 
!     CALL system (command(1:7+laenge)) 
      CALL EXECUTE_COMMAND_LINE (command(1:7+laenge), CMDSTAT=ier_num, CMDMSG=message) 
      IF (ier_num.eq.0) then 
         ier_typ = ER_NONE 
      ELSE 
         WRITE ( *, 1900) name
         WRITE ( *, 2000) ier_num 
         ier_num = - 5 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
 1900 FORMAT    ('Could not delete file: ', a)
 2000 FORMAT    (' ****SYST****Operating System/Shell Error Number:',i5,&
     &                  '****')                                         
      END SUBROUTINE do_del_file                    
!*****7*****************************************************************
      SUBROUTINE do_rename_file (nameold, namenew) 
!+                                                                      
!     Renames the file <nameold> to <namenew>                           
!-                                                                      
      USE blanks_mod
      USE errlist_mod 
USE lib_length
USE precision_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      CHARACTER (LEN=*) nameold, namenew 
      CHARACTER(80) command 
      CHARACTER(LEN=PREC_STRING) :: message
      INTEGER lo, ln 
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
!     CALL system (command(1:4+lo+ln)) 
      CALL EXECUTE_COMMAND_LINE (command(1:4+lo+ln), CMDSTAT=ier_num, CMDMSG=message) 
      IF (ier_num.eq.0) then 
         ier_typ = ER_NONE 
      ELSE 
         WRITE ( *, 2000) ier_num 
         ier_num = - 5 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
 2000 FORMAT    ('****SYST****Operating System/Shell Error Number:',i5, &
     &                  '****')                                         
      END SUBROUTINE do_rename_file                 
!*****7*****************************************************************
      REAL FUNCTION seknds (s) 
!-                                                                      
!     returns the elapsed user time since the last call or start (s=0)  
!     im seconds                                                        
!                                                                       
      IMPLICIT none 
!                                                                       
      REAL s 
!     REAL time (2), res 
      REAL :: current
!                                                                       
!     CALL etime (time, res) 
!     seknds = time (1) - s
      CALL CPU_TIME(current) 
      seknds = current - s
!
      END FUNCTION seknds                           
!*****7*****************************************************************
      SUBROUTINE open_def (idef) 
!                                                                       
      USE envir_mod 
      USE errlist_mod 
USE precision_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      CHARACTER(LEN=PREC_STRING) :: dffile 
      INTEGER idef, l1, l2 
      LOGICAL lread 
!                                                                       
      DATA dffile / ' ' / 
!                                                                       
!     try in present directory                                          
!                                                                       
      dffile = deffile 
      lread = .true. 
      CALL oeffne (idef, dffile, 'old') 
      IF (ier_num.ne.0) then 
!                                                                       
!     not found, try in home directory                                  
!                                                                       
         dffile (1:home_dir_l) = home_dir 
         l1 = home_dir_l + 1 
         l2 = home_dir_l + 1 + deffile_l + 1 
         dffile (l1:l2) = '/'//deffile 
         CALL oeffne (idef, dffile, 'old') 
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
            CALL oeffne (idef, dffile, 'old') 
         ENDIF 
      ENDIF 
      END SUBROUTINE open_def                       
!*****7***************************************************************  
      SUBROUTINE oeffne (inum, datei, stat) 
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
      IMPLICIT none 
!                                                                       
!                                                                       
      CHARACTER (LEN=* )    :: datei, stat 
      CHARACTER (LEN=PREC_STRING)  :: line 
      CHARACTER (LEN=PREC_STRING)  :: message
      INTEGER inum, ios 
      INTEGER l_datei 
      LOGICAL lda
!                                                                       
!                                                                       
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
!            line = datei(1:LEN_TRIM(datei))
!            line = ' '
!            line = 'test'
            OPEN (UNIT=inum, FILE = datei , STATUS = stat, IOSTAT =&
            ios, IOMSG=message)
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
      END SUBROUTINE oeffne                         
!*****7***************************************************************  
      SUBROUTINE oeffne_append (inum, datei, stat) 
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
      IMPLICIT none 
!                                                                       
!                                                                       
      CHARACTER (LEN=*) datei, stat 
      CHARACTER(LEN=PREC_STRING) :: line 
      CHARACTER(LEN=PREC_STRING) :: message 
      INTEGER inum, ios 
      INTEGER l_datei 
      LOGICAL lda
!                                                                       
!                                                                       
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
!*****7***************************************************************  
      SUBROUTINE do_sleep (seconds) 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER seconds 
!                                                                       
      CALL sleep (seconds) 
!                                                                       
      END SUBROUTINE do_sleep                       
!*****7***************************************************************  
      SUBROUTINE do_fexist (zeile, lp, lout) 
!                                                                       
      USE build_name_mod
      USE errlist_mod 
      USE get_params_mod
      USE param_mod 
USE precision_mod
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 10) 
!                                                                       
      CHARACTER (LEN=*), INTENT(INOUT) :: zeile 
      INTEGER , INTENT(INOUT)::lp 
      LOGICAL , INTENT(IN)   ::lout
!
      CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
      REAL(KIND=PREC_DP) :: werte (maxw) 
      INTEGER lpara (maxw)
      INTEGER ianz 
      LOGICAL lexist 
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
            res_para (0) = 1 
            res_para (1) = 0 
            IF(lout) WRITE (output_io, 1100) cpara (1) (1:lpara (1) ) 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
 1000 FORMAT     (' ------ > File ',a,' exists ...') 
 1100 FORMAT     (' ------ > File ',a,' does NOT exist ...') 
      END SUBROUTINE do_fexist                      
!*****7***********************************************************      
!
LOGICAL FUNCTION is_nan(x)
!
REAL, INTENT(IN) :: x
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
!*****7***********************************************************      
END MODULE sys_compiler
