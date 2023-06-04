MODULE appl_env_mod
!
private
!
public appl_env
public color_set_scheme
public lib_f90_init_updates
public lib_f90_test_updates
public lib_f90_update_discus
PUBLIC lib_f90_update_ubuntu
public program_files
public write_appl_env
public start_pgxwin
!
!*****7***************************************************************  
!                                                                       
CONTAINS
!
!*******************************************************************************!
!
SUBROUTINE appl_env (standalone) !, local_mpi_myid)
!-                                                                      
!     Reads environment variables, sets path for helpfile               
!     UNIX version ..                                                   
!+                                                                      
USE blanks_mod
USE errlist_mod
USE envir_mod 
use lib_config
USE lib_errlist_func
USE lib_length
USE terminal_mod
USE param_mod
USE precision_mod
USE prompt_mod 
USE string_convert_mod
USE support_mod
use sys_compiler
USE lib_do_operating_mod
!
IMPLICIT none 
!                                                                       
LOGICAL, INTENT(IN) :: standalone
!INTEGER, INTENT(IN) :: local_mpi_myid
!                                                                       
INTEGER, PARAMETER :: idef = 68
!
character(len=PREC_STRING) :: get_screen_size
character(len=PREC_STRING) :: get_screen_file
CHARACTER(len=255) :: cdummy
CHARACTER(LEN=12), DIMENSION(7), PARAMETER :: tmp_test = (/                     &
'/private/tmp', '/tmp        ', '/TMP        ',                                 &
'/var/tmp    ', '/Var/tmp    ', '/var/TMP    ', '/Var/TMP    ' /)
CHARACTER(LEN=PREC_LSTRING) :: line
CHARACTER(LEN=PREC_LSTRING) :: pathfile
CHARACTER(LEN=PREC_LSTRING) :: ufile
CHARACTER(LEN=15), PARAMETER :: etc_os_release = '/etc/os-release'
character(len=PREC_STRING) :: message
integer             :: exit_msg
integer             :: ier_cmd
INTEGER ico, ice, iii, i, j
INTEGER :: length
INTEGER :: ios ! I/O status
INTEGER pname_l 
!INTEGER, PARAMETER                        :: MAXW = 1
!CHARACTER (LEN=PREC_STRING), DIMENSION(1) :: cpara
!INTEGER                    , DIMENSION(1) :: lpara
LOGICAL :: lpresent
!logical :: lexist
!INTEGER :: lib_f90_getpid
!
IF(envir_done) RETURN
!
get_screen_size = ' '
!                                                                       
! Get the PID of the DISCUS Process
!
PID = lib_f90_getpid()
!
! Determine a temporary directory
!
tmp_dir   = '.'
tmp_dir_l = 1
find_tmp: DO i=1, 6
   line= tmp_test(i)
   length = LEN_TRIM(line)
   CALL do_fexist(line, length, .FALSE.)
   lpresent = res_para(1) == 1
!  INQUIRE(DIRECTORY=line, EXIST=lpresent)
   IF(lpresent) THEN
      tmp_dir = tmp_test(i)
      tmp_dir_l = LEN_TRIM(tmp_dir)
      EXIT find_tmp
   ENDIF
ENDDO find_tmp
!
!CALL lib_f90_init_updates
!
line = ''
CALL GET_COMMAND(line)                    ! Get complete command line
IF(line(1:23)=='discus_suite_noparallel') THEN
  start_line = 'discus_suite_noparallel'
  start_line_l = 23
ELSE
  start_line = 'discus_suite'
  start_line_l = 12
ENDIF
operating = ' ' 
color_theme = THEME_DEFAULT
INQUIRE(FILE='/proc/version', EXIST=lpresent)
IF(lpresent) THEN               ! /proc/version exists, Linux type OS
   OPEN(UNIT=idef, FILE='/proc/version', ACTION='read')
   READ(idef,'(a)') line
   CLOSE(UNIT=idef)
   CALL do_cap(line)
   IF(INDEX(line,'LINUX') > 0) THEN
      operating = OS_LINUX                  ! Assume native Linux
      color_theme = THEME_DEFAULT
      IF(INDEX(line,'MICROSOFT') > 0) THEN
         operating = OS_LINUX_WSL           ! Linux as Windows APP
         color_theme = THEME_DEFAULT
      ENDIF
         home_dir = ' '  ! Home dir needs work if MS-WSL icon in the future
         CALL get_environment_variable ('HOME', home_dir) 
         IF (home_dir.eq.' ') then 
            home_dir = '.' 
         ENDIF 
      home_dir_l = len_str (home_dir) 
!
      OPEN(idef,FILE='/etc/os-release', STATUS='old', IOSTAT=ios)
      get_linux_loop: DO
         READ(idef,'(a)', IOSTAT=ios) line
         IF(IS_IOSTAT_END(ios)) EXIT get_linux_loop
         j= LEN_TRIM(line)
         CALL rem_bl(line,j)
         j = INDEX(line,'=')
         CALL do_cap(line(1:j-1))
         IF(line(1:5) == 'NAME=') THEN
            operating_name=line(6:LEN_TRIM(line))
            j= LEN_TRIM(operating_name)
            DO i=1,j
               IF(operating_name(i:i)=='"') operating_name(i:i)=' '
            ENDDO
            CALL rem_bl(operating_name,j)
            EXIT get_linux_loop
         ENDIF
      ENDDO get_linux_loop
      CLOSE(idef)
   ELSEIF(INDEX(line,'CYGWIN') > 0) THEN
      operating   = OS_WINDOWS              ! Windows based Cygwin
      color_theme = THEME_LGHTYL_BLACK
!
!     Distinguish run within Cygwin or from DISCUS Icon at Windows
!
      INQUIRE(FILE='/Cygwin.bat',EXIST=lpresent)
      IF(lpresent) THEN
         OPEN(idef,FILE='/Cygwin.bat',STATUS='old')
         READ(idef,'(a)', IOSTAT=ios) line
         find_oper: DO WHILE(.NOT. IS_IOSTAT_END(ios))
            IF(line(1:5)=='chdir') THEN    ! FOUND correct line
              operat_top = line( 7:LEN_TRIM(LINE)-4)
              operating  = line(10:LEN_TRIM(LINE)-4) !'cygwin' or cygwin64
              EXIT find_oper
            ENDIF
            READ(idef,'(a)', IOSTAT=ios) line
         ENDDO find_oper
         CLOSE(idef)
         home_dir = ' ' 
         CALL get_environment_variable ('HOME', home_dir) 
         IF (home_dir.eq.' ') then 
            home_dir = '.' 
         ENDIF 
         home_dir_l = len_str (home_dir) 
         operating_name = 'CygwinNative'
      ELSE
         INQUIRE(FILE='/suite.sh',EXIST=lpresent)
         IF(lpresent) THEN
            operating = OS_WINDOWS  ! Cygwin from DISCUS Icon
            operating_name = 'CygwinIcon'
         ELSE
            operating_name = 'CygwinNative'
         ENDIF
         user_profile = ' '
         CALL get_environment_variable ('USERPROFILE', user_profile)
         home_dir = user_profile
         home_dir_l = len_str (home_dir) 
      ENDIF
   ELSEIF(INDEX(line,'DARWIN') > 0) THEN
      operating   = OS_MACOSX               ! MAC OS X
      color_theme = THEME_DEFAULT
      home_dir = ' ' 
      CALL get_environment_variable ('HOME', home_dir) 
      IF (home_dir.eq.' ') then 
         home_dir = '.' 
      ENDIF 
      home_dir_l = len_str (home_dir) 
      operating_name = 'MacOSX'
   ENDIF
ELSE   !  /proc/version does not exist , likely a MAC OS X 
!  Read OS from uname
   WRITE(ufile,'(a,a,i10.10)') tmp_dir(1:tmp_dir_l), '/DISCUS_SUITE_UNAME.' , PID
   WRITE(line,'(a,a)') 'uname -av > ', ufile(1:LEN_TRIM(ufile))
!   line = 'uname -av > '//tmp_dir(1:tmp_dir_l)//'/DISCUS_SUITE_UNAME'
   length = LEN_TRIM(line)
   CALL do_operating(line, length)
   line = ' '
!  WRITE(line,'(a,a,i10.10)')                                           &
!         tmp_dir(1:tmp_dir_l), '/DISCUS_SUITE_UNAME' , PID
!  INQUIRE(FILE=line                     ,EXIST=lpresent)
   OPEN(UNIT=idef, FILE=ufile)
   READ(IDEF,'(a)') line
   CLOSE(IDEF)
   CALL do_cap(line)
   IF(INDEX(line,'DARWIN') > 0) THEN
      operating   = OS_MACOSX               ! MAC OS X
      color_theme = THEME_DEFAULT
      home_dir = ' '
      CALL get_environment_variable ('HOME', home_dir)
      IF (home_dir.eq.' ') then
         home_dir = '.'
      ENDIF
      home_dir_l = len_str (home_dir) 
      operating_name = 'MacOSX'
   ENDIF
!  Remove temporary file
   WRITE(line,'(a,a)') 'rm -f ', ufile(1:LEN_TRIM(ufile))
   length = LEN_TRIM(line)
   CALL do_operating(line, length)
ENDIF
!
IF(operating==OS_WINDOWS) THEN
   IF(home_dir(home_dir_l:home_dir_l) == '/') THEN
      home_dir(home_dir_l:home_dir_l) = '\'
   ELSEIF(home_dir(home_dir_l:home_dir_l) /= '\') THEN
      home_dir(home_dir_l+1:home_dir_l+1) = '\'
      home_dir_l = home_dir_l + 1
   ENDIF
ELSE
   IF(home_dir(home_dir_l:home_dir_l) /= '/') THEN
      home_dir(home_dir_l+1:home_dir_l+1) = '/'
      home_dir_l = home_dir_l + 1
   ENDIF
ENDIF
! All remaining test for OS should be obsolete
   IF(operating == ' ') THEN
      CALL get_environment_variable ('OS', operating) 
      IF(operating == '') THEN
         CALL get_environment_variable ('OSTYPE', operating) 
      ENDIF
!
      color_theme = THEME_DEFAULT
      INQUIRE(FILE=  etc_os_release ,EXIST=lpresent)
      IF(lpresent) THEN
         CALL oeffne(idef,   etc_os_release , 'old')
         i=0
         name_search: DO
            READ(idef,'(a)',IOSTAT = ios) line 
            IF(IS_IOSTAT_END(ios)) EXIT name_search
            i=i+1
            IF(i.gt.100) EXIT name_search
            CALL do_cap(line)
            IF(INDEX(line,'SUSE') > 0) THEN
               operating   = 'Linux'
               color_theme = THEME_LGHTYL_BLACK
               EXIT name_search
            ELSEIF(INDEX(line,'UBUNTU') > 0) THEN
               operating   = 'Linux'
               color_theme = THEME_DEFAULT
               EXIT name_search
            ELSEIF(INDEX(line,'DEBIAN') > 0) THEN
               operating   = 'Linux'
               color_theme = THEME_DEFAULT
               EXIT name_search
            ELSEIF(INDEX(line,'CYGWIN') > 0) THEN
               operating   = 'Windows'
               color_theme = THEME_LGHTYL_BLACK
               EXIT name_search
            ELSEIF(INDEX(line,'DARWIN') > 0) THEN
               operating   = 'darwin18'
               color_theme = THEME_DEFAULT
               EXIT name_search
            ENDIF
         ENDDO name_search
         CLOSE(idef)
      ENDIF
   ENDIF
!
      IF(operating==' ') THEN    ! Still not found try MAC specifics
         CALL get_environment_variable ('DISPLAY', line)
         IF(INDEX(line, 'macos') > 0) THEN
            operating = 'darwin18'
         ENDIF
      ENDIF
      pname_l = len_str (pname) 
      line = ' ' 
      lines = 42 
      CALL get_environment_variable ('LINES', line) 
      IF (line.ne.' ') then 
         READ (line, *, end = 10) lines 
   10    CONTINUE 
      ELSE 
         CALL get_environment_variable ('TERMCAP', line) 
         ico = index (line, 'co') + 3 
         ice = index (line (ico:256) , ':') + ico - 2 
         IF (ice.gt.ico) then 
            READ (line (ico:ice), *, end = 20, err = 20) lines 
         ENDIF 
   20    CONTINUE 
      ENDIF 
      lines = lines - 2 
!
!     Try to find the correct share/ folder
!         Start from from process name
!         Else assume /usr/local/bin
!         Else assume $HOME/bin
IF(standalone) THEN
!
!              Process folder
!
   CALL GET_ENVIRONMENT_VARIABLE ('_', cdummy) 
   iii = INDEX(cdummy,pname,.true.)
   i   = INDEX(cdummy,'/',.TRUE.)
   j   = INDEX(cdummy,'\',.TRUE.)
   iii = MAX(i,j)
   appl_dir = ' '
   appl_dir(1:MAX(1,MIN(LEN(cdummy),iii)))   = cdummy(1:iii  )
   appl_dir_l = len_str (appl_dir) 
   i=INDEX(appl_dir(1:len_trim(appl_dir)-1),'/', .TRUE.)
   share_dir = appl_dir(1:i) // 'share/'
   hlpfile = share_dir(1:share_dir_l) // pname(1:pname_l)//'.hlp'
   hlpfile_l = LEN(TRIM (hlpfile) )
   CALL do_fexist(hlpfile,hlpfile_l,.FALSE.)
   IF(res_para(1)==0) THEN   ! hlpdir does not exist
!
!           Try /usr/local/bin
!
      share_dir = '/usr/local/share/'
      share_dir_l = 17
      hlpfile = share_dir(1:share_dir_l) // pname(1:pname_l)//'.hlp'
      hlpfile_l = LEN(TRIM (hlpfile) )
      CALL do_fexist(hlpfile,hlpfile_l,.FALSE.)
      IF(res_para(1)==0) THEN   ! hlpdir does not exist
!
!                 last resort, take home dir
!
         share_dir   = home_dir(1:home_dir_l) // 'share/'
         share_dir_l = home_dir_l + 6
      ENDIF
   ENDIF
ENDIF
!
!                                                                       
      deffile = '.'//pname (1:pname_l) 
      deffile_l = len_str (deffile) 
!                                                                       
      nullfile = '/dev/null' 
!
      CALL program_files
!
      IF(index(operating, 'Windows') /= 0) THEN  ! We got a Windows
         CALL do_cwd (start_dir, start_dir_l) 
         IF(start_dir == '/' .AND. start_dir_l==1) THEN
!           This is started from the icon, set start directory
!           to user home 
!           otherwise use current directory
         start_dir = ' '
         IF(user_profile == ' ') THEN
            start_dir   = 'C:\Users'
            start_dir_l = 8
         ELSE
            start_dir   = user_profile
            start_dir_l = len_str(start_dir)
         ENDIF
         IF(start_dir(start_dir_l:start_dir_l) /= '/') THEN
            start_dir   = start_dir(1:start_dir_l) // '/'
            start_dir_l = start_dir_l + 1
         ENDIF
         ELSE
!           start dir is not '/', started via system <pname>
!           Use current start directory
!
            IF(start_dir(start_dir_l:start_dir_l) /= '/') THEN
               start_dir   = start_dir(1:start_dir_l) // '/'
               start_dir_l = start_dir_l + 1
            ENDIF
         ENDIF
         CALL do_chdir ( start_dir, start_dir_l, .FALSE.)
      ELSEIF(operating==OS_LINUX_WSL) THEN
         CALL do_cwd (start_dir, start_dir_l) 
!
!           Started from Icon set start directory to
!           User HOME
            line = 'echo $PATH > '//tmp_dir(1:tmp_dir_l)//'/discus_suite_path.txt'
            length = LEN_TRIM(line)
            CALL do_operating(line, length)
            pathfile = tmp_dir(1:tmp_dir_l)//'/discus_suite_path.txt'
            OPEN(UNIT=idef,FILE=pathfile,ACTION='READ')
            READ(idef,'(a)') line
            CLOSE(UNIT=idef)
            i = INDEX(line,'/mnt/c/Users')
            IF(i/=0) THEN
               j = INDEX(line(i+13:),'/')
               IF(j/=0) THEN
                  user_profile = ' '
                  user_profile = line(i+13:i+13+j-2)
               ENDIF
            ENDIF
         IF(start_dir == '/mnt/c/Users' ) THEN
            user_name = user_profile(1:len(user_name))
                  IF(start_dir(start_dir_l:start_dir_l)=='/') THEN
                     start_dir_l = start_dir_l - 1
                  ENDIF
                  start_dir = start_dir(1:start_dir_l) // '/' // &
                              user_name(1:LEN_TRIM(user_name))
                  start_dir_l = LEN_TRIM(start_dir)
                  CALL do_chdir(start_dir, start_dir_l, .FALSE.) 
                  start_dir_l = LEN_TRIM(start_dir)
         ENDIF
         IF(start_dir(start_dir_l:start_dir_l) /= '/') THEN
            start_dir   = start_dir(1:start_dir_l) // '/'
            start_dir_l = start_dir_l + 1
         ENDIF
         get_screen_file = tmp_dir(1:tmp_dir_l)//'/discus_suite_screen.txt'
         get_screen_size = 'xdpyinfo | grep dimensions > '//get_screen_file(1:len_trim(get_screen_file))
      ELSE
!                                                                       
         CALL do_cwd (start_dir, start_dir_l) 
         IF(start_dir(start_dir_l:start_dir_l) /= '/') THEN
            start_dir   = start_dir(1:start_dir_l) // '/'
            start_dir_l = start_dir_l + 1
         ENDIF
         if(operating==OS_LINUX) then
           get_screen_file = tmp_dir(1:tmp_dir_l)//'/discus_suite_screen.txt'
           get_screen_size = 'xdpyinfo | grep dimensions > '//get_screen_file(1:len_trim(get_screen_file))
         elseif(operating==OS_MACOSX) then
           get_screen_file = tmp_dir(1:tmp_dir_l)//'/discus_suite_screen.txt'
           get_screen_size = 'xrandr -q | grep "\*" | awk ''{print " dimensions: " $1 " pixels"; }'' >' &
                             //get_screen_file(1:len_trim(get_screen_file))
         endif
      ENDIF
      current_dir   = start_dir
      current_dir_l = start_dir_l
!
!     Define terminal color scheme
!
!     CALL color_set_scheme (standalone, local_mpi_myid)
!
CALL lib_f90_getpname(PID,PID,process_name)
! Get Parent Process ID, required knowledge of operating
!
PPID = lib_f90_getppid(PID)
CALL lib_f90_getpname(PID,PPID, parent_name)
!
envir_done = .TRUE.
!write(*,*) ' HOME >',home_dir(1:home_dir_l),'<'
!
! test for and create "${HOME}/.DISCUS"
!
discus_dir = home_dir(1:home_dir_l) // '.DISCUS/'
discus_dir_l = home_dir_l + 8
!cpara(1) = discus_dir
!lpara(1) = discus_dir_l
!call sys_inquire_directory(MAXW, cpara, lpara,lexist)
!if(.not.lexist) then
!   line = 'mkdir -p ' // discus_dir(1:discus_dir_l)
!   call execute_command_line(line)
!!else
!!  write(*,*) ' DISCUS DIR :>', discus_dir(1:discus_dir_l),'< exists'
!endif
!
call generic_read_config(discus_dir, discus_dir_l)   ! Read generic preferences
!
screen_size = 300    !  A very rough default
call execute_command_line(get_screen_size(1:len_trim(get_screen_size)) , WAIT=.true., &
     cmdstat=ier_cmd, cmdmsg=message, exitstat=exit_msg)
open(unit=idef,file=get_screen_file,action='READ')
read(idef, '(a)', iostat=ios) line
close(idef)
if(ios==0) then
   i = index(line,':')
   j = index(line(i+1:),'x')
   read(line(i+1:i+j-1), *, iostat=ios) screen_size(1)
   iii=index(line(i+j+1:),' pix')
   read(line(i+j+1:i+j+1+iii-1), *, iostat=ios) screen_size(2)
endif
!                                                                       
END SUBROUTINE appl_env                       
!
!*******************************************************************************!
!
subroutine start_pgxwin
!-
! Start the pgxwin server 
!+
!
use precision_mod
use envir_mod
use errlist_mod
!
integer, parameter :: IDEF = 68
character(len=PREC_STRING) :: line
character(len=PREC_STRING) :: message
integer             :: exit_msg
integer             :: ier_cmd
integer :: ios
!
if(operating==OS_LINUX_WSL) THEN
   line = 'ps aux| fgrep -v grep | fgrep pgxwin_server > /tmp/DISCUS.PGXWIN'
   call execute_command_line(line(1:len_trim(line)) , WAIT=.true., &
        CMDSTAT=ier_cmd, CMDMSG=message, EXITSTAT=exit_msg)
   open(unit=idef, file='/tmp/DISCUS.PGXWIN', iostat=ios, status='old')
   line = ' '
   if(ios==0) then
      read(idef, '(a)', iostat=ios) line
   endif
   close(idef)
   if(index(line,'pgxwin_server')==0.or. is_iostat_end(ios)) then
      open(unit=idef, file='/tmp/DISCUS.PGXWIN', iostat=ios, status='unknown')
      write(idef,'(a)') '#!/bin/bash'
      write(idef,'(a)') 'pgxwin_server & '
      close(idef)
      line   = 'chmod ugo+rwx ' // '/tmp/DISCUS.PGXWIN'
      call execute_command_line(line(1:LEN_TRIM(line)), WAIT=.true.,  &
           CMDSTAT=ier_cmd, CMDMSG=message, EXITSTAT=exit_msg)
      line = '/tmp/DISCUS.PGXWIN &'
      call execute_command_line(line(1:len_trim(line)), WAIT=.false., &
           CMDSTAT=ier_cmd, CMDMSG=message, EXITSTAT=exit_msg)
   endif
endif
!                                                                       
END SUBROUTINE start_pgxwin                       
!
!*******************************************************************************!
!
SUBROUTINE write_appl_env (standalone, local_mpi_myid)
!-
!  Writes the path for help manual etc to the welcome screen
!+
USE envir_mod
USE errlist_mod
use gen_mpi_mod
USE lib_errlist_func
use lib_config
USE precision_mod
USE prompt_mod
USE support_mod
USE terminal_mod
!
IMPLICIT NONE
!                                                                       
LOGICAL, INTENT(IN) :: standalone
INTEGER, INTENT(IN) :: local_mpi_myid
INTEGER             :: old_version
INTEGER             :: new_version
INTEGER             :: since_update
CHARACTER(LEN=10)   :: cversion
CHARACTER(LEN=PREC_STRING) :: line
!
IF(local_mpi_myid/=0) RETURN
!
!  Analyse if new version is avalable at GIThub
!
if(.not.gen_mpi_active) then
   CALL lib_f90_test_updates(old_version, new_version, cversion, since_update)
endif
old_version = 0
new_version = 0
!
!
IF(standalone .AND. local_mpi_myid==0) THEN
   IF(term_scheme_exists) THEN
!     WRITE ( *, 1900) TRIM(color_bg),TRIM(color_info), man_dir (1:LEN_TRIM(man_dir)) ,TRIM(color_fg)
!     WRITE ( *, 2000) TRIM(color_bg),TRIM(color_info),umac_dir (1:LEN_TRIM(umac_dir)),TRIM(color_fg)
!     WRITE ( *, 2100) TRIM(color_bg),TRIM(color_info),mac_dir (1:mac_dir_l),     TRIM(color_fg)
      line =  start_dir
      CALL display_dir(line)
      WRITE ( *, 2200) TRIM(color_bg),TRIM(color_info),line (1:LEN_TRIM(line)) ,TRIM(color_fg)
      WRITE ( *, 2300) TRIM(color_bg),TRIM(color_info),TRIM(color_fg)
      WRITE ( *, 2400) TRIM(color_bg),TRIM(color_info),TRIM(color_fg)
      WRITE ( *, 2500) TRIM(color_bg),TRIM(color_info),TRIM(color_fg)
      WRITE ( *, 2600) TRIM(color_bg),TRIM(color_info),TRIM(color_fg)
      IF(operating == OS_LINUX_WSL) THEN
         WRITE ( *, 2700) TRIM(color_bg),TRIM(color_info),TRIM(color_fg)
         if(generic_get_interval()>0) then
         IF(since_update>generic_get_interval()) THEN
            WRITE(*,*)
            WRITE(*,2750) TRIM(color_bg),TRIM(color_fg)
         ENDIF
         endif
      ENDIF
      IF(new_version > old_version ) THEN
         WRITE(*,*)
         WRITE(*,2800) TRIM(color_bg),TRIM(color_info),cversion, TRIM(color_fg)
      ENDIF
   ELSE
!     WRITE ( *,  900)  man_dir (1:LEN_TRIM(man_dir)) 
!     WRITE ( *, 1000) umac_dir (1:LEN_TRIM(umac_dir)) 
!     WRITE ( *, 1100) mac_dir (1:mac_dir_l) 
      line =  start_dir
      CALL display_dir(line)
      WRITE ( *, 1200) line (1:LEN_TRIM(line))
      WRITE ( *, 1300)
      WRITE ( *, 1400)
      WRITE ( *, 1500)
      WRITE ( *, 1600)
      IF(operating == OS_LINUX_WSL) THEN
         WRITE ( *, 1700)
         if(generic_get_interval()>0) then
         IF(since_update>generic_get_interval()) THEN
            WRITE(*,*)
            WRITE(*,1750) 
         ENDIF
         endif
      ENDIF
      IF(new_version > old_version ) THEN
         WRITE(*,*)
         WRITE(*,1800) cversion
      ENDIF
   ENDIF
ENDIF
CALL lib_f90_findterminal
!                                                                       
  900 FORMAT     (1x,'Manual files in  : ',a) 
 1000 FORMAT     (1x,'User macros in   : ',a) 
 1100 FORMAT     (1x,'System macros in : ',a) 
 1200 FORMAT     (1x,'Start directory  : ',a) 
 1300 FORMAT     (1x,'Access manuals at each section with   : ','manual')
 1400 FORMAT     (1x,'Access help at each section/menu with : ','help  ')
 1500 FORMAT     (1x,'News at each section/Command_lang in  : ','help News')
 1600 FORMAT     (1x,'Change font size with                 : ','CTRL + /CTRL -')
 1700 FORMAT     (1x,'Preferences: Right click + Preferences: ' )
 1750 FORMAT     (1x,'Recommend to update Ubuntu at DISCUS exit : ')
 1800 FORMAT     (1x,'New DISCUS version available at GIThub: ',a)
 1900 FORMAT     (1x,a,'Manual files in  : ',a,a,a) 
 2000 FORMAT     (1x,a,'User macros in   : ',a,a,a) 
 2100 FORMAT     (1x,a,'System macros in : ',a,a,a) 
 2200 FORMAT     (1x,a,'Start directory  : ',a,a,a) 
 2300 FORMAT     (1x,a,'Access manuals at each section with   : ',a,'manual',a)
 2400 FORMAT     (1x,a,'Access help at each section/menu with : ',a,'help  ',a)
 2500 FORMAT     (1x,a,'News at each section/Command_lang in  : ',a,'help News',a)
 2600 FORMAT     (1x,a,'Change font size with                 : ',a,'CTRL + /CTRL -',a)
 2700 FORMAT     (1x,a,'Preferences:',a,' Right click + Preferences',a,':')
 2750 FORMAT     (1x,a,'Recommend to update Ubuntu at DISCUS exit : ',a)
 2800 FORMAT     (1x,a,'New DISCUS version available at GIThub: ',a,a,a)
!
END SUBROUTINE write_appl_env                       
!
!
!*******************************************************************************!
!
SUBROUTINE color_set_scheme (standalone, local_mpi_myid)
!-
!  Set the terminal color scheme 
!  If the file share/discus.term.scheme exists it is used
!  Else we try a default color scheme according to the
!  Operating system
!+
!
USE blanks_mod
USE terminal_mod
USE envir_mod
USE param_mod
USE precision_mod
USE prompt_mod
USE string_convert_mod
USE support_mod
!
IMPLICIT NONE
!
LOGICAL, INTENT(IN) :: standalone
INTEGER, INTENT(IN) :: local_mpi_myid
!
INTEGER, PARAMETER :: idef = 68
CHARACTER(LEN=PREC_STRING) :: line, color
INTEGER :: ios ! I/O status
INTEGER :: i, icolon, length
!
term_scheme_file = appl_dir(1:appl_dir_l)//'../share/discus.term.scheme'
term_scheme_file_l = LEN_TRIM(term_scheme_file)
CALL do_fexist(term_scheme_file,term_scheme_file_l,.FALSE.)
IF(res_para(1)==0) THEN   ! discus.term.scheme does not exist
   term_scheme_exists = .FALSE.
   IF(color_theme == THEME_DEFAULT) THEN  ! Use default foreground and background
      color_fg   = COLOR_FG_DEFAULT
      color_bg   = COLOR_BG_DEFAULT
      color_err  = COLOR_FG_RED   ! Lets hope the background is NOT red...
      color_info = COLOR_FG_BLUE  ! Lets hope the background is NOT blue...
   ELSEIF(color_theme == THEME_LGHTYL_BLACK) THEN  ! Use default foreground and background
      color_fg   = COLOR_FG_BLACK
      color_bg   = COLOR_BG_LIGHT_YELLOW
      color_err  = COLOR_FG_RED   ! Lets hope the background is NOT red...
      color_info = COLOR_FG_BLUE  ! Lets hope the background is NOT blue...
   ELSEIF(color_theme == THEME_BLACK_WHITE ) THEN  ! Use default foreground and background
      color_fg   = COLOR_FG_BLACK
      color_bg   = COLOR_BG_WHITE
      color_err  = COLOR_FG_RED   ! Lets hope the background is NOT red...
      color_info = COLOR_FG_BLUE  ! Lets hope the background is NOT blue...
   ENDIF
ELSE
   color_fg   = COLOR_FG_DEFAULT  ! Set a sensible default anyway
   color_bg   = COLOR_BG_DEFAULT  ! Set a sensible default anyway
   color_err  = COLOR_FG_RED      ! Lets hope the background is NOT red...
   color_info = COLOR_FG_BLUE     ! Lets hope the background is NOT blue...
   color_theme = THEME_USER
   term_scheme_exists = .TRUE.
   CALL oeffne(idef, term_scheme_file, 'old')
   i=0
   color_search: DO
      READ(idef,'(a)',IOSTAT = ios) line
      IF(IS_IOSTAT_END(ios)) EXIT color_search
      i=i+1
      IF(i.gt.100) EXIT color_search
      CALL do_cap(line)
      icolon = INDEX(line, ':')
      color  = line(icolon+1:LEN_TRIM(line))
      length = LEN_TRIM(line)
      CALL rem_leading_bl(color,length)
      CALL color_set_fg(line,icolon,'FOREGROUND:', color, color_fg)
      CALL color_set_fg(line,icolon,'ERROR:',      color, color_err)
      CALL color_set_fg(line,icolon,'HIGHLIGHT:',  color, color_high)
      CALL color_set_fg(line,icolon,'INFORMATION:',color, color_info)
      CALL color_set_bg(line,icolon,'BACKGROUND:', color, color_bg)
   END DO color_search
   CLOSE(idef)
ENDIF
IF(standalone.AND.local_mpi_myid==0) THEN  ! set standard Background and Foreground
   WRITE(*,'(a)') TRIM(color_bg)
   WRITE(*,'(a)') TRIM(color_fg)
!
DO i=1, LEN(line)
  line(i:i) = '.'
ENDDO
WRITE(*,*) line(1:MIN(256,LEN_TRIM(line)))
WRITE(*,*) 
WRITE(*,*) 
WRITE(*,*) 
WRITE(*,*) 
WRITE(*,*) 
ENDIF
!
END SUBROUTINE color_set_scheme
!
!
!*******************************************************************************!
!
SUBROUTINE color_set_fg(line,icolon,color_type , color_string, color_color)
!-
!  Sets colors for foreground
!+
!
USE ber_params_mod
USE errlist_mod
USE precision_mod
USE terminal_mod
IMPLICIT NONE
!
CHARACTER (LEN=*) , INTENT(IN)  :: line          ! == 'FORGROUND:' or similar
INTEGER           , INTENT(IN)  :: icolon        ! Position of colon
CHARACTER (LEN=*) , INTENT(IN)  :: color_type    ! String for comparison
CHARACTER (LEN=*) , INTENT(IN)  :: color_string  ! User color name/value
CHARACTER (LEN=*) , INTENT(OUT) :: color_color   ! Intended output color string
!
CHARACTER (LEN=   3) :: zeile
INTEGER              :: ianz
INTEGER, PARAMETER                 :: MAXW = 1
CHARACTER (LEN=PREC_STRING), DIMENSION(1) :: cpara
INTEGER             , DIMENSION(1) :: lpara
REAL(KIND=PREC_DP)  , DIMENSION(1) :: werte
!
IF(line(1:icolon) == color_type) THEN
   IF(color_string=='BLACK') THEN
      color_color = COLOR_FG_BLACK
   ELSEIF(color_string=='RED') THEN
      color_color = COLOR_FG_RED
   ELSEIF(color_string=='GREEN') THEN
      color_color = COLOR_FG_GREEN
   ELSEIF(color_string=='YELLOW') THEN
      color_color = COLOR_FG_YELLOW
   ELSEIF(color_string=='BLUE') THEN
      color_color = COLOR_FG_BLUE
   ELSEIF(color_string=='MAGENTA') THEN
      color_color = COLOR_FG_MAGENTA
   ELSEIF(color_string=='CYAN') THEN
      color_color = COLOR_FG_CYAN
   ELSEIF(color_string=='WHITE') THEN
      color_color = COLOR_FG_WHITE
   ELSE
      cpara(1) = color_string
      lpara(1) = LEN_TRIM(color_string)
      ianz   = 1
      CALL ber_params ( ianz, cpara, lpara, werte, MAXW)
      IF(ier_num == 0 .AND. 0<=NINT(werte(1)) .AND. NINT(WERTE(1))<=255) THEN
         WRITE(zeile,'(I3.3)') NINT(werte(1))
         color_color = ESC//'[38;5;'//zeile(1:3)//'m'
      ELSE
         color_color = COLOR_FG_DEFAULT
      ENDIF
   ENDIF
ENDIF
END SUBROUTINE color_set_fg
!
!
!*******************************************************************************!
!
SUBROUTINE color_set_bg(line,icolon,color_type , color_string, color_color)
!-
!  Sets colors for background
!+
!
USE ber_params_mod
USE errlist_mod
USE precision_mod
USE terminal_mod
IMPLICIT NONE
!
CHARACTER (LEN=*) , INTENT(IN)  :: line          ! == 'FORGROUND:' or similar
INTEGER           , INTENT(IN)  :: icolon        ! Position of colon
CHARACTER (LEN=*) , INTENT(IN)  :: color_type    ! String for comparison
CHARACTER (LEN=*) , INTENT(IN)  :: color_string  ! User color name/value
CHARACTER (LEN=*) , INTENT(OUT) :: color_color   ! Intended output color string
!
CHARACTER (LEN=   3) :: zeile
INTEGER              :: ianz
INTEGER, PARAMETER                 :: MAXW = 1
CHARACTER (LEN=PREC_STRING), DIMENSION(1) :: cpara
INTEGER             , DIMENSION(1) :: lpara
REAL(KIND=PREC_DP)  , DIMENSION(1) :: werte
!
IF(line(1:icolon) == color_type) THEN
   IF(color_string=='BLACK') THEN
      color_color = COLOR_FG_BLACK
   ELSEIF(color_string=='RED') THEN
      color_color = COLOR_BG_RED
   ELSEIF(color_string=='GREEN') THEN
      color_color = COLOR_BG_LIGHT_YELLOW
   ELSEIF(color_string=='LIGHT_YELLOW') THEN
      color_color = COLOR_BG_LIGHT_YELLOW
   ELSEIF(color_string=='YELLOW') THEN
      color_color = COLOR_BG_YELLOW
   ELSEIF(color_string=='BLUE') THEN
      color_color = COLOR_BG_BLUE
   ELSEIF(color_string=='WHITE') THEN
      color_color = COLOR_BG_WHITE
   ELSE
      cpara(1) = color_string
      lpara(1) = LEN_TRIM(color_string)
      ianz   = 1
      CALL ber_params ( ianz, cpara, lpara, werte, MAXW)
      IF(ier_num == 0 .AND. 0<=NINT(werte(1)) .AND. NINT(WERTE(1))<=255) THEN
         WRITE(zeile,'(I3.3)') NINT(werte(1))
         color_color = ESC//'[48;5;'//zeile(1:3)//'m'
      ELSE
         color_color = COLOR_BG_DEFAULT
      ENDIF
   ENDIF
ENDIF
!
END SUBROUTINE color_set_bg
!
!*******************************************************************************!
!
SUBROUTINE  program_files
!-                                                                      
!     Sets path for helpfile , mac directories
!     UNIX version ..                                                   
!+                                                                      
USE envir_mod 
USE errlist_mod
USE precision_mod
USE prompt_mod 
USE support_mod
!
IMPLICIT none 
!                                                                       
!INTEGER, PARAMETER :: ird = 34
!
!CHARACTER(LEN=PREC_STRING) :: line
CHARACTER(LEN=7)    :: progname
INTEGER :: progname_l
INTEGER :: i
!LOGICAL :: l_exist
!
IF(lstandalone) THEN
   progname = pname          ! Make a local copy of the program name pname
ELSE
   progname = 'suite'
ENDIF
!
progname_l = LEN(TRIM(progname))
!
!
mac_dir = ' ' 
mac_dir = share_dir // pname (1:LEN_TRIM(pname)) //'/'
mac_dir_l = LEN_TRIM (mac_dir) 
!                                                                       
umac_dir = home_dir(1:home_dir_l)//'mac/'//progname(1:progname_l) //'/'
umac_dir_l = LEN(TRIM (umac_dir) )
!                                                                       
hlpfile = ' ' 
hlpdir  = ' ' 
hlpdir  = share_dir(1:MIN(LEN(hlpdir),LEN_TRIM(share_dir)))
hlp_dir_l = share_dir_l
hlpfile   = hlpdir(1:hlp_dir_l)//progname(1:progname_l)//'.hlp'
hlpfile_l = LEN(TRIM (hlpfile) )
!                                                                       
colorfile = ' ' 
colorfile = share_dir(1:share_dir_l) // 'color.map'
colorfile_l = LEN(TRIM (colorfile) )
!
!   Search for the Installation file "DiscusSuite.txt" to obtain paths
!
man_dir = hlpdir
!
i=LEN_TRIM(man_dir)
IF(operating(1:7)=='Windows') THEN
   IF(man_dir(i:i) /='\') THEN
      man_dir(i+1:i+1) = '\'
   ENDIF
ELSE
   IF(man_dir(i:i) /='/') THEN
      man_dir(i+1:i+1) = '/'
   ENDIF
ENDIF
!
END SUBROUTINE  program_files
!
!*******************************************************************************
!
INTEGER FUNCTION lib_f90_getppid(cpid)
!-
! Determine Parent Process ID for current PID cpid
!+
!
USE envir_mod
USE errlist_mod
USE precision_mod
USE support_mod
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: cpid  ! current PID whose Parent Process ID is to be found
!
INTEGER, PARAMETER  :: ITMP     = 79  ! temporary unit number
!
CHARACTER(LEN=PREC_STRING) :: line
CHARACTER(LEN=PREC_STRING) :: temp_file
character(len=PREC_STRING) :: message
integer             :: exit_msg
integer             :: ier_cmd
INTEGER             :: tpid           ! temporary pid for current process for verification
INTEGER             :: tppid          ! temporary ppid for current process for verification
INTEGER             :: ios            ! I/O status 
!
tppid = 0   ! Default value 
!
!  Make a temp_file to write the system(ps) into
WRITE(temp_file, '(a,I10.10)') '/tmp/getppid.', cpid
!
IF(operating(1:6)=='darwin') THEN
   WRITE(line,'(a,i10,a,a)') 'ps j | grep ',PID,' | grep -v grep | awk ''{print $2, $3}'' >> ', &
       temp_file(1:LEN_TRIM(temp_file))
   call execute_command_line(line(1:len_trim(line)) , WAIT=.true., &
        cmdstat=ier_cmd, cmdmsg=message, exitstat=exit_msg)
ELSEIF(operating==OS_WINDOWS .OR. operating==OS_CYGWIN64 .OR. operating==OS_CYGWIN32) THEN
   WRITE(line,'(a,i10,a,a)') 'ps j | grep ',PID,' | grep discus  | awk ''{print $1, $2}'' >> ', &
       temp_file(1:LEN_TRIM(temp_file))
   call execute_command_line(line(1:len_trim(line)) , WAIT=.true., &
        CMDSTAT=ier_cmd, CMDMSG=message, EXITSTAT=exit_msg)
ELSE
   WRITE(line,'(a,i10,a,a)') 'ps j | grep ',PID,' | grep -v grep | awk ''{print $2, $1}'' >> ', &
       temp_file(1:LEN_TRIM(temp_file))
   call execute_command_line(line(1:len_trim(line)) , WAIT=.true., &
        CMDSTAT=ier_cmd, CMDMSG=message, EXITSTAT=exit_msg)
ENDIF
!
CALL oeffne( ITMP, temp_file, 'old')
IF(ier_num==0) THEN
   READ(ITMP,*,IOSTAT=ios) tpid, tppid
   loop: DO WHILE (.NOT.IS_IOSTAT_END(ios))
      IF(tpid==cpid) EXIT LOOP        ! got the correct line with pid and then tppid
      READ(ITMP,*,IOSTAT=ios) tpid, tppid
   ENDDO loop
ENDIF
CLOSE(ITMP)
line = 'rm -f ' // temp_file(1:len_trim(temp_file))
call execute_command_line(line(1:len_trim(line)) , WAIT=.true., &
     CMDSTAT=ier_cmd, CMDMSG=message, EXITSTAT=exit_msg)
!
lib_f90_getppid = tppid
!
END FUNCTION lib_f90_getppid
!
!*******************************************************************************
!
SUBROUTINE lib_f90_getpname(cpid, tpid,tpname)
!-
! Get the process name of the process with PID tpid
!+
USE envir_mod
USE precision_mod
!
IMPLICIT NONE
!
INTEGER         , INTENT(IN)  :: cpid  ! Current Process id
INTEGER         , INTENT(IN)  :: tpid  ! Process id
CHARACTER(LEN=*), INTENT(OUT) :: tpname  ! process name
!
INTEGER, PARAMETER  :: ITMP     = 79  ! temporary unit number
!
CHARACTER(LEN=PREC_STRING) :: line           ! Dummy line
CHARACTER(LEN=PREC_STRING) :: temp_file      ! temporary file
character(len=PREC_STRING) :: message
integer             :: exit_msg
integer             :: ier_cmd
INTEGER             :: ios            ! I/O status 
INTEGER             :: islash         ! I/O status 
!
IF(tpid == 0) THEN
  tpname = 'NULL'
  RETURN
ENDIF
!  Make a temp_file to write the system(ps) into
WRITE(temp_file, '(a,I10.10)') '/tmp/getpname.', cpid
!
IF(operating(1:6)=='darwin') THEN
   WRITE(line,'(a,I10,a,a)') 'ps -p ',tpid, ' -o comm= > ',temp_file(1:LEN_TRIM(temp_file))
ELSEIF(operating==OS_WINDOWS .OR. operating==OS_CYGWIN64 .OR. operating==OS_CYGWIN32) THEN
   WRITE(line,'(a,i10,a,i10,a,a)') 'ps -p ',tpid, ' | grep ', tpid, '| grep -v grep > ', &
   temp_file(1:LEN_TRIM(temp_file))
ELSEIF(operating==OS_LINUX) THEN
   WRITE(line,'(a,I10,a,a)') 'ps -p ',tpid, ' -o comm= > ',temp_file(1:LEN_TRIM(temp_file))
ELSEIF(operating==OS_LINUX_WSL) THEN
   WRITE(line,'(a,I10,a,a)') 'ps -p ',tpid, ' -o comm= > ',temp_file(1:LEN_TRIM(temp_file))
ENDIF
call execute_command_line(line(1:len_trim(line)) , WAIT=.true., &
     cmdstat=ier_cmd, cmdmsg=message, exitstat=exit_msg)
!
tpname = ' '
OPEN(UNIT=ITMP,FILE=temp_file,STATUS='old')
IF(operating==OS_WINDOWS .OR. operating==OS_CYGWIN64 .OR. operating==OS_CYGWIN32) THEN
   READ(ITMP,'(a)', IOSTAT=ios) line
   READ(LINE(65:),'(a)', IOSTAT=ios) tpname
   islash =INDEX(tpname,'/', .TRUE.)
   READ(tpname(islash+1:),'(a)', IOSTAT=ios) tpname
ELSE
   READ(ITMP,'(a)', IOSTAT=ios) tpname
   IF(IS_IOSTAT_END(ios)) THEN
      tpname = 'NULL'
   ENDIF
ENDIF
CLOSE(ITMP)
line = 'rm -f ' // temp_file(1:len_trim(temp_file))
call execute_command_line(line(1:len_trim(line)) , WAIT=.true., &
     cmdstat=ier_cmd, cmdmsg=message, exitstat=exit_msg)

END SUBROUTINE lib_f90_getpname
!
!*******************************************************************************
!
SUBROUTINE lib_f90_init_updates
!
USE envir_mod
USE errlist_mod
USE lib_errlist_func
use gen_mpi_mod
USE precision_mod
!
IMPLICIT NONE
!
integer, parameter :: IWR = 68
CHARACTER(LEN=128        ) :: cfile
character(len=PREC_STRING) :: ufile
CHARACTER(LEN=PREC_STRING) :: line
CHARACTER(LEN=PREC_STRING) :: message
INTEGER             :: exit_msg
INTEGER             :: ier_cmd
LOGICAL             :: lda
!
if(gen_mpi_active) return               ! No update checks while mpi is active
if(.not. check_github() ) return    ! no www network access 
WRITE(cfile,'(a,a,i10.10)') tmp_dir(1:len_trim(tmp_dir)),'/DISCUS_CURRENT.', PID ! Initiate search for new version
!
INQUIRE(FILE=cfile, EXIST=lda)
IF(lda) THEN 
   line = "rm -f /" // cfile(1:LEN_TRIM(cfile))
   CALL execute_command_line (line(1:LEN_TRIM(line)), wait=.true.,  &
        CMDSTAT=ier_cmd, CMDMSG=message, EXITSTAT=exit_msg)
ENDIF
! 140.82.121.3 = github.com
line = "curl -k --silent --location --connect-timeout 3 --max-time 3 ""https://github.com/tproffen/DiffuseCode/releases/latest"" > " &
       // cfile(1:LEN_TRIM(cfile))
write(ufile,'(a,a,i10.10)') tmp_dir(1:len_trim(tmp_dir)), '/DISCUS.UFILE.', PID
open(unit=IWR, file=ufile, status='unknown')
write(IWR, '(a)') '#!/bin/bash'
write(IWR, '(a)') line(1:len_trim(line))
close(IWR)
line   = 'chmod ugo+rwx ' // ufile(1:len_trim(ufile))
CALL execute_command_line (line(1:LEN_TRIM(line)), wait=.true.,  &
     CMDSTAT=ier_cmd, CMDMSG=message, EXITSTAT=exit_msg)
line = ufile(1:len_trim(ufile)) // ' & '
CALL execute_command_line (line(1:LEN_TRIM(line)), WAIT=.false., &
     CMDSTAT=ier_cmd, CMDMSG=message, EXITSTAT=exit_msg)
!call execute_command_line ('cat /tmp/DISCUS_CURRENT*')
CALL no_error                         ! Checked and printed in ==> write_appl_env
END SUBROUTINE lib_f90_init_updates
!
!*******************************************************************************
!
SUBROUTINE lib_f90_test_updates(old_version, new_version, cversion, since_update)
!
use charact_mod
USE envir_mod
USE errlist_mod
USE lib_errlist_func
use gen_mpi_mod
USE precision_mod
USE prompt_mod
!
IMPLICIT NONE
!
INTEGER          , INTENT(OUT) :: old_version
INTEGER          , INTENT(OUT) :: new_version
CHARACTER(LEN=10), INTENT(OUT) :: cversion
INTEGER          , INTENT(OUT) :: since_update
!
INTEGER, PARAMETER  :: IRD = 69
!
CHARACTER(LEN=128        ) :: cfile
CHARACTER(LEN=128        ) :: discus_update
CHARACTER(LEN=PREC_STRING) :: string
CHARACTER(LEN=PREC_STRING) :: message
INTEGER             :: exit_msg
INTEGER             :: ier_cmd
INTEGER             :: idot1, idot2
INTEGER             :: iquote
INTEGER             :: itag
INTEGER             :: icode
INTEGER             :: idate
INTEGER             :: i, j, k
INTEGER             :: ios
!
LOGICAL :: lda
logical :: lonline
CHARACTER(LEN=8) :: date
CHARACTER(LEN=10):: time
CHARACTER(LEN=5) :: zone
INTEGER, DIMENSION(8) :: values
!integer, dimension(12) :: days
integer :: iyear
!
!data days/ 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334/
!                                                                       
if(gen_mpi_active) return               ! No update checks while mpi is active
CALL DATE_AND_TIME (date, time, zone, values)
!
!  Analyse if new version is available at GIThub
!
WRITE(cfile,        '(a,a,i10.10)') tmp_dir(1:len_trim(tmp_dir)),'/DISCUS_CURRENT.', PID ! Initiate search for new version
WRITE(discus_update,'(a,a       )') tmp_dir(1:len_trim(tmp_dir)),'/DISCUS_UPDATE'        ! Initiate search for new version
old_version = 0
new_version = 0
!
lonline = .FALSE.                 ! Assume offline
INQUIRE(FILE=cfile, EXIST=lda)
!
ON_LINE: IF(lda) THEN
   OPEN(UNIT=IRD, FILE=cfile, STATUS='OLD')
   LOOP_SEARCH: do
      READ(IRD,'(a)', iostat=ios) string
      if(is_iostat_end(ios)) exit ON_LINE
      lonline = .TRUE.                 ! We are online
      itag = INDEX(string,'tag')
      icode = index(string,'DiffuseCode')
      IF(itag>0 .and. icode>0) THEN
         iquote = itag + INDEX(string(itag:LEN_TRIM(string)),'"') - 1
         cversion = string (itag+6:iquote-1)
         idot1 = INDEX(cversion, '.')
         idot2 = INDEX(cversion, '.', .TRUE.)
         READ(cversion(1:idot1-1), *) i
         READ(cversion(idot1+1:idot2-1), *) j
         READ(cversion(idot2+1:LEN_TRIM(cversion)), *) k
         new_version = i*10000 + j*100 + k
         exit ON_LINE
      ENDIF
   enddo LOOP_SEARCH
ENDIF ON_LINE
CLOSE(UNIT=IRD)
!
do i=1, len_trim(version)
   if(iachar(version(i:i))<0 .or. iachar(version(i:i))>nine) version(i:i) = ' ' ! Eliminate characters
enddo
!
idot1 = INDEX(version, '.')
idot2 = INDEX(version, '.', .TRUE.)
READ(version(1:idot1-1), *) i
READ(version(idot1+1:idot2-1), *) j
READ(version(idot2+1:LEN_TRIM(version)), *) k
old_version = i*10000 + j*100 + k
!
string = "rm -f " // cfile(1:LEN_TRIM(cfile))
CALL execute_command_line (string(1:LEN_TRIM(string)), CMDSTAT=ier_cmd, CMDMSG=message, EXITSTAT=exit_msg)
CALL no_error
if(.not. lonline) return
!
! Test for updates of operating system
!
last_update = 0
since_update = 0
!
cond_update: IF(operating==OS_LINUX_WSL) THEN
   string = 'stat /var/cache/apt/pkgcache.bin > '// discus_update(1:len_trim(discus_update))
   call execute_command_line(string, wait=.true., &
        cmdstat=ier_cmd, cmdmsg=message, exitstat=exit_msg)
   call no_error
!   call execute_command_line('cat /tmp/DISCUS_UPDATE')
   INQUIRE(FILE=discus_update,EXIST=lda)
   IF(lda) THEN
      OPEN(UNIT=IRD, FILE=discus_update, STATUS='old')
      do i=1, 6
         READ(IRD,'(a)', IOSTAT=ios) string
         if(ios/=0) exit cond_update
      enddo
      iyear=index(string,'20')
      READ(string(iyear:iyear+9),'(i4,1x,i2,1x,i2)') i,j,k
      last_update = days_since(i,j,k)
      CLOSE(IRD)
   ENDIF
   READ(date,'(i4,i2,i2)') i,j,k
   idate = days_since(i,j,k)
   since_update = idate-last_update
ENDIF cond_update
!
END SUBROUTINE lib_f90_test_updates
!
!*******************************************************************************
!
SUBROUTINE lib_f90_update_ubuntu
!
USE envir_mod
USE errlist_mod
USE lib_errlist_func
use gen_mpi_mod
USE precision_mod
USE prompt_mod
!
implicit none
!
INTEGER, PARAMETER  :: IRD = 69
!
CHARACTER(LEN=PREC_STRING) :: string
CHARACTER(LEN=PREC_STRING) :: cfile
CHARACTER(LEN=PREC_STRING) :: message
INTEGER             :: exit_msg
integer             :: ier_cmd
integer             :: ios
logical             :: lonline
logical             :: lda
!
if(gen_mpi_active) return               ! No update checks while mpi is active
if(.not. check_github() ) return    ! no www network access 
lonline = .TRUE.
WRITE(cfile,        '(a,a,i10.10)') tmp_dir(1:len_trim(tmp_dir)),'/DISCUS_UBUNTU.', PID ! Test file to see if we are on-line
string = "curl -k --silent --location ""https://github.com/tproffen/DiffuseCode/releases/latest"" > " &
       // cfile(1:LEN_TRIM(cfile))                                                      ! The actual command
call execute_command_line(string, wait=.true., &
        CMDSTAT=ier_cmd, CMDMSG=message, EXITSTAT=exit_msg)
INQUIRE(FILE=cfile,EXIST=lda)
IF_LDA: IF(lda) THEN
   open(unit=IRD, file=cfile, status='OLD')
   string = ' '
   read(IRD,'(a)', iostat=ios) string
   if(is_iostat_end(ios)) lonline = .FALSE.    ! Apparently offline
!  if(string(1:88) /= '<html><body>You are being <a href="https://github.com/tproffen/DiffuseCode/releases/tag/') &
!     lonline = .FALSE.                        ! Apparently offline
   on_line: if(lonline) then
!
   IF(operating == OS_LINUX .OR. operating == OS_LINUX_WSL) THEN
      WRITE(output_io, '(a)') 
      WRITE(output_io, '(a)')  ' Your Ubuntu system will be updated '
      string = 'sudo apt update'
      CALL execute_command_line (string(1:LEN_TRIM(string)), wait=.true., &
           CMDSTAT=ier_cmd, CMDMSG=message, EXITSTAT=exit_msg)
      string = 'sudo apt upgrade -y'
      CALL execute_command_line (string(1:LEN_TRIM(string)), wait=.true., &
           CMDSTAT=ier_cmd, CMDMSG=message, EXITSTAT=exit_msg)
      string = 'date --iso-8601 > '// tmp_dir(1:LEN_TRIM(tmp_dir)) // '/DISCUS_UPDATE'
      CALL execute_command_line (string(1:LEN_TRIM(string)), wait=.true., &
           CMDSTAT=ier_cmd, CMDMSG=message, EXITSTAT=exit_msg)
      string = 'chmod ugo+rw ' // tmp_dir(1:len_trim(tmp_dir)) // '/DISCUS_UPDATE'
      CALL execute_command_line (string(1:LEN_TRIM(string)), wait=.true., &
           CMDSTAT=ier_cmd, CMDMSG=message, EXITSTAT=exit_msg)
      WRITE(output_io, '(a)') 
   ENDIF
   endif on_line
!
   string = "rm -f " // cfile(1:LEN_TRIM(cfile))
!
   CALL execute_command_line (string(1:LEN_TRIM(string)), wait=.true., &
        CMDSTAT=ier_cmd, CMDMSG=message, EXITSTAT=exit_msg)
   CALL no_error
endif if_lda             ! test file is present
!
END SUBROUTINE lib_f90_update_ubuntu
!
!*******************************************************************************
!
SUBROUTINE lib_f90_update_discus(zeile, lp)
!
USE envir_mod
USE errlist_mod
use exit_para_mod
USE get_params_mod
use gen_mpi_mod
USE precision_mod
USE prompt_mod
USE string_convert_mod
USE support_mod
USE take_param_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: zeile
INTEGER         , INTENT(INOUT) :: lp
!
INTEGER, PARAMETER :: IRD = 69
INTEGER, PARAMETER :: IWR = 68
INTEGER, PARAMETER :: MAXW = 2
CHARACTER(LEN=PREC_STRING) :: string
CHARACTER(LEN=PREC_STRING) :: line
CHARACTER(LEN=PREC_STRING) :: cdir              ! previous current directory
CHARACTER(LEN=PREC_STRING) :: code_str          ! optional "code=" string
CHARACTER(LEN=PREC_STRING) :: inst_str          ! optional "install=" string
CHARACTER(LEN=PREC_STRING) :: prep_str          ! optional "prepare=" string
CHARACTER(LEN=128)          :: discus_version
CHARACTER(LEN=128)         :: discus_power
CHARACTER(LEN=32)          :: grep
CHARACTER(LEN=40)          :: script
CHARACTER(LEN=20)          :: verstring
CHARACTER(LEN= 8)          :: flag
CHARACTER(LEN=PREC_STRING) :: command
CHARACTER(LEN=PREC_STRING) :: message
INTEGER             :: exit_msg
integer             :: ier_cmd
integer             :: ios                      ! I/O status
INTEGER             :: length
!
CHARACTER(LEN=PREC_STRING), DIMENSION(MAXW) :: cpara
INTEGER           , DIMENSION(MAXW) :: lpara
INTEGER                             :: ianz
INTEGER, PARAMETER :: NOPTIONAL = 3
INTEGER, PARAMETER :: O_CODE    = 1
INTEGER, PARAMETER :: O_INSTALL = 2
INTEGER, PARAMETER :: O_PREPARE = 3
CHARACTER(LEN=1024), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=1024), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent!opt. para is present
REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 0 ! Number of values to calculate 
!
DATA oname  / 'code', 'install', 'prepare'/
DATA loname /  4    ,  7       ,  7       /
opara  =  (/ 'pre      ', 'fetch    '  , 'libraries'/)   ! Always provide fresh default values
lopara =  (/  3,           5           ,  9         /)
owerte =  (/  0.0,         0.0         ,  0.0       /)
!
if(gen_mpi_active) return               ! No update checks while mpi is active
!                                                                       
CALL get_params (zeile, ianz, cpara, lpara, maxw, lp)
IF (ier_num.ne.0) RETURN
!
CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
code_str = 'code=pre'
inst_str = 'install=fetch'
prep_str = 'prepare=libraries'
CALL do_low(opara(O_CODE))
CALL do_low(opara(O_INSTALL))
CALL do_low(opara(O_PREPARE))
IF(opara(O_CODE) == 'pre') THEN
   code_str = 'code=pre'
ELSEIF(opara(O_CODE) == 'git') THEN
   code_str = 'code=git'
ELSEIF(opara(O_CODE) == 'current') THEN
   code_str = 'code=current'
ELSEIF(opara(O_CODE)(lopara(O_CODE)-5:lopara(O_CODE)) == 'tar.gz') THEN
   code_str = 'code=' // opara(O_CODE)(1:lopara(O_CODE))
ELSE
   ier_num = -6
   ier_typ = ER_FORT
   ier_msg(1) = 'Optional parameter ''code'' must be ''pre'', ''git'' '
   ier_msg(2) = '  ''current'' or a file ending in tar.gz'
   RETURN
ENDIF
IF(opara(O_INSTALL) == 'fetch') THEN
   inst_str = 'install=fetch'
ELSEIF(opara(O_INSTALL) == 'local') THEN
   inst_str = 'install=local'
ELSEIF(opara(O_INSTALL)(lopara(O_INSTALL)-5:lopara(O_INSTALL)) == 'tar.gz') THEN
   inst_str = 'install=' // opara(O_INSTALL)(1:lopara(O_INSTALL))
ELSE
   ier_num = -6
   ier_typ = ER_FORT
   ier_msg(1) = 'Optional parameter ''install'' must be ''fetch'' or ''local'' '
   ier_msg(2) = 'or a file ending in tar.gz'
   RETURN
ENDIF
!
IF(opara(O_PREPARE) == 'libraries') THEN
   prep_str = 'prepare=libraries'
ELSEIF(opara(O_PREPARE) == 'none') THEN
   prep_str = 'prepare=none'
ELSE
   ier_num = -6
   ier_typ = ER_FORT
   ier_msg(1) = 'Optional parameter ''prepare'' must be ''libraries'' or ''none'' '
   RETURN
ENDIF
!
discus_power  ='/tmp/DISCUS_POWER'
WRITE(discus_version,'(a,a)') tmp_dir(1:len_trim(tmp_dir)),'/DISCUS_VERSION' ! Initiate search for new version
!
IF(terminal_wrp /= ' ') THEN
   OPEN(UNIT=IRD, FILE=terminal_wrp, STATUS='unknown')
   WRITE(IRD, '(a)') '#/bin/bash'
   WRITE(IRD, '(a)') 'if [ $# -eq 0 ]; then "${SHELL:-sh}"; else "$@"; fi '
   WRITE(IRD, '(a)') 'echo "The DISCUS_UPDATE exited with status $?. Press any key to close the terminal." '
   WRITE(IRD, '(a)') 'stty -icanon; dd ibs=1 count=1 >/dev/null 2>&1 '
   CLOSE(UNIT=IRD)
   line = 'chmod ugo+rwx ' // terminal_wrp(1:LEN_TRIM(terminal_wrp))
   CALL execute_command_line(line(1:LEN_TRIM(line)), wait=.true., &
        CMDSTAT=ier_cmd, CMDMSG=message, EXITSTAT=exit_msg)
ENDIF
!
IF(operating == OS_LINUX) then
   grep    = 'grep -m 1 -Poe'
   script  = 'bbb_install_script.sh'
   command = terminal_emu(1:LEN_TRIM(terminal_emu)) // ' '// &
             terminal_exe(1:LEN_TRIM(terminal_exe)) // ' '// &
             terminal_wrp(1:LEN_TRIM(terminal_wrp)) // ' '// &
             ' $HOME/' // script(1:LEN_TRIM(script)) //      &
             ' started=native ' //                           &
             code_str(1:LEN_TRIM(code_str)) // ' ' //        &
             prep_str(1:LEN_TRIM(prep_str)) // ' ' //        &
             inst_str(1:LEN_TRIM(inst_str))
!
elseif(operating == OS_LINUX_WSL) then
   grep    = 'grep -m 1 -Poe'
   script  = 'bbb_install_script.sh'
   command = terminal_emu(1:LEN_TRIM(terminal_emu)) // ' '// &
             terminal_exe(1:LEN_TRIM(terminal_exe)) // ' '// &
             terminal_wrp(1:LEN_TRIM(terminal_wrp)) // ' '// &
             ' $HOME/' // script(1:LEN_TRIM(script)) //      &
             ' started=powershell ' //                       &
             code_str(1:LEN_TRIM(code_str)) // ' ' //        &
             prep_str(1:LEN_TRIM(prep_str)) // ' ' //        &
             inst_str(1:LEN_TRIM(inst_str))
   command = 'wt.exe ''' // 'C:\Users\' //  &
              user_profile(1:LEN_TRIM(user_profile))  &
              // '\DISCUS_INSTALLATION\ccc_install_suite_Windows10_WSL.bat'' '
!
ELSEIF(operating == OS_MACOSX) THEN
   WRITE(discus_version,'(a,a)') tmp_dir(1:len_trim(tmp_dir)),'/DISCUS_VERSION' ! Initiate search for new version
   grep    = 'grep -m 1 -oe '
   script  = 'bbb_install_script_mac.sh'
!  command ='$HOME./' // script(1:LEN_TRIM(script)) // ' started=native'
!  command = terminal_emu(1:LEN_TRIM(terminal_emu)) // ' '// &
!            terminal_exe(1:LEN_TRIM(terminal_exe)) // ' '// &
!            terminal_wrp(1:LEN_TRIM(terminal_wrp)) // ' '// &
   OPEN(UNIT=IWR,FILE='/tmp/bbb.sh', STATUS='unknown')
   command = ' $HOME/' // script(1:LEN_TRIM(script)) //      &
             ' started=native ' //                           &
             code_str(1:LEN_TRIM(code_str)) // ' ' //        &
             prep_str(1:LEN_TRIM(prep_str)) // ' ' //        &
             inst_str(1:LEN_TRIM(inst_str))
   WRITE(IWR, '(a)' ) '#!/bin/zsh'
   WRITE(IWR, '(a)' ) command(1:LEN_TRIM(command))
   CLOSE(UNIT=IWR)
   command = 'chmod 700 /tmp/bbb.sh'
   CALL execute_command_line(command(1:LEN_TRIM(command)), wait=.true., &
        CMDSTAT=ier_cmd, CMDMSG=message, EXITSTAT=exit_msg)
   command = 'open -b com.apple.terminal /tmp/bbb.sh' 
ENDIF
!
! Get latest DISCUS Version
!
string = 'curl -k --silent --location https://github.com/tproffen/DiffuseCode/releases/latest ' // &
         '| grep "Release " ' // &
         '| ' // grep(1:LEN_TRIM(grep)) // ' ''v\.[0-9]*\.[0-9]*\.[0-9]*'' > ' //  discus_version
CALL execute_command_line(string(1:LEN_TRIM(string)), wait=.true., &
        CMDSTAT=ier_cmd, CMDMSG=message, EXITSTAT=exit_msg)
!
OPEN(UNIT=IRD,FILE=discus_version, STATUS='old')
READ(IRD,'(a)', iostat=ios) verstring
CLOSE(UNIT=IRD)
if(ios/=0) then
   ier_num =  -18
   ier_typ = ER_COMM
   ier_msg(1) = 'No internet connection ?'
   return
endif
!
! Download latest installation script
!
string = 'curl -k -o $HOME/' // script(1:LEN_TRIM(script)) //                       &
         ' -fSL https://github.com/tproffen/DiffuseCode/releases/download/' //   &
         verstring(1:LEN_TRIM(verstring)) // '/' // script(1:LEN_TRIM(script))
CALL execute_command_line(string(1:LEN_TRIM(string)), wait=.true., &
        CMDSTAT=ier_cmd, CMDMSG=message, EXITSTAT=exit_msg)
!
string = 'chmod 700 $HOME/' // script(1:LEN_TRIM(script))
call execute_command_line(string(1:LEN_TRIM(string)), wait=.true., &
        CMDSTAT=ier_cmd, CMDMSG=message, EXITSTAT=exit_msg)
!string = 'ls -l $HOME/' // script(1:LEN_TRIM(script))
!CALL EXECUTE_COMMAND_LINE(string(1:LEN_TRIM(string)), CMDSTAT=ier_num, CMDMSG=message, EXITSTAT=exit_msg)
!
! For LINUX_WSL we need to step into 'C:\Users\...\DISCUS_INSTALLATION\'
IF(operating == OS_LINUX_WSL) THEN
!
! Download latest installation script
!
   script  = 'bbb_install_suite_Windows10_WSL.ps1'
   string = 'curl -k -o $HOME/' // script(1:LEN_TRIM(script)) //                       &
            ' -fSL https://github.com/tproffen/DiffuseCode/releases/download/' //   &
            verstring(1:LEN_TRIM(verstring)) // '/' // script(1:LEN_TRIM(script))
   CALL execute_command_line(string(1:LEN_TRIM(string)), wait=.true., &
        CMDSTAT=ier_cmd, CMDMSG=message, EXITSTAT=exit_msg)
   string = 'cp $HOME/' // script(1:LEN_TRIM(script)) // ' /mnt/c/Users/' //          &
            user_profile(1:LEN_TRIM(user_profile)) // '/DISCUS_INSTALLATION/'  
   CALL execute_command_line(string(1:LEN_TRIM(string)), wait=.true., &
        CMDSTAT=ier_cmd, CMDMSG=message, EXITSTAT=exit_msg)
   script  = 'ccc_install_suite_Windows10_WSL.bat'
   string = 'curl -k -o $HOME/' // script(1:LEN_TRIM(script)) //                       &
            ' -fSL https://github.com/tproffen/DiffuseCode/releases/download/' //   &
            verstring(1:LEN_TRIM(verstring)) // '/' // script(1:LEN_TRIM(script))
   string = 'curl -k -o $HOME/' // script(1:LEN_TRIM(script)) //                       &
            ' -fSL https://github.com/rneder/DiffuseSuplement/releases/download/v.1.0.0/'  &
            // script(1:LEN_TRIM(script))
   CALL execute_command_line(string(1:LEN_TRIM(string)), wait=.true., &
        CMDSTAT=ier_cmd, CMDMSG=message, EXITSTAT=exit_msg)
   string = 'cp $HOME/' // script(1:LEN_TRIM(script)) // ' /mnt/c/Users/' //          &
            user_profile(1:LEN_TRIM(user_profile)) // '/DISCUS_INSTALLATION/'  
   CALL execute_command_line(string(1:LEN_TRIM(string)), wait=.true., &
        CMDSTAT=ier_cmd, CMDMSG=message, EXITSTAT=exit_msg)
!   string = '/mnt/c/Users/' // user_profile(1:LEN_TRIM(user_profile)) // '/DISCUS_INSTALLATION\'
!   length = LEN_TRIM(string)
!   CALL do_chdir (string, length, .FALSE.)
ENDIF
!write(*,*) 'POW : >', command(1:len_trim(command)),'<<'
!
! Check for other DISCUS instances / pgxwin_server
!
call lib_terminate(flag)
ex_do_exit = .false.
!
if(flag=='CONTINUE') then
!
! Finally run the DISCUS update via the installation script
!
   cdir = current_dir
   if(operating == OS_LINUX_WSL) then
      line = '/mnt/c/Users/' // user_profile(1:LEN_TRIM(user_profile)) // '/DISCUS_INSTALLATION'
   else
      line = home_dir
   endif
   length = len_trim(line)
   CALL do_chdir(line,length,.FALSE.)     ! Go to home dir
   CALL execute_command_line(command(1:LEN_TRIM(command)), WAIT=.FALSE.) !,                  CMDMSG=message, EXITSTAT=exit_msg)
!
! As installation script updated the operating system, touch update file
!
!   IF(operating == OS_LINUX .OR. operating == OS_LINUX_WSL) THEN
!      string = 'date --iso-8601 > ' // tmp_dir(1:LEN_TRIM(tmp_dir)) // '/DISCUS_UPDATE'
!   ELSEIF(operating == OS_MACOSX) THEN
!      string = 'date +%Y-%m-%d > ' // tmp_dir(1:LEN_TRIM(tmp_dir)) // '/DISCUS_UPDATE'
!   ENDIF
!   CALL EXECUTE_COMMAND_LINE (string(1:LEN_TRIM(string)), CMDSTAT=ier_num, CMDMSG=message, EXITSTAT=exit_msg)
!   string = 'chmod ugo+rw ' // tmp_dir(1:len_trim(tmp_dir)) // '/DISCUS_UPDATE'
!   CALL EXECUTE_COMMAND_LINE (string(1:LEN_TRIM(string)), CMDSTAT=ier_num, CMDMSG=message, EXITSTAT=exit_msg)
   IF(operating == OS_LINUX_WSL) THEN
      WRITE(output_io,'(a)') ' '
      WRITE(output_io,'(a)') ' This DISCUS Window will stop now. '
      WRITE(output_io,'(a)') ' The new version will be started '
      WRITE(output_io,'(a)') ' once the update is finished. '
      WRITE(output_io,'(a)') ' '
      close(output_io)
      zeile = 'stop'
      ex_do_exit = .true.
   ELSE
!    IF(operating_name=="ManjaroLinux" ) THEN
      WRITE(output_io,'(a)') ' '
      WRITE(output_io,'(a)') 'This DISCUS Window will stop now. '
      WRITE(output_io,'(a)') 'You might need an ENTER to recover the terminal prompt '
      WRITE(output_io,'(a)') ' '
      close(output_io)
      ex_do_exit = .true.
!    ELSE
!    WRITE(output_io,*) ' started script'
!    ENDIF
   ENDIF
!  line = cdir                           ! Return to old current directory
!  length = len_trim(line)
!  CALL do_chdir(line,length,.FALSE.)
endif
!
END SUBROUTINE lib_f90_update_discus
!
!*******************************************************************************
!
subroutine lib_terminate(flag)
!-
!  Checks if other discus_instances are running and asks to close these 
!  respectively asks if its OK to kill
!  If all are closed, pgxwin_server is closed
!+
!
use blanks_mod
use envir_mod
use precision_mod
use prompt_mod
use string_convert_mod
use str_comp_mod
!
implicit none
character(len=*), intent(out) :: flag
!
integer, parameter :: IRD = 37
!
character(len=PREC_STRING) :: pid_file
character(len=PREC_STRING) :: string
character(len=PREC_STRING) :: line
character(len=PREC_STRING) :: message
integer :: exit_msg
integer :: ier_cmd
integer :: j
integer :: ll
integer :: ios
logical :: lsingle
!
pid_file = tmp_dir(1:len_trim(tmp_dir)) // '/DISCUS.PIDS'
flag = 'CONTINUE'
!
loop_check: do
   lsingle = .true.
   if(operating == OS_LINUX_WSL) then
      string = 'ps --cols 1024 aux | fgrep discus_suite | fgrep -v grep '       &
               // '| grep -v ubuntu | grep -v terminal | awk ''{print $2} '' > ' &
               // pid_file(1:len_trim(pid_file))
   else
      string = 'ps --cols 1024 aux | fgrep discus_suite | fgrep -v grep '       &
               //                                    ' | awk ''{print $2} '' > ' &
               // pid_file(1:len_trim(pid_file))
   endif
   call execute_command_line(string, wait=.true., &
        CMDSTAT=ier_cmd, CMDMSG=message, EXITSTAT=exit_msg)
!  call execute_command_line('cat /tmp/DISCUS.PIDS')
   open(unit=IRD, file=pid_file, status='old')
   loop_test: do
      read(IRD,'(a)', iostat=ios) string
      if(is_iostat_end(ios)) exit loop_test
      read(string,*) j
      if(j /= PID ) then
         lsingle = .false.
         exit loop_test
      endif
   enddo loop_test
   close(unit=IRD)
   if(lsingle) then
      flag = 'CONTINUE'
      exit loop_check
   endif
!
   flag = 'CANCEL'
   line = 'ENTER'
   write(output_io,'(a)') ' '
   write(output_io,'(a)') ' Other DISCUS instances are running, which need to be closed.'
   write(output_io,'(a)') ' Close all other windows and hit ENTER key in this window'
   write(output_io,'(a)') ' or enter ''kill'' to have other instances stopped without save!'
   write(output_io,'(a)') ' or enter ''cancel'' to cancel the update'
   write(output_io,'(a)', advance='NO') ' Enter / Kill / Cancel '
   read(*, '(a)', iostat=ios) line
   call do_cap(line)
   ll = len_trim(line)
   call rem_leading_bl(line, ll)
   if(str_comp(line, 'CANCEL', 1, len_trim(line), 6)) then
      flag = 'CANCEL'
      exit loop_check
   elseif(str_comp(line, 'KILL', 1, len_trim(line), 4)) then
      open(unit=IRD, file=pid_file, status='old')
      loop_kill: do
         read(IRD,'(a)', iostat=ios) string
         if(is_iostat_end(ios)) exit loop_kill
         read(string,*) j
         if(j /= PID ) then
            write(line,'(a,i8)') 'kill -9 ',j
            call execute_command_line(line, wait=.true., &
                 CMDSTAT=ier_cmd, CMDMSG=message, EXITSTAT=exit_msg)
         endif
      enddo loop_kill
      close(IRD)
      flag = 'CONTINUE'
      exit loop_check
   endif
enddo loop_check
!
!  Close down pgxwin_server
!
if(flag=='CONTINUE') then       ! KIll pgxwin_server
   string = 'ps --cols  256 aux | fgrep pgxwin_server | fgrep -v grep | awk ''{print $2} '' > ' &
            // pid_file(1:len_trim(pid_file))
   call execute_command_line(string, wait=.true., &
        CMDSTAT=ier_cmd, CMDMSG=message, EXITSTAT=exit_msg)
   open(unit=IRD, file=pid_file, status='old')
   loop_pgx: do                 ! In case the server runs multiple times
      string = ' '
      read(IRD,'(a)', iostat=ios) string
      if(is_iostat_end(ios)) exit loop_pgx
      read(string,*) j
      if(j /= PID ) then
         write(line,'(a,i8)') 'kill -9 ',j
         call execute_command_line(line, wait=.true., &
              CMDSTAT=ier_cmd, CMDMSG=message, EXITSTAT=exit_msg)
      endif
   enddo loop_pgx
   close(IRD)
endif
!
end subroutine lib_terminate
!
!*******************************************************************************
!
SUBROUTINE lib_f90_findterminal
!-
! Try to find a terminal emulator installed in this system
!+
USE envir_mod
USE errlist_mod
USE precision_mod
USE lib_do_operating_mod
!
IMPLICIT NONE
!
INTEGER, PARAMETER :: MAXTERM=4
INTEGER, PARAMETER :: IDEF = 68
!
CHARACTER(LEN=16), DIMENSION(0:MAXTERM-1) :: terminals
CHARACTER(LEN=16), DIMENSION(0:MAXTERM-1) :: terminal_execute
CHARACTER(LEN=128        ) :: cfile
CHARACTER(LEN=PREC_STRING) :: line
CHARACTER(LEN=PREC_STRING) :: message
INTEGER             :: exit_msg
integer             :: ier_cmd   ! execute_command_line message
INTEGER :: ios
INTEGER :: length
INTEGER :: i,j
INTEGER :: istart
DATA terminals /'gnome-terminal  ',  &
                'terminator      ',  &
                'xterm           ',  &
                'konsole         '   &
               /
DATA terminal_execute /' --             ', &
                       ' --execute      ', &
                       ' -e             ', &
                       ' --hold -e      '  &
                      /
!
j = 0
IF(operating==OS_LINUX) THEN
  IF(operating_name=="ManjaroLinux" ) THEN
     istart = MAXTERM-1
  ELSE
     istart = 0
  ENDIF
ELSEIF(operating==OS_LINUX_WSL) THEN
  istart = 1
ELSEIF(OPERATING==OS_MACOSX) THEN
  istart = 2
ELSE
  istart = 0
ENDIF
!
WRITE(cfile,'(a,a,i10.10)') tmp_dir(1:LEN_TRIM(tmp_dir)), '/DISCUS_TERM.', PID
main_loop: DO i=istart, istart+MAXTERM
   j = MOD(i,MAXTERM)
   line = 'which ' // terminals(j) // ' > ' // cfile
   CALL execute_command_line(line, wait=.true., &
        CMDSTAT=ier_cmd, CMDMSG=message, EXITSTAT=exit_msg)
   OPEN(UNIT=IDEF, FILE=cfile, STATUS='OLD', IOSTAT=ios)
   IF(ios/=0) CYCLE main_loop
   line = ' '
   READ(IDEF, '(a)', IOSTAT=ios) line
   IF(IS_IOSTAT_END(ios)) CYCLE main_loop
   IF(LINE /= ' ') EXIT main_loop
ENDDO main_loop
CLOSE(UNIT=IDEF)
IF(line== ' ') THEN
   terminal_emu = 'NONE'
   terminal_exe = ' '
ELSE
   terminal_emu = line(1:len_trim(line))
   terminal_exe = terminal_execute(j)
   IF(j==MAXTERM-1) terminaL_wrp  = ' '     ! No terminal_wrapper for konsole
ENDIF 
WRITE(line,'(a,a)') 'rm -f ', cfile(1:LEN_TRIM(cfile))
length = LEN_TRIM(line)
CALL do_operating(line, length)
!
END SUBROUTINE lib_f90_findterminal
!
!*******************************************************************************
!
logical function check_github() result(lnetz)
!-
!  Check if a ping to github.com = 140.82.121.3 is successfull with less than 
!  100% packet loss
!  analyses lines like: 
!  "1 packets transmitted, 0 packets received, 100% packet loss"
!  "1 packets transmitted, 1 packets received, 0% packet loss"
!  "1 packets transmitted, 0 packets received, 100.0% packet loss"   ! Mac Version
!  "1 packets transmitted, 1 packets received, 0.0% packet loss"     ! Mac Version
!  Results are language dependent!, search for number prior to "%"
!+
!
use errlist_mod
use envir_mod
use precision_mod
!
implicit none
!
integer, parameter         :: IRD = 33
!
character(len=PREC_STRING) :: string    ! Generic string
character(len=PREC_STRING) :: string_e  ! Generic string
character(len=PREC_STRING) :: cfile     ! temporary file
CHARACTER(LEN=PREC_STRING) :: message   ! execute_command_line message
integer                    :: exit_msg  ! execute_command_line message
integer                    :: ier_cmd   ! execute_command_line message
integer                    :: ios       ! read error number
!integer                    :: iper      ! location of percent sign
!integer                    :: icom      ! prior location of ','
!integer                    :: idec      ! location of decimal point needed for MacOS!
integer                    :: i         ! packet loss percentage Linux WSL version
real(kind=PREC_SP)         :: r         ! packet loss percentage MacOs     version
!
lnetz = .false.
r = 0.0
i = 0
cfile = tmp_dir(1:tmp_dir_l) // '/check_github.txt'
string = 'ping -c 1 -W 3 140.82.121.3 2> /dev/null > '// cfile(1:len_trim(cfile)) &
         // '; echo $? >> ' // cfile(1:len_trim(cfile)) 
! Include a wait, to make reasonably sure that the temp file is written
call execute_command_line(string, CMDSTAT=ier_cmd, CMDMSG=message, EXITSTAT=exit_msg, wait=.true.)
!
string_e = '2'
open(unit=ird, file=cfile, status='old')
loop_check: do                          ! Loop to find the line
   read(ird, '(a)', iostat=ios) string
   if(is_iostat_end(ios)) exit loop_check
   string_e = string
!  iper = index(string, '%')            ! Locate percent sign
!  if(iper>0) then 
!     icom = index(string(1:iper-1),',', back=.TRUE.)  ! Locate prior comma
!     idec = index(string(icom+1:iper-1),'.')          ! Mac has decimal point!
!     if(idec>0) then
!        read(string(icom+1:iper-1),*) r               ! Real value for Mac
!        i = nint(r)
!     else
!        read(string(icom+1:iper-1),*) i               ! Integernumber for Linux / WSL
!     endif
!     if(i<100) then                                   ! Less than 100% loss, found network
!        lnetz = .true.
!     endif
!     exit loop_check
!  endif
enddo loop_check
!
close(unit=ird)
i = 2
read(string_e,'(i1)', iostat=ios) i
lnetz = i==0
!write(*,*) ' ERROR CODE ', string_e(1:len_trim(string_e)), i, lnetz
string = "rm -f " // cfile(1:LEN_TRIM(cfile))          ! Remove temporary file
call execute_command_line(string, wait=.true., &
        CMDSTAT=ier_cmd, CMDMSG=message, EXITSTAT=exit_msg)
!
end function check_github
!
!*******************************************************************************
!
END MODULE appl_env_mod
