MODULE appl_env_mod
!*****7***************************************************************  
!                                                                       
CONTAINS
!
!*******************************************************************************!
!
SUBROUTINE appl_env (standalone, local_mpi_myid)
!-                                                                      
!     Reads environment variables, sets path for helpfile               
!     UNIX version ..                                                   
!+                                                                      
USE errlist_mod
USE envir_mod 
USE terminal_mod
USE param_mod
USE prompt_mod 
USE string_convert_mod
!
IMPLICIT none 
!                                                                       
LOGICAL, INTENT(IN) :: standalone
INTEGER, INTENT(IN) :: local_mpi_myid
!                                                                       
INTEGER, PARAMETER :: idef = 68
!
CHARACTER(255) cdummy
CHARACTER(LEN=8), DIMENSION(6), PARAMETER :: tmp_test = (/'/tmp    ','/TMP    ', &
      '/var/tmp', '/Var/tmp', '/var/TMP', '/Var/TMP' /)
CHARACTER(LEN=2048) :: line
CHARACTER(LEN=2048) :: pathfile
INTEGER ico, ice, iii, i, j
INTEGER :: length
INTEGER :: ios ! I/O status
INTEGER len_str 
INTEGER pname_l 
LOGICAL lpresent
INTEGER :: lib_f90_getpid
!
IF(envir_done) RETURN
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
!     ENDIF
      home_dir_l = len_str (home_dir) 
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
      ELSE
         INQUIRE(FILE='/suite.sh',EXIST=lpresent)
         IF(lpresent) THEN
            operating = OS_WINDOWS  ! Cygwin from DISCUS Icon
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
   ENDIF
ELSE   !  /proc/version does not exist , likely a MAC OS X 
!  Read OS from uname
   line = 'uname -av > '//tmp_dir(1:tmp_dir_l)//'/DISCUS_SUITE_UNAME'
   CALL do_operating_comm(line)
   OPEN(UNIT=idef, FILE='/tmp/DISCUS_SUITE_UNAME')
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
   ENDIF
ENDIF
!
! All remaining test for OS should be obsolete
   IF(operating == ' ') THEN
      CALL get_environment_variable ('OS', operating) 
      IF(operating == '') THEN
         CALL get_environment_variable ('OSTYPE', operating) 
      ENDIF
!
      color_theme = THEME_DEFAULT
      INQUIRE(FILE='/etc/os-release',EXIST=lpresent)
      IF(lpresent) THEN
         CALL oeffne(idef, '/etc/os-release', 'old')
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
!     IF(index(operating, 'Windows') /= 0) THEN  ! We got a Windows
!        operating = 'Windows'
!        color_theme = THEME_LGHTYL_BLACK
!        INQUIRE(FILE='/Cygwin.bat',EXIST=lpresent)
!        IF(lpresent) THEN
!           OPEN(idef,FILE='/Cygwin.bat',STATUS='old')
!           READ(idef,'(a)', IOSTAT=ios) line
!           find_op: DO WHILE(.NOT. IS_IOSTAT_END(ios))
!              IF(line(1:5)=='chdir') THEN    ! FOUND correct line
!                operat_top = line( 7:LEN_TRIM(LINE)-4)
!                operating  = line(10:LEN_TRIM(LINE)-4) !'Cygwin'
!                EXIT find_op
!              ENDIF
!              READ(idef,'(a)', IOSTAT=ios) line
!           ENDDO find_op
!           CLOSE(idef)
!        ELSE
!           INQUIRE(FILE='/suite.sh',EXIST=lpresent)
!           IF(lpresent) THEN
!              operating = 'Windows'
!           ENDIF
!        ENDIF
!        user_profile = ' '
!        CALL get_environment_variable ('USERPROFILE', user_profile)
!        home_dir = user_profile
!     ELSE
!                                                                       
!        home_dir = ' ' 
!        CALL get_environment_variable ('HOME', home_dir) 
!        IF (home_dir.eq.' ') then 
!           home_dir = '.' 
!        ENDIF 
!     ENDIF
!     home_dir_l = len_str (home_dir) 
!                                                                       
!     appl_dir = ' ' 
!     CALL get_environment_variable (pname_cap, appl_dir) 
!     IF (appl_dir.eq.' ') then 
!        appl_dir = '.' 
!     ENDIF 
!
!     Try to find the correct share/ folder
!         Start with environment variable PNAME_CAP
!         Else try from process name
!         Else assume /usr/local/bin
      IF(standalone) THEN
!
!        Environment variable
!
         appl_dir = ' ' 
         CALL get_environment_variable (pname_cap, appl_dir) 
         IF (appl_dir.eq.' ') then 
            appl_dir = '.' 
         ENDIF 
         appl_dir_l = LEN_TRIM(appl_dir)
         IF(appl_dir(appl_dir_l:appl_dir_l) /= '/') THEN
            appl_dir   = appl_dir(1:appl_dir_l)//'/'
            appl_dir_l = appl_dir_l+ 1
         ENDIF
         hlpdir    = ' ' 
         hlpdir(1:appl_dir_l+9)=appl_dir(1:appl_dir_l)//'../share/'
         hlp_dir_l = appl_dir_l+9
         hlpfile   = hlpdir(1:hlp_dir_l)//pname(1:pname_l)//'.hlp'
         hlpfile_l = LEN(TRIM (hlpfile) )
         CALL do_fexist(hlpfile,hlpfile_l,.FALSE.)
         IF(res_para(1)==0) THEN   ! hlpdir does not exist
!
!           Try /usr/local/bin
!
            appl_dir   = '/usr/local/bin/'
            appl_dir_l = 15
            hlpdir    = ' ' 
            hlpdir(1:appl_dir_l+9) = appl_dir(1:appl_dir_l)//'../share/'
            hlp_dir_l = appl_dir_l+9
            hlpfile   = hlpdir(1:hlp_dir_l)//pname(1:pname_l)//'.hlp'
            hlpfile_l = LEN(TRIM (hlpfile) )
            CALL do_fexist(hlpfile,hlpfile_l,.FALSE.)
            IF(res_para(1)==0) THEN   ! hlpfile does not exist
!
!              Process folder
!
               CALL GET_ENVIRONMENT_VARIABLE ('_', cdummy) 
               iii = INDEX(cdummy,pname,.true.)
               i   = INDEX(cdummy,'/',.TRUE.)
               j   = INDEX(cdummy,'\',.TRUE.)
               iii = MAX(i,j)
               appl_dir   = cdummy(1:iii  )
               appl_dir_l = len_str (appl_dir) 
               hlpdir    = ' ' 
               hlpdir(1:appl_dir_l+9) = appl_dir(1:appl_dir_l)//'../share/'
               hlp_dir_l = appl_dir_l+9
               hlpfile   = hlpdir(1:hlp_dir_l)//pname(1:pname_l)//'.hlp'
               hlpfile_l = LEN(TRIM (hlpfile) )
               CALL do_fexist(hlpfile,hlpfile_l,.FALSE.)
               IF(res_para(1)==0) THEN   ! hlpfile does not exist
!
!                 last resort, take home dir
!
                  appl_dir   = home_dir
                  appl_dir_l = home_dir_l
               ENDIF
            ENDIF
         ENDIF
      ENDIF
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
         CALL do_chdir ( start_dir, start_dir_l, .false.)
      ELSEIF(operating==OS_LINUX_WSL) THEN
         CALL do_cwd (start_dir, start_dir_l) 
         IF(start_dir == '/mnt/c/Users' ) THEN
!
!           Started from Icon set start directory to
!           User HOME
            line = 'echo $PATH > '//tmp_dir(1:tmp_dir_l)//'/discus_suite_path.txt'
            CALL do_operating_comm(line)
            pathfile = tmp_dir(1:tmp_dir_l)//'/discus_suite_path.txt'
            OPEN(UNIT=idef,FILE=pathfile,ACTION='READ')
            READ(idef,'(a)') line
            CLOSE(UNIT=idef)
            i = INDEX(line,'/mnt/c/Users')
            IF(i/=0) THEN
               j = INDEX(line(i+13:),'/')
               IF(j/=0) THEN
                  user_name = ' '
                  user_name = line(i+13:i+13+j-2)
                  IF(start_dir(start_dir_l:start_dir_l)=='/') THEN
                     start_dir_l = start_dir_l - 1
                  ENDIF
                  start_dir = start_dir(1:start_dir_l) // '/' // &
                              user_name(1:LEN_TRIM(user_name))
                  start_dir_l = LEN_TRIM(start_dir)
                  CALL do_chdir (start_dir, start_dir_l, .TRUE.) 
               ENDIF
            ENDIF
         ENDIF
         IF(start_dir(start_dir_l:start_dir_l) /= '/') THEN
            start_dir   = start_dir(1:start_dir_l) // '/'
            start_dir_l = start_dir_l + 1
         ENDIF
      ELSE
!                                                                       
         CALL do_cwd (start_dir, start_dir_l) 
         IF(start_dir(start_dir_l:start_dir_l) /= '/') THEN
            start_dir   = start_dir(1:start_dir_l) // '/'
            start_dir_l = start_dir_l + 1
         ENDIF
      ENDIF
      current_dir   = start_dir
      current_dir_l = start_dir_l
!
!     Define terminal color scheme
!
      CALL color_set_scheme (standalone, local_mpi_myid)
!
! Get Parent Process ID, required knowledge of operating
!
PPID = lib_f90_getppid(PID)
!
envir_done = .TRUE.
!                                                                       
END SUBROUTINE appl_env                       
!
!*******************************************************************************!
!
SUBROUTINE write_appl_env (standalone, local_mpi_myid)
!-
!  Writes the path for help manual etc to the welcome screen
!+
USE envir_mod
USE terminal_mod
IMPLICIT NONE
!                                                                       
LOGICAL, INTENT(IN) :: standalone
INTEGER, INTENT(IN) :: local_mpi_myid
!
IF(standalone .AND. local_mpi_myid==0) THEN
   IF(term_scheme_exists) THEN
      WRITE ( *, 1900) TRIM(color_bg),TRIM(color_info), man_dir (1:LEN_TRIM(man_dir)) ,TRIM(color_fg)
      WRITE ( *, 2000) TRIM(color_bg),TRIM(color_info),umac_dir (1:LEN_TRIM(umac_dir)),TRIM(color_fg)
      WRITE ( *, 2100) TRIM(color_bg),TRIM(color_info),mac_dir (1:mac_dir_l),     TRIM(color_fg)
      WRITE ( *, 2200) TRIM(color_bg),TRIM(color_info),start_dir (1:start_dir_l) ,TRIM(color_fg)
      WRITE ( *, 2300) TRIM(color_bg),TRIM(color_info),TRIM(color_fg)
      WRITE ( *, 2400) TRIM(color_bg),TRIM(color_info),TRIM(color_fg)
      WRITE ( *, 2500) TRIM(color_bg),TRIM(color_info),TRIM(color_fg)
      WRITE ( *, 2600) TRIM(color_bg),TRIM(color_info),TRIM(color_fg),TRIM(color_info),TRIM(color_fg)
   ELSE
      WRITE ( *,  900)  man_dir (1:LEN_TRIM(man_dir)) 
      WRITE ( *, 1000) umac_dir (1:LEN_TRIM(umac_dir)) 
      WRITE ( *, 1100) mac_dir (1:mac_dir_l) 
      WRITE ( *, 1200) start_dir (1:start_dir_l) 
      WRITE ( *, 1300)
      WRITE ( *, 1400)
      WRITE ( *, 1500)
      WRITE ( *, 1600)
   ENDIF
ENDIF
!                                                                       
  900 FORMAT     (1x,'Manual files in  : ',a) 
 1000 FORMAT     (1x,'User macros in   : ',a) 
 1100 FORMAT     (1x,'System macros in : ',a) 
 1200 FORMAT     (1x,'Start directory  : ',a) 
 1300 FORMAT     (1x,'Access manuals at each section with   : ','manual')
 1400 FORMAT     (1x,'Access help at each section/menu with : ','help  ')
 1500 FORMAT     (1x,'News at each section/Command_lang in  : ','help News')
 1600 FORMAT     (1x,'Change font size with VT Fonts menu or: ','CTRL + /CTRL -')
 1900 FORMAT     (1x,a,'Manual files in  : ',a,a,a) 
 2000 FORMAT     (1x,a,'User macros in   : ',a,a,a) 
 2100 FORMAT     (1x,a,'System macros in : ',a,a,a) 
 2200 FORMAT     (1x,a,'Start directory  : ',a,a,a) 
 2300 FORMAT     (1x,a,'Access manuals at each section with   : ',a,'manual',a)
 2400 FORMAT     (1x,a,'Access help at each section/menu with : ',a,'help  ',a)
 2500 FORMAT     (1x,a,'News at each section/Command_lang in  : ',a,'help News',a)
 2600 FORMAT     (1x,a,'Change font size with ',a,'VT Fonts',a,' menu or: ',a,'CTRL + /CTRL -',a)
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
USE prompt_mod
USE string_convert_mod
!
IMPLICIT NONE
!
LOGICAL, INTENT(IN) :: standalone
INTEGER, INTENT(IN) :: local_mpi_myid
!
INTEGER, PARAMETER :: idef = 68
CHARACTER(LEN=1024) :: line, color
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
CHARACTER (LEN=1024), DIMENSION(1) :: cpara
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
CHARACTER (LEN=1024), DIMENSION(1) :: cpara
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
USE prompt_mod 
!
IMPLICIT none 
!                                                                       
INTEGER, PARAMETER :: ird = 34
!
CHARACTER(LEN=1024) :: line
CHARACTER(LEN=7)    :: progname
INTEGER :: progname_l
INTEGER :: i
LOGICAL :: l_exist
!
IF(lstandalone) THEN
  progname = pname
ELSE
  progname = 'suite'
ENDIF
!
progname_l = LEN(TRIM(progname))
!
mac_dir = ' ' 
mac_dir (1:appl_dir_l) = appl_dir 
mac_dir (appl_dir_l + 1:appl_dir_l + progname_l + 10) = '../share/'//     &
progname (1:progname_l) //'/'                                           
mac_dir_l = LEN(TRIM (mac_dir) )
!                                                                       
umac_dir = home_dir(1:home_dir_l)//'/mac/'//progname(1:progname_l) //'/'
umac_dir_l = LEN(TRIM (umac_dir) )
!                                                                       
hlpfile = ' ' 
hlpdir  = ' ' 
hlpdir  (1:appl_dir_l+9) = appl_dir(1:appl_dir_l) //'../share/'
hlp_dir_l = appl_dir_l+9
hlpfile   = hlpdir(1:hlp_dir_l)//progname(1:progname_l)//'.hlp'
hlpfile_l = LEN(TRIM (hlpfile) )
!                                                                       
colorfile = ' ' 
colorfile (1:appl_dir_l) = appl_dir 
colorfile (appl_dir_l + 1:appl_dir_l + 19) = '../share/color.map' 
colorfile_l = LEN(TRIM (colorfile) )
!
!   Search for the Installation file "DiscusSuite.txt" to obtain paths
!
man_dir = ' '
inst_file = appl_dir(1:LEN_TRIM(appl_dir)) // '../share/DiscusSuite.txt'
INQUIRE(FILE=inst_file,EXIST=l_exist)
IF(.NOT.l_exist) THEN
   inst_file = '/usr/local/share/DiscusSuite.txt'
   INQUIRE(FILE=inst_file,EXIST=l_exist)
   IF(.NOT.l_exist) THEN
            inst_file = '/usr/share/DiscusSuite.txt'
      INQUIRE(FILE=inst_file,EXIST=l_exist)
      IF(.NOT.l_exist) THEN
         inst_file = '/share/DiscusSuite.txt'
      ENDIF
   ENDIF
ENDIF
IF(l_exist) THEN
   CALL oeffne(IRD, inst_file,'old')
   IF(ier_num==0) THEN
      READ(IRD,'(a)') line
      READ(IRD,'(a)') line
      READ(IRD,'(a)') line
      IF(line(1:13)=='Manual      :') THEN
         READ(line(15:LEN_TRIM(line)),'(a)') man_dir
      ENDIF
   ENDIF
   CLOSE(IRD)
ELSE
   man_dir = appl_dir(1:LEN_TRIM(appl_dir)) // '../share/'
ENDIF
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
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: cpid  ! current PID whose Parent Process ID is to be found
!
INTEGER, PARAMETER  :: ITMP     = 79  ! temporary unit number
!
CHARACTER(LEN=1024) :: line
CHARACTER(LEN=1024) :: temp_file
INTEGER             :: tpid           ! temporary pid for current process for verification
INTEGER             :: tppid          ! temporary ppid for current process for verification
INTEGER             :: ios            ! I/O status 
!
tppid = 0   ! Default value 
!
!  Make a temp_file to write the system(ps) into
WRITE(temp_file, '(a,I5.5)') '/tmp/getppid.', cpid
!
IF(operating(1:6)=='darwin') THEN
   WRITE(line,'(a,i8,a,a)') 'ps j | grep ',PID,' | grep -v grep | awk ''{print $2, $3}'' >> ', &
       temp_file(1:LEN_TRIM(temp_file))
   CALL system(line)
ELSEIF(operating==OS_WINDOWS .OR. operating==OS_CYGWIN64 .OR. operating==OS_CYGWIN32) THEN
   WRITE(line,'(a,i8,a,a)') 'ps j | grep ',PID,' | grep discus  | awk ''{print $1, $2}'' >> ', &
       temp_file(1:LEN_TRIM(temp_file))
   CALL system(line)
ELSE
   WRITE(line,'(a,i8,a,a)') 'ps j | grep ',PID,' | grep -v grep | awk ''{print $2, $1}'' >> ', &
       temp_file(1:LEN_TRIM(temp_file))
   CALL system(line)
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
CALL SYSTEM(line)                    ! remove temporary file
!
lib_f90_getppid = tppid
!
END FUNCTION lib_f90_getppid
!
!*******************************************************************************
!
END MODULE appl_env_mod
