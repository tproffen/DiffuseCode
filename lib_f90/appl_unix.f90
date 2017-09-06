!*****7***************************************************************  
!                                                                       
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
      IMPLICIT none 
!                                                                       
      LOGICAL, INTENT(IN) :: standalone
      INTEGER, INTENT(IN) :: local_mpi_myid
!                                                                       
      INTEGER, PARAMETER :: idef = 68
!
      CHARACTER(255) cdummy
      CHARACTER(LEN=1024) :: line
      INTEGER ico, ice, iii, i, j
      INTEGER :: ios ! I/O status
      INTEGER len_str 
      INTEGER pname_l 
      LOGICAL lpresent
!                                                                       
      operating = ' ' 
      CALL get_environment_variable ('OS', operating) 
      IF(operating == '') THEN
         CALL get_environment_variable ('OSTYPE', operating) 
      ENDIF
!
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
            ENDIF
         ENDDO name_search
         CLOSE(idef)
      ENDIF
      pname_l = len_str (pname) 
      home_dir = ' ' 
      lines = 42 
      CALL get_environment_variable ('LINES', home_dir) 
      IF (home_dir.ne.' ') then 
         READ (home_dir, *, end = 10) lines 
   10    CONTINUE 
      ELSE 
         CALL get_environment_variable ('TERMCAP', home_dir) 
         ico = index (home_dir, 'co') + 3 
         ice = index (home_dir (ico:256) , ':') + ico - 2 
         IF (ice.gt.ico) then 
            READ (home_dir (ico:ice), *, end = 20, err = 20) lines 
         ENDIF 
   20    CONTINUE 
      ENDIF 
      lines = lines - 2 
!
      IF(index(operating, 'Windows') /= 0) THEN  ! We got a Windows
         user_profile = ' '
         CALL get_environment_variable ('USERPROFILE', user_profile)
         home_dir = user_profile
      ELSE
!                                                                       
         home_dir = ' ' 
         CALL get_environment_variable ('HOME', home_dir) 
         IF (home_dir.eq.' ') then 
            home_dir = '.' 
         ENDIF 
      ENDIF
      home_dir_l = len_str (home_dir) 
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
      END SUBROUTINE appl_env                       
!
SUBROUTINE write_appl_env (standalone, local_mpi_myid)
!
USE envir_mod
USE terminal_mod
IMPLICIT NONE
!                                                                       
LOGICAL, INTENT(IN) :: standalone
INTEGER, INTENT(IN) :: local_mpi_myid
!
IF(standalone .AND. local_mpi_myid==0) THEN
   IF(term_scheme_exists) THEN
      WRITE ( *, 2000) TRIM(color_bg),TRIM(color_info),umac_dir (1:umac_dir_l),   TRIM(color_fg)
      WRITE ( *, 2100)                TRIM(color_info),mac_dir (1:mac_dir_l),     TRIM(color_fg)
      WRITE ( *, 2200)                TRIM(color_info),start_dir (1:start_dir_l) ,TRIM(color_fg)
   ELSE
      WRITE ( *, 1000) umac_dir (1:umac_dir_l) 
      WRITE ( *, 1100) mac_dir (1:mac_dir_l) 
      WRITE ( *, 1200) start_dir (1:start_dir_l) 
   ENDIF
ENDIF
!                                                                       
 1000 FORMAT     (1x,'User macros in   : ',a) 
 1100 FORMAT     (1x,'System macros in : ',a) 
 1200 FORMAT     (1x,'Start directory  : ',a) 
 2000 FORMAT     (1x,a,'User macros in   : ',a,a,a) 
 2100 FORMAT     (1x,  'System macros in : ',a,a,a) 
 2200 FORMAT     (1x,  'Start directory  : ',a,a,a) 
!
END SUBROUTINE write_appl_env                       
!
SUBROUTINE color_set_scheme (standalone, local_mpi_myid)
!
!  Set the terminal color scheme 
!  If the file share/discus.term.scheme exists it is used
!  Else we try a default color scheme according to the
!  Operating system
!
USE terminal_mod
USE envir_mod
USE param_mod
USE prompt_mod
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
ENDIF
END SUBROUTINE color_set_scheme
!
SUBROUTINE color_set_fg(line,icolon,color_type , color_string, color_color)
!
USE errlist_mod
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
REAL                , DIMENSION(1) :: werte
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
SUBROUTINE color_set_bg(line,icolon,color_type , color_string, color_color)
!
USE errlist_mod
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
REAL                , DIMENSION(1) :: werte
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
      SUBROUTINE  program_files
!-                                                                      
!     Sets path for helpfile , mac directories
!     UNIX version ..                                                   
!+                                                                      
      USE envir_mod 
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
      INTEGER :: pname_l
!
      pname_l = LEN(TRIM(pname))
!
      mac_dir = ' ' 
      mac_dir (1:appl_dir_l) = appl_dir 
      mac_dir (appl_dir_l + 1:appl_dir_l + pname_l + 10) = '../share/'//     &
      pname (1:pname_l) //'/'                                           
      mac_dir_l = LEN(TRIM (mac_dir) )
!                                                                       
      umac_dir = home_dir(1:home_dir_l)//'/mac/'//pname(1:pname_l) //'/'
      umac_dir_l = LEN(TRIM (umac_dir) )
!                                                                       
      hlpfile = ' ' 
      hlpdir  = ' ' 
      hlpdir  (1:appl_dir_l+9) = appl_dir(1:appl_dir_l) //'../share/'
      hlp_dir_l = appl_dir_l+9
      hlpfile   = hlpdir(1:hlp_dir_l)//pname(1:pname_l)//'.hlp'
      hlpfile_l = LEN(TRIM (hlpfile) )
!                                                                       
      colorfile = ' ' 
      colorfile (1:appl_dir_l) = appl_dir 
      colorfile (appl_dir_l + 1:appl_dir_l + 19) = '../share/color.map' 
      colorfile_l = LEN(TRIM (colorfile) )
!
      END SUBROUTINE  program_files
