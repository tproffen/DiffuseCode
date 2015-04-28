!*****7***************************************************************  
!                                                                       
!                                                                       
      SUBROUTINE appl_env 
!-                                                                      
!     Reads environment variables, sets path for helpfile               
!     UNIX version ..                                                   
!+                                                                      
      USE envir_mod 
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      CHARACTER(255) cdummy
      INTEGER ico, ice, iii
      INTEGER len_str 
      INTEGER pname_l 
!                                                                       
      operating = ' ' 
      CALL get_environment_variable ('OS', operating) 
!
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
      appl_dir = ' ' 
      CALL get_environment_variable (pname_cap, appl_dir) 
      IF (appl_dir.eq.' ') then 
         appl_dir = '.' 
      ENDIF 
!
      CALL get_environment_variable ('_', cdummy) 
      iii=index(cdummy,pname,.true.)
      appl_dir=cdummy(1:iii-1)
      appl_dir_l = len_str (appl_dir) 
!                                                                       
      deffile = '.'//pname (1:pname_l) 
      deffile_l = len_str (deffile) 
!                                                                       
      mac_dir = ' ' 
      mac_dir (1:appl_dir_l) = appl_dir 
      mac_dir (appl_dir_l + 1:appl_dir_l + pname_l + 10) = '../share/'//     &
      pname (1:pname_l) //'/'                                           
      mac_dir_l = len_str (mac_dir) 
!                                                                       
      umac_dir = home_dir(1:home_dir_l)//'/mac/'//pname(1:pname_l) //'/'
      umac_dir_l = len_str (umac_dir) 
!                                                                       
      nullfile = '/dev/null' 
!                                                                       
      hlpfile = ' ' 
      hlpfile (1:appl_dir_l) = appl_dir 
      hlpfile (appl_dir_l + 1:appl_dir_l + 1 + pname_l + 13) = '../share/'//     &
      pname (1:pname_l) //'.hlp'                                        
      hlpfile_l = len_str (hlpfile) 
!                                                                       
      colorfile = ' ' 
      colorfile (1:appl_dir_l) = appl_dir 
      colorfile (appl_dir_l + 1:appl_dir_l + 19) = '../share/color.map' 
      colorfile_l = len_str (colorfile) 
!
      IF(index(operating, 'Windows') /= 0) THEN  ! We got a Windows
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
      WRITE ( *, 1000) umac_dir (1:umac_dir_l) 
      WRITE ( *, 1100) mac_dir (1:mac_dir_l) 
      WRITE ( *, 1200) start_dir (1:start_dir_l) 
!                                                                       
 1000 FORMAT     (1x,'User macros in   : ',a) 
 1100 FORMAT     (1x,'System macros in : ',a) 
 1200 FORMAT     (1x,'Start directory  : ',a) 
      END SUBROUTINE appl_env                       
