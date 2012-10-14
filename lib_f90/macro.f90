!*****7**************************************************************** 
!                                                                       
      SUBROUTINE file_kdo (line, i) 
!-                                                                      
!     Opens a new macro file.                                           
!+                                                                      
      IMPLICIT none 
!                                                                       
      include'errlist.inc' 
      include'envir.inc' 
      include'macro.inc' 
      include'prompt.inc' 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = MAC_MAX_PARA + 1) 
!                                                                       
      CHARACTER ( * ) line 
      CHARACTER(1024) filename, string, cpara (maxw) 
      CHARACTER (LEN=1024) :: ldir  ! current local directory
!                                                                       
      INTEGER lpara (maxw) 
      INTEGER              :: lslash          ! position os slash in filename
      INTEGER              :: ldir_length     ! length of directory string
      INTEGER              :: filename_length ! length of filename string
      INTEGER ip, ianz, i, lc 
      LOGICAL fileda 
      REAL werte (maxw) 
!
      INTEGER   len_str
!
!---- Get current working directory
!
      CALL do_cwd ( ldir, ldir_length )
      If(ldir(ldir_length:ldir_length) /= '/') THEN
        ldir_length = ldir_length + 1
        ldir(ldir_length:ldir_length) = '/'
      ENDIF
!                                                                       
!     If there is place for another macro                               
!                                                                       
      IF (mac_level.lt.MAC_MAX_LEVEL) then 
!                                                                       
!     --Get filename from command line and string for parameters        
!                                                                       
         string = line 
         ip = index (string, ' ') 
         string (ip:ip) = ',' 
         CALL get_params (string, ianz, cpara, lpara, maxw, i) 
!                                                                       
!     --Try to build filename                                           
!                                                                       
         CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
         IF (ier_num.eq.0) then 
!                                                                       
!       Try to open the macro file                                      
!                                                                       
            IF (ip.gt.1) then 
!                                                                       
!-----      Try filename as is                                          
!                                                                       
               filename = cpara (1) (1:lpara (1) ) 
               INQUIRE (file = filename, exist = fileda) 
!                                                                       
!-----      Try filename as is with appended '.mac'                     
!                                                                       
               IF (.not.fileda) then 
                  filename = cpara (1) (1:lpara (1) ) //'.mac' 
                  INQUIRE (file = filename, exist = fileda) 
               ENDIF
               filename_length = len_str(filename)
!
!-----      If file is found, convert to absolute path
!
               IF ( fileda ) THEN 
                 lslash = index ( filename , '/' )
                 IF ( lslash == 0 ) THEN       ! No slash in filename
                   filename = ldir(1:ldir_length) // filename(1:filename_length)
                 ENDIF
               ENDIF
!                                                                       
!-----      Try local directory with appended '.mac'                    
!                                                                       
               IF (.not.fileda) then 
                  string = ' ' 
                  lc = 1 
                  CALL do_chdir (string, lc, .false.) 
                  IF (string (lc:lc) .ne.'/') then 
                     lc = lc + 1 
                     string (lc:lc) = '/' 
                  ENDIF 
                  filename = string (1:lc) //cpara (1) (1:lpara (1) ) //&
                  '.mac'                                                
                  INQUIRE (file = filename, exist = fileda) 
!                                                                       
!-----      Try local directory                                         
!                                                                       
                  IF (.not.fileda) then 
                     filename = string (1:lc) //cpara (1) (1:lpara (1) ) 
                     INQUIRE (file = filename, exist = fileda) 
                  ENDIF 
               ENDIF 
!                                                                       
!-----      Try users system directory with appended '.mac'             
!                                                                       
               IF (.not.fileda) then 
                  filename = umac_dir (1:umac_dir_l) //cpara (1)        &
                  (1:lpara (1) ) //'.mac'                               
                  INQUIRE (file = filename, exist = fileda) 
               ENDIF 
!                                                                       
!-----      Try users system directory                                  
!                                                                       
               IF (.not.fileda) then 
                  filename = umac_dir (1:umac_dir_l) //cpara (1)        &
                  (1:lpara (1) )                                        
                  INQUIRE (file = filename, exist = fileda) 
               ENDIF 
!                                                                       
!-----      Try system directory with appended '.mac'                   
!                                                                       
               IF (.not.fileda) then 
                  filename = mac_dir (1:mac_dir_l) //cpara (1) (1:lpara &
                  (1) ) //'.mac'                                        
                  INQUIRE (file = filename, exist = fileda) 
               ENDIF 
!                                                                       
!-----      Try system directory                                        
!                                                                       
               IF (.not.fileda) then 
                  filename = mac_dir (1:mac_dir_l) //cpara (1) (1:lpara &
                  (1) )                                                 
                  INQUIRE (file = filename, exist = fileda) 
               ENDIF 
!                                                                       
               IF (fileda) then 
!                                                                       
!       File could be opened., append to maco-file listing              
!                                                                       
                  mac_level = mac_level + 1 
                  mac_name (mac_level) = filename 
                  mac_line (mac_level) = 0 
!                                                                       
!     ----Try to get parameters                                         
!                                                                       
                  IF (ip + 1.le.i) then 
                     IF (ianz.ge.2) then 
                        IF (ianz - 1.le.MAC_MAX_PARA) then 
!                                                                       
!     ----------Copy parameters to macro parameter list.                
!                                                                       
                           DO i = 2, ianz 
                           mac_para (i - 1, mac_level) = cpara (i) 
                           mac_leng (i - 1, mac_level) = lpara (i) 
                           ENDDO 
                           mac_n_par (mac_level) = ianz - 1 
                           WRITE (mac_para (0, mac_level), 1000) ianz - &
                           1                                            
                           mac_leng (0, mac_level) = 4 
                        ELSE 
                           ier_num = - 1 
                           ier_typ = ER_MAC 
                        ENDIF 
                     ENDIF 
                  ELSE 
                     mac_n_par (mac_level) = 0 
                     mac_para (0, mac_level) = '0' 
                     mac_leng (0, mac_level) = 1 
                  ENDIF 
!                                                                       
!       --if macro file is already active, or currently debugged,       
!         close current file                                            
!                                                                       
                  i = max (1, mac_level - 1) 
                  IF (lmakro.or.lmacro_dbg (i) ) then 
                     CLOSE (19) 
                  ENDIF 
!                                                                       
                  lmakro = .true. 
                  IF (prompt.ne.'macro ') oprompt = prompt 
      prompt = 'macro  ' 
                  OPEN (unit = 19, file = filename, status = 'old') 
               ELSE 
!       MAKRO nicht gefunden                                            
                  ier_num = - 12 
                  ier_typ = ER_MAC 
               ENDIF 
            ELSE 
!       ungueltiges MAKRO                                               
               ier_num = - 13 
               ier_typ = ER_MAC 
            ENDIF 
         ENDIF 
      ELSE 
!     Too many macro levels                                             
         ier_num = - 35 
         ier_typ = ER_MAC 
      ENDIF 
!                                                                       
      IF (ier_num.NE.0) THEN 
         CALL ERRLIST 
         ier_num = - 11 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
 1000 FORMAT    (i4) 
!                                                                       
      END SUBROUTINE file_kdo                       
!*****7*****************************************************************
      SUBROUTINE macro_read (line, laenge) 
!-                                                                      
!     Reads a single line from the current macro file                   
!+                                                                      
      IMPLICIT none 
!                                                                       
      include'charact.inc' 
      include'doact.inc' 
      include'errlist.inc' 
      include'macro.inc' 
      include'prompt.inc' 
!                                                                       
      CHARACTER ( * ) line 
!                                                                       
      CHARACTER(1) cdummy 
      CHARACTER(1024) zeile 
      CHARACTER(1024) string 
      INTEGER ndol 
      INTEGER lpar 
      INTEGER n_par 
      INTEGER laenge 
      INTEGER lx, nx, x, lll, sdol 
      INTEGER i, il 
      INTEGER lstring 
      LOGICAL lnum 
      REAL r_par 
!                                                                       
      INTEGER len_str 
      LOGICAL str_comp 
!                                                                       
      ier_num = 0 
      ier_typ = ER_NONE 
      mac_line (mac_level) = mac_line (mac_level) + 1 
      READ (19, 1000, end = 20) line 
      laenge = len_str (line) 
      IF (laenge.eq.0) then 
         line = ' ' 
         laenge = 0 
         RETURN 
      ENDIF 
!                                                                       
      IF (line (1:1) .ne.'#'.and.line (1:1) .ne.'!'.and.Laenge.ne.0)    &
      then                                                              
!                                                                       
!     If inside a macro, check for parameter substitution               
!                                                                       
         zeile = line 
         ndol = index (zeile (1:laenge) , '$') 
         IF (ndol.gt.laenge) then 
            ndol = 0 
         ENDIF 
!                                                                       
!     Replace all '$'-parameters                                        
!                                                                       
         DO while (ndol.ne.0) 
!                                                                       
!------ --Determine length of the numerical parameter string            
!       i.e.: '1', '12'                                                 
!                                                                       
         lx = 0 
         nx = ndol + lx + 1 
         x = iachar (zeile (nx:nx) ) 
         lnum = zero.le.x.and.x.le.nine 
         DO while (lnum.and.nx.lt.laenge) 
         lx = lx + 1 
         nx = ndol + lx + 1 
         x = iachar (zeile (nx:nx) ) 
         lnum = zero.le.x.and.x.le.nine 
         ENDDO 
         IF (nx.lt.laenge) then 
            nx = nx - 1 
         ELSEIF (nx.eq.laenge.and..not.lnum) then 
            nx = nx - 1 
         ENDIF 
!                                                                       
!     --Read parameter number and substitute rest of string             
!                                                                       
         IF (nx.ge.ndol + 1) then 
            string = zeile (ndol + 1:nx) 
            lstring = nx - (ndol + 1) + 1 
            CALL ber_params (1, string, lstring, r_par, 1) 
         ELSE 
            ier_num = - 12 
            ier_typ = ER_FORT 
         ENDIF 
!DBG        read(zeile(ndol+1:nx),*) n_par                              
         IF (ier_num.eq.0) then 
            n_par = nint (r_par) 
            IF (n_par.le.mac_n_par (mac_level) ) then 
               lpar = mac_leng (n_par, mac_level) 
               line = ' ' 
               IF (ndol.gt.1) then 
                  line (1:ndol - 1) = zeile (1:ndol - 1) 
               ENDIF 
!                                                                       
               line (ndol:ndol + lpar - 1) = mac_para (n_par, mac_level)&
               (1:lpar)                                                 
               lll = ndol + lpar - 1 
               IF (nx.lt.laenge) then 
                  line (ndol + lpar:ndol + lpar + laenge-nx) = zeile (  &
                  nx + 1:laenge)                                        
                  lll = lll + laenge-nx 
               ENDIF 
               zeile = ' ' 
               zeile = line 
               laenge = lll 
               sdol = ndol + 1 
               IF (sdol.gt.laenge) then 
                  ndol = 0 
               ELSE 
                  ndol = index (zeile (sdol:laenge) , '$') 
                  IF (ndol.gt.laenge) then 
                     ndol = 0 
                  ELSEIF (ndol.gt.0) then 
                     ndol = ndol + sdol - 1 
                  ENDIF 
               ENDIF 
            ELSE 
               il = len_str (line) 
               WRITE (output_io, 1000) line (1:il) 
               ier_num = - 41 
               ier_typ = ER_MAC 
               RETURN 
            ENDIF 
         ELSE 
            ier_num = 0 
            ier_typ = ER_NONE 
            sdol = ndol + 1 
            IF (sdol.gt.laenge) then 
               ndol = 0 
            ELSE 
               ndol = index (zeile (sdol:laenge) , '$') 
               IF (ndol.gt.laenge) then 
                  ndol = 0 
               ELSEIF (ndol.gt.0) then 
                  ndol = ndol + sdol - 1 
               ENDIF 
            ENDIF 
         ENDIF 
         ENDDO 
      ENDIF 
!                                                                       
!     line read, return to calling routine                              
!                                                                       
      il = len_str (line) 
      IF (prompt_status.ne.PROMPT_OFF) then 
         IF (il.gt.0) then 
            WRITE (output_io, 1000) line (1:il) 
         ELSE 
            WRITE (output_io, 1000) ' ' 
         ENDIF 
      ENDIF 
!                                                                       
!     Only if the string has length longer than zero                    
!                                                                       
      IF (il.gt.0) then 
!                                                                       
!     Check for 'stop' command, unless we are reading a block structure 
!                                                                       
         IF (str_comp (line, 'stop', 4, il, 4) .and..not.lblock_read)   &
         then                                                           
            WRITE ( *, 2000) char (7) 
            lmakro = .false. 
            lmacro_dbg (mac_level) = .true. 
            line = '#' 
            il = 1 
         ENDIF 
      ENDIF 
      RETURN 
!                                                                       
!     End of macro file, if nested macro file go back to previous       
!     else set lmakro to .false.                                        
!                                                                       
   20 CONTINUE 
      IF (prompt_status.ne.PROMPT_OFF) then 
         WRITE (output_io, 1000) ' ' 
      ENDIF 
!                                                                       
      CLOSE (19) 
      mac_line (mac_level) = 0 
      mac_name (mac_level) = ' ' 
      lmacro_dbg (mac_level) = .false. 
      mac_level = mac_level - 1 
!                                                                       
      IF (mac_level.gt.0) then 
         OPEN (UNIT = 19, FILE = mac_name (mac_level) , STATUS = 'OLD') 
         DO i = 1, mac_line (mac_level) 
         READ (19, 1000, end = 30) cdummy 
         ENDDO 
!                                                                       
!     --This level has been stopped                                     
!                                                                       
         IF (lmacro_dbg (mac_level) ) then 
            lmakro = .false. 
            line = '#' 
            il = 1 
         ENDIF 
      ELSE 
!                                                                       
!     --This was the last macro, reset macro status and return          
!       program prompt                                                  
!                                                                       
         CALL macro_close 
      ENDIF 
!                                                                       
      line = '#' 
      laenge = 1 
!                                                                       
      RETURN 
!                                                                       
   30 CONTINUE 
!                                                                       
!     unexpected EOF in previous macro file                             
!                                                                       
      line = '#' 
      laenge = 1 
      CALL macro_close 
      ier_num = - 36 
      ier_typ = ER_MAC 
!                                                                       
!                                                                       
 1000 FORMAT     (a) 
 2000 FORMAT     (' ------ > Macro halted, continue with cont ...',a1) 
      END SUBROUTINE macro_read                     
!*****7*****************************************************************
      SUBROUTINE macro_close 
!-                                                                      
!     Closes the macro file, switches macro status off and sets the     
!     macro level back to zero.                                         
!+                                                                      
      IMPLICIT none 
!                                                                       
      include'macro.inc' 
      include'prompt.inc' 
!                                                                       
      INTEGER i 
!                                                                       
      CLOSE (19) 
      lmakro = .false. 
      prompt = oprompt 
      mac_level = 0 
      DO i = 1, MAC_MAX_LEVEL 
      mac_line (i) = 0 
      mac_name (i) = ' ' 
      lmacro_dbg (i) = .false. 
      ENDDO 
!                                                                       
      END SUBROUTINE macro_close                    
!*****7*****************************************************************
      SUBROUTINE macro_continue (zeile, lcomm) 
!-                                                                      
!     Continues the macro file, that had been interupted for            
!     debugging purposes                                                
!+                                                                      
      IMPLICIT none 
!                                                                       
      include'doact.inc' 
      include'errlist.inc' 
      include'macro.inc' 
      include'prompt.inc' 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 1) 
      CHARACTER ( * ) zeile 
      CHARACTER(1024) cpara (maxw) 
      INTEGER lpara (maxw) 
      INTEGER ianz, lcomm 
!                                                                       
      INTEGER len_str 
      LOGICAL str_comp 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lcomm) 
!                                                                       
!     --No parameter for macro or block structures                      
!                                                                       
      IF (ianz.eq.0) then 
!                                                                       
!     --Only active inside macros stopped by a previous 'stop' command  
!                                                                       
         IF (mac_level.gt.0) then 
            lmakro = .true. 
            IF (prompt.ne.'macro ') oprompt = prompt 
            prompt = 'macro ' 
            lmacro_dbg (mac_level) = .false. 
         ENDIF 
!                                                                       
!     Only active inside macros stopped by a previous 'stop' command    
!                                                                       
         IF (lblock_dbg) then 
            lblock_dbg = .false. 
            lblock = .true. 
         ENDIF 
!                                                                       
!     --One parameter and '<pname>' command                             
!                                                                       
      ELSEIF (ianz.eq.1) then 
         IF (str_comp (cpara (1), pname, 3, lpara (1), len_str (pname) )&
         ) then                                                         
            CALL macro_close 
            lblock_dbg = .false. 
            lblock = .false. 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
      END SUBROUTINE macro_continue                 
