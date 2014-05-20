!+                                                                      
!     Handling of blockstructures and all math functions                
!     and parameter updating and referencing is done by                 
!     routines in this file.                                            
!                                                                       
!*****7*****************************************************************
      SUBROUTINE do_loop (line, lend, length) 
!+                                                                      
!     All commands included in a block structure are read and stored in 
!     character array. Executable commands are parsed to mach_kdo.      
!     If an illegal nesting of block structures is detected, an         
!     error flag is returned and the block structure is not executed    
!     at all.                                                           
!-                                                                      
      USE doact_mod 
      USE doexec_mod 
      USE doloop_mod 
      USE errlist_mod 
      USE learn_mod 
      USE class_macro_internal 
      USE prompt_mod 
      USE set_sub_generic_mod
!                                                                       
      IMPLICIT none 
!
      CHARACTER(1024) line 
      CHARACTER(1024) zeile 
      CHARACTER(20) prom 
      CHARACTER(4) befehl 
      CHARACTER(3) cprom (0:3) 
      INTEGER jlevel (0:maxlev) 
      INTEGER i, length, lp, lbef 
      LOGICAL lend, lreg 
!                                                                       
      INTEGER len_str 
      LOGICAL str_comp 
!                                                                       
      DATA cprom / '/do', '/if', '/do', '/do' / 
!                                                                       
!-----      read first block structure command                          
!                                                                       
      IF (line (1:2) .eq.'do'.and.index (line, '=') .ne.0) then 
         jlevel (0) = 0 
         i = length - 3 
         CALL rem_bl (line (4:length), i) 
         length = i + 3 
      ELSEIF (line (1:2) .eq.'do'.and.index (line, 'while') .ne.0) then 
         jlevel (0) = 2 
         i = length - 3 
         CALL rem_bl (line (4:length), i) 
         length = i + 3 
      ELSEIF (line (1:2) .eq.'do'.and.length.eq.2) then 
         jlevel (0) = 3 
      ELSEIF (line (1:2) .eq.'if'.and.index (line, 'then') .ne.0) then 
         jlevel (0) = 1 
         CALL rem_bl (line, length) 
      ELSE 
         ier_num = - 31 
         ier_typ = ER_FORT 
         GOTO 999 
      ENDIF 
      prom = prompt (1:len_str (prompt) ) //cprom (jlevel (0) ) 
      DO i = 0, maxlev 
      nlevel (i) = - 1 
      ENDDO 
      level = 0 
      nlevel (level) = 0 
      do_comm (0, 0) = line 
      do_leng (0, 0) = length 
!                                                                       
!.....read all commands                                                 
!                                                                       
      DO while (level.gt. - 1) 
      lblock_read = .true. 
      nlevel (level) = nlevel (level) + 1 
      IF (nlevel (level) .gt.maxcom) then 
         ier_num = - 15 
         ier_typ = ER_FORT 
         GOTO 999 
      ENDIF 
!                                                                       
   10 CONTINUE 
!                                                                       
      prom = prompt (1:len_str (prompt) ) //cprom (jlevel (level) ) 
      CALL get_cmd (line, length, befehl, lbef, zeile, lp, prom) 
      IF (line (1:2) .eq.'do'.and.index (line, '=') .ne.0) then 
         i = length - 3 
         CALL rem_bl (line (4:length), i) 
         length = i + 3 
      ELSEIF (line (1:2) .eq.'do'.and.index (line, 'while') .ne.0) then 
         i = length - 3 
         CALL rem_bl (line (4:length), i) 
         length = i + 3 
      ELSEIF (line (1:2) .eq.'if') then 
         CALL rem_bl (line, length) 
      ENDIF 
!                                                                       
!------ execute a macro file                                            
!                                                                       
      IF (line (1:1) .eq.'@') then 
         ier_num = 0 
         ier_typ = ER_NONE 
         IF (length.ge.2) then 
            CALL file_kdo (line (2:length), length - 1) 
         ELSE 
            ier_num = - 13 
            ier_typ = ER_MAC 
         ENDIF 
         GOTO 10 
      ENDIF 
!                                                                       
                                                                        
      do_comm (nlevel (level), level) = line 
      do_leng (nlevel (level), level) = length 
      ldostart (level) = .true. 
!                                                                       
!     sort commands into proper levels                                  
!                                                                       
      IF (line (1:3) .eq.'do '.or.line (1:2) .eq.'if') then 
         IF (level.lt.maxlev) then 
            do_comm (nlevel (level) , level) = '&&' 
            do_leng (nlevel (level), level) = 6 
            WRITE (do_comm (nlevel (level) , level) (3:6) , '(i4)')     &
            nlevel (level + 1) + 1                                      
            level = level + 1 
            nlevel (level) = nlevel (level) + 1 
            do_comm (nlevel (level), level) = line 
            do_leng (nlevel (level), level) = length 
            IF (line (1:2) .eq.'do'.and.index (line, '=') .ne.0) then 
               jlevel (level) = 0 
            ELSEIF (line(1:2).eq.'do'.and.index(line,'while').ne.0) then
               jlevel (level) = 2 
            ELSEIF (line (1:2) .eq.'do'.and.length.eq.2) then 
               jlevel (level) = 3 
            ELSEIF (line (1:2).eq.'if'.and.index(line,'then').ne.0) then
               jlevel (level) = 1 
            ELSE 
               ier_num = - 31 
               ier_typ = ER_FORT 
               GOTO 999 
            ENDIF 
         ELSE 
            ier_num = - 16 
            ier_typ = ER_FORT 
            GOTO 999 
         ENDIF 
      ELSEIF (line (1:5) .eq.'enddo') then 
         IF (jlevel (level) .eq.0.and.length.eq.5) then 
            level = level - 1 
         ELSEIF (jlevel (level) .eq.2.and.length.eq.5) then 
            level = level - 1 
         ELSEIF (jlevel (level) .eq.3.and.index (line, 'until') .ne.0)  then
            level = level - 1 
         ELSE 
            ier_num = - 19 
            ier_typ = ER_FORT 
            GOTO 999 
         ENDIF 
      ELSEIF (line (1:5) .eq.'endif') then 
         IF (jlevel (level) .eq.1) then 
            level = level - 1 
         ELSE 
            ier_num = - 19 
            ier_typ = ER_FORT 
            GOTO 999 
         ENDIF 
      ENDIF 
      ENDDO 
      lblock_read = .false. 
      ier_num = 0 
      ier_typ = ER_NONE 
!                                                                       
!-----      execute the block structure                                 
!                                                                       
      DO i = 0, maxlev 
      ilevel (i) = - 1 
      ltest (i) = .false. 
      ldostart (i) = .true. 
      ENDDO 
      level = 0 
      lblock = .true. 
!                                                                       
!-----      as long as there are commands and no error ...              
!                                                                       
      DO while (level.gt. - 1.and. (                                    &
      ier_num.eq.0.or.ier_num.ne.0.and.ier_sta.eq.ER_S_LIVE) )          
!                                                                       
!     increment the command array, follow up with block commands        
!                                                                       
      CALL do_execute (lreg, line, length) 
      IF (ier_num.ne.0.and.ier_sta.ne.ER_S_LIVE) then 
         GOTO 999 
      ELSEIF (ier_num.ne.0.and.ier_sta.eq.ER_S_LIVE) then 
         CALL errlist 
      ENDIF 
!                                                                       
!     Regular command, call mach_kdo or file_kdo                        
!                                                                       
      IF (lreg) then 
         IF (.not. (level.eq.0.and.ilevel (level) .eq.0) ) then 
            IF (str_comp (line (1:4) , 'stop', 4, length, 4) ) then 
               WRITE (output_io, 2000) achar (7) 
               lblock_dbg = .true. 
               lblock = .false. 
               line = '#' 
               length = 1 
!                                                                       
!     ----Continuous loop until debug mode is switched off              
!                                                                       
               DO while (lblock_dbg) 
               CALL get_cmd (line, length, befehl, lbef, zeile, lp, prom)
               IF (line (1:1) .eq.'@') then 
                  CALL file_kdo (line (2:length), length - 1) 
               ELSE 
                  CALL p_mache_kdo (line, lend, length) 
               ENDIF 
               IF (ier_num.ne.0.and.ier_sta.ne.ER_S_LIVE) then 
                  GOTO 999 
               ELSEIF (ier_num.ne.0.and.ier_sta.eq.ER_S_LIVE) then 
                  CALL errlist 
               ENDIF 
               ENDDO 
               IF (.not.lblock) then 
                  GOTO 999 
               ENDIF 
               line = '#' 
               length = 1 
            ENDIF 
            IF (line (1:1) .eq.'@') then 
               CALL file_kdo (line (2:length), length - 1) 
            ELSE 
               CALL p_mache_kdo (line, lend, length) 
            ENDIF 
            IF (ier_num.ne.0.and.ier_sta.ne.ER_S_LIVE) then 
               GOTO 999 
            ELSEIF (ier_num.ne.0.and.ier_sta.eq.ER_S_LIVE) then 
               CALL errlist 
            ENDIF 
         ENDIF 
      ENDIF 
      ENDDO 
!                                                                       
  999 CONTINUE 
      lblock = .false. 
      IF (ier_num.ne.0) then 
         WRITE (ier_msg (1), 3000) 
         WRITE (ier_msg (2), 3100) line (1:41) 
         IF(lmakro) THEN
            CALL macro_close
         ENDIF 
      ENDIF 
!                                                                       
!                                                                       
 2000 FORMAT    (a1) 
 3000 FORMAT    ('Erroneous line in block structure') 
 3100 FORMAT    (a41) 
      END SUBROUTINE do_loop                        
!*****7**************************************************************** 
      SUBROUTINE do_execute (lreg, line, laenge) 
!-                                                                      
!     This subroutine increments the commands in the block structure.   
!     Block structure commands are executed in this subroutine,         
!     regular commands are returned to the calling subroutine.          
!+                                                                      
      USE doexec_mod 
      USE doloop_mod 
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 2) 
!                                                                       
      CHARACTER ( * ) line 
      CHARACTER(1024) string, cpara (maxw) 
      INTEGER lpara (maxw) 
      INTEGER i, ianz, laenge 
      INTEGER ikl, ithen 
      LOGICAL lreg 
      REAL werte (maxw) 
!                                                                       
      LOGICAL if_test 
!                                                                       
      lreg = .false. 
!                                                                       
      ilevel (level) = ilevel (level) + 1 
      line = do_comm (ilevel (level), level) 
      laenge = do_leng (ilevel (level), level) 
      IF (line (1:2) .eq.'&&') then 
         level = level + 1 
         READ (line (3:6), * ) jump (level) 
         ilevel (level) = jump (level) 
         line = do_comm (ilevel (level), level) 
         laenge = do_leng (ilevel (level), level) 
      ENDIF 
!                                                                       
!     if do-loop command, evaluate counter                              
!                                                                       
      IF (line (1:3) .eq.'do ') then 
         IF (level.eq.0) then 
            jump (level) = ilevel (level) 
         ENDIF 
         CALL do_do (line, level, laenge) 
         IF (ier_num.ne.0) then 
            RETURN 
         ENDIF 
         IF (ldostart (level) ) then 
!     ----do loop is done                                               
            level = level - 1 
         ENDIF 
!                                                                       
!     enddo command, check for enddo until                              
!                                                                       
      ELSEIF (line (1:5) .eq.'enddo') then 
         CALL do_end (line, level, laenge) 
         IF (ier_num.ne.0) then 
            RETURN 
         ENDIF 
         IF (ldostart (level) ) then 
            ilevel (level) = ilevel (level) + 1 
            level = level - 1 
         ELSE 
            DO while (do_comm (ilevel (level) , level) (1:3) .ne.'do ') 
            ilevel (level) = ilevel (level) - 1 
            ENDDO 
            ilevel (level) = jump (level) - 1 
         ENDIF 
!                                                                       
!     if or elseif command                                              
!                                                                       
      ELSEIF (line (1:2) .eq.'if'.or.line (1:6) .eq.'elseif') then 
         IF (index (line, 'then') .eq.0) then 
            ier_num = - 31 
            ier_typ = ER_FORT 
            RETURN 
         ENDIF 
         IF (ltest (level) ) then 
!........... previous block of the current if-elseif had been true      
            ltest (level) = .false. 
            level = level - 1 
         ELSE 
!...........This is the first if statement thats true                   
            ikl = index (line, '(') 
            ithen = index (line, 'then') - 1 
            string = line (ikl:ithen) 
            ltest (level) = if_test (string, ithen - ikl + 1) 
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
            IF (ier_num.ne.0) goto 999 
            IF (.not.ltest (level) ) then 
!.............If statement is false, ignore the commands that follow    
               DO while(do_comm(ilevel(level)+1,level)(1:4).ne.'else'   &
                   .and.do_comm(ilevel(level)+1,level)(1:5).ne.'endif')
               ilevel (level) = ilevel (level) + 1 
               ENDDO 
            ENDIF 
         ENDIF 
!                                                                       
!     else command                                                      
!                                                                       
      ELSEIF (line (1:5) .eq.'else ') then 
         IF (ltest (level) ) then 
!...........A previous block of the currrent if-elseif had been true    
            ltest (level) = .false. 
            level = level - 1 
         ENDIF 
!                                                                       
!     elseif command                                                    
!                                                                       
      ELSEIF (line.eq.'endif') then 
         ltest (level) = .false. 
         level = level - 1 
!                                                                       
!     break command                                                     
!                                                                       
      ELSEIF (line (1:5) .eq.'break') then 
         CALL get_params (line (6:laenge), ianz, cpara, lpara, maxw, laenge-5)                                                      
         IF (ier_num.eq.0) then 
            IF (ianz.eq.1) then 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               IF (ier_num.eq.0) then 
                  IF (level.ge.nint (werte (1) ) - 1) then 
                     DO i = level, level - nint (werte (1) ) + 1, - 1
                        ldostart (i) = .true. 
                        ltest (i) = .false. 
                     ENDDO 
                     level = level - nint (werte (1) ) 
                  ELSE 
                     ier_num = - 28 
                     ier_typ = ER_FORT 
                  ENDIF 
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_FORT 
            ENDIF 
         ENDIF 
!                                                                       
!     regular command                                                   
!                                                                       
      ELSE 
         lreg = .true. 
      ENDIF 
!                                                                       
  999 CONTINUE 
!                                                                       
      END SUBROUTINE do_execute                     
!*****7**************************************************************** 
      SUBROUTINE do_do (line, level, laenge) 
!-                                                                      
!     reads the 'do' command and evaluates the corresponding counter    
!+                                                                      
      USE doloop_mod 
      USE errlist_mod 
      USE set_sub_generic_mod
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 3) 
!                                                                       
      CHARACTER ( * ) line 
      CHARACTER(1024) zeile, cpara (maxw) 
      INTEGER lpara (maxw) 
      INTEGER ipos, ikp, ianz, level, laenge, lll 
      INTEGER ianz_d, i 
      LOGICAL if_test 
      LOGICAL l_var 
      REAL werte (maxw), wert 
!                                                                       
      ier_num = - 6 
      ier_typ = ER_FORT 
!                                                                       
!     search for argument separator on the do-loop command line         
!                                                                       
      ipos = index (line, '=') 
      ikp = index (line, '[') 
!                                                                       
!     Do-loop of type: do counter = start,end[,increment]               
!                                                                       
      IF (ipos.gt.0) then 
!                                                                       
!       The counter is a user defined variable name ?                   
!                                                                       
         l_var = (ikp.eq.0.or.ipos.lt.ikp) 
         lll = ipos - ikp - 2 
         IF (.not.l_var) then 
            IF (lll.gt.0) then 
               zeile = ' ' 
               zeile (1:lll) = line (ikp + 1:ipos - 2) 
               CALL get_params (zeile, ianz_d, cpara, lpara, maxw, lll) 
               IF (ier_num.eq.0) then 
                  CALL ber_params (ianz_d, cpara, lpara, werte, maxw) 
                  IF (ier_num.ne.0) then 
                     RETURN 
                  ENDIF 
                  IF (ianz_d.ge.1.and.ianz_d.le.2) then 
                     DO i = 1, ianz_d 
                     do_kpara (i) = nint (werte (i) ) 
                     ENDDO 
                  ELSE 
                     ier_num = - 17 
                     ier_typ = ER_FORT 
                     RETURN 
                  ENDIF 
               ELSE 
                  RETURN 
               ENDIF 
            ELSE 
               ier_num = - 14 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
         ENDIF 
         ier_num = 0 
         ier_typ = ER_NONE 
         IF (ldostart (level) ) then 
!.........Do loop needs to be started                                   
            IF (ipos.ne.0) then 
               lll = laenge- (ipos + 1) + 1 
               CALL get_params (line (ipos + 1:laenge), ianz, cpara,    &
                                lpara, maxw, lll)                                        
               IF (ier_num.eq.0.and.ianz.ge.2) then 
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  IF (ier_num.ne.0) then 
                     RETURN 
                  ENDIF 
                  glow (level) = werte (1) 
                  ghigh (level) = werte (2) 
                  IF (ianz.eq.2) then 
                     ginc (level) = 1.0 
                  ELSE 
                     ginc (level) = werte (3) 
                  ENDIF 
                  nloop (level) = max (int ( (ghigh (level) - glow (    &
                  level) ) / ginc (level) ) + 1, 0)                     
                  iloop (level) = 0 
                  ldostart (level) = .false. 
                  ier_num = 0 
                  ier_typ = ER_NONE 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_FORT 
                  GOTO 999 
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_FORT 
               GOTO 999 
            ENDIF 
         ENDIF 
         iloop (level) = iloop (level) + 1 
         IF (iloop (level) .le.nloop (level) ) then 
!.........Do loop is running                                            
            wert = glow (level) + (iloop (level) - 1) * ginc (level) 
            IF (l_var) then 
               CALL upd_variable (line (4:ipos - 1), ipos - 4, wert,    &
               cpara (1), lpara (1) )
            ELSE 
               CALL p_upd_para (line (4:ikp - 1), do_kpara, 1, wert, ianz_d)
            ENDIF 
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
         ELSE 
!.........Do loop is done                                               
            ldostart (level) = .true. 
         ENDIF 
!                                                                       
!     do...until loop                                                   
!                                                                       
      ELSEIF (line (4:4) .eq.'    ') then 
         ier_num = 0 
         ier_typ = ER_NONE 
         ldostart (level) = .false. 
!                                                                       
!     DO WHILE loop                                                     
!                                                                       
      ELSEIF (laenge.gt.4.and.index(line (4:laenge),'while') .ne.0) then
         ier_num = 0 
         ier_typ = ER_NONE 
         ipos = index (line, '(') 
         zeile = line (ipos:laenge) 
         ldostart (level) = .not.if_test (zeile, laenge-ipos + 1) 
         IF (ier_num.ne.0) then 
            RETURN 
         ENDIF 
!                                                                       
!     error, wrong parameter on the do loop command line                
!                                                                       
      ELSE 
         ier_num = - 8 
         ier_typ = ER_COMM 
         ldostart (level) = .true. 
      ENDIF 
!                                                                       
  999 CONTINUE 
!                                                                       
      END SUBROUTINE do_do                          
!*****7**************************************************************** 
      SUBROUTINE do_end (line, level, laenge) 
!-                                                                      
!     reads the 'enddo'. If the 'until' is found, the expression is     
!     evaluated                                                         
!+                                                                      
      USE doloop_mod 
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      CHARACTER ( * ) line 
      CHARACTER(1024) zeile 
      INTEGER laenge, ipos, level 
      LOGICAL if_test 
!                                                                       
      ipos = index (line, 'until') 
      IF (ipos.ne.0) then 
         zeile = line (ipos + 5:laenge) 
         ldostart (level) = if_test (zeile, laenge- (ipos + 5) + 1) 
         IF (ier_num.ne.0) then 
            RETURN 
         ENDIF 
      ENDIF 
      END SUBROUTINE do_end                         
!*****7**************************************************************** 
      LOGICAL FUNCTION if_test (string, laenge) 
!-                                                                      
!     Tests the logical condition given in 'string'                     
!+                                                                      
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER MAXW 
      PARAMETER (MAXW = 20) 
!                                                                       
      CHARACTER ( * ) string 
      CHARACTER(2) comp 
      CHARACTER(1024) zeile, line, oldstr 
      CHARACTER(1024) string1, string2 
      CHARACTER(1024) cpara (MAXW) 
      INTEGER lpara (MAXW) 
      INTEGER ianz 
      INTEGER suche_vor, suche_nach 
      INTEGER laenge, icom, iz1, iz2 
      INTEGER ikl, iklz, ikla, ikla1, ikla2 
      INTEGER ll, i, lcom, inot, lll 
      INTEGER istring1, istring2 
      INTEGER istring1_len 
      INTEGER istring2_len 
      LOGICAL lscr, lscr1 
      LOGICAL lstring1, lstring2 
      REAL werte (MAXW) 
      REAL w1, w2 
      REAL berechne 
!                                                                       
      ier_num = 0 
      ier_typ = ER_NONE 
      oldstr = string 
      CALL rem_bl (string, laenge) 
      if_test = .false. 
      IF (laenge.eq.0.or.string.eq.' '.or.ier_num.ne.0) then 
         CONTINUE 
      ELSE 
         icom = max (index (string, '.lt.') , index (string, '.le.') ,  &
                     index (string, '.gt.') , index (string, '.ge.') ,  &
                     index (string, '.eq.') , index (string, '.ne.') )
         DO while (icom.ne.0) 
!                                                                       
!     --Found an operator, search for numbers before and after          
!                                                                       
         comp = string (icom + 1:icom + 2) 
         lll = icom - 1 
         iz1 = suche_vor (string (1:icom - 1), lll) 
         zeile = '('//string (iz1:icom - 1) //')' 
         ll = icom - 1 - iz1 + 3 
         istring1 = index (zeile, '''') 
         lstring1 = .false. 
         IF (istring1.gt.1) then 
!                                                                       
!     ----found a string variable                                       
!                                                                       
            istring2 = index (zeile (istring1 + 1:ll) , '''') + istring1                                                    
            IF (istring2.eq.istring1) then 
!                                                                       
!     ----Missing second ', first ' is there                            
!                                                                       
               ier_num = - 21 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
            string1 = zeile (istring1 + 1:istring2 - 1) 
            istring1_len = istring2 - istring1 - 1 
            lstring1 = .true. 
         ELSE 
            w1 = berechne (zeile, ll) 
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
         ENDIF 
         lll = (laenge) - (icom + 4) + 1 
         iz2 = suche_nach (string (icom + 4:laenge), lll) 
         zeile = '('//string (icom + 4:icom + 4 + iz2 - 1) //')' 
         ll = icom + 4 + iz2 - 1 - (icom + 4) + 3 
         istring1 = index (zeile, '''') 
         lstring2 = .false. 
         IF (istring1.gt.1) then 
!                                                                       
!     ----found a string variable                                       
!                                                                       
            istring2 = index (zeile (istring1 + 1:ll) , '''') + istring1                                                    
            IF (istring2.eq.istring1) then 
!                                                                       
!     ----Missing second ', first ' is there                            
!                                                                       
               ier_num = - 21 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
            string2 = zeile (istring1 + 1:istring2 - 1) 
            istring2_len = istring2 - istring1 - 1 
            lstring2 = .true. 
         ELSE 
            w2 = berechne (zeile, ll) 
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
         ENDIF 
         IF (lstring1.and.lstring2) then 
            CALL get_params (string1 (1:istring1_len), ianz, cpara,     &
            lpara, maxw, istring1_len)                                  
            CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
            string1      = cpara (1) (1:lpara (1) ) 
            istring1_len = lpara (1) 
            CALL get_params (string2 (1:istring2_len), ianz, cpara,     &
            lpara, maxw, istring2_len)                                  
            CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
            string2      = cpara (1) (1:lpara (1) ) 
            istring2_len = lpara (1) 
!                                                                       
            IF (comp.eq.'eq') then 
               lscr = string1 (1:istring1_len) .eq.string2 (1:          &
               istring2_len)                                            
            ELSEIF (comp.eq.'ne') then 
               lscr = string1 (1:istring1_len) .ne.string2 (1:          &
               istring2_len)                                            
            ELSEIF (comp.eq.'lt') then 
               lscr = string1 (1:istring1_len) .lt.string2 (1:          &
               istring2_len)                                            
            ELSEIF (comp.eq.'le') then 
               lscr = string1 (1:istring1_len) .le.string2 (1:          &
               istring2_len)                                            
            ELSEIF (comp.eq.'gt') then 
               lscr = string1 (1:istring1_len) .gt.string2 (1:          &
               istring2_len)                                            
            ELSEIF (comp.eq.'ge') then 
               lscr = string1 (1:istring1_len) .ge.string2 (1:          &
               istring2_len)                                            
            ENDIF 
         ELSE 
            IF (comp.eq.'lt') then 
               lscr = w1.lt.w2 
            ELSEIF (comp.eq.'le') then 
               lscr = w1.le.w2 
            ELSEIF (comp.eq.'gt') then 
               lscr = w1.gt.w2 
            ELSEIF (comp.eq.'ge') then 
               lscr = w1.ge.w2 
            ELSEIF (comp.eq.'eq') then 
               lscr = w1.eq.w2 
            ELSEIF (comp.eq.'ne') then 
               lscr = w1.NE.w2 
            ENDIF 
         ENDIF 
         lcom = 4 
         CALL ersetz_log (string, iz1, iz2, icom, laenge, lcom, lscr) 
         icom = max (index (string, '.lt.') , index (string, '.le.') ,  &
                     index (string, '.gt.') , index (string, '.ge.') ,  &
                     index (string, '.eq.') , index (string, '.ne.') )
         ENDDO 
         ikla = index (string, '(') 
         DO while (ikla.ne.0) 
         iklz = index (string (ikla + 1:laenge) , ')') + ikla 
         IF (iklz.eq.ikla) then 
            ier_num = - 9 
            ier_typ = ER_FORT 
            if_test = .false. 
            RETURN 
         ENDIF 
         ikla2 = index (string (ikla + 1:iklz) , '(') + ikla 
         ikla1 = ikla 
         DO while (ikla2.lt.iklz.and.ikla2.gt.ikla1) 
         ikla1 = ikla2 
         ikla2 = index (string (ikla1 + 1:iklz) , '(') + ikla1 
         ENDDO 
         ikl = max (ikla1, ikla2) 
         IF (ikl.ne.0) then 
            line = string (ikl + 1:iklz - 1) 
            ll = iklz - ikl - 1 
!                                                                       
!     ----Found a set of corresponding brackets                         
!                                                                       
!         Evaluate any '.not.'                                          
!                                                                       
            inot = index (line, '.not.') 
            DO while (inot.ne.0) 
            READ (line (inot + 5:inot + 5) , '(l1)') lscr 
            lscr = .not.lscr 
            iz1 = max (inot - 1, 1) 
            iz2 = 1 
            lcom = 5 
            CALL ersetz_log (line, iz1, iz2, inot, ll, lcom, lscr) 
            inot = index (line, '.not.') 
            ENDDO 
!                                                                       
!         Evaluate any '.and.'                                          
!                                                                       
            inot = index (line, '.and.') 
            DO while (inot.ne.0) 
            READ (line (inot - 1:inot - 1) , '(l1)') lscr1 
            READ (line (inot + 5:inot + 5) , '(l1)') lscr 
            lscr = lscr1.and.lscr 
            iz1 = inot - 1 
            iz2 = 1 
            lcom = 5 
            CALL ersetz_log (line, iz1, iz2, inot, ll, lcom, lscr) 
            inot = index (line, '.and.') 
            ENDDO 
!                                                                       
!         Evaluate any '.eqv.'                                          
!                                                                       
            inot = index (line, '.eqv.') 
            DO while (inot.ne.0) 
            READ (line (inot - 1:inot - 1) , '(l1)') lscr1 
            READ (line (inot + 5:inot + 5) , '(l1)') lscr 
            lscr = (lscr1.and.lscr) .or. (.not.lscr1.and..not.lscr) 
            iz1 = inot - 1 
            iz2 = 1 
            lcom = 5 
            CALL ersetz_log (line, iz1, iz2, inot, ll, lcom, lscr) 
            inot = index (line, '.eqv.') 
            ENDDO 
!                                                                       
!         Evaluate any '.xor.'                                          
!                                                                       
            inot = index (line, '.xor.') 
            DO while (inot.ne.0) 
            READ (line (inot - 1:inot - 1) , '(l1)') lscr1 
            READ (line (inot + 5:inot + 5) , '(l1)') lscr 
            lscr = (lscr1.and..not.lscr) .or. (.not.lscr1.and.lscr) 
            iz1 = inot - 1 
            iz2 = 1 
            lcom = 5 
            CALL ersetz_log (line, iz1, iz2, inot, ll, lcom, lscr) 
            inot = index (line, '.xor.') 
            ENDDO 
!                                                                       
!         Evaluate any '.or.'                                           
!                                                                       
            inot = index (line, '.or.') 
            DO while (inot.ne.0) 
            READ (line (inot - 1:inot - 1) , '(l1)') lscr1 
            READ (line (inot + 4:inot + 4) , '(l1)') lscr 
            lscr = lscr1.or.lscr 
            iz1 = inot - 1 
            iz2 = 1 
            lcom = 4 
            CALL ersetz_log (line, iz1, iz2, inot, ll, lcom, lscr) 
            inot = index (line, '.or.') 
            ENDDO 
            zeile = ' ' 
            IF (ikl.gt.1) zeile (1:ikl - 1) = string (1:ikl - 1) 
            zeile (ikl:ikl + ll - 1) = line (1:ll) 
            lll = ikl + ll - 1 
            IF (iklz + 1.le.laenge) then 
               zeile (ikl + ll:ikl + ll + laenge-iklz) =  &
                    string (iklz +  1:laenge)
               lll = ikl + ll + laenge-iklz 
            ENDIF 
            string = zeile 
            laenge = lll 
         ENDIF 
         ikla = index (string, '(') 
         ENDDO 
         i = 1 
         DO while (string (i:i) .eq.' ') 
         i = i + 1 
         ENDDO 
         zeile = string (i:laenge) 
         string = zeile 
         IF (string.eq.'T'.or.string.eq.'t') then 
            if_test = .true. 
         ELSEIF (string.eq.'F'.or.string.eq.'f') then 
            if_test = .false. 
         ELSE 
            ier_num = - 18 
            ier_typ = ER_FORT 
            if_test = .false. 
         ENDIF 
      ENDIF 
!                                                                       
      END FUNCTION if_test                          
!*****7**************************************************************** 
      SUBROUTINE ersetz_log (string, iz1, iz2, icom, laenge, lcom, lscr) 
!-                                                                      
!       Replaces the result of a logical operation within the           
!       string                                                          
!                                                                       
!       Version  : 0.0                                                  
!                                                                       
!       Date     : September 1992                                       
!                                                                       
!       Author   : r.n.                                                 
!                                                                       
!       modified :                                                      
!                                                                       
!+                                                                      
      IMPLICIT none 
!                                                                       
      CHARACTER ( * ) string 
      CHARACTER(1024) zeile 
      INTEGER iz1, iz2, icom, laenge, lcom, ll 
      LOGICAL lscr 
!                                                                       
      zeile = ' ' 
      IF (iz1.gt.1) zeile (1:iz1 - 1) = string (1:iz1 - 1) 
      WRITE (zeile (iz1:iz1) , '(L1)') lscr 
      ll = laenge- (lcom + iz2 + icom - iz1) + 1 
      IF (icom + lcom + iz2.le.laenge) then 
         zeile (iz1 + 1:ll) = string (icom + lcom + iz2:laenge) 
      ENDIF 
      string = zeile 
      laenge = ll 
!                                                                       
      END SUBROUTINE ersetz_log                     
!*****7**************************************************************** 
      INTEGER FUNCTION suche_vor (string, i) 
!-                                                                      
!     searches for the last '(' or '.') in the string                   
!+                                                                      
      IMPLICIT none 
!                                                                       
      CHARACTER ( * ) string 
      INTEGER i, level 
      LOGICAL l_hyp 
!                                                                       
      l_hyp = .false. 
      level = 0 
      DO while (i.gt.0.and.level.ge.0) 
      IF (string (i:i) .eq.')'.or.string (i:i) .eq.']') then 
         level = level + 1 
      ELSEIF (string (i:i) .eq.''''.and..not.l_hyp) then 
         level = level + 1 
         l_hyp = .true. 
      ELSEIF (string (i:i) .eq.''''.and.l_hyp) then 
         level = level - 1 
         l_hyp = .false. 
      ELSEIF (string (i:i) .eq.'('.or.string (i:i) .eq.'[') then 
         level = level - 1 
      ELSEIF (level.eq.0.and.string (i:i) .eq.',') then 
         level = - 1 
      ELSEIF (i.gt.4) then 
         IF (string (i - 4:i) .eq.'.not.'.or.string (i - 4:i) .eq.'.and.') then
            level = - 1 
         ELSEIF (i.gt.3) then 
            IF (string (i - 3:i) .eq.'.lt.'.or. &
                string (i - 3:i) .eq.'.le.'.or. &
                string (i - 3:i) .eq.'.gt.'.or. &
                string (i - 3:i) .eq.'.ge.'.or. &
                string (i - 3:i) .eq.'.eq.'.or. &
                string (i - 3:i) .eq.'.ne.'.or. &
                string (i - 3:i) .eq.'.or.'     ) then        
               level = - 1 
            ENDIF 
         ENDIF 
      ELSEIF (i.gt.3) then 
         IF (string (i - 3:i) .eq.'.lt.'.or. &
             string (i - 3:i) .eq.'.le.'.or. &
             string (i - 3:i) .eq.'.gt.'.or. &
             string (i - 3:i) .eq.'.ge.'.or. &
             string (i - 3:i) .eq.'.eq.'.or. &
             string (i - 3:i) .eq.'.ne.'.or. &
             string (i - 3:i) .eq.'.or.'     ) then                
            level = - 1 
         ENDIF 
      ENDIF 
      i = i - 1 
      ENDDO 
      suche_vor = i + 2 
      END FUNCTION suche_vor                        
!*****7**************************************************************** 
      INTEGER FUNCTION suche_nach (string, laenge) 
!-                                                                      
!     searches for the first '(' or '.') in the string                  
!+                                                                      
      IMPLICIT none 
!                                                                       
      CHARACTER ( * ) string 
      INTEGER laenge, i, level 
!                                                                       
      i = 1 
      level = 0 
      DO while (i.le.laenge.and.level.ge.0) 
      IF (string (i:i) .eq.')'.or.string (i:i) .eq.']') then 
         level = level - 1 
      ELSEIF (string (i:i) .eq.'('.or.string (i:i) .eq.'[') then 
         level = level + 1 
      ELSEIF (level.eq.0.and.string (i:i) .eq.',') then 
         level = - 1 
      ELSEIF (i.lt.laenge-4) then 
         IF (string (i:i + 4) .eq.'.not.'.or.string (i:i + 4) .eq.'.and.') then
            level = - 1 
         ELSEIF (i.lt.laenge-3) then 
            IF (string (i:i + 3) .eq.'.lt.'.or.   &
                string (i:i + 3) .eq.'.le.'.or.   &
                string (i:i + 3) .eq.'.gt.'.or.   &
                string (i:i + 3) .eq.'.ge.'.or.   &
                string (i:i + 3) .eq.'.eq.'.or.   &
                string (i:i + 3) .eq.'.ne.'.or.   &
                string (i:i + 3) .eq.'.or.') then        
               level = - 1 
            ENDIF 
         ENDIF 
      ELSEIF (i.lt.laenge-3) then 
         IF (string (i:i + 3) .eq.'.lt.'.or.   &
             string (i:i + 3) .eq.'.le.'.or.   &
             string (i:i + 3) .eq.'.gt.'.or.   &
             string (i:i + 3) .eq.'.ge.'.or.   &
             string (i:i + 3) .eq.'.eq.'.or.   &
             string (i:i + 3) .eq.'.ne.'.or.   &
             string (i:i + 3) .eq.'.or.') then                
            level = - 1 
         ENDIF 
      ENDIF 
      i = i + 1 
      ENDDO 
      suche_nach = i - 2 
      END FUNCTION suche_nach                       
!*****7**************************************************************** 
      INTEGER FUNCTION suche_vor2 (string, i) 
!-                                                                      
!     Searches for the last '(' or '*' or '+' or '-'                    
!     in the string                                                     
!+                                                                      
      IMPLICIT none 
!                                                                       
      CHARACTER ( * ) string 
      LOGICAL lcont 
      INTEGER i, j 
!                                                                       
      suche_vor2 = 0 
      lcont = .true. 
      DO while (i.gt.0.and.lcont) 
      IF (string (i:i) .eq.'('.or.string (i:i) .eq.'*'.or.string (i:i)  .eq.'/') then
!     &                 string(i:i).eq.'/'               .or.           
!     &                (string(i:i).eq.'-' .and. i.eq.1) .or.           
!     &                (string(i:i).eq.'+' .and. i.eq.1)                
         lcont = .false. 
      ELSEIF (i.gt.1) then 
         IF ( (string (i:i) .eq.'-'.or.string (i:i) .eq.'+') .and. &
              (string (i - 1:i - 1) .ne.'E'.and.string (i - 1:i - 1) .ne.'e') ) then                                                         
            lcont = .false. 
         ENDIF 
      ENDIF 
      IF (lcont) then 
         i = i - 1 
      ENDIF 
      ENDDO 
!DBG                                                                    
!DBG      if(j.eq.-1) then                                              
!DBG      do while(i.gt.0       .and. string(i:i).ne.'(' .and.          
!DBG     &   string(i:i).ne.'*'   .and.                                 
!DBG     &   string(i:i).ne.'/'   .and.                                 
!DBG     &   (string(i:i).ne.'-'   .or. (i.gt.1             .and.       
!DBG     &   string(i:i).eq.'-'   .and.                                 
!DBG     &   (string(i-1:i-1).eq.'E'.or.string(i-1:i-1).eq.'e') ) .or.  
!DBG     &   (string(i:i).eq.'-'   .and. i.eq.1)               ) .and.  
!DBG     &   (string(i:i).ne.'+'   .or. (i.gt.1             .and.       
!DBG     &   string(i:i).eq.'+'   .and.                                 
!DBG     &   (string(i-1:i-1).eq.'E'.or.string(i-1:i-1).eq.'e') ) .or.  
!DBG     &   (string(i:i).eq.'+'   .and. i.eq.1                  ) ))   
!DBG        i=i-1                                                       
!DBG      write(*,*) i                                                  
!DBG      ENDDO                                                         
!DBG      endif                                                         
      suche_vor2 = i + 1 
      IF (i.gt.0) then 
         IF (i.gt.1.and. (string (i:i) .eq.'-'.or.string (i:i) .eq.'+') ) then                                                         
            IF (i.eq.1) then 
               suche_vor2 = i 
            ELSE 
               j = i - 1 
               DO while (j.gt.1.and.string (j:j) .eq.' ') 
               j = j - 1 
               ENDDO 
               IF (string (j:j) .eq.'*'.or. &
                   string (j:j) .eq.'/'.or. &
                   string (j:j) .eq.'+'.or. &
                   string (j:j) .eq.'-'.or. &
                   string (j:j) .eq.'(') then
                  suche_vor2 = i 
               ENDIF 
            ENDIF 
         ENDIF 
      ENDIF 
      END FUNCTION suche_vor2                       
!*****7**************************************************************** 
      INTEGER FUNCTION suche_nach2 (string, laenge) 
!-                                                                      
!     Searches for the first '(' or '*' or '/' or '+' or '-'            
!     in the string                                                     
!       im string                                                       
!+                                                                      
      IMPLICIT none 
!                                                                       
      CHARACTER ( * ) string 
      LOGICAL lcont 
      INTEGER laenge 
      INTEGER i 
!                                                                       
      suche_nach2 = 0 
      i = 1 
      IF (laenge.eq.0) return 
      DO while (i.le.laenge.and.string (i:i) .eq.' ') 
      i = i + 1 
      ENDDO 
      IF (string (i:i) .eq.'-'.or.string (i:i) .eq.'+') then 
         i = i + 1 
      ENDIF 
      IF (laenge.gt.1) then 
         lcont = .true. 
         DO while (i.le.laenge.and.lcont) 
            IF (string (i:i) .eq.')'.or.string (i:i) .eq.'*'.or. &
                string (i:i) .eq.'/'.or.                         &
               (string (i:i) .eq.'-'.and.i.eq.1) .or.            &
               (string (i:i) .eq.'+'.and.i.eq.1) ) then                                
            lcont = .false. 
            ELSEIF (i.gt.1) then 
               IF ( (string (i:i) .eq.'-'.or.          &
                     string (i:i) .eq.'+'     )  .and. &
                    (string (i - 1:i - 1) .ne.'E'.and. &
                     string (i - 1:i - 1) .ne.'e') ) then                                             
               lcont = .false. 
               ENDIF 
            ENDIF 
            IF (lcont) then 
               i = i + 1 
            ENDIF 
         ENDDO 
!DBG                                                                    
         IF (i.eq. - 1111) then 
            DO while (i.le.laenge.and.                                &
               string (i:i) .ne.')' .and.  string (i:i) .ne.'*' .and. &
               string (i:i) .ne.'/' .and.                             &
              (string (i:i) .ne.'-' .or.                              &
               (i.gt.1             .and.                              &
                string (i:i) .eq.'-'.and.                             &
               (string (i - 1:i - 1) .eq.'E'.or.                      &
                string (i - 1:i - 1) .eq.'e') ) ) .and.               &                
              (string (i:i) .ne.'+'.or.                               &
               (i.gt.1             .and.                              &
                string (i:i) .eq.'+' .and.                            &
                (string (i - 1:i - 1) .eq.'E'.or.                     &
                 string (i - 1: i - 1) .eq.'e') ) ) )                                       
            i = i + 1 
            ENDDO 
         ENDIF 
         suche_nach2 = i - 1 
      ELSE 
         suche_nach2 = i 
      ENDIF 
      END FUNCTION suche_nach2                      
!****7***************************************************************** 
      SUBROUTINE do_math (line, indxg, length) 
!-                                                                      
!     Calculates the value of an expression and stores the result       
!     in the proper variable                                            
!                                                                       
      USE errlist_mod 
      USE set_sub_generic_mod
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 10) 
!                                                                       
      CHARACTER ( * ) line 
      CHARACTER(1024) zeile, cpara (maxw) 
!                                                                       
      INTEGER lpara (maxw) 
      INTEGER indxg, i, ikk, iii (maxw), ianz, lll 
      INTEGER length 
!                                                                       
      REAL berechne, wert, werte (maxw) 
!                                                                       
!     String substitution???                                            
!                                                                       
      IF (index (line, '"') .gt.0.or.index (line, '''') .gt.0) then 
         CALL do_string_alloc (line, indxg, length) 
         RETURN 
      ENDIF 
!                                                                       
!     Get the expression                                                
!                                                                       
      lll = length - (indxg + 1) + 1 
      CALL get_params (line (indxg + 1:length), ianz, cpara, lpara, maxw, lll)
      IF (ier_num.ne.0) then 
         RETURN 
      ELSEIF (ianz.eq.0) then 
         ier_num = - 6 
         ier_typ = ER_COMM 
         RETURN 
      ENDIF 
      i = lpara (1) 
      zeile = '('//cpara (1) (1:i) //')' 
      i = i + 2 
!                                                                       
!     Calculate the expression                                          
!                                                                       
      wert = berechne (zeile, i) 
      IF (ier_num.eq.0) then 
!                                                                       
!-----evaluate the index of the variable                                
!                                                                       
         lll = indxg - 1 
         CALL get_params (line (1:indxg - 1), ianz, cpara, lpara, maxw, lll)
         IF (ier_num.eq.0) then 
            line = cpara (1) 
            i = lpara (1) 
            ikk = index (line, '[') 
            IF (ikk.lt.i.and.ikk.gt.0) then 
               IF (line (i:i) .eq.']') then 
                  IF (i.gt.ikk + 1) then 
                     zeile = ' ' 
                     zeile (1:i - ikk - 1) = line (ikk + 1:i - 1) 
                     lll = i - ikk - 1 
                     CALL get_params (zeile, ianz, cpara, lpara, maxw, lll)
                     IF (ier_num.eq.0) then 
                        IF (ianz.ge.1.or.ianz.le.3) then 
                           CALL ber_params (ianz, cpara, lpara, werte, maxw)
                           IF (ier_num.eq.0) then 
                              DO i = 1, ianz 
                              iii (i) = nint (werte (i) ) 
                              ENDDO 
!                                                                       
!     ------------Store result in the variable                          
!                                                                       
                              CALL p_upd_para (line (1:ikk - 1), iii,ianz, wert, ianz)
                           ENDIF 
                        ELSE 
                           ier_num = - 6 
                           ier_typ = ER_FORT 
                        ENDIF 
                     ENDIF 
                  ELSE 
                     ier_num = - 10 
                     ier_typ = ER_FORT 
                  ENDIF 
               ELSE 
                  ier_num = - 9 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSEIF (ikk.eq.0) then 
               CALL upd_variable (line (1:i), i, wert, cpara (1), lpara (1) )
            ELSE 
               ier_num = - 6 
               ier_typ = ER_FORT 
            ENDIF 
         ENDIF 
      ELSE 
         zeile = line (1:indxg) //'"%c",'//line (indxg + 1:length) 
         length = length + 5 
         CALL do_string_alloc (zeile, indxg, length) 
      ENDIF 
!                                                                       
      END SUBROUTINE do_math                        
!****7***************************************************************** 
      SUBROUTINE do_string_alloc (line, indxg, length) 
!-                                                                      
!     Evaluates the parameters and stores the result                    
!     in the proper string variable                                     
!                                                                       
      USE errlist_mod 
      USE param_mod 
      USE set_sub_generic_mod
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 10) 
!                                                                       
      CHARACTER ( * ) line 
      CHARACTER(1024) zeile, cpara (maxw) 
      CHARACTER(1024) string 
!                                                                       
      INTEGER lpara (maxw) 
      INTEGER indxg, i, ikk, iii (maxw), ianz, lll 
      INTEGER ising 
      INTEGER length, l_string 
!                                                                       
      REAL wert, werte (maxw) 
!                                                                       
!     for flexibility                                                   
!                                                                       
      ising = index (line, '''') 
      DO while (ising.gt.0) 
      line (ising:ising) = '"' 
      ising = index (line, '''') 
      ENDDO 
!                                                                       
!     Get the expression                                                
!                                                                       
      lll = - (length - (indxg + 1) + 1) 
      CALL get_params (line (indxg + 1:length), ianz, cpara, lpara, maxw, lll)
      IF (ier_num.ne.0) then 
         RETURN 
      ELSEIF (ianz.eq.0) then 
         ier_num = - 6 
         ier_typ = ER_COMM 
         RETURN 
      ENDIF 
!                                                                       
!     Construct the string                                              
!                                                                       
      CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
      string = cpara (1) (1:lpara (1) ) 
      l_string = lpara (1) 
      res_para (0) = 1 
      res_para (1) = lpara (1) 
      IF (ier_num.eq.0) then 
!                                                                       
!-----evaluate the index of the variable                                
!                                                                       
         lll = indxg - 1 
         CALL get_params (line (1:indxg - 1), ianz, cpara, lpara, maxw, lll)
         IF (ier_num.eq.0) then 
            line = cpara (1) 
            i = lpara (1) 
            ikk = index (line, '[') 
            IF (ikk.lt.i.and.ikk.gt.0) then 
               IF (line (i:i) .eq.']') then 
                  IF (i.gt.ikk + 1) then 
                     zeile = ' ' 
                     zeile (1:i - ikk - 1) = line (ikk + 1:i - 1) 
                     lll = i - ikk - 1 
                     CALL get_params (zeile, ianz, cpara, lpara, maxw, lll)
                     IF (ier_num.eq.0) then 
                        IF (ianz.ge.1.or.ianz.le.3) then 
                           CALL ber_params (ianz, cpara, lpara, werte, maxw)
                           IF (ier_num.eq.0) then 
                              DO i = 1, ianz 
                              iii (i) = nint (werte (i) ) 
                              ENDDO 
!                                                                       
!     ------------Store result in the variable                          
!                                                                       
                              CALL p_upd_para (line (1:ikk - 1), iii, ianz, wert, ianz)
                           ENDIF 
                        ELSE 
                           ier_num = - 6 
                           ier_typ = ER_FORT 
                        ENDIF 
                     ENDIF 
                  ELSE 
                     ier_num = - 10 
                     ier_typ = ER_FORT 
                  ENDIF 
               ELSE 
                  ier_num = - 9 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSEIF (ikk.eq.0) then 
               CALL upd_variable (line (1:i), i, wert, string, l_string) 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_FORT 
            ENDIF 
         ENDIF 
      ENDIF 
!                                                                       
      END SUBROUTINE do_string_alloc                
!****7***************************************************************** 
      SUBROUTINE do_eval (line, i) 
!-                                                                      
!     evaluates the expression stored in line                           
!                                                                       
      USE errlist_mod 
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 10) 
!                                                                       
      CHARACTER(1024) line
      CHARACTER(1024) cpara (maxw) 
      CHARACTER(1024) cstr 
      INTEGER lpara (maxw) 
      INTEGER i, ianz, il 
      INTEGER length 
      REAL werte (maxw) 
!                                                                       
      INTEGER len_str 
!                                                                       
      IF (line.eq.' ') then 
         ier_num = - 6 
         ier_typ = ER_COMM 
         RETURN 
      ELSE 
!                                                                       
!     --Non blank line. Make length negative to avoid removing blanks   
!                                                                       
         i = - i 
         CALL get_params (line, ianz, cpara, lpara, maxw, i) 
         IF (ier_num.eq.0) then 
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
            IF (ier_num.eq.0) then 
               length = lpara (1) 
               DO i = 1, ianz 
               length = max (length, lpara (i) ) 
               ENDDO 
               DO i = 1, ianz 
               WRITE ( *, 2222) cpara (i) (1:length), werte (i) 
               IF (output_status.eq.OUTPUT_FILE) then 
                  WRITE (output_io, 2222) cpara (i) (1:length), werte ( i)
               ENDIF 
               IF (lconn.and.lsocket.and.i.eq.1) then 
                  WRITE (cstr, 2222) cpara (i) (1:lpara (i) ), werte (i) 
                  il = len_str (cstr) 
                  CALL socket_send (s_conid, cstr, il) 
               ENDIF 
               ENDDO 
            ENDIF 
         ENDIF 
      ENDIF 
!                                                                       
 2222 FORMAT    (' Value of ',a,' = ',g15.8) 
      END SUBROUTINE do_eval                        
!****7***************************************************************** 
      REAL FUNCTION berechne (string, laenge) 
!-                                                                      
!     Calculates the value of the expression stored in string           
!+                                                                      
      USE charact_mod
      USE errlist_mod 
      USE set_sub_generic_mod
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 3) 
!                                                                       
      CHARACTER ( * ) string 
      CHARACTER(1024) zeile, line, cpara (maxw) 
      INTEGER lpara (maxw) 
      INTEGER max 
      INTEGER laenge, ikla, iklz, ikla1, ikla2, ikl, ll, lll, ie 
      INTEGER ikpa, ikpa1, ikpa2, ikp, ikpz, lp, ianz, i, ikom 
      REAL werte (maxw) 
      REAL r 
!                                                                       
      REAL do_read_number 
!                                                                       
      ier_num = 0 
      ier_typ = ER_NONE 
      berechne = 0 
      IF (laenge.eq.0.or.string.eq.' '.or.ier_num.ne.0) then 
         CONTINUE 
      ELSE 
         CALL ersetz_variable (string, laenge) 
         ie = index (string, 'E') 
         DO while (ie.ne.0) 
         string (ie:ie) = 'e' 
         ie = index (string, 'E') 
         ENDDO 
         ikla = index (string, '(') 
         DO while (ikla.ne.0) 
         iklz = index (string (ikla + 1:laenge) , ')') + ikla 
         IF (iklz.gt.ikla + 1) then 
            ikla2 = index (string (ikla + 1:iklz) , '(') + ikla 
         ELSE 
            ikla2 = ikla 
         ENDIF 
         ikla1 = ikla 
         DO while (ikla2.lt.iklz.and.ikla2.gt.ikla1) 
         ikla1 = ikla2 
         IF (iklz.gt.ikla + 1) then 
            ikla2 = index (string (ikla1 + 1:iklz) , '(') + ikla1 
         ELSE 
            ikla2 = ikla1 
         ENDIF 
         ENDDO 
         ikl = max (ikla1, ikla2) 
         IF (ikl.ne.0) then 
!                                                                       
!     ----Print error messages if Parentheses are missing               
!                                                                       
            IF (iklz.eq.ikl) then 
               ier_num = - 11 
               ier_typ = ER_FORT 
               RETURN 
            ELSEIF (iklz.eq.ikl + 1) then 
               ier_num = - 12 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
            line = string (ikl + 1:iklz - 1) 
            ll = iklz - ikl - 1 
!                                                                       
!     ----Found a pair of parentheses, search for inner set             
!                                                                       
            ikpa = index (line, '[') 
            DO while (ikpa.ne.0) 
            ikpz = index (line (ikpa + 1:ll) , ']') + ikpa 
            IF (ikpz.gt.ikpa + 1) then 
               ikpa2 = index (line (ikpa + 1:ikpz) , '[') + ikpa 
            ELSE 
               ikpa2 = ikpa 
            ENDIF 
            ikpa1 = ikpa 
            DO while (ikpa2.lt.ikpz.and.ikpa2.gt.ikpa1) 
            ikpa1 = ikpa2 
            IF (ikpz.gt.ikpa + 1) then 
               ikpa2 = index (line (ikpa1 + 1:ikpz) , '[') + ikpa1 
            ELSE 
               ikpa2 = ikpa1 
            ENDIF 
            ENDDO 
            ikp = max (ikpa1, ikpa2) 
            IF (ikp.ne.0) then 
               IF (ikpz.gt.ikp + 1) then 
                  zeile = line (ikp + 1:ikpz - 1) 
                  lp = ikpz - ikp - 1 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     DO i = 1, ianz 
                     CALL eval (cpara (i), lpara (i) ) 
                     IF (ier_num.ne.0) then 
                        RETURN 
                     ENDIF 
                     werte (i) = do_read_number (cpara (i), lpara (i) ) 
                     IF (ier_num.ne.0) then 
                        GOTO 999 
                     ENDIF 
                     ENDDO 
                     CALL p_ersetz_para (ikp, ikpz, line, ll, werte, maxw, ianz)
                     IF (ier_num.ne.0) then 
                        RETURN 
                     ENDIF 
                  ELSE 
                     RETURN 
                  ENDIF 
               ELSEIF (ikpz.eq.ikp + 1) then 
                  ier_num = - 10 
                  ier_typ = ER_FORT 
                  RETURN 
               ELSE 
                  ier_num = - 9 
                  ier_typ = ER_FORT 
                  RETURN 
               ENDIF 
            ENDIF 
            ikpa = index (line, '[') 
            ENDDO 
!                                                                       
!                                                                       
            ikom = index (line, ',') 
            IF (ikom.eq.0) then 
               CALL eval (line, ll) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
            ENDIF 
            IF (ikl.ge.3) then 
               CALL calc_intr (string, line, ikl, iklz, laenge, ll) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
            ELSE 
               zeile = ' ' 
               IF (ikl.gt.1) zeile (1:ikl - 1) = string (1:ikl - 1) 
               zeile (ikl:ikl + ll - 1) = line (1:ll) 
               lll = ikl + ll - 1 
               IF (iklz + 1.le.laenge) then 
                  zeile (ikl + ll:ikl + ll + laenge-iklz - 1) =  &
                       string (iklz + 1:laenge)
                  lll = ikl + ll + laenge-iklz - 1 
               ENDIF 
               string = zeile 
               laenge = lll 
            ENDIF 
         ENDIF 
         ikla = index (string, '(') 
         ENDDO 
      ENDIF 
      r = do_read_number (string, laenge) 
      berechne = do_read_number (string, laenge) 
  999 CONTINUE 
!                                                                       
      END FUNCTION berechne                         
!****7***************************************************************** 
      SUBROUTINE berechne_char (string, laenge) 
!-                                                                      
!     Calculates the value of the character expression stored in string 
!+                                                                      
      USE charact_mod
      USE errlist_mod 
      USE set_sub_generic_mod
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 3) 
!                                                                       
      CHARACTER ( * ) string 
      CHARACTER(1024) zeile, line, cpara (maxw) 
      CHARACTER(1024) substring 
      INTEGER lpara (maxw) 
      INTEGER max 
      INTEGER laenge, ikla, iklz, ikla1, ikla2, ikl, ll, lll
      INTEGER ikpa, ikpa1, ikpa2, ikp, ikpz, lp, ianz, i, ikom 
      INTEGER lsub 
      INTEGER icol 
      INTEGER iapo 
      INTEGER j (2) 
      REAL werte (maxw) 
!                                                                       
      REAL berechne 
      REAL do_read_number 
!                                                                       
      lll = 1
      ier_num = 0 
      ier_typ = ER_NONE 
      IF (laenge.eq.0.or.string.eq.' '.or.ier_num.ne.0) then 
         CONTINUE 
      ELSE 
         CALL ersetz_variable (string, laenge) 
         ikla = index (string, '(') 
         DO while (ikla.ne.0) 
         iklz = index (string (ikla + 1:laenge) , ')') + ikla 
         IF (iklz.gt.ikla + 1) then 
            ikla2 = index (string (ikla + 1:iklz) , '(') + ikla 
         ELSE 
            ikla2 = ikla 
         ENDIF 
         ikla1 = ikla 
         DO while (ikla2.lt.iklz.and.ikla2.gt.ikla1) 
         ikla1 = ikla2 
         IF (iklz.gt.ikla + 1) then 
            ikla2 = index (string (ikla1 + 1:iklz) , '(') + ikla1 
         ELSE 
            ikla2 = ikla1 
         ENDIF 
         ENDDO 
         ikl = max (ikla1, ikla2) 
         IF (ikl.ne.0) then 
!                                                                       
!     ----Print error messages if Parentheses are missing               
!                                                                       
            IF (iklz.eq.ikl) then 
               ier_num = - 11 
               ier_typ = ER_FORT 
               RETURN 
            ELSEIF (iklz.eq.ikl + 1) then 
               ier_num = - 12 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
            line = string (ikl + 1:iklz - 1) 
            ll = iklz - ikl - 1 
!                                                                       
!     ----Found a pair of parentheses, search for inner set             
!                                                                       
            ikpa = index (line, '[') 
            DO while (ikpa.ne.0) 
            ikpz = index (line (ikpa + 1:ll) , ']') + ikpa 
            IF (ikpz.gt.ikpa + 1) then 
               ikpa2 = index (line (ikpa + 1:ikpz) , '[') + ikpa 
            ELSE 
               ikpa2 = ikpa 
            ENDIF 
            ikpa1 = ikpa 
            DO while (ikpa2.lt.ikpz.and.ikpa2.gt.ikpa1) 
            ikpa1 = ikpa2 
            IF (ikpz.gt.ikpa + 1) then 
               ikpa2 = index (line (ikpa1 + 1:ikpz) , '[') + ikpa1 
            ELSE 
               ikpa2 = ikpa1 
            ENDIF 
            ENDDO 
            ikp = max (ikpa1, ikpa2) 
            IF (ikp.ne.0) then 
               IF (ikpz.gt.ikp + 1) then 
                  zeile = line (ikp + 1:ikpz - 1) 
                  lp = ikpz - ikp - 1 
                  CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                  IF (ier_num.eq.0) then 
                     DO i = 1, ianz 
                     CALL eval (cpara (i), lpara (i) ) 
                     IF (ier_num.ne.0) then 
                        RETURN 
                     ENDIF 
                     werte (i) = do_read_number (cpara (i), lpara (i) ) 
                     IF (ier_num.ne.0) then 
                        GOTO 999 
                     ENDIF 
                     ENDDO 
                     CALL p_ersetz_para (ikp, ikpz, line, ll, werte, maxw, ianz)
                     IF (ier_num.ne.0) then 
                        RETURN 
                     ENDIF 
                  ELSE 
                     RETURN 
                  ENDIF 
               ELSEIF (ikpz.eq.ikp + 1) then 
                  ier_num = - 10 
                  ier_typ = ER_FORT 
                  RETURN 
               ELSE 
                  ier_num = - 9 
                  ier_typ = ER_FORT 
                  RETURN 
               ENDIF 
            ENDIF 
            ikpa = index (line, '[') 
            ENDDO 
!                                                                       
!                                                                       
            ikom = index (line, ',') 
            icol = index (line, ':') 
            iapo = index (line, '''') 
            IF (ikom.eq.0.and.icol.eq.0.and.iapo.eq.0) then 
               CALL eval (line, ll) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
            ENDIF 
            IF (ikl.ge.3.and.icol.eq.0) then 
               CALL calc_intr (string, line, ikl, iklz, laenge, ll) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
            ELSE 
               IF (icol.ge.1.and.ikl.gt.1) then 
                  cpara (1) = line (1:icol - 1) 
                  lpara (1) = icol - 1 
                  cpara (2) = line (icol + 1:ll) 
                  lpara (2) = ll - icol 
                  ianz = 2 
                  DO i = 1, ianz 
                  substring = '('//cpara (i) (1:lpara (i) ) //')' 
                  lsub = lpara (i) + 2 
                  j (i) = nint (berechne (substring, lsub) ) 
                  IF (ier_num.ne.0) then 
                     RETURN 
                  ENDIF 
                  ENDDO 
                  IF (ikl.gt.1) then 
                     IF (ikl.gt.2) then 
                        IF ( (string (2:2) .eq.''''.or.             &
                              string (2:2) .eq.'"'      ) .and.     &
                             (string (ikl - 1:ikl - 1) .eq.''''.or. &
                              string (ikl - 1:ikl - 1) .eq.'"') ) then
                           j (1) = j (1) + 1 
                           j (2) = j (2) + 1 
                        ENDIF 
                     ENDIF 
                     IF (j (1) .le.j(2) .and.1.le.j (1) .and.      &
                         j (1) .lt.ikl  .and.1.le.j (2) .and.      &
                         j (2) .lt.ikl                       ) then    
                        IF (string (1:1) .eq.'(') then 
                           zeile = string (j (1) + 1:j (2) + 1) 
                           lll = j (2) - j (1) + 1 
                        ELSE 
                           zeile = string (j (1) :j (2) ) 
                           lll = j (2) - j (1) + 1 
                        ENDIF 
                     ELSE 
                        ier_num = - 29 
                        ier_typ = ER_FORT 
                        string = ' ' 
                        laenge = 1 
                        RETURN 
                     ENDIF 
                  ENDIF 
                  string = zeile 
                  laenge = lll 
               ELSE 
                  zeile = ' ' 
                  IF (ikl.gt.1) zeile (1:ikl - 1) = string (1:ikl - 1) 
                  zeile (ikl:ikl + ll - 1) = line (1:ll) 
                  lll = ikl + ll - 1 
                  IF (iklz + 1.le.laenge) then 
                     zeile (ikl + ll:ikl + ll + laenge-iklz - 1)        &
                     = string (iklz + 1:laenge)                         
                     lll = ikl + ll + laenge-iklz - 1 
                  ENDIF 
                  string = zeile 
                  laenge = lll 
               ENDIF 
            ENDIF 
         ENDIF 
         ikla = index (string, '(') 
         ENDDO 
      ENDIF 
  999 CONTINUE 
!                                                                       
      END SUBROUTINE berechne_char                  
!*****7**************************************************************** 
      SUBROUTINE eval (line, ll) 
!-                                                                      
!       evaluates a line that has only the basic arithmetics            
!+                                                                      
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      CHARACTER(1024) line 
      INTEGER ipot, iz1, iz2, ll, imal, idiv, iverk 
      INTEGER iexpo, iplus, iminus, ipp 
      LOGICAL lreal 
      REAL w1, w2, ww 
!                                                                       
!...........Evaluate the exponentiation  '**'                           
      ipot = index (line, '**') 
      DO while (ipot.ne.0) 
      CALL get_w1_w2 (w1, w2, line, ipot, iz1, iz2, ll, 2, lreal) 
      IF (ier_num.ne.0) then 
         RETURN 
      ENDIF 
      ww = w1**w2 
      CALL ersetz (line, iz1, iz2, ww, ipot, ll, 2, lreal) 
      ipot = index (line, '**') 
      ENDDO 
!...... Multiplication, division                                        
      imal = index (line, '*') 
      idiv = index (line, '/') 
      IF (imal.gt.0.and.idiv.gt.0) then 
         iverk = min (imal, idiv) 
      ELSEIF (imal.gt.0.and.idiv.eq.0) then 
         iverk = imal 
      ELSEIF (imal.eq.0.and.idiv.gt.0) then 
         iverk = idiv 
      ELSEIF (imal.eq.0.and.idiv.eq.0) then 
         iverk = 0 
      ENDIF 
      DO while (iverk.ne.0) 
      CALL get_w1_w2 (w1, w2, line, iverk, iz1, iz2, ll, 1, lreal) 
      IF (ier_num.ne.0) then 
         RETURN 
      ENDIF 
      IF (line (iverk:iverk) .eq.'*') then 
         ww = w1 * w2 
      ELSEIF (line (iverk:iverk) .eq.'/') then 
         IF (w2.ne.0) then 
            IF (lreal) then 
               ww = w1 / w2 
            ELSE 
               ww = float (int (w1) / int (w2) ) 
            ENDIF 
         ELSE 
            ier_num = - 4 
            ier_typ = ER_FORT 
            RETURN 
         ENDIF 
      ENDIF 
      CALL ersetz (line, iz1, iz2, ww, iverk, ll, 1, lreal) 
      imal = index (line, '*') 
      idiv = index (line, '/') 
      IF (imal.gt.0.and.idiv.gt.0) then 
         iverk = min (imal, idiv) 
      ELSEIF (imal.gt.0.and.idiv.eq.0) then 
         iverk = imal 
      ELSEIF (imal.eq.0.and.idiv.gt.0) then 
         iverk = idiv 
      ELSEIF (imal.eq.0.and.idiv.eq.0) then 
         iverk = 0 
      ENDIF 
      ENDDO 
!...... addition,subtraction                                            
!...........search for all '+' that are not part of  xxxE+yyy           
      IF (line (1:1) .eq.'+') then 
         iplus = index (line (2:ll) , '+') 
         IF (iplus.ne.0) iplus = iplus + 1 
      ELSE 
         iplus = index (line, '+') 
      ENDIF 
      IF (iplus.gt.1) then 
         iexpo = index (line, 'e') 
         DO while (iplus.gt.1.and.iexpo.gt.1.and.iexpo + 1.eq.iplus) 
         ipp = index (line (iplus + 1:ll) , '+') 
         IF (ipp.eq.0) then 
            iplus = 0 
         ELSE 
            iplus = iplus + ipp 
            iexpo = iexpo + index (line (iexpo + 1:ll) , 'e') 
         ENDIF 
         ENDDO 
      ENDIF 
      IF (line (1:1) .eq.'-') then 
         iminus = index (line (2:ll) , '-') 
         IF (iminus.ne.0) iminus = iminus + 1 
      ELSE 
         iminus = index (line, '-') 
      ENDIF 
      IF (iminus.gt.1) then 
         iexpo = index (line, 'e') 
         DO while (iminus.gt.1.and.iexpo.gt.1.and.iexpo + 1.eq.iminus) 
         ipp = index (line (iminus + 1:ll) , '-') 
         IF (ipp.eq.0) then 
            iminus = 0 
         ELSE 
            iminus = iminus + ipp 
            iexpo = iexpo + index (line (iexpo + 1:ll) , 'e') 
         ENDIF 
         ENDDO 
      ENDIF 
      IF (iplus.gt.0.and.iminus.gt.0) then 
         iverk = min (iplus, iminus) 
      ELSEIF (iplus.gt.0.and.iminus.eq.0) then 
         iverk = iplus 
      ELSEIF (iplus.eq.0.and.iminus.gt.0) then 
         iverk = iminus 
      ELSEIF (iplus.eq.0.and.iminus.eq.0) then 
         iverk = 0 
      ENDIF 
      DO while (iverk.gt.1) 
      CALL get_w1_w2 (w1, w2, line, iverk, iz1, iz2, ll, 1, lreal) 
      IF (ier_num.ne.0) then 
         RETURN 
      ENDIF 
      IF (line (iverk:iverk) .eq.'+') then 
         ww = w1 + w2 
      ELSE 
         ww = w1 - w2 
      ENDIF 
      CALL ersetz (line, iz1, iz2, ww, iverk, ll, 1, lreal) 
      IF (line (1:1) .eq.'+') then 
         iplus = index (line (2:ll) , '+') 
         IF (iplus.ne.0) iplus = iplus + 1 
      ELSE 
         iplus = index (line, '+') 
      ENDIF 
      IF (iplus.gt.1) then 
         iexpo = index (line, 'e') 
         DO while (iplus.gt.1.and.iexpo.gt.1.and.iexpo + 1.eq.iplus) 
         ipp = index (line (iplus + 1:ll) , '+') 
         IF (ipp.eq.0) then 
            iplus = 0 
         ELSE 
            iplus = iplus + ipp 
            iexpo = iexpo + index (line (iexpo + 1:ll) , 'e') 
         ENDIF 
         ENDDO 
      ENDIF 
      IF (line (1:1) .eq.'-') then 
         iminus = index (line (2:ll) , '-') 
         IF (iminus.ne.0) iminus = iminus + 1 
      ELSE 
         iminus = index (line, '-') 
      ENDIF 
      IF (iminus.gt.1) then 
         iexpo = index (line, 'e') 
         DO while (iminus.gt.1.and.iexpo.gt.1.and.iexpo + 1.eq.iminus) 
         ipp = index (line (iminus + 1:ll) , '-') 
         IF (ipp.eq.0) then 
            iminus = 0 
         ELSE 
            iminus = iminus + ipp 
            iexpo = iexpo + index (line (iexpo + 1:ll) , 'e') 
         ENDIF 
         ENDDO 
      ENDIF 
      IF (iplus.gt.0.and.iminus.gt.0) then 
         iverk = min (iplus, iminus) 
      ELSEIF (iplus.gt.0.and.iminus.eq.0) then 
         iverk = iplus 
      ELSEIF (iplus.eq.0.and.iminus.gt.0) then 
         iverk = iminus 
      ELSEIF (iplus.eq.0.and.iminus.eq.0) then 
         iverk = 0 
      ENDIF 
      ENDDO 
      END SUBROUTINE eval                           
!*****7**************************************************************** 
      SUBROUTINE calc_intr (string, line, ikl, iklz, lll, lp) 
!                                                                       
!     Evaluate all intrinsic functions                                  
!                                                                       
      USE errlist_mod 
      USE random_mod
      USE wink_mod
      USE times_mod
      USE set_sub_generic_mod
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 9) 
!                                                                       
      CHARACTER ( * ) string, line 
      CHARACTER(1024) cpara (maxw) 
      CHARACTER(1024) zeile 
      CHARACTER(1024) answer 
      INTEGER lpara (maxw) 
      INTEGER lp, lll, ikom, ikl, iklz, i, ianz 
      INTEGER lcom 
      INTEGER ihyp 
      INTEGER dummy 
      REAL fl1, fl2, fl3, gbox_k, gbox_w, gbox_x 
      REAL werte (maxw), ww, a 
      REAL sind, cosd, tand, asind, acosd, atand 
      REAL atan2, atan2d 
!                                                                       
      INTEGER length_com 
      INTEGER len_str 
      REAL gasdev, ran1, poidev 
      REAL do_read_number 
!                                                                       
      ier_num = - 1 
      ier_typ = ER_FORT 
      ikom = index (line, ',') 
      ihyp = max (index (line, '''') , index (line, '"') ) 
      IF (ikom.eq.0.and.ihyp.eq.0) then 
         ww = do_read_number (line, lp) 
         IF (ier_num.ne.0) then 
            RETURN 
         ENDIF 
      ENDIF 
      ier_num = 0 
      ier_typ = ER_NONE 
      lcom = length_com (string, ikl) 
      IF (lcom.eq.0) then 
         CALL ersetz2 (string, ikl, iklz, ww, 0, lll) 
      ELSEIF (lcom.eq.6) then 
         IF (string (ikl - 6:ikl - 1) .eq.'getcwd') then 
            CALL holecwd (zeile, dummy) 
            i = len_str (zeile) 
            CALL ersetzc (string, ikl, iklz, zeile, i, 6, lll) 
         ELSEIF (string (ikl - 6:ikl - 1) .eq.'getenv') then 
            zeile = line (2:lp - 1) 
            CALL holeenv (zeile, answer) 
            i = len_str (answer) 
            CALL ersetzc (string, ikl, iklz, answer, i, 6, lll) 
         ELSE 
            CALL p_calc_intr_spec (string, line, ikl, iklz, ww, lll, lp) 
         ENDIF 
      ELSEIF (lcom.eq.5) then 
         IF (string (ikl - 5:ikl - 1) .eq.'asind') then 
            IF (abs (ww) .le.1.0) then 
               ww = asind (ww) 
               CALL ersetz2 (string, ikl, iklz, ww, 5, lll) 
            ELSE 
               ier_num = - 7 
               ier_typ = ER_FORT 
            ENDIF 
         ELSEIF (string (ikl - 5:ikl - 1) .eq.'acosd') then 
            IF (abs (ww) .le.1.0) then 
               ww = acosd (ww) 
               CALL ersetz2 (string, ikl, iklz, ww, 5, lll) 
            ELSE 
               ier_num = - 7 
               ier_typ = ER_FORT 
            ENDIF 
         ELSEIF (string (ikl - 5:ikl - 1) .eq.'atand') then 
            CALL get_params (line, ianz, cpara, lpara, 2, lp) 
            IF (ianz.eq.1) then 
               ww = atand (ww) 
            ELSEIF (ianz.eq.2) then 
               werte (1) = do_read_number (cpara (1), lpara (1) ) 
               werte (2) = do_read_number (cpara (2), lpara (2) ) 
               ww = atan2d (werte (1), werte (2) ) 
            ELSE 
               ier_num = - 27 
               ier_typ = ER_FORT 
            ENDIF 
            CALL ersetz2 (string, ikl, iklz, ww, 5, lll) 
         ELSEIF (string (ikl - 5:ikl - 1) .eq.'fmodt') then 
            CALL ersetzc (string, ikl, iklz, f_modt, 24, 5, lll) 
         ELSEIF (string (ikl - 5:ikl - 1) .eq.'fdate') then 
            CALL datum 
            CALL ersetzc (string, ikl, iklz, f_date, 24, 5, lll) 
         ELSE 
            CALL p_calc_intr_spec (string, line, ikl, iklz, ww, lll, lp) 
         ENDIF 
      ELSEIF (lcom.eq.4) then 
         IF (string (ikl - 4:ikl - 1) .eq.'asin') then 
            IF (abs (ww) .le.1.0) then 
               ww = asin (ww) 
               CALL ersetz2 (string, ikl, iklz, ww, 4, lll) 
            ELSE 
               ier_num = - 7 
               ier_typ = ER_FORT 
            ENDIF 
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'acos') then 
            IF (abs (ww) .le.1.0) then 
               ww = acos (ww) 
               CALL ersetz2 (string, ikl, iklz, ww, 4, lll) 
            ELSE 
               ier_num = - 7 
               ier_typ = ER_FORT 
            ENDIF 
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'atan') then 
            CALL get_params (line, ianz, cpara, lpara, 2, lp) 
            IF (ianz.eq.1) then 
               ww = atan (ww) 
            ELSEIF (ianz.eq.2) then 
               werte (1) = do_read_number (cpara (1), lpara (1) ) 
               werte (2) = do_read_number (cpara (2), lpara (2) ) 
               ww = atan2 (werte (1), werte (2) ) 
            ELSE 
               ier_num = - 27 
               ier_typ = ER_FORT 
            ENDIF 
            CALL ersetz2 (string, ikl, iklz, ww, 4, lll) 
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'sind') then 
            ww = sind (ww) 
            CALL ersetz2 (string, ikl, iklz, ww, 4, lll) 
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'cosd') then 
            ww = cosd (ww) 
            CALL ersetz2 (string, ikl, iklz, ww, 4, lll) 
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'tand') then 
            ww = tand (ww) 
            CALL ersetz2 (string, ikl, iklz, ww, 4, lll) 
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'tanh') then 
            ww = tanh (ww) 
            CALL ersetz2 (string, ikl, iklz, ww, 4, lll) 
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'sinh') then 
            ww = sinh (ww) 
            CALL ersetz2 (string, ikl, iklz, ww, 4, lll) 
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'cosh') then 
            ww = cosh (ww) 
            CALL ersetz2 (string, ikl, iklz, ww, 4, lll) 
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'sqrt') then 
            IF (ww.ge.0) then 
               ww = sqrt (ww) 
               CALL ersetz2 (string, ikl, iklz, ww, 4, lll) 
            ELSE 
               ier_num = - 5 
               ier_typ = ER_FORT 
            ENDIF 
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'nint') then 
            ww = float (nint (ww) ) 
            CALL ersetz2 (string, ikl, iklz, ww, 4, lll) 
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'frac') then 
            ww = ww - float (int (ww) ) 
            CALL ersetz2 (string, ikl, iklz, ww, 4, lll) 
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'gbox') then 
            CALL get_params (line, ianz, cpara, lpara, 3, lp) 
            IF (ianz.eq.3) then 
               DO i = 1, 3 
               CALL eval (cpara (i), lpara (i) ) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
               werte (i) = do_read_number (cpara (i), lpara (i) ) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
               IF (werte (i) .lt.0.0) then 
                  ier_num = - 28 
                  ier_typ = ER_APPL 
                  RETURN 
               ENDIF 
               ENDDO 
!           Determine relative size of exponentials and box             
               fl1 = 0.5 * werte (1) * sqrt (zpi) 
               fl2 = werte (2) 
               fl3 = 0.5 * werte (3) * sqrt (zpi) 
!           Normalize to 1                                              
               gbox_k = 1. / (fl1 + fl2 + fl3) 
               gbox_w = 1. - (fl1 + fl3) * gbox_k 
!           Get random number                                           
               gbox_x = ran1 (idum) 
!                                                                       
               IF (gbox_x.lt.fl1 * gbox_k) then 
                  ww = - werte (2) * 0.5 - abs (gasdev (werte (1) ) ) 
               ELSEIF (gbox_x.lt. (fl1 + fl2) * gbox_k) then 
                  ww = - werte (2) * 0.5 + (gbox_x - fl1 * gbox_k)      &
                       * werte (2) / gbox_w
               ELSE 
                  ww = werte (2) * 0.5 + abs (gasdev (werte (3) ) ) 
               ENDIF 
               CALL ersetz2 (string, ikl, iklz, ww, 4, lll) 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'gran') then 
            CALL get_params (line, ianz, cpara, lpara, 2, lp) 
            IF (ianz.ge.1) then 
               CALL eval (cpara (1), lpara (1) ) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
               werte (1) = do_read_number (cpara (1), lpara (1) ) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
               IF (werte (1) .lt.0.0) then 
                  ier_num = - 28 
                  ier_typ = ER_APPL 
                  RETURN 
               ENDIF 
!                                                                       
               IF (ianz.eq.1.or.cpara (2) .eq.'s') then 
                  a = werte (1) 
               ELSEIF (ianz.eq.2.and.cpara (2) .eq.'f') then 
                  a = werte (1) / sqrt (8. * log (2.) ) 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_FORT 
                  RETURN 
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
            ww = gasdev (a) 
            CALL ersetz2 (string, ikl, iklz, ww, 4, lll) 
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'logn') then 
            CALL get_params (line, ianz, cpara, lpara, 2, lp) 
            IF (ianz.ge.2) then 
               CALL eval (cpara (1), lpara (1) ) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
               werte (1) = do_read_number (cpara (1), lpara (1) ) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
               IF (werte (1) .lt.0.0) then 
                  ier_num = - 28 
                  ier_typ = ER_APPL 
                  RETURN 
               ENDIF 
!                                                                       
               CALL eval (cpara (2), lpara (2) ) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
               werte (2) = do_read_number (cpara (2), lpara (2) ) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
               IF (werte (2) .lt.0.0) then 
                  ier_num = - 28 
                  ier_typ = ER_APPL 
                  RETURN 
               ENDIF 
!                                                                       
               IF (ianz.eq.2.or.cpara (3) .eq.'s') then 
                  a = werte (2) 
               ELSEIF (ianz.eq.2.and.cpara (3) .eq.'f') then 
                  a = werte (2) / sqrt (8. * log (2.) ) 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_FORT 
                  RETURN 
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
            ww = exp (log (werte (1) ) + gasdev (a) ) 
            CALL ersetz2 (string, ikl, iklz, ww, 4, lll) 
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'pois') then 
            ww = poidev (ww, idum) 
            CALL ersetz2 (string, ikl, iklz, ww, 4, lll) 
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'date') then 
            CALL datum_intrinsic 
            CALL ersetzc (string, ikl, iklz, f_date, 24, 4, lll) 
         ELSE 
            CALL p_calc_intr_spec (string, line, ikl, iklz, ww, lll, lp) 
         ENDIF 
      ELSEIF (lcom.eq.3) then 
         IF (string (ikl - 3:ikl - 1) .eq.'sin') then 
            ww = sin (ww) 
            CALL ersetz2 (string, ikl, iklz, ww, 3, lll) 
         ELSEIF (string (ikl - 3:ikl - 1) .eq.'cos') then 
            ww = cos (ww) 
            CALL ersetz2 (string, ikl, iklz, ww, 3, lll) 
         ELSEIF (string (ikl - 3:ikl - 1) .eq.'tan') then 
            ww = tan (ww) 
            CALL ersetz2 (string, ikl, iklz, ww, 3, lll) 
         ELSEIF (string (ikl - 3:ikl - 1) .eq.'exp') then 
            ww = exp (ww) 
            CALL ersetz2 (string, ikl, iklz, ww, 3, lll) 
         ELSEIF (string (ikl - 3:ikl - 1) .eq.'abs') then 
            ww = abs (ww) 
            CALL ersetz2 (string, ikl, iklz, ww, 3, lll) 
         ELSEIF (string (ikl - 3:ikl - 1) .eq.'int') then 
            ww = float (int (ww) ) 
            CALL ersetz2 (string, ikl, iklz, ww, 3, lll) 
         ELSEIF (string (ikl - 3:ikl - 1) .eq.'max') then 
            CALL get_params (line, ianz, cpara, lpara, 2, lp) 
            IF (ianz.eq.2) then 
               DO i = 1, ianz 
               CALL eval (cpara (i), lpara (i) ) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
               werte (i) = do_read_number (cpara (i), lpara (i) ) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
               ENDDO 
               ww = max (werte (1), werte (2) ) 
               CALL ersetz2 (string, ikl, iklz, ww, 3, lll) 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
         ELSEIF (string (ikl - 3:ikl - 1) .eq.'min') then 
            CALL get_params (line, ianz, cpara, lpara, 2, lp) 
            IF (ianz.eq.2) then 
               DO i = 1, ianz 
               CALL eval (cpara (i), lpara (i) ) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
               werte (i) = do_read_number (cpara (i), lpara (i) ) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
               ENDDO 
               ww = min (werte (1), werte (2) ) 
               CALL ersetz2 (string, ikl, iklz, ww, 3, lll) 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
         ELSEIF (string (ikl - 3:ikl - 1) .eq.'mod') then 
            CALL get_params (line, ianz, cpara, lpara, 2, lp) 
            IF (ianz.eq.2) then 
               DO i = 1, ianz 
               CALL eval (cpara (i), lpara (i) ) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
               werte (i) = do_read_number (cpara (i), lpara (i) ) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
               ENDDO 
               ww = amod (werte (1), werte (2) ) 
               CALL ersetz2 (string, ikl, iklz, ww, 3, lll) 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
         ELSEIF (string (ikl - 3:ikl - 1) .eq.'ran') then 
            ww = ran1 (idum) 
            CALL ersetz2 (string, ikl, iklz, ww, 3, lll) 
         ELSE 
            CALL p_calc_intr_spec (string, line, ikl, iklz, ww, lll, lp) 
         ENDIF 
      ELSEIF (lcom.eq.2) then 
         IF (string (ikl - 2:ikl - 1) .eq.'ln') then 
            IF (ww.gt.0.0) then 
               ww = log (ww) 
               CALL ersetz2 (string, ikl, iklz, ww, 2, lll) 
            ELSE 
               ier_num = - 20 
               ier_typ = ER_FORT 
            ENDIF 
         ELSE 
            CALL p_calc_intr_spec (string, line, ikl, iklz, ww, lll, lp) 
         ENDIF 
      ELSE 
         CALL p_calc_intr_spec (string, line, ikl, iklz, ww, lll, lp) 
      ENDIF 
!                                                                       
      END SUBROUTINE calc_intr                      
!*****7**************************************************************** 
      SUBROUTINE get_w1_w2 (w1, w2, line, iverk, iz1, iz2, ll, lverk,   &
      lreal)                                                            
!                                                                       
!     Get the two numbers in front and after an operator                
!                                                                       
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      CHARACTER ( * ) line 
      CHARACTER(1024) zeile 
      INTEGER suche_vor2, suche_nach2 
      INTEGER iverk, iz1, iz2, ll, lverk, lll 
      LOGICAL lreal 
      REAL w1, w2 
!                                                                       
      REAL do_read_number 
!                                                                       
      lll = iverk - 1 
      iz1 = suche_vor2 (line (1:iverk - 1), lll) 
      zeile = line (iz1:iverk - 1) 
      w1 = do_read_number (zeile, iverk - iz1) 
      IF (ier_num.ne.0) then 
         RETURN 
      ENDIF 
      lreal = index (line (iz1:iverk - 1) , '.') .gt.0 
      lll = ll - (iverk + lverk) + 1 
      iz2 = suche_nach2 (line (iverk + lverk:ll), lll) 
      zeile = line (iverk + lverk:iverk + lverk + iz2 - 1) 
      w2 = do_read_number (zeile, iz2) 
      IF (ier_num.ne.0) then 
         RETURN 
      ENDIF 
      lreal = lreal.or.index (line (iverk + lverk:iverk + lverk + iz2 - &
      1) , '.') .gt.0                                                   
!                                                                       
      END SUBROUTINE get_w1_w2                      
!*****7**************************************************************** 
      SUBROUTINE ersetz (line, iz1, iz2, ww, iverk, ll, lverk, lreal) 
!+                                                                      
!     Replaces the corresponding part of line by the number ww          
!     Different format is used for real and integer numbers             
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER ( * ) line 
      CHARACTER(5) form 
      CHARACTER(1024) zeile 
      INTEGER iz1, iz2, iverk, ll, lverk, lw 
      INTEGER ltot 
      LOGICAL lreal 
      REAL ww 
!                                                                       
      zeile = ' ' 
      IF (iz1.gt.1) zeile (1:iz1 - 1) = line (1:iz1 - 1) 
      IF (lreal) then 
         lw = 15 
         WRITE (zeile (iz1:iz1 + lw - 1) , '(e15.8e2)') ww 
         zeile (iz1 + 11:iz1 + 11) = 'e' 
      ELSE 
         lw = int (alog (abs (ww) + 1.) ) + 2 
         WRITE (form, 1000) lw 
         WRITE (zeile (iz1:iz1 + lw - 1), form) int (ww) 
      ENDIF 
      IF (iverk + lverk + iz2.le.ll) then 
         ltot = (iz1 + lw) + (ll - iverk - lverk - iz2 + 1) - 1 
         IF (ltot.le.len (zeile) ) then 
            zeile (iz1 + lw:ltot) = line (iverk + lverk + iz2:ll) 
         ENDIF 
      ENDIF 
      line = zeile 
      ll = ll + lw - (lverk + iz2 + iverk - iz1 - 1) 
      CALL rem_bl (line, ll) 
!                                                                       
 1000 FORMAT    ('(i',i2,')') 
      END SUBROUTINE ersetz                         
!*****7**************************************************************** 
      SUBROUTINE ersetz2 (string, ikl, iklz, ww, lfunk, lll) 
!                                                                       
!     Replaces the intrinsic function and its argument by the           
!     corresponding value ww                                            
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER ( * ) string 
      CHARACTER(1024) zeile 
      INTEGER ikl, iklz, laenge, lfunk, lll 
      INTEGER ltot 
      REAL ww 
!                                                                       
      laenge = lll 
      zeile = ' ' 
      IF (ikl.gt.1) zeile (1:ikl - 1 - lfunk) = string (1:ikl - 1 - lfunk)
      WRITE (zeile (ikl - lfunk:ikl - lfunk + 14) , '(e15.8e2)') ww 
      zeile (ikl - lfunk + 11:ikl - lfunk + 11) = 'e' 
      lll = ikl - lfunk + 14 
      IF (iklz + 1.le.laenge) then 
         ltot = (ikl - lfunk + 15) + (laenge-iklz - 1 + 1) - 1 
         IF (ltot.le.len (zeile) ) then 
            zeile (ikl - lfunk + 15:ltot) = string (iklz + 1:laenge) 
            lll = lll + laenge- (iklz + 1) + 1 
         ENDIF 
      ENDIF 
      string = zeile 
      CALL rem_bl (string, lll) 
      END SUBROUTINE ersetz2                        
!*****7**************************************************************** 
      SUBROUTINE ersetzc (string, ikl, iklz, result, l_result, lfunk, lll)
!                                                                       
!     Replaces the intrinsic character function and its argument by the 
!     corresponding value line                                          
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER ( * ) string 
      CHARACTER ( * ) result 
      CHARACTER(1024) zeile 
      INTEGER ikl, iklz, laenge, lfunk, lll 
      INTEGER ltot 
      INTEGER l_result 
!                                                                       
      laenge = lll 
      zeile = ' ' 
      IF (ikl.gt.1) zeile (1:ikl - 1 - lfunk) = string (1:ikl - 1 - lfunk)
      zeile (ikl - lfunk:ikl - lfunk + l_result - 1) = result 
      lll = ikl - lfunk + l_result - 1 
      IF (iklz + 1.le.laenge) then 
         ltot = (ikl - lfunk + l_result) + (laenge-iklz - 1 + 1) - 1
         IF (ltot.le.len (zeile) ) then 
            zeile (ikl - lfunk + l_result:ltot) = string (iklz + 1: laenge)
         ENDIF 
         lll = lll + laenge- (iklz + 1) + 1 
      ENDIF 
      IF (zeile (1:1) .eq.'('.and.zeile (lll:lll) .eq.')') then 
         string = zeile (2:lll - 1) 
         lll = lll - 2 
      ELSE 
         string = zeile 
      ENDIF 
      END SUBROUTINE ersetzc                        
!*****7**************************************************************** 
      INTEGER FUNCTION bindex (string, char) 
!                                                                       
!     Searches for 'char' in 'string' from the back                     
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER ( * ) string 
      CHARACTER(1) char 
!                                                                       
      INTEGER il, ik 
      INTEGER len_str 
!                                                                       
      ik = 0 
      il = len_str (string) 
!                                                                       
      DO while (il.ge.1.and.ik.eq.0) 
      IF (string (il:il) .eq.char) then 
         ik = il 
      ELSE 
         il = il - 1 
      ENDIF 
      ENDDO 
!                                                                       
      bindex = ik 
!                                                                       
      END FUNCTION bindex                           
!*****7**************************************************************** 
      SUBROUTINE rem_bl (line, ll) 
!                                                                       
!     Removes all blanks from a string                                  
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER ( * ) line 
      CHARACTER(1024) zeile 
      INTEGER ll, ibl 
!                                                                       
      ibl = index (line (1:ll) , ' ') 
      DO while (ibl.gt.0) 
      zeile = ' ' 
      IF (ibl.gt.1) zeile (1:ibl - 1) = line (1:ibl - 1) 
      IF (ibl.lt.ll) then 
         zeile (ibl:ll - 1) = line (ibl + 1:ll) 
      ENDIF 
      ll = ll - 1 
      line = zeile 
      ibl = index (line (1:ll) , ' ') 
      ENDDO 
      END SUBROUTINE rem_bl                         
!*****7**************************************************************** 
      SUBROUTINE rem_leading_bl (line, ll) 
!                                                                       
!     Removes all leading blanks from a string                                  
!                                                                       
      USE charact_mod
      IMPLICIT none 
!                                                                       
      CHARACTER ( LEN=* ), INTENT(INOUT) :: line 
      INTEGER            , INTENT(INOUT) :: ll
      CHARACTER(LEN=1024)  :: zeile 
      INTEGER              :: i,j
!                                                                       
      zeile = ' '
      j     = 1
main: DO i = 1, ll
         j = i
         IF( .NOT. (line(i:i)==' ' .or. line(i:i)==tab)) THEN
            EXIT main
         ENDIF
      ENDDO main
      zeile = line(j:ll)
      ll    = ll - j + 1
      line  = zeile
      END SUBROUTINE rem_leading_bl                         
!*****7***************************************************************  
INTEGER FUNCTION len_str ( string )
!
! Returns the length of the string without trailing blanks, TABS, CR or LF
!
IMPLICIT NONE
!
CHARACTER (LEN=*), INTENT(IN) :: string
CHARACTER (LEN=3), PARAMETER  :: WHITE = achar(9)//achar(10)//achar(13)
INTEGER                       :: itab
INTEGER                       :: laenge
!
laenge = LEN_TRIM(string)               ! intrinsic LEN_TRIM, length without trailing blanks
itab   = SCAN(string(1:laenge),WHITE,.true.) ! to LEN_TRIM a TAB,CR, LF is non-white !@!%!?!
!
DO
  IF ( itab<laenge ) EXIT                ! TAB,CR, LF is not last character thats it
  IF ( laenge==0   ) EXIT                ! empty string
!
!  IF ( laenge==1 ) THEN                 ! length one AND tab is last character
!    laenge = 0                          ! No need to test this, as Fortran 2003 yields
!    EXIT                                ! a string of zero length for 1:laenge-1
!  ENDIF                                 ! with laenge == 1

  laenge = LEN_TRIM(string(1:laenge-1))
  itab   = SCAN(string(1:laenge),WHITE,.true.) ! .true. means scan backwards
ENDDO
!
len_str = laenge
!
END FUNCTION len_str
!*****7***************************************************************  
      SUBROUTINE del_params (ndel, ianz, cpara, lpara, nwerte) 
!-                                                                      
!     Deletes the first <ndel> parameters from the list and             
!     moves the rest to the beginning of the arrays.                    
!+                                                                      
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER nwerte 
      CHARACTER(1024) cpara (nwerte) 
      INTEGER lpara (nwerte) 
      INTEGER ndel, ianz, i 
!                                                                       
      IF (ndel.lt.ianz) then 
         DO i = ndel + 1, ianz 
         cpara (i - ndel) = cpara (i) 
         lpara (i - ndel) = lpara (i) 
         ENDDO 
         ianz = ianz - ndel 
      ENDIF 
!                                                                       
      END SUBROUTINE del_params                     
!*****7***************************************************************  
      SUBROUTINE get_params (string, ianz, cpara, lpara, nwerte, laenge) 
!-                                                                      
!     Reads parameters that have to be separated by ",". Each           
!       expression between "," is evaluated to give the value of        
!       the parameter. If the routine is called with a negative         
!     length, leading SPACE will NOT be removed !                       
!+                                                                      
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER nwerte 
      CHARACTER ( * ) string 
      CHARACTER(1024) cpara (nwerte) 
      INTEGER lpara (nwerte) 
      INTEGER i, lll 
      INTEGER laenge, ipos, ianz
      INTEGER level 
      LOGICAL quote 
      LOGICAL search 
      LOGICAL rmblk 
!                                                                       
      ipos    = 1 
      ianz    = 0 
      ier_num = 0 
      ier_typ = ER_NONE 
      DO i = 1, nwerte 
         cpara (i) = ' ' 
         lpara (i) = 1 
      ENDDO
!
      IF (laenge.eq.0) return 
!                                                                       
      IF (laenge.lt.0) then 
         laenge = - laenge 
         rmblk = .false. 
      ELSE 
         rmblk = .true. 
      ENDIF 
!
      IF ( len_trim(string).eq.0) THEN 
        laenge = 0
        RETURN
      ENDIF
!                                                                       
      search = .true. 
      level = 0 
      ianz = 1 
      lpara (ianz) = 0 
      DO i = 1, laenge 
      quote = string (i:i) .eq.'"'.or.string (i:i) .eq.'''' 
      IF (quote) then 
         search = .not.search 
      ENDIF 
      IF (search) then 
         IF (level.eq.0.and.string (i:i) .eq.',') then 
            IF (ianz.lt.nwerte) then 
               ianz = ianz + 1 
               lpara (ianz) = 0 
            ELSE 
               ier_num = - 17 
               ier_typ = ER_COMM 
               RETURN 
            ENDIF 
         ELSE 
            lpara (ianz) = lpara (ianz) + 1 
            cpara (ianz) (lpara (ianz) :lpara (ianz) ) = string (i:i) 
            IF (string (i:i) .eq.'('.or.string (i:i) .eq.'[') then 
               level = level + 1 
            ENDIF 
            IF (string (i:i) .eq.')'.or.string (i:i) .eq.']') then 
               level = level - 1 
            ENDIF 
         ENDIF 
      ELSE 
         lpara (ianz) = lpara (ianz) + 1 
         cpara (ianz) (lpara (ianz) :lpara (ianz) ) = string (i:i) 
         IF (string (i:i) .eq.'('.or.string (i:i) .eq.'[') then 
            level = level + 1 
         ENDIF 
         IF (string (i:i) .eq.')'.or.string (i:i) .eq.']') then 
            level = level - 1 
         ENDIF 
      ENDIF 
      ENDDO 
!                                                                       
!     remove leading blanks if length was given positive                
!                                                                       
      IF (rmblk) then 
         DO i = 1, ianz 
         lll = lpara (i) 
         CALL rem_bl (cpara (i), lll) 
         lpara (i) = lll 
         ENDDO 
      ENDIF 
!                                                                       
!      Check parameter length to catch things like ...,,...             
!      ianz is decremented until ALL get_parameter calls do             
!      an error check...                                                
!                                                                       
      DO i = 1, ianz 
      IF (lpara (i) .eq.0) then 
         ier_num = - 2 
         ier_typ = ER_COMM 
         ianz = max (0, ianz - 1) 
         RETURN 
      ENDIF 
      ENDDO 
!                                                                       
!      Check if we are still in search mode, ie. missing ""s            
!      ianz is decremented until ALL get_parameter calls do             
!      an error check...                                                
!                                                                       
      IF (.not.search) then 
         ier_num = - 30 
         ier_typ = ER_FORT 
         ianz = max (0, ianz - 1) 
      ENDIF 
!                                                                       
      IF (level.ne.0) then 
         ier_num = - 9 
         ier_typ = ER_FORT 
         ianz = max (0, ianz - 1) 
      ENDIF 
!                                                                       
      RETURN 
!CCC                                                                    
!     OLD CODE, KEPT UNTIL ERROR CHECK IS COMPLETED!                    
!                                                                       
!DBG      iko=index(string(ipos:laenge),',')                            
!DBG      if(iko.ne.0) then                                             
!DBG        lll = laenge - ipos + 1                                     
!DBG        ikk=suche_nach(string(ipos:laenge),lll)+1                   
!DBG        if(ikk.eq.laenge) then                                      
!DBG          ikk=ikk+1                                                 
!DBG        endif                                                       
!DBG        do while(ikk.gt.1 .and. ianz.lt.nwerte-1)                   
!DBG          ianz       =ianz+1                                        
!DBG          cpara(ianz)=string(ipos:ipos+ikk-2)                       
!DBG          lpara(ianz)= (ipos+ikk-2) - ipos + 1                      
!DBG          ipos       =ipos+ikk                                      
!DBG          if(ipos.le.laenge) then                                   
!DBG            iko=index(string(ipos:laenge),',')                      
!DBG            if(iko.ne.0) then                                       
!DBG              lll = laenge - ipos + 1                               
!DBG              ikk=suche_nach(string(ipos:laenge),lll)+1             
!DBG              if(ikk+ipos.gt.laenge) goto 10                        
!DBG            else                                                    
!DBG              ikk=0                                                 
!DBG            endif                                                   
!DBG          else                                                      
!DBG            goto 10                                                 
!DBG          endif                                                     
!DBG        ENDDO                                                       
!DBG      endif                                                         
!DBG10      continue                                                    
!DBG      if(ipos.le.laenge .and. string(ipos:laenge).ne. ' ') then     
!DBG        ianz       =ianz+1                                          
!DBG        cpara(ianz)=string(ipos:laenge)                             
!DBG        lpara(ianz)=laenge - ipos + 1                               
!DBG      endif                                                         
!DBG      ier_num = 0                                                   
!DBG      ier_typ = ER_NONE                                             
!DBGc                                                                   
!DBGc     remove leading blanks if length was given positive            
!DBGc                                                                   
!DBG      if (rmblk) then                                               
!DBG        do i=1,ianz                                                 
!DBG          lll      = lpara(i)                                       
!DBG          call rem_bl(cpara(i),lll)                                 
!DBG          lpara(i) = lll                                            
!DBG        ENDDO                                                       
!DBG      endif                                                         
      END SUBROUTINE get_params                     
!*****7***************************************************************  
      SUBROUTINE ber_params (ianz, cpara, lpara, werte, maxpara) 
!-                                                                      
!     Calculated the value of all expressions stored in cpara           
!+                                                                      
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER maxpara 
      CHARACTER(1024) cpara (maxpara) 
      INTEGER lpara (maxpara) 
      CHARACTER(1024) line 
      INTEGER ll, ianz, i 
      REAL werte (maxpara), wert 
      REAL berechne 
!                                                                       
!                                                                       
      DO i = 1, ianz 
      ll = lpara (i) 
      line = ' ' 
      line = '('//cpara (i) (1:ll) //')' 
      ll = ll + 2 
      wert = berechne (line, ll) 
      IF (ier_num.ne.0) then 
         GOTO 9999 
      ENDIF 
      werte (i) = wert 
      ENDDO 
 9999 CONTINUE 
      END SUBROUTINE ber_params                     
!*****7***************************************************************  
      INTEGER FUNCTION length_com (string, ikl) 
!-                                                                      
!     Determines the length of the variable or intrinsic function       
!     by searching backwards from the bracket to the first non          
!     character and non '_' character.                                  
!+                                                                      
      USE charact_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      CHARACTER ( * ) string 
!                                                                       
      INTEGER ikl, i, c 
      LOGICAL lchar 
!                                                                       
      i = ikl 
      lchar = .true. 
      DO while (lchar.and.i.gt.1) 
      i = i - 1 
      c = iachar (string (i:i) ) 
      lchar = a.le.c.and.c.le.z.or.c.eq.u 
      ENDDO 
!                                                                       
      IF (i.eq.1.and.lchar) then 
         length_com = ikl - 1 
      ELSE 
         length_com = ikl - i - 1 
      ENDIF 
!                                                                       
      END FUNCTION length_com                       
!*****7**************************************************************** 
      REAL FUNCTION do_read_number (string, laenge) 
!-                                                                      
!     Reads a numerical value from the string. This is the VMS version. 
!     In order to recognize erroneous nonmunerical data in the last     
!     character of string, the format is one character longer than the  
!     actual; length of the parameter. Otherwise unix does not          
!     recognize an integer like e.g. "3a" as an error.                  
!                                                                       
!     The VMS version properly detects an error for the (more sensible) 
!     length equal to the actual length of the parameter                
!+                                                                      
      USE errlist_mod 
      USE charact_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      CHARACTER ( * ) string 
      INTEGER laenge 
!                                                                       
      CHARACTER(1024) line 
      INTEGER i 
!                                                                       
      ier_num = - 1 
      ier_typ = ER_FORT 
      line = string 
!                                                                       
      IF (laenge.gt.0) then 
         i = iachar (string (laenge:laenge) ) 
         IF (zero.le.i.and.i.le.nine.or.i.eq.period.or.i.eq.blank1) then
            READ (line, *, end = 999, err = 999) do_read_number 
!                                                                       
            ier_num = 0 
            ier_typ = ER_NONE 
         ELSE 
            do_read_number = 0. 
         ENDIF 
      ENDIF 
!                                                                       
  999 CONTINUE 
!                                                                       
      END FUNCTION do_read_number                   
!*****7***************************************************************  
      SUBROUTINE do_build_name (ianz, cpara, lpara, werte, MAXW, fpara) 
!-                                                                      
!     Creates a filename from several parameters. the first parameter   
!     is interpreted as format.                                         
!+                                                                      
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER FORM_INTEGER 
      INTEGER FORM_INT_ZEROS 
      INTEGER FORM_REAL 
      INTEGER FORM_REAL_ZEROS 
      INTEGER FORM_CHARACTER 
      INTEGER FORM_CHAR_ZEROS 
      PARAMETER (FORM_INTEGER = 1) 
      PARAMETER (FORM_INT_ZEROS = 2) 
      PARAMETER (FORM_REAL = 3) 
      PARAMETER (FORM_REAL_ZEROS = 4) 
      PARAMETER (FORM_CHARACTER = 5) 
      PARAMETER (FORM_CHAR_ZEROS = 6) 
      INTEGER NUM_FORM 
      PARAMETER (NUM_FORM = 6) 
!                                                                       
      INTEGER MAXW 
!                                                                       
      CHARACTER(1024) cpara (MAXW) 
      CHARACTER(1024) string 
      INTEGER ianz 
      INTEGER fpara 
      INTEGER lpara (MAXW) 
      REAL werte (MAXW) 
!                                                                       
      CHARACTER(1) c_form (NUM_FORM) 
      CHARACTER(1024) fstring 
      CHARACTER(1024) line, number 
      INTEGER ftyp 
      INTEGER ind_ql, ind_qr 
      INTEGER ind_d, ind_i, ind_f 
      INTEGER ind_p, itot, idec 
      INTEGER ii, lp 
      INTEGER laenge 
      INTEGER lwert 
      INTEGER lsub 
      INTEGER npara 
      INTEGER pos 
      INTEGER ll, i 
      INTEGER fstring_l 
      INTEGER ising 
      REAL wert 
      LOGICAL lfloat 
      LOGICAL l_hyp 
      LOGICAL l_sin 
!                                                                       
      REAL berechne 
      INTEGER len_str 
!                                                                       
      DATA c_form / 'd', 'D', 'f', 'F', 'c', 'C' / 
!                                                                       
      string = ' ' 
      lwert = 0
      wert = 0.
      fstring_l = 0
      IF (lpara (fpara) .le.1) then 
!                                                                       
!     --If first parameter is too short for a string, return            
!                                                                       
         RETURN 
      ENDIF 
!                                                                       
!     Remove leading blanks                                             
!                                                                       
      DO i = 1, ianz 
      ll = 1 
      line (1:1) = cpara (i) (ll:ll) 
      DO while (line (1:1) .eq.' ') 
      ll = ll + 1 
      line (1:1) = cpara (i) (ll:ll) 
      IF (ll.gt.lpara (i) ) goto 7575 
      ENDDO 
 7575 CONTINUE 
      string = cpara (i) (ll:lpara (i) ) 
      cpara (i) = string (1:lpara (i) - ll + 1) 
      lpara (i) = lpara (i) - ll + 1 
      ENDDO 
!                                                                       
!     Get position of left quotation mark into ind_ql                   
!                                                                       
      ind_ql = index (cpara (fpara) (1:lpara (fpara) - 1) , '"') 
      IF (ind_ql.lt.1) then 
!                                                                       
!     --No quotation mark found, return                                 
!                                                                       
         RETURN 
      ENDIF 
!                                                                       
!     Get position of right quotation mark into ind_qr                  
!                                                                       
      ind_qr = ind_ql + index(cpara(fpara)(ind_ql + 1:lpara(fpara)),'"')
      IF (ind_qr.lt.ind_ql + 1) then 
         ier_num = - 30 
         ier_typ = ER_FORT 
         cpara (fpara) = ' ' 
         lpara (fpara) = 1 
         RETURN 
      ENDIF 
!                                                                       
!     Find position of first format specifiers                          
!                                                                       
      ind_d = index (cpara (fpara)(ind_ql + 1:ind_qr - 1) , '%') + ind_ql
!                                                                       
      laenge = 0 
      pos = ind_ql + 1 
      npara = fpara + 1 
!                                                                       
!     While there are any format specifiers left                        
!                                                                       
      DO while (ind_d.ge.pos) 
!                                                                       
!     --If there is space left for a format specifier                   
!                                                                       
      IF (ind_d+1.le.ind_qr - 1) then 
!                                                                       
!     ----Find first 'd' or 'f'                                         
!                                                                       
         ind_f = 9999 
         ftyp = 0 
         DO i = 1, NUM_FORM 
         ind_i = index (cpara (fpara) (ind_d+1:ind_qr - 1), c_form (i) ) 
         IF (ind_i.gt.0.and.ind_i.lt.ind_f) then 
            ind_f = ind_i 
            ftyp = i 
         ENDIF 
         ENDDO 
         IF (ftyp.eq.0) then 
            ier_num = - 21 
            ier_typ = ER_FORT 
            cpara (fpara) = ' ' 
            lpara (fpara) = 1 
            RETURN 
         ENDIF 
         ind_i = ind_f 
!                                                                       
!     --We have an integer format                                       
!       The string between '%' and 'd' or 'D' is evaluated and          
!       written into the format variable fstring. Allows error checking 
!       and flexible formats                                            
!                                                                       
         IF (ftyp.eq.FORM_INTEGER.or.ftyp.eq.FORM_INT_ZEROS) then 
            ll = ind_i - 1 
            IF (ll.gt.0) then 
      line = '('//cpara (fpara)  (ind_d+1:ind_d+ind_i - 1) //')' 
               ll = ll + 2 
               CALL rem_bl (line (1:ll), ll) 
               ier_num = 0 
               wert = berechne (line, ll) 
               IF (ier_num.ne.0) then 
                  ier_msg (1) = 'An error occurred while calculating' 
      ier_msg (2)  = 'the value of a format specifier    ' 
                  cpara (fpara) = ' ' 
                  lpara (fpara) = 1 
                  RETURN 
               ENDIF 
               WRITE (fstring, 3000) int (wert) 
               ll = 15 
               CALL rem_bl (fstring (1:ll), ll) 
            ELSE 
               fstring = '*' 
               ll = 1 
            ENDIF 
!                                                                       
!     --We have a real format                                           
!       The string between '%' and 'f' or 'F' is evaluated and          
!       written into the format variable fstring. Allows error checking 
!       and flexible formats                                            
!                                                                       
         ELSEIF (ftyp.eq.FORM_REAL.or.ftyp.eq.FORM_REAL_ZEROS) then 
            ll = ind_f - 1 
            IF (ll.gt.0) then 
      line = '('//cpara (fpara)  (ind_d+1:ind_d+ind_f - 1) //')' 
               ll = ll + 2 
               CALL rem_bl (line (1:ll), ll) 
               ind_p = index (line (1:ll) , '.') 
               IF (ind_p.gt.0) then 
                  number = line (1:ind_p - 1) //')' 
                  lp = ind_p 
                  itot = nint (berechne (number, lp) ) 
                  IF (ier_num.ne.0) then 
      ier_msg (1)  = 'An error occurred while calculating' 
      ier_msg (2)  = 'the value of a format specifier    ' 
                     cpara (fpara) = ' ' 
                     lpara (fpara) = 1 
                     RETURN 
                  ENDIF 
                  IF (ind_p.lt.ll - 1) then 
                     number = '('//line (ind_p + 1:ll) 
                     lp = ll - ind_p + 1 
                     idec = nint (berechne (number, lp) ) 
                     IF (ier_num.ne.0) then 
      ier_msg (1)  = 'An error occurred while calculating' 
      ier_msg (2)  = 'the value of a format specifier    ' 
                        cpara (fpara) = ' ' 
                        lpara (fpara) = 1 
                        RETURN 
                     ENDIF 
                  ELSE 
                     idec = 0 
                  ENDIF 
               ELSE 
                  itot = nint (berechne (line (1:ll), ll) ) 
                  IF (ier_num.ne.0) then 
      ier_msg (1)  = 'An error occurred while calculating' 
      ier_msg (2)  = 'the value of a format specifier    ' 
                     cpara (fpara) = ' ' 
                     lpara (fpara) = 1 
                     RETURN 
                  ENDIF 
                  idec = 0 
               ENDIF 
               IF (ier_num.ne.0) return 
               WRITE (fstring, 3100) itot, idec 
               ll = 24 
               CALL rem_bl (fstring (1:ll), ll) 
            ELSE 
               fstring = '*' 
               ll = 1 
            ENDIF 
!                                                                       
!     --We have a character format                                      
!       The string between '%' and 'c' or 'C' is evaluated and          
!       written into the format variable fstring. Allows error checking 
!       and flexible formats                                            
!                                                                       
         ELSEIF (ftyp.eq.FORM_CHARACTER.or.ftyp.eq.FORM_CHAR_ZEROS) then
            ll = ind_i - 1 
            IF (ll.gt.0) then 
      line = '('//cpara (fpara)  (ind_d+1:ind_d+ind_i - 1) //')' 
               ll = ll + 2 
               CALL rem_bl (line (1:ll), ll) 
               ier_num = 0 
               wert = berechne (line, ll) 
               IF (ier_num.ne.0) then 
                  ier_msg (1) = 'An error occurred while calculating' 
      ier_msg (2)  = 'the value of a format specifier    ' 
                  cpara (fpara) = ' ' 
                  lpara (fpara) = 1 
                  RETURN 
               ENDIF 
               WRITE (fstring, 3200) int (wert) 
               fstring_l = int (wert) 
               ll = 15 
               CALL rem_bl (fstring (1:ll), ll) 
            ELSE 
               fstring = 'a' 
               ll = 1 
            ENDIF 
         ENDIF 
      ELSE 
         ier_num = - 21 
         ier_typ = ER_FORT 
         cpara (fpara) = ' ' 
         lpara (fpara) = 1 
         RETURN 
      ENDIF 
!                                                                       
!     --Copy strings left before the current format specifier.          
!                                                                       
      IF (ind_d.gt.pos) then 
         lsub = ind_d-pos 
         string (laenge+1:laenge+lsub) = cpara (fpara) (pos:ind_d-1) 
         laenge = laenge+lsub 
      ENDIF 
!                                                                       
!     --Write numerical value of parameter into string                  
!                                                                       
      IF (ianz.ge.npara) then 
         IF (.not. (ftyp.eq.FORM_CHARACTER.or.ftyp.eq.FORM_CHAR_ZEROS) )&
         then                                                           
!                                                                       
!     ------Calculate value of current parameter                        
!                                                                       
            ll = lpara (npara) 
            line = ' ' 
            line = '('//cpara (npara) (1:ll) //')' 
            ll = ll + 2 
            CALL rem_bl (line (1:ll), ll) 
            wert = berechne (line, ll) 
            IF (ier_num.ne.0) then 
               cpara (fpara) = ' ' 
               lpara (fpara) = 1 
               RETURN 
            ENDIF 
            werte (npara) = wert 
         ELSE 
!                                                                       
!     ------Calculate value of current character parameter              
!           if enclosed in " " take it as is, otherwise calculate       
!                                                                       
            ll = lpara (npara) 
            l_hyp = .false. 
            IF (cpara (npara) (1:1) .eq.'"') then 
               l_hyp = .true. 
               l_sin = .false. 
            ENDIF 
            IF (cpara (npara) (1:1) .eq.'''') then 
               l_hyp = .true. 
               l_sin = .true. 
            ENDIF 
            IF (l_hyp) then 
               IF (l_sin) then 
                  ising = index (cpara (npara) (2:ll) , '''') + 1 
               ELSE 
                  ising = index (cpara (npara) (2:ll) , '"') + 1 
               ENDIF 
               IF (ising.eq.ll) then 
                  cpara (npara) = cpara (npara) (2:ll - 1) 
                  ll = ll - 2 
                  lpara (npara) = ll 
               ELSEIF (ising.gt.1) then 
                  line (1:1) = '(' 
                  line (2:ising - 1) = cpara (npara) (2:ising - 1) 
                  line (ising:ll - 1) = cpara (npara) (ising + 1:ll) 
                  line (ll:ll) = ')' 
                  CALL berechne_char (line, ll) 
                  IF (line (1:1) .eq.''''.and.line (ll:ll) .eq.'''') then
                     cpara (npara) = line (2:ll) 
                     ll = ll - 2 
                  ELSE 
                     cpara (npara) = line (1:ll) 
                  ENDIF 
                  IF (ier_num.ne.0) then 
                     cpara (fpara) = ' ' 
                     lpara (fpara) = 1 
                     RETURN 
                  ENDIF 
               ELSE 
                  ier_num = - 30 
                  ier_typ = ER_FORT 
                  cpara (fpara) = ' ' 
                  lpara (fpara) = 1 
                  RETURN 
               ENDIF 
            ELSE 
               line = ' ' 
               line = '('//cpara (npara) (1:ll) //')' 
               ll = ll + 2 
               CALL berechne_char (line, ll) 
               IF (line (1:1) .eq.''''.and.line (ll:ll) .eq.'''') then 
                  cpara (npara) = line (2:ll) 
                  ll = ll - 2 
               ELSE 
                  cpara (npara) = line (1:ll) 
               ENDIF 
               lpara (npara) = ll 
               IF (ier_num.ne.0) then 
                  cpara (fpara) = ' ' 
                  lpara (fpara) = 1 
                  RETURN 
               ENDIF 
            ENDIF 
         ENDIF 
!                                                                       
         IF (ftyp.eq.FORM_REAL) then 
            pos = ind_d+ind_f + 1 
            IF (fstring.eq.'*') then 
               WRITE (number, * ) werte (npara) 
               ii = len_str (number) 
               CALL rem_bl (number, ii) 
            ELSE 
               WRITE (number, fstring) werte (npara) 
            ENDIF 
            ii = len_str (number) 
            lwert = len_str (number) 
         ELSEIF (ftyp.eq.FORM_REAL_ZEROS) then 
            pos = ind_d+ind_f + 1 
            IF (fstring.eq.'*') then 
               WRITE (number, * ) werte (npara) 
               ii = len_str (number) 
               CALL rem_bl (number, ii) 
            ELSE 
               WRITE (number, fstring) werte (npara) 
            ENDIF 
            ii = len_str (number) 
            DO i = 1, ii 
            IF (number (i:i) .eq.' ') then 
               number (i:i) = '0' 
            ENDIF 
            ENDDO 
            DO i = 1, ii 
            IF (number (i:i) .eq.'-') then 
               number (i:i) = '0' 
               number (1:1) = '-' 
            ENDIF 
            ENDDO 
            lwert = len_str (number) 
         ELSEIF (ftyp.eq.FORM_INTEGER) then 
            pos = ind_d+ind_i + 1 
            ii = nint (werte (npara) ) 
            IF (fstring.eq.'*') then 
               WRITE (number, * ) ii 
               ii = len_str (number) 
               CALL rem_bl (number, ii) 
            ELSE 
               WRITE (number, fstring) ii 
            ENDIF 
            ii = len_str (number) 
            lwert = len_str (number) 
         ELSEIF (ftyp.eq.FORM_INT_ZEROS) then 
            pos = ind_d+ind_i + 1 
            ii = nint (werte (npara) ) 
            IF (fstring.eq.'*') then 
               WRITE (number, * ) ii 
               ii = len_str (number) 
               CALL rem_bl (number, ii) 
            ELSE 
               WRITE (number, fstring) ii 
            ENDIF 
            ii = len_str (number) 
            DO i = 1, ii 
            IF (number (i:i) .eq.' ') then 
               number (i:i) = '0' 
            ENDIF 
            ENDDO 
            DO i = 1, ii 
            IF (number (i:i) .eq.'-') then 
               number (i:i) = '0' 
               number (1:1) = '-' 
            ENDIF 
            ENDDO 
            lwert = len_str (number) 
         ELSEIF (ftyp.eq.FORM_CHARACTER) then 
            pos = ind_d+ind_i + 1 
            IF (fstring.eq.'a') then 
               number = cpara (npara) (1:lpara (npara) ) 
               lwert = lpara (npara) 
            ELSE 
               WRITE (number, fstring) cpara (npara) 
               lwert = fstring_l 
            ENDIF 
         ENDIF 
         IF (lwert.eq.0) then 
            ier_num = - 34 
            ier_typ = ER_FORT 
            RETURN 
         ENDIF 
         IF ( (number (1:1) .eq.''''.or.number (1:1) .eq.'"') .and.  &
              (number (lwert:lwert) .eq.''''.or.                     &
               number (lwert:lwert) .eq.'"') ) then
            string (laenge+1:laenge+lwert - 2) = number (2:lwert - 1) 
            laenge = laenge+lwert - 2 
         ELSE 
            string (laenge+1:laenge+lwert) = number (1:lwert) 
            laenge = laenge+lwert 
         ENDIF 
!                                                                       
!     ----Eliminate parameter from list                                 
!                                                                       
         DO i = npara, ianz - 1 
         cpara (i) = cpara (i + 1) 
         lpara (i) = lpara (i + 1) 
         ENDDO 
         cpara (ianz) = ' ' 
         lpara (ianz) = 0 
         ianz = ianz - 1 
      ELSE 
         ier_num = - 43 
         ier_typ = ER_APPL 
         cpara (fpara) = ' ' 
         lpara (fpara) = 1 
         RETURN 
      ENDIF 
!                                                                       
!     --Any format specifiers left ?                                    
!                                                                       
      IF (pos.le.ind_qr - 1) then 
         ind_d = index (cpara (fpara) (pos:ind_qr - 1) , '%') + pos - 1 
         lfloat = (cpara (fpara) (ind_d+1:ind_d+1) .eq.'f') 
      ELSE 
         ind_d = 0 
      ENDIF 
      ENDDO 
!                                                                       
!     Append rest of format specifying string to filename               
!                                                                       
      IF (ind_qr - 1.ge.pos) then 
         string (laenge+1:laenge+ind_qr - pos) = &
              cpara (fpara) (pos: ind_qr - 1)
         laenge = laenge+ind_qr - pos 
      ENDIF 
!                                                                       
      cpara (fpara) = string 
      lpara (fpara) = laenge 
!                                                                       
 3000 FORMAT    ('(i',i10,')') 
 3100 FORMAT    ('(f',i10,'.',i10,')') 
 3200 FORMAT    ('(a',i10,')') 
      END SUBROUTINE do_build_name                  
!*****7***************************************************************  
      SUBROUTINE count_col (zeile, ianz) 
!+                                                                      
!     This subroutine counts the number of columns of string 'zeile'    
!-                                                                      
      USE prompt_mod 
!
      IMPLICIT NONE
!                                                                       
      CHARACTER ( * ) zeile 
      INTEGER ianz, i 
      LOGICAL ein 
!                                                                       
      INTEGER len_str 
!                                                                       
      ianz = 0 
      ein = .false. 
!                                                                       
      DO i = 1, len_str (zeile) 
      IF (zeile (i:i) .ne.' ') then 
         ein = .true. 
      ELSEIF (zeile (i:i) .eq.' '.and.ein) then 
         ein = .false. 
         ianz = ianz + 1 
      ENDIF 
      ENDDO 
!                                                                       
      IF (ein) ianz = ianz + 1 
      END SUBROUTINE count_col                      
