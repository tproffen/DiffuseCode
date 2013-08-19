!******7****************************************************************
!     Here are menu and mouse related commands                          
!******7****************************************************************
      SUBROUTINE do_mouse (zeile, lp) 
!+                                                                      
!     Mouse command ..                                                  
!-                                                                      
      IMPLICIT none 
!                                                                       
      include'config.inc' 
      include'kuplot.inc' 
      include'errlist.inc' 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 4) 
!                                                                       
      CHARACTER ( * ) zeile 
      INTEGER lp 
!                                                                       
      CHARACTER(1024) cpara (maxw) 
      INTEGER lpara (maxw), ianz 
      REAL werte (maxw) 
!                                                                       
      LOGICAL str_comp 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, - lp) 
      IF (ier_num.ne.0) return 
!                                                                       
      IF (ianz.eq.0) then 
         lgui (iwin) = .true. 
         CALL do_menu 
         lgui (iwin) = .false. 
      ELSE 
         IF (str_comp (cpara (1) , 'point', 1, lpara (1) , 5) ) then 
            CALL do_koor_nogui 
         ELSEIF (str_comp (cpara (1) , 'line', 1, lpara (1) , 4) ) then 
            CALL do_region_nogui (1) 
         ELSEIF (str_comp (cpara (1) , 'rect', 1, lpara (1) , 4) ) then 
            CALL do_region_nogui (2) 
         ELSEIF (str_comp (cpara (1) , 'yrange', 1, lpara (1) , 6) )    &
         then                                                           
            CALL do_region_nogui (3) 
         ELSEIF (str_comp (cpara (1) , 'xrange', 1, lpara (1) , 6) )    &
         then                                                           
            CALL do_region_nogui (4) 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
      ENDIF 
!                                                                       
      END SUBROUTINE do_mouse                       
!******7****************************************************************
      SUBROUTINE do_menu 
!                                                                       
!     Activate KUPLOT menu                                              
!                                                                       
      IMPLICIT none 
!                                                                       
      include'config.inc' 
      include'prompt.inc' 
      include'kuplot.inc' 
      include'learn.inc' 
      include'errlist.inc' 
!                                                                       
      CHARACTER(1) ch 
      REAL xb, xt, yt, mx, my 
      INTEGER ini 
      LOGICAL lw, lend 
!                                                                       
      LOGICAL hit_button, n_in_f 
!                                                                       
      lw = .true. 
!                                                                       
      CALL do_plot (.true.) 
      IF (ier_num.ne.0) return 
      IF (iaf (iwin) .gt.1) call border_frame (iframe, 1, 0.02) 
!                                                                       
      WRITE (output_io, 1000) 
!                                                                       
!------ Here is the menu input loop                                     
!                                                                       
   10 CONTINUE 
      CALL frame_menu 
      CALL PGCURS (mx, my, ch) 
      IF (ch.eq.butt_l) then 
!                                                                       
!------ --- Exit selected                                               
!                                                                       
         IF (hit_button (1, mx, my) ) then 
            CALL do_plot (.false.) 
            IF (llearn) then 
               WRITE (33, 2000) ex (iwin, iframe, 1), ex (iwin, iframe, &
               2), ey (iwin, iframe, 1), ey (iwin, iframe, 2)           
            ENDIF 
            RETURN 
!                                                                       
!------ --- Coordinate command                                          
!                                                                       
         ELSEIF (hit_button (2, mx, my) ) then 
            CALL draw_button (2, .false., .true.) 
            CALL do_koor 
            CALL draw_button (2, .false., .false.) 
!                                                                       
!------ --- Select region command                                       
!                                                                       
         ELSEIF (hit_button (3, mx, my) ) then 
            CALL draw_button (3, .false., .true.) 
            CALL do_region 
            CALL draw_button (3, .false., .false.) 
!                                                                       
!------ --- Zoom                                                        
!                                                                       
         ELSEIF (hit_button (4, mx, my) ) then 
            CALL draw_button (4, .false., .true.) 
            CALL do_zoom (.true.) 
            CALL draw_button (4, .false., .false.) 
!                                                                       
!------ --- Move                                                        
!                                                                       
         ELSEIF (hit_button (5, mx, my) ) then 
            CALL draw_button (5, .false., .true.) 
            CALL do_zoom (.false.) 
            CALL draw_button (5, .false., .false.) 
!                                                                       
!------ --- Reset region                                                
!                                                                       
         ELSEIF (hit_button (6, mx, my) ) then 
            CALL draw_button (6, .false., .true.) 
            ex (iwin, iframe, 1) = - 9999. 
            ey (iwin, iframe, 1) = - 9999. 
            t (iwin, iframe, 1) = - 9999. 
            t (iwin, iframe, 2) = - 9999. 
            CALL draw_frame (iframe, lw) 
            CALL frame_menu 
            CALL draw_tframe (' ', ' ', 'Region reset ...') 
            CALL draw_button (6, .false., .false.) 
!                                                                       
!------ --- Select frame command                                        
!                                                                       
         ELSEIF (hit_button (7, mx, my) .and.iaf (iwin) .gt.1) then 
            CALL draw_button (7, .false., .true.) 
            CALL do_fselect 
            CALL draw_button (7, .false., .false.) 
!                                                                       
!------ --- Select distance command                                     
!                                                                       
         ELSEIF (hit_button (8, mx, my) ) then 
            CALL draw_button (8, .false., .true.) 
            CALL do_distance 
            CALL draw_button (8, .false., .false.) 
!                                                                       
!------ --- Select ZMIN                                                 
!                                                                       
         ELSEIF (hit_button (9, mx, my) .and.n_in_f (ini) ) then 
            CALL draw_button (9, .false., .true.) 
            CALL do_zscale (.true.) 
            CALL draw_button (9, .false., .false.) 
!                                                                       
!------ --- Select ZMAX                                                 
!                                                                       
         ELSEIF (hit_button (10, mx, my) .and.n_in_f (ini) ) then 
            CALL draw_button (10, .false., .true.) 
            CALL do_zscale (.false.) 
            CALL draw_button (10, .false., .false.) 
!                                                                       
!------ --- Enter KUPLOT command                                        
!                                                                       
         ELSEIF (hit_button (11, mx, my) ) then 
            CALL draw_button (11, .false., .true.) 
            CALL do_enter_command (lend) 
            IF (lend) return 
            CALL draw_button (11, .false., .false.) 
         ENDIF 
      ENDIF 
      IF (iaf (iwin) .gt.1) call border_frame (iframe, 1, 0.02) 
      GOTO 10 
!                                                                       
!------ Finish off                                                      
!                                                                       
 9999 CONTINUE 
                                                                        
!                                                                       
 1000 FORMAT    (' ------ > Click on EXIT MENU to return to ',          &
     &                            'command mode ...')                   
 2000 FORMAT     ('skal ',3(G12.6,','),G12.6) 
      END SUBROUTINE do_menu                        
!******7****************************************************************
      SUBROUTINE frame_menu 
!                                                                       
!     This sets the menu viewport and draws borders ...                 
!                                                                       
      IMPLICIT none 
!                                                                       
      include'config.inc' 
      include'kuplot.inc' 
!                                                                       
      REAL xpl (2), ypl (2) 
!                                                                       
      CALL PGSVP (0.0, 0.999, 0.0, 0.999) 
      CALL PGSWIN (0.0, 1.0, 0.0, 1.0) 
!                                                                       
      CALL PGSFS (2) 
      CALL PGSCI (6) 
      CALL PGRECT (0.0, 1.0, 0.0, 1.0) 
!                                                                       
      END SUBROUTINE frame_menu                     
!******7****************************************************************
      SUBROUTINE draw_menu 
!                                                                       
!     Draw the different menu parts                                     
!                                                                       
      IMPLICIT none 
!                                                                       
      include'config.inc' 
      include'kuplot.inc' 
      include'prompt.inc' 
      include'errlist.inc' 
!                                                                       
      CHARACTER(40) text 
      REAL xh, yh, xt, yt 
      INTEGER ini 
!                                                                       
      INTEGER len_str 
      LOGICAL n_in_f 
!                                                                       
!------ Draw borders and filling of menu region                         
!                                                                       
      CALL PGSVP (0.0, 0.999, 0.0, 0.999) 
      CALL PGSWIN (0.0, 1.0, 0.0, 1.0) 
      CALL PGSCI (6) 
      CALL PGSFS (1) 
      CALL PGRECT (1.0 - dev_draw (iwin, 1), 1.0, 0.0, 1.0) 
      CALL PGRECT (0.0, 1.0 - dev_draw (iwin, 1), 0.0, dev_draw (iwin,  &
      2) )                                                              
!                                                                       
!------ Plot KUPLOT logo                                                
!                                                                       
      CALL PGSCF (1) 
      CALL PGSCH (1.0) 
      CALL PGSCI (1) 
      xt = 1.0 - 0.50 * dev_draw (iwin, 1) 
      yt = 0.0 + 0.60 * dev_draw (iwin, 2) 
      CALL PGQCS (4, xh, yh) 
      CALL PGPTXT (xt, yt, 0.0, 0.5, 'KUPLOT') 
      yt = yt - 1.05 * yh 
      text = 'Version '//version (1:len_str (version) ) 
      CALL PGPTXT (xt, yt, 0.0, 0.5, text) 
!                                                                       
!------ Add menu buttons                                                
!                                                                       
      xt = 1.0 - 0.50 * dev_draw (iwin, 1) 
      yt = dev_draw (iwin, 1) * 0.9 
!                                                                       
      CALL def_button (2, 'Coordinates', xt, 0.95, yt, 0.04) 
      CALL def_button (8, 'Distances', xt, 0.90, yt, 0.04) 
      CALL def_button (3, 'Select region', xt, 0.80, yt, 0.04) 
      CALL def_button (4, 'Zoom', xt, 0.75, yt, 0.04) 
      CALL def_button (5, 'Move', xt, 0.70, yt, 0.04) 
      CALL def_button (6, 'Reset region', xt, 0.65, yt, 0.04) 
      CALL def_button (9, 'Select zmin', xt, 0.55, yt, 0.04) 
      CALL def_button (10, 'Select zmax', xt, 0.50, yt, 0.04) 
      CALL def_button (7, 'Select frame', xt, 0.40, yt, 0.04) 
      CALL def_button (11, 'Enter command', xt, 0.25, yt, 0.04) 
      CALL def_button (1, 'Exit menu', xt, 0.20, yt, 0.04) 
!                                                                       
      CALL draw_button (1, .false., .false.) 
      CALL draw_button (2, .false., .false.) 
      CALL draw_button (3, .false., .false.) 
      CALL draw_button (4, .false., .false.) 
      CALL draw_button (5, .false., .false.) 
      CALL draw_button (6, .false., .false.) 
      CALL draw_button (7, (iaf (iwin) .eq.1), .false.) 
      CALL draw_button (8, .false., .false.) 
      CALL draw_button (9, .not.n_in_f (ini), .false.) 
      CALL draw_button (10, .not.n_in_f (ini), .false.) 
      CALL draw_button (11, .false., .false.) 
!                                                                       
      CALL draw_tframe (' ', ' ', 'Use the EXIT button to return to comm&
     &and mode ...')                                                    
!                                                                       
      END SUBROUTINE draw_menu                      
!******7****************************************************************
      SUBROUTINE draw_tframe (tt1, tt2, tt3) 
!                                                                       
!     This routine writes the strings tt1-tt3 below plot                
!                                                                       
      IMPLICIT none 
!                                                                       
      include'config.inc' 
      include'kuplot.inc' 
!                                                                       
      CHARACTER ( * ) tt1, tt2, tt3 
      CHARACTER(80) t1, t2, t3 
      REAL xt, yt, xh, yh 
      INTEGER ltext 
!                                                                       
      INTEGER len_str 
!                                                                       
      t1 = tt1 
      t2 = tt2 
      t3 = tt3 
!                                                                       
!------ Remove old text                                                 
!                                                                       
      CALL PGSCI (6) 
      CALL PGSFS (1) 
      CALL PGRECT (0.0, 1.0 - dev_draw (iwin, 1), 0.0, dev_draw (iwin,  &
      2) )                                                              
      CALL frame_menu 
!                                                                       
!------ Plot new text                                                   
!                                                                       
      CALL PGSCF (1) 
      CALL PGSCI (0) 
      CALL PGSCH (1.25) 
      CALL PGQCS (4, xh, yh) 
      xt = 0.5 * (1.0 - dev_draw (iwin, 1) ) 
      yt = dev_draw (iwin, 2) - 1.5 * yh 
      ltext = len_str (t1) 
      IF (ltext.gt.0) call PGPTXT (xt, yt, 0.0, 0.5, t1 (1:ltext) ) 
      yt = yt - 1.1 * yh 
      ltext = len_str (t2) 
      IF (ltext.gt.0) call PGPTXT (xt, yt, 0.0, 0.5, t2 (1:ltext) ) 
      CALL PGSCH (1.0) 
      CALL PGQCS (4, xh, yh) 
      yt = yt - 1.5 * yh 
      ltext = len_str (t3) 
      IF (ltext.gt.0) call PGPTXT (xt, yt, 0.0, 0.5, t3 (1:ltext) ) 
!                                                                       
      END SUBROUTINE draw_tframe                    
!******7****************************************************************
      SUBROUTINE draw_command (cmd, befehl, ic, ib) 
!                                                                       
!     Get command from GUI area ..                                      
!                                                                       
      IMPLICIT none 
!                                                                       
      include'config.inc' 
      include'kuplot.inc' 
!                                                                       
      CHARACTER ( * ) cmd, befehl 
      INTEGER ic, ib 
      REAL xt, yt, xh, yh, x1, x2, y1, y2 
!                                                                       
      ic = 0 
      cmd = ' ' 
!                                                                       
      CALL PGSCI (6) 
      CALL PGSFS (1) 
      CALL PGRECT (0.0, 1.0 - dev_draw (iwin, 1), 0.0, dev_draw (iwin,  &
      2) )                                                              
!                                                                       
      CALL frame_menu 
!                                                                       
!------ Plot prompt                                                     
!                                                                       
      CALL PGSCF (1) 
      CALL PGSCH (1.0) 
      CALL PGQCS (4, xh, yh) 
      CALL PGSCI (1) 
!                                                                       
      x1 = 0.0 + 0.2 * xh 
      x2 = 1.0 - dev_draw (iwin, 1) 
      y1 = 1.80 * yh 
      y2 = 0.20 * yh 
      xt = x1 + 0.50 * xh 
      yt = y2 + 0.40 * yh 
!                                                                       
      CALL PGPTXT (xt, yt, 0.0, 0.0, 'kuplot >') 
      CALL PGRSTR (xt + 5.0 * xh, yt, 0.0, 0.0, cmd, ic, 6) 
      CALL extr_befehl (cmd, befehl, ic, ib) 
!                                                                       
      END SUBROUTINE draw_command                   
!******7****************************************************************
      SUBROUTINE extr_befehl (cmd, befehl, ic, ib) 
!                                                                       
!     Extract 'befehl' string                                           
!                                                                       
      IMPLICIT none 
!                                                                       
      include'config.inc' 
      include'kuplot.inc' 
!                                                                       
      CHARACTER ( * ) cmd, befehl 
      INTEGER ic, ib, il, ir, len_str 
!                                                                       
      ic = len_str (cmd) 
      il = 1 
      DO while (cmd (il:il) .eq.' '.and.il.lt.ic) 
      il = il + 1 
      ENDDO 
      ir = index (cmd (il:ic) , ' ') - 1 
      IF (ir.lt.0) ir = ic 
!                                                                       
      befehl = cmd (il:ir) 
      ib = len_str (befehl) 
!                                                                       
      END SUBROUTINE extr_befehl                    
!******7****************************************************************
      LOGICAL function hit_button (ib, mx, my) 
!                                                                       
!     Checks if button was pressed                                      
!                                                                       
      IMPLICIT none 
!                                                                       
      include'config.inc' 
      include'kuplot.inc' 
!                                                                       
      REAL mx, my, w2, h2 
      INTEGER ib 
!                                                                       
      w2 = 0.5 * bw (ib) 
      h2 = 0.5 * bh (ib) 
!                                                                       
      hit_button = mx.ge. (bx (ib) - w2) .and.mx.le. (bx (ib) + w2)     &
      .and.my.ge. (by (ib) - h2) .and.my.le. (by (ib) + h2)             
!                                                                       
      END FUNCTION hit_button                       
!******7****************************************************************
      LOGICAL function n_in_f (ini) 
!                                                                       
!     Checks for 2D files in current frame                              
!                                                                       
      IMPLICIT none 
!                                                                       
      include'config.inc' 
      include'kuplot.inc' 
!                                                                       
      INTEGER i, ini 
      LOGICAL k_in_f 
!                                                                       
      n_in_f = .false. 
      DO i = 1, iz - 1 
      n_in_f = n_in_f.or. (k_in_f (i) .and.lni (i) ) 
      IF (n_in_f) ini = i 
      ENDDO 
!                                                                       
      END FUNCTION n_in_f                           
!******7****************************************************************
      SUBROUTINE def_button (ib, text, xb, yb, wb, hb) 
!                                                                       
!     Defines a button for the menu                                     
!                                                                       
      IMPLICIT none 
!                                                                       
      include'config.inc' 
      include'kuplot.inc' 
!                                                                       
      CHARACTER ( * ) text 
      REAL xb, yb, wb, hb 
      INTEGER ib 
!                                                                       
      btext (ib) = text 
      bx (ib) = xb 
      by (ib) = yb 
      bw (ib) = wb 
      bh (ib) = hb 
!                                                                       
      END SUBROUTINE def_button                     
!******7****************************************************************
      SUBROUTINE draw_button (i, ldisabled, lactive) 
!                                                                       
!     This routine draws a button on menu panel                         
!                                                                       
      IMPLICIT none 
!                                                                       
      include'config.inc' 
      include'kuplot.inc' 
!                                                                       
      REAL xh, yh, w2, h2 
      INTEGER i, ib, it 
      LOGICAL lactive, ldisabled 
!                                                                       
      IF (ldisabled) then 
         ib = 6 
         it = 6 
      ELSE 
         IF (lactive) then 
            ib = 6 
            it = 0 
         ELSE 
            ib = 0 
            it = 6 
         ENDIF 
      ENDIF 
!                                                                       
      w2 = 0.5 * bw (i) 
      h2 = 0.5 * bh (i) 
!                                                                       
      CALL PGSCI (ib) 
      CALL PGSFS (1) 
      CALL PGRECT (bx (i) - w2, bx (i) + w2, by (i) - h2, by (i)        &
      + h2)                                                             
      CALL PGSFS (2) 
      CALL PGSCI (it) 
      CALL PGQCS (4, xh, yh) 
      CALL PGPTXT (bx (i), by (i) - yh * 0.3, 0.0, 0.5, btext (i) ) 
!                                                                       
      END SUBROUTINE draw_button                    
!******7****************************************************************
      SUBROUTINE do_enter_command (lend) 
!                                                                       
!     Allows the user to enter KUPLOT command                           
!                                                                       
      IMPLICIT none 
!                                                                       
      include'config.inc' 
      include'kuplot.inc' 
      include'doact.inc' 
      include'macro.inc' 
      include'errlist.inc' 
!                                                                       
      CHARACTER(1024) command, befehl 
      INTEGER ic, ib, mac_level_old 
      LOGICAL lend, lreg 
!                                                                       
      LOGICAL str_comp 
!                                                                       
      lend = .false. 
!                                                                       
!------ Now enter, check and execute command                            
!                                                                       
      mac_level_old = mac_level 
      CALL draw_command (command, befehl, ic, ib) 
      IF (ic.gt.0) then 
!                                                                       
!------ - Check for disabled commands                                   
!                                                                       
         IF (str_comp (befehl, 'help', 2, ib, 4) .or.str_comp (befehl,  &
         'do', 2, ib, 2) .or.str_comp (befehl, 'if', 2, ib, 2) ) then   
            GOTO 10 
         ELSE 
            CALL mache_kdo (command, lend, ic) 
            IF (ier_num.ne.0) goto 10 
!                                                                       
!------ --- Macro entered ?                                             
!                                                                       
            DO while (mac_level.gt.mac_level_old) 
            IF (lblock) then 
               CALL do_execute (lreg, command, ic) 
               IF (ier_num.ne.0.or..not.lreg) goto 10 
            ELSE 
               CALL macro_read (command, ic) 
               CALL extr_befehl (command, befehl, ic, ib) 
               IF (str_comp (befehl, 'do', 2, ib, 4) .or.str_comp (     &
               befehl, 'if', 2, ib, 2) ) then                           
                  CALL do_loop (command, lend, ic) 
               ELSE 
                  CALL mache_kdo (command, lend, ic) 
               ENDIF 
               IF (ier_num.ne.0) goto 10 
            ENDIF 
            ENDDO 
         ENDIF 
      ENDIF 
!                                                                       
!------ Redraw and return                                               
!                                                                       
   10 CONTINUE 
!                                                                       
      CALL do_plot (.not.lend) 
      IF (ier_num.ne.0) then 
         CALL errlist 
         IF (.not.lend) call draw_tframe ('Last command', command,      &
         ier_out)                                                       
      ENDIF 
!                                                                       
      END SUBROUTINE do_enter_command               
!******7****************************************************************
      SUBROUTINE do_zscale (lmin) 
!                                                                       
!     Set zmin,zmax via mouse                                           
!     LEFT : up by 5% / MIDDLE : down by 5% / RIGHT : back to menu      
!                                                                       
      IMPLICIT none 
!                                                                       
      include'config.inc' 
      include'kuplot.inc' 
      include'errlist.inc' 
!                                                                       
      CHARACTER(80) zeile 
      CHARACTER(1) key 
      REAL wx, wy, hub 
      INTEGER ini 
      LOGICAL lw, lmin, l2d, n_in_f 
!                                                                       
      lw = .true. 
      l2d = n_in_f (ini) 
!                                                                       
      CALL draw_tframe (' ', 'LEFT button: up 5% / MIDDLE button: down 5&
     &%', 'RIGHT botton: back')                                         
!                                                                       
!------ First point                                                     
!                                                                       
   10 CONTINUE 
      CALL open_viewport 
      CALL PGSCI (ilinecol (iwin, iframe, 0) ) 
      CALL PGBAND (0, 0, 0.0, 0.0, wx, wy, key) 
      hub = nz (iwin, iframe, 1) * z_inc (iwin, iframe, 1) 
      IF (key.eq.butt_r) then 
         GOTO 9999 
      ELSEIF (key.eq.butt_l) then 
         IF (lmin) then 
            z_min (iwin, iframe, 1) = z_min (iwin, iframe, 1) + 0.05 *  &
            hub                                                         
         ELSE 
            z_inc (iwin, iframe, 1) = (1.00 + 0.05) * hub / float (nz ( &
            iwin, iframe, 1) )                                          
         ENDIF 
      ELSEIF (key.eq.butt_m) then 
         IF (lmin) then 
            z_min (iwin, iframe, 1) = z_min (iwin, iframe, 1) - 0.05 *  &
            hub                                                         
         ELSE 
            z_inc (iwin, iframe, 1) = (1.00 - 0.05) * hub / float (nz ( &
            iwin, iframe, 1) )                                          
         ENDIF 
      ENDIF 
      CALL draw_frame (iframe, lw) 
      GOTO 10 
!                                                                       
 9999 CONTINUE 
!                                                                       
      CALL draw_frame (iframe, lw) 
      CALL frame_menu 
      WRITE (zeile, 1000) z_min (iwin, iframe, 1), nz (iwin, iframe, 1) &
      * z_inc (iwin, iframe, 1)                                         
!                                                                       
      CALL draw_tframe (zeile, ' ', ' ') 
!                                                                       
 1000 FORMAT     ('New range Z: ',G12.6,' to ',G12.6) 
!                                                                       
      END SUBROUTINE do_zscale                      
!******7****************************************************************
      SUBROUTINE do_fselect 
!                                                                       
!     Select active frame via mouse                                     
!     LEFT : select frame                                               
!                                                                       
      IMPLICIT none 
!                                                                       
      include'config.inc' 
      include'kuplot.inc' 
      include'errlist.inc' 
!                                                                       
      CHARACTER(80) zeile 
      CHARACTER(1) key 
      REAL wx, wy, fx, fy 
      INTEGER i, ifn 
      LOGICAL lw 
!                                                                       
      lw = .true. 
!                                                                       
      CALL draw_tframe (' ', ' ', 'LEFT button: select active frame') 
!                                                                       
!------ Get point                                                       
!                                                                       
      CALL PGSCI (ilinecol (iwin, iframe, 0) ) 
      CALL PGBAND (0, 0, 0.0, 0.0, wx, wy, key) 
      IF (key.ne.butt_l) return 
!                                                                       
!------ Determine frame                                                 
!                                                                       
      fx = wx / (1.0 - dev_draw (iwin, 1) ) 
      fy = wy / (1.0 - dev_draw (iwin, 2) ) - dev_draw (iwin, 2) 
!                                                                       
      ier_num = - 35 
      ier_typ = ER_APPL 
!                                                                       
      DO i = 1, iaf (iwin) 
      IF (frame (iwin, i, 1) .le.fx.and.frame (iwin, i, 3)              &
      .ge.fx.and.frame (iwin, i, 2) .le.fy.and.frame (iwin, i, 4)       &
      .ge.fy) then                                                      
         ifn = i 
         ier_num = 0 
         ier_typ = ER_NONE 
      ENDIF 
      ENDDO 
!                                                                       
      IF (ier_num.eq.0) then 
         CALL draw_frame (iframe, lw) 
         CALL draw_frame (ifn, lw) 
         CALL border_frame (ifn, 1, 0.02) 
         iframe = ifn 
      ENDIF 
!                                                                       
      WRITE (zeile, 1000) iframe 
      CALL frame_menu 
      CALL draw_tframe (zeile, ' ', ' ') 
!                                                                       
 1000 FORMAT     ('Selected active frame :',i3) 
!                                                                       
      END SUBROUTINE do_fselect                     
!******7****************************************************************
      SUBROUTINE do_region 
!                                                                       
!     Set plotting region via mouse                                     
!     LEFT : select region / RIGHT : back to showing ALL data           
!                                                                       
      IMPLICIT none 
!                                                                       
      include'config.inc' 
      include'kuplot.inc' 
      include'errlist.inc' 
!                                                                       
      CHARACTER(80) zeile (2) 
      CHARACTER(1) key 
      REAL px (2), py (2) 
      REAL wx1, wy1, wx2, wy2 
      LOGICAL lw 
!                                                                       
      lw = .true. 
!                                                                       
      CALL draw_tframe (' ', ' ', 'LEFT button: mark corners / RIGHT bot&
     &ton: reset')                                                      
!                                                                       
!------ First point                                                     
!                                                                       
      CALL open_viewport 
      CALL PGSCI (ilinecol (iwin, iframe, 0) ) 
      CALL PGBAND (0, 0, 0.0, 0.0, wx1, wy1, key) 
      IF (key.eq.butt_r) then 
         ex (iwin, iframe, 1) = - 9999.0 
         ey (iwin, iframe, 1) = - 9999.0 
         t (iwin, iframe, 1) = - 9999.0 
         t (iwin, iframe, 2) = - 9999.0 
      ENDIF 
      IF (key.ne.butt_l) goto 9999 
!                                                                       
!------ Draw cross at first point                                       
!                                                                       
      CALL PGSLS (4) 
      CALL trans_koor (wx1, wy1) 
!                                                                       
      px (1) = ex (iwin, iframe, 1) 
      py (1) = wy1 
      px (2) = ex (iwin, iframe, 2) 
      py (2) = wy1 
      CALL PGLINE (2, px, py) 
!                                                                       
      px (1) = wx1 
      py (1) = ey (iwin, iframe, 1) 
      px (2) = wx1 
      py (2) = ey (iwin, iframe, 2) 
      CALL PGLINE (2, px, py) 
!                                                                       
!------ Second point                                                    
!                                                                       
      CALL PGBAND (2, 0, wx1, wy1, wx2, wy2, key) 
      IF (key.ne.butt_l) goto 9999 
!                                                                       
!------ Set new scale                                                   
!                                                                       
      CALL trans_koor (wx2, wy2) 
!                                                                       
      ex (iwin, iframe, 1) = min (wx1, wx2) 
      ex (iwin, iframe, 2) = max (wx1, wx2) 
      ey (iwin, iframe, 1) = min (wy1, wy2) 
      ey (iwin, iframe, 2) = max (wy1, wy2) 
!                                                                       
      t (iwin, iframe, 1) = - 9999.0 
      t (iwin, iframe, 2) = - 9999.0 
!                                                                       
 9999 CONTINUE 
!                                                                       
      CALL draw_frame (iframe, lw) 
      CALL frame_menu 
      WRITE (zeile (1), 1000) ex (iwin, iframe, 1), ex (iwin, iframe, 2) 
      WRITE (zeile (2), 1010) ey (iwin, iframe, 1), ey (iwin, iframe, 2) 
!                                                                       
      CALL draw_tframe (zeile (1) , zeile (2) , ' ') 
!                                                                       
 1000 FORMAT     ('New range X: ',G12.6,' to ',G12.6) 
 1010 FORMAT     ('New range Y: ',G12.6,' to ',G12.6) 
!                                                                       
      END SUBROUTINE do_region                      
!******7****************************************************************
      SUBROUTINE do_zoom (lzoom) 
!                                                                       
!     Zoom in/out of picture                                            
!                                                                       
      IMPLICIT none 
!                                                                       
      include'config.inc' 
      include'kuplot.inc' 
      include'learn.inc' 
      include'errlist.inc' 
!                                                                       
      CHARACTER(80) zeile (2), thelp 
      CHARACTER(1) key 
      REAL wx, wy, w2, h2 
      LOGICAL lzoom, lw 
!                                                                       
      lw = .true. 
!                                                                       
      IF (lzoom) then 
      thelp = 'LEFT button: zoom in / MIDDLE button: zoom out'//' / RIGH&
     &T botton: exit'                                                   
      ELSE 
         thelp = 'LEFT button: new center / RIGHT botton: exit' 
      ENDIF 
      CALL draw_tframe (' ', ' ', thelp) 
!                                                                       
!------ Check mouse                                                     
!                                                                       
   10 CONTINUE 
      CALL open_viewport 
      CALL PGSCI (ilinecol (iwin, iframe, 0) ) 
      CALL PGBAND (0, 0, 0.0, 0.0, wx, wy, key) 
      IF (key.eq.butt_r) goto 20 
!                                                                       
!------ - Compute new boundaries                                        
!                                                                       
      CALL trans_koor (wx, wy) 
!                                                                       
      w2 = 0.5 * (ex (iwin, iframe, 2) - ex (iwin, iframe, 1) ) 
      h2 = 0.5 * (ey (iwin, iframe, 2) - ey (iwin, iframe, 1) ) 
!                                                                       
      IF (key.eq.butt_l.and.lzoom) then 
         w2 = 0.9 * w2 
         h2 = 0.9 * h2 
      ENDIF 
      IF (key.eq.butt_m.and.lzoom) then 
         w2 = 1.1 * w2 
         h2 = 1.1 * h2 
      ENDIF 
!                                                                       
      ex (iwin, iframe, 1) = wx - w2 
      ex (iwin, iframe, 2) = wx + w2 
      ey (iwin, iframe, 1) = wy - h2 
      ey (iwin, iframe, 2) = wy + h2 
      t (iwin, iframe, 1) = - 9999.0 
      t (iwin, iframe, 2) = - 9999.0 
!                                                                       
!------ - Redraw frame                                                  
!                                                                       
      CALL draw_frame (iframe, lw) 
!                                                                       
      CALL frame_menu 
      WRITE (zeile (1), 1000) ex (iwin, iframe, 1), ex (iwin, iframe, 2) 
      WRITE (zeile (2), 1010) ey (iwin, iframe, 1), ey (iwin, iframe, 2) 
      CALL draw_tframe (zeile (1), zeile (2), thelp) 
      GOTO 10 
   20 CONTINUE 
!                                                                       
      CALL frame_menu 
      WRITE (zeile (1), 1000) ex (iwin, iframe, 1), ex (iwin, iframe, 2) 
      WRITE (zeile (2), 1010) ey (iwin, iframe, 1), ey (iwin, iframe, 2) 
      CALL draw_tframe (zeile (1) , zeile (2) , ' ') 
!                                                                       
 1000 FORMAT     ('New range X: ',G12.6,' to ',G12.6) 
 1010 FORMAT     ('New range Y: ',G12.6,' to ',G12.6) 
!                                                                       
      END SUBROUTINE do_zoom                        
!******7****************************************************************
      SUBROUTINE do_distance 
!                                                                       
!     Shows distance between two points ..                              
!                                                                       
      IMPLICIT none 
!                                                                       
      include'config.inc' 
      include'kuplot.inc' 
      include'param.inc' 
      include'errlist.inc' 
!                                                                       
      CHARACTER(80) zeile 
      CHARACTER(1) key 
      REAL x1, y1, x2, y2, dist 
!                                                                       
      zeile = ' ' 
!                                                                       
      res_para (0) = 0 
   10 CONTINUE 
      CALL frame_menu 
      CALL draw_tframe (zeile, ' ', 'LEFT button: select points / RIGHT &
     &button: exit')                                                    
      CALL open_viewport 
      CALL PGSCI (ilinecol (iwin, iframe, 0) ) 
!                                                                       
      CALL PGBAND (7, 0, 0.0, 0.0, x1, y1, key) 
      IF (key.eq.butt_r) goto 20 
      CALL PGBAND (1, 0, x1, y1, x2, y2, key) 
      IF (key.eq.butt_r) goto 20 
!                                                                       
      y1 = yskal (iwin, iframe) * y1 
      y2 = yskal (iwin, iframe) * y2 
      dist = sqrt ( (x2 - x1) **2 + (y2 - y1) **2) 
!                                                                       
      WRITE (zeile, 2200) dist 
      IF (res_para (0) - 1.le.maxpar_res) then 
         res_para (0) = res_para (0) + 1 
         res_para (NINT(res_para (0)) - 1) = dist 
      ENDIF 
      GOTO 10 
!                                                                       
   20 CONTINUE 
!                                                                       
      CALL frame_menu 
!                                                                       
 2200 FORMAT    ('Distance = ',g12.6) 
!                                                                       
      END SUBROUTINE do_distance                    
!******7****************************************************************
      SUBROUTINE do_koor 
!                                                                       
!     Show coordinates of mouse pointer position                        
!                                                                       
      IMPLICIT none 
!                                                                       
      include'config.inc' 
      include'kuplot.inc' 
      include'param.inc' 
      include'errlist.inc' 
!                                                                       
      CHARACTER(80) zeile 
      CHARACTER(1) key 
      REAL wx, wy, wz, dxx, dyy 
      INTEGER i, ini, nxx, nyy 
      LOGICAL l2d, k_in_f 
!                                                                       
      zeile = ' ' 
!                                                                       
      l2d = .false. 
      DO i = 1, iz - 1 
      l2d = l2d.or. (k_in_f (i) .and.lni (i) ) 
      IF (l2d) ini = i 
      ENDDO 
!                                                                       
      res_para (0) = 0 
   10 CONTINUE 
      CALL frame_menu 
      CALL draw_tframe (zeile, ' ', 'LEFT button: select point / RIGHT b&
     &utton: exit')                                                     
      CALL open_viewport 
      CALL PGSCI (ilinecol (iwin, iframe, 0) ) 
      CALL PGBAND (7, 0, 0.0, 0.0, wx, wy, key) 
      IF (key.eq.butt_r) goto 20 
!                                                                       
      CALL trans_koor (wx, wy) 
!                                                                       
      IF (.not.l2d) then 
         WRITE (zeile, 2200) wx, wy 
         IF (res_para (0) - 2.le.maxpar_res) then 
            res_para (0) = res_para (0) + 2 
            res_para (NINT(res_para (0)) - 1) = wx 
            res_para (NINT(res_para (0)) ) = wy 
         ENDIF 
!                                                                       
      ELSE 
         dxx = (xmax (ini) - xmin (ini) ) / float (nx (ini) - 1) 
         dyy = (ymax (ini) - ymin (ini) ) / float (ny (ini) - 1) 
         nxx = nint ( ( (wx - xmin (ini) ) / dxx) ) + 1 
         nyy = nint ( ( (wy - ymin (ini) ) / dyy) ) + 1 
         wz = z (offz (ini - 1) + (nxx - 1) * ny (ini) + nyy) 
         WRITE (zeile, 2210) wx, wy, wz 
         IF (res_para (0) - 3.le.maxpar_res) then 
            res_para (0) = res_para (0) + 3 
            res_para (NINT(res_para (0)) - 2) = wx 
            res_para (NINT(res_para (0)) - 1) = wy 
            res_para (NINT(res_para (0)) ) = wz 
         ENDIF 
      ENDIF 
      GOTO 10 
!                                                                       
   20 CONTINUE 
!                                                                       
      CALL frame_menu 
      CALL draw_tframe (' ', ' ', ' ') 
!                                                                       
 1200 FORMAT    (1x,2(g12.6,1x)) 
 1210 FORMAT    (1x,3(g12.6,1x)) 
 2200 FORMAT    ('(x,y) = ',2(g12.6,1x)) 
 2210 FORMAT    ('(x,y,z) = ',3(g12.6,1x)) 
!                                                                       
      END SUBROUTINE do_koor                        
!******7****************************************************************
      SUBROUTINE do_koor_nogui 
!                                                                       
!     Get coordinates of mouse pointer position outside GUI window      
!                                                                       
      IMPLICIT none 
!                                                                       
      include'config.inc' 
      include'kuplot.inc' 
      include'prompt.inc' 
      include'param.inc' 
      include'errlist.inc' 
!                                                                       
      CHARACTER(1) key 
      REAL wx, wy, wz 
      REAL dxx, dyy 
      INTEGER i, ini, ikey 
      INTEGER nxx, nyy 
      LOGICAL l2d, k_in_f 
!                                                                       
      CALL do_plot (lgui (iwin) ) 
      IF (ier_num.ne.0) return 
!                                                                       
      WRITE (output_io, 1000) 
!                                                                       
      l2d = .false. 
      DO i = 1, iz - 1 
      l2d = l2d.or. (k_in_f (i) .and.lni (i) ) 
      IF (l2d) ini = i 
      ENDDO 
!                                                                       
      res_para (0) = 0 
      CALL open_viewport 
      CALL PGSCI (ilinecol (iwin, iframe, 0) ) 
      CALL PGBAND (7, 0, 0.0, 0.0, wx, wy, key) 
      IF (key.eq.butt_l) ikey = 1 
      IF (key.eq.butt_m) ikey = 2 
      IF (key.eq.butt_r) ikey = 3 
!                                                                       
      CALL trans_koor (wx, wy) 
!                                                                       
      IF (.not.l2d) then 
         res_para (0) = 3 
         res_para (1) = wx 
         res_para (2) = wy 
         res_para (3) = float (ikey) 
      ELSE 
         dxx = (xmax (ini) - xmin (ini) ) / float (nx (ini) - 1) 
         dyy = (ymax (ini) - ymin (ini) ) / float (ny (ini) - 1) 
         nxx = nint ( ( (wx - xmin (ini) ) / dxx) ) + 1 
         nyy = nint ( ( (wy - ymin (ini) ) / dyy) ) + 1 
         wz = z (offz (ini - 1) + (nxx - 1) * ny (ini) + nyy) 
         res_para (0) = 4 
         res_para (1) = wx 
         res_para (2) = wy 
         res_para (3) = wz 
         res_para (4) = float (ikey) 
      ENDIF 
!                                                                       
 1000 FORMAT     (' ------ > Select desired point using the mouse ...') 
      END SUBROUTINE do_koor_nogui                  
!******7****************************************************************
      SUBROUTINE do_region_nogui (imode) 
!                                                                       
!     Set region via mouse - without menu                               
!                                                                       
      IMPLICIT none 
!                                                                       
      include'config.inc' 
      include'prompt.inc' 
      include'kuplot.inc' 
      include'param.inc' 
      include'errlist.inc' 
!                                                                       
      CHARACTER(80) zeile (2) 
      CHARACTER(1) key 
      REAL wx1, wy1, wx2, wy2 
      INTEGER ikey, imode 
!                                                                       
      CALL do_plot (lgui (iwin) ) 
      IF (ier_num.ne.0) return 
!                                                                       
      WRITE (output_io, 1000) 
      ikey = - 1 
!                                                                       
!------ First point                                                     
!                                                                       
      CALL open_viewport 
      CALL PGSCI (ilinecol (iwin, iframe, 0) ) 
      CALL PGBAND (0, 0, 0.0, 0.0, wx1, wy1, key) 
      IF (key.ne.butt_l) goto 9999 
!                                                                       
!------ Second point                                                    
!                                                                       
      CALL PGBAND (imode, 0, wx1, wy1, wx2, wy2, key) 
      IF (key.eq.butt_l) ikey = 1 
      IF (key.eq.butt_m) ikey = 2 
      IF (key.eq.butt_r) ikey = 3 
!                                                                       
!------ Set new scale                                                   
!                                                                       
      CALL trans_koor (wx1, wy1) 
      CALL trans_koor (wx2, wy2) 
!                                                                       
      res_para (0) = 5 
      res_para (1) = wx1 
      res_para (2) = wy1 
      res_para (3) = wx2 
      res_para (4) = wy2 
      res_para (5) = float (ikey) 
!                                                                       
 9999 CONTINUE 
!                                                                       
 1000 FORMAT     (' ------ > Select desired points using the mouse ...') 
      END SUBROUTINE do_region_nogui                
!******7****************************************************************
      SUBROUTINE trans_koor (wx, wy) 
!                                                                       
!     Convert screen coordinates to world coordinates for ANGL not 90.  
!     and for log. axis ..                                              
!                                                                       
      IMPLICIT none 
!                                                                       
      include'config.inc' 
      include'kuplot.inc' 
      include'wink.inc' 
!                                                                       
      REAL wx, wy 
!                                                                       
      IF (shear (iwin, iframe) .ne.90.0) then 
         wx = wx - (wy - ey (iwin, iframe, 1) ) * yskal (iwin, iframe)  &
         / tan (rad * shear (iwin, iframe) )                            
      ENDIF 
!                                                                       
      IF (lachse (iwin, iframe, 1) ) wx = exp (wx * log (10.0) ) 
      IF (lachse (iwin, iframe, 2) ) wy = exp (wy * log (10.0) ) 
!                                                                       
      END SUBROUTINE trans_koor                     
!******7****************************************************************
      SUBROUTINE PGRSTR (X, Y, ANGLE, FJUST, TEXT, LSTR, BCI) 
      REAL X, Y, ANGLE, FJUST 
      CHARACTER ( * ) TEXT 
      INTEGER LSTR, BCI 
!-----------------------------------------------------------------------
      CHARACTER CH 
      INTEGER CI 
      REAL XCUR, YCUR, XBOX (4), YBOX (4) 
!                                                                       
      CALL PGQCI (CI) 
!                                                                       
   10 CONTINUE 
!     -- Draw current string                                            
      IF (LSTR.GT.0) THEN 
         CALL PGPTXT (X, Y, ANGLE, FJUST, TEXT (1:LSTR) ) 
         CALL PGQTXT (X, Y, ANGLE, FJUST, TEXT (1:LSTR), XBOX, YBOX) 
         XCUR = XBOX (4) 
         YCUR = YBOX (4) 
      ELSE 
         XCUR = X 
         YCUR = Y 
      ENDIF 
!         -- Read a character                                           
      CALL PGBAND (0, 1, XCUR, YCUR, XCUR, YCUR, CH) 
!         -- Erase old string                                           
      CALL PGSCI (BCI) 
      IF (LSTR.GT.0) CALL PGPTXT (X, Y, ANGLE, FJUST, TEXT (1:LSTR) ) 
      CALL PGSCI (CI) 
!         -- Avoid problem with PGPLOT escape character                 
      IF (CH.EQ.CHAR (92) ) CH = '*' 
!         -- Backspace (ctrl H) or delete removes last character        
      IF (ICHAR (CH) .EQ.8.OR.ICHAR (CH) .EQ.127) THEN 
         IF (LSTR.GT.0) TEXT (LSTR:LSTR) = ' ' 
         IF (LSTR.GT.0) LSTR = LSTR - 1 
!         -- Ctrl U removes entire string                               
      ELSEIF (ICHAR (CH) .EQ.21) THEN 
         TEXT (1:LSTR) = ' ' 
         LSTR = 0 
!         -- Any other non-printing character terminates input          
      ELSEIF (ICHAR (CH) .LT.32) THEN 
         IF (LSTR.GT.0) CALL PGPTXT (X, Y, ANGLE, FJUST, TEXT (1:LSTR) ) 
         GOTO 20 
!         -- Otherwise, add character to string if there is room        
      ELSEIF (LSTR.LT.LEN (TEXT) ) THEN 
         LSTR = LSTR + 1 
         TEXT (LSTR:LSTR) = CH 
      ENDIF 
      GOTO 10 
!                                                                       
   20 RETURN 
      END SUBROUTINE PGRSTR                         
