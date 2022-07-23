module kuplot_draw_mod
!
!******7****************************************************************
!     Here are menu and mouse related commands                          
!******7****************************************************************
!
contains
!
!******7****************************************************************
!
SUBROUTINE do_mouse (zeile, lp) 
!+                                                                      
!     Mouse command ..                                                  
!-                                                                      
USE errlist_mod 
USE get_params_mod
USE kuplot_config 
USE kuplot_mod 
!
USE precision_mod
USE str_comp_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER, PARAMETER :: maxw = 4 
!                                                                       
CHARACTER(len=*), intent(inout) :: zeile 
INTEGER         , intent(inout) :: lp 
!                                                                       
CHARACTER(LEN=PREC_STRING), dimension(MAXW) :: cpara ! (maxw) 
INTEGER                   , dimension(MAXW) :: lpara ! (maxw), ianz 
integer :: ianz
!                                                                       
!                                                                       
      lp = -lp
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
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
      USE errlist_mod 
      USE learn_mod 
      USE prompt_mod 
      USE kuplot_config 
      USE kuplot_mod 
use kuplot_draw_low_mod
use kuplot_draw_tframe_mod
use kuplot_low_mod
use kuplot_plot_mod
use kuplot_plot_low_mod

!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER(1) ch 
      REAL mx, my 
      INTEGER ini 
      LOGICAL lw, lend 
!                                                                       
!      LOGICAL hit_button, n_in_f 
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
      IF (ch.eq.butt_l.OR.ch==keyb_l) then 
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
!                                                                       
!------ --- Select layer
!
         ELSEIF (hit_button (12, mx, my) ) THEN
            CALL draw_button (12, .false., .true.) 
            CALL do_layer(.true.)
continue
            CALL draw_button (12, .false., .false.) 
         ENDIF 
      ENDIF 
      IF (iaf (iwin) .gt.1) call border_frame (iframe, 1, 0.02) 
      GOTO 10 
!                                                                       
!------ Finish off                                                      
!                                                                       
!                                                                       
 1000 FORMAT    (' ------ > Click on EXIT MENU to return to ',          &
     &                            'command mode ...')                   
 2000 FORMAT     ('skal ',3(G13.6,','),G13.6) 
      END SUBROUTINE do_menu                        
!
!*******************************************************************************
!
!******7****************************************************************
      SUBROUTINE draw_command (cmd, befehl, ic, ib) 
!                                                                       
!     Get command from GUI area ..                                      
!                                                                       
      USE kuplot_config 
      USE kuplot_mod 
use kuplot_draw_tframe_mod
!                                                                       
      IMPLICIT none 
!                                                                       
CHARACTER(len=*) , intent(out) :: cmd
CHARACTER(len=*) , intent(out) :: befehl 
INTEGER , intent(out) :: ic, ib 
integer :: i6
!
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
      xt = xt + 5.0 * xh
      x1 = 0.0
      y1 = 0.0
      i6 = 6
      CALL PGRSTR (xt           , yt, x1 , y1 , cmd, ic, i6) 
      CALL extr_befehl (cmd, befehl, ic, ib) 
!                                                                       
      END SUBROUTINE draw_command                   
!******7****************************************************************
      SUBROUTINE extr_befehl (cmd, befehl, ic, ib) 
!                                                                       
!     Extract 'befehl' string                                           
!                                                                       
      USE kuplot_config 
      USE kuplot_mod 
USE lib_length
!                                                                       
      IMPLICIT none 
!                                                                       
CHARACTER(len=*) , intent(out) :: cmd
CHARACTER(len=*) , intent(out) :: befehl 
INTEGER , intent(out) :: ic, ib 
!
      INTEGER ::      il, ir
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
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER , intent(in) :: ib 
      REAL , intent(in) :: mx, my
      REAL ::      w2, h2 
!                                                                       
      w2 = 0.5 * bw (ib) 
      h2 = 0.5 * bh (ib) 
!                                                                       
      hit_button = mx.ge. (bx (ib) - w2) .and.mx.le. (bx (ib) + w2)     &
      .and.my.ge. (by (ib) - h2) .and.my.le. (by (ib) + h2)             
!                                                                       
      END FUNCTION hit_button                       
!******7****************************************************************
!******7****************************************************************
      SUBROUTINE do_enter_command (lend) 
!                                                                       
!     Allows the user to enter KUPLOT command                           
!                                                                       
      USE do_if_mod
      USE doact_mod 
      USE do_execute_mod
      USE errlist_mod 
      USE class_macro_internal
      USE kuplot_config 
      USE kuplot_mod 
use kuplot_draw_tframe_mod
use kuplot_plot_mod
!
USE lib_errlist_func
USE lib_macro_func
USE precision_mod
USE str_comp_mod
!                                                                       
      IMPLICIT none 
!                                                                       
logical , intent(out) :: lend
      CHARACTER(LEN=PREC_STRING) :: command, befehl 
      INTEGER ic, ib, macro_level_old 
      LOGICAL ::    lreg 
!                                                                       
!                                                                       
      lend = .false. 
!                                                                       
!------ Now enter, check and execute command                            
!                                                                       
      macro_level_old = macro_level 
      CALL draw_command (command, befehl, ic, ib) 
      IF (ic.gt.0) then 
!                                                                       
!------ - Check for disabled commands                                   
!                                                                       
         IF (str_comp (befehl, 'help', 2, ib, 4) .or.str_comp (befehl,  &
         'do', 2, ib, 2) .or.str_comp (befehl, 'if', 2, ib, 2) ) then   
            GOTO 10 
         ELSE 
            CALL kuplot_mache_kdo (command, lend, ic) 
            IF (ier_num.ne.0) goto 10 
!                                                                       
!------ --- Macro entered ?                                             
!                                                                       
            DO while (macro_level.gt.macro_level_old) 
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
                  CALL kuplot_mache_kdo (command, lend, ic) 
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
!     LEFT : up by 5% / RIGHT : down by 5% / MIDDLE : back to menu      
!                                                                       
      USE errlist_mod 
      USE kuplot_config 
      USE kuplot_mod 
use kuplot_draw_low_mod
use kuplot_draw_tframe_mod
use kuplot_low_mod
use kuplot_plot_low_mod
!                                                                       
      IMPLICIT none 
!                                                                       
logical , intent(in) :: lmin
!
      CHARACTER(80) zeile 
      CHARACTER(1) key 
      REAL wx, wy, hub 
      INTEGER ini 
      LOGICAL lw, l2d! , n_in_f 
!                                                                       
      lw = .true. 
      l2d = n_in_f (ini) 
!                                                                       
      CALL draw_tframe ('LEFT button/Keyboard l: down 5%',    &
                        'RIGHT button/Keyboard r: up 5%', &
                        'MIDDLE botton/Keyboard m: back'   )                                         
!                                                                       
!------ First point                                                     
!                                                                       
   10 CONTINUE 
      CALL open_viewport 
      CALL PGSCI (ilinecol (iwin, iframe, 0) ) 
      CALL PGBAND (0, 0, 0.0, 0.0, wx, wy, key) 
      hub = nz (iwin, iframe, 1) * z_inc (iwin, iframe, 1) 
      IF (key.eq.butt_m .OR. key==keyb_m) then 
         GOTO 9999 
      ELSEIF (key.eq.butt_r .OR. key==keyb_r) then 
         IF (lmin) then 
            z_min (iwin, iframe, 1) = z_min (iwin, iframe, 1) + 0.05 *  &
            hub                                                         
         ELSE 
            z_inc (iwin, iframe, 1) = (1.00 + 0.05) * hub / REAL(nz ( &
            iwin, iframe, 1) )                                          
         ENDIF 
      ELSEIF (key.eq.butt_l .OR. key==keyb_l) then 
         IF (lmin) then 
            z_min (iwin, iframe, 1) = z_min (iwin, iframe, 1) - 0.05 *  &
            hub                                                         
         ELSE 
            z_inc (iwin, iframe, 1) = (1.00 - 0.05) * hub / REAL(nz ( &
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
 1000 FORMAT     ('New range Z: ',G13.6,' to ',G13.6) 
!                                                                       
      END SUBROUTINE do_zscale                      
!******7****************************************************************
SUBROUTINE do_layer  (lmin) 
!                                                                       
!     Set H5 layer increment ++1 by mouse                               
!     LEFT : up by 5% / RIGHT : down by 5% / MIDDLE : back to menu      
!                                                                       
USE errlist_mod 
USE kuplot_config 
USE kuplot_mod 
use kuplot_draw_low_mod
use kuplot_draw_tframe_mod
USE kuplot_place
use lib_data_struc_h5
use kuplot_low_mod
use kuplot_plot_low_mod
USE param_mod
USE prompt_mod
!                                                                       
IMPLICIT none 
!                                                                       
LOGICAL, INTENT(IN) :: lmin
!
CHARACTER(LEN=80) :: zeile, string 
CHARACTER(LEN=1)  :: key 
REAL              :: wx, wy!, hub 
REAL              :: zz     ! Current height
INTEGER ini 
integer, dimension(3) :: h5_dims
INTEGER :: n_layer
LOGICAL :: is_direct
LOGICAL :: lw, l2d!, n_in_f 
!                                                                       
lw = .true. 
l2d = n_in_f (ini) 
!                                                                       
call hdf5_get_dims(0, h5_dims)
main_loop: DO
   n_layer   = hdf5_get_layer()
   zz        = hdf5_get_height()
   is_direct = hdf5_get_direct()
   IF(is_direct) THEN
      WRITE(zeile,'(''Currently at layer:'',i7,2x,'' w = '',F10.4)') n_layer, zz
   ELSE
      WRITE(zeile,'(''Currently at layer:'',i7,2x,'' l = '',F10.4)') n_layer, zz
   ENDIF
   CALL frame_menu 
   CALL draw_tframe (zeile(1:LEN_TRIM(zeile)), &
      'LEFT button/Keyboard l: down  / RIGHT button/Keyboard r: up' &
     &  , 'MIDDLE button/Keyboard m: back')                                         
!                                                                       
!------ First point                                                     
!                                                                       
!   10 CONTINUE 
   CALL open_viewport 
   CALL PGSCI (ilinecol (iwin, iframe, 0) ) 
   CALL PGBAND (0, 0, 0.0, 0.0, wx, wy, key) 
!  hub = nz (iwin, iframe, 1) * z_inc (iwin, iframe, 1) 
   IF (key == butt_m .OR. key==keyb_m) then 
      EXIT main_loop
   ELSEIF (key == butt_l .OR. key==keyb_l) then 
      CALL place_kuplot(h5_dims, -1, .FALSE., .FALSE., .FALSE.,                &
   MAXARRAY, MAXKURVTOT, fname, iz, x, y, z, nx, ny, xmin, xmax, ymin, ymax,     &
   offxy, offz, lni, lh5, ku_ndims, lenc, ier_num, ier_typ, output_io)
   ELSEIF (key == butt_r .OR. key==keyb_r) then 
      CALL place_kuplot(h5_dims,  1, .FALSE., .FALSE., .FALSE.,                &
   MAXARRAY, MAXKURVTOT, fname, iz, x, y, z, nx, ny, xmin, xmax, ymin, ymax,     &
   offxy, offz, lni, lh5, ku_ndims, lenc, ier_num, ier_typ, output_io)
   ENDIF 
   CALL draw_frame (iframe, lw) 
ENDDO main_loop
!                                                                       
CALL draw_frame (iframe, lw) 
CALL frame_menu 
!
n_layer   = hdf5_get_layer()
zz        = hdf5_get_height()
is_direct = hdf5_get_direct()
IF(is_direct) THEN
   WRITE (zeile, 1000)  n_layer, zz
   string = 'layer in res[1], height in res[2], direct==1 in res[3]'
ELSE
   WRITE (zeile, 1100)  n_layer, zz
   string = 'layer in res[1], height in res[2], reciprocal==0 in res[3]'
ENDIF
!                                                                       
CALL draw_tframe (zeile, ' ', string)
res_para(1) = REAL(n_layer)
res_para(2) = zz
IF(is_direct) THEN
   res_para(3) = 1.0
ELSE
   res_para(3) = 0.0
ENDIF
res_para(0) = 3.0
!                                                                       
 1000 FORMAT     ('New layer  : ',I7,2x, 'w = ', F10.4) 
 1100 FORMAT     ('New layer  : ',I7,2x, 'l = ', F10.4) 
!                                                                       
END SUBROUTINE do_layer                      
!
!******7****************************************************************
      SUBROUTINE do_fselect 
!                                                                       
!     Select active frame via mouse                                     
!     LEFT : select frame                                               
!                                                                       
      USE errlist_mod 
      USE kuplot_config 
      USE kuplot_mod 
use kuplot_draw_low_mod
use kuplot_draw_tframe_mod
use kuplot_plot_low_mod
!                                                                       
      IMPLICIT none 
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
      IF (.NOT. (key==butt_l .OR. key==keyb_l)) return 
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
      USE errlist_mod 
      USE kuplot_config 
      USE kuplot_mod 
use kuplot_draw_low_mod
use kuplot_draw_tframe_mod
use kuplot_plot_low_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER(80) zeile (2) 
      CHARACTER(1) key 
      REAL px (2), py (2) 
      REAL wx1, wy1, wx2, wy2 
      LOGICAL lw 
!                                                                       
      lw = .true. 
!                                                                       
      CALL draw_tframe (' ', ' ', 'LEFT button: mark corners / RIGHT botton: reset')                                                      
!                                                                       
!------ First point                                                     
!                                                                       
      CALL open_viewport 
      CALL PGSCI (ilinecol (iwin, iframe, 0) ) 
      CALL PGBAND (0, 0, 0.0, 0.0, wx1, wy1, key) 
      IF (key.eq.butt_r .OR. key==keyb_r) then 
         ex (iwin, iframe, 1) = - 9999.0 
         ey (iwin, iframe, 1) = - 9999.0 
         t (iwin, iframe, 1) = - 9999.0 
         t (iwin, iframe, 2) = - 9999.0 
      ENDIF 
      IF (.NOT. (key==butt_l .OR. key==keyb_l)) goto 9999 
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
      IF (.NOT. (key==butt_l .OR. key==keyb_l)) goto 9999
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
 1000 FORMAT     ('New range X: ',G13.6,' to ',G13.6) 
 1010 FORMAT     ('New range Y: ',G13.6,' to ',G13.6) 
!                                                                       
      END SUBROUTINE do_region                      
!******7****************************************************************
      SUBROUTINE do_zoom (lzoom) 
!                                                                       
!     Zoom in/out of picture                                            
!                                                                       
      USE errlist_mod 
      USE learn_mod 
      USE kuplot_config 
      USE kuplot_mod 
use kuplot_draw_low_mod
use kuplot_draw_tframe_mod
use kuplot_plot_low_mod
!                                                                       
      IMPLICIT none 
!                                                                       
LOGICAL , intent(in) :: lzoom
      CHARACTER(80) zeile (2), thelp 
      CHARACTER(1) key 
      REAL wx, wy, w2, h2 
      LOGICAL ::     lw 
!                                                                       
      lw = .true. 
!                                                                       
      IF (lzoom) then 
         thelp = 'LEFT button: zoom in / RIGHT button: zoom out'//' / MIDDLE&
     & botton: exit'                                                   
      ELSE 
         thelp = 'LEFT button: new center / MIDDLE botton: exit' 
      ENDIF 
      CALL draw_tframe (' ', ' ', thelp) 
!                                                                       
!------ Check mouse                                                     
!                                                                       
   10 CONTINUE 
      CALL open_viewport 
      CALL PGSCI (ilinecol (iwin, iframe, 0) ) 
      CALL PGBAND (0, 0, 0.0, 0.0, wx, wy, key) 
      IF (key.eq.butt_m.OR. key==keyb_m) goto 20 
!                                                                       
!------ - Compute new boundaries                                        
!                                                                       
      CALL trans_koor (wx, wy) 
!                                                                       
      w2 = 0.5 * (ex (iwin, iframe, 2) - ex (iwin, iframe, 1) ) 
      h2 = 0.5 * (ey (iwin, iframe, 2) - ey (iwin, iframe, 1) ) 
!                                                                       
      IF ((key.eq.butt_l.OR.key==keyb_l).and.lzoom) then 
         w2 = 0.9 * w2 
         h2 = 0.9 * h2 
      ENDIF 
      IF ((key.eq.butt_r.OR.key==keyb_r).and.lzoom) then 
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
 1000 FORMAT     ('New range X: ',G13.6,' to ',G13.6) 
 1010 FORMAT     ('New range Y: ',G13.6,' to ',G13.6) 
!                                                                       
      END SUBROUTINE do_zoom                        
!******7****************************************************************
      SUBROUTINE do_distance 
!                                                                       
!     Shows distance between two points ..                              
!                                                                       
      USE errlist_mod 
      USE param_mod 
      USE kuplot_config 
      USE kuplot_mod 
use kuplot_draw_low_mod
use kuplot_draw_tframe_mod
use kuplot_plot_low_mod
!                                                                       
      IMPLICIT none 
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
 2200 FORMAT    ('Distance = ',g13.6) 
!                                                                       
      END SUBROUTINE do_distance                    
!******7****************************************************************
      SUBROUTINE do_koor 
!                                                                       
!     Show coordinates of mouse pointer position                        
!                                                                       
      USE errlist_mod 
      USE param_mod 
      USE kuplot_config 
      USE kuplot_mod 
use kuplot_draw_low_mod
use kuplot_draw_tframe_mod
use kuplot_low_mod 
use kuplot_plot_low_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER(80) zeile 
      CHARACTER(1) key 
      REAL wx, wy, wz, dxx, dyy 
      INTEGER i, ini, nxx, nyy 
      LOGICAL l2d! , k_in_f 
!                                                                       
      zeile = ' ' 
      ini = 1
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
      IF (key.eq.butt_r .OR. key==keyb_r) goto 20 
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
         dxx = (xmax (ini) - xmin (ini) ) / REAL(nx (ini) - 1) 
         dyy = (ymax (ini) - ymin (ini) ) / REAL(ny (ini) - 1) 
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
 2200 FORMAT    ('(x,y) = ',2(g13.6,1x)) 
 2210 FORMAT    ('(x,y,z) = ',3(g13.6,1x)) 
!                                                                       
      END SUBROUTINE do_koor                        
!******7****************************************************************
      SUBROUTINE do_koor_nogui 
!                                                                       
!     Get coordinates of mouse pointer position outside GUI window      
!                                                                       
      USE errlist_mod 
      USE param_mod 
      USE prompt_mod 
      USE kuplot_config 
      USE kuplot_mod 
use kuplot_low_mod 
use kuplot_plot_mod
use kuplot_plot_low_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER(1) key 
      REAL wx, wy, wz 
      REAL dxx, dyy 
      INTEGER i, ini, ikey 
      INTEGER nxx, nyy 
      LOGICAL l2d! , k_in_f 
!                                                                       
      ini = 1
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
      IF (key.eq.butt_l.OR.key==keyb_l) ikey = 1 
      IF (key.eq.butt_m.OR.key==keyb_m) ikey = 2 
      IF (key.eq.butt_r.OR.key==keyb_r) ikey = 3 
!                                                                       
      CALL trans_koor (wx, wy) 
!                                                                       
      IF (.not.l2d) then 
         res_para (0) = 3 
         res_para (1) = wx 
         res_para (2) = wy 
         res_para (3) = REAL(ikey) 
      ELSE 
         dxx = (xmax (ini) - xmin (ini) ) / REAL(nx (ini) - 1) 
         dyy = (ymax (ini) - ymin (ini) ) / REAL(ny (ini) - 1) 
         nxx = nint ( ( (wx - xmin (ini) ) / dxx) ) + 1 
         nyy = nint ( ( (wy - ymin (ini) ) / dyy) ) + 1 
         wz = z (offz (ini - 1) + (nxx - 1) * ny (ini) + nyy) 
         res_para (0) = 4 
         res_para (1) = wx 
         res_para (2) = wy 
         res_para (3) = wz 
         res_para (4) = REAL(ikey) 
      ENDIF 
!                                                                       
 1000 FORMAT     (' ------ > Select desired point using the mouse ...') 
      END SUBROUTINE do_koor_nogui                  
!******7****************************************************************
      SUBROUTINE do_region_nogui (imode) 
!                                                                       
!     Set region via mouse - without menu                               
!                                                                       
      USE errlist_mod 
      USE param_mod 
      USE prompt_mod 
      USE kuplot_config 
      USE kuplot_mod 
use kuplot_plot_mod
use kuplot_plot_low_mod
!                                                                       
      IMPLICIT none 
!                                                                       
integer, intent(in) :: imode
!
      CHARACTER(1) key 
      REAL wx1, wy1, wx2, wy2 
      INTEGER ikey
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
      IF (key.eq.butt_l.OR.key==keyb_l) ikey = 1 
      IF (key.eq.butt_m.OR.key==keyb_l) ikey = 2 
      IF (key.eq.butt_r.OR.key==keyb_l) ikey = 3 
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
      res_para (5) = REAL(ikey) 
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
      USE wink_mod
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
REAL , intent(inout) :: wx, wy 
!                                                                       
      IF (shear (iwin, iframe) .ne.90.0) then 
         wx = wx - (wy - ey (iwin, iframe, 1) ) * yskal (iwin, iframe)  &
         / tan (REAL(rad) * shear (iwin, iframe) )
      ENDIF 
!                                                                       
      IF (lachse (iwin, iframe, 1) ) wx = exp (wx * log (10.0) ) 
      IF (lachse (iwin, iframe, 2) ) wy = exp (wy * log (10.0) ) 
!                                                                       
      END SUBROUTINE trans_koor                     
!******7****************************************************************
      SUBROUTINE PGRSTR (X, Y, ANGLE, FJUST, TEXT, LSTR, BCI) 
!
REAL , intent(inout) :: X, Y, ANGLE, FJUST 
CHARACTER (len= * ) , intent(inout) :: TEXT 
INTEGER , intent(inout) :: LSTR, BCI 
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
!
end module kuplot_draw_mod
