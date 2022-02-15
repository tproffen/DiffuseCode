module kuplot_draw_low_mod
!
!-
!  Low level mouse and menu related routines that do not call other routines
!  besides PGPLOT
!+
!
contains
!
!*******************************************************************************
!
      SUBROUTINE draw_menu 
!                                                                       
!     Draw the different menu parts                                     
!                                                                       
      USE errlist_mod 
      USE prompt_mod 
      USE kuplot_config 
      USE kuplot_mod 
use kuplot_draw_tframe_mod
use kuplot_low_mod
USE lib_length
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER(40) text 
      REAL xh, yh, xt, yt 
      INTEGER ini 
!                                                                       
!      LOGICAL n_in_f 
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
      CALL def_button (2, 'Coordinates',    xt, 0.95, yt, 0.04) 
      CALL def_button (8, 'Distances',      xt, 0.90, yt, 0.04) 
      CALL def_button (3, 'Select region',  xt, 0.80, yt, 0.04) 
      CALL def_button (4, 'Zoom',           xt, 0.75, yt, 0.04) 
      CALL def_button (5, 'Move',           xt, 0.70, yt, 0.04) 
      CALL def_button (6, 'Reset region',   xt, 0.65, yt, 0.04) 
      CALL def_button (9, 'Select zmin',    xt, 0.55, yt, 0.04) 
      CALL def_button (10, 'Select zmax',   xt, 0.50, yt, 0.04) 
      CALL def_button (12, 'Select layer',  xt, 0.60, yt, 0.04) 
      CALL def_button (7, 'Select frame',   xt, 0.40, yt, 0.04) 
      CALL def_button (11, 'Enter command', xt, 0.25, yt, 0.04) 
      CALL def_button (1, 'Exit menu',      xt, 0.20, yt, 0.04) 
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
      CALL draw_button (12, .not.(n_in_f (ini).AND.lh5(iz-1)), .false.) 
!                                                                       
      CALL draw_tframe (' ', ' ', 'Use the EXIT button to return to comm&
     &and mode ...')                                                    
!                                                                       
      END SUBROUTINE draw_menu                      
      SUBROUTINE def_button (ib, text, xb, yb, wb, hb) 
!                                                                       
!     Defines a button for the menu                                     
!                                                                       
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
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
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
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
!
!*******************************************************************************
!
end module kuplot_draw_low_mod
