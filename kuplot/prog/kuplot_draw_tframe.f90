module kuplot_draw_tframe_mod
!
contains
!
!
!*******************************************************************************
!
      SUBROUTINE draw_tframe (tt1, tt2, tt3) 
!                                                                       
!     This routine writes the strings tt1-tt3 below plot                
!                                                                       
      USE kuplot_config 
      USE kuplot_mod 
USE lib_length
!                                                                       
      IMPLICIT none 
!                                                                       
CHARACTER(len=*), intent(in) :: tt1, tt2, tt3 
!
      CHARACTER(80) t1, t2, t3 
      REAL xt, yt, xh, yh 
      INTEGER ltext 
!                                                                       
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
!
!*******************************************************************************
!
!
      SUBROUTINE frame_menu
!                                                                       
!     This sets the menu viewport and draws borders ...                 
!                                                                       
      USE kuplot_config
      USE kuplot_mod
!                                                                       
      IMPLICIT none
!                                                                       
      CALL PGSVP (0.0, 0.999, 0.0, 0.999)
      CALL PGSWIN (0.0, 1.0, 0.0, 1.0)
!                                                                       
      CALL PGSFS (2)
      CALL PGSCI (6)
      CALL PGRECT (0.0, 1.0, 0.0, 1.0)
!                                                                       
      END SUBROUTINE frame_menu
!
!*******************************************************************************
!
end module kuplot_draw_tframe_mod
