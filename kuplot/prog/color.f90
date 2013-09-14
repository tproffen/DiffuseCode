!*****7*****************************************************************
!     Here are all the routines for color handling for the              
!     KUPLOT bitmaps.                                                   
!*****7*****************************************************************
      SUBROUTINE set_color (zeile, lp) 
!+                                                                      
!     Set colours for background & pens                                 
!-                                                                      
      USE errlist_mod 
      USE config_mod 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 5) 
!                                                                       
      CHARACTER ( * ) zeile 
      INTEGER lp 
!                                                                       
      CHARACTER(1024) cpara (maxw) 
      INTEGER lpara (maxw) 
      INTEGER ianz, icol 
      REAL werte (maxw) 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
!                                                                       
      IF (ianz.eq.0) then 
         CALL show_color 
      ELSEIF (ianz.eq.4) then 
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         IF (ier_num.ne.0) return 
         icol = nint (werte (1) ) 
         IF (icol.lt.0.or.icol.gt.18) then 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ELSE 
            IF (werte (2) .le.1.0.and.werte (2) .ge.0.0.and.werte (3)   &
            .le.1.0.and.werte (3) .ge.0.0.and.werte (4)                 &
            .le.1.0.and.werte (4) .ge.0.0) then                         
               colour (iwin, icol, 1) = werte (2) 
               colour (iwin, icol, 2) = werte (3) 
               colour (iwin, icol, 3) = werte (4) 
            ELSE 
               ier_num = - 18 
               ier_typ = ER_APPL 
            ENDIF 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
      END SUBROUTINE set_color                      
!*****7*****************************************************************
      SUBROUTINE set_cmap (zeile, lp) 
!                                                                       
!     Set colour map                                                    
!                                                                       
      USE errlist_mod 
      USE config_mod 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 5) 
!                                                                       
      CHARACTER ( * ) zeile 
      CHARACTER(1024) cpara (maxw) 
      REAL werte (maxw) 
      INTEGER lpara (maxw) 
      INTEGER ianz, lp 
      LOGICAL str_comp 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
!                                                                       
      IF (ianz.ge.1) then 
         IF (str_comp (cpara (1) , 'gray', 3, lpara (1) , 4) ) then 
            CALL cmap_gray (.true.) 
         ELSEIF (str_comp (cpara (1) , 'fire', 3, lpara (1) , 4) ) then 
            CALL cmap_fire (.true.) 
         ELSEIF (str_comp (cpara (1) , 'kupl', 3, lpara (1) , 4) ) then 
            CALL cmap_kupl (.true.) 
         ELSEIF (str_comp (cpara (1) , 'invert', 3, lpara (1) , 6) )    &
         then                                                           
            CALL cmap_invert (.true.) 
         ELSEIF (str_comp (cpara (1) , 'read', 3, lpara (1) , 4) ) then 
            CALL do_build_name (ianz, cpara, lpara, werte, maxw, 2) 
            IF (ier_num.eq.0) call cmap_read (cpara (2) ) 
         ELSEIF (str_comp (cpara (1) , 'write', 3, lpara (1) , 5) )     &
         then                                                           
            CALL do_build_name (ianz, cpara, lpara, werte, maxw, 2) 
            IF (ier_num.eq.0) call cmap_write (cpara (2) ) 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
      END SUBROUTINE set_cmap                       
!*****7*****************************************************************
      SUBROUTINE cmap_invert (lout) 
!                                                                       
!     Setting colormap to gray                                          
!                                                                       
      USE prompt_mod 
      USE config_mod 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      REAL dummy (maxcol, 3) 
      INTEGER i 
      LOGICAL lout 
!                                                                       
      IF (lout) write (output_io, 1000) 
!                                                                       
      DO i = 1, maxcol 
      dummy (i, 1) = col_map (iwin, i, 1) 
      dummy (i, 2) = col_map (iwin, i, 2) 
      dummy (i, 3) = col_map (iwin, i, 3) 
      ENDDO 
!                                                                       
      DO i = 1, maxcol 
      col_map (iwin, i, 1) = dummy (maxcol - i + 1, 1) 
      col_map (iwin, i, 2) = dummy (maxcol - i + 1, 2) 
      col_map (iwin, i, 3) = dummy (maxcol - i + 1, 3) 
      ENDDO 
!                                                                       
 1000 FORMAT     (' ------ > Inverting current colour map ... ') 
!                                                                       
      END SUBROUTINE cmap_invert                    
!*****7*****************************************************************
      SUBROUTINE cmap_gray (lout) 
!                                                                       
!     Setting colormap to gray                                          
!                                                                       
      USE prompt_mod 
      USE config_mod 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER i 
      LOGICAL lout 
!                                                                       
      IF (lout) write (output_io, 1000) 
!                                                                       
      DO i = 1, maxcol 
      col_map (iwin, i, 1) = float (maxcol - i + 1) / float (maxcol) 
      col_map (iwin, i, 2) = float (maxcol - i + 1) / float (maxcol) 
      col_map (iwin, i, 3) = float (maxcol - i + 1) / float (maxcol) 
      ENDDO 
!                                                                       
 1000 FORMAT     (' ------ > Setting colour map to : gray ') 
      END SUBROUTINE cmap_gray                      
!*****7*****************************************************************
      SUBROUTINE cmap_kupl (lout) 
!                                                                       
!     Setting colormap to original KUPL map                             
!                                                                       
      USE prompt_mod 
      USE config_mod 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      REAL rh, rf, rq, rp, rt 
      INTEGER i, ifarb 
      LOGICAL lout 
!                                                                       
      IF (lout) write (output_io, 1000) 
!                                                                       
      DO ifarb = 1, maxcol 
      rh = 0.1 + float (ifarb - 1) / 284.0 
      rh = 6.0 * rh 
      i = int (rh) 
      rf = rh - float (i) 
      rp = 0.0 
      rq = 1.0 - rf 
      rt = (1.0 - (1.0 - rf) ) 
!                                                                       
      IF (rt.gt.1.0) rt = 1.0 
      IF (rp.gt.1.0) rp = 1.0 
      IF (rq.gt.1.0) rq = 1.0 
!                                                                       
      IF (i.eq.3) then 
         col_map (iwin, ifarb, 1) = 1. 
         col_map (iwin, ifarb, 2) = rt 
         col_map (iwin, ifarb, 3) = rp 
      ELSEIF (i.eq.4) then 
         col_map (iwin, ifarb, 1) = rq 
         col_map (iwin, ifarb, 2) = 1. 
         col_map (iwin, ifarb, 3) = rp 
      ELSEIF (i.eq.5) then 
         col_map (iwin, ifarb, 1) = rp 
         col_map (iwin, ifarb, 2) = 1. 
         col_map (iwin, ifarb, 3) = rt 
      ELSEIF (i.eq.0) then 
         col_map (iwin, ifarb, 1) = rp 
         col_map (iwin, ifarb, 2) = rq 
         col_map (iwin, ifarb, 3) = 1. 
      ELSEIF (i.eq.1) then 
         col_map (iwin, ifarb, 1) = rt 
         col_map (iwin, ifarb, 2) = rp 
         col_map (iwin, ifarb, 3) = 1. 
      ELSEIF (i.eq.2) then 
         col_map (iwin, ifarb, 1) = 1. 
         col_map (iwin, ifarb, 2) = rp 
         col_map (iwin, ifarb, 3) = rq 
      ENDIF 
      ENDDO 
!                                                                       
 1000 FORMAT     (' ------ > Setting colour map to : kupl ') 
      END SUBROUTINE cmap_kupl                      
!*****7*****************************************************************
      SUBROUTINE cmap_fire (lout) 
!                                                                       
!     Setting colormap to fire                                          
!                                                                       
      USE prompt_mod 
      USE config_mod 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER i, ii, m 
      LOGICAL lout 
!                                                                       
      IF (lout) write (output_io, 1000) 
!                                                                       
      m = maxcol / 3 
      ii = 1 
!                                                                       
      DO i = 1, m 
      col_map (iwin, ii, 1) = float (3 * i) / float (maxcol) 
      col_map (iwin, ii, 2) = 0.0 
      col_map (iwin, ii, 3) = 0.0 
      ii = ii + 1 
      ENDDO 
!                                                                       
      DO i = 1, m 
      col_map (iwin, ii, 1) = 1.0 
      col_map (iwin, ii, 2) = float (3 * i) / float (maxcol) 
      col_map (iwin, ii, 3) = 0.0 
      ii = ii + 1 
      ENDDO 
!                                                                       
      DO i = 1, m 
      col_map (iwin, ii, 1) = 1.0 
      col_map (iwin, ii, 2) = 1.0 
      col_map (iwin, ii, 3) = float (3 * i) / float (maxcol) 
      ii = ii + 1 
      ENDDO 
!                                                                       
 1000 FORMAT     (' ------ > Setting colour map to : fire ') 
      END SUBROUTINE cmap_fire                      
!*****7*****************************************************************
      SUBROUTINE cmap_read (filname) 
!                                                                       
!     Reading colormap from file ..                                     
!                                                                       
      USE errlist_mod 
      USE prompt_mod 
      USE config_mod 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER ( * ) filname 
      INTEGER i, ir, ig, ib 
!                                                                       
      CALL oeffne (33, filname, 'old', .false.) 
      IF (ier_num.eq.0) then 
         WRITE (output_io, 1000) filname 
!                                                                       
         DO i = 1, maxcol 
         READ (33, 2000, err = 100, end = 200) ir, ig, ib 
         col_map (iwin, i, 1) = ir / 255. 
         col_map (iwin, i, 2) = ig / 255. 
         col_map (iwin, i, 3) = ib / 255. 
         ENDDO 
!                                                                       
         CLOSE (33) 
      ENDIF 
      RETURN 
!                                                                       
  100 ier_num = - 3 
      ier_typ = ER_IO 
      CLOSE (33) 
      RETURN 
!                                                                       
  200 ier_num = - 6 
      ier_typ = ER_IO 
      CLOSE (33) 
      RETURN 
!                                                                       
 1000 FORMAT     (' ------ > Reading colour map from file : ',A30) 
 2000 FORMAT     (1x,3z2) 
      END SUBROUTINE cmap_read                      
!*****7*****************************************************************
      SUBROUTINE cmap_write (filname) 
!                                                                       
!     Writing colormap to file ..                                       
!                                                                       
      USE errlist_mod 
      USE prompt_mod 
      USE config_mod 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER ( * ) filname 
      INTEGER i, ir, ig, ib 
!                                                                       
      CALL oeffne (33, filname, 'unknown', .false.) 
      IF (ier_num.eq.0) then 
         WRITE (output_io, 1000) filname 
!                                                                       
         DO i = 1, maxcol 
         ir = int (col_map (iwin, i, 1) * 255.) 
         ig = int (col_map (iwin, i, 2) * 255.) 
         ib = int (col_map (iwin, i, 3) * 255.) 
         WRITE (33, 2000) ir, ig, ib 
         ENDDO 
!                                                                       
         CLOSE (33) 
      ENDIF 
!                                                                       
 1000 FORMAT     (' ------ > Writing colour map to file : ',A30) 
 2000 FORMAT     ('#',3z2.2) 
      END SUBROUTINE cmap_write                     
