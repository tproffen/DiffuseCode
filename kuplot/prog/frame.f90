!*****7*****************************************************************
!     These routines set the frame parameters for KUPLOT. Each frame    
!     can be seen as individual plotting area.                          
!*****7*****************************************************************
      SUBROUTINE set_window (zeile, lp) 
!+                                                                      
!     Sets active window for plotting                                   
!-                                                                      
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
      USE prompt_mod 
      USE kuplot_config 
      USE kuplot_mod 
USE precision_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 2) 
!                                                                       
      CHARACTER ( * ) zeile 
      INTEGER lp 
!                                                                       
      CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
      INTEGER lpara (maxw) 
      INTEGER ianz, iw 
      REAL(KIND=PREC_DP) :: werte (maxw) 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
!                                                                       
      IF (ianz.eq.0) then 
         WRITE (output_io, 1000) iwin 
      ELSEIF (ianz.eq.1.or.ianz.eq.2) then 
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         IF (ier_num.ne.0) return 
         iw = nint (werte (1) ) 
         IF (iw.lt.1.or.iw.gt.MAXWIN) then 
            ier_num = - 42 
            ier_typ = ER_APPL 
         ELSE 
            iwin = iw 
            IF (ianz.eq.2) dev_sf (iwin, x11) = werte (2) 
            IF (dev_id (iw, x11) .eq. - 1) call copy_window (1, iw) 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
 1000 FORMAT     (' ------ > Active window : ',i3) 
      END SUBROUTINE set_window                     
!*****7*****************************************************************
      SUBROUTINE set_afra (zeile, lp) 
!+                                                                      
!     Sets active frame for input parameters                            
!-                                                                      
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
      USE prompt_mod 
      USE kuplot_config 
      USE kuplot_mod 
USE precision_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 2) 
!                                                                       
      CHARACTER ( * ) zeile 
      INTEGER lp 
!                                                                       
      CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
      INTEGER lpara (maxw) 
      INTEGER ianz, ifr 
      REAL(KIND=PREC_DP) :: werte (maxw) 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
!                                                                       
      IF (ianz.eq.0) then 
         WRITE (output_io, 1000) iframe 
      ELSEIF (ianz.eq.1) then 
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         IF (ier_num.ne.0) return 
         ifr = nint (werte (1) ) 
         IF (ifr.lt.1.or.ifr.gt.iaf (iwin) ) then 
            ier_num = - 15 
            ier_typ = ER_APPL 
         ELSE 
            iframe = ifr 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
 1000 FORMAT     (' ------ > Active frame : ',i3) 
      END SUBROUTINE set_afra                       
!*****7*****************************************************************
      SUBROUTINE set_bfra (zeile, lp) 
!+                                                                      
!     Set backgound color for frames                                    
!-                                                                      
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
      USE kuplot_config 
      USE kuplot_mod 
USE precision_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 5) 
!                                                                       
      CHARACTER ( * ) zeile 
      INTEGER lp 
!                                                                       
      CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
      INTEGER lpara (maxw) 
      INTEGER ianz, ifr 
      REAL(KIND=PREC_DP) :: werte (maxw) 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
!                                                                       
      IF (ianz.eq.0) then 
         CALL show_frames 
      ELSEIF (ianz.eq.4) then 
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         IF (ier_num.ne.0) return 
         ifr = nint (werte (1) ) 
         IF (ifr.lt.1.or.ifr.gt.iaf (iwin) ) then 
            ier_num = - 15 
            ier_typ = ER_APPL 
         ELSE 
            IF (werte (2) .le.1.0.and.werte (2) .ge.0.0.and.werte (3)   &
            .le.1.0.and.werte (3) .ge.0.0.and.werte (4)                 &
            .le.1.0.and.werte (4) .ge.0.0) then                         
               frback (iwin, ifr, 1) = werte (2) 
               frback (iwin, ifr, 2) = werte (3) 
               frback (iwin, ifr, 3) = werte (4) 
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
      END SUBROUTINE set_bfra                       
!*****7*****************************************************************
      SUBROUTINE set_cfra (zeile, lp) 
!+                                                                      
!     Copy plot parameters of one frame to other frame(s)               
!-                                                                      
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
      USE kuplot_config 
      USE kuplot_mod 
USE precision_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = maxframe+2) 
!                                                                       
      CHARACTER ( * ) zeile 
      INTEGER lp 
!                                                                       
      CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
      INTEGER lpara (maxw) 
      INTEGER ianz, ifrom, ito, i 
      REAL(KIND=PREC_DP) :: werte (maxw) 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
!                                                                       
      IF (ianz.ge.2) then 
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         IF (ier_num.ne.0) return 
         ifrom = nint (werte (1) ) 
         DO i = 2, ianz 
         ito = nint (werte (i) ) 
         IF (ifrom.lt.1.or.ifrom.gt.iaf (iwin)                          &
         .or.ito.lt.1.or.ito.gt.iaf (iwin) .or.ifrom.eq.ito) then       
            ier_num = - 15 
            ier_typ = ER_APPL 
            RETURN 
         ELSE 
            CALL copy_frame (ifrom, ito) 
         ENDIF 
         ENDDO 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
      END SUBROUTINE set_cfra                       
!*****7*****************************************************************
      SUBROUTINE set_sfra (zeile, lp) 
!+                                                                      
!     Sets frame sizes and locations                                    
!-                                                                      
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
      USE kuplot_config 
      USE kuplot_mod 
USE precision_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 6) 
!                                                                       
      CHARACTER ( * ) zeile 
      INTEGER lp 
!                                                                       
      CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
      INTEGER lpara (maxw) 
      INTEGER ianz, ifr 
      REAL(KIND=PREC_DP) :: werte (maxw) 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
!                                                                       
      IF (ianz.eq.0) then 
         CALL show_frames 
      ELSEIF (ianz.eq.5) then 
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         IF (ier_num.ne.0) return 
!                                                                       
         ifr = nint (werte (1) ) 
         IF (ifr.lt.1.or.ifr.gt.iaf (iwin) ) then 
            ier_num = - 15 
            ier_typ = ER_APPL 
         ELSE 
            IF (werte (2) .lt.werte (4) .and.werte (3) .lt.werte (5) )  &
            then                                                        
               frame (iwin, ifr, 1) = werte (2) 
               frame (iwin, ifr, 2) = werte (3) 
               frame (iwin, ifr, 3) = werte (4) 
               frame (iwin, ifr, 4) = werte (5) 
            ELSE 
               ier_num = - 7 
               ier_typ = ER_APPL 
            ENDIF 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
      END SUBROUTINE set_sfra                       
!*****7*****************************************************************
      SUBROUTINE set_kfra (zeile, lp) 
!+                                                                      
!     Set datasets to be displayed in each frame                        
!     or give textfile for text frames.                                 
!-                                                                      
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
      USE kuplot_config 
      USE kuplot_mod 
!
      USE build_name_mod
USE precision_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = maxkurvtot + 2) 
!                                                                       
      CHARACTER ( * ) zeile 
      INTEGER lp 
!                                                                       
      CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
      INTEGER lpara (maxw) 
      INTEGER ianz, ifr, ikk, i 
      REAL(KIND=PREC_DP) :: werte (maxw) 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
!                                                                       
      IF (ianz.eq.0) then 
         CALL show_frames 
      ELSEIF (ianz.ge.2) then 
         CALL ber_params (1, cpara, lpara, werte, maxw) 
         ifr = nint (werte (1) ) 
         IF (ier_num.ne.0) return 
         IF (ifr.lt.1.or.ifr.gt.iaf (iwin) ) then 
            ier_num = - 15 
            ier_typ = ER_APPL 
         ELSE 
            CALL del_params (1, ianz, cpara, lpara, maxw) 
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
            IF (ier_num.eq.0) then 
               DO i = 1, maxkurvtot 
               infra (iwin, ifr, i) = 0 
               ENDDO 
               DO i = 1, ianz 
               ikk = nint (werte (i) ) 
               IF (ikk.lt.1.or.ikk.gt. (iz - 1) ) then 
                  ier_num = - 4 
                  ier_typ = ER_APPL 
               ELSE 
                  infra (iwin, ifr, i) = ikk 
               ENDIF 
               ENDDO 
            ELSE 
               CALL no_error 
               CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
               IF (ier_num.ne.0) return 
               infra (iwin, ifr, 1) = - 1 
               ftext (iwin, ifr) = cpara (1) (1:MIN(40,LEN_TRIM(cpara(1))))
            ENDIF 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
      END SUBROUTINE set_kfra                       
!*****7*****************************************************************
      SUBROUTINE set_nfra (zeile, lp) 
!+                                                                      
!     Set mumber of frames to be used                                   
!-                                                                      
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
      USE prompt_mod 
      USE kuplot_config 
      USE kuplot_mod 
USE precision_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 2) 
!                                                                       
      CHARACTER ( * ) zeile 
      INTEGER lp 
!                                                                       
      CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
      INTEGER lpara (maxw) 
      INTEGER ianz, ifr, idummy 
      REAL(KIND=PREC_DP) :: werte (maxw) 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
!                                                                       
      IF (ianz.eq.0) then 
         CALL show_frames 
      ELSEIF (ianz.eq.1) then 
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         idummy = nint (werte (1) ) 
         IF (ier_num.ne.0) return 
!                                                                       
         IF (idummy.le.0.or.idummy.gt.maxframe) then 
            ier_num = - 17 
            ier_typ = ER_APPL 
         ELSE 
            iaf (iwin) = idummy 
            iframe = 1 
            DO ifr = 1, iaf (iwin) 
            CALL copy_frame (1, ifr) 
            ENDDO 
!                                                                       
!-----------default frames                                              
!                                                                       
            IF (iaf (iwin) .eq.1) then 
               frame (iwin, 1, 1) = 0.0 
               frame (iwin, 1, 2) = 0.0 
               frame (iwin, 1, 3) = 1.0 
               frame (iwin, 1, 4) = 1.0 
!                                                                       
            ELSEIF (iaf (iwin) .eq.2) then 
               frame (iwin, 1, 1) = 0.0 
               frame (iwin, 1, 2) = 0.0 
               frame (iwin, 1, 3) = 0.5 
               frame (iwin, 1, 4) = 1.0 
!                                                                       
               frame (iwin, 2, 1) = 0.5 
               frame (iwin, 2, 2) = 0.0 
               frame (iwin, 2, 3) = 1.0 
               frame (iwin, 2, 4) = 1.0 
!                                                                       
            ELSEIF (iaf (iwin) .eq.3) then 
               frame (iwin, 1, 1) = 0.0 
               frame (iwin, 1, 2) = 0.0 
               frame (iwin, 1, 3) = 0.333 
               frame (iwin, 1, 4) = 1.0 
!                                                                       
               frame (iwin, 2, 1) = 0.333 
               frame (iwin, 2, 2) = 0.0 
               frame (iwin, 2, 3) = 0.666 
               frame (iwin, 2, 4) = 1.0 
!                                                                       
               frame (iwin, 3, 1) = 0.666 
               frame (iwin, 3, 2) = 0.0 
               frame (iwin, 3, 3) = 1.0 
               frame (iwin, 3, 4) = 1.0 
!                                                                       
            ELSEIF (iaf (iwin) .eq.4) then 
               frame (iwin, 1, 1) = 0.0 
               frame (iwin, 1, 2) = 0.0 
               frame (iwin, 1, 3) = 0.5 
               frame (iwin, 1, 4) = 0.5 
!                                                                       
               frame (iwin, 2, 1) = 0.5 
               frame (iwin, 2, 2) = 0.0 
               frame (iwin, 2, 3) = 1.0 
               frame (iwin, 2, 4) = 0.5 
!                                                                       
               frame (iwin, 3, 1) = 0.0 
               frame (iwin, 3, 2) = 0.5 
               frame (iwin, 3, 3) = 0.5 
               frame (iwin, 3, 4) = 1.0 
!                                                                       
               frame (iwin, 4, 1) = 0.5 
               frame (iwin, 4, 2) = 0.5 
               frame (iwin, 4, 3) = 1.0 
               frame (iwin, 4, 4) = 1.0 
            ELSEIF (iaf (iwin) .eq.5.or.iaf (iwin) .eq.6) then 
               frame (iwin, 1, 1) = 0.0 
               frame (iwin, 1, 2) = 0.0 
               frame (iwin, 1, 3) = 0.333 
               frame (iwin, 1, 4) = 0.5 
!                                                                       
               frame (iwin, 2, 1) = 0.333 
               frame (iwin, 2, 2) = 0.0 
               frame (iwin, 2, 3) = 0.666 
               frame (iwin, 2, 4) = 0.5 
!                                                                       
               frame (iwin, 3, 1) = 0.666 
               frame (iwin, 3, 2) = 0.0 
               frame (iwin, 3, 3) = 1.0 
               frame (iwin, 3, 4) = 0.5 
!                                                                       
               frame (iwin, 4, 1) = 0.0 
               frame (iwin, 4, 2) = 0.5 
               frame (iwin, 4, 3) = 0.333 
               frame (iwin, 4, 4) = 1.0 
!                                                                       
               frame (iwin, 5, 1) = 0.333 
               frame (iwin, 5, 2) = 0.5 
               frame (iwin, 5, 3) = 0.666 
               frame (iwin, 5, 4) = 1.0 
!                                                                       
               frame (iwin, 6, 1) = 0.666 
               frame (iwin, 6, 2) = 0.5 
               frame (iwin, 6, 3) = 1.0 
               frame (iwin, 6, 4) = 1.0 
            ELSE 
               WRITE (output_io, 2000) iaf (iwin) 
            ENDIF 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
 2000 FORMAT     (' ****WARN****No default for ',i2,                    &
     &                   ' frames, use sfra ..')                        
      END SUBROUTINE set_nfra                       
!*****7*****************************************************************
      SUBROUTINE copy_frame (if1, if2) 
!+                                                                      
!     Copy parameters associated with frame if1 to frame if2            
!-                                                                      
      USE kuplot_config 
      USE kuplot_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER if1, if2 
      INTEGER i, j, ihl 
!                                                                       
      achse (iwin, if2, 1) = achse (iwin, if1, 1) 
      achse (iwin, if2, 2) = achse (iwin, if1, 2) 
      achse (iwin, if2, 3) = achse (iwin, if1, 3) 
      lachse (iwin, if2, 1) = lachse (iwin, if1, 1) 
      lachse (iwin, if2, 2) = lachse (iwin, if1, 2) 
      lachse (iwin, if2, 3) = lachse (iwin, if1, 3) 
      titel (iwin, if2, 1) = titel (iwin, if1, 1) 
      titel (iwin, if2, 2) = titel (iwin, if1, 2) 
      ftext (iwin, if2) = ftext (iwin, if1) 
      yskal_u (iwin, if2) = yskal_u (iwin, if1) 
      lyskal (iwin, if2) = lyskal (iwin, if1) 
      ex (iwin, if2, 1) = ex (iwin, if1, 1) 
      ex (iwin, if2, 2) = ex (iwin, if1, 2) 
      ey (iwin, if2, 1) = ey (iwin, if1, 1) 
      ey (iwin, if2, 2) = ey (iwin, if1, 2) 
      t (iwin, if2, 1) = t (iwin, if1, 1) 
      t (iwin, if2, 2) = t (iwin, if1, 2) 
      shear (iwin, if2) = shear (iwin, if1) 
      sfl (iwin, if2) = sfl (iwin, if1) 
      igrid (iwin, if2) = igrid (iwin, if1) 
      ibox (iwin, if2) = ibox (iwin, if1) 
      ifname (iwin, if2) = ifname (iwin, if1) 
      frjust (iwin, if2) = frjust (iwin, if1) 
      frback (iwin, if2, 1) = frback (iwin, if1, 1) 
      frback (iwin, if2, 2) = frback (iwin, if1, 2) 
      frback (iwin, if2, 3) = frback (iwin, if1, 3) 
      fonscal (iwin, if2) = fonscal (iwin, if1) 
!                                                                       
      lab_angle (iwin, if2, 1) = lab_angle (iwin, if1, 1) 
      tick_ma_h (iwin, if2, 1) = tick_ma_h (iwin, if1, 1) 
      tick_mi_h (iwin, if2, 1) = tick_mi_h (iwin, if1, 1) 
      tick_nsub (iwin, if2, 1) = tick_nsub (iwin, if1, 1) 
!                                                                       
      DO j = 1, maxbond 
      bond_rad (iwin, if2, j) = bond_rad (iwin, if1, j) 
      bond_sig (iwin, if2, j) = bond_sig (iwin, if1, j) 
      bond_lwid (iwin, if2, j) = bond_lwid (iwin, if1, j) 
      bond_lcol (iwin, if2, j) = bond_lcol (iwin, if1, j) 
      bond_ltyp (iwin, if2, j) = bond_ltyp (iwin, if1, j) 
      ENDDO 
!                                                                       
      DO j = 1, 4 
      ibuf (iwin, if2, j) = ibuf (iwin, if1, j) 
      ENDDO 
!                                                                       
      DO j = 1, nfon 
      fonsize (iwin, if2, j) = fonsize (iwin, if1, j) 
      fon_id (iwin, if2, j) = fon_id (iwin, if1, j) 
      foncol (iwin, if2, j) = foncol (iwin, if1, j) 
      ENDDO 
!                                                                       
      DO j = 1, maxan 
      antext (iwin, if2, j) = antext (iwin, if1, j) 
      anx (iwin, if2, j) = anx (iwin, if1, j) 
      any (iwin, if2, j) = any (iwin, if1, j) 
      antx (iwin, if2, j) = antx (iwin, if1, j) 
      anty (iwin, if2, j) = anty (iwin, if1, j) 
      anjust (iwin, if2, j) = anjust (iwin, if1, j) 
      anangle (iwin, if2, j) = anangle (iwin, if1, j) 
      ENDDO 
!                                                                       
      iho (iwin, if2) = iho (iwin, if1) 
      DO ihl = 1, maxhl 
      nz (iwin, if2, ihl) = nz (iwin, if1, ihl) 
      z_min (iwin, if2, ihl) = z_min (iwin, if1, ihl) 
      z_inc (iwin, if2, ihl) = z_inc (iwin, if1, ihl) 
      ENDDO 
!                                                                       
      linewid (iwin, if2, 0) = linewid (iwin, if1, 0) 
      ilinecol (iwin, if2, 0) = ilinecol (iwin, if1, 0) 
      DO i = 1, maxkurvtot 
      clegend (iwin, if2, i) = clegend (iwin, if1, i) 
      ilegend (iwin, if2, i) = ilegend (iwin, if1, i) 
      sizemark (iwin, if2, i) = sizemark (iwin, if1, i) 
      linewid (iwin, if2, i) = linewid (iwin, if1, i) 
      ilinecol (iwin, if2, i) = ilinecol (iwin, if1, i) 
      ilinetyp (iwin, if2, i) = ilinetyp (iwin, if1, i) 
      ifillcol (iwin, if2, i) = ifillcol (iwin, if1, i) 
      ifilltyp (iwin, if2, i) = ifilltyp (iwin, if1, i) 
      ilineart (iwin, if2, i) = ilineart (iwin, if1, i) 
      imarkcol (iwin, if2, i) = imarkcol (iwin, if1, i) 
      imarktyp (iwin, if2, i) = imarktyp (iwin, if1, i) 
      imarkmax (iwin, if2, i) = imarkmax (iwin, if1, i) 
      ierr (iwin, if2, i) = ierr (iwin, if1, i) 
      ierrcol (iwin, if2, i) = ierrcol (iwin, if1, i) 
      hlineart (iwin, if2, i) = hlineart (iwin, if1, i) 
      hlabel (iwin, if2, i) = hlabel (iwin, if1, i) 
      DO ihl = 1, maxhl 
      hlinecol (iwin, if2, i, ihl) = hlinecol (iwin, if1, i, ihl) 
      hlinetyp (iwin, if2, i, ihl) = hlinetyp (iwin, if1, i, ihl) 
      ENDDO 
      DO ihl = 1, 4 
      fillrange (iwin, if2, i, ihl) = fillrange (iwin, if1, i, ihl) 
      ENDDO 
      ENDDO 
!                                                                       
      info_orig (iwin, if2, 1) = info_orig (iwin, if1, 1) 
      info_orig (iwin, if2, 2) = info_orig (iwin, if1, 2) 
!                                                                       
      END SUBROUTINE copy_frame                     
!*****7*****************************************************************
      SUBROUTINE copy_window (iw1, iw2) 
!+                                                                      
!     Copy parameters associated with window iw1 to frame iw2           
!-                                                                      
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER iw1, iw2 
      INTEGER i, j, k 
!                                                                       
      ipara (iw2) = ipara (iw1) 
      fpara (iw2) = fpara (iw1) 
      orient (iw2) = orient (iw1) 
      tot_frame (iw2) = tot_frame (iw1) 
      iaf (iw2) = iaf (iw1) 
      icol_old (iw2) = icol_old (iw1) 
      lgui (iw2) = lgui (iw1) 
!                                                                       
      DO i = 1, maxframe 
      ibuf (iw2, i, 1) = ibuf (iw1, i, 1) 
      ibuf (iw2, i, 2) = ibuf (iw1, i, 2) 
      ibuf (iw2, i, 3) = ibuf (iw1, i, 3) 
      ibuf (iw2, i, 4) = ibuf (iw1, i, 4) 
      linewid (iw2, i, 0) = linewid (iw1, i, 0) 
      ilinecol (iw2, i, 0) = ilinecol (iw1, i, 0) 
      iho (iw2, i) = iho (iw1, i) 
!                                                                       
      info_orig (iw2, i, 1) = info_orig (iw1, i, 1) 
      info_orig (iw2, i, 2) = info_orig (iw1, i, 2) 
!                                                                       
      DO k = 1, maxbond 
      bond_rad (iw2, i, k) = bond_rad (iw1, i, k) 
      bond_sig (iw2, i, k) = bond_sig (iw1, i, k) 
      bond_lwid (iw2, i, k) = bond_lwid (iw1, i, k) 
      bond_lcol (iw2, i, k) = bond_lcol (iw1, i, k) 
      bond_ltyp (iw2, i, k) = bond_ltyp (iw1, i, k) 
      ENDDO 
!                                                                       
      DO j = 1, maxkurvtot 
      DO k = 1, maxhl 
      hlinecol (iw2, i, j, k) = hlinecol (iw1, i, j, k) 
      hlinetyp (iw2, i, j, k) = hlinetyp (iw1, i, j, k) 
      ENDDO 
      DO k = 1, 4 
      fillrange (iw2, i, j, k) = fillrange (iw1, i, j, k) 
      ENDDO 
      hlabel (iw2, i, j) = hlabel (iw1, i, j) 
      linewid (iw2, i, j) = linewid (iw1, i, j) 
      ilegend (iw2, i, j) = ilegend (iw1, i, j) 
      ilinecol (iw2, i, j) = ilinecol (iw1, i, j) 
      ifillcol (iw2, i, j) = ifillcol (iw1, i, j) 
      imarkcol (iw2, i, j) = imarkcol (iw1, i, j) 
      ierrcol (iw2, i, j) = ierrcol (iw1, i, j) 
      ierr (iw2, i, j) = ierr (iw1, i, j) 
      ilinetyp (iw2, i, j) = ilinetyp (iw1, i, j) 
      ifilltyp (iw2, i, j) = ifilltyp (iw1, i, j) 
      imarktyp (iw2, i, j) = imarktyp (iw1, i, j) 
      imarkmax (iw2, i, j) = imarkmax (iw1, i, j) 
      ilineart (iw2, i, j) = ilineart (iw1, i, j) 
      hlineart (iw2, i, j) = hlineart (iw1, i, j) 
      sizemark (iw2, i, j) = sizemark (iw1, i, j) 
      infra (iw2, i, j) = infra (iw1, i, j) 
      ENDDO 
                                                                        
      DO j = 1, maxan 
      antext (iw2, i, j) = antext (iw1, i, j) 
      anx (iw2, i, j) = anx (iw1, i, j) 
      any (iw2, i, j) = any (iw1, i, j) 
      antx (iw2, i, j) = antx (iw1, i, j) 
      anty (iw2, i, j) = anty (iw1, i, j) 
      anjust (iw2, i, j) = anjust (iw1, i, j) 
      anangle (iw2, i, j) = anangle (iw1, i, j) 
      ENDDO 
!                                                                       
      DO j = 1, maxhl 
      z_min (iw2, i, j) = z_min (iw1, i, j) 
      z_inc (iw2, i, j) = z_inc (iw1, i, j) 
      nz (iw2, i, j) = nz (iw1, i, j) 
      ENDDO 
!                                                                       
      frjust (iw2, i) = frjust (iw1, i) 
      frback (iw2, i, 1) = frback (iw1, i, 1) 
      frback (iw2, i, 2) = frback (iw1, i, 2) 
      frback (iw2, i, 3) = frback (iw1, i, 3) 
!                                                                       
      fonscal (iw2, i) = fonscal (iw1, i) 
!                                                                       
      DO j = 1, 6 
      fonsize (iw2, i, j) = fonsize (iw1, i, j) 
      fon_id (iw2, i, j) = fon_id (iw1, i, j) 
      foncol (iw2, i, j) = foncol (iw1, i, j) 
      ENDDO 
!                                                                       
      DO j = 0, 20 
      colour (iw2, j, 1) = colour (iw1, j, 1) 
      colour (iw2, j, 2) = colour (iw1, j, 2) 
      colour (iw2, j, 3) = colour (iw1, j, 3) 
      ENDDO 
!                                                                       
      DO j = 1, maxcol 
      col_map (iw2, j, 1) = col_map (iw1, j, 1) 
      col_map (iw2, j, 2) = col_map (iw1, j, 2) 
      col_map (iw2, j, 3) = col_map (iw1, j, 3) 
      ENDDO 
!                                                                       
      ifname (iw2, i) = ifname (iw1, i) 
      igrid (iw2, i) = igrid (iw1, i) 
      ibox (iw2, i) = ibox (iw1, i) 
      ftext (iw2, i) = ftext (iw1, i) 
      ex (iw2, i, 1) = ex (iw1, i, 1) 
      ex (iw2, i, 2) = ex (iw1, i, 2) 
      ey (iw2, i, 1) = ey (iw1, i, 1) 
      ey (iw2, i, 2) = ey (iw1, i, 2) 
      t (iw2, i, 1) = t (iw1, i, 1) 
      t (iw2, i, 2) = t (iw1, i, 2) 
      yskal_u (iw2, i) = yskal_u (iw1, i) 
      lyskal (iw2, i) = lyskal (iw1, i) 
!                                                                       
      lab_angle (iw2, i, 1) = lab_angle (iw1, i, 1) 
      tick_ma_h (iw2, i, 1) = tick_ma_h (iw1, i, 1) 
      tick_mi_h (iw2, i, 1) = tick_mi_h (iw1, i, 1) 
      tick_nsub (iw2, i, 1) = tick_nsub (iw1, i, 1) 
!                                                                       
      titel (iw2, i, 1) = titel (iw1, i, 1) 
      titel (iw2, i, 2) = titel (iw1, i, 2) 
      achse (iw2, i, 1) = achse (iw1, i, 1) 
      achse (iw2, i, 2) = achse (iw1, i, 2) 
      achse (iw2, i, 3) = achse (iw1, i, 3) 
      lachse (iw2, i, 1) = lachse (iw1, i, 1) 
      lachse (iw2, i, 2) = lachse (iw1, i, 2) 
      lachse (iw2, i, 3) = lachse (iw1, i, 3) 
      shear (iw2, i) = shear (iw1, i) 
      sfl (iw2, i) = sfl (iw1, i) 
      frame (iw2, i, 1) = frame (iw1, i, 1) 
      frame (iw2, i, 2) = frame (iw1, i, 2) 
      frame (iw2, i, 3) = frame (iw1, i, 3) 
      frame (iw2, i, 4) = frame (iw1, i, 4) 
      ENDDO 
!                                                                       
      END SUBROUTINE copy_window                    
