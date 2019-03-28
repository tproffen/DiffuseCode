!*****7*****************************************************************
!                                                                       
      SUBROUTINE kuplot_initarrays 
!+                                                                      
!     This is the initialisation routine called at startup time         
!     of KUPLOT. In the first section there are a few definitions       
!     that might need adjustment.                                       
!+                                                                      
      USE nexus_kuplot
!
      USE envir_mod 
      USE errlist_mod 
      USE learn_mod
      USE prompt_mod 
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER i, j, k, iw 
      INTEGER stift (0:14) 
!                                                                       
      DATA stift / 6, 3, 1, 2, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 6 / 
!                                                                       
      CALL do_nxinit 
!                                                                       
!------ Sockt defaults                                                  
!                                                                       
      s_port = 3331 
      s_ipallowed = "localhost" 
!                                                                       
!------ Device definitions for PGPLOT                                   
!                                                                       
      CALL init_devices 
!                                                                       
      DO iw = 1, maxwin 
      dev_width (iw) = 10.0 
      dev_height (iw) = 7.5 
!                                                                       
      dev_sf (iw, x11) = 0.7 
      dev_id (iw, x11) = - 1 
      dev_sf (iw, ps) = 1.0 
      dev_id (iw, ps) = - 1 
      dev_sf (iw, vps) = 1.0 
      dev_id (iw, vps) = - 1 
      dev_sf (iw, pic) = 0.7 
      dev_id (iw, pic) = - 1 
      dev_sf (iw, vpic) = 0.7 
      dev_id (iw, vpic) = - 1 
      dev_sf (iw, png ) = 1.0 
      dev_id (iw, png ) = - 1 
      dev_sf (iw, lat) = 1.0 
      dev_id (iw, lat) = - 1 
      ENDDO 
!                                                                       
!------ common /errors/                                                 
!                                                                       
      ier_num = 0 
      ier_typ = ER_NONE 
      ier_sta = ER_S_CONT 
!                                                                       
      DO i = 1, 3 
      ier_msg (i) = ' ' 
      ENDDO 
!                                                                       
!------ setting default color map for bitmaps                           
!                                                                       
      DO iwin = 1, maxwin 
      CALL cmap_gray (.false.) 
      ENDDO 
!                                                                       
!------ setting default line colours                                    
!------ Background                                                      
!                                                                       
      DO iw = 1, maxwin 
      colour (iw, 0, 1) = 1.0 
      colour (iw, 0, 2) = 1.0 
      colour (iw, 0, 3) = 1.0 
!                                                                       
!------ - Colours set to PGPLOT defaults                                
!                                                                       
      colour (iw, 1, 1) = 1.0 
      colour (iw, 1, 2) = 0.0 
      colour (iw, 1, 3) = 0.0 
!                                                                       
      colour (iw, 2, 1) = 0.0 
      colour (iw, 2, 2) = 1.0 
      colour (iw, 2, 3) = 0.0 
!                                                                       
      colour (iw, 3, 1) = 0.0 
      colour (iw, 3, 2) = 0.0 
      colour (iw, 3, 3) = 1.0 
!                                                                       
      colour (iw, 4, 1) = 1.0 
      colour (iw, 4, 2) = 0.0 
      colour (iw, 4, 3) = 1.0 
!                                                                       
      colour (iw, 5, 1) = 1.0 
      colour (iw, 5, 2) = 1.0 
      colour (iw, 5, 3) = 0.0 
!                                                                       
      colour (iw, 6, 1) = 0.0 
      colour (iw, 6, 2) = 0.0 
      colour (iw, 6, 3) = 0.0 
!                                                                       
      colour (iw, 7, 1) = 0.6 
      colour (iw, 7, 2) = 0.0 
      colour (iw, 7, 3) = 0.0 
!                                                                       
      colour (iw, 8, 1) = 0.0 
      colour (iw, 8, 2) = 0.6 
      colour (iw, 8, 3) = 0.0 
!                                                                       
      colour (iw, 9, 1) = 0.0 
      colour (iw, 9, 2) = 0.0 
      colour (iw, 9, 3) = 0.6 
!                                                                       
      colour (iw, 10, 1) = 0.6 
      colour (iw, 10, 2) = 0.0 
      colour (iw, 10, 3) = 0.6 
!                                                                       
      colour (iw, 11, 1) = 0.6 
      colour (iw, 11, 2) = 0.6 
      colour (iw, 11, 3) = 0.0 
!                                                                       
      colour (iw, 12, 1) = 0.5 
      colour (iw, 12, 2) = 0.5 
      colour (iw, 12, 3) = 0.5 
!                                                                       
      colour (iw, 13, 1) = 0.0 
      colour (iw, 13, 2) = 1.0 
      colour (iw, 13, 3) = 1.0 
!                                                                       
      colour (iw, 14, 1) = 0.0 
      colour (iw, 14, 2) = 0.6 
      colour (iw, 14, 3) = 0.6 
!                                                                       
      colour (iw, 15, 1) = 1.0 
      colour (iw, 15, 2) = 1.0 
      colour (iw, 15, 3) = 1.0 
!                                                                       
      ENDDO 
!                                                                       
!------ here start all the defaults for the variables                   
!                                                                       
      ftyp = 'NONE' 
      fstart = .true. 
      wtyp = 'INV' 
      frall = .false. 
      urf = 0.1 
      ncycle = 30 
      npara = 2 
      np1 = 2 
      np2 = 0 
      np3 = 0 
      fit_ifen = 5 
      exclude9999 = 0 
      ifen = 9 
!                                                                       
      offxy (0) = 0 
      offz (0) = 0 
!                                                                       
      DO i = 1, maxpara 
      p (i) = 0.0 
      pinc (i) = 1.0 
      ENDDO 
      p_origin = 0 
!                                                                       
      DO i = 1, maxkurvtot 
      lni (i) = .false. 
      ikfirst (i) = .true. 
      ENDDO 
!                                                                       
      iz = 1 
      iframe = 1 
      iwin = 1 
!                                                                       
      DO iw = 1, maxwin 
      ipara (iw) = 2 
      fpara (iw) = 'kupl.par' 
      orient (iw) = .true. 
      tot_frame (iw) = .false. 
      iaf (iw) = 1 
      icol_old (iw) = - 1 
!                                                                       
      lgui (iw) = .false. 
      iden (iw) = .false. 
!                                                                       
      DO i = 1, maxframe 
      ax_d (iw, i, 1) = 1.5 
      ax_d (iw, i, 2) = 1.5 
      lab_angle (iw, i, 1) = 0.0 
      lab_angle (iw, i, 2) = 0.0 
      lab_d (iw, i, 1) = 0.3 
      lab_d (iw, i, 2) = 0.3 
      tick_ma_h (iw, i, 1) = 0.5 
      tick_ma_h (iw, i, 2) = 0.5 
      tick_mi_h (iw, i, 1) = 0.5 
      tick_mi_h (iw, i, 2) = 0.5 
      tick_nsub (iw, i, 1) = 2 
      tick_nsub (iw, i, 2) = 2 
!                                                                       
      ibuf (iw, i, 1) = 0.1 
      ibuf (iw, i, 2) = 0.1 
      ibuf (iw, i, 3) = 0.1 
      ibuf (iw, i, 4) = 0.1 
      linewid (iw, i, 0) = 0.3 
      ilinecol (iw, i, 0) = 6 
      iho (iw, i) = 1 
!                                                                       
      info_orig (iw, i, 1) = - 9999. 
      info_orig (iw, i, 2) = - 9999. 
!                                                                       
      DO k = 1, maxbond 
      bond_rad (iw, i, k) = 0.0 
      bond_sig (iw, i, k) = 0.2 
      bond_lwid (iw, i, k) = 0.5 
      bond_lcol (iw, i, k) = stift (mod (k, 15) ) 
      bond_ltyp (iw, i, k) = 1 
      ENDDO 
!                                                                       
      DO j = 1, maxkurvtot 
      DO k = 1, maxhl 
      hlinecol (iw, i, j, k) = stift (mod (j, 15) ) 
      hlinetyp (iw, i, j, k) = 1 
      ENDDO 
      DO k = 1, 4 
      fillrange (iw, i, j, k) = - 9999.0 
      ENDDO 
      hlabel (iw, i, j) = 0 
      linewid (iw, i, j) = 0.2 
      ilegend (iw, i, j) = 0 
      ilinecol (iw, i, j) = stift (mod (j, 15) ) 
      ifillcol (iw, i, j) = stift (mod (j, 15) ) 
      imarkcol (iw, i, j) = stift (mod (j, 15) ) 
      ierrcol (iw, i, j) = stift (mod (j, 15) ) 
      ierr (iw, i, j) = 0 
      ilinetyp (iw, i, j) = 1 
      ifilltyp (iw, i, j) = 0 
      imarktyp (iw, i, j) = 0 
      imarkmax (iw, i, j) = 0 
      ilineart (iw, i, j) = 1 
      hlineart (iw, i, j) = 1 
      sizemark (iw, i, j) = 0.2 
      rel_mark (iw, i, j) = 0
      infra (iw, i, j) = j 
      ENDDO 
                                                                        
      DO j = 1, maxan 
      antext (iw, i, j) = 'OFF' 
      anx (iw, i, j) = 0.0 
      any (iw, i, j) = 0.0 
      antx (iw, i, j) = 0.0 
      anty (iw, i, j) = 0.0 
      anjust (iw, i, j) = if_left 
      anangle (iw, i, j) = 0.0 
      ENDDO 
!                                                                       
      DO j = 1, maxhl 
      z_min (iw, i, j) = 50.0 
      z_inc (iw, i, j) = 50.0 
      nz (iw, i, j) = 10 
      ENDDO 
!                                                                       
      frjust (iw, i) = if_left 
      frback (iw, i, 1) = colour (iw, 0, 1) 
      frback (iw, i, 2) = colour (iw, 0, 2) 
      frback (iw, i, 3) = colour (iw, 0, 3) 
!                                                                       
      fonscal (iw, i) = 1.0 
!                                                                       
      fonsize (iw, i, 1) = 160.0 
      fonsize (iw, i, 2) = 140.0 
      fonsize (iw, i, 3) = 120.0 
      fonsize (iw, i, 4) = 120.0 
      fonsize (iw, i, 5) = 120.0 
      fonsize (iw, i, 6) = 120.0 
!                                                                       
      fon_id (iw, i, 1) = 2 
      fon_id (iw, i, 2) = 2 
      fon_id (iw, i, 3) = 2 
      fon_id (iw, i, 4) = 2 
      fon_id (iw, i, 5) = 1 
      fon_id (iw, i, 6) = 2 
!                                                                       
      foncol (iw, i, 1) = 6 
      foncol (iw, i, 2) = 6 
      foncol (iw, i, 3) = 6 
      foncol (iw, i, 4) = 6 
      foncol (iw, i, 5) = 6 
      foncol (iw, i, 6) = 6 
!                                                                       
      ifname (iw, i) = .true. 
      igrid (iw, i) = .false. 
      ibox (iw, i) = 3 
      ftext (iw, i) = 'kupl.par' 
      ex (iw, i, 1) = - 9999. 
      ex (iw, i, 2) = - 9999. 
      ey (iw, i, 1) = - 9999. 
      ey (iw, i, 2) = - 9999. 
      t (iw, i, 1) = - 9999. 
      t (iw, i, 2) = - 9999. 
      yskal_u (iw, i) = - 9999. 
      lyskal (iw, i) = .false. 
!                                                                       
      titel (iw, i, 1) = ' ' 
      titel (iw, i, 2) = ' ' 
      achse (iw, i, 1) = ' x ' 
      achse (iw, i, 2) = ' y ' 
      achse (iw, i, 3) = 'OFF' 
      lachse (iw, i, 1) = .false. 
      lachse (iw, i, 2) = .false. 
      lachse (iw, i, 3) = .false. 
      shear (iw, i) = 90.0 
      sfl (iw, i) = .true. 
      frame (iw, i, 1) = 0.0 
      frame (iw, i, 2) = 0.0 
      frame (iw, i, 3) = 1.0 
      frame (iw, i, 4) = 1.0 
      ENDDO 
      ENDDO 
!                                                                       
      END SUBROUTINE kuplot_initarrays                     
!*****7*****************************************************************
      SUBROUTINE do_read_def (zeile, lp) 
!+                                                                      
!     Reading defaults                                                  
!-                                                                      
      USE build_name_mod
      USE errlist_mod 
      USE get_params_mod
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER maxw, idef 
      PARAMETER (maxw = 5) 
      PARAMETER (idef = 9) 
!                                                                       
      CHARACTER ( * ) zeile 
      CHARACTER(1024) cpara (maxw) 
      CHARACTER(80) cdummy 
      INTEGER lp, ianz 
      INTEGER lpara (maxw) 
      REAL werte (maxw) 
!                                                                       
      INTEGER len_str 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
      IF (ianz.eq.0) then 
         cdummy = 'kuplot.def' 
      ELSE 
         CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
         IF (ier_num.ne.0) return 
         cdummy = cpara (1)(1:MIN(80,LEN_TRIM(cpara(1))))
      ENDIF 
!                                                                       
      CALL oeffne (idef, cdummy, 'old') 
      IF (ier_num.eq.0) then 
         CALL read_def (idef) 
         WRITE (output_io, 1000) cdummy (1:len_str (cdummy) ) 
         CLOSE (idef) 
      ENDIF 
!                                                                       
 1000 FORMAT     (' ------ > Reading defaults from file ',a,' ..') 
      END SUBROUTINE do_read_def                    
!*****7*****************************************************************
      SUBROUTINE do_write_def (zeile, lp) 
!+                                                                      
!     Writing defaults                                                  
!-                                                                      
      USE build_name_mod
      USE errlist_mod 
      USE get_params_mod
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER maxw, idef 
      PARAMETER (maxw = 5) 
      PARAMETER (idef = 9) 
!                                                                       
      CHARACTER ( * ) zeile 
      CHARACTER(1024) cpara (maxw) 
      CHARACTER(80) cdummy 
      INTEGER lpara (maxw) 
      INTEGER lp, ianz 
      REAL werte (maxw) 
!                                                                       
      INTEGER len_str 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
      IF (ianz.eq.0) then 
         cdummy = 'kuplot.def' 
      ELSE 
         CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
         IF (ier_num.ne.0) return 
         cdummy = cpara (1)(1:MIN(80,LEN_TRIM(cpara(1))))
      ENDIF 
!                                                                       
      CALL oeffne (idef, cdummy, 'unknown') 
      IF (ier_num.eq.0) then 
         CALL write_def (idef) 
         WRITE (output_io, 1000) cdummy (1:len_str (cdummy) ) 
         CLOSE (idef) 
      ENDIF 
!                                                                       
 1000 FORMAT     (' ------ > Writing defaults to file ',a,' ..') 
      END SUBROUTINE do_write_def                   
!*****7*****************************************************************
      SUBROUTINE write_def (idef) 
!+                                                                      
!       This routine saves the current settings to a file               
!-                                                                      
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER(1) q 
      INTEGER idef, i, j, k 
!                                                                       
      INTEGER len_str 
!                                                                       
      q = char (39) 
!                                                                       
      WRITE (idef, 1000) 
      WRITE (idef, 9998) 'MAXKURVTOT,MAXFRAME,MAXBOND,MAXHL,MAXWIN' 
      WRITE (idef, 1100) MAXKURVTOT, MAXFRAME, MAXBOND, MAXHL, MAXWIN 
      WRITE (idef, 9998) 'Title line 1 and 2 for each frame' 
      WRITE (idef, 9997) ( ( (q, titel (k, i, j) (1:len_str (titel (k,  &
      i, j) ) ), q, k = 1, MAXWIN), i = 1, MAXFRAME), j = 1, 2)         
      WRITE (idef, 9998) 'Axes lables x,y and z for each frame' 
      WRITE (idef, 9997) ( ( (q, achse (k, i, j) (1:len_str (achse (k,  &
      i, j) ) ), q, k = 1, MAXWIN), i = 1, MAXFRAME), j = 1, 3)         
      WRITE (idef, 9998) 'Buffer space for each frame (l,r,b,t)' 
      WRITE (idef, 1600) ibuf 
      WRITE (idef, 9998) 'Code for box around plot (-> buf)' 
      WRITE (idef, 1800) ibox 
      WRITE (idef, 9998) 'Code (T/F) for grid plotting (-> grid)' 
      WRITE (idef, 1900) igrid 
      WRITE (idef, 9998) 'Code (T/F) for paper orientation (-> orient)' 
      WRITE (idef, 1900) orient 
      WRITE (idef, 9998) 'Marker size for each data set and each frame' 
      WRITE (idef, 1600) sizemark 
      WRITE (idef, 9998) 'Line width ...' 
      WRITE (idef, 1600) linewid 
      WRITE (idef, 9998) 'Line colour ...' 
      WRITE (idef, 1800) ilinecol 
      WRITE (idef, 9998) 'Marker colour ...' 
      WRITE (idef, 1800) imarkcol 
      WRITE (idef, 9998) 'Line type ...' 
      WRITE (idef, 1800) ilinetyp 
      WRITE (idef, 9998) 'Marker type ...' 
      WRITE (idef, 1800) imarktyp 
      WRITE (idef, 9998) 'Line style ...' 
      WRITE (idef, 1800) ilineart 
      WRITE (idef, 9998) 'Maximum marker type ... ' 
      WRITE (idef, 1800) imarkmax 
      WRITE (idef, 9998) 'Error bar type ... ' 
      WRITE (idef, 1800) ierr 
      WRITE (idef, 9998) 'Frame corners for each frame' 
      WRITE (idef, 1500) frame 
      WRITE (idef, 9998) 'Frame background colour for each frame' 
      WRITE (idef, 1600) frback 
      WRITE (idef, 9998) 'Font sizes ...' 
      WRITE (idef, 1700) fonsize 
      WRITE (idef, 9998) 'Font overall scaling ...' 
      WRITE (idef, 1700) fonscal 
      WRITE (idef, 9998) 'Font ids ...' 
      WRITE (idef, 1800) fon_id 
      WRITE (idef, 9998) 'Font colour ...' 
      WRITE (idef, 1800) foncol 
      WRITE (idef, 9998) 'Code text justification ...' 
      WRITE (idef, 1800) frjust 
      WRITE (idef, 9998) 'Code filename plotting ...' 
      WRITE (idef, 1900) ifname 
      WRITE (idef, 9998) 'Code border around frames ...' 
      WRITE (idef, 1900) tot_frame 
      WRITE (idef, 9998) 'Number of contour sets' 
      WRITE (idef, 1800) iho 
      WRITE (idef, 9998) 'Lowest contour per set per frame' 
      WRITE (idef, 1400) z_min 
      WRITE (idef, 9998) 'Contour increment ...' 
      WRITE (idef, 1400) z_inc 
      WRITE (idef, 9998) 'Contour line colour ...' 
      WRITE (idef, 1800) hlinecol 
      WRITE (idef, 9998) 'Contour line type ...' 
      WRITE (idef, 1800) hlinetyp 
      WRITE (idef, 9998) 'Contour line style ...' 
      WRITE (idef, 1800) hlineart 
      WRITE (idef, 9998) 'Contour line labelling ...' 
      WRITE (idef, 1800) hlabel 
      WRITE (idef, 9998) 'Process -9999. ...' 
      WRITE (idef, 1800) exclude9999 
      WRITE (idef, 9998) 'Bond length ...' 
      WRITE (idef, 1400) bond_rad 
      WRITE (idef, 9998) 'Bond length sigma ...' 
      WRITE (idef, 1400) bond_sig 
      WRITE (idef, 9998) 'Bond line width ...' 
      WRITE (idef, 1600) bond_lwid 
      WRITE (idef, 9998) 'Bond line colour ...' 
      WRITE (idef, 1800) bond_lcol 
      WRITE (idef, 9998) 'Bond line type ...' 
      WRITE (idef, 1800) bond_ltyp 
      WRITE (idef, 9998) 'Line colours ...' 
      WRITE (idef, 1600) colour 
      WRITE (idef, 9998) 'Scale factors of graphics devices ...' 
      WRITE (idef, 1600) dev_sf 
      WRITE (idef, 9998) 'Colours for error bars' 
      WRITE (idef, 1800) ierrcol 
      WRITE (idef, 9998) 'Angle for text annotations' 
      WRITE (idef, 1800) anangle 
      WRITE (idef, 9998) 'Colour for filling' 
      WRITE (idef, 1800) ifillcol 
      WRITE (idef, 9998) 'Style for filling' 
      WRITE (idef, 1800) ifilltyp 
      WRITE (idef, 9998) 'Range for filling' 
      WRITE (idef, 1400) fillrange 
      WRITE (idef, 9998) 'Angle of axes labels' 
      WRITE (idef, 1600) lab_angle 
      WRITE (idef, 9998) 'Distance numbers from axes' 
      WRITE (idef, 1600) lab_d 
      WRITE (idef, 9998) 'Distance labels from axes' 
      WRITE (idef, 1600) ax_d 
      WRITE (idef, 9998) 'Height major tick marks' 
      WRITE (idef, 1500) tick_ma_h 
      WRITE (idef, 9998) 'Height minor tick marks' 
      WRITE (idef, 1500) tick_mi_h 
      WRITE (idef, 9998) 'Number of ticks between major labels' 
      WRITE (idef, 1800) tick_nsub 
      WRITE (idef, 9998) 'Plot date/time/username under plot' 
      WRITE (idef, 1900) iden 
      WRITE (idef, 9998) 'Log. axes flag' 
      WRITE (idef, 1900) lachse 
!                                                                       
 1000 FORMAT    (32('#'),' KUPLOT Defaults ',31('#')) 
 1100 FORMAT    (5(i8,1x)) 
 1400 FORMAT    (6(g11.5,1x)) 
 1500 FORMAT    (10(f7.4,1x)) 
 1600 FORMAT    (10(f7.2,1x)) 
 1700 FORMAT    (10(f7.0,1x)) 
 1800 FORMAT    (20(i3,1x)) 
 1900 FORMAT    (40(l1,1x)) 
 9997 FORMAT    (10(a1,a,a1,2x)) 
 9998 FORMAT    ('# ',a) 
      END SUBROUTINE write_def                      
!*****7*****************************************************************
      SUBROUTINE read_def (idef) 
!+                                                                      
!       This routine reads defaults from a file                         
!-                                                                      
      USE errlist_mod 
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER(80) cdummy 
      INTEGER idef, i1, i2, i3, i4, i5 
!                                                                       
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) i1, i2, i3, i4, i5 
!                                                                       
!------ Check is file has correct array sizes                           
!                                                                       
      IF (i1.ne.MAXKURVTOT.or.i2.ne.MAXFRAME.or.i3.ne.MAXBOND.or.i4.ne.M&
     &AXHL.or.i5.ne.MAXWIN) then                                        
         ier_num = - 36 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
!                                                                       
!------ All seems ok, so we read the defaults                           
!                                                                       
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) titel 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) achse 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) ibuf 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) ibox 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) igrid 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) orient 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) sizemark 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) linewid 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) ilinecol 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) imarkcol 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) ilinetyp 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) imarktyp 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) ilineart 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) imarkmax 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) ierr 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) frame 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) frback 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) fonsize 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) fonscal 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) fon_id 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) foncol 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) frjust 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) ifname 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) tot_frame 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) iho 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) z_min 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) z_inc 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) hlinecol 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) hlinetyp 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) hlineart 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) hlabel 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) exclude9999 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) bond_rad 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) bond_sig 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) bond_lwid 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) bond_lcol 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) bond_ltyp 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) colour 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) dev_sf 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) ierrcol 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) anangle 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) ifillcol 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) ifilltyp 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) fillrange 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) lab_angle 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) lab_d 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) ax_d 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) tick_ma_h 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) tick_mi_h 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) tick_nsub 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) iden 
      READ (idef, 1, err = 98, end = 99) cdummy 
      READ (idef, *, err = 98, end = 99) lachse 
!                                                                       
      RETURN 
!                                                                       
   98 CONTINUE 
      ier_num = - 3 
      ier_typ = ER_IO 
      RETURN 
!                                                                       
   99 CONTINUE 
      ier_num = - 6 
      ier_typ = ER_IO 
      RETURN 
!                                                                       
    1 FORMAT  (a) 
      END SUBROUTINE read_def                       
!*****7*****************************************************************
      SUBROUTINE kuplot_auto_def 
!+                                                                      
!     If the default setup file is there load it ..                     
!-                                                                      
      USE envir_mod 
      USE errlist_mod 
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
      INTEGER idef 
!                                                                       
!                                                                       
      idef = 33 
      ier_num = 0 
      ier_typ = ER_NONE 
!                                                                       
      CALL open_def (idef) 
      IF (ier_num.eq.0) then 
         WRITE (output_io, 1000) deffile 
         CALL read_def (idef) 
      ELSE 
         WRITE (output_io, * ) 
         ier_num = 0 
         ier_typ = ER_NONE 
      ENDIF 
      CLOSE (idef) 
!                                                                       
 1000 FORMAT   (1x,'Reading defaults : ',a) 
      END SUBROUTINE kuplot_auto_def                       
