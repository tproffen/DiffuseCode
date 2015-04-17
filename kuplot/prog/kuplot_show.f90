!*****7*****************************************************************
!     This routine shows all settings of KUPLOT parameters. Some        
!     parameter settings can also be displayed by entering the          
!     corresponding command without parameters.                         
!*****7*****************************************************************
      SUBROUTINE kuplot_do_show (line, lp) 
!                                                                       
!     Main show menu                                                    
!                                                                       
      USE errlist_mod 
      USE prompt_mod 
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 2) 
!                                                                       
      CHARACTER ( * ) line 
      CHARACTER(1024) cpara (maxw) 
      REAL werte (maxw) 
      INTEGER lpara (maxw), lp 
      INTEGER ianz, iku 
      LOGICAL str_comp 
!                                                                       
      CALL get_params (line, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
!                                                                       
      IF (ianz.eq.0) then 
!                                                                       
!------ - Output all information                                        
!                                                                       
         CALL show_data ( - 1) 
         CALL show_plot 
         CALL show_style 
         CALL show_hlin 
         CALL show_annotation 
         CALL show_frames 
         CALL show_font 
!                                                                       
      ELSEIF (ianz.eq.1.or.ianz.eq.2) then 
!                                                                       
!------ - Contour line information                                      
!                                                                       
         IF (str_comp (cpara (1) , 'hlin', 2, lpara (1) , 4) ) then 
            CALL show_hlin 
!                                                                       
!------ - Data information                                              
!                                                                       
         ELSEIF (str_comp (cpara (1) , 'data', 2, lpara (1) , 4) ) then 
            IF (ianz.eq.1) then 
               CALL show_data ( - 1) 
            ELSE 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               iku = nint (werte (1) ) 
               CALL show_data (iku) 
            ENDIF 
!                                                                       
!------ - Plot information                                              
!                                                                       
         ELSEIF (str_comp (cpara (1) , 'plot', 2, lpara (1) , 4) ) then 
            CALL show_plot 
!                                                                       
!------ - Style information (colors, ...)                               
!                                                                       
         ELSEIF (str_comp (cpara (1) , 'style', 2, lpara (1) , 5) )     &
         then                                                           
            CALL show_style 
!                                                                       
!------ - Show annotations                                              
!                                                                       
         ELSEIF (str_comp (cpara (1) , 'annot', 2, lpara (1) , 5) )     &
         then                                                           
            CALL show_annotation 
!                                                                       
!------ - Current color map                                             
!                                                                       
         ELSEIF (str_comp (cpara (1) , 'cmap', 3, lpara (1) , 4) ) then 
            CALL show_cmap 
!                                                                       
!------ - Current color settings                                        
!                                                                       
         ELSEIF (str_comp (cpara (1) , 'color', 3, lpara (1) , 5) )     &
         then                                                           
            CALL show_color 
!                                                                       
!------ - Current frame setting                                         
!                                                                       
         ELSEIF (str_comp (cpara (1) , 'frame', 2, lpara (1) , 5) )     &
         then                                                           
            CALL show_frames 
!                                                                       
!------ - Configuration of KUPLOT                                       
!                                                                       
         ELSEIF (str_comp (cpara (1) , 'config', 3, lpara (1) , 5) )    &
         then                                                           
            CALL show_config 
!                                                                       
!------ - Current font settings                                         
!                                                                       
         ELSEIF (str_comp (cpara (1) , 'font', 2, lpara (1) , 4) ) then 
            CALL show_font 
!                                                                       
         ELSE 
            CALL do_show_generic (cpara, lpara, maxw) 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
      END SUBROUTINE kuplot_do_show                        
!*****7*****************************************************************
      SUBROUTINE show_config 
!+                                                                      
!     Show KUPLOT configuration                                         
!-                                                                      
      USE prompt_mod 
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      WRITE (output_io, 1000) 
      WRITE (output_io, 2000) MAXARRAY 
      WRITE (output_io, 2010) MAXKURVTOT 
      WRITE (output_io, 2020) MAXZ 
      WRITE (output_io, 2030) MAXHL 
      WRITE (output_io, 2040) MAXFRAME 
      WRITE (output_io, 2045) MAXWIN 
      WRITE (output_io, 2050) MAXSP 
      WRITE (output_io, 2060) MAXAN 
      WRITE (output_io, 2070) MAXBOND 
      WRITE (output_io, 2080) MAXPARA 
!                                                                       
 1000 FORMAT     (' Current KUPLOT configuration :',/) 
 2000 FORMAT     ('   Maximum total number of data points      : ',i8) 
 2010 FORMAT     ('   Maximum number of data sets              : ',i8) 
 2020 FORMAT     ('   Maximum bitmap size (maxz,maxz)          : ',i8) 
 2030 FORMAT     ('   Maximum number of contour line sets      : ',i8) 
 2040 FORMAT     ('   Maximum number of frames                 : ',i8) 
 2045 FORMAT     ('   Maximum number of windows                : ',i8) 
 2050 FORMAT     ('   Maximum number of points for step/spline : ',i8) 
 2060 FORMAT     ('   Maximum number of annotations            : ',i8) 
 2070 FORMAT     ('   Maximum number of different bonds        : ',i8) 
 2080 FORMAT     ('   Maximum number of fit parameters         : ',i8,/) 
      END SUBROUTINE show_config                    
!*****7*****************************************************************
      SUBROUTINE show_hlin 
!+                                                                      
!     Show contour line information                                     
!-                                                                      
      USE errlist_mod 
      USE prompt_mod 
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER j 
!                                                                       
      WRITE (output_io, 1000) iwin, iframe 
      WRITE (output_io, 1010) 
      DO j = 1, iho (iwin, iframe) 
      WRITE (output_io, 1020) j, z_min (iwin, iframe, j), z_inc (iwin,  &
      iframe, j), nz (iwin, iframe, j)                                  
      ENDDO 
      WRITE (output_io, * ) 
!                                                                       
 1000 FORMAT     (' Contour line settings for window ',i3,              &
     &                   ' / frame ',i3,/)                              
 1010 FORMAT     (3x,'Cont.set    minimum       ',                      &
     &                   'increment   # of contours',/,3x,56('-'))      
 1020 FORMAT     (4x,i3,7x,g11.5,3x,g11.5,6x,i4) 
      END SUBROUTINE show_hlin                      
!*****7*****************************************************************
      SUBROUTINE show_data (iku) 
!+                                                                      
!     Print information about loaded data                               
!-                                                                      
      USE errlist_mod 
      USE prompt_mod 
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER iku, ia, ie, ik 
!                                                                       
      INTEGER len_str 
!                                                                       
!------ check if we want ALL data sets                                  
!                                                                       
      IF (iku.eq. - 1) then 
         ia = 1 
         ie = iz - 1 
      ELSE 
         IF (iku.gt.iz - 1.or.iku.lt.1) then 
            ier_num = - 4 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
         ia = iku 
         ie = iku 
      ENDIF 
!                                                                       
!------ Print information                                               
!                                                                       
      DO ik = ia, ie 
      IF (lni (ik) ) then 
         CALL get_extrema_xy (x, ik, nx (ik), xmin, xmax) 
         CALL get_extrema_xy (y, ik, ny (ik), ymin, ymax) 
         CALL get_extrema_z (z, ik, nx (ik), ny (ik), zmin, zmax) 
         WRITE (output_io, 1000) ik, fname (ik) (1:len_str (fname (ik) )&
         ), nx (ik), ny (ik), xmin (ik), xmax (ik), ymin (ik), ymax (ik)&
         , zmin (ik), zmax (ik)                                         
      ELSE 
         CALL get_extrema_xy (x, ik, len (ik), xmin, xmax) 
         CALL get_extrema_xy (y, ik, len (ik), ymin, ymax) 
         WRITE (output_io, 1010) ik, fname (ik) (1:len_str (fname (ik) )&
         ), len (ik), xmin (ik), xmax (ik), ymin (ik), ymax (ik)        
      ENDIF 
      ENDDO 
!                                                                       
 1000 FORMAT     (' Data set no.: ',i3,//,                              &
     &            3x,'Filename : ',a,/                                  &
     &                     3x,'Size     : ',i6,' x ',i6,' points',//    &
     &            3x,'Range  x : ',g12.4,' to ',g12.4/                  &
     &            3x,'Range  y : ',g12.4,' to ',g12.4/                  &
     &            3x,'Range  z : ',g12.4,' to ',g12.4,/)                
 1010 FORMAT     (' Data set no.: ',i3,//,                              &
     &            3x,'Filename : ',a,/                                  &
     &                     3x,'Size     : ',i6,' points',//             &
     &            3x,'Range  x : ',g12.4,' to ',g12.4/                  &
     &            3x,'Range  y : ',g12.4,' to ',g12.4,/)                
!                                                                       
      END SUBROUTINE show_data                      
!*****7**************************************************************** 
      SUBROUTINE show_plot 
!+                                                                      
!     Show information about plot settings                              
!-                                                                      
      USE errlist_mod 
      USE prompt_mod 
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      WRITE (output_io, 1000) iwin, iframe, yskal_u (iwin, iframe),     &
      shear (iwin, iframe), ilinecol (iwin, iframe, 0), linewid (iwin,  &
      iframe, 0), ibuf (iwin, iframe, 1), ibuf (iwin, iframe, 2),       &
      ibuf (iwin, iframe, 3), ibuf (iwin, iframe, 4)                    
!                                                                       
      IF (ibox (iwin, iframe) .eq.0) then 
         WRITE (output_io, 1005) 'no frame' 
      ELSEIF (abs (ibox (iwin, iframe) ) .eq.1) then 
         WRITE (output_io, 1005) 'only box around plot' 
      ELSEIF (abs (ibox (iwin, iframe) ) .eq.2) then 
         WRITE (output_io, 1005) 'box and tick marks / labels' 
      ELSEIF (abs (ibox (iwin, iframe) ) .eq.3) then 
         WRITE (output_io, 1005) 'box, tick marks and x,y-axis' 
      ENDIF 
!                                                                       
      IF (ibox (iwin, iframe) .ge.0) then 
         WRITE (output_io, 1007) 'LEFT' 
      ELSE 
         WRITE (output_io, 1007) 'RIGHT' 
      ENDIF 
!                                                                       
      IF (ifname (iwin, iframe) ) then 
         WRITE (output_io, 1010) 'on ' 
      ELSE 
         WRITE (output_io, 1010) 'off' 
      ENDIF 
      IF (tot_frame (iwin) ) then 
         WRITE (output_io, 1013) 'on ' 
      ELSE 
         WRITE (output_io, 1013) 'off' 
      ENDIF 
      IF (igrid (iwin, iframe) ) then 
         WRITE (output_io, 1015) 'on ' 
      ELSE 
         WRITE (output_io, 1015) 'off' 
      ENDIF 
      IF (iden (iwin) ) then 
         WRITE (output_io, 1025) 'on ' 
      ELSE 
         WRITE (output_io, 1025) 'off' 
      ENDIF 
!                                                                       
 1000 FORMAT    (' Parameter for window ',i3,' / frame ',i3,//          &
     &                   3x,'Ratio y/x-axis           : ',g11.5/        &
     &                   3x,'Angle between x&y axis   : ',f7.2/         &
     &                   3x,'Color of grid lines      : ',i7/           &
     &                   3x,'Line width of frame [cm] : ',f7.2/         &
     &                   3x,'Buffer around plot       : ',4(f4.2,1x))   
 1005 FORMAT    ( 3x,'Type of frame            : ',a) 
 1007 FORMAT    ( 3x,'Location of y-labels     : ',a) 
 1010 FORMAT    ( 3x,'Plot filenames           : ',a) 
 1013 FORMAT    ( 3x,'Plot frame               : ',a) 
 1015 FORMAT    ( 3x,'Plot grid                : ',a) 
 1025 FORMAT    ( 3x,'Plot user/date info.     : ',a/) 
      END SUBROUTINE show_plot                      
!*****7**************************************************************** 
      SUBROUTINE show_frames 
!+                                                                      
!     Show current frame settings                                       
!-                                                                      
      USE prompt_mod 
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER i, j 
!                                                                       
      WRITE (output_io, 1000) iwin, iframe 
      WRITE (output_io, 2000) 
      DO i = 1, iaf (iwin) 
      WRITE (output_io, 2500) i, (frame (iwin, i, j), j = 1, 4),        &
      (frback (iwin, i, j), j = 1, 3)                                   
      ENDDO 
!                                                                       
      WRITE (output_io, 3000) 
      DO i = 1, iaf (iwin) 
      IF (infra (iwin, i, 1) .ne. - 1) then 
         WRITE (output_io, 3500) i, (infra (iwin, i, j), j = 1,         &
         maxkurvtot)                                                    
      ELSE 
         WRITE (output_io, 3600) i, ftext (iwin, i) 
      ENDIF 
      ENDDO 
      WRITE (output_io, * ) 
!                                                                       
 1000 FORMAT (' Frame settings (window ',i3,' / frame ',i3,') :',/) 
 2000 FORMAT (3x,'Frame  xmin ymin xmax ymax  Backgr.  r   g   b',/     &
     &        3x,'------------------------------------------------')    
 2500 FORMAT (3x,i3,4x,4(f4.2,1x),9x,3(f3.1,1x)) 
 3000 FORMAT (/,3x,'Frame  data sets/text file',/,3x,67('-')) 
 3500 FORMAT (3x,i3,4x,20i3) 
 3600 FORMAT (3x,i3,4x,a40) 
      END SUBROUTINE show_frames                    
!*****7**************************************************************** 
      SUBROUTINE show_font 
!+                                                                      
!     Show current font settings                                        
!-                                                                      
      USE prompt_mod 
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER(20) fon_name (4) 
!                                                                       
      DATA fon_name / 'Standard', 'Roman', 'Italic', 'Script' / 
!                                                                       
      WRITE (output_io, 1000) iwin, iframe, fonscal (iwin, iframe) 
      WRITE (output_io, 1100) 1, 'Main title line    ', fonsize (iwin, i&
     &frame, 1)  / 10., fon_name (fon_id (iwin, iframe, 1) ) , fon_id (i&
     &win, iframe, 1) , foncol (iwin, iframe, 1)                        
      WRITE (output_io, 1100) 2, 'Subtitle line      ', fonsize (iwin, i&
     &frame, 2)  / 10., fon_name (fon_id (iwin, iframe, 2) ) , fon_id (i&
     &win, iframe, 2) , foncol (iwin, iframe, 2)                        
      WRITE (output_io, 1100) 3, 'Axis labels        ', fonsize (iwin, i&
     &frame, 3)  / 10., fon_name (fon_id (iwin, iframe, 3) ) , fon_id (i&
     &win, iframe, 3) , foncol (iwin, iframe, 3)                        
      WRITE (output_io, 1100) 4, 'Numbers at axis    ', fonsize (iwin, i&
     &frame, 4)  / 10., fon_name (fon_id (iwin, iframe, 4) ) , fon_id (i&
     &win, iframe, 4) , foncol (iwin, iframe, 4)                        
      WRITE (output_io, 1100) 5, 'Text in text frame ', fonsize (iwin,  &
      iframe, 5) / 10., fon_name (fon_id (iwin, iframe, 5) ) , fon_id ( &
      iwin, iframe, 5) , foncol (iwin, iframe, 5)                       
      WRITE (output_io, 1100) 6, 'Filename & caption ', fonsize (iwin,  &
      iframe, 6) / 10., fon_name (fon_id (iwin, iframe, 6) ) , fon_id ( &
      iwin, iframe, 6) , foncol (iwin, iframe, 6)                       
      WRITE (output_io, * ) 
!                                                                       
 1000 FORMAT (' Font settings for window ',i3,' / frame ',i3,/,         &
     &      3x,'Overall font scaling factor  : ',f6.2,//,               &
     &      3x,'Font  where ?              size   fontname',            &
     &      12x,'font-id   color',/,3x,69('-'))                         
 1100 FORMAT (3x,i3,3x,a,2x,f4.0,2x,a,4x,i2,7x,i1) 
      END SUBROUTINE show_font                      
!*****7**************************************************************** 
      SUBROUTINE show_annotation 
!+                                                                      
!     Show annotations                                                  
!-                                                                      
      USE prompt_mod 
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER(1) just 
      INTEGER i, il, len_str 
      LOGICAL noann 
!                                                                       
      noann = .true. 
!                                                                       
      WRITE (output_io, 1000) iwin, iframe 
      DO i = 1, maxan 
      IF (antext (iwin, iframe, i) (1:3) .ne.'OFF') then 
         IF (anjust (iwin, iframe, i) .eq.if_left) just = 'L' 
         IF (anjust (iwin, iframe, i) .eq.if_right) just = 'R' 
         IF (anjust (iwin, iframe, i) .eq.if_centre) just = 'C' 
         il = min (len_str (antext (iwin, iframe, i) ), 15) 
         WRITE (output_io, 2000) i, anx (iwin, iframe, i), any (iwin,   &
         iframe, i), anangle (iwin, iframe, i), just, antext (iwin,     &
         iframe, i) (1:il)                                              
         IF (anx (iwin, iframe, i) .ne.antx (iwin, iframe, i) .or.any ( &
         iwin, iframe, i) .ne.anty (iwin, iframe, i) ) then             
            WRITE (output_io, 2100) antx (iwin, iframe, i), anty (iwin, &
            iframe, i)                                                  
         ENDIF 
         noann = .false. 
      ENDIF 
      ENDDO 
!                                                                       
      IF (noann) write (output_io, 3000) 
      WRITE (output_io, * ) 
!                                                                       
 1000 FORMAT     (' Annotation for window ',i3,' / frame  ',i3,//       &
     &      3x,'  #',6x,'x',13x,'y',8x,'angle  just  text',/,           &
     &      3x,69('-'))                                                 
 2000 FORMAT     (3x,i3,2(2x,g12.6),2x,f5.1,2x,a1,4x,a) 
 2100 FORMAT     (6x,2(2x,g12.6),14x,'Text position ..') 
 3000 FORMAT     ('   *** no annotations set ***') 
      END SUBROUTINE show_annotation                
!*****7**************************************************************** 
      SUBROUTINE show_cmap 
!+                                                                      
!     Show current color map                                            
!-                                                                      
      USE prompt_mod 
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER i, j 
!                                                                       
      WRITE (output_io, 1000) 
      WRITE (output_io, 2000) (i, (col_map (iwin, i, j), j = 1, 3),     &
      i = 1, maxcol)                                                    
!                                                                       
 1000 FORMAT     (' Current colour map values (RGB) : ',/) 
 2000 FORMAT     (4(3x,i3,':',3(1x,f3.1))) 
      END SUBROUTINE show_cmap                      
!*****7**************************************************************** 
      SUBROUTINE show_color 
!+                                                                      
!     Show current colour settings                                      
!-                                                                      
      USE prompt_mod 
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER i, j 
!                                                                       
      WRITE (output_io, 1000) 
      WRITE (output_io, 2000) (colour (iwin, 0, j), j = 1, 3) 
      WRITE (output_io, 2050) (colour (iwin, 6, j), j = 1, 3) 
!                                                                       
      DO i = 1, 15 
      WRITE (output_io, 2100) i, (colour (iwin, i, j), j = 1, 3) 
      ENDDO 
!                                                                       
 1000 FORMAT     (  ' Current colour settings (RGB) : ',/) 
 2000 FORMAT     (  '   Background  :',3(1x,f4.2)) 
 2050 FORMAT     (  '   Foreground  :',3(1x,f4.2),' -> Pen 6 !'/) 
 2100 FORMAT     (  '   Pen ',i2,'      :',3(1x,f4.2)) 
      END SUBROUTINE show_color                     
!*****7**************************************************************** 
      SUBROUTINE show_style 
!+                                                                      
!     Show information about plot styles (color, ..)                    
!-                                                                      
      USE errlist_mod 
      USE prompt_mod 
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER(40) cdummy 
      INTEGER i, j, ik 
      LOGICAL d2, d3 
!                                                                       
      INTEGER len_str 
      LOGICAL k_in_f 
!                                                                       
      d3 = .false. 
      d2 = .false. 
      DO i = 1, iz - 1 
      d3 = d3.or. (lni (i) .and.k_in_f (i) ) 
      d2 = d2.or. (.not.lni (i) .and.k_in_f (i) ) 
      ENDDO 
!                                                                       
!------ 2D data set parameters                                          
!                                                                       
      IF (d2) then 
         WRITE (output_io, 1000) iwin, iframe 
         DO ik = 1, iz - 1 
         IF (.not.lni (ik) .and.k_in_f (ik) ) then 
            cdummy = 'off' 
            IF (ilegend (iwin, iframe, ik) .eq.1) cdummy = clegend (    &
            iwin, iframe, ik)                                           
            WRITE (output_io, 1100) ik, ilinecol (iwin, iframe, ik),    &
            ilinetyp (iwin, iframe, ik), ilineart (iwin, iframe, ik),   &
            imarkcol (iwin, iframe, ik), imarktyp (iwin, iframe, ik),   &
            imarkmax (iwin, iframe, ik), ierr (iwin, iframe, ik),       &
            ierrcol (iwin, iframe, ik), linewid (iwin, iframe, ik),     &
            sizemark (iwin, iframe, ik), cdummy (1:min (len_str (cdummy)&
            , 7) )                                                      
         ENDIF 
         ENDDO 
         WRITE (output_io, * ) 
!                                                                       
         WRITE (output_io, 1500) iwin, iframe 
         DO ik = 1, iz - 1 
         IF (.not.lni (ik) .and.k_in_f (ik) ) then 
            IF (fillrange (iwin, iframe, ik, 1) .eq. - 9999) then 
               WRITE (output_io, 1610) ik, ifillcol (iwin, iframe, ik), &
               ifilltyp (iwin, iframe, ik)                              
            ELSEIF (fillrange (iwin, iframe, ik, 3) .eq. - 9999) then 
               WRITE (output_io, 1620) ik, ifillcol (iwin, iframe, ik), &
               ifilltyp (iwin, iframe, ik), (fillrange (iwin, iframe,   &
               ik, i), i = 1, 2)                                        
            ELSE 
               WRITE (output_io, 1600) ik, ifillcol (iwin, iframe, ik), &
               ifilltyp (iwin, iframe, ik), (fillrange (iwin, iframe,   &
               ik, i), i = 1, 2)                                        
            ENDIF 
         ENDIF 
         ENDDO 
         WRITE (output_io, * ) 
      ENDIF 
!                                                                       
!------ 3D data set parameters                                          
!                                                                       
      IF (d3) then 
         WRITE (output_io, 2000) iwin, iframe 
         DO ik = 1, iz - 1 
         IF (lni (ik) .and.k_in_f (ik) ) then 
            cdummy = 'off' 
            IF (ilegend (iwin, iframe, ik) .eq.1) cdummy = clegend (    &
            iwin, iframe, ik)                                           
            WRITE (output_io, 2200) ik, hlineart (iwin, iframe, ik),    &
            (hlinecol (iwin, iframe, ik, j), j = 1, maxhl), (hlinetyp ( &
            iwin, iframe, ik, j), j = 1, maxhl), linewid (iwin, iframe, &
            ik), cdummy (1:min (len_str (cdummy), 10) )                 
         ENDIF 
         ENDDO 
         WRITE (output_io, * ) 
      ENDIF 
!                                                                       
 1000 FORMAT(' Plot style for 2D data sets (Window ',i3,                &
     &      ' / frame ',i3,') : ',//,                                   &
     &      3x,'d.set  lcol  ltyp  lart  mcol  mtyp  ptyp  etyp  ecol', &
     &      '   lwid   msiz   caption',/,3x,77('-'))                    
 1100 FORMAT( 3x,i3,1x,8(4x,i2),3x,f5.2,2x,f5.2,3x,a) 
 1500 FORMAT(' Fill style for 2D data sets (Window ',i3,                &
     &      ' / frame ',i3,') : ',//,                                   &
     &      3x,'d.set  fcol  ftyp',10x,'range-x',19x,'range y',/,       &
     &      3x,77('-'))                                                 
 1600 FORMAT( 3x,i3,1x,2(4x,i2),3x,4(g12.4,1x)) 
 1610 FORMAT( 3x,i3,1x,2(4x,i2),11x,'not set ..') 
 1620 FORMAT( 3x,i3,1x,2(4x,i2),3x,2(g12.4,1x),8x,'not set ..') 
 2000 FORMAT(' Plot style for 3D data sets (Window ',i3,                &
     &      ' / frame ',i3,') : ',//,                                   &
     &     3x,'d.set',2x,'hart',5x,'hcol',9x,'htyp',8x,'lwid',          &
     &      4x,'caption',/,3x,76('-'))                                  
 2200 FORMAT(3x,i3,5x,i2,4x,4(i1,','),i1,4x,4(i1,','),i1,3x,f5.2,5x,a) 
!                                                                       
      END SUBROUTINE show_style                     
