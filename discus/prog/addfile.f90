!*****7*****************************************************************
!                                                                       
      SUBROUTINE do_addfile (string, laenge) 
!-                                                                      
!     Adds two files. the names plus a scale factor are read from       
!     string. Routine uses arrays 'dsi' and 'tcsf' => results of        
!     Fourier might be destroyed !                                      
!+                                                                      
      USE config_mod 
      USE diffuse_mod 
      USE errlist_mod 
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER maxw, if1, if2, if3 
      PARAMETER (maxw = 15, if1 = 1, if2 = 2, if3 = 3) 
!                                                                       
      CHARACTER ( * ) string 
      CHARACTER(3) ftyp 
      CHARACTER(7) stat 
      CHARACTER(1024) datei, line, fname 
      CHARACTER(1024) cpara (maxw), ccpara (maxw) 
      INTEGER lpara (maxw) 
      INTEGER nx (2), ny (2), i, j, ianz, laenge 
      INTEGER line_no 
      LOGICAL lcont, lread 
      REAL xmin (2), xmax (2), ymin (2), ymax (2), scale 
      REAL x1, x2, y1, y2, z1, z2, zz1, zz2 
      REAL werte (maxw) 
!                                                                       
!                                                                       
      ier_num = 0 
      ier_typ = ER_NONE 
      line_no = 0 
!                                                                       
      CALL get_params (string, ianz, cpara, lpara, maxw, laenge) 
      IF (ier_num.ne.0) then 
         RETURN 
      ELSEIF (ianz.ge.4) then 
!                                                                       
!     --Build filenames                                                 
!                                                                       
!       First parameter (plus format parameters) is filename 1          
!                                                                       
         CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
         IF (ier_num.ne.0) then 
            RETURN 
         ENDIF 
         ccpara (1) = cpara (1) 
!                                                                       
         DO i = 1, ianz - 1 
         cpara (i) = cpara (i + 1) 
         lpara (i) = lpara (i + 1) 
         ENDDO 
         cpara (ianz) = ' ' 
         lpara (ianz) = 0 
         ianz = ianz - 1 
!                                                                       
!       Second parameter (plus format parameters) is filename 2         
!                                                                       
         CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
         IF (ier_num.ne.0) then 
            RETURN 
         ENDIF 
         ccpara (2) = cpara (1) 
!                                                                       
         DO i = 1, ianz - 1 
         cpara (i) = cpara (i + 1) 
         lpara (i) = lpara (i + 1) 
         ENDDO 
         cpara (ianz) = ' ' 
         lpara (ianz) = 0 
         ianz = ianz - 1 
!                                                                       
!       Third parameter is multiplier                                   
!                                                                       
         CALL ber_params (1, cpara, lpara, werte, maxw) 
         IF (ier_num.ne.0) then 
            RETURN 
         ENDIF 
         scale = werte (1) 
!                                                                       
         DO i = 1, ianz - 1 
         cpara (i) = cpara (i + 1) 
         lpara (i) = lpara (i + 1) 
         ENDDO 
         cpara (ianz) = ' ' 
         lpara (ianz) = 0 
         ianz = ianz - 1 
!                                                                       
!       Fourth parameter (plus format parameters) is filename 3         
!                                                                       
         CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
         IF (ier_num.ne.0) then 
            RETURN 
         ENDIF 
         ccpara (4) = cpara (1) 
!                                                                       
         lcont = ccpara (1) .eq.ccpara (2) .or.ccpara (1) .eq.ccpara (4) 
!                                                                       
         IF (lcont) then 
            datei = 'discus.scr' 
            fname = ccpara (1) 
         ELSE 
            datei = ccpara (1) 
         ENDIF 
         stat = 'unknown' 
         lread = .false. 
         CALL oeffne (if1, datei, stat, lread) 
         IF (ier_num.ne.0) then 
            RETURN 
         ENDIF 
!                                                                       
         datei = ccpara (2) 
         stat = 'old' 
         lread = .true. 
         CALL oeffne (if2, datei, stat, lread) 
         IF (ier_num.ne.0) then 
            RETURN 
         ENDIF 
!                                                                       
         datei = ccpara (4) 
         CALL oeffne (if3, datei, stat, lread) 
         IF (ier_num.ne.0) then 
            RETURN 
         ENDIF 
!                                                                       
!     --There are one or two paramters left, depending on the           
!       presence or absence of file type specifier                      
!                                                                       
         IF (ianz.eq.2) then 
            ftyp = cpara (2) 
            laenge = lpara (2) 
            CALL rem_bl (ftyp, laenge) 
         ELSE 
            ftyp = 'ni' 
         ENDIF 
!                                                                       
         IF (ier_num.eq.0) then 
!                                                                       
            IF (ftyp.eq.'ni') then 
               ier_num = - 3 
               ier_typ = ER_IO 
               READ (if2, *, err = 900) nx (1), ny (1) 
               READ (if2, *, err = 900) xmin (1), xmax (1), ymin (1),   &
               ymax (1)                                                 
               READ (if3, *, err = 900) nx (2), ny (2) 
               READ (if3, *, err = 900) xmin (2), xmax (2), ymin (2),   &
               ymax (2)                                                 
               ier_num = 0 
               ier_typ = ER_NONE 
               IF (nx (1) .eq.nx (2) .and.ny (1) .eq.ny (2) .and.xmin ( &
               1) .eq.xmin (2) .and.ymin (1) .eq.ymin (2) .and.xmax (1) &
               .eq.xmax (2) .and.ymax (1) .eq.ymax (2) ) then           
                  WRITE (if1, * ) nx (1), ny (1) 
                  WRITE (if1, * ) xmin (1), xmax (1), ymin (1), ymax (1) 
                  DO i = 1, ny (1) 
                  ier_num = - 3 
                  ier_typ = ER_IO 
                  READ (if2, *, err = 900) (dsi (j), j = 1, nx (1) ) 
                  DO j = 1, nx (1) 
                  tcsf (j) = cmplx (dble (dsi (j) ), 0.0d0) 
                  ENDDO 
                  READ (if3, *, err = 900) (dsi (j), j = 1, nx (1) ) 
                  ier_num = 0 
                  ier_typ = ER_NONE 
                  DO j = 1, nx (1) 
                  dsi (j) = real (tcsf (j) ) + scale * dsi (j) 
                  ENDDO 
                  WRITE (if1, 3010) (dsi (j), j = 1, nx (1) ) 
                  WRITE (if1, * ) 
                  ENDDO 
               ELSE 
                  ier_num = - 22 
                  ier_typ = ER_APPL 
               ENDIF 
            ELSEIF (ftyp.eq.'gnu') then 
               ier_num = - 3 
               ier_typ = ER_IO 
   10          CONTINUE 
               line_no = line_no + 1 
               READ (if2, 1000, end = 800, err = 900) line 
               IF (line.eq.' '.or.line.eq.char (13) ) then 
                  READ (if3, 1000, end = 900, err = 900) line 
                  IF (line.eq.' '.or.line.eq.char (13) ) then 
                     WRITE (if1, * ) 
                  ELSE 
                     ier_num = - 23 
                     ier_typ = ER_APPL 
                     GOTO 900 
                  ENDIF 
               ELSE 
                  READ (line, *, end = 800, err = 900) x1, y1, z1, zz1 
                  READ (if3, *, end = 900, err = 900) x2, y2, z2, zz2 
                  IF (x1.eq.x2.and.y1.eq.y2.and.zz1.eq.zz2) then 
                     WRITE (if1, 3020) x1, y1, z1 + scale * z2, zz1 
                  ELSE 
                     ier_num = - 24 
                     ier_typ = ER_APPL 
                     GOTO 900 
                  ENDIF 
               ENDIF 
               GOTO 10 
            ELSEIF (ftyp.eq.'1d') then 
               ier_num = - 3 
               ier_typ = ER_IO 
   20          CONTINUE 
               line_no = line_no + 1 
               READ (if2, *, end = 800, err = 900) x1, y1 
               READ (if3, *, end = 900, err = 900) x2, y2 
               IF (x1.eq.x2) then 
                  WRITE (if1, 3030) x1, y1 + scale * y2 
               ELSE 
                  ier_num = - 25 
                  ier_typ = ER_APPL 
                  GOTO 900 
               ENDIF 
               GOTO 20 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
               GOTO 900 
            ENDIF 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
         GOTO 900 
      ENDIF 
  800 CONTINUE 
      ier_num = 0 
      ier_typ = ER_NONE 
  900 CONTINUE 
      CLOSE (if1) 
      CLOSE (if2) 
      CLOSE (if3) 
      IF (lcont) then 
         CALL do_rename_file ('discus.scr', fname) 
      ENDIF 
!                                                                       
      IF (ier_num.ne.0) then 
         IF (ftyp.eq.'ni') then 
            WRITE (output_io, * ) ' Index y: ', i, ' Index x: ', j 
         ELSE 
            WRITE (output_io, * ) ' Line number: ', line_no 
         ENDIF 
      ENDIF 
!                                                                       
 1000 FORMAT    (a) 
 3000 FORMAT    (1x,2i4/1x,4e12.5e3) 
 3010 FORMAT    (5(1x,e11.5)) 
 3020 FORMAT    (4(1x,e11.5)) 
 3030 FORMAT    (2(1x,e11.5)) 
!                                                                       
      END SUBROUTINE do_addfile                     
