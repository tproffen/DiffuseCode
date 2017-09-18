      SUBROUTINE kuplot_ersetz_para (ikl, iklz, string, ll, ww, maxw, ianz)
!       Replaces a substring in an expression by the value of the       
!       appropriate parameter. Modified version for KUPLOT.             
!+                                                                      
      USE errlist_mod 
      USE param_mod 
      USE kuplot_config 
      USE kuplot_mod 
      USE variable_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER,                    INTENT(IN   ) :: ikl
      INTEGER,                    INTENT(IN   ) :: iklz
      CHARACTER (LEN=*),          INTENT(OUT  ) :: string 
      INTEGER,                    INTENT(OUT  ) :: ll
      INTEGER,                    INTENT(IN   ) :: maxw
      REAL   , DIMENSION(1:maxw), INTENT(IN   ) :: ww
      INTEGER,                    INTENT(IN   ) :: ianz
!
      INTEGER mmaxw
      PARAMETER (mmaxw = 3) 
!                                                                       
      CHARACTER(1024) zeile 
      INTEGER i, laenge
      INTEGER lcomm, ltyp 
      INTEGER idummy 
      INTEGER kpara (mmaxw) 
!                                                                       
      INTEGER length_com 
!                                                                       
      laenge = ll 
      ltyp = 1 
      zeile = ' ' 
      DO i = 1, maxw 
      kpara (i) = nint (ww (i) ) 
      ENDDO 
!                                                                       
      lcomm = length_com (string, ikl) 
!                                                                       
      IF (lcomm.eq.1) then 
         IF (ikl.gt.2) zeile (1:ikl - 1) = string (1:ikl - 2) 
!                                                                       
!-------- general integer parameter i[n]                                
!                                                                       
         IF (string (ikl - 1:ikl - 1) .eq.'i') then 
            IF (ianz.eq.1) then 
               IF (0.le.kpara (1) .and.kpara (1) .le.MAXPAR) then 
                  WRITE (zeile (ikl - 1:ikl + 13) , '(i15)') inpara (   &
                  kpara (1) )                                           
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
            ENDIF 
!                                                                       
!-------- program information n[n]                                      
!                                                                       
         ELSEIF (string (ikl - 1:ikl - 1) .eq.'n') then 
            IF (ianz.eq.1) then 
               IF (1.le.kpara (1) .and.kpara (1) .le.5) then 
                  IF (kpara (1) .eq.1) then 
                     idummy = iz - 1 
                     WRITE (zeile (ikl - 1:ikl + 13) , '(i15)') idummy 
                  ELSEIF (kpara (1) .eq.2) then 
                     WRITE (zeile (ikl - 1:ikl + 13) , '(i15)')         &
                     maxkurvtot                                         
                  ELSEIF (kpara (1) .eq.3) then 
                     WRITE (zeile (ikl - 1:ikl + 13) , '(i15)')         &
                     maxframe                                           
                  ELSEIF (kpara (1) .eq.4) then 
                     WRITE (zeile (ikl - 1:ikl + 13) , '(i15)') maxan 
                  ELSEIF (kpara (1) .eq.5) then 
                     WRITE (zeile (ikl - 1:ikl + 13) , '(i15)') maxbond 
                  ENDIF 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
            ENDIF 
!                                                                       
!-------- general real parameter r[n]                                   
!                                                                       
         ELSEIF (string (ikl - 1:ikl - 1) .eq.'r') then 
            IF (ianz.eq.1) then 
               IF (0.le.kpara (1) .and.kpara (1) .le.MAXPAR) then 
                  WRITE (zeile (ikl - 1:ikl + 13) , '(e15.8e2)') rpara (&
                  kpara (1) )                                           
                  zeile (ikl + 10:ikl + 10) = 'e' 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
            ENDIF 
!                                                                       
!-------- Fit parameter p[n] and error s[n]                             
!                                                                       
         ELSEIF (string (ikl - 1:ikl - 1) .eq.'p') then 
            IF (ianz.eq.1) then 
               IF (1.le.kpara (1) .and.kpara (1) .le.MAXPARA) then 
                  WRITE (zeile (ikl - 1:ikl + 13) , '(e15.8e2)') p (    &
                  kpara (1) )                                           
                  zeile (ikl + 10:ikl + 10) = 'e' 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
            ENDIF 
!                                                                       
         ELSEIF (string (ikl - 1:ikl - 1) .eq.'s') then 
            IF (ianz.eq.1) then 
               IF (1.le.kpara (1) .and.kpara (1) .le.MAXPARA) then 
                  WRITE (zeile (ikl - 1:ikl + 13) , '(e15.8e2)') dp (   &
                  kpara (1) )                                           
                  zeile (ikl + 10:ikl + 10) = 'e' 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
            ENDIF 
!                                                                       
!-------- data arrays x[n,m], y[n,m], z[n,mx,my]                        
!                                                                       
         ELSEIF (string (ikl - 1:ikl - 1) .eq.'x') then 
            IF (ianz.eq.2) then 
               IF (1.le.kpara (1) .and.kpara (1) .le. (iz - 1)          &
               .and.1.le.kpara (2) .and.kpara (2) .le.len (kpara (1) ) )&
               then                                                     
                  WRITE (zeile (ikl - 1:ikl + 13) , '(e15.8e2)') x (    &
                  offxy (kpara (1) - 1) + kpara (2) )                   
                  zeile (ikl + 10:ikl + 10) = 'e' 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
            ENDIF 
!                                                                       
         ELSEIF (string (ikl - 1:ikl - 1) .eq.'y') then 
            IF (ianz.eq.2) then 
               IF (1.le.kpara (1) .and.kpara (1) .le. (iz - 1)          &
               .and.1.le.kpara (2) .and.kpara (2) .le.len (kpara (1) ) )&
               then                                                     
                  WRITE (zeile (ikl - 1:ikl + 13) , '(e15.8e2)') y (    &
                  offxy (kpara (1) - 1) + kpara (2) )                   
                  zeile (ikl + 10:ikl + 10) = 'e' 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
            ENDIF 
!                                                                       
         ELSEIF (string (ikl - 1:ikl - 1) .eq.'z') then 
            IF (ianz.eq.3) then 
               IF (1.le.kpara (1) .and.kpara (1) .le. (iz - 1)          &
               .and.1.le.kpara (2) .and.kpara (2) .le.nx (kpara (1) )   &
               .and.1.le.kpara (3) .and.kpara (3) .le.ny (kpara (1) ) ) &
               then                                                     
                  WRITE (zeile (ikl - 1:ikl + 13) , '(e15.8e2)') z (    &
                  offz (kpara (1) - 1) + (kpara (2) - 1) * ny (kpara (1)&
                  ) + kpara (3) )                                       
                  zeile (ikl + 10:ikl + 10) = 'e' 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
            ENDIF 
!                                                                       
         ELSE 
            ier_num = - 2 
            ier_typ = ER_FORT 
         ENDIF 
      ELSEIF (lcomm.eq.2) then 
         IF (ikl.gt.3) zeile (1:ikl - 2) = string (1:ikl - 3) 
!                                                                       
!-------- data arrays dx[n,m], dy[n,m]                                  
!                                                                       
         IF (string (ikl - 2:ikl - 1) .eq.'dx') then 
            IF (ianz.eq.2) then 
               IF (1.le.kpara (1) .and.kpara (1) .le. (iz - 1)          &
               .and.1.le.kpara (2) .and.kpara (2) .le.len (kpara (1) ) )&
               then                                                     
                  WRITE (zeile (ikl - 2:ikl + 13) , '(e15.8e2)') dx (   &
                  offxy (kpara (1) - 1) + kpara (2) )                   
                  zeile (ikl + 9:ikl + 9) = 'e' 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
            ENDIF 
!                                                                       
         ELSEIF (string (ikl - 2:ikl - 1) .eq.'dy') then 
            IF (ianz.eq.2) then 
               IF (1.le.kpara (1) .and.kpara (1) .le. (iz - 1)          &
               .and.1.le.kpara (2) .and.kpara (2) .le.len (kpara (1) ) )&
               then                                                     
                  WRITE (zeile (ikl - 2:ikl + 13) , '(e15.8e2)') dy (   &
                  offxy (kpara (1) - 1) + kpara (2) )                   
                  zeile (ikl + 9:ikl + 9) = 'e' 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
            ENDIF 
!                                                                       
!-------- data set dimensions nx[n], ny[n], np[n], ni[n]                
!                                                                       
         ELSEIF (string (ikl - 2:ikl - 1) .eq.'nx') then 
            IF (ianz.eq.1) then 
               IF (1.le.kpara (1) .and.kpara (1) .le. (iz - 1) ) then 
                  IF (lni (kpara (1) ) ) then 
                     WRITE(zeile(ikl-2:ikl+13),'(i15)') nx(kpara(1))
                  ELSE
                     ier_num = -62 
                     ier_typ = ER_APPL 
                  ENDIF 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
            ENDIF 
!                                                                       
         ELSEIF (string (ikl - 2:ikl - 1) .eq.'ny') then 
            IF (ianz.eq.1) then 
               IF (1.le.kpara (1) .and.kpara (1) .le. (iz - 1) ) then 
                  IF (lni (kpara (1) ) ) then 
                     WRITE(zeile(ikl-2:ikl+13),'(i15)') ny(kpara(1))
                  ELSE
                     ier_num = -62 
                     ier_typ = ER_APPL 
                  ENDIF 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
            ENDIF 
!                                                                       
         ELSEIF (string (ikl - 2:ikl - 1) .eq.'np') then 
            IF (ianz.eq.1) then 
               IF (1.le.kpara (1) .and.kpara (1) .le. (iz - 1) ) then 
                  IF (lni (kpara (1) ) ) then 
                     WRITE(zeile(ikl-2:ikl+13),'(i15)') ny(kpara(1))*ny(kpara(1))
                  ELSE
                     WRITE(zeile(ikl-2:ikl+13),'(i15)') len(kpara(1))
                  ENDIF 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
            ENDIF 
!                                                                       
         ELSEIF (string (ikl - 2:ikl - 1) .eq.'ni') then 
            IF (ianz.eq.1) then 
               IF (1.le.kpara (1) .and.kpara (1) .le. (iz - 1) ) then 
                  IF (lni (kpara (1) ) ) then 
                     WRITE (zeile (ikl - 2:ikl + 13) , '(i15)') 1 
                  ELSE 
                     WRITE (zeile (ikl - 2:ikl + 13) , '(i15)') 0 
                  ENDIF 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
            ENDIF 
         ELSE 
            ier_num = - 2 
            ier_typ = ER_FORT 
         ENDIF 
!                                                                       
      ELSEIF (lcomm.eq.3) then 
         IF (ikl.gt.4) zeile (1:ikl - 3) = string (1:ikl - 4) 
!                                                                       
!-------- result array res[n]                                           
!                                                                       
         IF (string (ikl - 3:ikl - 1) .eq.'res') then 
            IF (ikl.gt.4) zeile (1:ikl - 3) = string (1:ikl - 4) 
            IF (ianz.eq.1) then 
               IF (0.le.kpara (1) .and.kpara (1) .le.MAXPAR_RES) then 
                  WRITE (zeile (ikl - 3:ikl + 13) , '(e15.8e2)')        &
                  res_para (kpara (1) )                                 
                  zeile (ikl + 8:ikl + 8) = 'e' 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
            ENDIF 
         ELSE 
            ier_num = - 2 
            ier_typ = ER_FORT 
         ENDIF 
!                                                                       
      ELSEIF (lcomm.eq.4) then 
         IF (ikl.gt.5) zeile (1:ikl - 4) = string (1:ikl - 5) 
!                                                                       
!-------- data set dimensions xmin[n], xmax[n], ..                      
!                                                                       
         IF (string (ikl - 4:ikl - 1) .eq.'xmin') then 
            IF (ianz.eq.1) then 
               IF (1.le.kpara (1) .and.kpara (1) .le. (iz - 1) ) then 
                  CALL get_extrema 
                  WRITE (zeile (ikl - 4:ikl + 13) , '(e15.8e2)') xmin ( &
                  kpara (1) )                                           
                  zeile (ikl + 7:ikl + 7) = 'e' 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
            ENDIF 
!                                                                       
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'xmax') then 
            IF (ianz.eq.1) then 
               IF (1.le.kpara (1) .and.kpara (1) .le. (iz - 1) ) then 
                  CALL get_extrema 
                  WRITE (zeile (ikl - 4:ikl + 13) , '(e15.8e2)') xmax ( &
                  kpara (1) )                                           
                  zeile (ikl + 7:ikl + 7) = 'e' 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
            ENDIF 
!                                                                       
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'ymin') then 
            IF (ianz.eq.1) then 
               IF (1.le.kpara (1) .and.kpara (1) .le. (iz - 1) ) then 
                  CALL get_extrema 
                  WRITE (zeile (ikl - 4:ikl + 13) , '(e15.8e2)') ymin ( &
                  kpara (1) )                                           
                  zeile (ikl + 7:ikl + 7) = 'e' 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
            ENDIF 
!                                                                       
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'ymax') then 
            IF (ianz.eq.1) then 
               IF (1.le.kpara (1) .and.kpara (1) .le. (iz - 1) ) then 
                  CALL get_extrema 
                  WRITE (zeile (ikl - 4:ikl + 13) , '(e15.8e2)') ymax ( &
                  kpara (1) )                                           
                  zeile (ikl + 7:ikl + 7) = 'e' 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
            ENDIF 
!                                                                       
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'zmin') then 
            IF (ianz.eq.1) then 
               IF (1.le.kpara (1) .and.kpara (1) .le. (iz - 1) .and.lni &
               (kpara (1) ) ) then                                      
                  CALL get_extrema 
                  WRITE (zeile (ikl - 4:ikl + 13) , '(e15.8e2)') zmin ( &
                  kpara (1) )                                           
                  zeile (ikl + 7:ikl + 7) = 'e' 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
            ENDIF 
!                                                                       
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'zmax') then 
            IF (ianz.eq.1) then 
               IF (1.le.kpara (1) .and.kpara (1) .le. (iz - 1) .and.lni &
               (kpara (1) ) ) then                                      
                  CALL get_extrema 
                  WRITE (zeile (ikl - 4:ikl + 13) , '(e15.8e2)') zmax ( &
                  kpara (1) )                                           
                  zeile (ikl + 7:ikl + 7) = 'e' 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
            ENDIF 
!                                                                       
!------ - Scale factors for graphic devices                             
!                                                                       
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'size') then 
            IF (ianz.eq.1) then 
               IF (1.eq.kpara (1) ) then 
                  WRITE (zeile (ikl - 4:ikl + 13) , '(e15.8e2)') dev_sf &
                  (iwin, x11)                                           
                  zeile (ikl + 7:ikl + 7) = 'e' 
               ELSEIF (2.eq.kpara (1) ) then 
                  WRITE (zeile (ikl - 4:ikl + 13) , '(e15.8e2)') dev_sf &
                  (iwin, ps)                                            
                  zeile (ikl + 7:ikl + 7) = 'e' 
               ELSEIF (3.eq.kpara (1) ) then 
                  WRITE (zeile (ikl - 4:ikl + 13) , '(e15.8e2)') dev_sf &
                  (iwin, pic)                                           
                  zeile (ikl + 7:ikl + 7) = 'e' 
               ELSEIF (4.eq.kpara (1) ) then 
                  WRITE (zeile (ikl - 4:ikl + 13) , '(e15.8e2)') dev_sf &
                  (iwin, png)                                           
                  zeile (ikl + 7:ikl + 7) = 'e' 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
            ENDIF 
!                                                                       
!------ - Current plot window dimensions                                
!                                                                       
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'pwin') then 
            CALL skalieren 
            IF (ianz.eq.1) then 
               IF (1.eq.kpara (1) ) then 
                  WRITE (zeile (ikl - 4:ikl + 13) , '(e15.8e2)') pex (  &
                  iwin, iframe, 1)                                      
                  zeile (ikl + 7:ikl + 7) = 'e' 
               ELSEIF (2.eq.kpara (1) ) then 
                  WRITE (zeile (ikl - 4:ikl + 13) , '(e15.8e2)') pex (  &
                  iwin, iframe, 2)                                      
                  zeile (ikl + 7:ikl + 7) = 'e' 
               ELSEIF (3.eq.kpara (1) ) then 
                  WRITE (zeile (ikl - 4:ikl + 13) , '(e15.8e2)') pey (  &
                  iwin, iframe, 1)                                      
                  zeile (ikl + 7:ikl + 7) = 'e' 
               ELSEIF (4.eq.kpara (1) ) then 
                  WRITE (zeile (ikl - 4:ikl + 13) , '(e15.8e2)') pey (  &
                  iwin, iframe, 2)                                      
                  zeile (ikl + 7:ikl + 7) = 'e' 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
            ENDIF 
!                                                                       
!------ - color map variables cmap[n,m], cmax[n]                        
!                                                                       
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'cmap') then 
            IF (ianz.eq.2) then 
               IF (1.le.kpara (1) .and.kpara (1)                        &
               .le.maxcol.and.1.le.kpara (2) .and.kpara (2) .le.3) then 
                  WRITE (zeile (ikl - 4:ikl + 13) , '(e15.8e2)')        &
                  col_map (iwin, kpara (1) , kpara (2) )                
                  zeile (ikl + 7:ikl + 7) = 'e' 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
            ENDIF 
!                                                                       
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'cmax') then 
            IF (ianz.eq.1) then 
               IF (kpara (1) .eq.1) then 
                  WRITE (zeile (ikl - 4:ikl + 13) , '(i15)') maxcol 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
            ENDIF 
!                                                                       
!------ - Drawing of axes relates variables                             
!                                                                       
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'axis') then 
            IF (ianz.eq.2) then 
               IF (1.le.kpara (1) .and.kpara (1) .le.6.and.1.le.kpara ( &
               2) .and.kpara (2) .le.2) then                            
!                                                                       
                  IF (kpara (1) .eq.1) then 
                     WRITE (zeile (ikl - 4:ikl + 13) , '(e15.8e2)')     &
                     lab_angle (iwin, iframe, kpara (2) )               
                     zeile (ikl + 7:ikl + 7) = 'e' 
                  ELSEIF (kpara (1) .eq.2) then 
                     WRITE (zeile (ikl - 4:ikl + 13) , '(e15.8e2)')     &
                     tick_ma_h (iwin, iframe, kpara (2) )               
                     zeile (ikl + 7:ikl + 7) = 'e' 
                  ELSEIF (kpara (1) .eq.3) then 
                     WRITE (zeile (ikl - 4:ikl + 13) , '(e15.8e2)')     &
                     tick_mi_h (iwin, iframe, kpara (2) )               
                     zeile (ikl + 7:ikl + 7) = 'e' 
                  ELSEIF (kpara (1) .eq.4) then 
                     WRITE (zeile (ikl - 4:ikl + 13) , '(i15)')         &
                     tick_nsub (iwin, iframe, kpara (2) )               
                  ELSEIF (kpara (1) .eq.5) then 
                     WRITE (zeile (ikl - 4:ikl + 13) , '(e15.8e2)')     &
                     lab_d (iwin, iframe, kpara (2) )                   
                     zeile (ikl + 7:ikl + 7) = 'e' 
                  ELSEIF (kpara (1) .eq.6) then 
                     WRITE (zeile (ikl - 4:ikl + 13) , '(e15.8e2)')     &
                     ax_d (iwin, iframe, kpara (2) )                    
                     zeile (ikl + 7:ikl + 7) = 'e' 
                  ENDIF 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
            ENDIF 
!                                                                       
         ELSE 
            ier_num = - 2 
            ier_typ = ER_FORT 
         ENDIF 
!
      ELSEIF(lcomm==8) THEN                                                                  
!                                                                       
         IF (string (ikl - 8:ikl - 1) .eq.'ref_para') then 
            IF (ianz.eq.1) then 
               IF (ikl.gt.lcomm + 1) zeile (1:ikl - lcomm - 1) = string &
               (1:ikl - lcomm - 1)                                      
               IF (0.lt.kpara(1).and.kpara(1).le.MAXPAR_REF   ) then 
                  WRITE (zeile (ikl - 8:ikl + 13) , '(e15.8e2)')        &
                  ref_para (kpara(1))
                  zeile (ikl + 3:ikl + 3) = 'e' 
               ELSE 
                  ier_num = -133 
                  ier_typ = ER_APPL 
                  RETURN 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
         ELSE 
            ier_num = - 2 
            ier_typ = ER_FORT 
         ENDIF 
!
      ELSE 
         ier_num = - 2 
         ier_typ = ER_FORT 
      ENDIF 
      IF (ier_num.eq.0) then 
         ll = laenge+15 - ltyp - (iklz - ikl + 1) 
         IF (iklz + 1.le.laenge) zeile (ikl + 14:ll) = string (iklz + 1:&
         laenge)                                                        
         string = zeile 
         CALL rem_bl (string, ll) 
      ENDIF 
      END SUBROUTINE kuplot_ersetz_para                    
!*****7*****************************************************************
      SUBROUTINE kuplot_upd_para (ctype, ww, maxw, wert, ianz) 
!-                                                                      
!       updates the parameter spezified by ctype, index ww  to the      
!       new value of wert                                               
!+                                                                      
      USE errlist_mod 
      USE param_mod 
      USE prompt_mod 
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER (LEN=*),          INTENT(IN) :: ctype 
      INTEGER,                    INTENT(IN) :: maxw
      INTEGER,                    INTENT(IN) :: ianz 
      INTEGER, DIMENSION(1:MAXW), INTENT(IN) :: ww
      REAL   ,                    INTENT(IN) :: wert 
!
      INTEGER idummy 
!                                                                       
!------ integer parameter i[n]                                          
!                                                                       
      IF (ctype.eq.'i') then 
         IF (ianz.eq.1) then 
            IF (0.le.ww (1) .and.ww (1) .le.MAXPAR) then 
               inpara (ww (1) ) = int (wert) 
            ELSE 
               ier_num = - 8 
               ier_typ = ER_FORT 
            ENDIF 
         ELSE 
            ier_num = - 13 
            ier_typ = ER_FORT 
         ENDIF 
!                                                                       
!------ real parameter r[n]                                             
!                                                                       
      ELSEIF (ctype.eq.'r') then 
         IF (ianz.eq.1) then 
            IF (0.le.ww (1) .and.ww (1) .le.MAXPAR) then 
               rpara (ww (1) ) = wert 
            ELSE 
               ier_num = - 8 
               ier_typ = ER_FORT 
            ENDIF 
         ELSE 
            ier_num = - 13 
            ier_typ = ER_FORT 
         ENDIF 
!                                                                       
!------ fit parameter p[n]                                              
!                                                                       
      ELSEIF (ctype.eq.'p') then 
         IF (ianz.eq.1) then 
            IF (1.le.ww (1) .and.ww (1) .le.MAXPARA) then 
               p (ww (1) ) = wert 
            ELSE 
               ier_num = - 8 
               ier_typ = ER_FORT 
            ENDIF 
         ELSE 
            ier_num = - 13 
            ier_typ = ER_FORT 
         ENDIF 
!                                                                       
!------ Number of loaded data sets n[1] - to clear some                 
!                                                                       
      ELSEIF (ctype.eq.'n') then 
         IF (ianz.eq.1) then 
            IF (ww (1) .eq.1) then 
               idummy = nint (wert) 
               IF (idummy.gt.0.and.idummy.le. (iz - 1) ) then 
                  WRITE (output_io, 2000) iz - idummy - 1 
                  iz = idummy + 1 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
            ELSE 
               ier_num = - 8 
               ier_typ = ER_FORT 
            ENDIF 
         ELSE 
            ier_num = - 13 
            ier_typ = ER_FORT 
         ENDIF 
!                                                                       
!------ data arrays x[n,m], y[n,m], z[n,nx,ny], dx[n,m], dy[n,m]        
!                                                                       
      ELSEIF (ctype.eq.'x') then 
         IF (ianz.eq.2) then 
            IF (1.le.ww (1) .and.ww (1) .le. (iz - 1) .and.1.le.ww (2)  &
            .and.ww (2) .le.len (ww (1) ) ) then                        
               x (offxy (ww (1) - 1) + ww (2) ) = wert 
            ELSE 
               ier_num = - 8 
               ier_typ = ER_FORT 
            ENDIF 
         ELSE 
            ier_num = - 13 
            ier_typ = ER_FORT 
         ENDIF 
!                                                                       
      ELSEIF (ctype.eq.'dx') then 
         IF (ianz.eq.2) then 
            IF (1.le.ww (1) .and.ww (1) .le. (iz - 1) .and.1.le.ww (2)  &
            .and.ww (2) .le.len (ww (1) ) ) then                        
               dx (offxy (ww (1) - 1) + ww (2) ) = wert 
            ELSE 
               ier_num = - 8 
               ier_typ = ER_FORT 
            ENDIF 
         ELSE 
            ier_num = - 13 
            ier_typ = ER_FORT 
         ENDIF 
!                                                                       
      ELSEIF (ctype.eq.'y') then 
         IF (ianz.eq.2) then 
            IF (1.le.ww (1) .and.ww (1) .le. (iz - 1) .and.1.le.ww (2)  &
            .and.ww (2) .le.len (ww (1) ) ) then                        
               y (offxy (ww (1) - 1) + ww (2) ) = wert 
            ELSE 
               ier_num = - 8 
               ier_typ = ER_FORT 
            ENDIF 
         ELSE 
            ier_num = - 13 
            ier_typ = ER_FORT 
         ENDIF 
!                                                                       
      ELSEIF (ctype.eq.'dy') then 
         IF (ianz.eq.2) then 
            IF (1.le.ww (1) .and.ww (1) .le. (iz - 1) .and.1.le.ww (2)  &
            .and.ww (2) .le.len (ww (1) ) ) then                        
               dy (offxy (ww (1) - 1) + ww (2) ) = wert 
            ELSE 
               ier_num = - 8 
               ier_typ = ER_FORT 
            ENDIF 
         ELSE 
            ier_num = - 13 
            ier_typ = ER_FORT 
         ENDIF 
!                                                                       
      ELSEIF (ctype.eq.'z') then 
         IF (ianz.eq.3) then 
            IF (1.le.ww (1) .and.ww (1) .le. (iz - 1) .and.1.le.ww (2)  &
            .and.ww (2) .le.nx (ww (1) ) .and.1.le.ww (3) .and.ww (3)   &
            .le.ny (ww (1) ) ) then                                     
               z (offz (ww (1) - 1) + (ww (2) - 1) * ny (ww (1) )       &
               + ww (3) ) = wert                                        
            ELSE 
               ier_num = - 8 
               ier_typ = ER_FORT 
            ENDIF 
         ELSE 
            ier_num = - 13 
            ier_typ = ER_FORT 
         ENDIF 
!                                                                       
!------ scale factors for graphic devices                               
!                                                                       
      ELSEIF (ctype.eq.'size') then 
         IF (ianz.eq.1) then 
            IF (1.eq.ww (1) ) then 
               dev_sf (iwin, x11) = wert 
            ELSEIF (2.eq.ww (1) ) then 
               dev_sf (iwin, ps) = wert 
               dev_sf (iwin, vps) = wert 
            ELSEIF (3.eq.ww (1) ) then 
               dev_sf (iwin, pic) = wert 
               dev_sf (iwin, vpic) = wert 
            ELSEIF (4.eq.ww (1) ) then 
               dev_sf (iwin, png) = wert 
            ELSE 
               ier_num = - 8 
               ier_typ = ER_FORT 
            ENDIF 
         ELSE 
            ier_num = - 13 
            ier_typ = ER_FORT 
         ENDIF 
!                                                                       
!------ color map values cmap[], cmax[]                                 
!                                                                       
      ELSEIF (ctype.eq.'cmap') then 
         IF (ianz.eq.2) then 
            IF (1.le.ww (1) .and.ww (1) .le.maxcol.and.1.le.ww (2)      &
            .and.ww (2) .le.3) then                                     
               IF (wert.ge.0.0.and.wert.le.1.0) then 
                  col_map (iwin, ww (1), ww (2) ) = wert 
               ELSE 
                  ier_num = - 18 
                  ier_typ = ER_APPL 
               ENDIF 
            ELSE 
               ier_num = - 8 
               ier_typ = ER_FORT 
            ENDIF 
         ELSE 
            ier_num = - 13 
            ier_typ = ER_FORT 
         ENDIF 
!                                                                       
!------ axis drawing parameters                                         
!                                                                       
      ELSEIF (ctype.eq.'axis') then 
         IF (ianz.eq.2) then 
            IF (1.le.ww (1) .and.ww (1) .le.6.and.1.le.ww (2) .and.ww ( &
            2) .le.2) then                                              
               IF (ww (1) .eq.1) then 
                  lab_angle (iwin, iframe, ww (2) ) = wert 
               ELSEIF (ww (1) .eq.2) then 
                  tick_ma_h (iwin, iframe, ww (2) ) = wert 
               ELSEIF (ww (1) .eq.3) then 
                  tick_ma_h (iwin, iframe, ww (2) ) = wert 
               ELSEIF (ww (1) .eq.4) then 
                  tick_nsub (iwin, iframe, ww (2) ) = nint (wert) 
               ELSEIF (ww (1) .eq.5) then 
                  lab_d (iwin, iframe, ww (2) ) = wert 
               ELSEIF (ww (1) .eq.6) then 
                  ax_d (iwin, iframe, ww (2) ) = wert 
               ENDIF 
            ELSE 
               ier_num = - 8 
               ier_typ = ER_FORT 
            ENDIF 
         ELSE 
            ier_num = - 13 
            ier_typ = ER_FORT 
         ENDIF 
!                                                                       
!------ result array parameter res[n]                                   
!                                                                       
      ELSEIF (ctype.eq.'res') then 
         IF (ianz.eq.1) then 
            IF (0.le.ww (1) .and.ww (1) .le.MAXPAR_RES) then 
               res_para (ww (1) ) = wert 
            ELSE 
               ier_num = - 8 
               ier_typ = ER_FORT 
            ENDIF 
         ELSE 
            ier_num = - 13 
            ier_typ = ER_FORT 
         ENDIF 
      ELSEIF (ctype.eq.'ref_para') THEN
         IF (ianz.eq.1) then 
            IF (0.le.ww (1) .and.ww (1) .le.MAXPAR_REF) then 
               ref_para (ww (1) ) = wert 
            ELSE 
               ier_num = - 8 
               ier_typ = ER_FORT 
            ENDIF 
         ELSE 
            ier_num = - 13 
            ier_typ = ER_FORT 
            RETURN 
         ENDIF 
      ELSE 
         ier_num = - 2 
         ier_typ = ER_FORT 
      ENDIF 
!                                                                       
 2000 FORMAT     (' ------ > Last ',i3,' data set(s) deleted ...') 
      END SUBROUTINE kuplot_upd_para                       
!*****7***************************************************************  
      SUBROUTINE kuplot_calc_intr_spec (string, line, ikl, iklz, ww, laenge,   &
      lp)                                                               
!-                                                                      
!     These are special intrinsic function for the KUPLOT. Any          
!     intrinsic function that references KUPLOT specific values         
!     is found in this subroutine.                                      
!+                                                                      
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
      CHARACTER (LEN=*), INTENT(IN   ) :: string
      CHARACTER (LEN=*), INTENT(INOUT) :: line 
      INTEGER,           INTENT(IN)    :: ikl
      INTEGER,           INTENT(IN)    :: iklz
      INTEGER,           INTENT(IN)    :: laenge
      INTEGER,           INTENT(IN)    :: lp
      REAL   ,           INTENT(OUT)   :: ww
!                                                                       
      INTEGER il 
!                                                                       
      INTEGER len_str 
!                                                                       
!------ so far no special functions :-)                                 
!------ so we should not end up here ...                                
!                                                                       
      ier_num = - 3 
      ier_typ = ER_FORT 
!                                                                       
      il = max (len (ier_msg (1) ) - 11, len_str (string) ) 
      ier_msg (1) = 'Function : '//string (1:30) 
!                                                                       
      END SUBROUTINE kuplot_calc_intr_spec                 
!*****7**************************************************************** 
      SUBROUTINE kuplot_validate_var_spec (zeile, lp) 
!-                                                                      
!       checks whether the variable name is legal, KUPLOT specific part 
!                                                                       
!       Version : 1.0                                                   
!       Date    : 22 Oct 03                                             
!                                                                       
!       Author  : R.B. Neder  (reinhard.neder@krist.uni-erlangen.de)    
!+                                                                      
      USE errlist_mod 
      USE kuplot_config 
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER (LEN=*), INTENT(IN) :: zeile 
      INTEGER,           INTENT(IN) :: lp 
!                                                                       
      INTEGER reserved_n 
      PARAMETER (reserved_n = 22) 
                                                                        
      CHARACTER(12) reserved (reserved_n) 
      INTEGER i 
!                                                                       
      LOGICAL str_comp 
      DATA reserved / 'x', 'y', 'z', 'dx', 'dy', 'nx', 'ny', 'ni', 'np',&
      'xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zman', 'n', 'p', 'axis', &
      'pwin', 'cmap', 'cmax', 'size' /                                  
!                                                                       
      ier_num = 0 
      ier_typ = ER_NONE 
!                                                                       
      DO i = 1, reserved_n 
      IF (index (reserved (i), zeile (1:lp) ) .ne.0) then 
         ier_num = - 25 
         ier_typ = ER_FORT 
      ENDIF 
      ENDDO 
!                                                                       
      END SUBROUTINE kuplot_validate_var_spec              
