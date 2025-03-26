module kuplot_update_mod
!
contains
!
      SUBROUTINE kuplot_ersetz_para (ikl, iklz, string, ll, ww, maxw, ianz)
!       Replaces a substring in an expression by the value of the       
!       appropriate parameter. Modified version for KUPLOT.             
!+                                                                      
      USE blanks_mod
      USE errlist_mod 
      USE param_mod 
      USE kuplot_config 
      USE kuplot_mod 
use kuplot_extrema_mod
use kuplot_para_mod
!
      USE variable_mod
      USE lib_upd_mod
USE lib_length
USE lib_errlist_func
      USE precision_mod
      USE precision_command_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER,                    INTENT(IN   ) :: ikl
      INTEGER,                    INTENT(IN   ) :: iklz
      CHARACTER (LEN=*),          INTENT(INOUT) :: string 
      INTEGER,                    INTENT(INOUT) :: ll
      INTEGER,                    INTENT(IN   ) :: maxw
      REAL(KIND=PREC_DP)   , DIMENSION(1:maxw), INTENT(IN   ) :: ww
      INTEGER,                    INTENT(IN   ) :: ianz
!
      INTEGER mmaxw
      PARAMETER (mmaxw = 3) 
!                                                                       
      CHARACTER(LEN=PREC_STRING) :: zeile 
      INTEGER i, laenge
      INTEGER lcomm, ltyp 
      INTEGER idummy 
      INTEGER kpara (mmaxw) 
!                                                                       
!
CALL lib_ersetz_para (ikl, iklz, string, ll, ww, maxw, ianz)
IF(ier_num == 0) RETURN
CALL no_error
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
!-------- program information n[n]                                      
!                                                                       
         IF (string (ikl - 1:ikl - 1) .eq.'n') then 
            IF (ianz.eq.1) then 
               IF (1.le.kpara (1) .and.kpara (1) .le.5) then 
                  IF (kpara (1) .eq.1) then 
                     idummy = iz - 1 
                     WRITE (zeile (ikl - 1:ikl + PREC_WIDTH-2) , PREC_F_INTE) idummy 
                  ELSEIF (kpara (1) .eq.2) then 
                     WRITE (zeile (ikl - 1:ikl + PREC_WIDTH-2) , PREC_F_INTE)         &
                     maxkurvtot                                         
                  ELSEIF (kpara (1) .eq.3) then 
                     WRITE (zeile (ikl - 1:ikl + PREC_WIDTH-2) , PREC_F_INTE)         &
                     maxframe                                           
                  ELSEIF (kpara (1) .eq.4) then 
                     WRITE (zeile (ikl - 1:ikl + PREC_WIDTH-2) , PREC_F_INTE) maxan 
                  ELSEIF (kpara (1) .eq.5) then 
                     WRITE (zeile (ikl - 1:ikl + PREC_WIDTH-2) , PREC_F_INTE) maxbond 
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
!-------- Fit parameter p[n] and error s[n]                             
!                                                                       
         ELSEIF (string (ikl - 1:ikl - 1) .eq.'p') then 
            IF (ianz.eq.1) then 
               IF (1.le.kpara (1) .and.kpara (1) .le.MAXPARA) then 
                  WRITE (zeile (ikl - 1:ikl + PREC_WIDTH-2) , PREC_F_REAL) p (    &
                  kpara (1) )                                           
                  zeile (ikl + PREC_MANTIS-lcomm:ikl + PREC_MANTIS-lcomm) = 'd' 
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
                  WRITE (zeile (ikl - 1:ikl + PREC_WIDTH-2) , PREC_F_REAL) dp (   &
                  kpara (1) )                                           
                  zeile (ikl + PREC_MANTIS-lcomm:ikl + PREC_MANTIS-lcomm) = 'd' 
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
               .and.1.le.kpara (2) .and.kpara (2) .le.lenc(kpara (1) ) )&
               then                                                     
                  WRITE (zeile (ikl - 1:ikl + PREC_WIDTH-2) , PREC_F_REAL) x (    &
                  offxy (kpara (1) - 1) + kpara (2) )                   
                  zeile (ikl + PREC_MANTIS-lcomm:ikl + PREC_MANTIS-lcomm) = 'd' 
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
               .and.1.le.kpara (2) .and.kpara (2) .le.lenc(kpara (1) ) )&
               then                                                     
                  WRITE (zeile (ikl - 1:ikl + PREC_WIDTH-2) , PREC_F_REAL) y (    &
                  offxy (kpara (1) - 1) + kpara (2) )                   
                  zeile (ikl + PREC_MANTIS-lcomm:ikl + PREC_MANTIS-lcomm) = 'd' 
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
                  WRITE (zeile (ikl - 1:ikl + PREC_WIDTH-2) , PREC_F_REAL) z (    &
                  offz (kpara (1) - 1) + (kpara (2) - 1) * ny (kpara (1)&
                  ) + kpara (3) )                                       
                  zeile (ikl + PREC_MANTIS-lcomm:ikl + PREC_MANTIS-lcomm) = 'd' 
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
               .and.1.le.kpara (2) .and.kpara (2) .le.lenc(kpara (1) ) )&
               then                                                     
                  WRITE (zeile (ikl - 2:ikl + PREC_WIDTH-2) , PREC_F_REAL) dx (   &
                  offxy (kpara (1) - 1) + kpara (2) )                   
                  zeile (ikl + PREC_MANTIS-lcomm:ikl + PREC_MANTIS-lcomm) = 'd' 
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
               .and.1.le.kpara (2) .and.kpara (2) .le.lenc(kpara (1) ) )&
               then                                                     
                  WRITE (zeile (ikl - 2:ikl + PREC_WIDTH-2) , PREC_F_REAL) dy (   &
                  offxy (kpara (1) - 1) + kpara (2) )                   
                  zeile (ikl + PREC_MANTIS-lcomm:ikl + PREC_MANTIS-lcomm) = 'd' 
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
                     WRITE(zeile(ikl-2:ikl+PREC_WIDTH-2),PREC_F_INTE) nx(kpara(1))
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
                     WRITE(zeile(ikl-2:ikl+PREC_WIDTH-2),PREC_F_INTE) ny(kpara(1))
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
                     WRITE(zeile(ikl-2:ikl+PREC_WIDTH-2),PREC_F_INTE) ny(kpara(1))*ny(kpara(1))
                  ELSE
                     WRITE(zeile(ikl-2:ikl+PREC_WIDTH-2),PREC_F_INTE) lenc(kpara(1))
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
                     WRITE (zeile (ikl - 2:ikl + PREC_WIDTH-2) , PREC_F_INTE) 1 
                  ELSE 
                     WRITE (zeile (ikl - 2:ikl + PREC_WIDTH-2) , PREC_F_INTE) 0 
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
      ELSEIF (lcomm.eq.4) then 
         IF (ikl.gt.5) zeile (1:ikl - 4) = string (1:ikl - 5) 
!                                                                       
!-------- data set dimensions xmin[n], xmax[n], ..                      
!                                                                       
         IF (string (ikl - 4:ikl - 1) .eq.'xmin') then 
            IF (ianz.eq.1) then 
               IF (1.le.kpara (1) .and.kpara (1) .le. (iz - 1) ) then 
                  CALL get_extrema 
                  WRITE (zeile (ikl - 4:ikl + PREC_WIDTH-2) , PREC_F_REAL) xmin ( &
                  kpara (1) )                                           
                  zeile (ikl + PREC_MANTIS-lcomm:ikl + PREC_MANTIS-lcomm) = 'd' 
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
                  WRITE (zeile (ikl - 4:ikl + PREC_WIDTH-2) , PREC_F_REAL) xmax ( &
                  kpara (1) )                                           
                  zeile (ikl + PREC_MANTIS-lcomm:ikl + PREC_MANTIS-lcomm) = 'd' 
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
                  WRITE (zeile (ikl - 4:ikl + PREC_WIDTH-2) , PREC_F_REAL) ymin ( &
                  kpara (1) )                                           
                  zeile (ikl + PREC_MANTIS-lcomm:ikl + PREC_MANTIS-lcomm) = 'd' 
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
                  WRITE (zeile (ikl - 4:ikl + PREC_WIDTH-2) , PREC_F_REAL) ymax ( &
                  kpara (1) )                                           
                  zeile (ikl + PREC_MANTIS-lcomm:ikl + PREC_MANTIS-lcomm) = 'd' 
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
                  WRITE (zeile (ikl - 4:ikl + PREC_WIDTH-2) , PREC_F_REAL) zmin ( &
                  kpara (1) )                                           
                  zeile (ikl + PREC_MANTIS-lcomm:ikl + PREC_MANTIS-lcomm) = 'd' 
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
                  WRITE (zeile (ikl - 4:ikl + PREC_WIDTH-2) , PREC_F_REAL) zmax ( &
                  kpara (1) )                                           
                  zeile (ikl + PREC_MANTIS-lcomm:ikl + PREC_MANTIS-lcomm) = 'd' 
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
                  WRITE (zeile (ikl - 4:ikl + PREC_WIDTH-2) , PREC_F_REAL) dev_sf &
                  (iwin, x11)                                           
                  zeile (ikl + PREC_MANTIS-lcomm:ikl + PREC_MANTIS-lcomm) = 'd' 
               ELSEIF (2.eq.kpara (1) ) then 
                  WRITE (zeile (ikl - 4:ikl + PREC_WIDTH-2) , PREC_F_REAL) dev_sf &
                  (iwin, ps)                                            
                  zeile (ikl + PREC_MANTIS-lcomm:ikl + PREC_MANTIS-lcomm) = 'd' 
               ELSEIF (3.eq.kpara (1) ) then 
                  WRITE (zeile (ikl - 4:ikl + PREC_WIDTH-2) , PREC_F_REAL) dev_sf &
                  (iwin, pic)                                           
                  zeile (ikl + PREC_MANTIS-lcomm:ikl + PREC_MANTIS-lcomm) = 'd' 
               ELSEIF (4.eq.kpara (1) ) then 
                  WRITE (zeile (ikl - 4:ikl + PREC_WIDTH-2) , '(e15.8e2)') dev_sf &
                  (iwin, png)                                           
                  zeile (ikl + PREC_MANTIS-lcomm:ikl + PREC_MANTIS-lcomm) = 'd' 
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
                  WRITE (zeile (ikl - 4:ikl + PREC_WIDTH-2) , PREC_F_REAL) pex (  &
                  iwin, iframe, 1)                                      
                  zeile (ikl + PREC_MANTIS-lcomm:ikl + PREC_MANTIS-lcomm) = 'd' 
               ELSEIF (2.eq.kpara (1) ) then 
                  WRITE (zeile (ikl - 4:ikl + PREC_WIDTH-2) , PREC_F_REAL) pex (  &
                  iwin, iframe, 2)                                      
                  zeile (ikl + PREC_MANTIS-lcomm:ikl + PREC_MANTIS-lcomm) = 'd' 
               ELSEIF (3.eq.kpara (1) ) then 
                  WRITE (zeile (ikl - 4:ikl + PREC_WIDTH-2) , PREC_F_REAL) pey (  &
                  iwin, iframe, 1)                                      
                  zeile (ikl + PREC_MANTIS-lcomm:ikl + PREC_MANTIS-lcomm) = 'd' 
               ELSEIF (4.eq.kpara (1) ) then 
                  WRITE (zeile (ikl - 4:ikl + PREC_WIDTH-2) , PREC_F_REAL) pey (  &
                  iwin, iframe, 2)                                      
                  zeile (ikl + PREC_MANTIS-lcomm:ikl + PREC_MANTIS-lcomm) = 'd' 
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
                  WRITE (zeile (ikl - 4:ikl + PREC_WIDTH-2) , PREC_F_REAL)        &
                  col_map (iwin, kpara (1) , kpara (2) )                
                  zeile (ikl + PREC_MANTIS-lcomm:ikl + PREC_MANTIS-lcomm) = 'd' 
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
                  WRITE (zeile (ikl - 4:ikl + PREC_WIDTH-2) , PREC_F_INTE) maxcol 
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
                     WRITE (zeile (ikl - 4:ikl + PREC_WIDTH-2) , PREC_F_REAL)     &
                     lab_angle (iwin, iframe, kpara (2) )               
                     zeile (ikl + PREC_MANTIS-lcomm:ikl + PREC_MANTIS-lcomm) = 'd' 
                  ELSEIF (kpara (1) .eq.2) then 
                     WRITE (zeile (ikl - 4:ikl + PREC_WIDTH-2) , PREC_F_REAL)     &
                     tick_ma_h (iwin, iframe, kpara (2) )               
                     zeile (ikl + PREC_MANTIS-lcomm:ikl + PREC_MANTIS-lcomm) = 'd' 
                  ELSEIF (kpara (1) .eq.3) then 
                     WRITE (zeile (ikl - 4:ikl + PREC_WIDTH-2) , PREC_F_REAL)     &
                     tick_mi_h (iwin, iframe, kpara (2) )               
                     zeile (ikl + PREC_MANTIS-lcomm:ikl + PREC_MANTIS-lcomm) = 'd' 
                  ELSEIF (kpara (1) .eq.4) then 
                     WRITE (zeile (ikl - 4:ikl + PREC_WIDTH-2) , PREC_F_INTE)         &
                     tick_nsub (iwin, iframe, kpara (2) )               
                  ELSEIF (kpara (1) .eq.5) then 
                     WRITE (zeile (ikl - 4:ikl + PREC_WIDTH-2) , PREC_F_REAL)     &
                     lab_d (iwin, iframe, kpara (2) )                   
                     zeile (ikl + PREC_MANTIS-lcomm:ikl + PREC_MANTIS-lcomm) = 'd' 
                  ELSEIF (kpara (1) .eq.6) then 
                     WRITE (zeile (ikl - 4:ikl + PREC_WIDTH-2) , PREC_F_REAL)     &
                     ax_d (iwin, iframe, kpara (2) )                    
                     zeile (ikl + 7:ikl + 7) = 'd' 
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
      ELSE 
         ier_num = - 2 
         ier_typ = ER_FORT 
      ENDIF 
      IF (ier_num.eq.0) then 
         ll = laenge+PREC_WIDTH - ltyp - (iklz - ikl + 1) 
         IF (iklz + 1.le.laenge) zeile (ikl + PREC_WIDTH-1:ll) = string (iklz + 1:&
         laenge)                                                        
         string = zeile 
!        CALL rem_bl (string, ll) 
      ELSE
         ll = min (40, laenge)
         WRITE (ier_msg (1), '(a)') string (1:ll)
      ENDIF 
      ll = LEN_TRIM(string)
!
      END SUBROUTINE kuplot_ersetz_para                    
!
!*****7*****************************************************************
!
SUBROUTINE kuplot_upd_para (ctype, ww, maxw, wert, ianz, cstring, substr) 
!-                                                                      
!       updates the parameter spezified by ctype, index ww  to the      
!       new value of wert                                               
!+                                                                      
      USE errlist_mod 
      USE param_mod 
      USE prompt_mod 
      USE kuplot_config 
      USE kuplot_mod 
USE lib_errlist_func
      USE lib_upd_mod
USE precision_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER (LEN=*),          INTENT(IN) :: ctype 
      INTEGER,                    INTENT(IN) :: maxw
      INTEGER,                    INTENT(IN) :: ianz 
      INTEGER, DIMENSION(1:MAXW), INTENT(IN) :: ww
      REAL(KIND=PREC_DP)        , INTENT(IN) :: wert 
      CHARACTER (LEN=*),          INTENT(IN) :: cstring
INTEGER, DIMENSION(2), INTENT(IN)    :: substr ! Indices of substring
!
      INTEGER idummy 
CALL lib_upd_para (ctype, ww, maxw, wert, ianz, cstring, substr)
IF(ier_num==0 .OR. (ier_num==-40 .AND. ier_typ==ER_FORT)) RETURN
CALL no_error
!                                                                       
!------ fit parameter p[n]                                              
!                                                                       
      IF (ctype.eq.'p') then 
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
            .and.ww (2) .le.lenc (ww (1) ) ) then                        
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
            .and.ww (2) .le.lenc (ww (1) ) ) then                        
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
            .and.ww (2) .le.lenc (ww (1) ) ) then                        
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
            .and.ww (2) .le.lenc (ww (1) ) ) then                        
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
               dev_sf (iwin, x11) = real(wert , kind=PREC_SP)
            ELSEIF (2.eq.ww (1) ) then 
               dev_sf (iwin, ps) = real(wert , kind=PREC_SP)
               dev_sf (iwin, vps) = real(wert , kind=PREC_SP)
            ELSEIF (3.eq.ww (1) ) then 
               dev_sf (iwin, pic) = real(wert , kind=PREC_SP)
               dev_sf (iwin, vpic) = real(wert , kind=PREC_SP)
            ELSEIF (4.eq.ww (1) ) then 
               dev_sf (iwin, png) = real(wert , kind=PREC_SP)
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
                  col_map (iwin, ww (1), ww (2) ) = real(wert , kind=PREC_SP)
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
                  lab_angle (iwin, iframe, ww (2) ) = real(wert , kind=PREC_SP)
               ELSEIF (ww (1) .eq.2) then 
                  tick_ma_h (iwin, iframe, ww (2) ) = real(wert , kind=PREC_SP)
               ELSEIF (ww (1) .eq.3) then 
                  tick_ma_h (iwin, iframe, ww (2) ) = real(wert , kind=PREC_SP)
               ELSEIF (ww (1) .eq.4) then 
                  tick_nsub (iwin, iframe, ww (2) ) = nint (wert) 
               ELSEIF (ww (1) .eq.5) then 
                  lab_d (iwin, iframe, ww (2) ) = real(wert , kind=PREC_SP)
               ELSEIF (ww (1) .eq.6) then 
                  ax_d (iwin, iframe, ww (2) ) = real(wert , kind=PREC_SP)
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
         WRITE (ier_msg (1), '(a)') ctype
      ENDIF 
!                                                                       
 2000 FORMAT     (' ------ > Last ',i3,' data set(s) deleted ...') 
      END SUBROUTINE kuplot_upd_para                       
!*****7***************************************************************  
SUBROUTINE kuplot_calc_intr_spec(string, line, ikl, iklz, ww, laenge, lp)
!-                                                                      
!     These are special intrinsic function for the KUPLOT. Any          
!     intrinsic function that references KUPLOT specific values         
!     is found in this subroutine.                                      
!+                                                                      
USE errlist_mod 
USE lib_length
USE precision_mod
!
      IMPLICIT none 
!                                                                       
      CHARACTER (LEN=*), INTENT(INOUT) :: string
      CHARACTER (LEN=*), INTENT(INOUT) :: line 
      INTEGER,           INTENT(IN)    :: ikl
      INTEGER,           INTENT(IN)    :: iklz
      INTEGER,           INTENT(INOUT) :: laenge
      INTEGER,           INTENT(INOUT) :: lp
      REAL(KIND=PREC_DP),INTENT(INOUT) :: ww
!                                                                       
      INTEGER il 
!                                                                       
!                                                                       
!------ so far no special functions :-)                                 
!------ so we should not end up here ...                                
!                                                                       
      ier_num = - 3 
      ier_typ = ER_FORT 
!                                                                       
      il = max (LEN(ier_msg (1) ) - 11, len_str(string) ) 
      ier_msg (1) = 'Function : '//string (1:30) 
!                                                                       
      END SUBROUTINE kuplot_calc_intr_spec                 
!
!*****7**************************************************************** 
!
SUBROUTINE kuplot_calc_intr_log_spec(string, length)
!
IMPLICIT NONE
CHARACTER(LEN=*) , INTENT(INOUT) :: string
INTEGER          , INTENT(INOUT) :: length
!
END SUBROUTINE kuplot_calc_intr_log_spec
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
USE reserved_mod
USE errlist_mod 
USE kuplot_config 
!                                                                       
IMPLICIT none 
!                                                                       
CHARACTER (LEN=*), INTENT(IN) :: zeile 
INTEGER,           INTENT(IN) :: lp 
!                                                                       
INTEGER :: i , length
!                                                                       
ier_num = 0 
ier_typ = ER_NONE 
!                                                                       
main: DO i = 1, kuplot_reserved_n 
!  IF (index (kuplot_reserved (i), zeile (1:lp) ) .ne.0) THEN 
   length = MAX(LEN_TRIM(kuplot_reserved(i)), LEN_TRIM(zeile(1:lp)))
   length = MIN(length, LEN(kuplot_reserved), LEN(zeile))
   IF(kuplot_reserved (i)(1:length)== zeile(1:length) ) THEN           
      ier_num = - 25 
      ier_typ = ER_FORT 
      EXIT main
   ENDIF 
ENDDO  main
!                                                                       
END SUBROUTINE kuplot_validate_var_spec              
!
!*******************************************************************************
!
SUBROUTINE kuplot_get_var_type(line,length, var_is_type)
!
! Returns the variable type : INTEGER, REAL, CHARACTER, and Scalar versus field
!
USE constants_mod
USE lib_get_var
USE variable_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*)     , INTENT(IN)  :: line
INTEGER              , INTENT(IN)  :: length
INTEGER, DIMENSION(3), INTENT(OUT) :: var_is_type
!
INTEGER, PARAMETER :: MAXPAR = 23
CHARACTER(LEN=16), DIMENSION(MAXPAR) :: kuplot_names
INTEGER          , DIMENSION(MAXPAR) :: kuplot_type
INTEGER          , DIMENSION(MAXPAR) :: kuplot_dim
LOGICAL          , DIMENSION(MAXPAR) :: kuplot_ro 
INTEGER :: i
!
DATA kuplot_names  &
    /'zmax    ', 'zmin    ', 'ymax    ', 'ymin    ', 'xmax    ', &
     'xmin    ', 'size    ', 'pwin    ', 'cmax    ', 'cmap    ', &
     'axis    ', 'ny      ', 'nx      ', 'np      ', 'ni      ', &
     'dy      ', 'dx      ', 'z       ', 'y       ', 'x       ', &
     's       ', 'p       ', 'n       '                          &
    /
DATA kuplot_type &
    /  IS_REAL ,   IS_REAL ,   IS_REAL ,   IS_REAL ,   IS_REAL , &
       IS_REAL ,   IS_REAL ,   IS_REAL ,   IS_REAL ,   IS_REAL , &
       IS_REAL ,   IS_INTE ,   IS_INTE ,   IS_INTE ,   IS_INTE , &
       IS_REAL ,   IS_REAL ,   IS_REAL ,   IS_REAL ,   IS_REAL , &
       IS_REAL ,   IS_REAL ,   IS_INTE                           &
    /
DATA kuplot_dim  &
    /  IS_VEC  ,   IS_VEC  ,   IS_VEC  ,   IS_VEC  ,   IS_ARR  , &
       IS_VEC  ,   IS_VEC  ,   IS_VEC  ,   IS_VEC  ,   IS_VEC  , &
       IS_ARR  ,   IS_VEC  ,   IS_VEC  ,   IS_ARR  ,   IS_VEC  , &
       IS_ARR  ,   IS_ARR  ,   IS_ARR  ,   IS_ARR  ,   IS_ARR  , &
       IS_VEC  ,   IS_VEC  ,   IS_VEC                            &
    /
DATA kuplot_ro  &
    /  .TRUE.  ,   .TRUE.  ,   .TRUE.  ,   .TRUE.  ,   .TRUE.  , &
       .TRUE.  ,   .FALSE. ,   .TRUE.  ,   .FALSE. ,   .FALSE. , &
       .FALSE. ,   .TRUE.  ,   .TRUE.  ,   .TRUE.  ,   .TRUE.  , &
       .FALSE. ,   .FALSE. ,   .FALSE. ,   .FALSE. ,   .FALSE. , &
       .TRUE.  ,   .FALSE. ,   .FALSE.                           &
    /
!
var_is_type(:) = IS_UNKNOWN
!
main: DO i=1, MAXPAR
   IF(line(1:length) == kuplot_names(i)(1:LEN_TRIM(kuplot_names(i)))) THEN
      var_is_type(1) = kuplot_type(i)
      var_is_type(2) = kuplot_dim (i)
      IF(kuplot_ro(i)) THEN
         var_is_type(3) = IS_READ
      ELSE
         var_is_type(3) = IS_WRITE
      ENDIF
      RETURN
   ENDIF
ENDDO main
!
CALL lib_get_var_type(line, length, var_is_type)
!
!
END SUBROUTINE kuplot_get_var_typE
!
end module kuplot_update_mod
