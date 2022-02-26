module kuplot_color_mod
!
!*****7*****************************************************************
!     Here are all the routines for color handling for the              
!     KUPLOT bitmaps.                                                   
!*****7*****************************************************************
!
contains
!
SUBROUTINE set_color (zeile, lp) 
!+                                                                      
!     Set colours for background & pens                                 
!-                                                                      
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
      USE kuplot_config 
      USE kuplot_mod 
use kuplot_show_mod
USE precision_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 5) 
!                                                                       
CHARACTER(len=*), intent(inout) :: zeile 
INTEGER, intent(inout) :: lp 
!                                                                       
      CHARACTER(LEN=PREC_STRING) cpara (maxw) 
      INTEGER lpara (maxw) 
      INTEGER ianz, icol 
      REAL(KIND=PREC_DP) werte (maxw) 
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
      USE get_params_mod
      USE kuplot_config 
      USE kuplot_mod
! 
      USE build_name_mod
USE precision_mod
USE str_comp_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 5) 
!                                                                       
CHARACTER(len=*), intent(inout) :: zeile 
INTEGER, intent(inout) :: lp 
!                                                                       
      CHARACTER(LEN=PREC_STRING) cpara (maxw) 
      REAL(KIND=PREC_DP) werte (maxw) 
REAL :: zzmin, zzmax
      INTEGER lpara (maxw) 
      INTEGER ianz
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
!                                                                       
      IF (ianz.ge.1) then 
         IF (str_comp (cpara (1) , 'gray', 3, lpara (1) , 4) ) then 
            CALL cmap_gray (.true.) 
            col_map_type(iwin,iframe) = COL_MAP_GRAY
         ELSEIF (str_comp (cpara (1) , 'ice', 3, lpara (1) , 3) ) then 
            CALL cmap_ice (.true.) 
            col_map_type(iwin,iframe) = COL_MAP_ICE
         ELSEIF (str_comp (cpara (1) , 'fire', 3, lpara (1) , 4) ) then 
            CALL cmap_fire (.true.) 
            col_map_type(iwin,iframe) = COL_MAP_FIRE
         ELSEIF (str_comp (cpara (1) , 'kupl', 3, lpara (1) , 4) ) then 
            CALL cmap_kupl (.true.) 
            col_map_type(iwin,iframe) = COL_MAP_KUPL
         ELSEIF (str_comp (cpara (1) , 'thermal', 3, lpara (1) , 7) ) then 
            zzmin = z_min (iwin, iframe, 1) 
            zzmax = nz(iwin, iframe, 1) * z_inc(iwin, iframe, 1) + z_min(iwin, iframe, 1)
            CALL cmap_thermal (zzmin, zzmax, .true.) 
            col_map_type(iwin,iframe) = COL_MAP_THER
         ELSEIF (str_comp (cpara (1) , 'pdf', 3, lpara (1) , 3) ) then 
            zzmin = z_min (iwin, iframe, 1) 
            zzmax = nz(iwin, iframe, 1) * z_inc(iwin, iframe, 1) + z_min(iwin, iframe, 1)
            CALL cmap_pdf (zzmin, zzmax, .true.) 
            col_map_type(iwin,iframe) = COL_MAP_PDF
         ELSEIF (str_comp (cpara (1) , 'invert', 3, lpara (1) , 6) )    &
         then                                                           
            CALL cmap_invert (.true.) 
            col_map_type(iwin,iframe) = -col_map_type(iwin,iframe)
         ELSEIF (str_comp (cpara (1) , 'read', 3, lpara (1) , 4) ) then 
            CALL do_build_name (ianz, cpara, lpara, werte, maxw, 2) 
            IF (ier_num.eq.0) call cmap_read (cpara (2) ) 
            col_map_type(iwin,iframe) = COL_MAP_READ
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
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      LOGICAL, intent(in) :: lout 
!
      REAL dummy (maxcol, 3) 
      INTEGER i 
!                                                                       
      IF (lout) WRITE(output_io, 1000) 
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
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
LOGICAL, intent(in) :: lout 
!
      INTEGER i 
!                                                                       
      IF (lout) WRITE (output_io, 1000) 
!                                                                       
      DO i = 1, maxcol 
      col_map (iwin, i, 1) = REAL(maxcol - i + 1) / REAL(maxcol) 
      col_map (iwin, i, 2) = REAL(maxcol - i + 1) / REAL(maxcol) 
      col_map (iwin, i, 3) = REAL(maxcol - i + 1) / REAL(maxcol) 
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
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
LOGICAL, intent(in) :: lout 
!
      REAL rh, rf, rq, rp, rt 
      INTEGER i, ifarb 
!                                                                       
      IF (lout) WRITE (output_io, 1000) 
!                                                                       
      DO ifarb = 1, maxcol 
      rh = 0.1 + REAL(ifarb - 1) / 284.0 
      rh = 6.0 * rh 
      i = int (rh) 
      rf = rh - REAL(i) 
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
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
LOGICAL, intent(in) :: lout 
!
      INTEGER i, ii, m 
!                                                                       
      IF (lout) WRITE (output_io, 1000) 
!                                                                       
      m = maxcol / 3 
      ii = 1 
!                                                                       
      DO i = 1, m 
      col_map (iwin, ii, 1) = REAL(3 * i) / REAL(maxcol) 
      col_map (iwin, ii, 2) = 0.0 
      col_map (iwin, ii, 3) = 0.0 
      ii = ii + 1 
      ENDDO 
!                                                                       
      DO i = 1, m 
      col_map (iwin, ii, 1) = 1.0 
      col_map (iwin, ii, 2) = REAL(3 * i) / REAL(maxcol) 
      col_map (iwin, ii, 3) = 0.0 
      ii = ii + 1 
      ENDDO 
!                                                                       
      DO i = 1, m 
      col_map (iwin, ii, 1) = 1.0 
      col_map (iwin, ii, 2) = 1.0 
      col_map (iwin, ii, 3) = REAL(3 * i) / REAL(maxcol) 
      ii = ii + 1 
      ENDDO 
!                                                                       
 1000 FORMAT     (' ------ > Setting colour map to : fire ') 
      END SUBROUTINE cmap_fire                      
!
!*****7*****************************************************************
!
      SUBROUTINE cmap_ice (lout) 
!                                                                       
!     Setting colormap to ice
!                                                                       
      USE prompt_mod 
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
LOGICAL, intent(in) :: lout 
!
      INTEGER i, ii, m 
!     LOGICAL lout 
!                                                                       
      IF (lout) WRITE (output_io, 1000) 
!                                                                       
      m = maxcol / 3 
      ii = 1 
!                                                                       
!     DO i = 1, m 
!     col_map (iwin, ii, 1) = 0.0 
!     col_map (iwin, ii, 2) = 0.0 
!     col_map (iwin, ii, 3) = REAL(3 * i) / REAL(maxcol)  /3.
!     ii = ii + 1 
!     ENDDO 
!                                                                       
!     DO i = 1, m 
!     col_map (iwin, ii, 1) = 0.0 
!     col_map (iwin, ii, 2) = REAL(3 * i) / REAL(maxcol)    /6.
!     col_map (iwin, ii, 3) = REAL(3 * i) / REAL(maxcol) *2./3.
!     ii = ii + 1 
!     ENDDO 
!                                                                       
!     DO i = 1, m 
!     col_map (iwin, ii, 1) = REAL(3 * i) / REAL(maxcol) *1./6.
!     col_map (iwin, ii, 2) = REAL(3 * i) / REAL(maxcol) *1./3.
!     col_map (iwin, ii, 3) = REAL(3 * i) / REAL(maxcol) 
!     ii = ii + 1 
!     ENDDO 
DO i=1, maxcol
   col_map (iwin, i, 3) = REAL(i)/REAL(maxcol)
ENDDO
DO i=1, maxcol
   col_map (iwin, i, 2) = REAL(i)/REAL(maxcol) * 0.800000
ENDDO
DO i=1, maxcol
   col_map (iwin, i, 1) = REAL(i)/REAL(maxcol) * 0.250000
ENDDO
!                                                                       
 1000 FORMAT     (' ------ > Setting colour map to : ice ') 
      END SUBROUTINE cmap_ice                       
!*****7*****************************************************************
!
SUBROUTINE cmap_thermal(zzmin, zzmax, lout)
!
USE prompt_mod 
USE kuplot_config 
USE kuplot_mod 
!
IMPLICIT NONE
!
REAl, INTENT(IN) :: zzmin
REAl, INTENT(IN) :: zzmax
LOGICAL, INTENT(IN) :: lout
!
INTEGER :: i, m
INTEGER :: ii
INTEGER :: istart
REAL    :: scalef
REAL    :: red, green!, blue
!
IF(zzmax > 0.0 .AND. zzmin > 0.0 ) THEN  ! positive only
   CALL cmap_fire(lout)
ELSEIF(zzmax < 0.0 .AND. zzmin < 0.0 ) THEN  ! negative only
   CALL cmap_ice(lout)
ELSE
!
!  Calculate istart , keep in window [1:maxcol]
   istart = MIN(maxcol,MAX(1,NINT(maxcol * (1. - zzmax/(zzmax-zzmin)))))
!                                                                       
   scalef = 1./(MAX(ABS(zzmax),ABS(zzmin))/(zzmax-zzmin))
   m  = INT(maxcol / 3 / scalef)
   col_map(iwin,istart,1:3) = 0.0   ! Zero level at black
   ii = istart+1 
   positive: IF(ii<=maxcol) THEN 
!                                                                       
      DO i = 1, m 
         col_map (iwin, ii, 1) = MIN(1.0,REAL(3 * i) / REAL(maxcol) *scalef)
         col_map (iwin, ii, 2) = 0.0 
         col_map (iwin, ii, 3) = 0.0 
         ii = ii + 1 
         IF(ii>maxcol) EXIT positive
      ENDDO 
!                                                                       
      red = col_map (iwin, ii-1,1)
      DO i = 1, m 
         col_map (iwin, ii, 1) = red 
         col_map (iwin, ii, 2) = MIN(1.,REAL(3 * i) / REAL(maxcol) *scalef)
         col_map (iwin, ii, 3) = 0.0 
         ii = ii + 1 
         IF(ii>maxcol) EXIT positive
      ENDDO 
!                                                                       
      green = col_map (iwin, ii-1,2)
      DO i = 1, maxcol - ii + 1
         col_map (iwin, ii, 1) = red 
         col_map (iwin, ii, 2) = green
         col_map (iwin, ii, 3) = MIN(1.0, REAL(3 * i) / REAL(maxcol) *scalef)
         ii = ii + 1 
         IF(ii>maxcol) EXIT positive
      ENDDO 
   ENDIF positive
!
!  Negative colors like ICE
!
   DO i=istart-1,1,-1
      col_map (iwin, i, 3) = REAL(istart-i)/REAL(maxcol) * scalef
   ENDDO
   DO i=istart-1,1,-1
      col_map (iwin, i, 2) = REAL(istart-i)/REAL(maxcol) * 0.800000 * scalef
   ENDDO
   DO i=istart-1,1,-1
      col_map (iwin, i, 1) = REAL(istart-i)/REAL(maxcol) * 0.250000 * scalef
   ENDDO
ENDIF
!
IF(lout) WRITE(output_io, 1000)
1000 FORMAT     (' ------ > Setting colour map to : thermal ') 
!
END SUBROUTINE cmap_thermal
!
!*****7*****************************************************************
!
SUBROUTINE cmap_pdf(zzmin, zzmax, lout)
!
USE prompt_mod 
USE kuplot_config 
USE kuplot_mod 
!
IMPLICIT NONE
!
REAl, INTENT(IN) :: zzmin
REAl, INTENT(IN) :: zzmax
LOGICAL, INTENT(IN) :: lout
!
INTEGER :: i, m, n
INTEGER :: ii
INTEGER :: istart
REAL    :: scalef
!REAL    :: red, green!, blue
!
!  CALL cmap_fire(lout)
!  CALL cmap_invert(lout)
!LSEIF(zzmax < 0.0 .AND. zzmin < 0.0 ) THEN  ! negative only
!  CALL cmap_ice(lout)
!  CALL cmap_invert(lout)
!LSE
!
!  Calculate istart , keep in window [1:maxcol]
!  istart = MIN(maxcol,MAX(1,NINT(maxcol * (1. - zzmax/(zzmax-zzmin)))))
!                                                                       
!  scalef = 1./(MAX(ABS(zzmax),ABS(zzmin))/(zzmax-zzmin))
!  m  = INT(maxcol / 3 / scalef)
!  col_map(iwin,istart,1:3) = 1.0   ! Zero level at white
!                                                                       
IF(zzmax > 0.0 .AND. zzmin > 0.0 ) THEN  ! positive only
   istart = 1
ELSEIF(zzmax < 0.0 .AND. zzmin < 0.0 ) THEN  ! negative only
   istart = maxcol
ELSE
   istart = (maxcol+1)/2
ENDIF
!
   m = (maxcol-istart) / 3 
n = MOD((maxcol-istart),3)
   ii = istart
!                                                                       
   DO i = 1, m 
      col_map (iwin, ii, 1) = 1.0 
      col_map (iwin, ii, 2) = 1.0 
      col_map (iwin, ii, 3) = 1.0 - REAL(3 * i) / REAL(maxcol) 
      ii = ii + 1 
   ENDDO 
!                                                                       
   DO i = 1, m 
      col_map (iwin, ii, 1) = 1.0 
      col_map (iwin, ii, 2) = 1.0 - REAL(3 * i) / REAL(maxcol) 
      col_map (iwin, ii, 3) = 0.0 
      ii = ii + 1 
   ENDDO 
!                                                                       
final:   DO i = 1, m + n + 1
      if(ii>maxcol) EXIT final
      col_map (iwin, ii, 1) = 1.0 - REAL(3 * i) / REAL(maxcol) 
      col_map (iwin, ii, 2) = 0.0 
      col_map (iwin, ii, 3) = 0.0 
      ii = ii + 1 
   ENDDO  final
!                                                                       
!
!  Negative colors like ICE
!
scalef = 1.
   DO i=istart-1,1, -1
      col_map (iwin, i, 3) = MAX(0.0, 1.0 - REAL(istart-i)/REAL(istart) * 0.250000 * scalef)
   ENDDO
   DO i=istart-1,1, -1
      col_map (iwin, i, 2) = MAX(0.0, 1.0 - REAL(istart-i)/REAL(istart) * 0.800000 * scalef)
   ENDDO
   DO i=istart-1,1, -1
      col_map (iwin, i, 1) = MAX(0.0, 1.0 - REAL(istart-i)/REAL(istart) * 1.000000 * scalef)
   ENDDO
!ENDIF
!
IF(lout) WRITE(output_io, 1000)
1000 FORMAT     (' ------ > Setting colour map to : pdf ') 
!
END SUBROUTINE cmap_pdf
!
!*****7*****************************************************************
!
      SUBROUTINE cmap_read (filname) 
!                                                                       
!     Reading colormap from file ..                                     
!                                                                       
      USE errlist_mod 
      USE prompt_mod 
      USE kuplot_config 
      USE kuplot_mod 
USE support_mod
!                                                                       
      IMPLICIT none 
!                                                                       
CHARACTER(len=*), intent(in) :: filname 
!
INTEGER i, ir, ig, ib 
!                                                                       
      CALL oeffne (33, filname, 'old') 
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
      USE kuplot_config 
      USE kuplot_mod 
USE support_mod
!                                                                       
      IMPLICIT none 
!                                                                       
CHARACTER(len=*), intent(in) :: filname 
!
      INTEGER i, ir, ig, ib 
!                                                                       
      CALL oeffne (33, filname, 'unknown') 
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
end module kuplot_color_mod
