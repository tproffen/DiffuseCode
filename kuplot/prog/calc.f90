module kuplot_calc_mod
!
contains
!
SUBROUTINE do_calc (zeile, lp) 
!*****7*****************************************************************
!                                                                       
!     These routines perform data manipulation - which is much          
!     faster than using the build in interpreter - but less             
!     flexible ..                                                       
!                                                                       
!*****7*****************************************************************
!                                                                       
USE ber_params_mod
USE errlist_mod 
USE get_params_mod
USE kuplot_config 
USE kuplot_mod 
use kuplot_extrema_mod
use kuplot_show_mod
USE string_convert_mod
USE precision_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER , PARAMETER :: maxw = 4
!                                                                       
CHARACTER(len=*), intent(inout) :: zeile
integer         , intent(inout) :: lp
! 
CHARACTER(LEN=PREC_STRING), dimension(maxw) :: cpara ! (maxw) 
CHARACTER(len=3)                            :: oper 
CHARACTER(len=2)                            :: unt 
REAL(KIND=PREC_DP)        , dimension(MAXW) :: werte ! (maxw) 
INTEGER                   , dimension(MAXW) :: lpara ! (maxw)
INTEGER :: ianz, ilen, ik 
!                                                                       
CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
IF (ier_num.ne.0) return 
!                                                                       
IF (ianz.lt.3.or.ianz.gt.4) then 
   ier_num = - 6 
   ier_typ = ER_COMM 
   RETURN 
ENDIF 
!                                                                       
oper = cpara (1) (1:3) 
unt = cpara (2) (1:2) 
CALL do_cap (oper) 
CALL do_cap (unt) 
!                                                                       
CALL del_params (2, ianz, cpara, lpara, maxw) 
CALL ber_params (ianz, cpara, lpara, werte, maxw) 
IF (ier_num.ne.0) return 
!                                                                       
ik = nint (werte (1) ) 
IF (ik.gt.0.and.ik.lt.iz) then 
   ilen = lenc (ik) 
   IF (unt.eq.'WX') then 
      IF (lni (ik) ) ilen = nx (ik) 
      if(oper=='SIG') then
         CALL calc_xy_d (x,dx, ilen, ik) 
      else
         CALL calc_xy (x, ilen, ik, oper, werte, maxw, ianz) 
      endif
      IF (ier_num.ne.0) return 
      CALL get_extrema_xy (x, ik, ilen, xmin, xmax) 
      CALL show_data (ik) 
   ELSEIF (unt.eq.'WY') then 
      IF (lni (ik) ) ilen = ny (ik) 
      if(oper=='SIG') then
         CALL calc_xy_d (y,dy, ilen, ik)
      else
         CALL calc_xy (y, ilen, ik, oper, werte, maxw, ianz) 
      endif
      IF (ier_num.ne.0) return 
      CALL get_extrema_xy (y, ik, ilen, ymin, ymax) 
      CALL show_data (ik) 
   ELSEIF (unt.eq.'DX') then 
      CALL calc_xy (dx, lenc (ik), ik, oper, werte, maxw, ianz) 
   ELSEIF (unt.eq.'DY') then 
      CALL calc_xy (dy, lenc (ik), ik, oper, werte, maxw, ianz) 
   ELSEIF (unt.eq.'WZ') then 
      IF (lni (ik) ) then 
         CALL calc_z(z, nx(ik), ny(ik), ik, oper, werte, maxw, ianz)
         IF (ier_num.ne.0) return 
         CALL get_extrema_z (z, ik, nx (ik), ny (ik), zmin, zmax) 
         CALL show_data (ik) 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
   ELSE 
      ier_num = - 6 
      ier_typ = ER_COMM 
   ENDIF 
ELSE 
   ier_num = - 4 
   ier_typ = ER_APPL 
ENDIF 
!                                                                  
END SUBROUTINE do_calc                        
!
!*****7*****************************************************************
!
SUBROUTINE calc_xy (a, ilen, ik, op, werte, maxw, ianz) 
!+                                                                      
!     berechnungen fuer x und y feld                                    
!-                                                                      
USE errlist_mod 
USE kuplot_config 
USE kuplot_mod 
USE precision_mod
!                                                                       
IMPLICIT none 
!                                                                       
REAL(kind=PREC_DP), dimension(MAXARRAY), intent(inout) :: a ! (maxarray) 
integer                            , intent(in)    :: ilen
integer                            , intent(in)    :: ik
CHARACTER(len=*)                   , intent(in)    :: op 
integer                            , intent(in)    :: MAXW
integer                            , intent(in)    :: ianz
REAL(KIND=PREC_DP), DIMENSION(MAXW), intent(in)    :: werte (maxw)
      INTEGER :: i 
REAL(kind=PREC_DP) :: summand, faktor , thresh
!                                                                       
IF(op.eq.'INV') then 
   DO i = 1, ilen 
      IF (a (offxy (ik - 1) + i) .ne.0.0) then 
         a (offxy (ik - 1) + i) = 1.0 / a (offxy (ik - 1) + i) 
      ELSE 
         a (offxy (ik - 1) + i) = 0.0 
      ENDIF 
   ENDDO 
ELSEIF (op.eq.'LOG') then 
   DO i = 1, ilen 
      IF (a (offxy (ik - 1) + i) .gt.0.0) then 
         a (offxy (ik - 1) + i) = log (a (offxy (ik - 1) + i) ) 
      ELSE 
         a (offxy (ik - 1) + i) = 0.0 
      ENDIF 
   ENDDO 
ELSEIF (op.eq.'EXP') then 
   DO i = 1, ilen 
      a (offxy (ik - 1) + i) = exp (a (offxy (ik - 1) + i) ) 
   ENDDO 
ELSEIF (op.eq.'SQU') then 
   DO i = 1, ilen 
      a(offxy(ik - 1) + i) = a(offxy(ik - 1) + i) * a(offxy(ik - 1) + i)
   ENDDO 
ELSEIF (op.eq.'SQR') then 
   DO i = 1, ilen 
      IF (a (offxy (ik - 1) + i) .ge.0.0) then 
         a (offxy (ik - 1) + i) = sqrt (a (offxy (ik - 1) + i) ) 
      ELSE 
         a (offxy (ik - 1) + i) = 0.0 
      ENDIF 
   ENDDO 
ELSEIF (op.eq.'THR') then 
   thresh = -9999.D0
   IF (ianz.eq.2) thresh = werte (2) 
   DO i = 1, ilen 
      IF( a(offxy (ik - 1) + i)<= thresh .AND. a(offxy (ik - 1) + i)/=-9999.) THEN
         a (offxy (ik - 1) + i) = -9999.0
      ENDIF
   ENDDO 
ELSEIF (op.eq.'ADD') then 
   summand = 0.0D0 
   IF (ianz.eq.2) summand = werte (2) 
   DO i = 1, ilen 
      a (offxy (ik - 1) + i) = a (offxy (ik - 1) + i) + summand 
   ENDDO 
ELSEIF (op.eq.'MUL') then 
   faktor = 1.0D0 
   IF (ianz.eq.2) faktor = werte (2) 
   DO i = 1, ilen 
      a (offxy (ik - 1) + i) = a (offxy (ik - 1) + i) * faktor 
   ENDDO 
ELSEIF (op.eq.'ABS') then 
   DO i = 1, ilen 
      a (offxy (ik - 1) + i) = abs (a (offxy (ik - 1) + i) ) 
      ENDDO 
ELSE 
   ier_num = - 6 
   ier_typ = ER_COMM 
ENDIF 
!                                                                       
END SUBROUTINE calc_xy                        
!
!*****7*****************************************************************
!
SUBROUTINE calc_xy_d (a, da, ilen, ik) 
!+                                                                      
!     berechnungen fuer x und y feld                                    
!-                                                                      
USE errlist_mod 
USE kuplot_config 
USE kuplot_mod 
use lib_random_func
USE precision_mod
!                                                                       
IMPLICIT none 
!                                                                       
REAL(kind=PREC_DP), dimension(MAXARRAY), intent(inout) :: a ! (maxarray) 
REAL(kind=PREC_DP), dimension(MAXARRAY), intent(in   ) :: da! (maxarray) 
integer                                , intent(in)    :: ilen
integer                                , intent(in)    :: ik
INTEGER :: i 
!                                                                       
do i = 1, ilen
  a(offxy(ik - 1) + i) = a(offxy(ik - 1) + i) + gasdev(da(offxy(ik - 1) + i))
enddo
!                                                                       
END SUBROUTINE calc_xy_d
!
!*****7*****************************************************************
!
SUBROUTINE calc_z (a, nxx, nyy, ik, op, werte, maxw, ianz) 
!+                                                                      
!     berechnungen fuer z feld                                          
!-                                                                      
USE errlist_mod 
USE kuplot_config 
USE kuplot_mod 
USE precision_mod
!                                                                       
IMPLICIT none 
!                                                                       
REAL(kind=PREC_DP), dimension(MAXARRAY), intent(inout) :: a ! (maxarray) 
integer                                , intent(in)    :: nxx
integer                                , intent(in)    :: nyy
integer                                , intent(in)    :: ik
CHARACTER(len=*)                       , intent(in)    :: op 
integer                                , intent(in)    :: MAXW
REAL(KIND=PREC_DP), dimension(MAXW)    , intent(in)    :: werte ! (maxw)
integer                                , intent(in)    :: ianz

INTEGER :: i, j, ikk 
REAL(kind=PREC_DP) :: faktor, summand , thresh
!                                                                       
IF (op.eq.'INV') then 
   DO i = 1, nxx 
      DO j = 1, nyy 
         ikk = offz (ik - 1) + (i - 1) * ny (ik) + j 
         IF (a (ikk) .ne.0.0.and.a (ikk) .ne. - 9999.0) then 
            a (ikk) = 1.0 / a (ikk) 
         ENDIF 
      ENDDO 
   ENDDO 
ELSEIF (op.eq.'LOG') then 
   DO i = 1, nxx 
      DO j = 1, nyy 
         ikk = offz (ik - 1) + (i - 1) * ny (ik) + j 
         IF (a (ikk) .ne.0.0.and.a (ikk) .ne. - 9999.0) then 
            a (ikk) = log (a (ikk) ) 
         ENDIF 
      ENDDO 
   ENDDO 
ELSEIF (op.eq.'EXP') then 
   DO i = 1, nxx 
      DO j = 1, nyy 
         ikk = offz (ik - 1) + (i - 1) * ny (ik) + j 
         IF (a (ikk) .ne. - 9999.0) then 
            a (ikk) = exp (a (ikk) ) 
         ENDIF 
      ENDDO 
   ENDDO 
ELSEIF (op.eq.'SQU') then 
   DO i = 1, nxx 
      DO j = 1, nyy 
         ikk = offz (ik - 1) + (i - 1) * ny (ik) + j 
         IF (a (ikk) .ne. - 9999.0) then 
            a (ikk) = a (ikk) **2 
         ENDIF 
      ENDDO 
   ENDDO 
ELSEIF (op.eq.'SQR') then 
   DO i = 1, nxx 
      DO j = 1, nyy 
         ikk = offz (ik - 1) + (i - 1) * ny (ik) + j 
         IF (a (ikk) .ge.0.0.and.a (ikk) .ne. - 9999.0) then 
            a (ikk) = sqrt (a (ikk) ) 
         ENDIF 
      ENDDO 
   ENDDO 
ELSEIF (op.eq.'THR') then 
   thresh = -9999.00
   IF (ianz.eq.2) thresh = werte (2) 
   DO i = 1, nxx 
      DO j = 1, nyy 
         ikk = offz (ik - 1) + (i - 1) * ny (ik) + j 
         IF (a (ikk) <= thresh.and.a (ikk) .ne. - 9999.0) then 
            a (ikk) =  -9999.0
         ENDIF 
      ENDDO 
   ENDDO 
ELSEIF (op.eq.'ADD') then 
   summand = 0.0 
   IF (ianz.eq.2) summand = werte (2) 
   DO i = 1, nxx 
      DO j = 1, nyy 
         ikk = offz (ik - 1) + (i - 1) * ny (ik) + j 
         IF (a (ikk) .ne. - 9999.0) then 
            a (ikk) = a (ikk) + summand 
         ENDIF 
      ENDDO 
   ENDDO 
ELSEIF (op.eq.'MUL') then 
   faktor = 1.0 
   IF (ianz.eq.2) faktor = werte (2) 
   DO i = 1, nxx 
      DO j = 1, nyy 
         ikk = offz (ik - 1) + (i - 1) * ny (ik) + j 
         IF (a (ikk) .ne. - 9999.0) then 
            a (ikk) = a (ikk) * faktor 
         ENDIF 
      ENDDO 
   ENDDO 
ELSE 
   ier_num = - 6 
   ier_typ = ER_COMM 
ENDIF 
!                                                                       
END SUBROUTINE calc_z                         
!
!*****7*****************************************************************
!
SUBROUTINE do_exclude (zeile, lp) 
!+                                                                      
!     This routine allows the user to exclude data regions              
!     and replace by linear interpolation of end points.                
!-                                                                      
USE ber_params_mod
USE errlist_mod 
USE get_params_mod
USE kuplot_config 
USE kuplot_mod 
use kuplot_extrema_mod
USE precision_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER maxw 
PARAMETER (maxw = 3) 
!                                                                       
CHARACTER(len=*), intent(inout) :: zeile 
integer         , intent(inout) :: lp
!
CHARACTER(LEN=PREC_STRING), dimension(MAXW) :: cpara ! (maxw) 
REAL(KIND=PREC_DP) :: werte (maxw) 
REAL(kind=PREC_DP) :: a, b, x1, x2, y1, y2, xleft, xright 
INTEGER, dimension(MAXW) :: lpara !(maxw)
INTEGER :: i, ipkt, ianz, ik, ileft, iright
!                                                                       
!      INTEGER closest_point 
!                                                                       
CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
IF (ier_num.ne.0) return 
!                                                                       
IF (ianz.ne.3) then 
   ier_num = - 6 
   ier_typ = ER_COMM 
   RETURN 
ENDIF 
!                                                                       
CALL ber_params (ianz, cpara, lpara, werte, maxw) 
IF (ier_num.ne.0) return 
!                                                                       
ik = nint (werte (1) ) 
xleft = werte (2) 
xright = werte (3) 
!                                                                       
IF (ik.gt.0.and.ik.lt.iz.and..not.lni (ik) ) then 
   ileft = closest_point (ik, xleft) 
   iright = closest_point (ik, xright) 
   x1 = x (offxy (ik - 1) + ileft) 
   y1 = y (offxy (ik - 1) + ileft) 
   x2 = x (offxy (ik - 1) + iright) 
   y2 = y (offxy (ik - 1) + iright) 
!                                                                       
   IF (x2.le.x1) return 
!                                                                       
!-------- region is in the middle                                       
!                                                                       
   IF (xleft.gt.xmin (ik) .and.xright.lt.xmax (ik) ) then 
      a = (y2 - y1) / (x2 - x1) 
      b = y1 - a * x1 
      DO i = 1, lenc (ik) 
         ipkt = offxy (ik - 1) + i 
         IF (x (ipkt) .gt.xleft.and.x (ipkt) .lt.xright) then 
            y (ipkt) = a * x (ipkt) + b 
         ENDIF 
      ENDDO 
!                                                                       
!-------- left is outside                                               
!                                                                       
   ELSEIF (xleft.lt.xmin (ik) .and.xright.lt.xmax (ik) ) then 
      IF(ileft<iright) THEN   ! normal dataset
         DO i = iright, lenc (ik) 
            ipkt = offxy (ik - 1) + i 
            x (ipkt - iright + 1) = x (ipkt) 
            dx(ipkt - iright + 1) = dx(ipkt) 
            y (ipkt - iright + 1) = y (ipkt) 
            dy(ipkt - iright + 1) = dy(ipkt) 
         ENDDO 
         lenc (ik) = lenc (ik) - iright + 1 
      ELSE     ! Data set in inverse sequence from + to -
         lenc(ik) = iright
      ENDIF
      CALL get_extrema_xy (x, ik, lenc(ik), xmin, xmax) 
!                                                                       
!-------- right is outside                                              
!                                                                       
   ELSEIF (xleft.gt.xmin (ik) .and.xright.gt.xmax (ik) ) then 
      IF(ileft<iright) THEN   ! normal dataset
         lenc (ik) = ileft 
      ELSE     ! Data set in inverse sequence from + to -
         DO i=1, lenc(ik) - ileft + 1
            ipkt = offxy (ik - 1) + i 
            x (ipkt) = x (ipkt + ileft -1)
            dx(ipkt) = dx(ipkt + ileft -1)
            y (ipkt) = y (ipkt + ileft -1)
            dy(ipkt) = dy(ipkt + ileft -1)
         ENDDO
         lenc (ik) = lenc (ik) - ileft + 1 
      ENDIF
      CALL get_extrema_xy (x, ik, lenc(ik), xmin, xmax) 
!                                                                       
   ELSE 
      ier_num = - 7 
      ier_typ = ER_APPL 
   ENDIF 
ELSE 
   ier_num = - 4 
   ier_typ = ER_APPL 
ENDIF 
!                                                                       
END SUBROUTINE do_exclude                     
!
!*****7*****************************************************************
!
INTEGER function closest_point (ik, xvalue) 
!+                                                                      
!     Find the point closest to xvalue;                                 
!-                                                                      
USE kuplot_config 
USE kuplot_mod 
!
use precision_mod
!                                                                       
IMPLICIT none 
!                                                                       
integer, intent(in) :: ik
REAL(kind=PREC_DP)   , intent(in) :: xvalue
!
INTEGER :: i, ipkt, iclose 
REAL(kind=PREC_DP) :: xdist_old, xdist_new 
!                                                                       
xdist_old = abs (xvalue-x (offxy (ik - 1) + 1) ) 
iclose = 1 
DO i = 2, lenc (ik) 
   ipkt = offxy (ik - 1) + i 
   xdist_new = abs (xvalue-x (offxy (ik - 1) + i) ) 
   IF (xdist_new.lt.xdist_old) then 
      iclose = i 
      xdist_old = xdist_new 
   ENDIF 
ENDDO 
!                                                                       
closest_point = iclose 
!                                                                       
END FUNCTION closest_point                    
!
!*****7*****************************************************************
!
SUBROUTINE do_kmath (zeile, lp) 
!+                                                                      
!     Arithmetic between complete data sets                             
!-                                                                      
USE ber_params_mod
USE errlist_mod 
USE get_params_mod
USE prompt_mod 
USE kuplot_config 
USE kuplot_mod 
use kuplot_extrema_mod
use kuplot_show_mod
USE precision_mod
USE str_comp_mod
USE string_convert_mod
!                                                                       
IMPLICIT none 
!                                                                       
CHARACTER(len=*), intent(inout) :: zeile 
integer         , intent(inout) :: lp
!
INTEGER maxw 
PARAMETER (maxw = 5) 
INTEGER EXCL, NEUT, IGNO 
PARAMETER (EXCL = 0, NEUT = 1, IGNO = 2) 
!                                                                       
CHARACTER(LEN=PREC_STRING), dimension(MAXW) :: cpara (maxw) 
CHARACTER(len=3)                            :: oper 
REAL(KIND=PREC_DP) :: werte (maxw) 
INTEGER, dimension(MAXW) :: lpara !(maxw)
INTEGER :: imin 
INTEGER :: ianz, ik1, ik2, ik3, ikk1, ikk2, ikk3, ix, iy, i, ii 
LOGICAL :: loverwrite 
!                                                                       
!------ get parameters                                                  
!                                                                       
CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
IF (ier_num.ne.0) return 
loverwrite = str_comp (cpara (ianz) , 'overwrite', 2, lpara (ianz)&
      , 9)                                                              
IF (loverwrite) then 
   ianz = ianz - 1 
ENDIF 
!                                                                       
!------ check arguments                                                 
!                                                                       
IF (ianz.eq.4) then 
   IF (str_comp (cpara (4) , 'excl', 2, lpara (4) , 4) ) then 
      exclude9999 = EXCL 
   ELSEIF (str_comp (cpara (4) , 'neut', 2, lpara (4) , 4) ) then 
      exclude9999 = NEUT 
   ELSEIF (str_comp (cpara (4) , 'igno', 2, lpara (4) , 4) ) then 
      exclude9999 = IGNO 
   ELSE 
      ier_num = - 6 
      ier_typ = ER_COMM 
      RETURN 
   ENDIF 
   ianz = 3 
ELSEIF (ianz.eq.3) then 
   exclude9999 = EXCL 
ENDIF 
!                                                                       
IF (ianz.ne.3) then 
   ier_num = - 6 
   ier_typ = ER_COMM 
   RETURN 
ENDIF 
!                                                                       
oper = cpara (1) (1:3) 
CALL do_cap (oper) 
IF (oper.ne.'ADD'.and.oper.ne.'MUL'.and.oper.ne.'SUB'.and.oper.ne.&
     &'DIV') then                                                       
   ier_num = - 6 
   ier_typ = ER_COMM 
   RETURN 
ENDIF 
!                                                                       
CALL del_params (1, ianz, cpara, lpara, maxw) 
CALL ber_params (ianz, cpara, lpara, werte, maxw) 
IF (ier_num.ne.0) return 
!                                                                       
ik1 = nint (werte (1) ) 
ik2 = nint (werte (2) ) 
!                                                                       
IF (ik1.lt.1.or.ik1.gt. (iz - 1) .or.ik2.lt.1.or.ik2.gt. (iz - 1) &
      ) then                                                            
   ier_num = - 4 
   ier_typ = ER_APPL 
   RETURN 
ENDIF 
!                                                                       
!------ space left for new data set ??                                  
!                                                                       
IF (iz.gt.maxkurvtot.and..not.loverwrite) then 
   ier_num = - 1 
   ier_typ = ER_APPL 
   RETURN 
ENDIF 
!                                                                       
IF (loverwrite) then 
   ik3 = ik1 
   WRITE (output_io, 5000) ik3 
ELSE 
   ik3 = iz 
   WRITE (output_io, 5100) ik3 
ENDIF 
!                                                                       
 5000 FORMAT    (1x,'Result overwrites dataset ',i3,' ..') 
 5100 FORMAT    (1x,'Result creates new dataset ',i3,' ..') 
!                                                                       
!------ if we have a 3d data set                                        
!                                                                       
IF (lni (ik1) .and.lni (ik2) .and. (nx (ik1) .eq.nx (ik2) ) .and. &
    (ny (ik1) .eq.ny (ik2) ) ) then                                   
   DO ix = 1, nx (ik1) 
      DO iy = 1, ny (ik1) 
         ikk1 = offz (ik1 - 1) + (ix - 1) * ny (ik1) + iy 
         ikk2 = offz (ik2 - 1) + (ix - 1) * ny (ik1) + iy 
         ikk3 = offz (ik3 - 1) + (ix - 1) * ny (ik1) + iy 
         IF (exclude9999.eq.EXCL) then 
            IF (z (ikk1) .eq. - 9999..or.z (ikk2) .eq. - 9999.) then 
               z (ikk3) = - 9999.0 
            ELSE 
               IF (oper.eq.'MUL') z (ikk3) = z (ikk1) * z (ikk2) 
               IF (oper.eq.'SUB') z (ikk3) = z (ikk1) - z (ikk2) 
               IF (oper.eq.'ADD') z (ikk3) = z (ikk1) + z (ikk2) 
               IF (oper.eq.'DIV') then 
                  IF (z (ikk2) .ne.0.0) then 
                     z (ikk3) = z (ikk1) / z (ikk2) 
                  ELSE 
                     z (ikk3) = 0.0 
                  ENDIF 
               ENDIF 
            ENDIF 
         ELSEIF (exclude9999.eq.NEUT) then 
            IF (z (ikk1) .eq. - 9999..and.z (ikk2) .eq. - 9999.) then 
               z (ikk3) = - 9999.0 
            ELSEIF (z (ikk1) .eq. - 9999.) then 
               z (ikk3) = z (ikk2) 
            ELSEIF (z (ikk2) .eq. - 9999.) then 
               z (ikk3) = z (ikk1) 
            ELSE 
               IF (oper.eq.'MUL') z (ikk3) = z (ikk1) * z (ikk2) 
               IF (oper.eq.'SUB') z (ikk3) = z (ikk1) - z (ikk2) 
               IF (oper.eq.'ADD') z (ikk3) = z (ikk1) + z (ikk2) 
               IF (oper.eq.'DIV') then 
                  IF (z (ikk3) .ne.0.0) then 
                     z (ikk3) = z (ikk1) / z (ikk2) 
                  ELSE 
                     z (ikk3) = 0.0 
                  ENDIF 
               ENDIF 
            ENDIF 
         ELSEIF (exclude9999.eq.IGNO) then 
            IF (oper.eq.'MUL') z (ikk3) = z (ikk1) * z (ikk2) 
            IF (oper.eq.'SUB') z (ikk3) = z (ikk1) - z (ikk2) 
            IF (oper.eq.'ADD') z (ikk3) = z (ikk1) + z (ikk2) 
            IF (oper.eq.'DIV') then 
               IF (z (ikk3) .ne.0.0) then 
                  z (ikk3) = z (ikk1) / z (ikk2) 
               ELSE 
                  z (ikk3) = 0.0 
               ENDIF 
            ENDIF 
         ENDIF 
      ENDDO 
   ENDDO 
   DO i = 1, nx (ik1) 
      x (offxy (ik3 - 1) + i) = x (offxy (ik1 - 1) + i) 
   ENDDO 
   DO i = 1, ny (ik1) 
      y (offxy (ik3 - 1) + i) = y (offxy (ik1 - 1) + i) 
   ENDDO 
!                                                                       
   IF (.not.loverwrite) then 
      nx (iz) = nx (ik1) 
      ny (iz) = ny (ik1) 
      lenc (iz) = lenc (ik1) 
      lni (iz) = .true. 
      fname (iz) = 'result.dat' 
      fform (iz) = fform (ik1) 
      offxy (iz) = offxy (iz - 1) + lenc (iz) 
      offz (iz) = offz (iz - 1) + nx (ik1) * ny (ik1) 
      iz = iz + 1 
   ENDIF 
   CALL get_extrema_xy (x, ik3, nx (ik3), xmin, xmax) 
   CALL get_extrema_xy (y, ik3, ny (ik3), ymin, ymax) 
   CALL get_extrema_z (z, ik3, nx (ik3), ny (ik3), zmin, zmax) 
   CALL show_data (ik3) 
!                                                                       
!------ here for a 2d data set                                          
!                                                                       
ELSEIF (.not.lni (ik1) .and..not.lni (ik2) ) then 
        imin = min (lenc (ik1), lenc (ik2) ) 
   DO ii = 1, imin 
      ikk1 = offxy (ik1 - 1) + ii 
      ikk2 = offxy (ik2 - 1) + ii 
      ikk3 = offxy (ik3 - 1) + ii 
      x (ikk3) = x (ikk1) 
      IF (oper.eq.'MUL') then 
         y (ikk3) = y (ikk1) * y (ikk2) 
         dy (ikk3) = sqrt ( (y (ikk1) * dy (ikk2) ) **2 + (y (ikk2)  &
            * dy (ikk1) ) **2)                                          
      ELSEIF (oper.eq.'ADD') then 
         y (ikk3) = y (ikk1) + y (ikk2) 
         dy (ikk3) = sqrt (dy (ikk1) **2 + dy (ikk2) **2) 
      ELSEIF (oper.eq.'SUB') then 
         y (ikk3) = y (ikk1) - y (ikk2) 
         dy (ikk3) = sqrt (dy (ikk1) **2 + dy (ikk2) **2) 
      ELSEIF (oper.eq.'DIV') then 
         IF (y (ikk2) .ne.0.0) then 
            y (ikk3) = y (ikk1) / y (ikk2) 
            dy (ikk3) = sqrt ( (dy (ikk1) / y (ikk2) ) **2 + (dy (   &
            ikk2) * y (ikk1) / y (ikk2) **2) **2)                    
         ELSE 
            y (ikk3) = 0.0 
         ENDIF 
      ENDIF 
   ENDDO 
!                                                                       
   IF (.not.loverwrite) then 
      lenc (iz) = imin 
      fform (iz) = fform (ik1) 
      fname (iz) = 'result.dat' 
      offxy (iz) = offxy (iz - 1) + lenc (iz) 
      iz = iz + 1 
   ENDIF 
   CALL get_extrema_xy (x, ik3, lenc (ik3), xmin, xmax) 
   CALL get_extrema_xy (y, ik3, lenc (ik3), ymin, ymax) 
   CALL show_data (ik3) 
!                                                                       
!------ incompatible data sets                                          
!                                                                       
ELSE 
   ier_num = - 23 
   ier_typ = ER_APPL 
ENDIF 
!                                                                       
END SUBROUTINE do_kmath                       
!
!*****7*****************************************************************
!
SUBROUTINE do_merge (zeile, lp) 
!+                                                                      
!     Merge different data sets                                         
!-                                                                      
USE ber_params_mod
USE errlist_mod 
USE get_params_mod
USE kuplot_config 
USE kuplot_mod 
use kuplot_extrema_mod
use kuplot_show_mod
USE precision_mod
USE str_comp_mod
!                                                                       
IMPLICIT none 
!                                                                       
CHARACTER (LEN=*), INTENT(INOUT) :: zeile 
INTEGER          , INTENT(INOUT) :: lp
!
INTEGER maxw 
PARAMETER (maxw = MAXKURVTOT) 
!                                                                       
CHARACTER(LEN=PREC_STRING), dimension(:), allocatable :: cpara !(maxw) 
REAL(KIND=PREC_DP) :: werte (maxw) 
REAL(kind=PREC_DP) :: mdelta, mmin, mmax 
REAL(kind=PREC_DP) ::  mdeltay, mminy, mmaxy 
INTEGER :: lpara (maxw)
INTEGER :: ianz, ik, ibin, i, ip, ntot, maxpp 
INTEGER :: ntoty
INTEGER :: maxzz
INTEGER :: ix, iy       ! Dummy loop indices
INTEGER :: ind_in, ind_out     ! Indices in inout and output 2D data
INTEGER, DIMENSION(:), ALLOCATABLE :: npt ! number of data points
LOGICAL :: lvalid, ladd, lall 
!                                                                       
!                                                                       
!------ space left for new data set ??                                  
!                                                                       
IF (iz.gt.maxkurvtot) then 
   ier_num = - 1 
   ier_typ = ER_APPL 
   RETURN 
ENDIF 
!
allocate(cpara(maxw))
!                                                                       
!------ get parameters                                                  
!                                                                       
CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
IF (ier_num.ne.0) then
   deallocate(cpara)
   return 
endif
!                                                                       
ladd = str_comp (cpara (ianz) , 'add', 1, lpara (ianz) , 3) 
IF (ladd) ianz = ianz - 1 
!                                                                       
lall = str_comp (cpara (ianz) , 'all', 1, lpara (ianz) , 3) 
IF (lall) then 
   werte (1) = - 1 
   deallocate(cpara)
ELSE 
!                                                                       
   CALL ber_params (ianz, cpara, lpara, werte, maxw) 
   deallocate(cpara)
   IF (ier_num.ne.0) then
      return 
   endif
ENDIF 
!                                                                       
!------ check arguments                                                 
!                                                                       
IF (nint (werte (1) ) .ne. - 1.and.ianz.lt.2) then 
   ier_num = - 6 
   ier_typ = ER_COMM 
   RETURN 
ENDIF 
!                                                                       
IF (ianz.eq.1) then 
   IF (nint (werte (1) ) .eq. - 1) then 
      DO i = 1, iz - 1 
         werte (i) = i 
      ENDDO 
      ianz = iz - 1 
   ENDIF 
ENDIF 
!                                                                       
!------ First we get extend of new data set (first one defines DELTA)   
!                                                                       
ik = nint (werte (1) ) 
IF(.NOT.lni(ik)) THEN                    ! 1-D DATA SETS
!                                                                       
   lvalid = .true. 
   DO i = 1, ianz 
      ik = nint (werte (i) ) 
      lvalid = lvalid.and.ik.ge.1.and.ik.lt.iz 
      IF (lvalid) then 
         lvalid = lvalid.and..not.lni (ik) 
      ENDIF 
   ENDDO 
   IF (.not.lvalid) then 
      ier_num = - 4 
      ier_typ = ER_APPL 
      RETURN 
   ENDIF 
!
   mdelta = x(offxy(ik - 1) + lenc(ik)) - x(offxy(ik - 1) + 1)
   mdelta = mdelta / REAL(lenc (ik) - 1) 
   mmin = xmin (nint (werte (1) ) ) 
   mmax = xmax (nint (werte (1) ) ) 
   DO i = 2, ianz 
      ik = nint (werte (i) ) 
      mmax = max (mmax, xmax (ik) ) 
      mmin = min (mmin, xmin (ik) ) 
   ENDDO 
!                                                                       
   maxpp = maxarray - offxy (iz - 1) 
   ntot = nint ( (mmax - mmin) / mdelta) + 1 
!                                                                       
   IF (ntot.gt.maxpp) then 
      ier_num = - 6 
      ier_typ = ER_APPL 
      RETURN 
   ENDIF 
!                                                                       
!------ Now we merge the data one by one ...                            
!                                                                       
   DO i = 1, ntot 
      ibin = offxy (iz - 1) + i 
      x (ibin) = 0.0 
      y (ibin) = 0.0 
      dx (ibin) = 0.0 
      dy (ibin) = 0.0 
   ENDDO 
!                                                                       
   DO i = 1, ianz 
      ik = nint (werte (i) ) 
      DO ip = 1, lenc (ik) 
         ibin = offxy(iz - 1) + nint((x(offxy(ik - 1) + ip) - mmin) / mdelta) + 1
         y (ibin) = y (ibin) + y (offxy (ik - 1) + ip) 
         dy (ibin) = dy (ibin) + dy (offxy (ik - 1) + ip) **2 
         dx (ibin) = dx (ibin) + 1.0 
      ENDDO 
   ENDDO 
!                                                                       
   DO i = 1, ntot 
      ibin = offxy (iz - 1) + i 
      IF (.not.ladd.and.dx (ibin) .gt.0.0) then 
         y (ibin) = y (ibin) / dx (ibin) 
         dy (ibin) = sqrt (dy (ibin) ) / dx (ibin) 
      ENDIF 
      x (ibin) = mmin + (i - 1) * mdelta 
      dx (ibin) = 0.0 
   ENDDO 
!                                                                       
      lenc (iz) = ntot 
      fform (iz) = fform (nint (werte (1) ) ) 
      offxy (iz) = offxy (iz - 1) + lenc (iz) 
!
   fname(iz) = 'merge.dat' 
   iz = iz + 1 
   CALL get_extrema_xy (x, iz - 1, lenc (iz), xmin, xmax) 
   CALL get_extrema_xy (y, iz - 1, lenc (iz), ymin, ymax) 
   CALL show_data (iz - 1) 
ELSE                               ! 2D Data sets
!                                                                       
   lvalid = .true. 
   DO i = 1, ianz 
      ik = NINT(werte(i)) 
      lvalid = lvalid .AND. ik >= 1 .AND. ik < iz  .AND. lni(ik)
   ENDDO 
   IF (.NOT.lvalid) THEN 
      ier_num = - 4 
      ier_typ = ER_APPL 
      RETURN 
   ENDIF 
!
   mdelta  = (xmax(ik) - xmin(ik)) / REAL(nx(ik) - 1)   ! Step in x
   mdeltay = (ymax(ik) - ymin(ik)) / REAL(ny(ik) - 1)   ! Step in y
   mmin    = xmin(NINT(werte(1)))                        ! Min/max in x/y
   mmax    = xmax(NINT(werte(1))) 
   mminy   = ymin(NINT(werte(1))) 
   mmaxy   = ymax(NINT(werte(1))) 
   DO i = 2, ianz 
      ik = NINT(werte(i)) 
      mmax  = MAX(mmax, xmax(ik)) 
      mmin  = MIN(mmin, xmin(ik)) 
      mmaxy = MAX(mmaxy, ymax(ik)) 
      mminy = MIN(mminy, ymin(ik)) 
   ENDDO 
   ntot  = NINT((mmax - mmin) / mdelta) + 1             ! Data points along x
   ntoty = NINT((mmaxy - mminy) / mdeltay) + 1          ! Data points along y 
!
   maxpp  = maxarray - offxy (iz - 1)
   maxzz  = maxarray - offz (iz - 1)
   DO i=1, ntot
      x(offxy(iz-1)+i) = mmin + (i-1)*mdelta
   ENDDO 
   DO i=1, ntoty
      y(offxy(iz-1)+i) = mminy + (i-1)*mdeltay
   ENDDO 
   DO i=1,ntot*ntoty           ! Zero data in new field
      z(offz(iz-1)+i) = 0.0
   ENDDO
!
   ALLOCATE(npt(ntot*ntoty))
   npt(:) = 0
!
   DO i = 1, ianz
      ik = NINT(werte(i))
         DO iy=1, ny(ik)
      DO ix=1, nx(ik)
            ind_in = offz(ik-1)+(ix-1)*ny(ik) + iy
            ind_out=  (NINT((x(offxy(ik-1)+ix)-mmin )/mdelta )  )*ntoty &
                     +(NINT((y(offxy(ik-1)+iy)-mminy)/mdeltay)+1)
!write(*,'(7i3,3f6.1,i3)') i,ix, iy, ind_in, ind_out, (NINT((x(offxy(ik-1)+ix)-mmin )/mdelta )  )*ntoty &
! , (NINT((y(offxy(ik-1)+iy)-mminy)/mdeltay)+1), x(offxy(ik-1)+ix), mmin,mdelta, ntoty
            npt(ind_out) = npt(ind_out) + 1             ! Increment number of contribution data
            ind_out      = ind_out + offz(iz-1)
            z(ind_out)   = z(ind_out) + z(ind_in)
         ENDDO
      ENDDO
   ENDDO
!
   DO ibin=1, ntot*ntoty
      IF(.NOT.ladd .AND. npt(ibin) > 0) THEN
         ind_out = offz(iz-1) + ibin
         z(ind_out) = z(ind_out) / npt(ibin)
      ENDIF
   ENDDO
!
   DEALLOCATE(npt)
!
   DO ix=1,ntot
      x(offxy(iz-1)+ix) = mmin + (ix-1)*mdelta
   ENDDO
   DO iy=1,ntot
      y(offxy(iz-1)+iy) = mminy + (iy-1)*mdeltay
   ENDDO
!
   lni(iz)   = .TRUE.
   fname(iz) = 'merge.dat' 
   fform(iz) = fform(NINT(werte(1)))
   nx(iz)    = ntot 
   ny(iz)    = ntoty
   lenc(iz)   = MAX(nx(iz), ny(iz))
   offxy(iz) = offxy(iz-1) + lenc(iz)
   offz (iz) = offz (iz-1) + nx(iz)*ny(iz)
   iz = iz + 1 
   CALL get_extrema_xy (x, iz - 1, lenc (iz-1), xmin, xmax) 
   CALL get_extrema_xy (y, iz - 1, lenc (iz-1), ymin, ymax) 
   CALL get_extrema_z  (z, iz - 1, nx(iz-1), ny(iz-1), zmin, zmax) 
   CALL show_data(iz - 1) 
ENDIF
!
!                                                                       
END SUBROUTINE do_merge                       
!
!*****7*****************************************************************
!
SUBROUTINE do_rebin (zeile, lp) 
!+                                                                      
!     Rebin data to given grid (only 2D - otherwise save/read GNU)      
!-                                                                      
USE ber_params_mod
USE errlist_mod 
USE get_params_mod
USE kuplot_config 
USE kuplot_mod 
use kuplot_extrema_mod
use kuplot_show_mod
USE precision_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER maxw 
PARAMETER (maxw = 25) 
!                                                                       
CHARACTER (LEN=*), INTENT(INOUT) :: zeile 
INTEGER          , INTENT(INOUT) :: lp
!
CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
REAL(KIND=PREC_DP) :: werte (maxw)
REAL(kind=PREC_DP) :: bdelta, fraction 
INTEGER :: lpara (maxw), maxpp, ntot 
INTEGER :: ianz, ik, ibin, jbin, ip, i, j, k, iupper, iubin 
!                                                                       
!------ space left for new data set ??                                  
!                                                                       
IF (iz.gt.maxkurvtot) then 
   ier_num = - 1 
   ier_typ = ER_APPL 
   RETURN 
ENDIF 
!                                                                       
!------ get parameters                                                  
!                                                                       
CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
IF (ier_num.ne.0) return 
!                                                                       
!------ check arguments                                                 
!                                                                       
IF (ianz.ne.2) then 
   ier_num = - 6 
   ier_typ = ER_COMM 
   RETURN 
ENDIF 
!                                                                       
CALL ber_params (ianz, cpara, lpara, werte, maxw) 
IF (ier_num.ne.0) return 
!                                                                       
ik = nint (werte (1) ) 
bdelta = werte (2) 
!                                                                       
IF (ik.lt.1.or.ik.gt.iz.or.lni (ik) ) then 
   ier_num = - 4 
   ier_typ = ER_APPL 
   RETURN 
ENDIF 
!                                                                       
maxpp = maxarray - offxy (iz - 1) 
ntot = nint ( (xmax (ik) - xmin (ik) ) / bdelta) + 1 
!                                                                       
IF (ntot.gt.maxpp) then 
   ier_num = - 6 
   ier_typ = ER_APPL 
   RETURN 
ENDIF 
!                                                                       
!------ Now we rebin the data                                           
!                                                                       
DO i = 1, ntot 
   ibin = offxy (iz - 1) + i 
   x (ibin) = 0.0 
   y (ibin) = 0.0 
   dx (ibin) = 0.0 
   dy (ibin) = 0.0 
ENDDO 
!                                                                       
!------ Find closest bin and add data point to it                       
!                                                                       
DO ip = 1, lenc (ik) 
   ibin = offxy(iz - 1) + nint((x(offxy(ik - 1) + ip) - xmin(ik) ) / bdelta) + 1
   y (ibin) = y (ibin) + y (offxy (ik - 1) + ip) 
   dy (ibin) = dy (ibin) + dy (offxy (ik - 1) + ip) **2 
   dx (ibin) = dx (ibin) + 1.0 
ENDDO 
!                                                                       
!------ Normalize new data points                                       
!                                                                       
DO i = 1, ntot 
   ibin = offxy (iz - 1) + i 
   IF (dx (ibin) .gt.0.0) then 
      y (ibin) = y (ibin) / dx (ibin) 
      dy (ibin) = sqrt (dy (ibin) ) / dx (ibin) 
   ENDIF 
ENDDO 
!                                                                       
!------ Set x arrays                                                    
!                                                                       
DO i = 1, ntot 
   ibin = offxy (iz - 1) + i 
   x (ibin) = xmin (ik) + (i - 1) * bdelta 
ENDDO 
!                                                                       
!------ Find empty bins and assign average value from neighbours        
!------ needs to be done                                                
!                                                                       
DO i = 2, ntot - 1 
   ibin = offxy (iz - 1) + i 
   IF (dx (ibin) .eq.0.0) then 
      iupper = 0 
      loop_empty: DO k = i + 1, ntot 
         jbin = offxy (iz - 1) + k 
         IF (dx (jbin) .gt.0.0) then 
            iupper = k 
            exit loop_empty
         ENDIF 
      ENDDO  loop_empty
!  71    CONTINUE 
      continue
      iubin = offxy (iz - 1) + iupper 
      fraction = (y(iubin) - y(ibin - 1) ) / (x(iubin) - x(ibin - 1) )
      DO j = i, iupper - 1 
         jbin = offxy (iz - 1) + j 
         y (jbin) = y (ibin - 1) + fraction * (x (jbin) - x (ibin - 1) ) 
         dy (jbin) = dy (ibin - 1) 
         dx (jbin) = 1.0 
      ENDDO 
   ENDIF 
ENDDO 
!                                                                       
!------ Set dx arrays                                                   
!                                                                       
DO i = 1, ntot 
   ibin = offxy (iz - 1) + i 
   dx (ibin) = 0.0 
ENDDO 
!                                                                       
lenc (iz) = ntot 
fform (iz) = fform (ik) 
fname (iz) = 'rebin.dat' 
offxy (iz) = offxy (iz - 1) + lenc (iz) 
iz = iz + 1 
CALL get_extrema_xy (x, iz - 1, lenc (iz), xmin, xmax) 
CALL get_extrema_xy (y, iz - 1, lenc (iz), ymin, ymax) 
CALL show_data (iz - 1) 
!                                                                       
END SUBROUTINE do_rebin                       
!
!*****7*****************************************************************
!
SUBROUTINE do_interpolate(zeile, lp) 
!+                                                                      
!     Interpolate data set <ik> on grid of <ig> using spline            
!-                                                                      
USE ber_params_mod
USE errlist_mod 
USE get_params_mod
USE kuplot_config 
USE kuplot_mod 
use kuplot_extrema_mod
use kuplot_math_mod
use kuplot_show_mod
!
USE precision_mod
use spline_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER maxw 
PARAMETER (maxw = 2) 
!                                                                       
CHARACTER (LEN=*), INTENT(INOUT) :: zeile 
INTEGER          , INTENT(INOUT) :: lp
!
CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
REAL(KIND=PREC_DP) :: werte (maxw)
REAL(kind=PREC_DP) :: x_igg 
REAL(kind=PREC_DP) :: yyy 
REAL(kind=PREC_DP),dimension(:), allocatable :: xpl !(maxarray)
REAL(kind=PREC_DP),dimension(:), allocatable :: ypl !(maxarray)
REAL(kind=PREC_DP),dimension(:), allocatable :: y2a !(maxarray)
!REAL(kind=PREC_DP) :: y2a (maxarray) 
INTEGER :: lpara (maxw), maxpp, ntot 
INTEGER :: ianz, ibin, igg, ik, ig, i , j, k
!logical :: ascend        ! x-data increase=true or decrease=false
!                                                                       
!------ space left for new data set ??                                  
!                                                                       
IF (iz.gt.maxkurvtot) then 
   ier_num = - 1 
   ier_typ = ER_APPL 
   RETURN 
ENDIF 
!                                                                       
!------ get parameters                                                  
!                                                                       
CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
IF (ier_num.ne.0) return 
!                                                                       
!------ check arguments                                                 
!                                                                       
IF (ianz.ne.2) then 
   ier_num = - 6 
   ier_typ = ER_COMM 
   RETURN 
ENDIF 
!                                                                       
CALL ber_params (ianz, cpara, lpara, werte, maxw) 
IF (ier_num.ne.0) return 
!                                                                       
ik = nint (werte (1) ) 
ig = nint (werte (2) ) 
!                                                                       
IF (   ik.lt.1.or.ik.gt.iz.or.lni (ik)                               &
   .or.ig.lt.1.or.ig.gt.iz.or.lni (ig) .or.ik.eq.ig) then            
   ier_num = - 4 
   ier_typ = ER_APPL 
   RETURN 
ENDIF 
!                                                                       
maxpp = maxarray - offxy (iz - 1) 
ntot = lenc (ig) 
!                                                                       
IF (ntot.gt.maxpp) then 
   ier_num = - 6 
   ier_typ = ER_APPL 
   RETURN 
ENDIF 
!
allocate(xpl(1:lenc(ik)))
allocate(ypl(1:lenc(ik)))
allocate(y2a(1:lenc(ik)))
!                                                                       
!------ Now we interpolate the data                                     
!                                                                       
!ascend = x (offxy (ik - 1) +1) < x (offxy (ik - 1) + lenc(ik))
!write(*,*) 'x(1), x(E) ', x (offxy (ik - 1) +1), x (offxy (ik - 1) + lenc(ik))
!write(*,*) 'min , max K', xmin(ik), xmax(ik), ascend
!write(*,*) 'min , max G', xmin(ig), xmax(ig)
j = 0
loop_copy: DO i = 1, lenc (ik) 
!   if(ascend) then
      if(x (offxy (ik - 1) + i) < xmin(ig)              ) cycle loop_copy
      if(x (offxy (ik - 1) + i) > xmax(ig)              ) exit  loop_copy
!   else
!      if(x (offxy (ik - 1) + i) < xmax(ig)              ) cycle loop_copy
!      if(x (offxy (ik - 1) + i) > xmin(ig)              ) exit  loop_copy
!   endif
   j = j + 1
   xpl (j) = x (offxy (ik - 1) + i) 
   ypl (j) = y (offxy (ik - 1) + i) 
ENDDO loop_copy
!write(*,*) 'lengths ', lenc(ik), lenc(ig), j
!write(*,*) 'x(1), x(E) ', x (offxy (ig - 1) +1), x (offxy (ig - 1) + lenc(ig))
!write(*,*) 'xpl , xpl  ', xpl(1), xpl(j)
!                                                                       
!CALL spline(xpl, ypl, lenc (ik), 1D30, 1D30, y2a) 
!CALL spline(lenc(ik), xpl, ypl, 1D30, 1D30, y2a) 
CALL spline(j       , xpl, ypl, 1D30, 1D30, y2a) 
!                                                                       
k = 0
loop_spline:DO i = 1, ntot 
!   if(ascend) then
!     if(x(offxy (ig - 1) + i)<minval(xpl)) cycle loop_spline
!     if(x(offxy (ig - 1) + i)>maxval(xpl)) exit  loop_spline
      if(x(offxy (ig - 1) + i)<xpl(1)     ) cycle loop_spline
      if(x(offxy (ig - 1) + i)>xpl(j)     ) exit  loop_spline
!   else
!      if(x(offxy (ig - 1) + i)<xpl(j)     ) cycle loop_spline
!      if(x(offxy (ig - 1) + i)>xpl(1)     ) exit  loop_spline
!   endif
   k = k + 1
!  ibin = offxy (iz - 1) + i 
   ibin = offxy (iz - 1) + k
   igg = offxy (ig - 1) + i 
!  CALL splint(xpl, ypl, y2a, lenc (ik), x (igg), yyy,ier_num)
   x_igg = real(x(igg), kind=PREC_DP)
!  CALL splint(lenc(ik), xpl, ypl, y2a, x_igg, yyy,ier_num)
   CALL splint(j       , xpl, ypl, y2a, x_igg, yyy,ier_num)
   IF(ier_num /= 0) THEN
         ier_typ = ER_APPL
         RETURN
   ENDIF 
!                                                                       
   x (ibin) = x(igg) 
   y (ibin) = yyy 
   dx(ibin) = dx(igg)
   dy(ibin) = dy(igg)
ENDDO loop_spline
!write(*,*) ' LENC  ', k
!
deallocate(xpl)
deallocate(ypl)
deallocate(y2a)
!                                                                       
!lenc (iz) = ntot 
lenc (iz) = k 
fform (iz) = fform (ig) 
fname (iz) = 'interpolate.dat' 
offxy (iz) = offxy (iz - 1) + lenc (iz) 
iz = iz + 1 
CALL get_extrema_xy (x, iz - 1, lenc (iz), xmin, xmax) 
CALL get_extrema_xy (y, iz - 1, lenc (iz), ymin, ymax) 
CALL show_data (iz - 1) 
!                                                                       
END SUBROUTINE do_interpolate                 
!
!*****7*****************************************************************
!
SUBROUTINE do_convolute (zeile, lp) 
!+                                                                      
!     Convoluting two data sets                                         
!-                                                                      
USE ber_params_mod
USE errlist_mod 
USE get_params_mod
USE kuplot_config 
USE kuplot_mod 
use kuplot_extrema_mod
use kuplot_show_mod
USE precision_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER maxw 
PARAMETER (maxw = 2) 
!                                                                       
CHARACTER (LEN=*), INTENT(INOUT) :: zeile 
INTEGER          , INTENT(INOUT) :: lp
!
CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
REAL(KIND=PREC_DP) :: werte (maxw) 
REAL(kind=PREC_DP) :: dr1, dr2, ysum, start
INTEGER :: lpara (maxw), maxpp 
INTEGER :: ianz, ik1, ik2, ii, i, j 
!                                                                       
!------ space left for new data set ??                                  
!                                                                       
IF (iz.gt.maxkurvtot) then 
   ier_num = - 1 
   ier_typ = ER_APPL 
   RETURN 
ENDIF 
!                                                                       
!------ get parameters                                                  
!                                                                       
CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
IF (ier_num.ne.0) return 
!                                                                       
!------ check arguments                                                 
!                                                                       
IF (ianz.ne.2) then 
   ier_num = - 6 
   ier_typ = ER_COMM 
   RETURN 
ENDIF 
!                                                                       
CALL ber_params (ianz, cpara, lpara, werte, maxw) 
IF (ier_num.ne.0) return 
!                                                                       
ik1 = nint (werte (1) ) 
ik2 = nint (werte (2) ) 
!                                                                       
IF (ik1.lt.1.or.ik1.gt.iz.or.lni (ik1)                            &
      .or.ik2.lt.1.or.ik2.gt.iz.or.lni (ik2) .or.ik1.eq.ik2) then       
   ier_num = - 4 
   ier_typ = ER_APPL 
   RETURN 
ENDIF 
!                                                                       
dr1 = (xmax (ik1) - xmin (ik1) ) / REAL(lenc (ik1) - 1) 
dr2 = (xmax (ik2) - xmin (ik2) ) / REAL(lenc (ik2) - 1) 
!                                                                       
IF (abs (dr1 - dr2) .gt.1e-8) then 
   ier_num = - 45 
   ier_typ = ER_APPL 
ENDIF 
!                                                                       
maxpp = maxarray - offxy (iz - 1) 
IF (lenc (ik1) .gt.maxpp) then 
   ier_num = - 6 
   ier_typ = ER_APPL 
ENDIF 
IF (ier_num.ne.0) return 
!                                                                       
!------ Now do convolution                                              
!                                                                       
start = x (offxy (ik2 - 1) + 1) 
!                                                                       
DO i = 1, lenc (ik1) 
   ysum = 0.0 
   DO j = 1, lenc (ik1) 
      ii = (i - j) - nint (start / dr1) + 1 
      IF (ii.ge.1.and.ii.le.lenc (ik2) ) then 
         ysum = ysum + y (offxy (ik1 - 1) + j) * y (offxy (ik2 - 1)     &
         + ii)                                                          
      ENDIF 
   ENDDO 
!                                                                       
   x (offxy (iz - 1) + i) = x (offxy (ik1 - 1) + i) 
   y (offxy (iz - 1) + i) = ysum * dr1 
!                                                                       
   dx (offxy (iz - 1) + i) = 0.0 
   dy (offxy (iz - 1) + i) = 0.0 
ENDDO 
!                                                                       
!------ Do other settings                                               
!                                                                       
lenc (iz) = lenc (ik1) 
fform (iz) = fform (ik1) 
fname (iz) = 'convolute.dat' 
offxy (iz) = offxy (iz - 1) + lenc (iz) 
iz = iz + 1 
CALL get_extrema_xy (x, iz - 1, lenc (iz), xmin, xmax) 
CALL get_extrema_xy (y, iz - 1, lenc (iz), ymin, ymax) 
CALL show_data (iz - 1) 
!                                                                       
END SUBROUTINE do_convolute                   
!
!*****7*****************************************************************
!
SUBROUTINE do_derivative (zeile, lp) 
!+                                                                      
!     Calculate n-th derivative of given data set.                      
!-                                                                      
USE ber_params_mod
USE errlist_mod 
USE get_params_mod
USE kuplot_config 
USE kuplot_mod 
use kuplot_extrema_mod
use kuplot_math_mod
use kuplot_show_mod
USE precision_mod
USE prompt_mod 
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER :: maxw, maxsm 
!                                                                       
PARAMETER (maxw = 2) 
PARAMETER (maxsm = 51) 
!                                                                       
CHARACTER (LEN=*), INTENT(INOUT) :: zeile 
INTEGER          , INTENT(INOUT) :: lp
!
CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
REAL(KIND=PREC_DP) :: werte(maxw) 
REAL(KIND=PREC_SP) :: cc(maxsm) 
REAL(KIND=PREC_DP) :: deltax, fnorm 
INTEGER :: lpara (maxw)
INTEGER :: ianz, i, ik, id, ip, im, inew, iold, maxpp 
!                                                                       
!------ space left for new data set ??                                  
!                                                                       
IF (iz.gt.maxkurvtot) then 
   ier_num = - 1 
   ier_typ = ER_APPL 
   RETURN 
ENDIF 
!                                                                       
!------ get parameters                                                  
!                                                                       
CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
IF (ier_num.ne.0) return 
!                                                                       
!------ check arguments                                                 
!                                                                       
IF (ianz.ne.1.and.ianz.ne.2) then 
   ier_num = - 6 
   ier_typ = ER_COMM 
   RETURN 
ENDIF 
!                                                                       
CALL ber_params (ianz, cpara, lpara, werte, maxw) 
IF (ier_num.ne.0) return 
!                                                                       
id = 1 
ik = nint (werte (1) ) 
IF (ianz.eq.2) id = nint (werte (2) ) 
!                                                                       
IF (id.lt.1.or.id.gt.4) then 
   ier_num = - 40 
   ier_typ = ER_APPL 
   RETURN 
ENDIF 
!                                                                       
IF (ik.lt.1.or.ik.gt.iz.or.lni (ik) ) then 
   ier_num = - 4 
   ier_typ = ER_APPL 
   RETURN 
ENDIF 
!                                                                       
maxpp = maxarray - offxy (iz - 1) 
IF (lenc (ik) .gt.maxpp) then 
   ier_num = - 6 
   ier_typ = ER_APPL 
   RETURN 
ENDIF 
!                                                                       
WRITE ( output_io, 1000) id, ik 
 1000 FORMAT     (1x,'------ > Calculating derivative ',i2,             &
     &                      ' for data set ',i3,' ...')                 
!                                                                       
!------ Some setup                                                      
!                                                                       
im = 6 
ip = 2 * im + 1 
!                                                                       
!------ First we copy the data to new data set                          
!                                                                       
lenc (iz) = lenc (ik) 
fform (iz) = fform (ik) 
fname (iz) = 'deriv.dat' 
!                                                                       
DO i = 1, lenc (ik) 
   inew = offxy (iz - 1) + i 
   iold = offxy (ik - 1) + i 
   x (inew) = x (iold) 
   y (inew) = y (iold) 
   dx (inew) = dx (iold) 
   dy (inew) = dy (iold) 
ENDDO 
!                                                                       
!------ Set up smoothing to return derivative                           
!                                                                       
CALL SAVGOL(cc, ip, ip / 2, ip / 2, id, im) 
IF (ier_num.ne.0) return 
CALL do_glatt_y (iz, ip, lenc (iz), cc, maxsm, .true.) 
IF (ier_num.ne.0) return 
!                                                                       
!------ Divide by the step size                                         
!                                                                       
CALL get_extrema_xy (x, iz, lenc (iz), xmin, xmax) 
deltax = (xmax (iz) - xmin (iz) ) / REAL(lenc (iz) - 1) 
!                                                                       
fnorm = 1.0 
DO i = 2, id 
   fnorm = fnorm * REAL(i) 
ENDDO 
!                                                                       
fnorm = deltax**id / fnorm 
!                                                                       
DO i = 1 + im, lenc (iz) - im 
   inew = offxy (iz - 1) + i - im 
   iold = offxy (iz - 1) + i 
   x (inew) = x (iold) 
   y (inew) = y (iold) / fnorm 
   dx (inew) = dx (iold) 
   dy (inew) = dy (iold) 
ENDDO 
!                                                                       
lenc (iz) = lenc (iz) - 2 * im 
offxy (iz) = offxy (iz - 1) + lenc (iz) 
!                                                                       
!------ Finish off                                                      
!                                                                       
iz = iz + 1 
!                                                                       
CALL get_extrema_xy (y, iz - 1, lenc (iz - 1), ymin, ymax) 
CALL show_data (iz - 1) 
!                                                                       
END SUBROUTINE do_derivative                  
!
!*****7*****************************************************************
!
SUBROUTINE do_match (zeile, lp) 
!+                                                                      
!     Calculates scaling and offset to give best match between          
!     two data sets.                                                    
!-                                                                      
USE ber_params_mod
USE errlist_mod 
USE get_params_mod
USE param_mod 
USE prompt_mod 
USE kuplot_config 
USE kuplot_mod 
use kuplot_extrema_mod
use kuplot_show_mod
use kuplot_wichtung_mod
USE string_convert_mod
USE precision_mod
USE str_comp_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER maxw 
PARAMETER (maxw = 25) 
!                                                                       
CHARACTER (LEN=*), INTENT(INOUT) :: zeile 
INTEGER          , INTENT(INOUT) :: lp
!
CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
REAL(KIND=PREC_DP) :: werte (maxw) 
REAL(kind=PREC_DP) :: sk, ba, wtot, ee, cc, ce, c, e, wi 
INTEGER :: lpara (maxw)
INTEGER :: ianz, iko, ikc, ikko, ikkc, ip 
INTEGER :: iianz 
INTEGER :: iweight 
INTEGER :: i, j 
LOGICAL :: lback, lskal, lmix 
real(kind=PREC_DP) :: bck_k = 0.0D0
!                                                                       
integer, PARAMETER :: W_ONE = 0
integer, PARAMETER :: W_SQUA = 1 
integer, PARAMETER :: W_SQRT = 2 
integer, PARAMETER :: W_INV = 3 
integer, PARAMETER :: W_LOG = 4 
integer, PARAMETER :: W_ISQ = 5 
integer, PARAMETER :: W_LIN = 6 
integer, PARAMETER :: W_DAT = 7 
integer, PARAMETER :: W_BCK = 8 
!                                                                       
!                                                                       
!     REAL r_wichtung 
!                                                                       
!------ get parameters                                                  
!                                                                       
CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
IF (ier_num.ne.0) return 
!                                                                       
!------ check arguments                                                 
!                                                                       
IF (ianz.lt.3) then 
   ier_num = - 6 
   ier_typ = ER_COMM 
   RETURN 
ENDIF 
!                                                                       
!------ Getting arguments                                               
!                                                                       
lskal = str_comp (cpara (1) , 'scal', 1, lpara (1) , 4)  .or.     &
        str_comp (cpara (1) , 'all', 1, lpara (1) , 3)                
lback = str_comp (cpara (1) , 'offset', 1, lpara (1) , 6) .or.    &
        str_comp (cpara (1) , 'all', 1, lpara (1) , 3)                
lmix = str_comp (cpara (1) , 'mix', 1, lpara (1) , 3) 
!                                                                       
IF (lmix) then 
   CALL do_match_mix (ianz, cpara, lpara, maxw) 
   RETURN 
ENDIF 
!                                                                       
CALL del_params (1, ianz, cpara, lpara, maxw) 
iianz = min (ianz, 2) 
CALL ber_params (ianz, cpara, lpara, werte, maxw) 
IF (ier_num.ne.0) return 
!                                                                       
iko = nint (werte (1) ) 
ikc = nint (werte (2) ) 
!                                                                       
CALL del_params (2, ianz, cpara, lpara, maxw) 
CALL do_low(cpara(1)(1:lpara(1)))
iweight = W_ONE 
IF (ianz.eq.0.or.str_comp(cpara(1) , 'one', 1, lpara(1), 3)) then
   iweight = W_ONE 
ELSEIF (str_comp (cpara (1) , 'squa', 3, lpara (1) , 4) ) then 
   iweight = W_SQUA 
ELSEIF (str_comp (cpara (1) , 'sqrt', 3, lpara (1) , 4) ) then 
   iweight = W_SQRT 
ELSEIF (str_comp (cpara (1) , 'inv', 1, lpara (1) , 3) ) then 
   iweight = W_INV 
ELSEIF (str_comp (cpara (1) , 'log', 2, lpara (1) , 3) ) then 
   iweight = W_LOG 
ELSEIF (str_comp (cpara (1) , 'isq', 1, lpara (1) , 3) ) then 
   iweight = W_ISQ 
ELSEIF (str_comp (cpara (1) , 'lin', 2, lpara (1) , 3) ) then 
   iweight = W_LIN 
ELSEIF (str_comp (cpara (1) , 'dat', 1, lpara (1) , 3) ) then 
   iweight = W_DAT 
ELSEIF (str_comp (cpara (1) , 'bck', 1, lpara (1) , 3) ) then 
   iweight = W_BCK 
   iweight = W_DAT 
ELSE
   ier_num = -27
   ier_typ = ER_APPL 
   RETURN 
ENDIF 
!                                                                       
IF (iko.lt.1.or.iko.gt.iz.or.ikc.lt.1.or.ikc.gt.iz.or.iko.eq.ikc) then
   ier_num = - 4 
   ier_typ = ER_APPL 
   RETURN 
ENDIF 
!                                                                       
!------ Now we calculate scale and offset                               
!                                                                       
wtot = 0.0 
e = 0.0 
c = 0.0 
ee = 0.0 
cc = 0.0 
ce = 0.0 
!                                                                       
wi = 1.0 
!                                                                       
!     xy-data                                                           
!                                                                       
IF (.not.lni (iko) .and..not.lni (ikc) ) then 
   DO ip = 1, lenc (ikc) 
      ikko = offxy (iko - 1) + ip 
      ikkc = offxy (ikc - 1) + ip 
      wi = r_wichtung (y (ikko), dy (ikko), iweight, z(ikc), bck_k) 
      wtot = wtot + wi 
      e = e+wi * y (ikko) 
      ee = ee+wi * y (ikko) **2 
      c = c + wi * y (ikkc) 
      cc = cc + wi * y (ikkc) **2 
      ce = ce+wi * y (ikkc) * y (ikko) 
   ENDDO 
!                                                                       
!      z-data                                                           
!                                                                       
ELSEIF (lni (iko) .and.lni (ikc) ) then 
   DO i = 1, nx (iko) 
      DO j = 1, ny (iko) 
         ikko = offz (iko - 1) + (i - 1) * ny (iko) + j 
         ikkc = offz (ikc - 1) + (i - 1) * ny (ikc) + j 
         IF (z (ikko) .ne. - 9999.0.and.z (ikko) .ne. - 9999.0.and.z (  &
         ikkc) .ne. - 9999.0.and.z (ikkc) .ne. - 9999.0) then           
            wi = r_wichtung (z (ikko), 1.0, iweight, z(ikc), bck_k) 
            wtot = wtot + wi 
            e = e+wi * z (ikko) 
            ee = ee+wi * z (ikko) **2 
            c = c + wi * z (ikkc) 
            cc = cc + wi * z (ikkc) **2 
            ce = ce+wi * z (ikkc) * z (ikko) 
         ENDIF 
      ENDDO 
   ENDDO 
ENDIF 
!                                                                       
IF (lskal.and.lback) then 
   sk = (wtot * ce-e * c) / (wtot * cc - c * c) 
   ba = (e-sk * c) / wtot 
ELSEIF (lskal.and..not.lback) then 
   sk = ce / cc 
   ba = 0.0 
ELSEIF (.not.lskal.and.lback) then 
   sk = 1.0 
   ba = (e-c) / wtot 
ELSE 
   sk = 1.0 
   ba = 0.0 
ENDIF 
!                                                                       
WRITE (output_io, 1000) sk, ba 
!                                                                       
!------ Modifying data set ikc                                          
!                                                                       
IF (.not.lni (iko) .and..not.lni (ikc) ) then 
   DO ip = 1, lenc (ikc) 
      ikkc = offxy (ikc - 1) + ip 
      y (ikkc) = sk * y (ikkc) + ba 
   ENDDO 
!                                                                       
   CALL get_extrema_xy (x, ikc, lenc (ikc), xmin, xmax) 
   CALL get_extrema_xy (y, ikc, lenc (ikc), ymin, ymax) 
ELSEIF (lni (iko) .and.lni (ikc) ) then 
   DO i = 1, nx (ikc) 
      DO j = 1, ny (ikc) 
         ikkc = offz (ikc - 1) + (i - 1) * ny (ikc) + j 
         z (ikkc) = sk * z (ikkc) + ba 
      ENDDO 
   ENDDO 
!                                                                       
   CALL get_extrema_z (z, ikc, nx (ikc), ny (ikc), zmin, zmax) 
ENDIF 
CALL show_data (ikc) 
res_para (0) = 2 
res_para (1) = sk 
res_para (2) = ba 
!                                                                       
 1000 FORMAT     (' ------ > Scale = ',g12.5,' Offset = ',g12.5,/) 
!                                                                       
END SUBROUTINE do_match                       
!
!*****7*****************************************************************
!
SUBROUTINE do_match_mix (ianz, cpara, lpara, maxw) 
!+                                                                      
!     Calculates fraction to match x*data_1 + (1-x)*data_2              
!     to given data set.                                                
!-                                                                      
USE ber_params_mod
USE errlist_mod 
USE get_params_mod
USE prompt_mod 
USE kuplot_config 
USE kuplot_mod 
use kuplot_extrema_mod
use kuplot_show_mod
USE lib_length
USE precision_mod
!                                                                       
IMPLICIT none 
!                                                                       
integer         , intent(inout) :: ianz
integer         , intent(in)    :: MAXW
CHARACTER(len=*), dimension(MAXW), intent(inout) :: cpara ! (maxw) 
INTEGER         , dimension(MAXW), intent(inout) :: lpara (maxw)
!
REAL(KIND=PREC_DP), dimension(MAXW) :: werte ! (maxw) 
!                                                                       
REAL(kind=PREC_DP) :: a, b, aa, bb, ab, ao, bo, fra 
INTEGER :: ip, iko, ik1, ik2, iik, ii1, ii2, iio 
!                                                                       
!                                                                       
!------ space left for new data set ??                                  
!                                                                       
IF (iz.gt.maxkurvtot) then 
   ier_num = - 1 
   ier_typ = ER_APPL 
   RETURN 
ENDIF 
!                                                                       
CALL del_params (1, ianz, cpara, lpara, maxw) 
CALL ber_params (ianz, cpara, lpara, werte, maxw) 
IF (ier_num.ne.0) return 
!                                                                       
IF (ianz.ne.3) then 
   ier_num = - 6 
   ier_typ = ER_COMM 
   RETURN 
ENDIF 
!                                                                       
iko = nint (werte (1) ) 
ik1 = nint (werte (2) ) 
ik2 = nint (werte (3) ) 
!                                                                       
IF (iko.lt.1.or.iko.gt.iz.or.lni (iko) .or.                          &
    ik1.lt.1.or.ik1.gt.iz.or.lni (ik1) .or.                           &
    ik2.lt.1.or.ik2.gt.iz.or.lni (ik2) .or.                           &
    iko.eq.ik1.or.iko.eq.ik2.or.ik1.eq.ik2) then                  
   ier_num = - 4 
   ier_typ = ER_APPL 
   RETURN 
ENDIF 
!                                                                       
IF (lenc (iko) .ne.lenc (ik1) .or.lenc (iko) .ne.lenc (ik2) .or.lenc ( &
      ik1) .ne.lenc (ik2) .or.xmin (ik1) .ne.xmin (iko) .or.xmin (ik2)   &
      .ne.xmin (iko) .or.xmin (ik1) .ne.xmin (ik2) .or.xmax (ik1)       &
      .ne.xmax (iko) .or.xmax (ik2) .ne.xmax (iko) .or.xmax (ik1)       &
      .ne.xmax (ik2) ) then                                             
   ier_num = - 44 
   ier_typ = ER_APPL 
   RETURN 
ENDIF 
!                                                                       
!------ Now we calculate mixing fraction                                
!                                                                       
a = 0.0 
b = 0.0 
aa = 0.0 
ao = 0.0 
bo = 0.0 
ab = 0.0 
bb = 0.0 
!                                                                       
DO ip = 1, lenc (iko) 
   iik = offxy (iz - 1) + ip 
   ii1 = offxy (ik1 - 1) + ip 
   ii2 = offxy (ik2 - 1) + ip 
   iio = offxy (iko - 1) + ip 
!                                                                       
   a = a + y (ii1) 
   b = b + y (ii2) 
   aa = aa + y (ii1) **2 
   bb = bb + y (ii2) **2 
   ao = ao + y (ii1) * y (iio) 
   bo = bo + y (ii2) * y (iio) 
   ab = ab + y (ii1) * y (ii2) 
ENDDO 
!                                                                       
fra = (bb + ao - ab - bo) / (aa - 2.0 * ab + bb) 
!                                                                       
WRITE (output_io, 1000) fra, fname (ik1) (1:len_str (fname (ik1))), &
       (1.0 - fra), fname (ik2) (1:len_str (fname (ik2) ) )           
!                                                                       
!------ Create new data set                                             
!                                                                       
DO ip = 1, lenc (iko) 
   iik = offxy (iz - 1) + ip 
   ii1 = offxy (ik1 - 1) + ip 
   ii2 = offxy (ik2 - 1) + ip 
   iio = offxy (iko - 1) + ip 
   x (iik) = x (iio) 
   y (iik) = fra * y (ii1) + (1.0 - fra) * y (ii2) 
ENDDO 
!                                                                       
lenc (iz) = lenc (iko) 
fform (iz) = fform (iko) 
fname (iz) = 'mix.dat' 
offxy (iz) = offxy (iz - 1) + lenc (iz) 
iz = iz + 1 
CALL get_extrema_xy (x, iz - 1, lenc (iz), xmin, xmax) 
CALL get_extrema_xy (y, iz - 1, lenc (iz), ymin, ymax) 
CALL show_data (iz - 1) 
!                                                                       
 1000 FORMAT     (' ------ > Mix = ',f12.5,' * ',a,' + ',               &
     &                              f12.5,' * ',a/)                     
!                                                                       
END SUBROUTINE do_match_mix                   
!
!*****7*****************************************************************
!
SUBROUTINE do_rvalue(zeile, lp, lout) 
!                                                                       
!     Calculates the residual of two data sets.                         
!                                                                       
!     Version : 4.1                                                     
!     Date    : 04. August 2004                                         
!                                                                       
!     Author  : R.B. Neder (reinhar.neder@krist.uni-erlangen.de)        
!                                                                       
!*****7*****************************************************************
!                                                                       
USE ber_params_mod
USE errlist_mod 
USE get_params_mod
USE param_mod 
USE prompt_mod 
USE kuplot_config 
USE kuplot_mod 
USE take_param_mod
USE string_convert_mod
USE precision_mod
USE str_comp_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER maxw 
PARAMETER (maxw = 5) 
!                                                                       
CHARACTER (LEN=*), INTENT(INOUT) :: zeile 
INTEGER          , INTENT(INOUT) :: lp
LOGICAL          , INTENT(IN)    :: lout
!
CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
REAL(KIND=PREC_DP)         :: werte (maxw) 
INTEGER                    :: lpara (maxw)
INTEGER :: ianz, ik, il 
INTEGER :: iianz 
INTEGER :: iweight 
REAL(kind=PREC_DP) :: rval 
REAL(kind=PREC_DP) :: wrval 
real(kind=PREC_DP) :: bck_k
!                                                                       
integer, PARAMETER :: W_ONE = 0 
integer, PARAMETER :: W_SQUA = 1 
integer, PARAMETER :: W_SQRT = 2 
integer, PARAMETER :: W_INV = 3 
integer, PARAMETER :: W_LOG = 4 
integer, PARAMETER :: W_ISQ = 5 
integer, PARAMETER :: W_LIN = 6
integer, PARAMETER :: W_DAT = 7
integer, PARAMETER :: W_BCK = 8 
!                                                                       
!                                                                       
INTEGER, PARAMETER :: NOPTIONAL = 1
CHARACTER(LEN=   7), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=PREC_STRING), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent  !opt. para present
REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 1 ! Number of values to calculate 
INTEGER :: irvalue
INTEGER :: i
!
DATA oname  / 'partial'/
DATA loname /  7        /
!
opara  =  (/ '0.0000'   /) ! Always provide fresh default values
lopara =  (/  6         /)
owerte =  (/  0.00      /)
bck_k = 0.0
!
CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
IF (ier_num.ne.0) return 
!                                                                       
IF (ianz.lt.2) then 
   ier_num = - 6 
   ier_typ = ER_COMM 
   RETURN 
ENDIF 
!                                                                       
CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
IF (ier_num.ne.0) RETURN 
irvalue = NINT(owerte(1))
!
iianz = 2 
CALL ber_params (iianz, cpara, lpara, werte, maxw) 
CALL del_params (2, ianz, cpara, lpara, maxw) 
!                                                                       
IF (ier_num.ne.0) return 
!                                                                       
ik = nint (werte (1) ) 
il = nint (werte (2) ) 
!                                                                       
!      CALL del_params (1, ianz, cpara, lpara, maxw) 
CALL do_low(cpara(1)(1:lpara(1)))
iweight = W_ONE 
IF(ianz.eq.0.or.str_comp(cpara(1), 'one', 1, lpara(1), 3) ) then
   iweight = W_ONE 
ELSEIF (str_comp (cpara (1) , 'squa', 3, lpara (1) , 4) ) then 
   iweight = W_SQUA 
ELSEIF (str_comp (cpara (1) , 'sqrt', 3, lpara (1) , 4) ) then 
   iweight = W_SQRT 
ELSEIF (str_comp (cpara (1) , 'inv', 1, lpara (1) , 3) ) then 
   iweight = W_INV 
ELSEIF (str_comp (cpara (1) , 'log', 2, lpara (1) , 3) ) then 
   iweight = W_LOG 
ELSEIF (str_comp (cpara (1) , 'isq', 1, lpara (1) , 3) ) then 
   iweight = W_ISQ 
ELSEIF (str_comp (cpara (1) , 'lin', 2, lpara (1) , 3) ) then 
   iweight = W_LIN 
ELSEIF (str_comp (cpara (1) , 'dat', 1, lpara (1) , 3) ) then 
   iweight = W_DAT 
ELSEIF (str_comp (cpara (1) , 'bck', 1, lpara (1) , 3) ) then 
   CALL del_params (1, ianz, cpara, lpara, maxw) 
   iianz = 1 
   CALL ber_params (iianz, cpara, lpara, werte, maxw) 
   bck_k = werte(1)
   iweight = W_BCK 
ELSE
   ier_num = -27
   ier_typ = ER_APPL 
   RETURN 
ENDIF 
!                                                                       
IF (ik.gt.0.and.ik.lt.iz.and.il.gt.0.and.il.lt.iz) then 
   IF (lni (ik) .and.lni (il) ) then 
      IF (nx (ik) .eq.nx (il) .and.ny (ik) .eq.ny (il) ) then 
         CALL rvalue_z (z, nx (ik), ny (ik), ik, il, iweight,     &
         rval, wrval, bck_k)                                             
         res_para (0) = 2 
         res_para (1) = rval 
         res_para (2) = wrval 
         rvalues  (1, irvalue) = rval 
         rvalues  (2, irvalue) = wrval 
         rvalue_yes   = .true.
      ELSE 
         ier_num = - 23 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
   ELSEIF (.not.lni (ik) .and..not.lni (il) ) then 
      IF (lenc (ik) .eq.lenc (il) ) then 
         CALL rvalue_y (y, dy, lenc(ik), ik, il,        &
         iweight, rval, wrval, bck_k)                                    
         res_para (0) = 2 
         res_para (1) = rval 
         res_para (2) = wrval 
         rvalues  (1, irvalue) = rval 
         rvalues  (2, irvalue) = wrval 
         rvalue_yes   = .true.
      ELSE 
         ier_num = - 23 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
   ELSE 
      ier_num = - 6 
      ier_typ = ER_COMM 
      RETURN 
   ENDIF 
ELSE 
   ier_num = - 6 
   ier_typ = ER_COMM 
   RETURN 
ENDIF 
!
IF(irvalue > 1) THEN
   rvalues(1,0) = 0.0
   rvalues(2,0) = 0.0
   nrvalues = MAX(irvalue,nrvalues)
   DO i=1, MAX(irvalue,nrvalues)
      rvalues(1,0) = rvalues(1,0) + rvalues(1,i)
      rvalues(2,0) = rvalues(2,0) + rvalues(2,i)
   ENDDO
   rvalues(1,0) = rvalues(1,0) /nrvalues
   rvalues(2,0) = rvalues(2,0) /nrvalues
ELSE
   nrvalues = MAX(1, nrvalues)
ENDIF
!                                                                       
IF(lout) THEN
   WRITE (output_io, 1000) rval 
   WRITE (output_io, 2000) wrval 
ENDIF
!                                                                       
 1000 FORMAT    (' R-value          : ',f8.4) 
 2000 FORMAT    (' weighted R-value : ',f8.4) 
!
END SUBROUTINE do_rvalue                      
!
!*****7*****************************************************************
!
SUBROUTINE rvalue_z (a, nxx, nyy, ik, il, iweight, rval, wrval, bck_k) 
!+                                                                      
!     berechnungen fuer z feld                                          
!-                                                                      
USE errlist_mod 
USE kuplot_config 
USE kuplot_mod 
use kuplot_wichtung_mod
!
use precision_mod
!                                                                       
IMPLICIT none 
!
REAL(kind=PREC_DP), dimension(MAXARRAY), intent(in) :: a
integer                                , intent(in) :: nxx
integer                                , intent(in) :: nyy
integer                                , intent(in) :: ik
integer                                , intent(in) :: il
integer                                , intent(in) :: iweight
REAL(kind=PREC_DP)                     , intent(out) :: rval 
REAL(kind=PREC_DP)                     , intent(out) :: wrval 
REAL(kind=PREC_DP)                     , intent(in) :: bck_k
!                                                                       
INTEGER :: i, j, ikk, ikl 
REAL(kind=PREC_DP) :: sumrz, sumrn 
REAL(kind=PREC_DP) :: sumwrz, sumwrn 
REAL(kind=PREC_DP) :: wght 
!                                                                       
!     REAL r_wichtung 
!                                                                       
sumrz = 0.0 
sumrn = 0.0 
sumwrz = 0.0 
sumwrn = 0.0 
DO i = 1, nxx 
   DO j = 1, nyy 
      ikk = offz (ik - 1) + (i - 1) * ny (ik) + j 
      ikl = offz (il - 1) + (i - 1) * ny (il) + j 
      IF (a (ikl) .ne. - 9999.0.and.a (ikk) .ne. - 9999.0) then 
         sumrz = sumrz + abs (a (ikk) - a (ikl) ) 
         sumrn = sumrn + abs (a (ikk) ) 
         wght = r_wichtung (a (ikk), 1.0, iweight, a(ikl), bck_k) 
         sumwrz = sumwrz + wght * (a (ikk) - a (ikl) ) **2 
         sumwrn = sumwrn + wght * (a (ikk) ) **2 
      ENDIF 
   ENDDO 
ENDDO 
!                                                                       
rval = sumrz / sumrn 
wrval = sqrt (sumwrz / sumwrn) 
!                                                                       
END SUBROUTINE rvalue_z                       
!
!*****7*****************************************************************
!
SUBROUTINE rvalue_y(a, da, klen, ik, il, iweight, rval, wrval, bck_k)
!+                                                                      
!     berechnungen fuer y feld                                          
!-                                                                      
USE errlist_mod 
USE kuplot_config 
USE kuplot_mod 
use kuplot_wichtung_mod
!                                                                       
IMPLICIT none 
!
REAL(kind=PREC_DP), dimension(MAXARRAY), intent(in)  :: a
REAL(kind=PREC_DP), dimension(MAXARRAY), intent(in)  :: da
INTEGER                                , intent(in)  :: klen
INTEGER                                , intent(in)  :: ik
INTEGER                                , intent(in)  :: il
INTEGER                                , intent(in)  :: iweight
REAL(kind=PREC_DP)                     , intent(out) :: rval 
REAL(kind=PREC_DP)                     , intent(out) :: wrval 
REAL(kind=PREC_DP)                     , intent(in)  :: bck_k
!                                                                       
INTEGER :: i, ikk, ikl 
REAL(kind=PREC_DP) :: sumrz, sumrn 
REAL(kind=PREC_DP) :: sumwrz, sumwrn 
REAL(kind=PREC_DP) :: wght 
!                                                                       
!      REAL r_wichtung 
!                                                                       
sumrz = 0.0 
sumrn = 0.0 
sumwrz = 0.0 
sumwrn = 0.0 
DO i = 1, klen 
   ikk    = offxy (ik - 1) + i 
   ikl    = offxy (il - 1) + i 
   sumrz  = sumrz + abs (a (ikk) - a (ikl) ) **2
   sumrn  = sumrn + abs (a (ikk) ) **2
   wght   = r_wichtung (a (ikk), da (ikk), iweight, a(ikl), bck_k) 
   sumwrz = sumwrz + wght * (a (ikk) - a (ikl) ) **2 
   sumwrn = sumwrn + wght * (a (ikk) ) **2 
ENDDO 
!                                                                       
rval  = sqrt (sumrz / sumrn) 
wrval = sqrt (sumwrz / sumwrn) 
                                                                        
END SUBROUTINE rvalue_y                       
!
!*****7*****************************************************************
!
!REAL(kind=PREC_DP) function r_wichtung (z, dz, iweight, zc, bck_k) 
!!                                                                       
!!     calculates the weight                                             
!!                                                                       
!use precision_mod
!!
!IMPLICIT none 
!!                                                                       
!REAL(kind=PREC_DP), intent(in) :: z
!REAL(kind=PREC_SP), intent(in) :: dz 
!INTEGER           , intent(in) :: iweight 
!REAL(kind=PREC_SP), intent(in) :: zc 
!REAL(kind=PREC_DP), intent(in) :: bck_k
!!                                                                       
!INTEGER, parameter :: W_ONE  = 0
!INTEGER, parameter :: W_SQUA = 1
!INTEGER, parameter :: W_SQRT = 2
!INTEGER, parameter :: W_INV  = 3
!INTEGER, parameter :: W_LOG  = 4
!INTEGER, parameter :: W_ISQ  = 5
!INTEGER, parameter :: W_LIN  = 6
!INTEGER, parameter :: W_DAT  = 7
!INTEGER , parameter :: W_BCK = 8          ! Currently not used
!!     PARAMETER (W_ONE = 0) 
!!     PARAMETER (W_SQUA = 1) 
!!     PARAMETER (W_SQRT = 2) 
!!     PARAMETER (W_INV = 3) 
!!     PARAMETER (W_LOG = 4) 
!!     PARAMETER (W_ISQ = 5) 
!!     PARAMETER (W_LIN = 6) 
!!     PARAMETER (W_DAT = 7) 
!!     PARAMETER (W_BCK = 8) 
!!                                                                       
!         
!r_wichtung = 0.0 
!IF (iweight.eq.W_ONE) then 
!   r_wichtung = 1.0 
!ELSEIF (iweight.eq.W_SQUA) then 
!   r_wichtung = z**2 
!ELSEIF (iweight.eq.W_SQRT) then 
!   IF (z.ge.0) then 
!      r_wichtung = sqrt (z) 
!   ELSE 
!      r_wichtung = 0.0 
!   ENDIF 
!ELSEIF (iweight.eq.W_INV) then 
!   IF (z.ne.0) then 
!      r_wichtung = 1.0 / abs (z) 
!   ELSE 
!      r_wichtung = 0.0 
!   ENDIF 
!ELSEIF (iweight.eq.W_LOG) then 
!   IF (z.gt.0) then 
!      r_wichtung = log (z) 
!   ELSE 
!      r_wichtung = 0.0 
!   ENDIF 
!ELSEIF (iweight.eq.W_ISQ) then 
!   IF (z.gt.0) then 
!      r_wichtung = 1. / sqrt (z) 
!   ELSE 
!      r_wichtung = 0.0 
!   ENDIF 
!ELSEIF (iweight.eq.W_LIN) then 
!   r_wichtung = z 
!ELSEIF (iweight.eq.W_DAT) then 
!   IF (dz.ne.0) then 
!      r_wichtung = 1.0 / dz**2
!   ELSE 
!      r_wichtung = 1.0
!   ENDIF 
!elseif(iweight==W_BCK) then
!   if((z-zc)>0.0) then
!      r_wichtung = dz*(0.8*bck_k)
!   else
!      r_wichtung = dz/(0.8*bck_k)
!   endif
!!           r_wichtung = exp(-abs(bck_k)*(z-zc))
!ENDIF 
!!                                                                       
!END FUNCTION r_wichtung                       
!
end module kuplot_calc_mod
