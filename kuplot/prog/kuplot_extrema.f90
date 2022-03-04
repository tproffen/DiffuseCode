module kuplot_extrema_mod
!
!  Routines to find extrema in data sets
!+
!
contains
!
!*******************************************************************************
!
SUBROUTINE get_extrema 
!+                                                                      
!     get max/min values of files                                       
!-                                                                      
USE kuplot_config 
USE kuplot_mod 
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER :: ik 
!                                                                       
DO ik = 1, iz - 1 
   IF (lni (ik) ) then 
      CALL get_extrema_xy (x, ik, nx (ik), xmin, xmax) 
      CALL get_extrema_xy (y, ik, ny (ik), ymin, ymax) 
      CALL get_extrema_z (z, ik, nx (ik), ny (ik), zmin, zmax) 
   ELSE 
      CALL get_extrema_xy (x, ik, lenc (ik), xmin, xmax) 
      CALL get_extrema_xy (y, ik, lenc (ik), ymin, ymax) 
   ENDIF 
ENDDO 
!                                                                       
END SUBROUTINE get_extrema                    
!
!*******************************************************************************
!
SUBROUTINE get_extrema_xy_local (i, ymi, yma) 
!+                                                                      
!     extrema der kurve ik local bestimmen                              
!-                                                                      
USE kuplot_config 
USE kuplot_mod 
!
use precision_mod
!                                                                       
IMPLICIT none 
!                                                                       
integer, intent(in) :: i
REAL(kind=PREC_DP), intent(out) :: ymi
REAL(kind=PREC_DP), intent(out) :: yma 
!
INTEGER :: ip 
!                                                                       
yma = - 1e38 
ymi = 1e38 
DO ip = 1, lenc(i) 
   IF (x (offxy (i - 1) + ip) .ge.ex (iwin, iframe, 1) .and.x (offxy &
      (i - 1) + ip) .le.ex (iwin, iframe, 2) ) then                     
      yma = max (yma, real(y (offxy (i - 1) + ip), kind=PREC_DP) ) 
      ymi = min (ymi, real(y (offxy (i - 1) + ip), kind=PREC_DP) ) 
   ENDIF 
ENDDO 
!                                                                       
END SUBROUTINE get_extrema_xy_local           
!
!**7*****************************************************************   
!
SUBROUTINE get_extrema_xy (a, ik, ilen, amin, amax) 
!+                                                                      
!     extrema der kurve ik bestimmen                                    
!-                                                                      
USE kuplot_config 
USE kuplot_mod 
!                                                                       
IMPLICIT none 
!                                                                       
REAL    , intent(in) :: a (maxarray) 
integer , intent(in) :: ik
integer , intent(in) :: ilen
REAL(kind=PREC_DP), dimension(MAXKURVTOT), intent(inout) :: amax (maxkurvtot), amin (maxkurvtot) 
!
INTEGER :: ip
!                                                                       
amax (ik) = a (offxy (ik - 1) + 1) 
amin (ik) = a (offxy (ik - 1) + 1) 
DO ip = 2, ilen 
   amax (ik) = max (amax (ik), real(a (offxy (ik - 1) + ip), kind=PREC_DP) ) 
   amin (ik) = min (amin (ik), real(a (offxy (ik - 1) + ip), kind=PREC_DP) ) 
ENDDO 
!                                                                       
END SUBROUTINE get_extrema_xy                 
!
!**7*****************************************************************   
      SUBROUTINE get_extrema_z (a, ik, nxx, nyy, amin, amax) 
!+                                                                      
!     extrema der kurve ik bestimmen                                    
!-                                                                      
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
REAL , intent(in) :: a (maxarray) 
integer , intent(in) :: nxx
integer , intent(in) :: nyy
REAL , intent(out) :: amax (maxkurvtot), amin (maxkurvtot) 
!
      INTEGER :: ik, ikk, i, j 
!                                                                       
      amax (ik) = a (offz (ik - 1) + 1) 
      amin (ik) = a (offz (ik - 1) + 1) 
      IF (amin (ik) .eq. - 9999.) amin (ik) = 1e+30 
      DO i = 1, nxx 
      DO j = 2, nyy 
      ikk = offz (ik - 1) + (i - 1) * ny (ik) + j 
      amax (ik) = max (amax (ik), a (ikk) ) 
      IF (a (ikk) .ne. - 9999.) amin (ik) = min (amin (ik), a (ikk) ) 
      ENDDO 
      ENDDO 
!                                                                       
      END SUBROUTINE get_extrema_z                  
!
!**7******************************************************************* 
!
SUBROUTINE do_fmax_xy (ik, wmax, ixm, maxmax, ima) 
!+                                                                      
!     maxima bestimmen fuer xy-files                                    
!-                                                                      
USE errlist_mod 
USE kuplot_config 
USE kuplot_mod 
!                                                                       
use precision_mod
!
IMPLICIT none 
!                                                                       
integer                              , intent(in)  :: ik
integer                              , intent(in)  :: maxmax 
integer                              , intent(out) :: ima
REAL(kind=PREC_DP), dimension(MAXMAX), intent(out) :: wmax! (maxmax) 
integer           , dimension(MAXMAX), intent(out) :: ixm
!
INTEGER :: ix, jf 
LOGICAL :: max_da 
!                                                                       
ima = 1 
DO ix = 1 + ifen, lenc(ik) - ifen 
   IF (x (offxy (ik - 1) + ix) .ge.ex (iwin, iframe, 1) .and.x (     &
      offxy (ik - 1) + ix) .le.ex (iwin, iframe, 2) ) then              
      max_da = .true. 
      DO jf = - ifen+1, - 1 
         max_da = max_da.and. (y (offxy (ik - 1) + ix) .ge.    &
                               y (offxy (ik - 1) + ix + jf) )  &
                        .AND. (y (offxy (ik - 1) + ix) .gt.    &
                               y (offxy (ik - 1) + ix + jf-1) )                                              
      ENDDO 
      DO jf = 1, ifen -1
         max_da = max_da.and. (y (offxy (ik - 1) + ix) .ge.    &
                               y (offxy (ik - 1) + ix + jf) )  &
                        .AND. (y (offxy (ik - 1) + ix) .gt.    &
                               y (offxy (ik - 1) + ix + jf+1) )                                              
      ENDDO 
!                                                                       
      IF (max_da) then 
         wmax (ima) = y (offxy (ik - 1) + ix) 
         ixm (ima) = ix 
         ima = ima + 1 
         IF (ima.gt.maxmax) then 
            ier_num = - 21 
            ier_typ = ER_APPL 
            ima = maxmax 
            RETURN 
         ENDIF 
      ENDIF 
   ENDIF 
ENDDO 
ima = ima - 1 
!                                                                       
IF (ima.eq.0) then 
   ier_num = - 22 
   ier_typ = ER_APPL 
ENDIF 
!                                                                       
END SUBROUTINE do_fmax_xy                     
!
!**7******************************************************************* 
!
SUBROUTINE do_fmax_z (ik, wmax, ixm, iym, maxmax, ima) 
!+                                                                      
!     maxima bestimmen aus xyz files                                    
!-                                                                      
USE errlist_mod 
USE kuplot_config 
USE kuplot_mod 
!
use precision_mod
!                                                                       
IMPLICIT none 
!                                                                       
integer, intent(in) :: ik
INTEGER, intent(in) :: maxmax 
integer, intent(out) :: ima
REAL(kind=PREC_DP), dimension(MAXMAX), intent(out) :: wmax! (maxmax) 
integer           , dimension(MAXMAX), intent(out) :: ixm
integer           , dimension(MAXMAX), intent(out) :: iym
!
INTEGER :: ikk, iix, ix, iy, jf, kf 
LOGICAL :: max_da 
!                                                                       
      ima = 1 
      ikk = offz (ik - 1) 
      DO ix = 1 + ifen, nx (ik) - ifen 
      DO iy = 1 + ifen, ny (ik) - ifen 
      iix = ix - 1 
      IF (z (ikk + iix * ny (ik) + iy) .eq. - 9999.0) goto 20 
      max_da = .true. 
      DO jf = - ifen, ifen 
      DO kf = - ifen, ifen 
      IF (jf.ne.0.or.kf.ne.0) then 
         max_da = max_da.and. (z (ikk + iix * ny (ik) + iy) ) .ge.z (   &
         ikk + (iix + jf) * ny (ik) + iy + kf)                          
         IF (.not.max_da) goto 20 
      ENDIF 
      ENDDO 
      ENDDO 
      IF (max_da) then 
         wmax (ima) = z (ikk + iix * ny (ik) + iy) 
         ixm (ima) = ix 
         iym (ima) = iy 
         ima = ima + 1 
         IF (ima.gt.maxmax) then 
            ier_num = - 21 
            ier_typ = ER_APPL 
            ima = maxmax 
            RETURN 
         ENDIF 
      ENDIF 
   20 CONTINUE 
      ENDDO 
      ENDDO 
      ima = ima - 1 
!                                                                       
      IF (ima.eq.0) then 
         ier_num = - 22 
         ier_typ = ER_APPL 
      ENDIF 
!                                                                       
END SUBROUTINE do_fmax_z                      
!
!*******************************************************************************
!
end module kuplot_extrema_mod
