MODULE koordinate_mod
!
PRIVATE
PUBLIC koor_shear, koor_log
!
INTERFACE koor_shear
   MODULE PROCEDURE koor_shear_s, koor_shear_1
END INTERFACE koor_shear
!
INTERFACE koor_log
   MODULE PROCEDURE koor_log_s, koor_log_1
END INTERFACE koor_log
!
CONTAINS
!***********************************************************************
      SUBROUTINE koor_shear_s (npkt, xa, ya) 
!+                                                                      
!     transformation beim shearen Scalar Version
!-                                                                      
      IMPLICIT none 
!                                                                       
      INTEGER npkt, i 
      REAL xa       , ya
!                                                                       
      include'config.inc' 
      include'kuplot.inc' 
      include'wink.inc' 
!                                                                       
      IF (shear (iwin, iframe) .ne.90.0.and.sfl (iwin, iframe) ) then 
!        DO i = 1, npkt 
         IF (ya     .gt.pey (iwin, iframe, 1) ) then 
            xa     = xa     + (ya     - pey (iwin, iframe, 1) ) * yskal &
            (iwin, iframe) / tan (rad * shear (iwin, iframe) )          
         ELSE 
            xa     = xa     
         ENDIF 
         ya     = ya     
!        ENDDO 
      ENDIF 
!                                                                       
      END SUBROUTINE koor_shear_s
!***********************************************************************
      SUBROUTINE koor_shear_1 (npkt, xa, ya) 
!+                                                                      
!     transformation beim shearen                                       
!-                                                                      
      IMPLICIT none 
!                                                                       
      INTEGER npkt, i 
      REAL xa (npkt), ya (npkt) 
!                                                                       
      include'config.inc' 
      include'kuplot.inc' 
      include'wink.inc' 
!                                                                       
      IF (shear (iwin, iframe) .ne.90.0.and.sfl (iwin, iframe) ) then 
         DO i = 1, npkt 
         IF (ya (i) .gt.pey (iwin, iframe, 1) ) then 
            xa (i) = xa (i) + (ya (i) - pey (iwin, iframe, 1) ) * yskal &
            (iwin, iframe) / tan (rad * shear (iwin, iframe) )          
         ELSE 
            xa (i) = xa (i) 
         ENDIF 
         ya (i) = ya (i) 
         ENDDO 
      ENDIF 
!                                                                       
      END SUBROUTINE koor_shear_1
!***********************************************************************
      SUBROUTINE koor_log_s (npkt, xa, ya) 
!+                                                                      
!     transformation for LOG axes , scalar version
!-                                                                      
      IMPLICIT none 
!                                                                       
      INTEGER npkt, i 
      REAL xa       , ya
      REAL log10 
!                                                                       
      include'config.inc' 
      include'kuplot.inc' 
      include'wink.inc' 
!                                                                       
      log10 = log (10.0) 
!                                                                       
      IF (lachse (iwin, iframe, 1) ) then 
!        DO i = 1, npkt 
         IF (xa     .ne.0.0) then 
            xa     = log (abs (xa     ) ) / log10 
         ELSE 
            xa     = - 9999. 
         ENDIF 
!        ENDDO 
      ENDIF 
!                                                                       
      IF (lachse (iwin, iframe, 2) ) then 
!        DO i = 1, npkt 
         IF (ya     .ne.0.0) then 
            ya     = log (abs (ya     ) ) / log10 
         ELSE 
            ya     = - 9999. 
         ENDIF 
!        ENDDO 
      ENDIF 
!                                                                       
      END SUBROUTINE koor_log_s
!***********************************************************************
      SUBROUTINE koor_log_1 (npkt, xa, ya) 
!+                                                                      
!     transformation for LOG axes , 1D version
!-                                                                      
      IMPLICIT none 
!                                                                       
      INTEGER npkt, i 
      REAL xa (npkt), ya (npkt) 
      REAL log10 
!                                                                       
      include'config.inc' 
      include'kuplot.inc' 
      include'wink.inc' 
!                                                                       
      log10 = log (10.0) 
!                                                                       
      IF (lachse (iwin, iframe, 1) ) then 
         DO i = 1, npkt 
         IF (xa (i) .ne.0.0) then 
            xa (i) = log (abs (xa (i) ) ) / log10 
         ELSE 
            xa (i) = - 9999. 
         ENDIF 
         ENDDO 
      ENDIF 
!                                                                       
      IF (lachse (iwin, iframe, 2) ) then 
         DO i = 1, npkt 
         IF (ya (i) .ne.0.0) then 
            ya (i) = log (abs (ya (i) ) ) / log10 
         ELSE 
            ya (i) = - 9999. 
         ENDIF 
         ENDDO 
      ENDIF 
!                                                                       
      END SUBROUTINE koor_log_1                       
END MODULE koordinate_mod
