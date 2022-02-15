module kuplot_low_mod
!
! Low level functions
!
contains
!
!******7****************************************************************
      LOGICAL function n_in_f (ini) 
!                                                                       
!     Checks for 2D files in current frame                              
!                                                                       
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER i, ini 
!      LOGICAL k_in_f 
!                                                                       
      n_in_f = .false. 
      DO i = 1, iz - 1 
      n_in_f = n_in_f.or. (k_in_f (i) .and.lni (i) ) 
      IF (n_in_f) ini = i 
      ENDDO 
!                                                                       
      END FUNCTION n_in_f                           
!*****7*****************************************************************
      LOGICAL function k_in_f (ik) 
!+                                                                      
!     Check if data set ik should be plotted in actual frame            
!-                                                                      
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER ik, i 
!                                                                       
      k_in_f = .false. 
      DO i = 1, maxkurvtot 
      k_in_f = k_in_f.or. (infra (iwin, iframe, i) .eq.ik) 
      ENDDO 
!                                                                       
      END FUNCTION k_in_f                           
!*****7*****************************************************************
      LOGICAL function inrect (ex, ey, x, y) 
!-                                                                      
!       Checks if point x,y is within rectangle ex(2),ey(2).            
!+                                                                      
      REAL x, y, ex (2), ey (2) 
!                                                                       
      inrect = (x.ge.ex (1) .and.x.le.ex (2) ) .and. (y.ge.ey (1)       &
      .and.y.le.ey (2) )                                                
      END FUNCTION inrect                           
!***********************************************************************
!
end module kuplot_low_mod
