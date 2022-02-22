module kuplot_low_mod
!
! Low level functions
!
contains
!
!******7****************************************************************
!
LOGICAL function n_in_f (ini) 
!                                                                       
!     Checks for 2D files in current frame                              
!                                                                       
USE kuplot_config 
USE kuplot_mod 
!                                                                       
IMPLICIT none 
!                                                                       
integer, intent(out) :: ini
!
INTEGER :: i
!                                                                       
ini = 0
n_in_f = .false. 
DO i = 1, iz - 1 
   n_in_f = n_in_f.or. (k_in_f (i) .and.lni (i) ) 
   IF(n_in_f) ini = i 
ENDDO 
!                                                                       
END FUNCTION n_in_f                           
!
!*****7*****************************************************************
!
LOGICAL function k_in_f (ik) 
!+                                                                      
!     Check if data set ik should be plotted in actual frame            
!-                                                                      
USE kuplot_config 
USE kuplot_mod 
!                                                                       
IMPLICIT none 
!                                                                       
integer, intent(in) :: ik
!
INTEGER :: i
!                                                                       
k_in_f = .false. 
DO i = 1, maxkurvtot 
   k_in_f = k_in_f.or. (infra (iwin, iframe, i) == ik) 
ENDDO 
!                                                                       
END FUNCTION k_in_f                           
!
!*****7*****************************************************************
!
LOGICAL function inrect (ex, ey, x, y) 
!-                                                                      
!       Checks if point x,y is within rectangle ex(2),ey(2).            
!+                                                                      
!
implicit none
!
REAL, intent(in) :: x
REAL, intent(in) :: y
REAL,  dimension(2), intent(in) :: ex ! (2)
REAL,  dimension(2), intent(in) :: ey ! (2) 
!                                                                       
inrect = (x >= ex(1) .and. x <= ex(2)) .and.    &
         (y >= ey(1) .and. y <= ey(2)        )                                                
!
END FUNCTION inrect                           
!
!***********************************************************************
!
end module kuplot_low_mod
