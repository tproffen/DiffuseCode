MODULE refine_constraint_mod
!
!  Rotines related to parameter constraints
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE refine_constrain_auto
!-
! Apply and set automatic constraints for special parameter names
! These are:
! P_lat...
! P_biso...
! P_dia...
!+
USE refine_params_mod
!
IMPLICIT NONE
!
INTEGER :: i    ! Dummy loop variable
!
DO i=1,refine_par_n
   IF(refine_params(i)(1:5)=='P_lat'  .OR.                                      &
      refine_params(i)(1:6)=='P_biso' .OR.                                      &
      refine_params(i)(1:5)=='P_dia'                                            &
                                           ) THEN
      IF(refine_range(i,1)<refine_range(i,2)) THEN  ! User did set a range
         refine_range(i,1) = MAX(0.0, refine_range(i,1))
      ELSE
         refine_range(i,1) = 0.0
         refine_range(i,2) = HUGE(0.0)
      ENDIF
   ENDIF
ENDDO
!
END SUBROUTINE refine_constrain_auto
!
!*******************************************************************************
!
END MODULE refine_constraint_mod
