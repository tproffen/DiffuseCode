MODULE kuplot_3dm_draw
!
IMPLICIT NONE
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE kuplot_draw_3d_static(ik)
!
! Draws a 3D rendering of a bitmap
!
USE kuplot_config
USE kuplot_mod
USE kuplot_3dm_mod
!
IMPLICIT NONE
!
INTEGER, INTENT(in) :: ik
!
!REAL :: zzmin, zzmax    ! min max z values as of hlin command
!
END SUBROUTINE kuplot_draw_3d_static
!
!*******************************************************************************
!
END MODULE kuplot_3dm_draw
