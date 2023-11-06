MODULE crystal_task_mod
!-
!  Collection of basic tasks for the crystal
!+
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE crystal_calc_mass
!
!  Calculate crystal mass
!
USE crystal_mod
USE element_data_mod
!
USE precision_mod
!
IMPLICIT NONE
!
INTEGER :: iscat
INTEGER :: i
!
cr_mass  = 0.0E0
cr_nreal = 0.0E0
!
!     Prepare and calculate average atom numbers
!
cr_niscat = 0
DO i=1,cr_natoms
   cr_niscat(cr_iscat(i,1)) = cr_niscat(cr_iscat(i,1)) + 1
ENDDO
cr_nreal = SUM(NINT(cr_niscat(1:cr_nscat)*cr_occ(1:cr_nscat)))  ! Add real atom numbers
!
DO iscat=1,cr_nscat
   cr_mass = cr_mass + cr_niscat(iscat)*cr_occ(iscat)*get_mass(cr_at_lis(iscat))
ENDDO
!
END SUBROUTINE crystal_calc_mass
!
!*******************************************************************************
!

END MODULE crystal_task_mod
