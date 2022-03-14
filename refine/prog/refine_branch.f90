SUBROUTINE refine_branch(zeile, length, lreset, lloop)
!
!  Specific REFINE Version of a branch subroutine
!  Call DISCUS/KUPLOT via system
!  Currently this gives an error message in a stand alone 
!  program
!
USE errlist_mod
!
IMPLICIT NONE
!
CHARACTER (LEN=*), INTENT(IN) :: zeile
INTEGER          , INTENT(IN) :: length
LOGICAL          , INTENT(IN) :: lreset
integer          , INTENT(IN) :: lloop
!
ier_num = -7
ier_typ = ER_COMM
!
END SUBROUTINE refine_branch
