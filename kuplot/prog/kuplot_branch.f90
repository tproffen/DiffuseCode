SUBROUTINE kuplot_branch(zeile, length, lreset)
!
!  Specific KUPLOT Version of a branch subroutine
!  Call DISCUS via system
!
USE errlist_mod
!
IMPLICIT NONE
!
CHARACTER (LEN=*), INTENT(IN) :: zeile
INTEGER          , INTENT(IN) :: length
LOGICAL          , INTENT(IN) :: lreset
!
ier_num = -7
ier_typ = ER_COMM
!
END SUBROUTINE kuplot_branch
