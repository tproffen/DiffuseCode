SUBROUTINE discus_branch(zeile, length, lreset, lloop)
!
!  Specific DISCUS Version of a branch subroutine
!  Call KUPLOT via system
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
END SUBROUTINE discus_branch
