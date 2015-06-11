SUBROUTINE diffev_branch(zeile, length)
!
!  Specific DIFFEV Version of a branch subroutine
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
!
ier_num = -7
ier_typ = ER_COMM
!
END SUBROUTINE diffev_branch
