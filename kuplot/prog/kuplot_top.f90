SUBROUTINE kuplot_top(zeile)
!
!  Specific KUPLOT Version of a branch subroutine
!  Call DISCUS via system
!
USE errlist_mod
!
IMPLICIT NONE
!
CHARACTER (LEN=*), INTENT(IN) :: zeile
!
ier_num = -7
ier_typ = ER_COMM
!
END SUBROUTINE kuplot_top
