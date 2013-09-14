!
! routines specific for MAC_OS
!
LOGICAL FUNCTION IS_IOSTAT_END ( status_flag )
!
! Is actually FORTRAN2003 standard, should not be needed???
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: status_flag
!
IS_IOSTAT_END = status_flag == -1
!
END FUNCTION IS_IOSTAT_END
