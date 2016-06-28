MODULE operating_mod
!
CONTAINS
!
SUBROUTINE operating_exit
!
USE prompt_mod
IMPLICIT NONE
!
! Currently no need for specific exit
!IF(.NOT. lstandalone .OR. pname_CAP == 'KUPLOT') THEN
!   WRITE(*,*) ' '
!   WRITE(*,*) pname_cap,' is finished'
!   WRITE(*,*) 'Close PGPLOT window '
!   WRITE(*,*) 'If this is the primary window close pgplot_server as well'
!   WRITE(*,*) 'Finally close this window as well'
!   WRITE(*,*) ' '
!ENDIF
!
END SUBROUTINE operating_exit
END MODULE operating_mod
