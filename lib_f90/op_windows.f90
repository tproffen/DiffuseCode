MODULE operating_mod
!
INTEGER, PARAMETER :: OP_LINUX   = 1
INTEGER, PARAMETER :: OP_MAC     = 2
INTEGER, PARAMETER :: OP_WINDOWS = 3
!
INTEGER, PARAMETER :: OP_SYSTEM = OP_WINDOWS
!
CONTAINS
!
SUBROUTINE get_mpi_path(mpi_path)
!
! Find path to 'mpiexec' command
! For Windows version this is fixed
!
IMPLICIT NONE
!
CHARACTER(LEN=1024), INTENT(OUT) :: mpi_path
!
mpi_path = '/bin'
!
END SUBROUTINE get_mpi_path
!
SUBROUTINE get_discus_path(discus_path, discus_name)
!
! Find path to 'discus_suite' command
! For Windows version this is fixed
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(OUT) :: discus_path
CHARACTER(LEN=*), INTENT(OUT) :: discus_name
!
discus_path = '/bin/'
discus_name = 'discus_suite_parallel.exe'
!
END SUBROUTINE get_discus_path
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
