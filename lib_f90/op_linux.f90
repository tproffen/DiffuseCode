MODULE operating_mod
!
! Operating system dependent code
! LINUX Version
!
INTEGER, PARAMETER :: OP_LINUX   = 1
INTEGER, PARAMETER :: OP_MAC     = 2
INTEGER, PARAMETER :: OP_WINDOWS = 3
!
INTEGER, PARAMETER :: OP_SYSTEM = OP_LINUX
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
mpi_path = ' '  ! Turn MPI path empty
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
discus_path = '/usr/local/bin/'
discus_name = 'discus_suite'
!
END SUBROUTINE get_discus_path
!
SUBROUTINE operating_exit
!
USE prompt_mod
IMPLICIT NONE
!
! Currently no need for specifics
!
END SUBROUTINE operating_exit
END MODULE operating_mod
