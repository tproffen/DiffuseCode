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
USE envir_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=1024), INTENT(OUT) :: mpi_path
!
CHARACTER(LEN=1024)  :: mpi_file
LOGICAL              :: lda
!
IF(start_line(1:start_line_l) == 'discus_suite_noparallel') THEN
   mpi_file = '/bin/mpiexec'
   INQUIRE(FILE=mpi_file, EXIST=lda)
   IF(lda) THEN
      mpi_path = '/bin'
   ELSE
      mpi_file = '/usr/bin/mpiexec'
      INQUIRE(FILE=mpi_file, EXIST=lda)
      IF(lda) THEN
         mpi_path = '/usr/bin'
      ELSE
         mpi_path = ' '
      ENDIF
   ENDIF
ELSE
   mpi_path = ' '  ! Turn MPI path empty
ENDIF
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
USE envir_mod
USE prompt_mod
!
IMPLICIT NONE
!
! Currently no need for specifics
!
IF(operating==OS_LINUX_WSL) THEN
   WRITE(output_io,'(a)') 'DISCUS_SUITE is finished, please close this window'
ENDIF
!
END SUBROUTINE operating_exit
!
END MODULE operating_mod
