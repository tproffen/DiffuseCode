MODULE parallel_mod
!-
!   Variables related to parallel processing using OMP
!+
LOGICAL   :: par_omp_use = .TRUE.      ! User wants to use OMP
INTEGER   :: par_omp_maxthreads = -1   ! Maximum number of threads to use, -1==use all available
INTEGER   :: par_omp_phys = 1          ! Number of physical cores
INTEGER   :: par_omp_logi = 1          ! Number of logical  cores in case of hyperthreading
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE get_cores()
!
USE envir_mod
USE errlist_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=1024), DIMENSION(2) :: string
CHARACTER(LEN=1024) :: ofile
INTEGER             :: ianz
!
IF(operating==OS_LINUX .OR. operating==OS_LINUX_WSL) THEN
   ofile = tmp_dir(1:LEN_TRIM(tmp_dir)) // '/discus.phys_cores'
   WRITE(string(1),'(a,a)') 'cat /proc/cpuinfo | grep ''cpu cores'' | uniq | awk ''{print $4}'' > ', &
            ofile(1:LEN_TRIM(ofile))
   WRITE(string(2),'(a,a)') 'cat /proc/cpuinfo | grep ''physical id'' | sort -u | wc -l >>', &
            ofile(1:LEN_TRIM(ofile))
   ianz = 2
   par_omp_phys =  read_phys(string, ofile, ianz)
!
   ofile = tmp_dir(1:LEN_TRIM(tmp_dir)) // '/discus.logi_cores'
   WRITE(string(1),'(a,a)') 'nproc > ', ofile(1:LEN_TRIM(ofile))
   ianz = 1
   par_omp_logi =  read_phys(string, ofile, ianz)
ELSEIF(operating==OS_MACOSX) THEN
   ofile = tmp_dir(1:LEN_TRIM(tmp_dir)) // '/discus.phys_cores'
   WRITE(string(1),'(a,a)') 'sysctl hw.physicalcpu | awk ''{print $2}'' > ', ofile(1:LEN_TRIM(ofile))
   ianz = 1
   par_omp_phys =  read_phys(string, ofile, ianz)
!
   ofile = tmp_dir(1:LEN_TRIM(tmp_dir)) // '/discus.logi_cores'
   WRITE(string(1),'(a,a)') 'sysctl hw.logicalcpu | awk ''{print $2}'' > ', ofile(1:LEN_TRIM(ofile))
   ianz = 1
   par_omp_logi =  read_phys(string, ofile, ianz)
ENDIF

END SUBROUTINE get_cores
!
!*******************************************************************************
!
INTEGER FUNCTION read_phys(string, ofile, ianz)
!
use precision_mod
USE support_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), DIMENSION(2), INTENT(INOUT) :: string
CHARACTER(LEN=*), INTENT(IN) :: ofile
INTEGER         , INTENT(IN) :: ianz
!
INTEGER, PARAMETER :: IRD = 88
!
character(len=prec_STRING) :: message
integer :: ier_cmd, exit_msg
INTEGER :: j
INTEGER :: k
INTEGER :: ios
!
read_phys = 1
DO k=1,ianz
   call execute_command_line(string(k), wait=.true.,                             &
        cmdstat=ier_cmd, cmdmsg=message, exitstat=exit_msg)
ENDDO
CALL oeffne(IRD, ofile, 'old')
!
DO k=1,ianz
   READ(IRD, *, IOSTAT=ios) j
   IF(ios==0) THEN
      read_phys = read_phys * j
   ENDIF
ENDDO
CLOSE(IRD)
WRITE(string(1),'(a,a)') 'rm -f ', ofile(1:LEN_TRIM(ofile))
call execute_command_line(string(1), wait=.false.,                               &
     cmdstat=ier_cmd, cmdmsg=message, exitstat=exit_msg)
!
END FUNCTION read_phys
!
!*******************************************************************************
!
END MODULE parallel_mod
