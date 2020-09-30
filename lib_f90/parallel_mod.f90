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
CHARACTER(LEN=1024) :: string
CHARACTER(LEN=1024) :: ofile
!
IF(operating==OS_LINUX .OR. operating==OS_LINUX_WSL) THEN
   ofile = tmp_dir(1:LEN_TRIM(tmp_dir)) // '/discus.phys_cores'
   WRITE(string,'(a,a)') 'cat /proc/cpuinfo | grep ''cpu cores'' | uniq | awk ''{print $4}'' > ', &
            ofile(1:LEN_TRIM(ofile))
   par_omp_phys =  read_phys(string, ofile)
!
   WRITE(string,'(a,a)') 'nproc > ', ofile(1:LEN_TRIM(ofile))
   par_omp_logi =  read_phys(string, ofile)
ELSEIF(operating==OS_MACOSX) THEN
   ofile = tmp_dir(1:LEN_TRIM(tmp_dir)) // '/discus.phys_cores'
   WRITE(string,'(a,a)') 'sysctl hw.physicalcpu | awk ''{print $2}'' > ', ofile(1:LEN_TRIM(ofile))
   par_omp_phys =  read_phys(string, ofile)
!
   WRITE(string,'(a,a)') 'sysctl hw.logicalcpu | awk ''{print $2}'' > ', ofile(1:LEN_TRIM(ofile))
   par_omp_logi =  read_phys(string, ofile)
ENDIF

END SUBROUTINE get_cores
!
!*******************************************************************************
!
INTEGER FUNCTION read_phys(string, ofile)
!
USE support_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: string
CHARACTER(LEN=*), INTENT(IN) :: ofile
!
INTEGER, PARAMETER :: IRD = 88
!
INTEGER :: j
INTEGER :: ios
!
read_phys = 1
CALL EXECUTE_COMMAND_LINE(string)
CALL oeffne(IRD, ofile, 'old')
READ(IRD, *, IOSTAT=ios) j
CLOSE(IRD)
IF(ios==0) THEN
   read_phys = j
ENDIF
WRITE(string,'(a,a)') 'rm -f ', ofile(1:LEN_TRIM(ofile))
CALL EXECUTE_COMMAND_LINE(string)
!
END FUNCTION read_phys
!
!*******************************************************************************
!
END MODULE parallel_mod
