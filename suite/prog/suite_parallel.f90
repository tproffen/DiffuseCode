MODULE suite_parallel_mod
!
IMPLICIT NONE
!
CONTAINS
!
SUBROUTINE suite_do_parallel ( zeile, length)
!
USE ber_params_mod
USE errlist_mod
USE get_params_mod
USE operating_mod
USE prompt_mod
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(IN) :: zeile
INTEGER         , INTENT(INOUT) :: length
!
INTEGER, PARAMETER :: MAXW = 30
CHARACTER(LEN=14)                    :: cdummy
CHARACTER(LEN= 6)                    :: cnumproc
CHARACTER(LEN=20)  , PARAMETER       :: env_num_proc='NUMBER_OF_PROCESSORS'
CHARACTER(LEN=1024)                  :: mfile
CHARACTER(LEN=1024)                  :: line
CHARACTER(LEN=1024), DIMENSION(MAXW) :: cpara
INTEGER            , DIMENSION(MAXW) :: lpara
REAL               , DIMENSION(MAXW) :: werte
INTEGER                              :: numproc
INTEGER                              :: ianz
INTEGER                              :: iianz
INTEGER                              :: i
INTEGER                              :: ios
LOGICAL                              :: lexist
CHARACTER(LEN=1024) :: mpi_path    = ' '
CHARACTER(LEN=1024) :: discus_path = ' '
CHARACTER(LEN=1024) :: discus_name = ' '
!
CALL get_mpi_path(mpi_path)
IF(mpi_path==' ') THEN
   ier_num = -8
   ier_typ = ER_COMM
   ier_msg(1) = 'The parallel command is not avaialable for '
   ier_msg(2) = 'your Operating system. Use an explicit command:'
   ier_msg(3) = 'mpiexec -n PROCESSES discus_suite -macro name.mac'
   RETURN
ENDIF
WRITE(output_io,*) 'MPI startup at WINDOWS might be very slow...'
CALL get_discus_path(discus_path, discus_name)
CALL holeenv(env_num_proc, cnumproc)
READ(cnumproc,*    ,IOSTAT=ios) numproc
IF(ios/=0) numproc = 3  ! Error, set to default
!
CALL get_params(zeile, ianz, cpara, lpara, MAXW, length)
IF(ier_num == 0) THEN
   IF(ianz == 0) THEN   ! Do auto-help
      cdummy = 'suite parallel'
      CALL do_hel(cdummy, LEN_TRIM(cdummy))
   ELSE
      IF(IANZ > 1) THEN  ! More than one parameter, check nprocessors
         iianz = 1
         CALL ber_params(iianz, cpara, lpara, werte, MAXW)
         IF(ier_num== 0) THEN  ! OK, 1st parameter is macro name
            numproc = NINT(werte(1))
            CALL del_params(1, ianz, cpara, lpara, MAXW)
         ENDIF
         CALL no_error   ! Switch error status off
         IF(numproc < 1) THEN
            ier_num = -6
            ier_typ = ER_COMM
            ier_msg(1) = ' Number of processors must be 1 or larger!'
            ier_msg(2) = ' Check value of first parameter'
         ENDIF
      ENDIF
      IF(ier_num == 0 .AND. IANZ >= 1) THEN  ! At least macro file, maybe params
         mfile = cpara(1)
         INQUIRE(FILE=mfile, EXIST=lexist)
         IF(lexist) THEN
            WRITE(line, '(3a,i6,4a)')  'mpiexec --prefix ',mpi_path(1:LEN_TRIM(mpi_path)), &
            ' -n ',numproc, ' ',discus_path(1:LEN_TRIM(discus_path)), &
            discus_name(1:LEN_TRIM(discus_name)),' -macro'
            DO I=1, ianz
               line = line(1:LEN_TRIM(line))//' '//cpara(i)(1:lpara(i))
            ENDDO
            CALL do_operating_comm(line)
         ELSE
            ier_num = -12
            ier_typ = ER_MAC
            ier_msg(1) ='Parallel command did not find the macro file:'
            ier_msg(2) = mfile(1:MIN(LEN(ier_msg),LEN_TRIM(mfile)))
            ier_msg(3) = 'Check the value/string, see help'
         ENDIF
      ENDIF
   ENDIF
ENDIF
!
END SUBROUTINE suite_do_parallel
!
END MODULE suite_parallel_mod
