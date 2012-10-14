!*****7***************************************************************
SUBROUTINE RUN_MPI_INIT 
!
! MPI Version
!
USE mpi
USE run_mpi_mod
!
IMPLICIT none
INCLUDE 'errlist.inc'
!
CALL MPI_INIT (ier_num)
!
IF ( ier_num /= 0 ) THEN
   ier_msg(1) = 'MPI SYSTEM could not be initialized'
   WRITE(ier_msg(2),3000) ier_num
   ier_num = -22
   ier_typ = ER_APPL
   RETURN
ENDIF 
!
CALL MPI_COMM_RANK (MPI_COMM_WORLD, run_mpi_myid,     ier_num)
!
IF ( ier_num /= 0 ) THEN
   ier_msg(1) = 'MPI SYSTEM did not return RANK'
   WRITE(ier_msg(2),3000) ier_num
   ier_num = -22
   ier_typ = ER_APPL
   RETURN
ENDIF 
!
CALL MPI_COMM_SIZE (MPI_COMM_WORLD, run_mpi_numprocs, ier_num)
!
IF ( ier_num /= 0 ) THEN
   ier_msg(1) = 'MPI SYSTEM did not return SIZE'
   WRITE(ier_msg(2),3000) ier_num
   ier_num = -22
   ier_typ = ER_APPL
   RETURN
ENDIF 
!
! For future use with MPI_TYPE_...
!
!run_mpi_offsets(0)     = 0
!run_mpi_oldtypes(0)    = MPI_LOGICAL
!run_mpi_blockcounts(0) = 1
!
!CALL MPI_TYPE_EXTENT ( MPI_LOGICAL, run_mpi_extent, ier_num )
!
!run_mpi_offsets(1)     = run_mpi_extent
!run_mpi_oldtypes(1)    = MPI_INTEGER
!run_mpi_blockcounts(1) = RUN_MPI_COUNT_INTEGER
!
!CALL MPI_TYPE_EXTENT ( MPI_INTEGER, run_mpi_extent, ier_num )
!
!run_mpi_offsets(2)     = run_mpi_offsets(1) + run_mpi_blockcounts(0)*run_mpi_extent
!run_mpi_oldtypes(2)    = MPI_CHARACTER
!run_mpi_blockcounts(2) = 4*2048
!
!CALL MPI_TYPE_CREATE_STRUCT ( 3, run_mpi_blockcounts, run_mpi_offsets,  &
!     run_mpi_oldtypes, run_mpi_data_type, ier_num )
!CALL MPI_TYPE_COMMIT ( run_mpi_data_type, ier_num)
!!!!write(*,*) '############## data type , myid', run_mpi_data_type,run_mpi_myid, ier_num
!
WRITE(*,4000)
!
3000 FORMAT('MPI system returned error no. ',i8)
4000 FORMAT(1x,'MPI initilization successful ..')
!
END SUBROUTINE RUN_MPI_INIT
!
!*****7***************************************************************
SUBROUTINE RUN_MPI_MASTER 
!
USE mpi
USE population
USE run_mpi_mod
!
IMPLICIT none
INCLUDE 'errlist.inc'
!
INTEGER, DIMENSION(1:MPI_STATUS_SIZE) :: run_mpi_status
!
CHARACTER (LEN=2048)  :: send_direc
INTEGER               :: send_direc_l
INTEGER               :: sender
INTEGER               :: i
!DBG
INTEGER               :: ierr
!
INTEGER  :: len_str
INTEGER  :: system
!
CALL do_cwd ( send_direc, send_direc_l )    ! Get current working directory
run_mpi_senddata%direc_l = send_direc_l
run_mpi_senddata%direc   = send_direc
!
run_mpi_numsent = 0
run_mpi_numjobs = MIN ( run_mpi_numprocs - 1, pop_c * run_mpi_senddata%nindiv )
!
! SEND OUT INITIAL JOBS
!
!
!  COPY permanent part into INTEGER array send_data ! Just while structure does not work
!
run_mpi_send_data     = 0
IF (run_mpi_senddata%repeat) THEN
   run_mpi_send_data( 1) = 1
ENDIF
run_mpi_send_data( 2) = run_mpi_senddata%nindiv
run_mpi_send_data( 5) = 0
run_mpi_send_data( 6) = run_mpi_senddata%direc_l
run_mpi_send_data( 7) = run_mpi_senddata%prog_l
run_mpi_send_data( 8) = run_mpi_senddata%mac_l
run_mpi_send_data( 9) = run_mpi_senddata%out_l
DO i = 1,run_mpi_senddata%direc_l
   run_mpi_send_data(  9+i) = IACHAR(run_mpi_senddata%direc(i:i))
ENDDO
DO i = 1,run_mpi_senddata%prog_l
   run_mpi_send_data( 89+i) = IACHAR(run_mpi_senddata%prog (i:i))
ENDDO
DO i = 1,run_mpi_senddata%mac_l
   run_mpi_send_data(169+i) = IACHAR(run_mpi_senddata%mac  (i:i))
ENDDO
DO i = 1,run_mpi_senddata%out_l
   run_mpi_send_data(249+i) = IACHAR(run_mpi_senddata%out  (i:i))
ENDDO
!
!  Start initial jobs
!
DO i = 1, run_mpi_numjobs                   !  Start the intial jobs
   ier_num = 0
   run_mpi_senddata%member = mod( run_mpi_numsent,  pop_c) + 1
   run_mpi_senddata%kid    =      run_mpi_numsent / pop_c  + 1
!
!  Temporarily, while structure does not work
!
   run_mpi_send_data( 3) = run_mpi_senddata%member
   run_mpi_send_data( 4) = run_mpi_senddata%kid
!
!write(*,*) ' M sending Job ',i, run_mpi_senddata%member, run_mpi_senddata%kid
   CALL MPI_SEND ( run_mpi_send_data, 329, MPI_INTEGER, i,i, MPI_COMM_WORLD, ierr)
!====   CALL MPI_SEND ( run_mpi_senddata, 1, run_mpi_data_type, i, i, &
!====                   MPI_COMM_WORLD, ier_num )
!
   run_mpi_numsent = run_mpi_numsent + 1
ENDDO
!
!------       --Receive results and hand out new jobs
!
DO i = 1, pop_c * run_mpi_senddata%nindiv
   CALL MPI_RECV ( run_mpi_send_data, 329, MPI_INTEGER, MPI_ANY_SOURCE, MPI_ANY_TAG, &
               MPI_COMM_WORLD, run_mpi_status, ierr)
!====   CALL MPI_RECV ( run_mpi_senddata, 1, run_mpi_data_type, MPI_ANY_SOURCE, &
!====                   MPI_ANY_TAG, MPI_COMM_WORLD, run_mpi_status, ier_num )
!
   sender = run_mpi_status(MPI_SOURCE)
   IF ( run_mpi_numsent < pop_c*run_mpi_senddata%nindiv ) THEN
      run_mpi_senddata%member = mod( run_mpi_numsent,  pop_c) + 1
      run_mpi_senddata%kid    =      run_mpi_numsent / pop_c  + 1
!
!     Temporarily, while structure does not work
!
      run_mpi_send_data( 3) = run_mpi_senddata%member
      run_mpi_send_data( 4) = run_mpi_senddata%kid
!
      CALL MPI_SEND ( run_mpi_send_data, 329, MPI_INTEGER, sender,i, MPI_COMM_WORLD, ierr)
!====      CALL MPI_SEND ( run_mpi_senddata, 1, run_mpi_data_type, sender, i, &
!====                      MPI_COMM_WORLD, ier_num )
      run_mpi_numsent = run_mpi_numsent + 1
   ENDIF
ENDDO
!
!------       --End of loop over all members in the population
!
!
!
END SUBROUTINE RUN_MPI_MASTER
!
!*****7***************************************************************
!
SUBROUTINE RUN_MPI_SLAVE
!
USE mpi
USE run_mpi_mod
!
IMPLICIT none
INCLUDE 'errlist.inc'
!
CHARACTER (LEN=2048)   :: line
CHARACTER (LEN=2048)   :: output
!
INTEGER, DIMENSION(1:MPI_STATUS_SIZE) :: run_mpi_status
!
INTEGER  :: i
INTEGER  :: ierr
INTEGER  :: job_l
INTEGER  :: output_l
!DBG
!
INTEGER  :: len_str
INTEGER  :: system
!
ierr = 0
!
! Infinite loop, as long as new jobs come in, terminated my TAG=0
!
slave: DO
   CALL MPI_RECV ( run_mpi_send_data, 329, MPI_INTEGER, MPI_ANY_SOURCE, MPI_ANY_TAG, &
               MPI_COMM_WORLD, run_mpi_status, ierr)
!write(*,*) ' S received task ', run_mpi_send_data(1),run_mpi_send_data(2)
!====   CALL MPI_RECV ( run_mpi_senddata, 1, run_mpi_data_type, MPI_ANY_SOURCE, &
!====                   MPI_ANY_TAG, MPI_COMM_WORLD, run_mpi_status, ier_num )
!
!  temporarily, while TYPE does not work
!
run_mpi_senddata%repeat   = run_mpi_send_data ( 1)  == 1
run_mpi_senddata%nindiv   = run_mpi_send_data ( 2)
run_mpi_senddata%member   = run_mpi_send_data ( 3)
run_mpi_senddata%kid      = run_mpi_send_data ( 4)
run_mpi_senddata%direc_l  = run_mpi_send_data ( 6)
run_mpi_senddata%prog_l   = run_mpi_send_data ( 7)
run_mpi_senddata%mac_l    = run_mpi_send_data ( 8)
run_mpi_senddata%out_l    = run_mpi_send_data ( 9)
DO i = 1,run_mpi_senddata%direc_l
   run_mpi_senddata%direc(i:i) = ACHAR(run_mpi_send_data(  9+i))
ENDDO
DO i = 1,run_mpi_senddata%prog_l
   run_mpi_senddata%prog (i:i) = ACHAR(run_mpi_send_data( 89+i))
ENDDO
DO i = 1,run_mpi_senddata%mac_l
   run_mpi_senddata%mac  (i:i) = ACHAR(run_mpi_send_data(169+i))
ENDDO
DO i = 1,run_mpi_senddata%out_l
   run_mpi_senddata%out  (i:i) = ACHAR(run_mpi_send_data(249+i))
ENDDO
!write(*,*) ' S run_mpi_senddata%repeat',run_mpi_senddata%repeat, run_mpi_senddata%nindiv
!
!------       --- If TAG is zero exit, else compute
!
   IF (run_mpi_status (MPI_TAG) .eq.0) then
      EXIT slave
   ENDIF
!
!------       ----- Here we will call DISCUS/KUPLOT 
!
   repeat: IF ( run_mpi_senddata%repeat ) THEN
      IF ( run_mpi_senddata%out == '/dev/null' ) THEN
         output   = '/dev/null'
         output_l = 9
      ELSE
         WRITE(output,1000) run_mpi_senddata%out(1:run_mpi_senddata%out_l),&
                       run_mpi_senddata%member, run_mpi_senddata%kid
         output_l = run_mpi_senddata%out_l + 10
      ENDIF
      WRITE(line,2000) run_mpi_senddata%prog (1:run_mpi_senddata%prog_l ), &
                       run_mpi_senddata%mac  (1:run_mpi_senddata%mac_l  ), &
                       run_mpi_senddata%direc(1:run_mpi_senddata%direc_l), &
                       run_mpi_senddata%member, run_mpi_senddata%kid     , &
                       output(1:output_l)
   ELSE repeat                                      ! no repatition required
      IF ( run_mpi_senddata%out == '/dev/null' ) THEN
         output   = '/dev/null'
         output_l = 9
      ELSE
         WRITE(output,1100) run_mpi_senddata%out(1:run_mpi_senddata%out_l),&
                       run_mpi_senddata%member
         output_l = run_mpi_senddata%out_l + 5
      ENDIF
      WRITE(line,2100) run_mpi_senddata%prog (1:run_mpi_senddata%prog_l ), &
                       run_mpi_senddata%mac  (1:run_mpi_senddata%mac_l  ), &
                       run_mpi_senddata%direc(1:run_mpi_senddata%direc_l), &
                       run_mpi_senddata%member                           , &
                       output(1:output_l)
!write(*,*) 'line', line(1:len_str(line))
   ENDIF repeat

   job_l = len_str(line)
  ierr =  system ( line(1:job_l))
!
!  Answer back to master
!
   run_mpi_senddata%ierr = ierr
!
   CALL MPI_SEND ( run_mpi_send_data, 329, MPI_INTEGER, 0,0, MPI_COMM_WORLD, ierr)
!====   CALL MPI_SEND ( run_mpi_senddata, 1, run_mpi_data_type, 0, 0, &
!====                   MPI_COMM_WORLD, ier_num )
ENDDO slave
!
1000 FORMAT ( a,'.',i4.4,'.',i4.4)
1100 FORMAT ( a,'.',i4.4)
2000 FORMAT ( a,' -macro ',a,' ',a,i8,2x,i8,' > ',a) 
2100 FORMAT ( a,' -macro ',a,' ',a,i8,      ' > ',a) 
!
END SUBROUTINE RUN_MPI_SLAVE
!
!*****7***************************************************************
!
SUBROUTINE RUN_MPI_FINALIZE
!
USE mpi
USE run_mpi_mod
IMPLICIT none
INCLUDE 'errlist.inc'
!
INTEGER :: i
!
!------       -- Send termination signal to all slave processes
!
IF ( run_mpi_myid == 0 ) THEN
   run_mpi_send_data = 0
   DO i = 1, run_mpi_numprocs - 1
      ier_num = 0
      CALL MPI_SEND ( run_mpi_send_data, 329, MPI_INTEGER, i,0, MPI_COMM_WORLD, ier_num)
!====   CALL MPI_SEND ( run_mpi_senddata, 1, run_mpi_data_type, i, 0, &
!====                   MPI_COMM_WORLD, ier_num )
   ENDDO
ENDIF
!
CALL MPI_FINALIZE ( ier_num )
!
END SUBROUTINE RUN_MPI_FINALIZE
