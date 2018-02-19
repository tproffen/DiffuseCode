MODULE DIFFEV_MPI_MOD
!
PUBLIC :: run_mpi_init, run_mpi_master, run_mpi_slave, run_mpi_finalize
CONTAINS
!*****7***************************************************************
SUBROUTINE run_mpi_init 
!
! MPI Version
!
! Runs a refinement in parallel using MPI
!
! two options exist, one to start the slave program via a system call, 
! the other starting the slave via a socket.
!
! In the system mode, the program is executed and closed upon each
! individual calculation. The slave program must read GENERATION and 
! the trial files. The program is started with the command:
! discus -macro <pwd>, <kid>, <indiv> > <logfile>
!
! In the socket mode, the program is started at the start of the distributed 
! calculation and remains active until the end, i.e. throughout one entire
! generation. The variables in GENERATION, the variable NINDIV, and the trial
! parameters are passed to the slave program. For each kid/indiv combination
! the slave macro is started and should not expect any macro parameters.
!
USE mpi
USE run_mpi_mod
USE population
USE times_mod
!
USE errlist_mod
USE mpi_slave_mod
USE prompt_mod
!
IMPLICIT none
!
INTEGER                        :: run_mpi_integer_extent
INTEGER                        :: run_mpi_logical_extent
INTEGER                        :: run_mpi_charact_extent
INTEGER                        :: run_mpi_real_extent
!
CHARACTER (LEN=MPI_MAX_PROCESSOR_NAME) :: node_name   = ' '
INTEGER, DIMENSION(1:MPI_STATUS_SIZE) :: run_mpi_status
INTEGER               :: sender        ! Id of slave that answered
INTEGER               :: local_id      ! BUG Patch MPI_ID messy with structure
INTEGER               :: i, j
INTEGER               :: length
INTEGER               :: ierr
LOGICAL               :: success
!
CALL MPI_INIT (ier_num)       ! initialize the MPI system
!
IF ( ier_num /= 0 ) THEN
   ier_msg(1) = 'MPI SYSTEM could not be initialized'
   WRITE(ier_msg(2),3000) ier_num
   ier_num = -22
   ier_typ = ER_APPL
   RETURN
ENDIF 
!
!  Get the rank of this process store in run_mpi_myid
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
!  Get the size of the distribution, store in run_mpi_numprocs
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
IF ( run_mpi_numprocs < 2 ) THEN
!  ier_msg(1) = 'MPI SYSTEM returned one CPU   '
!  ier_msg(2) = 'DIFFEV must be started with   '
!  ier_msg(3) = 'mpiexec -n X diffev; x >= 2   '
!  WRITE(ier_msg(2),3000) ier_num
!
!  MPI Does not seem to be active, quietly turn off
   ier_num = 0
   ier_typ = ER_NONE
   RETURN
ENDIF
!
! For future use with MPI_TYPE_...
!
CALL MPI_TYPE_EXTENT ( MPI_INTEGER,   run_mpi_integer_extent, ier_num )
CALL MPI_TYPE_EXTENT ( MPI_LOGICAL,   run_mpi_logical_extent, ier_num )
CALL MPI_TYPE_EXTENT ( MPI_CHARACTER, run_mpi_charact_extent, ier_num )
CALL MPI_TYPE_EXTENT ( MPI_REAL     , run_mpi_real_extent   , ier_num )
!
run_mpi_offsets(0)     = 0
run_mpi_oldtypes(0)    = MPI_INTEGER
run_mpi_blockcounts(0) = RUN_MPI_COUNT_INTEGER*run_mpi_integer_extent
!
run_mpi_offsets(1)     = run_mpi_offsets(0) +  run_mpi_integer_extent*RUN_MPI_COUNT_INTEGER
run_mpi_oldtypes(1)    = MPI_LOGICAL
run_mpi_blockcounts(1) = RUN_MPI_COUNT_LOGICAL*run_mpi_logical_extent
!
run_mpi_offsets(2)     = run_mpi_offsets(1) +  run_mpi_logical_extent*RUN_MPI_COUNT_LOGICAL
run_mpi_oldtypes(2)    = MPI_CHARACTER
run_mpi_blockcounts(2) = RUN_MPI_COUNT_CHARACTER*run_mpi_charact_extent
!
run_mpi_offsets(3)     = run_mpi_offsets(2) +  run_mpi_charact_extent*RUN_MPI_COUNT_CHARACTER
run_mpi_oldtypes(3)    = MPI_REAL
run_mpi_blockcounts(3) = RUN_MPI_COUNT_REAL   *run_mpi_real_extent
!
run_mpi_offsets(4)     = run_mpi_offsets(3) +  run_mpi_real_extent*RUN_MPI_COUNT_REAL
run_mpi_oldtypes(4)    = MPI_REAL
run_mpi_blockcounts(4) = RUN_MPI_COUNT_TRIAL  *run_mpi_real_extent
!
CALL MPI_TYPE_STRUCT ( 5, run_mpi_blockcounts, run_mpi_offsets,  &
     run_mpi_oldtypes, run_mpi_data_type, ier_num )
CALL MPI_TYPE_COMMIT ( run_mpi_data_type, ier_num)
!write(*,*) '############## data type , myid', run_mpi_data_type,run_mpi_myid, ier_num
!
!WRITE(*,4000)
!
run_mpi_active = .true.
!
socket_status = PROMPT_OFF  ! Turn off socket responses
!
mpi_active = .true.
!
!  Build a list of nodes, node names, and assignment of the processes onto these nodes
!
local_id = run_mpi_myid
IF(run_mpi_myid==0)  THEN  !   MASTER
   IF(ALLOCATED(node_names)) DEALLOCATE(node_names)
   IF(ALLOCATED(node_names)) DEALLOCATE(slave_on_node)
   IF(ALLOCATED(node_names)) DEALLOCATE(slave_is_node)
   IF(ALLOCATED(node_names)) DEALLOCATE(slave_per_node)
   ALLOCATE(node_names    (1:run_mpi_numprocs))
   ALLOCATE(slave_on_node (1:run_mpi_numprocs))
   ALLOCATE(slave_is_node (1:run_mpi_numprocs))
   ALLOCATE(slave_per_node(1:run_mpi_numprocs))
   node_names    (:) = ' '
   slave_on_node (:) = ' '
   slave_is_node (:) = -1
   slave_per_node(:) = 0
   NUM_NODE = 0
   sdl_length = 1! 580! + 200
   DO i=1, run_mpi_numprocs-1
      run_mpi_senddata%kid = -i
      CALL MPI_SEND ( run_mpi_senddata, 1, run_mpi_data_type, i, i, &
                      MPI_COMM_WORLD, ier_num )
   ENDDO
   DO i=1, run_mpi_numprocs-1
      CALL MPI_RECV ( run_mpi_senddata, 1, run_mpi_data_type, MPI_ANY_SOURCE, &
                      MPI_ANY_TAG, MPI_COMM_WORLD, run_mpi_status, ier_num )
      run_mpi_myid = local_id ! BUG PATCH MPI_ID gets messed up by receive with structure?????
      sender = run_mpi_status(MPI_SOURCE)     ! Identify the slave
      success = .FALSE.
      search: DO j=1,NUM_NODE
         IF(run_mpi_senddata%direc == node_names(j)) THEN
            success = .TRUE.
            slave_is_node(sender) = j
            slave_per_node(j)     = slave_per_node(j) + 1
            EXIT search
         ENDIF
      ENDDO search
      IF(.NOT.success) THEN
         NUM_NODE = NUM_NODE + 1
         node_names(NUM_NODE) = run_mpi_senddata%direc(1:run_mpi_senddata%direc_l)
         slave_is_node(sender) = NUM_NODE
         slave_per_node(NUM_NODE) = 1
      ENDIF
      slave_on_node(sender) = run_mpi_senddata%direc
   ENDDO
   run_mpi_max_slaves = MAXVAL(slave_per_node,1)
   IF(ALLOCATED(node_has_slaves)) DEALLOCATE(node_has_slaves)
   ALLOCATE(node_has_slaves(1:NUM_NODE,0:run_mpi_max_slaves))
   node_has_slaves(:,:) = 0
   DO j=1, run_mpi_numprocs-1
      i = slave_is_node(j)
      node_has_slaves(i,0) = node_has_slaves(i,0) + 1
      node_has_slaves(i,node_has_slaves(i,0)) = j
   ENDDO
!
ELSE
   CALL MPI_PROBE( MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, run_mpi_status, ierr) ! Querry incomming size
   CALL MPI_GET_COUNT(run_mpi_status, run_mpi_data_type, sdl_length, ierr)                  ! Determine size
!
!  Now receive incomming message
!
   CALL MPI_RECV ( run_mpi_senddata, sdl_length, run_mpi_data_type, MPI_ANY_SOURCE, &
                   MPI_ANY_TAG, MPI_COMM_WORLD, run_mpi_status, ier_num )
   run_mpi_myid = local_id ! BUG PATCH MPI_ID gets messed up by receive with structure?????
   CALL MPI_GET_PROCESSOR_NAME ( node_name, length, ier_num)
   run_mpi_senddata%direc   = node_name
   run_mpi_senddata%direc_l = length

   CALL MPI_SEND ( run_mpi_senddata, 1, run_mpi_data_type, 0, 0, &
                   MPI_COMM_WORLD, ier_num )
ENDIF
!
3000 FORMAT('MPI system returned error no. ',i8)
!4000 FORMAT(1x,'MPI initilization successful ..')
!
END SUBROUTINE run_mpi_init
!
!*****7***************************************************************
SUBROUTINE run_mpi_master 
!
!  Main routine of the distribution. packs the relevant parameters and
!  sends them to run_mpi_numprocs -1 slaves
!
USE diffev_allocate_appl
USE diffev_random
USE mpi
USE population
USE run_mpi_mod
USE errlist_mod
USE prompt_mod
!
IMPLICIT none
!
INTEGER, DIMENSION(1:MPI_STATUS_SIZE) :: run_mpi_status
!
CHARACTER (LEN=2048)  :: send_direc    ! working directory
INTEGER               :: send_direc_l  ! working directory length
INTEGER               :: sender        ! Id of slave that answered
INTEGER               :: i,j
INTEGER               :: nseeds        ! Number of seeds for randum numbers
INTEGER               :: run_mpi_numsent  ! Number of jobs sent out
INTEGER               :: run_mpi_numjobs  ! Number of initial jobs
INTEGER               :: nprog         ! number of different program/mac combinations
LOGICAL               :: prog_start    ! External program needs to be started
!
LOGICAL               :: prog_exist    ! program/macro combination exists in data base
INTEGER               :: local_id      ! BUG Patch MPI_ID messy with structure
INTEGER               :: inode         ! Current slave process
INTEGER               :: min_kid       ! Current kid with least completed INDIVs
INTEGER               :: lastgen = -1  ! At the previous call to Rrun_mpi_master we were inthis GENERATION
!
local_id = run_mpi_myid  ! BUG PATCH MPI_ID gets messed up be receive with structure
!
sdl_length = 1 !580! + 200
!
! Test if program / macro combination exist /not exists ; sockets only
!
IF (run_mpi_senddata%use_socket) THEN
   prog_exist = .false.
   find_prog_entry: DO i=1,run_mpi_nprog        ! Search existing entries
      IF(prog_entry(i) == run_mpi_senddata%prog(1:run_mpi_senddata%prog_l) &
                          // ' '//                                         &
                          run_mpi_senddata%mac(1:run_mpi_senddata%mac_l)) THEN
         prog_exist = .true.
         prog_start = .false.                   ! Program should already be running
         run_mpi_senddata%prog_num = i          ! This is program/macro identifier
         run_mpi_senddata%s_remote = socket_id(i) ! This is program/macro identifier
         EXIT find_prog_entry
      ENDIF
   ENDDO find_prog_entry
!
   IF(.NOT. prog_exist) THEN                    ! New program to be started by slaves
      IF(run_mpi_nprog==RUN_MPI_MAXPROG) THEN   ! Need more space
         nprog = RUN_MPI_MAXPROG + 2            ! increment by two programs
         CALL alloc_socket_nprogs ( nprog, run_mpi_numprocs)
      ELSE                                      ! Sufficient space
         run_mpi_nprog = run_mpi_nprog + 1      ! Increment no of known prog/mac entries
         prog_entry(run_mpi_nprog) = run_mpi_senddata%prog(1:run_mpi_senddata%prog_l) &
                                     // ' '//                                         &
                                     run_mpi_senddata%mac(1:run_mpi_senddata%mac_l)
      ENDIF
      run_mpi_senddata%prog_num = run_mpi_nprog ! This is program/macro identifier
      prog_start = .true.                       ! Program needs to be started
      port_id  (run_mpi_nprog,:) = 0            ! Initial port number
      socket_id(run_mpi_nprog)   = 0            ! Initial socket ID
   ENDIF
ELSE
   prog_start = .true.                          ! Always start non socket program
ENDIF
!
run_mpi_senddata%prog_start = prog_start        ! Copy start flag into send structure
!
CALL do_cwd ( send_direc, send_direc_l )        ! Get current working directory
run_mpi_senddata%direc_l = send_direc_l         ! Copy directory into send structure
run_mpi_senddata%direc   = send_direc(1:MIN(send_direc_l,200))
!
run_mpi_numsent = 0                             ! No jobs sent yet
run_mpi_numjobs = MIN ( run_mpi_numprocs - 1, pop_c * run_mpi_senddata%nindiv )
!
run_mpi_senddata%l_get_state = l_get_random_state  ! Inquire random number status
!
node_has_kids(:,:,2) =  0                   ! No individuals are done
IF(pop_gen > lastgen) THEN                 ! New GENERATION , new job distributio
   node_has_kids(:,:,:) = 0
   kid_on_node(:)       = 0
ENDIF
!
!  Start initial jobs
!
initial:DO i = 1, run_mpi_numjobs                   !  Start the intial jobs
   ier_num = 0
   run_mpi_senddata%kid    = mod( run_mpi_numsent,  pop_c) + 1
   run_mpi_senddata%indiv  =      run_mpi_numsent / pop_c  + 1
   IF(pop_gen > lastgen) THEN                 ! New GENERATION , new job distribution
      kid_on_node(run_mpi_senddata%kid) = slave_is_node(i)
      kid_on_node(0                   ) = kid_on_node  (0) + 1
      node_has_kids(slave_is_node(i),0,1) = node_has_kids(slave_is_node(i),0,1) + 1 !Increment number of kids on this node
      node_has_kids(slave_is_node(i),node_has_kids(slave_is_node(i),0,1),1) = run_mpi_senddata%kid
      node_has_kids(slave_is_node(i),node_has_kids(slave_is_node(i),0,1),2) = run_mpi_senddata%indiv
   ELSE                                     ! Same GENERATION use old Job DISTRIBUTION
      inode = slave_is_node(i)             ! We are on this node
      min_kid = MINLOC(node_has_kids(inode,1:node_has_kids(inode,0,1),2),1)
      IF( node_has_kids(inode,min_kid,2) < run_mpi_senddata%nindiv) THEN   ! More indivs are needed
         node_has_kids(inode,min_kid,2) = node_has_kids(inode,min_kid,2) + 1
         run_mpi_senddata%kid    = node_has_kids(inode, min_kid,1)
         run_mpi_senddata%indiv  = node_has_kids(inode, min_kid,2)
      ELSE
         CYCLE initial
      ENDIF
   ENDIF
   IF(run_mpi_senddata%prog_start) THEN     ! Program needs to be started, increment port no
      run_mpi_senddata%port     = 2000 + MOD(run_mpi_numprocs*run_mpi_senddata%prog_num + i,3600)
   ELSE                                     ! Program is running use old port no
      run_mpi_senddata%port     = port_id  (run_mpi_senddata%prog_num,i)
   ENDIF
   DO j=1,pop_dimx                          ! Encode current trial values
      run_mpi_senddata%trial_names (j) = pop_name(j                  ) ! Takes value for kid that answered
      run_mpi_senddata%trial_values(j) = pop_t(j,run_mpi_senddata%kid) ! Takes value for kid that answered
   ENDDO
!
   sdl_length = 1! 580! + 200
   CALL MPI_SEND ( run_mpi_senddata, 1, run_mpi_data_type, i, i, &
                   MPI_COMM_WORLD, ier_num )
!
   run_mpi_numsent = run_mpi_numsent + 1
ENDDO initial
!
!------       --Receive results and hand out new jobs
!
rec_hand: DO i = 1, pop_c * run_mpi_senddata%nindiv
   CALL MPI_RECV ( run_mpi_senddata, 1, run_mpi_data_type, MPI_ANY_SOURCE, &
                   MPI_ANY_TAG, MPI_COMM_WORLD, run_mpi_status, ier_num )
   run_mpi_myid = local_id ! BUG PATCH MPI_ID gets messed up by receive with structure?????
!
   sender = run_mpi_status(MPI_SOURCE)     ! Identify the slave
   IF(run_mpi_senddata%ierr /=0 ) THEN
      ier_msg(1) = 'A slave program exited with error message'
      WRITE(ier_msg(2), 2000)  sender, ier_num
      ier_num = -26
      ier_typ = ER_APPL
2000 FORMAT('Error message from slave ',I4,' is ',i8)
      EXIT rec_hand
   ENDIF
   IF(run_mpi_senddata%use_socket ) THEN   ! SOCKET option is active
       socket_id(run_mpi_senddata%prog_num) = run_mpi_senddata%s_remote  ! Store socket info in data base
       port_id  (run_mpi_senddata%prog_num,sender) = run_mpi_senddata%port  ! Store port info in data base
   ENDIF
!
   IF(l_get_random_state) THEN             ! Update random number status
      pop_random(:,run_mpi_senddata%kid) = 0
      nseeds = run_mpi_senddata%nseeds
      j = MIN(nseeds, RUN_MPI_NSEEDS)
      pop_random(1:j,run_mpi_senddata%kid) = run_mpi_senddata%seeds(1:j)
      pop_random(0  ,run_mpi_senddata%kid) = j
   ENDIF
   IF(run_mpi_senddata%l_rvalue) THEN      ! R-value is returned
      j = run_mpi_senddata%n_rvalue_o
      n_rvalue_o = j
      trial_val(run_mpi_senddata%kid,0:j) = run_mpi_senddata%rvalue(0:j)
   ENDIF
!
   IF ( run_mpi_numsent < pop_c*run_mpi_senddata%nindiv ) THEN  ! There are more jobs to do
      run_mpi_senddata%kid    = mod( run_mpi_numsent,  pop_c) + 1
      run_mpi_senddata%indiv  =      run_mpi_numsent / pop_c  + 1
      inode = slave_is_node(sender)         ! We are on this node
      IF(pop_gen > lastgen .AND. kid_on_node(0)<pop_c .AND.  &
         node_has_kids(inode,0,1) < slave_per_node(inode)*run_mpi_kid_per_core ) THEN    ! Kids need to be distributed
         kid_on_node(0) = kid_on_node(0) + 1
         kid_on_node(kid_on_node(0)) = inode  ! Place next kid onto node of this slave
         node_has_kids(inode,0,1) = node_has_kids(inode,0,1) + 1 !Increment number of kids on this node
         node_has_kids(inode,node_has_kids(inode,0,1),1) = kid_on_node(0)
         node_has_kids(inode,node_has_kids(inode,0,1),2) = 1
         run_mpi_senddata%kid    = kid_on_node(0)
         run_mpi_senddata%indiv  = 1
         min_kid = node_has_kids(inode,0,1)
      ELSE                             ! All kids are distributed, increment INDIV
         min_kid = MINLOC(node_has_kids(inode,1:node_has_kids(inode,0,1),2),1)
         IF( node_has_kids(inode,min_kid,2) < run_mpi_senddata%nindiv) THEN   ! More indivs are needed
            node_has_kids(inode,min_kid,2) = node_has_kids(inode,min_kid,2) + 1
            run_mpi_senddata%kid    = node_has_kids(inode, min_kid,1)
            run_mpi_senddata%indiv  = node_has_kids(inode, min_kid,2)
         ELSE
            CYCLE rec_hand
         ENDIF
      ENDIF
      DO j=1,pop_dimx                          ! Encode current trial values
         run_mpi_senddata%trial_names (j) = pop_name(j                  ) ! Takes value for kid that answered
         run_mpi_senddata%trial_values(j) = pop_t(j,run_mpi_senddata%kid) ! Takes value for kid that answered
      ENDDO
!
      CALL MPI_SEND ( run_mpi_senddata, 1, run_mpi_data_type, sender, i, &
                      MPI_COMM_WORLD, ier_num )
!
      run_mpi_numsent = run_mpi_numsent + 1
   ENDIF
ENDDO rec_hand
IF(ier_num == -26) THEN   ! Fatal error occured, wait for remaining jobs
   DO i=1, run_mpi_numjobs-1
         CALL MPI_RECV ( run_mpi_senddata, 1, run_mpi_data_type, MPI_ANY_SOURCE, &
                   MPI_ANY_TAG, MPI_COMM_WORLD, run_mpi_status, ier_num )
   ENDDO
   ier_num = -26 ! Reinstate the error message
ENDIF
run_mpi_myid = local_id  ! BUG PATCH MPI_ID gets messed up by receive with structure?????
!
!------       --End of loop over all kids in the population
!
lastgen = pop_gen
!
END SUBROUTINE run_mpi_master
!
!*****7***************************************************************
!
SUBROUTINE run_mpi_slave
!
!  The slaves run in an indefinite loop until they are closed
!
!  The application program is started as a "system program" call
!  or via a socket
!
USE mpi
USE run_mpi_mod
USE diffev_setup_mod
!
USE set_sub_generic_mod
USE errlist_mod
USE mpi_slave_mod
USE prompt_mod
USE variable_mod
!
IMPLICIT none
!
CHARACTER (LEN=2048)   :: line
CHARACTER (LEN=2048)   :: output
!
INTEGER, DIMENSION(1:MPI_STATUS_SIZE) :: run_mpi_status
!
INTEGER  :: j
INTEGER  :: il
INTEGER  :: ierr
INTEGER  :: job_l
INTEGER  :: output_l
INTEGER  :: port  = 2000
!
INTEGER, PARAMETER :: seconds = 1       ! Wait time for application to start up
!
INTEGER  :: len_str
INTEGER  :: system
INTEGER  :: socket_close
INTEGER  :: socket_connect
INTEGER  :: socket_get
INTEGER  :: socket_send
!
INTEGER  :: local_id  ! BUG PATCH MPI_ID gets messed up by receive with structure?????
!
ierr = 0
!
local_id = run_mpi_myid ! BUG PATCH MPI_ID gets messed up by receive with structure?????
!
! Infinite loop, as long as new jobs come in, terminated my TAG=0
!
slave: DO
   CALL MPI_PROBE( MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, run_mpi_status, ierr) ! Querry incomming size
   CALL MPI_GET_COUNT(run_mpi_status, run_mpi_data_type, sdl_length, ierr)                  ! Determine size
!
!  Now receive incomming message
!
   CALL MPI_RECV ( run_mpi_senddata, sdl_length, run_mpi_data_type, MPI_ANY_SOURCE, &
                   MPI_ANY_TAG, MPI_COMM_WORLD, run_mpi_status, ier_num )
   run_mpi_myid = local_id ! BUG PATCH MPI_ID gets messed up by receive with structure?????
!
   s_remote = run_mpi_senddata%s_remote
   port     = run_mpi_senddata%port
!
!------       --- If TAG is zero exit, else compute
!
   tag_exit: IF (run_mpi_status (MPI_TAG) == 0) then   ! End of diffev program
      IF ( run_mpi_senddata%use_socket ) THEN          ! SOCKET option is active
         line  = 'exit'
         job_l = 4
         ierr = socket_send    (s_remote, line, job_l) ! End application program
         CALL sleep(seconds)                           ! wait for application to close
         ierr = socket_close   (s_remote)              ! Close socket
         lremote = .false.
         IF(run_mpi_senddata%prog_num == -1) THEN         ! This is the last program to close
            lremote = .false.
            EXIT slave                                    ! now we can end slave
         ENDIF
      ELSE
         lremote = .false.
         EXIT slave                                    ! now we can end slave
      ENDIF
   ELSE  tag_exit                                      ! Normal operation
!
!------       Define output redirection
!
      out_null: IF ( run_mpi_senddata%out(1:run_mpi_senddata%out_l) == '/dev/null' ) THEN   ! no output desired
         output   = '/dev/null'
         output_l = 9
      ELSE out_null                                          ! redirect to a log file
         out_socket: IF ( run_mpi_senddata%use_socket ) THEN ! SOCKET option is active
               WRITE(output,1100) run_mpi_senddata%out(1:run_mpi_senddata%out_l),&
                             run_mpi_myid
               output_l = run_mpi_senddata%out_l + 5
         ELSE out_socket                                     ! No socket active
            out_repeat:IF ( run_mpi_senddata%repeat ) THEN   ! call with kid AND indiv
               WRITE(output,1000) run_mpi_senddata%out(1:run_mpi_senddata%out_l),&
                             run_mpi_senddata%kid   , run_mpi_senddata%indiv
               output_l = run_mpi_senddata%out_l + 10
            ELSE out_repeat                                  ! call with kid only
               WRITE(output,1100) run_mpi_senddata%out(1:run_mpi_senddata%out_l),&
                             run_mpi_senddata%kid   
               output_l = run_mpi_senddata%out_l + 5
            ENDIF out_repeat
         ENDIF out_socket
      ENDIF out_null
!
!------       ----- Here we will call DISCUS/KUPLOT 
!
      use_socket: IF ( run_mpi_senddata%use_socket ) THEN   ! SOCKET option is active
         IF ( run_mpi_senddata%prog_start ) THEN            ! server needs to be started
            WRITE(line,4000) run_mpi_senddata%prog (1:run_mpi_senddata%prog_l ), &
                 port, output(1:output_l)
            job_l = run_mpi_senddata%prog_l + 42 + output_l
            ierr =  system ( line(1:job_l))                 ! Start with -remote ... optios
            CALL sleep(seconds)                             ! wait for application to start
            ierr = socket_connect (s_remote, '127.0.0.1',  9, port) ! Connect to socket
            ierr = socket_get     (s_remote, line, il)      ! Get socket answer
!
!           Send standard variable names and values
!
            ierr = socket_send    (s_remote, 'variable integer, port      ', 28)
            CALL   socket_wait
            ierr = socket_send    (s_remote, 'variable integer, myid      ', 28)
            CALL   socket_wait
            ierr = socket_send    (s_remote, 'variable integer, generation', 28)
            CALL   socket_wait
            ierr = socket_send    (s_remote, 'variable integer, member    ', 28)
            CALL   socket_wait
            ierr = socket_send    (s_remote, 'variable integer, children  ', 28)
            CALL   socket_wait
            ierr = socket_send    (s_remote, 'variable integer, parameters', 28)
            CALL   socket_wait
            ierr = socket_send    (s_remote, 'variable integer, nindiv    ', 28)
            CALL   socket_wait
            ierr = socket_send    (s_remote, 'variable integer, indiv     ', 28)
            CALL   socket_wait
            ierr = socket_send    (s_remote, 'variable integer, kid       ', 28)
            CALL   socket_wait
            lremote = .true.         ! remote Socket program is running
         ENDIF
!
!        Socket should  be active at this point
         WRITE(line, 4010) 'port      ',port                        ! Current member size
         ierr = socket_send    (s_remote, line, 21)
         CALL   socket_wait
         WRITE(line, 4010) 'myid      ',run_mpi_myid                ! Current member size
         ierr = socket_send    (s_remote, line, 21)
         CALL   socket_wait
         WRITE(line, 4010) 'member    ',run_mpi_senddata%member     ! Current member size
         ierr = socket_send    (s_remote, line, 21)
         CALL   socket_wait
         WRITE(line, 4010) 'children  ',run_mpi_senddata%children   ! Current children size
         ierr = socket_send    (s_remote, line, 21)
         CALL   socket_wait
         WRITE(line, 4010) 'parameters',run_mpi_senddata%parameters ! Current parameter size
         ierr = socket_send    (s_remote, line, 21)
         CALL   socket_wait
         WRITE(line, 4010) 'nindiv    ',run_mpi_senddata%nindiv     ! Current member size
         ierr = socket_send    (s_remote, line, 21)
         CALL   socket_wait
         WRITE(line, 4020) run_mpi_senddata%direc(1:run_mpi_senddata%direc_l) ! Current directory
         job_l = len_str(line)
         ierr = socket_send    (s_remote, line, job_l)
         CALL   socket_wait
         WRITE(line, 4010) 'generation',run_mpi_senddata%generation !Current Generation
         ierr = socket_send    (s_remote, line, 21)
         CALL   socket_wait
         WRITE(line, 4010) 'indiv     ',run_mpi_senddata%indiv      ! Current indiv numner
         ierr = socket_send    (s_remote, line, 21)
         CALL   socket_wait
         WRITE(line, 4010) 'kid       ',run_mpi_senddata%kid        ! Current kid numner
         ierr = socket_send    (s_remote, line, 21)
         CALL   socket_wait
         DO j=1,run_mpi_senddata%parameters                         ! Current trial values
            WRITE(line, 4040) 200+j,run_mpi_senddata%trial_values(j)
!            line = ' ' !WORK WORK WORK 
            ierr = socket_send    (s_remote, line, 29)
            CALL   socket_wait
         ENDDO
!        Now send the macro
         WRITE(line, 4030) run_mpi_senddata%mac  (1:run_mpi_senddata%mac_l  )
         job_l = len_str(line)
         ierr = socket_send    (s_remote, line, job_l)     ! Send macro
         CALL   socket_wait
      ELSE use_socket                                      ! explicitely start DISCUS/KUPLOT
         repeat: IF ( run_mpi_senddata%repeat ) THEN       ! NINDIV calculations needed
            WRITE(line,2000) run_mpi_senddata%prog (1:run_mpi_senddata%prog_l ), &
                             run_mpi_senddata%mac  (1:run_mpi_senddata%mac_l  ), &
                             run_mpi_senddata%direc(1:run_mpi_senddata%direc_l), &
                             run_mpi_senddata%kid   , run_mpi_senddata%indiv   , &
                             output(1:output_l)
      ELSE repeat                                          ! no repetition required
         WRITE(line,2100) run_mpi_senddata%prog (1:run_mpi_senddata%prog_l ), &
                          run_mpi_senddata%mac  (1:run_mpi_senddata%mac_l  ), &
                          run_mpi_senddata%direc(1:run_mpi_senddata%direc_l), &
                          run_mpi_senddata%kid                              , &
                          output(1:output_l)
      ENDIF repeat
      job_l = len_str(line)
      mpi_is_slave = .true.
      mpi_slave_error = 0
!             Execute the "generic" cost function calculation
      CALL p_execute_cost( run_mpi_senddata%repeat,                          &
                           LEN(run_mpi_senddata%prog),                       &
                           run_mpi_senddata%prog , run_mpi_senddata%prog_l , &
                           LEN(run_mpi_senddata%mac),                        &
                           run_mpi_senddata%mac  , run_mpi_senddata%mac_l  , &
                           LEN(run_mpi_senddata%direc),                      &
                           run_mpi_senddata%direc, run_mpi_senddata%direc_l, &
                           run_mpi_senddata%kid  , run_mpi_senddata%indiv  , &
                           run_mpi_senddata%n_rvalue_i, run_mpi_senddata%n_rvalue_o, &
                           RUN_MPI_MAXRVALUE,     &
                           run_mpi_senddata%rvalue, run_mpi_senddata%l_rvalue,     &
                           LEN(output),                                      &
                           output                , output_l ,                      &
                           run_mpi_senddata%generation, run_mpi_senddata%member,   &
                           run_mpi_senddata%children, run_mpi_senddata%parameters, &
                                                   run_mpi_senddata%nindiv  , &
                           run_mpi_senddata%trial_names ,                          &
                           run_mpi_senddata%trial_values, RUN_MPI_COUNT_TRIAL,     &
                           run_mpi_senddata%l_get_state,                           &
                           run_mpi_senddata%nseeds, run_mpi_senddata%seeds,        &
                           ierr )
   ENDIF use_socket
!
!  Answer back to master)
!
   run_mpi_senddata%ierr     = ierr
   run_mpi_senddata%s_remote = s_remote
   run_mpi_senddata%port     = port
!
   CALL MPI_SEND ( run_mpi_senddata, 1, run_mpi_data_type, 0, 0, &
                   MPI_COMM_WORLD, ier_num )
   ENDIF tag_exit
ENDDO slave
run_mpi_myid = local_id ! BUG PATCH MPI_ID gest messed up by receive with structure?????
!
1000 FORMAT ( a,'.',i4.4,'.',i4.4)
1100 FORMAT ( a,'.',i4.4)
2000 FORMAT ( a,' -macro ',a,' ',a,i8,2x,i8,' > ',a) 
2100 FORMAT ( a,' -macro ',a,' ',a,i8,      ' > ',a) 
4000 FORMAT ( a, ' -remote -access=127.0.0.1 -port=',I4.4,' > ',a,' &')
4010 FORMAT ( a,' = ',i8)
4020 FORMAT ( 'cd ',a)
4030 FORMAT ( '@',a)
4040 FORMAT ( 'r[',i3.3,'] = ',e20.10)
!
END SUBROUTINE run_mpi_slave
!
!*****7***************************************************************
!
RECURSIVE SUBROUTINE run_mpi_finalize
!
USE mpi
USE run_mpi_mod
USE population
USE errlist_mod
IMPLICIT none
!
INTEGER :: i,j
!
!------       -- Send termination signal to all slave processes
!
IF ( run_mpi_myid == 0 ) THEN
   IF (run_mpi_senddata%use_socket) THEN   ! Closing down, socket variation
         DO i = 1, run_mpi_numprocs - 1
      DO j=1,run_mpi_nprog
         IF(j==run_mpi_nprog) THEN   ! Last program to close, shutdown slave
            run_mpi_senddata%prog_num = -1
         ENDIF
         ier_num = 0
         CALL MPI_SEND ( run_mpi_senddata, 1, run_mpi_data_type, i, 0, &
                         MPI_COMM_WORLD, ier_num )
         ENDDO
      ENDDO
   ELSE                                     ! Closing down, system variation
      DO i = 1, run_mpi_numprocs - 1
   CALL MPI_SEND ( run_mpi_senddata, 1, run_mpi_data_type, i, 0, &
                   MPI_COMM_WORLD, ier_num )
      ENDDO
   ENDIF
ENDIF
!
CALL MPI_FINALIZE ( ier_num )
!
IF(ALLOCATED(node_names)) DEALLOCATE(node_names)
IF(ALLOCATED(node_names)) DEALLOCATE(slave_on_node)
IF(ALLOCATED(node_names)) DEALLOCATE(slave_is_node)
!
END SUBROUTINE run_mpi_finalize
!
END MODULE DIFFEV_MPI_MOD
