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
!
USE times_mod
!
USE blanks_mod
USE envir_mod
USE errlist_mod
USE gen_mpi_mod
USE mpi_slave_mod
USE prompt_mod
USE variable_mod
!
IMPLICIT none
!
INTEGER                        :: run_mpi_integer_extent
INTEGER                        :: run_mpi_logical_extent
INTEGER                        :: run_mpi_charact_extent
INTEGER                        :: run_mpi_real_extent
!
!                                      Make a local data type to identify node names
INTEGER, PARAMETER       :: GEN_MPI_COUNT_INTEGER   =   1
INTEGER, PARAMETER       :: GEN_MPI_COUNT_CHARACTER =   1
TYPE gen_mpi_type                      ! 
   INTEGER               :: node_l     !  
   CHARACTER (LEN=240)   :: node_name  ! 
END TYPE gen_mpi_type
TYPE ( gen_mpi_type)     :: gen_mpi_senddata
INTEGER                  :: gen_sdl_length
INTEGER, DIMENSION(0:4)  :: gen_mpi_oldtypes
INTEGER, DIMENSION(0:4)  :: gen_mpi_blockcounts
INTEGER, DIMENSION(0:4)  :: gen_mpi_offsets
INTEGER                  :: gen_mpi_data_type
!
CHARACTER (LEN=MPI_MAX_PROCESSOR_NAME) :: node_name   = ' '
INTEGER, DIMENSION(1:MPI_STATUS_SIZE)  :: gen_mpi_status
INTEGER               :: sender        ! Id of slave that answered
INTEGER               :: local_id      ! BUG Patch MPI_ID messy with structure
INTEGER               :: i, j
INTEGER               :: length
INTEGER               :: ierr
LOGICAL               :: success
!
!
gen_mpi_active = .FALSE.
mpi_active = .FALSE.
!
IF(.NOT. (parent_name=='mpiexec' .OR. parent_name=='mpirun' .OR.   &
          parent_name=='orterun')) THEN
!  MPI Does not seem to be active, quietly turn off
   RETURN
ENDIF
!
!PID = lib_f90_getpid()                   ! Get PID and parent processes
!!  pstree command is bugged may create error "/proc/xxxx no such file or directory"
!WRITE(line, '(a,i10, a, i10.10)') 'pstree -s -l -p ', PID, ' > /tmp/discus_suite.',PID
!CALL EXECUTE_COMMAND_LINE(line)
!WRITE(line, '(a, i10.10)') '/tmp/discus_suite.',PID
!CALL oeffne(idef, line, 'old')
!READ(idef, '(a)') string                 ! Get result of pstree command
!CLOSE(UNIT=idef)
!WRITE(line, '(a, i10.10)') 'rm -f /tmp/discus_suite.',PID  ! Remove temporary file
!CALL EXECUTE_COMMAND_LINE(line)
!WRITE(line, '(i10)') PID
!i = 10
!CALL rem_bl(line,i)                               ! remove leading blanks
!ind_pid = INDEX(string, line(1:len_trim(line)))   ! Locate PID
!ind_mpi = INDEX(string, 'mpiexec')                ! Locate mpiexec command
!IF(ind_mpi==0 .OR. ind_pid < ind_mpi) THEN        ! No mpiexec as parent process
!   gen_mpi_active = .FALSE.
!   mpi_active = .FALSE.
!!  MPI Does not seem to be active, quietly turn off
!   RETURN
!ENDIF
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
!  Get the rank of this process, store in gen_mpi_myid
!
CALL MPI_COMM_RANK (MPI_COMM_WORLD, gen_mpi_myid,     ier_num)
!
IF ( ier_num /= 0 ) THEN
   ier_msg(1) = 'MPI SYSTEM did not return RANK'
   WRITE(ier_msg(2),3000) ier_num
   ier_num = -22
   ier_typ = ER_APPL
   RETURN
ENDIF 
!
!  Get the size of the distribution, store in gen_mpi_numprocs
!
CALL MPI_COMM_SIZE (MPI_COMM_WORLD, gen_mpi_numprocs, ier_num)
!
IF ( ier_num /= 0 ) THEN
   ier_msg(1) = 'MPI SYSTEM did not return SIZE'
   WRITE(ier_msg(2),3000) ier_num
   ier_num = -22
   ier_typ = ER_APPL
   RETURN
ENDIF
!
!  TEST for MPIEXEC was used as parent process but insufficient CPUs
!
IF ( gen_mpi_numprocs < 2 ) THEN
!   IF(ind_mpi==0 .OR. ind_pid < ind_mpi) THEN        ! No mpiexec as parent process
!      gen_mpi_active = .FALSE.
!      mpi_active = .FALSE.
!!     MPI Does not seem to be active, quietly turn off
!      RETURN
!   ELSE
!     MPI !,0)was used to start with only one process
      CALL MPI_FINALIZE(ier_num)
      gen_mpi_active = .FALSE.
      mpi_active = .FALSE.
!
      ier_msg(1) = 'MPI SYSTEM returned one CPU   '
      ier_msg(2) = 'SUITE  must be started with   '
      ier_msg(3) = 'mpiexec -n X diffev; X >= 2   '
!
      ier_num = -22
      ier_typ = ER_APPL
      RETURN
!   ENDIF
ENDIF
!
! For future use with MPI_TYPE_...
!
CALL MPI_TYPE_EXTENT ( MPI_INTEGER,   run_mpi_integer_extent, ier_num )
CALL MPI_TYPE_EXTENT ( MPI_LOGICAL,   run_mpi_logical_extent, ier_num )
CALL MPI_TYPE_EXTENT ( MPI_CHARACTER, run_mpi_charact_extent, ier_num )
CALL MPI_TYPE_EXTENT ( MPI_REAL     , run_mpi_real_extent   , ier_num )
!
! Build a local mpi data structure, just to get node names
!
gen_mpi_offsets(0)     = 0
gen_mpi_oldtypes(0)    = MPI_INTEGER
gen_mpi_blockcounts(0) = GEN_MPI_COUNT_INTEGER*run_mpi_integer_extent
!
gen_mpi_offsets(1)     = gen_mpi_offsets(0) +  run_mpi_integer_extent*GEN_MPI_COUNT_INTEGER
gen_mpi_oldtypes(1)    = MPI_CHARACTER
gen_mpi_blockcounts(1) = GEN_MPI_COUNT_CHARACTER*run_mpi_charact_extent
!
CALL MPI_TYPE_STRUCT ( 2, gen_mpi_blockcounts, gen_mpi_offsets,  &
     gen_mpi_oldtypes, gen_mpi_data_type, ier_num )
CALL MPI_TYPE_COMMIT ( gen_mpi_data_type, ier_num)
!
!
! For future use with MPI_TYPE_...
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
!write(*,*) '############## data type , myid', run_mpi_data_type,gen_mpi_myid, ier_num
!
!WRITE(*,4000)
!
gen_mpi_active = .true.
!
!socket_status = PROMPT_OFF  ! Turn off socket responses
!
mpi_active = .true.
!
!  Build a list of nodes, node names, and assignment of the processes onto these nodes
!
local_id = gen_mpi_myid
IF(gen_mpi_myid==0)  THEN  !   MASTER
   IF(ALLOCATED(node_names)) DEALLOCATE(node_names)
   IF(ALLOCATED(node_names)) DEALLOCATE(slave_is_node)
   ALLOCATE(node_names    (1:gen_mpi_numprocs))
   ALLOCATE(slave_is_node (1:gen_mpi_numprocs))
   node_names    (:) = ' '
   slave_is_node (:) = -1
   NUM_NODE = 0
   gen_sdl_length = 1
   DO i=1, gen_mpi_numprocs-1
!
      CALL MPI_SEND ( gen_mpi_senddata, 1, gen_mpi_data_type, i, i, &
                      MPI_COMM_WORLD, ier_num )
   ENDDO
   DO i=1, gen_mpi_numprocs-1
      CALL MPI_RECV ( gen_mpi_senddata, 1, gen_mpi_data_type, MPI_ANY_SOURCE, &
                      MPI_ANY_TAG, MPI_COMM_WORLD, gen_mpi_status, ier_num )
      gen_mpi_myid = local_id ! BUG PATCH MPI_ID gets messed up by receive with structure?????
      sender = gen_mpi_status(MPI_SOURCE)     ! Identify the slave
      success = .FALSE.
      search: DO j=1,NUM_NODE
         IF(gen_mpi_senddata%node_name == node_names(j)) THEN
            success = .TRUE.
            slave_is_node(sender) = j
!           slave_per_node(j)     = slave_per_node(j) + 1
            EXIT search
         ENDIF
      ENDDO search
      IF(.NOT.success) THEN
         NUM_NODE = NUM_NODE + 1
         node_names(NUM_NODE) = ' '
         node_names(NUM_NODE)(1:gen_mpi_senddata%node_l) = &
              gen_mpi_senddata%node_name(1:gen_mpi_senddata%node_l)
         slave_is_node(sender) = NUM_NODE
      ENDIF
   ENDDO
   var_val(VAR_NUM_NODES) = NUM_NODE
!
ELSE
   CALL MPI_PROBE( MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, gen_mpi_status, ierr) ! Querry incomming size
   CALL MPI_GET_COUNT(gen_mpi_status, gen_mpi_data_type, gen_sdl_length, ierr)                  ! Determine size
!
!  Now receive incomming message
!
   CALL MPI_RECV ( gen_mpi_senddata, gen_sdl_length, gen_mpi_data_type, MPI_ANY_SOURCE, &
                   MPI_ANY_TAG, MPI_COMM_WORLD, gen_mpi_status, ier_num )
   gen_mpi_myid = local_id ! BUG PATCH MPI_ID gets messed up by receive with structure?????
   CALL MPI_GET_PROCESSOR_NAME ( node_name, length, ier_num)
   gen_mpi_senddata%node_name  = node_name(1:MIN(LEN(gen_mpi_senddata%node_name),LEN_TRIM(node_name)))
   gen_mpi_senddata%node_l     = length

   CALL MPI_SEND ( gen_mpi_senddata, 1, gen_mpi_data_type, 0, 0, &
                   MPI_COMM_WORLD, ier_num )
   var_val(VAR_NUM_NODES) = 0
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
!  sends them to gen_mpi_numprocs -1 slaves
!
USE diffev_allocate_appl
USE diffev_random
USE diffev_distrib_mod
USE mpi
USE population
USE run_mpi_mod
USE gen_mpi_mod
USE errlist_mod
USE lib_errlist_func
USE precision_mod
USE prompt_mod
USE support_mod
!
IMPLICIT none
!
INTEGER, PARAMETER    :: one_indiv = 1 ! Send NINDIV=1 to distrib if repeat==false
!
INTEGER, DIMENSION(1:MPI_STATUS_SIZE) :: run_mpi_status
!
CHARACTER (LEN=PREC_LSTRING)  :: send_direc    ! working directory
INTEGER               :: send_direc_l  ! working directory length
INTEGER               :: sender        ! Id of slave that answered
INTEGER               :: i,j,kid
INTEGER               :: nseeds        ! Number of seeds for randum numbers
INTEGER               :: run_mpi_numsent  ! Number of jobs sent out
INTEGER               :: run_mpi_numjobs  ! Number of initial jobs
!INTEGER               :: nprog         ! number of different program/mac combinations
!LOGICAL               :: prog_start    ! External program needs to be started
!
!LOGICAL               :: prog_exist    ! program/macro combination exists in data base
INTEGER               :: local_id      ! BUG Patch MPI_ID messy with structure
INTEGER               :: inode         ! Current slave process
INTEGER               :: slave         ! Current slave
INTEGER               :: num_hand = 0  ! Number of loops handing out jobs
INTEGER               :: numrec   = 0  ! Number of RJobs received
INTEGER               :: numtasks = 0  ! Number of loops handing out jobs
!
local_id = gen_mpi_myid  ! BUG PATCH MPI_ID gets messed up be receive with structure
!
sdl_length = 1 !580! + 200
!
! Test if program / macro combination exist /not exists ; sockets only
!
!IF (run_mpi_senddata%use_socket) THEN
!   prog_exist = .false.
!   find_prog_entry: DO i=1,run_mpi_nprog        ! Search existing entries
!      IF(prog_entry(i) == run_mpi_senddata%prog(1:run_mpi_senddata%prog_l) &
!                          // ' '//                                         &
!                          run_mpi_senddata%mac(1:run_mpi_senddata%mac_l)) THEN
!         prog_exist = .true.
!         prog_start = .false.                   ! Program should already be running
!         run_mpi_senddata%prog_num = i          ! This is program/macro identifier
!         run_mpi_senddata%s_remote = socket_id(i) ! This is program/macro identifier
!         EXIT find_prog_entry
!      ENDIF
!   ENDDO find_prog_entry
!!
!   IF(.NOT. prog_exist) THEN                    ! New program to be started by slaves
!      IF(run_mpi_nprog==RUN_MPI_MAXPROG) THEN   ! Need more space
!         nprog = RUN_MPI_MAXPROG + 2            ! increment by two programs
!         CALL alloc_socket_nprogs ( nprog, gen_mpi_numprocs)
!      ELSE                                      ! Sufficient space
!         run_mpi_nprog = run_mpi_nprog + 1      ! Increment no of known prog/mac entries
!         prog_entry(run_mpi_nprog) = run_mpi_senddata%prog(1:run_mpi_senddata%prog_l) &
!                                     // ' '//                                         &
!                                     run_mpi_senddata%mac(1:run_mpi_senddata%mac_l)
!      ENDIF
!      run_mpi_senddata%prog_num = run_mpi_nprog ! This is program/macro identifier
!      prog_start = .true.                       ! Program needs to be started
!      port_id  (run_mpi_nprog,:) = 0            ! Initial port number
!      socket_id(run_mpi_nprog)   = 0            ! Initial socket ID
!   ENDIF
!ELSE
!   prog_start = .true.                          ! Always start non socket program
!ENDIF
!
!run_mpi_senddata%prog_start = prog_start        ! Copy start flag into send structure
!
CALL do_cwd ( send_direc, send_direc_l )        ! Get current working directory
run_mpi_senddata%direc_l = send_direc_l         ! Copy directory into send structure
run_mpi_senddata%direc   = send_direc(1:MIN(send_direc_l,200))
run_mpi_senddata%ierr_msg_l = 80
run_mpi_senddata%ierr_msg_n =  7
!
!
IF(pop_gen /= lastgen) THEN                   ! New GENERATION , new job distribution
   IF(ALLOCATED(kid_at_node  )) DEALLOCATE(kid_at_node)
   IF(ALLOCATED(node_has_kids)) DEALLOCATE(node_has_kids)
   IF(ALLOCATED(node_max_kids)) DEALLOCATE(node_max_kids)
   ALLOCATE(kid_at_node (1:pop_c))
   ALLOCATE(node_has_kids(1:NUM_NODE,0:pop_c))
   ALLOCATE(node_max_kids(1:NUM_NODE))
   kid_at_node (:)    = 0
   node_has_kids(:,:) = 0
   node_max_kids(:)   = 0
!
!     Create even load of kids per node
!
   kid = 0
   spread: DO
      DO j=NUM_NODE,1,-1
         kid = kid + 1
         node_max_kids(j) = node_max_kids(j) + 1
         IF(kid==pop_c) EXIT spread
      ENDDO
   ENDDO spread
   run_mpi_senddata%l_first_job = .TRUE.
ELSE
   IF(.NOT.ALLOCATED(kid_at_node  )) THEN
      ALLOCATE(kid_at_node (1:pop_c))
      kid_at_node (:)    = 0
   ENDIF
   IF(.NOT.ALLOCATED(node_has_kids  )) THEN
      ALLOCATE(node_has_kids(1:NUM_NODE,0:pop_c))
      node_has_kids (:,:)  = 0
   ENDIF
   IF(.NOT.ALLOCATED(node_max_kids)) THEN
      ALLOCATE(node_max_kids(1:NUM_NODE))
      node_max_kids(:)    = 0
   ENDIF
   run_mpi_senddata%l_first_job = .FALSE.
ENDIF
!
IF(ALLOCATED(kid_at_indiv) ) DEALLOCATE(kid_at_indiv)
IF(ALLOCATED(node_finished) ) DEALLOCATE(node_finished)
ALLOCATE(kid_at_indiv(1:pop_c))
ALLOCATE(node_finished(1:NUM_NODE))
kid_at_indiv(: ) = 0
node_finished(:) = .FALSE.
!
run_mpi_numsent = 0                             ! No jobs sent yet
IF(run_mpi_senddata%repeat) THEN
   run_mpi_numjobs = MIN ( gen_mpi_numprocs - 1, pop_c * run_mpi_senddata%nindiv )
   numtasks = pop_c * run_mpi_senddata%nindiv
ELSE
   run_mpi_numjobs = MIN ( gen_mpi_numprocs - 1, pop_c                           )
   numtasks = pop_c
ENDIF
!
run_mpi_senddata%l_get_state = l_get_random_state  ! Inquire random number status
!
!  Start initial jobs
!
slave = 1
!
initial:DO                                           !  Start the intial jobs
   ier_num = 0
!
   run_mpi_senddata%kid    = 0                      ! Will be > 0, if distrib is OK
   run_mpi_senddata%indiv  = 0
   inode = slave_is_node(slave)                     ! Current slave is on this node
      IF(run_mpi_senddata%repeat) THEN              ! Parallel computing of indivs requested
         IF(pop_gen /= lastgen) THEN                ! New GENERATION , new job distribution
   IF(.NOT.node_finished(inode)) THEN               ! This node has kids/indivs available
            CALL distrib_even(run_mpi_senddata%kid, run_mpi_senddata%indiv, &
                 NUM_NODE,  pop_c, run_mpi_senddata%nindiv, inode,          &
                 kid_at_indiv, kid_at_node,                                 &
                 node_has_kids, node_max_kids, node_finished)
         ENDIF
   IF(ALL(node_finished)) THEN                     ! All nodes have been populated
      EXIT initial
   ENDIF
         ELSE                                       ! Same generation place kid onto previous node
            CALL distrib_preve(run_mpi_senddata%kid, run_mpi_senddata%indiv, &
                 NUM_NODE,  pop_c, run_mpi_senddata%nindiv, inode,           &
                 kid_at_indiv, node_has_kids, node_finished)
         ENDIF
      ELSE                                          ! Serial distribution of indivs
         IF(pop_gen /= lastgen) THEN                ! New GENERATION , new job distribution
            CALL distrib_sequential(run_mpi_senddata%kid, run_mpi_senddata%indiv, &
                 NUM_NODE, pop_c, run_mpi_senddata%nindiv, numtasks, inode,       &
                 run_mpi_numsent, kid_at_indiv, kid_at_node, node_has_kids)
         ELSE                                       ! Same generation place kid onto previous node
            CALL distrib_preve(run_mpi_senddata%kid, run_mpi_senddata%indiv, &
                 NUM_NODE,  pop_c, one_indiv, inode,           &
                 kid_at_indiv, node_has_kids, node_finished)
      ENDIF
   ENDIF
   IF(run_mpi_senddata%kid> 0) THEN   ! Proper assignment start the job
!      IF(run_mpi_senddata%prog_start) THEN     ! Program needs to be started, increment port no
!         run_mpi_senddata%port     = 2000 + MOD(gen_mpi_numprocs*run_mpi_senddata%prog_num + slave,3600)
!      ELSE                                     ! Program is running use old port no
!         run_mpi_senddata%port     = port_id  (run_mpi_senddata%prog_num,slave)
!      ENDIF
      DO j=1,pop_dimx                          ! Encode current trial values
         run_mpi_senddata%trial_names (j) = pop_name(j                  ) ! Takes value for kid that answered
         run_mpi_senddata%trial_values(j) = pop_t(j,run_mpi_senddata%kid) ! Takes value for kid that answered
      ENDDO
!
      sdl_length = 1! 580! + 200
      CALL MPI_SEND ( run_mpi_senddata, 1, run_mpi_data_type, slave, slave, &
                      MPI_COMM_WORLD, ier_num )
!
      run_mpi_numsent = run_mpi_numsent + 1
!
      slave = slave + 1
      IF(run_mpi_numsent >  run_mpi_numjobs .OR. slave > gen_mpi_numprocs-1) EXIT initial
   ELSE
      WRITE(output_io,'(a)') 'DISTRIBUTION failed'
      WRITE(output_io,'(a)') 'Please document and report to author '
      STOP
   ENDIF
   IF(run_mpi_numsent == numtasks) EXIT initial
ENDDO initial
!
!------       --Receive results and hand out new jobs
!
num_hand = 0
numrec   = 0
rec_hand: DO   ! i = 1, num_hand
   CALL MPI_RECV ( run_mpi_senddata, 1, run_mpi_data_type, MPI_ANY_SOURCE, &
                   MPI_ANY_TAG, MPI_COMM_WORLD, run_mpi_status, ier_num )
   gen_mpi_myid = local_id ! BUG PATCH MPI_ID gets messed up by receive with structure?????
!
   numrec = numrec + 1
!
   sender = run_mpi_status(MPI_SOURCE)     ! Identify the slave
   IF(run_mpi_senddata%ierr /=0 ) THEN
      ier_num = run_mpi_senddata%ierr
      ier_typ = run_mpi_senddata%ierr_typ
      ier_msg = run_mpi_senddata%ierr_msg
      call errlist
      ier_num = 0
      ier_typ = 0
      CALL diffev_error_macro              ! Write a recovery macro
      ier_msg(1) = 'A slave program exited with error message'
      WRITE(ier_msg(2), 2000)  sender, run_mpi_senddata%ierr, run_mpi_senddata%ierr_typ
      ier_num = -26
      ier_typ = ER_APPL
2000 FORMAT('Error message from slave ',I4,' is ',i4,' ',i4)
      EXIT rec_hand
   ENDIF
!   IF(run_mpi_senddata%use_socket ) THEN   ! SOCKET option is active
!       socket_id(run_mpi_senddata%prog_num) = run_mpi_senddata%s_remote  ! Store socket info in data base
!       port_id  (run_mpi_senddata%prog_num,sender) = run_mpi_senddata%port  ! Store port info in data base
!   ENDIF
!
   IF(l_get_random_state==-1) THEN             ! Update random number status
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
   slave = sender
   inode = slave_is_node(slave)
   IF(run_mpi_senddata%repeat)   THEN                                  ! Even distribution, use a kid on this node
      IF(pop_gen/=lastgen) THEN                           ! New generation, create new distribution
            CALL distrib_even(run_mpi_senddata%kid, run_mpi_senddata%indiv, &
                 NUM_NODE,  pop_c, run_mpi_senddata%nindiv, inode,         &
                 kid_at_indiv, kid_at_node,                               &
                 node_has_kids, node_max_kids, node_finished)
         ELSE                                               ! Same Generation use previous distribution
            CALL distrib_preve(run_mpi_senddata%kid, run_mpi_senddata%indiv, &
                 NUM_NODE,  pop_c, run_mpi_senddata%nindiv, inode,          &
                 kid_at_indiv, node_has_kids, node_finished)
         ENDIF
   ELSE                                                  ! Sequential distribution, take next job
      IF(pop_gen/=lastgen) THEN                           ! New generation, create new distribution
         CALL distrib_sequential(run_mpi_senddata%kid, run_mpi_senddata%indiv, &
                 NUM_NODE, pop_c, run_mpi_senddata%nindiv, numtasks, inode,       &
                 run_mpi_numsent, kid_at_indiv, kid_at_node, node_has_kids)
!        CALL distrib_sequential(run_mpi_senddata%kid, run_mpi_senddata%indiv, &
!             pop_c, numtasks, inode, run_mpi_numsent, kid_at_indiv, kid_at_node)
      ELSE                                               ! Same Generation use previous distribution
            CALL distrib_preve(run_mpi_senddata%kid, run_mpi_senddata%indiv, &
                 NUM_NODE,  pop_c, one_indiv, inode,          &
                 kid_at_indiv, node_has_kids, node_finished)
      ENDIF
   ENDIF
!
      IF(run_mpi_senddata%kid>0) THEN
         DO j=1,pop_dimx                          ! Encode current trial values
            run_mpi_senddata%trial_names (j) = pop_name(j                  ) ! Takes value for kid that answered
            run_mpi_senddata%trial_values(j) = pop_t(j,run_mpi_senddata%kid) ! Takes value for kid that answered
         ENDDO
!
      CALL MPI_SEND ( run_mpi_senddata, 1, run_mpi_data_type, sender, slave, &
                      MPI_COMM_WORLD, ier_num )
!
      run_mpi_numsent = run_mpi_numsent + 1
!
   ENDIF
   IF(numrec == numtasks) EXIT rec_hand
ENDDO rec_hand
!
IF(ier_num == -26) THEN   ! Fatal error occured, wait for remaining jobs
   DO i=1, run_mpi_numjobs-1
         CALL MPI_RECV ( run_mpi_senddata, 1, run_mpi_data_type, MPI_ANY_SOURCE, &
                   MPI_ANY_TAG, MPI_COMM_WORLD, run_mpi_status, ier_num )
   ENDDO
   ier_num = -26 ! Reinstate the error message
ENDIF
gen_mpi_myid = local_id  ! BUG PATCH MPI_ID gets messed up by receive with structure?????
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
USE diffev_setup_cost_mod
!use diffev_allocate_appl
!
USE gen_mpi_mod
USE set_sub_generic_mod
USE errlist_mod
USE mpi_slave_mod
USE lib_length
USE prompt_mod
!USE sockets_mod
USE precision_mod
USE variable_mod
!
IMPLICIT none
!
CHARACTER (LEN=PREC_LSTRING)   :: line
CHARACTER (LEN=PREC_LSTRING)   :: output
!
INTEGER, DIMENSION(1:MPI_STATUS_SIZE) :: run_mpi_status
!
!INTEGER  :: j
!INTEGER  :: il
INTEGER  :: ierr
INTEGER  :: job_l
INTEGER  :: output_l
INTEGER  :: port  = 2000
!
!INTEGER, PARAMETER :: seconds = 1       ! Wait time for application to start up
!
!INTEGER  :: system
!INTEGER  :: socket_close
!INTEGER  :: socket_connect
!INTEGER  :: socket_get
!INTEGER  :: socket_send
!
!real(kind=PREC_DP), dimension(:), allocatable :: trial_values
INTEGER  :: local_id  ! BUG PATCH MPI_ID gets messed up by receive with structure?????
!
ierr = 0
!call alloc_senddata(14, 1)
!
local_id = gen_mpi_myid ! BUG PATCH MPI_ID gets messed up by receive with structure?????
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
   gen_mpi_myid = local_id ! BUG PATCH MPI_ID gets messed up by receive with structure?????
!
   s_remote = run_mpi_senddata%s_remote
   port     = run_mpi_senddata%port
!
!------       --- If TAG is zero exit, else compute
!
   tag_exit: IF (run_mpi_status (MPI_TAG) == 0) then   ! End of diffev program
!      IF ( run_mpi_senddata%use_socket ) THEN          ! SOCKET option is active
!         line  = 'exit'
!         job_l = 4
!         ierr = socket_send    (s_remote, line, job_l) ! End application program
!         CALL sleep(seconds)                           ! wait for application to close
!         ierr = socket_close   (s_remote)              ! Close socket
!         lremote = .false.
!         IF(run_mpi_senddata%prog_num == -1) THEN         ! This is the last program to close
!            lremote = .false.
!            EXIT slave                                    ! now we can end slave
!         ENDIF
!      ELSE
         lremote = .false.
         EXIT slave                                    ! now we can end slave
!      ENDIF
   ELSE  tag_exit                                      ! Normal operation
!
!------       Define output redirection
!
      out_null: IF ( run_mpi_senddata%out(1:run_mpi_senddata%out_l) == '/dev/null' ) THEN   ! no output desired
         output   = '/dev/null'
         output_l = 9
      ELSE out_null                                          ! redirect to a log file
!         out_socket: IF ( run_mpi_senddata%use_socket ) THEN ! SOCKET option is active
!               WRITE(output,1100) run_mpi_senddata%out(1:run_mpi_senddata%out_l),&
!                             gen_mpi_myid
!               output_l = run_mpi_senddata%out_l + 5
!         ELSE out_socket                                     ! No socket active
            out_repeat:IF ( run_mpi_senddata%repeat ) THEN   ! call with kid AND indiv
               WRITE(output,1000) run_mpi_senddata%out(1:run_mpi_senddata%out_l),&
                             run_mpi_senddata%kid   , run_mpi_senddata%indiv
               output_l = run_mpi_senddata%out_l + 10
            ELSE out_repeat                                  ! call with kid only
               WRITE(output,1100) run_mpi_senddata%out(1:run_mpi_senddata%out_l),&
                             run_mpi_senddata%kid   
               output_l = run_mpi_senddata%out_l + 5
            ENDIF out_repeat
!         ENDIF out_socket
      ENDIF out_null
!
!------       ----- Here we will call DISCUS/KUPLOT 
!
!      use_socket: IF ( run_mpi_senddata%use_socket ) THEN   ! SOCKET option is active
!         IF ( run_mpi_senddata%prog_start ) THEN            ! server needs to be started
!            WRITE(line,4000) run_mpi_senddata%prog (1:run_mpi_senddata%prog_l ), &
!                 port, output(1:output_l)
!            job_l = run_mpi_senddata%prog_l + 42 + output_l
!            ierr =  system ( line(1:job_l))                 ! Start with -remote ... optios
!            CALL sleep(seconds)                             ! wait for application to start
!            ierr = socket_connect (s_remote, '127.0.0.1',  9, port) ! Connect to socket
!            ierr = socket_get     (s_remote, line, il)      ! Get socket answer
!!
!!           Send standard variable names and values
!!
!            ierr = socket_send    (s_remote, 'variable integer, port      ', 28)
!            CALL   socket_wait
!            ierr = socket_send    (s_remote, 'variable integer, myid      ', 28)
!            CALL   socket_wait
!            ierr = socket_send    (s_remote, 'variable integer, generation', 28)
!            CALL   socket_wait
!            ierr = socket_send    (s_remote, 'variable integer, member    ', 28)
!            CALL   socket_wait
!            ierr = socket_send    (s_remote, 'variable integer, children  ', 28)
!            CALL   socket_wait
!            ierr = socket_send    (s_remote, 'variable integer, parameters', 28)
!            CALL   socket_wait
!            ierr = socket_send    (s_remote, 'variable integer, nindiv    ', 28)
!            CALL   socket_wait
!            ierr = socket_send    (s_remote, 'variable integer, indiv     ', 28)
!            CALL   socket_wait
!            ierr = socket_send    (s_remote, 'variable integer, kid       ', 28)
!            CALL   socket_wait
!            lremote = .true.         ! remote Socket program is running
!         ENDIF
!!
!!        Socket should  be active at this point
!         WRITE(line, 4010) 'port      ',port                        ! Current member size
!         ierr = socket_send    (s_remote, line, 21)
!         CALL   socket_wait
!         WRITE(line, 4010) 'myid      ',gen_mpi_myid                ! Current member size
!         ierr = socket_send    (s_remote, line, 21)
!         CALL   socket_wait
!         WRITE(line, 4010) 'member    ',run_mpi_senddata%member     ! Current member size
!         ierr = socket_send    (s_remote, line, 21)
!         CALL   socket_wait
!         WRITE(line, 4010) 'children  ',run_mpi_senddata%children   ! Current children size
!         ierr = socket_send    (s_remote, line, 21)
!         CALL   socket_wait
!         WRITE(line, 4010) 'parameters',run_mpi_senddata%parameters ! Current parameter size
!         ierr = socket_send    (s_remote, line, 21)
!         CALL   socket_wait
!         WRITE(line, 4010) 'nindiv    ',run_mpi_senddata%nindiv     ! Current member size
!         ierr = socket_send    (s_remote, line, 21)
!         CALL   socket_wait
!         WRITE(line, 4020) run_mpi_senddata%direc(1:run_mpi_senddata%direc_l) ! Current directory
!         job_l = len_str(line)
!         ierr = socket_send    (s_remote, line, job_l)
!         CALL   socket_wait
!         WRITE(line, 4010) 'generation',run_mpi_senddata%generation !Current Generation
!         ierr = socket_send    (s_remote, line, 21)
!         CALL   socket_wait
!         WRITE(line, 4010) 'indiv     ',run_mpi_senddata%indiv      ! Current indiv numner
!         ierr = socket_send    (s_remote, line, 21)
!         CALL   socket_wait
!         WRITE(line, 4010) 'kid       ',run_mpi_senddata%kid        ! Current kid numner
!         ierr = socket_send    (s_remote, line, 21)
!         CALL   socket_wait
!         DO j=1,run_mpi_senddata%parameters                         ! Current trial values
!            WRITE(line, 4040) 200+j,run_mpi_senddata%trial_values(j)
!!            line = ' ' !WORK WORK WORK 
!            ierr = socket_send    (s_remote, line, 29)
!            CALL   socket_wait
!         ENDDO
!!        Now send the macro
!         WRITE(line, 4030) run_mpi_senddata%mac  (1:run_mpi_senddata%mac_l  )
!         job_l = len_str(line)
!!         ierr = socket_send    (s_remote, line, job_l)     ! Send macro
!         CALL   socket_wait
!      ELSE use_socket                                      ! explicitely start DISCUS/KUPLOT
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
!deallocate(run_mpi_senddata%trial_values)
!allocate(run_mpi_senddata%trial_values(1:ubound(run_mpi_senddata%trial_values,1)))
!write(*,*) ' SLAVE ', gen_mpi_myid, allocated(run_mpi_senddata%trial_values)
!write(*,*) ' SLAVE ', gen_mpi_myid, ubound(run_mpi_senddata%trial_values)
!write(*,*) ' SLAVE ', 'allocated'
!trial_values = run_mpi_senddata%trial_values
!write(*,*) ' SLAVE1', gen_mpi_myid,                         trial_values(1)
!write(*,*) ' SLAVE2', gen_mpi_myid,        run_mpi_senddata%trial_values(1)
      CALL p_execute_cost( run_mpi_senddata%repeat,                          &
                           LEN(run_mpi_senddata%prog),                       &
                           run_mpi_senddata%prog , run_mpi_senddata%prog_l , &
                           LEN(run_mpi_senddata%mac),                        &
                           run_mpi_senddata%mac  , run_mpi_senddata%mac_l  , &
                           LEN(run_mpi_senddata%direc),                      &
                           run_mpi_senddata%direc, run_mpi_senddata%direc_l, &
                           run_mpi_senddata%kid  , run_mpi_senddata%indiv  , &
                           run_mpi_senddata%n_rvalue_i, run_mpi_senddata%n_rvalue_o, &
                           run_mpi_senddata%RUN_MPI_MAXRVALUE,     &
                           run_mpi_senddata%rvalue, run_mpi_senddata%l_rvalue,     &
                           LEN(output),                                      &
                           output                , output_l ,                      &
                           run_mpi_senddata%RUN_MPI_MAX_FLAGS, run_mpi_senddata%global_flags, &
                           run_mpi_senddata%generation, run_mpi_senddata%member,   &
                           run_mpi_senddata%children, run_mpi_senddata%parameters, &
                                                   run_mpi_senddata%nindiv  , &
                           run_mpi_senddata%trial_names ,                          &
                           run_mpi_senddata%trial_values, run_mpi_senddata%RUN_MPI_COUNT_TRIAL,     &
                           run_mpi_senddata%l_get_state,                           &
                           run_mpi_senddata%nseeds, run_mpi_senddata%seeds,        &
                           run_mpi_senddata%l_first_job,                           &
                           run_mpi_senddata%ierr, run_mpi_senddata%ierr_typ,        &
                           run_mpi_senddata%ierr_msg_l, run_mpi_senddata%ierr_msg_n,&
                           run_mpi_senddata%ierr_msg )
!   ENDIF use_socket
!
!  Answer back to master)
!
!  run_mpi_senddata%ierr     = ierr
!  run_mpi_senddata%ierr_typ = ierr_typ
   run_mpi_senddata%s_remote = s_remote
   run_mpi_senddata%port     = port
!
   CALL MPI_SEND ( run_mpi_senddata, 1, run_mpi_data_type, 0, 0, &
                   MPI_COMM_WORLD, ier_num )
   ENDIF tag_exit
ENDDO slave
gen_mpi_myid = local_id ! BUG PATCH MPI_ID gest messed up by receive with structure?????
!
1000 FORMAT ( a,'.',i4.4,'.',i4.4)
1100 FORMAT ( a,'.',i4.4)
2000 FORMAT ( a,' -macro ',a,' ',a,i8,2x,i8,' > ',a) 
2100 FORMAT ( a,' -macro ',a,' ',a,i8,      ' > ',a) 
!4000 FORMAT ( a, ' -remote -access=127.0.0.1 -port=',I4.4,' > ',a,' &')
!4010 FORMAT ( a,' = ',i8)
!4020 FORMAT ( 'cd ',a)
!4030 FORMAT ( '@',a)
!4040 FORMAT ( 'r[',i3.3,'] = ',e20.10)
!
END SUBROUTINE run_mpi_slave
!
!*****7***************************************************************
!
RECURSIVE SUBROUTINE run_mpi_finalize
!
USE mpi
USE run_mpi_mod
USE mpi_slave_mod
USE gen_mpi_mod
!
USE errlist_mod
IMPLICIT none
!
INTEGER :: i
!
IF(mpi_active) THEN
!
!------       -- Send termination signal to all slave processes
!
IF ( gen_mpi_myid == 0 ) THEN
!   IF (run_mpi_senddata%use_socket) THEN   ! Closing down, socket variation
!         DO i = 1, gen_mpi_numprocs - 1
!      DO j=1,run_mpi_nprog
!         IF(j==run_mpi_nprog) THEN   ! Last program to close, shutdown slave
!            run_mpi_senddata%prog_num = -1
!         ENDIF
!         ier_num = 0
!         CALL MPI_SEND ( run_mpi_senddata, 1, run_mpi_data_type, i, 0, &
!                         MPI_COMM_WORLD, ier_num )
!         ENDDO
!      ENDDO
!   ELSE                                     ! Closing down, system variation
      DO i = 1, gen_mpi_numprocs - 1
   CALL MPI_SEND ( run_mpi_senddata, 1, run_mpi_data_type, i, 0, &
                   MPI_COMM_WORLD, ier_num )
      ENDDO
!   ENDIF
ENDIF
!
CALL MPI_FINALIZE ( ier_num )
ENDIF
!
IF(ALLOCATED(node_names)) DEALLOCATE(node_names)
!
IF(ALLOCATED(node_names)) DEALLOCATE(slave_is_node)
!
END SUBROUTINE run_mpi_finalize
!
END MODULE DIFFEV_MPI_MOD
