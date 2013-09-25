!
!
MODULE run_mpi_mod
!
IMPLICIT none
!
INTEGER                  :: run_mpi_myid
INTEGER                  :: run_mpi_numprocs
!
INTEGER                  :: run_mpi_numjobs
INTEGER                  :: run_mpi_numsent
!
INTEGER, PARAMETER       :: RUN_MPI_COUNT_LOGICAL   = 1
INTEGER, PARAMETER       :: RUN_MPI_COUNT_INTEGER   = 8
INTEGER, PARAMETER       :: RUN_MPI_COUNT_CHARACTER = 4
INTEGER                  :: RUN_MPI_MAXPROG         = 1
INTEGER, DIMENSION(0:2)  :: run_mpi_oldtypes
INTEGER, DIMENSION(0:2)  :: run_mpi_blockcounts
INTEGER, DIMENSION(0:2)  :: run_mpi_offsets
INTEGER                  :: run_mpi_extent
INTEGER                  :: run_mpi_data_type
INTEGER                  :: run_mpi_data_type2
INTEGER                  :: run_mpi_port_base   = 0
INTEGER                  :: run_mpi_nprog       = 0
!
LOGICAL                  :: run_mpi_active      = .false. ! With/without MPI?
!
TYPE run_mpi_type                          ! MPI with types does not work yet
   SEQUENCE
   LOGICAL               :: repeat     !  1
   INTEGER               :: generation !  2
   INTEGER               :: member     !  3
   INTEGER               :: children   !  4
   INTEGER               :: parameters !  5
   INTEGER               :: nindiv     !  6
   INTEGER               :: kid        !  7
   INTEGER               :: indiv      !  8
   INTEGER               :: ierr       !  9
   INTEGER               :: direc_l    ! 10
   INTEGER               :: prog_l     ! 11
   INTEGER               :: mac_l      ! 12
   INTEGER               :: out_l      ! 13
   LOGICAL               :: use_socket ! 14
   INTEGER               :: prog_num   ! 15
   LOGICAL               :: prog_start ! 16
   INTEGER               :: s_remote   ! 17
   INTEGER               :: port       ! 18
   CHARACTER (LEN=240)   :: direc      ! 21 : 260
   CHARACTER (LEN=100)   :: prog       !261 : 360
   CHARACTER (LEN=100)   :: mac        !361 : 460
   CHARACTER (LEN=100)   :: out        !461 : 560
   REAL    ,DIMENSION(:)  , ALLOCATABLE :: trial_values   !  (MAXPOP)
END TYPE run_mpi_type
!
TYPE ( run_mpi_type)     :: run_mpi_senddata
INTEGER, DIMENSION(:), ALLOCATABLE :: run_mpi_send_data
INTEGER                         :: sdl_length
!
CHARACTER(LEN=200), DIMENSION(:), ALLOCATABLE :: prog_entry
INTEGER           , DIMENSION(:), ALLOCATABLE :: socket_id
INTEGER           , DIMENSION(:,:), ALLOCATABLE :: port_id

!
INTEGER                  :: progs_size_of = 0
!
END MODULE run_mpi_mod
