!
!
MODULE run_mpi_mod
!
IMPLICIT none
!
INTEGER                  :: run_mpi_myid
INTEGER                  :: run_mpi_numprocs
!
INTEGER, PARAMETER       :: RUN_MPI_COUNT_INTEGER   =  16
INTEGER, PARAMETER       :: RUN_MPI_COUNT_LOGICAL   =   4
INTEGER, PARAMETER       :: RUN_MPI_COUNT_CHARACTER = 540
INTEGER, PARAMETER       :: RUN_MPI_COUNT_REAL      =   2
INTEGER, PARAMETER       :: RUN_MPI_COUNT_TRIAL     = 200
INTEGER                  :: RUN_MPI_MAXPROG         = 1
INTEGER, DIMENSION(0:4)  :: run_mpi_oldtypes
INTEGER, DIMENSION(0:4)  :: run_mpi_blockcounts
INTEGER, DIMENSION(0:4)  :: run_mpi_offsets
INTEGER                  :: run_mpi_extent
INTEGER                  :: run_mpi_data_type
INTEGER                  :: run_mpi_data_type2
INTEGER                  :: run_mpi_port_base   = 0
INTEGER                  :: run_mpi_nprog       = 0
!
LOGICAL                  :: run_mpi_active      = .false. ! With/without MPI?
!
TYPE run_mpi_type                          ! MPI with types does not work yet
   INTEGER               :: generation !  2  !  1
   INTEGER               :: member     !  3  !  2
   INTEGER               :: children   !  4  !  3
   INTEGER               :: parameters !  5  !  4
   INTEGER               :: nindiv     !  6  !  5
   INTEGER               :: kid        !  7  !  6
   INTEGER               :: indiv      !  8  !  7
   INTEGER               :: ierr       !  9  !  8
   INTEGER               :: direc_l    ! 10  !  9
   INTEGER               :: prog_l     ! 11  ! 10
   INTEGER               :: mac_l      ! 12  ! 11
   INTEGER               :: out_l      ! 13  ! 12
   INTEGER               :: prog_num   ! 15  ! 13
   INTEGER               :: s_remote   ! 17  ! 14
   INTEGER               :: port       ! 18  ! 15
   INTEGER               :: spacer     ! 18  ! 15
   LOGICAL               :: repeat     !  1  ! 16
   LOGICAL               :: use_socket ! 14  ! 17
   LOGICAL               :: prog_start ! 16  ! 18
   LOGICAL               :: l_rvalue   ! 19  ! 19
   CHARACTER (LEN=240)   :: direc      ! 21 : 260
   CHARACTER (LEN=100)   :: prog       !261 : 360
   CHARACTER (LEN=100)   :: mac        !361 : 460
   CHARACTER (LEN=100)   :: out        !461 : 560
   REAL                  :: rvalue     !561 : 580
   REAL                  :: rvalue2    !561 : 580
   REAL,DIMENSION(1:RUN_MPI_COUNT_TRIAL) :: trial_values   !  (MAXPOP)
END TYPE run_mpi_type
!
TYPE ( run_mpi_type)     :: run_mpi_senddata
INTEGER                  :: sdl_length
!
CHARACTER(LEN=200), DIMENSION(:), ALLOCATABLE :: prog_entry
INTEGER           , DIMENSION(:), ALLOCATABLE :: socket_id
INTEGER           , DIMENSION(:,:), ALLOCATABLE :: port_id

!
INTEGER                  :: progs_size_of = 0
!
END MODULE run_mpi_mod
