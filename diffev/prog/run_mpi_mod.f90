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
INTEGER, PARAMETER       :: RUN_MPI_NSEEDS          =  12 !Number of seeds for random
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
   INTEGER               :: generation !  1  !  1
   INTEGER               :: member     !  2  !  2
   INTEGER               :: children   !  3  !  3
   INTEGER               :: parameters !  4  !  4
   INTEGER               :: nindiv     !  5  !  5
   INTEGER               :: kid        !  6  !  6
   INTEGER               :: indiv      !  7  !  7
   INTEGER               :: ierr       !  8  !  8
   INTEGER               :: direc_l    !  9  !  9
   INTEGER               :: prog_l     ! 10  ! 10
   INTEGER               :: mac_l      ! 11  ! 11
   INTEGER               :: out_l      ! 12  ! 12
   INTEGER               :: prog_num   ! 13  ! 13
   INTEGER               :: s_remote   ! 14  ! 14
   INTEGER               :: port       ! 15  ! 15
   INTEGER               :: nseeds     ! 16  ! 16
   INTEGER,DIMENSION(1:RUN_MPI_NSEEDS) :: seeds   !  (nseeds)
   LOGICAL               :: repeat     ! 21  ! 21
   LOGICAL               :: use_socket ! 22  ! 22
   LOGICAL               :: prog_start ! 23  ! 23
   LOGICAL               :: l_rvalue   ! 24  ! 24
   LOGICAL               :: l_get_state! 25  ! 25
   LOGICAL               :: spacer1    ! 26  ! 26
   LOGICAL               :: spacer2    ! 27  ! 27
   LOGICAL               :: spacer3    ! 28  ! 28
   CHARACTER (LEN=240)   :: direc      ! 29 : 268
   CHARACTER (LEN=100)   :: prog       !269 : 368
   CHARACTER (LEN=100)   :: mac        !369 : 468
   CHARACTER (LEN=100)   :: out        !469 : 568
   REAL                  :: rvalue     !569 : 569
   REAL                  :: rvalue2    !570 : 570
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
