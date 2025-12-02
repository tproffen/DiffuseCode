MODULE run_mpi_mod
!
use precision_mod
use lib_global_flags_mod
!
IMPLICIT none
!
INTEGER, PARAMETER       :: RUN_MPI_COUNT_INTEGER   =  16
INTEGER, PARAMETER       :: RUN_MPI_COUNT_LOGICAL   =   4
INTEGER, PARAMETER       :: RUN_MPI_COUNT_CHARACTER = 540
INTEGER, PARAMETER       :: RUN_MPI_COUNT_REAL      =   2
INTEGER, PARAMETER       :: RUN_MPI_COUNT_TRIAL     = 1000
INTEGER, PARAMETER       :: RUN_MPI_MAXRVALUE       =  15
INTEGER, PARAMETER       :: RUN_MPI_MAX_FLAGS       =  LIB_GLOBAL_FLAGS_MAX
!INTEGER                  :: RUN_MPI_COUNT_TRIAL     = 200
!INTEGER                  :: RUN_MPI_MAXRVALUE       =  15
INTEGER, PARAMETER       :: RUN_MPI_NSEEDS          =  64 !Number of seeds for random
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
TYPE run_mpi_type                          ! MPI with types does not work yet
   INTEGER               :: generation !  1  !  1
   INTEGER               :: member     !  2  !  2
   INTEGER               :: children   !  3  !  3
   INTEGER               :: parameters !  4  !  4
   INTEGER               :: nindiv     !  5  !  5
   INTEGER               :: kid        !  6  !  6
   INTEGER               :: indiv      !  7  !  7
   INTEGER               :: ierr       !  8  !  8
   INTEGER               :: ierr_typ   !  8  !  8   ! SEQUENCE ?!?
   INTEGER               :: ierr_msg_l     !  8  !  8   ! SEQUENCE ?!?
   INTEGER               :: ierr_msg_n     !  8  !  8   ! SEQUENCE ?!?
   INTEGER               :: direc_l    !  9  !  9
   INTEGER               :: prog_l     ! 10  ! 10
   INTEGER               :: mac_l      ! 11  ! 11
   INTEGER               :: out_l      ! 12  ! 12
   INTEGER               :: prog_num   ! 13  ! 13
   INTEGER               :: s_remote   ! 14  ! 14
   INTEGER               :: port       ! 15  ! 15
   INTEGER               :: nseeds     ! 16  ! 16
   INTEGER,DIMENSION(1:RUN_MPI_NSEEDS) :: seeds   !  (nseeds)
   INTEGER               :: n_rvalue_i ! 16  ! 16
   INTEGER               :: n_rvalue_o ! 16  ! 16
   INTEGER               :: RUN_MPI_COUNT_TRIAL
   INTEGER               :: RUN_MPI_MAXRVALUE
   INTEGER               :: RUN_MPI_MAX_FLAGS
   integer(KIND=PREC_INT_BYTE), dimension(RUN_MPI_MAX_FLAGS) :: global_flags
   LOGICAL               :: repeat     ! 21  ! 21
!   LOGICAL               :: use_socket ! 22  ! 22
   LOGICAL               :: prog_start ! 23  ! 23
   LOGICAL               :: l_rvalue   ! 24  ! 24
   integer               :: l_get_state! 25  ! 25
   LOGICAL               :: l_first_job! 26  ! 26
   LOGICAL               :: spacer2    ! 27  ! 27
   LOGICAL               :: spacer3    ! 28  ! 28
   character(len=80), dimension(7) :: ierr_msg
   CHARACTER (LEN=240)   :: direc      ! 29 : 268
   CHARACTER (LEN=100)   :: prog       !269 : 368
   CHARACTER (LEN=100)   :: mac        !369 : 468
   CHARACTER (LEN=100)   :: out        !469 : 568
   real(kind=PREC_DP), dimension(3,3)                  :: data_ext       ! F_xmin, etc. for the 3 dimensions
   REAL(kind=PREC_DP),DIMENSION(0:RUN_MPI_MAXRVALUE  ) :: rvalue         !  (MAXPOP)
   REAL(kind=PREC_DP),DIMENSION(1:RUN_MPI_COUNT_TRIAL) :: trial_values   !  (MAXPOP)
   CHARACTER(LEN=16),DIMENSION(1:RUN_MPI_COUNT_TRIAL) :: trial_names   !  (MAXPOP)
!  REAL(kind=PREC_DP),DIMENSION(:), allocatable        :: rvalue         !  (MAXPOP)
!  REAL(kind=PREC_DP),DIMENSION(:), allocatable        :: trial_values   !  (MAXPOP)
!  CHARACTER(LEN=16),DIMENSION(:), allocatable        :: trial_names   !  (MAXPOP)
END TYPE run_mpi_type
!
TYPE ( run_mpi_type)     :: run_mpi_senddata
INTEGER                  :: sdl_length
INTEGER                  :: NUM_NODE = 0
INTEGER                  :: run_mpi_max_slaves = 0                !Maximum slaves on any node
INTEGER                  :: run_mpi_kid_per_core
!
!
CHARACTER(LEN=200), DIMENSION(:), ALLOCATABLE :: prog_entry
!INTEGER           , DIMENSION(:), ALLOCATABLE :: socket_id
INTEGER           , DIMENSION(:,:), ALLOCATABLE :: port_id
!
!
contains
!
!*******************************************************************************
!
subroutine run_mpi_senddata_init
!-
! Initialize all entrie in run_mpi_senddata
!
run_mpi_senddata%generation = 0
run_mpi_senddata%member     = 0
run_mpi_senddata%children   = 0
run_mpi_senddata%parameters = 0
run_mpi_senddata%nindiv     = 0
run_mpi_senddata%kid        = 0
run_mpi_senddata%indiv      = 0
run_mpi_senddata%ierr       = 0
run_mpi_senddata%ierr_typ   = 0
run_mpi_senddata%ierr_msg_l = 0
run_mpi_senddata%ierr_msg_n = 0
run_mpi_senddata%direc_l    = 1
run_mpi_senddata%prog_l     = 6
run_mpi_senddata%mac_l      = 8
run_mpi_senddata%out_l      = 9
run_mpi_senddata%prog_num   = 0
run_mpi_senddata%s_remote   = 0
run_mpi_senddata%port       = 0
run_mpi_senddata%nseeds     = RUN_MPI_NSEEDS
run_mpi_senddata%seeds      = 0
run_mpi_senddata%n_rvalue_i = 0
run_mpi_senddata%n_rvalue_o = 0
run_mpi_senddata%RUN_MPI_COUNT_TRIAL = RUN_MPI_COUNT_TRIAL
run_mpi_senddata%RUN_MPI_MAXRVALUE   = RUN_MPI_MAXRVALUE
run_mpi_senddata%RUN_MPI_MAX_FLAGS   = RUN_MPI_MAX_FLAGS
run_mpi_senddata%global_flags = lib_global_flags
run_mpi_senddata%repeat       = .false.
run_mpi_senddata%prog_start   = .true.
run_mpi_senddata%l_rvalue     = .false.
run_mpi_senddata%l_get_state  = 0
run_mpi_senddata%l_first_job  = .true.
run_mpi_senddata%spacer2      = .true.
run_mpi_senddata%spacer3      = .true.
run_mpi_senddata%ierr_msg     = ' '
run_mpi_senddata%direc        = '.'
run_mpi_senddata%prog         = 'discus'
run_mpi_senddata%mac          = 'none.mac'
run_mpi_senddata%out          = '/dev/null'
run_mpi_senddata%trial_values = 0.0_PREC_DP
run_mpi_senddata%rvalue       = 0.0_PREC_DP
run_mpi_senddata%trial_names  = ' '
!
!
end subroutine run_mpi_senddata_init
!
END MODULE run_mpi_mod
