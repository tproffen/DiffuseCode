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
INTEGER, DIMENSION(0:2)  :: run_mpi_oldtypes
INTEGER, DIMENSION(0:2)  :: run_mpi_blockcounts
INTEGER, DIMENSION(0:2)  :: run_mpi_offsets
INTEGER                  :: run_mpi_extent
INTEGER                  :: run_mpi_data_type
INTEGER                  :: run_mpi_data_type2
!
TYPE run_mpi_type                          ! MPI with types does not work yet
   SEQUENCE
   LOGICAL               :: repeat   !  1
   INTEGER               :: nindiv   !  2
   INTEGER               :: member   !  3
   INTEGER               :: kid      !  4
   INTEGER               :: ierr     !  5
   INTEGER               :: direc_l  !  6
   INTEGER               :: prog_l   !  7
   INTEGER               :: mac_l    !  8
   INTEGER               :: out_l    !  9
   CHARACTER (LEN=240)   :: direc    ! 10 : 249
   CHARACTER (LEN=100)   :: prog     !250 : 349
   CHARACTER (LEN=100)   :: mac      !350 : 449
   CHARACTER (LEN=100)   :: out      !450 : 549
END TYPE run_mpi_type
!
TYPE ( run_mpi_type)   :: run_mpi_senddata
INTEGER, DIMENSION(1:549) :: run_mpi_send_data
!
END MODULE run_mpi_mod
