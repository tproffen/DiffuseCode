MODULE gen_mpi_mod
!
! Generic MPI related variables
!
IMPLICIT NONE
!
INTEGER                  :: gen_mpi_myid        ! Id of master or slave
INTEGER                  :: gen_mpi_numprocs    ! Number of available processors 
!
LOGICAL                  :: gen_mpi_active      = .false. ! With/without MPI?
!
CHARACTER(LEN=1024), DIMENSION(:)  , ALLOCATABLE :: node_names     ! Identifiers for each node
INTEGER            , DIMENSION(:)  , ALLOCATABLE :: slave_is_node  ! Slave is on this node number
INTEGER            , DIMENSION(:)  , ALLOCATABLE :: kid_at_indiv   ! Kid is currently at this child
INTEGER            , DIMENSION(:)  , ALLOCATABLE :: kid_at_node    ! Kid is placed onto this node
INTEGER            , DIMENSION(:,:), ALLOCATABLE :: node_has_kids  ! Which kids are at this node
INTEGER            , DIMENSION(:)  , ALLOCATABLE :: node_max_kids  ! Maximum kids to be placed onto this node
LOGICAL            , DIMENSION(:)  , ALLOCATABLE :: node_finished  ! This node has done all its jobs
!
END MODULE gen_mpi_mod
