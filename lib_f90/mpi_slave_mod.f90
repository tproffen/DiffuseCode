MODULE mpi_slave_mod
!
LOGICAL  :: mpi_active   = .false.
LOGICAL  :: mpi_is_slave = .false.
INTEGER  :: mpi_slave_error = 0
INTEGER  :: mpi_slave_err_typ = 0
character(len=80), dimension(7) :: mpi_slave_msg = ' '
!
END MODULE mpi_slave_mod
