SUBROUTINE refine_loop_mpi(prog_n, prog_l, mac_n, mac_l, out_n, out_l, repeat, nindiv)
!
!  This routine is called if a "run_mpi" command occurs within a do loop or
!  if block and we are running discus_suite without MPI active.
!
!  Copies the run_mpi parameters detected from a "run_mpi" within a loop
!  into the senddata structure.
!
USE run_mpi_mod
!
IMPLICIT NONE
!
CHARACTER (LEN=*), INTENT(IN) :: prog_n
CHARACTER (LEN=*), INTENT(IN) :: mac_n
CHARACTER (LEN=*), INTENT(IN) :: out_n
INTEGER          , INTENT(IN) :: prog_l
INTEGER          , INTENT(IN) :: mac_l
INTEGER          , INTENT(IN) :: out_l
LOGICAL          , INTENT(IN) :: repeat
INTEGER          , INTENT(IN) :: nindiv
!
!run_mpi_senddata%prog   = prog_n
!run_mpi_senddata%prog_l = prog_l
!run_mpi_senddata%mac    = mac_n
!run_mpi_senddata%mac_l  = mac_l
!run_mpi_senddata%out    = out_n
!run_mpi_senddata%out_l  = out_l
!run_mpi_senddata%repeat = repeat
!run_mpi_senddata%nindiv = nindiv
!
!run_mpi_senddata%generation = pop_gen    ! Current GENERATION no
!run_mpi_senddata%member     = pop_n      ! Number of members
!run_mpi_senddata%children   = pop_c      ! Number of children
!run_mpi_senddata%parameters = pop_dimx   ! Number of parameters
!run_mpi_senddata%use_socket = .false.
!
!run_mpi_senddata%l_get_state = l_get_random_state
!
END SUBROUTINE refine_loop_mpi
