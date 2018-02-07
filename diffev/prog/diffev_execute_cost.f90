SUBROUTINE diffev_execute_cost( repeat,    &
                         prog_len,         &
                         prog  ,  prog_l , &
                         mac_len,          &
                         mac   ,  mac_l  , &
                         direc_len,        &
                         direc ,  direc_l, &
                         kid   ,  indiv  , &
                         n_rvalue_i, n_rvalue_o, NRVAL,  &
                         rvalue, l_rvalue, &
                         output_len,       &
                         output, output_l, &
                         generation, member, &
                         children, parameters, &
                                 nindiv  , &
                         trial_names ,         &
                         trial_values, NTRIAL, &
                         l_get_random_state,     &
                         rd_nseeds,rd_seeds,     &
                         ierr )
!
! specific funtion to execute the cost function from diffev
!
!  This version runs a "system" call, just as in run_mpi_slave,
!  This version is called if diffev runs as standalone and 
!  MPI is not active.
!
IMPLICIT NONE
LOGICAL                , INTENT(IN) :: repeat
INTEGER                , INTENT(IN) :: prog_len
INTEGER                , INTENT(IN) :: prog_l
INTEGER                , INTENT(IN) :: mac_len
INTEGER                , INTENT(IN) :: mac_l
INTEGER                , INTENT(IN) :: direc_len
INTEGER                , INTENT(IN) :: direc_l
INTEGER                , INTENT(IN) :: output_len
INTEGER                , INTENT(IN) :: output_l
CHARACTER(LEN=prog_len  ), INTENT(IN) :: prog  
CHARACTER(LEN=mac_len   ), INTENT(IN) :: mac   
CHARACTER(LEN=direc_len ), INTENT(IN) :: direc 
INTEGER                , INTENT(IN) :: kid   
INTEGER                , INTENT(IN) :: indiv
INTEGER                , INTENT(IN ):: n_rvalue_i
INTEGER                , INTENT(OUT):: n_rvalue_o
INTEGER                , INTENT(IN ):: NRVAL
REAL, DIMENSION(0:NRVAL   ),INTENT(OUT):: rvalue
LOGICAL                , INTENT(OUT):: l_rvalue
CHARACTER(LEN=output_len), INTENT(IN) :: output
INTEGER                , INTENT(IN) :: generation
INTEGER                , INTENT(IN) :: member
INTEGER                , INTENT(IN) :: children
INTEGER                , INTENT(IN) :: parameters
INTEGER                , INTENT(IN) :: nindiv
INTEGER                , INTENT(OUT):: ierr 
INTEGER                , INTENT(IN) :: NTRIAL
LOGICAL                , INTENT(IN)  :: l_get_random_state
INTEGER                , INTENT(OUT) :: rd_nseeds
INTEGER, DIMENSION(64) , INTENT(OUT) :: rd_seeds
CHARACTER(LEN=16),DIMENSION(1:NTRIAL),INTENT(IN) :: trial_names
REAL,DIMENSION(1:NTRIAL),INTENT(IN) :: trial_values
!
CHARACTER(LEN=2048) :: line
INTEGER             :: job_l
!
INTEGER             :: len_str
INTEGER             :: system
!
! If instructed, get state of random number generator
! As there is no connection to slave program set all to zero
!
IF(l_get_random_state) THEN
   rd_nseeds   = 1
   rd_seeds(:) = 0
ENDIF
!
IF ( repeat ) THEN       ! NINDIV calculations needed
   WRITE(line,2000) prog (1:prog_l ), &
                    mac  (1:mac_l  ), &
                    direc(1:direc_l), &
                    kid   , indiv   , &
                    output(1:output_l)
ELSE
   WRITE(line,2100) prog (1:prog_l ), &
                    mac  (1:mac_l  ), &
                    direc(1:direc_l), &
                    kid                              , &
                    output(1:output_l)
ENDIF
job_l = len_str(line)
ierr  = system( line(1:job_l))
!
2000 FORMAT ( a,' -macro ',a,' ',a,i8,2x,i8,' > ',a)
2100 FORMAT ( a,' -macro ',a,' ',a,i8,      ' > ',a)
!
END SUBROUTINE diffev_execute_cost
