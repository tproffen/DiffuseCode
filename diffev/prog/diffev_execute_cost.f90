SUBROUTINE diffev_execute_cost( repeat,    &
                         prog_len,         &
                         prog  ,  prog_l , &
                         mac_len,          &
                         mac   ,  mac_l  , &
                         direc_len,        &
                         direc ,  direc_l, &
                         kid   ,  indiv  , &
                         rvalue, l_rvalue, &
                         output_len,       &
                         output, output_l, &
                         generation, member, &
                         children, parameters, &
                                 nindiv  , &
                         trial_values, NTRIAL, &
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
REAL                   , INTENT(OUT):: rvalue
LOGICAL                , INTENT(OUT):: l_rvalue
CHARACTER(LEN=output_len), INTENT(IN) :: output
INTEGER                , INTENT(IN) :: generation
INTEGER                , INTENT(IN) :: member
INTEGER                , INTENT(IN) :: children
INTEGER                , INTENT(IN) :: parameters
INTEGER                , INTENT(IN) :: nindiv
INTEGER                , INTENT(OUT):: ierr 
INTEGER                , INTENT(IN) :: NTRIAL
REAL,DIMENSION(1:NTRIAL),INTENT(IN) :: trial_values
!
CHARACTER(LEN=2048) :: line
INTEGER             :: job_l
!
INTEGER             :: len_str
INTEGER             :: system
!
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
