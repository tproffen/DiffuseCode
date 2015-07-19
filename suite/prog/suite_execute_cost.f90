SUBROUTINE suite_execute_cost( repeat,           &
                         prog  ,  prog_l , &
                         mac   ,  mac_l  , &
                         direc ,  direc_l, &
                         kid   ,  indiv  , &
                         rvalue, l_rvalue, &
                         output, output_l, &
                         generation, member, &
                         children, parameters, &
                         trial_v, NTRIAL, &
                         ierr )
!
USE discus_setup_mod
USE discus_loop_mod
USE kuplot_setup_mod
USE kuplot_loop_mod
USE suite_setup_mod
!
USE errlist_mod
USE prompt_mod
USE param_mod
!
IMPLICIT NONE
LOGICAL                , INTENT(IN) :: repeat
INTEGER                , INTENT(IN) :: prog_l
INTEGER                , INTENT(IN) :: mac_l
INTEGER                , INTENT(IN) :: direc_l
INTEGER                , INTENT(IN) :: output_l
CHARACTER(LEN=prog_l  ), INTENT(IN) :: prog  
CHARACTER(LEN=mac_l   ), INTENT(IN) :: mac   
CHARACTER(LEN=direc_l ), INTENT(IN) :: direc 
INTEGER                , INTENT(IN) :: kid   
INTEGER                , INTENT(IN) :: indiv
CHARACTER(LEN=output_l), INTENT(IN) :: output
REAL                   , INTENT(OUT):: rvalue
LOGICAL                , INTENT(OUT):: l_rvalue
INTEGER                , INTENT(IN) :: generation
INTEGER                , INTENT(IN) :: member
INTEGER                , INTENT(IN) :: children
INTEGER                , INTENT(IN) :: parameters
INTEGER                , INTENT(IN) :: NTRIAL
REAL,DIMENSION(1:NTRIAL),INTENT(IN) :: trial_v
INTEGER                , INTENT(OUT):: ierr
!
CHARACTER(LEN=2048) :: line
CHARACTER(LEN=2048) :: logfile
LOGICAL, SAVE       :: l_discus_init =.false.
LOGICAL, SAVE       :: l_kuplot_init =.false.
INTEGER             :: i
INTEGER             :: job_l
LOGICAL :: str_comp
INTEGER :: len_str
!
IF(str_comp(output(1:output_l), '/dev/null', 9, output_l, 9)) THEN
   line  = 'prompt, redirect, off'
   job_l =  21
   CALL do_set(line, job_l)
ELSEIF(str_comp(output(1:output_l),'on',2,output_l,2)) THEN
   prompt_status = PROMPT_REDIRECT
   socket_status = PROMPT_REDIRECT
   IF(output_IO==37) THEN
      CLOSE(OUTPUT_IO)
   ENDIF
   output_io = 6
ELSE
   prompt_status = PROMPT_REDIRECT
   socket_status = PROMPT_REDIRECT
   output_io = 37
   logfile   = output(1:output_l)
   OPEN(UNIT=output_IO, FILE=logfile, STATUS='unknown')
ENDIF
!
DO i=1,parameters
   rpara(200+i) = trial_v(i)
ENDDO
inpara(200) = 4
inpara(201) = generation
inpara(202) = member
inpara(203) = children
inpara(204) = parameters
!
!  build macro line
!
IF ( repeat ) THEN       ! NINDIV calculations needed
   WRITE(line,2000) mac  (1:mac_l  ), &
                    direc(1:direc_l), &
                    kid   , indiv
ELSE
   WRITE(line,2100) mac  (1:mac_l  ), &
                    direc(1:direc_l), &
                    kid
ENDIF
job_l = len_str(line)
CALL file_kdo(line, job_l)
!
rvalue_yes = .false.   ! Reset the global R-value flag
l_rvalue   = .false.
rvalue     = 0.0
IF(str_comp(prog, 'discus', 6, prog_l, 6)) THEN
   IF(.NOT. l_discus_init) THEN
      CALL discus_setup   (lstandalone)
      l_discus_init = .true.
   ENDIF
   CALL discus_set_sub ()
   CALL suite_set_sub_branch
   CALL discus_loop ()
ELSEIF(str_comp(prog, 'kuplot', 6, prog_l, 6)) THEN
   IF(.NOT. l_kuplot_init) THEN
      CALL kuplot_setup   (lstandalone)
      l_kuplot_init = .true.
   ENDIF
   CALL kuplot_set_sub ()
   CALL suite_set_sub_branch
   CALL kuplot_loop ()
ENDIF
!
IF(rvalue_yes) THEN    ! We got an r-value
   l_rvalue = .true.
   rvalue   = rvalues(2)
ENDIF
!
IF(output_IO==37) THEN
   CLOSE(OUTPUT_IO)
ENDIF
!
ierr = 0

2000 FORMAT ( a,' ',a,',',i8,',',i8)
2100 FORMAT ( a,' ',a,',',i8       )
!
END SUBROUTINE suite_execute_cost
