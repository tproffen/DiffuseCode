SUBROUTINE suite_execute_cost( repeat,           &
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
                         trial_v, NTRIAL, &
                         l_get_random_state,     &
                         rd_idum,rd_iff,         &
                         rd_ix1, rd_ix2, rd_ix3, &
                         ierr )
!
USE diffev_setup_mod
USE discus_setup_mod
USE discus_loop_mod
USE kuplot_setup_mod
USE kuplot_loop_mod
USE suite_setup_mod
!
USE errlist_mod
USE mpi_slave_mod
USE prompt_mod
USE param_mod
USE random_state_mod
USE variable_mod
USE lib_f90_allocate_mod
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
CHARACTER(LEN=prog_len), INTENT(IN) :: prog  
CHARACTER(LEN=mac_len ), INTENT(IN) :: mac   
CHARACTER(LEN=direc_len), INTENT(IN) :: direc 
INTEGER                , INTENT(IN) :: kid   
INTEGER                , INTENT(IN) :: indiv
CHARACTER(LEN=output_len), INTENT(IN) :: output
REAL                   , INTENT(OUT):: rvalue
LOGICAL                , INTENT(OUT):: l_rvalue
INTEGER                , INTENT(IN) :: generation
INTEGER                , INTENT(IN) :: member
INTEGER                , INTENT(IN) :: children
INTEGER                , INTENT(IN) :: parameters
INTEGER                , INTENT(IN) :: nindiv
INTEGER                , INTENT(IN) :: NTRIAL
REAL,DIMENSION(1:NTRIAL),INTENT(IN) :: trial_v
LOGICAL                , INTENT(IN)  :: l_get_random_state
INTEGER                , INTENT(OUT) :: rd_idum
INTEGER                , INTENT(OUT) :: rd_iff
INTEGER                , INTENT(OUT) :: rd_ix1
INTEGER                , INTENT(OUT) :: rd_ix2
INTEGER                , INTENT(OUT) :: rd_ix3
INTEGER                , INTENT(OUT):: ierr
!
CHARACTER(LEN=2048) :: line
CHARACTER(LEN=2048) :: logfile
CHARACTER(LEN= 7   ) :: ct_pname_old,ct_pname_cap_old
INTEGER              :: ct_prompt_status_old

LOGICAL, SAVE       :: l_discus_init =.false.
LOGICAL, SAVE       :: l_kuplot_init =.false.
INTEGER             :: i
INTEGER             :: job_l
INTEGER             :: ios
LOGICAL :: str_comp
INTEGER :: len_str
!
CALL do_chdir(direc,direc_l,.FALSE.)    ! Set current directeory as passed from master
!
! If instructed, get state of random number generator
!
IF(l_get_random_state) THEN
   CALL random_current(rd_idum, rd_iff, rd_ix1, rd_ix2, rd_ix3)
ENDIF
!
! Store old program name and prompt status
!
ct_pname_old         = pname
ct_pname_cap_old     = pname_cap
ct_prompt_status_old = prompt_status

IF(str_comp(output(1:output_l), '/dev/null', 9, output_l, 9)) THEN
   line  = 'prompt, redirect, off'
   job_l =  21
   CALL do_set(line, job_l)
ELSEIF(str_comp(output(1:output_l),'on',2,output_l,2)) THEN
!   prompt_status = PROMPT_REDIRECT
!   socket_status = PROMPT_REDIRECT
   IF(output_IO==37) THEN
      CLOSE(OUTPUT_IO)
   ENDIF
   output_io = 6
ELSE
!   prompt_status = PROMPT_REDIRECT
!   socket_status = PROMPT_REDIRECT
   output_io = 37
   logfile   = output(1:output_l)
   OPEN(UNIT=output_IO, FILE=logfile, STATUS='unknown', IOSTAT=ios)
   IF(ios/=0) THEN
      ier_num = -2
      ier_typ = ER_IO
      ier_msg(1) = logfile(1:80)
      ier_msg(2) = 'Check spelling, existence of directories'
      RETURN
   ENDIF
ENDIF
!
! Interim solution up to 5.6.4 will be phased out
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
! Long term solution copy into specialized variables
!
IF(UBOUND(ref_para, 1) < parameters) THEN
   CALL alloc_ref_para(parameters)
   MAXPAR_REF = parameters
ENDIF
DO i=1,parameters            ! Trial parameters for this kid
   ref_para(i) = trial_v(i)
ENDDO
!
var_val( var_ref+0) = generation
var_val( var_ref+1) = member
var_val( var_ref+2) = children
var_val( var_ref+3) = parameters
var_val( var_ref+4) = kid
var_val( var_ref+5) = indiv
var_val( var_ref+6) = nindiv
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
IF(ier_num == 0 ) THEN  ! Defined macro with no error
!
   rvalue_yes = .false.   ! Reset the global R-value flag
   l_rvalue   = .false.
   rvalue     = 0.0
   IF(str_comp(prog, 'discus', 6, prog_l, 6)) THEN
      IF(.NOT. l_discus_init) THEN
         CALL discus_setup   (lstandalone)
         l_discus_init = .true.
      ENDIF
      pname     = 'discus'
      pname_cap = 'DISCUS'
      prompt    = pname
      oprompt   = pname
      CALL discus_set_sub ()
      CALL suite_set_sub_branch
      CALL discus_loop ()
   ELSEIF(str_comp(prog, 'kuplot', 6, prog_l, 6)) THEN
      IF(.NOT. l_kuplot_init) THEN
         CALL kuplot_setup   (lstandalone)
         l_kuplot_init = .true.
      ENDIF
      pname     = 'kuplot'
      pname_cap = 'KUPLOT'
      prompt    = pname
      oprompt   = pname
      CALL kuplot_set_sub ()
      CALL suite_set_sub_branch
      CALL kuplot_loop ()
   ENDIF
!
   CALL macro_close_mpi(mac, mac_l)
!
   IF(rvalue_yes) THEN    ! We got an r-value
      l_rvalue = .true.
      rvalue   = rvalues(2)
   ENDIF
ELSE
   mpi_slave_error = ier_num
ENDIF
!
IF(output_IO==37) THEN
   CLOSE(OUTPUT_IO)
ENDIF
!
ierr = mpi_slave_error
!
! Restore old program secion settings, i.e. always DIFFEV...
!
pname         = ct_pname_old
pname_cap     = ct_pname_cap_old
prompt_status = ct_prompt_status_old
prompt        = pname
oprompt       = pname
!
!
IF(pname=='discus') THEN      ! Return to DISCUS branch
      CALL discus_set_sub ()
ELSEIF(pname=='diffev') THEN  ! Return to DIFFEV branch
      CALL diffev_set_sub ()
ELSEIF(pname=='kuplot') THEN  ! Return to KUPLOT branch
      CALL kuplot_set_sub ()
ENDIF
CALL suite_set_sub_branch ()
CALL program_files ()
!
2000 FORMAT ( a,' ',a,',',i8,',',i8)
2100 FORMAT ( a,' ',a,',',i8       )
!
END SUBROUTINE suite_execute_cost
