SUBROUTINE suite_execute_cost( repeat,           &
                         prog_len,         &
                         prog  ,  prog_l , &
                         mac_len,          &
                         mac   ,  mac_l  , &
                         direc_len,        &
                         direc ,  direc_l, &
                         kid   ,  indiv  , &
                         n_rvalue_i, n_rvalue_o, NRVAL,        &
                         rvalue, l_rvalue, &
                         output_len,       &
                         output, output_l, &
                         generation, member, &
                         children, parameters, &
                                 nindiv  , &
                         trial_n,         &
                         trial_v, NTRIAL, &
                         l_get_random_state,     &
                         rd_nseeds,rd_seeds,     &
                         l_first_job,            &
                         ierr )
!
USE class_internal

USE diffev_setup_mod
USE discus_setup_mod
USE discus_loop_mod
USE structur, ONLY: rese_cr
USE kuplot_setup_mod
USE kuplot_loop_mod
USE suite_setup_mod
USE suite_set_sub_mod
!
USE appl_env_mod
USe define_variable_mod
USE do_set_mod
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
INTEGER                , INTENT(IN)  :: n_rvalue_i
INTEGER                , INTENT(OUT) :: n_rvalue_o
INTEGER                , INTENT(IN)  :: NRVAL
REAL, DIMENSION(0:NRVAL), INTENT(OUT) :: rvalue
!REAL                   , INTENT(OUT):: rvalue
LOGICAL                , INTENT(OUT):: l_rvalue
INTEGER                , INTENT(IN) :: generation
INTEGER                , INTENT(IN) :: member
INTEGER                , INTENT(IN) :: children
INTEGER                , INTENT(IN) :: parameters
INTEGER                , INTENT(IN) :: nindiv
INTEGER                , INTENT(IN) :: NTRIAL
CHARACTER(LEN=16),DIMENSION(1:NTRIAL),INTENT(IN) :: trial_n
REAL,DIMENSION(1:NTRIAL),INTENT(IN) :: trial_v
LOGICAL                , INTENT(IN)  :: l_get_random_state
INTEGER                , INTENT(OUT) :: rd_nseeds
INTEGER, DIMENSION(64) , INTENT(OUT) :: rd_seeds
LOGICAL                , INTENT(IN ) :: l_first_job
INTEGER                , INTENT(OUT):: ierr
!
LOGICAL, PARAMETER :: IS_DIFFEV = .TRUE.
CHARACTER(LEN=2048) :: line
CHARACTER(LEN=1024) :: empty = ' '
CHARACTER(LEN=2048) :: logfile
CHARACTER(LEN= 7   ) :: ct_pname_old,ct_pname_cap_old
CHARACTER(LEN= 9+LEN(trial_n))   :: string
INTEGER              :: ct_prompt_status_old

LOGICAL, SAVE       :: l_discus_init =.false.
LOGICAL, SAVE       :: l_kuplot_init =.false.
INTEGER             :: i,j
INTEGER             :: length = 1
INTEGER             :: job_l
INTEGER             :: ios
INTEGER             :: ct_ostate    ! Original program state
INTEGER             :: ct_oprogr    ! Original program
INTEGER             :: ct_ompifirst ! Original mpi_first state
LOGICAL :: str_comp
INTEGER :: len_str
!
ier_mpi = .TRUE.
!
CALL do_chdir(direc,direc_l,.FALSE.)    ! Set current directory as passed from master
!
! Store old program name and prompt status
!
ct_pname_old         = pname
ct_pname_cap_old     = pname_cap
ct_prompt_status_old = prompt_status
ct_ostate            = NINT(var_val(VAR_STATE))
ct_oprogr            = NINT(var_val(VAR_PROGRAM))
ct_ompifirst         = NINT(var_val(VAR_MPI_FIRST))

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
loop_par: DO i=1,parameters            ! Trial parameters for this kid
   ref_para(i) = trial_v(i)
   string = ' '
   string = 'real, '//trial_n(i)
   length = LEN(string)
   CALL define_variable(string,length, IS_DIFFEV)
   DO j=1,var_num
      IF(var_name(j)==trial_n(i)) THEN
         var_val(j) = trial_v(i)
         CYCLE loop_par
      ENDIF
   ENDDO
ENDDO loop_par
!
var_val( var_ref+0) = generation
var_val( var_ref+1) = member
var_val( var_ref+2) = children
var_val( var_ref+3) = parameters
var_val( var_ref+4) = kid
var_val( var_ref+5) = indiv
var_val( var_ref+6) = nindiv
IF(repeat) THEN
   var_val( var_ref+7) = var_val(VAR_TRUE)
ELSE
   var_val( var_ref+7) = var_val(VAR_FALSE)
ENDIF
IF(l_first_job) THEN
   var_val(VAR_MPI_FIRST) = var_val(VAR_TRUE)
ELSE
   var_val(VAR_MPI_FIRST) = var_val(VAR_FALSE)
ENDIF
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
   rvalue(:)  = HUGE(0.0E0)
!
!  Reset KUPLOT
!
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
   CALL do_rese (empty, length)
!
!  Reset DISCUS
!
   IF(.NOT. l_discus_init) THEN
      CALL discus_setup   (lstandalone)
      l_discus_init = .true.
   ENDIF
!
! If instructed, get state of random number generator
!
   IF(l_get_random_state) THEN
      CALL random_current(rd_nseeds, rd_seeds)
   ENDIF
!
   IF(str_comp(prog, 'discus', 6, prog_l, 6)) THEN
      CALL store_remove_all(store_root)    ! Allways do a DISCUS internal storage reset
      pname     = 'discus'
      pname_cap = 'DISCUS'
      prompt    = pname
      oprompt   = pname
      var_val(VAR_PROGRAM)= var_val(VAR_DISCUS)   ! Set program to DISCUS
      var_val(VAR_STATE)  = var_val(VAR_IS_SECTION)
      CALL discus_set_sub ()
      CALL suite_set_sub_branch
      CALL rese_cr
      CALL discus_loop ()
   ELSEIF(str_comp(prog, 'kuplot', 6, prog_l, 6)) THEN
      pname     = 'kuplot'
      pname_cap = 'KUPLOT'
      prompt    = pname
      oprompt   = pname
      var_val(VAR_PROGRAM)= var_val(VAR_KUPLOT)   ! Set program to KUPLOT
      var_val(VAR_STATE)  = var_val(VAR_IS_SECTION)
      CALL kuplot_set_sub ()
      CALL suite_set_sub_branch
      CALL kuplot_loop ()
   ENDIF
!
   mpi_slave_error = ier_num
   CALL macro_close_mpi(mac, mac_l)
   IF(ier_num /= 0) mpi_slave_error = ier_num
!
   IF(rvalue_yes) THEN    ! We got an r-value
      l_rvalue = .true.
      n_rvalue_o = 1
      rvalue(0)= rvalues(2,0)
      IF(ISNAN(rvalues(2,0))) THEN
         ier_num = -49
         ier_typ = ER_FORT
      ENDIF
      n_rvalue_o = 0
      IF(nrvalues == n_rvalue_i) THEN
         rvalue(1:n_rvalue_i) = rvalues(2,1:n_rvalue_i)
         n_rvalue_o = nrvalues
      ENDIF
      DO i =0, n_rvalue_i
         IF(ISNAN(rvalue(i))) THEN
            ier_num = -49
            ier_typ = ER_FORT
         ENDIF
      ENDDO
      IF(ier_num /= 0) mpi_slave_error = ier_num
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
var_val(VAR_STATE)     = ct_ostate
var_val(VAR_PROGRAM)   = ct_oprogr
var_val(VAR_MPI_FIRST) = ct_ompifirst
!
!
IF(pname=='discus') THEN      ! Return to DISCUS branch
      CALL discus_set_sub ()
      var_val(VAR_PROGRAM)= var_val(VAR_DISCUS)   ! Set program to DISCUS
ELSEIF(pname=='diffev') THEN  ! Return to DIFFEV branch
      CALL diffev_set_sub ()
      var_val(VAR_PROGRAM)= var_val(VAR_DIFFEV)   ! Set program to DIFFEV
ELSEIF(pname=='kuplot') THEN  ! Return to KUPLOT branch
      CALL kuplot_set_sub ()
      var_val(VAR_PROGRAM)= var_val(VAR_KUPLOT)   ! Set program to KUPLOT
ENDIF
CALL suite_set_sub_branch ()
CALL program_files ()
ier_mpi = .FALSE.
!
2000 FORMAT ( a,' ',a,',',i8,',',i8)
2100 FORMAT ( a,' ',a,',',i8       )
!
END SUBROUTINE suite_execute_cost
