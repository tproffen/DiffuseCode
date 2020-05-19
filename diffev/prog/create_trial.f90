MODULE create_trial_mod
!
CONTAINS
!*****7**************************************************************** 
!                                                                       
      SUBROUTINE create_trial 
!                                                                       
!     creates a new generation.                                         
!     Each parameter is calculated as                                   
!     CHILD_P(1) = P(1) + f*(P(2)-P(3))      or:                        
!     CHILD_P(1) = P(1),                                                
!     where P(i) are randomly chosen members of the parent population.  
!     The crossover constant diff_cr determines whether a cross over    
!     takes place or whether the parameter is inherited from the        
!                                                                       
!     to avoid getting stuck in a local minimum or at the edge of       
!     an allowed region we do:                                          
!     - the minimum shift is gaussian distributed,                      
!            applied if (P(2)-P(3))=0                                   
!     - If the shift leads outside the allowed region,                  
!            the new parameter                                          
!       is gaussian distributed on the inside of the allowed region.    
!     The sigma is supplied by the user                                 
!                                                                       
!       For those parameters that are not affected by the cross-over, a 
!     local shift may be applied using the user supplied local sigma.   
!                                                                       
      USE population
      USE diff_evol
      USE random_mod
!                                                                       
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER j, jt 
!                                                                       
      REAL ran1 
!                                                                       
      CALL adapt_sigma 
!                                                                       
      IF (diff_sel_mode.eq.SEL_COMP) then 
!                                                                       
!     --Each child is compared to its parent, generate straight from    
!       parent                                                          
!                                                                       
         DO jt = 1, pop_n 
            j = jt 
            CALL create_single (j, jt) 
            CALL write_trial (jt) 
         ENDDO 
      ELSE 
!                                                                       
!     --Selection mode is "BEST".                                       
!       If there are equal or more children than parents, then          
!         create pop_n children straigth from each parent               
!         create additional children from random parents                
!                                                                       
         IF (pop_n.le.pop_c) then 
            DO jt = 1, pop_n 
               j = jt 
               CALL create_single (j, jt) 
               CALL write_trial (jt) 
            ENDDO 
            DO jt = pop_n + 1, pop_c 
               j = int (pop_n * ran1 (idum) ) + 1 
               CALL create_single (j, jt) 
               CALL write_trial (jt) 
            ENDDO 
         ELSE 
!     --There are less children than parents, create all children       
!       from randomly chosen parents                                    
!                                                                       
            DO jt = 1, pop_c 
               j = int (pop_n * ran1 (idum) ) + 1 
               CALL create_single (j, jt) 
               CALL write_trial (jt) 
            ENDDO 
         ENDIF 
      ENDIF 
!                                                                       
      END SUBROUTINE create_trial                   
!*****7**************************************************************** 
SUBROUTINE create_single (j, jt) 
!                                                                       
USE diffev_config
USE constraint
USE diff_evol
USE population
USE triple_perm
!USE do_if_mod
USE random_mod
USE errlist_mod 
USE do_execute_mod
USE precision_mod
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER,         INTENT(IN   ) :: j
INTEGER,         INTENT(IN   ) :: jt
!
CHARACTER (LEN=PREC_STRING)    :: line 
INTEGER                        :: i, ii, k, l
INTEGER                        :: n_tried 
INTEGER                        :: length 
INTEGER                        :: j1, j2, j3 
LOGICAL                        :: l_ok 
LOGICAL                        :: l_unchanged 
REAL                           :: shift 
REAL                           :: value 
REAL                           :: w 
!                                                                       
REAL                           :: ran1 
REAL(KIND=PREC_DP), EXTERNAL   :: gasdev 
!                                                                       
w = 0.0 
n_tried = 0 
l_ok = .false. 
main:DO while (.not.l_ok.and.n_tried.lt.MAX_CONSTR_TRIAL) 
   n_tried = n_tried+1 
local:  IF (ran1 (idum) .gt.diff_local) then 
      CALL do_triple_perm (j, j1, j2, j3, pop_n) 
      l_unchanged = .true. 
      ii = int (ran1 (idum) * pop_dimx) 
kdimx1:    DO k = 1, pop_dimx 
         i = mod (k + ii - 1, pop_dimx) + 1 
refine1: IF (pop_refine (i) ) then 
cross1:     IF (ran1 (idum) .lt.diff_cr.or.l_unchanged) then 
               l_unchanged = .false. 
!                                                                 
!-----------Determine crossover parameter shift                   
!                                                                 
               shift = diff_f * (pop_x (i, j1) - pop_x (i, j2) ) 
!                                                                 
!-----------Zero shift correction                                 
!                                                                 
 zero1:        IF (shift.eq.0) then 
                  IF (abs (pop_sigma (i) ) .gt.0.0) then 
                     shift = gasdev (DBLE(pop_sigma (i))) 
                  ELSE 
                     shift = 0.0 
                  ENDIF 
               ENDIF zero1
!                                                                 
!     ------Determine trial value                                 
!                                                                 
               value = pop_x (i, pop_best) ! sensible default value, will never be used...
 trial1:       IF (diff_donor_mode.eq.ADD_TO_RANDOM) then 
                  value = diff_k * pop_x (i, j3) + (1. - diff_k)        &
                  * pop_x (i, j) + shift                                
               ELSEIF (diff_donor_mode.eq.ADD_TO_BEST) then 
                  value = diff_k * pop_x (i, pop_best) + (1. - diff_k)  &
                  * pop_x (i, j) + shift                                
               ENDIF trial1
!                                                                 
!-----------Upper/Lower limit correction                          
!                                                                 
 limit1:       IF (value.gt.pop_xmax (i) ) then 
                  value = pop_xmax (i) - abs (gasdev (DBLE(pop_sigma (i) ))) 
               ELSEIF (value.lt.pop_xmin (i) ) then 
                  value = pop_xmin (i) + abs (gasdev (DBLE(pop_sigma (i) ))) 
               ENDIF limit1
!                                                                 
!-----------Apply final shift                                     
!                                                                 
 shift1:       IF (pop_type (i) .eq.POP_INTEGER) THEN 
                  pop_t (i, jt) = nint (max (min (value, pop_xmax (i) ),&
                  pop_xmin (i) ) )                                      
               ELSE 
                  pop_t (i, jt) = max (min (value, pop_xmax (i) ),      &
                  pop_xmin (i) )                                        
               ENDIF shift1
            ELSE 
!                                                                 
!     --------No Cross-Over, apply local shift only               
!                                                                 
 lshift1:      IF (abs (pop_lsig (i) ) .gt.0.0) then 
                  shift = gasdev (DBLE(pop_lsig (i))) 
               ELSE 
                  shift = 0.0 
               ENDIF lshift1
               value = pop_x (i, j) + shift 
!                                                                 
!-----------Upper/Lower limit correction                          
!                                                                 
 limit2:       IF (value.gt.pop_xmax (i) ) then 
                  value = pop_xmax (i) - abs (gasdev (DBLE(pop_lsig (i) ))) 
               ELSEIF (value.lt.pop_xmin (i) ) then 
                  value = pop_xmin (i) + abs (gasdev (DBLE(pop_lsig (i) ))) 
               ENDIF limit2
!                                                                 
!-----------Apply final shift                                     
!                                                                 
 shift2:       IF (pop_type (i) .eq.POP_INTEGER) THEN 
                  pop_t (i, jt) = nint (max (min (value, pop_xmax (i) ),&
                  pop_xmin (i) ) )                                      
               ELSE 
                  pop_t (i, jt) = max (min (value, pop_xmax (i) ),      &
                  pop_xmin (i) )                                        
               ENDIF shift2
            ENDIF cross1
            pop_para (i) = pop_t (i, jt) 
         ELSE 
            pop_t (i, jt) = pop_x (i, j) 
         ENDIF refine1
      ENDDO kdimx1
   ELSE local
  idimx2:   DO i = 1, pop_dimx 
 refine2:IF (pop_refine (i) ) then 
!                                                                 
!     --------No Cross-Over, apply local shift only               
!                                                                 
 lshift2:   IF (abs (pop_lsig (i) ) .gt.0.0) then 
               shift = gasdev (DBLE(pop_lsig (i))) 
            ELSE 
               shift = 0.0 
            ENDIF lshift2
            value = pop_x (i, j) + shift 
!                                                                 
!-----------Upper/Lower limit correction                          
!                                                                 
 limit3:    IF (value.gt.pop_xmax (i) ) then 
               value = pop_xmax (i) - abs (gasdev (DBLE(pop_lsig (i) ))) 
            ELSEIF (value.lt.pop_xmin (i) ) then 
               value = pop_xmin (i) + abs (gasdev (DBLE(pop_lsig (i) ))) 
            ENDIF limit3
!                                                                 
!-----------Apply final shift                                     
!                                                                 
 shift3:    IF (pop_type (i) .eq.POP_INTEGER) THEN 
               pop_t (i, jt) = nint (max (min (value, pop_xmax (i) ),   &
               pop_xmin (i) ) )
            ELSE 
               pop_t (i, jt)  = max (min (value, pop_xmax (i) ),  pop_xmin (i) ) 
            ENDIF shift3
            pop_para (i) = pop_t (i, jt) 
         ELSE refine2
            pop_t (i, jt) = pop_x (i, j) 
         ENDIF refine2
      ENDDO idimx2
   ENDIF local
   l_ok = .true. 
   DO l = 1, constr_number 
      line = ' ' 
      line = '('//constr_line (l) (1:constr_length (l) ) //')' 
      length = constr_length (l) + 2 
      l_ok = l_ok.and.if_test (line, length) 
   ENDDO 
ENDDO main
nlok: IF (.not.l_ok) then 
!                                                                 
!     ----Try a local modification, before giving up              
!                                                                 
   DO i = 1, pop_dimx 
      IF (pop_refine (i) ) then 
!                                                                 
!     --------No Cross-Over, apply local shift only               
!                                                                 
         IF (abs (pop_lsig (i) ) .gt.0.0) then 
            shift = gasdev (DBLE(pop_lsig (i))) 
         ELSE 
            shift = 0.0 
         ENDIF 
         value = pop_x (i, j) + shift 
!                                                                 
!-----------Upper/Lower limit correction                          
!                                                                 
         IF (value.gt.pop_xmax (i) ) then 
            value = pop_xmax (i) - abs (gasdev (DBLE(pop_lsig (i) ))) 
         ELSEIF (value.lt.pop_xmin (i) ) then 
            value = pop_xmin (i) + abs (gasdev (DBLE(pop_lsig (i) ))) 
         ENDIF 
!                                                                 
!-----------Apply final shift                                     
!                                                                 
         IF (pop_type (i) .eq.POP_INTEGER) THEN 
            pop_t (i, jt) = nint (max (min (value, pop_xmax (i) ),   &
            pop_xmin (i) ) )                                         
         ELSE 
            pop_t (i, jt)  = max (min (value, pop_xmax (i) ),  pop_xmin (i) ) 
         ENDIF 
         pop_para (i) = pop_t (i, jt) 
      ELSE 
         pop_t (i, jt) = pop_x (i, j) 
      ENDIF 
   ENDDO 
   l_ok = .true. 
   DO l = 1, constr_number 
      line = ' ' 
      line = '('//constr_line (l) (1:constr_length (l) ) //')' 
      length = constr_length (l) + 2 
      l_ok = l_ok.and.if_test (line, length) 
   ENDDO 
   IF ( l_ok ) THEN
      pop_current_trial = .true.
   ELSE
      ier_num = - 8 
      ier_typ = ER_APPL 
      RETURN 
   ENDIF 
ENDIF nlok
!                                                                       
END SUBROUTINE create_single                  
!*****7**************************************************************** 
SUBROUTINE write_trial (jt) 
!
! Writes the trial file of just one kid, and updates the ref_para array
!                                                                       
USE diff_evol
USE population
USE support_diffev_mod
USE param_mod
USE variable_mod
USE sys_compiler
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER, PARAMETER             :: iwr = 7
!                                                                       
INTEGER,           INTENT(IN)  :: jt 
!
CHARACTER (LEN=7)              :: stat  = 'unknown'
!                                                                       
INTEGER                        :: i
INTEGER                        :: len_file 
!                                                                       
DO i = 1, pop_dimx 
   ref_para(i) = pop_t (i, jt)
ENDDO 
var_val( var_ref+4) = jt
!
IF(pop_trial_file_wrt) THEN
   len_file = pop_ltrialfile 
   CALL make_file (pop_trialfile, len_file, 4, jt) 
   CALL do_del_file (pop_trialfile) 
   CALL oeffne (iwr, pop_trialfile, stat) 
!                                                                       
   WRITE (iwr, 2000) pop_gen, pop_n, pop_c, pop_dimx 
   WRITE (iwr, 2100) jt 
   WRITE (iwr, 2200) 
   DO i = 1, pop_dimx 
      WRITE (iwr, 3000) pop_t (i, jt), i, pop_name (i) 
   ENDDO 
   CLOSE (iwr) 
ENDIF
!                                                                       
 2000 FORMAT ('# generation members children parameters',/ 4(i8,2x)) 
 2100 FORMAT ('# current member ',/i5) 
 2200 FORMAT ('# parameter list') 
 3000 FORMAT (2x,e18.10,2x,i4,'  !  ',a) 
!                                                                       
END SUBROUTINE write_trial                    
!*****7**************************************************************** 
SUBROUTINE write_genfile 
!                                                                       
USE population
USE diff_evol
USE diffev_random
USE run_mpi_mod
!
USE random_state_mod
USE sys_compiler
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER, PARAMETER  :: iwr = 7
!                                                                       
CHARACTER (LEN=7)   :: stat = 'unknown'
CHARACTER (LEN=16)  :: string
INTEGER             :: i, j
INTEGER :: nseed_run    ! Actual number of seed used by compiler
!                                                                       
CALL do_del_file (pop_genfile) 
CALL oeffne (iwr, pop_genfile, stat) 
WRITE (iwr, 1000) pop_gen, pop_n, pop_c, pop_dimx 
WRITE (iwr, 1100) pop_trialfile (1:LEN_TRIM(pop_trialfile)) 
WRITE (iwr, 1200) trial_results (1:LEN_TRIM(trial_results))
WRITE (iwr, 1300) parent_results (1:LEN_TRIM(parent_results)) 
WRITE (iwr, 1400) parent_summary (1:LEN_TRIM(parent_summary))
WRITE (iwr, 1500) parent_current (1:LEN_TRIM(parent_current))
WRITE (iwr, 1510) pop_back_number
DO i= 1, pop_back_number
   WRITE (iwr, 1520) pop_back_fil(i)(1:LEN_TRIM(pop_back_fil(i)))
   WRITE (iwr, 1525) pop_back_ext(i)(1:LEN_TRIM(pop_back_ext(i)))
   WRITE (iwr, 1530) pop_back_trg(i)(1:LEN_TRIM(pop_back_trg(i)))
ENDDO
WRITE (iwr, 2000) diff_donor_mode, diff_sel_mode
WRITE (iwr, 2100) diff_cr, diff_f, diff_local, diff_k
!
nseed_run = random_nseeds()
random_nseed   = MIN(RUN_MPI_NSEEDS, nseed_run)  !  to be debugged depend on compiler ???
!
WRITE(IWR,3000) random_nseed
DO i = 0, (random_nseed-1)/6
   WRITE(IWR,3100) (random_best(j+1),j = i*6,MIN(i*6+5,random_nseed-1))
ENDDO
!
WRITE (IWR, 1640) 'XMIN','XMAX','SMIN','SMAX'
string='Member'
WRITE (IWR, 1800) string(1:16),-1, 1.,REAL(pop_dimx), &
                                   1., REAL(pop_dimx),&
                                   'I', .FALSE.
string='Rvalue'
WRITE (IWR, 1800) string(1:16), 0, parent_val(pop_best,0),&
                               parent_val(pop_worst,0),&
                               child_val(pop_best,0),&
                               child_val(pop_worst,0),&
 'I', .FALSE.
DO i = 1, pop_dimx
   string = pop_name(i)(1:LEN_TRIM(pop_name(i)))
!  WRITE(IWR, 1800) pop_name(i)(1:LEN_TRIM(pop_name(i))), &
!  WRITE(IWR, 1800) string(1:16)                        , &
   WRITE(IWR, 1800) pop_name(i)(1:16)                   , &
                  i,pop_xmin(i), pop_xmax(i),             &
                    pop_smin(i), pop_smax(i), 'R',        &
                    pop_refine(i)
ENDDO
CLOSE (iwr) 
!                                                                       
 1000 FORMAT ('# generation members children parameters',/ 4(i8,2x)) 
 1100 FORMAT ('# trial file'/a) 
 1200 FORMAT ('# result file'/a) 
 1300 FORMAT ('# log file'/a) 
 1400 FORMAT ('# summary file'/a) 
 1500 FORMAT ('# last    file'/a) 
 1510 FORMAT ('# backup number'/1x,i4) 
 1520 FORMAT ('# backup source'/a) 
 1525 FORMAT ('# backup extension'/a) 
 1530 FORMAT ('# backup target'/a) 
 1640 FORMAT ('# Parameter      no.',4(1x,a18))
 1800 FORMAT (a16,1x,i3,4(1x,E18.10),1x,a1,1x,l1)
 3000 FORMAT ('# random seeds '/i5) 
 3100 FORMAT (5(1x,i11),1x,i11)
!
 2000 FORMAT ('# donor mode, selection mode'/2(2x,i2))
 2100 FORMAT ('# Cross_over, Factor, Local, position'/,4(2x, E18.10)) 
!
END SUBROUTINE write_genfile                  
!*****7**************************************************************** 
SUBROUTINE read_genfile 
!                                                                       
USE diffev_allocate_appl
USE diff_evol
USE population
USE diffev_random
!
USE errlist_mod
USE lib_f90_allocate_mod
USE variable_mod
USE sys_compiler
!                                                                       
IMPLICIT none 
!
!                                                                       
INTEGER, PARAMETER  :: iwr = 7
!                                                                       
CHARACTER (LEN=7)   :: stat = 'old'
CHARACTER (LEN=128) :: line
INTEGER             :: io_status
INTEGER             :: i, j
INTEGER             :: n_pop, n_dim
INTEGER             :: m_pop_gen, m_pop_n, m_pop_c, m_pop_dimx
LOGICAL             :: l_write
REAL                :: r1, r2, r3, r4 
REAL                :: highest_r
!                                                                       
CALL oeffne (iwr, pop_genfile, stat) 
IF ( ier_num/=0) THEN
   ier_msg(1) = 'Could not open the GENERATION file'
   CLOSE ( iwr)
   RETURN
ENDIF
ier_num = -31                                       ! Turn error on for read
ier_typ = ER_APPL
!
READ (iwr, *   ,iostat=IO_status)                   ! Read 1st header line
IF(IS_IOSTAT_END(IO_STATUS)) THEN
   CLOSE(iwr)
   RETURN
ENDIF
!
READ (iwr, *   ,iostat=IO_status) r1, r2, r3, r4    ! Read generation values
IF(IS_IOSTAT_END(IO_STATUS)) THEN
   CLOSE(iwr)
   RETURN
ENDIF
!
READ (iwr, *   ,iostat=IO_status)                   ! Read trialfile header line
IF(IS_IOSTAT_END(IO_STATUS)) THEN
   CLOSE(iwr)
   RETURN
ENDIF
!
READ (iwr, 1000,iostat=IO_status) pop_trialfile     ! Read trialfile name
pop_ltrialfile = LEN_TRIM(pop_trialfile)
!
READ (iwr, *   ,iostat=IO_status)                   ! read reaultfile header line
IF(IS_IOSTAT_END(IO_STATUS)) THEN
   CLOSE(iwr)
   RETURN
ENDIF
!
READ (iwr, 1000,iostat=IO_status) trial_results     ! Read resultfile name 
IF(IS_IOSTAT_END(IO_STATUS)) THEN
   CLOSE(iwr)
   RETURN
ENDIF
ltrial_results = LEN_TRIM(trial_results)
!
READ (iwr, *   ,iostat=IO_status)                   ! read logfile header line
IF(IS_IOSTAT_END(IO_STATUS)) THEN
   CLOSE(iwr)
   RETURN
ENDIF
!
READ (iwr, 1000,iostat=IO_status) parent_results    ! Read parent log file name 
IF(IS_IOSTAT_END(IO_STATUS)) THEN
   CLOSE(iwr)
   RETURN
ENDIF
lparent_results = LEN_TRIM(parent_results)
!
READ (iwr, *   ,iostat=IO_status)                   ! read summary header line
IF(IS_IOSTAT_END(IO_STATUS)) THEN
   CLOSE(iwr)
   RETURN
ENDIF
!
READ (iwr, 1000,iostat=IO_status) parent_summary    ! Read summary file name 
IF(IS_IOSTAT_END(IO_STATUS)) THEN
   CLOSE(iwr)
   RETURN
ENDIF
lparent_summary     = LEN_TRIM(parent_summary)
!
ier_num = 0                                         ! Turn error message off as 
ier_typ = 0                                         ! old file finish here
READ (iwr, *   ,iostat=IO_status)                   ! read current header line
IF(IS_IOSTAT_END(IO_STATUS)) THEN
   CLOSE(iwr)
   RETURN
ENDIF
!
ier_num = -31                                       ! Turn error on for read
ier_typ = ER_APPL
READ (iwr, 1000,iostat=IO_status) parent_current    ! Read current file name
IF(IS_IOSTAT_END(IO_STATUS)) THEN
   CLOSE(iwr)
   RETURN
ENDIF
lparent_current = LEN_TRIM(parent_current)
!
READ (iwr, 1000,iostat=IO_status) line              ! read backup header line
IF(IS_IOSTAT_END(IO_STATUS)) THEN
   CLOSE(iwr)
   RETURN
ENDIF
!
READ (iwr, *   ,iostat=IO_status) pop_back_number  ! read number of backups
IF(IS_IOSTAT_END(IO_STATUS)) THEN
   CLOSE(iwr)
   RETURN
ENDIF
!
DO i=1, pop_back_number
   READ (iwr, *   ,iostat=IO_status) 
   IF(IS_IOSTAT_END(IO_STATUS)) THEN
      CLOSE(iwr)
      RETURN
   ENDIF
!
   READ (iwr, 1000,iostat=IO_status) pop_back_fil(i) 
   IF(IS_IOSTAT_END(IO_STATUS)) THEN
      CLOSE(iwr)
      RETURN
   ENDIF
!
   READ (iwr, *   ,iostat=IO_status) 
   IF(IS_IOSTAT_END(IO_STATUS)) THEN
      CLOSE(iwr)
      RETURN
   ENDIF
!
   READ (iwr, 1000,iostat=IO_status) pop_back_ext(i) 
   IF(IS_IOSTAT_END(IO_STATUS)) THEN
      CLOSE(iwr)
      RETURN
   ENDIF
!
   READ (iwr, *   ,iostat=IO_status) 
   IF(IS_IOSTAT_END(IO_STATUS)) THEN
      CLOSE(iwr)
      RETURN
   ENDIF
!
   READ (iwr, 1000,iostat=IO_status) pop_back_trg(i) 
   IF(IS_IOSTAT_END(IO_STATUS)) THEN
      CLOSE(iwr)
      RETURN
   ENDIF
!
   pop_back_fil_l(i) = LEN_TRIM(pop_back_fil(i))
   pop_back_ext_l(i) = LEN_TRIM(pop_back_ext(i))
   pop_back_trg_l(i) = LEN_TRIM(pop_back_trg(i))
ENDDO
!
IF(pop_back_number>0) pop_backup = .TRUE.
READ (iwr, *   ,iostat=IO_status)                  ! Read donor  header line
IF(IS_IOSTAT_END(IO_STATUS)) THEN
   CLOSE(iwr)
   RETURN
ENDIF
!
READ (iwr, *, IOSTAT=io_status   ) diff_donor_mode, diff_sel_mode
READ (iwr, *   ,iostat=IO_status)                  ! Read CR, F  header line
IF(IS_IOSTAT_END(IO_STATUS)) THEN
   CLOSE(iwr)
   RETURN
ENDIF
!
READ (iwr, *, IOSTAT=io_status   ) diff_cr, diff_f, diff_local, diff_k
IF(IS_IOSTAT_END(IO_STATUS)) THEN
   CLOSE(iwr)
   RETURN
ENDIF
!
! Read random status of best element
!
line = ' '
j = 0
READ(iwr,'(a)', IOSTAT=io_status) line
IF(IS_IOSTAT_END(IO_STATUS)) THEN
   CLOSE(iwr)
   RETURN
ENDIF
!
!
READ (iwr, *, IOSTAT=io_status   ) random_nseed
IF(IS_IOSTAT_END(IO_STATUS)) THEN
   CLOSE(iwr)
   RETURN
ENDIF
!
!
DO i = 0, (random_nseed-1)/6
   READ(iwr,'(a)', IOSTAT=io_status) line
   IF(IS_IOSTAT_END(IO_STATUS)) THEN
      CLOSE(iwr)
      RETURN
   ENDIF
!
   READ(line,*, IOSTAT=IO_status) (random_best(j+1),j = i*6,MIN(i*6+5,random_nseed-1))
   IF(IS_IOSTAT_END(IO_STATUS)) THEN
      CLOSE(iwr)
      RETURN
   ENDIF
!
ENDDO
!READ (iwr, *   ,iostat=IO_status) random_best(1:random_nseed)
!
READ (iwr, *   ,iostat=IO_status) 
IF(IS_IOSTAT_END(IO_STATUS)) THEN
   CLOSE(iwr)
   RETURN
ENDIF
!
READ (iwr, *   ,iostat=IO_status)
IF(IS_IOSTAT_END(IO_STATUS)) THEN
   CLOSE(iwr)
   RETURN
ENDIF
!
READ (iwr, *   ,iostat=IO_status)
IF(IS_IOSTAT_END(IO_STATUS)) THEN
   CLOSE(iwr)
   RETURN
ENDIF
!
l_write = .FALSE.              ! Generation file does not have to be written
IF(.NOT.pop_current)   THEN    ! Need to update the population etc
   m_pop_gen  = MAX(pop_gen , NINT(r1))
   m_pop_n    = MAX(pop_n   , NINT(r2))
   m_pop_c    = MAX(pop_c   , NINT(r3))
   m_pop_dimx = MAX(pop_dimx, NINT(r4))
   IF(m_pop_n>MAXPOP .OR. m_pop_c>MAXPOP .OR. m_pop_dimx> MAXDIMX) THEN   ! Need to allocate
      n_pop = MAX(m_pop_n, m_pop_c, MAXPOP)
      n_dim = MAX(m_pop_dimx, MAXDIMX)
      CALL alloc_population( n_pop, n_dim)
      CALL alloc_ref_para(  n_dim)
      IF(ier_num < 0) THEN
         CLOSE(iwr)
         RETURN
      ENDIF
   ENDIF
   IF ( pop_gen > 0 ) THEN
      pop_dimx_new = .true.
   ENDIF
   IF ( pop_gen == 0 ) THEN
      pop_dimx_new = .false.
   ENDIF
   pop_dimx_init = .true.      ! The dimension has been initialized in this run
   DO i=1,NINT(r4)
      READ(iwr,'(a)', IOSTAT=io_status) line
      IF(IS_IOSTAT_END(IO_STATUS)) THEN
         ier_num = -17
         ier_typ = ER_APPL
         ier_msg(1) = 'Number of parameters differs from current value'
         ier_msg(2) = 'Check value of pop_dimx[1], and adjust         '
         ier_msg(3) = 'Conflicting values within GENERATION file      '
         CLOSE(iwr)
         RETURN
      ENDIF
!
      READ(line(1:16), '(a)', IOSTAT=io_status) pop_name(i)
      IF(IS_IOSTAT_END(IO_STATUS)) THEN
         ier_num = -17
         ier_typ = ER_APPL
         ier_msg(1) = 'Error reading parameter name in GENERATION file'
         ier_msg(3) = 'Conflicting values within GENERATION file      '
         CLOSE(iwr)
         RETURN
      ENDIF
!
      READ(line(21:96),*,IOSTAT=io_status) pop_xmin(i), pop_xmax(i), &
                                           pop_smin(i), pop_smax(i)
      IF(IS_IOSTAT_END(IO_STATUS)) THEN
         ier_num = -17
         ier_typ = ER_APPL
         ier_msg(1) = 'Error reading parameter limits in GENERATION file'
         ier_msg(3) = 'Conflicting values within GENERATION file      '
         CLOSE(iwr)
         RETURN
      ENDIF
!
      pop_type(i) = POP_REAL
      IF(line(98:98) == 'I') pop_type(i) = POP_INTEGER
      READ(line(100:100),'(L1)', IOSTAT=io_status) pop_refine(i)
      IF(IS_IOSTAT_END(IO_STATUS)) THEN
         ier_num = -17
         ier_typ = ER_APPL
         ier_msg(1) = 'Error reading parameter status in GENERATION file'
         ier_msg(3) = 'Conflicting values within GENERATION file      '
         CLOSE(iwr)
         RETURN
      ENDIF
!
   ENDDO
!
!  Test for changes in population size
!
   IF(pop_not_first .AND. NINT(r1)>0 .AND. NINT(r1)>= pop_gen) THEN   ! We are in a GENERATION > zero
      pop_gen = NINT(r1)
      IF(NINT(r2)<pop_n) THEN   ! Population size increased  by pop_n[] command
               highest_r = MAXVAL( parent_val(1:NINT(r1),0))
               parent_val(pop_n+1:pop_n,0) = 10.*highest_r
               FORALL(i=pop_n+1:pop_n) pop_x (:,i)     = pop_x (:,1)
               r2 = REAL(pop_n) 
               var_val(var_ref+1) = pop_n     ! Update global user variable
               l_write = .TRUE.
      ELSEIF(NINT(r2)>pop_n) THEN   ! Population size decreased  by pop_n[] command
               r2 = REAL(pop_n) 
               var_val(var_ref+1) = pop_n     ! Update global user variable
               l_write = .TRUE.
      ENDIF
      IF(NINT(r3)/=pop_c) THEN   ! Children size changed by pop_c[} command
               r3 = REAL(pop_c)
               var_val(var_ref+2) = pop_c     ! Update global user variable
               CALL create_trial                 ! Make a new set
               l_write = .TRUE.
      ENDIF
      IF(NINT(r4)<pop_dimx) THEN   ! Dimension has  increased by pop_dimx[] command
               r4 = REAL(pop_dimx)
               l_write = .TRUE.
      ENDIF
   ELSE
      pop_gen  = NINT(r1)
      pop_n    = NINT(r2)
      pop_c    = NINT(r3)
      pop_dimx = NINT(r4)
!      CALL alloc_ref_para(MAXDIMX)
   ENDIF
ENDIF
CLOSE (iwr) 
IF(l_write) THEN
   CALL write_genfile                ! Write the "GENERATION" file
   l_write = .FALSE.
ENDIF
! In case io_status is zero, reset the error message
IF(io_status==0) THEN
   ier_num = 0
   ier_typ = 0
ENDIF
!
IF(io_status /= 0) THEN
   ier_num = -19
   ier_typ = ER_APPL
   ier_msg(1) = 'Could not read the GENERATION file'
   RETURN
ENDIF
!                                                                       
IF ( pop_n >  NINT (r2) ) THEN
   ier_num = -16
   ier_typ = ER_APPL
   ier_msg(1) = 'Number of members differs from current value'
   ier_msg(2) = 'Check value of pop_n[1], and adjust         '
   RETURN
ENDIF
IF ( pop_c >  NINT (r3) ) THEN
   ier_num = -16
   ier_typ = ER_APPL
   ier_msg(1) = 'Number of children differs from current value'
   ier_msg(2) = 'Check value of pop_c[1], and adjust          '
   RETURN
ENDIF
IF ( pop_dimx /= NINT (r4) ) THEN
   ier_num = -17
   ier_typ = ER_APPL
   ier_msg(1) = 'Number of parameters differs from current value'
   ier_msg(2) = 'Check value of pop_dimx[1], and adjust         '
   RETURN
ENDIF
pop_gen  = nint (r1) 
var_val(var_ref+0) = pop_gen  ! Update global user variable
!pop_n    = nint (r2) 
!pop_c    = nint (r3) 
!pop_dimx = nint (r4) 
                                                                        
!                                                                       
 1000 FORMAT    (a) 
!                                                                       
END SUBROUTINE read_genfile                   
!*****7**************************************************************** 
SUBROUTINE adapt_sigma 
!                                                                       
USE diff_evol
USE population
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER i, j 
!                                                                       
DO j = 1, pop_dimx 
   pop_pmin (j) = pop_x (j, 1) 
   pop_pmax (j) = pop_x (j, 1) 
   DO i = 2, pop_n 
      pop_pmin (j) = min (pop_pmin (j), pop_x (j, i) ) 
      pop_pmax (j) = max (pop_pmax (j), pop_x (j, i) ) 
   ENDDO 
ENDDO 
!                                                                       
IF (pop_gen.gt.0) then 
   DO j = 1, pop_dimx 
      IF (pop_ad_sigma (j) ) then 
         pop_sigma (j) = pop_sig_ad (j) * (abs (pop_pmax (j) - pop_pmin (j) ) )                                            
      ENDIF 
   ENDDO 
   DO j = 1, pop_dimx 
      IF (pop_ad_lsigma (j) ) then 
         pop_lsig (j) = pop_lsig_ad (j) * (abs (pop_pmax (j) - pop_pmin (j) ) )                                            
      ENDIF 
   ENDDO 
ENDIF 
!                                                                       
END SUBROUTINE adapt_sigma                    
END MODULE create_trial_mod
