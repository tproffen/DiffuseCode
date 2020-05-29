MODULE compare
!
USE errlist_mod
USE sorting_mod
!
PRIVATE
PUBLIC  :: do_compare
PUBLIC  :: do_dismiss
PUBLIC  :: do_read_values
PUBLIC  :: patch_para
PUBLIC  :: read_par_values
PUBLIC  :: write_current
PUBLIC  :: write_parents
!
CONTAINS
!
!*****7**************************************************************** 
!                                                                       
   SUBROUTINE do_compare 
!
! Main comparison routine, branches into the two modes
! compare to parent
! compare to best of all (membeers + children)
!
   USE diff_evol
   USE diffev_random
   USE population
!                                                                       
   IMPLICIT none 
!                                                                       
   IF (diff_sel_mode.eq.SEL_COMP) THEN 
      CALL compare_toparent 
   ELSEIF (diff_sel_mode.eq.SEL_BEST_ALL) THEN 
      CALL compare_best_all 
   ENDIF 
   CALL diffev_best_macro
   lastgen = pop_gen - 1
   END SUBROUTINE do_compare                        
!*****7**************************************************************** 
   SUBROUTINE compare_toparent 
!
! compares all children to their immediate parent
!
   USE diff_evol
   USE diffev_random
   USE population
USE sys_compiler
!                                                                       
   IMPLICIT none 
!
   CHARACTER(LEN=600)             :: string 
   INTEGER                        :: i, j , nb
   REAL                           ::  best, worst 
!                                                                       
   CALL do_read_values(.TRUE.)
   IF ( ier_num /=0) RETURN
!                                                                       
!     Selection mode: Compare each child to its parent                  
!                                                                       
   during: IF (pop_gen.gt.0) THEN 
      DO j = 1, pop_n 
         IF (trial_val (j,0) .lt.parent_val (j,0) ) THEN 
            DO i = 1, pop_dimx 
               child (i, j) = trial (i, j) 
               child_val (j,0:n_rvalue_i) = trial_val (j,0:n_rvalue_i)
            ENDDO 
            bck_during: IF(pop_backup) THEN    ! copy current best calculations into backup 
                                               ! This effectively replaces kup.backup.mac
               IF(j<=pop_c) THEN
               DO nb = 1, pop_back_number
                  WRITE(string,1000) pop_back_fil(nb)(1:pop_back_fil_l(nb)), j, &
                                     pop_back_ext(nb)(1:pop_back_ext_l(nb)),    &
                                     pop_back_trg(nb)(1:pop_back_trg_l(nb)), j, &
                                     pop_back_ext(nb)(1:pop_back_ext_l(nb))
                  CALL do_operating_comm(string)
               ENDDO
               ENDIF
            ENDIF bck_during
         ELSE 
            DO i = 1, pop_dimx 
               child (i, j) = pop_x (i, j) 
               child_val (j,0:n_rvalue_i) = parent_val (j,0:n_rvalue_i) 
            ENDDO 
         ENDIF 
      ENDDO 
   ELSE during
      DO j = 1, pop_n 
         DO i = 1, pop_dimx 
            child (i, j) = trial (i, j) 
            child_val (j,0:n_rvalue_i) = trial_val (j,0:n_rvalue_i) 
         ENDDO 
         bck_prior: IF(pop_backup) THEN    ! copy current best calculations into backup 
                                            ! This effectively replaces kup.backup.mac
            IF(j<=pop_c) THEN
            DO nb = 1, pop_back_number
               WRITE(string,1000) pop_back_fil(nb)(1:pop_back_fil_l(nb)), j, &
                                  pop_back_ext(nb)(1:pop_back_ext_l(nb)),    &
                                  pop_back_trg(nb)(1:pop_back_trg_l(nb)), j, &
                                  pop_back_ext(nb)(1:pop_back_ext_l(nb))
               CALL do_operating_comm(string)
               IF ( ier_num /=0) RETURN
            ENDDO
            ENDIF
         ENDIF bck_prior
      ENDDO 
   ENDIF during
!                                                                       
!     --Determine best and worst member                                 
!                                                                       
   best = child_val (1,0) + 1.0 
   worst = child_val (1,0) - 1.0 
   DO j = 1, pop_n 
      IF (child_val (j,0) .lt.best) THEN 
         best     = child_val(j,0)
         pop_best = j 
      ENDIF 
      IF (child_val (j,0) .gt.worst) THEN 
         worst    = child_val(j,0)
         pop_worst = j 
      ENDIF 
   ENDDO 
!
!  Archive random state for best member
!
   CALL diffev_random_save(pop_random(:,pop_best))
!                                                                       
!------ write the parameters and the results for the current generation 
!------ Copy current child into parent parameters and create new trial  
!     parameters                                                        
!                                                                       
   CALL write_parents 
!                                                                       
1000 FORMAT('cp ',a,i4.4,a,' ',a,i4.4,a)
!
   END SUBROUTINE compare_toparent               
!*****7**************************************************************** 
   SUBROUTINE compare_best_all 
!
! compares all children to the combined group of (parent + children)
!
   USE diff_evol
   USE diffev_random
   USE population
   USE random_mod
USE sys_compiler
!
   IMPLICIT none 
!                                                                       
!                                                                       
   CHARACTER(LEN=600)             :: string 
   INTEGER                        :: list_number 
   INTEGER                        :: list_index (2 * MAXPOP) 
   REAL                           :: list_val (2 * MAXPOP) 
   REAL                           :: best, worst
!                                                                       
   INTEGER                        :: i, j, k, ii , nb
!
!                                                                       
   CALL do_read_values(.TRUE.)
   IF ( ier_num /=0) RETURN
!                                                                       
!     Selection mode: use the pop_n best of parents and children        
!                                                                       
!                                                                       
!     Create complete list of all r-values                              
!                                                                       
   list_val(:) = 0.0
   during: IF (pop_gen.gt.0) THEN 
      DO j = 1, pop_n 
         list_val (j) = parent_val (j,0) 
      ENDDO 
      best  = trial_val(1,0)
      worst = trial_val(1,0)
      DO j = 1, pop_c 
         list_val (pop_n + j) = trial_val (j,0) 
         best = MIN(best , trial_val(j,0))
         worst= MAX(worst, trial_val(j,0))
      ENDDO 
      DO j = pop_c+1, pop_n
         list_val(j) = worst + (worst-best)
      ENDDO
      list_number = pop_n + pop_c 
   ELSE during
      best  = trial_val(1,0)
      worst = trial_val(1,0)
      DO j = 1, pop_c 
         list_val (j) = trial_val (j,0) 
         best = MIN(best , trial_val(j,0))
         worst= MAX(worst, trial_val(j,0))
      ENDDO 
      DO j = pop_c+1, pop_n
         list_val(j) = worst + (worst-best)
      ENDDO
      list_number = pop_c 
   ENDIF during
!                                                                       
!     heapsort index array  on r-values                                 
!                                                                       
list_index(:) = 0
   CALL indexx (list_number, list_val, list_index) 
!                                                                       
!     copy the pop_n best into the child variables                      
!                                                                       
   copy: IF (pop_gen.gt.0) THEN 
      bck_during: IF(pop_backup) THEN    ! copy current best calculations into backup 
                                         ! This effectively replaces kup.backup.mac
         DO k = pop_n ,1, -1
            IF (list_index (k) .le.pop_n) THEN 
!                                                                       
!     ------ an old parent                                                 
!                                                                       
               ii = list_index (k) 
               IF(ii /= k) THEN
                  IF(ii> 0 .AND. ii<=pop_c) THEN
                  DO nb = 1, pop_back_number
                     WRITE(string,1000) pop_back_trg(nb)(1:pop_back_trg_l(nb)),ii, &
                                        pop_back_ext(nb)(1:pop_back_ext_l(nb)),    &
                                        pop_back_trg(nb)(1:pop_back_trg_l(nb)), k, &
                                        pop_back_ext(nb)(1:pop_back_ext_l(nb))
                     CALL do_operating_comm(string)
                     IF ( ier_num /=0) RETURN
                  ENDDO
                  ENDIF
               ENDIF
            ELSE 
!                                                                       
!     ------ a child                                                       
!                                                                       
               ii = list_index (k) - pop_n 
               IF(ii> 0 .AND. ii<=pop_c) THEN
               DO nb = 1, pop_back_number
                  WRITE(string,1000) pop_back_fil(nb)(1:pop_back_fil_l(nb)),ii, &
                                     pop_back_ext(nb)(1:pop_back_ext_l(nb)),    &
                                     pop_back_trg(nb)(1:pop_back_trg_l(nb)), k, &
                                     pop_back_ext(nb)(1:pop_back_ext_l(nb))
                  CALL do_operating_comm(string)
               ENDDO
               ENDIF
            ENDIF
         ENDDO
      ENDIF bck_during
      DO k = 1, pop_n 
         IF (list_index (k) .le.pop_n) THEN 
!                                                                       
!     --- an old parent                                                 
!                                                                       
            ii = list_index (k) 
            DO i = 1, pop_dimx 
               child (i, k) = pop_x (i, ii) 
               child_val (k,0:n_rvalue_i) = parent_val (ii,0:n_rvalue_i) 
            ENDDO 
         ELSE 
!                                                                       
!     --- a child                                                       
!                                                                       
            ii = list_index (k) - pop_n 
            DO i = 1, pop_dimx 
               child (i, k) = trial (i, ii) 
               child_val (k,0:n_rvalue_i) = trial_val (ii,0:n_rvalue_i) 
            ENDDO 
         ENDIF 
      ENDDO 
      k = 1
      IF (.NOT.(list_index (k) .le.pop_n)) THEN 
!
!     This is a new "best", update random state
!
         ii = list_index (k) - pop_n 
         CALL diffev_random_save(pop_random(:,ii))
      ENDIF
   ELSE copy
      bck_prior: IF(pop_backup) THEN    ! copy current best calculations into backup 
         DO k = pop_n ,1, -1
            ii = list_index (k) 
            IF(k>pop_c) ii = list_index(list_number)  ! Not enough children copy worst
            IF(ii> 0 .AND. ii<=pop_c) THEN
            DO nb = 1, pop_back_number
               WRITE(string,1000) pop_back_fil(nb)(1:pop_back_fil_l(nb)),ii, &
                                  pop_back_ext(nb)(1:pop_back_ext_l(nb)),    &
                                  pop_back_trg(nb)(1:pop_back_trg_l(nb)), k, &
                                  pop_back_ext(nb)(1:pop_back_ext_l(nb))
               CALL do_operating_comm(string)
            ENDDO
            ENDIF
         ENDDO
      ENDIF bck_prior
      DO k = 1, pop_n 
         ii = list_index (k) 
         DO i = 1, pop_dimx 
            child (i, k) = trial (i, ii) 
            child_val (k,0:n_rvalue_i) = trial_val (ii,0:n_rvalue_i) 
         ENDDO 
      ENDDO 
      k = 1
      ii = list_index (1)
      CALL diffev_random_save(pop_random(:,ii))
   ENDIF copy
!                                                                       
!     determine best/worst member                                       
!                                                                       
   pop_best = 1 
   pop_worst = pop_n 
!                                                                       
!     write parents, copy current children into parents and             
!     create new trial file                                             
!                                                                       
   CALL write_parents 
!
1000 FORMAT('cp ',a,i4.4,a,' ',a,i4.4,a)
!                                                                       
   END SUBROUTINE compare_best_all               
!*****7**************************************************************** 
   SUBROUTINE read_obj_values 
!
! Read result of the cost function evaluation
!
   USE population
   USE support_diffev_mod
USE sys_compiler
!                                                                       
   IMPLICIT none 
!
!                                                                       
   INTEGER, PARAMETER             :: iwr = 7
!                                                                       
   CHARACTER (LEN=7)              :: stat = 'unknown'
   INTEGER                        :: j 
   INTEGER                        :: iostatus
   INTEGER                        :: len_file 
   REAL                           :: r 
!                                                                       
   silent: IF(pop_result_file_rd) THEN
   DO j = 1, pop_c 
      len_file = ltrial_results 
      CALL make_file (trial_results, len_file, 4, j) 
      CALL oeffne (iwr, trial_results, stat) 
      READ (iwr, * ,iostat=iostatus) r, trial_val (j,0) 
      CLOSE (iwr) 
      IF ( iostatus /= 0) THEN
         ier_num = -11
         ier_typ = ER_APPL
         ier_msg(1) = 'Error while reading'
         WRITE (ier_msg(2),2000) j
         RETURN
      ENDIF
   ENDDO 
   ENDIF silent
!
   2000 FORMAT ('Child No. ',i4)
!                                                                       
   END SUBROUTINE read_obj_values                
!*****7**************************************************************** 
   SUBROUTINE read_par_values 
!
! Read the logfile
!
   USE population
   USE create_trial_mod
   USE support_diffev_mod
!
USe lib_length
   USE precision_mod
   USE prompt_mod
   USE random_mod
   USE terminal_mod
USE sys_compiler
!
!                                                                       
   IMPLICIT none 
!
!                                                                       
   INTEGER, PARAMETER             :: iwr = 7 
!                                                                       
   CHARACTER (LEN=7)              :: stat = 'unknown'
   CHARACTER (LEN=PREC_STRING)    :: line 
   CHARACTER (LEN=PREC_STRING)    :: fname
   INTEGER                        :: j, i, ii , k
   INTEGER                        :: len_file,length 
   INTEGER                        :: pop_dimx_old
   INTEGER                        :: iostatus, isuccess
   INTEGER                        :: ieq
   LOGICAL                        :: istda, lcurrent
   REAL                           :: best, worst 
   REAL                           :: r
   REAL                           :: temp_val_min, temp_val_max
   REAL                           :: temp_pop_min, temp_pop_max
!                                                                       
!
   iostatus = 0
   temp_val_min = 0.0
   temp_val_max = 0.0
!                                                                       
! Read old Parent value, if not yet initialized                     
!                                                                       
   init: IF (.not.pop_current.and.pop_gen.gt.0) THEN 
!  init: IF (                     pop_gen.gt.0) THEN 
!
!     loop over dimension to find old dimension
!
      length       = len_str(parent_current)
      IF(length == 0 ) THEN                    ! Current file not defined
         length       = len_str(parent_results)
         lcurrent = .false.
      ELSE
         lcurrent = .true.
      ENDIF
      pop_dimx_old = 0
      search0: DO i = 1, pop_dimx
         IF(pop_name(i)=='PARA0000') CYCLE search0
         IF(lcurrent) THEN
            WRITE(fname, 950) parent_current(1:length), pop_name(i)(1:LEN_TRIM(pop_name(i)))   ! Use current file
         ELSE
            WRITE(fname, 950) parent_results(1:length), pop_name(i)(1:LEN_TRIM(pop_name(i)))   ! Use full parameter file
         ENDIF
         INQUIRE ( FILE=fname, EXIST=istda )   ! does file exist?
         IF ( .NOT. istda) EXIT                ! if not exit loop
         pop_dimx_old = i                      ! current old dimension
      ENDDO search0
      search_old: DO i = 1, pop_dimx_old                   ! Loop over old dimension
         IF(pop_name(i)=='PARA0000') CYCLE search_old
         IF(lcurrent) THEN
            WRITE(fname, 950) parent_current(1:length), pop_name(i)(1:LEN_TRIM(pop_name(i)))   ! Use current file
         ELSE
            WRITE(fname, 950) parent_results(1:length), pop_name(i)(1:LEN_TRIM(pop_name(i)))   ! Use full parameter file
         ENDIF
         CALL oeffne (iwr, fname, stat) 
         ii = - 1 
         DO while (ii.ne.pop_gen - 1)          ! Loop over all previous generations
            READ (iwr, '(a)' ,END=30,ERR=30,iostat=iostatus) line 
            DO while (line (2:2) .ne.'S') 
               READ (iwr, '(a)' ,END=30,ERR=30,iostat=iostatus) line 
            ENDDO 
            ieq = INDEX(line,'=')
            READ (line (3:ieq-2), *  ,END=30,ERR=30,iostat=iostatus) ii 
         ENDDO 
         READ (iwr, '(a)' ,END=30,ERR=30,iostat=iostatus) line 
         isuccess =  0 
         get_params: DO j = 1, pop_c                       ! Read all old parameters
            READ (iwr, *     ,iostat=iostatus) &
              &   ii, parent_val (j,0), pop_x (i, j)
            child_val(j,0) = parent_val(j,0)
            IF(IS_IOSTAT_END(iostatus)) THEN 
               isuccess = -13
               EXIT get_params
            ELSEIF(iostatus/=0 .OR. ISNAN(parent_val (j,0)) .OR. ISNAN(pop_x (i, j)) ) THEN
               isuccess= -13
               ier_num = -13
               ier_typ = ER_APPL
               WRITE(error_io,'(a,a)') ' ***APPL*** ',fname(1:LEN_TRIM(fname))
               WRITE(error_io,'(a)')   ' ***APPL*** Generation,    member,     child'
               WRITE(error_io,'(a  ,3(2x,I9))') ' ***APPL***',pop_gen-1,i,j
               WRITE(error_io,'(a,2(2x,G16.10E2))') ' ***APPL*** ',parent_val(j,0), pop_x(i,j)
               CLOSE(iwr)
               RETURN
            ENDIF
         ENDDO get_params
         IF(isuccess == -13) THEN   ! Premature END of FILE assume population has increased
            temp_val_min = MINVAL(parent_val(:,0))
            temp_val_max = MAXVAL(parent_val(:,0))
            temp_pop_min = pop_x(i,1)
            temp_pop_max = pop_x(i,1)
            DO ii = 1, j-1
               temp_pop_min = MIN(temp_pop_min,pop_x(i,ii))
               temp_pop_max = MAX(temp_pop_min,pop_x(i,ii))
            ENDDO
            DO ii = j, pop_c        ! Create new parameters
               parent_val(ii,0) = temp_val_max + temp_val_min
               CALL RANDOM_NUMBER(r)
               pop_x(i,ii)    = temp_pop_min + (temp_pop_max-temp_pop_min)*r
            ENDDO
         ENDIF
      ENDDO search_old
      IF(n_rvalue_i>1) THEN     !Read partial Rvalues
         i= 0  
      search_partial: DO k = 1, n_rvalue_i                   ! Loop over partial Rvalues
         IF(pop_name(i)=='PARA0000') CYCLE search_partial
         IF(lcurrent) THEN
            WRITE(fname, 970) parent_current(1:length), pop_name(i)(1:LEN_TRIM(pop_name(i))),k   ! Use current file
         ELSE
            WRITE(fname, 970) parent_results(1:length), pop_name(i)(1:LEN_TRIM(pop_name(i))),k   ! Use full parameter file
         ENDIF
         CALL oeffne (iwr, fname, stat) 
         ii = - 1 
         DO while (ii.ne.pop_gen - 1)          ! Loop over all previous generations
            READ (iwr, '(a)' ,END=30,ERR=30,iostat=iostatus) line 
            DO while (line (2:2) .ne.'S') 
               READ (iwr, '(a)' ,END=30,ERR=30,iostat=iostatus) line 
            ENDDO 
            ieq = INDEX(line,'=')
            READ (line (3:ieq-2), *  ,END=30,ERR=30,iostat=iostatus) ii 
         ENDDO 
         READ (iwr, '(a)' ,END=30,ERR=30,iostat=iostatus) line 
         isuccess =  0 
         get_params_k: DO j = 1, pop_c                       ! Read all old parameters
            READ (iwr, *     ,iostat=iostatus) &
              &   ii, parent_val (j,k)
            IF(IS_IOSTAT_END(iostatus)) THEN 
               isuccess = -13
               EXIT get_params_k
            ELSEIF(iostatus/=0 .OR. ISNAN(parent_val (j,k))) THEN
               ier_num = -37
               ier_typ = ER_APPL
               WRITE(error_io,'(a,a)') ' ***APPL*** ',fname(1:LEN_TRIM(fname))
               WRITE(error_io,'(a)')   ' ***APPL*** Generation,    member,     child'
               WRITE(error_io,'(a  ,3(2x,I9))') ' ***APPL***',pop_gen-1,i,j
               WRITE(error_io,'(a, (2x,G16.10E2))') ' ***APPL*** ',parent_val(j,k)
               CLOSE(iwr)
               RETURN
            ENDIF
         ENDDO get_params_k
         IF(isuccess == -13) THEN   ! Premature END of FILE assume population has increased
            DO ii = j, pop_c        ! Create new parameters
               parent_val(ii,k) = temp_val_max + temp_val_min
            ENDDO
         ENDIF
      ENDDO search_partial
      ENDIF
!                                                                       
!     --Determine best and worst member                                 
!                                                                       
      best  = parent_val (1,0) + 1.0 
      worst = parent_val (1,0) - 1.0 
               WRITE(error_io,'(10x,3(2x,I9))') pop_gen-1,i,j
      DO j = 1, pop_n 
         IF (parent_val (j,0) .lt.best) THEN 
            best     = parent_val(j,0)
            pop_best = j 
         ENDIF 
         IF (parent_val (j,0) .gt.worst) THEN 
            worst    = parent_val(j,0)
            pop_worst = j 
         ENDIF 
      ENDDO 
!
!     -- set values of children, needed for a 'fix <i>,best'
!     -- set values of pop_para, needed for a  init <i>
!
      DO i=1,pop_dimx
         DO j=1,pop_c
            child(i,j) = pop_x(i,j)
         ENDDO
      ENDDO
      pop_current = .true. 
   ENDIF init
   iostatus = 0
!                                                                       
! Read old trial value, if not yet initialized, and a stand alone program                     
!                                                                       
  is_alone: IF (lstandalone) THEN
  init_trial: IF (.not.pop_current_trial.and.pop_gen.gt.0) THEN 
      DO j = 1, pop_c 
        len_file = pop_ltrialfile 
        CALL make_file (pop_trialfile, len_file, 4, j) 
        CALL oeffne (iwr, pop_trialfile, stat) 
        READ (iwr, * ,END=20,ERR=20,iostat=iostatus) 
        READ (iwr, * ,END=20,ERR=20,iostat=iostatus) 
        READ (iwr, * ,END=20,ERR=20,iostat=iostatus) 
        READ (iwr, * ,END=20,ERR=20,iostat=iostatus) 
        READ (iwr, * ,END=20,ERR=20,iostat=iostatus) 
        DO i = 1, pop_dimx 
           READ (iwr, * ,END=20,ERR=20,iostat=iostatus) trial (i, j) 
        ENDDO 
      20 CONTINUE
      CLOSE (iwr) 
         IF ( iostatus /= 0) THEN
            ier_num = -12
               WRITE(error_io,'(10x,3(2x,I9))') pop_gen-1,i,j
            ier_typ = ER_APPL
            ier_msg(1) = 'Error while reading'
            WRITE (ier_msg(2),2000) j
            RETURN
         ENDIF
!
      ENDDO 
      pop_t = trial  ! (i,j) 
      pop_current_trial = .true.
   ELSE
      trial = pop_t
      pop_current_trial = .true.
     ENDIF init_trial
  ELSE  is_alone
     init_slave: IF (.not.pop_current_trial.and.pop_gen.gt.0) THEN 
        CALL create_trial
        pop_current_trial = .true.
      ELSE init_slave
         trial = pop_t
         pop_current_trial = .true.
     ENDIF init_slave
  ENDIF is_alone
!
! Read error  parent result file
!
   30 CONTINUE
   CLOSE (iwr) 
   IF ( iostatus /= 0) THEN
      ier_num = -13
      ier_typ = ER_APPL
!      ier_msg(1) = fname
      WRITE(output_io,'(a,a,a,a)') TRIM(color_err),' ***COMP*** Error while reading ',  &
         fname(1:LEN_TRIM(fname)),TRIM(color_fg)
      WRITE(output_io,'(a,a,a  )') TRIM(color_err), &
         ' ***COMP*** Check generation number in GENERATIONS and logfile',TRIM(color_fg)
      RETURN
   ENDIF
!
    950 FORMAT (A,'.',A   )
    970 FORMAT (A,'.',A,'.',I4.4   )
   2000 FORMAT ('Child No. ',i4)
!                                                                       
   END SUBROUTINE read_par_values                
!*****7**************************************************************** 
   SUBROUTINE write_parents 
!                                                                       
! Writes the logfile and summary file
!
   USE diffev_allocate_appl
   USE create_trial_mod
   USE diff_evol
   USE population
   USE lib_f90_allocate_mod
USE lib_length
   USE precision_mod
   USE variable_mod
USE sys_compiler
!
   IMPLICIT none 
!                                                                       
!                                                                       
   INTEGER, PARAMETER             :: iwr = 7
!                                                                       
   INTEGER                        :: i, j,k , kmax
   INTEGER                        :: i1, i2 
   INTEGER                        :: length
!                                                                       
!   CHARACTER (LEN=7)              :: stat  = 'append'
   CHARACTER (LEN=PREC_LSTRING)   :: line 
   CHARACTER (LEN=PREC_STRING)    :: fname
!   LOGICAL                        :: lread = .false.
!                                                                       
   REAL                           :: pave, pmin, pmax, psig 
   REAL                           :: sx, sx2, arg 
   REAL                           :: sw, wg
!
!                                                                       
   changed: IF ( pop_dimx_new ) THEN      ! Dimension has changed, patch parameter and summary file
      IF(pop_dimx.gt.MAXDIMX) THEN
         CALL alloc_population(pop_c, pop_dimx)  ! Local DIFFEV Variables
         CALL alloc_ref_para(pop_dimx)           ! Global refinement paramater variable
      ENDIF
      CALL patch_para
      pop_dimx_new = .false.
   ENDIF changed
!                                                                       
!------ write the parameters and the results for the current generation 
!                                                                       
   length = len_str(parent_results)
   kmax = 0
   IF(n_rvalue_i>1) kmax = n_rvalue_i
   i      = 0                                     ! 0 is the R-value
   partial_p: DO k=0,kmax
      IF(k==0) THEN
         WRITE (fname, 950) parent_results(1:length), pop_name(i)(1:LEN_TRIM(pop_name(i)))
      ELSE
         WRITE (fname, 970) parent_results(1:length), pop_name(i)(1:LEN_TRIM(pop_name(i))), k
      ENDIF
      CALL oeffne_append (iwr, fname, 'unknown' )
      IF (ier_num.ne.0) THEN 
         RETURN 
      ENDIF 
!                                                                       
!  write current generation as scan number                           
!                                                                       
      WRITE (iwr, 1100) pop_gen 
!                                                                       
!  write titles                                                      
!                                                                       
      WRITE (iwr, 1250) '#L Member Rvalue Rvalue '
!                                                                       
!  write the parameters of the individual members                    
!                                                                       
      DO j = 1, pop_n 
         line = ' ' 
         WRITE (iwr, 1300) j, child_val (j,k), child_val (j,k) 
      ENDDO 
      CLOSE (iwr) 
   ENDDO partial_p
!!!!!!!!!!!!!!!
!
!  Loop over all parameters pop_dimx
!
   params: DO i = 1, pop_dimx
!
      fname = ' '
      WRITE (fname, 950) parent_results(1:length), pop_name(i)(1:LEN_TRIM(pop_name(i)))
      IF(fname=='PARA0000') CYCLE params
      CALL oeffne_append (iwr, fname, 'unknown')
      IF (ier_num.ne.0) THEN 
         RETURN 
      ENDIF 
!                                                                       
!     write current generation as scan number                           
!                                                                       
      WRITE (iwr, 1100) pop_gen 
!                                                                       
!     write titles                                                      
!                                                                       
      WRITE (iwr, 1250) '#L Member Rvalue '//pop_name (i) (1:pop_lname (i) ) 
!                                                                       
!     write the parameters of the individual members                    
!                                                                       
      DO j = 1, pop_n 
         WRITE (iwr, 1300) j, child_val (j,0), child (i, j) 
      ENDDO 
      CLOSE (iwr) 
   ENDDO params
!
   CALL write_current
!                                                                       
!     Write the Summary files
!                                                                       
   length = len_str(parent_summary)
   partial_s: DO k=0,kmax
      i    = 0
      IF(k==0) THEN
         WRITE (fname, 950) parent_summary(1:length), pop_name(i)(1:LEN_TRIM(pop_name(i)))
      ELSE
         WRITE (fname, 970) parent_summary(1:length), pop_name(i)(1:LEN_TRIM(pop_name(i))), k
      ENDIF
      CALL oeffne_append (iwr, fname, 'unknown') 
      IF (ier_num.ne.0) THEN 
         RETURN 
      ENDIF 
!
      line = ' ' 
      WRITE (line (1:7), 4000) pop_gen 
      i    = 1 
      pave = 0.0 
      psig = 0.0 
      pmax = child_val (1,k) 
      pmin = child_val (1,k) 
      sx   = 0.0 
      sx2  = 0.0 
!
      DO j = 1, pop_n 
      sx = sx + child_val (j,k) 
         pmax = max (pmax, child_val (j,k) ) 
         pmin = min (pmin, child_val (j,k) ) 
      ENDDO 
!
      pave = sx / pop_n 
!
      DO j = 1, pop_n 
         sx2 = sx2 + (child_val (j,k)-pave) **2 
      ENDDO 
      arg  = sx2 / (pop_n - 1) 
      IF (arg.lt.0.0) THEN 
            psig = 0.0 
      ELSE 
            psig = sqrt (abs (arg) ) 
      ENDIF 
!
      i1 = 8 + (i - 1) * 72 
      i2 = 8 + (i - 1) * 72 + 71 
      WRITE (line (i1:i2), 4100) pave, pmin, pmax, psig
      WRITE (iwr, 4200) line (1:i2) 
      CLOSE (IWR)
   ENDDO partial_s
!
   summary: DO i = 1, pop_dimx 
      WRITE (fname, 950) parent_summary(1:length), pop_name(i)(1:LEN_TRIM(pop_name(i)))
      IF(fname=='PARA0000') CYCLE summary
      CALL oeffne_append (iwr, fname, 'unknown') 
      IF (ier_num.ne.0) return
      pave = 0.0 
      psig = 0.0 
      pmax = child (i, 1) 
      pmin = child (i, 1) 
      sx   = 0.0 
      sx2  = 0.0 
      sw   = 0.0
      DO j = 1, pop_n 
         wg   = 1./(1.+ABS(child_val(j,0) - child_val(pop_best,0)))
         sw   = sw + wg
         sx   = sx + wg*child (i, j) 
         pmax = max (pmax, child (i, j) ) 
         pmin = min (pmin, child (i, j) ) 
      ENDDO 
      pave = sx / pop_n / ABS(sw/pop_n)
      sw   = 0.0
      DO j = 1, pop_n 
         wg   = 1./(1.+(child_val(j,0) - child_val(pop_best,0))**2)
         sw   = sw + wg
         sx2  = sx2 + wg*(child (i, j)-pave) **2 
      ENDDO 
      arg  = sx2 / (pop_n - 1) / ABS(sw/pop_n)
      IF (arg.lt.0.0) THEN 
         psig = 0.0 
      ELSE 
         psig = sqrt (abs (arg) ) 
      ENDIF 
!
      line = ' ' 
      WRITE (line (1:7), 4000) pop_gen 
      i1 =  8
      i2 =  8 + 71 
      WRITE (line (i1:i2), 4100) pave, pmin, pmax, psig 
      WRITE (iwr, 4200) line (1:i2) 
      CLOSE (iwr) 
   ENDDO summary
!
!                                                                       
!  CLOSE (iwr) 
!                                                                       
!------ Copy current child into parent parameters and create new trial  
!     parameters                                                        
!                                                                       
   DO j = 1, pop_n 
      DO i = 1, pop_dimx 
         pop_x (i, j) = child (i, j) 
      ENDDO 
      parent_val (j,0:n_rvalue_i) = child_val (j,0:n_rvalue_i) 
   ENDDO 
!                                                                       
   pop_gen = pop_gen + 1 
   var_val(var_ref+0) = pop_gen ! Update global user variable
!                                                                       
   CALL create_trial 
   CALL write_genfile 
!                                                                       
!   RETURN 
!     999 CONTINUE 
!   WRITE ( * , * ) ' Error opening file' 
!                                                                       
     950 FORMAT (A,'.',a   )
     970 FORMAT (A,'.',a,'.',I4.4   )
    1100 FORMAT ('#S ',i7,' = Generation Number ') 
!    1200 FORMAT (a10) 
    1250 FORMAT (a) 
    1300 FORMAT (i5,2(2x,e18.10))
!    1310 FORMAT (2x,e18.10) 
!    2100 FORMAT (i5,' = Member     Number ',i5) 
!    3000 FORMAT (2x,e18.10) 
!    5000 FORMAT (2x,i5,2x,e18.10) 
    4000 FORMAT (i7) 
    4100 FORMAT (4(1x,e17.10)) 
    4200 FORMAT (a) 
!                                                                       
   END SUBROUTINE write_parents                  
!
!*******************************************************************************
!
   SUBROUTINE write_current
!
!   USE diffev_allocate_appl
   USE create_trial_mod
   USE diff_evol
USE lib_length
   USE population
   USE precision_mod
USE sys_compiler
!
   IMPLICIT NONE
!
   INTEGER, PARAMETER             :: iwr = 7
!                                                                       
   INTEGER                        :: i, j , k, kmax
   INTEGER                        :: length
!
   CHARACTER (LEN=PREC_STRING)    :: fname
   CHARACTER (LEN=PREC_STRING)    :: line
!                                                                       
!                                                                       
!------ write the parameters and the results for the current generation 
!                                                                       
   length = len_str(parent_current)
   IF(length  > 0 ) THEN                    ! Current file is defined
      kmax = 0
      IF(n_rvalue_i>1) kmax = n_rvalue_i
      i      = 0                                     ! 0 is the R-value
      partial_c: DO k=0,kmax
         IF(k==0) THEN
            WRITE (fname, 950) parent_current(1:length), pop_name(i)(1:LEN_TRIM(pop_name(i)))
         ELSE
            WRITE (fname, 970) parent_current(1:length), pop_name(i)(1:LEN_TRIM(pop_name(i))),k
         ENDIF
         CALL oeffne(iwr, fname, 'unknown' )
         IF (ier_num.ne.0) THEN 
            RETURN 
         ENDIF 
!     write header
         WRITE(iwr, 1000)
!                                                                       
!     write current generation as scan number                           
!                                                                       
         WRITE (iwr, 1100) pop_gen 
!                                                                       
!     write titles                                                      
!                                                                       
         WRITE (iwr, 1250) '#L Member Rvalue Rvalue '
!                                                                       
!     write the parameters of the individual members                    
!                                                                       
         DO j = 1, pop_n 
            line = ' ' 
            WRITE (iwr, 1300) j, child_val (j,k), child_val (j,k) 
         ENDDO 
         CLOSE (iwr) 
      ENDDO partial_c
!
!     Loop over all parameters pop_dimx
!
      current:DO i = 1, pop_dimx
!
         fname = ' '
         WRITE (fname, 950) parent_current(1:length),pop_name(i)(1:LEN_TRIM(pop_name(i)))
         CALL oeffne(iwr, fname, 'unknown')
         IF (ier_num.ne.0) THEN 
            RETURN 
         ENDIF 
!        write header
         WRITE(iwr, 1000)
!                                                                       
!        write current generation as scan number                           
!                                                                       
         WRITE (iwr, 1100) pop_gen 
!                                                                       
!        write titles                                                      
!                                                                       
         WRITE (iwr, 1250) '#L Member Rvalue '//pop_name (i) (1:pop_lname (i) ) 
!                                                                       
!        write the parameters of the individual members                    
!                                                                       
         DO j = 1, pop_n 
            WRITE (iwr, 1300) j, child_val (j,0), child (i, j) 
         ENDDO 
         CLOSE (iwr) 
      ENDDO current
   ENDIF
!
     950 FORMAT (A,'.',A   )
     970 FORMAT (A,'.',A,'.',I4.4   )
    1000 FORMAT ('#C Current file by DIFFEV')
    1100 FORMAT ('#S ',i7,' = Generation Number ') 
    1250 FORMAT (a) 
    1300 FORMAT (i5,2(2x,e18.10))
!
   END SUBROUTINE write_current
!
!*****7**************************************************************** 
   SUBROUTINE patch_para
!
! If the number of parameters i.e. pop_dimx has increased, 
! new Parameter(=logfile)
! and Summary files need to be written accordingly.
!
   USE population
USE lib_length
   USE precision_mod
USE sys_compiler
!
   IMPLICIT none
!
!
   INTEGER, PARAMETER    :: ird = 8
   INTEGER, PARAMETER    :: iwr = 7
!
   CHARACTER (LEN=   7)  :: stat
   CHARACTER (LEN=8*PREC_STRING)  :: line
   CHARACTER (LEN=PREC_STRING)  :: fname
   CHARACTER (LEN=PREC_STRING)  :: oname
!
   INTEGER               :: length
   INTEGER               :: i,j
   INTEGER               :: pop_dimx_old
   INTEGER               :: ios           ! I/O status variable
   LOGICAL               :: lread
   LOGICAL               :: istda
!
   LOGICAL               :: IS_IOSTAT_END
!
   stat    = 'unknown'
   lread   = .false.
!
!  loop over dimension to find old dimension
!
   length       = len_str(parent_results)
   pop_dimx_old = 0
   DO i = 1, pop_dimx
      WRITE(fname, 950) parent_results(1:length), pop_name(i)(1:LEN_TRIM(pop_name(i)))
      INQUIRE ( FILE=fname, EXIST=istda )   ! does file exist?
      IF ( .NOT. istda) EXIT                ! if not exit loop
      pop_dimx_old = i                      ! current old dimension
   ENDDO
!
!  create new Parameter files, use .0000 as template
!
   i = 0
   WRITE(fname, 950) parent_results(1:length), pop_name(i)(1:LEN_TRIM(pop_name(i)))
!
!  Loop over new files
!
   DO i = pop_dimx_old+1, pop_dimx
     WRITE(oname, 950) parent_results(1:length), pop_name(i)(1:LEN_TRIM(pop_name(i)))
     CALL OEFFNE(ird,fname,stat)              ! open Rvalue as template
     CALL OEFFNE(iwr,oname,stat)              ! open new Parameter file
     input: DO                                      ! loop over all lines
        READ ( ird, 1000, IOSTAT=ios) line
        IF ( IS_IOSTAT_END(ios)) EXIT input         ! EOF finish reading
        IF ( line(1:2)=='#C') THEN
           WRITE(IWR,1000) line(1:len_str(line))
        ELSEIF ( line(1:2)=='#S') THEN
           WRITE(IWR,1000) line(1:len_str(line))
        ELSEIF ( line(1:2)=='#L') THEN
           WRITE(IWR,1000) line(1:16)//' '//pop_name(i)(1:pop_lname (i) )
        ELSE
           WRITE(IWR,1000) line(1:25)//'    0.0000000000E+00'
        ENDIF
     ENDDO input                                    ! end loop over all lines
     CLOSE ( IRD )
     CLOSE ( IWR )
   ENDDO
!
!  create new LASTFILE  files, use .0000 as template
!
   IF(parent_current/= ' ') THEN    ! A file name exists, update
   i = 0
   length = LEN_TRIM(parent_current)
   WRITE(fname, 950) parent_current(1:length), pop_name(i)(1:LEN_TRIM(pop_name(i)))
!
!  Loop over new files
!
   DO i = pop_dimx_old+1, pop_dimx
     WRITE(oname, 950) parent_results(1:length), pop_name(i)(1:LEN_TRIM(pop_name(i)))
     CALL OEFFNE(ird,fname,stat)                    ! open Rvalue as template
     CALL OEFFNE(iwr,oname,stat)                    ! open new Parameter file
     input_c: DO                                      ! loop over all lines
        READ ( ird, 1000, IOSTAT=ios) line
        IF ( IS_IOSTAT_END(ios)) EXIT input_c       ! EOF finish reading
        IF ( line(1:2)=='#C') THEN
           WRITE(IWR,1000) line(1:len_str(line))
        ELSEIF ( line(1:2)=='#S') THEN
           WRITE(IWR,1000) line(1:len_str(line))
        ELSEIF ( line(1:2)=='#L') THEN
           WRITE(IWR,1000) line(1:16)//' '//pop_name(i)(1:pop_lname (i) )
        ELSE
           WRITE(IWR,1000) line(1:25)//'    0.0000000000E+00'
        ENDIF
     ENDDO input_c                                  ! end loop over all lines
     CLOSE ( IRD )
     CLOSE ( IWR )
   ENDDO
   ENDIF
!
!  create new Summary files, use .0000 as template
!
   length       = len_str(parent_summary)
   i = 0
   WRITE(fname, 950) parent_summary(1:length), pop_name(i)(1:LEN_TRIM(pop_name(i)))
!
!  Loop over new files
!
   DO i = pop_dimx_old+1, pop_dimx
      WRITE(oname, 950) parent_summary(1:length), pop_name(i)(1:LEN_TRIM(pop_name(i)))
      CALL OEFFNE(ird,fname,stat)                    ! open Rvalue as template
      IF(ier_num /= 0) THEN
         RETURN
      ENDIF
      CALL OEFFNE(iwr,oname,stat)                    ! open new Parameter file
      READ ( ird, 1000, IOSTAT=ios) line             ! read the three header lines
      READ ( ird, 1000, IOSTAT=ios) line
      READ ( ird, 1000, IOSTAT=ios) line
      IF ( IS_IOSTAT_END(ios)) THEN                  ! EOF finish reading
         ier_num = -13
         ier_typ = ER_APPL
         ier_msg(1) = 'Error reading the header lines of summary'
         ier_msg(2) = 'Summary files for R-values '
         RETURN
      ENDIF
      WRITE ( IWR, 1100) i                           ! Write Header line 1
      WRITE ( IWR, 1000) '#S 1'                      ! Write Header line 2
      line = '#L GEN '                              ! Prepare Header line 3
      WRITE (line ( 8:), 1200) pop_name(i)(1:LEN_TRIM(pop_name(i))), &
                                 pop_name(i)(1:LEN_TRIM(pop_name(i))), &
                                 pop_name(i)(1:LEN_TRIM(pop_name(i))), &
                                 pop_name(i)(1:LEN_TRIM(pop_name(i)))
!     DO j =  9, 13
!        IF (line (j:j) .eq.' ') line (j:j) = '_'
!     ENDDO
!     DO j = 27, 31
!        IF (line (j:j) .eq.' ') line (j:j) = '_'
!     ENDDO
!     DO j = 45, 49
!        IF (line (j:j) .eq.' ') line (j:j) = '_'
!     ENDDO
!     DO j = 63, 67
!        IF (line (j:j) .eq.' ') line (j:j) = '_'
!     ENDDO
      WRITE (iwr, 1000) line(1:LEN_TRIM(line))        ! Write Header line 1
      input2: DO                                      ! loop over all lines
         READ ( ird, *   , IOSTAT=ios) j
         IF ( IS_IOSTAT_END(ios)) EXIT input2        ! EOF finish reading
         WRITE ( IWR, 1300) j,0.0, 0.0, 0.0, 0.0
      ENDDO input2                                   ! end loop over all lines
      CLOSE ( IRD )
      CLOSE ( IWR )
   ENDDO
!
 950 FORMAT ( A,'.',A   )
1000 FORMAT ( A )
1100 FORMAT ( '#C Summaryfile by DIFFEV, Parameter no. ',I4.4)
1200 FORMAT ( ' ',a,'_AVE',2x,a,'_MIN',2x,a,'_MAX ',2x,a,'_SIG ')
1300 FORMAT ( I4, 4(1x,E17.10))
!
   END SUBROUTINE patch_para
!*****7**************************************************************** 
   SUBROUTINE do_dismiss ( lb, ub)
!                                                                       
!  Sets the R-values of the worst parents to a very high value.
!  The purpose is to ensure replacement in the next generation
!
   USE population
!
   IMPLICIT none 
!
!
   INTEGER, INTENT(IN) :: lb
   INTEGER, INTENT(IN) :: ub
!
   INTEGER             :: list_index(MAXPOP)
!
   INTEGER             :: j
   REAL                :: shift
!                                                                       
   CALL do_read_values(.TRUE.)       ! If necessary read parameter values from logfile
   IF ( ier_num /=0) RETURN
!                                                                       
!  heapsort index array  on r-values                                 
!                                                                       
   CALL indexx (MAXPOP, parent_val(:,0), list_index) 
   shift = int(parent_val(list_index(pop_n),0)) + 1.0E10
!                                                                       
!     copy the pop_n best into the child variables                      
!
   DO j = lb, ub
      parent_val (list_index(j),: ) = parent_val (list_index(j),:) + shift
   ENDDO 
!
   END SUBROUTINE do_dismiss
!*****7**************************************************************** 
   SUBROUTINE do_read_values(forced)
!
!  Reads the parameter values from the log file
!                                                                       
   USE create_trial_mod
   USE population
!
   IMPLICIT none 
!
   LOGICAL, INTENT(IN) :: forced
!                                                                       
   init: IF (pop_gen.gt.0 .OR. forced) THEN 
      CALL read_genfile 
      IF ( ier_num /=0) THEN
         ier_msg(1) = 'check existence of GENERATION'
         ier_msg(2) = 'has population been properly initialized?'
         RETURN
      ENDIF
      CALL read_obj_values 
      IF ( ier_num /=0) THEN
         ier_msg(1) = 'check for errors in result files'
         ier_msg(2) = 'has a calculation of the R-values'
         ier_msg(3) = 'been performed?'
         RETURN
      ENDIF
      CALL read_par_values 
      IF ( ier_num /=0) THEN
         ier_msg(1) = 'check existence of logfile/summary'
         ier_msg(2) = 'has population been properly initialized?'
         RETURN
      ENDIF
      pop_not_first = .TRUE.                     ! This is no longer the first time
   ENDIF init
!
   END SUBROUTINE do_read_values
!*****7**************************************************************** 
END MODULE compare
