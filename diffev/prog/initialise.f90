!*****7**************************************************************** 
!                                                                       
!
MODULE initialise
!
USE errlist_mod 
!
PRIVATE
PUBLIC  :: do_initialise
PUBLIC  :: init_x
PUBLIC  :: do_fix
!
!
CONTAINS
!
   SUBROUTINE do_initialise (l_init_x)
!-                                                                      
!     These routines creates the initial two generations.               
!     The first generation is purely random, the second is derived      
!     from the first by the standard procedure for generating a         
!     trial generation                                                  
!                                                                       
!     Version : 1.0                                                     
!     Date    : 14 Nov 2000                                             
!                                                                       
!     Authors : R.B. Neder  (reinhard.neder@krist.uni-erelangen.de)      
!-                                                                      
   USE diffev_allocate_appl
   USE constraint
   USE create_trial_mod
   USE population
   USE random_mod
!
   IMPLICIT none 
!                                                                       
   LOGICAL, INTENT(IN)  :: l_init_x     ! Shall parameter be initialised?
!                                                                       
   INTEGER, PARAMETER   :: iwr = 7
!                                                                       
   CHARACTER (LEN=2048) :: line 
   CHARACTER (LEN=1024) :: fname 
   INTEGER              :: i
   INTEGER              :: i1, i2 
   INTEGER              :: length
!                                                                       
   CHARACTER (LEN=7)    :: stat  = 'unknown'
!
   INTEGER, EXTERNAL    :: len_str
!                                                                       
!
   IF ( pop_c > MAXPOP) THEN
      ier_num = -9
      ier_typ = ER_APPL
      ier_msg(1) = 'The population size has been reduced !'
      ier_msg(2) = 'Most likely a deallocation not followed'
      ier_msg(3) = ' by a new definition of the population'
      RETURN
   ENDIF
!
   IF ( pop_dimx > MAXDIMX) THEN
      ier_num = -14
      ier_typ = ER_APPL
      ier_msg(1) = 'The parameter number has been reduced !'
      ier_msg(2) = 'Most likely a deallocation not followed'
      ier_msg(3) = ' by a new definition of the population'
      RETURN
   ENDIF
!
! should never be needed!!!   CALL diffev_alloc_population(pop_c , pop_dimx)
   IF(l_init_x) THEN 
      CALL init_x ( 1, pop_dimx)   ! Initialise parameters
      IF ( ier_num /= 0) RETURN
!                                                                 
      pop_current = .true. 
      pop_current_trial = .true. 
   ENDIF
!                                                                 
   CALL write_genfile 
!                                                                 
!    Initialise Logfiles
!                                                                 
   length = len_str(parent_results)
   DO i = 0, pop_dimx
      WRITE (fname, 950) parent_results(1:length), pop_name(i)(1:LEN_TRIM(pop_name(i)))
      CALL oeffne (iwr, fname, stat) 
      IF (ier_num.ne.0) return 
      WRITE (iwr, 1000) 
      CLOSE (iwr) 
   ENDDO
!
   IF(n_rvalue_i > 1 ) THEN    ! partial R-values are requested
      DO i=1,n_rvalue_i
         WRITE (fname, 970) parent_results(1:length), pop_name(0)(1:LEN_TRIM(pop_name(0))), i
         CALL oeffne (iwr, fname, stat)
         IF (ier_num.ne.0) return
         WRITE (iwr, 1000)
         CLOSE (iwr)
      ENDDO
   ENDIF
!                                                                 
!    Initialise Summary Files                                     
!                                                                 
   length = len_str(parent_summary)
!
   i    = 0
   WRITE (fname, 950) parent_summary(1:length), pop_name(i)(1:LEN_TRIM(pop_name(i)))
   CALL oeffne (iwr, fname, stat) 
   IF (ier_num.ne.0) return 
   WRITE (iwr, 2000) i
   WRITE (iwr, 2100) 
   line = '#L GEN '
   i1 =  8
   i2 = 71
   WRITE (line (i1:i2), 2150) 'RAVE ','RMIN ','RMAX ','RSIG ' 
   WRITE (iwr, 3000) line (1:LEN_TRIM(line)) 
   CLOSE (iwr) 
!
   IF(n_rvalue_i > 1 ) THEN    ! partial R-values are requested
      DO i=1,n_rvalue_i
         WRITE (fname, 970) parent_summary(1:length), pop_name(0)(1:LEN_TRIM(pop_name(0))), i
         CALL oeffne (iwr, fname, stat) 
         IF (ier_num.ne.0) return 
         WRITE (iwr, 2050) i
         WRITE (iwr, 2100) 
         line = '#L GEN '
         i1 =  8
         i2 = 71
         WRITE (line (i1:i2), 2150) 'RAVE ','RMIN ','RMAX ','RSIG ' 
         WRITE (iwr, 3000) line (1:LEN_TRIM(line)) 
         CLOSE (iwr) 
      ENDDO
   ENDIF
!
   DO i = 1, pop_dimx 
      WRITE (fname, 950) parent_summary(1:length), pop_name(i)(1:LEN_TRIM(pop_name(i)))
      CALL oeffne (iwr, fname, stat) 
      IF (ier_num.ne.0) return 
      WRITE (iwr, 2000) i
      WRITE (iwr, 2100) 
      line = '#L GEN '
      i1 =  8 
      i2 = 71
      WRITE (line (i1:), 2200) pop_name(i)(1:LEN_TRIM(pop_name(i))), &
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
      WRITE (iwr, 3000) line (1:LEN_TRIM(line)) 
      CLOSE (iwr) 
   ENDDO 
!
!  Create file with parameter names
!
   length = len_str(parent_results)
   fname  = parent_results(1:length) // '.name'
   CALL oeffne (iwr, fname, stat)
   IF (ier_num.ne.0) return
   WRITE(IWR, 2240) 'Member','XMIN','XMAX','SMIN','SMAX'
   WRITE(IWR, 2250) 'Rvalue', REAL(pop_gen), REAL(pop_n), &
                    REAL(pop_c), REAL(pop_dimx)
   DO i = 1, pop_dimx
      WRITE(IWR, 2250) pop_name(i)(1:len_str(pop_name(i))), &
                       pop_xmin(i), pop_xmax(i),            &
                       pop_smin(i), pop_smax(i)
   ENDDO
   CLOSE (IWR)
!                                                                       
     950 FORMAT (A,'.',A   )
     970 FORMAT (A,'.',A,'.',I4.4   )
    1000 FORMAT ('#C Logfile by DIFFEV') 
    2000 FORMAT ('#C Summaryfile by DIFFEV, Parameter no. ',i4.4) 
    2050 FORMAT ('#C Summaryfile by DIFFEV, Partial R-value no. ',i4.4) 
    2100 FORMAT ('#S 1') 
    2150 FORMAT (' ',a5,      13x,a5,      13x,a5,      13x,a5       ) 
    2200 FORMAT (' ',a,'_AVE',2x,a,'_MIN',2x,a,'_MAX',2x,a,'_SIG') 
    2240 FORMAT (a16,4(2x,a18)) 
    2250 FORMAT (a16,4(2x,E18.10)) 
    3000 FORMAT (a)
!                                                                       
   END SUBROUTINE do_initialise                     
!*****7**************************************************************** 
   SUBROUTINE init_x (lb, ub)
!
   USE diffev_config
   USE population
   USE constraint
   USE create_trial_mod
!  USE do_if_mod
   USE do_execute_mod
   USE random_mod
!
   IMPLICIT none
!
   INTEGER, INTENT(IN) :: lb,ub
!                                                                       
!
   INTEGER, PARAMETER   :: iwr = 7
!
   CHARACTER (LEN=2048) :: line 
   CHARACTER (LEN=1024) :: fname 
   INTEGER              :: n_tried 
   INTEGER              :: length 
   INTEGER              :: i,j,l
   LOGICAL              :: l_ok
   REAL                 :: w 
!                                                                       
   CHARACTER (LEN=7)    :: stat  = 'unknown'
!
   INTEGER, EXTERNAL    :: len_str
!
   REAL   , EXTERNAL    :: ran1 
!
   children: DO j = 1, pop_c 
      n_tried = 0 
      l_ok = .false. 
      attempt: DO while (.not.l_ok.and.n_tried.lt.MAX_CONSTR_TRIAL) 
         n_tried = n_tried+1 
         w = 0.00 
         DO i=1,pop_dimx ! set default values
            pop_para(i) = pop_t (i,j)
         ENDDO
         params: DO i = lb, ub
            IF(pop_smin(i) < pop_xmin(i) .or. pop_xmax(i) < pop_smax(i)) THEN
               ier_num = -25
               ier_typ = ER_APPL
               ier_msg(1) = 'Conflicting values for boundaries'
               ier_msg(2) = 'Check values for pop_xmin/max; pop_smin/smax'
               write(ier_msg(3),'(a,i5)') 'Parameter ',i
               RETURN
            ENDIF
            IF(.NOT. pop_refine(i)) THEN  ! parameter is fixed, check pop_xmin/max
               IF(pop_xmin(i) /= pop_xmax(i) .OR. &
                  MINVAL(pop_x(i,:)) /= MAXVAL(pop_x(i,:))) THEN
                  ier_num = -28
                  ier_typ = ER_APPL
                  write(ier_msg(1),'(a,a,a)') 'Parameter: ',pop_name(i),' is fixed but'
                  ier_msg(2) = 'but limits pop_xmin/pop_xmax are not identical'
                  ier_msg(3) = 'or trial parameters are not identical'
                  RETURN
               ENDIF
            ENDIF
            IF (pop_type (i) .eq.POP_INTEGER) THEN 
               pop_t (i, j) = nint (pop_smin (i) + ran1 (idum) * (pop_smax (i) - pop_smin (i) ) )
            ELSE 
               pop_t (i, j) = pop_smin (i) + ran1 (idum) * (pop_smax (i) - pop_smin (i) )
            ENDIF 
            pop_para (i) = pop_t (i, j) 
         ENDDO  params
         l_ok = .true. 
         DO l = 1, constr_number 
            line = ' ' 
            line = '('//constr_line (l) (1:constr_length (l) ) //')' 
            length = constr_length (l) + 2 
            l_ok = l_ok.and.if_test (line, length) 
         ENDDO
      ENDDO  attempt
      IF (.not.l_ok) then 
         ier_num = - 8 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
   ENDDO children
!                                                                 
   DO j = 1, pop_c 
      DO i = lb,ub 
         IF( pop_gen == 0 ) THEN
            pop_x(i, j) = pop_t(i, j) 
            child(i, j) = pop_t(i,j)
         ENDIF
      ENDDO 
      CALL write_trial (j) 
   ENDDO 
!
   DO i = lb,ub 
      pop_refine(i) = .true.
      IF(pop_sig_ad(i) == 0.0) THEN
         pop_sigma (i) = (pop_smax(i)-pop_smin(i))*0.2
         pop_sig_ad(i) = 0.2
      ELSE
         pop_sigma (i) = (pop_smax(i)-pop_smin(i))*pop_sig_ad(i)
      ENDIF
   ENDDO
!
!  Create file with parameter names
!
   length = len_str(parent_results)
   fname  = parent_results(1:length) // '.name'
   CALL oeffne (iwr, fname, stat)
   IF (ier_num.ne.0) return
   WRITE(IWR, 2250) 'Member'
   WRITE(IWR, 2250) 'Rvalue'
   DO i = 1, pop_dimx
      WRITE(IWR, 2250) pop_name(i)(1:len_str(pop_name(i)))
   ENDDO
   CLOSE (IWR)
!
    2250 FORMAT (a) 
!
   END SUBROUTINE init_x
!*****7**************************************************************** 
   SUBROUTINE do_fix ( lb, value )
!
!  This subroutine fixes the values of parameter <lb> for all children 
!  and trials to <value>. Refinement of this parameter is switched off.
!
   USE create_trial_mod
   USE population
!
   IMPLICIT NONE
!
!
   INTEGER, INTENT(IN)  :: lb
   REAL   , INTENT(IN)  :: value
!
   INTEGER              :: j
!
   IF ( pop_xmin(lb) <= value .and. value <= pop_xmax(lb) ) THEN
      IF(value == child(lb, pop_best)) THEN       ! Fix to value of current best member
      fix_best: DO j=1, pop_n
         IF(j == pop_best) CYCLE fix_best         ! Keep R-value(s) of best
         parent_val(j,0:n_rvalue_i) = 1.E12      ! Invalidate al other parents
      ENDDO fix_best
      ELSE
         parent_val(:,0:n_rvalue_i) = 1.E12      ! Invalidate al other parents
      ENDIF
      DO j=1,pop_n
         child  (lb,j) = value
         pop_x  (lb,j) = value
         pop_t  (lb,j) = value
         CALL write_trial (j)
      ENDDO
      pop_xmin(lb) = value
      pop_xmax(lb) = value
      pop_smin(lb) = value
      pop_smax(lb) = value
!
      pop_refine(lb) = .false.
      pop_sigma (lb) = 0.0
      pop_lsig  (lb) = 0.0
   ELSE
      ier_num = -18
      ier_typ = ER_APPL
      ier_msg(1) = 'Attempting to fix a parameter to a value'
      ier_msg(2) = 'that is outside the limits pop_xmin/max'
      ier_msg(3) = 'Adjust pop_xmin/pop_xmax.'
   ENDIF
!
   END SUBROUTINE do_fix
!
END MODULE initialise
