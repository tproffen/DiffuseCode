!module diffev_mache_kdo_mod
!
!contains
!
!*****7*****************************************************************
!                                                                       
SUBROUTINE diffev_mache_kdo (line, lend, length) 
!+                                                                      
!     This is the main routine for command interpretation, each         
!     command is identified here and the corresponding subroutine       
!     executed. A leading @ indicates a macro.                          
!-                                                                      
USE constraint
USE charact_mod 
USE add_param_mod
USE diff_evol
USE diffev_mpi_mod
USE diffev_reset
USE population
USE diffev_allocate_appl
USE create_trial_mod
USE compare
USE initialise
USE run_mpi_mod
USE diffev_show_mod
USE diffev_refine
USE diffev_random
USE diffev_release_mod
USE diffev_set_gen_mod
!
USE ber_params_mod
USE blanks_mod
USE build_name_mod
use calc_expr_mod
USE define_variable_mod
USE errlist_mod 
use exit_para_mod
USE gen_mpi_mod
USE get_params_mod
USE kdo_all_mod
USE learn_mod 
USE lib_errlist_func
USE lib_macro_func
USE macro_mod
USE prompt_mod
USE set_sub_generic_mod
USE precision_mod
USE str_comp_mod
USE take_param_mod
USE variable_mod
IMPLICIT none 
!                                                                       
!                                                                       
LOGICAL, PARAMETER :: IS_DIFFEV = .TRUE.
INTEGER, PARAMETER :: maxw = 20
!                                                                       
CHARACTER (LEN= *  ), INTENT(INOUT) :: line 
LOGICAL             , INTENT(  OUT) :: lend 
INTEGER             , INTENT(INOUT) :: length 
!
CHARACTER (LEN=PREC_STRING)                  :: zeile   = ' '
CHARACTER (LEN=PREC_STRING), DIMENSION(MAXW) :: cpara   = ' '
CHARACTER (LEN=  11)                  :: befehl  = ' '
CHARACTER (LEN=PREC_STRING)                  :: string  = ' '
INTEGER                               :: indxb, indxg, lcomm, lbef, indxt 
INTEGER                               :: i, j, k, ii , nb
INTEGER                               :: n_pop  ! dummy for allocation
INTEGER                               :: lb,ub
INTEGER                               :: kid, indiv, nindiv
INTEGER                               :: ianz 
INTEGER                               :: iianz 
INTEGER                               :: str_length
INTEGER             , DIMENSION(MAXW) :: lpara = 0
!INTEGER, SAVE                         :: lastgen = -1
LOGICAL                               :: back_new
LOGICAL                               :: lexist
LOGICAL                               :: lbest
LOGICAL                               :: l_init_x = .true.
!                                                                       
REAL(KIND=PREC_DP)  , DIMENSION(MAXW) :: werte = 0.0
REAL(KIND=PREC_DP)                    :: value
!                                                                       
INTEGER, PARAMETER :: NOPTIONAL = 4
INTEGER, PARAMETER :: O_PARTIAL = 1
INTEGER, PARAMETER :: O_REPEAT  = 2
INTEGER, PARAMETER :: O_LOGFILE = 3
INTEGER, PARAMETER :: O_COMPUTE = 4
CHARACTER(LEN=   7), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=PREC_STRING), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent  !opt. para present
REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 2 ! Number of values to calculate 
!
DATA oname  / 'partial', 'repeat' , 'logfile', 'compute'  /
DATA loname /  7       ,  6       ,  7       ,  7  /
opara  =  (/ '0.000000', '1.000000', 'none    ', 'parallel' /)   ! Always provide fresh default values
lopara =  (/  8        ,  8        ,  8        ,  8         /)
owerte =  (/  0.0      ,  1.0      ,  0.0      ,  1.0       /)
!
!                                                                       
CALL no_error 
!                                                                 
!-------If a commentary return immediately                        
!                                                                 
IF (line (1:1) .EQ.' ' .or. line (1:1) .eq.'#'.or. &
    line (1:1) .eq.'!' .or. length.eq.0           ) RETURN                                     
!                                                                 
!     Only the first 5 characters are significant. The command    
!     consists of the four nonblank characters                    
!                                                                 
befehl = '    ' 
indxt  = INDEX (line, tab)       ! find a tabulator
IF(indxt==0) indxt = length + 1
indxb  = index (line, ' ') 
IF(indxb==0) indxb = length + 1
indxb  = MIN(indxb,indxt)
lbef   = min (indxb - 1, 9) 
befehl = line (1:lbef) 
!                                                                 
!------ command parameters start at the first character following 
!     the blank                                                   
!                                                                 
zeile = ' ' 
lcomm = 0 
IF (indxb + 1.le.length) THEN 
   zeile = line (indxb + 1:length) 
   lcomm = length - indxb 
   call rem_leading_bl(zeile, lcomm)
ENDIF 
!                                                                 
!-------Suche nach einem "="                                      
!                                                                 
indxg = index (line, '=') 
IF (indxg.ne.0.and.                                              &
    &    .not. (str_comp (befehl, 'echo',    2, lbef, 4) ) .and.   &
    &    .not. (str_comp (befehl, 'system',  2, lbef, 6) ) .and.   &
    &    .not. (str_comp (befehl, 'fput',    2, lbef, 4) ) .and.   &
    &    .not. (str_comp (befehl, 'help',    2, lbef, 4) .or.      &
    &     str_comp (befehl, '?   ',  2, lbef, 4) )       .AND.   &
          INDEX(line,'==') == 0                               ) THEN      
!    &    .not. (str_comp (befehl, 'socket',2, lbef, 5) ) .and.   &
!                                                                 
!-------Zuweisung eines Funktionswertes                           
!                                                                 
   CALL do_math (line, indxg, length) 
ELSE 
!
   cpara(:) = ' '
   lpara(:) = 0
   werte(:) = 0.0
!                                                                 
!     --execute a macro file                                      
!                                                                 
   IF (befehl (1:1) .eq.'@') THEN 
      IF (length.ge.2) THEN 
          line = line(2:length)
          length = length -1
          CALL file_kdo(line, length)
!         CALL file_kdo (line (2:length), length - 1) 
      ELSE 
         ier_num = - 13 
         ier_typ = ER_MAC 
      ENDIF 
!                                                                 
!-------Terminate DIFFEV 'exit'                                   
!                                                                 
   ELSEIF (str_comp (befehl, 'exit', 2, lbef, 4) ) THEN 
      LEND = .TRUE. 
!                                                                 
!     ----Start of DIFFEV special commands                        
!                                                                 

!     -- Allocate array sizes
!
   ELSEIF (str_comp (befehl, 'allocate', 3, lbef,  8) ) THEN
      CALL diffev_do_allocate_appl (zeile, lcomm)
!
!     -- Deallocate array sizes
!
   ELSEIF (str_comp (befehl, 'deallocate', 3, lbef, 10) ) THEN
      CALL diffev_do_deallocate_appl (zeile, lcomm)
!                                                                 
!     -- define adaptation of sigmas 'adapt'                      
!                                                                 
   ELSEIF (str_comp (befehl, 'adapt', 3, lbef, 5) ) THEN 
      CALL get_params (zeile, ianz, cpara, lpara, maxw, length) 
      IF (ier_num.eq.0) THEN 
         IF (ianz.eq.2.or.ianz.eq.3) THEN 
            IF (str_comp (cpara(1),'sigma',3,lpara (1) , 5) ) THEN
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               iianz = 1 
               CALL ber_params (iianz, cpara, lpara, werte, maxw) 
               ii = nint (werte (1) ) 
               IF(ii==0 .or. ii>pop_dimx .or. ii>MAXDIMX) THEN
                 ier_num = -14
                 ier_typ = ER_APPL
                 RETURN
               ENDIF
               IF (ianz.eq.1.or.str_comp (cpara (2) , 'yes', 1, lpara (1) , 3) ) THEN
                  pop_ad_sigma (ii) = .true. 
               ELSEIF (ianz.eq.2.and.str_comp (cpara (2) , 'no', 2, lpara (2) , 2) ) THEN
                  pop_ad_sigma (ii) = .false. 
               ELSE 
                  CALL del_params (1, ianz, cpara, lpara, maxw) 
                  CALL ber_params (ianz, cpara, lpara, werte,maxw)
                  IF (ier_num.eq.0) THEN 
                     pop_ad_sigma (ii) = .true. 
                     pop_sig_ad (ii) = werte (1) 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ENDIF 
            ELSEIF (str_comp (cpara(1),'lsig', 3, lpara(1),4)) THEN
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               iianz = 1 
               CALL ber_params (iianz, cpara, lpara, werte, maxw) 
               ii = nint (werte (1) ) 
               IF(ii==0 .or. ii>pop_dimx .or. ii>MAXDIMX) THEN
                 ier_num = -14
                 ier_typ = ER_APPL
                 RETURN
               ENDIF
               IF (ianz.eq.1.or.str_comp (cpara (2) , 'yes', 1, lpara (1) , 3) ) THEN
                  pop_ad_lsigma (ii) = .true. 
               ELSEIF (ianz.eq.2.and.str_comp (cpara (2) , 'no', 2, lpara (2) , 2) ) THEN
                  pop_ad_lsigma (ii) = .false. 
               ELSE 
                  CALL del_params (1, ianz, cpara, lpara, maxw) 
                  CALL ber_params (ianz, cpara, lpara, werte,maxw)
                  IF (ier_num.eq.0) THEN 
                     pop_ad_lsigma (ii) = .true. 
                     pop_lsig_ad (ii) = werte (1) 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
      ENDIF 
!                                                                 
!     -- set backup option 
!                                                                 
   ELSEIF (str_comp (befehl, 'backup', 3, lbef, 6) ) THEN 
      CALL get_params (zeile, ianz, cpara, lpara, maxw, length) 
      IF (ier_num.eq.0) THEN 
         IF(cpara(1)=='NONE') THEN
            pop_backup = .false.
            pop_back_number = 0
         ELSE
            IF(pop_back_number==MAXBACK) THEN
               CALL alloc_backup(n_pop)
            ENDIF
            back_new = .TRUE.
            loop_back: DO i=1,pop_back_number
               IF(ianz==3) THEN
                  IF((pop_back_fil(i)                                == cpara(1) .OR.    &
                      pop_back_fil(i)(1:LEN_TRIM(pop_back_fil(i))-1) == cpara(1)) .AND.  &
                     (pop_back_ext(i)                                == cpara(2) .OR.    &
                      pop_back_ext(i)(2:LEN_TRIM(pop_back_ext(i))  ) == cpara(2)) .AND.  &
                     (pop_back_trg(i)                                == cpara(3) .OR.    &
                      pop_back_trg(i)(1:LEN_TRIM(pop_back_trg(i))-1) == cpara(3))       ) THEN
                     back_new = .FALSE.
                     EXIT loop_back
                  ENDIF
               ELSE IF(ianz==2) THEN
                  IF((pop_back_fil(i)                                == cpara(1) .OR.    &
                      pop_back_fil(i)(1:LEN_TRIM(pop_back_fil(i))-1) == cpara(1)) .AND.  &
                     (pop_back_trg(i)                                == cpara(2) .OR.    &
                      pop_back_trg(i)(1:LEN_TRIM(pop_back_trg(i))-1) == cpara(2))       ) THEN
                     back_new = .FALSE.
                     EXIT loop_back
                  ENDIF
               ENDIF
            ENDDO loop_back
            IF(back_new) THEN
            pop_back_number = pop_back_number + 1
            nb              = pop_back_number
            IF(cpara(ianz)(lpara(ianz):lpara(ianz))=='.') THEN
               pop_back_trg  (nb) = cpara(ianz)(1:lpara(ianz))
               pop_back_trg_l(nb) = lpara(ianz)
            ELSE
               pop_back_trg  (nb) = cpara(ianz)(1:lpara(ianz)) // '.'
               pop_back_trg_l(nb) = lpara(ianz) + 1
            ENDIF
            ianz = ianz - 1
            IF(ianz==2) THEN
               IF(cpara(ianz)(1:1)=='.') THEN
                  pop_back_ext  (nb) = cpara(ianz)(1:lpara(ianz))
                  pop_back_ext_l(nb) = lpara(ianz)
               ELSE
                  pop_back_ext  (nb) = '.' // cpara(ianz)(1:lpara(ianz))
                  pop_back_ext_l(nb) = lpara(ianz)+1
               ENDIF
            ELSE
               pop_back_ext   = ' '
               pop_back_ext_l = 1
            ENDIF
            IF(cpara(1)(lpara(1):lpara(1))=='.') THEN
               pop_back_fil  (nb) = cpara(1)(1:lpara(1))
               pop_back_fil_l(nb) = lpara(1)
            ELSE
               pop_back_fil  (nb) = cpara(1)(1:lpara(1)) // '.'
               pop_back_fil_l(nb) = lpara(1) + 1
            ENDIF 
            pop_backup = .true.
            ENDIF    ! if(back_new)
         ENDIF 
      ENDIF 
!                                                                 
!     -- compare the trial results to the last generation         
!                                                                 
   ELSEIF (str_comp (befehl, 'comp', 3, lbef, 4) ) THEN 
      CALL get_params (zeile, ianz, cpara, lpara, maxw, length) 
      IF (ier_num.eq.0) THEN 
         IF (ianz.eq.0) THEN
            pop_result_file_rd = .true.
         ELSEIF (ianz.eq.1) THEN
            IF(str_comp (cpara(1), 'silent',6,lpara(1),6)) THEN
               IF(lstandalone) THEN 
                  ier_num = -27
                  ier_typ = ER_APPL
                  RETURN 
               ELSE
                  pop_result_file_rd = .false.
               ENDIF
            ENDIF
         ELSE
            ier_num = -6
            ier_typ = ER_COMM
            RETURN 
         ENDIF
      ELSE
         ier_num = -6
         ier_typ = ER_COMM
         RETURN 
      ENDIF
      CALL do_compare 
!
!     Turn random state log ON
!
      CALL diffev_random_on
!                                                                 
!     -- Define a hard constraint 'constraint'                    
!                                                                 
   ELSEIF (str_comp (befehl, 'constraint', 3, lbef, 10) ) THEN 
      IF(constr_number >= MAX_CONSTR) THEN
         CALL alloc_constraint ( constr_number + 1)
         IF(ier_num < 0) THEN
            RETURN
         ENDIF
      ENDIF
      IF (constr_number.lt.MAX_CONSTR) THEN 
         constr_number = constr_number + 1 
         constr_line (constr_number) = zeile (1:length) 
         constr_length (constr_number) = length 
      ELSE 
         ier_num = - 7 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
!                                                                 
!     -- initialise the sequence                                  
!                                                                 
   ELSEIF (str_comp (befehl, 'dismiss', 3, lbef, 7) ) THEN 
      IF (pop_n.gt.3) THEN
         CALL get_params (zeile, ianz, cpara, lpara, maxw, length) 
         IF (ier_num.eq.0) THEN 
            IF (ianz.eq.0 .or.  &
                str_comp (cpara(1), 'all',3, lpara(1), 3)) THEN 
               lb = 1
               ub = pop_n
               CALL do_dismiss  (lb, ub)
            ELSE
               CALL ber_params (ianz, cpara, lpara, werte, maxw)
               IF (ier_num.eq.0) THEN
                  IF (ianz.eq.1) THEN 
                     lb = pop_n + 1 - nint(werte(1))
                     ub = pop_n
                  ELSE
                     ier_num = -6
                     ier_typ = ER_COMM
                     RETURN 
                  ENDIF 
                  IF ( 0<lb .and. lb<=ub .and. ub<=pop_n) THEN
                     CALL do_dismiss ( lb,ub)
                  ELSE
                     ier_num = -6
                     ier_typ = ER_COMM
                     RETURN 
                  ENDIF 
               ENDIF 
            ENDIF 
         ENDIF 
      ELSE 
         ier_num = - 3 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
!                                                                 
!     -- set the donor type                                       
!                                                                 
   ELSEIF (str_comp (befehl, 'donor', 4, lbef, 4) ) THEN 
      CALL get_params (zeile, ianz, cpara, lpara, maxw, length) 
      IF (ier_num.eq.0) THEN 
         IF (ianz.eq.1) THEN 
            IF (str_comp (cpara(1),'best',3,lpara(1),4) ) THEN
               diff_donor_mode = ADD_TO_BEST 
            ELSEIF (str_comp(cpara(1),'random',3,lpara(1),6)) THEN
               diff_donor_mode = ADD_TO_RANDOM 
            ENDIF 
         ENDIF 
      ENDIF 
!                                                                 
!     -- fix a parameter at a value                               
!                                                                 
   ELSEIF (str_comp (befehl, 'fix', 3, lbef, 3) ) THEN 
      INQUIRE(FILE='GENERATION', EXIST=lexist)
      IF(lexist) THEN                ! A GENERATION FILE EXISTS
         CALL do_read_values(.TRUE.) ! We need to read values as this can be the first command after a continue
      ENDIF
      IF (pop_n.gt.3) THEN
         IF((pop_gen==0 .AND. .NOT. pop_initialized) .OR. .NOT.lexist) THEN   ! Population was not yet initialized
            IF(pop_trialfile == ' ') pop_trial_file_wrt= .FALSE.
            CALL do_initialise (l_init_x)
            pop_initialized = .TRUE.
         ENDIF
         IF(ier_num/=0) THEN
            ier_msg(1) = 'Check the GENERATION, the parameter and'
            ier_msg(2) = 'the lastfile for conflicting generation values'
            RETURN
         ENDIF
         CALL get_params (zeile, ianz, cpara, lpara, maxw, length) 
         IF (ier_num == 0) THEN 
            IF (ianz == 1) THEN 
               cpara(2) = 'best'
               lpara(2) = 4
               ianz = 2
            ENDIF
!           
            IF (ianz == 2) THEN 
               lbest = .false.
!              Check if parameter names were provided
!              DO i = 1, ianz
               fix_t: DO k=1,pop_dimx
                     IF(cpara(1)==pop_name(k)) THEN
                        WRITE(cpara(1),'(I4)') k
                        EXIT fix_t
                     ENDIF
                  ENDDO fix_t
!              ENDDO
               IF (str_comp (cpara(2), 'best', 3, lpara(2), 4) ) THEN 
                  lbest = .true.
                  ianz  = ianz - 1
                  CALL ber_params (ianz, cpara, lpara, werte, maxw)
                  IF (ier_num.eq.0) THEN
                     lb = nint(werte(1))
                  ELSE
                     ier_num = -6
                     ier_typ = ER_COMM
                     ier_msg(1) = 'Parameter name unknown or not a proper number'
                     RETURN 
                  ENDIF 
               ELSE
                  CALL ber_params (ianz, cpara, lpara, werte, maxw)
                  IF (ier_num.eq.0) THEN
                     lb = nint(werte(1))
                  ELSE
                     ier_num = -6
                     ier_typ = ER_COMM
                     ier_msg(1) = 'Parameter name unknown or not a proper number'
                     RETURN 
                  ENDIF 
               ENDIF 
               IF ( 0<lb .and. lb<=pop_dimx) THEN
!                 CALL do_read_values(.FALSE.)   ! Read values from logfile
                  IF(lbest) THEN
                     value = child(lb,pop_best)
                  ELSE
                     value = werte(2)
                  ENDIF
                  CALL do_fix ( lb, value)
               ELSE
                  ier_num = -6
                  ier_typ = ER_COMM
                  RETURN 
               ENDIF 
            ELSE
               ier_num = -6
               ier_typ = ER_COMM
               RETURN
            ENDIF 
         ENDIF 
      ELSE 
         ier_num = - 3 
         ier_typ = ER_APPL 
         ier_msg(1) = 'To fix a parameter a proper population size'
         ier_msg(2) = 'must have been defined ==> pop_n[1], pop_c[1]'
         RETURN 
      ENDIF 
!                                                                 
!     -- initialise the sequence                                  
!                                                                 
   ELSEIF (str_comp (befehl, 'init', 3, lbef, 4) ) THEN 
      CALL get_params (zeile, ianz, cpara, lpara, maxw, length) 
      IF (ier_num.eq.0) THEN 
         IF (ianz.eq.0) THEN 
            pop_trial_file_wrt = .true.
            l_init_x = .true.
            CALL do_initialise (l_init_x)
            pop_initialized = .TRUE.
         ELSE
            pop3: IF (pop_n.gt.3.or.                                                 &
                str_comp (cpara(ianz),'logfile',3, lpara(ianz), 7).AND.        &
                ianz==1 ) THEN
               IF(str_comp (cpara(ianz),'silent',6, lpara(ianz), 6).AND. ianz==1) THEN
                  IF(lstandalone) THEN 
                     ier_num = -27
                     ier_typ = ER_APPL
                     RETURN 
                  ELSE
                     pop_trial_file_wrt = .false.
                     ianz = ianz - 1
                     l_init_x = .true.
                     CALL do_initialise (l_init_x)
                     pop_initialized = .TRUE.
                  ENDIF
               ELSEIF(str_comp (cpara(ianz),'logfile',3, lpara(ianz), 7).AND. ianz==1) THEN
                  IF(.NOT.pop_current) THEN       ! Population has not been read/calculated
                     INQUIRE(FILE='GENERATION', EXIST=lexist)
                     IF(lexist) THEN                     ! A GENERATION FILE EXISTS
                        CALL do_read_values(.TRUE.)         ! We need to read values as this can be the first command after a continue
                     ELSE
                        ier_num = -35
                        ier_typ = ER_APPL
                        RETURN
                     ENDIF
                  ENDIF
                  IF(pop_gen>0) THEN
                  l_init_x = .false.
                  CALL do_initialise (l_init_x)     ! Write empty log files
                  pop_gen = 0
                  CALL write_parents             ! Add the current scan to the parameter files
                  ELSE
                     ier_num = -20
                     ier_typ = ER_APPL
                     RETURN
                  ENDIF
               ELSE
!
!                 If last parameter is 'silent' turn trial files off, else leave current status
!
                  IF(str_comp (cpara(ianz),'silent',6, lpara(ianz), 6)) THEN
                     pop_trial_file_wrt = .false.
                     ianz = ianz -1
                  ENDIF
                  IF (ianz == 1 .OR. ianz ==2) THEN 
!               CALL read_genfile ! init,<i> works only in generations > 0
!               IF(ier_num /= 0) RETURN
                     IF (pop_gen <= 0 ) THEN
                        ier_num = -20
                        ier_typ = ER_APPL
                        ier_msg(1) = 'A single parameter can only be '
                        ier_msg(2) = 'initialized in generations > 0 '
                        ier_msg(3) = 'Run at least one compare       '
                        RETURN
                     ENDIF
!                    Check if parameter names were provided
                     DO i = 1, ianz
                        DO k=1,pop_dimx
                           IF(cpara(i)==pop_name(k)) THEN
                              WRITE(cpara(i),'(I4)') k
                           ENDIF
                        ENDDO
                     ENDDO
                     CALL ber_params (ianz, cpara, lpara, werte, maxw)
                     IF (ier_num == 0) THEN
                        IF (ianz == 1) THEN 
                           lb = NINT(werte(1))
                           ub = NINT(werte(1))
                        ELSEIF (ianz.eq.2 ) THEN
                           lb = MIN(NINT(werte(1)), NINT(werte(2)))
                           ub = MAX(NINT(werte(2)), NINT(werte(1)))
                        ELSE
                           ier_num = -6
                           ier_typ = ER_COMM
                           RETURN 
                        ENDIF 
                        IF ( 0<lb .and. lb<=ub .and. ub<=pop_dimx) THEN
                           CALL do_read_values(.FALSE.)   ! Read values from logfile
                           IF( ier_num /= 0) RETURN
                           CALL init_x ( lb,ub)  ! Initialize parameter range
                           IF( ier_num /= 0) RETURN
                           CALL write_genfile    ! Write the "GENERATION" file
                        ELSE
                           ier_num = -6
                           ier_typ = ER_COMM
                           RETURN 
                        ENDIF 
                     ENDIF 
                  ENDIF 
               ENDIF
            ELSE  pop3
               ier_num = - 3 
               ier_typ = ER_APPL 
               RETURN 
            ENDIF  pop3
         ENDIF 
!
!           Turn random state log ON
!
            CALL diffev_random_on
      ENDIF 
!                                                                 
!     -- set the logfile file                                     
!                                                                 
   ELSEIF (str_comp (befehl, 'logfile', 3, lbef, 7) ) THEN 
      CALL get_params (zeile, ianz, cpara, lpara, maxw, length) 
      IF (ier_num.eq.0) THEN 
         CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
         IF (ier_num.eq.0) THEN 
            parent_results = cpara (1)(1:lpara(1)) 
            lparent_results = lpara (1) 
         ENDIF 
      ENDIF 
!                                                                 
!     -- set the current file                                     
!                                                                 
   ELSEIF (str_comp (befehl, 'lastfile', 3, lbef, 8) ) THEN 
      CALL get_params (zeile, ianz, cpara, lpara, maxw, length) 
      IF (ier_num.eq.0) THEN 
         CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
         IF (ier_num.eq.0) THEN 
            parent_current = cpara (1)(1:lpara(1)) 
            lparent_current = lpara (1) 
         ENDIF 
      ENDIF 
!                                                                 
!     -- add a new parameter to the dimension                     
!                                                                 
   ELSEIF (str_comp (befehl, 'newparam', 3, lbef, 8) ) THEN 
      CALL add_param(zeile, length)
!                                                                 
!     -- set the name of a refinement parameter                   
!                                                                 
   ELSEIF (str_comp (befehl, 'pop_name', 3, lbef, 8) ) THEN 
      CALL get_params (zeile, ianz, cpara, lpara, maxw, length) 
      IF (ier_num.eq.0) THEN 
         iianz = 1 
         CALL ber_params (iianz, cpara, lpara, werte, maxw) 
         IF (ier_num.eq.0) THEN 
            IF (nint(werte(1)) >  0       .and.         &
                nint(werte(1)) <= MAXDIMX .and.         &
                nint(werte(1)) <= pop_dimx     ) THEN
                IF(lpara(2)>LEN(pop_name)) THEN
                   ier_num = -31
                   ier_typ = ER_APPL
                   WRITE(ier_msg(1),'(a,i2)') 'Parameter length must be <= ',LEN(pop_name)
                ELSE
                   pop_lname (nint (werte(1))) = min(lpara(2), LEN(pop_name))
                   pop_name (nint(werte(1))) = cpara(2)(1:pop_lname(nint(werte(1)))) 
                   DO j=1,var_num
                      IF(pop_name(NINT(werte(1))) == var_name(j)) THEN
                        ier_num = -29
                        ier_typ = ER_APPL
                        ier_msg(1) = 'The parameter names must be unique. Check'
                        ier_msg(2) = 'with >>''variable show'' for existing names.'
                        RETURN
                     ENDIF
                  ENDDO
                  pop_dimx_init = .TRUE.      ! The dimension has been initialized in this run
                  string     = 'real, '//pop_name(NINT(werte(1)))
                  str_length = 6+pop_lname(NINT(werte(1)))
                  CALL define_variable(string, str_length, IS_DIFFEV)
               ENDIF 
            ELSE
               ier_num = -14
               ier_typ = ER_APPL
            ENDIF 
         ENDIF 
      ENDIF 
   ELSEIF (str_comp (befehl, 'read', 4, lbef, 4) ) THEN
      pop_current = .FALSE.
      CALL do_read_values(.TRUE.)         ! Always try to read parameters as instructed by user
   ELSEIF (str_comp (befehl, 'run_mpi', 7, lbef, 7) ) THEN
      INQUIRE(FILE='GENERATION', EXIST=lexist)
      IF(lexist) THEN                     ! A GENERATION FILE EXISTS
         CALL do_read_values(.TRUE.)         ! We need to read values as this can be the first command after a continue
      ENDIF
      IF((pop_gen==0 .AND. .NOT. pop_initialized) .OR. .NOT.lexist) THEN   ! Population was not yet initialized
         IF(pop_trialfile == ' ') pop_trial_file_wrt= .FALSE.
         CALL do_initialise (l_init_x)
         pop_initialized = .TRUE.
      ENDIF
      IF(ier_num/=0) THEN
         ier_msg(1) = 'Check the GENERATION, the parameter and'
         ier_msg(2) = 'the lastfile for conflicting generation values'
         RETURN
      ENDIF
!     IF ( .not. gen_mpi_active .or. gen_mpi_numprocs < 2 ) THEN
      IF (       gen_mpi_active .and. gen_mpi_numprocs < 2 ) THEN
         ier_num =  -24
         ier_typ = ER_APPL
         ier_msg(1) = 'To run an MPI distributed application requires'
         ier_msg(2) = 'that diffev has been compiled with MPI'
         ier_msg(3) = 'and been started with mpiexec -option diffev'
         RETURN
      ENDIF
      IF ( pop_c < 4 .or. pop_n < 4 ) THEN
         ier_num = -3
         ier_typ = ER_APPL
         ier_msg(1) = 'To run an MPI distributed application requires'
         ier_msg(2) = 'a properly defined population'
         ier_msg(3) = 'Check population size, dimension etc'
         RETURN
      ENDIF
!
      call alloc_senddata(max(1,pop_dimx), max(1,n_rvalue_i))
!
      DO i=1, pop_dimx
         IF(.NOT. pop_refine(i)) THEN  ! parameter is fixed, check pop_xmin/max
            IF(pop_xmin(i) /= pop_xmax(i) .OR.         &
               MINVAL(pop_t(i,1:MIN(pop_n,pop_c))) /=  &
               MAXVAL(pop_t(i,1:MIN(pop_n,pop_c)))     ) THEN
               ier_num = -28
               ier_typ = ER_APPL
               WRITE(ier_msg(1),'(a,i4,a)') 'Parameter no.: ',i,' is fixed but'
               ier_msg(2) = 'Limits pop_xmin/pop_xmax are not identical'
               ier_msg(3) = 'or trial parameters are not identical'
               RETURN
            ENDIF
         ENDIF
      ENDDO
      CALL get_params (zeile, ianz, cpara, lpara, maxw, length) 
      IF (ier_num.eq.0) THEN 
         CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                           oname, loname, opara, lopara, lpresent, owerte)
         IF (ier_num == 0) THEN 
         run_mpi_senddata%l_get_state = l_get_random_state   ! Log random number state
         IF(cpara(3) == 'DOLOOP') THEN          ! Special signal set if MPI not active
                                                ! and 'run_mpi' within a do loop
            run_mpi_senddata%generation = pop_gen    ! Current GENERATION no
            run_mpi_senddata%member     = pop_n      ! Number of members
            run_mpi_senddata%children   = pop_c      ! Number of children
            run_mpi_senddata%parameters = pop_dimx   ! Number of parameters
            run_mpi_senddata%prog   = cpara(1)(1:lpara(1))
            run_mpi_senddata%prog_l = lpara(1)
            run_mpi_senddata%mac    = cpara(2)(1:lpara(2))
            run_mpi_senddata%mac_l  = lpara(2)
            run_mpi_senddata%n_rvalue_i = n_rvalue_i ! We expect this many Rvalues
            CALL refine_no_mpi(.true.)
         ELSE
!         IF ( ianz == 5 ) THEN
!            run_mpi_senddata%use_socket = str_comp(cpara(5), 'socket', 6, lpara(5), 6)
!            ianz = 4
!         ELSE
!            run_mpi_senddata%use_socket = .false.
!         ENDIF
         IF ( ianz >= 2 .and. ianz <= 4 ) THEN
            run_mpi_senddata%generation = pop_gen    ! Current GENERATION no
            run_mpi_senddata%member     = pop_n      ! Number of members
            run_mpi_senddata%children   = pop_c      ! Number of children
            run_mpi_senddata%parameters = pop_dimx   ! Number of parameters
            run_mpi_senddata%prog   = cpara(1)(1:100) ! Program to run
            run_mpi_senddata%prog_l = lpara(1)
            run_mpi_senddata%mac    = cpara(2)(1:100) ! Macro to run
            run_mpi_senddata%mac_l  = lpara(2)
            run_mpi_senddata%n_rvalue_i = n_rvalue_i ! We expect this many Rvalues
            IF ( ianz == 4 ) THEN
               run_mpi_senddata%out   = cpara(4)(1:100)! Target for program output 
               run_mpi_senddata%out_l = lpara(4)
            ELSE
               IF(opara(O_LOGFILE)=='none') THEN
                  run_mpi_senddata%out   = '/dev/null'   ! Default output
                  run_mpi_senddata%out_l = 9
               ELSE
                  run_mpi_senddata%out   = opara(O_LOGFILE)(1:lopara(O_LOGFILE))
                  run_mpi_senddata%out_l = lopara(O_LOGFILE)
               ENDIF
            ENDIF 
            run_mpi_senddata%repeat = .false.! repeat = .false. means no repetition
            run_mpi_senddata%nindiv = 1      ! nindiv is at least 1
            IF ( ianz >  2 ) THEN            ! Number of local repeations
               CALL del_params (2, ianz, cpara, lpara, maxw) 
               ianz = 1
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               IF (ier_num.eq.0) THEN 
                  IF ( nint(werte(1)) > 1 ) THEN
                     run_mpi_senddata%repeat = .true.             ! repeat = false means no repetition
                     IF (str_comp (opara(O_COMPUTE), 'parallel', 4, lopara(O_COMPUTE), 8) ) THEN
                        run_mpi_senddata%repeat = .TRUE.
                     ELSE
                        run_mpi_senddata%repeat = .FALSE.
                     ENDIF
                  ELSE !
                     run_mpi_senddata%repeat = .FALSE.
                  ENDIF 
                  run_mpi_senddata%nindiv = max(1,nint(werte(1))) ! nindiv is at least 1
               ENDIF 
            ELSE 
                     IF (str_comp (opara(O_COMPUTE), 'parallel', 4, lopara(O_COMPUTE), 8) ) THEN
                        run_mpi_senddata%repeat = .TRUE.
                     ELSE
                        run_mpi_senddata%repeat = .FALSE.
                     ENDIF
                  run_mpi_senddata%nindiv = max(1,nint(owerte(O_REPEAT))) ! nindiv is at least 1
            ENDIF 
            IF(gen_mpi_active .AND. pop_gen>lastgen)  THEN ! Flag errors if new generation
               IF(run_mpi_senddata%repeat) THEN            ! parallel refinement of indivs
!                 IF(NUM_NODE>pop_c+1) THEN
!                    ier_num = -32
!                    ier_typ = ER_APPL
!                    ier_msg(1) = 'For compute:parallel refinement of indivs'
!                    ier_msg(2) = 'Node number must be <= population size'
!                    ier_msg(3) = 'Node*core must be <= pop_c*REF_NINDIV+1'
!                    RETURN
!                 ENDIF
                  IF(gen_mpi_numprocs>pop_c*run_mpi_senddata%nindiv+2) THEN
                     ier_num = -33
                     ier_typ = ER_APPL
                     ier_msg(1) = 'For compute:parallel refinement of indivs'
                     ier_msg(2) = 'Node number must be <= population size+1'
                     ier_msg(3) = 'Node*core must be <= pop_c*REF_NINDIV+1'
                     RETURN
                  ENDIF
               ELSE                                       ! Serial refinement of indivs
                  IF(NUM_NODE>pop_c+1) THEN
                     ier_num = -33
                     ier_typ = ER_APPL
                     ier_msg(1) = 'For compute:serial   refinement of indivs'
                     ier_msg(2) = 'Node number must be <= population size+1'
                     ier_msg(3) = ' '
                     RETURN
                  ENDIF
               ENDIF
            ENDIF
            IF ( ier_num == 0) THEN
                  IF(.NOT.lstandalone) THEN
                     init_slave: IF (.not.pop_current_trial.and.pop_gen.gt.0) THEN
                        CALL read_par_values
                     ENDIF init_slave
                  ENDIF
               IF(gen_mpi_active) THEN   !Parallel processing with MPI
!                 run_mpi_kid_per_core = INT(pop_c/(run_mpi_max_slaves*NUM_NODE))+1
!                 IF(.NOT.ALLOCATED(kid_on_node)) THEN
!                    ALLOCATE(kid_on_node(0:pop_c))
!                    kid_on_node(:) = 0
!                 ENDIF
!                 IF(.NOT.ALLOCATED(node_has_kids)) THEN
!                    ALLOCATE(node_has_kids(1:NUM_NODE,0:run_mpi_max_slaves*run_mpi_kid_per_core,2))
!                    node_has_kids(:,:,:) = 0
!                 ENDIF
                  CALL run_mpi_master 
!                 lastgen = pop_gen     ! We finished a refinement in this generation
!                 DEALLOCATE(kid_on_node)
!                 DEALLOCATE(node_has_kids)
               ELSE    ! gen_mpi_active
!                 MPI is not active, refine with non parallel algorithm
                  CALL refine_no_mpi(.false.)
               ENDIF   ! gen_mpi_active
            ENDIF 
         ELSE
            ier_num =  -6
            ier_typ =  -ER_COMM
         ENDIF 
         ENDIF 
      ENDIF 
      ENDIF 
!
!     Turn random state log OFF, and documentation ON
!
!     IF(diffev_random_status()) 
      CALL diffev_random_write_on( &
                                    run_mpi_senddata%prog,    &
                                    run_mpi_senddata%prog_l,  &
                                    run_mpi_senddata%mac,     &
                                    run_mpi_senddata%mac_l,   &
                                    run_mpi_senddata%repeat  )
      CALL diffev_random_off
!                                                                 
!-------  Show parameters 'show'                                  
!                                                                 
   ELSEIF (str_comp (befehl, 'show', 3, lbef, 4) ) THEN 
      CALL diffev_do_show (zeile, lcomm) 
!                                                                 
!     -- set the average information file                         
!                                                                 
   ELSEIF (str_comp (befehl, 'summaryfile', 3, lbef, 11) ) THEN 
      CALL get_params (zeile, ianz, cpara, lpara, maxw, length) 
      IF (ier_num.eq.0) THEN 
         CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
         IF (ier_num.eq.0) THEN 
            parent_summary = cpara (1)(1:lpara(1)) 
            lparent_summary = lpara (1) 
         ENDIF 
      ENDIF 
!                                                                 
!     -- define the refine status of each gene                    
!                                                                 
   ELSEIF (str_comp (befehl, 'refine', 3, lbef, 6) ) THEN 
      CALL get_params (zeile, ianz, cpara, lpara, maxw, length) 
      IF (ier_num.eq.0) THEN 
         IF (ianz.ge.1) THEN 
            IF (str_comp (cpara(1),'all',1, lpara (1) , 3) ) THEN
               DO i = 1, MAXDIMX 
               pop_refine (i) = .true. 
               ENDDO 
            ELSEIF (str_comp (cpara(1),'none',1,lpara(1),3) ) THEN
               DO i = 1, MAXDIMX 
               pop_refine (i) = .false. 
               ENDDO 
            ELSE 
!              Check if parameter names were provided
               DO i = 1, ianz
                  DO k=1,pop_dimx
                     IF(cpara(i)==pop_name(k)) THEN
                        WRITE(cpara(i),'(I4)') k
                        lpara(i) = 4
                     ENDIF
                  ENDDO
               ENDDO
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               IF (ier_num.eq.0) THEN 
                  DO i = 1, ianz 
                  j = nint (werte (i) ) 
                  IF(j==0 .or. j>pop_dimx .or. j>MAXDIMX) THEN
                    ier_num = -14
                    ier_typ = ER_APPL
                    RETURN
                  ENDIF
                  IF (j.gt.0) THEN 
                     pop_refine (j) = .true. 
                  ELSEIF (j.lt.0) THEN 
                     pop_refine (abs (j) ) = .false. 
                  ENDIF 
                  ENDDO 
               ENDIF 
            ENDIF 
         ENDIF 
      ELSE
         IF(ier_num == -17) THEN  ! too many numbers in this line
            ier_msg(1) = 'The refine command can only take 20 values.'
            ier_msg(2) = 'Please use several refine command lines for'
            ier_msg(3) = 'all these parameters.'
         ENDIF
      ENDIF 
!
!     -- Release a previously fixed parameter 'release'          
!
   ELSEIF (str_comp (befehl, 'release', 3, lbef, 7) ) THEN 
      IF (pop_gen <= 0 ) THEN
         ier_num = -20
         ier_typ = ER_APPL
         ier_msg(1) = 'A single parameter can only be '
         ier_msg(2) = 'released in generations > 0 '
         ier_msg(3) = 'Run at least one compare       '
         RETURN
      ENDIF
      CALL do_release(zeile, length)
!                                                                 
!     -- set the selection mode 'selection'                       
!                                                                 
   ELSEIF (str_comp (befehl, 'selection', 2, lbef, 9) ) THEN 
      CALL get_params (zeile, ianz, cpara, lpara, maxw, length) 
      IF (ier_num.eq.0) THEN 
         IF (ianz.ge.1) THEN 
            IF (str_comp (cpara(1),'best',3,lpara(1),4)) THEN
               IF (str_comp (cpara(2),'all',3,lpara(2), 3) ) THEN
                  diff_sel_mode = SEL_BEST_ALL 
               ELSEIF (str_comp(cpara(2),'children',3,lpara(2), 8) ) THEN
                  diff_sel_mode = SEL_BEST_CHILD 
               ENDIF 
            ELSEIF (str_comp(cpara(1),'compare',3,lpara(1),7)) THEN
               diff_sel_mode = SEL_COMP 
            ENDIF 
         ENDIF 
      ENDIF 
!                                                                 
!     -- set the trial file                                       
!                                                                 
   ELSEIF (str_comp (befehl, 'trialfile', 3, lbef, 9) ) THEN 
      CALL get_params (zeile, ianz, cpara, lpara, maxw, length) 
      IF (ier_num.eq.0) THEN 
         IF(str_comp (cpara(1), 'silent',6,lpara(1),6)) THEN
            IF(lstandalone) THEN 
               ier_num = -27
               ier_typ = ER_APPL
               RETURN 
            ELSE
               pop_trial_file_wrt = .false.
               pop_trialfile = ' '
               pop_ltrialfile = 1
            ENDIF 
         ELSE
            pop_trial_file_wrt = .true.
            CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
            IF (ier_num.eq.0) THEN 
               pop_trialfile = cpara (1)(1:lpara(1)) 
               pop_ltrialfile = lpara (1) 
            ENDIF 
         ENDIF 
      ENDIF 
!                                                                 
!     -- set the parameter type                                   
!                                                                 
   ELSEIF (str_comp (befehl, 'type', 4, lbef, 4) ) THEN 
      CALL get_params (zeile, ianz, cpara, lpara, maxw, length) 
      IF (ier_num.eq.0) THEN 
         IF (ianz.eq.2) THEN 
            IF (str_comp(cpara(1),'integer', 3, lpara (1) , 7)) THEN
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               IF (ier_num.eq.0) THEN 
                  pop_type (nint (werte (1) ) ) = POP_INTEGER 
               ENDIF 
            ELSEIF (str_comp(cpara(1),'real',3,lpara(1),4) ) THEN
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               IF (ier_num.eq.0) THEN 
                  pop_type (nint (werte (1) ) ) = POP_REAL 
               ENDIF 
            ENDIF 
         ENDIF 
      ENDIF 
!                                                                 
!     -- Restart in a previous generation                         
!                                                                 
   ELSEIF (str_comp (befehl, 'restart' , 5, lbef, 7) ) THEN 
      CALL get_params (zeile, ianz, cpara, lpara, maxw, length) 
      IF (ier_num == 0) THEN 
         IF (ianz == 1) THEN 
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
            IF (ier_num.eq.0) THEN 
               CALL set_gen(NINT(werte(1)))
            ENDIF 
         ELSE
            ier_num = -6
            ier_typ = ER_COMM
         ENDIF 
      ENDIF 
!                                                                 
!     -- set the result file                                      
!                                                                 
   ELSEIF (str_comp (befehl, 'reset'   , 5, lbef, 5) ) THEN 
      CALL diffev_do_reset
   ELSEIF (str_comp (befehl, 'restrial', 5, lbef, 8) ) THEN 
      CALL get_params (zeile, ianz, cpara, lpara, maxw, length) 
      IF (ier_num.eq.0) THEN 
         CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                           oname, loname, opara, lopara, lpresent, owerte)
         IF(ier_num==0) THEN
!
         IF(str_comp (cpara(1), 'silent',6,lpara(1),6)) THEN
            IF(lstandalone) THEN 
               ier_num = -27
               ier_typ = ER_APPL
               RETURN 
            ELSE
               pop_result_file_rd = .false.
               trial_results = ' '
               ltrial_results = 1
            ENDIF 
         ELSE
            CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
            IF (ier_num.eq.0) THEN 
               trial_results = cpara (1)(1:lpara(1)) 
               ltrial_results = lpara (1) 
            ENDIF 
         ENDIF 
         n_rvalue_i = NINT(owerte(O_PARTIAL))
!        IF(ianz>1) THEN
!           CALL del_params (1, ianz, cpara, lpara, maxw) 
!           CALL ber_params (ianz, cpara, lpara, werte, maxw) 
!           IF(ier_num/=0) RETURN
!           n_rvalue_i = NINT(werte(1))
!        ENDIF
      ENDIF 
      ENDIF 
!
!     -- Write a new set of children for the current generation
!
   ELSEIF (str_comp (befehl, 'write', 3, lbef, 5) ) THEN
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lcomm)
      IF (ier_num.eq.0) THEN
         IF (str_comp (cpara (1) , 'children', 3, lpara (1) , 8)) THEN
            CALL read_par_values              ! Make sure parent values are set
            CALL create_trial                 ! Make a new set
            CALL write_genfile                ! Write the "GENERATION" file
            CALL write_current                ! Update the Current parameter file
         ELSE IF (str_comp (cpara (1) , 'generation', 3, lpara (1) ,10)) THEN
            CALL write_genfile                ! Write the "GENERATION" file
         ELSE IF (str_comp (cpara (1) , 'kid', 3, lpara (1) ,3)) THEN
            cpara(1) = '0'
            lpara(1) = 1
            kid    = 1
            indiv  = 1
            nindiv = 1
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
            kid = NINT(werte(2))
            IF(ianz==3) THEN
               indiv  = NINT(werte(3))
            ENDIF
            IF(ianz==4) THEN
               nindiv = NINT(werte(4))
            ENDIF
            var_val( var_ref+0) = pop_gen
            var_val( var_ref+1) = pop_n
            var_val( var_ref+2) = pop_c
            var_val( var_ref+3) = pop_dimx
            var_val( var_ref+4) = kid
            var_val( var_ref+5) = indiv
            var_val( var_ref+6) = nindiv
            CALL write_trial(kid)             ! Update the Current parameter file
         ELSE
            ier_num = - 6
            ier_typ = ER_COMM
         ENDIF
      ENDIF
!                                                                       
!       Branch to DISCUS/ KUPLOT (standalone call system, suite do branch)
!                                                                       
         ELSEIF (str_comp (befehl, 'branch', 2, lbef, 6) ) THEN
            CALL p_branch (zeile, lcomm, .FALSE., 0     )
!                                                                 
!------   Try general commands                                    
!                                                                 
   ELSE 
      CALL kdo_all (befehl, lbef, zeile, lcomm) 
      IF(zeile == 'EXIT') THEN ! kdo_all detected "continue suite"
         lend = .TRUE.
      ENDIF 
   ENDIF 
ENDIF 
!
if(ex_do_exit) lend = .true.   ! A global exit was flagged
!                                                                       
END SUBROUTINE diffev_mache_kdo                      
!
!end module diffev_mache_kdo_mod
