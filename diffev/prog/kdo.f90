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
USE diff_evol
USE diffev_mpi_mod
USE population
USE diffev_allocate_appl
USE create_trial_mod
USE compare
USE initialise
USE run_mpi_mod
use show_mod
!
USE errlist_mod 
USE learn_mod 
USE macro_mod
IMPLICIT none 
!                                                                       
!                                                                       
INTEGER, PARAMETER   :: maxw = 20
!                                                                       
CHARACTER (LEN= *  ), INTENT(INOUT) :: line 
LOGICAL             , INTENT(  OUT) :: lend 
INTEGER             , INTENT(INOUT) :: length 
!
CHARACTER (LEN=1024)                  :: zeile 
CHARACTER (LEN=1024), DIMENSION(MAXW) :: cpara
CHARACTER (LEN=   9)                  :: befehl 
INTEGER                               :: indxb, indxg, lcomm, lbef, indxt 
INTEGER                               :: i, j, ii , nb
INTEGER                               :: n_pop  ! dummy for allocation
INTEGER                               :: lb,ub
INTEGER                               :: ianz 
INTEGER                               :: iianz 
INTEGER                               :: lpara (maxw) 
LOGICAL                               :: lbest
!                                                                       
REAL                                  :: werte (maxw) 
REAL                                  :: value
LOGICAL, EXTERNAL                     :: str_comp 
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
IF (indxb + 1.le.length) then 
   zeile = line (indxb + 1:length) 
   lcomm = length - indxb 
   call rem_leading_bl(zeile, lcomm)
ENDIF 
!                                                                 
!-------Suche nach einem "="                                      
!                                                                 
indxg = index (line, '=') 
IF (indxg.ne.0.and.                                              &
    &    .not. (str_comp (befehl, 'echo',  2, lbef, 4) ) .and.   &
    &    .not. (str_comp (befehl, 'syst',  2, lbef, 4) ) .and.   &
    &    .not. (str_comp (befehl, 'fput',  2, lbef, 4) ) .and.   &
    &    .not. (str_comp (befehl, 'socket',2, lbef, 5) ) .and.   &
    &    .not. (str_comp (befehl, 'help',  2, lbef, 4) .or.      &
    &     str_comp (befehl, '?   ',  2, lbef, 4) )      ) then      
!                                                                 
!-------Zuweisung eines Funktionswertes                           
!                                                                 
   CALL do_math (line, indxg, length) 
ELSE 
!                                                                 
!     --execute a macro file                                      
!                                                                 
   IF (befehl (1:1) .eq.'@') then 
      IF (length.ge.2) then 
          CALL file_kdo (line(2:length), length -1)
!         CALL file_kdo (line (2:length), length - 1) 
      ELSE 
         ier_num = - 13 
         ier_typ = ER_MAC 
      ENDIF 
!                                                                 
!-------Terminate DIFFEV 'exit'                                   
!                                                                 
   ELSEIF (str_comp (befehl, 'exit', 2, lbef, 4) ) then 
      LEND = .TRUE. 
!                                                                 
!     ----Start of DIFFEV special commands                        
!                                                                 

!     -- Allocate array sizes
!
   ELSEIF (str_comp (befehl, 'allocate', 3, lbef,  8) ) then
      CALL diffev_do_allocate_appl (zeile, lcomm)
!
!     -- Deallocate array sizes
!
   ELSEIF (str_comp (befehl, 'deallocate', 3, lbef, 10) ) then
      CALL diffev_do_deallocate_appl (zeile, lcomm)
!                                                                 
!     -- define adaptation of sigmas 'adapt'                      
!                                                                 
   ELSEIF (str_comp (befehl, 'adapt', 3, lbef, 5) ) then 
      CALL get_params (zeile, ianz, cpara, lpara, maxw, length) 
      IF (ier_num.eq.0) then 
         IF (ianz.eq.2.or.ianz.eq.3) then 
            IF (str_comp (cpara(1),'sigma',3,lpara (1) , 5) ) then
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               iianz = 1 
               CALL ber_params (iianz, cpara, lpara, werte, maxw) 
               ii = nint (werte (1) ) 
               IF(ii==0 .or. ii>pop_dimx .or. ii>MAXDIMX) THEN
                 ier_num = -14
                 ier_typ = ER_APPL
                 RETURN
               ENDIF
               IF (ianz.eq.1.or.str_comp (cpara (2) , 'yes', 1, lpara (1) , 3) ) then
                  pop_ad_sigma (ii) = .true. 
               ELSEIF (ianz.eq.2.and.str_comp (cpara (2) , 'no', 2, lpara (2) , 2) ) then
                  pop_ad_sigma (ii) = .false. 
               ELSE 
                  CALL del_params (1, ianz, cpara, lpara, maxw) 
                  CALL ber_params (ianz, cpara, lpara, werte,maxw)
                  IF (ier_num.eq.0) then 
                     pop_ad_sigma (ii) = .true. 
                     pop_sig_ad (ii) = werte (1) 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ENDIF 
            ELSEIF (str_comp (cpara(1),'lsig', 3, lpara(1),4))then
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               iianz = 1 
               CALL ber_params (iianz, cpara, lpara, werte, maxw) 
               ii = nint (werte (1) ) 
               IF(ii==0 .or. ii>pop_dimx .or. ii>MAXDIMX) THEN
                 ier_num = -14
                 ier_typ = ER_APPL
                 RETURN
               ENDIF
               IF (ianz.eq.1.or.str_comp (cpara (2) , 'yes', 1, lpara (1) , 3) ) then
                  pop_ad_lsigma (ii) = .true. 
               ELSEIF (ianz.eq.2.and.str_comp (cpara (2) , 'no', 2, lpara (2) , 2) ) then
                  pop_ad_lsigma (ii) = .false. 
               ELSE 
                  CALL del_params (1, ianz, cpara, lpara, maxw) 
                  CALL ber_params (ianz, cpara, lpara, werte,maxw)
                  IF (ier_num.eq.0) then 
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
   ELSEIF (str_comp (befehl, 'backup', 3, lbef, 6) ) then 
      CALL get_params (zeile, ianz, cpara, lpara, maxw, length) 
      IF (ier_num.eq.0) then 
         IF(cpara(1)=='NONE') THEN
            pop_backup = .false.
         ELSE
            IF(pop_back_number==MAXBACK) THEN
               CALL alloc_backup(n_pop)
            ENDIF
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
         ENDIF 
      ENDIF 
!                                                                 
!     -- compare the trial results to the last generation         
!                                                                 
   ELSEIF (str_comp (befehl, 'comp', 3, lbef, 4) ) then 
      CALL get_params (zeile, ianz, cpara, lpara, maxw, length) 
      IF (ier_num.eq.0) then 
         IF (ianz.eq.0) THEN
            pop_result_file_rd = .true.
         ELSEIF (ianz.eq.1) THEN
            IF(str_comp (cpara(1), 'silent',6,lpara(1),6)) THEN
               pop_result_file_rd = .false.
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
!     -- Define a hard constraint 'constraint'                    
!                                                                 
   ELSEIF (str_comp (befehl, 'const', 3, lbef, 5) ) then 
      IF(constr_number >= MAX_CONSTR) THEN
         CALL alloc_constraint ( constr_number + 1)
         IF(ier_num < 0) then
            RETURN
         ENDIF
      ENDIF
      IF (constr_number.lt.MAX_CONSTR) then 
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
   ELSEIF (str_comp (befehl, 'dismiss', 3, lbef, 7) ) then 
      IF (pop_n.gt.3) then
         CALL get_params (zeile, ianz, cpara, lpara, maxw, length) 
         IF (ier_num.eq.0) then 
            IF (ianz.eq.0 .or.  &
                str_comp (cpara(1), 'all',3, lpara(1), 3)) then 
               lb = 1
               ub = pop_n
               CALL do_dismiss  (lb, ub)
            ELSE
               CALL ber_params (ianz, cpara, lpara, werte, maxw)
               IF (ier_num.eq.0) then
                  IF (ianz.eq.1) then 
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
   ELSEIF (str_comp (befehl, 'donor', 4, lbef, 4) ) then 
      CALL get_params (zeile, ianz, cpara, lpara, maxw, length) 
      IF (ier_num.eq.0) then 
         IF (ianz.eq.1) then 
            IF (str_comp (cpara(1),'best',3,lpara(1),4) ) then
               diff_donor_mode = ADD_TO_BEST 
            ELSEIF (str_comp(cpara(1),'random',3,lpara(1),6)) then
               diff_donor_mode = ADD_TO_RANDOM 
            ENDIF 
         ENDIF 
      ENDIF 
!                                                                 
!     -- fix a parameter at a value                               
!                                                                 
   ELSEIF (str_comp (befehl, 'fix', 3, lbef, 3) ) then 
      IF (pop_n.gt.3) then
         CALL get_params (zeile, ianz, cpara, lpara, maxw, length) 
         IF (ier_num == 0) then 
            IF (ianz == 2) then 
               lbest = .false.
               IF (str_comp (cpara(2), 'best', 3, lpara(2), 4) ) then 
                  lbest = .true.
                  ianz  = ianz - 1
                  CALL ber_params (ianz, cpara, lpara, werte, maxw)
                  IF (ier_num.eq.0) then
                     lb = nint(werte(1))
                  ELSE
                     ier_num = -6
                     ier_typ = ER_COMM
                     RETURN 
                  ENDIF 
               ELSE
                  CALL ber_params (ianz, cpara, lpara, werte, maxw)
                  IF (ier_num.eq.0) then
                     lb = nint(werte(1))
                  ELSE
                     ier_num = -6
                     ier_typ = ER_COMM
                     RETURN 
                  ENDIF 
               ENDIF 
               IF ( 0<lb .and. lb<=pop_dimx) THEN
                  CALL do_read_values   ! Read values from logfile
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
            ENDIF 
         ENDIF 
      ELSE 
         ier_num = - 3 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
!                                                                 
!     -- initialise the sequence                                  
!                                                                 
   ELSEIF (str_comp (befehl, 'init', 3, lbef, 4) ) then 
      IF (pop_n.gt.3) then
         CALL get_params (zeile, ianz, cpara, lpara, maxw, length) 
         IF (ier_num.eq.0) then 
            IF (ianz.eq.0) then 
               CALL do_initialise 
            ELSE
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
               CALL ber_params (ianz, cpara, lpara, werte, maxw)
               IF (ier_num.eq.0) then
                  IF (ianz.eq.1) then 
                     lb = nint(werte(1))
                     ub = nint(werte(1))
                  ELSEIF (ianz.eq.2 ) THEN
                     lb = nint(werte(1))
                     ub = nint(werte(2))
                  ELSE
                     ier_num = -6
                     ier_typ = ER_COMM
                     RETURN 
                  ENDIF 
                  IF ( 0<lb .and. lb<=ub .and. ub<=pop_dimx) THEN
                     CALL do_read_values   ! Read values from logfile
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
      ELSE 
         ier_num = - 3 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
!                                                                 
!     -- set the logfile file                                     
!                                                                 
   ELSEIF (str_comp (befehl, 'logfile', 3, lbef, 7) ) then 
      CALL get_params (zeile, ianz, cpara, lpara, maxw, length) 
      IF (ier_num.eq.0) then 
         CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
         IF (ier_num.eq.0) then 
            parent_results = cpara (1)(1:lpara(1)) 
            lparent_results = lpara (1) 
         ENDIF 
      ENDIF 
!                                                                 
!     -- set the name of a refinement parameter                   
!                                                                 
   ELSEIF (str_comp (befehl, 'pop_name', 3, lbef, 8) ) then 
      CALL get_params (zeile, ianz, cpara, lpara, maxw, length) 
      IF (ier_num.eq.0) then 
         iianz = 1 
         CALL ber_params (iianz, cpara, lpara, werte, maxw) 
         IF (ier_num.eq.0) then 
            IF (nint(werte(1)) >  0       .and.         &
                nint(werte(1)) <= MAXDIMX .and.         &
                nint(werte(1)) <= pop_dimx     ) THEN
            pop_lname (nint (werte (1) ) ) = min (lpara (2), 8)
            pop_name (nint (werte (1) ) ) = cpara (2)(1:pop_lname(nint (werte (1) ) )) 
            ELSE
               ier_num = -14
               ier_typ = ER_APPL
            ENDIF 
         ENDIF 
      ENDIF 
   ELSEIF (str_comp (befehl, 'run_mpi', 7, lbef, 7) ) then 
      IF ( .not. run_mpi_active .or. run_mpi_numprocs < 2 ) THEN
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
      CALL get_params (zeile, ianz, cpara, lpara, maxw, length) 
      IF (ier_num.eq.0) then 
         IF ( ianz == 5 ) THEN
            run_mpi_senddata%use_socket = str_comp(cpara(5), 'socket', 6, lpara(5), 6)
            ianz = 4
         ELSE
            run_mpi_senddata%use_socket = .false.
         ENDIF
         IF ( ianz >= 2 .and. ianz <= 4 ) THEN
            run_mpi_senddata%generation = pop_gen    ! Current GENERATION no
            run_mpi_senddata%member     = pop_n      ! Number of members
            run_mpi_senddata%children   = pop_c      ! Number of children
            run_mpi_senddata%parameters = pop_dimx   ! Number of parameters
            run_mpi_senddata%prog   = cpara(1)(1:100) ! Program to run
            run_mpi_senddata%prog_l = lpara(1)
            run_mpi_senddata%mac    = cpara(2)(1:100) ! Macro to run
            run_mpi_senddata%mac_l  = lpara(2)
            IF ( ianz == 4 ) THEN
               run_mpi_senddata%out   = cpara(4)(1:100)! Target for program output 
               run_mpi_senddata%out_l = lpara(4)
            ELSE
               run_mpi_senddata%out   = '/dev/null'   ! Default output
               run_mpi_senddata%out_l = 9
            ENDIF 
            run_mpi_senddata%repeat = .false.! repeat = .false. means no repetition
            run_mpi_senddata%nindiv = 1      ! nindiv is at least 1
            IF ( ianz >  2 ) THEN            ! Number of local repeations
               CALL del_params (2, ianz, cpara, lpara, maxw) 
               ianz = 1
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               IF (ier_num.eq.0) then 
                  IF ( nint(werte(1)) > 0 ) THEN
                     run_mpi_senddata%repeat = .true.             ! repeat = false means no repetition
                  ENDIF 
                  run_mpi_senddata%nindiv = max(1,nint(werte(1))) ! nindiv is at least 1
               ENDIF 
            ENDIF 
            IF ( ier_num == 0) THEN
               call run_mpi_master 
            ENDIF 
         ELSE
            ier_num =  -6
            ier_typ =  -ER_COMM
         ENDIF 
      ENDIF 
!                                                                 
!-------  Show parameters 'show'                                  
!                                                                 
   ELSEIF (str_comp (befehl, 'show', 3, lbef, 4) ) then 
      CALL diffev_do_show (zeile, lcomm) 
!                                                                 
!     -- set the average information file                         
!                                                                 
   ELSEIF (str_comp (befehl, 'summaryfile', 3, lbef, 11) ) then 
      CALL get_params (zeile, ianz, cpara, lpara, maxw, length) 
      IF (ier_num.eq.0) then 
         CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
         IF (ier_num.eq.0) then 
            parent_summary = cpara (1)(1:lpara(1)) 
            lparent_summary = lpara (1) 
         ENDIF 
      ENDIF 
!                                                                 
!     -- define the refine status of each gene                    
!                                                                 
   ELSEIF (str_comp (befehl, 'refine', 2, lbef, 6) ) then 
      CALL get_params (zeile, ianz, cpara, lpara, maxw, length) 
      IF (ier_num.eq.0) then 
         IF (ianz.ge.1) then 
            IF (str_comp (cpara(1),'all',1, lpara (1) , 3) ) then
               DO i = 1, MAXDIMX 
               pop_refine (i) = .true. 
               ENDDO 
            ELSEIF (str_comp (cpara(1),'none',1,lpara(1),3) ) then
               DO i = 1, MAXDIMX 
               pop_refine (i) = .false. 
               ENDDO 
            ELSE 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               IF (ier_num.eq.0) then 
                  DO i = 1, ianz 
                  j = nint (werte (i) ) 
                  IF(j==0 .or. j>pop_dimx .or. j>MAXDIMX) THEN
                    ier_num = -14
                    ier_typ = ER_APPL
                    RETURN
                  ENDIF
                  IF (j.gt.0) then 
                     pop_refine (j) = .true. 
                  ELSEIF (j.lt.0) then 
                     pop_refine (abs (j) ) = .false. 
                  ENDIF 
                  ENDDO 
               ENDIF 
            ENDIF 
         ENDIF 
      ENDIF 
!                                                                 
!     -- set the selection mode 'selection'                       
!                                                                 
   ELSEIF (str_comp (befehl, 'selection', 2, lbef, 9) ) then 
      CALL get_params (zeile, ianz, cpara, lpara, maxw, length) 
      IF (ier_num.eq.0) then 
         IF (ianz.ge.1) then 
            IF (str_comp (cpara(1),'best',3,lpara(1),4)) then
               IF (str_comp (cpara(2),'all',3,lpara(2), 3) ) then
                  diff_sel_mode = SEL_BEST_ALL 
               ELSEIF (str_comp(cpara(2),'children',3,lpara(2), 8) ) then
                  diff_sel_mode = SEL_BEST_CHILD 
               ENDIF 
            ELSEIF (str_comp(cpara(1),'compare',3,lpara(1),7)) then
               diff_sel_mode = SEL_COMP 
            ENDIF 
         ENDIF 
      ENDIF 
!                                                                 
!     -- set the trial file                                       
!                                                                 
   ELSEIF (str_comp (befehl, 'trialfile', 3, lbef, 9) ) then 
      CALL get_params (zeile, ianz, cpara, lpara, maxw, length) 
      IF (ier_num.eq.0) then 
         CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
         IF (ier_num.eq.0) then 
            pop_trialfile = cpara (1)(1:lpara(1)) 
            pop_ltrialfile = lpara (1) 
         ENDIF 
      ENDIF 
!                                                                 
!     -- set the parameter type                                   
!                                                                 
   ELSEIF (str_comp (befehl, 'type', 4, lbef, 4) ) then 
      CALL get_params (zeile, ianz, cpara, lpara, maxw, length) 
      IF (ier_num.eq.0) then 
         IF (ianz.eq.2) then 
            IF (str_comp(cpara(1),'integer', 3, lpara (1) , 7)) then
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               IF (ier_num.eq.0) then 
                  pop_type (nint (werte (1) ) ) = POP_INTEGER 
               ENDIF 
            ELSEIF (str_comp(cpara(1),'real',3,lpara(1),4) ) then
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               IF (ier_num.eq.0) then 
                  pop_type (nint (werte (1) ) ) = POP_REAL 
               ENDIF 
            ENDIF 
         ENDIF 
      ENDIF 
!                                                                 
!     -- set the result file                                      
!                                                                 
   ELSEIF (str_comp (befehl, 'restrial', 3, lbef, 8) ) then 
      CALL get_params (zeile, ianz, cpara, lpara, maxw, length) 
      IF (ier_num.eq.0) then 
         CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
         IF (ier_num.eq.0) then 
            trial_results = cpara (1)(1:lpara(1)) 
            ltrial_results = lpara (1) 
         ENDIF 
      ENDIF 
!
!     -- Write a new set of children for the current generation
!
   ELSEIF (str_comp (befehl, 'write', 3, lbef, 5) ) then
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lcomm)
      IF (ier_num.eq.0) then
         IF (str_comp (cpara (1) , 'children', 3, lpara (1) , 8)) then
            CALL create_trial                 ! Make a new set
            CALL write_genfile                ! Write the "GENERATION" file
         ELSE IF (str_comp (cpara (1) , 'generation', 3, lpara (1) ,10)) then
            CALL write_genfile                ! Write the "GENERATION" file
         ELSE
            ier_num = - 6
            ier_typ = ER_COMM
         ENDIF
      ENDIF
!                                                                 
!------   Try general commands                                    
!                                                                 
   ELSE 
      CALL kdo_all (befehl, lbef, zeile, lcomm) 
   ENDIF 
ENDIF 
!                                                                       
END SUBROUTINE diffev_mache_kdo                      
