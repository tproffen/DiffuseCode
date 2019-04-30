!*****7*****************************************************************
!                                                                       
SUBROUTINE refine_mache_kdo (line, lend, length) 
!+                                                                      
!     This is the main routine for command interpretation, each         
!     command is identified here and the corresponding subroutine       
!     executed. A leading @ indicates a macro.                          
!-                                                                      
USE diffev_mpi_mod
!
USE refine_reset
USE refine_add_param_mod
USE refine_allocate_appl
USE refine_fix_mod
USE refine_load_mod
USE refine_run_mod
USE refine_set_mod
USE refine_show_mod
!
USE ber_params_mod
USE blanks_mod
USE build_name_mod
use calc_expr_mod
USE charact_mod 
USE define_variable_mod
USE errlist_mod 
USE gen_mpi_mod
USE get_params_mod
USE kdo_all_mod
USE learn_mod 
USE macro_mod
USE prompt_mod
USE set_sub_generic_mod
USE take_param_mod
USE variable_mod
!
IMPLICIT none 
!                                                                       
!                                                                       
INTEGER, PARAMETER :: MAXW = 20
LOGICAL, PARAMETER :: LDATA = .TRUE.
LOGICAL, PARAMETER :: LSIGMA = .FALSE.
!                                                                       
CHARACTER (LEN= *  ), INTENT(INOUT) :: line 
LOGICAL             , INTENT(  OUT) :: lend 
INTEGER             , INTENT(INOUT) :: length 
!
CHARACTER (LEN=1024)                  :: zeile   = ' '
CHARACTER (LEN=1024), DIMENSION(MAXW) :: cpara   = ' '
CHARACTER (LEN=   9)                  :: befehl  = ' '
INTEGER                               :: indxb, indxg, lcomm, lbef, indxt 
!INTEGER                               :: n_pop  ! dummy for allocation
!INTEGER                               :: kid, indiv, nindiv
INTEGER             , DIMENSION(MAXW) :: lpara = 0
!INTEGER, SAVE                         :: lastgen = -1
!LOGICAL                               :: back_new
!LOGICAL                               :: lexist
!LOGICAL                               :: lbest
!                                                                       
REAL                , DIMENSION(MAXW) :: werte = 0.0
LOGICAL, EXTERNAL                     :: str_comp 
!                                                                       
!INTEGER, PARAMETER :: NOPTIONAL = 4
!CHARACTER(LEN=1024), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
!CHARACTER(LEN=1024), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
!INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
!INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
!LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent  !opt. para present
!REAL               , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
!INTEGER, PARAMETER                        :: ncalc = 2 ! Number of values to calculate 
!
!DATA oname  / 'partial', 'repeat' , 'logfile', 'compute'  /
!DATA loname /  7       ,  6       ,  7       ,  7  /
!opara  =  (/ '0.000000', '1.000000', 'none    ', 'parallel' /)   ! Always provide fresh default values
!lopara =  (/  8        ,  8        ,  8        ,  8         /)
!owerte =  (/  0.0      ,  1.0      ,  0.0      ,  1.0       /)
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
is_math: IF (indxg.ne.0.and.                                              &
    &    .not. (str_comp (befehl, 'echo',  2, lbef, 4) ) .and.   &
    &    .not. (str_comp (befehl, 'syst',  2, lbef, 4) ) .and.   &
    &    .not. (str_comp (befehl, 'fput',  2, lbef, 4) ) .and.   &
    &    .not. (str_comp (befehl, 'socket',2, lbef, 5) ) .and.   &
    &    .not. (str_comp (befehl, 'help',  2, lbef, 4) .or.      &
    &     str_comp (befehl, '?   ',  2, lbef, 4) )       .AND.   &
          INDEX(line,'==') == 0                               ) THEN      
!                                                                 
!-------Zuweisung eines Funktionswertes                           
!                                                                 
   CALL do_math (line, indxg, length) 
ELSE  is_math
!
   cpara(:) = ' '
   lpara(:) = 0
   werte(:) = 0.0
!                                                                 
!     --execute a macro file                                      
!                                                                 
   is_befehl: IF (befehl (1:1) .eq.'@') THEN 
      IF (length.ge.2) THEN 
          CALL file_kdo (line(2:length), length -1)
!         CALL file_kdo (line (2:length), length - 1) 
      ELSE 
         ier_num = - 13 
         ier_typ = ER_MAC 
      ENDIF 
!                                                                 
!-------Terminate REFINE 'exit'                                   
!                                                                 
   ELSEIF (str_comp (befehl, 'exit', 2, lbef, 4) ) THEN is_befehl
      LEND = .TRUE. 
!                                                                 
!     ----Start of REFINE special commands                        
!                                                                 
!
!     -- Allocate array sizes
!
   ELSEIF (str_comp (befehl, 'allocate', 3, lbef,  8) ) THEN is_befehl
      CALL refine_do_allocate_appl (zeile, lcomm)
!
!     -- Deallocate array sizes
!
   ELSEIF (str_comp (befehl, 'deallocate', 3, lbef, 10) ) THEN is_befehl
      CALL refine_do_deallocate_appl (zeile, lcomm)
!                                                                 
!     -- set 
!                                                                 
   ELSEIF (str_comp (befehl, 'set', 3, lbef, 3) ) THEN  is_befehl
      CALL refine_set(zeile, length)
!
!     -- define data set
!
   ELSEIF (str_comp (befehl, 'data', 3, lbef, 4) ) THEN  is_befehl
      CALL refine_load(LDATA, zeile, length)
!
!     -- define data set
!
   ELSEIF (str_comp (befehl, 'sigma', 3, lbef, 5) ) THEN  is_befehl
      CALL refine_load(LSIGMA, zeile, length)
!                                                                 
!     -- Finish command will be ignored
!                                                                 
   ELSEIF (str_comp (befehl, 'finished', 3, lbef, 8) ) THEN  is_befehl
      CONTINUE
!                                                                 
!     -- Fix a parameter
!                                                                 
   ELSEIF (str_comp (befehl, 'fix', 3, lbef, 3) ) THEN  is_befehl
      CALL refine_fix(zeile, length)
!                                                                 
!     -- add a new parameter to the dimension                     
!                                                                 
   ELSEIF (str_comp (befehl, 'newparam', 3, lbef, 8) ) THEN  is_befehl
      CALL refine_add_param(zeile, length)
!                                                                 
!     -- set the name of a refinement parameter                   
!                                                                 
   ELSEIF (str_comp (befehl, 'run', 3, lbef, 3) ) THEN is_befehl
!
      CALL refine_run(zeile, length)
!                                                                 
!-------  Show parameters 'show'                                  
!                                                                 
   ELSEIF (str_comp (befehl, 'show', 3, lbef, 4) ) THEN  is_befehl
      CALL refine_do_show (zeile, lcomm) 
!                                                                 
!                                                                       
!       Branch to DISCUS/ KUPLOT (standalone call system, suite do branch)
!                                                                       
   ELSEIF (str_comp (befehl, 'branch', 2, lbef, 6) ) THEN is_befehl
      CALL p_branch (zeile, lcomm, .FALSE.)
!                                                                 
!------   Try general commands                                    
!                                                                 
   ELSE  is_befehl
      CALL kdo_all (befehl, lbef, zeile, lcomm) 
      IF(zeile == 'EXIT') THEN ! kdo_all detected "continue suite"
         lend = .TRUE.
      ENDIF 
   ENDIF  is_befehl
ENDIF  is_math
!                                                                       
END SUBROUTINE refine_mache_kdo                      
