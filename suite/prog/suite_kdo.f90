!module suite_mache_kdo_mod
!
!contains
!
!*****7*****************************************************************
!                                                                       
SUBROUTINE suite_mache_kdo (line, lend, length) 
!+                                                                      
!     This is the main routine for command interpretation, each         
!     command is identified here and the corresponding subroutine       
!     executed. A leading @ indicates a macro.                          
!-                                                                      
!USE discus,  discus_interactive => interactive
USE diffev_setup_mod
USE diffev_setup_sub_mod
USE diffev_loop_mod
USE diffev_reset
!
USE discus_setup_mod
USE discus_setup_sub_mod
USE discus_loop_mod
USE discus_reset_all_mod
!
use experi_setup_mod
use experi_setup_sub_mod
use experi_loop_mod
!
USE kuplot_setup_mod
USE kuplot_setup_sub_mod
USE kuplot_loop_mod
use kuplot_reset_mod
!
USE refine_reset
USE refine_setup_mod
USE refine_setup_sub_mod
USE refine_loop_mod
!
USE suite_init_mod
USE suite_parallel_mod
USE suite_setup_mod
USE suite_set_sub_mod
!
USE charact_mod 
!USE allocate_appl
!
USE blanks_mod
USE appl_env_mod
USE doact_mod
USE calc_expr_mod
USE errlist_mod 
use exit_para_mod
USE class_macro_internal
USE kdo_all_mod
USE learn_mod 
USE lib_errlist_func
USE lib_macro_func
USE prompt_mod
USE envir_mod
USE str_comp_mod
USE variable_mod
!USE set_sub_generic_mod
!
IMPLICIT none 
!                                                                       
CHARACTER (LEN= *  ), INTENT(INOUT) :: line 
LOGICAL             , INTENT(  OUT) :: lend 
INTEGER             , INTENT(INOUT) :: length 
!
CHARACTER (LEN=MAX(PREC_STRING,LEN(line))) :: zeile 
CHARACTER (LEN=   9)                  :: befehl 
INTEGER                               :: indxb, indxg, lcomm, lbef, indxt 
!                                                                       
!                                                                       
lstandalone = .false.  ! Switch to slave mode for DIFFEV/DISCUS/KUPLOT
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
         .not. (str_comp (befehl, 'echo',    2, lbef, 4) ) .and.   &
         .not. (str_comp (befehl, 'system',  2, lbef, 6) ) .and.   &
         .not. (str_comp (befehl, 'fput',    2, lbef, 4) ) .and.   &
         .not. (str_comp (befehl, 'help',    2, lbef, 4) .or.      &
          str_comp (befehl, '?   ',  2, lbef, 4) )       .AND.   &
         INDEX(line,'==') == 0                                ) THEN      
!         .not. (str_comp (befehl, 'socket',2, lbef, 5) ) .and.   &
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
         line = line(2:length)
         length = length - 1
         CALL file_kdo (line, length)
      ELSE 
         ier_num = - 13 
         ier_typ = ER_MAC 
      ENDIF 
!                                                                 
!-------Terminate DISCUS_SUITE 'exit'                                   
!                                                                 
   ELSEIF (str_comp (befehl, 'exit', 2, lbef, 4) ) then 
      LEND = .TRUE. 
   ELSEIF (str_comp (befehl, 'finished', 2, lbef, 8) ) then 
      LEND = .TRUE. 
!                                                                 
!     ----Start of DISCUS_SUITE special commands                        
!                                                                 

!                                                                 
!     -- branch to DIFFEV
!
   ELSEIF ((linteractive.OR.lblock.OR.lmakro) .AND. str_comp (befehl, 'diffev', 3, lbef, 6) ) then
       IF(suite_diffev_init) then
          pname     = 'diffev'
          pname_cap = 'DIFFEV'
          prompt    = pname
          oprompt   = pname
          CALL program_files ()
       ELSE
         CALL diffev_setup   (lstandalone)
         suite_diffev_init = .TRUE.
       ENDIF
       var_val(VAR_STATE)   = var_val(VAR_IS_SECTION)
       var_val(VAR_PROGRAM) = var_val(VAR_DIFFEV)
       CALL diffev_set_sub ()
       CALL suite_set_sub_branch
       CALL diffev_loop    ()
       lend      = .false.
       pname     = 'suite'
       pname_cap = 'SUITE'
       prompt    = pname
       oprompt   = pname
       var_val(VAR_STATE)   = var_val(VAR_IS_TOP)
       var_val(VAR_PROGRAM) = var_val(VAR_SUITE)
       CALL suite_set_sub
       IF(ier_num == -9 .AND. ier_typ == 1) THEN
          CALL program_files ()
          ier_num = -9
          ier_typ = ER_COMM
       ELSE
          CALL program_files ()
       ENDIF
!                                                                 
!     -- branch to DISCUS
!
   ELSEIF ((linteractive.OR.lblock.OR.lmakro) .AND. str_comp (befehl, 'discus', 3, lbef, 6) ) then
       IF(suite_discus_init) then
          pname     = 'discus'
          pname_cap = 'DISCUS'
          prompt    = pname
          oprompt   = pname
          CALL program_files ()
       ELSE
         CALL discus_setup   (lstandalone)
         suite_discus_init = .TRUE.
       ENDIF
       var_val(VAR_STATE)   = var_val(VAR_IS_SECTION)
       var_val(VAR_PROGRAM) = var_val(VAR_DISCUS)
       CALL discus_set_sub ()
       CALL suite_set_sub_branch
       CALL discus_loop    ()
       lend = .false.
       pname  = 'suite'
       pname_cap = 'SUITE'
       prompt = pname
       oprompt   = pname
       var_val(VAR_STATE)   = var_val(VAR_IS_TOP)
       var_val(VAR_PROGRAM) = var_val(VAR_SUITE)
       CALL suite_set_sub
       IF(ier_num == -9 .AND. ier_typ == 1) THEN
          CALL program_files ()
          ier_num = -9
          ier_typ = ER_COMM
       ELSE
          CALL program_files ()
       ENDIF
!                                                                 
!     -- branch to EXPERI
!
   ELSEIF ((linteractive.OR.lblock.OR.lmakro) .AND. str_comp (befehl, 'experi', 3, lbef, 6) ) then
       IF(suite_experi_init) then
          pname     = 'experi'
          pname_cap = 'EXPERI'
          prompt    = pname
          oprompt   = pname
          CALL program_files ()
       ELSE
         CALL experi_setup
         suite_refine_init = .TRUE.
       ENDIF
       var_val(VAR_STATE)   = var_val(VAR_IS_SECTION)
       var_val(VAR_PROGRAM) = var_val(VAR_EXPERI)
       CALL experi_set_sub ()
       CALL suite_set_sub_branch
       CALL experi_loop    ()
       lend      = .false.
       pname     = 'suite'
       pname_cap = 'SUITE'
       prompt    = pname
       oprompt   = pname
       var_val(VAR_STATE)   = var_val(VAR_IS_TOP)
       var_val(VAR_PROGRAM) = var_val(VAR_SUITE)
       CALL suite_set_sub
       IF(ier_num == -9 .AND. ier_typ == 1) THEN
          CALL program_files ()
          ier_num = -9
          ier_typ = ER_COMM
       ELSE
          CALL program_files ()
       ENDIF
!                                                                 
!     -- branch to KUPLOT
!
   ELSEIF ((linteractive.OR.lblock.OR.lmakro) .AND. str_comp (befehl, 'kuplot', 3, lbef, 6) ) then
       IF(suite_kuplot_init) then
          pname     = 'kuplot'
          pname_cap = 'KUPLOT'
          prompt    = pname
          oprompt   = pname
          CALL program_files ()
       ELSE
         CALL kuplot_setup   (lstandalone)
         suite_kuplot_init = .TRUE.
       ENDIF
       var_val(VAR_STATE)   = var_val(VAR_IS_SECTION)
       var_val(VAR_PROGRAM) = var_val(VAR_KUPLOT)
       CALL kuplot_set_sub ()
       CALL suite_set_sub_branch
       CALL kuplot_loop    ()
       lend   = .false.
       pname  = 'suite'
       pname_cap = 'SUITE'
       prompt = pname
       oprompt   = pname
       var_val(VAR_STATE)   = var_val(VAR_IS_TOP)
       var_val(VAR_PROGRAM) = var_val(VAR_SUITE)
       CALL suite_set_sub
       IF(ier_num == -9 .AND. ier_typ == 1) THEN
          CALL program_files ()
          ier_num = -9
          ier_typ = ER_COMM
       ELSE
          CALL program_files ()
       ENDIF
!                                                                 
!     -- branch to REFINE
!
   ELSEIF ((linteractive.OR.lblock.OR.lmakro) .AND. str_comp (befehl, 'refine', 3, lbef, 6) ) then
       IF(suite_refine_init) then
          pname     = 'refine'
          pname_cap = 'REFINE'
          prompt    = pname
          oprompt   = pname
          CALL program_files ()
       ELSE
         CALL refine_setup   (lstandalone)
         suite_refine_init = .TRUE.
       ENDIF
       var_val(VAR_STATE)   = var_val(VAR_IS_SECTION)
       var_val(VAR_PROGRAM) = var_val(VAR_REFINE)
       CALL refine_set_sub ()
       CALL suite_set_sub_branch
       CALL refine_loop    ()
       lend      = .false.
       pname     = 'suite'
       pname_cap = 'SUITE'
       prompt    = pname
       oprompt   = pname
       var_val(VAR_STATE)   = var_val(VAR_IS_TOP)
       var_val(VAR_PROGRAM) = var_val(VAR_SUITE)
       CALL suite_set_sub
       IF(ier_num == -9 .AND. ier_typ == 1) THEN
          CALL program_files ()
          ier_num = -9
          ier_typ = ER_COMM
       ELSE
          CALL program_files ()
       ENDIF
!                                                                 
!     -- Run a parallel version of discus_suite
!
   ELSEIF (str_comp (befehl, 'parallel', 3, lbef, 8) ) then
         CALL suite_do_parallel (zeile, lcomm)
!!!      CALL do_deallocate_appl (zeile, lcomm)
!
!     -- Test a simple forpy application
!
   ELSEIF ((linteractive.OR.lblock.OR.lmakro) .AND. str_comp (befehl, 'testforpy', 9, lbef, 9) ) then
      call suite_test_forpy! (zeile, lcomm, lend)
!
!     -- Test a simple jmol plot
!
   ELSEIF ((linteractive.OR.lblock.OR.lmakro) .AND. str_comp (befehl, 'testjmol', 8, lbef, 8) ) then
      call suite_test_jmol(lend)
!
!     -- Test a simple plot
!
   ELSEIF ((linteractive.OR.lblock.OR.lmakro) .AND. str_comp (befehl, 'testplot', 8, lbef, 8) ) then
      call suite_test_plot(zeile, lcomm, lend)
!                                                                 
!------   Reset the entire suite
!                                                                 
   ELSEIF(str_comp(befehl, 'reset', 3, lbef, 5)) THEN
      CALL kuplot_do_reset( 'all' , 3)
      CALL discus_reset_all
      CALL diffev_do_reset
      CALL refine_do_reset
!                                                                 
!------   Try general commands                                    
!                                                                 
   ELSE 
      CALL kdo_all (befehl, lbef, zeile, lcomm) 
      if(zeile=='stop') lend = .true.
   ENDIF 
ENDIF 
!
if(ex_do_exit) lend = .true.
!                                                                       
END SUBROUTINE suite_mache_kdo
!
!*******************************************************************************
!
subroutine suite_test_plot(zeile, lcomm, lend)
!-
! Do an automatic plot , serves to test if plot works
!+
!
use kuplot_setup_mod
use kuplot_setup_sub_mod
use kuplot_loop_mod
use kuplot_reset_mod
use kuplot_load_mod
use kuplot_plot_mod
use kuplot_mod
use kuplot_para_mod
!
use suite_init_mod
use suite_parallel_mod
use suite_setup_mod
use suite_set_sub_mod
!
use appl_env_mod
use errlist_mod
use learn_mod 
use prompt_mod
use precision_mod
use variable_mod
!
implicit none
character(len=*), intent(in)    :: zeile
integer         , intent(in)    :: lcomm
logical         , intent(  out) :: lend 
!
character(len=PREC_STRING) :: line
integer                    :: lp
!
if(suite_kuplot_init) then
   pname     = 'kuplot'
   pname_cap = 'KUPLOT'
   prompt    = pname
   oprompt   = pname
   call program_files ()
else
   call kuplot_setup   (lstandalone)
   suite_kuplot_init = .TRUE.
endif
var_val(VAR_STATE)   = var_val(VAR_IS_SECTION)
var_val(VAR_PROGRAM) = var_val(VAR_KUPLOT)
call kuplot_set_sub ()
call suite_set_sub_branch
if(zeile /= ' ') then
   line = zeile(1:lcomm)
   lp = lcomm
else
   line = 'sind(r[0]), 0,720,1'
   lp = 19
endif
call do_func(line, lp)
line = ' '
lp = 1
call set_skal(line, lp)
call set_mark(line, lp)
call do_plot (.false.)
iz = iz - 1
!      call kuplot_loop    ()
lend   = .false.
pname  = 'suite'
pname_cap = 'SUITE'
prompt = pname
oprompt   = pname
var_val(VAR_STATE)   = var_val(VAR_IS_TOP)
var_val(VAR_PROGRAM) = var_val(VAR_SUITE)
call suite_set_sub
if(ier_num == -9 .and. ier_typ == 1) then
   call program_files ()
   ier_num = -9
   ier_typ = ER_COMM
else
   call program_files ()
endif
!
end subroutine suite_test_plot
!
!
!*******************************************************************************
!
subroutine suite_test_jmol(lend)
!-
! Do an automatic jmol plot , serves to test if jmol works
!+
!
use discus_setup_mod
use discus_setup_sub_mod
use discus_loop_mod
use discus_reset_all_mod
!
use crystal_mod
use discus_plot_mod
use discus_plot_menu
use modify_mod
USE spcgr_apply
use structur
!use kuplot_load_mod
!use kuplot_plot_mod
!use kuplot_mod
!use kuplot_para_mod
!
use suite_init_mod
use suite_parallel_mod
use suite_setup_mod
use suite_set_sub_mod
!
use appl_env_mod
use errlist_mod
use learn_mod 
use prompt_mod
use precision_mod
use variable_mod
!
implicit none
!character(len=*), intent(in)    :: zeile
!integer         , intent(in)    :: lcomm
logical         , intent(  out) :: lend 
!
character(len=PREC_STRING) :: line
integer                    :: lp
!
if(suite_discus_init) then
   pname     = 'discus'
   pname_cap = 'DISCUS'
   prompt    = pname
   oprompt   = pname
   call program_files ()
else
   call discus_setup   (lstandalone)
   suite_discus_init = .TRUE.
endif
var_val(VAR_STATE)   = var_val(VAR_IS_SECTION)
var_val(VAR_PROGRAM) = var_val(VAR_DISCUS)
call discus_set_sub ()
call suite_set_sub_branch
!
CALL rese_cr
cr_name = 'freely created structure'
cr_spcgr (1:1)  = 'P'
cr_spcgr (2:2)  = '1'
cr_spcgr (3:16) = '              '
cr_set          = 'abc'
cr_spcgrno = 1
cr_syst = 1
spcgr_para = 1
CALL get_symmetry_matrices    ! Define matrices for default P1
cr_a0 = 3.5000D0
cr_win = 90.0D0
cr_icc (1) = 1
cr_icc (2) = 1
cr_icc (3) = 1
cr_natoms = 0
cr_ncatoms = 1
cr_ncreal  = 1
cr_nscat = 0
as_natoms = 0
line = 'Ni, 0.00, 0.00, 0.00, 0.5'
lp = 32
call do_ins(line, lp)
line = 'Ni, 0.50, 0.50, 0.00, 0.5'
lp = 32
call do_ins(line, lp)
line = 'Ni, 1.00, 0.00, 0.00, 0.5'
lp = 32
call do_ins(line, lp)
line = 'Ni, 0.00, 1.00, 0.00, 0.5'
lp = 32
call do_ins(line, lp)
line = 'Ni, 1.00, 1.00, 0.00, 0.5'
lp = 32
call do_ins(line, lp)
line = 'Ni, 0.50, 0.00, 0.50, 0.5'
lp = 32
call do_ins(line, lp)
line = 'Ni, 0.00, 0.50, 0.50, 0.5'
lp = 32
call do_ins(line, lp)
line = 'Ni, 0.50, 1.00, 0.50, 0.5'
lp = 32
call do_ins(line, lp)
line = 'Ni, 1.00, 0.50, 0.50, 0.5'
lp = 32
call do_ins(line, lp)
line = 'Ni, 0.00, 0.00, 1.00, 0.5'
lp = 32
call do_ins(line, lp)
line = 'Ni, 0.50, 0.50, 1.00, 0.5'
lp = 32
call do_ins(line, lp)
line = 'Ni, 1.00, 0.00, 1.00, 0.5'
lp = 32
call do_ins(line, lp)
line = 'Ni, 0.00, 1.00, 1.00, 0.5'
lp = 32
call do_ins(line, lp)
line = 'Ni, 1.00, 1.00, 1.00, 0.5'
lp = 32
call do_ins(line, lp)
!
call  plot_reset
!
pl_out = 'dummy_plot.cif'
line = 'all'
lp = 3
CALL atom_select (line, lp, 0, PL_MAXSCAT, pl_latom, &
                  pl_lsite, 0, PL_MAXSITE, &
                  pl_sel_atom, .true., .true.)
!
pl_ext_all = .true.
pl_prog = 'jmol'
call plot_test_jmol(pl_prog, pl_jmol)
!
line = 'plot:inter, kill:yes'
lp = 20
!write(*,*) ' RUN JMOL '
call plot_run(line, lp)
!
call plot_reset
call rese_cr
!
lend   = .false.
pname  = 'suite'
pname_cap = 'SUITE'
prompt = pname
oprompt   = pname
var_val(VAR_STATE)   = var_val(VAR_IS_TOP)
var_val(VAR_PROGRAM) = var_val(VAR_SUITE)
call suite_set_sub
if(ier_num == -9 .and. ier_typ == 1) then
   call program_files ()
   ier_num = -9
   ier_typ = ER_COMM
else
   call program_files ()
endif
!
end subroutine suite_test_jmol
!
!*******************************************************************************
!
subroutine suite_test_forpy
!-
!  Test a quick access to forpy library
!+
!
use errlist_mod
use precision_mod
!
use lib_forpython_mod
use forpy_mod
!
implicit none
!
type(ndarray)   :: p_arr   ! A python array
real(kind=PREC_DP), dimension(2,3) :: matrix
!
matrix(1,1) = 1.0
matrix(1,2) = 2.0
matrix(1,3) = 3.0
matrix(2,1) = 4.0
matrix(2,2) = 5.0
matrix(2,3) = 6.0
!
call forpy_start(ier_num)
ier_num = ndarray_create_ones(p_arr, [2, 2], dtype="float64", order="F")
ier_num = print_py(p_arr)
ier_num = ndarray_create(p_arr, matrix)
ier_num = print_py(p_arr)
call p_arr%destroy
!
end subroutine suite_test_forpy
!
!*******************************************************************************
!end module suite_mache_kdo_mod
