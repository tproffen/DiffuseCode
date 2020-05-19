MODULE refine_blk_appl
!                                                                       
!
PRIVATE
PUBLIC :: refine_autodef 
PUBLIC :: refine_initarrays
PUBLIC :: init_refine
!
CONTAINS
!
SUBROUTINE refine_initarrays
!
!       Startvalues for all important arrays.                           
!                                                                       
USE prompt_mod 
!
IMPLICIT NONE
!                                                                 
!------ Sockt defaults                                            
!                                                                 
s_port = 3333 
s_ipallowed = "localhost" 
!
!   CALL init_refine
!
END SUBROUTINE refine_initarrays
!
SUBROUTINE init_refine 
!                                                                       
!  Initializes the variable for refine 
!
IMPLICIT none 
!                                                                 
!     /refine/                                             
!                                                                 
!  diff_cr         = 0.9 
!  diff_f          = 0.81 
!  diff_k          = 1.0 
!  diff_local      = 0.0 
!  diff_donor_mode = ADD_TO_RANDOM 
!  diff_sel_mode   = SEL_COMP 
!                                                                 
!  pop_genfile    = 'GENERATION'
!  pop_trialfile  = ' '
!  trial_results  = ' '
!  parent_results = 'LOGFILE'
!  parent_summary = 'AVERAGES'
!
!                                                                       
END SUBROUTINE init_refine                     
!
!*****7*****************************************************************
!
SUBROUTINE refine_autodef 
!-                                                                      
!     Tries to open a default file for the integer and real variables   
!+                                                                      
USE envir_mod 
USE errlist_mod 
USE param_mod 
USE sys_compiler
!
IMPLICIT NONE
!
!                                                                       
INTEGER              :: idef  = 34
!                                                                       
ier_num = 0 
ier_typ = ER_NONE 
CALL open_def (idef) 
IF (ier_num.eq.0) then 
   ier_num = - 3 
   ier_typ = ER_IO 
   READ (idef, *, end = 998, err = 998) inpara 
   READ (idef, *, end = 998, err = 998) rpara 
ENDIF 
ier_num = 0 
ier_typ = ER_NONE 
998 CONTINUE 
IF (ier_num.ne.0) then 
   CALL errlist 
   ier_num = 0 
   ier_typ = ER_NONE 
ENDIF 
!
CLOSE (idef) 
!
!
END SUBROUTINE refine_autodef                        
!
END MODULE refine_blk_appl
