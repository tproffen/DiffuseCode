MODULE blk_appl
!                                                                       
USE population
USE diff_evol
!
PRIVATE
PUBLIC :: autodef 
PUBLIC :: initarrays
PUBLIC :: init_diffev
PUBLIC :: init_population
!
CONTAINS
!
   SUBROUTINE initarrays
!
!       Startvalues for all important arrays.                           
!                                                                       
   IMPLICIT NONE
!
   include'prompt.inc' 
!                                                                 
!------ Sockt defaults                                            
!                                                                 
   s_port = 3333 
   s_ipallowed = "localhost" 
!
   CALL init_diffev
   CALL init_population
!
   END SUBROUTINE initarrays
!
   SUBROUTINE init_diffev 
!                                                                       
!  Initializes the variable for diffev 
!
   IMPLICIT none 
!                                                                 
!     /diffev/                                             
!                                                                 
   diff_cr         = 0.9 
   diff_f          = 0.81 
   diff_k          = 1.0 
   diff_local      = 0.0 
   diff_donor_mode = ADD_TO_RANDOM 
   diff_sel_mode   = SEL_COMP 
!                                                                 
   pop_genfile    = 'GENERATION'
   pop_trialfile  = 'TRIALFILE'
   trial_results  = 'RESULTFILE'
   parent_results = 'LOGFILE'
   parent_summary = 'AVERAGES'
!
!                                                                       
   END SUBROUTINE init_diffev                     
!
   SUBROUTINE init_population 
!                                                                       
!+                                                                      
!   Initializes the population
!-                                                                      
!
   IMPLICIT none 
!                                                                 
   INTEGER        :: i, j 
!                                                                 
!                                                                 
   pop_n     = 1 
   pop_c     = pop_n 
   pop_dimx  = 1 
   pop_best  = 1 
   pop_worst = 1 
!                                                                 
   pop_current = .false. 
!                                                                 
! Initialize array of size MAXDIMX
!
   pop_ad_sigma   = .true. 
   pop_ad_lsigma  = .true. 
   pop_sigma      = 1.0 
   pop_lsig       = 1.0 
   pop_sig_ad     = 0.5 
   pop_lsig_ad    = 0.05 
   pop_xmin       = 0.0 
   pop_xmax       = 0.0 
!                                                                 
! Initialize array of size MAXPOP
!
   pop_x  = 0.0 
   pop_t  = 0.0 
   child  = 0.0 
   trial  = 0.0 
!                                                                 
   DO i = 1, MAXDIMX 
      WRITE (pop_name (i), 1000) i 
      DO j = 1, 8 
         IF (pop_name (i) (j:j) .eq.' ') pop_name (i) (j:j) = '0' 
      ENDDO 
   ENDDO 
!                                                                 
! Initialize array of size MAXPOP
!
   child_val  = 0.0 
   trial_val  = 0.0 
   parent_val = 0.0 
!                                                                       
    1000 FORMAT    ('PARA',i4) 
!                                                                       
   END SUBROUTINE init_population                     
!*****7*****************************************************************
   SUBROUTINE autodef 
!-                                                                      
!     Tries to open a default file for the integer and real variables   
!+                                                                      
      include'envir.inc' 
      include'errlist.inc' 
      include'param.inc' 
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
   CLOSE (idef) 
!
    1000 FORMAT    (a) 
!
   END SUBROUTINE autodef                        
!
END MODULE blk_appl
