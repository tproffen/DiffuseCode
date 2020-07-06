MODULE diffev_blk_appl
!                                                                       
USE population
USE diff_evol
!
PRIVATE
PUBLIC :: diffev_autodef 
PUBLIC :: diffev_initarrays
PUBLIC :: init_diffev
PUBLIC :: init_population
!
CONTAINS
!
   SUBROUTINE diffev_initarrays
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
   CALL init_diffev
   CALL init_population
!
   END SUBROUTINE diffev_initarrays
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
   pop_trialfile  = ' '
   trial_results  = ' '
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
   USE variable_mod
   IMPLICIT none 
!                                                                 
   INTEGER        :: i, j 
!                                                                 
!                                                                 
   pop_n     = 1 
   pop_c     = pop_n 
   lastgen   = -1
   pop_dimx  = 1 
   pop_best  = 1 
   pop_worst = 1 
   var_val(var_ref+1) = pop_n     ! Update the global user variable
   var_val(var_ref+2) = pop_c     ! Update the global user variable
   var_val(var_ref+3) = pop_dimx  ! Update the global user variable
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
   pop_name(0) = 'Rvalue'
   DO i = 1, MAXDIMX 
      WRITE (pop_name (i), 1000) i 
      DO j = 1, 8 
         IF (pop_name (i) (j:j) .eq.' ') pop_name (i) (j:j) = '0' 
      ENDDO 
   ENDDO 
!                                                                 
! Initialize array of size MAXPOP
!
   child_val(:,:)  = 0.0 
   trial_val(:,:)  = 0.0 
   parent_val(:,:) = 0.0 
   n_rvalue_i      = 1
   n_rvalue_o      = 1
!                                                                       
    1000 FORMAT    ('PARA',i4.4) 
!                                                                       
   END SUBROUTINE init_population                     
!*****7*****************************************************************
   SUBROUTINE diffev_autodef 
!-                                                                      
!     Tries to open a default file for the integer and real variables   
!+                                                                      
      USE envir_mod 
      USE errlist_mod 
USE lib_errlist_func
      USE param_mod 
USE support_mod
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
   CLOSE (idef) 
!
!
   END SUBROUTINE diffev_autodef                        
!
END MODULE diffev_blk_appl
