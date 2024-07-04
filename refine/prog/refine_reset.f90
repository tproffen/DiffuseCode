MODULE refine_reset
!
CONTAINS
!
SUBROUTINE refine_do_reset
!-
!   Reset REFINE to system start
!+
!
USE refine_allocate_appl
USE refine_blk_appl
USE refine_control_mod
USE refine_data_mod
USE refine_fit_erg
use refine_head_mod
use refine_log_mod
USE refine_params_mod
USE do_variable_mod
USE precision_mod
!
IMPLICIT NONE
!
CHARACTER (LEN=PREC_STRING)                   :: zeile
!
INTEGER              :: lcomm
!
! Remove all parameter names from the variable entry
!
zeile ='default'
lcomm = 7
!
CALL refine_do_allocate_appl(zeile,lcomm)
CALL refine_initarrays
!
! refine_params_mod
!
refine_autoconstr =.TRUE.   ! Do automatic constraints
refine_fwhm       = .FALSE. ! u,v,w, are refined=TRUE or fixed=FALSE
refine_fwhm_ind   = HUGE(0)
refine_eta        = .FALSE. ! eta    are refined=TRUE or fixed=FALSE
refine_eta_ind    = HUGE(0)
refine_par_n    = 0 ! number of parameters
refine_fix_n    = 0 
ref_dim         = 1
!
! refine_fit_erg
!
refine_chisqr   = -1.0  ! Sum of dev^2 in previous cycle
refine_conf     = -1.0  ! sum of dev^2
refine_lamda    = -1.0  ! Final convergence criterion
refine_rval     = -1.0  ! weighted R-value
refine_rexp     = -1.0  ! unweighted R-value
!
!refine_control_mod
!
refine_cycles   = 1     ! Maximum number of cycles defined by user
conv_status     = .TRUE.
conv_dp_sig     = 0.005 ! Maximum DeltaP/sigma for convergence
conv_dchi2      = 0.5   ! Maximum Chi^2 change for convergence
conv_chi2       = 0.5   ! Minimum Chi^2 value  for convergence
conv_conf       = 0.01  ! Minimum confidence level
conv_lambda     = 1.00E10 ! Maximum lambda     level
lconvergence    = .FALSE. ! Convergence has been reached
lconv           = .FALSE. ! Convergence has been reached
refine_lamda_s  = 0.001   ! MRQ lamda start
refine_lamda_u  = 16.00   ! MRQ lamda Increase (up)
refine_lamda_d  = 0.500   ! MRQ lamda Decrease (down)
refine_init     = .TRUE.  ! Initialize MRQ
!
! refine_data_mod
!
ref_load_u  = ' '  ! Load string data
ref_csigma_u= ' '  ! Load string Sigma's
ref_load    = ' '  ! Load string data
ref_csigma  = ' '  ! Load string Sigma's
ref_kload   = 0    ! Data set within KUPLOT
ref_ksigma  = 0    ! Sigma set within KUPLOT
ref_kupl    = 0    ! Data set within KUPLOT that needs to be kept
!
ref_run_u   = ' '
conv_dp_sig_u  = .FALSE. ! Maximum DeltaP/sigma for convergence; user defined
conv_dchi2_u   = .FALSE. ! Maximum Chi^2 change for convergence; user defined
conv_chi2_u    = .FALSE. ! Minimum Chi^2 value  for convergence; user defined
conv_conf_u    = .FALSE. ! Minimum confidence level; user defined
conv_lamb_u    = .FALSE. ! Maxmimum lambda    level; user defined
!
refine_lamda_s_u = .FALSE.!MRQ lamda start
refine_lamda_u_u = .FALSE.!MRQ lamda Increase (up)
refine_lamda_d_u = .FALSE.!MRQ lamda Decrease (down)
!
refine_log = .false.
!
call reset_header      ! Clear header lines for "refine_best", "refine_new.res"
!
END SUBROUTINE refine_do_reset
!
END MODULE refine_reset
