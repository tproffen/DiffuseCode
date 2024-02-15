MODULE refine_control_mod
!
use precision_mod
IMPLICIT NONE
!
character(len=PREC_STRING)        :: ref_run_u      = ' '     ! User command on run line
INTEGER                           :: refine_cycles  = 1       ! Maximum number of cycles defined by user
LOGICAL                           :: conv_status    = .TRUE.  ! Apply convergence criteria
REAL(kind=PREC_DP)                :: conv_dp_sig    = 0.005D0 ! Maximum DeltaP/sigma for convergence
REAL(kind=PREC_DP)                :: conv_dchi2     = 0.5D0   ! Maximum Chi^2 change for convergence
REAL(kind=PREC_DP)                :: conv_chi2      = 0.5D0   ! Minimum Chi^2 value  for convergence
REAL(kind=PREC_DP)                :: conv_conf      = 0.01D0  ! Minimum confidence level
logical                           :: conv_dp_sig_u  = .FALSE. ! Maximum DeltaP/sigma for convergence; user defined
logical                           :: conv_dchi2_u   = .FALSE. ! Maximum Chi^2 change for convergence; user defined
logical                           :: conv_chi2_u    = .FALSE. ! Minimum Chi^2 value  for convergence; user defined
logical                           :: conv_conf_u    = .FALSE. ! Minimum confidence level; user defined
LOGICAL                           :: lconvergence   = .FALSE. ! Convergence has been reached
logical           , dimension(3)  :: lconv          = .false. ! Which convergence criterium has been met
REAL(kind=PREC_DP)                :: refine_lamda_s = 0.001D0 ! MRQ lamda start
REAL(kind=PREC_DP)                :: refine_lamda_u = 4.000D0 ! MRQ lamda Increase (up)
REAL(kind=PREC_DP)                :: refine_lamda_d = 0.500D0 ! MRQ lamda Decrease (down)
logical                           :: refine_lamda_s_u = .FALSE.!MRQ lamda start
logical                           :: refine_lamda_u_u = .FALSE.!MRQ lamda Increase (up)
logical                           :: refine_lamda_d_u = .FALSE.!MRQ lamda Decrease (down)
LOGICAL                           :: refine_init    = .TRUE.  ! Initialize 'run'
!
END MODULE refine_control_mod
