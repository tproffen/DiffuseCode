MODULE refine_control_mod
!
IMPLICIT NONE
!
INTEGER :: refine_cycles = 1     ! Maximum number of cycles defined by user
REAL    :: conv_dp_sig   = 0.005 ! Maximum DeltaP/sigma for convergence
REAL    :: conv_dchi2    = 0.5   ! Maximum Chi^2 change for convergence
REAL    :: conv_chi2     = 0.5   ! Minimum Chi^2 value  for convergence
REAL    :: conv_conf     = 0.01  ! Minimum confidence level
LOGICAL :: lconvergence  = .FALSE. ! Convergence has been reached
!
END MODULE refine_control_mod
