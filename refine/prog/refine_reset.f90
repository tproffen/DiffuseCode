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
USE refine_params_mod
USE do_variable_mod
!
IMPLICIT NONE
!
CHARACTER (LEN=1024)                   :: zeile
!INTEGER            , PARAMETER         :: MAXW=2
!INTEGER                                :: ianz
!CHARACTER(LEN=1024), DIMENSION(1:MAXW) :: cpara
!INTEGER            , DIMENSION(1:MAXW) :: lpara
!
!LOGICAL, PARAMETER   :: is_refine = .TRUE.
INTEGER              :: lcomm
!INTEGER              :: i
!
! Remove all parameter names from the variable entry
!
zeile ='default'
lcomm = 7
!
CALL refine_do_allocate_appl(zeile,lcomm)
CALL refine_initarrays
!
refine_par_n    = 0 ! number of parameters
refine_fix_n    = 0 
ref_dim(1)      = 1
ref_dim(2)      = 1
!
! refine_fit_erg
!
refine_chisqr   = -1.0  ! Sum of dev^2 in previous cycle
refine_conf     = -1.0  ! sum of dev^2
refine_lamda    = -1.0  ! Final convergence criterion
refine_rval     = -1.0  ! weighted R-value
refine_rexp     = -1.0  ! unweighted R-value
!
!refine_contol_mod
!
refine_cycles   = 1     ! Maximum number of cycles defined by user
conv_dp_sig     = 0.005 ! Maximum DeltaP/sigma for convergence
conv_dchi2      = 0.5   ! Maximum Chi^2 change for convergence
conv_chi2       = 0.5   ! Minimum Chi^2 value  for convergence
conv_conf       = 0.01  ! Minimum confidence level
lconvergence    = .FALSE. ! Convergence has been reached
refine_lamda_s  = 0.001   ! MRQ lamda start
refine_lamda_u  = 4.000   ! MRQ lamda Increase (up)
refine_lamda_d  = 0.500   ! MRQ lamda Decrease (down)
!
! refine_data_mod
!
ref_load    = ' '  ! Load string data
ref_csigma  = ' '  ! Load string Sigma's
ref_kload   = 0    ! Data set within KUPLOT
ref_ksigma  = 0    ! Sigma set within KUPLOT
ref_kupl    = 0    ! Data set within KUPLOT that needs to be kept

!
END SUBROUTINE refine_do_reset
!
END MODULE refine_reset
