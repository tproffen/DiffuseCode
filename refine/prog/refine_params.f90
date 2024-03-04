MODULE refine_params_mod
!
!  Variables related to the refined parameters 
!
use precision_mod
!
IMPLICIT NONE
!
INTEGER                                        :: REF_MAXPARAM    = 0 ! maximum number of parameters
INTEGER                                        :: REF_MAXPARAM_FIX= 0 ! maximum number of fixed parameters
INTEGER, parameter                             :: REF_MAXPARAM_SPC=  9 ! maximum number of special parameters
INTEGER                                        :: refine_par_n    = 0 ! number of parameters
INTEGER                                        :: refine_fix_n    = 0 ! number of fixed parameters
LOGICAL                                        :: refine_autoconstr =.TRUE. ! Do automatic constraints
LOGICAL          , DIMENSION(0:3)              :: refine_fwhm       = .FALSE. ! u,v,w, are refined=TRUE or fixed=FALSE
INTEGER          , DIMENSION(3)                :: refine_fwhm_ind = HUGE(0) ! u,v,w, are at these locations in _p/ _f
LOGICAL          , DIMENSION(0:3)              :: refine_eta        = .FALSE. ! eta    are refined=TRUE or fixed=FALSE
INTEGER          , DIMENSION(3)                :: refine_eta_ind = HUGE(0)  ! eta    are at these locations in _p/ _f
CHARACTER(LEN=16), DIMENSION(:)  , ALLOCATABLE :: refine_params       ! parameter names
CHARACTER(LEN=16), DIMENSION(:)  , ALLOCATABLE :: refine_fixed        ! parameter names, fixed values
character(len=16), dimension(REF_MAXPARAM_SPC) :: refine_spc_name     ! parameters with special meaning
real(kind=PREC_DP),dimension(REF_MAXPARAM_SPC) :: refine_spc_delta    ! Fixed shifts for parameters with special meaning
real(kind=PREC_DP),dimension(REF_MAXPARAM_SPC) :: refine_spc_shift    ! Fixed relative shifts for parameters with special meaning
integer           ,dimension(REF_MAXPARAM_SPC) :: refine_spc_nderiv   ! Fixed relative shifts for parameters with special meaning
REAL(kind=PREC_DP),DIMENSION(:,:), ALLOCATABLE :: refine_range        ! allowed parameter range
REAL(kind=PREC_DP),DIMENSION(:)  , ALLOCATABLE :: refine_p            ! Current parameter value
REAL(kind=PREC_DP),DIMENSION(:)  , ALLOCATABLE :: refine_f            ! Current fixed parameter value
REAL(kind=PREC_DP),DIMENSION(:)  , ALLOCATABLE :: refine_dp           ! Current parameter sigma
REAL(kind=PREC_DP),DIMENSION(:,:), ALLOCATABLE :: refine_cl           ! Correlation matrix
REAL(kind=PREC_DP),DIMENSION(:,:), ALLOCATABLE :: refine_alpha        ! temporary Correlation matrix
REAL(kind=PREC_DP),DIMENSION(:  ), ALLOCATABLE :: refine_beta         ! temporary parameter shift
REAL(kind=PREC_DP),DIMENSION(:)  , ALLOCATABLE :: refine_shift        ! P*shift gives shift to calc derivative
INTEGER          , DIMENSION(:)  , ALLOCATABLE :: refine_nderiv       ! Number of p+n*DELTA to calc derivative
INTEGER          , DIMENSION(:)  , ALLOCATABLE :: refine_kderiv       ! KUPLOT data set that has derivative
!
data refine_spc_name  /   &
     'P_eta_l         ',  &
     'P_eta_q         ',  &
     'P_biso          ',  &
     'P_lat           ',  &
     'P_dia           ',  &
     'P_eta           ',  &
     'P_u             ',  &
     'P_v             ',  &
     'P_w             '   &
    /
data refine_spc_delta /   &      ! Default shifts if Parameter value is zero
     1.0D-5,   &   ! P_eta_l
     1.0D-6,   &   ! P_eta_q
     1.0D-3,   &   ! P_biso
     1.0D-3,   &   ! P_lat
     1.0D-1,   &   ! P_dia
     1.0D-3,   &   ! P_eta
     1.0D-6,   &   ! P_u
     1.0D-5,   &   ! P_v
     1.0D-5 /      ! P_w
data refine_spc_shift /   &      ! Default relative shifts if "shift:" is omitted
     3.0D-4,   &   ! P_eta_l
     3.0D-5,   &   ! P_eta_q
     3.0D-3,   &   ! P_biso
     3.0D-3,   &   ! P_lat
     1.0D-1,   &   ! P_dia
     3.0D-3,   &   ! P_eta
     3.0D-3,   &   ! P_u
     3.0D-3,   &   ! P_v
     3.0D-3 /      ! P_w
data refine_spc_nderiv /   &      ! Default derivative points if "points:" is omitted
     3     ,   &   ! P_eta_l
     3     ,   &   ! P_eta_q
     3     ,   &   ! P_biso
     3     ,   &   ! P_lat
     5     ,   &   ! P_dia
     3     ,   &   ! P_eta
     5     ,   &   ! P_u
     5     ,   &   ! P_v
     5      /      ! P_w
!
END MODULE refine_params_mod
