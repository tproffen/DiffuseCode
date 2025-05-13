MODULE stack_mod
!+
!
!     Variables needed for the generalized stacking faults
!-
USE discus_config_mod
use precision_mod
!
SAVE
!
INTEGER, PARAMETER ::       ST_F_ALL = 0
INTEGER, PARAMETER ::       ST_F_AVE = 1
INTEGER, PARAMETER ::       ST_F_SUB = 2
!     
INTEGER, PARAMETER ::       ST_DIST_MATRIX = 0
INTEGER, PARAMETER ::       ST_DIST_FILE   = 1
INTEGER, PARAMETER ::       ST_DIST_LIST   = 2
!     
integer, parameter ::       ST_ATOM_OFF    =  0
integer, parameter ::       ST_ATOM_ON     =  1
integer, parameter ::       ST_ATOM_STRICT = -1
!
INTEGER            ::  st_layer_increment  = 5
!
!INTEGER       BL_ST_MAXTYPE
!PARAMETER    (BL_ST_MAXTYPE   = ST_MAXTYPE*ST_MAXTYPE)
!INTEGER       BL_ST_MAXTYPE3
!PARAMETER    (BL_ST_MAXTYPE3  = ST_MAXTYPE*ST_MAXTYPE*3)
!INTEGER       BL_ST_MAXLAYER3
!PARAMETER    (BL_ST_MAXLAYER3 = ST_MAXLAYER*3)
!
INTEGER, dimension(3) ::  ST_MAXQXY    = 1
INTEGER ::  ST_MAXLAYER  = 1
INTEGER ::  ST_MAXTYPE   = 1
!
CHARACTER (LEN=1024), DIMENSION(:)    , ALLOCATABLE :: st_layer        ! (  ST_MAXTYPE)
CHARACTER (LEN=1024), DIMENSION(:)    , ALLOCATABLE :: st_layer_c      ! (  ST_MAXTYPE)
INTEGER            , DIMENSION(:)    , ALLOCATABLE :: st_llayer       ! (  ST_MAXTYPE)
INTEGER            , DIMENSION(:)    , ALLOCATABLE :: st_number       ! (  ST_MAXTYPE)
INTEGER            , DIMENSION(:)    , ALLOCATABLE :: st_ndisp        ! (  ST_MAXTYPE)
INTEGER            , DIMENSION(:)    , ALLOCATABLE :: st_chem         ! (  ST_MAXTYPE)
REAL(kind=PREC_DP) , DIMENSION(:,:)  , ALLOCATABLE :: st_disp         ! (3,ST_MAXTYPE)
REAL(kind=PREC_DP) , DIMENSION(:,:)  , ALLOCATABLE :: st_corr         ! (  ST_MAXTYPE, ST_MAXTYPE)
REAL(kind=PREC_DP) , DIMENSION(:,:,:), ALLOCATABLE :: st_sigma        ! (  ST_MAXTYPE, ST_MAXTYPE,3)
REAL(kind=PREC_DP) , DIMENSION(:,:,:), ALLOCATABLE :: st_trans        ! (  ST_MAXTYPE, ST_MAXTYPE,3)
!
INTEGER            , DIMENSION(:)    , ALLOCATABLE :: st_type         ! (  ST_MAXLAYER)
LOGICAL            , DIMENSION(:)    , ALLOCATABLE :: st_internal     ! (  ST_MAXTYPE )
REAL(kind=PREC_DP) , DIMENSION(:,:)  , ALLOCATABLE :: st_origin       ! (3,ST_MAXLAYER)
REAL(kind=PREC_DP) , DIMENSION(:)    , ALLOCATABLE :: st_rot_ang_no   ! (  ST_MAXLAYER)
REAL(kind=PREC_DP) , DIMENSION(:)    , ALLOCATABLE :: st_rot_ang_m1   ! (  ST_MAXLAYER)
REAL(kind=PREC_DP) , DIMENSION(:)    , ALLOCATABLE :: st_rot_ang_m2   ! (  ST_MAXLAYER)
!
! COMPLEX (KIND=KIND(0.0D0)), DIMENSION(:)    , ALLOCATABLE :: st_csf          ! (ST_MAXQXY)                                      ! Who is st_csf? Should it be 2D or 3D? ( Neder's original code )
COMPLEX (KIND=KIND(0.0D0)), DIMENSION(:,:,:), ALLOCATABLE :: st_csf            ! (ST_MAXQXY, ST_MAXTYPE, ST_MAXTYPE)              ! My code
!
CHARACTER (LEN=1024)                       :: st_infile        = ' '
!
LOGICAL                                    :: st_new_form      = .TRUE.    ! New formfactors need to be copied
INTEGER                                    :: st_distr         = ST_DIST_MATRIX
INTEGER                                    :: st_infile_l      = 1
INTEGER                                    :: st_nlayer        = 0
INTEGER                                    :: st_ntypes        = 0
INTEGER                                    :: st_nchem         = 0
INTEGER                                    :: st_first         = 0
integer                                    :: st_ncunit        = 1  ! Number of layers per unit translation
LOGICAL                                    :: st_mod_sta       = .false.
integer                                    :: st_mod_atom      = 0
LOGICAL                                    :: st_tra_aver      = .false.
LOGICAL                                    :: st_rot_mode      = .false.
LOGICAL                                    :: st_rot_status    = .false.
LOGICAL                                    :: st_rot_no_lspace = .true.
LOGICAL                                    :: st_rot_m1_lspace = .true.
LOGICAL                                    :: st_rot_m2_lspace = .true.
LOGICAL                                    :: st_cr_magnetic  = .FALSE.   ! Crystal is magnetic YES / NO
logical                                    :: st_is_anis      = .FALSE.   ! Crystal has anisotropic ADPs 
logical                                    :: st_is_asym      = .TRUE.    ! Crystal is an asymmetric unit
logical                                    :: st_is_homo      = .TRUE.    ! Crystal is homogeneous
logical                                    :: st_is_stack     = .false.   ! Crystal is build by stacking faults
REAL(kind=PREC_DP)                         :: st_aver          = 0.0D0
REAL(kind=PREC_DP)                         :: st_prob          = 0.0D0
REAL(kind=PREC_DP)   , DIMENSION(3,3)      :: st_mod           = 0.0D0
REAL(kind=PREC_DP)   , DIMENSION(3,3)      :: st_inv           = 0.0D0
REAL(kind=PREC_DP)   , DIMENSION(3)        :: st_t_aver        = 0.0D0
REAL(kind=PREC_DP)   , DIMENSION(3)        :: st_off           = 0.0D0
REAL(kind=PREC_DP)   , DIMENSION(3)        :: st_sigma_off     = 0.0D0
REAL(kind=PREC_DP)   , DIMENSION(3)        :: st_rot_no        = 0.0D0
REAL(kind=PREC_DP)   , DIMENSION(3)        :: st_rot_m1        = 0.0D0
REAL(kind=PREC_DP)   , DIMENSION(3)        :: st_rot_m2        = 0.0D0
REAL(kind=PREC_DP)                         :: st_rot_si_no     = 0.0D0
REAL(kind=PREC_DP)                         :: st_rot_si_m1     = 0.0D0
REAL(kind=PREC_DP)                         :: st_rot_si_m2     = 0.0D0
!
!
END MODULE stack_mod
