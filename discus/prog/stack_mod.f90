MODULE stack_mod
!+
!
!     Variables needed for the generalized stacking faults
!-
USE config_mod
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
INTEGER            ::  st_layer_increment  = 5
!
!INTEGER       BL_ST_MAXTYPE
!PARAMETER    (BL_ST_MAXTYPE   = ST_MAXTYPE*ST_MAXTYPE)
!INTEGER       BL_ST_MAXTYPE3
!PARAMETER    (BL_ST_MAXTYPE3  = ST_MAXTYPE*ST_MAXTYPE*3)
!INTEGER       BL_ST_MAXLAYER3
!PARAMETER    (BL_ST_MAXLAYER3 = ST_MAXLAYER*3)
!
INTEGER ::  ST_MAXQXY    = 1
INTEGER ::  ST_MAXLAYER  = 1
INTEGER ::  ST_MAXTYPE   = 1
!
CHARACTER (LEN=200), DIMENSION(:)    , ALLOCATABLE :: st_layer        ! (  ST_MAXTYPE)
CHARACTER (LEN=200), DIMENSION(:)    , ALLOCATABLE :: st_layer_c      ! (  ST_MAXTYPE)
INTEGER            , DIMENSION(:)    , ALLOCATABLE :: st_llayer       ! (  ST_MAXTYPE)
INTEGER            , DIMENSION(:)    , ALLOCATABLE :: st_number       ! (  ST_MAXTYPE)
INTEGER            , DIMENSION(:)    , ALLOCATABLE :: st_ndisp        ! (  ST_MAXTYPE)
INTEGER            , DIMENSION(:)    , ALLOCATABLE :: st_chem         ! (  ST_MAXTYPE)
REAL               , DIMENSION(:,:)  , ALLOCATABLE :: st_disp         ! (3,ST_MAXTYPE)
REAL               , DIMENSION(:,:)  , ALLOCATABLE :: st_corr         ! (  ST_MAXTYPE, ST_MAXTYPE)
REAL               , DIMENSION(:,:,:), ALLOCATABLE :: st_sigma        ! (  ST_MAXTYPE, ST_MAXTYPE,3)
REAL               , DIMENSION(:,:,:), ALLOCATABLE :: st_trans        ! (  ST_MAXTYPE, ST_MAXTYPE,3)
!
INTEGER            , DIMENSION(:)    , ALLOCATABLE :: st_type         ! (  ST_MAXLAYER)
LOGICAL            , DIMENSION(:)    , ALLOCATABLE :: st_internal     ! (  ST_MAXTYPE )
REAL               , DIMENSION(:,:)  , ALLOCATABLE :: st_origin       ! (3,ST_MAXLAYER)
REAL               , DIMENSION(:)    , ALLOCATABLE :: st_rot_ang_no   ! (  ST_MAXLAYER)
REAL               , DIMENSION(:)    , ALLOCATABLE :: st_rot_ang_m1   ! (  ST_MAXLAYER)
REAL               , DIMENSION(:)    , ALLOCATABLE :: st_rot_ang_m2   ! (  ST_MAXLAYER)
!
COMPLEX            , DIMENSION(:)    , ALLOCATABLE :: st_csf          ! (ST_MAXQXY)
!
CHARACTER (LEN=200)                        :: st_infile        = ' '
!
INTEGER                                    :: st_distr         = ST_DIST_MATRIX
INTEGER                                    :: st_infile_l      = 1
INTEGER                                    :: st_nlayer        = 0
INTEGER                                    :: st_ntypes        = 0
INTEGER                                    :: st_nchem         = 0
LOGICAL                                    :: st_mod_sta       = .false.
LOGICAL                                    :: st_tra_aver      = .false.
LOGICAL                                    :: st_rot_mode      = .false.
LOGICAL                                    :: st_rot_status    = .false.
LOGICAL                                    :: st_rot_no_lspace = .true.
LOGICAL                                    :: st_rot_m1_lspace = .true.
LOGICAL                                    :: st_rot_m2_lspace = .true.
REAL                                       :: st_aver          = 0.0
REAL                                       :: st_prob          = 0.0
REAL   , DIMENSION(3,3)                    :: st_mod           = 0.0
REAL   , DIMENSION(3,3)                    :: st_inv           = 0.0
REAL   , DIMENSION(3)                      :: st_t_aver        = 0.0
REAL   , DIMENSION(3)                      :: st_off           = 0.0
REAL   , DIMENSION(3)                      :: st_sigma_off     = 0.0
REAL   , DIMENSION(3)                      :: st_rot_no        = 0.0
REAL   , DIMENSION(3)                      :: st_rot_m1        = 0.0
REAL   , DIMENSION(3)                      :: st_rot_m2        = 0.0
REAL                                       :: st_rot_si_no     = 0.0
REAL                                       :: st_rot_si_m1     = 0.0
REAL                                       :: st_rot_si_m2     = 0.0
!
INTEGER                                    :: st_size_of       = 0
!
END MODULE stack_mod
