MODULE population
!-
!
!      Variables needed to define the population
!+
!
use precision_mod
!
SAVE
PUBLIC
!
INTEGER, PARAMETER      ::      POP_INTEGER = 0
INTEGER, PARAMETER      ::      POP_REAL    = 1
!
CHARACTER (LEN=200)     ::      pop_genfile
CHARACTER (LEN=200)     ::      pop_trialfile  = ' '
!
CHARACTER (LEN=200)     ::      trial_results  = ' '
CHARACTER (LEN=200)     ::      parent_results = ' '
CHARACTER (LEN=200)     ::      parent_current = ' '
CHARACTER (LEN=200)     ::      parent_summary = ' '
!
CHARACTER (LEN= 16) ,DIMENSION(:) , ALLOCATABLE :: pop_name ! (MAXDIMX)
!
INTEGER                 :: MAXPOP                      ! Maximum population size
INTEGER                 :: MAXDIMX                     ! Maximum parameter number
INTEGER                 :: MAXBACK                     ! Maximum backup number
INTEGER                 :: pop_n
INTEGER                 :: pop_c
INTEGER                 :: pop_dimx
LOGICAL                 :: pop_dimx_new
LOGICAL                 :: pop_dimx_init = .FALSE.
INTEGER                 :: pop_gen
INTEGER                 :: lastgen
INTEGER                 :: pop_best
INTEGER                 :: pop_worst
INTEGER                 :: pop_lgenfile
INTEGER                 :: pop_ltrialfile
INTEGER                 :: pop_trialfile_ext
INTEGER                 :: ltrial_results
INTEGER                 :: trial_results_ext
INTEGER                 :: n_rvalue_i
INTEGER                 :: n_rvalue_o
INTEGER                 :: lparent_current = 0
INTEGER                 :: lparent_results = 0
INTEGER                 :: lparent_summary = 0
INTEGER                 :: pop_back_number = 0
LOGICAL                 :: pop_backup = .false.
INTEGER ,DIMENSION(:)  , ALLOCATABLE :: pop_type       !  (MAXDIMX)
INTEGER ,DIMENSION(:)  , ALLOCATABLE :: pop_lname      !  (MAXDIMX)
INTEGER ,DIMENSION(:,:), ALLOCATABLE :: pop_random     !  (nseeds:MAXPOP )
!
LOGICAL                              :: pop_not_first     = .FALSE.
LOGICAL                              :: pop_current       = .FALSE.
LOGICAL                              :: pop_current_trial = .FALSE.
LOGICAL                              :: pop_initialized   = .FALSE.
LOGICAL, PARAMETER                   :: pop_result_file_rd= .FALSE.
LOGICAL, PARAMETER                   :: pop_trial_file_wrt= .FALSE.
LOGICAL ,DIMENSION(:)  , ALLOCATABLE :: pop_refine     !  (MAXDIMX)
!
LOGICAL ,DIMENSION(:)  , ALLOCATABLE :: pop_ad_sigma   !  (MAXDIMX)
LOGICAL ,DIMENSION(:)  , ALLOCATABLE :: pop_ad_lsigma  !  (MAXDIMX)
!
REAL(kind=PREC_DP)    ,DIMENSION(:,:), ALLOCATABLE :: pop_x          !  (MAXDIMX,MAXPOP) Parent parameters
REAL(kind=PREC_DP)    ,DIMENSION(:,:), ALLOCATABLE :: pop_t          !  (MAXDIMX,MAXPOP) Trial  parameters
REAL(kind=PREC_DP)    ,DIMENSION(:)  , ALLOCATABLE :: pop_para       !  (MAXDIMX)
REAL(kind=PREC_DP)    ,DIMENSION(:)  , ALLOCATABLE :: pop_xmin       !  (MAXDIMX) Absolute lower boundary
REAL(kind=PREC_DP)    ,DIMENSION(:)  , ALLOCATABLE :: pop_xmax       !  (MAXDIMX) Absolute upper boundary
REAL(kind=PREC_DP)    ,DIMENSION(:)  , ALLOCATABLE :: pop_pmin       !  (MAXDIMX)
REAL(kind=PREC_DP)    ,DIMENSION(:)  , ALLOCATABLE :: pop_pmax       !  (MAXDIMX)
REAL(kind=PREC_DP)    ,DIMENSION(:)  , ALLOCATABLE :: pop_smin       !  (MAXDIMX) Start    lower boundary
REAL(kind=PREC_DP)    ,DIMENSION(:)  , ALLOCATABLE :: pop_smax       !  (MAXDIMX) Start    upper boundary
REAL(kind=PREC_DP)    ,DIMENSION(:)  , ALLOCATABLE :: pop_sigma      !  (MAXDIMX) Sigma for edge correction
REAL(kind=PREC_DP)    ,DIMENSION(:)  , ALLOCATABLE :: pop_lsig       !  (MAXDIMX) Sigma for local search
REAL(kind=PREC_DP)    ,DIMENSION(:)  , ALLOCATABLE :: pop_sig_ad     !  (MAXDIMX)
REAL(kind=PREC_DP)    ,DIMENSION(:)  , ALLOCATABLE :: pop_lsig_ad    !  (MAXDIMX)
!
REAL(kind=PREC_DP)    ,DIMENSION(:,:), ALLOCATABLE :: child          !  (MAXDIMX,MAXPOP) Child parameters
REAL(kind=PREC_DP)    ,DIMENSION(:,:), ALLOCATABLE :: trial          !  (MAXDIMX,MAXPOP) Trial parameters
!
REAL(kind=PREC_DP)    ,DIMENSION(:,:), ALLOCATABLE :: child_val      !  (MAXPOP)
REAL(kind=PREC_DP)    ,DIMENSION(:,:), ALLOCATABLE :: trial_val      !  (MAXPOP)
REAL(kind=PREC_DP)    ,DIMENSION(:,:), ALLOCATABLE :: parent_val     !  (MAXPOP)
!
CHARACTER (LEN=200),DIMENSION(:)  , ALLOCATABLE :: pop_back_fil ! (MAXBACK)
CHARACTER (LEN=200),DIMENSION(:)  , ALLOCATABLE :: pop_back_ext ! (MAXBACK)
CHARACTER (LEN=200),DIMENSION(:)  , ALLOCATABLE :: pop_back_trg ! (MAXBACK)
INTEGER            ,DIMENSION(:)  , ALLOCATABLE :: pop_back_fil_l ! (MAXBACK)
INTEGER            ,DIMENSION(:)  , ALLOCATABLE :: pop_back_ext_l ! (MAXBACK)
INTEGER            ,DIMENSION(:)  , ALLOCATABLE :: pop_back_trg_l ! (MAXBACK)
!
!
END MODULE population
