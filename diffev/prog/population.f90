MODULE population
!-
!
!      Variables needed to define the population
!+
!
SAVE
PUBLIC
!
INTEGER, PARAMETER      ::      POP_INTEGER = 0
INTEGER, PARAMETER      ::      POP_REAL    = 1
!
CHARACTER (LEN=200)     ::      pop_genfile
CHARACTER (LEN=200)     ::      pop_trialfile
!
CHARACTER (LEN=200)     ::      trial_results
CHARACTER (LEN=200)     ::      parent_results = ' '
CHARACTER (LEN=200)     ::      parent_current = ' '
CHARACTER (LEN=200)     ::      parent_summary = ' '
!
CHARACTER (LEN=  8) ,DIMENSION(:) , ALLOCATABLE :: pop_name ! (MAXDIMX)
!
INTEGER                 :: MAXPOP                      ! Maximum population size
INTEGER                 :: MAXDIMX                     ! Maximum parameter number
INTEGER                 :: MAXBACK                     ! Maximum backup number
INTEGER                 :: pop_n
INTEGER                 :: pop_c
INTEGER                 :: pop_dimx
LOGICAL                 :: pop_dimx_new
INTEGER                 :: pop_gen
INTEGER                 :: pop_best
INTEGER                 :: pop_worst
INTEGER                 :: pop_lgenfile
INTEGER                 :: pop_ltrialfile
INTEGER                 :: pop_trialfile_ext
INTEGER                 :: ltrial_results
INTEGER                 :: trial_results_ext
INTEGER                 :: lparent_current = 0
INTEGER                 :: lparent_results = 0
INTEGER                 :: lparent_summary = 0
INTEGER                 :: pop_back_number = 0
LOGICAL                 :: pop_backup = .false.
INTEGER ,DIMENSION(:)  , ALLOCATABLE :: pop_type       !  (MAXDIMX)
INTEGER ,DIMENSION(:)  , ALLOCATABLE :: pop_lname      !  (MAXDIMX)
INTEGER ,DIMENSION(:,:), ALLOCATABLE :: pop_random     !  (nseeds:MAXPOP )
!
LOGICAL                              :: pop_current       = .false.
LOGICAL                              :: pop_current_trial = .false.
LOGICAL                              :: pop_result_file_rd= .true.
LOGICAL                              :: pop_trial_file_wrt= .true.
LOGICAL ,DIMENSION(:)  , ALLOCATABLE :: pop_refine     !  (MAXDIMX)
!
LOGICAL ,DIMENSION(:)  , ALLOCATABLE :: pop_ad_sigma   !  (MAXDIMX)
LOGICAL ,DIMENSION(:)  , ALLOCATABLE :: pop_ad_lsigma  !  (MAXDIMX)
!
REAL    ,DIMENSION(:,:), ALLOCATABLE :: pop_x          !  (MAXDIMX,MAXPOP)
REAL    ,DIMENSION(:,:), ALLOCATABLE :: pop_t          !  (MAXDIMX,MAXPOP)
REAL    ,DIMENSION(:)  , ALLOCATABLE :: pop_para       !  (MAXDIMX)
REAL    ,DIMENSION(:)  , ALLOCATABLE :: pop_xmin       !  (MAXDIMX)
REAL    ,DIMENSION(:)  , ALLOCATABLE :: pop_xmax       !  (MAXDIMX)
REAL    ,DIMENSION(:)  , ALLOCATABLE :: pop_pmin       !  (MAXDIMX)
REAL    ,DIMENSION(:)  , ALLOCATABLE :: pop_pmax       !  (MAXDIMX)
REAL    ,DIMENSION(:)  , ALLOCATABLE :: pop_smin       !  (MAXDIMX)
REAL    ,DIMENSION(:)  , ALLOCATABLE :: pop_smax       !  (MAXDIMX)
REAL    ,DIMENSION(:)  , ALLOCATABLE :: pop_sigma      !  (MAXDIMX)
REAL    ,DIMENSION(:)  , ALLOCATABLE :: pop_lsig       !  (MAXDIMX)
REAL    ,DIMENSION(:)  , ALLOCATABLE :: pop_sig_ad     !  (MAXDIMX)
REAL    ,DIMENSION(:)  , ALLOCATABLE :: pop_lsig_ad    !  (MAXDIMX)
!
REAL    ,DIMENSION(:,:), ALLOCATABLE :: child          !  (MAXDIMX,MAXPOP)
REAL    ,DIMENSION(:,:), ALLOCATABLE :: trial          !  (MAXDIMX,MAXPOP)
!
REAL    ,DIMENSION(:)  , ALLOCATABLE :: child_val      !  (MAXPOP)
REAL    ,DIMENSION(:)  , ALLOCATABLE :: trial_val      !  (MAXPOP)
REAL    ,DIMENSION(:)  , ALLOCATABLE :: parent_val     !  (MAXPOP)
!
CHARACTER (LEN=200),DIMENSION(:)  , ALLOCATABLE :: pop_back_fil ! (MAXBACK)
CHARACTER (LEN=200),DIMENSION(:)  , ALLOCATABLE :: pop_back_ext ! (MAXBACK)
CHARACTER (LEN=200),DIMENSION(:)  , ALLOCATABLE :: pop_back_trg ! (MAXBACK)
INTEGER            ,DIMENSION(:)  , ALLOCATABLE :: pop_back_fil_l ! (MAXBACK)
INTEGER            ,DIMENSION(:)  , ALLOCATABLE :: pop_back_ext_l ! (MAXBACK)
INTEGER            ,DIMENSION(:)  , ALLOCATABLE :: pop_back_trg_l ! (MAXBACK)
!
INTEGER                 :: pop_size_of  ! Bytes allocated for population
!
END MODULE population
