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
CHARACTER (LEN=200)     ::      parent_results
CHARACTER (LEN=200)     ::      parent_summary
!
CHARACTER (LEN=  8) ,DIMENSION(:) , ALLOCATABLE :: pop_name ! (MAXDIMX)
!
INTEGER                 :: MAXPOP                      ! Maximum population size
INTEGER                 :: MAXDIMX                     ! Maximum parameter number
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
INTEGER                 :: lparent_results
INTEGER                 :: lparent_summary
INTEGER ,DIMENSION(:)  , ALLOCATABLE :: pop_type       !  (MAXDIMX)
INTEGER ,DIMENSION(:)  , ALLOCATABLE :: pop_lname      !  (MAXDIMX)
!
LOGICAL                              :: pop_current       = .false.
LOGICAL                              :: pop_current_trial = .false.
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
INTEGER                 :: pop_size_of  ! Bytes allocated for population
!
!      common /pop/      pop_n,pop_c,pop_dimx,pop_gen,pop_best,          &
!     &                      pop_worst,                                  &
!     &                      pop_lgenfile,pop_genfile,                   &
!     &                      pop_ltrialfile,pop_trialfile,               &
!     &                  pop_x,pop_t,pop_para,pop_xmin,pop_xmax,         &
!     &                  pop_smin,pop_smax,pop_sigma,pop_lsig,           &
!     &                      pop_sig_ad,pop_lsig_ad,                     &
!     &                      pop_ad_sigma,pop_ad_lsigma,                 &
!     &                  pop_type,pop_refine,pop_current,                &
!     &                  pop_name,pop_lname
!RBN     &                      pop_lparfile,pop_parfile,

!      common /values/ ltrial_results,trial_results,                     &
!     &                  lparent_results,parent_results,                 &
!     &                  lparent_summary,parent_summary,                 &
!     &                  child_val,trial_val,parent_val

!      common /kinder/ child,trial
!
END MODULE population
