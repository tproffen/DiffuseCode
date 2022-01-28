MODULE phases_mod
!+
!     variables needed for multiple phases
!-
USE discus_config_mod
USE precision_mod
!
SAVE
!
INTEGER :: PHA_MAXPHA  = 1
INTEGER :: PHA_MAXPTS  = 1
INTEGER :: PHA_MAXSCAT = 1
!
! Section for multiple phases
!
LOGICAL                                         :: pha_multi = .FALSE.   ! Multiple phases yes / no
INTEGER                                         :: pha_n     = 1         ! number of phases
INTEGER                                         :: pha_curr  = 1         ! curent    phase
INTEGER           , DIMENSION(:)  , ALLOCATABLE :: pha_nscat   ! No of atom types for each phase (      i)
INTEGER           , DIMENSION(:)  , ALLOCATABLE :: pha_calc    ! Calc mode Comp/Debye each phase (      i)
REAL(KIND=PREC_DP), DIMENSION(:)  , ALLOCATABLE :: pha_frac    ! weight fraction User intent     (      i)
REAL(KIND=PREC_DP), DIMENSION(:)  , ALLOCATABLE :: pha_weight  ! Weight          for each phase   (      i)
REAL(KIND=PREC_DP), DIMENSION(:)  , ALLOCATABLE :: pha_scale   ! Scale temp>frac for each phase   (      i)
REAL(KIND=PREC_DP), DIMENSION(:)  , ALLOCATABLE :: pha_nreal   ! Number real atoms at     phase   (      i)
REAL(KIND=PREC_DP), DIMENSION(:)  , ALLOCATABLE :: pha_ncreal  ! Number real atoms /unit call at phase   (      i)
REAL(KIND=PREC_DP), DIMENSION(:,:), ALLOCATABLE :: pha_powder  ! Powder pattern for each phase   (q,    i)
REAL(KIND=PREC_DP), DIMENSION(:,:,:), ALLOCATABLE :: pha_form    ! Form factors   for each phase   (q,is,i)
REAL(KIND=PREC_DP), DIMENSION(:,:), ALLOCATABLE :: pha_adp     ! B-values ADP   for each phase   (is,   i)
REAL(KIND=PREC_DP), DIMENSION(:,:), ALLOCATABLE :: pha_occ     ! Occupancies    for each phase   (is ,  i)
INTEGER           , DIMENSION(:,:), ALLOCATABLE :: pha_niscat  ! Atom numbers of iscat at phasof iscat at phase
!
!
END MODULE phases_mod
