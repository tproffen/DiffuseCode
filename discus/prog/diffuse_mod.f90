MODULE diffuse_mod
!+
!     Contains all variables for Fourier transform
!-
USE discus_config_mod
USE precision_mod
!
SAVE
!
INTEGER , PARAMETER  :: LOT_OFF  = 1
INTEGER , PARAMETER  :: LOT_BOX  = 2
INTEGER , PARAMETER  :: LOT_ELI  = 3
!
INTEGER , PARAMETER  :: POWDER   = 0
INTEGER , PARAMETER  :: RECIPR   = 1
!
INTEGER , PARAMETER  :: EXTERNAL = 0
INTEGER , PARAMETER  :: INTERNAL = 1
!
INTEGER(KIND=PREC_INT_LARGE) , PARAMETER  :: I2PI     = 2**16
INTEGER(KIND=PREC_INT_LARGE) , PARAMETER  :: MASK     = I2PI-1
!
INTEGER , PARAMETER  :: RAD_XRAY = 1
INTEGER , PARAMETER  :: RAD_NEUT = 2
INTEGER , PARAMETER  :: RAD_ELEC = 3
INTEGER , PARAMETER  :: RAD_INTER= 0
INTEGER , PARAMETER  :: RAD_WAAS = 1
INTEGER , PARAMETER  :: RAD_DISC = 2
!
INTEGER , PARAMETER  :: FOUR_ZL  = -4  ! Zone Axis with lots
INTEGER , PARAMETER  :: FOUR_3L  = -3  ! 3D Fourier with lots
INTEGER , PARAMETER  :: FOUR_2L  = -2  ! 2D Fourier with lots
INTEGER , PARAMETER  :: FOUR_1L  = -1  ! 1D Fourier with lots
INTEGER , PARAMETER  :: FOUR_NN  =  0  ! None calculated
INTEGER , PARAMETER  :: FOUR_1D  =  1
INTEGER , PARAMETER  :: FOUR_2D  =  2
INTEGER , PARAMETER  :: FOUR_3D  =  3
INTEGER , PARAMETER  :: FOUR_ZA  =  4
INTEGER , PARAMETER  :: PATT_1D  =  5
INTEGER , PARAMETER  :: PATT_2D  =  6
INTEGER , PARAMETER  :: PATT_3D  =  7
INTEGER , PARAMETER  :: POWD_CO  =  8  ! Powder pattern complete 
INTEGER , PARAMETER  :: POWD_DY  =  9  ! Powder pattern Debye Algorithm
!
integer , parameter  :: FOUR_TURBO = 0 ! Standard Turbo mode
integer , parameter  :: FOUR_NUFFT = 1 ! RNon uniform FFT mode
!
integer                                 :: four_tech = FOUR_TURBO    ! Fourier technique
!
INTEGER                                 ::  DIF_MAXAT    ! current size of array at
INTEGER                                 ::  DIF_MAXSCAT  ! current size of array at
COMPLEX (KIND=PREC_DP    ) , DIMENSION(:, :), ALLOCATABLE  ::  cfact        ! (0:CFPKT, 1:MAXSCAT)
COMPLEX (KIND=PREC_DP    ) , DIMENSION(:, :), ALLOCATABLE  ::  cfact_pure   ! (0:CFPKT, 1:MAXSCAT)
!
COMPLEX (KIND=PREC_DP    ) , DIMENSION(:,:,:)   , ALLOCATABLE  ::  csf       ! (1:num(1),1:num(2),1:num(3))   
COMPLEX (KIND=PREC_DP    ) , DIMENSION(:,:,:)   , ALLOCATABLE  ::  tcsf      ! (1:num(1),1:num(2),1:num(3))
COMPLEX (KIND=PREC_DP    ) , DIMENSION(:,:,:)   , ALLOCATABLE  ::  acsf      ! (1:num(1),1:num(2),1:num(3))    
!
REAL    (KIND=PREC_DP    ) , DIMENSION(:,:,:)   , ALLOCATABLE  ::  rpdf         ! (1:MAXQXY)
COMPLEX (KIND=PREC_DP    ) , DIMENSION(0:MASK)             ::  cex       = (0.0D0,0.0D0)
REAL    (KIND=PREC_DP    ) , DIMENSION(:,:,:)   , ALLOCATABLE  ::  dsi          ! (1:MAXQXY)
REAL    (KIND=PREC_DP    ) , DIMENSION(:,:,:)   , ALLOCATABLE  ::  dsi3d        ! (1:MAXQXY)
!
COMPLEX (KIND=PREC_DP    ) , DIMENSION(:,:,:)   , ALLOCATABLE  ::  csf_sum      ! (1:MAXQXY)
REAL    (KIND=PREC_DP    ) , DIMENSION(:,:,:)   , ALLOCATABLE  ::  dsi_sum      ! (1:MAXQXY)
REAL    (KIND=PREC_DP    ) , DIMENSION(:, :), ALLOCATABLE  ::  xat          ! (1:NMAX, 1:3)
REAL    (KIND=PREC_DP    ), DIMENSION(1:3)                ::  xm        = 0.0D0
REAL    (KIND=PREC_DP    ), DIMENSION(1:3)                ::  win       = 0.0D0
REAL    (KIND=PREC_DP    ), DIMENSION(1:3)                ::  vin       = 0.0D0
REAL    (KIND=PREC_DP    ), DIMENSION(1:3)                ::  uin       = 0.0D0
REAL    (KIND=PREC_DP    )              ::  fave      = 0.0
REAL    (KIND=PREC_DP    )              ::  fave_sca  = 1.0
INTEGER , DIMENSION(:,:,:), ALLOCATABLE ::  istl         ! (1:MAXQXY)   Sine(theta)/lambda
INTEGER , DIMENSION(1:3)                ::  num       = 1
INTEGER                                 ::  nlots     = 1
INTEGER                                 ::  ilots     = LOT_OFF
INTEGER , DIMENSION(1:3)                ::  ls_xyz    = 5
INTEGER                                 ::  nxat      = 1
INTEGER                                 ::  four_mode = INTERNAL
LOGICAL                                 ::  lot_all   = .false.
LOGICAL                                 ::  ffour     = .false.
LOGICAL                                 ::  lperiod   = .true.
LOGICAL                                 ::  four_log  = .false.
LOGICAL                                 ::  four_was_run  = .false. ! TRUE if a fourier has been calculated
!
REAL(kind=PREC_DP)                      ::  braggmax = 0.0
REAL(kind=PREC_DP)                      ::  braggmin = 0.0
REAL(kind=PREC_DP)                      ::  diffuave = 0.0
REAL(kind=PREC_DP)                      ::  diffusig = 0.0
REAL(kind=PREC_DP)                      ::  diffumax = 0.0
REAL(kind=PREC_DP)                      ::  diffumin = 0.0
REAL(kind=PREC_DP)                      ::  ps_low   = 1.20
REAL(kind=PREC_DP)                      ::  ps_high  = 0.01
REAL(kind=PREC_DP)                      ::  zmin     = 0.0
REAL(kind=PREC_DP)                      ::  zmax     = 0.0
!
CHARACTER(LEN=4)                        ::  lambda   = 'MOA1'
INTEGER                                 ::  four_exp = 0
INTEGER , DIMENSION(1:3)                ::  diff_inc_hkl  = (/ 1, 1,  1 /)   ! Data points for <f> 
INTEGER , DIMENSION(1:3)                ::  inc      = (/ 121, 121,  1 /)
INTEGER , DIMENSION(1:6)                ::  lmn      = 0
LOGICAL                                 ::  ano      = .false.
LOGICAL                                 ::  ldbw     = .false.
LOGICAL                                 ::  lxray    = .true.
LOGICAL                                 ::  diff_lsingle  = .true.
INTEGER                                 ::  diff_radiation = RAD_XRAY
INTEGER                                 ::  diff_table     = RAD_INTER
INTEGER                                 ::  diff_power     = 4
character(len=PREC_STRING)              ::  diff_file      = ' '
character(len=PREC_STRING)              ::  diff_trust     = ' '
real(kind=prec_dp)                      ::  diff_exti = 0.0
REAL(kind=PREC_DP), DIMENSION(1:3, 1:4) ::  eck      = reshape((/ 0.0, 0.0,  0.0, &    ! (hkl, corner_number)
                                                                  5.0, 0.0,  0.0, &
                                                                  0.0, 5.0,  0.0, &
                                                                  0.0, 0.0,  0.0/),shape(eck))
REAL(kind=PREC_DP), DIMENSION(1:3, 1:3) ::  vi       = reshape((/0.05, 0.00, 0.00, &   ! (hkl, number_abs_od_top)
                                                                 0.0 , 0.05, 0.00, &
                                                                 0.00, 0.00, 0.00/),shape(vi))
REAL(kind=PREC_DP), DIMENSION(1:3, 1:4) ::  diff_eck_hkl = reshape((/ 0.0, 0.0,  0.0, &
                                                                  5.0, 0.0,  0.0, &
                                                                  0.0, 5.0,  0.0, &
                                                                  0.0, 0.0,  0.0/),shape(diff_eck_hkl))
REAL(kind=PREC_DP), DIMENSION(1:3, 1:3) ::  diff_vi_hkl  = reshape((/1.00, 0.00, 0.00, &
                                                                 0.0 , 1.00, 0.00, &
                                                                 0.00, 0.00, 0.00/),shape(diff_vi_hkl))
REAL(kind=PREC_DP), DIMENSION(1:3, 1:3) ::  off_shift= 0.00              ! (hkl, abs_ord_top)
REAL(kind=PREC_DP)                      ::  renergy  = 17.480782
REAL(kind=PREC_DP)                      ::  rlambda  =  0.709260
LOGICAL                                 ::  l_energy = .false.
REAL(kind=PREC_DP), DIMENSION(1:3,1:4)             ::  diff_eck_u  = 0.0   ! User supplied corners
REAL(kind=PREC_DP), DIMENSION(1:3, 1:3)           ::  diff_vi_u   = 0.0
INTEGER , DIMENSION(1:3)                ::  diff_inc_u  = 1
!
LOGICAL                                 ::  diff_l_friedel = .FALSE.     ! Use Fridels law to reduce calculation time
LOGICAL, DIMENSION(0:3)                 ::  diff_l_even    = .FALSE.     ! User inc is even along dimension
LOGICAL, DIMENSION(1:3)                 ::  diff_l_all     = .FALSE.     ! User inc is even along all dimension
INTEGER                                 ::  diff_idim      = 0           ! Which dimension has been cut in half
!
LOGICAL                                 ::  l_zone       = .false.
REAL(kind=PREC_DP)   , DIMENSION(1:3)   ::  zone_uvw     = (/ 0.0, 0.0, 1.0/) 
REAL(kind=PREC_DP)   , DIMENSION(1:3)   ::  zone_ewald   = 0.0
REAL(kind=PREC_DP)                      ::  zone_res     = 0.0
REAL(kind=PREC_DP)                      ::  zone_delta_d = 0.015
!
REAL(KIND=PREC_DP), DIMENSION(4,3)      ::  diff_res             ! Resolution sigma and vectors
REAL(KIND=PREC_DP), DIMENSION(3,3)      ::  diff_tr              ! Resolution transformation matrix
real(kind=PREC_DP), dimension(:,:,:,:), allocatable :: four_dbw  ! Debye-Waller map
!
INTEGER                                 ::  four_last = FOUR_NN  ! No Fourier calculated yet
!
INTEGER, PARAMETER :: FOUR_ACCUM_INIT     = -1
INTEGER, PARAMETER :: FOUR_ACCUM_SINGLE   =  0
INTEGER, PARAMETER :: FOUR_ACCUM_ACCUM    =  1
INTEGER, PARAMETER :: FOUR_ACCUM_FINISHED =  2
INTEGER                                 :: four_accum = 0        ! Run a single Fourier (-1==init, 0==single, 1==add, 2==finished)
LOGICAL                                 :: four_symm  = .FALSE.  ! Apply reciprocal space symmetry averaging
logical                                 :: four_friedel = .true. ! Reduce computation time via Friedel, if possible
!
! Filter parameter
!
integer, parameter :: FOUR_FILTER_OFF     = 0
integer, parameter :: FOUR_FILTER_LANCZOS = 1
integer            :: four_filter
integer            :: four_nscale
integer            :: four_width
real(kind=PREC_DP) :: four_rscale
real(kind=PREC_DP) :: four_damp
!
END MODULE diffuse_mod
