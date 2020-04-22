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
INTEGER                                 ::  DIF_MAXAT    ! current size of array at
INTEGER                                 ::  DIF_MAXSCAT  ! current size of array at
COMPLEX (KIND=KIND(0.0D0)) , DIMENSION(:, :), ALLOCATABLE  ::  cfact        ! (0:CFPKT, 1:MAXSCAT)
COMPLEX (KIND=KIND(0.0D0)) , DIMENSION(:, :), ALLOCATABLE  ::  cfact_pure   ! (0:CFPKT, 1:MAXSCAT)
COMPLEX (KIND=KIND(0.0D0)) , DIMENSION(:)   , ALLOCATABLE  ::  csf          ! (1:MAXQXY)
COMPLEX (KIND=KIND(0.0D0)) , DIMENSION(:)   , ALLOCATABLE  ::  tcsf         ! (1:MAXQXY)
COMPLEX (KIND=KIND(0.0D0)) , DIMENSION(:)   , ALLOCATABLE  ::  acsf         ! (1:MAXQXY)
REAL    (KIND=KIND(0.0D0)) , DIMENSION(:)   , ALLOCATABLE  ::  rpdf         ! (1:MAXQXY)
COMPLEX (KIND=KIND(0.0D0)) , DIMENSION(0:MASK)             ::  cex       = (0.0D0,0.0D0)
REAL    (KIND=KIND(0.0D0)) , DIMENSION(:)   , ALLOCATABLE  ::  dsi          ! (1:MAXQXY)
REAL    , DIMENSION(:, :), ALLOCATABLE  ::  xat          ! (1:NMAX, 1:3)
REAL    (KIND=KIND(0.0D0)), DIMENSION(1:3)                ::  xm        = 0.0D0
REAL    (KIND=KIND(0.0D0)), DIMENSION(1:3)                ::  win       = 0.0D0
REAL    (KIND=KIND(0.0D0)), DIMENSION(1:3)                ::  vin       = 0.0D0
REAL    (KIND=KIND(0.0D0)), DIMENSION(1:3)                ::  uin       = 0.0D0
REAL                                    ::  fave      = 0.0
INTEGER , DIMENSION(:)   , ALLOCATABLE  ::  istl         ! (1:MAXQXY)
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
REAL                                    ::  braggmax = 0.0
REAL                                    ::  braggmin = 0.0
REAL                                    ::  diffuave = 0.0
REAL                                    ::  diffusig = 0.0
REAL                                    ::  diffumax = 0.0
REAL                                    ::  diffumin = 0.0
REAL                                    ::  ps_low   = 1.20
REAL                                    ::  ps_high  = 0.01
REAL                                    ::  zmin     = 0.0
REAL                                    ::  zmax     = 0.0
!
CHARACTER(LEN=4)                        ::  lambda   = 'MOA1'
INTEGER                                 ::  four_exp = 0
INTEGER , DIMENSION(1:3)                ::  inc      = (/ 121, 121,  1 /)
INTEGER , DIMENSION(1:6)                ::  lmn      = 0
LOGICAL                                 ::  ano      = .false.
LOGICAL                                 ::  ldbw     = .false.
LOGICAL                                 ::  lxray    = .true.
INTEGER                                 ::  diff_radiation = RAD_XRAY
INTEGER                                 ::  diff_power     = 4
REAL    , DIMENSION(1:3, 1:4)           ::  eck      = reshape((/ 0.0, 0.0,  0.0, &
                                                                  5.0, 0.0,  0.0, &
                                                                  0.0, 5.0,  0.0, &
                                                                  0.0, 0.0,  0.0/),shape(eck))
REAL    , DIMENSION(1:3, 1:3)           ::  vi       = reshape((/0.05, 0.00, 0.00, &
                                                                 0.0 , 0.05, 0.00, &
                                                                 0.00, 0.00, 0.00/),shape(vi))
REAL    , DIMENSION(1:3, 1:3)           ::  off_shift= 0.00
REAL                                    ::  renergy  = 17.480782
REAL                                    ::  rlambda  =  0.709260
LOGICAL                                 ::  l_energy = .false.
INTEGER                                 ::  dif_size_of = 0.0
!
LOGICAL                                 ::  l_zone       = .false.
REAL   , DIMENSION(1:3)                 ::  zone_uvw     = (/ 0.0, 0.0, 1.0/) 
REAL   , DIMENSION(1:3)                 ::  zone_ewald   = 0.0
REAL                                    ::  zone_res     = 0.0
REAL                                    ::  zone_delta_d = 0.015
!
REAL(KIND=PREC_DP), DIMENSION(4,3)      ::  diff_res             ! Resolution sigma and vectors
REAL(KIND=PREC_DP), DIMENSION(3,3)      ::  diff_tr              ! Resolution transformation matrix
!
INTEGER                                 ::  four_last = FOUR_NN  ! No Fourier calculated yet
!
END MODULE diffuse_mod
