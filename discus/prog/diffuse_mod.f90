MODULE diffuse_mod
!+
!     Contains all variables for Fourier transform
!-
USE config_mod
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
INTEGER , PARAMETER  :: I2PI     = 2**16
INTEGER , PARAMETER  :: MASK     = I2PI-1
!
INTEGER , PARAMETER  :: RAD_XRAY = 1
INTEGER , PARAMETER  :: RAD_NEUT = 2
INTEGER , PARAMETER  :: RAD_ELEC = 3
!
INTEGER                                 ::  DIF_MAXAT    ! current size of array at
INTEGER                                 ::  DIF_MAXSCAT  ! current size of array at
!COMPLEX , DIMENSION(0:CFPKT, 1:MAXSCAT) ::  cfact        ! (0:CFPKT, 1:MAXSCAT)
!COMPLEX , DIMENSION(1:MAXQXY)           ::  csf          ! (1:MAXQXY)
!COMPLEX , DIMENSION(1:MAXQXY)           ::  tcsf         ! (1:MAXQXY)
!COMPLEX , DIMENSION(1:MAXQXY)           ::  acsf         ! (1:MAXQXY)
COMPLEX , DIMENSION(:, :), ALLOCATABLE  ::  cfact        ! (0:CFPKT, 1:MAXSCAT)
COMPLEX , DIMENSION(:)   , ALLOCATABLE  ::  csf          ! (1:MAXQXY)
COMPLEX , DIMENSION(:)   , ALLOCATABLE  ::  tcsf         ! (1:MAXQXY)
COMPLEX , DIMENSION(:)   , ALLOCATABLE  ::  acsf         ! (1:MAXQXY)
COMPLEX , DIMENSION(0:MASK)             ::  cex       = (0.0,0.0)
!REAL    , DIMENSION(1:MAXQXY)           ::  dsi          ! (1:MAXQXY)
!REAL    , DIMENSION(1:NMAX, 1:3)        ::  xat          ! (1:NMAX, 1:3)
REAL    , DIMENSION(:)   , ALLOCATABLE  ::  dsi          ! (1:MAXQXY)
REAL    , DIMENSION(:, :), ALLOCATABLE  ::  xat          ! (1:NMAX, 1:3)
REAL    , DIMENSION(1:3)                ::  xm        = 0.0
REAL    , DIMENSION(1:3)                ::  vin       = 0.0
REAL    , DIMENSION(1:3)                ::  uin       = 0.0
REAL                                    ::  fave      = 0.0
!INTEGER , DIMENSION(1:MAXQXY)           ::  istl         ! (1:MAXQXY)
INTEGER , DIMENSION(:)   , ALLOCATABLE  ::  istl         ! (1:MAXQXY)
INTEGER , DIMENSION(1:2)                ::  num       = 1
INTEGER                                 ::  nlots     = 1
INTEGER                                 ::  ilots     = LOT_OFF
INTEGER , DIMENSION(1:3)                ::  ls_xyz    = 5
INTEGER                                 ::  nxat      = 1
INTEGER                                 ::  four_mode = INTERNAL
LOGICAL                                 ::  ffour     = .false.
LOGICAL                                 ::  lperiod   = .true.
LOGICAL                                 ::  four_log  = .false.
!
REAL                                    ::  braggmax
REAL                                    ::  braggmin
REAL                                    ::  diffuave
REAL                                    ::  diffusig
REAL                                    ::  diffumax
REAL                                    ::  diffumin
REAL                                    ::  ps_low   = 1.20
REAL                                    ::  ps_high  = 0.01
REAL                                    ::  zmin
REAL                                    ::  zmax
!
CHARACTER(LEN=4)                        ::  lambda   = 'MOA1'
INTEGER                                 ::  four_exp = 0
INTEGER , DIMENSION(1:2)                ::  inc      = 121
LOGICAL                                 ::  ano      = .false.
LOGICAL                                 ::  ldbw     = .false.
LOGICAL                                 ::  lxray    = .true.
INTEGER                                 ::  diff_radiation = RAD_XRAY
INTEGER                                 ::  diff_power     = 4
REAL    , DIMENSION(1:3, 1:3)           ::  eck      = reshape((/-0.5,-0.5,0.0, &
                                                                  5.5,-0.5,0.0, &
                                                                  -0.5,5.5,0.0/),shape(eck))
REAL    , DIMENSION(1:3, 1:2)           ::  vi       = reshape((/0.05,0.00,0.00, &
                                                                 0.00,0.05,0.00/),shape(vi))
REAL                                    ::  rlambda  = 0.709260
INTEGER                                 ::  dif_size_of
!
!     COMMON /diffuse/     cfact,csf,cex,istl,ffour,nlots,ls_xyz,       &
!    &                     ilots,lperiod,fave,tcsf,xat,nxat,            &
!    &                     acsf,dsi,num,xm,vin,uin,four_log,            &
!    &                     four_mode
!     COMMON /qspace/      braggmax,braggmin,diffumax,diffumin,         &
!    &                     diffuave,diffusig,ps_low,ps_high,zmin,zmax
!     COMMON /fourier_cmm/ lambda,inc,ano,ldbw,lxray,                   &
!    &                     eck,vi,rlambda,four_exp
!
END MODULE diffuse_mod
