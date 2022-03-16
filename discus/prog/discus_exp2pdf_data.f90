MODULE exp2pdf_data_mod
!
USE precision_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=PREC_STRING)                                :: exp_load  = ' '  ! Load string data
CHARACTER(LEN=PREC_STRING)                                :: exp_cback  = ' ' ! Load string Background
CHARACTER(LEN=PREC_STRING)                                :: exp_csigma = ' ' ! Load string Sigma's
CHARACTER(LEN=PREC_STRING)                                :: exp_radiation = 'xray' ! Radiation 'xray', 'neutron', 'electron'
character(len=4)          , DIMENSION(:    ), ALLOCATABLE :: exp_atname       ! Atom names
REAL(kind=PREC_DP)        , DIMENSION(:    ), ALLOCATABLE :: exp_atocc        ! Occupancies
integer                                                   :: exp_natom = 0    ! Number of atom types
INTEGER                                                   :: exp_kload = 0    ! Data set within KUPLOT
INTEGER                                                   :: exp_kback = 0    ! Data set within KUPLOT
INTEGER                                                   :: exp_ksigma= 0    ! Sigma set within KUPLOT
INTEGER                                                   :: exp_kfq   = 0    ! Sigma set within KUPLOT
INTEGER                                                   :: exp_kupl  = 0    ! Data set within KUPLOT that needs to be kept
INTEGER                   , DIMENSION(3)                  :: exp_dim    ! Dimensions of data set
INTEGER                   , DIMENSION(3)                  :: exp_dim_b  ! Dimensions of data set
integer                                                   :: exp_nstep  ! Data points on equdistant grid
integer                                                   :: exp_nstep_f! Data points on equdistant grid
integer                                                   :: exp_npoly = 7    ! Order of polynomial through F(Q)
logical                                                   :: exp_comp_current = .true.  ! Composition is that of current structure
logical                                                   :: exp_inter        = .true.  ! Show intermittent plots 
REAL(kind=PREC_DP)        , DIMENSION(:,:,:), ALLOCATABLE :: exp_data   ! the actual data set
REAL(kind=PREC_DP)        , DIMENSION(:,:,:), ALLOCATABLE :: exp_back   ! the actual data set
REAL(kind=PREC_DP)        , DIMENSION(:,:,:), ALLOCATABLE :: exp_sigma  ! sigma at each data point
REAL(kind=PREC_DP)        , DIMENSION(:,:,:), ALLOCATABLE :: exp_sigmab ! sigma at each background point
REAL(kind=PREC_DP)        , DIMENSION(:    ), ALLOCATABLE :: exp_x      ! x-values of data set
REAL(kind=PREC_DP)        , DIMENSION(:    ), ALLOCATABLE :: exp_y      ! y-values of data set
REAL(kind=PREC_DP)        , DIMENSION(:    ), ALLOCATABLE :: exp_xb     ! x-values of background set
REAL(kind=PREC_DP)        , DIMENSION(:    ), ALLOCATABLE :: exp_faver2 ! Calculated <f>**2
REAL(kind=PREC_DP)        , DIMENSION(:    ), ALLOCATABLE :: exp_temp_x ! the actual data set
REAL(kind=PREC_DP)        , DIMENSION(:    ), ALLOCATABLE :: exp_temp_y ! the actual data set
REAL(kind=PREC_DP)        , DIMENSION(:    ), ALLOCATABLE :: exp_temp_dy ! the actual data set sigmas
!REAL(kind=PREC_DP)        , DIMENSION(:    ), ALLOCATABLE :: exp_temp_b ! the actual background
real(kind=PREC_DP)                                        :: exp_bscale = 1.0D0 ! Background scale
real(kind=PREC_DP)                                        :: exp_qmin   ! Internal      Qmin
real(kind=PREC_DP)                                        :: exp_qmax   ! Internal      Qmax
real(kind=PREC_DP)                                        :: exp_qmin_i = 0.0D0   ! User supplied "instrumental" Qmin
real(kind=PREC_DP)                                        :: exp_qmin_f = 0.0D0   ! User supplied Qmin for Fourier
real(kind=PREC_DP)                                        :: exp_qmax_f = 1.0D9 ! User supplied Qmax for Fourier
real(kind=PREC_DP)                                        :: exp_qmax_u = 1.0D9 ! User supplied Qmax for adhoc correction
logical                                                   :: exp_qmin_il= .false. ! User supplied Qmin for adhoc correction
logical                                                   :: exp_qmin_fl= .false. ! User supplied Qmin for Fourier
logical                                                   :: exp_qmax_fl= .false. ! User supplied Qmax for Fourier
logical                                                   :: exp_qmax_ul= .false. ! User supplied Qmax for adhoc correction
logical                                                   :: exp_qfirst_l =  .false.  ! Qvalue at first maximum
real(kind=PREC_DP)                                        :: exp_qfirst_o =  0.000D0  ! Qvalue at first maximum
real(kind=PREC_DP)                                        :: exp_qfirst_c =  0.000D0  ! Qvalue at first maximum
real(kind=PREC_DP)                                        :: exp_qscale   =  0.000D0  ! Qvalue at first maximum
real(kind=PREC_DP)                                        :: exp_qstep =   0.001D0  ! Internal Q-step usually 0.001
real(kind=PREC_DP)                                        :: exp_rmin  =   0.01D0   ! PDF Rmin
!
CHARACTER(LEN=PREC_STRING)                                :: exp_outgr   = 'discus.grobs' ! Write GROBS
CHARACTER(LEN=PREC_STRING)                                :: exp_outiq   = 'discus.iqobs' ! Radiation 'xray', 'neutron', 'electron'
CHARACTER(LEN=PREC_STRING)                                :: exp_outfq   = 'discus.fqobs' ! Radiation 'xray', 'neutron', 'electron'
CHARACTER(LEN=PREC_STRING)                                :: exp_outsq   = 'discus.sqobs' ! Radiation 'xray', 'neutron', 'electron'
logical                                                   :: exp_outgr_l = .true.         ! Radiation 'xray', 'neutron', 'electron'
logical                                                   :: exp_outiq_l = .false.        ! Write IQ
logical                                                   :: exp_outfq_l = .false.        ! Write FQ
logical                                                   :: exp_outsq_l = .false.        ! Write SQ
real(kind=PREC_DP)                                        :: exp_rmax  = 100.01D0   ! PDF Rmax
real(kind=PREC_DP)                                        :: exp_rstep =   0.01D0   ! PDF Rstep
integer                                                   :: exp_npdf   ! Number data points in PDF
REAL(kind=PREC_DP)        , DIMENSION(:    ), ALLOCATABLE :: exp_pdf_x  ! the actual PDF set
REAL(kind=PREC_DP)        , DIMENSION(:    ), ALLOCATABLE :: exp_pdf_y  ! the actual PDF set
!
END MODULE exp2pdf_data_mod
