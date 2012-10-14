MODULE rmc_mod
!+
!     Variables for REVERSE-MONTE-CARLO level
!-
USE config_mod
!
SAVE
!
INTEGER                 :: RMC_MAXSCAT       = 1
!
INTEGER, PARAMETER      :: rmc_mode_shift    = 1
INTEGER, PARAMETER      :: rmc_mode_swchem   = 2
INTEGER, PARAMETER      :: rmc_mode_swdisp   = 3
INTEGER, PARAMETER      :: rmc_mode_external = 4
!
INTEGER, PARAMETER      :: rmc_data_nipl = 1
INTEGER, PARAMETER      :: rmc_data_pgm  = 2
!
INTEGER, PARAMETER      :: rmc_wic_eins = 1
INTEGER, PARAMETER      :: rmc_wic_sqrt = 2
INTEGER, PARAMETER      :: rmc_wic_log  = 3
INTEGER, PARAMETER      :: rmc_wic_lin  = 4
INTEGER, PARAMETER      :: rmc_wic_qua  = 5
INTEGER, PARAMETER      :: rmc_wic_inv  = 6
INTEGER, PARAMETER      :: rmc_wic_isq  = 7
INTEGER, PARAMETER      :: rmc_wic_dat  = 8
!
INTEGER, PARAMETER      :: rmc_local_all     = 1
INTEGER, PARAMETER      :: rmc_local_loc     = 2
INTEGER, PARAMETER      :: rmc_local_locsite = 3
INTEGER, PARAMETER      :: rmc_local_site    = 4
!
CHARACTER (LEN=80), DIMENSION(RMC_MAX_PLANES)      :: rmc_fname
CHARACTER (LEN=80)                                 :: rmc_lname
CHARACTER (LEN= 4), DIMENSION(RMC_MAX_PLANES)      :: rmc_lambda
!
complex, DIMENSION(RMC_MAX_SQ,RMC_MAX_LOTS)        :: rmc_csf
complex, DIMENSION(RMC_MAX_SQ,RMC_MAX_LOTS)        :: rmc_csf_new
complex, DIMENSION(0:CFPKT,DEF_MAXSCAT,RMC_MAX_PLANES) :: rcfact
!
REAL   , DIMENSION(3,3,RMC_MAX_SYM,RMC_MAX_PLANES) :: rmc_eck
REAL   , DIMENSION(3,2,RMC_MAX_SYM,RMC_MAX_PLANES) :: rmc_vi
REAL   , DIMENSION(RMC_MAX_Q)                      :: rmc_int
REAL   , DIMENSION(RMC_MAX_Q)                      :: rmc_wic
REAL   , DIMENSION(4,RMC_MAX_PLANES)               :: rmc_xy
REAL   , DIMENSION(RMC_MAX_PLANES)                 :: rmc_rlambda
REAL   , DIMENSION(RMC_MAX_PLANES)                 :: rmc_skal
REAL   , DIMENSION(RMC_MAX_PLANES)                 :: rmc_back
REAL   , DIMENSION(RMC_MAX_PLANES)                 :: rmc_chi2
REAL   , DIMENSION(RMC_MAX_PLANES)                 :: rmc_wtot
REAL   , DIMENSION(3,0:DEF_MAXSCAT   )                 :: rmc_maxmove
REAL   , DIMENSION(DEF_MAXSCAT,DEF_MAXSCAT)                :: rmc_mindist
REAL                                               :: rmc_mindist_max
REAL                                               :: rmc_ave
REAL                                               :: rmc_sigma
REAL                                               :: rmc_qmin,rmc_qmax
REAL                                               :: rmc_llim,rmc_ulim
!
INTEGER                                            :: rmc_nplane
INTEGER                                            :: rmc_data
INTEGER, DIMENSION(RMC_MAX_PLANES+1,RMC_MAX_SYM)   :: offsq
INTEGER, DIMENSION(RMC_MAX_PLANES+1)               :: offq
INTEGER, DIMENSION(RMC_MAX_PLANES)                 :: rmc_wic_typ
INTEGER, DIMENSION(RMC_MAX_SQ)                     :: ristl
INTEGER, DIMENSION(RMC_MAX_PLANES)                 :: rmc_nsym
INTEGER, DIMENSION(RMC_MAX_PLANES)                 :: rmc_constrain
INTEGER, DIMENSION(2,RMC_MAX_PLANES)               :: rmc_num
INTEGER, DIMENSION(3,RMC_MAX_LOTS)                 :: rmc_lots_orig
INTEGER, DIMENSION(3)                              :: rmc_csize
INTEGER                                            :: rmc_nlots,rmc_ilots
INTEGER                                            :: rmc_maxcyc,rmc_display
INTEGER                                            :: rmc_mode,rmc_local
INTEGER, DIMENSION(0:1)                            :: rmc_sel_prop
!
LOGICAL, DIMENSION(RMC_MAX_PLANES)                 :: rmc_lxray
LOGICAL, DIMENSION(RMC_MAX_PLANES)                 :: rmc_ano 
LOGICAL, DIMENSION(RMC_MAX_PLANES)                 :: rmc_ldbw 
LOGICAL, DIMENSION(:), ALLOCATABLE                 :: rmc_allowed  ! (0:RMC_MAXSCAT)
LOGICAL                                            :: rmc_doskal,rmc_doback
LOGICAL                                            :: rmc_calc_f,rmc_log,rmc_dosym,rmc_nosym
LOGICAL                                            :: rmc_ranloc,rmc_sel_atom
!
INTEGER                                            :: rmc_size_of  ! Bytes allocates for rmc
!     
!     COMMON /rmbl/ rmc_fname,rmc_lname,rmc_lambda,                     &
!    &              rmc_wic_typ,rmc_maxmove,                            &
!    &              rmc_maxcyc,rmc_nplane,rmc_int,rmc_wic,rmc_csf,      &
!    &              rmc_csf_new,rmc_sigma,rmc_calc_f,rmc_mindist,       &
!    &              rmc_display,rmc_mode,rmc_allowed,rmc_mindist_max,   &
!    &              rmc_lxray,rmc_ano,rmc_ldbw,rmc_skal,rmc_back,       &
!    &              rmc_wtot,rmc_chi2,rmc_doskal,rmc_ave,rmc_constrain, &
!    &              rmc_doback,rmc_rlambda,rmc_log,rmc_nsym,rmc_dosym,  &
!    &              rmc_qmin,rmc_qmax,rmc_llim,rmc_ulim,rmc_eck,rmc_vi, &
!    &              rmc_num,rmc_data,rmc_xy,ristl,rcfact,rmc_local,     &
!    &              rmc_lots_orig,rmc_nlots,rmc_ilots,rmc_csize,        &
!    &              rmc_nosym,rmc_ranloc,rmc_sel_atom,offq,offsq,       &
!    &              rmc_sel_prop
!
END MODULE rmc_mod
