MODULE fourier_reset_mod
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE fourier_reset
!
USE discus_allocate_appl_mod
USE diffuse_mod
IMPLICIT NONE
!+
!     Contains all variables for Fourier transform
!-
!USE precision_mod
!
CALL alloc_diffuse(1,1,1)
!
DIF_MAXAT   = 1 ! current size of array at
DIF_MAXSCAT = 1 ! current size of array at
IF(ALLOCATED(cfact))      cfact(:,:)      = (0.0D0,0.0D0)  ! (0:CFPKT, 1:MAXSCAT)
IF(ALLOCATED(cfact_pure)) cfact_pure(:,:) = (0.0D0,0.0D0)  ! (0:CFPKT, 1:MAXSCAT)
IF(ALLOCATED(csf))        csf(:)          = (0.0D0,0.0D0)  ! (1:MAXQXY)
IF(ALLOCATED(tcsf))       tcsf(:)         = (0.0D0,0.0D0)  ! (1:MAXQXY)
IF(ALLOCATED(acsf))       acsf(:)         = (0.0D0,0.0D0)  ! (1:MAXQXY)
IF(ALLOCATED(xat))        xat(:,:)        = 0.0            ! (1:NMAX, 1:3)
IF(ALLOCATED(istl))       istl(:)         = 0              ! (1:MAXQXY)
cex(:)    = (0.0D0,0.0D0)
dsi(:)    = 0.0D0          ! (1:MAXQXY)
xm(:)     = 0.0D0
win(:)    = 0.0D0
vin(:)    = 0.0D0
uin(:)    = 0.0D0
fave      = 0.0
num(:)    = 1
nlots     = 1
ilots     = LOT_OFF
ls_xyz(:) = 5
nxat      = 1
four_mode = INTERNAL
lot_all   = .false.
ffour     = .false.
lperiod   = .true.
four_log  = .false.
four_was_run  = .false. ! TRUE if a fourier has been calculated
!
braggmax = 0.0
braggmin = 0.0
diffuave = 0.0
diffusig = 0.0
diffumax = 0.0
diffumin = 0.0
ps_low   = 1.20
ps_high  = 0.01
zmin     = 0.0
zmax     = 0.0
!
lambda   = 'MOA1'
four_exp = 0
inc(:)   = (/ 121, 121,  1 /)
lmn(:)   = 0
ano      = .false.
ldbw     = .false.
lxray    = .true.
diff_radiation = RAD_XRAY
diff_power     = 4
eck      = reshape((/ 0.0, 0.0,  0.0, &
                      5.0, 0.0,  0.0, &
                      0.0, 5.0,  0.0, &
                      0.0, 0.0,  0.0/),shape(eck))
vi       = reshape((/0.05, 0.00, 0.00, &
                     0.0 , 0.05, 0.00, &
                     0.00, 0.00, 0.00/),shape(vi))
off_shift= 0.00
renergy  = 17.480782
rlambda  =  0.709260
l_energy = .false.
!
l_zone = .false.
zone_uvw(:)   = (/0.0, 0.0, 1.0/)
zone_ewald(:) = 0.0
zone_res      = 0.0
zone_delta_d  = 0.015
!
four_last = FOUR_NN  ! No Fourier calculated yet
!
four_accum = 0
four_symm  = .FALSE.
!
END SUBROUTINE fourier_reset
!
!*******************************************************************************
!
END MODULE fourier_reset_mod
