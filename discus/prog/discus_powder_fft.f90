MODULE powder_fft_mod
!
CONTAINS
!
!*******************************************************************************
!
SUBROUTINE  fft_fq(npkt_wrt, xwrt, ywrt, qmin, qmax, deltaq, rmin, rmax, rstep, &
                   npkt_fft, npkt_pdf, xfour, yfour)
!
USE fast_fourier_mod
USE spline_mod
USE wink_mod
use errlist_mod
use precision_mod
!
implicit none
!
INTEGER                       , INTENT(IN)  :: npkt_wrt
REAL(kind=PREC_DP)   , DIMENSION(0:npkt_wrt), INTENT(IN)  :: xwrt
REAL(kind=PREC_DP)   , DIMENSION(0:npkt_wrt), INTENT(IN)  :: ywrt
REAL(KIND=PREC_DP)            , INTENT(IN)  :: qmin, qmax, deltaq
REAL(KIND=PREC_DP)            , INTENT(IN)  :: rmin, rmax, rstep
INTEGER                       , INTENT(IN)  :: npkt_fft !points in powder pattern for Fast Fourier = 2**16
INTEGER                       , INTENT(IN)  :: npkt_pdf
REAL(kind=PREC_DP)   , DIMENSION(0:npkt_pdf), INTENT(OUT) :: xfour
REAL(kind=PREC_DP)   , DIMENSION(0:npkt_pdf), INTENT(OUT) :: yfour
!
INTEGER   :: i
INTEGER   :: lensav      ! Array size for ip
INTEGER   :: lenwrk      ! Array size for w
INTEGER   :: iqmin
INTEGER   :: irmin
INTEGER   :: nlow = 0
REAL(KIND=PREC_DP) :: dq
REAL(KIND=PREC_DP) :: qmax_l
REAL(KIND=PREC_DP), DIMENSION(:), ALLOCATABLE :: temp   ! Temporary intensities for FFT
INTEGER           , DIMENSION(:), ALLOCATABLE :: ip
REAL(KIND=PREC_DP), DIMENSION(:), ALLOCATABLE :: w
REAL(KIND=PREC_DP), DIMENSION(:), ALLOCATABLE :: xfft   ! Temporary array for FFT result
REAL(KIND=PREC_DP), DIMENSION(:), ALLOCATABLE :: yfft   ! Temporary array for FFT result
!
!write(*,*) ' PDF OUT ', npkt_wrt, npkt_fft, npkt_pdf, rmin, rmax, (rmax-rmin) / REAL((npkt_pdf-1), KIND=PREC_DP)
!write(*,*) ' PDF ste ', (rmax-rmin) / REAL((npkt_pdf-1), KIND=PREC_DP), rstep
!write(*,*) ' PDF q   ', qmin, qmax, deltaq
!write(*,*) ' PDFxwrt ', xwrt(1), xwrt(2), xwrt(npkt_wrt)
!write(*,*) ' PDFywrt ', ywrt(1), ywrt(2), ywrt(npkt_wrt)
!n
nlow = 0
dq = deltaq !  (xwrt(npkt_wrt)-xwrt(1))/(npkt_wrt-1)
qmax_l   =  PI/rstep           ! by using 2*PI/rstep, the step size in direct space is halved
!
iqmin = MAX(0,NINT(qmin/deltaq))
lensav= 4+INT(SQRT(FLOAT(npkt_fft)/2))
lenwrk= npkt_fft*5/4-1
ALLOCATE(temp(0:npkt_fft+1))
ALLOCATE(ip(0:lensav))
ALLOCATE(w (0:lenwrk))
ALLOCATE(xfft(0:npkt_fft+1))
ALLOCATE(yfft(0:npkt_fft+1))
temp = 0.0D0
ip   = 0
w    = 0.0D0
!
DO i=0,iqmin-1                       ! Augment straight line to q=0
   temp(i) = REAL(i,KIND=PREC_DP)*dq*ywrt(1)/xwrt(1)
ENDDO
!
temp(iqmin:iqmin+npkt_wrt-1) = ywrt(1:npkt_wrt) ! Add actual powder pattern
!open(77,file='POWDER/augm_fft.FQ',status='unknown')
!DO i=0,iqmin+npkt_wrt-1
!  write(77,'(2(2x,G17.7E3))') (i)*dq,temp(i)
!enddo
!close(77)
!
CALL ddst(npkt_fft, 1, temp, ip, w)
!
irmin = nint(rmin/rstep)
!write(*,*) 'PDF OUT ', rmin, rstep, irmin, npkt_pdf
!write(*,*) 'PDF_OUT ', dq, qmax, (npkt_fft-1)*dq
qmax_l = (npkt_fft-1)*dq
!write(*,*) 'PDF_OUT ', qmax, rstep, PI/qmax
!open(77,file='POWDER/fft.PDF',status='unknown')
!DO i=0,npkt_pdf+irmin
!  write(77,'(2(2x,G17.7E3))') ( 0.5+i)*PI/QMAX_l,temp(i)*2/PI*dq  ! -0.5
!enddo
!close(77)
!open(77,file='POWDER/fft_2.PDF',status='unknown')
DO i=0,npkt_fft+1
  xfft (i) = (i+0.50)*PI/qmax_l
  yfft (i) = temp(i  )*2/PI*dq
!  write(77,'(2(2x,G17.7E3))') xfft(i), yfft(i)
ENDDO
!write(*,*) ' temp ', allocated(temp)
!write(*,*) ' temp ', lbound(temp), ubound(temp), minval(temp), maxval(temp)
!write(*,*) ' xfft ', lbound(xfft), ubound(xfft), minval(xfft), maxval(xfft)
!write(*,*) ' yfft ', lbound(yfft), ubound(yfft), minval(yfft), maxval(yfft)
!close(77)
!write(*,*) 'DO SPLINE ', REAL(rmin), REAL(rmax), REAL(rstep), npkt_pdf, nlow, npkt_fft+1
!write(*,*) ' DO SPLINE ', nlow, npkt_fft+1, npkt_pdf
CALL spline_prep(nlow, npkt_fft+1, xfft, yfft, rmin, rmax, rstep, npkt_pdf, xfour, yfour)
!write(*,*) ' IER ', ier_num, ier_typ
!write(*,*) ' xfour ', lbound(xfour), ubound(xfour), minval(xfour), maxval(xfour), xfour(ubound(xfour,1))
!write(*,*) ' yfour ', lbound(yfour), ubound(yfour), minval(yfour), maxval(yfour), yfour(ubound(yfour,1))
!
!write(*,*) ' temp ', allocated(temp), allocated(ip)
DEALLOCATE(temp)
DEALLOCATE(ip)
DEALLOCATE(w)
DEALLOCATE(xfft)
DEALLOCATE(yfft)
!
END SUBROUTINE  fft_fq
!
!*****7*****************************************************************
!
SUBROUTINE powder_conv_corrlin   (dat, tthmin, tthmax, dtth, sigma2, &
      corrlin, corrquad, rcut, pow_width, POW_MAXPKT)
!-
!     Convolute PDF with Gaussian function 
!     FWHM = (sigma**2 - corrlin/r)
!+
USE discus_config_mod 
!
USE gauss_lorentz_pseudo_mod
USE trig_degree_mod
USE wink_mod
!
IMPLICIT none 
!
INTEGER, INTENT(IN) :: POW_MAXPKT
REAL(kind=PREC_DP)   , DIMENSION(0:POW_MAXPKT), INTENT(INOUT) :: dat       ! Data to be convoluted
REAL(kind=PREC_DP)                            , INTENT(IN)    :: tthmin    ! 2Theta min
REAL(kind=PREC_DP)                            , INTENT(IN)    :: tthmax    ! 2Theta max
REAL(kind=PREC_DP)                            , INTENT(IN)    :: dtth      ! 2Theta step
REAL(kind=PREC_DP)                            , INTENT(IN)    :: sigma2    ! Gaussian Sigma^2
REAL(kind=PREC_DP)                            , INTENT(IN)    :: corrlin   ! 1/r deppendend width correction
REAL(kind=PREC_DP)                            , INTENT(IN)    :: corrquad  ! 1/r^2 deppendend width correction
REAL(KIND=PREC_DP)              , INTENT(IN)    :: rcut      ! minimum  distance for clin/(r-rmin)
REAL(kind=PREC_DP)                            , INTENT(IN)    :: pow_width ! Number of FWHM's to calculate
!
REAL(KIND=PREC_DP), PARAMETER :: eightln2  = 2.772588722239781237669D0 * 2.0D0
!
REAL(KIND=PREC_DP)            :: fwhm     ! Current FWHM at Theta
REAL(kind=PREC_DP), DIMENSION(0:POW_MAXPKT) :: dummy    ! temporary data (0:POW_MAXPKT) 
REAL(KIND=PREC_DP)            :: tth      ! Theta within convolution, main data set
REAL(KIND=PREC_DP)    :: sigmasq          ! actual scaled local sigma**2
REAL(KIND=PREC_DP)    :: sigmamin         ! minimum       local sigma**2
REAL(KIND=PREC_DP)    :: sigmaminsq       ! minimum       local sigma
REAL(KIND=PREC_DP)    :: eta              ! actual eta at current 2Theta
REAL(KIND=PREC_DP)    :: dist_min  ! minimum  distance for clin/(r-rmin)
INTEGER :: imax, i, j, ii,jj  ! Dummy loop indices
INTEGER :: i1, i2          ! Pseudo Voigt lookup indices
REAL(kind=PREC_DP)    :: pseudo          ! scale factor for lookup table
INTEGER :: max_ps 
!                                                                       
!------ Now convolute                                                   
!                                                                       
imax = INT( (tthmax - tthmin) / dtth )
!write(*,*) 'CONV ', tthmin, tthmax, dtth, imax
sigmasq = sigma2!*SQRT(eightln2)
sigmamin = sigmasq * 0.10D0
!sigmamin = sigmasq * 0.050D0
sigmaminsq = sqrt(sigmamin)
dist_min = REAL(rcut, KIND=PREC_DP)
!
eta = 0.0     ! Gaussian function
dummy = 0.0   ! dummy(:)
!write(*,*) ' IN CONV ', imax, tthmin, tthmax, dtth
!write(*,*) ' CQ      ', corrquad, sigmasq, dist_min
!tth = 1.810
!fwhm = SQRT(MAX(sigmasq - corrlin/(tth-dist_min) - corrquad/(tth-dist_min)**2, sigmamin))
!write(*,7777) ' FWHM ',tth, fwhm , sigmasq, sigmasq - corrlin/(tth-dist_min) - corrquad/(tth-dist_min)**2, sigmamin
!7777 format(a,f5.1, f9.5, f9.5, f12.5, f , POW_MAXPKT)9.5)
!tth = 2.300
!fwhm = SQRT(MAX(sigmasq - corrlin/(tth-dist_min) - corrquad/(tth-dist_min)**2, sigmamin))
!write(*,7777) ' FWHM ',tth, fwhm , sigmasq, sigmasq - corrlin/(tth-dist_min) - corrquad/(tth-dist_min)**2, sigmamin
!tth = 3.700
!fwhm = SQRT(MAX(sigmasq - corrlin/(tth-dist_min) - corrquad/(tth-dist_min)**2, sigmamin))
!write(*,7777) ' FWHM ',tth, fwhm , sigmasq, sigmasq - corrlin/(tth-dist_min) - corrquad/(tth-dist_min)**2, sigmamin
!tth = 10.000
!fwhm = SQRT(MAX(sigmasq - corrlin/(tth-dist_min) - corrquad/(tth-dist_min)**2, sigmamin))
!write(*,7777) ' FWHM ',tth, fwhm , sigmasq, sigmasq - corrlin/(tth-dist_min) - corrquad/(tth-dist_min)**2, sigmamin
!write(*,'(4f10.5)') sigmasq, sigmamin, sqrt(sigmasq*sqrt(eightln2)), sqrt(sigmamin*sqrt(eightln2)) 
!write(*,'(4f10.5)') corrlin, corrquad

!tth = 2.320
!fwhm = SQRT(MAX(sigmasq - corrlin/tth - corrquad/tth**2, sigmamin))
!write(*,'(a, 4f10.5)') ' FWHM ',tth, fwhm , sigmasq, sigmasq - corrlin/tth - corrquad/tth**2
!tth = 4.8000
!fwhm = SQRT(MAX(sigmasq - corrlin/tth - corrquad/tth**2, sigmamin))
!write(*,'(a, 4f10.5)') ' FWHM ',tth, fwhm , sigmasq, sigmasq - corrlin/tth - corrquad/tth**2
!tth = 10.000
!fwhm = SQRT(MAX(sigmasq - corrlin/tth - corrquad/tth**2, sigmamin))
!write(*,'(a, 4f10.5)') ' FWHM ',tth, fwhm , sigmasq, sigmasq - corrlin/tth - corrquad/tth**2
!tth = tthmax
!fwhm = SQRT(MAX(sigmasq - corrlin/tth - corrquad/tth**2, sigmamin))
!write(*,'(a, 4f10.5)') ' FWHM ',tth, fwhm , sigmasq, sigmasq - corrlin/tth - corrquad/tth**2
!open(45,file='POWDER/new.pdf', status='unknown')
!write(*,*) ' DISTMIN ', dist_min, tthmin, dtth
!open(45,file='POWDER/fwhm.pdf', status='unknown')
main_pts: DO i = 0, imax 
   tth = tthmin + i * dtth 
   IF(tth<0.01) THEN
      dummy(i) = dat(i)
      cycle main_pts
   ENDIF
!
!  fwhm = SQRT(MAX(sigmasq - corrlin/tth -corrquad/tth**2, 0.00001))
   fwhm = sigmaminsq
   if(tth>dist_min) &
   fwhm = SQRT(MAX((sigmasq - corrlin/(tth-dist_min) -corrquad/(tth-dist_min)**2)*sqrt(eightln2), sigmamin))
!
!if(abs(tth-2.2)<0.005) then
!if( abs(tth-nint(tth))<0.005 .and. tth < 30.) then
!  write(*,'(a,7f10.5)') ' TTH FWHM ', tth, tth-dist_min, sigmasq, corrquad/(tth-dist_min), &
!  corrlin/(tth-dist_min), &
!  fwhm, sqrt(sigmasq*sqrt(eightln2))
!!read(*,*) max_ps
!endif
!
   max_ps = max(21,INT((pow_width * fwhm) / dtth ))
   pseudo =     dtth/fwhm*glp_npt  ! scale factor for look up table
!
   ii = max(i-1-max_ps+1,1)
   jj = min(i+1+max_ps-1,imax, POW_MAXPKT)
!  do j=ii,jj
!     i1 = abs(j-i)
!     dummy(i) = dummy(i) + dat(j) *( glp_pseud_indx(i1, eta, fwhm))
!  enddo
!     cycle main_pts
!
   i1 = 0                               ! == 0 * dtth
   i2 = MIN(INT(2*i*pseudo), GLP_MAX)   ! == 2*i*dtth
!write(45,'(f8.3,2x, g20.8e3, f5.1, f7.0)') tth, fwhm, 0.0, real(max_ps)
   dummy (i) = dat (i) * ( glp_pseud_indx(i1, eta, fwhm)  &
                          -glp_pseud_indx(i2, eta, fwhm))                             
!
   ii = MAX (i - 1 - max_ps + 1, 1)
   first:DO j = ii, i - 1 
      i1 = MIN(    INT((i - j) * pseudo)      , GLP_MAX)  ! == tth1 = (i - j) * dtth
      i2 = MIN(max(INT((i + j) * pseudo),int((i1+1)*pseudo)), GLP_MAX)  ! == tth2 = (i + j) * dtth
      dummy(i) = dummy(i) + dat(j) *( glp_pseud_indx(i1, eta, fwhm)  &
                                     -glp_pseud_indx(i2, eta, fwhm))                    
   ENDDO first
!
   ii = MIN(i + 1 + max_ps - 1, imax, POW_MAXPKT)
   secnd: DO j = i + 1, ii 
      i1 = MIN(    INT((j - i) * pseudo)      , GLP_MAX)  ! == tth1 = (j - i) * dtth
      i2 = MIN(max(INT((j + i) * pseudo),int((i1+1)*pseudo)), GLP_MAX)  ! == tth1 = (j + i) * dtth
      dummy(i) = dummy(i) + dat(j) *( glp_pseud_indx(i1, eta, fwhm)  &
                                     -glp_pseud_indx(i2, eta, fwhm))
   ENDDO secnd
ENDDO main_pts
!                                                                       
DO i = 0, imax 
   dat (i) = dummy (i) * dtth   ! scale with stepwidth
!!   write(45,'(2f18.8)') i*dtth+tthmin, dat(i)
ENDDO 
!close(45)
!                                                                       
END SUBROUTINE powder_conv_corrlin          
!
!*****7*****************************************************************
!
!
!*******************************************************************************
END MODULE powder_fft_mod
