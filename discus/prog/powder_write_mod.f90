module powder_period_temp_mod
!
use precision_mod
!
implicit none
!
INTEGER, save   :: npkt_pdf_temp    ! number of points in powder pattern ready to write
REAL(kind=PREC_DP), DIMENSION(:), ALLOCATABLE :: xwrt_pdf_temp ! y-values of powder pattern ready for output
REAL(kind=PREC_DP), DIMENSION(:), ALLOCATABLE :: ywrt_pdf_temp ! y-values of powder pattern ready for output
!
end module powder_period_temp_mod
!
!
!
MODULE powder_write_mod
!-
! Write the powder pattern in the selected format onto disk
!
! Multiple phases are averaged in phases_average. This 
! produces specialized output for the different output formats.
! Currently correlated motion is treated here, as a special case 
! for the PDF. Should eventually be moved into phases as well. 
!+
!
USE errlist_mod 
use precision_mod                                                                       
!
IMPLICIT NONE
!
PUBLIC
!
CONTAINS
!
SUBROUTINE powder_out (value, ltemp)
!-                                                                      
!     Write the powder pattern                                          
!+                                                                      
USE crystal_mod
USE discus_config_mod 
USE debye_mod 
USE diffuse_mod 
use discus_output_powder_mod
USE output_mod 
USE phases_mod
USE phases_set_mod
USE powder_mod 
USE powder_fft_mod
USE powder_tables_mod
USE pdf_mod
use powder_period_temp_mod
!
USE spline_mod
USE wink_mod
USE precision_mod
USE trig_degree_mod
USE qval_mod
!
IMPLICIT none 
!                                                                       
!     INTEGER, PARAMETER :: iff = 2 
!                                                                       
INTEGER, INTENT(IN) :: value ! Type of output
logical, intent(in) :: ltemp ! Prepare output in case of periodic boundary conditions
!                            ! True only for periodic boundary and PDF output
!                                                                       
INTEGER   :: POW_WR_MAXPKT   ! Maximum number of data points for write
INTEGER   :: ii, j , iii
integer   :: istart      ! Index at which user qmin starts
INTEGER   :: all_status  ! Allocation status
INTEGER   :: npkt        ! number of points in powder pattern; actual calculation
INTEGER   :: npkt_u      ! number of points in powder pattern; User input in 'powder'
INTEGER   :: npkt_equi   ! number of points in equidistant powder pattern
INTEGER   :: npkt_wrt    ! number of points in powder pattern ready to write
INTEGER   :: npkt_fft    ! number of points in powder pattern for Fast Fourier
LOGICAL   :: lread 
REAL(kind=PREC_DP), DIMENSION(:), ALLOCATABLE :: xpl  ! x-values of calculated powder pattern
REAL(kind=PREC_DP), DIMENSION(:), ALLOCATABLE :: ypl  ! y-values of calculated powder pattern
REAL(kind=PREC_DP), DIMENSION(:), ALLOCATABLE :: lpv  ! Values of LP correction versus Q/Theta 
REAL(kind=PREC_DP), DIMENSION(:), ALLOCATABLE :: y2a  ! y-values of splined    powder pattern
REAL(kind=PREC_DP), DIMENSION(:), ALLOCATABLE :: xwrt ! x-values of powder pattern ready for output
REAL(kind=PREC_DP), DIMENSION(:), ALLOCATABLE :: ywrt ! y-values of powder pattern ready for output
REAL(kind=PREC_DP) :: ttheta
REAL(kind=PREC_DP) :: lpscale      ! Scale factor introduced by LP correction
REAL(kind=PREC_DP) :: q=0.0
REAL(kind=PREC_DP)      :: normalizer
REAL(kind=PREC_DP) :: xmin, xmax, xdel
REAL(kind=PREC_DP)      :: xequ    ! x-position of equdistant curve
REAL(kind=PREC_DP)      :: yequ    ! y-value    of equdistant curve
REAL(kind=PREC_DP)      :: tthmin  ! minimum for equdistant curve
REAL(kind=PREC_DP)      :: tthmax  ! minimum for equdistant curve
REAL(KIND=PREC_DP)      ::   qmin  ! minimum for equdistant curve
REAL(KIND=PREC_DP)      ::   qmax  ! maximum for equdistant curve
REAL(KIND=PREC_DP)      :: deltaq  ! step    for equdistant curve
REAL(kind=PREC_DP)      :: arg
real(kind=PREC_DP)      :: back_zero
!REAL(KIND=PREC_DP) :: u2aver_scale = 2.00   ! Scale to multiply <u^2> if conversion
!                     ! with corrlin_corrquad is needed. This increases the calculated
!                     ! intensity actually by more than the damping by <u^2>, in order
!                     ! to sharpen the distances sufficiently for later broadening
!
real(kind=PREC_DP) :: rmax
! REAL(KIND=PREC_DP) :: rmin, rmax, rstep
!REAL(KIND=PREC_DP) :: rminf, rmaxf, rstepf
!INTEGER            :: npkt_pdf, npkt_pdff
!integer            :: npkt_back
!integer            :: npkt_ppdf
!real(kind=PREC_DP) :: qmin_back, qmax_back, qstp_back
!real(kind=PREC_DP) :: rmin_back, rmax_back, rstp_back
!REAL(kind=PREC_DP), DIMENSION(:), ALLOCATABLE :: xfour
!REAL(kind=PREC_DP), DIMENSION(:), ALLOCATABLE :: yfour
!REAL(kind=PREC_DP), DIMENSION(:), ALLOCATABLE :: xback
!REAL(kind=PREC_DP), DIMENSION(:), ALLOCATABLE :: yback
!REAL(kind=PREC_DP), DIMENSION(:), ALLOCATABLE :: rr_back
!REAL(kind=PREC_DP), DIMENSION(:), ALLOCATABLE :: gr_back
!
!REAL(kind=PREC_DP) :: sigma   ! sigma for PDF Corrlin correction
!
!
npkt_fft = 2**18
!
IF(.NOT. (value == val_inten  .OR. value == val_sq      .OR. &
          value == val_fq     .OR. value == val_iq      .OR. &
          value == val_f2aver .OR. value == val_faver2  .OR. &
                                   value == val_pdf            )) THEN
   ier_msg(1) = ' Powder output is defined only for:'
   ier_msg(2) = ' Intensity, S(Q), F(Q), <f>^2, <f^2>'
   ier_msg(3) = ' Intensity/N, PDF'
   ier_num = -124
   ier_typ = ER_APPL
   RETURN
ENDIF
!
!write(*,*) out_user_limits .and. .not. (cpow_form.eq.'tth'.or.cpow_form=='r'), &
!out_user_limits, .not. cpow_form.eq.'tth', .not. cpow_form.eq.'r', cpow_form
if(out_user_limits .and. .not. (cpow_form.eq.'tth'.or.cpow_form=='r')) then     ! Limits in 'output' menu
!  xmin = out_user_values(1)
!  xmax = out_user_values(2)
!  xdel = out_user_values(3)
   npkt_u = out_user_inc(1)
else                                                        ! Limits from powder menu
!  xmin   = pow_qmin_u
!  xmax   = pow_qmax_u
!  xdel   = pow_deltaq_u
!  npkt_u = NINT((xmax+xdel-xmin)/xdel) + 0
   npkt_u = pow_npkt_u
endif
!                                                                       
xmin = pow_qmin_c
xmax = pow_qmax_c
xdel = pow_deltaq_c
npkt = NINT((xmax+xdel-xmin)/xdel) + 0            ! Use automatic maximum for write until spline
!write(*,*) ' INITIAL P  ', npkt, pow_qmin_u, pow_qmax_u, pow_deltaq_u
!write(*,*) ' INITIAL C  ', npkt, pow_qmin_c, pow_qmax_c, pow_deltaq_c
!write(*,*) ' INITIAL    ', npkt, pow_qmin  , pow_qmax  , pow_deltaq  
!write(*,*) ' INITIAL x  ', npkt, xmin, xmax, xdel
!write(*,*) ' initial u  ', npkt_u, out_user_values(1:3), out_user_inc(1)
!
POW_WR_MAXPKT = MAX(npkt, npkt_u, out_user_inc(1), POW_MAXPKT)
!
ALLOCATE(xpl(0:POW_WR_MAXPKT),stat = all_status)  ! Allocate array for calculated powder pattern
ALLOCATE(ypl(0:POW_WR_MAXPKT),stat = all_status)  ! Allocate array for calculated powder pattern
ALLOCATE(lpv(0:POW_WR_MAXPKT),stat = all_status)  ! Allocate array for LP correction
xpl     = 0.0   ! (:)
ypl     = 0.0   ! (:)
lpv     = 0.0   ! (:)
!
CALL phases_average(xmin, xdel, npkt)          ! Calculate average: powder, f2aver, faver2, fu
!
lread = .false. 
IF (ier_num /= 0) THEN 
   DEALLOCATE(xpl    ,stat = all_status)  ! DeAllocate array for calculated powder pattern
   DEALLOCATE(ypl    ,stat = all_status)  ! DeAllocate array for calculated powder pattern
   DEALLOCATE(lpv    ,stat = all_status)  ! DeAllocate array for LP correction
   RETURN
ENDIF
!
IF(value == val_f2aver) THEN          ! Output is f^2 aver
   DO j = 0, npkt
          ypl (j) = (pow_f2aver(j))
   ENDDO
ELSEIF(value == val_faver2) THEN     ! Output is faver^2 
   DO j = 0, npkt
          ypl (j) = (pow_faver2(j))
   ENDDO
ELSEIF(value == val_sq .OR. value == val_fq) THEN
   DO j = 0, npkt
          ypl (j) = pow_sq  (j)   ! copy from convoluted pattern
                                  ! After phases_average
   ENDDO
!ELSEIF(value == val_pdf) THEN
!!
!!  Read S(Q) for PDF
!   DO j = 0, npkt
!          ypl (j) = pow_sq  (j)   ! copy from convoluted pattern
!                                  ! After phases_average
!   ENDDO
ELSE                              ! All other output is  inte
   DO j = 0, npkt
          ypl (j) = pow_conv(j)   ! copy from convoluted pattern
                                  ! After phases_average
   ENDDO
ENDIF           ! Output is if_block if(value==val_f2aver)
!
!DBGCQUAD
!open(77,file='POWDER/INITIAL.inte',status='unknown')
!DO ii=0,npkt
!  write(77,'(2(2x,G17.7E3))') xmin+(ii)*xdel,     ypl(ii)
!enddo
!close(77)
!
!------ copy the powder pattern into output array, if necessary this will be put on
!       equidistant scale
!                                                                       
lpv    = 0.0
lpv(0) = 1.0
!
!write(*,*) ' STEP INTO K12 ', xmin, xdel, ubound(xpl)
CALL pow_k12(npkt, POW_WR_MAXPKT, pow_ka21, pow_ka21_u, xmin, xdel, xpl, ypl)
!
! Adjust zero point along x-scale
!
!write(*,*) ' PRAE  ZERO ', npkt, xpl(1), xpl(npkt-1), xpl(npkt), pow_tthzero
!write(*,*) ' PRAE       ', pow_npkt_u, out_user_values(1), xmin, xdel, npkt_u
!write(*,*) ' PRAE TTHmin', 2.D0*asind(pow_qmin_u/2.0D0/zpi*rlambda), pow_qmin_u
!write(*,*) ' PRAE TTHmax', 2.D0*asind(pow_qmax_u/2.0D0/zpi*rlambda), pow_qmax_u
!write(*,*) ' PRAE tth0  ', pow_tthzero
!write(*,*) ' PRAE m,d   ', xmin, xdel
istart = 0
j      = 0  ! Flag for tth > 180
loop_zero: DO ii = 0, npkt
   q    = (ii) * xdel + xmin 
   arg = (q / 2.D0 /zpi *rlambda )
   if(arg>1.00D0) then
      j = ii
      exit loop_zero
   endif
   ttheta = 2.D0*asind ( (q / 2.D0 /zpi *rlambda ))
   if(q<pow_qmin_u) istart = ii
   lpv(ii) = polarisation (ttheta)
!      lpv(ii ) = lp                         ! For debye to get I(Q) ???
   IF (cpow_form.eq.'tth') THEN 
      xpl(ii) = ttheta - pow_tthzero
   ELSEIF (cpow_form.eq.'q  ') THEN 
      xpl(ii) = q - pow_qzero
   ELSEIF (cpow_form.eq.'r  ') THEN 
      xpl(ii) = q - pow_qzero                ! Initially the x-axis is in Q
   ENDIF 
ENDDO loop_zero 
!
if(j>0) then     ! 2Theta > 180°
   npkt = j -1
endif
!write(*,*) ' aft  m,d   ', xmin, xdel
!write(*,*) ' AFTER ZERO ', npkt, xpl(1), xpl(npkt-1), xpl(npkt), q, npkt_u, ttheta
!
lpscale = 1.0
!
!     Prepare S(Q) or F(Q)
!
normalizer = 1.0D0
prsq: IF(value == val_sq .or. value == val_fq   .OR. value == val_inten  .OR. & 
         value == val_iq .OR.                        value == val_pdf         ) THEN
!
   valq: IF(value == val_sq) THEN                   ! Calc S(Q)
      CONTINUE
   ELSEIF(value == val_fq) THEN  valq               ! Calc F(Q)IF
      ypl = ypl - 1.0
      DO j = 0, npkt   
         q = ((j)*xdel + xmin)
         ypl(j) = ypl(j) * q
      ENDDO
!
   ELSEIF(value == val_iq) THEN   valq              ! Calc F(Q)IF
         lpscale = 1.0
         IF(deb_conv .OR. .NOT.ldbw) THEN              ! DEBYE was done with convolution of ADP
            DO j = 0, npkt   
               q = (j)*xdel + xmin
               ypl(j) =  (ypl(j)                    /normalizer    &
                         +         lpv(j) * (                      &
                         + pow_faver2(j) - & !exp(-q**2*pow_u2aver*0.25)      * &
                          pow_f2aver(j)              ) )
            ENDDO
         ELSE
            DO j = 0, npkt   
               q = (j)*xdel + xmin
               ypl(j) =  (ypl(j)                    /normalizer   &
                         +         lpv(j) * (                      &
                         + pow_faver2(j) - exp(-q**2*pow_u2aver)   * &
                          pow_f2aver(j)              ))
            ENDDO
         ENDIF
   ELSEIF(value == val_inten) THEN   valq              ! Calc Intensity
      CONTINUE
   ENDIF valq
ENDIF prsq  !Prepare S(Q), F(Q)
!
!DBGCQUAD
!write(*,*) ' normalizer ', normalizer
!open(77,file='POWDER/normalized.inte',status='unknown')
!DO ii=0,npkt
!!               q = ((ii-1)*xdel + xmin)
!write(77,'(2(2x,G17.7E3))') xpl(ii)         , ypl(ii)
!enddo
!close(77)
!

!write(*,*) ' PRAE  K12  ', npkt, xpl(1), xpl(npkt-1), xpl(npkt), npkt_u
rmax     = out_user_values(2)
!write(*,*) ' AFTER K12I ', npkt, xpl(1), xpl(npkt-1), xpl(npkt), npkt_u
!DBGCQUAD
!open(77,file='POWDER/normalized.k12',status='unknown')
!DO ii=0,npkt
!!!               q = ((ii-1)*xdel + xmin)
!write(77,'(2(2x,G17.7E3))') xpl(ii)         , ypl(ii)
!enddo
!close(77)
!
IF( cpow_form == 'tth' ) THEN
!
!write(*,*) ' TTH ', out_user_limits
   IF(out_user_limits .AND. value /= val_pdf) THEN ! User provided values, not for PDF
      pow_tthmin   = out_user_values(1)
      pow_tthmax   = out_user_values(2)
      pow_deltatth = out_user_values(3)
      npkt_u       = out_user_inc(1)
!write(*,*) ' PRE   EQUI ', npkt_u, pow_tthmin, pow_tthmax, pow_deltatth
   ELSE                                          ! Convert q limits
!!write(*,*) ' NO USER LIMITS ', pow_qmin_u, xmin
!write(*,*) ' NO user        ', xpl(0), xpl(npkt)
      arg        = pow_qmin_u/(zpi) * rlambda / 2.d0 ! Directly with arg in asind()
      pow_tthmin = nint(2.*asind(arg)*1.0D5)/1.0D5     ! results in error ??????????
      arg        = MIN(1.0D0,pow_qmax_u/zpi * rlambda / 2.)
      if(arg>1.0D0) then
         pow_tthmax = 180.0D0
      else
         pow_tthmax = nint(2.*asind(arg)*1.0D5)/1.0D5
      endif
      pow_deltatth = xpl(istart+2)-xpl(istart+1)
      npkt_u     = int((pow_tthmax-pow_tthmin)/pow_deltatth) + 1
   ENDIF
!write(*,*) ' PRE   equi ', 0        , pow_tthmin, pow_tthmax, pow_deltatth, npkt_u
   IF(pow_tthmin < xpl(0) ) THEN           ! User lower limit too low!
      tthmin = nint(xpl(0)*1.0D5)/1.0D5
      npkt_u = INT((pow_tthmax-tthmin)/pow_deltatth) + 1             
   ELSE
      tthmin = pow_tthmin
   ENDIF
!write(*,*) ' PRE   equi ', 1        ,     tthmin, pow_tthmax, pow_deltatth, npkt_u
   IF(pow_tthmax > xpl(npkt) ) THEN              ! User upper limit too high!
      tthmax = nint(xpl(npkt)*1.0D5)/1.0D5
      npkt_u = INT((tthmax-tthmin)/pow_deltatth) + 1             
   ELSE
      tthmax = pow_tthmax
   ENDIF
!write(*,*) ' PRE   equi ', 2        ,     tthmin,     tthmax, pow_deltatth, npkt_u
   xmin = nint(tthmin*1.0D5)/1.0D5               ! Adjust limits needed later to cut 
   xmax = nint(tthmax*1.0D5)/1.0D5               ! off rounding errors
   npkt_equi =     INT((tthmax-tthmin)/pow_deltatth) + 1             
!write(*,*) ' PRE   EQUI ', npkt_equi, tthmin, tthmax, pow_deltatth, npkt_u
   ALLOCATE(y2a (0:POW_WR_MAXPKT),stat = all_status) ! Allocate array for calculated powder pattern
   ALLOCATE(xwrt(0:npkt_equi),stat = all_status)  ! Allocate array for powder pattern ready to write
   ALLOCATE(ywrt(0:npkt_equi),stat = all_status)  ! Allocate array for powder pattern ready to write
   xwrt = 0.0
   ywrt = 0.0
   y2a  = 0.0
   CALL spline (npkt, xpl, ypl, 1.d31, 1.d31, y2a)
   DO ii = 0, npkt_equi
      xequ = tthmin + (ii)*pow_deltatth
      CALL splint (npkt, xpl, ypl, y2a, xequ, yequ, ier_num)
      IF(ier_num/=0) THEN
         DEALLOCATE( xpl, stat = all_status)
         DEALLOCATE( ypl, stat = all_status)
         DEALLOCATE( lpv, stat = all_status)
         DEALLOCATE( y2a, stat = all_status)
         DEALLOCATE( xwrt, stat = all_status)
         DEALLOCATE( ywrt, stat = all_status)
         RETURN
      ENDIF
      xwrt(ii) = xequ
      ywrt(ii) = yequ
   ENDDO
   npkt_wrt = npkt_equi
   DEALLOCATE(y2a, stat = all_status)
!
ELSEIF( cpow_form == 'q') THEN        ! axis is Q
!        IF ( pow_axis      == POW_AXIS_TTH  .or.  &        ! Non matching form, spline onto equidistant steps
!            ((pow_four_type == POW_COMPL) .AND. value == val_pdf) .OR. &
!             pow_four_type == POW_DEBYE              ) THEN ! DEBYE, always spline
!write(*,*) ' OUT_USE_LIMITS', out_user_limits, out_user_values
   IF(out_user_limits) THEN                     ! User provided values
      qmin   = out_user_values(1)
      qmax   = out_user_values(2)
      deltaq = out_user_values(3)
   ELSEIF(npkt_u>npkt) THEN
      qmin   = pow_qmin_u
      qmax   = pow_qmax_u
      deltaq = pow_deltaq_u
   ELSE                                          ! Convert q limits
      qmin = pow_qmin
      qmax = pow_qmax
      deltaq = pow_deltaq
   ENDIF
!write(*,*) 'qmin, xpl)(0) ', qmin, xpl(0), qmin < xpl(0)
   IF(qmin < xpl(0) ) THEN                     ! User lower limit too low!
      qmin =            (INT( (xpl(1)            )/deltaq) + 1)*deltaq
      npkt_u =     NINT((qmax-qmin)/deltaq) + 1             
   ENDIF
!
!write(*,*) ' QLIMIT *? ', qmax, xpl(npkt), npkt
!write(*,*) ' QLIMIT *? ', qmin, qmax, deltaq , pow_deltaq_u
   IF(qmax > xpl(npkt) ) THEN                  ! User upper limit too high!
               qmax =            (INT( (         xpl(npkt))/deltaq) - 1)*deltaq
      npkt_u =     NINT((qmax-qmin)/deltaq) + 1             
   ENDIF
   xmin =   nint(qmin/deltaq)*deltaq              ! Adjust limits needed later to cut 
   xmax =   nint(qmax/deltaq)*deltaq              ! off rounding errors
   npkt_equi =     NINT((qmax-qmin)/deltaq) + 1             
   ALLOCATE(y2a (0:POW_WR_MAXPKT),stat = all_status) ! Allocate array for calculated powder pattern
   ALLOCATE(xwrt(0:npkt_equi),stat = all_status)  ! Allocate array for powder pattern ready to write
   ALLOCATE(ywrt(0:npkt_equi),stat = all_status)  ! Allocate array for powder pattern ready to write
   xwrt = 0.0
   ywrt = 0.0
   y2a  = 0.0
   CALL spline (npkt+1, xpl, ypl, 1.d31, 1.d31, y2a)
   DO ii = 0, npkt_equi
      xequ = qmin + (ii)*deltaq
      CALL splint (npkt+1, xpl, ypl, y2a, xequ, yequ, ier_num)
      IF(ier_num/=0) THEN
         DEALLOCATE( xpl, stat = all_status)
         DEALLOCATE( ypl, stat = all_status)
         DEALLOCATE( lpv, stat = all_status)
         DEALLOCATE( y2a, stat = all_status)
         DEALLOCATE( xwrt, stat = all_status)
         DEALLOCATE( ywrt, stat = all_status)
         RETURN
      ENDIF
      xwrt(ii) = xequ
      ywrt(ii) = yequ
   ENDDO
   npkt_wrt = npkt_equi
   DEALLOCATE(y2a, stat = all_status)
ELSE                    ! cpow_form == 
   ier_num = - 6
   ier_typ = ER_COMM
   DEALLOCATE( xpl, stat = all_status)
   DEALLOCATE( ypl, stat = all_status)
   DEALLOCATE( lpv, stat = all_status)
   RETURN
ENDIF                   ! cpow_form == 
!write(*,*) ' AFTER EQUI ', npkt_u, npkt_equi, npkt_wrt, lbound(xwrt), ubound(xwrt), lbound(ywrt), ubound(ywrt)
!write(*,*) ' AFTER EQUI ', npkt_wrt, xwrt(0), xwrt(npkt_wrt-1), xwrt(npkt_wrt)
npkt_wrt = npkt_u -1    ! xwrt is in range [0, npkt_u-1]
!write(*,*) ' AFTER EQUI ', npkt_wrt, xwrt(0), xwrt(npkt_wrt-1), xwrt(npkt_wrt)
!!
!DBGCQUAD
!open(77,file='POWDER/doequistep.inte',status='unknown')
!DO ii=1,npkt_wrt
!!!               q = ((ii-1)*xdel + xmin)
!write(77,'(2(2x,G17.7E3))') xwrt(ii)         , ywrt(ii)
!enddo
!close(77)
!
!     Scale intensity and add a background
!
place_ywrt: IF(value==val_inten) THEN
   IF( cpow_form == 'tth' ) THEN           ! "theta scale, convert x back to Q for Background parameters
      if(pow_back_q) then
         back_zero = 2.0*asind(pow_back(-1)/4.0_PREC_DP/PI*rlambda)
      else
         back_zero = pow_back(-1)
      endif
      DO ii=0,npkt_wrt
         q = 4.0D0*PI*sind(xwrt(ii))/rlambda
         ywrt(ii) = pow_scale*ywrt(ii)
         DO iii=0,pow_nback
            ywrt(ii) = ywrt(ii) + pow_back(iii)*(xwrt(ii)-back_zero)**iii
         ENDDO
      ENDDO
   elseif( cpow_form == 'q' ) THEN
      if(pow_back_q) then
         back_zero = pow_back(-1)
      else
         back_zero = 4.0D0*PI*sind(0.50_PREC_DP*pow_back(-1))/rlambda
      endif
      DO ii=0,npkt_wrt
         ywrt(ii) = pow_scale*ywrt(ii)
         DO iii=0,pow_nback
            ywrt(ii) = ywrt(ii) + pow_back(iii)*(xwrt(ii)-back_zero)**iii
         ENDDO
      ENDDO
   endif
ELSEIF(value==val_fq) THEN place_ywrt
   DO ii=0,npkt_wrt
      ywrt(ii) = pow_scale*ywrt(ii)
   ENDDO
!
ENDIF place_ywrt
!
!open(77,file='POWDER/final.write',status='unknown')
!DO ii=1,npkt_wrt
!write(77,'(2(2x,G17.7E3))') xwrt(ii), ywrt(ii)
!enddo
!close(77)
!read(*,*) ii
!     Finally write the pattern
!
if(ltemp) then         ! copy pdf into temporary array
   if(pow_lperiod .and. value==val_pdf) then
      if(allocated(xwrt_pdf_temp)) deallocate(xwrt_pdf_temp)
      if(allocated(ywrt_pdf_temp)) deallocate(ywrt_pdf_temp)
      allocate(xwrt_pdf_temp(0:npkt_wrt))
      allocate(ywrt_pdf_temp(0:npkt_wrt))
      xwrt_pdf_temp = xwrt
      ywrt_pdf_temp = ywrt
      npkt_pdf_temp = npkt_wrt
   endif
else                   ! normal write 
!   if(pow_lperiod .and. value==val_pdf) then
!!     Determine user output limits for PDF
!!     The PDF prepared for periodic boundary condition is splied onto the user grid.
!!
!      IF(out_user_limits) THEN
!         rmin     = out_user_values(1)
!         rmax     = out_user_values(2)
!         rstep    = out_user_values(3)
!         npkt_pdf = NINT((out_user_values(2)-out_user_values(1))/out_user_values(3))! + 1
!      ELSE
!         rmin     = pdf_rminu
!         rmax     = pdf_rmaxu
!         rstep    = pdf_deltaru
!         npkt_pdf = NINT((rmax-rmin)/pdf_deltaru)! + 1
!      ENDIF
!!     rstep = REAL((rmax-rmin)/(npkt_pdf-1), KIND=PREC_DP)
!!
!      ALLOCATE(y2a(0:POW_WR_MAXPKT),stat = all_status) ! Allocate array for calculated powder pattern
!      call spline(npkt_pdf_temp+1, xwrt_pdf_temp, ywrt_pdf_temp, 1.d31, 1.d31, y2a)
!      npkt_wrt = npkt_pdf
!!
!      if(allocated(xwrt)) deallocate(xwrt)
!      if(allocated(ywrt)) deallocate(ywrt)
!      allocate(xwrt(0:npkt_wrt))
!      allocate(ywrt(0:npkt_wrt))
!      DO ii = 0, npkt_wrt
!         xequ = rmin + (ii)*rstep
!         CALL splint (npkt_pdf+1, xwrt_pdf_temp, ywrt_pdf_temp, y2a, xequ, yequ, ier_num)
!         IF(ier_num/=0) THEN
!            DEALLOCATE( xpl, stat = all_status)
!            DEALLOCATE( ypl, stat = all_status)
!            DEALLOCATE( lpv, stat = all_status)
!            DEALLOCATE( y2a, stat = all_status)
!            DEALLOCATE( xwrt, stat = all_status)
!            DEALLOCATE( ywrt, stat = all_status)
!            RETURN
!         ENDIF
!         xwrt(ii) = xequ
!         ywrt(ii) = yequ
!      ENDDO
!      DEALLOCATE(y2a, stat = all_status)
!   endif
!
!write(*,*) ' FINAL WRITE'
!open(77,file='POWDER/final2.write',status='unknown')
!DO ii=1,npkt_wrt
!write(77,'(2(2x,G17.7E3))') xwrt(ii), ywrt(ii)
!enddo
!close(77)
!read(*,*) ii
!write(*,*) ' FINAL      ', npkt_wrt, xwrt(0), xwrt(npkt_wrt-1), xwrt(npkt_wrt)
   CALL powder_do_write (outfile, npkt_wrt, xwrt, ywrt)
endif
!
DEALLOCATE( ypl, stat = all_status)
DEALLOCATE( lpv, stat = all_status)
DEALLOCATE( xpl, stat = all_status)
DEALLOCATE( xwrt, stat = all_status)
DEALLOCATE( ywrt, stat = all_status)
!                                                                       
END SUBROUTINE powder_out                     
!
!*******************************************************************************
!
SUBROUTINE powder_out_pdf (value, ltemp)
!-                                                                      
!     Convert the powder pattern to a PDF
!+                                                                      
USE crystal_mod
USE discus_config_mod 
USE debye_mod 
USE diffuse_mod 
use discus_output_powder_mod
USE output_mod 
USE phases_mod
USE phases_set_mod
USE powder_mod 
USE powder_fft_mod
USE powder_tables_mod
USE pdf_mod
use powder_period_temp_mod
!
USE spline_mod
USE wink_mod
USE precision_mod
USE trig_degree_mod
USE qval_mod
!
IMPLICIT none 
!                                                                       
!     INTEGER, PARAMETER :: iff = 2 
!                                                                       
INTEGER, INTENT(IN) :: value ! Type of output
logical, intent(in) :: ltemp ! Prepare output in case of periodic boundary conditions
!                            ! True only for periodic boundary and PDF output
!                                                                       
INTEGER   :: POW_WR_MAXPKT   ! Maximum number of data points for write
INTEGER   :: ii, j
integer   :: istart      ! Index at which user qmin starts
INTEGER   :: all_status  ! Allocation status
INTEGER   :: npkt        ! number of points in powder pattern; actual calculation
INTEGER   :: npkt_u      ! number of points in powder pattern; User input in 'powder'
INTEGER   :: npkt_equi   ! number of points in equidistant powder pattern
INTEGER   :: npkt_wrt    ! number of points in powder pattern ready to write
INTEGER   :: npkt_fft    ! number of points in powder pattern for Fast Fourier
LOGICAL   :: lread 
REAL(kind=PREC_DP), DIMENSION(:), ALLOCATABLE :: xpl  ! x-values of calculated powder pattern
REAL(kind=PREC_DP), DIMENSION(:), ALLOCATABLE :: ypl  ! y-values of calculated powder pattern
REAL(kind=PREC_DP), DIMENSION(:), ALLOCATABLE :: lpv  ! Values of LP correction versus Q/Theta 
REAL(kind=PREC_DP), DIMENSION(:), ALLOCATABLE :: y2a  ! y-values of splined    powder pattern
REAL(kind=PREC_DP), DIMENSION(:), ALLOCATABLE :: xwrt ! x-values of powder pattern ready for output
REAL(kind=PREC_DP), DIMENSION(:), ALLOCATABLE :: ywrt ! y-values of powder pattern ready for output
REAL(kind=PREC_DP) :: ttheta
REAL(kind=PREC_DP) :: lpscale      ! Scale factor introduced by LP correction
REAL(kind=PREC_DP) :: q=0.0
REAL(kind=PREC_DP)      :: normalizer
REAL(kind=PREC_DP) :: xmin, xmax, xdel
REAL(kind=PREC_DP)      :: xequ    ! x-position of equdistant curve
REAL(kind=PREC_DP)      :: yequ    ! y-value    of equdistant curve
REAL(KIND=PREC_DP)      ::   qmin  ! minimum for equdistant curve
REAL(KIND=PREC_DP)      ::   qmax  ! maximum for equdistant curve
REAL(KIND=PREC_DP)      :: deltaq  ! step    for equdistant curve
REAL(kind=PREC_DP)      :: arg
REAL(KIND=PREC_DP) :: u2aver_scale = 2.00   ! Scale to multiply <u^2> if conversion
!                     ! with corrlin_corrquad is needed. This increases the calculated
!                     ! intensity actually by more than the damping by <u^2>, in order
!                     ! to sharpen the distances sufficiently for later broadening
!
REAL(KIND=PREC_DP) :: rmin, rmax, rstep
REAL(KIND=PREC_DP) :: rminf, rmaxf, rstepf
INTEGER            :: npkt_pdf, npkt_pdff
integer            :: npkt_back
integer            :: npkt_ppdf
real(kind=PREC_DP) :: qmin_back, qmax_back, qstp_back
real(kind=PREC_DP) :: rmin_back, rmax_back, rstp_back
REAL(kind=PREC_DP), DIMENSION(:), ALLOCATABLE :: xfour
REAL(kind=PREC_DP), DIMENSION(:), ALLOCATABLE :: yfour
REAL(kind=PREC_DP), DIMENSION(:), ALLOCATABLE :: xback
REAL(kind=PREC_DP), DIMENSION(:), ALLOCATABLE :: yback
REAL(kind=PREC_DP), DIMENSION(:), ALLOCATABLE :: rr_back
REAL(kind=PREC_DP), DIMENSION(:), ALLOCATABLE :: gr_back
!
REAL(kind=PREC_DP) :: sigma   ! sigma for PDF Corrlin correction
!
real(kind=PREC_DP) :: sigma_damp
real(kind=PREC_DP) :: yyy
real(kind=PREC_DP) :: qqq
!
npkt_fft = 2**18
!
IF(.NOT. (value == val_inten  .OR. value == val_sq      .OR. &
          value == val_fq     .OR. value == val_iq      .OR. &
          value == val_f2aver .OR. value == val_faver2  .OR. &
                                   value == val_pdf            )) THEN
   ier_msg(1) = ' Powder output is defined only for:'
   ier_msg(2) = ' Intensity, S(Q), F(Q), <f>^2, <f^2>'
   ier_msg(3) = ' Intensity/N, PDF'
   ier_num = -124
   ier_typ = ER_APPL
   RETURN
ENDIF
!
!write(*,*) out_user_limits .and. .not. (cpow_form.eq.'tth'.or.cpow_form=='r'), &
!out_user_limits, .not. cpow_form.eq.'tth', .not. cpow_form.eq.'r', cpow_form
if(out_user_limits .and. .not. (cpow_form.eq.'tth'.or.cpow_form=='r')) then     ! Limits in 'output' menu
!  xmin = out_user_values(1)
!  xmax = out_user_values(2)
!  xdel = out_user_values(3)
   npkt_u = out_user_inc(1)
else                                                        ! Limits from powder menu
!  xmin   = pow_qmin_u
!  xmax   = pow_qmax_u
!  xdel   = pow_deltaq_u
!  npkt_u = NINT((xmax+xdel-xmin)/xdel) + 0
   npkt_u = pow_npkt_u
endif
!                                                                       
xmin = pow_qmin_c
xmax = pow_qmax_c
xdel = pow_deltaq_c
npkt = NINT((xmax+xdel-xmin)/xdel) + 0            ! Use automatic maximum for write until spline
!write(*,*) ' INITIAL P  ', npkt, pow_qmin_u, pow_qmax_u, pow_deltaq_u
!write(*,*) ' INITIAL C  ', npkt, pow_qmin_c, pow_qmax_c, pow_deltaq_c
!write(*,*) ' INITIAL    ', npkt, pow_qmin  , pow_qmax  , pow_deltaq  
!write(*,*) ' INITIAL x  ', npkt, xmin, xmax, xdel
!write(*,*) ' initial u  ', npkt_u, out_user_values(1:3), out_user_inc(1)
!
POW_WR_MAXPKT = MAX(npkt, npkt_u, out_user_inc(1), POW_MAXPKT)
!
ALLOCATE(xpl(0:POW_WR_MAXPKT),stat = all_status)  ! Allocate array for calculated powder pattern
ALLOCATE(ypl(0:POW_WR_MAXPKT),stat = all_status)  ! Allocate array for calculated powder pattern
ALLOCATE(lpv(0:POW_WR_MAXPKT),stat = all_status)  ! Allocate array for LP correction
xpl     = 0.0   ! (:)
ypl     = 0.0   ! (:)
lpv     = 0.0   ! (:)
!
CALL phases_average(xmin, xdel, npkt)          ! Calculate average: powder, f2aver, faver2, fu
!
lread = .false. 
IF (ier_num /= 0) THEN 
   DEALLOCATE(xpl    ,stat = all_status)  ! DeAllocate array for calculated powder pattern
   DEALLOCATE(ypl    ,stat = all_status)  ! DeAllocate array for calculated powder pattern
   DEALLOCATE(lpv    ,stat = all_status)  ! DeAllocate array for LP correction
   RETURN
ENDIF
!
!  Read S(Q) for PDF
DO j = 0, npkt
   ypl (j) = pow_sq  (j)   ! copy from convoluted pattern
                           ! After phases_average
ENDDO
!
!DBGCQUAD
!open(77,file='POWDER/INITIAL.inte',status='unknown')
!DO ii=0,npkt
!  write(77,'(2(2x,G17.7E3))') xmin+(ii)*xdel,     ypl(ii)
!enddo
!close(77)
!
!------ copy the powder pattern into output array, if necessary this will be put on
!       equidistant scale
!                                                                       
lpv    = 0.0
lpv(0) = 1.0
!
!write(*,*) ' STEP INTO K12 ', xmin, xdel, ubound(xpl)
CALL pow_k12(npkt, POW_WR_MAXPKT, pow_ka21, pow_ka21_u, xmin, xdel, xpl, ypl)
!
! Adjust zero point along x-scale
!
!write(*,*) ' PRAE  ZERO ', npkt, xpl(1), xpl(npkt-1), xpl(npkt), pow_tthzero
!write(*,*) ' PRAE       ', pow_npkt_u, out_user_values(1), xmin, xdel, npkt_u
!write(*,*) ' PRAE TTHmin', 2.D0*asind(pow_qmin_u/2.0D0/zpi*rlambda), pow_qmin_u
!write(*,*) ' PRAE TTHmax', 2.D0*asind(pow_qmax_u/2.0D0/zpi*rlambda), pow_qmax_u
!write(*,*) ' PRAE tth0  ', pow_tthzero
!write(*,*) ' PRAE m,d   ', xmin, xdel
istart = 0
j      = 0  ! Flag for tth > 180
loop_zero: DO ii = 0, npkt
   q    = (ii) * xdel + xmin 
   arg = (q / 2.D0 /zpi *rlambda )
   if(arg>1.00D0) then
      j = ii
      exit loop_zero
   endif
   ttheta = 2.D0*asind ( (q / 2.D0 /zpi *rlambda ))
   if(q<pow_qmin_u) istart = ii
   lpv(ii) = polarisation (ttheta)
!      lpv(ii ) = lp                         ! For debye to get I(Q) ???
   IF (cpow_form.eq.'tth') THEN 
      xpl(ii) = ttheta - pow_tthzero
   ELSEIF (cpow_form.eq.'q  ') THEN 
      xpl(ii) = q - pow_qzero
   ELSEIF (cpow_form.eq.'r  ') THEN 
      xpl(ii) = q - pow_qzero                ! Initially the x-axis is in Q
   ENDIF 
ENDDO loop_zero 
!
if(j>0) then     ! 2Theta > 180°
   npkt = j -1
endif
!write(*,*) ' aft  m,d   ', xmin, xdel
!write(*,*) ' AFTER ZERO ', npkt, xpl(1), xpl(npkt-1), xpl(npkt), q, npkt_u, ttheta
!
lpscale = 1.0
!
! Prepare F(Q)
!
normalizer = 1.0D0
!
!open(77,file='POWDER/fq.inte',status='unknown')
ypl = ypl - 1.0
DO j = 0, npkt   
   q = (j)*xdel + xmin
   ypl(j) = ypl(j) * q
!write(77,'(2(2x,G17.7E3))') q               , ypl(j)
ENDDO
!close(77)
!write(*,*) ' F(Q) PREP ', .NOT. (deb_conv .OR. .NOT.ldbw), deb_conv, ldbw, cr_natoms
!
!
! If corrlin or corrquad, we need to divide by average Debye-Waller term
IF(pdf_clin_a/=0.0 .OR. pdf_cquad_a/=0.0) THEN
!
!  Divide by Debye-Waller
!
   DO ii=0, npkt
      ypl(ii) = ypl(ii)/EXP(-0.5*(pow_u2aver*u2aver_scale)*xpl(ii)**2)
   ENDDO
!
endif
!
rmax     = out_user_values(2)
!
!write(*,*) ' POWDER L  ', pow_qmin, pow_qmax, pow_deltaq
!write(*,*) ' USER LIM  ', pow_qmin_u, pow_qmax_u, pow_deltaq_u
!write(*,*) ' VALUE     ', value == val_pdf, value, out_user_limits
!write(*,*) ' NPKT_u    ', npkt_u
if(pdf_clin_a>0.0D0 .or. pdf_cquad_a>0.0) then
   qmin   = xpl(0)
   qmax   = xpl(npkt)
   deltaq = (qmax-qmin)/real(npkt-1,kind=PREC_DP)
else
   qmin = pow_qmin
   qmax = pow_qmax
   deltaq = pow_deltaq
endif
IF(qmin < xpl(0) ) THEN                     ! User lower limit too low!
   qmin =   (INT( (xpl(1)            )/deltaq) + 1)*deltaq
   npkt_u = NINT((qmax-qmin)/deltaq) + 1             
ENDIF
!
!write(*,*) ' QLIMIT *? ', qmax, xpl(npkt), npkt
!write(*,*) ' QLIMIT *? ', qmin, qmax, deltaq , pow_deltaq_u
IF(qmax > xpl(npkt) ) THEN                  ! User upper limit too high!
   qmax =   (INT( (         xpl(npkt))/deltaq) - 1)*deltaq
   npkt_u = NINT((qmax-qmin)/deltaq) + 1             
ENDIF
!
xmin =   nint(qmin/deltaq)*deltaq              ! Adjust limits needed later to cut 
xmax =   nint(qmax/deltaq)*deltaq              ! off rounding errors
npkt_equi =     NINT((qmax-qmin)/deltaq) + 1             
!write(*,*) ' NPKT_EQUI ', npkt_equi
!   ALLOCATE(y2a (0:POW_WR_MAXPKT),stat = all_status) ! Allocate array for calculated powder pattern
ALLOCATE(xwrt(0:npkt_equi),stat = all_status)  ! Allocate array for powder pattern ready to write
ALLOCATE(ywrt(0:npkt_equi),stat = all_status)  ! Allocate array for powder pattern ready to write
xwrt = 0.0
ywrt = 0.0
!NEW  For the PDF Q-data are on equidistant scale already, xwrt is equal to xpl 
xwrt = xpl
ywrt = ypl
npkt_wrt = npkt_equi
!
if(pdf_clin_a>0.0 .OR. pdf_cquad_a>0.0) then
!
!  Dampen last 10 A^-1
   qqq = (qmax -pow_qmax_u) * 0.800
   sigma_damp = qqq /deltaq*0.1250
   j = npkt_wrt-nint(qqq /deltaq)
!write(*,'(a, 3i8)')    ' DAMP FQ ', j, npkt_wrt, npkt_wrt-j
!write(*, '(a,5f12.5)') ' RANGE ', qqq, qmax, pow_qmax_u, sigma_damp, xwrt(j)
!   do ii= npkt_wrt-nint(qqq/deltaq), npkt_wrt
   do ii= j                        , npkt_wrt
      ywrt(ii) = ywrt(ii)*exp(-0.5D0*(real(ii-j,kind=PREC_DP)/sigma_damp)**2)
   enddo
endif
!
place_ywrt: IF(value==val_pdf) THEN  ! Transform F(Q) into PDF
   if(.not.ltemp .and. pow_lperiod) exit place_ywrt ! PDF has been prepared for periodic
!
! treat for correlated motion effects
!
!write(*,*) ' PDF CORRLIN ? ', pdf_clin_a>0.0 .OR. pdf_cquad_a>0.0
!write(*,*) ' NPKT_WRT ', npkt_wrt, xwrt(npkt_wrt)
!write(*,*) ' PDF_C*   ', pdf_clin_a, pdf_cquad_a
!write(*,*) ' CORRQUAD ', pdf_clin_a>0.0 .OR. pdf_cquad_a>0.0
   corr_if: IF(pdf_clin_a>0.0 .OR. pdf_cquad_a>0.0) THEN
      IF(out_user_limits) THEN
         rmin     = 0.01D0 ! out_user_values(1)
         rmax     = max(200.00_PREC_DP,out_user_values(2)*1.25D0)
         rstep    = min(  0.01_PREC_DP,out_user_values(3))
      ELSE
         rmin     = 0.01D0 !pdf_rminu
         rmax     = max(200.00_PREC_DP, pdf_rmaxu*1.25D0)
         rstep    = min(  0.01_PREC_DP, pdf_deltaru)
      ENDIF
      npkt_pdf = NINT((rmax-rmin)/rstep) + 1
      rstep = (rmax-rmin)/REAL((npkt_pdf-1), KIND=PREC_DP)
! write(*,*) ' PDF npkt', npkt_pdf, rmin, rmax, rstep
      ALLOCATE(xfour(0:npkt_pdf))
      ALLOCATE(yfour(0:npkt_pdf))
      zero_last1: DO ii = npkt_wrt,2,-1
         IF(ywrt(ii)*ywrt(ii-1)<0) THEN   ! Different signs , found zero point
            ywrt(ii) = 0.0
            EXIT zero_last1
         ELSEIF(ywrt(ii)==0.0) THEN       ! Found exact zero point
            EXIT zero_last1
         ELSE
            ywrt(ii) = 0.0
         ENDIF
      ENDDO zero_last1
      npkt_pdf = NINT((rmax-rmin)/pdf_deltaru) + 1
!
!open(77,file='POWDER/prae_corrlin.FQ',status='unknown')
!DO ii=1,npkt_wrt
!write(77,'(2(2x,G17.7E3))') xwrt(ii), ywrt(ii)
!enddo
!close(77)
!write(*,*) ' ABOUT TO DO 1st FFT INTO PDF ', npkt_fft, npkt_pdf
!write(*,*) ' ABOUT TO DO 1st FFT INTO PDF ', rmin, rmax, rstep 
!write(*,*) ' ABOUT TO DO 1st FFT INTO PDF ', npkt_wrt, ubound(xwrt)
      CALL fft_fq(npkt_wrt, xwrt, ywrt, qmin, qmax, deltaq, rmin, rmax, rstep, &
                  npkt_fft, npkt_pdf, xfour, yfour)
!write(*,*) ' after       1st FFT INTO PDF ', npkt_pdf
!
!open(77,file='POWDER/prae_corrlin.PDF',status='unknown')
!DO ii=1,npkt_pdf
!write(77,'(2(2x,G17.7E3))') xfour(ii), yfour(ii)
!enddo
!close(77)
!write(*,*) ' <u^2> , B', pow_u2aver, pow_u2aver*8.*3.1415**2
      sigma = 2.0*(pow_u2aver*u2aver_scale)              ! TO BE REPLACED BY ATOMIC B VALUE
      CALL powder_conv_corrlin(yfour, rmin,rmax, rstep,   &
                               sigma, pdf_clin_a, pdf_cquad_a, pdf_rcut, pow_width,   &
                               npkt_pdf)
!yfour(1:151) = 0.0D0
!     Dampen initial PDF for small distances to rough -4PI rho_0 r 
      xfour(npkt_pdf) = xfour(npkt_pdf-1) + rstep
      yfour(npkt_pdf) = yfour(npkt_pdf-1)
      j = int((1.25 - xfour(1))/rstep)
      sigma_damp = real(j,kind=PREC_DP)*0.0150_PREC_DP
      yyy = sum(yfour(j-30:j+30))/61.0D0
      do ii=1, j
! write(*,'(a, 2i4,3f12.5)') ' DAMPING ', ii, j-ii, exp(-0.50D0* ((j-ii)/sigma_damp)*+2), yfour(ii), &
!yfour(ii)*exp(-0.50D0* ((j-ii)/sigma_damp)*+2)
!  write(*,'(a, 2i4,3f12.5)') ' DAMPING ', ii, j-ii, exp(-0.50D0* ((j-ii)/sigma_damp)*+2), yfour(ii), &
!yfour(ii)*exp(-0.50D0* ((j-ii)/sigma_damp)*+2)
         yfour(ii) = xfour(ii)/xfour(j)*yyy      + &
                    (yfour(ii) - xfour(ii)/xfour(j)*yyy     ) * exp(-0.50D0* ((j-ii)/sigma_damp)**2)
!  yfour(ii) = yfour(ii)*exp(-0.50D0* ((j-ii)/sigma_damp)*+2)
      enddo
!
!     Dampen long distacne part of PDF, since the convolution does not do this part well
!
      j = int((5.0)/rstep)
      sigma_damp = j*0.0625
      do ii = npkt_pdf-j, npkt_pdf
!        write(*,'(a, 2i6,3f12.5)') ' d pdf   ', ii, ii-npkt_pdf+j, exp(-0.50D0* ((ii-npkt_pdf+j)/sigma_damp)*+2), yfour(ii), &
!        yfour(ii)*exp(-0.50D0* ((ii-npkt_pdf+j)/sigma_damp)*+2)
         yfour(ii) = yfour(ii)*exp(-0.50D0* ((ii-npkt_pdf+j)/sigma_damp)**2)
      enddo
!
      yfour = yfour*0.5_PREC_DP*pi
!open(77,file='POWDER/post_corrlin.PDF',status='unknown')
!DO ii=1,npkt_pdf
!write(77,'(2(2x,G17.7E3))') xfour(ii), yfour(ii)
!enddo
!close(77)
      qmin_back = pow_qmin_u ! 0.555D0
      qmax_back = pow_qmax_u !30.000D0
      qstp_back = pow_deltaq_u ! 0.005D0
      npkt_back = nint((qmax_back-qmin_back)/qstp_back) + 1
      allocate(xback(0:npkt_back))
      allocate(yback(0:npkt_back))
      CALL fft_fq(npkt_pdf, xfour, yfour, rmin, rmax, rstep, qmin_back, qmax_back, qstp_back, &
                  npkt_fft, npkt_back, xback, yback)
      xback(npkt_back) = xback(npkt_back-1) + qstp_back
      yback(npkt_back) = yback(npkt_back-1) 
!open(77,file='POWDER/post_corrlin.FQ',status='unknown')
!DO ii=1,npkt_back
!write(77,'(2(2x,G17.7E3))') xback(ii), yback(ii)
!enddo
!close(77)
!   rmin_back =   0.01D0
!   rmax_back = 100.00D0
!   rstp_back =   0.01D0
      IF(out_user_limits) THEN
         rmin_back     = 0.01D0 ! out_user_values(1)
         rmin_back     = out_user_values(1)
         rmax_back     = out_user_values(2)
         rstp_back     = out_user_values(3)
      ELSE
         rmin_back     = pdf_rminu
         rmax_back     = pdf_rmaxu
         rstp_back     = pdf_deltaru
      ENDIF
      npkt_ppdf = nint((rmax_back-rmin_back)/rstp_back) + 1
      allocate(rr_back(0:npkt_ppdf))
      allocate(gr_back(0:npkt_ppdf))
      CALL fft_fq(npkt_back, xback, yback, qmin_back, qmax_back, qstp_back, rmin_back, rmax_back, rstp_back, &
                  npkt_fft, npkt_ppdf, rr_back, gr_back)
!open(77,file='POWDER/post_corrlin.PPDF',status='unknown')
!DO ii=1,npkt_ppdf
!write(77,'(2(2x,G17.7E3))') rr_back(ii), gr_back(ii)
!enddo
!close(77)
      deallocate(xback)
      deallocate(yback)
      deallocate(xfour)
      deallocate(yfour)
      allocate(xfour(0:npkt_ppdf))
      allocate(yfour(0:npkt_ppdf))
      xfour = rr_back
      yfour = gr_back
      npkt_pdf = npkt_ppdf
!
!  The final limits to be written need to be adjusted, as corrlin convolution sets a 
!  rmin at 0.5, respectively a rmax at rmax*1.25
!
      IF(out_user_limits) THEN
         rminf    = out_user_values(1)
         rmaxf    = out_user_values(2)
         rstepf   = out_user_values(3)
         npkt_pdff= NINT((out_user_values(2)-out_user_values(1))/out_user_values(3)) + 1
!        npkt_pdf = NINT((out_user_values(2)*1.25-out_user_values(1))/out_user_values(3)) + 1
      ELSE
         rminf    = pdf_rminu
         rmaxf    = pdf_rmaxu
         rstepf   = pdf_deltaru
         npkt_pdff= NINT((rmax-rmin)/pdf_deltaru) + 1
      ENDIF
!
      DEALLOCATE(xwrt)
      DEALLOCATE(ywrt)
!
!RMIN   IF(rminf<rmin) THEN       !rminuser < rmin; rmin is set to 0.5 if User_values are present
      IF(rminf<rmin_back) THEN       !rminuser < rmin; rmin is set to 0.5 if User_values are present
!RMIN      j = NINT((rmin-rminf)/rstepf)
         j = NINT((rmin_back-rminf)/rstepf)
         ALLOCATE(xwrt(0:npkt_pdff))
         ALLOCATE(ywrt(0:npkt_pdff))
         DO ii = 0, j
            xwrt(ii) = rminf + (ii)*rstepf
            ywrt(ii) = 0.0
         ENDDO
         DO ii = j+1, npkt_pdff
            xwrt(ii) = xfour(ii-j)
            ywrt(ii) = yfour(ii-j)
         ENDDO
         npkt_wrt = npkt_pdff-1   ! Finally set correct points for write
!RMIN   ELSEIF(rminf>rmin) THEN  ! rminuser > rmin; rmin is set to 0.5 if User_values are present
      ELSEIF(rminf>rmin_back) THEN  ! rminuser > rmin; rmin is set to 0.5 if User_values are present
         ALLOCATE(xwrt(0:npkt_pdff))
         ALLOCATE(ywrt(0:npkt_pdff))
!RMIN         j = NINT((rmin-rminf)/rstepf)     ! j will be < 0
         j = NINT((rmin_back-rminf)/rstepf)     ! j will be < 0
         DO ii = 0, npkt_pdff-1
            xwrt(ii) = xfour(ii-j)
            ywrt(ii) = yfour(ii-j)
         ENDDO
         npkt_wrt = npkt_pdff-1   ! Finally set correct points for write
      ELSE                     ! rminuser == rmin
         ALLOCATE(xwrt(0:npkt_pdf+1))
         ALLOCATE(ywrt(0:npkt_pdf+1))
         j = 0
         DO ii=0,npkt_pdf - 1
            xwrt(ii) =xfour(ii+1)
            ywrt(ii) =yfour(ii+1)
         ENDDO
         npkt_wrt = npkt_pdf -1   ! Finally set correct points for write
      ENDIF
      DEALLOCATE(xfour)
      DEALLOCATE(yfour)
   ELSE corr_if       ! No correlated motion
      IF(out_user_limits) THEN
         rmin     = out_user_values(1)
         rmax     = out_user_values(2)
         rstep    = out_user_values(3)
         npkt_pdf = NINT((out_user_values(2)-out_user_values(1))/out_user_values(3)) + 1
      ELSE
         rmin     = pdf_rminu
         rmax     = pdf_rmaxu
         rstep    = pdf_deltaru
         npkt_pdf = NINT((rmax-rmin)/pdf_deltaru) + 1
      ENDIF
      rstep = REAL((rmax-rmin)/(npkt_pdf-1), KIND=PREC_DP)
!
      ALLOCATE(xfour(0:npkt_pdf))
      ALLOCATE(yfour(0:npkt_pdf))
!
      zero_last: DO ii = npkt_wrt,2,-1
         IF(ywrt(ii)*ywrt(ii-1)<0) THEN   ! Different signs , found zero point
            ywrt(ii) = 0.0
            EXIT zero_last
         ELSEIF(ywrt(ii)==0.0) THEN       ! Found exact zero point
            EXIT zero_last
         ELSE
            ywrt(ii) = 0.0
         ENDIF
      ENDDO zero_last
!
!write(*,*) ' PRAE_FFT.FQ'
!open(77,file='POWDER/prae_fft.FQ',status='unknown')
!DO ii=1,npkt_wrt
!write(77,'(2(2x,G17.7E3))') xwrt(ii), ywrt(ii)
!enddo
!close(77)
!read(*,*) ii
!     npkt_pdf = npkt_pdf + 1
      CALL fft_fq(npkt_wrt, xwrt, ywrt, qmin, qmax, deltaq, rmin, rmax, rstep, &
                  npkt_fft, npkt_pdf, xfour, yfour)
      yfour = yfour * PI*0.50_PREC_DP
!write(*,*) ' RLIMITS ', rmin, rmax, rstep, npkt_pdf
!open(77,file='POWDER/ende_corrlin.PDF',status='unknown')
!DO ii=1,npkt_pdf
!write(77,'(2(2x,G17.7E3))') xfour(ii), yfour(ii)
!enddo
!close(77)
      DEALLOCATE(xwrt)
      DEALLOCATE(ywrt)
      ALLOCATE(xwrt(0:npkt_pdf  ))
      ALLOCATE(ywrt(0:npkt_pdf  ))
      DO ii=0,npkt_pdf  
         xwrt(ii) =xfour(ii+1)
         ywrt(ii) =yfour(ii+1)
      ENDDO
      DEALLOCATE(xfour)
      DEALLOCATE(yfour)
      npkt_wrt = npkt_pdf-1
   ENDIF corr_if
!
   IF(pow_lperiod) THEN      ! Make the PDF periodic
      DO ii=0, npkt_wrt
         IF(xwrt(ii)<pow_period) THEN
            ywrt(ii) =ywrt(ii)/(1.0-1.5*(xwrt(ii)/pow_period*pow_period_cut)   &
                                   +0.5*(xwrt(ii)/pow_period*pow_period_cut)**3)
         ELSE
            ywrt(ii) = 0.00
         ENDIF
      ENDDO
   ENDIF
!
   DO ii=0,npkt_wrt
      ywrt(ii) = pow_scale*ywrt(ii)
   ENDDO
!
ENDIF place_ywrt
!
!open(77,file='POWDER/final.write',status='unknown')
!DO ii=1,npkt_wrt
!write(77,'(2(2x,G17.7E3))') xwrt(ii), ywrt(ii)
!enddo
!close(77)
!write(*,*) ' FINAL 1    ', npkt_wrt, xwrt(0), xwrt(npkt_wrt-1), xwrt(npkt_wrt)
!read(*,*) ii
!     Finally write the pattern
!
if(ltemp) then         ! copy pdf into temporary array
   if(pow_lperiod .and. value==val_pdf) then
      if(allocated(xwrt_pdf_temp)) deallocate(xwrt_pdf_temp)
      if(allocated(ywrt_pdf_temp)) deallocate(ywrt_pdf_temp)
      allocate(xwrt_pdf_temp(0:npkt_wrt))
      allocate(ywrt_pdf_temp(0:npkt_wrt))
      xwrt_pdf_temp = xwrt
      ywrt_pdf_temp = ywrt
      npkt_pdf_temp = npkt_wrt
   endif
else                   ! normal write 
   if(pow_lperiod .and. value==val_pdf) then
!     Determine user output limits for PDF
!     The PDF prepared for periodic boundary condition is splined onto the user grid.
!
      IF(out_user_limits) THEN
         rmin     = out_user_values(1)
         rmax     = out_user_values(2)
         rstep    = out_user_values(3)
         npkt_pdf = NINT((out_user_values(2)-out_user_values(1))/out_user_values(3))! + 1
      ELSE
         rmin     = pdf_rminu
         rmax     = pdf_rmaxu
         rstep    = pdf_deltaru
         npkt_pdf = NINT((rmax-rmin)/pdf_deltaru)! + 1
      ENDIF
!     rstep = REAL((rmax-rmin)/(npkt_pdf-1), KIND=PREC_DP)
!
      ALLOCATE(y2a(0:POW_WR_MAXPKT),stat = all_status) ! Allocate array for calculated powder pattern
      call spline(npkt_pdf_temp+1, xwrt_pdf_temp, ywrt_pdf_temp, 1.d31, 1.d31, y2a)
      npkt_wrt = npkt_pdf
!
      if(allocated(xwrt)) deallocate(xwrt)
      if(allocated(ywrt)) deallocate(ywrt)
      allocate(xwrt(0:npkt_wrt))
      allocate(ywrt(0:npkt_wrt))
      DO ii = 0, npkt_wrt
         xequ = rmin + (ii)*rstep
         CALL splint (npkt_pdf+1, xwrt_pdf_temp, ywrt_pdf_temp, y2a, xequ, yequ, ier_num)
         IF(ier_num/=0) THEN
            DEALLOCATE( xpl, stat = all_status)
            DEALLOCATE( ypl, stat = all_status)
            DEALLOCATE( lpv, stat = all_status)
            DEALLOCATE( y2a, stat = all_status)
            DEALLOCATE( xwrt, stat = all_status)
            DEALLOCATE( ywrt, stat = all_status)
            RETURN
         ENDIF
         xwrt(ii) = xequ
         ywrt(ii) = yequ
      ENDDO
      DEALLOCATE(y2a, stat = all_status)
   endif
!
!write(*,*) ' FINAL WRITE'
!open(77,file='POWDER/final2.write',status='unknown')
!DO ii=1,npkt_wrt
!write(77,'(2(2x,G17.7E3))') xwrt(ii), ywrt(ii)
!enddo
!close(77)
!read(*,*) ii
!write(*,*) ' FINAL      ', npkt_wrt, xwrt(0), xwrt(npkt_wrt-1), xwrt(npkt_wrt)
   CALL powder_do_write (outfile, npkt_wrt, xwrt, ywrt)
endif
!
DEALLOCATE( ypl, stat = all_status)
DEALLOCATE( lpv, stat = all_status)
DEALLOCATE( xpl, stat = all_status)
DEALLOCATE( xwrt, stat = all_status)
DEALLOCATE( ywrt, stat = all_status)
!                                                                       
END SUBROUTINE powder_out_pdf
!
!*******************************************************************************
!
subroutine powder_out_tof (value, ltemp)
!-
! Write the powder pattern as a TOF data set
!+
!                                                                       
use discus_output_powder_mod
use powder_mod
use output_mod 
!
use precision_mod
use wink_mod
!
implicit none
!
INTEGER, INTENT(IN) :: value ! Type of output
logical, intent(in) :: ltemp ! Prepare output in case of periodic boundary conditions
!                            ! True only for periodic boundary and PDF output
integer :: i
integer :: npkt_wrt
integer :: all_status
real(kind=PREC_DP) :: q      ! Q 
real(kind=PREC_DP) :: d      ! d-value 
real(kind=PREC_DP) :: t      ! time
REAL(kind=PREC_DP), DIMENSION(:), ALLOCATABLE :: xwrt ! x-values of powder pattern ready for output
REAL(kind=PREC_DP), DIMENSION(:), ALLOCATABLE :: ywrt ! y-values of powder pattern ready for output
!
ALLOCATE(xwrt(0:pow_npkt_u),stat = all_status)  ! Allocate array for calculated powder pattern
ALLOCATE(ywrt(0:pow_npkt_u),stat = all_status)  ! Allocate array for calculated powder pattern
!
do i=0, pow_npkt_u
   q = (real(i,kind=PREC_DP))*pow_deltaq + pow_qmin
   d = zpi/q
   t = pow_difc*d + pow_difa*d**2 + pow_tzero 
   if(d>0.0D0) t = t + pow_difb/d
   xwrt(pow_npkt_u-i) = t
   ywrt(pow_npkt_u-i) = pow_conv(i)
enddo
npkt_wrt = pow_npkt_u
CALL powder_do_write (outfile, npkt_wrt, xwrt, ywrt)
!
end subroutine powder_out_tof
!
!*******************************************************************************
!
SUBROUTINE powder_convolute                     
!-
!  Convolutes the powder pattern with profile function
!+
USE debye_mod
USE diffuse_mod
USE powder_mod
!
USE precision_mod
USE prompt_mod
USE support_mod
!
IMPLICIT none 
!                                                                       
INTEGER :: npkt
INTEGER :: j
INTEGER :: nzero      ! number of zeros in intensity
!
REAL(KIND=PREC_DP) :: xmin, xmax, xdel, xxmax
REAL(KIND=PREC_DP) :: scalef
REAL(KIND=PREC_DP) :: pow_tmp_sum
REAL(KIND=PREC_DP) :: pow_uuu_sum
!
REAL(KIND=PREC_DP)           :: ss       ! time
!
!if(pow_profile==0) then
!WRITE (output_io, * ) ' Starting convolution'!, pow_eta, pow_eta_l, pow_eta_q, pow_u, pow_v, pow_w
!endif
ss = seknds (0.0)
!
xmin  = 0.0 
xmax  = 0.0 
xxmax = 0.0 
xdel  = 0.0 
IF (pow_four_type.eq.POW_COMPL .or. pow_four_type==POW_NUFFT .or. pow_four_type==POW_GRID) THEN 
!  IF (pow_axis.eq.POW_AXIS_Q) THEN 
      xmin = pow_qmin 
      xmax = pow_qmax 
      xdel = pow_deltaq 
!  ELSEIF (pow_axis.eq.POW_AXIS_TTH) THEN 
!     xmin = pow_tthmin 
!     xmax = pow_tthmax 
!     xdel = pow_deltatth 
!  ELSE 
!     ier_num = - 104 
!     ier_typ = ER_APPL 
!     ier_msg (1) = 'Use command ==> set axis,{"tth"|"q"}' 
!     ier_msg (2) = 'within the powder menu to define the axis' 
!     ier_msg (3) = ' ' 
!     RETURN 
!  ENDIF 
   npkt = MIN(NINT((xmax+xdel-xmin)/xdel) + 2, POW_MAXPKT)
ELSEIF (pow_four_type.eq.POW_DEBYE ) THEN
!  IF (pow_axis.eq.POW_AXIS_Q) THEN 
      xmin = pow_qmin 
      xmax = pow_qmax 
      xdel = (pow_qmax - pow_qmin) / (num (1) ) 
!  ELSEIF (pow_axis.eq.POW_AXIS_TTH) THEN 
!     xmin = pow_tthmin 
!     xmax = pow_tthmax 
!     xdel = (pow_tthmax - pow_tthmin) / (num (1) ) 
!  ELSE 
!     ier_num = - 104 
!     ier_typ = ER_APPL 
!     ier_msg (1) = 'Use command ==> set axis,{"tth"|"q"}' 
!     ier_msg (2) = 'within the powder menu to define the axis' 
!     ier_msg (3) = ' ' 
!     RETURN 
!  ENDIF 
   npkt = MIN(num(1), POW_MAXPKT)
ENDIF 
!
IF (.not.(pow_four_type == POW_COMPL .or. pow_four_type==POW_NUFFT .or. pow_four_type==POW_GRID)) THEN 
!                                                              
!     This is a Debye calculation, copy rsf or csf into pow_conv
!                                                              
   IF (npkt    .le.POW_MAXPKT) THEN 
      DO j = 1, npkt    
         pow_conv (j) = (rsf (j) )     ! Double precision no longer needed
      ENDDO 
   ENDIF 
ELSE
   pow_conv(:) = 2.0D0*(pow_qsp(:))   ! Double precision no longer needed, Correct for half space
ENDIF 
!open(77,file='POWDER/prae_conv.inte',status='unknown')
!DO j=0,npkt
!  write(77,'(2(2x,G17.7E3))')     xmin+j *xdel, pow_conv(j)
!enddo
!close(77)
if(pow_profile==0) return
!                                                              
!- -Does the powder pattern have to be convoluted by a profile function?
!                                                              
!        Determine integral to scale after convolution
pow_tmp_sum = 0.0
nzero = 0
DO j=1,npkt
   pow_tmp_sum = pow_tmp_sum + pow_conv(j)
   if(pow_conv(j) == 0.0) nzero = nzero + 1
ENDDO
!lconv = .FALSE.
IF (pow_profile == POW_PROFILE_GAUSS) THEN 
   IF (pow_delta.gt.0.0) THEN 
      xxmax = xmax + xdel
      CALL powder_conv_res(pow_conv, xmin,xxmax, xdel,         &
                           pow_delta, POW_MAXPKT)
   ENDIF 
!  lconv = .TRUE.
elseif (pow_profile == POW_PROFILE_PEARS) then 
   
      xxmax = xmax + xdel
         CALL powder_conv_pears_uvw_asym(pow_conv, xmin,xxmax, xdel,   &
         pow_eta, pow_eta_l, pow_eta_q, pow_u, pow_v, pow_w, pow_asym,  &
         pow_width, POW_MAXPKT, pow_four_type, pow_axis, rlambda, pow_pr_fwhm)
ELSEIF (pow_profile == POW_PROFILE_PSVGT) THEN 
!   IF(pow_u/=0.0 .OR. pow_v/=0.0 .OR. pow_eta_l/=0.0 .OR. pow_eta_q/=0.0) THEN
      xxmax = xmax + xdel
!     IF(pow_asym(1,1)/=0.0 .OR. pow_asym(2,1)/=0.0 .OR.                 &
!        pow_asym(3,1)/=0.0 .OR. pow_asym(3,1)/=0.0            ) THEN       
   if(pow_u/=0.0 .or. pow_v/=0.0 .or. pow_w/=0.0) then
      if(maxval(pow_asym)>0.0 .or. minval(pow_asym)<0.0) then
         CALL powder_conv_psvgt_uvw_asym(pow_conv, xmin,xxmax, xdel,   &
         pow_eta, pow_eta_l, pow_eta_q, pow_u, pow_v, pow_w, pow_asym,  &
         pow_width, POW_MAXPKT, pow_four_type, pow_axis, rlambda, pow_pr_fwhm)
      ELSE          ! Symmetric case
         CALL powder_conv_psvgt_uvw(pow_conv, xmin,xxmax, xdel,   &
         pow_eta, pow_eta_l, pow_eta_q, pow_u, pow_v, pow_w, pow_width,  &
         POW_MAXPKT, pow_four_type, pow_axis, rlambda, pow_pr_fwhm)
      ENDIF
   endif
!   ELSE 
!write(*,*) ' powder_conv_psvgt_fix'
!      xxmax = xmax + xdel
!      CALL powder_conv_psvgt_fix (pow_conv, xmin,xxmax, xdel,   &
!      pow_eta, pow_w, pow_width, POW_MAXPKT, pow_four_type)
!   ENDIF 
!  lconv = .TRUE.
elseif (pow_profile == POW_PROFILE_TOF) THEN
      xxmax = xmax + xdel
   call powder_conv_tof(pow_conv, xmin,xxmax, xdel, pow_width, POW_MAXPKT,      &
        pow_difa, pow_difb, pow_difc, pow_tzero, pow_tof_a0, pow_tof_a1,        &
        pow_tof_b0, pow_tof_b1, pow_tof_bq, pow_tof_s0, pow_tof_s1, pow_tof_s2, &
        pow_tof_sq, pow_tof_z, pow_tof_y, pow_tof_x, pow_tof_siz, pow_tof_str,  &
        pow_four_type, POW_COMPL, POW_NUFFT, POW_GRID, POW_DEBYE)
ENDIF 
!open(77,file='POWDER/inte_conv.FQ',status='unknown')
!DO j=0,npkt
!  write(77,'(2(2x,G17.7E3))')     xmin+j *xdel, pow_conv(j)
!enddo
!close(77)
pow_uuu_sum = 0.0
DO j=1,npkt
   pow_uuu_sum = pow_uuu_sum + pow_conv(j)
ENDDO
scalef = 1.0
IF(pow_four_type.eq.POW_DEBYE) THEN
   scalef = pow_tmp_sum/pow_uuu_sum
ELSEIF(pow_four_type==POW_COMPL .or. pow_four_type==POW_NUFFT .or. pow_four_type==POW_GRID) THEN
   scalef = 1./xdel
ENDIF
pow_conv(:) = pow_conv(:) * scalef
ss = seknds (ss) 
WRITE(output_io, '('' Elapsed time Convolution : '',G13.6,'' sec'')') ss
!write(*,*) ' POST CONVOLUTE '
!open(77,file='POWDER/post_conv.FQ',status='unknown')
!DO j=0,npkt
!  write(77,'(2(2x,G17.7E3))')     xmin+j *xdel, pow_conv(j)
!enddo
!close(77)
!read(*,*) j
!
END SUBROUTINE powder_convolute                     
!
!*******************************************************************************
!
SUBROUTINE four_fq(npkt_wrt, xwrt, ywrt, rmin, rmax, npkt_pdf, xfour, yfour)
use precision_mod
!
USE wink_mod
!
INTEGER                       , INTENT(IN) :: npkt_wrt
REAL(kind=PREC_DP)   , DIMENSION(1:npkt_wrt), INTENT(IN) :: xwrt
REAL(kind=PREC_DP)   , DIMENSION(1:npkt_wrt), INTENT(IN) :: ywrt
REAL(KIND=PREC_DP)            , INTENT(IN) :: rmin, rmax
INTEGER                       , INTENT(IN) :: npkt_pdf
REAL(kind=PREC_DP)   , DIMENSION(1:npkt_pdf), INTENT(OUT) :: xfour
REAL(kind=PREC_DP)   , DIMENSION(1:npkt_pdf), INTENT(OUT) :: yfour
!
INTEGER :: i,j, k, kmax
REAL(KIND=PREC_DP) :: rr, dr
REAL(KIND=PREC_DP) :: dq
REAL(KIND=PREC_DP) :: qmin
REAL(KIND=PREC_DP), DIMENSION(:), ALLOCATABLE :: sine
INTEGER            :: iqmin
!
qmin = xwrt(1)
dr = (rmax-rmin)/(npkt_pdf-1)
dq = (xwrt(npkt_wrt)-xwrt(1))/(npkt_wrt-1)
iqmin = NINT(qmin/dq) - 1
!
! Create SINE Lookup
!
kmax = 10000 !NINT(ZPI*10000.D0)
ALLOCATE(sine(0:KMAX))
!
DO i=0,kmax
   sine(i) = SIN(ZPI*i/10000.0D0)
ENDDO
!
! Augment straight line from Q=0, F(0)=0 to qmin, F(qmin)
DO i = 1, npkt_pdf
   rr = rmin + (i-1)*dr
   xfour(i) = (rr)
   yfour(i) = 0.0
   DO j = 1, iqmin
      k = MOD(INT((j-1)*dq*rr*10000.0D0/ZPI),kmax)
      yfour(i) = yfour(i) + (j-1)*dq*ywrt(1)/qmin*sine(k)
   ENDDO
ENDDO
! Do rest of F(Q)
DO i = 1, npkt_pdf
   rr = rmin + (i-1)*dr
   DO j=1, npkt_wrt
      k = MOD(INT(xwrt(j)*rr*10000.0D0/ZPI),kmax)
      yfour(i) = yfour(i) + ywrt(j)*sine(k)
   ENDDO
   yfour(i) = yfour(i)*2/PI*dq
ENDDO
!
DEALLOCATE(sine)
!
END SUBROUTINE four_fq
!Q!
!Q!*******************************************************************************
!Q!
!Q      REAL(KIND=PREC_DP) function lorentz (ttheta, flag_fq) 
!Q!+                                                                      
!Q!-                                                                      
!Q      USE discus_config_mod 
!Quse precision_mod
!Q      USE powder_mod 
!Q      USE trig_degree_mod
!Q      IMPLICIT none 
!Q!                                                                       
!Q!                                                                       
!Q      REAL(KIND=PREC_DP)   , INTENT(IN) :: ttheta 
!Q      INTEGER, INTENT(IN) :: flag_fq
!Q!                                                                       
!Q!
!Q!Q      lorentz = 1.0
!Q      
!Q      IF (pow_four_type.eq.POW_DEBYE) THEN 
!Q         lorentz = 1.0 
!Q      ELSE 
!Q         IF(flag_fq==0) THEN
!Q         IF (pow_lp.eq.POW_LP_BRAGG) THEN 
!Q            lorentz = 0.5 / sind (0.5 * ttheta) / sind (ttheta) 
!Q         ELSEIF (pow_lp.eq.POW_LP_NEUT) THEN 
!Q            lorentz = 0.5 / sind (0.5 * ttheta) / sind (ttheta) 
!Q         ELSEIF (pow_lp.eq.POW_LP_NONE) THEN 
!Q            lorentz = 1.0 
!Q         ELSEIF (pow_lp.eq.POW_LP_SYNC) THEN 
!Q            lorentz = 0.5 / sind (0.5 * ttheta) / sind (ttheta) 
!Q         ENDIF 
!Q         ELSEIF(flag_fq==1) THEN 
!Q            lorentz = 0.5 / sind (0.5 * ttheta) / sind (ttheta) 
!Q         ENDIF 
!Q      ENDIF 
!Q!                                                                       
!Q      END FUNCTION lorentz                          
!*****7*****************************************************************
!Q      REAL(KIND=PREC_DP) FUNCTION polarisation (ttheta) 
!Q!+                                                                      
!Q!-                                                                      
!Q      USE discus_config_mod 
!Q      USE powder_mod 
!Quse precision_mod
!Q      USE trig_degree_mod
!Q      IMPLICIT none 
!Q!                                                                       
!Q!                                                                       
!Q      REAL(KIND=PREC_DP) :: ttheta 
!Q!
!Q!
!Q      polarisation = 1.0
!Q      
!Q      IF (pow_lp.eq.POW_LP_BRAGG) THEN 
!Q         polarisation = (1. + (cosd (ttheta) ) **2 * pow_lp_fac)        &
!Q         / (1. + pow_lp_fac)                                            
!Q      ELSEIF (pow_lp.eq.POW_LP_NEUT) THEN 
!Q         polarisation = 1.0 
!Q      ELSEIF (pow_lp.eq.POW_LP_NONE) THEN 
!Q         polarisation = 1.0 
!Q      ELSEIF (pow_lp.eq.POW_LP_SYNC) THEN 
!Q         polarisation = pow_lp_fac + (1. - pow_lp_fac) * (cosd (ttheta) &
!Q         ) **2 * pow_lp_cos                                             
!Q      ENDIF 
!Q!                                                                       
!Q      END FUNCTION polarisation                     
!Q!*****7*****************************************************************
!Q      REAL(KIND=PREC_DP) FUNCTION lorentz_pol (ttheta) 
!Q!+                                                                      
!Q!-                                                                      
!Q      USE discus_config_mod 
!Q      USE powder_mod 
!Q      USE trig_degree_mod
!Q      IMPLICIT none 
!Q!                                                                       
!Q!                                                                       
!Q      REAL(KIND=PREC_DP) :: ttheta 
!Q!                                                                       
!Q      IF (pow_four_type.eq.POW_DEBYE) THEN 
!Q         lorentz_pol = 1.0 
!Q      ELSE
!Q         lorentz_pol = (1-pow_lp_fac+pow_lp_fac*(cosd(pow_lp_ang))**2*(cosd(ttheta))**2)/ &
!Q                       (2.*(sind(0.5*ttheta))**2*cosd(0.5*ttheta))
!Q      ENDIF
!Q      END FUNCTION lorentz_pol
!*****7*****************************************************************
      SUBROUTINE powder_conv_res (dat, tthmin, tthmax, dtth, delta, POW_MAXPKT)
!-                                                                      
!     Convolute powder pattern with resolution function (Gaussian)      
!+                                                                      
      USE discus_config_mod 
      USE wink_mod
use precision_mod
!
      IMPLICIT none 
!                                                                       
!
      INTEGER, INTENT(IN) :: POW_MAXPKT
!                                                                       
      REAL(kind=PREC_DP) :: dat (0:POW_MAXPKT) 
      REAL(kind=PREC_DP) :: tthmin, tthmax, dtth, delta 
!                                                                       
      REAL(kind=PREC_DP) :: dummy (0:POW_MAXPKT) 
      REAL(kind=PREC_DP) :: gauss (0:2 * POW_MAXPKT) 
      REAL(kind=PREC_DP) :: tth
      INTEGER imax, i, j, ii 
      INTEGER max_ps 
!                                                                       
!------ Setup Gaussian                                                  
!                                                                       
      max_ps = int( (10.0 * delta) / dtth )
      DO i = 0, max_ps 
      tth = i * dtth 
      gauss (i) = 1.0 / sqrt ((pi)) / delta * exp ( - (tth**2 / delta**2))  
      ENDDO 
!                                                                       
      DO i = max_ps + 1, 2 * POW_MAXPKT 
      gauss (i) = 0.0 
      ENDDO 
!                                                                       
!------ Now convolute                                                   
!                                                                       
      imax = int( (tthmax - tthmin) / dtth )
      DO i = 0, imax 
      dummy (i) = dat (i) * (gauss (0) - gauss (2 * i) ) 
      ii = max (i - 1 - max_ps + 1, 0) 
      DO j = ii, i - 1 
      dummy (i) = dummy (i) + dat (j) * (gauss (i - j) - gauss (i + j) ) 
      ENDDO 
      ii = min (i + 1 + max_ps - 1, imax) 
      DO j = i + 1, ii 
      dummy (i) = dummy (i) + dat (j) * (gauss (j - i) - gauss (j + i) ) 
      ENDDO 
      ENDDO 
!                                                                       
      DO i = 0, imax 
      dat (i) = dummy (i) !* dtth 
      ENDDO 
!                                                                       
      END SUBROUTINE powder_conv_res                
!
!*****7*****************************************************************
!
SUBROUTINE powder_conv_psvgt_fix (dat, tthmin, tthmax, dtth, eta, &
      w, pow_width, POW_MAXPKT, pow_type)
!-                                                                      
!     Convolute powder pattern with resolution function (Pseudo-Voigt)  
!     Constant FWHM, Constant eta                                       
!+                                                                      
USE discus_config_mod 
!
USE gauss_lorentz_pseudo_mod
USE precision_mod
USE wink_mod
IMPLICIT none 
!                                                                       
!
INTEGER, INTENT(IN) :: POW_MAXPKT
!                                                                       
REAL(KIND=PREC_DP), DIMENSION(0:POW_MAXPKT), INTENT(INOUT) ::  dat !(0:POW_MAXPKT) 
REAL(KIND=PREC_DP)                         , INTENT(IN)    ::  tthmin, tthmax, dtth, eta
REAL(KIND=PREC_DP)                         , INTENT(IN)    ::  w 
REAL(KIND=PREC_DP)                         , INTENT(IN)    ::  pow_width 
INTEGER                      , INTENT(IN)    ::  pow_type   ! = 0==COMPl =1==DEBYE
!
INTEGER, PARAMETER  :: POW_COMPL = 0
INTEGER, PARAMETER  :: POW_DEBYE = 1
INTEGER, PARAMETER  :: POW_NUFFT = 2
INTEGER, PARAMETER  :: POW_GRID  = 3
!                                                                       
REAL(KIND=PREC_DP)                               :: fwhm
REAL(KIND=PREC_DP), DIMENSION(0:POW_MAXPKT)      :: dummy
REAL(KIND=PREC_DP), DIMENSION(0:2 * POW_MAXPKT)  :: psvgt
REAL(KIND=PREC_DP)                 :: eta_dp
REAL(KIND=PREC_DP)                 :: fwhm_dp
INTEGER                            :: imax, i, j, ii 
INTEGER                            :: max_ps 
!write(*,*) ' CONVOLUTION WITH PSEUDO      ', maxval(dat), pow_type==POW_COMPL, imax
!
if(w<1.E-8) return         ! FWHM == zero, no convolution
!                                                                       
!------ Setup Pseudo-Voigt                                              
!                                                                       
fwhm = sqrt (abs (w) ) 
max_ps = int( (pow_width * fwhm) / dtth )
psvgt = 0.0
eta_dp  = REAL(eta,KIND=PREC_DP)
fwhm_dp = REAL(fwhm,KIND=PREC_DP)
!write(*,'(8(f8.5,1x))') eta_dp, fwhm_dp
!write(*,*)
!open(77,file='pseudo.comp', status='unknown')
DO i = 0, max_ps 
   ii = MIN(INT(i*dtth/fwhm*GLP_NPT), GLP_MAX)
   psvgt (i) = glp_pseud_indx(ii  , eta_dp, fwhm_dp)
ENDDO 
!                                                                       
!------ Now convolute                                                   
!                                                                       
imax = int( (tthmax - tthmin) / dtth ) - 1
dummy = 0.0
!write(*,*) ' IMAX ', imax, MINVAL(dat), maxval(dat)
!write(*,*) ' TYPE ', pow_type, POW_COMPL, POW_DEBYE
!
!------ If COMPLETE, check for zeros, else do all points
!
IF(pow_type==POW_COMPL .or. pow_type==POW_NUFFT .or. pow_type==POW_GRID) THEN
   DO i = 0, imax 
!if(dat(i)>10.) write(*,'(i5,f7.3,f14.5)') i, tthmin + i*dtth, dat(i)
      dummy (i) = dat (i) * (psvgt (0) - psvgt (2 * i) ) 
      ii = max (i - 1 - max_ps + 1, 0  ) 
      first: DO j = ii, i - 1 
         IF(dat(j)==0.0) CYCLE first
         dummy (i) = dummy (i) + dat (j) * (psvgt (i - j) - psvgt (i + j) ) 
      ENDDO first
      ii = min (i + 1 + max_ps - 1, imax) 
      secnd: DO j = i + 1, ii 
         IF(dat(j)==0.0) CYCLE secnd
         dummy (i) = dummy (i) + dat (j) * (psvgt (j - i) - psvgt (j + i) ) 
      ENDDO secnd
   ENDDO 
ELSEIF(pow_type==POW_DEBYE) THEN       ! Debye calculation no checks fro zeros in DAT
   DO i = 0, imax 
      dummy (i) = dat (i) * (psvgt (0) - psvgt (2 * i) ) 
      ii = max (i - 1 - max_ps + 1, 0  ) 
      DO j = ii, i - 1 
         dummy (i) = dummy (i) + dat (j) * (psvgt (i - j) - psvgt (i + j) ) 
      ENDDO
      ii = min (i + 1 + max_ps - 1, imax) 
      DO j = i + 1, ii 
         dummy (i) = dummy (i) + dat (j) * (psvgt (j - i) - psvgt (j + i) ) 
      ENDDO
   ENDDO 
ENDIF
!write(*,*) ' IMAX ', imax, MINVAL(dummy), maxval(dummy)
!                                                                       
DO i = 0, imax 
   dat (i) = dummy (i) * dtth 
ENDDO 
!                                                                       
END SUBROUTINE powder_conv_psvgt_fix          
!
!*****7*****************************************************************
!
SUBROUTINE alt_powder_conv_psvgt_uvw (dat, tthmin, tthmax, dtth, eta0,&
      eta_l, eta_q, u, v, w, pow_width, POW_MAXPKT, pow_type, axis, &
      rlambda, pow_pr_fwhm )
!-
!     Convolute powder pattern with resolution function (Pseudo-Voigt)  
!     FWHM according to caglioti equation, Constant eta                 
!     FWHM = sqrt ( U*tan**2(Theta) + V*tan(Theta) + W)  IF 2Theta scale
!     FWHM = sqrt ( U*Q^2  +  V*Q  +  W               )  IF Q      scale
!     Symmetric case
!+
USE discus_config_mod 
!
USE gauss_lorentz_pseudo_mod
use precision_mod                                                                       
USE trig_degree_mod
USE wink_mod
!
IMPLICIT none 
!
INTEGER, INTENT(IN) :: POW_MAXPKT
REAL(KIND=PREC_DP)   , DIMENSION(0:POW_MAXPKT), INTENT(INOUT) :: dat   ! Data to be convoluted
REAL(KIND=PREC_DP)                            , INTENT(IN)    :: tthmin ! 2Theta min
REAL(KIND=PREC_DP)                            , INTENT(IN)    :: tthmax ! 2Theta max
REAL(KIND=PREC_DP)                            , INTENT(IN)    :: dtth   ! 2Theta step
REAL(KIND=PREC_DP)                            , INTENT(IN)    :: eta0   ! Lor/Gaus mix constant part
REAL(KIND=PREC_DP)                            , INTENT(IN)    :: eta_l   ! Lor/Gaus mix variable part
REAL(KIND=PREC_DP)                            , INTENT(IN)    :: eta_q   ! Lor/Gaus mix variable part
REAL(KIND=PREC_DP)                            , INTENT(IN)    :: u      ! u*tan^2(Theta)
REAL(KIND=PREC_DP)                            , INTENT(IN)    :: v      ! v*tan  (Theta)
REAL(KIND=PREC_DP)                            , INTENT(IN)    :: w      ! w
REAL(KIND=PREC_DP)                            , INTENT(IN)    :: pow_width ! Number of FWHM's to calculate
INTEGER                         , INTENT(IN)    :: pow_type  ! == 0 for COMPLETE, ==1 for Debye
INTEGER                         , INTENT(IN)    :: axis   ! == 2 for 2theta, == 1 for Q
REAL(KIND=PREC_DP)                            , INTENT(IN)    :: rlambda ! Wavelength
INTEGER                                       , INTENT(IN)    :: pow_pr_fwhm   ! == 1 for Cagliotti, 2 for area 
!
INTEGER, PARAMETER  :: POW_COMPL = 0
INTEGER, PARAMETER  :: POW_DEBYE = 1
INTEGER, PARAMETER  :: POW_NUFFT = 2
INTEGER, PARAMETER  :: POW_GRID  = 3
INTEGER, PARAMETER  :: POW_PROFILE_CAGLIOTTI  = 1
INTEGER, PARAMETER  :: POW_PROFILE_AREA       = 2
!
REAL(KIND=PREC_DP)            :: fwhm     ! Current FWHM at Theta
REAL(kind=PREC_DP), DIMENSION(0:POW_MAXPKT) :: dummy    ! temporary data (0:POW_MAXPKT) 
REAL(KIND=PREC_DP)            :: theta    ! Real angle Theta for Q-axis
REAL(KIND=PREC_DP)            :: tth      ! 2Theta within convolution, main data set
REAL(kind=PREC_DP)            :: tantth   ! tan(Theta)
REAL(kind=PREC_DP)            :: costth   ! cos(2Theta)
REAL(KIND=PREC_DP)    :: eta              ! actual eta at current 2Theta
REAL(KIND=PREC_DP)    :: ddtth            ! actual stepwidth  at current 2Theta
INTEGER :: imax, i, j, ii  ! Dummy loop indices
INTEGER :: i1, i2          ! Pseudo Voigt lookup indices
REAL(KIND=PREC_DP)    :: pseudo          ! scale factor for lookup table
INTEGER :: max_ps 
!write(*,*) ' CONVOLUTION WITH PSEUDO      ', maxval(dat), pow_type==POW_COMPL, imax
!                                                                       
!------ Now convolute                                                   
!                                                                       
ddtth = 0.001D0
imax = INT( (tthmax - tthmin) / dtth ) - 1
!
dummy = 0.0D0   ! dummy(:)

IF(pow_type==POW_COMPL .or. pow_type==POW_NUFFT .or. pow_type==POW_GRID) THEN     ! Complete, check for zeros in DAT
   main_pts: DO i = 0, imax 
!
      fwhm = 0.00001D0
!
      tth = tthmin + i * ddtth       ! This is tth(axis=2) or Q(axis=1)
      tantth = tand (tth * 0.5D0)    ! tan(Theta) = tan(2Theta/2)
      IF(axis==2 ) THEN              ! 2Theta axis
         if(pow_pr_fwhm==POW_PROFILE_CAGLIOTTI) then
            fwhm = SQRT (MAX (ABS (u * tantth**2 + v * tantth + w), 0.00001D0) ) 
         elseif(pow_pr_fwhm==POW_PROFILE_AREA) then
            costth = cosd(tth)
            fwhm = sqrt(max(abs(u*costth**4 + v*costth**2 + w), 0.00001D0))
         endif
      ELSEif(axis==1) THEN
         if(pow_pr_fwhm==POW_PROFILE_CAGLIOTTI) then
            theta  = asind(tth*rlambda*0.2500D0/PI) ! Theta = asin(Q*lamnda/4pi)
            tantth = tand (theta                  ) ! tan(Theta) = tan( asin(Q*lamnda/4pi)
            fwhm = SQRT( MAX( ABS( u*tantth**2 + v*tantth + w), 0.00001D0) ) * &
                   PI**2/90.0D0/rlambda*cosd(theta)
         elseif(pow_pr_fwhm==POW_PROFILE_AREA) then
            costth = cosd(tth)
            fwhm = sqrt(max(abs(u*costth**4 + v*costth**2 + w), 0.00001D0)) * &
                   PI**2/90.0D0/rlambda*cosd(theta)
         endif
      ENDIF
!
      max_ps = INT((pow_width * fwhm) / ddtth )
      eta = MIN(1.0D0, MAX(0.0D0, eta0 + eta_l * tth + eta_q*tth**2) ) 
!write(*,'(8(f8.5,1x))') eta, fwhm, eta0, eta_l, eta_q, u*tantth**2 + v*tantth + w, tth, theta
      pseudo =     ddtth/fwhm*glp_npt  ! scale factor for look up table
      i1 = 0                               ! == 0 * ddtth
      i2 = MIN(INT(2*i*pseudo), GLP_MAX)   ! == 2*i*ddtth
      dummy (i) = dat (i) * ( glp_pseud_indx(i1, eta, fwhm)  &
                             -glp_pseud_indx(i2, eta, fwhm))                             
!
      ii = MAX (i - 1 - max_ps + 1, 0) 
      first:DO j = ii, i !- 1 
         IF(dat(j)==0.0) CYCLE first
         i1 = MIN(INT((i - j) * pseudo), GLP_MAX)  ! == tth1 = (i - j) * ddtth
         i2 = MIN(INT((i + j) * pseudo), GLP_MAX)  ! == tth2 = (i + j) * ddtth
!write(*,'(8(f11.5,1x))') eta, fwhm, 1.*i1, 1.*i2, glp_pseud_indx(i1, eta, fwhm), glp_pseud_indx(i2, eta, fwhm)
         dummy(i) = dummy(i) + dat(j) *( glp_pseud_indx(i1, eta, fwhm)  &
                                        -glp_pseud_indx(i2, eta, fwhm))                    
      ENDDO first
!
      ii = MIN(i + 1 + max_ps - 1, imax) 
      secnd: DO j = i + 1, ii 
         IF(dat(j)==0.0) CYCLE secnd
         i1 = MIN(INT((j - i) * pseudo), GLP_MAX)  ! == tth1 = (j - i) * ddtth
         i2 = MIN(INT((j + i) * pseudo), GLP_MAX)  ! == tth1 = (j + i) * ddtth
         dummy(i) = dummy(i) + dat(j) *( glp_pseud_indx(i1, eta, fwhm)  &
                                        -glp_pseud_indx(i2, eta, fwhm))
      ENDDO secnd
   ENDDO main_pts
ELSEIF(pow_type==POW_DEBYE) THEN     ! DEBYE, do not check for zeros in DAT
   main_pts_deb: DO i = 0, imax 
!
      fwhm = 0.00001D0
!
      tth = tthmin + i * ddtth       ! This is tth(axis=2) or Q(axis=1)
      tantth = tand (tth * 0.5D0)    ! tan(Theta) = tan(2Theta/2)
      IF(axis==2 ) THEN              ! 2Theta axis
         if(pow_pr_fwhm==POW_PROFILE_CAGLIOTTI) then
            fwhm = SQRT (MAX (ABS (u * tantth**2 + v * tantth + w), 0.00001D0) ) 
         elseif(pow_pr_fwhm==POW_PROFILE_AREA) then
            costth = cosd(tth)
            fwhm = sqrt(max(abs(u*costth**4 + v*costth**2 + w), 0.00001D0))
         endif
      ELSEif(axis==1) THEN
         if(pow_pr_fwhm==POW_PROFILE_CAGLIOTTI) then
            theta  = asind(tth*rlambda*0.2500D0/PI) ! Theta = asin(Q*lamnda/4pi)
            tantth = tand (theta                  ) ! tan(Theta) = tan( asin(Q*lamnda/4pi)
            fwhm = SQRT( MAX( ABS( u*tantth**2 + v*tantth + w), 0.00001D0) ) * &
                   PI**2/90.0D0/rlambda*cosd(theta)
         elseif(pow_pr_fwhm==POW_PROFILE_AREA) then
            costth = cosd(tth)
            fwhm = sqrt(max(abs(u*costth**4 + v*costth**2 + w), 0.00001D0)) * &
                   PI**2/90.0D0/rlambda*cosd(theta)
         endif
      ENDIF
!      tth = tthmin + i * ddtth   ! This is tth(axis=2) or Q(axis=1)
!      tantth = tand (tth * 0.5D0) 
!!
!      fwhm = 0.00001D0
!      IF(axis==2 ) THEN       ! 2Theta axis
!         fwhm = SQRT (MAX (ABS (u * tantth**2 + v * tantth + w), 0.00001D0) ) 
!      ELSEif(axis==1) THEN
!         fwhm = SQRT( MAX( ABS( u*tth**2 + v*tth + w), 0.00001D0) )
!      ENDIF
!
      max_ps = INT((pow_width * fwhm) / ddtth )
      eta = MIN(1.0D0, MAX(0.0D0, eta0 + eta_l * tth + eta_q*tth**2) ) 
      pseudo =     ddtth/fwhm*glp_npt  ! scale factor for look up table
      i1 = 0                               ! == 0 * ddtth
      i2 = MIN(INT(2*i*pseudo), GLP_MAX)   ! == 2*i*ddtth
      dummy (i) = dat (i) * ( glp_pseud_indx(i1, eta, fwhm)  &
                             -glp_pseud_indx(i2, eta, fwhm))                             
!
      ii = MAX (i - 1 - max_ps + 1, 0) 
      first_deb:DO j = ii, i - 1 
         i1 = MIN(INT((i - j) * pseudo), GLP_MAX)  ! == tth1 = (i - j) * ddtth
         i2 = MIN(INT((i + j) * pseudo), GLP_MAX)  ! == tth2 = (i + j) * ddtth
         dummy(i) = dummy(i) + dat(j) *( glp_pseud_indx(i1, eta, fwhm)  &
                                        -glp_pseud_indx(i2, eta, fwhm))                    
      ENDDO first_deb
!
      ii = MIN(i + 1 + max_ps - 1, imax) 
      secnd_deb: DO j = i + 1, ii 
         i1 = MIN(INT((j - i) * pseudo), GLP_MAX)  ! == tth1 = (j - i) * ddtth
         i2 = MIN(INT((j + i) * pseudo), GLP_MAX)  ! == tth1 = (j + i) * ddtth
         dummy(i) = dummy(i) + dat(j) *( glp_pseud_indx(i1, eta, fwhm)  &
                                        -glp_pseud_indx(i2, eta, fwhm))
      ENDDO secnd_deb
   ENDDO main_pts_deb
ENDIF
!                                                                       
DO i = 0, imax 
   dat (i) = dummy (i) * ddtth   ! scale with stepwidth
ENDDO 
!                                                                       
END SUBROUTINE alt_powder_conv_psvgt_uvw          
!
!*****7*****************************************************************
!
SUBROUTINE powder_conv_psvgt_uvw (dat, tthmin, tthmax, dtth, eta0,&
      eta_l, eta_q, u, v, w, pow_width, POW_MAXPKT, pow_type, axis, &
      rlambda, pow_pr_fwhm )
!-
!     Convolute powder pattern with resolution function (Pseudo-Voigt)  
!     FWHM according to caglioti equation, Constant eta                 
!     FWHM = sqrt ( U*tan**2(Theta) + V*tan(Theta) + W)  IF 2Theta scale
!     FWHM = sqrt ( U*Q^2  +  V*Q  +  W               )  IF Q      scale
!     Symmetric case
!+
USE discus_config_mod 
!
USE gauss_lorentz_pseudo_mod
use precision_mod                                                                       
USE trig_degree_mod
USE wink_mod
!
IMPLICIT none 
!
INTEGER, INTENT(IN) :: POW_MAXPKT
REAL(KIND=PREC_DP)   , DIMENSION(0:POW_MAXPKT), INTENT(INOUT) :: dat   ! Data to be convoluted
REAL(KIND=PREC_DP)                            , INTENT(IN)    :: tthmin ! 2Theta min
REAL(KIND=PREC_DP)                            , INTENT(IN)    :: tthmax ! 2Theta max
REAL(KIND=PREC_DP)                            , INTENT(IN)    :: dtth   ! 2Theta step
REAL(KIND=PREC_DP)                            , INTENT(IN)    :: eta0   ! Lor/Gaus mix constant part
REAL(KIND=PREC_DP)                            , INTENT(IN)    :: eta_l   ! Lor/Gaus mix variable part
REAL(KIND=PREC_DP)                            , INTENT(IN)    :: eta_q   ! Lor/Gaus mix variable part
REAL(KIND=PREC_DP)                            , INTENT(IN)    :: u      ! u*tan^2(Theta)
REAL(KIND=PREC_DP)                            , INTENT(IN)    :: v      ! v*tan  (Theta)
REAL(KIND=PREC_DP)                            , INTENT(IN)    :: w      ! w
REAL(KIND=PREC_DP)                            , INTENT(IN)    :: pow_width ! Number of FWHM's to calculate
INTEGER                         , INTENT(IN)    :: pow_type  ! == 0 for COMPLETE, ==1 for Debye
INTEGER                         , INTENT(IN)    :: axis   ! == 2 for 2theta, == 1 for Q
REAL(KIND=PREC_DP)                            , INTENT(IN)    :: rlambda ! Wavelength
INTEGER                                       , INTENT(IN)    :: pow_pr_fwhm   ! == 1 for Cagliotti, 2 for area 
!
INTEGER, PARAMETER  :: POW_COMPL = 0
INTEGER, PARAMETER  :: POW_DEBYE = 1
INTEGER, PARAMETER  :: POW_NUFFT = 2
INTEGER, PARAMETER  :: POW_GRID  = 3
!INTEGER, PARAMETER  :: POW_PROFILE_CAGLIOTTI  = 1
!INTEGER, PARAMETER  :: POW_PROFILE_AREA       = 2
!
REAL(KIND=PREC_DP)            :: fwhm     ! Current FWHM at Theta
REAL(kind=PREC_DP), DIMENSION(0:POW_MAXPKT) :: dummy    ! temporary data (0:POW_MAXPKT) 
!REAL(KIND=PREC_DP)            :: theta    ! Real angle Theta for Q-axis
REAL(KIND=PREC_DP)            :: tth      ! 2Theta within convolution, main data set
!REAL(kind=PREC_DP)            :: tantth   ! tan(Theta)
!REAL(kind=PREC_DP)            :: costth   ! cos(2Theta)
REAL(KIND=PREC_DP)    :: eta              ! actual eta at current 2Theta
INTEGER :: imax, i, j, ii  ! Dummy loop indices
INTEGER :: i1, i2          ! Pseudo Voigt lookup indices
REAL(KIND=PREC_DP)    :: pseudo          ! scale factor for lookup table
INTEGER :: max_ps 
!write(*,*) ' CONVOLUTION WITH PSEUDO      ', maxval(dat), pow_type==POW_COMPL, imax
!                                                                       
!------ Now convolute                                                   
!                                                                       
imax = INT( (tthmax - tthmin) / dtth ) - 1
!
dummy = 0.0D0   ! dummy(:)

IF(pow_type==POW_COMPL .or. pow_type==POW_NUFFT .or. pow_type==POW_GRID) THEN     ! Complete, check for zeros in DAT
   main_pts: DO i = 0, imax 
      if(dat(i) == 0.0D0) cycle main_pts
!
      tth = tthmin + i * dtth        ! This is tth(axis=2) or Q(axis=1)
      fwhm = powder_calc_fwhm_symm(i, tth, axis, u,v,w, rlambda, pow_pr_fwhm)
      eta  = powder_calc_eta      (i, tth, axis, eta_q,eta_l,eta0, rlambda, pow_pr_fwhm, 0.0_PREC_DP, 1.0_PREC_DP)
!
!     eta    = MIN(1.0D0, MAX(0.0D0, eta0 + eta_l * tth + eta_q*tth**2) ) 
      pseudo =     dtth/fwhm*glp_npt  ! scale factor for look up table
      max_ps = min(INT((pow_width * fwhm) / dtth ), int(GLP_MAX/pseudo))
      i1 = max(0, i-max_ps)
      i2 = min(   i+max_ps, imax)
      first: do j = i1, i2
         ii = abs(j-i)*nint(pseudo)
         if(ii<=GLP_MAX) dummy(j) = dummy(j) + dat(i) * glp_pseud_indx(ii, eta, fwhm)
      enddo first
   ENDDO main_pts
ELSEIF(pow_type==POW_DEBYE) THEN     ! DEBYE, do not check for zeros in DAT
!write(*,*) ' powder_conv_psvgt_uvw', u,v,w, imax
i=imax 
i= 0
tth = tthmin + i * dtth
fwhm = powder_calc_fwhm_symm(i, tth, axis, u,v,w, rlambda, pow_pr_fwhm)
!write(*,*) ' tth ', tth, dtth, fwhm
!write(*,*) powder_calc_fwhm_symm(i, tth, axis, u,v,w, rlambda, pow_pr_fwhm)
      eta  = powder_calc_eta      (i, tth, axis, eta_q,eta_l,eta0, rlambda, pow_pr_fwhm, 0.0_PREC_DP, 1.0_PREC_DP)
      pseudo =     dtth/fwhm*glp_npt  ! scale factor for look up table
      max_ps = min(INT((pow_width * fwhm) / dtth ), int(GLP_MAX/pseudo))
      i1 = max(0, i-max_ps)
      i2 = min(   i+max_ps, imax)
!write(*,*) ' tth ', tth, eta, pseudo, max_ps, i1, i2
!write(*,*) ' DONE '
   main_pts_deb: DO i = 0, imax 
!
      tth = tthmin + i * dtth        ! This is tth(axis=2) or Q(axis=1)
      fwhm = powder_calc_fwhm_symm(i, tth, axis, u,v,w, rlambda, pow_pr_fwhm)
!
!     eta    = MIN(1.0D0, MAX(0.0D0, eta0 + eta_l * tth + eta_q*tth**2) ) 
      eta  = powder_calc_eta      (i, tth, axis, eta_q,eta_l,eta0, rlambda, pow_pr_fwhm, 0.0_PREC_DP, 1.0_PREC_DP)
      pseudo =     dtth/fwhm*glp_npt  ! scale factor for look up table
      max_ps = min(INT((pow_width * fwhm) / dtth ), int(GLP_MAX/pseudo))
      i1 = max(0, abs(i-max_ps + 1))
      i2 = min(   i+max_ps-1, imax)
      first_deb: do j = i1, i2
         ii = abs(j-i)*nint(pseudo)
!if(ii>40000) then
!write(*,*) ' i   ', i, tth, fwhm, eta, pseudo
!write(*,*) ' IND ', i1, i2, j, ii, j - i
!endif
         if(ii<=GLP_MAX) dummy(j) = dummy(j) + dat(i) * glp_pseud_indx(ii, eta, fwhm)
      enddo first_deb
   ENDDO main_pts_deb
ENDIF
!                                                                       
DO i = 0, imax 
   dat (i) = dummy (i) * dtth   ! scale with stepwidth
ENDDO 
!                                                                       
END SUBROUTINE powder_conv_psvgt_uvw          
!
!*****7*****************************************************************
!
real(kind=PREC_DP) function powder_calc_fwhm_symm(i, tth, axis, u,v,w,  &
     rlambda, pow_pr_fwhm) result(fwhm)
!-
!  Calculates the FWHM as Cagliotti or as are detector equation
!+
use gauss_lorentz_pseudo_mod
use precision_mod                                                                       
use trig_degree_mod
use wink_mod
!
integer           , intent(in) :: i
REAL(KIND=PREC_DP), INTENT(IN) :: tth    ! 2Theta / Q
INTEGER           , INTENT(IN) :: axis   ! == 2 for 2theta, == 1 for Q
REAL(KIND=PREC_DP), INTENT(IN) :: rlambda ! Wavelength
INTEGER           , INTENT(IN) :: pow_pr_fwhm   ! == 1 for Cagliotti, 2 for area 
REAL(KIND=PREC_DP), INTENT(IN) :: u      ! u*tan^2(Theta)
REAL(KIND=PREC_DP), INTENT(IN) :: v      ! v*tan  (Theta)
REAL(KIND=PREC_DP), INTENT(IN) :: w      ! w
!
INTEGER, PARAMETER  :: POW_PROFILE_CAGLIOTTI  = 1
INTEGER, PARAMETER  :: POW_PROFILE_AREA       = 2
INTEGER, PARAMETER  :: POW_PROFILE_POLY       = 3
!
REAL(KIND=PREC_DP)            :: theta    ! Real angle Theta for Q-axis
!REAL(KIND=PREC_DP)            :: tth      ! 2Theta within convolution, main data set
REAL(kind=PREC_DP)            :: tantth   ! tan(Theta)
REAL(kind=PREC_DP)            :: costth   ! cos(2Theta)
!
fwhm = 0.00001D0
!
tantth = tand (tth * 0.5D0)    ! tan(Theta) = tan(2Theta/2)
IF(axis==2 ) THEN              ! 2Theta axis
   if(pow_pr_fwhm==POW_PROFILE_CAGLIOTTI) then
      fwhm = SQRT (MAX (ABS (u * tantth**2 + v * tantth + w), 0.00001D0) ) 
   elseif(pow_pr_fwhm==POW_PROFILE_AREA) then
      costth = cosd(tth)
      fwhm = sqrt(max(abs(u*costth**4 + v*costth**2 + w), 0.00001D0))
   elseif(pow_pr_fwhm==POW_PROFILE_POLY) then
      fwhm = SQRT (MAX (ABS (u * tth**2 + v * tth + w), 0.00001D0) ) 
   endif
ELSEif(axis==1) THEN
   if(pow_pr_fwhm==POW_PROFILE_CAGLIOTTI) then
      theta  = asind(tth*rlambda*0.2500D0/PI) ! Theta = asin(Q*lambda/4pi)
      tantth = tand (theta                  ) ! tan(Theta) = tan( asin(Q*lamnda/4pi)
      fwhm = SQRT( MAX( ABS( u*tantth**2 + v*tantth + w), 0.00001D0) ) * &
             PI**2/90.0D0/rlambda*cosd(theta)
   elseif(pow_pr_fwhm==POW_PROFILE_AREA) then
      theta  = asind(tth*rlambda*0.2500D0/PI) ! Theta = asin(Q*lambda/4pi)
      costth = cosd(2.0D0*theta)
      fwhm = sqrt(max(abs(u*costth**4 + v*costth**2 + w), 0.00001D0)) * &
             PI**2/90.0D0/rlambda*cosd(theta)
   elseif(pow_pr_fwhm==POW_PROFILE_POLY) then
      fwhm = SQRT (MAX (ABS (u * tth**2 + v * tth + w), 0.00001D0) ) 
   endif
ENDIF
!
end function powder_calc_fwhm_symm
!
!*****7*****************************************************************
!
real(kind=PREC_DP) function powder_calc_eta      (i, tth, axis, eta_q,eta_l,eta_0,  &
     rlambda, pow_pr_fwhm, eta_min, eta_max) result(eta)
!-
!  Calculates the FWHM as Cagliotti or as are detector equation
!+
use gauss_lorentz_pseudo_mod
use precision_mod                                                                       
use trig_degree_mod
use wink_mod
!
integer           , intent(in) :: i
REAL(KIND=PREC_DP), INTENT(IN) :: tth    ! 2Theta / Q
INTEGER           , INTENT(IN) :: axis   ! == 2 for 2theta, == 1 for Q
REAL(KIND=PREC_DP), INTENT(IN) :: rlambda ! Wavelength
INTEGER           , INTENT(IN) :: pow_pr_fwhm   ! == 1 for Cagliotti, 2 for area 
REAL(KIND=PREC_DP), INTENT(IN) :: eta_q  ! eta proportional to tth**2
REAL(KIND=PREC_DP), INTENT(IN) :: eta_l  ! eta proportional to tth**1
REAL(KIND=PREC_DP), INTENT(IN) :: eta_0  ! eta proportional to tth**0
real(kind=PREC_DP), intent(in) :: eta_min  ! Minimum eta for PSVOIGT = r0 for Pearson 0.5 
real(kind=PREC_DP), intent(in) :: eta_max  ! Maximum eta for PSVOIGT = 1, for Pearson 1000
!
INTEGER, PARAMETER  :: POW_PROFILE_CAGLIOTTI  = 1
INTEGER, PARAMETER  :: POW_PROFILE_AREA       = 2
INTEGER, PARAMETER  :: POW_PROFILE_POLY       = 3
!
REAL(KIND=PREC_DP)            :: ttheta    ! Real angle Theta for Q-axis
!REAL(KIND=PREC_DP)            :: tth      ! 2Theta within convolution, main data set
!REAL(kind=PREC_DP)            :: tantth   ! tan(Theta)
!REAL(kind=PREC_DP)            :: costth   ! cos(2Theta)
!
eta = 0.00000D0
!
IF(axis==2 ) THEN              ! 2Theta axis
   if(pow_pr_fwhm==POW_PROFILE_CAGLIOTTI) then
      eta = min(MAX (ABS (eta_q * tth**2 + eta_l * tth + eta_0), eta_min  ) ,eta_max)
   elseif(pow_pr_fwhm==POW_PROFILE_AREA) then
      eta = min(MAX (ABS (eta_q * tth**2 + eta_l * tth + eta_0), eta_min  ) ,eta_max)
   elseif(pow_pr_fwhm==POW_PROFILE_POLY) then
      eta = min(MAX (ABS (eta_q * tth**2 + eta_l * tth + eta_0), eta_min  ) ,eta_max)
   endif
ELSEif(axis==1) THEN
   if(pow_pr_fwhm==POW_PROFILE_CAGLIOTTI) then
      ttheta  = 2.0D0*asind(tth*rlambda*0.2500D0/PI) ! Theta = asin(Q*lambda/4pi)
      eta = min(MAX (ABS (eta_q * ttheta**2 + eta_l * ttheta + eta_0), eta_min  ) ,eta_max)
   elseif(pow_pr_fwhm==POW_PROFILE_AREA) then
      ttheta  = 2.0D0*asind(tth*rlambda*0.2500D0/PI) ! Theta = asin(Q*lambda/4pi)
      eta = min(MAX (ABS (eta_q * ttheta**2 + eta_l * ttheta + eta_0), eta_min  ) ,eta_max)
   elseif(pow_pr_fwhm==POW_PROFILE_POLY) then
      eta = min(MAX (ABS (eta_q * tth**2 + eta_l * tth + eta_0), eta_min  ) ,eta_max)
   endif
ENDIF
!
end function powder_calc_eta
!
!*****7*****************************************************************
!
SUBROUTINE powder_conv_psvgt_uvw_asym (dat, tthmin, tthmax, dtth, eta0,&
      eta_l, eta_q, u, v, w, asym          , pow_width, POW_MAXPKT,    &
      pow_type, axis, rlambda, pow_pr_fwhm)
!-
!     Convolute powder pattern with resolution function (Pseudo-Voigt)  
!     FWHM according to caglioti equation, Constant eta                 
!     FWHM = sqrt ( U*tan**2(Theta) + V*tan(Theta) + W)                 
!     Symmetric case
!+
USE discus_config_mod 
!
USE gauss_lorentz_pseudo_mod
use precision_mod                                                                       
USE trig_degree_mod
USE wink_mod
!
IMPLICIT none 
!
INTEGER, INTENT(IN) :: POW_MAXPKT
REAL(KIND=PREC_DP)   , DIMENSION(0:POW_MAXPKT), INTENT(INOUT) :: dat   ! Data to be convoluted
REAL(kind=PREC_DP)              , INTENT(IN)    :: tthmin ! 2Theta min
REAL(kind=PREC_DP)              , INTENT(IN)    :: tthmax ! 2Theta max
REAL(kind=PREC_DP)              , INTENT(IN)    :: dtth   ! 2Theta step
REAL(kind=PREC_DP)              , INTENT(IN)    :: eta0   ! Lor/Gaus mix constant part
REAL(kind=PREC_DP)              , INTENT(IN)    :: eta_l   ! Lor/Gaus mix variable part
REAL(kind=PREC_DP)              , INTENT(IN)    :: eta_q   ! Lor/Gaus mix variable part
REAL(kind=PREC_DP)              , INTENT(IN)    :: u      ! u*tan^2(Theta)
REAL(kind=PREC_DP)              , INTENT(IN)    :: v      ! v*tan  (Theta)
REAL(kind=PREC_DP)              , INTENT(IN)    :: w      ! w
REAL(kind=PREC_DP)   , DIMENSION(2,-1:2)      , INTENT(IN)    :: asym   ! Asymmetry terms p1 to P4
REAL(kind=PREC_DP)              , INTENT(IN)    :: pow_width ! Number of FWHM's to calculate
INTEGER                         , INTENT(IN)    :: pow_type  ! == 0 for COMPLETE, ==1 for Debye
INTEGER                         , INTENT(IN)    :: axis   ! == 2 for 2theta, == 1 for Q
REAL(KIND=PREC_DP)              , INTENT(IN)    :: rlambda ! Wavelength
INTEGER                         , INTENT(IN)    :: pow_pr_fwhm   ! == 1 for Cagliotti, 2 for area 
!
INTEGER, PARAMETER  :: POW_COMPL = 0
INTEGER, PARAMETER  :: POW_DEBYE = 1
INTEGER, PARAMETER  :: POW_NUFFT = 2
INTEGER, PARAMETER  :: POW_GRID  = 3
!INTEGER, PARAMETER  :: POW_PROFILE_CAGLIOTTI  = 1
!INTEGER, PARAMETER  :: POW_PROFILE_AREA       = 2
!
REAL(KIND=PREC_DP)            :: fwhm     ! Current FWHM at Theta
REAL(kind=PREC_DP), DIMENSION(0:POW_MAXPKT) :: dummy    ! temporary data (0:POW_MAXPKT) 
REAL(KIND=PREC_DP)            :: tth      ! Theta within convolution, main data set
!REAL(KIND=PREC_DP)            :: theta    ! Real angle Theta for Q-axis
!REAL(kind=PREC_DP)            :: tantth   ! tan(Theta)
!REAL(kind=PREC_DP)            :: costth   ! cos(2Theta)
REAL(KIND=PREC_DP)    :: eta              ! actual eta at current 2Theta
!REAL(KIND=PREC_DP)    :: ddtth            ! actual stepwidth  at current 2Theta
REAL(KIND=PREC_DP)    :: tth1            ! Theta values for asymmetry
REAL(KIND=PREC_DP)    :: tth2            ! Theta values for asymmetry
REAL(kind=PREC_DP)    :: p1, p2          ! 2Theta dependen asymmety parameters
!REAL(kind=PREC_DP)    :: p1, p2, p3, p4  ! 2Theta dependen asymmety parameters
INTEGER :: imax, i, j, ii  ! Dummy loop indices
INTEGER :: i1, i2          ! Pseudo Voigt lookup indices
REAL(kind=PREC_DP)    :: pseudo          ! scale factor for lookup table
REAL(kind=PREC_DP)    :: pra1!, pra2      ! Asymmetry pre factors
INTEGER :: max_ps 
!                                                                       
!------ Now convolute                                                   
!                                                                       
imax = INT( (tthmax - tthmin) / dtth )
!write(*,*) ' IMAX ', imax, tthmin, tthmax, dtth
dummy = 0.0D0
!
IF(pow_type==POW_COMPL .or. pow_type==POW_NUFFT .or. pow_type==POW_GRID) THEN
   main_compl: DO i = 0, imax 
!
      if(dat(i)==0.0D0) cycle main_compl
!
      tth = tthmin + i * dtth        ! This is tth(axis=2) or Q(axis=1)
      fwhm = powder_calc_fwhm_symm(i, tth, axis, u,v,w, rlambda, pow_pr_fwhm)
!
      eta  = powder_calc_eta      (i, tth, axis, eta_q,eta_l,eta0, rlambda, pow_pr_fwhm, 0.0_PREC_DP, 1.0_PREC_DP)
      pseudo =     dtth/fwhm*glp_npt  ! scale factor for look up table
      max_ps = min(INT((pow_width * fwhm) / dtth ), int(GLP_MAX/pseudo))
      tth1 = 0 * dtth 
      tth2 = 2 * i * dtth 
      p1 = asym(1,0) + asym(1,1)*tth +     asym(1,2)*tth**2 + asym(1,-1)/tth
      p2 = asym(2,0) + asym(2,1)*tth +     asym(2,2)*tth**2 + asym(2,-1)/tth
!     p3 = asym(3,0) + asym(3,1)*tth +     asym(3,2)*tth**2 + asym(3,-1)/tth
!     p4 = asym(4,0) + asym(4,1)*tth +     asym(4,2)*tth**2 + asym(4,-1)/tth
      pra1 = profile_asymmetry (tth, tth1, fwhm, p1, p2) !, p3, p4) 
      i1 = max(0, i-max_ps)
      i2 = min(   i+max_ps, imax)
      first: do j = i1, i2
         ii = abs(j-i)*nint(pseudo)
         tth1 = (j-i-1)*dtth
         pra1 = profile_asymmetry (tth, tth1, fwhm, p1, p2) !, p3, p4) 
         dummy(j) = dummy(j) + dat(i) * glp_pseud_indx(ii, eta, fwhm)*pra1
      enddo first
   ENDDO main_compl
ELSEIF(pow_type==POW_DEBYE) THEN
   main_debye: DO i = 0, imax 
!
      tth = tthmin + i * dtth        ! This is tth(axis=2) or Q(axis=1)
      fwhm = powder_calc_fwhm_symm(i, tth, axis, u,v,w, rlambda, pow_pr_fwhm)
!
      eta  = powder_calc_eta      (i, tth, axis, eta_q,eta_l,eta0, rlambda, pow_pr_fwhm, 0.0_PREC_DP, 1.0_PREC_DP)
      pseudo =     dtth/fwhm*glp_npt  ! scale factor for look up table
      max_ps = min(INT((pow_width * fwhm) / dtth ), int(GLP_MAX/pseudo))
      tth1 = 0 * dtth 
      tth2 = 2 * i * dtth 
      p1 = asym(1,0) + asym(1,1)*tth +     asym(1,2)*tth**2 + asym(1,-1)/tth
      p2 = asym(2,0) + asym(2,1)*tth +     asym(2,2)*tth**2 + asym(2,-1)/tth
!     p3 = asym(3,0) + asym(3,1)*tth +     asym(3,2)*tth**2 + asym(3,-1)/tth
!     p4 = asym(4,0) + asym(4,1)*tth +     asym(4,2)*tth**2 + asym(4,-1)/tth
      i1 = max(0, i-max_ps)
      i2 = min(   i+max_ps, imax)
      first_deb: do j = i1, i2
         ii = abs(j-i)*nint(pseudo)
         tth1 = (j-i-1)*dtth
         pra1 = profile_asymmetry (tth, tth1, fwhm, p1, p2) !, p3, p4) 
         dummy(j) = dummy(j) + dat(i) * glp_pseud_indx(ii, eta, fwhm)*pra1
      enddo first_deb
   ENDDO main_debye
ENDIF
!                                                                       
DO i = 0, imax 
   dat (i) = dummy (i) * dtth   ! scale with stepwidth
ENDDO 
!                                                                       
END SUBROUTINE powder_conv_psvgt_uvw_asym
!
!*****7*****************************************************************
!
SUBROUTINE powder_conv_pears_uvw_asym (dat, tthmin, tthmax, dtth, eta0,&
      eta_l, eta_q, u, v, w, asym          , pow_width, POW_MAXPKT,    &
      pow_type, axis, rlambda, pow_pr_fwhm)
!-
!     Convolute powder pattern with resolution function (Pearson VII )  
!     FWHM according to caglioti equation, Constant eta                 
!     FWHM = sqrt ( U*tan**2(Theta) + V*tan(Theta) + W)                 
!     Symmetric case
!+
USE discus_config_mod 
!
use gamma_mod
USE gauss_lorentz_pseudo_mod
use precision_mod                                                                       
USE trig_degree_mod
USE wink_mod
!
IMPLICIT none 
!
INTEGER, INTENT(IN) :: POW_MAXPKT
REAL(KIND=PREC_DP)   , DIMENSION(0:POW_MAXPKT), INTENT(INOUT) :: dat   ! Data to be convoluted
REAL(kind=PREC_DP)              , INTENT(IN)    :: tthmin ! 2Theta min
REAL(kind=PREC_DP)              , INTENT(IN)    :: tthmax ! 2Theta max
REAL(kind=PREC_DP)              , INTENT(IN)    :: dtth   ! 2Theta step
REAL(kind=PREC_DP)              , INTENT(IN)    :: eta0   ! Lor/Gaus mix constant part
REAL(kind=PREC_DP)              , INTENT(IN)    :: eta_l   ! Lor/Gaus mix variable part
REAL(kind=PREC_DP)              , INTENT(IN)    :: eta_q   ! Lor/Gaus mix variable part
REAL(kind=PREC_DP)              , INTENT(IN)    :: u      ! u*tan^2(Theta)
REAL(kind=PREC_DP)              , INTENT(IN)    :: v      ! v*tan  (Theta)
REAL(kind=PREC_DP)              , INTENT(IN)    :: w      ! w
REAL(kind=PREC_DP)   , DIMENSION(2,-1:2)      , INTENT(IN)    :: asym   ! Asymmetry terms p1 to P4
REAL(kind=PREC_DP)              , INTENT(IN)    :: pow_width ! Number of FWHM's to calculate
INTEGER                         , INTENT(IN)    :: pow_type  ! == 0 for COMPLETE, ==1 for Debye
INTEGER                         , INTENT(IN)    :: axis   ! == 2 for 2theta, == 1 for Q
REAL(KIND=PREC_DP)              , INTENT(IN)    :: rlambda ! Wavelength
INTEGER                         , INTENT(IN)    :: pow_pr_fwhm   ! == 1 for Cagliotti, 2 for area 
!
INTEGER, PARAMETER  :: POW_COMPL = 0
INTEGER, PARAMETER  :: POW_DEBYE = 1
INTEGER, PARAMETER  :: POW_NUFFT = 2
INTEGER, PARAMETER  :: POW_GRID  = 3
!INTEGER, PARAMETER  :: POW_PROFILE_CAGLIOTTI  = 1
!INTEGER, PARAMETER  :: POW_PROFILE_AREA       = 2
!
REAL(KIND=PREC_DP)            :: fwhm     ! Current FWHM at Theta
REAL(kind=PREC_DP), DIMENSION(0:POW_MAXPKT) :: dummy    ! temporary data (0:POW_MAXPKT) 
REAL(KIND=PREC_DP)            :: tth      ! Theta within convolution, main data set
!REAL(KIND=PREC_DP)            :: theta    ! Real angle Theta for Q-axis
!REAL(kind=PREC_DP)            :: tantth   ! tan(Theta)
!REAL(kind=PREC_DP)            :: costth   ! cos(2Theta)
REAL(KIND=PREC_DP)    :: eta              ! actual eta at current 2Theta
!REAL(KIND=PREC_DP)    :: ddtth            ! actual stepwidth  at current 2Theta
REAL(KIND=PREC_DP)    :: tth1            ! Theta values for asymmetry
REAL(KIND=PREC_DP)    :: tth2            ! Theta values for asymmetry
REAL(kind=PREC_DP)    :: p1, p2          ! 2Theta dependen asymmety parameters
!REAL(kind=PREC_DP)    :: p1, p2, p3, p4  ! 2Theta dependen asymmety parameters
INTEGER :: imax, i, j, ii  ! Dummy loop indices
INTEGER :: i1, i2          ! Pseudo Voigt lookup indices
REAL(kind=PREC_DP)    :: pseudo          ! scale factor for lookup table
REAL(kind=PREC_DP)    :: pra1!, pra2      ! Asymmetry pre factors
real(kind=PREC_DP)    :: alpha
real(kind=PREC_DP)    :: divisor
INTEGER :: max_ps 
!                                                                       
!------ Now convolute                                                   
!                                                                       
imax = INT( (tthmax - tthmin) / dtth )
dummy = 0.0D0
!
IF(pow_type==POW_COMPL .or. pow_type==POW_NUFFT .or. pow_type==POW_GRID) THEN
   main_compl: DO i = 0, imax 
!
      if(dat(i)==0.0D0) cycle main_compl
!
      tth = tthmin + i * dtth        ! This is tth(axis=2) or Q(axis=1)
      fwhm = powder_calc_fwhm_symm(i, tth, axis, u,v,w, rlambda, pow_pr_fwhm)
!
      eta  = powder_calc_eta      (i, tth, axis, eta_q,eta_l,eta0, rlambda, pow_pr_fwhm, 0.5_PREC_DP, 1000.0_PREC_DP)
      eta   = max(0.50D0, eta)
      alpha = fwhm*0.50D0/sqrt(2.0D0**( 1.0D0/eta) - 1.0D0)
      max_ps =     INT((pow_width * fwhm) / dtth )
      tth1 = 0 * dtth 
      tth2 = 2 * i * dtth 
      p1 = asym(1,0) + asym(1,1)*tth +     asym(1,2)*tth**2 + asym(1,-1)/tth
      p2 = asym(2,0) + asym(2,1)*tth +     asym(2,2)*tth**2 + asym(2,-1)/tth
!     p3 = asym(3,0) + asym(3,1)*tth +     asym(3,2)*tth**2 + asym(3,-1)/tth
!     p4 = asym(4,0) + asym(4,1)*tth +     asym(4,2)*tth**2 + asym(4,-1)/tth
      i1 = max(0, i-max_ps)
      i2 = min(   i+max_ps, imax)
      divisor = (alpha*func_beta(eta-0.5D0, 0.5D0))
!
      first: do j = i1, i2
         ii = abs(j-i)*nint(pseudo)
         tth1 = (j-i-1)*dtth
         pra1 = profile_asymmetry (tth, tth1, fwhm, p1, p2) !, p3, p4) 
         dummy(j) = dummy(j) + dat(i) * pra1 * &
                    ( (1.0D0+((tth1)/alpha)**2)**(-eta) )/ divisor
      enddo first
   ENDDO main_compl
ELSEIF(pow_type==POW_DEBYE) THEN
   main_debye: DO i = 0, imax 
!
      tth = tthmin + i * dtth        ! This is tth(axis=2) or Q(axis=1)
      fwhm = powder_calc_fwhm_symm(i, tth, axis, u,v,w, rlambda, pow_pr_fwhm)
!
      eta  = powder_calc_eta      (i, tth, axis, eta_q,eta_l,eta0, rlambda, pow_pr_fwhm, 0.5_PREC_DP, 1000.0_PREC_DP)
      eta   = max(0.50D0, eta)
      alpha = fwhm*0.50D0/sqrt(2.0D0**( 1.0D0/eta) - 1.0D0)
      max_ps =     INT((pow_width * fwhm) / dtth )
      tth1 = 0 * dtth 
      tth2 = 2 * i * dtth 
      p1 = asym(1,0) + asym(1,1)*tth +     asym(1,2)*tth**2 + asym(1,-1)/tth
      p2 = asym(2,0) + asym(2,1)*tth +     asym(2,2)*tth**2 + asym(2,-1)/tth
!     p3 = asym(3,0) + asym(3,1)*tth +     asym(3,2)*tth**2 + asym(3,-1)/tth
!     p4 = asym(4,0) + asym(4,1)*tth +     asym(4,2)*tth**2 + asym(4,-1)/tth
      i1 = max(0, i-max_ps)
      i2 = min(   i+max_ps, imax)
      divisor = (alpha*func_beta(eta-0.5D0, 0.5D0))
      first_deb: do j = i1, i2
         ii = abs(j-i)*nint(pseudo)
         tth1 = (j-i-1)*dtth
         pra1 = profile_asymmetry(tth, tth1, fwhm, p1, p2) !, p3, p4) 
         dummy(j) = dummy(j) + dat(i) * pra1 * &
                    ( (1.0D0+((tth1)/alpha)**2)**(-eta) )/ divisor
      enddo first_deb
   ENDDO main_debye
ENDIF
!                                                                       
DO i = 0, imax 
   dat (i) = dummy (i) * dtth   ! scale with stepwidth
ENDDO 
!                                                                       
END SUBROUTINE powder_conv_pears_uvw_asym
!
!*****7*****************************************************************
!
SUBROUTINE powder_conv_psvgt_uvw_asymt(dat, tthmin, tthmax, dtth, eta0,&
      eta_l, eta_q, u, v, w, p1, p2, p3, p4, pow_width, POW_MAXPKT)
!-
!     Convolute powder pattern with resolution function (Pseudo-Voigt)  
!     FWHM according to caglioti equation, Constant eta                 
!     FWHM = sqrt ( U*tan**2(Theta) + V*tan(Theta) + W)                 
!+
USE discus_config_mod 
!
USE gauss_lorentz_pseudo_mod
use precision_mod                                                                       
USE trig_degree_mod
USE wink_mod
!
IMPLICIT none 
!
INTEGER, INTENT(IN) :: POW_MAXPKT
REAL(KIND=PREC_DP)   , DIMENSION(0:POW_MAXPKT), INTENT(INOUT) :: dat   ! Data to be convoluted
REAL(KIND=PREC_DP)                            , INTENT(IN)    :: tthmin ! 2Theta min
REAL(KIND=PREC_DP)                            , INTENT(IN)    :: tthmax ! 2Theta max
REAL(KIND=PREC_DP)                            , INTENT(IN)    :: dtth   ! 2Theta step
REAL(KIND=PREC_DP)                            , INTENT(IN)    :: eta0   ! Lor/Gaus mix constant part
REAL(KIND=PREC_DP)                            , INTENT(IN)    :: eta_l   ! Lor/Gaus mix variable part
REAL(KIND=PREC_DP)                            , INTENT(IN)    :: eta_q   ! Lor/Gaus mix variable part
REAL(KIND=PREC_DP)                            , INTENT(IN)    :: u      ! u*tan^2(Theta)
REAL(KIND=PREC_DP)                            , INTENT(IN)    :: v      ! v*tan  (Theta)
REAL(KIND=PREC_DP)                            , INTENT(IN)    :: w      ! w
REAL(KIND=PREC_DP)                            , INTENT(IN)    :: p1     ! Asymmetry terms p1 to P4
REAL(KIND=PREC_DP)                            , INTENT(IN)    :: p2     ! Asymmetry terms p1 to P4
REAL(KIND=PREC_DP)                            , INTENT(IN)    :: p3     ! Asymmetry terms p1 to P4
REAL(KIND=PREC_DP)                            , INTENT(IN)    :: p4     ! Asymmetry terms p1 to P4
REAL(KIND=PREC_DP)                            , INTENT(IN)    :: pow_width ! Number of FWHM's to calculate
!
REAL(KIND=PREC_DP)    :: fwhm            ! Current FWHM at Theta
REAL(kind=PREC_DP), DIMENSION(0:POW_MAXPKT) :: dummy  ! temporary data (0:POW_MAXPKT) 
REAL(KIND=PREC_DP)    :: tth             ! Theta within convolution, main data set
REAL(kind=PREC_DP)    :: tantth          ! tan(Theta)
REAL(KIND=PREC_DP)    :: tth1            ! Theta values for asymmetry
REAL(KIND=PREC_DP)    :: tth2            ! Theta values for asymmetry
REAL(KIND=PREC_DP)    :: eta             ! actual eta at current 2Theta
REAL(kind=PREC_DP)    :: pra1, pra2      ! Asymmetry pre factors
INTEGER :: imax, i, j, ii  ! Dummy loop indices
INTEGER :: max_ps 
!                                                                       
!write(*,*) ' IN CONV_VARIABLE '
pra1 = 1.0D0
pra2 = 1.0D0
!                                                                       
!------ Now convolute                                                   
!                                                                       
imax = INT( (tthmax - tthmin) / dtth )
DO i = 0, imax 
   tth = tthmin + i * dtth 
   tantth = tand (tth * 0.5D0) 
!     atheta = tth * 0.5 
!     atwoth = tth 
   fwhm = SQRT (MAX (ABS (u * tantth**2 + v * tantth + w), 0.00001D0) ) 
!      fwhm1 = fwhm 
   max_ps = INT((pow_width * fwhm) / dtth )
   eta = MIN(1.0D0, MAX(0.0D0, eta0 + eta_l * tth + eta_q*tth**2) ) 
   tth1 = 0 * dtth 
   tth2 = 2 * i * dtth 
   pra1 = profile_asymmetry (tth, tth1, fwhm, p1, p2) !, p3, p4) 
   pra2 = profile_asymmetry (tth, tth2, fwhm, p1, p2) !, p3, p4) 
!  dummy (i) = dat (i) * (pseudovoigt (tth1, eta, fwhm) * pra1 -     &
!                         pseudovoigt (tth2, eta, fwhm) * pra2)                             
   dummy (i) = dat (i) * (glp_pseud   (tth1, eta, fwhm) * pra1 -     &
                          glp_pseud   (tth2, eta, fwhm) * pra2)                             
!       do j=0,i-1                                                      
   ii = MAX (i - 1 - max_ps + 1, 0) 
   DO j = ii, i - 1 
      tth1 = (i - j) * dtth 
      tth2 = (i + j) * dtth 
      pra1 = profile_asymmetry(tth, tth1, fwhm, p1, p2) !, p3, p4) 
      pra2 = profile_asymmetry(tth, tth2, fwhm, p1, p2) !, p3, p4) 
!     dummy(i) = dummy(i) + dat(j) *( pseudovoigt(tth1, eta, fwhm) * pra1  &
!                                    -pseudovoigt(tth2, eta, fwhm) * pra2)                    
      dummy(i) = dummy(i) + dat(j) *( glp_pseud  (tth1, eta, fwhm) * pra1  &
                                     -glp_pseud  (tth2, eta, fwhm) * pra2)                    
   ENDDO 
!       do j=i+1,imax                                                   
   ii = MIN(i + 1 + max_ps - 1, imax) 
   DO j = i + 1, ii 
      tth1 = (j - i) * dtth 
      tth2 = (j + i) * dtth 
      pra1 = profile_asymmetry (tth, - tth1, fwhm, p1, p2) !, p3, p4) 
      pra2 = profile_asymmetry (tth, - tth2, fwhm, p1, p2) !, p3, p4) 
!     dummy(i) = dummy(i) + dat(j) *( pseudovoigt(tth1, eta, fwhm) * pra1  &
!                                    -pseudovoigt(tth2, eta, fwhm) * pra2)
      dummy(i) = dummy(i) + dat(j) *( glp_pseud  (tth1, eta, fwhm) * pra1  &
                                     -glp_pseud  (tth2, eta, fwhm) * pra2)
   ENDDO 
ENDDO 
!                                                                       
DO i = 0, imax 
   dat (i) = dummy (i) * dtth   ! scale with stepwidth
ENDDO 
!                                                                       
END SUBROUTINE powder_conv_psvgt_uvw_asymt          
!*****7*****************************************************************
      SUBROUTINE powder_conv_psvgt_uvw_Qscale (dat, tthmin, tthmax,     &
      dtth, eta0, eta_l, eta_q, u, v, w, p1, p2, p3, p4, pow_width, rlambda,    &
      pow_axis, POW_AXIS_Q, POW_MAXPKT)                                             
!-                                                                      
!     Convolute powder pattern with resolution function (Pseudo-Voigt)  
!     FWHM according to caglioti equation, Constant eta                 
!     FWHM = sqrt ( U*tan**2(Theta) + V*tan(Theta) + W)                 
!+                                                                      
      USE discus_config_mod 
      USE trig_degree_mod
      USE precision_mod
      USE wink_mod
      IMPLICIT none 
!                                                                       
!
      INTEGER, INTENT(IN) :: POW_MAXPKT
!                                                                       
      REAL(KIND=PREC_DP):: dat (0:POW_MAXPKT) 
      REAL(KIND=PREC_DP):: tthmin, tthmax, dtth, eta0, eta_l, eta_q
      REAL(KIND=PREC_DP):: u, v, w 
      REAL(KIND=PREC_DP):: p1, p2, p3, p4 
      REAL(KIND=PREC_DP):: pow_width 
      REAL(KIND=PREC_DP):: rlambda 
      INTEGER pow_axis 
      INTEGER POW_AXIS_Q 
!                                                                       
      REAL(KIND=PREC_DP):: dummy (0:POW_MAXPKT) 
      REAL(KIND=PREC_DP) tth
      REAL(KIND=PREC_DP):: tantth 
      REAL(KIND=PREC_DP) tth1 
      REAL(KIND=PREC_DP) tth2 
      REAL(KIND=PREC_DP):: atheta 
      REAL(KIND=PREC_DP):: atwoth 
      REAL(KIND=PREC_DP) fwhm1 , fwhm
      REAL(KIND=PREC_DP) eta 
      REAL(KIND=PREC_DP):: pra1, pra2 
      INTEGER imax, i, j, ii 
      INTEGER max_ps 
!                                                                       
!------ Now convolute                                                   
!                                                                       
      imax = int( (tthmax - tthmin) / dtth )
      DO i = 0, imax 
      tth = tthmin + i * dtth 
      tantth = tand (tth * 0.5D0) 
      atheta = tth * 0.5D0 
      atwoth = tth 
      fwhm = sqrt (max (abs (u * tantth**2 + v * tantth + w), 0.00001D0) ) 
      fwhm1 = fwhm 
      IF (pow_axis.eq.POW_AXIS_Q) THEN 
         atheta = asind (tth * rlambda / (fpi)) 
         tantth = tand (atheta) 
         fwhm1 = sqrt (max (abs (u * tantth**2 + v * tantth + w), 0.00001D0) )                                                     
         fwhm = 0.500D0 * ((fpi) * sind (atheta + 0.5D0 * fwhm1) / rlambda -  &
                fpi * sind (atheta - 0.5D0 * fwhm1) / rlambda)                   
      ENDIF 
      max_ps = int( (pow_width * fwhm) / dtth )
      eta = min (1.0D0, max (0.0D0, eta0 + eta_l * tth + eta_q*tth**2) ) 
      tth1 = 0 * dtth 
      tth2 = 2 * i * dtth 
      pra1 = profile_asymmetry (tth, tth1, fwhm, p1, p2) !, p3, p4) 
      pra2 = profile_asymmetry (tth, tth2, fwhm, p1, p2) !, p3, p4) 
!     dummy (i) = dat (i) * (pseudovoigt (tth1, eta, fwhm) * pra1 -     &
!     pseudovoigt (tth2, eta, fwhm) * pra2)                             
!       do j=0,i-1                                                      
      ii = max (i - 1 - max_ps + 1, 0) 
      DO j = ii, i - 1 
      tth1 = (i - j) * dtth 
      tth2 = (i + j) * dtth 
      pra1 = profile_asymmetry (tth, tth1, fwhm, p1, p2) !, p3, p4) 
      pra2 = profile_asymmetry (tth, tth2, fwhm, p1, p2) !, p3, p4) 
!     dummy (i) = dummy (i) + dat (j) * (pseudovoigt (tth1, eta, fwhm)  &
!     * pra1 - pseudovoigt (tth2, eta, fwhm) * pra2)                    
      ENDDO 
!       do j=i+1,imax                                                   
      ii = min (i + 1 + max_ps - 1, imax) 
      DO j = i + 1, ii 
      tth1 = (j - i) * dtth 
      tth2 = (j + i) * dtth 
      pra1 = profile_asymmetry (tth, - tth1, fwhm, p1, p2) !, p3, p4) 
      pra2 = profile_asymmetry (tth, - tth2, fwhm, p1, p2) !, p3, p4) 
!     dummy (i) = dummy (i) + dat (j) * (pseudovoigt (tth1, eta, fwhm)  &
!     * pra1 - pseudovoigt (tth2, eta, fwhm) * pra2)                    
      ENDDO 
      ENDDO 
!                                                                       
      DO i = 0, imax 
      dat (i) = dummy (i) ! * dtth 
      ENDDO 
!                                                                       
      END SUBROUTINE powder_conv_psvgt_uvw_Qscale   
!
!*****7*****************************************************************
!
SUBROUTINE powder_conv_corrlin_old(dat, tthmin, tthmax, dtth, sigma2, &
      corrlin, corrquad, rcut, pow_width, POW_MAXPKT)
!-
!     Convolute PDF with Gaussian function 
!     FWHM = (sigma**2 - corrlin/r)
!+
USE discus_config_mod 
!
USE gauss_lorentz_pseudo_mod
use precision_mod                                                                       
USE trig_degree_mod
USE wink_mod
!
IMPLICIT none 
!
INTEGER, INTENT(IN) :: POW_MAXPKT
REAL(KIND=PREC_DP)   , DIMENSION(0:POW_MAXPKT), INTENT(INOUT) :: dat       ! Data to be convoluted
REAL(KIND=PREC_DP)                            , INTENT(IN)    :: tthmin    ! 2Theta min
REAL(KIND=PREC_DP)                            , INTENT(IN)    :: tthmax    ! 2Theta max
REAL(KIND=PREC_DP)                            , INTENT(IN)    :: dtth      ! 2Theta step
REAL(KIND=PREC_DP)                            , INTENT(IN)    :: sigma2    ! Gaussian Sigma^2
REAL(KIND=PREC_DP)                            , INTENT(IN)    :: corrlin   ! 1/r deppendend width correction
REAL(KIND=PREC_DP)                            , INTENT(IN)    :: corrquad  ! 1/r^2 deppendend width correction
REAL(KIND=PREC_DP)              , INTENT(IN)    :: rcut      ! minimum  distance for clin/(r-rmin)
REAL(KIND=PREC_DP)                            , INTENT(IN)    :: pow_width ! Number of FWHM's to calculate
!
REAL(KIND=PREC_DP), PARAMETER :: eightln2  = 2.772588722239781237669D0 * 2.0D0
!
REAL(KIND=PREC_DP)            :: fwhm     ! Current FWHM at Theta
REAL(kind=PREC_DP), DIMENSION(0:POW_MAXPKT) :: dummy    ! temporary data (0:POW_MAXPKT) 
REAL(KIND=PREC_DP)            :: tth      ! Theta within convolution, main data set
REAL(KIND=PREC_DP)    :: sigmasq          ! actual scaled local sigma**2
REAL(KIND=PREC_DP)    :: sigmamin         ! minimum       local sigma**2
REAL(KIND=PREC_DP)    :: eta              ! actual eta at current 2Theta
REAL(KIND=PREC_DP)    :: dist_min  ! minimum  distance for clin/(r-rmin)
INTEGER :: imax, i, j, ii  ! Dummy loop indices
INTEGER :: i1, i2          ! Pseudo Voigt lookup indices
REAL(KIND=PREC_DP)    :: pseudo          ! scale factor for lookup table
INTEGER :: max_ps 
!                                                                       
!------ Now convolute                                                   
!                                                                       
imax = INT( (tthmax - tthmin) / dtth )
sigmasq = sigma2*SQRT(eightln2)
sigmamin = sigmasq * 0.020D0
dist_min = REAL(rcut, KIND=PREC_DP)
!
eta = 0.0D0     ! Gaussian function
dummy = 0.0D0   ! dummy(:)
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
!write(*,*) sigmasq, sigmamin
!write(*,*) corrlin, corrquad
!tth = 2.320
!fwhm = SQRT(MAX(sigmasq - corrlin/tth - corrquad/tth**2, sigmamin))
!write(*,*) ' FWHM ',tth, fwhm , sigmasq, sigmasq - corrlin/tth - corrquad/tth**2
!tth = 4.8000
!fwhm = SQRT(MAX(sigmasq - corrlin/tth - corrquad/tth**2, sigmamin))
!write(*,*) ' FWHM 2.500 ',tth, fwhm 
!tth = tthmax
!fwhm = SQRT(MAX(sigmasq - corrlin/tth - corrquad/tth**2, sigmamin))
!write(*,*) ' FWHM FINAL ',tth, fwhm 
main_pts: DO i = 0, imax 
   tth = tthmin + i * dtth 
   IF(tth<0.50) THEN
      dummy(i) = dat(i)
      cycle main_pts
   ENDIF
!
!  fwhm = SQRT(MAX(sigmasq - corrlin/tth -corrquad/tth**2, 0.00001))
   fwhm = SQRT(MAX(sigmasq - corrlin/(tth-dist_min) -corrquad/(tth-dist_min)**2, sigmamin))
!
   max_ps = INT((pow_width * fwhm) / dtth )
   pseudo =     dtth/fwhm*glp_npt  ! scale factor for look up table
   i1 = 0                               ! == 0 * dtth
   i2 = MIN(INT(2*i*pseudo), GLP_MAX)   ! == 2*i*dtth
   dummy (i) = dat (i) * ( glp_pseud_indx(i1, eta, fwhm)  &
                          -glp_pseud_indx(i2, eta, fwhm))                             
!
   ii = MAX (i - 1 - max_ps + 1, 1)
   first:DO j = ii, i - 1 
      i1 = MIN(INT((i - j) * pseudo), GLP_MAX)  ! == tth1 = (i - j) * dtth
      i2 = MIN(INT((i + j) * pseudo), GLP_MAX)  ! == tth2 = (i + j) * dtth
      dummy(i) = dummy(i) + dat(j) *( glp_pseud_indx(i1, eta, fwhm)  &
                                     -glp_pseud_indx(i2, eta, fwhm))                    
   ENDDO first
!
   ii = MIN(i + 1 + max_ps - 1, imax, POW_MAXPKT)
   secnd: DO j = i + 1, ii 
      i1 = MIN(INT((j - i) * pseudo), GLP_MAX)  ! == tth1 = (j - i) * dtth
      i2 = MIN(INT((j + i) * pseudo), GLP_MAX)  ! == tth1 = (j + i) * dtth
      dummy(i) = dummy(i) + dat(j) *( glp_pseud_indx(i1, eta, fwhm)  &
                                     -glp_pseud_indx(i2, eta, fwhm))
   ENDDO secnd
ENDDO main_pts
!                                                                       
DO i = 0, imax 
   dat (i) = dummy (i) * dtth   ! scale with stepwidth
ENDDO 
!                                                                       
END SUBROUTINE powder_conv_corrlin_old
!
!*****7*****************************************************************
!
subroutine powder_conv_tof(dat, xxmin,xxmax, xxdel, pow_width, MAXPKT,          &
   difa, difb, difc, zero, alpha0, alpha1, beta0, beta1, betaq, sigma02,        &
   sigma12, sigma22, sigmaq, gamma0, gamma1, gamma2, gsize, gstrain,            &
   pow_four_type, POW_COMPL, POW_NUFFT, POW_GRID, POW_DEBYE)
!-
!   Convolute with TOF profile function
!+
USE profile_tof_mod
use precision_mod                                                                       
USE trig_degree_mod
USE wink_mod
!
IMPLICIT none 
!
INTEGER, INTENT(IN) :: MAXPKT
REAL(KIND=PREC_DP)   , DIMENSION(0:MAXPKT), INTENT(INOUT) :: dat   ! Data to be convoluted
REAL(KIND=PREC_DP)                        , INTENT(IN)    :: xxmin  ! Q      min
REAL(KIND=PREC_DP)                        , INTENT(IN)    :: xxmax  ! Q      max
REAL(KIND=PREC_DP)                        , INTENT(IN)    :: xxdel  ! Q      step
REAL(KIND=PREC_DP)                        , INTENT(IN)    :: pow_width     ! Number of FWHM's to calculate
REAL(KIND=PREC_DP)                        , INTENT(IN)    :: difa          ! dval to time
REAL(KIND=PREC_DP)                        , INTENT(IN)    :: difb          ! dval to time
REAL(KIND=PREC_DP)                        , INTENT(IN)    :: difc          ! dval to time
REAL(KIND=PREC_DP)                        , INTENT(IN)    :: zero          ! dval to time
real(kind=PREC_DP)                        , intent(in)    :: alpha0     ! [         / mys]
real(kind=PREC_DP)                        , intent(in)    :: alpha1     ! [         / mys]
real(kind=PREC_DP)                        , intent(in)    :: beta0      ! [         / mys]
real(kind=PREC_DP)                        , intent(in)    :: beta1      ! [         / mys]
real(kind=PREC_DP)                        , intent(in)    :: betaq      ! [         / mys]
real(kind=PREC_DP)                        , intent(in)    :: sigma02   ! [             mys  ]
real(kind=PREC_DP)                        , intent(in)    :: sigma12   ! [             mys^2]
real(kind=PREC_DP)                        , intent(in)    :: sigma22   ! [             mys^2]
real(kind=PREC_DP)                        , intent(in)    :: sigmaq    ! [             mys^2]
real(kind=PREC_DP)                        , intent(in)    :: gamma0    ! [             mys  ]
real(kind=PREC_DP)                        , intent(in)    :: gamma1    ! [             mys  ]
real(kind=PREC_DP)                        , intent(in)    :: gamma2    ! [             mys  ]
real(kind=PREC_DP)                        , intent(in)    :: gsize     ! [             mys  ]
real(kind=PREC_DP)                        , intent(in)    :: gstrain   ! [             mys  ]
integer                                   , INTENT(IN)    :: pow_four_type ! Complete or Debye
integer                                   , INTENT(IN)    :: POW_COMPL     ! Complete or Debye
integer                                   , INTENT(IN)    :: POW_NUFFT     ! Complete or Debye
integer                                   , INTENT(IN)    :: POW_GRID      ! Complete or Debye
integer                                   , INTENT(IN)    :: POW_DEBYE     ! Debye
!
integer :: i, j
integer :: imax, jmax
real(kind=PREC_DP), dimension(:), allocatable :: conv
real(kind=PREC_DP)                            :: dval
!
real(kind=PREC_DP) :: alpha     ! [         / mys]
real(kind=PREC_DP) :: beta      ! [         / mys]
real(kind=PREC_DP) :: sg        ! [             mys  ] sqrt ( big_gamma/(8 ln(2))
real(kind=PREC_DP) :: sg2       ! [             mys^2]        big_gamma/(8 ln(2))
real(kind=PREC_DP) :: big_gamma    ! [             mys  ]
real(kind=PREC_DP) :: pre_g   ! Prefactor Gaussian   (1-eta)* NORM
real(kind=PREC_DP) :: pre_l   ! Prefactor Lorentzian  2*eta * NORM/PI
!
real(kind=PREC_DP) :: time_c    ! [             mys  ]
real(kind=PREC_DP) :: time_l    ! [             mys  ]
real(kind=PREC_DP) :: dt
!real(kind=PREC_DP) :: dv
real(kind=PREC_DP) :: ww
!integer :: bank
!
!bank = 3
!if(bank==3) then
!   alpha1 =  0.099198D0
!   beta0  =  3.340090D0
!   beta1  =  1.000000D-005
!   sigma02 =  0.000000D+000
!   sigma12 =   620.716000D0
!   sigma22 =   176.425000D0
!   gamma0 =   0.000000D+000
!   gamma1 =   22.268400D0
!   gamma2 =   0.000000D+000
!elseif(bank==4) then
!   alpha1 =   0.122511
!   beta0  =   0.058640
!   beta1  =   0.095311
!   sigma02 =   0.000000E+000
!   sigma12 =   43.864000
!   sigma22 =   122.734000
!   gamma0 =   1.812640
!   gamma1 =   21.436200
!   gamma2 =   4.326450
!elseif(bank==5) then
!   alpha1 =   0.094346
!   beta0  =   0.020882
!   beta1  =   0.034961
!   sigma02 =   2.971850
!   sigma12 =   10.000000
!   sigma22 =   12.938600
!   gamma0 =   0.000000E+000
!   gamma1 =   9.039320
!   gamma2 =   4.840550
!endif
!gsize   = 0.0
!gstrain = 0.0
!
allocate(conv(0:MAXPKT))
!
conv = 0.0D0
!
imax = int( ( xxmax -  xxmin) / xxdel )
jmax = int(0.5/xxdel)
!
if(pow_four_type==POW_COMPL .or. pow_four_type==POW_NUFFT .or. POW_FOUR_TYPE==POW_GRID) then        ! Complete calculation
   do i=0, imax
      if(dat(i) > 0.0D0) then            ! Non-zero data do convolution
         dval   = zpi/(xxmin + i*xxdel) ! D-value from Q
         time_c = difc*dval + difa*dval**2 + zero
         call tof3_para( dval, alpha0, alpha1, beta0, beta1, betaq,             &
                sigma02, sigma12, sigma22, sigmaq,                              &
                gamma0, gamma1, gamma2, gsize, gstrain, difc,                   &
                alpha, beta, sg, sg2, big_gamma, pre_g, pre_l)
         do j=max(0,i-jmax), min(i+jmax,imax)
            dval   = zpi/(xxmin + j*xxdel) ! D-value from Q
            time_l = difc*dval + difa*dval**2 + zero
            dt = time_l - time_c
            ww = profile_tof3(dt, dval, alpha, beta, sg, &
                 sg2, big_gamma, pre_g, pre_l)
            conv(j) = conv(j) + ww*dat(i)
         enddo
      endif
   enddo
else                                     ! Debye-Scattering-Equation
endif
!
dat(0:imax) = conv(0:imax)
!
deallocate(conv)
!
end subroutine powder_conv_tof
!
!*****7*****************************************************************
!
      REAL(KIND=PREC_DP) function pseudovoigt (dtth, eta, fwhm) 
!-                                                                      
!     calculates the value of a pseudo-voigt function at dtth off the   
!     central position                                                  
!
use precision_mod                                                                       
!
      IMPLICIT none 
      REAL(KIND=PREC_DP), INTENT(IN) :: dtth 
      REAL(KIND=PREC_DP), INTENT(IN) :: eta 
      REAL(KIND=PREC_DP), INTENT(IN) :: fwhm 
!                                                                       
      REAL(KIND=PREC_DP), PARAMETER :: four_ln2  = 2.772588722
      REAL(KIND=PREC_DP), PARAMETER :: sq4ln2_pi = 0.939437279
      REAL(KIND=PREC_DP), PARAMETER :: two_pi    = 0.636619772
      REAL(KIND=PREC_DP)            :: pref_g 
      REAL(KIND=PREC_DP)            :: pref_l 
!                                                                       
!                                                                       
      pref_g = sq4ln2_pi / fwhm 
      pref_l = two_pi * fwhm 
!                                                                       
      pseudovoigt = eta * pref_l / (fwhm**2 + 4. * dtth**2) + (1. - eta)&
      * pref_g * exp ( - four_ln2 * (dtth / fwhm) **2)                  
!                                                                       
      END FUNCTION pseudovoigt                      
!*****7*****************************************************************
!
REAL(KIND=PREC_DP) function profile_asymmetry (tth, dtth, fwhm, p1, p2)!, p3, p4) 
!-                                                                      
!     calculates the asymmetry parameter for the profile function       
!                                                                       
USE trig_degree_mod
USE precision_mod
!
IMPLICIT none 
!                                                                       
REAL(KIND=PREC_DP), INTENT(IN) :: tth 
REAL(KIND=PREC_DP), INTENT(IN) :: dtth 
REAL(KIND=PREC_DP), INTENT(IN) :: fwhm 
REAL(KIND=PREC_DP), INTENT(IN) :: p1, p2!, p3, p4 
!                                                                       
REAL(KIND=PREC_DP) :: zz 
REAL(KIND=PREC_DP) :: fa, fb 
!                                                                       
!                                                                       
zz = dtth / fwhm 
!                                                                       
fa = 2. * zz * exp ( - zz**2) 
fb = 2. * (2. * zz**2 - 3.) * fa 
!                                                                       
profile_asymmetry = 1.0 + (p1 * fa + p2 * fb)! / tand(0.5 * tth)  &
!profile_asymmetry = 1.0 + (p1 * fa + p2 * fb) / tand(0.5 * tth)  &
!                        + (p3 * fa + p4 * fb) / tand(tth)                                
!                                                                       
END FUNCTION profile_asymmetry                
!
!*****7*****************************************************************
!
SUBROUTINE powder_f2aver ( num1 )
!
!     This subroutine calculates the average atomic form factor
!     <f^2> and <f>^2
!
USE crystal_mod 
USE chem_mod
USE diffuse_mod 
USE powder_mod 
USE powder_tables_mod 
USE wink_mod
!
IMPLICIT NONE
!
INTEGER,                INTENT(IN) :: num1
!
INTEGER, DIMENSION(:), ALLOCATABLE :: natom
!
INTEGER :: iscat
INTEGER :: i
REAL( KIND(0.0D0))             :: signum
!!!
pow_f2aver(:)    = 0.0D0
pow_faver2(:)    = 0.0D0
pow_f2    (:,:)  = 0.0D0
pow_u2aver       = 0.0
pow_nreal        = 0
!
!     Prepare and calculate average atom numbers
!
ALLOCATE(natom(0:cr_nscat))
natom = 0
DO i=1,cr_natoms
   natom(cr_iscat(1,i)) = natom(cr_iscat(1,i)) + 1
ENDDO
pow_nreal = SUM(NINT(natom(1:cr_nscat)*cr_occ(1:cr_nscat)))  ! Add real atom numbers 
IF(chem_quick) THEN
   pow_ncreal = pow_nreal/(cr_icc(1)*cr_icc(2)*cr_icc(3))
ELSE
   pow_ncreal = cr_ncatoms
ENDIF
!write(*,*) ' IN POW_f2aver nreal, ', pow_nreal, pow_ncreal, cr_nscat
!
DO iscat = 1, cr_nscat
   signum = 1.0D0
   IF(REAL(cfact_pure(1, iscat))< 0.0D0) signum = -1.0D0
!write(*,*) ' ', SQRT(DBLE (       cfact_pure(powder_istl(1), iscat)  * &
!                       conjg (cfact_pure(powder_istl(1), iscat)))), &
!                natom(iscat), cr_occ(iscat), pow_nreal,  cr_dw(iscat)
   DO i = 0, num1
      pow_f2aver (i) = pow_f2aver (i)  + &
                 DBLE (       cfact_pure(powder_istl(i), iscat)  * &
                       conjg (cfact_pure(powder_istl(i), iscat)))  &
               * natom (iscat)/pow_nreal*cr_occ(iscat)
      pow_faver2 (i) = pow_faver2 (i) +  &
            SQRT(DBLE (       cfact_pure(powder_istl(i), iscat)  * &
                       conjg (cfact_pure(powder_istl(i), iscat)))) &
               * natom (iscat)/pow_nreal *cr_occ(iscat)            &
               * signum
      pow_f2(i,iscat)= pow_f2(i,iscat) + &
                 DBLE (       cfact_pure(powder_istl(i), iscat)  * &
                       conjg (cfact_pure(powder_istl(i), iscat)))  &
               * natom (iscat)/pow_nreal*cr_occ(iscat)
   ENDDO
   pow_u2aver = pow_u2aver + (cr_dw(iscat)) * natom (iscat)/pow_nreal*cr_occ(iscat)
ENDDO
pow_faver2(:) = pow_faver2(:)**2
pow_u2aver    = pow_u2aver /8./REAL(pi)**2
!write(*,*) ' f2aver ', pow_f2aver(0), pow_f2aver(1)
!write(*,*) ' faver2 ', pow_faver2(0), pow_faver2(1)
!write(*,*) ' U2aver ', pow_u2aver, pow_u2aver*8*PI**2
!write(*,*) ' BVAL   ', cr_dw(1:cr_nscat)
!write(*,*) ' Natom  ', natom(1:cr_nscat)
DEALLOCATE(natom)
!
!
END SUBROUTINE powder_f2aver
!
!*******************************************************************************
!
SUBROUTINE pow_k12(npkt, POW_MAXPKT, pow_ka21, pow_ka21_u, xmin, xdel, xpl, ypl)
!-
! Add Kalpha1/2
! at this point xpl is on Q-scale
! The powder pattern is scaled by q/len_ratio. The original data are added 
! to the two data points j, j` with weights according to the location of q2
! between these two indices.
!+
USE diffuse_mod
USE element_data_mod
!
USE precision_mod
use param_mod
!
IMPLICIT NONE
!
INTEGER                                    , INTENT(IN)    :: npkt
INTEGER                                    , INTENT(IN)    :: POW_MAXPKT
REAL(KIND=PREC_DP)                         , INTENT(IN)    :: pow_ka21
LOGICAL                                    , INTENT(IN)    :: pow_ka21_u
REAL(KIND=PREC_DP)                         , INTENT(IN)    :: xmin   ! Q min
REAL(KIND=PREC_DP)                         , INTENT(IN)    :: xdel   ! Q Step 
REAL(KIND=PREC_DP), DIMENSION(0:POW_MAXPKT), INTENT(IN)    :: xpl
REAL(KIND=PREC_DP), DIMENSION(0:POW_MAXPKT), INTENT(INOUT) :: ypl
!
CHARACTER(LEN=4) :: alpha            ! Kalpa1 name
INTEGER :: i,j,k,l
REAL(KIND=PREC_DP)    :: int_ratio   ! Ka2/Ka1 intensity ratio
REAL(KIND=PREC_DP)    :: len_ratio   ! Ka2/Ka1 intensity ratio
REAL(KIND=PREC_DP)    :: int_ratio_b ! Ka2/K   intensity ratio
REAL(KIND=PREC_DP)    :: len_ratio_b ! Ka2/Kb  intensity ratio
REAL(KIND=PREC_DP)    :: part1       ! Ka2/Ka1 intensity ratio
REAL(KIND=PREC_DP)    :: part2       ! Ka2/Ka1 intensity ratio
REAL(KIND=PREC_DP)    :: q1, q2, q3  ! Ka2/Ka1 intensity ratio
real(kind=PREC_DP), dimension(:), allocatable :: ytemp  ! Temporary data
real(kind=PREC_DP) :: alpha_len
real(kind=PREC_DP) :: white_abs
real(kind=PREC_DP) :: white_damp
real(kind=PREC_DP) :: white_inte
!
int_ratio   = 0.0
int_ratio_b = 0.0
len_ratio   = 1.0
!
allocate(ytemp(0:POW_MAXPKT))
ytemp = ypl   ! Copy Ka1 into accumulated pattern
!
cond_lambda: IF(lambda(3:4)=='12' .OR. lambda=='W12 ') THEN    ! Only for explicit a12
   l = get_wave_number(lambda)
!   call get_sym_length(l-4, alpha, alpha_len)
!   write(*,*) ' ALPHA     ', alpha_len
   IF(l==0) exit cond_lambda   ! Not a listed wave length
!write(*,*) ' WAVE LENGTH NUMBER Is ', l
!
   IF(pow_ka21_u) THEN
      int_ratio = pow_ka21
   ELSE
      int_ratio = get_ka21_inte(l)
   ENDIF
   len_ratio = get_ka12_len(l)
!   len_ratio = 1./rpara(499)
!
   DO i=npkt,1,-1
      q1 = (xmin + xdel * real(i,kind=PREC_DP))            ! Original Q
      q2 = (xmin + xdel * real(i,kind=PREC_DP))/len_ratio  ! Q at Kalpha2
      j = int((q2-xmin)/xdel)                              ! Lower pixel
      part1 = 1.0D0 - ((q2-xmin)/xdel - j)                 ! Weight for lower pixel
      part2 = 1.0D0 - part1                                ! Weight for upper pixel
      if(j < npkt) THEN
         ytemp(j) = ytemp(j) + int_ratio*ypl(i)*part1          ! Add to lower pixel
         j = j + 1
         if(j < npkt) THEN
            ytemp(j) = ytemp(j) + int_ratio*ypl(i)*part2       ! Add to upper pixel
         ENDIF
      ENDIF
   ENDDO
elseif(lambda(len_trim(lambda)-1:len_trim(lambda))=='NI') THEN    cond_lambda ! Only for explicit Nickel filter
   l = get_wave_number(lambda)
   IF(l==0) exit cond_lambda   ! Not a listed wave length
!
   IF(pow_ka21_u) THEN
      int_ratio = pow_ka21
   ELSE
      int_ratio = get_ka21_inte(l-1)    ! Ka12 is stores onw down from ..NI
   ENDIF
   len_ratio = get_ka12_len(l-1)
!write(*,*) ' WAVE LENGTH NUMBER Is ', l, len_ratio, int_ratio
!   len_ratio = 1./rpara(499)
!
   DO i=npkt,1,-1
      q1 = (xmin + xdel * real(i,kind=PREC_DP))            ! Original Q
      q2 = (xmin + xdel * real(i,kind=PREC_DP))/len_ratio  ! Q at Kalpha2
      j = int((q2-xmin)/xdel)                              ! Lower pixel
      part1 = 1.0D0 - ((q2-xmin)/xdel - j)                 ! Weight for lower pixel
      part2 = 1.0D0 - part1                                ! Weight for upper pixel
      if(j < npkt) THEN
         ytemp(j) = ytemp(j) + int_ratio*ypl(i)*part1          ! Add to lower pixel
         j = j + 1
         if(j < npkt) THEN
            ytemp(j) = ytemp(j) + int_ratio*ypl(i)*part2       ! Add to upper pixel
         ENDIF
      ENDIF
   ENDDO
!  Now do Kbeta
   call get_sym_length(l-5, alpha, alpha_len)
   int_ratio_b = get_kabe_inte(l)
   len_ratio_b = get_kabe_len(l)
!
!  white_abs    = 1.490
!  white_damp   = 1.560
!  white_inte   = 0.00016_PREC_DP
!
!  white_abs    = rpara(443)
!  white_damp   = rpara(444)
!  white_inte   = rpara(445)
!
   white_abs    = get_white_abs(l)
   white_damp   = get_white_damp(l)
   white_inte   = get_white_inte(l)
!write(*,*) ' ALPHA     ', alpha_len
!write(*,*) ' WHITE abs ', white_abs
!write(*,*) ' WHITE dam ', white_damp
!write(*,*) ' WHITE int ', white_inte
!
   DO i=npkt,1,-1
      q1 = (xmin + xdel * real(i,kind=PREC_DP))              ! Original Q
      q2 = (xmin + xdel * real(i,kind=PREC_DP))/len_ratio_b  ! Q at Kbeta  
      j = int((q2-xmin)/xdel)                                ! Lower pixel
      part1 = 1.0D0 - ((q2-xmin)/xdel - j)                   ! Weight for lower pixel
      part2 = 1.0D0 - part1                                  ! Weight for upper pixel
      if(j > 0   ) THEN
         ytemp(j) = ytemp(j) + int_ratio_b*ypl(i)*part1          ! Add to lower pixel
         j = j + 1
         if(j > 0   ) THEN
            ytemp(j) = ytemp(j) + int_ratio_b*ypl(i)*part2       ! Add to upper pixel
         ENDIF
      ENDIF
   ENDDO
!
! now do the white line
!
   do i=1, npkt
      if(ypl(i) > 10.0) then
!     i = (3.088 - xmin)/xdel
      q1 = (xmin + xdel * real(i,kind=PREC_DP))              ! Original Q
      q2 = (xmin + xdel * real(i,kind=PREC_DP))*(white_abs /alpha_len) ! Q at Ni Absorption edge
      q3 = (xmin + xdel * real(i,kind=PREC_DP))*(white_damp/alpha_len) -q2 ! Sigma
      j = int((q2-xmin)/xdel)                                ! Lower pixel
      k = int((q2+4.0_PREC_DP*q3-xmin)/xdel)                 ! Upper limit
      do l = j, k
         if(l>0 .and. l<npkt) then
            q1 = (xmin + xdel * real(l,kind=PREC_DP))              ! Original Q
            ytemp(l) = ytemp(l) + ypl(i)*white_inte    *exp(-0.5_PREC_DP*(q1-q2)**2/q3**2)
         endif
      enddo
     endif
  enddo
 
ENDIF cond_lambda
!
ypl = ytemp / (1.+int_ratio + int_ratio_b)

deallocate(ytemp)
!
END SUBROUTINE pow_k12
!
!*******************************************************************************
!
END MODULE powder_write_mod
