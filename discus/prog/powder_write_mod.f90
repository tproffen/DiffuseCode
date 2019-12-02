MODULE powder_write_mod
!
USE errlist_mod 
!
IMPLICIT NONE
!
PUBLIC
!
CONTAINS
!
SUBROUTINE powder_out (value)
!-                                                                      
!     Write the powder pattern                                          
!+                                                                      
USE crystal_mod
USE discus_config_mod 
USE debye_mod 
USE diffuse_mod 
USE output_mod 
USE powder_mod 
USE powder_tables_mod
USE pdf_mod
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
!                                                                       
INTEGER   :: ii, j , iii
INTEGER   :: all_status  ! Allocation status
INTEGER   :: npkt        ! number of points in powder pattern
INTEGER   :: npkt_equi   ! number of points in equidistant powder pattern
INTEGER   :: npkt_wrt    ! number of points in powder pattern ready to write
INTEGER   :: npkt_fft    ! number of points in powder pattern for Fast Fourier
LOGICAL   :: lread 
!!      LOGICAL   :: lconv = .FALSE.  ! Did convolution with profile
REAL, DIMENSION(:), ALLOCATABLE :: pow_tmp  ! Local temporary copy of intensities
REAL, DIMENSION(:), ALLOCATABLE :: xpl  ! x-values of calculated powder pattern
REAL, DIMENSION(:), ALLOCATABLE :: ypl  ! y-values of calculated powder pattern
REAL, DIMENSION(:), ALLOCATABLE :: lpv  ! Values of LP correction versus Q/Theta 
REAL, DIMENSION(:), ALLOCATABLE :: y2a  ! y-values of splined    powder pattern
REAL, DIMENSION(:), ALLOCATABLE :: xwrt ! x-values of powder pattern ready for output
REAL, DIMENSION(:), ALLOCATABLE :: ywrt ! y-values of powder pattern ready for output
REAL :: ttheta, lp=1.0
REAL :: lpscale      ! Scale factor introduced by LP correction
REAL ss, st 
REAL :: q=0.0, stl=0.0, dstar=0.0
REAL      :: normalizer
REAL xmin, xmax, xdel , xpos
REAL (PREC_DP)     :: xstart  ! qmin  for sin Theta / lambda calculation
REAL (PREC_DP)     :: xdelta  ! qstep for sin Theta / lambda calculation
REAL      :: xequ    ! x-position of equdistant curve
REAL      :: yequ    ! y-value    of equdistant curve
REAL      :: tthmin  ! minimum for equdistant curve
REAL      :: tthmax  ! minimum for equdistant curve
REAL      ::   qmin  ! minimum for equdistant curve
REAL      ::   qmax  ! maximum for equdistant curve
REAL      :: arg
REAL(KIND=PREC_DP) :: u2aver_scale = 2.00   ! Scale to multiply <u^2> if conversion
!                     ! with corrlin_corrquad is needed. This increases the calculated
!                     ! intensity actually by more than the damping by <u^2>, in oder
!                     ! to sharpen the distances sufficiently for later broadening
!
REAL(KIND=PREC_DP) :: rmin, rmax, rstep
INTEGER            :: npkt_pdf
REAL, DIMENSION(:), ALLOCATABLE :: xfour
REAL, DIMENSION(:), ALLOCATABLE :: yfour
!
REAL :: sigma   ! sigma for PDF COrrlin correction
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
ALLOCATE(pow_tmp(0:POW_MAXPKT),stat = all_status)  ! Allocate array for powder pattern copy
ALLOCATE(xpl(0:POW_MAXPKT),stat = all_status)  ! Allocate array for calculated powder pattern
ALLOCATE(ypl(0:POW_MAXPKT),stat = all_status)  ! Allocate array for calculated powder pattern
ALLOCATE(lpv(0:POW_MAXPKT),stat = all_status)  ! Allocate array for LP correction
pow_tmp = 0.0
xpl     = 0.0
ypl     = 0.0
xdel    = 0.0
!                                                                       
IF (pow_four_type.eq.POW_COMPL) THEN 
!  IF (pow_axis.eq.POW_AXIS_Q) THEN 
      xmin = pow_qmin_u 
      xmax = pow_qmax_u 
!     IF(value == val_pdf) THEN
!        IF(pdf_clin_a>0.0 .OR. pdf_cquad_a>0.0) THEN
!           xmax = MAX(30.0D0, pow_qmax_u*1.2D0)
!        ENDIF
!     ENDIF
      xdel = pow_deltaq_u
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
!     DEALLOCATE(pow_tmp,stat = all_status)  ! DeAllocate array for powder pattern copy
!     DEALLOCATE(xpl    ,stat = all_status)  ! DeAllocate array for calculated powder pattern
!     DEALLOCATE(ypl    ,stat = all_status)  ! DeAllocate array for calculated powder pattern
!     DEALLOCATE(lpv    ,stat = all_status)  ! DeAllocate array for LP correction
!     RETURN 
!  ENDIF 
   npkt = MIN(NINT((xmax+xdel-xmin)/xdel) + 2, POW_MAXPKT)
ELSEIF (pow_four_type.eq.POW_DEBYE) THEN
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
!     DEALLOCATE(pow_tmp,stat = all_status)  ! DeAllocate array for powder pattern copy
!     DEALLOCATE(xpl    ,stat = all_status)  ! DeAllocate array for calculated powder pattern
!     DEALLOCATE(ypl    ,stat = all_status)  ! DeAllocate array for calculated powder pattern
!     DEALLOCATE(lpv    ,stat = all_status)  ! DeAllocate array for LP correction
!     RETURN 
!  ENDIF 
   npkt = MIN(num(1), POW_MAXPKT)
ENDIF 
!
!     Prepare average form factors for S(Q) or F(Q), Normalized Intensity, PDF, or faver2, f2aver
!
IF(value == val_sq     .OR. value == val_fq     .OR. &
   value == val_iq     .OR. value == val_inten  .OR. &
                            value == val_pdf    .OR. &
   value == val_f2aver .OR. value == val_faver2      ) THEN
!  IF (pow_axis.eq.POW_AXIS_Q .OR. (pow_axis.eq.POW_AXIS_TTH .AND.value == val_inten)) THEN 
      IF(.NOT.(pow_four_mode == POW_STACK)) THEN  ! Stack did its own faver2
         IF (pow_four_type.eq.POW_COMPL) THEN     ! Need to initialize pow_istl
            xstart = pow_qmin  /zpi
            xdelta = pow_deltaq/zpi
            CALL powder_stltab(npkt,xstart,xdelta) ! Really only needed for <f^2> and <f>^2 for F(Q) and S(Q)
         ENDIF
         CALL powder_f2aver (npkt   )             ! Calculate average form factors <f>2 and <f^2>
      ENDIF
!  ELSE                                           ! F(Q) works for Q-axis only
!     ier_msg (1) = 'Use command ==> form, powder,q'
!     ier_msg (2) = 'within the output menu to define the axis' 
!     ier_num = -125
!     ier_typ = ER_APPL
!     DEALLOCATE(pow_tmp,stat = all_status)  ! DeAllocate array for powder pattern copy
!     DEALLOCATE(xpl    ,stat = all_status)  ! DeAllocate array for calculated powder pattern
!     DEALLOCATE(ypl    ,stat = all_status)  ! DeAllocate array for calculated powder pattern
!     DEALLOCATE(lpv    ,stat = all_status)  ! DeAllocate array for LP correction
!     RETURN
! ENDIF
ENDIF
!
lread = .false. 
IF (ier_num /= 0) THEN 
   DEALLOCATE(pow_tmp,stat = all_status)  ! DeAllocate array for powder pattern copy
   DEALLOCATE(xpl    ,stat = all_status)  ! DeAllocate array for calculated powder pattern
   DEALLOCATE(ypl    ,stat = all_status)  ! DeAllocate array for calculated powder pattern
   DEALLOCATE(lpv    ,stat = all_status)  ! DeAllocate array for LP correction
   RETURN
ENDIF
!
IF(value == val_f2aver) THEN          ! Output is f^2 aver
   DO j = 0, npkt
      pow_tmp (j) = REAL(pow_f2aver(j))
   ENDDO
ELSEIF(value == val_faver2) THEN     ! Output is faver^2 
   DO j = 0, npkt
      pow_tmp (j) = REAL(pow_faver2(j))
   ENDDO
ELSE                         ! All other output
   DO j = 0, npkt
      pow_tmp (j) = pow_conv(j)   ! copy from convoluted pattern
   ENDDO
ENDIF           ! Output is if_block if(value==val_f2aver)
!CORR open(77,file='POWDER/INITIAL.inte',status='unknown')
!CORR DO ii=1,npkt
!CORR   write(77,'(2(2x,G17.7E3))') xmin+(ii-1)*xdel, pow_tmp(ii)
!CORR enddo
!CORR close(77)
!
!------ copy the powder pattern into output array, if necessary this will be put on
!       equidistant scale
!                                                                       
IF (pow_four_type.eq.POW_COMPL) THEN
   pow_tmp = pow_tmp/(REAL(pow_nreal, KIND=PREC_DP))**2 * &
                      REAL(cr_ncatoms, KIND=PREC_DP)/REAL(cr_v, KIND=PREC_DP) * &
             20.00D0  ! Needs to be verified where this 20 comes from
ENDIF
!
lpv    = 0.0
lpv(0) = 1.0
ypl(0) = pow_tmp(0)
!
copy: IF (pow_four_type.eq.POW_COMPL) THEN
   DO ii = 1, npkt
      xpos = ii * xdel + xmin 
!     IF (pow_axis.eq.POW_AXIS_Q) THEN 
         q      = xpos
         dstar  = q / REAL(zpi)
         stl    = q / REAL(zpi) / 2.
         ttheta = 2.*asind ( REAL(q / 2. /REAL(zpi) *rlambda ))
!     ELSEIF (pow_axis.eq.POW_AXIS_TTH) THEN 
!        ttheta = xpos
!        stl    =                  sind (ttheta * 0.5) / rlambda 
!        dstar  = 2. *             sind (ttheta * 0.5) / rlambda 
!        q      = 2. * REAL(zpi) * sind (ttheta * 0.5) / rlambda 
!     ENDIF 
      IF(value == val_sq .OR. value == val_fq .OR. value == val_pdf)  THEN
! WRONG!         lp = lorentz(ttheta,1)
      lp = 1./q**2
      ELSEIF(value == val_f2aver .OR. value == val_faver2) THEN  ! f2aver or faver2
         lp = 1
      ELSE
!WRONG!         lp     = lorentz (ttheta,0) * polarisation (ttheta) 
      lp = 1./q**2 * polarisation(ttheta)
      ENDIF 
      lpv(ii) = polarisation (ttheta)
      IF (cpow_form.eq.'tth') THEN 
         xpl(ii) = ttheta
      ELSEIF (cpow_form.eq.'stl') THEN 
         xpl(ii) = stl
      ELSEIF (cpow_form.eq.'q  ') THEN 
         xpl(ii) = q
      ELSEIF (cpow_form.eq.'dst') THEN 
         xpl(ii) = dstar
      ELSEIF (cpow_form.eq.'lop') THEN 
         xpl(ii) = ttheta
      ENDIF 
      ypl(ii) = pow_tmp(ii) * lp
   ENDDO 
ELSEIF (pow_four_type.eq.POW_DEBYE) THEN  copy
!  IF (pow_axis.eq.POW_AXIS_DSTAR) THEN 
!     CONTINUE 
!  ELSEIF (pow_axis.eq.POW_AXIS_Q) THEN 
      xm(1)  = pow_qmin / REAL(zpi) 
      ss     = pow_qmax / REAL(zpi) 
      st     = (pow_qmax - pow_deltaq) / REAL(zpi) 
      uin(1) = pow_deltaq / REAL(zpi) 
!  ELSEIF (pow_axis.eq.POW_AXIS_TTH) THEN 
!     xm(1)  = 2 * sind (0.5 * pow_tthmin) / rlambda 
!     ss     = 2 * sind (0.5 * pow_tthmax) / rlambda 
!     st     = 2 * sind (0.5 * (pow_tthmax - pow_deltatth) ) / rlambda
!     uin(1) = (ss - st) / 2. 
!  ENDIF 
   DO ii = 1, npkt    
      dstar = REAL(xm (1) + (ii - 1) * uin (1) ) 
      stl = REAL(0.5D0 * (xm (1) + (ii - 1) * uin (1) ) )
      q = REAL(zpi) * REAL(xm (1) + (ii - 1) * uin (1) ) 
      ttheta = 2. * asind (dstar * rlambda / 2.) 
!
      IF(value == val_sq .or. value == val_fq .OR. value == val_pdf) THEN
         lp = 1                     ! For S(Q) and F(Q) nor Polarisation corr.
      ELSEIF(value == val_f2aver .or. value == val_faver2) THEN  ! f2aver or faver2
         lp = 1
      ELSE
         lp = polarisation (ttheta) 
      ENDIF
      lpv(ii ) = lp
      IF (cpow_form.eq.'tth') THEN 
         xpl(ii) = ttheta
      ELSEIF (cpow_form.eq.'stl') THEN 
         xpl(ii) = stl
      ELSEIF (cpow_form.eq.'q  ') THEN 
         xpl(ii) = q
      ELSEIF (cpow_form.eq.'dst') THEN 
         xpl(ii) = dstar
      ENDIF 
      ypl(ii) = pow_tmp(ii) * lp
   ENDDO 
ENDIF  copy
!
IF(value == val_pdf) THEN   ! Divide by average Debye-Waller term
   IF(pdf_clin_a/=0.0 .OR. pdf_cquad_a/=0.0) THEN
!open(77,file='POWDER/normal.inte',status='unknown')
!DO ii=1,npkt
!write(77,'(2(2x,G17.7E3))') xmin+(ii-1)*xdel, ypl(ii)
!enddo
!close(77)
!  IF(pow_axis.eq.POW_AXIS_Q) THEN
      DO ii=1, npkt
         ypl(ii) = ypl(ii)/EXP(-0.5*(pow_u2aver*u2aver_scale)*xpl(ii)**2)
      ENDDO
!  ENDIF
ENDIF
!open(77,file='POWDER/divided.inte',status='unknown')
!DO ii=1,npkt
!write(77,'(2(2x,G17.7E3))') xmin+(ii-1)*xdel, ypl(ii)
!enddo
!close(77)
ENDIF        ! Output is if_block if(value==val_f2aver)
lpscale = 1.0
!
!     Prepare S(Q) or F(Q)
!
normalizer = 1.0D0
prsq: IF(value == val_sq .or. value == val_fq   .OR. value == val_inten  .OR. & 
         value == val_iq .OR.                        value == val_pdf         ) THEN
!  axq: IF (pow_axis.eq.POW_AXIS_Q .OR. (pow_axis.eq.POW_AXIS_TTH .AND.value == val_inten)) THEN 
      IF (pow_four_type.eq.POW_COMPL) THEN
!        pow_tmp_sum = 0.0                           ! Determine normalizer, such that 
!        pow_uuu_sum = 0.0                           ! the average F(q) is 0.0
!           jstart = MAX(1,2-int(xmin/xdel))            ! Exclude q = 0
!        jstart = MAX(1,int((1-xmin)/xdel)+1)        ! Exclude q < 1.0
!        DO j = jstart, npkt
!           q = ((j-1)*xdel + xmin)
!           pow_tmp_sum = pow_tmp_sum + ypl(j)/REAL(pow_faver2(j))* q
!           pow_uuu_sum = pow_uuu_sum       + exp(-q**2*pow_u2aver)*q * &
!                         pow_f2aver(j)/pow_faver2(j)
!        ENDDO
!        normalizer = pow_tmp_sum/pow_uuu_sum
         normalizer = 1.0D0
      ELSE
         normalizer = REAL(pow_nreal)
      ENDIF
!
      valq: IF(value == val_sq) THEN                   ! Calc S(Q)
         IF(deb_conv .OR. .NOT.ldbw) THEN              ! DEBYE was done with convolution of ADP
            DO j = 1, npkt   
               q = ((j-1)*xdel + xmin)
               ypl(j) =  (ypl(j)/REAL(pow_faver2(j))/normalizer    &
                         + 1.0 - & !exp(-q**2*pow_u2aver*0.25)      * &
                          pow_f2aver(j)/pow_faver2(j))
            ENDDO
         ELSE
            DO j = 1, npkt   
               q = ((j-1)*xdel + xmin)
               ypl(j) =  (ypl(j)/REAL(pow_faver2(j))/normalizer   &
                         + 1.0 - exp(-q**2*pow_u2aver)   * &
                          pow_f2aver(j)/pow_faver2(j))
            ENDDO
         ENDIF
      ELSEIF(value == val_fq) THEN  valq                       ! Calc F(Q)IF
         IF(deb_conv .OR. .NOT.ldbw) THEN                      ! DEBYE was done with convolution of ADP
            DO j = 1, npkt   
               q = ((j-1)*xdel + xmin)
               ypl(j) =  (ypl(j)/REAL(pow_faver2(j))/normalizer   &
                          -pow_f2aver(j)/pow_faver2(j)) * q
            ENDDO
         ELSE
            DO j = 1, npkt   
               q = ((j-1)*xdel + xmin)
               ypl(j) =  (ypl(j)/REAL(pow_faver2(j))/normalizer   &
                               - exp(-q**2*pow_u2aver)   *  &
                          pow_f2aver(j)/pow_faver2(j)) * q
            ENDDO
         ENDIF
      ELSEIF(value == val_pdf) THEN  valq                      ! Calc PDF IF
         IF(deb_conv .OR. .NOT.ldbw) THEN                      ! DEBYE was done with convolution of ADP
            DO j = 1, npkt   
               q = ((j-1)*xdel + xmin)
               ypl(j) =  (ypl(j)/REAL(pow_faver2(j))/normalizer   &
                          -pow_f2aver(j)/pow_faver2(j)) * q
            ENDDO
         ELSE
            IF(pdf_clin_a/=0.0 .OR. pdf_cquad_a/=0.0) THEN
               DO j = 1, npkt   
                  q = ((j-1)*xdel + xmin)
                  ypl(j) =  (ypl(j)/REAL(pow_faver2(j))/normalizer   &
                             - exp(-0.00*q**2*(pow_u2aver/u2aver_scale))   *  &
                             pow_f2aver(j)/pow_faver2(j)) * q
!                            - 1.0 *                                 &
               ENDDO
            ELSE
               DO j = 1, npkt   
                  q = ((j-1)*xdel + xmin)
                  ypl(j) =  (ypl(j)/REAL(pow_faver2(j))/normalizer   &
                             - exp(-q**2*pow_u2aver)   *  &
                             pow_f2aver(j)/pow_faver2(j)) * q
               ENDDO
            ENDIF
         ENDIF
      ELSEIF(value == val_iq) THEN   valq              ! Calc F(Q)IF
         lpscale = 1.0
         IF(deb_conv .OR. .NOT.ldbw) THEN              ! DEBYE was done with convolution of ADP
            DO j = 1, npkt   
               q = ((j-1)*xdel + xmin)
               ypl(j) =  (ypl(j)                    /normalizer    &
                         +         lpv(j) * (                      &
                         + pow_faver2(j) - & !exp(-q**2*pow_u2aver*0.25)      * &
                          pow_f2aver(j)              ) )
            ENDDO
         ELSE
            DO j = 1, npkt   
               q = ((j-1)*xdel + xmin)
               ypl(j) =  (ypl(j)                    /normalizer   &
                         +         lpv(j) * (                      &
                         + pow_faver2(j) - exp(-q**2*pow_u2aver)   * &
                          pow_f2aver(j)              ))
            ENDDO
         ENDIF
      ELSEIF(value == val_inten) THEN   valq              ! Calc Intensity
         IF(deb_conv .OR. .NOT.ldbw) THEN              ! DEBYE was done with convolution of ADP
            DO j = 1, npkt   
               q = ((j-1)*xdel + xmin)
               ypl(j) =  (ypl(j)/normalizer    &
                         ) 
!              ypl(j) =  (ypl(j)                    /normalizer    &
!                        +         lpv(j) * (                      &
!                        + pow_faver2(j) - & !exp(-q**2*pow_u2aver*0.25)      * &
!                         pow_f2aver(j)              ) ) * REAL(pow_nreal)
            ENDDO
         ELSE
            DO j = 1, npkt   
               q = ((j-1)*xdel + xmin)
               ypl(j) =   ypl(j)/normalizer +                 &
                         (+ 1.0 - exp(-q**2*pow_u2aver))*pow_faver2(j) 
                         !* REAL(pow_nreal)
!              ypl(j) =  (ypl(j)                    /normalizer   &
!                        +         lpv(j) * (                      &
!                        + pow_faver2(j) - exp(-q**2*pow_u2aver)   * &
!                         pow_f2aver(j)              )) * REAL(pow_nreal)
            ENDDO
         ENDIF
      ENDIF valq
!  ELSE axq                                       ! F(Q) works for Q-axis only
!           Should never occur, as covered by "prepare S(Q)" section
!     ier_msg (1) = 'Use command ==> form, powder,q'
!     ier_msg (2) = 'within the output menu to define the axis' 
!     ier_num = -125
!     ier_typ = ER_APPL
!     DEALLOCATE(pow_tmp,stat = all_status)  ! DeAllocate array for powder pattern copy
!     DEALLOCATE(xpl    ,stat = all_status)  ! DeAllocate array for calculated powder pattern
!     DEALLOCATE(ypl    ,stat = all_status)  ! DeAllocate array for calculated powder pattern
!     DEALLOCATE(lpv    ,stat = all_status)  ! DeAllocate array for calculated powder pattern
!     RETURN
!  ENDIF axq               ! pow_axis      == ??
ENDIF prsq  !Prepare S(Q), F(Q)
!
rmax     = out_user_values(2)
!
IF( cpow_form == 'tth' ) THEN
!        IF ( pow_axis      == POW_AXIS_Q  .or.  &        ! Non matching form, spline onto equidistant steps
!             pow_four_type == POW_DEBYE            ) THEN ! DEBYE, always spline
            IF(out_user_limits .AND. value /= val_pdf) THEN ! User provided values, not for PDF
               pow_tthmin   = out_user_values(1)
               pow_tthmax   = out_user_values(2)
               pow_deltatth = out_user_values(3)
            ELSE                                          ! Convert q limits
!              IF ( pow_axis == POW_AXIS_Q) THEN             ! Convert q limits to 2Theta
                  arg        = xmin/REAL(zpi) * rlambda / 2. ! Directly with arg in asind()
                  pow_tthmin = 2.*asind(arg)                 ! results in error ??????????
                  arg        = MIN(1.0D0,xmax/REAL(zpi) * rlambda / 2.)
                  pow_tthmax = 2.*asind(arg)
                  pow_deltatth = xpl(2)-xpl(1)
!              ENDIF
            ENDIF
            arg = 2.*asind(xpl(1)/REAL(zpi) * rlambda / 2.)
            IF(pow_tthmin <          arg ) THEN           ! User lower limit too low!
!           IF(pow_tthmin < xpl(1) ) THEN                 ! User lower limit too low!
               tthmin =              (INT( (         arg        )/pow_deltatth) + 1)*pow_deltatth
            ELSE
               tthmin = pow_tthmin
            ENDIF
            arg = 2.*asind((MIN(1.0, xpl(npkt)/REAL(zpi) * rlambda / 2.)))
            IF(pow_tthmax > arg       ) THEN              ! User upper limit too high!
               tthmax =              (INT( ( arg                )/pow_deltatth) - 1)*pow_deltatth
            ELSE
               tthmax = pow_tthmax
            ENDIF
            xmin = tthmin                                  ! Adjust limits needed later to cut 
            xmax = tthmax                                  ! off rounding errors
            npkt_equi =     INT((tthmax-tthmin)/pow_deltatth) + 1             
            ALLOCATE(y2a (1:POW_MAXPKT),stat = all_status) ! Allocate array for calculated powder pattern
            ALLOCATE(xwrt(1:npkt_equi),stat = all_status)  ! Allocate array for powder pattern ready to write
            ALLOCATE(ywrt(1:npkt_equi),stat = all_status)  ! Allocate array for powder pattern ready to write
            xwrt = 0.0
            ywrt = 0.0
            y2a  = 0.0
            CALL spline (npkt, xpl, ypl, 1.e31, 1.e31, y2a)
            DO ii = 1, npkt_equi
               xequ = tthmin + (ii-1)*pow_deltatth
               CALL splint (npkt, xpl, ypl, y2a, xequ, yequ, ier_num)
               IF(ier_num/=0) THEN
                  DEALLOCATE( pow_tmp, stat = all_status)
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
!        ELSE                                              ! Matching form no spline needed
!           ALLOCATE(xwrt(1:npkt     ),stat = all_status)  ! Allocate array for powder pattern ready to write
!           ALLOCATE(ywrt(1:npkt     ),stat = all_status)  ! Allocate array for powder pattern ready to write
!           xwrt = 0.0
!           ywrt = 0.0
!           DO ii = 1,npkt
!              xwrt(ii) = xpl(ii)
!              ywrt(ii) = ypl(ii)
!           ENDDO
!           npkt_wrt = npkt
!        ENDIF                   ! pow_axis      == ??
ELSEIF( cpow_form == 'q' ) THEN                       ! axis is Q
!        IF ( pow_axis      == POW_AXIS_TTH  .or.  &        ! Non matching form, spline onto equidistant steps
!            ((pow_four_type == POW_COMPL) .AND. value == val_pdf) .OR. &
!             pow_four_type == POW_DEBYE              ) THEN ! DEBYE, always spline
            IF(value == val_pdf) THEN                       ! Set limits for PDF
               CONTINUE
!              IF(pow_axis      == POW_AXIS_TTH) THEN
!                 pow_qmin   = REAL(zpi)*2/rlambda*sind(xmin)
!                 pow_qmax   = REAL(zpi)*2/rlambda*sind(xmax)
!                 pow_deltaq = xpl(npkt)-xpl(npkt-1)
!              ENDIF
            ELSE                                            ! Set limits for powder pattern
               IF(out_user_limits) THEN                     ! User provided values
                  pow_qmin   = out_user_values(1)
                  pow_qmax   = out_user_values(2)
                  pow_deltaq = out_user_values(3)
               ELSE                                          ! Convert q limits
               CONTINUE
!                 IF(pow_axis      == POW_AXIS_TTH) THEN
!                    pow_qmin   = REAL(zpi)*2/rlambda*sind(xmin)
!                    pow_qmax   = REAL(zpi)*2/rlambda*sind(xmax)
!                    pow_deltaq = xpl(npkt)-xpl(npkt-1)
!                 ENDIF
               ENDIF
            ENDIF
            IF(pow_qmin < xpl(1) ) THEN                     ! User lower limit too low!
               qmin =            (INT( (xpl(1)            )/pow_deltaq) + 1)*pow_deltaq
            ELSE
               qmin = pow_qmin
            ENDIF
            IF(pow_qmax > xpl(npkt) ) THEN                  ! User upper limit too high!
               qmax =            (INT( (         xpl(npkt))/pow_deltaq) - 1)*pow_deltaq
            ELSE
               qmax = pow_qmax
!     IF(value == val_pdf) THEN
!        IF(pdf_clin_a>0.0 .OR. pdf_cquad_a>0.0) THEN
!           qmax = MAX(30.0D0, pow_qmax_u*1.2D0)
!        ENDIF
!     ENDIF
            ENDIF
            xmin =   qmin                                  ! Adjust limits needed later to cut 
            xmax =   qmax                                  ! off rounding errors
            npkt_equi =     NINT((qmax-qmin)/pow_deltaq) + 1             
            ALLOCATE(y2a (1:POW_MAXPKT),stat = all_status) ! Allocate array for calculated powder pattern
            ALLOCATE(xwrt(1:npkt_equi),stat = all_status)  ! Allocate array for powder pattern ready to write
            ALLOCATE(ywrt(1:npkt_equi),stat = all_status)  ! Allocate array for powder pattern ready to write
            xwrt = 0.0
            ywrt = 0.0
            y2a  = 0.0
            CALL spline (npkt, xpl, ypl, 1.e31, 1.e31, y2a)
            DO ii = 1, npkt_equi
               xequ = qmin + (ii-1)*pow_deltaq
               CALL splint (npkt, xpl, ypl, y2a, xequ, yequ, ier_num)
               IF(ier_num/=0) THEN
                  DEALLOCATE( pow_tmp, stat = all_status)
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
!        ELSE                                              ! Matching form no spline needed
!           ALLOCATE(xwrt(1:npkt     ),stat = all_status)  ! Allocate array for powder pattern ready to write
!           ALLOCATE(ywrt(1:npkt     ),stat = all_status)  ! Allocate array for powder pattern ready to write
!           DO ii = 1,npkt
!              xwrt(ii) = xpl(ii)
!              ywrt(ii) = ypl(ii)
!           ENDDO
!           npkt_wrt = npkt
!        ENDIF                   ! pow_axis      == ??
ELSE                    ! cpow_form == 
         DO ii = 1,npkt
            xwrt(ii) = xpl(ii)
            ywrt(ii) = ypl(ii)
         ENDDO
         npkt_wrt = npkt
ENDIF                   ! cpow_form == 
!
      cut: DO
         IF(xwrt(npkt_wrt) > xmax) THEN
            npkt_wrt = npkt_wrt-1  ! Truncate in case of rounding errors
         ELSEIF(xwrt(npkt_wrt) < xwrt(npkt_wrt-1)) THEN
            npkt_wrt = npkt_wrt-1  ! Truncate in case of rounding errors
         ELSE
           EXIT cut
         ENDIF
      ENDDO cut
!
!     Scale intensity and add a background
!
IF(value==val_inten) THEN
         DO ii=1,npkt_wrt
            ywrt(ii) = pow_scale*ywrt(ii)
            DO iii=0,pow_nback
               ywrt(ii) = ywrt(ii) + pow_back(iii)*xwrt(ii)**iii
            ENDDO
         ENDDO
ELSEIF(value==val_pdf) THEN    ! Transform F(Q) into PDF
!
! treat for correlated motion effects
!
   corr_if: IF(pdf_clin_a>0.0 .OR. pdf_cquad_a>0.0) THEN
      IF(out_user_limits) THEN
         rmin     = 0.50 ! out_user_values(1)
         rmax     = out_user_values(2)*1.25
         rstep    = out_user_values(3)
         npkt_pdf = NINT((out_user_values(2)*1.25-0.5000000000000000)/out_user_values(3)) + 1
!        npkt_pdf = NINT((out_user_values(2)*1.25-out_user_values(1))/out_user_values(3)) + 1
      ELSE
         rmin     = pdf_rminu
         rmax     = pdf_rmaxu*1.25
         rstep    = pdf_deltaru
         npkt_pdf = NINT((rmax-rmin)/pdf_deltaru) + 1
      ENDIF
      rstep = REAL((rmax-rmin)/(npkt_pdf-1), KIND=PREC_DP)
      ALLOCATE(xfour(1:npkt_pdf))
      ALLOCATE(yfour(1:npkt_pdf))
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
!     j = ii - 1
!      peak_last1: DO ii = j,2, -1
!         IF((ywrt(ii-1)-ywrt(ii))*(ywrt(ii-2)-ywrt(ii-1))<=0.0) EXIT peak_last1
!      ENDDO peak_last1
!      j = ii-1
!      j = npkt_wrt - 600
!      smooth: DO ii = j+1, npkt_wrt !+ (npkt_wrt-ii)
!         ywrt(ii) = ywrt(j)*exp(-(REAL((ii-j))/(0.25*REAL((npkt_wrt-j )))))
!      ENDDO smooth
!write(*,*) 'found peak in FQ at ', xwrt(j), ywrt(j)
!
!open(77,file='POWDER/prae_corrlin.FQ',status='unknown')
!DO ii=1,npkt_wrt
!write(77,'(2(2x,G17.7E3))') xwrt(ii), ywrt(ii)
!enddo
!close(77)
!CORR write(*,*) ' ABOUT TO DO 1st FFT INTO PDF ', npkt_fft, npkt_pdf
      CALL fft_fq(npkt_wrt, xwrt, ywrt, rmin, rmax, rstep, npkt_fft, npkt_pdf, xfour, yfour)
!write(*,*) ' after       1st FFT INTO PDF ', npkt_pdf
!open(77,file='POWDER/prae_corrlin.PDF',status='unknown')
!DO ii=1,npkt_pdf
!write(77,'(2(2x,G17.7E3))') xfour(ii), yfour(ii)
!enddo
!close(77)
!write(*,*) ' <u^2> , B', pow_u2aver, pow_u2aver*8.*3.1415**2
      sigma = 2.0*(pow_u2aver*u2aver_scale)              ! TO BE REPLACED BY ATOMIC B VALUE
      CALL powder_conv_corrlin(yfour, REAL(rmin),REAL(rmax), REAL(rstep),   &
                               sigma, pdf_clin_a, pdf_cquad_a, pdf_rcut, pow_width,   &
                               POW_MAXPKT)
!open(77,file='POWDER/post_corrlin.PDF',status='unknown')
!DO ii=1,npkt_pdf
!write(77,'(2(2x,G17.7E3))') xfour(ii), yfour(ii)
!enddo
!close(77)
      CALL fft_fq(npkt_pdf, xfour, yfour, pow_qmin_u, pow_qmax_u, pow_deltaq_u, npkt_fft, npkt_wrt, xwrt, ywrt)
!open(77,file='POWDER/post_corrlin.FQ',status='unknown')
!DO ii=1,npkt_wrt
!write(77,'(2(2x,G17.7E3))') xwrt(ii), ywrt(ii)
!enddo
!close(77)
!     CALL fft_fq(npkt_wrt, xwrt, ywrt, rmin, rmax, rstep, npkt_fft, npkt_pdf, xfour, yfour)
!      DEALLOCATE(xfour)
!      DEALLOCATE(yfour)
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
      ALLOCATE(xfour(1:npkt_pdf))
      ALLOCATE(yfour(1:npkt_pdf))
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
      CALL  fft_fq(npkt_wrt, xwrt, ywrt, rmin, rmax, rstep, npkt_fft, npkt_pdf, xfour, yfour)
!write(*,*) ' RLIMITS ', rmin, rmax, rstep, npkt_pdf
!open(77,file='POWDER/ende_corrlin.PDF',status='unknown')
!DO ii=1,npkt_pdf
!write(77,'(2(2x,G17.7E3))') xfour(ii), yfour(ii)
!enddo
!close(77)
!        sigma = 2.0*pow_u2aver                    ! TO BE REPLACED BY ATOMIC B VALUE
!        IF(pdf_clin_a>0.0 .OR. pdf_cquad_a>0.0) THEN
!           CALL powder_conv_corrlin(yfour, REAL(rmin),REAL(rmax), REAL(rstep),   &
!                                    sigma, pdf_clin_a, pdf_cquad_a, pow_width,   &
!                                    POW_MAXPKT)
!        ENDIF
   ENDIF corr_if
   DEALLOCATE(xwrt)
   DEALLOCATE(ywrt)
   ALLOCATE(xwrt(1:npkt_pdf))
   ALLOCATE(ywrt(1:npkt_pdf))
   DO ii=1,npkt_pdf
      xwrt(ii) =xfour(ii)
      ywrt(ii) =yfour(ii)
   ENDDO
   DEALLOCATE(xfour)
   DEALLOCATE(yfour)
   npkt_wrt = npkt_pdf
!
ENDIF
!
!     Finally write the pattern
!
CALL powder_do_write (outfile, npkt_wrt, xwrt, ywrt)
!
DEALLOCATE( pow_tmp, stat = all_status)
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
!
IMPLICIT none 
!                                                                       
INTEGER :: npkt
INTEGER :: j
INTEGER :: nzero      ! number of zeros in intensity
!
REAL(KIND=PREC_SP) :: xmin, xmax, xdel, xxmax
REAL(KIND=PREC_SP) :: scalef
REAL(KIND=PREC_SP) :: pow_tmp_sum
REAL(KIND=PREC_SP) :: pow_uuu_sum
!
REAL           :: ss       ! time
REAL, EXTERNAL :: seknds
!
WRITE (output_io, * ) ' Starting convolution'
ss = seknds (0.0)
!
xmin  = 0.0 
xmax  = 0.0 
xxmax = 0.0 
xdel  = 0.0 
IF (pow_four_type.eq.POW_COMPL) THEN 
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
IF (pow_four_type.ne.POW_COMPL) THEN 
!                                                              
!     This is a Debye calculation, copy rsf or csf into pow_conv
!                                                              
   IF (npkt    .le.POW_MAXPKT) THEN 
      DO j = 1, npkt    
         pow_conv (j) = REAL(rsf (j) )     ! Double precision no longer needed
      ENDDO 
   ENDIF 
ELSE
   pow_conv(:) = 2.0*REAL(pow_qsp(:))   ! Double precision no longer needed, Correct for half space
ENDIF 
!open(77,file='POWDER/prae_conv.FQ',status='unknown')
!DO j=0,npkt
!  write(77,'(2(2x,G17.7E3))')     xmin+j *xdel, pow_conv(j)
!enddo
!close(77)
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
ELSEIF (pow_profile == POW_PROFILE_PSVGT) THEN 
   IF(pow_u/=0.0 .OR. pow_v/=0.0 .OR. pow_eta_l/=0.0 .OR. pow_eta_q/=0.0) THEN
      xxmax = xmax + xdel
      IF(pow_p1/=0.0 .OR. pow_p2/=0.0 .OR.                 &
         pow_p3/=0.0 .OR. pow_p4/=0.0            ) THEN       
         CALL powder_conv_psvgt_uvw_asym(pow_conv, xmin,xxmax, xdel,   &
         pow_eta, pow_eta_l, pow_eta_q, pow_u, pow_v, pow_w, pow_p1, pow_p2,  &
         pow_p3, pow_p4, pow_width, POW_MAXPKT, pow_four_type, pow_axis)
      ELSE          ! Symmetric case
         CALL powder_conv_psvgt_uvw(pow_conv, xmin,xxmax, xdel,   &
         pow_eta, pow_eta_l, pow_eta_q, pow_u, pow_v, pow_w, pow_width,  &
         POW_MAXPKT, pow_four_type, pow_axis)
      ENDIF
   ELSE 
      xxmax = xmax + xdel
      CALL powder_conv_psvgt_fix (pow_conv, xmin,xxmax, xdel,   &
      pow_eta, pow_w, pow_width, POW_MAXPKT, pow_four_type)
   ENDIF 
!  lconv = .TRUE.
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
ELSEIF(pow_four_type==POW_COMPL) THEN
   scalef = 1./xdel
ENDIF
pow_conv(:) = pow_conv(:) * scalef
ss = seknds (ss) 
WRITE (output_io, '(/,'' Elapsed time    : '',G12.6,'' sec'')') ss
!open(77,file='POWDER/post_conv.FQ',status='unknown')
!DO j=0,npkt
!  write(77,'(2(2x,G17.7E3))')     xmin+j *xdel, pow_conv(j)
!enddo
!close(77)
!
END SUBROUTINE powder_convolute                     
!
!*******************************************************************************
!
SUBROUTINE four_fq(npkt_wrt, xwrt, ywrt, rmin, rmax, npkt_pdf, xfour, yfour)
!
USE wink_mod
!
INTEGER                       , INTENT(IN) :: npkt_wrt
REAL   , DIMENSION(1:npkt_wrt), INTENT(IN) :: xwrt
REAL   , DIMENSION(1:npkt_wrt), INTENT(IN) :: ywrt
REAL(KIND=PREC_DP)            , INTENT(IN) :: rmin, rmax
INTEGER                       , INTENT(IN) :: npkt_pdf
REAL   , DIMENSION(1:npkt_pdf), INTENT(OUT) :: xfour
REAL   , DIMENSION(1:npkt_pdf), INTENT(OUT) :: yfour
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
   xfour(i) = REAL(rr)
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
!
!*******************************************************************************
!
SUBROUTINE  fft_fq(npkt_wrt, xwrt, ywrt, rmin, rmax, rstep, npkt_fft, npkt_pdf, xfour, yfour)
!
USE fast_fourier_mod
USE spline_mod
USE wink_mod
!
INTEGER                       , INTENT(IN)  :: npkt_wrt
REAL   , DIMENSION(1:npkt_wrt), INTENT(IN)  :: xwrt
REAL   , DIMENSION(1:npkt_wrt), INTENT(IN)  :: ywrt
REAL(KIND=PREC_DP)            , INTENT(IN)  :: rmin, rmax, rstep
INTEGER                       , INTENT(IN)  :: npkt_fft !points in powder pattern for Fast Fourier = 2**16
INTEGER                       , INTENT(IN)  :: npkt_pdf
REAL   , DIMENSION(1:npkt_pdf), INTENT(OUT) :: xfour
REAL   , DIMENSION(1:npkt_pdf), INTENT(OUT) :: yfour
!
INTEGER   :: i
INTEGER   :: lensav      ! Array size for ip
INTEGER   :: lenwrk      ! Array size for w
INTEGER   :: iqmin
INTEGER   :: iqmax
INTEGER   :: irmin
!REAL(KIND=PREC_DP) :: rstep
REAL(KIND=PREC_DP) :: dq
!REAL(KIND=PREC_DP) :: qmin
REAL(KIND=PREC_DP) :: qmax
REAL(KIND=PREC_DP), DIMENSION(:), ALLOCATABLE :: temp   ! Temporary intensities for FFT
INTEGER           , DIMENSION(:), ALLOCATABLE :: ip
REAL(KIND=PREC_DP), DIMENSION(:), ALLOCATABLE :: w
REAL(KIND=PREC_SP), DIMENSION(:), ALLOCATABLE :: xfft   ! Temporary array for FFT result
REAL(KIND=PREC_SP), DIMENSION(:), ALLOCATABLE :: yfft   ! Temporary array for FFT result
!
!write(*,*) ' PDF OUT ', npkt_wrt, npkt_fft, npkt_pdf, rmin, rmax, (rmax-rmin) / REAL((npkt_pdf-1), KIND=PREC_DP)
!write(*,*) ' PDF ste ', (rmax-rmin) / REAL((npkt_pdf-1), KIND=PREC_DP), rstep
!write(*,*) ' PDF OUT ', xwrt(1), xwrt(2), xwrt(npkt_wrt)
!write(*,*) ' PDF OUT ', ywrt(1), ywrt(2), ywrt(npkt_wrt)
!n
dq = (xwrt(npkt_wrt)-xwrt(1))/(npkt_wrt-1)
!rstep    = (rmax-rmin) / REAL((npkt_pdf-1), KIND=PREC_DP)
qmax     =  PI/rstep           ! by using 2*PI/rstep, the step size in direct space is halved
!
iqmin = NINT(xwrt(1       )/dq)
iqmax = NINT(xwrt(npkt_wrt)/dq)
lensav= 4+INT(SQRT(FLOAT(npkt_fft)/2))
lenwrk= npkt_fft*5/4-1
ALLOCATE(temp(0:npkt_fft-1))
ALLOCATE(ip(0:lensav))
ALLOCATE(w (0:lenwrk))
ALLOCATE(xfft(1:npkt_fft+1))
ALLOCATE(yfft(1:npkt_fft+1))
temp = 0.0D0
ip   = 0
w    = 0.0D0
DO i=0,iqmin-1                       ! Augment straight line to q=0
   temp(i) = REAL(i,KIND=PREC_DP)*dq*ywrt(1)/xwrt(1)
ENDDO
temp(iqmin:iqmax) = ywrt(1:npkt_wrt) ! Add actual powder pattern
!open(77,file='POWDER/prae_fft.FQ',status='unknown')
!DO i=0,iqmax
!  write(77,'(2(2x,G17.7E3))') (i)*dq,temp(i)
!enddo
!close(77)
!
!write(*,*) ' DO  FOURIER '
CALL ddst(npkt_fft, 1, temp, ip, w)
!write(*,*) ' DID FOURIER '
!
irmin = nint(rmin/rstep)
!write(*,*) 'PDF OUT ', rmin, rstep, irmin, npkt_pdf
!write(*,*) 'PDF_OUT ', dq, qmax, (npkt_fft-1)*dq
qmax = (npkt_fft-1)*dq
!write(*,*) 'PDF_OUT ', qmax, rstep, PI/qmax
!open(77,file='POWDER/fft.PDF',status='unknown')
!DO i=0,npkt_pdf+irmin
!  write(77,'(2(2x,G17.7E3))') ( 0.5+i)*PI/QMAX,temp(i)*2/PI*dq  ! -0.5
!enddo
!close(77)
!open(77,file='POWDER/fft_2.PDF',status='unknown')
DO i=1,npkt_fft
  xfft (i) = (i-0.50)*PI/qmax
  yfft (i) = temp(i-1)*2/PI*dq
!  write(77,'(2(2x,G17.7E3))') xfft(i), yfft(i)
ENDDO
!close(77)
!write(*,*) 'DO SPLINE ', REAL(rmin), REAL(rmax), REAL(rstep), npkt_pdf
CALL spline_prep(npkt_fft, xfft, yfft, REAL(rmin), REAL(rmax), REAL(rstep), npkt_pdf, xfour, yfour)
!
DEALLOCATE(temp)
DEALLOCATE(ip)
DEALLOCATE(w)
DEALLOCATE(xfft)
DEALLOCATE(yfft)
!
END SUBROUTINE  fft_fq
!*****7*****************************************************************
      REAL function lorentz (ttheta, flag_fq) 
!+                                                                      
!-                                                                      
      USE discus_config_mod 
      USE powder_mod 
      USE trig_degree_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      REAL   , INTENT(IN) :: ttheta 
      INTEGER, INTENT(IN) :: flag_fq
!                                                                       
!
!     REAL sind
!                                                                       
      lorentz = 1.0
      
      IF (pow_four_type.eq.POW_DEBYE) THEN 
         lorentz = 1.0 
      ELSE 
         IF(flag_fq==0) THEN
         IF (pow_lp.eq.POW_LP_BRAGG) THEN 
            lorentz = 0.5 / sind (0.5 * ttheta) / sind (ttheta) 
         ELSEIF (pow_lp.eq.POW_LP_NEUT) THEN 
            lorentz = 0.5 / sind (0.5 * ttheta) / sind (ttheta) 
         ELSEIF (pow_lp.eq.POW_LP_NONE) THEN 
            lorentz = 1.0 
         ELSEIF (pow_lp.eq.POW_LP_SYNC) THEN 
            lorentz = 0.5 / sind (0.5 * ttheta) / sind (ttheta) 
         ENDIF 
         ELSEIF(flag_fq==1) THEN 
            lorentz = 0.5 / sind (0.5 * ttheta) / sind (ttheta) 
         ENDIF 
      ENDIF 
!                                                                       
      END FUNCTION lorentz                          
!*****7*****************************************************************
      REAL FUNCTION polarisation (ttheta) 
!+                                                                      
!-                                                                      
      USE discus_config_mod 
      USE powder_mod 
      USE trig_degree_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      REAL ttheta 
!
!
!     REAL cosd 
!                                                                       
      polarisation = 1.0
      
      IF (pow_lp.eq.POW_LP_BRAGG) THEN 
         polarisation = (1. + (cosd (ttheta) ) **2 * pow_lp_fac)        &
         / (1. + pow_lp_fac)                                            
      ELSEIF (pow_lp.eq.POW_LP_NEUT) THEN 
         polarisation = 1.0 
      ELSEIF (pow_lp.eq.POW_LP_NONE) THEN 
         polarisation = 1.0 
      ELSEIF (pow_lp.eq.POW_LP_SYNC) THEN 
         polarisation = pow_lp_fac + (1. - pow_lp_fac) * (cosd (ttheta) &
         ) **2 * pow_lp_cos                                             
      ENDIF 
!                                                                       
      END FUNCTION polarisation                     
!*****7*****************************************************************
      REAL FUNCTION lorentz_pol (ttheta) 
!+                                                                      
!-                                                                      
      USE discus_config_mod 
      USE powder_mod 
      USE trig_degree_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      REAL ttheta 
!                                                                       
!     REAL sind, cosd 
!                                                                       
      IF (pow_four_type.eq.POW_DEBYE) THEN 
         lorentz_pol = 1.0 
      ELSE
         lorentz_pol = (1-pow_lp_fac+pow_lp_fac*(cosd(pow_lp_ang))**2*(cosd(ttheta))**2)/ &
                       (2.*(sind(0.5*ttheta))**2*cosd(0.5*ttheta))
      ENDIF
      END FUNCTION lorentz_pol
!*****7*****************************************************************
      SUBROUTINE powder_conv_res (dat, tthmin, tthmax, dtth, delta, POW_MAXPKT)
!-                                                                      
!     Convolute powder pattern with resolution function (Gaussian)      
!+                                                                      
      USE discus_config_mod 
      USE wink_mod
      IMPLICIT none 
!                                                                       
!
      INTEGER, INTENT(IN) :: POW_MAXPKT
!                                                                       
      REAL dat (0:POW_MAXPKT) 
      REAL tthmin, tthmax, dtth, delta 
!                                                                       
      REAL dummy (0:POW_MAXPKT) 
      REAL gauss (0:2 * POW_MAXPKT) 
      REAL tth
      INTEGER imax, i, j, ii 
      INTEGER max_ps 
!                                                                       
!------ Setup Gaussian                                                  
!                                                                       
      max_ps = int( (10.0 * delta) / dtth )
      DO i = 0, max_ps 
      tth = i * dtth 
      gauss (i) = 1.0 / sqrt (REAL(pi)) / delta * exp ( - (tth**2 / delta**2))  
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
REAL, DIMENSION(0:POW_MAXPKT), INTENT(INOUT) ::  dat !(0:POW_MAXPKT) 
REAL                         , INTENT(IN)    ::  tthmin, tthmax, dtth, eta
REAL                         , INTENT(IN)    ::  w 
REAL                         , INTENT(IN)    ::  pow_width 
INTEGER                      , INTENT(IN)    ::  pow_type   ! = 0==COMPl =1==DEBYE
!
INTEGER, PARAMETER  :: POW_COMPL = 0
INTEGER, PARAMETER  :: POW_DEBYE = 1
!                                                                       
REAL                               :: fwhm
REAL, DIMENSION(0:POW_MAXPKT)      :: dummy
REAL, DIMENSION(0:2 * POW_MAXPKT)  :: psvgt
!REAL                               :: tth
!REAL(KIND=PREC_DP)                 :: tth_dp
REAL(KIND=PREC_DP)                 :: eta_dp
REAL(KIND=PREC_DP)                 :: fwhm_dp
INTEGER                            :: imax, i, j, ii 
INTEGER                            :: max_ps 
!                                                                       
!REAL pseudovoigt 
!                                                                       
!------ Setup Pseudo-Voigt                                              
!                                                                       
fwhm = sqrt (abs (w) ) 
max_ps = int( (pow_width * fwhm) / dtth )
psvgt = 0.0
eta_dp  = REAL(eta,KIND=PREC_DP)
fwhm_dp = REAL(fwhm,KIND=PREC_DP)
!open(77,file='pseudo.comp', status='unknown')
DO i = 0, max_ps 
   ii = MIN(INT(i*dtth/fwhm*GLP_NPT), GLP_MAX)
   psvgt (i) = glp_pseud_indx(ii  , eta_dp, fwhm_dp)
ENDDO 
!                                                                       
!------ Now convolute                                                   
!                                                                       
imax = int( (tthmax - tthmin) / dtth )
dummy = 0.0
!write(*,*) ' IMAX ', imax, MINVAL(dat), maxval(dat)
!write(*,*) ' TYPE ', pow_type, POW_COMPL, POW_DEBYE
!
!------ If COMPLETE, check for zeros, else do all points
!
IF(pow_type==POW_COMPL) THEN
   DO i = 0, imax 
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
!*****7*****************************************************************
SUBROUTINE powder_conv_psvgt_uvw (dat, tthmin, tthmax, dtth, eta0,&
      eta_l, eta_q, u, v, w, pow_width, POW_MAXPKT, pow_type, axis)
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
USE trig_degree_mod
USE wink_mod
!
IMPLICIT none 
!
INTEGER, INTENT(IN) :: POW_MAXPKT
REAL   , DIMENSION(0:POW_MAXPKT), INTENT(INOUT) :: dat   ! Data to be convoluted
REAL                            , INTENT(IN)    :: tthmin ! 2Theta min
REAL                            , INTENT(IN)    :: tthmax ! 2Theta max
REAL                            , INTENT(IN)    :: dtth   ! 2Theta step
REAL                            , INTENT(IN)    :: eta0   ! Lor/Gaus mix constant part
REAL                            , INTENT(IN)    :: eta_l   ! Lor/Gaus mix variable part
REAL                            , INTENT(IN)    :: eta_q   ! Lor/Gaus mix variable part
REAL                            , INTENT(IN)    :: u      ! u*tan^2(Theta)
REAL                            , INTENT(IN)    :: v      ! v*tan  (Theta)
REAL                            , INTENT(IN)    :: w      ! w
REAL                            , INTENT(IN)    :: pow_width ! Number of FWHM's to calculate
INTEGER                         , INTENT(IN)    :: pow_type  ! == 0 for COMPLETE, ==1 for Debye
INTEGER                         , INTENT(IN)    :: axis   ! == 2 for 2theta, == 1 for Q
!
INTEGER, PARAMETER  :: POW_COMPL = 0
INTEGER, PARAMETER  :: POW_DEBYE = 1
!
REAL(KIND=PREC_DP)            :: fwhm     ! Current FWHM at Theta
REAL, DIMENSION(0:POW_MAXPKT) :: dummy    ! temporary data (0:POW_MAXPKT) 
REAL(KIND=PREC_DP)            :: tth      ! Theta within convolution, main data set
REAL                          :: tantth   ! tan(Theta)
REAL(KIND=PREC_DP)    :: eta              ! actual eta at current 2Theta
INTEGER :: imax, i, j, ii  ! Dummy loop indices
INTEGER :: i1, i2          ! Pseudo Voigt lookup indices
REAL    :: pseudo          ! scale factor for lookup table
INTEGER :: max_ps 
!                                                                       
!------ Now convolute                                                   
!                                                                       
imax = INT( (tthmax - tthmin) / dtth )
!
dummy = 0.0   ! dummy(:)
IF(pow_type==POW_COMPL) THEN     ! Complete, check for zeros in DAT
   main_pts: DO i = 0, imax 
      tth = tthmin + i * dtth 
      tantth = tand (tth * 0.5) 
!
      fwhm = 0.00001
      IF(axis==2 ) THEN       ! 2Theta axis
         fwhm = SQRT (MAX (ABS (u * tantth**2 + v * tantth + w), 0.00001) ) 
      ELSEif(axis==1) THEN
         fwhm = SQRT( MAX( ABS( u*tth**2 + v*tth + w), 0.00001D0) )
      ENDIF
!
      max_ps = INT((pow_width * fwhm) / dtth )
      eta = MIN(1.0D0, MAX(0.0D0, eta0 + eta_l * tth + eta_q*tth**2) ) 
      pseudo =     dtth/fwhm*glp_npt  ! scale factor for look up table
      i1 = 0                               ! == 0 * dtth
      i2 = MIN(INT(2*i*pseudo), GLP_MAX)   ! == 2*i*dtth
      dummy (i) = dat (i) * ( glp_pseud_indx(i1, eta, fwhm)  &
                             -glp_pseud_indx(i2, eta, fwhm))                             
!
      ii = MAX (i - 1 - max_ps + 1, 0) 
      first:DO j = ii, i - 1 
         IF(dat(j)==0.0) CYCLE first
         i1 = MIN(INT((i - j) * pseudo), GLP_MAX)  ! == tth1 = (i - j) * dtth
         i2 = MIN(INT((i + j) * pseudo), GLP_MAX)  ! == tth2 = (i + j) * dtth
         dummy(i) = dummy(i) + dat(j) *( glp_pseud_indx(i1, eta, fwhm)  &
                                        -glp_pseud_indx(i2, eta, fwhm))                    
      ENDDO first
!
      ii = MIN(i + 1 + max_ps - 1, imax) 
      secnd: DO j = i + 1, ii 
         IF(dat(j)==0.0) CYCLE secnd
         i1 = MIN(INT((j - i) * pseudo), GLP_MAX)  ! == tth1 = (j - i) * dtth
         i2 = MIN(INT((j + i) * pseudo), GLP_MAX)  ! == tth1 = (j + i) * dtth
         dummy(i) = dummy(i) + dat(j) *( glp_pseud_indx(i1, eta, fwhm)  &
                                        -glp_pseud_indx(i2, eta, fwhm))
      ENDDO secnd
   ENDDO main_pts
ELSEIF(pow_type==POW_DEBYE) THEN     ! DEBYE, do not check for zeros in DAT
   main_pts_deb: DO i = 0, imax 
      tth = tthmin + i * dtth 
      tantth = tand (tth * 0.5) 
!
      fwhm = 0.00001
      IF(axis==2 ) THEN       ! 2Theta axis
         fwhm = SQRT (MAX (ABS (u * tantth**2 + v * tantth + w), 0.00001) ) 
      ELSEif(axis==1) THEN
         fwhm = SQRT( MAX( ABS( u*tth**2 + v*tth + w), 0.00001D0) )
      ENDIF
!
      max_ps = INT((pow_width * fwhm) / dtth )
      eta = MIN(1.0D0, MAX(0.0D0, eta0 + eta_l * tth + eta_q*tth**2) ) 
      pseudo =     dtth/fwhm*glp_npt  ! scale factor for look up table
      i1 = 0                               ! == 0 * dtth
      i2 = MIN(INT(2*i*pseudo), GLP_MAX)   ! == 2*i*dtth
      dummy (i) = dat (i) * ( glp_pseud_indx(i1, eta, fwhm)  &
                             -glp_pseud_indx(i2, eta, fwhm))                             
!
      ii = MAX (i - 1 - max_ps + 1, 0) 
      first_deb:DO j = ii, i - 1 
         i1 = MIN(INT((i - j) * pseudo), GLP_MAX)  ! == tth1 = (i - j) * dtth
         i2 = MIN(INT((i + j) * pseudo), GLP_MAX)  ! == tth2 = (i + j) * dtth
         dummy(i) = dummy(i) + dat(j) *( glp_pseud_indx(i1, eta, fwhm)  &
                                        -glp_pseud_indx(i2, eta, fwhm))                    
      ENDDO first_deb
!
      ii = MIN(i + 1 + max_ps - 1, imax) 
      secnd_deb: DO j = i + 1, ii 
         i1 = MIN(INT((j - i) * pseudo), GLP_MAX)  ! == tth1 = (j - i) * dtth
         i2 = MIN(INT((j + i) * pseudo), GLP_MAX)  ! == tth1 = (j + i) * dtth
         dummy(i) = dummy(i) + dat(j) *( glp_pseud_indx(i1, eta, fwhm)  &
                                        -glp_pseud_indx(i2, eta, fwhm))
      ENDDO secnd_deb
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
SUBROUTINE powder_conv_psvgt_uvw_asym (dat, tthmin, tthmax, dtth, eta0,&
      eta_l, eta_q, u, v, w, p1, p2, p3, p4, pow_width, POW_MAXPKT,    &
      pow_type, axis)
!-
!     Convolute powder pattern with resolution function (Pseudo-Voigt)  
!     FWHM according to caglioti equation, Constant eta                 
!     FWHM = sqrt ( U*tan**2(Theta) + V*tan(Theta) + W)                 
!     Symmetric case
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
REAL   , DIMENSION(0:POW_MAXPKT), INTENT(INOUT) :: dat   ! Data to be convoluted
REAL                            , INTENT(IN)    :: tthmin ! 2Theta min
REAL                            , INTENT(IN)    :: tthmax ! 2Theta max
REAL                            , INTENT(IN)    :: dtth   ! 2Theta step
REAL                            , INTENT(IN)    :: eta0   ! Lor/Gaus mix constant part
REAL                            , INTENT(IN)    :: eta_l   ! Lor/Gaus mix variable part
REAL                            , INTENT(IN)    :: eta_q   ! Lor/Gaus mix variable part
REAL                            , INTENT(IN)    :: u      ! u*tan^2(Theta)
REAL                            , INTENT(IN)    :: v      ! v*tan  (Theta)
REAL                            , INTENT(IN)    :: w      ! w
REAL                            , INTENT(IN)    :: p1     ! Asymmetry terms p1 to P4
REAL                            , INTENT(IN)    :: p2     ! Asymmetry terms p1 to P4
REAL                            , INTENT(IN)    :: p3     ! Asymmetry terms p1 to P4
REAL                            , INTENT(IN)    :: p4     ! Asymmetry terms p1 to P4
REAL                            , INTENT(IN)    :: pow_width ! Number of FWHM's to calculate
INTEGER                         , INTENT(IN)    :: pow_type  ! == 0 for COMPLETE, ==1 for Debye
INTEGER                         , INTENT(IN)    :: axis   ! == 2 for 2theta, == 1 for Q
!
INTEGER, PARAMETER            :: POW_COMPL = 0
INTEGER, PARAMETER            :: POW_DEBYE = 1
!
REAL(KIND=PREC_DP)            :: fwhm     ! Current FWHM at Theta
REAL, DIMENSION(0:POW_MAXPKT) :: dummy    ! temporary data (0:POW_MAXPKT) 
REAL(KIND=PREC_DP)            :: tth      ! Theta within convolution, main data set
REAL                          :: tantth   ! tan(Theta)
REAL(KIND=PREC_DP)    :: eta              ! actual eta at current 2Theta
REAL(KIND=PREC_DP)    :: tth1            ! Theta values for asymmetry
REAL(KIND=PREC_DP)    :: tth2            ! Theta values for asymmetry
INTEGER :: imax, i, j, ii  ! Dummy loop indices
INTEGER :: i1, i2          ! Pseudo Voigt lookup indices
REAL    :: pseudo          ! scale factor for lookup table
REAL    :: pra1, pra2      ! Asymmetry pre factors
INTEGER :: max_ps 
!                                                                       
!------ Now convolute                                                   
!                                                                       
imax = INT( (tthmax - tthmin) / dtth )
dummy = 0.0
!
IF(pow_type==POW_COMPL) THEN
   main_compl: DO i = 0, imax 
      tth = tthmin + i * dtth 
      tantth = tand (tth * 0.5) 
!
      IF(axis==2 ) THEN       ! 2Theta axis
         fwhm = SQRT (MAX (ABS (u * tantth**2 + v * tantth + w), 0.00001) ) 
      ELSEif(axis==1) THEN
         fwhm = SQRT( MAX( ABS( u*tth**2 + v*tth + w), 0.00001D0) )
      ENDIF
!
      max_ps = INT((pow_width * fwhm) / dtth )
      eta = MIN(1.0D0, MAX(0.0D0, eta0 + eta_l * tth + eta_q*tth**2) ) 
      pseudo =     dtth/fwhm*glp_npt  ! scale factor for look up table
      i1 = 0                               ! == 0 * dtth
      i2 = MIN(INT(2*i*pseudo), GLP_MAX)   ! == 2*i*dtth
      tth1 = 0 * dtth 
      tth2 = 2 * i * dtth 
      pra1 = profile_asymmetry (tth, tth1, fwhm, p1, p2, p3, p4) 
      pra2 = profile_asymmetry (tth, tth2, fwhm, p1, p2, p3, p4) 
      dummy (i) = dat (i) * ( glp_pseud_indx(i1, eta, fwhm)*pra1  &
                             -glp_pseud_indx(i2, eta, fwhm)*pra2)                             
!
      ii = MAX (i - 1 - max_ps + 1, 0) 
      first_compl: DO j = ii, i - 1 
         IF(dat(j)==0.0) CYCLE first_compl
            i1 = MIN(INT((i - j) * pseudo), GLP_MAX)  ! == tth1 = (i - j) * dtth
         i2 = MIN(INT((i + j) * pseudo), GLP_MAX)  ! == tth2 = (i + j) * dtth
         tth1 = (i - j) * dtth 
         tth2 = (i + j) * dtth 
         pra1 = profile_asymmetry(tth, tth1, fwhm, p1, p2, p3, p4) 
         pra2 = profile_asymmetry(tth, tth2, fwhm, p1, p2, p3, p4) 
         dummy(i) = dummy(i) + dat(j) *( glp_pseud_indx(i1, eta, fwhm)*pra1  &
                                        -glp_pseud_indx(i2, eta, fwhm)*pra2)                    
      ENDDO first_compl
!
      ii = MIN(i + 1 + max_ps - 1, imax) 
      secnd_compl: DO j = i + 1, ii 
         IF(dat(j)==0.0) CYCLE secnd_compl
         i1 = MIN(INT((j - i) * pseudo), GLP_MAX)  ! == tth1 = (j - i) * dtth
         i2 = MIN(INT((j + i) * pseudo), GLP_MAX)  ! == tth1 = (j + i) * dtth
         tth1 = (j - i) * dtth 
         tth2 = (j + i) * dtth 
         pra1 = profile_asymmetry (tth, - tth1, fwhm, p1, p2, p3, p4) 
         pra2 = profile_asymmetry (tth, - tth2, fwhm, p1, p2, p3, p4) 
         dummy(i) = dummy(i) + dat(j) *( glp_pseud_indx(i1, eta, fwhm)*pra1  &
                                        -glp_pseud_indx(i2, eta, fwhm)*pra2)
      ENDDO secnd_compl
   ENDDO main_compl
ELSEIF(pow_type==POW_DEBYE) THEN
   main_debye: DO i = 0, imax 
      tth = tthmin + i * dtth 
      tantth = tand (tth * 0.5) 
!
      IF(axis==2 ) THEN       ! 2Theta axis
         fwhm = SQRT (MAX (ABS (u * tantth**2 + v * tantth + w), 0.00001) ) 
      ELSEif(axis==1) THEN
         fwhm = SQRT( MAX( ABS( u*tth**2 + v*tth + w), 0.00001D0) )
      ENDIF
!
      max_ps = INT((pow_width * fwhm) / dtth )
      eta = MIN(1.0D0, MAX(0.0D0, eta0 + eta_l * tth + eta_q*tth**2) ) 
      pseudo =     dtth/fwhm*glp_npt  ! scale factor for look up table
      i1 = 0                               ! == 0 * dtth
      i2 = MIN(INT(2*i*pseudo), GLP_MAX)   ! == 2*i*dtth
      tth1 = 0 * dtth 
      tth2 = 2 * i * dtth 
      pra1 = profile_asymmetry (tth, tth1, fwhm, p1, p2, p3, p4) 
      pra2 = profile_asymmetry (tth, tth2, fwhm, p1, p2, p3, p4) 
      dummy (i) = dat (i) * ( glp_pseud_indx(i1, eta, fwhm)*pra1  &
                             -glp_pseud_indx(i2, eta, fwhm)*pra2)                             
!
      ii = MAX (i - 1 - max_ps + 1, 0) 
      first_debye: DO j = ii, i - 1 
         i1 = MIN(INT((i - j) * pseudo), GLP_MAX)  ! == tth1 = (i - j) * dtth
         i2 = MIN(INT((i + j) * pseudo), GLP_MAX)  ! == tth2 = (i + j) * dtth
         tth1 = (i - j) * dtth 
         tth2 = (i + j) * dtth 
         pra1 = profile_asymmetry(tth, tth1, fwhm, p1, p2, p3, p4) 
         pra2 = profile_asymmetry(tth, tth2, fwhm, p1, p2, p3, p4) 
         dummy(i) = dummy(i) + dat(j) *( glp_pseud_indx(i1, eta, fwhm)*pra1  &
                                        -glp_pseud_indx(i2, eta, fwhm)*pra2)                    
      ENDDO first_debye
!
      ii = MIN(i + 1 + max_ps - 1, imax) 
      secnd_debye: DO j = i + 1, ii 
         i1 = MIN(INT((j - i) * pseudo), GLP_MAX)  ! == tth1 = (j - i) * dtth
         i2 = MIN(INT((j + i) * pseudo), GLP_MAX)  ! == tth1 = (j + i) * dtth
         tth1 = (j - i) * dtth 
         tth2 = (j + i) * dtth 
         pra1 = profile_asymmetry (tth, - tth1, fwhm, p1, p2, p3, p4) 
         pra2 = profile_asymmetry (tth, - tth2, fwhm, p1, p2, p3, p4) 
         dummy(i) = dummy(i) + dat(j) *( glp_pseud_indx(i1, eta, fwhm)*pra1  &
                                        -glp_pseud_indx(i2, eta, fwhm)*pra2)
      ENDDO secnd_debye
   ENDDO main_debye
ENDIF
!                                                                       
DO i = 0, imax 
   dat (i) = dummy (i) * dtth   ! scale with stepwidth
ENDDO 
!                                                                       
END SUBROUTINE powder_conv_psvgt_uvw_asym
!*****7*****************************************************************
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
USE trig_degree_mod
USE wink_mod
!
IMPLICIT none 
!
INTEGER, INTENT(IN) :: POW_MAXPKT
REAL   , DIMENSION(0:POW_MAXPKT), INTENT(INOUT) :: dat   ! Data to be convoluted
REAL                            , INTENT(IN)    :: tthmin ! 2Theta min
REAL                            , INTENT(IN)    :: tthmax ! 2Theta max
REAL                            , INTENT(IN)    :: dtth   ! 2Theta step
REAL                            , INTENT(IN)    :: eta0   ! Lor/Gaus mix constant part
REAL                            , INTENT(IN)    :: eta_l   ! Lor/Gaus mix variable part
REAL                            , INTENT(IN)    :: eta_q   ! Lor/Gaus mix variable part
REAL                            , INTENT(IN)    :: u      ! u*tan^2(Theta)
REAL                            , INTENT(IN)    :: v      ! v*tan  (Theta)
REAL                            , INTENT(IN)    :: w      ! w
REAL                            , INTENT(IN)    :: p1     ! Asymmetry terms p1 to P4
REAL                            , INTENT(IN)    :: p2     ! Asymmetry terms p1 to P4
REAL                            , INTENT(IN)    :: p3     ! Asymmetry terms p1 to P4
REAL                            , INTENT(IN)    :: p4     ! Asymmetry terms p1 to P4
REAL                            , INTENT(IN)    :: pow_width ! Number of FWHM's to calculate
!
REAL(KIND=PREC_DP)    :: fwhm            ! Current FWHM at Theta
REAL, DIMENSION(0:POW_MAXPKT) :: dummy  ! temporary data (0:POW_MAXPKT) 
REAL(KIND=PREC_DP)    :: tth             ! Theta within convolution, main data set
REAL    :: tantth          ! tan(Theta)
REAL(KIND=PREC_DP)    :: tth1            ! Theta values for asymmetry
REAL(KIND=PREC_DP)    :: tth2            ! Theta values for asymmetry
!REAL    :: atheta 
!REAL    :: atwoth 
!REAL    :: fwhm1 
REAL(KIND=PREC_DP)    :: eta             ! actual eta at current 2Theta
REAL    :: pra1, pra2      ! Asymmetry pre factors
INTEGER :: imax, i, j, ii  ! Dummy loop indices
INTEGER :: max_ps 
!                                                                       
!      REAL pseudovoigt 
!      REAL profile_asymmetry 
!      REAL tand
!write(*,*) ' IN CONV_VARIABLE '
pra1 = 1.0D0
pra2 = 1.0D0
!                                                                       
!------ Now convolute                                                   
!                                                                       
imax = INT( (tthmax - tthmin) / dtth )
DO i = 0, imax 
   tth = tthmin + i * dtth 
   tantth = tand (tth * 0.5) 
!     atheta = tth * 0.5 
!     atwoth = tth 
   fwhm = SQRT (MAX (ABS (u * tantth**2 + v * tantth + w), 0.00001) ) 
!      fwhm1 = fwhm 
   max_ps = INT((pow_width * fwhm) / dtth )
   eta = MIN(1.0D0, MAX(0.0D0, eta0 + eta_l * tth + eta_q*tth**2) ) 
   tth1 = 0 * dtth 
   tth2 = 2 * i * dtth 
   pra1 = profile_asymmetry (tth, tth1, fwhm, p1, p2, p3, p4) 
   pra2 = profile_asymmetry (tth, tth2, fwhm, p1, p2, p3, p4) 
!  dummy (i) = dat (i) * (pseudovoigt (tth1, eta, fwhm) * pra1 -     &
!                         pseudovoigt (tth2, eta, fwhm) * pra2)                             
   dummy (i) = dat (i) * (glp_pseud   (tth1, eta, fwhm) * pra1 -     &
                          glp_pseud   (tth2, eta, fwhm) * pra2)                             
!       do j=0,i-1                                                      
   ii = MAX (i - 1 - max_ps + 1, 0) 
   DO j = ii, i - 1 
      tth1 = (i - j) * dtth 
      tth2 = (i + j) * dtth 
      pra1 = profile_asymmetry(tth, tth1, fwhm, p1, p2, p3, p4) 
      pra2 = profile_asymmetry(tth, tth2, fwhm, p1, p2, p3, p4) 
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
      pra1 = profile_asymmetry (tth, - tth1, fwhm, p1, p2, p3, p4) 
      pra2 = profile_asymmetry (tth, - tth2, fwhm, p1, p2, p3, p4) 
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
      REAL dat (0:POW_MAXPKT) 
      REAL tthmin, tthmax, dtth, eta0, eta_l, eta_q
      REAL u, v, w 
      REAL p1, p2, p3, p4 
      REAL pow_width 
      REAL rlambda 
      INTEGER pow_axis 
      INTEGER POW_AXIS_Q 
!                                                                       
      REAL dummy (0:POW_MAXPKT) 
      REAL(KIND=PREC_DP) tth
      REAL tantth 
      REAL(KIND=PREC_DP) tth1 
      REAL(KIND=PREC_DP) tth2 
      REAL atheta 
      REAL atwoth 
      REAL(KIND=PREC_DP) fwhm1 , fwhm
      REAL(KIND=PREC_DP) eta 
      REAL pra1, pra2 
      INTEGER imax, i, j, ii 
      INTEGER max_ps 
!                                                                       
!      REAL pseudovoigt 
!      REAL profile_asymmetry 
!     REAL tand, sind, asind 
!                                                                       
!------ Now convolute                                                   
!                                                                       
      imax = int( (tthmax - tthmin) / dtth )
      DO i = 0, imax 
      tth = tthmin + i * dtth 
      tantth = tand (tth * 0.5) 
      atheta = tth * 0.5 
      atwoth = tth 
      fwhm = sqrt (max (abs (u * tantth**2 + v * tantth + w), 0.00001) ) 
      fwhm1 = fwhm 
      IF (pow_axis.eq.POW_AXIS_Q) THEN 
         atheta = asind (tth * rlambda / REAL(fpi)) 
         tantth = tand (atheta) 
         fwhm1 = sqrt (max (abs (u * tantth**2 + v * tantth + w),       &
         0.00001) )                                                     
         fwhm = 0.500 * (REAL(fpi) * sind (atheta + 0.5 * fwhm1) / rlambda -  &
         REAL(fpi) * sind (atheta - 0.5 * fwhm1) / rlambda)                   
      ENDIF 
      max_ps = int( (pow_width * fwhm) / dtth )
      eta = min (1.0D0, max (0.0D0, eta0 + eta_l * tth + eta_q*tth**2) ) 
      tth1 = 0 * dtth 
      tth2 = 2 * i * dtth 
      pra1 = profile_asymmetry (tth, tth1, fwhm, p1, p2, p3, p4) 
      pra2 = profile_asymmetry (tth, tth2, fwhm, p1, p2, p3, p4) 
!     dummy (i) = dat (i) * (pseudovoigt (tth1, eta, fwhm) * pra1 -     &
!     pseudovoigt (tth2, eta, fwhm) * pra2)                             
!       do j=0,i-1                                                      
      ii = max (i - 1 - max_ps + 1, 0) 
      DO j = ii, i - 1 
      tth1 = (i - j) * dtth 
      tth2 = (i + j) * dtth 
      pra1 = profile_asymmetry (tth, tth1, fwhm, p1, p2, p3, p4) 
      pra2 = profile_asymmetry (tth, tth2, fwhm, p1, p2, p3, p4) 
!     dummy (i) = dummy (i) + dat (j) * (pseudovoigt (tth1, eta, fwhm)  &
!     * pra1 - pseudovoigt (tth2, eta, fwhm) * pra2)                    
      ENDDO 
!       do j=i+1,imax                                                   
      ii = min (i + 1 + max_ps - 1, imax) 
      DO j = i + 1, ii 
      tth1 = (j - i) * dtth 
      tth2 = (j + i) * dtth 
      pra1 = profile_asymmetry (tth, - tth1, fwhm, p1, p2, p3, p4) 
      pra2 = profile_asymmetry (tth, - tth2, fwhm, p1, p2, p3, p4) 
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
REAL   , DIMENSION(0:POW_MAXPKT), INTENT(INOUT) :: dat       ! Data to be convoluted
REAL                            , INTENT(IN)    :: tthmin    ! 2Theta min
REAL                            , INTENT(IN)    :: tthmax    ! 2Theta max
REAL                            , INTENT(IN)    :: dtth      ! 2Theta step
REAL                            , INTENT(IN)    :: sigma2    ! Gaussian Sigma^2
REAL                            , INTENT(IN)    :: corrlin   ! 1/r deppendend width correction
REAL                            , INTENT(IN)    :: corrquad  ! 1/r^2 deppendend width correction
REAL(KIND=PREC_SP)              , INTENT(IN)    :: rcut      ! minimum  distance for clin/(r-rmin)
REAL                            , INTENT(IN)    :: pow_width ! Number of FWHM's to calculate
!
!REAL(KIND=PREC_DP), PARAMETER :: four_ln2  = 2.772588722239781237669D0
REAL(KIND=PREC_DP), PARAMETER :: eightln2  = 2.772588722239781237669D0 * 2.0D0
!
REAL(KIND=PREC_DP)            :: fwhm     ! Current FWHM at Theta
REAL, DIMENSION(0:POW_MAXPKT) :: dummy    ! temporary data (0:POW_MAXPKT) 
REAL(KIND=PREC_DP)            :: tth      ! Theta within convolution, main data set
!REAL                          :: tantth   ! tan(Theta)
REAL(KIND=PREC_DP)    :: sigmasq          ! actual scaled local sigma**2
REAL(KIND=PREC_DP)    :: sigmamin         ! minimum       local sigma**2
REAL(KIND=PREC_DP)    :: eta              ! actual eta at current 2Theta
REAL(KIND=PREC_SP)    :: dist_min  ! minimum  distance for clin/(r-rmin)
INTEGER :: imax, i, j, ii  ! Dummy loop indices
INTEGER :: i1, i2          ! Pseudo Voigt lookup indices
REAL    :: pseudo          ! scale factor for lookup table
INTEGER :: max_ps 
!                                                                       
!------ Now convolute                                                   
!                                                                       
imax = INT( (tthmax - tthmin) / dtth )
sigmasq = sigma2*SQRT(eightln2)
sigmamin = sigmasq * 0.20D0
dist_min = REAL(rcut, KIND=PREC_DP)
!
eta = 0.0     ! Gaussian function
dummy = 0.0   ! dummy(:)
!tth = 1.810
!fwhm = SQRT(MAX(sigmasq - corrlin/(tth-dist_min) - corrquad/(tth-dist_min)**2, sigmamin))
!write(*,7777) ' FWHM ',tth, fwhm , sigmasq, sigmasq - corrlin/(tth-dist_min) - corrquad/(tth-dist_min)**2, sigmamin
!7777 format(a,f5.1, f9.5, f9.5, f12.5, f9.5)
!tth = 2.300
!fwhm = SQRT(MAX(sigmasq - corrlin/(tth-dist_min) - corrquad/(tth-dist_min)**2, sigmamin))
!write(*,7777) ' FWHM ',tth, fwhm , sigmasq, sigmasq - corrlin/(tth-dist_min) - corrquad/(tth-dist_min)**2, sigmamin
!tth = 3.700
!fwhm = SQRT(MAX(sigmasq - corrlin/(tth-dist_min) - corrquad/(tth-dist_min)**2, sigmamin))
!write(*,7777) ' FWHM ',tth, fwhm , sigmasq, sigmasq - corrlin/(tth-dist_min) - corrquad/(tth-dist_min)**2, sigmamin
!tth = 10.000
!fwhm = SQRT(MAX(sigmasq - corrlin/(tth-dist_min) - corrquad/(tth-dist_min)**2, sigmamin))
!write(*,7777) ' FWHM ',tth, fwhm , sigmasq, sigmasq - corrlin/(tth-dist_min) - corrquad/(tth-dist_min)**2, sigmamin
!tth = 2.320
!fwhm = SQRT(MAX(sigmasq - corrlin/tth - corrquad/tth**2, 0.00001))
!write(*,*) ' FWHM ',tth, fwhm , sigmasq, sigmasq - corrlin/tth - corrquad/tth**2
!tth = 4.8000
!fwhm = SQRT(MAX(sigmasq - corrlin/tth - corrquad/tth**2, 0.00001))
!write(*,*) ' FWHM 2.500 ',tth, fwhm 
!tth = tthmax
!fwhm = SQRT(MAX(sigmasq - corrlin/tth - corrquad/tth**2, 0.00001))
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
   ii = MIN(i + 1 + max_ps - 1, imax) 
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
END SUBROUTINE powder_conv_corrlin          
!
!*****7*****************************************************************
!
      REAL function pseudovoigt (dtth, eta, fwhm) 
!-                                                                      
!     calculates the value of a pseudo-voigt function at dtth off the   
!     central position                                                  
!                                                                       
      IMPLICIT none 
      REAL, INTENT(IN) :: dtth 
      REAL, INTENT(IN) :: eta 
      REAL, INTENT(IN) :: fwhm 
!                                                                       
      REAL, PARAMETER :: four_ln2  = 2.772588722
      REAL, PARAMETER :: sq4ln2_pi = 0.939437279
      REAL, PARAMETER :: two_pi    = 0.636619772
      REAL            :: pref_g 
      REAL            :: pref_l 
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
REAL function profile_asymmetry (tth, dtth, fwhm, p1, p2, p3, p4) 
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
REAL, INTENT(IN) :: p1, p2, p3, p4 
!                                                                       
REAL :: zz 
REAL :: fa, fb 
!                                                                       
!     REAL :: tand 
!                                                                       
zz = dtth / fwhm 
!                                                                       
fa = 2. * zz * exp ( - zz**2) 
fb = 2. * (2. * zz**2 - 3.) * fa 
!                                                                       
profile_asymmetry = 1.0 + (p1 * fa + p2 * fb) / tand(0.5 * tth)  &
                        + (p3 * fa + p4 * fb) / tand(tth)                                
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
pow_u2aver       = 0.0
pow_nreal        = 0
!
!     Prepare and calculate average atom numbers
!
ALLOCATE(natom(0:cr_nscat))
natom = 0
DO i=1,cr_natoms
   natom(cr_iscat(i)) = natom(cr_iscat(i)) + 1
ENDDO
pow_nreal = SUM(natom(1:))  ! Add real atom numbers 
!write(*,*) ' IN POW_f2aver nreal, ', pow_nreal
!
DO iscat = 1, cr_nscat
   signum = 1.0D0
   IF(REAL(cfact_pure(1, iscat))< 0.0D0) signum = -1.0D0
   DO i = 1, num1
      pow_f2aver (i) = pow_f2aver (i)  + &
                 DBLE (       cfact_pure(powder_istl(i), iscat)  * &
                       conjg (cfact_pure(powder_istl(i), iscat)))  &
               * natom (iscat)/pow_nreal
      pow_faver2 (i) = pow_faver2 (i) +  &
            SQRT(DBLE (       cfact_pure(powder_istl(i), iscat)  * &
                       conjg (cfact_pure(powder_istl(i), iscat)))) &
               * natom (iscat)/pow_nreal                           &
               * signum
   ENDDO
   pow_u2aver = pow_u2aver + (cr_dw(iscat)) * natom (iscat)/pow_nreal
ENDDO
pow_faver2(:) = pow_faver2(:)**2
pow_u2aver    = pow_u2aver /8./REAL(pi)**2
DEALLOCATE(natom)
!
!
END SUBROUTINE powder_f2aver
!
!*******************************************************************************
!
END MODULE powder_write_mod
