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
      SUBROUTINE powder_out 
!-                                                                      
!     Write the powder pattern                                          
!+                                                                      
      USE discus_config_mod 
      USE debye_mod 
      USE diffuse_mod 
      USE output_mod 
      USE powder_mod 
      USE wink_mod
      IMPLICIT none 
!                                                                       
      INTEGER iff 
      PARAMETER (iff = 2) 
!                                                                       
!                                                                       
      INTEGER ii, j , iii
      INTEGER   :: all_status  ! Allocation status
      INTEGER   :: npkt        ! number of points in powder pattern
      INTEGER   :: npkt_equi   ! number of points in equidistant powder pattern
      LOGICAL lread 
      REAL, DIMENSION(:), ALLOCATABLE :: pow_tmp  ! Local temporary copy of intensities
      REAL, DIMENSION(:), ALLOCATABLE :: xpl  ! x-values of calculated powder pattern
      REAL, DIMENSION(:), ALLOCATABLE :: ypl  ! y-values of calculated powder pattern
      REAL, DIMENSION(:), ALLOCATABLE :: y2a  ! y-values of splined    powder pattern
      REAL :: ttheta, lp=1.0
      REAL ss, st 
      REAL :: q=0.0, stl=0.0, dstar=0.0
      REAL xmin, xmax, xdel , xpos
      REAL      :: xequ    ! x-position of equdistant curve
      REAL      :: yequ    ! y-value    of equdistant curve
      REAL      :: tthmin  ! minimum for equdistant curve
      REAL      :: tthmax  ! minimum for equdistant curve
      REAL      ::   qmin  ! minimum for equdistant curve
      REAL      ::   qmax  ! minimum for equdistant curve
!                                                                       
!      REAL lorentz 
!      REAL polarisation 
      REAL sind, asind 
!
      ALLOCATE(pow_tmp(0:POW_MAXPKT),stat = all_status)  ! Allocate array for powder pattern copy
      ALLOCATE(xpl(0:POW_MAXPKT),stat = all_status)  ! Allocate array for calculated powder pattern
      ALLOCATE(ypl(0:POW_MAXPKT),stat = all_status)  ! Allocate array for calculated powder pattern
      pow_tmp = 0.0
      xpl     = 0.0
      ypl     = 0.0
!                                                                       
      IF (pow_four_type.eq.POW_COMPL.or.pow_four_type.eq.POW_NEW) then 
         IF (pow_axis.eq.POW_AXIS_Q) then 
            xmin = pow_qmin 
            xmax = pow_qmax 
            xdel = pow_deltaq 
         ELSEIF (pow_axis.eq.POW_AXIS_TTH) then 
            xmin = pow_tthmin 
            xmax = pow_tthmax 
            xdel = pow_deltatth 
         ELSE 
            ier_num = - 104 
            ier_typ = ER_APPL 
            ier_msg (1) = 'Use command ==> set axis,{"tth"|"q"}' 
            ier_msg (2) = 'within the powder menu to define the axis' 
            ier_msg (3) = ' ' 
            RETURN 
         ENDIF 
         npkt = NINT((xmax-xmin)/xdel) + 1
      ELSEIF (pow_four_type.eq.POW_HIST ) THEN
         IF (pow_axis.eq.POW_AXIS_Q) then 
            xmin = pow_qmin 
            xmax = pow_qmax 
            xdel = (pow_qmax - pow_qmin) / (num (1) ) 
         ELSEIF (pow_axis.eq.POW_AXIS_TTH) then 
            xmin = pow_tthmin 
            xmax = pow_tthmax 
            xdel = (pow_tthmax - pow_tthmin) / (num (1) ) 
         ELSE 
            ier_num = - 104 
            ier_typ = ER_APPL 
            ier_msg (1) = 'Use command ==> set axis,{"tth"|"q"}' 
            ier_msg (2) = 'within the powder menu to define the axis' 
            ier_msg (3) = ' ' 
            RETURN 
         ENDIF 
         npkt = num(1)
      ENDIF 
      lread = .false. 
      IF (ier_num.eq.0) then 
         IF (pow_four_type.ne.POW_COMPL) then 
!                                                                       
!     This is a Debye calculation, copy rsf or csf into pow_tmp         
!                                                                       
            IF (pow_four_type.eq.POW_DEBYE) then 
               IF (num (1) .lt.POW_MAXPKT) then 
                  DO j = 1, num (1) 
                  pow_tmp (j) = real (csf (j) ) 
                  ENDDO 
               ENDIF 
            ELSEIF (                                                    &
            pow_four_type.eq.POW_FAST.or.pow_four_type.eq.POW_HIST)     &
            then                                                        
               IF (num (1) .le.POW_MAXPKT) then 
                  DO j = 1, num (1) 
                  pow_tmp (j) = rsf (j) 
                  ENDDO 
               ENDIF 
            ENDIF 
         ELSE
            pow_tmp(:) = pow_qsp(:)
         ENDIF 
!                                                                       
!- -Does the powder pattern have to be convoluted by a profile function?
!                                                                       
         IF (pow_profile.eq.POW_PROFILE_GAUSS) then 
            IF (pow_delta.gt.0.0) then 
               CALL powder_conv_res (pow_tmp, xmin, xmax, xdel,         &
               pow_delta, POW_MAXPKT)                                    
            ENDIF 
         ELSEIF (pow_profile.eq.POW_PROFILE_PSVGT) then 
           IF (pow_u.ne.0.0.or.pow_v.ne.0.0.or.pow_etax.ne.0.0.or.      &
               pow_p1.ne.0.0.or.pow_p2.ne.0.0.or.pow_p3.ne.0.0.or.      &
               pow_p4.ne.0.0                                      ) THEN       
!DBG                                          .or.                      
!DBG     &                pow_axis.eq.POW_AXIS_Q                 ) then 
               CALL powder_conv_psvgt_uvw (pow_tmp, xmin, xmax, xdel,   &
               pow_eta, pow_etax, pow_u, pow_v, pow_w, pow_p1, pow_p2,  &
               pow_p3, pow_p4, pow_width, POW_MAXPKT)
            ELSE 
               CALL powder_conv_psvgt_fix (pow_tmp, xmin, xmax, xdel,   &
               pow_eta, pow_w, pow_width, POW_MAXPKT)
            ENDIF 
         ENDIF 
!                                                                       
!------ copy the powder pattern into output array, if necessary this will be put on
!       equidistant scale
!                                                                       
         IF (pow_four_type.eq.POW_COMPL.or.pow_four_type.eq.POW_NEW) THEN                                                           
            DO ii = 0, npkt - 1
               iii = ii + 1
               xpos = ii * xdel + xmin 
               IF (pow_axis.eq.POW_AXIS_Q) then 
                  q      = xpos
                  dstar  = q / zpi
                  stl    = q / zpi / 2.
                  ttheta = 2.*asind ( q / 2. /zpi *rlambda )
                  lp     = lorentz (ttheta) * polarisation (ttheta) 
               ELSEIF (pow_axis.eq.POW_AXIS_TTH) then 
                  ttheta = xpos
                  stl    =            sind (ttheta * 0.5) / rlambda 
                  dstar  = 2. *       sind (ttheta * 0.5) / rlambda 
                  q      = 2. * zpi * sind (ttheta * 0.5) / rlambda 
                  lp     = lorentz (ttheta) * polarisation (ttheta) 
               ENDIF 
               IF (cpow_form.eq.'tth') then 
                  xpl(iii) = ttheta
               ELSEIF (cpow_form.eq.'stl') then 
                  xpl(iii) = stl
               ELSEIF (cpow_form.eq.'q  ') then 
                  xpl(iii) = q
               ELSEIF (cpow_form.eq.'dst') then 
                  xpl(iii) = dstar
               ELSEIF (cpow_form.eq.'lop') then 
                  xpl(iii) = ttheta
               ENDIF 
               ypl(iii) = pow_tmp(ii) * lp
!              write(iff,*), xpl(iii),ypl(iii)
            ENDDO 
         ELSEIF (pow_four_type.eq.POW_HIST) then 
            IF (pow_axis.eq.POW_AXIS_DSTAR) then 
            ELSEIF (pow_axis.eq.POW_AXIS_Q) then 
               xm(1)  = pow_qmin / zpi 
               ss     = pow_qmax / zpi 
               st     = (pow_qmax - pow_deltaq) / zpi 
               uin(1) = pow_deltaq / zpi 
            ELSEIF (pow_axis.eq.POW_AXIS_TTH) then 
               xm(1)  = 2 * sind (0.5 * pow_tthmin) / rlambda 
               ss     = 2 * sind (0.5 * pow_tthmax) / rlambda 
               st     = 2 * sind (0.5 * (pow_tthmax - pow_deltatth) ) / rlambda
               uin(1) = (ss - st) / 2. 
            ENDIF 
            DO ii = 1, num (1) 
            dstar = (xm (1) + (ii - 1) * uin (1) ) 
            stl = .5 * (xm (1) + (ii - 1) * uin (1) ) 
            q = zpi * (xm (1) + (ii - 1) * uin (1) ) 
            ttheta = 2. * asind (dstar * rlambda / 2.) 
!DBG          lp     = lorentz(ttheta)*polarisation(ttheta)             
            lp = polarisation (ttheta) 
            IF (cpow_form.eq.'tth') then 
               xpl(ii) = ttheta
            ELSEIF (cpow_form.eq.'stl') then 
               xpl(ii) = stl
            ELSEIF (cpow_form.eq.'q  ') then 
               xpl(ii) = q
            ELSEIF (cpow_form.eq.'dst') then 
               xpl(ii) = dstar
            ENDIF 
               ypl(ii) = pow_tmp(ii) * lp
            ENDDO 
         ENDIF 
      ENDIF 
!
      CALL oeffne (iff, outfile, 'unknown') 
      IF( cpow_form == 'tth' ) THEN
         IF ( pow_axis      == POW_AXIS_Q  .or.  &        ! Non matching form, spline onto equidistant steps
              pow_four_type == POW_HIST            ) THEN ! DEBYE, always spline
            IF(pow_tthmin < xpl(1) ) THEN                 ! User lower limit too low!
               tthmin = pow_tthmin + (INT( (xpl(1)-pow_tthmin   )/pow_deltatth) + 1)*pow_deltatth
            ELSE
               tthmin = pow_tthmin
            ENDIF
            IF(pow_tthmax > xpl(npkt) ) THEN              ! User upper limit too high!
               tthmax = pow_tthmax - (INT( (pow_tthmax-xpl(npkt))/pow_deltatth) + 1)*pow_deltatth
            ELSE
               tthmax = pow_tthmax
            ENDIF
            ALLOCATE(y2a(1:POW_MAXPKT),stat = all_status)  ! Allocate array for calculated powder pattern
            y2a = 0.0
            CALL spline (npkt, xpl, ypl, 1e31, 1e31, y2a)
            npkt_equi = INT((tthmax-tthmin)/pow_deltatth) + 1
            DO ii = 1, npkt_equi
               xequ = tthmin + (ii-1)*pow_deltatth
               CALL splint (npkt, xpl, ypl, y2a, xequ, yequ)
               IF(ier_num/=0) THEN
                  CLOSE(iff)
                  DEALLOCATE( pow_tmp, stat = all_status)
                  DEALLOCATE( xpl, stat = all_status)
                  DEALLOCATE( ypl, stat = all_status)
                  DEALLOCATE( y2a, stat = all_status)
                  RETURN
               ENDIF
               WRITE( iff, *) xequ, yequ
            ENDDO
            DEALLOCATE(y2a, stat = all_status)
         ELSE
            DO ii = 1,npkt
               WRITE( iff, *) xpl(ii),ypl(ii)
            ENDDO
         ENDIF
      ELSEIF( cpow_form == 'q' ) THEN                       ! axis is Q
         IF ( pow_axis      == POW_AXIS_TTH  .or.  &        ! Non matching form, spline onto equidistant steps
              pow_four_type == POW_HIST              ) THEN ! DEBYE, always spline
            IF(pow_qmin < xpl(1) ) THEN                     ! User lower limit too low!
               qmin = pow_qmin + (INT( (xpl(1)-pow_qmin   )/pow_deltaq) + 1)*pow_deltaq
            ELSE
               qmin = pow_qmin
            ENDIF
            IF(pow_qmax > xpl(npkt) ) THEN                  ! User upper limit too high!
               qmax = pow_qmax - (INT( (pow_qmax-xpl(npkt))/pow_deltaq) + 1)*pow_deltaq
            ELSE
               qmax = pow_qmax
            ENDIF
            ALLOCATE(y2a(1:POW_MAXPKT),stat = all_status)  ! Allocate array for calculated powder pattern
            CALL spline (npkt, xpl, ypl, 1e31, 1e31, y2a)
            npkt_equi = NINT((qmax-qmin)/pow_deltaq) + 1
            DO ii = 1, npkt_equi
               xequ = qmin + (ii-1)*pow_deltaq
               CALL splint (npkt, xpl, ypl, y2a, xequ, yequ)
               IF(ier_num/=0) THEN
                  CLOSE(iff)
                  DEALLOCATE( pow_tmp, stat = all_status)
                  DEALLOCATE( xpl, stat = all_status)
                  DEALLOCATE( ypl, stat = all_status)
                  DEALLOCATE( y2a, stat = all_status)
                  RETURN
               ENDIF
               WRITE( iff, *) xequ, yequ
            ENDDO
            DEALLOCATE(y2a, stat = all_status)
         ELSE
            DO ii = 1,npkt
               WRITE( iff, *) xpl(ii),ypl(ii)
            ENDDO
         ENDIF
      ELSE
         DO ii = 1,npkt
            WRITE( iff, *) xpl(ii),ypl(ii)
         ENDDO
      ENDIF
!
      CLOSE(iff)
!
      DEALLOCATE( pow_tmp, stat = all_status)
      DEALLOCATE( xpl, stat = all_status)
      DEALLOCATE( ypl, stat = all_status)
!                                                                       
      END SUBROUTINE powder_out                     
!*****7*****************************************************************
      REAL function lorentz (ttheta) 
!+                                                                      
!-                                                                      
      USE discus_config_mod 
      USE powder_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      REAL ttheta 
!                                                                       
!
      REAL sind
!                                                                       
      lorentz = 1.0
      
      IF (pow_four_type.eq.POW_DEBYE) then 
         lorentz = 1.0 
      ELSE 
         IF (pow_lp.eq.POW_LP_BRAGG) then 
            lorentz = 0.5 / sind (0.5 * ttheta) / sind (ttheta) 
         ELSEIF (pow_lp.eq.POW_LP_NEUT) then 
            lorentz = 0.5 / sind (0.5 * ttheta) / sind (ttheta) 
         ELSEIF (pow_lp.eq.POW_LP_NONE) then 
            lorentz = 1.0 
         ELSEIF (pow_lp.eq.POW_LP_SYNC) then 
            lorentz = 0.5 / sind (0.5 * ttheta) / sind (ttheta) 
         ENDIF 
      ENDIF 
!                                                                       
      END FUNCTION lorentz                          
!*****7*****************************************************************
      REAL function polarisation (ttheta) 
!+                                                                      
!-                                                                      
      USE discus_config_mod 
      USE powder_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      REAL ttheta 
!
!
      REAL cosd 
!                                                                       
      polarisation = 1.0
      
      IF (pow_lp.eq.POW_LP_BRAGG) then 
         polarisation = (1. + (cosd (ttheta) ) **2 * pow_lp_fac)        &
         / (1. + pow_lp_fac)                                            
      ELSEIF (pow_lp.eq.POW_LP_NEUT) then 
         polarisation = 1.0 
      ELSEIF (pow_lp.eq.POW_LP_NONE) then 
         polarisation = 1.0 
      ELSEIF (pow_lp.eq.POW_LP_SYNC) then 
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
      IMPLICIT none 
!                                                                       
!                                                                       
      REAL ttheta 
!                                                                       
      REAL sind, cosd 
!                                                                       
      IF (pow_four_type.eq.POW_DEBYE) then 
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
      gauss (i) = 1.0 / sqrt (pi) / delta * exp ( - (tth**2 / delta**2) &
      )                                                                 
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
      dat (i) = dummy (i) * dtth 
      ENDDO 
!                                                                       
      END SUBROUTINE powder_conv_res                
!*****7*****************************************************************
SUBROUTINE powder_conv_psvgt_fix (dat, tthmin, tthmax, dtth, eta, &
      w, pow_width, POW_MAXPKT)
!-                                                                      
!     Convolute powder pattern with resolution function (Pseudo-Voigt)  
!     Constant FWHM, Constant eta                                       
!+                                                                      
USE discus_config_mod 
USE wink_mod
IMPLICIT none 
!                                                                       
!
INTEGER, INTENT(IN) :: POW_MAXPKT
!                                                                       
REAL dat (0:POW_MAXPKT) 
REAL tthmin, tthmax, dtth, fwhm, eta
REAL w 
REAL pow_width 
!                                                                       
REAL dummy (0:POW_MAXPKT) 
REAL psvgt (0:2 * POW_MAXPKT) 
REAL tth
INTEGER imax, i, j, ii 
INTEGER max_ps 
!                                                                       
!REAL pseudovoigt 
!                                                                       
!------ Setup Pseudo-Voigt                                              
!                                                                       
fwhm = sqrt (abs (w) ) 
max_ps = int( (pow_width * fwhm) / dtth )
psvgt = 0.0
DO i = 0, max_ps 
   tth = i * dtth 
   psvgt (i) = pseudovoigt (tth, eta, fwhm) 
ENDDO 
!                                                                       
!     DO i = max_ps + 1, 2 * POW_MAXPKT 
!     psvgt (i) = 0.0 
!     ENDDO 
!                                                                       
!------ Now convolute                                                   
!                                                                       
      imax = int( (tthmax - tthmin) / dtth )
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
!      IF (i + ii.le.imax) then 
!         WRITE ( * , * ) ' i,j, psvgt(j-i), psvgt(j+i)', i, ii, psvgt ( &
!         ii - i) , psvgt (ii + i)                                       
!      ENDIF 
ENDDO 
!                                                                       
DO i = 0, imax 
   dat (i) = dummy (i) * dtth 
ENDDO 
!                                                                       
END SUBROUTINE powder_conv_psvgt_fix          
!*****7*****************************************************************
      SUBROUTINE powder_conv_psvgt_uvw (dat, tthmin, tthmax, dtth, eta0,&
      etax, u, v, w, p1, p2, p3, p4, pow_width, POW_MAXPKT)
!-                                                                      
!     Convolute powder pattern with resolution function (Pseudo-Voigt)  
!     FWHM according to caglioti equation, Constant eta                 
!     FWHM = sqrt ( U*tan**2(Theta) + V*tan(Theta) + W)                 
!+                                                                      
      USE discus_config_mod 
      USE wink_mod
      IMPLICIT none 
!                                                                       
!
      INTEGER, INTENT(IN) :: POW_MAXPKT
!                                                                       
      REAL dat (0:POW_MAXPKT) 
      REAL tthmin, tthmax, dtth, fwhm, eta0, etax 
      REAL u, v, w 
      REAL p1, p2, p3, p4 
      REAL pow_width 
!                                                                       
      REAL dummy (0:POW_MAXPKT) 
      REAL tth
      REAL tantth 
      REAL tth1 
      REAL tth2 
      REAL atheta 
      REAL atwoth 
      REAL fwhm1 
      REAL eta 
      REAL pra1, pra2 
      INTEGER imax, i, j, ii 
      INTEGER max_ps 
!                                                                       
!      REAL pseudovoigt 
!      REAL profile_asymmetry 
      REAL tand
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
      max_ps = int( (pow_width * fwhm) / dtth )
      eta = min (1.0, max (0.0, eta0 + etax * tth) ) 
      tth1 = 0 * dtth 
      tth2 = 2 * i * dtth 
      pra1 = profile_asymmetry (tth, tth1, fwhm, p1, p2, p3, p4) 
      pra2 = profile_asymmetry (tth, tth2, fwhm, p1, p2, p3, p4) 
      dummy (i) = dat (i) * (pseudovoigt (tth1, eta, fwhm) * pra1 -     &
      pseudovoigt (tth2, eta, fwhm) * pra2)                             
!       do j=0,i-1                                                      
      ii = max (i - 1 - max_ps + 1, 0) 
      DO j = ii, i - 1 
      tth1 = (i - j) * dtth 
      tth2 = (i + j) * dtth 
      pra1 = profile_asymmetry (tth, tth1, fwhm, p1, p2, p3, p4) 
      pra2 = profile_asymmetry (tth, tth2, fwhm, p1, p2, p3, p4) 
      dummy (i) = dummy (i) + dat (j) * (pseudovoigt (tth1, eta, fwhm)  &
      * pra1 - pseudovoigt (tth2, eta, fwhm) * pra2)                    
      ENDDO 
!       do j=i+1,imax                                                   
      ii = min (i + 1 + max_ps - 1, imax) 
      DO j = i + 1, ii 
      tth1 = (j - i) * dtth 
      tth2 = (j + i) * dtth 
      pra1 = profile_asymmetry (tth, - tth1, fwhm, p1, p2, p3, p4) 
      pra2 = profile_asymmetry (tth, - tth2, fwhm, p1, p2, p3, p4) 
      dummy (i) = dummy (i) + dat (j) * (pseudovoigt (tth1, eta, fwhm)  &
      * pra1 - pseudovoigt (tth2, eta, fwhm) * pra2)                    
      ENDDO 
      ENDDO 
!                                                                       
      DO i = 0, imax 
      dat (i) = dummy (i) * dtth 
      ENDDO 
!                                                                       
      END SUBROUTINE powder_conv_psvgt_uvw          
!*****7*****************************************************************
      SUBROUTINE powder_conv_psvgt_uvw_Qscale (dat, tthmin, tthmax,     &
      dtth, eta0, etax, u, v, w, p1, p2, p3, p4, pow_width, rlambda,    &
      pow_axis, POW_AXIS_Q, POW_MAXPKT)                                             
!-                                                                      
!     Convolute powder pattern with resolution function (Pseudo-Voigt)  
!     FWHM according to caglioti equation, Constant eta                 
!     FWHM = sqrt ( U*tan**2(Theta) + V*tan(Theta) + W)                 
!+                                                                      
      USE discus_config_mod 
      USE wink_mod
      IMPLICIT none 
!                                                                       
!
      INTEGER, INTENT(IN) :: POW_MAXPKT
!                                                                       
      REAL dat (0:POW_MAXPKT) 
      REAL tthmin, tthmax, dtth, fwhm, eta0, etax 
      REAL u, v, w 
      REAL p1, p2, p3, p4 
      REAL pow_width 
      REAL rlambda 
      INTEGER pow_axis 
      INTEGER POW_AXIS_Q 
!                                                                       
      REAL dummy (0:POW_MAXPKT) 
      REAL tth
      REAL tantth 
      REAL tth1 
      REAL tth2 
      REAL atheta 
      REAL atwoth 
      REAL fwhm1 
      REAL eta 
      REAL pra1, pra2 
      INTEGER imax, i, j, ii 
      INTEGER max_ps 
!                                                                       
!      REAL pseudovoigt 
!      REAL profile_asymmetry 
      REAL tand, sind, asind 
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
      IF (pow_axis.eq.POW_AXIS_Q) then 
         atheta = asind (tth * rlambda / fpi) 
         tantth = tand (atheta) 
         fwhm1 = sqrt (max (abs (u * tantth**2 + v * tantth + w),       &
         0.00001) )                                                     
         fwhm = 0.500 * (fpi * sind (atheta + 0.5 * fwhm1) / rlambda -  &
         fpi * sind (atheta - 0.5 * fwhm1) / rlambda)                   
      ENDIF 
      max_ps = int( (pow_width * fwhm) / dtth )
      eta = min (1.0, max (0.0, eta0 + etax * tth) ) 
      tth1 = 0 * dtth 
      tth2 = 2 * i * dtth 
      pra1 = profile_asymmetry (tth, tth1, fwhm, p1, p2, p3, p4) 
      pra2 = profile_asymmetry (tth, tth2, fwhm, p1, p2, p3, p4) 
      dummy (i) = dat (i) * (pseudovoigt (tth1, eta, fwhm) * pra1 -     &
      pseudovoigt (tth2, eta, fwhm) * pra2)                             
!       do j=0,i-1                                                      
      ii = max (i - 1 - max_ps + 1, 0) 
      DO j = ii, i - 1 
      tth1 = (i - j) * dtth 
      tth2 = (i + j) * dtth 
      pra1 = profile_asymmetry (tth, tth1, fwhm, p1, p2, p3, p4) 
      pra2 = profile_asymmetry (tth, tth2, fwhm, p1, p2, p3, p4) 
      dummy (i) = dummy (i) + dat (j) * (pseudovoigt (tth1, eta, fwhm)  &
      * pra1 - pseudovoigt (tth2, eta, fwhm) * pra2)                    
      ENDDO 
!       do j=i+1,imax                                                   
      ii = min (i + 1 + max_ps - 1, imax) 
      DO j = i + 1, ii 
      tth1 = (j - i) * dtth 
      tth2 = (j + i) * dtth 
      pra1 = profile_asymmetry (tth, - tth1, fwhm, p1, p2, p3, p4) 
      pra2 = profile_asymmetry (tth, - tth2, fwhm, p1, p2, p3, p4) 
      dummy (i) = dummy (i) + dat (j) * (pseudovoigt (tth1, eta, fwhm)  &
      * pra1 - pseudovoigt (tth2, eta, fwhm) * pra2)                    
      ENDDO 
      ENDDO 
!                                                                       
      DO i = 0, imax 
      dat (i) = dummy (i) * dtth 
      ENDDO 
!                                                                       
      END SUBROUTINE powder_conv_psvgt_uvw_Qscale   
!*****7*****************************************************************
      REAL function pseudovoigt (dtth, eta, fwhm) 
!-                                                                      
!     calculates the value of a pseudo-voigt function at dtth off the   
!     central position                                                  
!                                                                       
      IMPLICIT none 
      REAL dtth 
      REAL eta 
      REAL fwhm 
!                                                                       
      REAL pi 
      REAL four_ln2 
      REAL sq4ln2_pi 
      REAL two_pi 
      REAL pref_g 
      REAL pref_l 
!                                                                       
      DATA pi / 3.141592654 / 
      DATA four_ln2 / 2.772588722 / 
      DATA sq4ln2_pi / 0.939437279 / 
      DATA two_pi / 0.636619772 / 
!                                                                       
      pref_g = sq4ln2_pi / fwhm 
      pref_l = two_pi * fwhm 
!                                                                       
      pseudovoigt = eta * pref_l / (fwhm**2 + 4. * dtth**2) + (1. - eta)&
      * pref_g * exp ( - four_ln2 * (dtth / fwhm) **2)                  
!                                                                       
      END FUNCTION pseudovoigt                      
!*****7*****************************************************************
      REAL function profile_asymmetry (tth, dtth, fwhm, p1, p2, p3, p4) 
!-                                                                      
!     calculates the asymmetry parameter for the profile function       
!                                                                       
      IMPLICIT none 
!                                                                       
      REAL tth 
      REAL dtth 
      REAL fwhm 
      REAL p1, p2, p3, p4 
!                                                                       
      REAL zz 
      REAL fa, fb 
!                                                                       
      REAL tand 
!                                                                       
      zz = dtth / fwhm 
!                                                                       
      fa = 2. * zz * exp ( - zz**2) 
      fb = 2. * (2. * zz**2 - 3.) * fa 
!                                                                       
      profile_asymmetry = 1.0 + (p1 * fa + p2 * fb) / tand (0.5 * tth)  &
      + (p3 * fa + p4 * fb) / tand (tth)                                
!                                                                       
      END FUNCTION profile_asymmetry                
!*****7*****************************************************************
!                                                                       
      SUBROUTINE spline (n, x, y, yp1, ypn, y2) 
!
!      PARAMETER (nmax = maxarray) 
!
      INTEGER,              INTENT(IN)  :: n
      REAL, DIMENSION(1:n), INTENT(IN)  :: x
      REAL, DIMENSION(1:n), INTENT(IN)  :: y
      REAL                , INTENT(IN)  :: yp1
      REAL                , INTENT(IN)  :: ypn
      REAL, DIMENSION(1:n), INTENT(OUT) :: y2
!
      INTEGER               :: i,k
      REAL, DIMENSION(1:n)  :: u
      REAL                  :: p, qn, sig, un
!
!     INTEGER  :: klo, khi, k
!     REAL     :: a,b
!     DIMENSION x (n), y (n), y2 (n), u (n) 
      IF (yp1.gt..99e30) then 
         y2 (1) = 0. 
         u (1) = 0. 
      ELSE 
         y2 (1) = - 0.5 
         u (1) = (3. / (x (2) - x (1) ) ) * ( (y (2) - y (1) ) /        &
         (x (2) - x (1) ) - yp1)                                        
      ENDIF 
      DO 11 i = 2, n - 1 
         sig = (x (i) - x (i - 1) ) / (x (i + 1) - x (i - 1) ) 
         p = sig * y2 (i - 1) + 2. 
         y2 (i) = (sig - 1.) / p 
         u (i) = (6. * ( (y (i + 1) - y (i) ) / (x (i + 1) - x (i) )    &
         - (y (i) - y (i - 1) ) / (x (i) - x (i - 1) ) ) / (x (i + 1)   &
         - x (i - 1) ) - sig * u (i - 1) ) / p                          
   11 END DO 
      IF (ypn.gt..99e30) then 
         qn = 0. 
         un = 0. 
      ELSE 
         qn = 0.5 
         un = (3. / (x (n) - x (n - 1) ) ) * (ypn - (y (n) - y (n - 1) )&
         / (x (n) - x (n - 1) ) )                                       
      ENDIF 
      y2 (n) = (un - qn * u (n - 1) ) / (qn * y2 (n - 1) + 1.) 
      DO 12 k = n - 1, 1, - 1 
         y2 (k) = y2 (k) * y2 (k + 1) + u (k) 
   12 END DO 
      RETURN 
      END SUBROUTINE spline                         
!                                                                       
      SUBROUTINE splint (n, xa, ya, y2a, x, y) 
!
      INTEGER,              INTENT(IN)  :: n
      REAL, DIMENSION(1:n), INTENT(IN)  :: xa
      REAL, DIMENSION(1:n), INTENT(IN)  :: ya
      REAL, DIMENSION(1:n), INTENT(IN)  :: y2a
      REAL                , INTENT(IN)  :: x
      REAL                , INTENT(OUT) :: y
      INTEGER  :: klo, khi, k
      REAL     :: a,b,h
!
      klo = 1 
      khi = n 
    1 IF (khi - klo.gt.1) then 
         k = (khi + klo) / 2 
         IF (xa (k) .gt.x) then 
            khi = k 
         ELSE 
            klo = k 
         ENDIF 
         GOTO 1 
      ENDIF 
      h = xa (khi) - xa (klo) 
      IF (h.eq.0.) THEN
         ier_num = -121
         ier_typ = ER_APPL
         WRITE(ier_msg(1),'(''x- pos: '',F10.4,2x, F10.4)') xa(khi), xa(klo)
         WRITE(ier_msg(2),'(''x- pos: '',I6   ,6x, I6  )' )    khi ,    klo 
         RETURN
      ENDIF
      a = (xa (khi) - x) / h 
      b = (x - xa (klo) ) / h 
      y = a * ya (klo) + b * ya (khi) + ( (a**3 - a) * y2a (klo)        &
      + (b**3 - b) * y2a (khi) ) * (h**2) / 6.                          
      RETURN 
      END SUBROUTINE splint                         
END MODULE powder_write_mod
