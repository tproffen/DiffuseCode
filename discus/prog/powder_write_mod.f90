MODULE powder_write_mod
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
      USE config_mod 
      USE debye_mod 
      USE diffuse_mod 
      USE output_mod 
      USE powder_mod 
      IMPLICIT none 
!                                                                       
      INTEGER iff 
      PARAMETER (iff = 2) 
!                                                                       
      include'errlist.inc' 
      include'wink.inc' 
!                                                                       
      INTEGER ii, j 
      LOGICAL lread 
      REAL ttheta, lp 
      REAL ss, st 
      REAL dst, q, stl, dstar 
      REAL value 
      REAL xmin, xmax, xdel 
!                                                                       
!      REAL lorentz 
!      REAL polarisation 
      REAL cosd, sind, asind 
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
      ELSEIF (pow_four_type.eq.POW_HIST.or.pow_four_type.eq.POW_FAST.or.&
     &pow_four_type.eq.POW_DEBYE) then                                  
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
      ENDIF 
      lread = .false. 
      CALL oeffne (iff, outfile, 'unknown', lread) 
      IF (ier_num.eq.0) then 
         IF (pow_four_type.ne.POW_COMPL) then 
!                                                                       
!     This is a Debye calculation, copy rsf or csf into pow_qsp         
!                                                                       
            IF (pow_four_type.eq.POW_DEBYE) then 
               IF (num (1) .lt.POW_MAXPKT) then 
                  DO j = 1, num (1) 
                  pow_qsp (j) = real (csf (j) ) 
                  ENDDO 
               ENDIF 
            ELSEIF (                                                    &
            pow_four_type.eq.POW_FAST.or.pow_four_type.eq.POW_HIST)     &
            then                                                        
               IF (num (1) .le.POW_MAXPKT) then 
                  DO j = 1, num (1) 
                  pow_qsp (j) = rsf (j) 
                  ENDDO 
               ENDIF 
            ENDIF 
         ENDIF 
!                                                                       
!- -Does the powder pattern have to be convoluted by a profile function?
!                                                                       
         IF (pow_profile.eq.POW_PROFILE_GAUSS) then 
            IF (pow_delta.gt.0.0) then 
               CALL powder_conv_res (pow_qsp, xmin, xmax, xdel,         &
               pow_delta, pow_width, POW_MAXPKT)                                    
            ENDIF 
         ELSEIF (pow_profile.eq.POW_PROFILE_PSVGT) then 
           IF (pow_u.ne.0.0.or.pow_v.ne.0.0.or.pow_etax.ne.0.0.or.      &
               pow_p1.ne.0.0.or.pow_p2.ne.0.0.or.pow_p3.ne.0.0.or.      &
               pow_p4.ne.0.0                                      ) THEN       
!DBG                                          .or.                      
!DBG     &                pow_axis.eq.POW_AXIS_Q                 ) then 
               CALL powder_conv_psvgt_uvw (pow_qsp, xmin, xmax, xdel,   &
               pow_eta, pow_etax, pow_u, pow_v, pow_w, pow_p1, pow_p2,  &
               pow_p3, pow_p4, pow_width, rlambda, pow_axis, POW_AXIS_Q,&
               POW_MAXPKT)
            ELSE 
               CALL powder_conv_psvgt_fix (pow_qsp, xmin, xmax, xdel,   &
               pow_eta, pow_etax, pow_u, pow_v, pow_w, pow_p1, pow_p2,  &
               pow_p3, pow_p4, pow_width, POW_MAXPKT)                               
            ENDIF 
         ENDIF 
!                                                                       
!------ write the powder pattern                                        
!                                                                       
         IF (pow_four_type.eq.POW_COMPL.or.pow_four_type.eq.POW_NEW)    &
         then                                                           
            j = (pow_tthmax - pow_tthmin) / pow_deltatth 
            DO ii = 0, j 
            ttheta = ii * pow_deltatth + pow_tthmin 
            stl = sind (ttheta * 0.5) / rlambda 
            dstar = 2. * sind (ttheta * 0.5) / rlambda 
            q = 2. * zpi * sind (ttheta * 0.5) / rlambda 
            lp = lorentz (ttheta) * polarisation (ttheta) 
!DBG_RBN          if(pow_lp_apply) then                                 
!DBG_RBN            if(ttheta.eq.0) then                                
!DBG_RBN              lp = 0.0                                          
!DBG_RBN            ELSE                                                
! ALT FALSCH        lp     = (1+(cosd(ttheta))**2)/(sind(0.5*ttheta))**2
!DBG_RBN              lp = 0.5/sind(0.5*ttheta)/sind(ttheta)            
!DBG_RBN            endif                                               
!DBG_RBN          ELSE                                                  
!DBG_RBN            lp = 1                                              
!DBG_RBN          endif                                                 
            IF (cpow_form.eq.'tth') then 
               WRITE (iff, * ) ttheta, pow_qsp (ii) * lp 
            ELSEIF (cpow_form.eq.'stl') then 
               WRITE (iff, * ) stl, pow_qsp (ii) * lp 
            ELSEIF (cpow_form.eq.'q  ') then 
               WRITE (iff, * ) q, pow_qsp (ii) * lp 
            ELSEIF (cpow_form.eq.'dst') then 
               WRITE (iff, * ) dstar, pow_qsp (ii) * lp 
            ELSEIF (cpow_form.eq.'lop') then 
               WRITE (iff, * ) ttheta, lp 
            ENDIF 
            ENDDO 
         ELSEIF (pow_four_type.eq.POW_DEBYE) then 
            xm (1) = 2 * sind (0.5 * pow_tthmin) / rlambda 
            ss = 2 * sind (0.5 * pow_tthmax) / rlambda 
            st = 2 * sind (0.5 * (pow_tthmax - pow_deltatth) ) /        &
            rlambda                                                     
            uin (1) = (ss - st) / 2. 
            DO ii = 1, num (1) 
            dstar = (xm (1) + (ii - 1) * uin (1) ) 
            stl = .5 * (xm (1) + (ii - 1) * uin (1) ) 
            q = zpi * (xm (1) + (ii - 1) * uin (1) ) 
            ttheta = 2. * asind (dstar * rlambda / 2.) 
!DBG          lp     = lorentz(ttheta)*polarisation(ttheta)             
            lp = polarisation (ttheta) 
            IF (cpow_form.eq.'tth') then 
               WRITE (iff, * ) ttheta, pow_qsp (ii) * lp 
            ELSEIF (cpow_form.eq.'stl') then 
               WRITE (iff, * ) stl, pow_qsp (ii) * lp 
            ELSEIF (cpow_form.eq.'q  ') then 
               WRITE (iff, * ) q, pow_qsp (ii) * lp 
            ELSEIF (cpow_form.eq.'dst') then 
               WRITE (iff, * ) dstar, pow_qsp (ii) * lp 
            ENDIF 
            ENDDO 
         ELSEIF (pow_four_type.eq.POW_FAST) then 
            xm (1) = 2 * sind (0.5 * pow_tthmin) / rlambda 
            ss = 2 * sind (0.5 * pow_tthmax) / rlambda 
            st = 2 * sind (0.5 * (pow_tthmax - pow_deltatth) ) /        &
            rlambda                                                     
            uin (1) = (ss - st) / 2. 
            DO ii = 1, num (1) 
            dstar = (xm (1) + (ii - 1) * uin (1) ) 
            stl = .5 * (xm (1) + (ii - 1) * uin (1) ) 
            q = zpi * (xm (1) + (ii - 1) * uin (1) ) 
            ttheta = 2. * asind (dstar * rlambda / 2.) 
!DBG          lp     = lorentz(ttheta)*polarisation(ttheta)             
            lp = polarisation (ttheta) 
            IF (cpow_form.eq.'tth') then 
               WRITE (iff, * ) ttheta, pow_qsp (ii) * lp 
            ELSEIF (cpow_form.eq.'stl') then 
               WRITE (iff, * ) stl, pow_qsp (ii) * lp 
            ELSEIF (cpow_form.eq.'q  ') then 
               WRITE (iff, * ) q, pow_qsp (ii) * lp 
            ELSEIF (cpow_form.eq.'dst') then 
               WRITE (iff, * ) dstar, pow_qsp (ii) * lp 
            ENDIF 
            ENDDO 
         ELSEIF (pow_four_type.eq.POW_HIST) then 
            IF (pow_axis.eq.POW_AXIS_DSTAR) then 
            ELSEIF (pow_axis.eq.POW_AXIS_Q) then 
               xm (1) = pow_qmin / zpi 
               ss = pow_qmax / zpi 
               st = (pow_qmax - pow_deltaq) / zpi 
               uin (1) = pow_deltaq / zpi 
            ELSEIF (pow_axis.eq.POW_AXIS_TTH) then 
               xm (1) = 2 * sind (0.5 * pow_tthmin) / rlambda 
               ss = 2 * sind (0.5 * pow_tthmax) / rlambda 
               st = 2 * sind (0.5 * (pow_tthmax - pow_deltatth) )       &
               / rlambda                                                
               uin (1) = (ss - st) / 2. 
            ENDIF 
            DO ii = 1, num (1) 
            dstar = (xm (1) + (ii - 1) * uin (1) ) 
            stl = .5 * (xm (1) + (ii - 1) * uin (1) ) 
            q = zpi * (xm (1) + (ii - 1) * uin (1) ) 
            ttheta = 2. * asind (dstar * rlambda / 2.) 
!DBG          lp     = lorentz(ttheta)*polarisation(ttheta)             
            lp = polarisation (ttheta) 
            value = rsf (ii) 
            IF (cpow_form.eq.'tth') then 
               WRITE (iff, * ) ttheta, pow_qsp (ii) * lp 
            ELSEIF (cpow_form.eq.'stl') then 
               WRITE (iff, * ) stl, pow_qsp (ii) * lp 
            ELSEIF (cpow_form.eq.'q  ') then 
               WRITE (iff, * ) q, pow_qsp (ii) * lp 
            ELSEIF (cpow_form.eq.'dst') then 
               WRITE (iff, * ) dstar, pow_qsp (ii) * lp 
            ENDIF 
            ENDDO 
         ENDIF 
         CLOSE (iff) 
      ENDIF 
!                                                                       
      END SUBROUTINE powder_out                     
!*****7*****************************************************************
      REAL function lorentz (ttheta) 
!+                                                                      
!-                                                                      
      USE config_mod 
      USE powder_mod 
      IMPLICIT none 
!                                                                       
      include'errlist.inc' 
!                                                                       
      REAL ttheta 
      REAL r 
!                                                                       
      REAL sind, cosd 
!                                                                       
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
      USE config_mod 
      USE powder_mod 
      IMPLICIT none 
!                                                                       
      include'errlist.inc' 
!                                                                       
      REAL ttheta 
      REAL r 
!                                                                       
      REAL sind, cosd 
!                                                                       
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
      SUBROUTINE powder_conv_res (dat, tthmin, tthmax, dtth, delta,     &
      pow_width, POW_MAXPKT)                                                        
!-                                                                      
!     Convolute powder pattern with resolution function (Gaussian)      
!+                                                                      
      USE config_mod 
      IMPLICIT none 
!                                                                       
      include'wink.inc' 
!
      INTEGER, INTENT(IN) :: POW_MAXPKT
!                                                                       
      REAL dat (0:POW_MAXPKT) 
      REAL tthmin, tthmax, dtth, delta 
      REAL pow_width 
!                                                                       
      REAL dummy (0:POW_MAXPKT) 
      REAL gauss (0:2 * POW_MAXPKT) 
      REAL tth, ysum 
      INTEGER imax, i, j, ii 
      INTEGER max_ps 
!                                                                       
!------ Setup Gaussian                                                  
!                                                                       
      max_ps = (10.0 * delta) / dtth 
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
      imax = (tthmax - tthmin) / dtth 
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
      etax, u, v, w, p1, p2, p3, p4, pow_width, POW_MAXPKT)
!-                                                                      
!     Convolute powder pattern with resolution function (Pseudo-Voigt)  
!     Constant FWHM, Constant eta                                       
!+                                                                      
USE config_mod 
IMPLICIT none 
!                                                                       
include'wink.inc' 
!
INTEGER, INTENT(IN) :: POW_MAXPKT
!                                                                       
REAL dat (0:POW_MAXPKT) 
REAL tthmin, tthmax, dtth, fwhm, eta, etax 
REAL u, v, w 
REAL p1, p2, p3, p4 
REAL pow_width 
!                                                                       
REAL dummy (0:POW_MAXPKT) 
REAL psvgt (0:2 * POW_MAXPKT) 
REAL tth, ysum 
INTEGER imax, i, j, ii 
INTEGER max_ps 
!                                                                       
!REAL pseudovoigt 
!                                                                       
!------ Setup Pseudo-Voigt                                              
!                                                                       
fwhm = sqrt (abs (w) ) 
max_ps = (pow_width * fwhm) / dtth 
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
      imax = (tthmax - tthmin) / dtth 
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
      etax, u, v, w, p1, p2, p3, p4, pow_width, rlambda, pow_axis,      &
      POW_AXIS_Q, POW_MAXPKT)
!-                                                                      
!     Convolute powder pattern with resolution function (Pseudo-Voigt)  
!     FWHM according to caglioti equation, Constant eta                 
!     FWHM = sqrt ( U*tan**2(Theta) + V*tan(Theta) + W)                 
!+                                                                      
      USE config_mod 
      IMPLICIT none 
!                                                                       
      include'wink.inc' 
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
      REAL psvgt (0:2 * POW_MAXPKT) 
      REAL tth, ysum 
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
      imax = (tthmax - tthmin) / dtth 
      DO i = 0, imax 
      tth = tthmin + i * dtth 
      tantth = tand (tth * 0.5) 
      atheta = tth * 0.5 
      atwoth = tth 
      fwhm = sqrt (max (abs (u * tantth**2 + v * tantth + w), 0.00001) ) 
      fwhm1 = fwhm 
      max_ps = (pow_width * fwhm) / dtth 
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
      USE config_mod 
      IMPLICIT none 
!                                                                       
      include'wink.inc' 
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
      REAL psvgt (0:2 * POW_MAXPKT) 
      REAL tth, ysum 
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
      imax = (tthmax - tthmin) / dtth 
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
      max_ps = (pow_width * fwhm) / dtth 
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
END MODULE powder_write_mod
