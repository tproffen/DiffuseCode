!*****7**************************************************************   
!     Here are all the 'math' routines for KUPLOT.                      
!*****7**************************************************************   
      SUBROUTINE do_sort (zeile, lp) 
!+                                                                      
!     sort xy file after increasing x values                            
!-                                                                      
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
      USE prompt_mod 
      USE kuplot_config 
      USE kuplot_mod 
USE precision_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 1) 
!                                                                       
      CHARACTER ( * ) zeile 
      CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
      INTEGER lpara (maxw), lp 
      INTEGER ianz, ik, ir, l, i, j 
      REAL xxx, yyy, dxx, dyy 
      REAL(KIND=PREC_DP) :: werte (maxw) 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      CALL ber_params (ianz, cpara, lpara, werte, maxw) 
      IF (ier_num.ne.0) return 
!                                                                       
      IF (ianz.ne.1) then 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ELSE 
         ik = nint (werte (1) ) 
         IF (ik.lt.iz.and.ik.gt.0) then 
            WRITE (output_io, 1000) ik 
            l = lenc(ik) / 2 + 1 
            ir = lenc(ik) 
   10       CONTINUE 
            IF (l.gt.1) then 
               l = l - 1 
               xxx = x (offxy (ik - 1) + l) 
               yyy = y (offxy (ik - 1) + l) 
               dxx = dx (offxy (ik - 1) + l) 
               dyy = dy (offxy (ik - 1) + l) 
            ELSE 
               xxx = x (offxy (ik - 1) + ir) 
               yyy = y (offxy (ik - 1) + ir) 
               dxx = dx (offxy (ik - 1) + ir) 
               dyy = dy (offxy (ik - 1) + ir) 
               x (offxy (ik - 1) + ir) = x (offxy (ik - 1) + 1) 
               y (offxy (ik - 1) + ir) = y (offxy (ik - 1) + 1) 
               dx (offxy (ik - 1) + ir) = dx (offxy (ik - 1) + 1) 
               dy (offxy (ik - 1) + ir) = dy (offxy (ik - 1) + 1) 
               ir = ir - 1 
               IF (ir.eq.1) then 
                  x (offxy (ik - 1) + 1) = xxx 
                  y (offxy (ik - 1) + 1) = yyy 
                  dx (offxy (ik - 1) + 1) = dxx 
                  dy (offxy (ik - 1) + 1) = dyy 
                  RETURN 
               ENDIF 
            ENDIF 
            i = l 
            j = l + l 
   20       CONTINUE 
            IF (j.le.ir) then 
               IF (j.lt.ir) then 
                  IF (x (offxy (ik - 1) + j) .lt.x (offxy (ik - 1)      &
                  + j + 1) ) j = j + 1                                  
               ENDIF 
               IF (xxx.lt.x (offxy (ik - 1) + j) ) then 
                  x (offxy (ik - 1) + i) = x (offxy (ik - 1) + j) 
                  y (offxy (ik - 1) + i) = y (offxy (ik - 1) + j) 
                  dx (offxy (ik - 1) + i) = dx (offxy (ik - 1) + j) 
                  dy (offxy (ik - 1) + i) = dy (offxy (ik - 1) + j) 
                  i = j 
                  j = j + j 
               ELSE 
                  j = ir + 1 
               ENDIF 
               GOTO 20 
            ENDIF 
            x (offxy (ik - 1) + i) = xxx 
            y (offxy (ik - 1) + i) = yyy 
            dx (offxy (ik - 1) + i) = dxx 
            dy (offxy (ik - 1) + i) = dyy 
            GOTO 10 
         ELSE 
            ier_num = - 4 
            ier_typ = ER_APPL 
         ENDIF 
      ENDIF 
!                                                                       
 1000 FORMAT     (' ------ > Sorting data set ',i3,' ..') 
      END SUBROUTINE do_sort                        
!*****7**************************************************************   
      SUBROUTINE do_four (zeile, lp) 
!+                                                                      
!     Calculates simple discrete Fourier Transform from                 
!     given data set. Range in Fourier space is determined              
!     by sample width in real space.                                    
!-                                                                      
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
      USE kuplot_config 
      USE kuplot_mod 
USE precision_mod
      USE string_convert_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 5) 
!                                                                       
      CHARACTER ( * ) zeile 
      INTEGER lp 
!                                                                       
      CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
      INTEGER lpara (maxw) 
      INTEGER ianz, ik 
      REAL(KIND=PREC_DP) :: werte (maxw) 
      LOGICAL fo (4), lpi 
!                                                                       
!------ get parameters                                                  
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
!                                                                       
      IF (ianz.eq.1.or.ianz.eq.2.or.ianz.eq.4.or.ianz.eq.5) then 
         IF (ianz.eq.2.or.ianz.eq.5) then 
            CALL do_cap (cpara (ianz) ) 
            fo (1) = (index (cpara (ianz) (1:lpara (ianz) ) , 'R')      &
            .ne.0)                                                      
            fo (2) = (index (cpara (ianz) (1:lpara (ianz) ) , 'I')      &
            .ne.0)                                                      
            fo (3) = (index (cpara (ianz) (1:lpara (ianz) ) , 'A')      &
            .ne.0)                                                      
            fo (4) = (index (cpara (ianz) (1:lpara (ianz) ) , 'P')      &
            .ne.0)                                                      
            lpi = (index (cpara (ianz) (1:lpara (ianz) ) , 'H') .ne.0) 
            ianz = ianz - 1 
!                                                                       
            IF (.not. (fo (1) .or.fo (2) .or.fo (3) .or.fo (4) ) ) then 
               fo (1) = .true. 
            ENDIF 
         ELSE 
            fo (1) = .true. 
            fo (2) = .false. 
            fo (3) = .false. 
            fo (4) = .false. 
            lpi = .false. 
         ENDIF 
!                                                                       
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         ik = nint (werte (1) ) 
         IF (ik.lt.iz.and.ik.gt.0) then 
            IF (lni (ik) ) then 
               IF (ianz.le.2) then 
                  CALL do_four_z (ik, fo, lpi) 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
            ELSE 
               CALL do_four_y (ik, fo, lpi, ianz, werte, maxw) 
            ENDIF 
         ELSE 
            ier_num = - 4 
            ier_typ = ER_APPL 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
      END SUBROUTINE do_four                        
!*****7**************************************************************** 
      SUBROUTINE do_four_y (ik, fo, lpi, ianz, werte, maxw) 
!                                                                       
!     This routine perfoms the Fourier Transform for 2D data            
!                                                                       
      USE errlist_mod 
      USE prompt_mod 
      USE wink_mod
      USE kuplot_config 
      USE kuplot_mod 
USE precision_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER(10) oname (4) 
      INTEGER ianz, maxw, i, j, ipkt, inew, ik, k, kk 
      INTEGER maxpkt 
      INTEGER ift (0:4) 
      REAL dxx, h, hmin, hmax, dh, fc, f (4) 
      REAL(KIND=PREC_DP) :: werte (maxw) 
      LOGICAL fo (4), lpi 
!                                                                       
      DATA oname / 'real.xy', 'imag.xy', 'ampl.xy', 'phase.xy' / 
!                                                                       
!------ Enough space for new data sets                                  
!                                                                       
      inew = 0 
      ift (0) = ik 
      DO i = 1, 4 
      IF (fo (i) ) then 
         inew = inew + 1 
         ift (i) = iz - 1 + inew 
      ENDIF 
      ENDDO 
!                                                                       
      IF (inew.eq.0) then 
         ier_num = - 6 
         ier_typ = ER_COMM 
         RETURN 
      ENDIF 
!                                                                       
      IF ( (iz + inew - 1) .gt.maxkurvtot) then 
         ier_num = - 1 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
!                                                                       
!------ Figure out range                                                
!                                                                       
      dxx = x (offxy (ik - 1) + 2) - x (offxy (ik - 1) + 1) 
      IF (ianz.eq.4) then 
         hmin = werte (2) 
         hmax = werte (3) 
         dh = werte (4) 
      ELSE 
         hmin = - lenc(ik) / dxx 
         hmax = lenc(ik) / dxx 
         dh = 1.0 / (lenc(ik) * dxx) 
      ENDIF 
      ipkt = nint ( (hmax - hmin) / dh) + 1 
!                                                                       
      maxpkt = maxarray - offxy (iz - 1) 
      IF (inew * ipkt.gt.maxpkt) then 
         ier_num = - 6 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
!                                                                       
!------ Setting properties of newly created data sets                   
!                                                                       
      DO i = 1, 4 
      IF (fo (i) ) then 
         fname (ift (i) ) = oname (i) 
         fform (ift (i) ) = 'XY' 
         lni (ift (i) ) = .false. 
         lenc(ift (i) ) = ipkt 
         offxy (ift (i) ) = offxy (ift (i) - 1) + lenc(ift (i) ) 
         offz (ift (i) ) = offz (ift (i) - 1) 
         iz = iz + 1 
      ENDIF 
      ENDDO 
!                                                                       
!------ Calculating the Fourier Transform                               
!                                                                       
      IF (lpi) then 
         fc = REAL(zpi) 
      ELSE 
         fc = 1. 
      ENDIF 
!                                                                       
      WRITE (output_io, 1000) ik, lpi 
!                                                                       
      kk = 0 
      DO k = 1, ipkt 
      h = hmin + dh * (k - 1) 
      f (1) = 0.0 
      f (2) = 0.0 
      DO j = 1, lenc(ik) 
      f (1) = f (1) + y (offxy (ik - 1) + j) * cos (fc * h * x (offxy ( &
      ik - 1) + j) )                                                    
      f (2) = f (2) + y (offxy (ik - 1) + j) * sin (fc * h * x (offxy ( &
      ik - 1) + j) )                                                    
      ENDDO 
      f (1) = dxx * f (1) 
      f (2) = dxx * f (2) 
      f (3) = sqrt (f (1) **2 + f (2) **2) 
      f (4) = atan2(f (2) , f (1) ) / REAL(rad) 
      kk = kk + 1 
!                                                                       
      DO i = 1, 4 
      IF (fo (i) ) then 
         x (offxy (ift (i) - 1) + kk) = h 
         y (offxy (ift (i) - 1) + kk) = f (i) 
         dx (offxy (ift (i) - 1) + kk) = 0.0 
         dy (offxy (ift (i) - 1) + kk) = 0.0 
      ENDIF 
      ENDDO 
      ENDDO 
!                                                                       
      CALL get_extrema 
!                                                                       
 1000 FORMAT     (' ------ > Calulating Fourier Transform of data set ',&
     &                              I2,' (2pi=',L1,') ...')             
      END SUBROUTINE do_four_y                      
!*****7**************************************************************** 
      SUBROUTINE do_four_z (ik, fo, lpi) 
!                                                                       
!     This routine perfoms the Fourier Transform for 3D data            
!                                                                       
      USE errlist_mod 
      USE prompt_mod 
      USE wink_mod
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      CHARACTER(12) oname (4) 
      INTEGER i, jx, jy, ipkt, jpkt, inew, ik, kkx, kky 
      INTEGER kx, ky, maxpkt, maxzz 
      INTEGER ift (0:4) 
      REAL dxx, dyy, hx, hy, fc, f (4) 
      LOGICAL fo (4), lpi 
!                                                                       
      DATA oname / 'real.nipl', 'imag.nipl', 'ampl.nipl', 'phase.nipl' / 
!                                                                       
!------ Enough space for new data sets                                  
!                                                                       
      inew = 0 
      ift (0) = ik 
      DO i = 1, 4 
      IF (fo (i) ) then 
         inew = inew + 1 
         ift (i) = iz - 1 + inew 
      ENDIF 
      ENDDO 
!                                                                       
      IF (inew.eq.0) then 
         ier_num = - 6 
         ier_typ = ER_COMM 
         RETURN 
      ENDIF 
!                                                                       
      IF ( (iz + inew - 1) .gt.maxkurvtot) then 
         ier_num = - 1 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
!                                                                       
      maxpkt = maxarray - offxy (iz - 1) 
      maxzz = maxarray - offz (iz - 1) 
!                                                                       
      ipkt = 2 * (nx (ik) / 2) + 1 
      jpkt = 2 * (ny (ik) / 2) + 1 
      IF (nx (ik) .le.0.or.ny (ik) .le.0.or.inew * nx (ik) * ny (ik)    &
      .gt.maxzz.or.inew * max (nx (ik), ny (ik) ) .gt.maxpkt) then      
         ier_num = - 6 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
!                                                                       
!------ Setting properties of newly created data sets                   
!                                                                       
      DO i = 1, 4 
      IF (fo (i) ) then 
         fname (ift (i) ) = oname (i) 
         fform (ift (i) ) = 'NI' 
         lni (ift (i) ) = .true. 
         nx (ift (i) ) = ipkt 
         ny (ift (i) ) = jpkt 
         lenc(ift (i) ) = max (ipkt, jpkt) 
         offxy (ift (i) ) = offxy (ift (i) - 1) + lenc(ift (i) ) 
         offz (ift (i) ) = offz (ift (i) - 1) + nx (ift (i) ) * ny (ift &
         (i) )                                                          
         iz = iz + 1 
      ENDIF 
      ENDDO 
!                                                                       
!------ Calculating the Fourier Transform                               
!                                                                       
      WRITE (output_io, 1000) ik 
!                                                                       
      IF (lpi) then 
         fc = REAL(zpi) 
      ELSE 
         fc = 1. 
      ENDIF 
!                                                                       
      kkx = 0 
      dxx = x (offxy (ik - 1) + 2) - x (offxy (ik - 1) + 1) 
      dyy = y (offxy (ik - 1) + 2) - y (offxy (ik - 1) + 1) 
      DO kx = - nx (ik) / 2, nx (ik) / 2 
      kkx = kkx + 1 
      kky = 0 
      hx = kx / (nx (ik) * dxx) 
      DO ky = - ny (ik) / 2, ny (ik) / 2 
      kky = kky + 1 
      hy = ky / (ny (ik) * dyy) 
      f (1) = 0.0 
      f (2) = 0.0 
      DO jx = 1, nx (ik) 
      DO jy = 1, ny (ik) 
      f (1) = f (1) + z (offz (ik - 1) + (jx - 1) * ny (ik) + jy)       &
      * cos (fc * (hx * x (offxy (ik - 1) + jx) + hy * y (offxy (ik - 1)&
      + jy) ) )                                                         
      f (2) = f (2) + z (offz (ik - 1) + (jx - 1) * ny (ik) + jy)       &
      * sin (fc * (hx * x (offxy (ik - 1) + jx) + hy * y (offxy (ik - 1)&
      + jy) ) )                                                         
      ENDDO 
      ENDDO 
      f (1) = dxx * dyy * f (1) 
      f (2) = dxx * dyy * f (2) 
      f (3) = sqrt (f (1) **2 + f (2) **2) 
      f (4) = atan2(f (2) , f (1) ) / REAL(rad) 
!                                                                       
      DO i = 1, 4 
      IF (fo (i) ) then 
         z (offz (ift (i) - 1) + (kkx - 1) * ny (ift (i) ) + kky)       &
         = f (i)                                                        
         x (offxy (ift (i) - 1) + kkx) = hx 
         y (offxy (ift (i) - 1) + kky) = hy 
         dx (offxy (ift (i) - 1) + kkx) = 0.0 
         dy (offxy (ift (i) - 1) + kky) = 0.0 
      ENDIF 
      ENDDO 
      ENDDO 
      ENDDO 
!                                                                       
      CALL get_extrema 
!                                                                       
 1000 FORMAT     (' ------ > Calulating Fourier Transform of data set ',&
     &                              I2,' ...')                          
      END SUBROUTINE do_four_z                      
!*****7**************************************************************** 
      SUBROUTINE do_glat (zeile, lp, lsmooth) 
!+                                                                      
!     Smooth data set                                                   
!-                                                                      
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
      USE prompt_mod 
      USE kuplot_config 
      USE kuplot_mod 
USE precision_mod
      USE string_convert_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw, maxsm 
      PARAMETER (maxw = 4) 
      PARAMETER (maxsm = 501) 
!                                                                       
      CHARACTER ( * ) zeile 
      INTEGER lp 
      LOGICAL lsmooth 
!                                                                       
      CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
      CHARACTER(2) dir 
      INTEGER lpara (maxw) 
      INTEGER ianz, ik, ip, im, il 
      REAL(KIND=PREC_DP) :: werte (maxw) 
      REAL cc (maxsm) 
      LOGICAL gx, gy 
!                                                                       
      IF (lsmooth) then 
         DO ip = 1, maxsm 
         cc (ip) = 0.0 
         ENDDO 
      ELSE 
         DO ip = 1, maxsm 
         cc (ip) = 1.0 
         ENDDO 
      ENDIF 
!                                                                       
!------ get parameters                                                  
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
!                                                                       
      IF (ianz.eq.2.or.ianz.eq.3.or.ianz.eq.4) then 
         IF (lsmooth.and.ianz.eq.4) then 
            CALL do_cap (cpara (4) ) 
            gx = (cpara (4) (1:1) .eq.'X') 
            gy = (cpara (4) (1:1) .eq.'Y') 
            cpara (4) = '0.0' 
         ELSEIF (.not.lsmooth.and.ianz.eq.3) then 
            CALL do_cap (cpara (3) ) 
            gx = (cpara (3) (1:1) .eq.'X') 
            gy = (cpara (3) (1:1) .eq.'Y') 
            cpara (3) = '0.0' 
         ELSE 
            gx = .true. 
            gy = .true. 
         ENDIF 
!                                                                       
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         ik = nint (werte (1) ) 
         ip = nint (werte (2) ) 
         IF (ik.lt.iz.and.ik.gt.0) then 
!                                                                       
            IF (gx) dir = ' x' 
            IF (gy) dir = ' y' 
            IF (gx.and.gy) dir = 'xy' 
!                                                                       
            IF (lsmooth) then 
               il = 0 
               im = 2 
               IF (ianz.gt.2) im = nint (werte (3) ) 
               CALL SAVGOL (cc, ip, ip / 2, ip / 2, il, im) 
               IF (ier_num.ne.0) return 
            ENDIF 
!                                                                       
            IF (lni (ik) ) then 
               IF (lsmooth) then 
                  WRITE (output_io, 2000) ik, ip, im 
               ELSE 
                  WRITE (output_io, 2010) ik, ip 
               ENDIF 
               WRITE (output_io, 2100) dir 
               CALL do_glatt_z (ik, ip, nx (ik), ny (ik), gx, gy, cc,   &
               maxsm, lsmooth)                                          
            ELSE 
               IF (lsmooth) then 
                  WRITE (output_io, 2000) ik, ip, im 
               ELSE 
                  WRITE (output_io, 2010) ik, ip 
               ENDIF 
               CALL do_glatt_y (ik, ip, lenc(ik), cc, maxsm, lsmooth) 
            ENDIF 
         ELSE 
            ier_num = - 4 
            ier_typ = ER_APPL 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
 2000 FORMAT (' Smoothing data set ',i3,' via Savitzky-Golay :',/,      &
     &        '   Number of points           : ',i3,/,                  &
     &        '   Order of smoothing polynom : ',i3)                    
 2010 FORMAT (' Smoothing data set ',i3,' via sliding average :',/,     &
     &        '   Number of points           : ',i3)                    
 2100 FORMAT ('   Direction                  : ',1x,a2) 
      END SUBROUTINE do_glat                        
!*****7*****************************************************************
      SUBROUTINE do_inte (zeile, lp) 
!+                                                                      
!     Compute integral for 2d/3d data sets                              
!-                                                                      
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
      USE param_mod 
      USE prompt_mod 
      USE kuplot_config 
      USE kuplot_mod 
USE precision_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 5) 
!                                                                       
      CHARACTER ( * ) zeile 
      CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
      INTEGER lpara (maxw), lp 
      INTEGER ianz, ipkt, ip, ik, ix, iy 
      REAL rint, drint, rdx, rdy, xa, ya, xs 
      REAL wx1, wx2, wy1, wy2 
      REAL(KIND=PREC_DP) :: werte (maxw) 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      CALL ber_params (ianz, cpara, lpara, werte, maxw) 
      IF (ianz.lt.1) then 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
      IF (ier_num.ne.0) return 
      CALL skalieren 
      ik = nint (werte (1) ) 
      wx1 = ex (iwin, iframe, 1) 
      wx2 = ex (iwin, iframe, 2) 
      wy1 = ey (iwin, iframe, 1) 
      wy2 = ey (iwin, iframe, 2) 
!                                                                       
      IF (ik.gt.0.and.ik.lt.iz) then 
         IF (ianz.eq.1) then 
            CONTINUE 
         ELSEIF (ianz.eq.3) then 
            wx1 = werte (2) 
            wx2 = werte (3) 
         ELSEIF (ianz.eq.5.and.lni (ik) ) then 
            wx1 = werte (2) 
            wx2 = werte (3) 
            wy1 = werte (4) 
            wy2 = werte (5) 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
      ELSE 
         ier_num = - 4 
         ier_typ = ER_APPL 
      ENDIF 
!                                                                       
      IF (ier_num.ne.0) return 
!                                                                       
      rint = 0.0 
      drint = 0.0 
      ipkt = 0 
!                                                                       
!------ 2D integration                                                  
!                                                                       
      IF (lni (ik) ) then 
         DO ix = 1, nx (ik) 
         DO iy = 1, ny (ik) 
         IF (x (offxy (ik - 1) + ix) .ge.wx1.and.y (offxy (ik - 1)      &
         + iy) .ge.wy1.and.x (offxy (ik - 1) + ix) .le.wx2.and.y (offxy &
         (ik - 1) + iy) .le.wy2) then                                   
            zwert = z (offz (ik - 1) + (ix - 1) * ny (ik) + iy) 
            IF (zwert.ne. - 9999.0) then 
               rint = rint + zwert 
               ipkt = ipkt + 1 
            ENDIF 
         ENDIF 
         ENDDO 
         ENDDO 
         rdx = (xmax (ik) - xmin (ik) ) / REAL(nx (ik) ) 
         rdy = (ymax (ik) - ymin (ik) ) / REAL(ny (ik) ) 
         WRITE (output_io, 2000) ik, wx1, wx2, wy1, wy2, rint * rdx *   &
         rdy, drint, ipkt                                               
         IF (maxpar_res.lt.7) then 
            ier_num = - 24 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
         res_para (0) = 7 
         res_para (1) = rint * rdx * rdy 
         res_para (2) = drint 
         res_para (3) = REAL(ipkt) 
         res_para (4) = wx1 
         res_para (5) = wx2 
         res_para (6) = wy1 
         res_para (7) = wy2 
!                                                                       
!------ 1D integration                                                  
!                                                                       
      ELSE 
         DO ix = 2, lenc(ik) 
         ip = offxy (ik - 1) + ix 
         xa = (x (ip) + x (ip - 1) ) / 2.0 
         ya = (y (ip) + y (ip - 1) ) / 2.0 
         xs = abs (x (ip) - x (ip - 1) ) 
!                                                                       
         IF (xa.ge.wx1.and.xa.le.wx2) then 
            rint = rint + xs * ya 
            drint = drint + xs**2 * (dy (ip) **2 + dy (ip - 1) **2)     &
            / 4.0                                                       
            ipkt = ipkt + 1 
         ENDIF 
         ENDDO 
         WRITE (output_io, 2100) ik, wx1, wx2, rint, sqrt (drint),      &
         ipkt                                                           
         IF (maxpar_res.lt.5) then 
            ier_num = - 24 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
         res_para (0) = 5 
         res_para (1) = rint 
         res_para (2) = sqrt (drint) 
         res_para (3) = REAL(ipkt) 
         res_para (4) = wx1 
         res_para (5) = wx2 
      ENDIF 
      WRITE (output_io, * ) 
!                                                                       
 2000 FORMAT    (' Integration result for data set ',i3,' :',/          &
     &                  3x,'x-range    : ',g12.4,' to ',g12.4,/         &
     &                  3x,'y-range    : ',g12.4,' to ',g12.4,/         &
     &         3x,'Integral   : ',g12.4,' +- ',g12.4,'  (',i7,' pkt)')  
 2100 FORMAT    (' Integration result for data set ',i3,' :',/          &
     &                  3x,'x-range    : ',g12.4,' to ',g12.4,/         &
     &         3x,'Integral   : ',g12.4,' +- ',g12.4,'  (',i7,' pkt)')  
!                                                                       
      END SUBROUTINE do_inte                        
!**7******************************************************************* 
      SUBROUTINE do_mean (zeile, lp, lout) 
!+                                                                      
!     calculates mean and standard deviation                            
!-                                                                      
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
      USE param_mod 
      USE prompt_mod 
      USE kuplot_config 
      USE kuplot_mod 
USE precision_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 5) 
!                                                                       
      CHARACTER ( * ) zeile 
      CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
      INTEGER lpara (maxw), lp 
      INTEGER ianz, ik, ix, iy, i 
      REAL(KIND=PREC_DP) :: werte (maxw) 
      REAL wx1, wx2, wy1, wy2 
      REAL sum (3), sum2 (3), wert (3), ave (3), sig (3), n 
      LOGICAL lout 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      CALL ber_params (ianz, cpara, lpara, werte, maxw) 
      IF (ianz.lt.1) then 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
      IF (ier_num.ne.0) return 
!                                                                       
      ik = nint (werte (1) ) 
      IF (ik.ge.iz.or.lenc(ik) .eq.0) then 
         ier_num = - 4 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
!                                                                       
      CALL skalieren 
      IF (ianz.eq.1) then 
         wx1 = ex (iwin, iframe, 1) 
         wx2 = ex (iwin, iframe, 2) 
         wy1 = ey (iwin, iframe, 1) 
         wy2 = ey (iwin, iframe, 2) 
      ELSEIF (ianz.eq.3) then 
         wx1 = werte (2) 
         wx2 = werte (3) 
         wy1 = ey (iwin, iframe, 1) 
         wy2 = ey (iwin, iframe, 2) 
      ELSEIF (ianz.eq.5.and.lni (ik) ) then 
         wx1 = werte (2) 
         wx2 = werte (3) 
         wy1 = werte (4) 
         wy2 = werte (5) 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
         RETURN 
      ENDIF 
!                                                                       
      DO i = 1, 3 
      sum (i) = 0.0 
      sum2 (i) = 0.0 
      ENDDO 
      n = 0 
!                                                                       
      IF (lni (ik) ) then 
         DO ix = 1, nx (ik) 
         DO iy = 1, ny (ik) 
         IF (x (offxy (ik - 1) + ix) .ge.wx1.and.y (offxy (ik - 1)      &
         + iy) .ge.wy1.and.x (offxy (ik - 1) + ix) .le.wx2.and.y (offxy &
         (ik - 1) + iy) .le.wy2) then                                   
            wert (1) = x (offxy (ik - 1) + ix) 
            wert (2) = y (offxy (ik - 1) + iy) 
            wert (3) = z (offz (ik - 1) + (ix - 1) * ny (ik) + iy) 
            IF (werte (3) .ne. - 9999.0) then 
               DO i = 1, 3 
               sum (i) = sum (i) + wert (i) 
               sum2 (i) = sum2 (i) + wert (i) **2 
               ENDDO 
               n = n + 1 
            ENDIF 
         ENDIF 
         ENDDO 
         ENDDO 
         IF (n.gt.0) then 
            DO i = 1, 3 
            ave (i) = sum (i) / n 
            sig (i) = sqrt (sum2 (i) / n - (sum (i) / n) **2) 
            ENDDO 
            IF (lout) write (output_io, 2000) ik, wx1, wx2, wy1, wy2,   &
            int (n), (ave (i), sig (i), i = 1, 3)                       
            IF (maxpar_res.lt.11) then 
               ier_num = - 24 
               ier_typ = ER_APPL 
               RETURN 
            ENDIF 
            res_para (0) = 11 
            res_para (1) = ave (1) 
            res_para (2) = ave (2) 
            res_para (3) = ave (3) 
            res_para (4) = sig (1) 
            res_para (5) = sig (2) 
            res_para (6) = sig (3) 
            res_para (7) = int (n) 
            res_para (8) = wx1 
            res_para (9) = wx2 
            res_para (10) = wy1 
            res_para (11) = wy2 
         ELSE 
            ier_num = - 20 
            ier_typ = ER_APPL 
         ENDIF 
      ELSE 
         DO ix = 1, lenc(ik) 
         IF (x (offxy (ik - 1) + ix) .ge.wx1.and.x (offxy (ik - 1)      &
         + ix) .le.wx2) then                                            
            wert (1) = x (offxy (ik - 1) + ix) 
            wert (2) = y (offxy (ik - 1) + ix) 
            DO i = 1, 2 
            sum (i) = sum (i) + wert (i) 
            sum2 (i) = sum2 (i) + wert (i) **2 
            ENDDO 
            n = n + 1 
         ENDIF 
         ENDDO 
         IF (n.gt.0) then 
            DO i = 1, 2 
            ave (i) = sum (i) / n 
            sig (i) = sqrt (sum2 (i) / n - (sum (i) / n) **2) 
            ENDDO 
            IF (lout) write (output_io, 2100) ik, wx1, wx2, int (n),    &
            (ave (i), sig (i), i = 1, 2)                                
            IF (maxpar_res.lt.7) then 
               ier_num = - 24 
               ier_typ = ER_APPL 
               RETURN 
            ENDIF 
            res_para (0) = 7 
            res_para (1) = ave (1) 
            res_para (2) = ave (2) 
            res_para (3) = sig (1) 
            res_para (4) = sig (2) 
            res_para (5) = int (n) 
            res_para (6) = wx1 
            res_para (7) = wx2 
         ELSE 
            ier_num = - 20 
            ier_typ = ER_APPL 
         ENDIF 
      ENDIF 
      IF (lout) write (output_io, * ) 
!                                                                       
 2000 FORMAT    (' Average and standard deviation data set ',i3,' :',/  &
     &                  3x,'x-range   : ',g12.4,' to ',g12.4,/          &
     &                  3x,'y-range   : ',g12.4,' to ',g12.4,/          &
     &                  3x,'# points  : ',i7//                          &
     &         3x,'Average x : ',g12.4,' sigma :  ',g12.4,/             &
     &         3x,'Average y : ',g12.4,' sigma :  ',g12.4,/             &
     &         3x,'Average z : ',g12.4,' sigma :  ',g12.4)              
 2100 FORMAT    (' Average and standard deviation data set ',i3,' :',/  &
     &                  3x,'x-range   : ',g12.4,' to ',g12.4,/          &
     &                  3x,'# points  : ',i7//                          &
     &         3x,'Average x : ',g12.4,' sigma :  ',g12.4,/             &
     &         3x,'Average y : ',g12.4,' sigma :  ',g12.4)              
!                                                                       
      END SUBROUTINE do_mean                        
!**7******************************************************************* 
      SUBROUTINE do_smax (zeile, lp) 
!+                                                                      
!     Search for maxima in data set                                     
!-                                                                      
      USE ber_params_mod
      USE build_name_mod
      USE errlist_mod 
      USE get_params_mod
      USE param_mod 
      USE prompt_mod 
      USE kuplot_config 
      USE kuplot_mod 
USE lib_length
USE precision_mod
USE support_mod
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxmax, maxw 
      PARAMETER (maxmax = 50) 
      PARAMETER (maxw = 10) 
!                                                                       
      CHARACTER ( * ) zeile 
      CHARACTER(LEN=PREC_STRING) :: cpara (maxw) 
      CHARACTER(80) cdummy 
      CHARACTER(6) cout 
      REAL(KIND=PREC_DP) :: werte (maxw) 
      REAL wmax (maxmax) 
      INTEGER lpara (maxw), lp 
      INTEGER ixm (maxmax), iym (maxmax) 
      INTEGER ianz, ifil, im, ima, ik 
!                                                                       
!                                                                       
      cout = 'maxima' 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
!                                                                       
      IF (ianz.ge.2) then 
         CALL ber_params (2, cpara, lpara, werte, maxw) 
         IF (ier_num.ne.0) return 
         ik = nint (werte (1) ) 
         ifen = nint (werte (2) ) 
         IF (ik.gt.0.and.ik.lt.iz) then 
            IF (ianz.gt.2) then 
               CALL del_params (2, ianz, cpara, lpara, maxw) 
               CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
               IF (ier_num.ne.0) return 
               ifil = 66 
               cdummy = cpara (1) (1:MIN(80,LEN_TRIM(cpara(1))))
               CALL oeffne (ifil, cdummy, 'unknown') 
               IF (ier_num.ne.0) return 
               WRITE (output_io, 900) cdummy (1:len_str (cdummy) ) 
            ELSE 
               ifil = output_io 
            ENDIF 
!                                                                       
!------ --- 3D data set                                                 
!                                                                       
            IF (lni (ik) ) then 
               CALL do_fmax_z (ik, wmax, ixm, iym, maxmax, ima) 
               IF (ier_num.eq.0) then 
                  WRITE (ifil, 1000) cout, ik, ifen 
                  IF (maxpar_res.lt.2 * ima) then 
                     ier_num = - 24 
                     ier_typ = ER_APPL 
                     RETURN 
                  ENDIF 
                  res_para (0) = 2.0 * ima 
                  DO im = 1, ima 
                  WRITE (ifil, 1100) im, x (offxy (ik - 1) + ixm (im) ),&
                  y (offxy (ik - 1) + iym (im) ), wmax (im)             
                  res_para (2 * im - 1) = REAL(ixm (im) ) 
                  res_para (2 * im) = REAL(iym (im) ) 
                  ENDDO 
               ENDIF 
!                                                                       
!------ --- 2D data set                                                 
!                                                                       
            ELSE 
               CALL do_fmax_xy (ik, wmax, ixm, maxmax, ima) 
               IF (ier_num.eq.0) then 
                  WRITE (ifil, 2000) cout, ik, ifen 
                  IF (maxpar_res.lt.ima) then 
                     ier_num = - 24 
                     ier_typ = ER_APPL 
                     RETURN 
                  ENDIF 
                  res_para (0) = ima 
                  DO im = 1, ima 
                  WRITE (ifil, 2100) im, x (offxy (ik - 1) + ixm (im) ),&
                  wmax (im)                                             
                  res_para (im) = REAL(ixm (im) ) 
                  ENDDO 
               ENDIF 
            ENDIF 
            IF (ifil.ne.6) close (ifil) 
         ELSE 
            ier_num = - 4 
            ier_typ = ER_APPL 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
  900 FORMAT (' ------ > Saving results to file ',a,' ..') 
 1000 FORMAT (' Found ',a,' data set ',i3,' (ifen = ',i3,') :',//,      &
     &        3x,'No.        pos. x       pos. y         value'/        &
     &        3x,'--------------------------------------------------')  
 1100 FORMAT (3x,i3,2x,f11.3,2x,f11.3,2x,f18.3) 
 2000 FORMAT (' Found ',a,' data set ',i3,' (ifen = ',i3,') :',//,      &
     &        3x,'No.        pos. x         value'/                     &
     &        3x,'-------------------------------------')               
 2100 FORMAT (3x,i3,2x,f11.3,2x,f18.3) 
      END SUBROUTINE do_smax                        
!**7*****************************************************************   
!     SUPPORT ROUTINES                                                  
!**7*****************************************************************   
      SUBROUTINE do_glatt_z (ik, np, nxx, nyy, gx, gy, cc, mc, lsm) 
!+                                                                      
!     glaetten des nipl-2d-arrays                                       
!-                                                                      
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER mc 
!                                                                       
      REAL a (maxarray), cc (mc) 
      REAL xnw, s 
      INTEGER ik, np, nxx, nyy 
      INTEGER nd, na, ne, ip, jp, im, jm, i, j, k 
      LOGICAL gx, gy, lsm 
!                                                                       
      nd = (np - 1) / 2 
      na = nd+1 
!                                                                       
!------ glaettung in x-richtung                                         
!                                                                       
      IF (.not.gx) goto 222 
      ne = nxx - nd 
      DO j = 1, nyy 
      DO i = na, ne 
      xnw = np 
      im = i - nd 
      ip = i + nd 
      s = 0.0 
      DO k = im, ip 
      IF (z (offz (ik - 1) + (k - 1) * ny (ik) + j) .ne. - 9999.0.and.z &
      (offz (ik - 1) + (k - 1) * ny (ik) + j) .ne.0.0) then             
         s = s + cc (abs (i - k) + 1) * z (offz (ik - 1) + (k - 1)      &
         * ny (ik) + j)                                                 
      ELSE 
         xnw = xnw - 1 
      ENDIF 
      ENDDO 
      IF (xnw.gt.0) then 
         IF (lsm) then 
            a (i) = s 
         ELSE 
            a (i) = s / xnw 
         ENDIF 
      ELSE 
         a (i) = - 9999.0 
      ENDIF 
      ENDDO 
      DO i = na, ne 
      z (offz (ik - 1) + (i - 1) * ny (ik) + j) = a (i) 
      ENDDO 
      ENDDO 
!                                                                       
!------ glaettung in y-richtung                                         
!                                                                       
  222 CONTINUE 
      IF (.not.gy) goto 333 
      ne = nyy - nd 
      DO i = 1, nxx 
      DO j = na, ne 
      xnw = np 
      jm = j - nd 
      jp = j + nd 
      s = 0.0 
      DO k = jm, jp 
      IF (z (offz (ik - 1) + (i - 1) * ny (ik) + k) .ne. - 9999.0.and.z &
      (offz (ik - 1) + (i - 1) * ny (ik) + k) .ne.0.0) then             
         s = s + cc (abs (j - k) + 1) * z (offz (ik - 1) + (i - 1)      &
         * ny (ik) + k)                                                 
      ELSE 
         xnw = xnw - 1 
      ENDIF 
      ENDDO 
      IF (xnw.gt.0) then 
         IF (lsm) then 
            a (j) = s 
         ELSE 
            a (j) = s / xnw 
         ENDIF 
      ELSE 
         a (j) = - 9999.0 
      ENDIF 
      ENDDO 
      DO j = na, ne 
      z (offz (ik - 1) + (i - 1) * ny (ik) + j) = a (j) 
      ENDDO 
      ENDDO 
!                                                                       
  333 CONTINUE 
      END SUBROUTINE do_glatt_z                     
!**7*****************************************************************   
      SUBROUTINE do_glatt_y (ik, ip, llen, cc, mc, lsm) 
!+                                                                      
!     glaetten der kurve ik ueber ip punkte                             
!-                                                                      
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER mc 
!                                                                       
      REAL ay (maxarray), ady (maxarray), cc (mc) 
      REAL s, ds, xnw 
      INTEGER ik, ip, llen 
      INTEGER na, ne, nd, ic, i, k, kk 
      LOGICAL lsm 
!                                                                       
      nd = (ip - 1) / 2 
      na = nd+1 
      ne = llen - nd 
!                                                                       
      DO i = na, ne 
      xnw = ip 
      s = 0.0 
      ds = 0.0 
!                                                                       
      DO k = - nd, 0 
      kk = k + i 
      ic = 1 - k 
      s = s + cc (ic) * y (offxy (ik - 1) + kk) 
      ds = ds + cc (ic) **2 * dy (offxy (ik - 1) + kk) **2 
      ENDDO 
!                                                                       
      DO k = 1, nd 
      kk = k + i 
      ic = 2 * nd+2 - k 
      s = s + cc (ic) * y (offxy (ik - 1) + kk) 
      ds = ds + cc (ic) **2 * dy (offxy (ik - 1) + kk) **2 
      ENDDO 
!                                                                       
      IF (lsm) then 
         ay (i) = s 
         ady (i) = sqrt (ds) 
      ELSE 
         ay (i) = s / REAL(ip) 
         ady (i) = sqrt (ds / REAL(ip) ) 
      ENDIF 
      ENDDO 
!                                                                       
      DO i = na, ne 
      y (offxy (ik - 1) + i) = ay (i) 
      dy (offxy (ik - 1) + i) = ady (i) 
      ENDDO 
!                                                                       
      END SUBROUTINE do_glatt_y                     
!**7*****************************************************************   
      SUBROUTINE get_extrema_xy_local (i, ymi, yma) 
!+                                                                      
!     extrema der kurve ik local bestimmen                              
!-                                                                      
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      REAL ymi, yma 
      INTEGER i, ip 
!                                                                       
      yma = - 1e38 
      ymi = 1e38 
      DO ip = 1, lenc(i) 
      IF (x (offxy (i - 1) + ip) .ge.ex (iwin, iframe, 1) .and.x (offxy &
      (i - 1) + ip) .le.ex (iwin, iframe, 2) ) then                     
         yma = max (yma, y (offxy (i - 1) + ip) ) 
         ymi = min (ymi, y (offxy (i - 1) + ip) ) 
      ENDIF 
      ENDDO 
!                                                                       
      END SUBROUTINE get_extrema_xy_local           
!**7*****************************************************************   
      SUBROUTINE get_extrema_xy (a, ik, ilen, amin, amax) 
!+                                                                      
!     extrema der kurve ik bestimmen                                    
!-                                                                      
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      REAL a (maxarray) 
      REAL amax (maxkurvtot), amin (maxkurvtot) 
      INTEGER ik, ip, ilen 
!                                                                       
      amax (ik) = a (offxy (ik - 1) + 1) 
      amin (ik) = a (offxy (ik - 1) + 1) 
      DO ip = 2, ilen 
      amax (ik) = max (amax (ik), a (offxy (ik - 1) + ip) ) 
      amin (ik) = min (amin (ik), a (offxy (ik - 1) + ip) ) 
      ENDDO 
!                                                                       
      END SUBROUTINE get_extrema_xy                 
!**7*****************************************************************   
      SUBROUTINE get_extrema_z (a, ik, nxx, nyy, amin, amax) 
!+                                                                      
!     extrema der kurve ik bestimmen                                    
!-                                                                      
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      REAL a (maxarray) 
      REAL amax (maxkurvtot), amin (maxkurvtot) 
      INTEGER ik, nxx, nyy, ikk, i, j 
!                                                                       
      amax (ik) = a (offz (ik - 1) + 1) 
      amin (ik) = a (offz (ik - 1) + 1) 
      IF (amin (ik) .eq. - 9999.) amin (ik) = 1e+30 
      DO i = 1, nxx 
      DO j = 2, nyy 
      ikk = offz (ik - 1) + (i - 1) * ny (ik) + j 
      amax (ik) = max (amax (ik), a (ikk) ) 
      IF (a (ikk) .ne. - 9999.) amin (ik) = min (amin (ik), a (ikk) ) 
      ENDDO 
      ENDDO 
!                                                                       
      END SUBROUTINE get_extrema_z                  
!**7******************************************************************* 
      SUBROUTINE do_fmax_xy (ik, wmax, ixm, maxmax, ima) 
!+                                                                      
!     maxima bestimmen fuer xy-files                                    
!-                                                                      
      USE errlist_mod 
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxmax 
      REAL wmax (maxmax) 
      INTEGER ixm (maxmax), ik, ima 
      INTEGER ix, jf 
      LOGICAL max_da 
!                                                                       
      ima = 1 
      DO ix = 1 + ifen, lenc(ik) - ifen 
      IF (x (offxy (ik - 1) + ix) .ge.ex (iwin, iframe, 1) .and.x (     &
      offxy (ik - 1) + ix) .le.ex (iwin, iframe, 2) ) then              
         max_da = .true. 
         DO jf = - ifen+1, - 1 
         max_da = max_da.and. (y (offxy (ik - 1) + ix) .ge.    &
                               y (offxy (ik - 1) + ix + jf) )  &
                        .AND. (y (offxy (ik - 1) + ix) .gt.    &
                               y (offxy (ik - 1) + ix + jf-1) )                                              
         ENDDO 
         DO jf = 1, ifen -1
         max_da = max_da.and. (y (offxy (ik - 1) + ix) .ge.    &
                               y (offxy (ik - 1) + ix + jf) )  &
                        .AND. (y (offxy (ik - 1) + ix) .gt.    &
                               y (offxy (ik - 1) + ix + jf+1) )                                              
         ENDDO 
!                                                                       
         IF (max_da) then 
            wmax (ima) = y (offxy (ik - 1) + ix) 
            ixm (ima) = ix 
            ima = ima + 1 
            IF (ima.gt.maxmax) then 
               ier_num = - 21 
               ier_typ = ER_APPL 
               ima = maxmax 
               RETURN 
            ENDIF 
         ENDIF 
      ENDIF 
      ENDDO 
      ima = ima - 1 
!                                                                       
      IF (ima.eq.0) then 
         ier_num = - 22 
         ier_typ = ER_APPL 
      ENDIF 
!                                                                       
      END SUBROUTINE do_fmax_xy                     
!**7******************************************************************* 
      SUBROUTINE do_fmax_z (ik, wmax, ixm, iym, maxmax, ima) 
!+                                                                      
!     maxima bestimmen aus xyz files                                    
!-                                                                      
      USE errlist_mod 
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER maxmax 
      REAL wmax (maxmax) 
      INTEGER ik, ima, ixm (maxmax), iym (maxmax) 
      INTEGER ikk, iix, ix, iy, jf, kf 
      LOGICAL max_da 
!                                                                       
      ima = 1 
      ikk = offz (ik - 1) 
      DO ix = 1 + ifen, nx (ik) - ifen 
      DO iy = 1 + ifen, ny (ik) - ifen 
      iix = ix - 1 
      IF (z (ikk + iix * ny (ik) + iy) .eq. - 9999.0) goto 20 
      max_da = .true. 
      DO jf = - ifen, ifen 
      DO kf = - ifen, ifen 
      IF (jf.ne.0.or.kf.ne.0) then 
         max_da = max_da.and. (z (ikk + iix * ny (ik) + iy) ) .ge.z (   &
         ikk + (iix + jf) * ny (ik) + iy + kf)                          
         IF (.not.max_da) goto 20 
      ENDIF 
      ENDDO 
      ENDDO 
      IF (max_da) then 
         wmax (ima) = z (ikk + iix * ny (ik) + iy) 
         ixm (ima) = ix 
         iym (ima) = iy 
         ima = ima + 1 
         IF (ima.gt.maxmax) then 
            ier_num = - 21 
            ier_typ = ER_APPL 
            ima = maxmax 
            RETURN 
         ENDIF 
      ENDIF 
   20 CONTINUE 
      ENDDO 
      ENDDO 
      ima = ima - 1 
!                                                                       
      IF (ima.eq.0) then 
         ier_num = - 22 
         ier_typ = ER_APPL 
      ENDIF 
!                                                                       
      END SUBROUTINE do_fmax_z                      
!**7******************************************************************* 
      SUBROUTINE extract_subarray (xf, yf, zf, xmit, ymit, nsize, ik,   &
      ie)                                                               
!+                                                                      
!     schneiden eines subarrays fuer interpolationsroutinen             
!-                                                                      
      USE kuplot_config 
      USE kuplot_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      INTEGER nsize 
      REAL xf (nsize), yf (nsize), zf (nsize, nsize) 
      REAL xmit, ymit 
      INTEGER ik, ie, ip, ii, jj, ixx, iyy 
!                                                                       
      ie = 0 
      ip = nsize / 2 + 1 
      ixx = nint ( (xmit - xmin (ik) ) / (xmax (ik) - xmin (ik) )       &
      * REAL(nx (ik) - 1) ) + 1                                       
      iyy = nint ( (ymit - ymin (ik) ) / (ymax (ik) - ymin (ik) )       &
      * REAL(ny (ik) - 1) ) + 1                                       
      IF (ixx.lt.ip.or. (ixx + ip) .gt.nx (ik) .or.iyy.lt.ip.or. (iyy + &
      ip) .gt.ny (ik) ) then                                            
         ie = - 1 
         RETURN 
      ENDIF 
      DO ii = 1, nsize 
      DO jj = 1, nsize 
      xf (ii) = x (offxy (ik - 1) + ii - ip + ixx) 
      yf (jj) = y (offxy (ik - 1) + jj - ip + iyy) 
      zf (ii, jj) = z (offz (ik - 1) + (ii - ip + ixx) * ny (ik)        &
      + jj - ip + iyy)                                                  
      IF (zf (ii, jj) .eq. - 9999.) then 
         ie = - 1 
         RETURN 
      ENDIF 
      ENDDO 
      ENDDO 
!                                                                       
      END SUBROUTINE extract_subarray               
!********************************************************************   
      SUBROUTINE polin2 (x1a, x2a, ya, m, n, x1, x2, y, dy,ier) 
IMPLICIT integer(i-n)
IMPLICIT REAL (a-h, o-z)
      PARAMETER (nmax = 50, mmax = 50) 
      DIMENSION x1a (m), x2a (n), ya (m, n), yntmp (nmax), ymtmp (mmax) 
      INTEGER :: ier
      ier = 0
      DO 12 j = 1, m 
         DO 11 k = 1, n 
            yntmp (k) = ya (j, k) 
   11    END DO 
         CALL polint (x2a, yntmp, n, x2, ymtmp (j), dy, ier) 
         IF(ier/=0) RETURN
   12 END DO 
      CALL polint (x1a, ymtmp, m, x1, y, dy, ier) 
      RETURN 
      END SUBROUTINE polin2                         
!                                                                       
      SUBROUTINE polint (xa, ya, n, x, y, dy, ier) 
IMPLICIT integer(i-n)
IMPLICIT REAL (a-h, o-z)
      PARAMETER (nmax = 50) 
      DIMENSION xa (n), ya (n), c (nmax), d (nmax) 
      INTEGER :: ier
      ier = 0
      ns = 1 
      dif = abs (x - xa (1) ) 
      DO 11 i = 1, n 
         dift = abs (x - xa (i) ) 
         IF (dift.lt.dif) then 
            ns = i 
            dif = dift 
         ENDIF 
         c (i) = ya (i) 
         d (i) = ya (i) 
   11 END DO 
      y = ya (ns) 
      ns = ns - 1 
      DO 13 m = 1, n - 1 
         DO 12 i = 1, n - m 
            ho = xa (i) - x 
            hp = xa (i + m) - x 
            w = c (i + 1) - d (i) 
            den = ho - hp 
            IF (den.eq.0.) THEN
               ier = -60
               RETURN
            ENDIF
            den = w / den 
            d (i) = hp * den 
            c (i) = ho * den 
   12    END DO 
         IF (2 * ns.lt.n - m) then 
            dy = c (ns + 1) 
         ELSE 
            dy = d (ns) 
            ns = ns - 1 
         ENDIF 
         y = y + dy 
   13 END DO 
      RETURN 
      END SUBROUTINE polint                         
!                                                                       
      SUBROUTINE spline (x, y, n, yp1, ypn, y2) 
      USE kuplot_config 
IMPLICIT integer(i-n)
IMPLICIT REAL (a-h, o-z)
      PARAMETER (nmax = maxarray) 
      DIMENSION x (n), y (n), y2 (n), u (nmax) 
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
      SUBROUTINE splint (xa, ya, y2a, n, x, y, ier) 
IMPLICIT integer(i-n)
IMPLICIT REAL (a-h, o-z)
      DIMENSION xa (n), ya (n), y2a (n) 
      INTEGER :: ier
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
         ier = -60
         RETURN
      ENDIF
      a = (xa (khi) - x) / h 
      b = (x - xa (klo) ) / h 
      y = a * ya (klo) + b * ya (khi) + ( (a**3 - a) * y2a (klo)        &
      + (b**3 - b) * y2a (khi) ) * (h**2) / 6.                          
      RETURN 
      END SUBROUTINE splint                         
                                                                        
      SUBROUTINE SAVGOL (c, np, nl, nr, ld, m) 
      USE errlist_mod 
IMPLICIT integer(i-n)
IMPLICIT REAL (a-h, o-z)
      INTEGER ld, m, nl, np, nr, MMAX 
      REAL c (np) 
      PARAMETER (MMAX = 6) 
      INTEGER imj, ipj, j, k, kk, mm, indx (MMAX + 1) 
      REAL d, fac, sum, a (MMAX + 1, MMAX + 1), b (MMAX + 1) 
!                                                                       
      IF (np.lt.nl + nr +                                               &
      1.or.nl.lt.0.or.nr.lt.0.or.ld.gt.m.or.m.gt.MMAX.or.nl +           &
      nr.lt.m.or.np.le.3) then                                          
         ier_num = - 39 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
      DO ipj = 0, 2 * m 
      sum = 0. 
      IF (ipj.eq.0) sum = 1. 
      DO k = 1, nr 
      sum = sum + REAL(k) **ipj 
      ENDDO 
      DO k = 1, nl 
      sum = sum + REAL( - k) **ipj 
      ENDDO 
      mm = min (ipj, 2 * m - ipj) 
      DO imj = - mm, mm, 2 
      a (1 + (ipj + imj) / 2, 1 + (ipj - imj) / 2) = sum 
      ENDDO 
      ENDDO 
      CALL ludcmp (a, m + 1, MMAX + 1, indx, d, ier_num)
      IF(ier_num /= 0) THEN
         ier_typ = ER_NONE
         RETURN
      ENDIF 
      DO j = 1, m + 1 
      b (j) = 0. 
      ENDDO 
      b (ld+1) = 1. 
      CALL lubksb (a, m + 1, MMAX + 1, indx, b) 
      DO kk = 1, np 
      c (kk) = 0. 
      ENDDO 
      DO k = - nl, nr 
      sum = b (1) 
      fac = 1. 
      DO mm = 1, m 
      fac = fac * k 
      sum = sum + b (mm + 1) * fac 
      ENDDO 
      kk = mod (np - k, np) + 1 
      c (kk) = sum 
      ENDDO 
      RETURN 
      END SUBROUTINE SAVGOL                         
                                                                        
      SUBROUTINE LUBKSB (A, N, NP, INDX, B) 
IMPLICIT integer(i-n)
IMPLICIT REAL (a-h, o-z)
      DIMENSION A (NP, NP), INDX (N), B (N) 
      II = 0 
      DO 12 I = 1, N 
         LL = INDX (I) 
         SUM = B (LL) 
         B (LL) = B (I) 
         IF (II.NE.0) THEN 
            DO 11 J = II, I - 1 
               SUM = SUM - A (I, J) * B (J) 
   11       END DO 
         ELSEIF (SUM.NE.0.) THEN 
            II = I 
         ENDIF 
         B (I) = SUM 
   12 END DO 
      DO 14 I = N, 1, - 1 
         SUM = B (I) 
         IF (I.LT.N) THEN 
            DO 13 J = I + 1, N 
               SUM = SUM - A (I, J) * B (J) 
   13       END DO 
         ENDIF 
         B (I) = SUM / A (I, I) 
   14 END DO 
      RETURN 
      END SUBROUTINE LUBKSB                         
      SUBROUTINE LUDCMP (A, N, NP, INDX, D, ier) 
IMPLICIT integer(i-n)
IMPLICIT REAL (a-h, o-z)
      PARAMETER (NMAX = 100, TINY = 1.0E-20) 
      DIMENSION A (NP, NP), INDX (N), VV (NMAX) 
      INTEGER  , INTENT(OUT) :: ier
      imax = 0
      D = 1. 
      DO 12 I = 1, N 
         AAMAX = 0. 
         DO 11 J = 1, N 
            IF (ABS (A (I, J) ) .GT.AAMAX) AAMAX = ABS (A (I, J) ) 
   11    END DO 
         IF (AAMAX.EQ.0.) THEN !! PAUSE 'Singular matrix.' 
            ier = -61
            RETURN
         ENDIF
         VV (I) = 1. / AAMAX 
   12 END DO 
      DO 19 J = 1, N 
         IF (J.GT.1) THEN 
            DO 14 I = 1, J - 1 
               SUM = A (I, J) 
               IF (I.GT.1) THEN 
                  DO 13 K = 1, I - 1 
                     SUM = SUM - A (I, K) * A (K, J) 
   13             END DO 
                  A (I, J) = SUM 
               ENDIF 
   14       END DO 
         ENDIF 
         AAMAX = 0. 
         DO 16 I = J, N 
            SUM = A (I, J) 
            IF (J.GT.1) THEN 
               DO 15 K = 1, J - 1 
                  SUM = SUM - A (I, K) * A (K, J) 
   15          END DO 
               A (I, J) = SUM 
            ENDIF 
            DUM = VV (I) * ABS (SUM) 
            IF (DUM.GE.AAMAX) THEN 
               IMAX = I 
               AAMAX = DUM 
            ENDIF 
   16    END DO 
         IF (J.NE.IMAX) THEN 
            DO 17 K = 1, N 
               DUM = A (IMAX, K) 
               A (IMAX, K) = A (J, K) 
               A (J, K) = DUM 
   17       END DO 
            D = - D 
            VV (IMAX) = VV (J) 
         ENDIF 
         INDX (J) = IMAX 
         IF (J.NE.N) THEN 
            IF (A (J, J) .EQ.0.) A (J, J) = TINY 
            DUM = 1. / A (J, J) 
            DO 18 I = J + 1, N 
               A (I, J) = A (I, J) * DUM 
   18       END DO 
         ENDIF 
   19 END DO 
      IF (A (N, N) .EQ.0.) A (N, N) = TINY 
      RETURN 
      END SUBROUTINE LUDCMP                         
