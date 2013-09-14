!*****7*****************************************************************
      SUBROUTINE theory (xx, f, df, iwert) 
!                                                                       
      USE config_mod 
      USE mixscat_mod 
!                                                                       
      IMPLICIT none 
!                                                                       
      REAL xx, f, df (maxpara) 
      INTEGER iwert, ind 
!                                                                       
      DO ind = 1, npara 
      df (ind) = 0.0 
      ENDDO 
!                                                                       
      f = p (1) 
      DO ind = 1, npara - 1 
      IF (xx.ne.0) f = f + p (ind+1) * (xx**ind) 
      ENDDO 
!                                                                       
!-------Derivatives                                                     
!                                                                       
      IF (iwert.gt.0) then 
         DO ind = 0, npara - 1 
         IF (pinc (ind+1) .ne.0) then 
            df (ind+1) = (xx**ind) 
         ENDIF 
         ENDDO 
      ENDIF 
!                                                                       
      END SUBROUTINE theory                         
!*****7*****************************************************************
!     Altered subroutine FITTE for MIXSCAT use                          
!*****7*****************************************************************
      SUBROUTINE fit_poly 
!+                                                                      
!     LS routine (FITTE)                                                
!-                                                                      
      USE param_mod 
      USE prompt_mod 
!
      USE config_mod 
      USE mixscat_mod 
!
      IMPLICIT NONE
!                                                                       
      INTEGER :: i, j, k,l, ll, l1, nl   ! Loop indices
      INTEGER :: m, n, nf
      REAL :: f
      REAL :: fend
      REAL :: g
      REAL :: h
      REAL :: s
      REAL :: sqsum
      REAL :: sum3
      REAL :: sum4
      REAL :: syi
      REAL :: ttt
      REAL :: xx
      REAL :: zalt
      REAL :: zdif
      REAL :: zl
      REAL :: zwert
      REAL df (maxpara), dl (maxpara), dz (maxpara), pl (maxpara) 
!                                                                       
      n = npara 
      m = fit_np 
      j = 0 
   10 CONTINUE 
      zwert = 0.0 
      DO k = 1, n 
      dz (k) = 0.0 
      DO l = k, n 
      cl (l, k) = 0.0 
      ENDDO 
      ENDDO 
      DO i = 1, m 
      xx = fit_dat_x (i) 
      CALL theory (xx, f, df, i) 
      s = fit_dat_y (i) - f 
      zwert = zwert + fit_dat_w (i) * s * s 
      DO k = 1, n 
      IF (abs (df (k) ) .lt.0.0000001) df (k) = 0.0 
      ttt = fit_dat_w (i) * df (k) 
      dz (k) = dz (k) - s * ttt 
      DO l = k, n 
      cl (l, k) = cl (l, k) + ttt * df (l) 
      ENDDO 
      ENDDO 
      ENDDO 
      f = 1.0 
      IF (urf) 60, 60, 40 
   40 CONTINUE 
      h = 0.0 
      DO k = 1, n 
      IF (dz (k) .ne.0.0) h = h + (dz (k) * dz (k) ) / (zwert * cl (k,  &
      k) )                                                              
      ENDDO 
      f = 1.000001 + urf * h 
   60 CONTINUE 
      zalt = zl 
      zdif = zwert - zl 
      fend = f 
      IF (j) 80, 80, 70 
   70 CONTINUE 
      IF (zwert - zl) 80, 300, 300 
   80 CONTINUE 
      zl = zwert 
      DO k = 1, n 
      pl (k) = p (k) 
      ENDDO 
      DO k = 1, n 
      DO l = 1, k 
      l1 = l - 1 
      s = cl (k, l) 
      IF (l - k) 100, 140, 100 
  100 IF (l1) 130, 130, 110 
  110 DO i = 1, l1 
      s = s - cl (i, l) * cl (i, k) 
      ENDDO 
  130 cl (l, k) = s 
      ENDDO 
  140 IF (s) 145, 142, 145 
  142 dl (k) = 0.0 
      GOTO 190 
  145 s = s * f 
      IF (l1) 170, 170, 150 
  150 DO i = 1, l1 
      ttt = cl (i, k) 
      cl (i, k) = ttt * dl (i) 
      s = s - ttt * cl (i, k) 
      ENDDO 
  170 IF (s) 190, 190, 180 
  180 dl (k) = 1.0 / s 
  190 CONTINUE 
      ENDDO 
      IF (j - ncycle) 200, 200, 300 
  200 j = j + 1 
      IF (n - 1) 230, 230, 210 
  210 DO l = 2, n 
      l1 = l - 1 
      DO k = 1, l1 
      dz (l) = dz (l) - cl (k, l) * dz (k) 
      ENDDO 
      ENDDO 
  230 dz (n) = dz (n) * dl (n) 
      IF (n - 1) 260, 260, 240 
  240 DO nl = 2, n 
      l = n - nl + 1 
      l1 = l + 1 
      dz (l) = dz (l) * dl (l) 
      DO k = l1, n 
      dz (l) = dz (l) - cl (l, k) * dz (k) 
      ENDDO 
      ENDDO 
  260 DO k = 1, n 
      pinc (k) = 0.01 * dz (k) 
      p (k) = p (k) - dz (k) 
      ENDDO 
!                                                                       
      GOTO 10 
!                                                                       
  300 g = 1.0 
      h = float (m - n) 
      IF (h * zl.gt.0.0) g = zl / h 
      DO k = 1, n 
      dl (k) = dl (k) * g 
      cl (k, k) = 1.0 
      ENDDO 
      IF (n - 1) 350, 350, 320 
  320 DO l = 2, n 
      l1 = l - 1 
      DO k = 1, l1 
      s = 0.0 
      DO i = k, l1 
      s = s - cl (i, l) * cl (i, k) 
      ENDDO 
      cl (l, k) = s 
      ENDDO 
      ENDDO 
  350 DO k = 1, n 
      DO l = k, n 
      s = 0.0 
      IF (l - k) 380, 360, 380 
  360 DO i = l, n 
      ttt = cl (i, k) 
      cl (i, k) = ttt * dl (i) 
      s = s + ttt * cl (i, k) 
      ENDDO 
      GOTO 400 
  380 DO i = l, n 
      s = s + cl (i, l) * cl (i, k) 
      ENDDO 
  400 CONTINUE 
      cl (l, k) = s 
      ENDDO 
      ENDDO 
      DO l = 1, n 
      IF (dl (l) ) 502, 501, 502 
  501 dp (l) = 1.0 
      dp (l) = 0.0 
      GOTO 570 
  502 f = sqrt (cl (l, l) ) 
      dp (l) = f 
      DO k = 1, l 
      dz (k) = cl (l, k) / f 
      ENDDO 
      DO k = l, n 
      dz (k) = cl (k, l) / f 
      ENDDO 
      DO k = 1, l 
      IF (dp (k) .eq.0.0) then 
         cl (l, k) = dz (k) 
      ELSE 
         cl (l, k) = dz (k) / dp (k) 
      ENDIF 
      ENDDO 
  570 CONTINUE 
      ENDDO 
!                                                                       
      DO k = 1, n 
      p (k) = pl (k) 
      ENDDO 
      sqsum = 0.0 
      syi = 0.0 
      sum3 = 0.0 
      sum4 = 0.0 
      nf = n 
      DO i = 1, n 
      IF (df (i) .eq.0.) nf = nf - 1 
      ENDDO 
      DO i = 1, m 
      xx = fit_dat_x (i) 
      CALL theory (xx, f, df, - i) 
      h = fit_dat_y (i) - f 
      sum3 = sum3 + fit_dat_w (i) * h * h 
      sum4 = sum4 + fit_dat_w (i) * fit_dat_y (i) * fit_dat_y (i) 
      ENDDO 
      r4 = sqrt (sum3 / sum4) 
      rexp = sqrt ( (m - nf) / sum4) 
!                                                                       
      END SUBROUTINE fit_poly                       
