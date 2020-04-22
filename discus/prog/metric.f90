MODULE metric_mod
!
CONTAINS
!*****7*****************************************************************
      REAL function do_blen (lspace, u, v) 
!-                                                                      
!     Calculates the length of the vector v-u in the space defined by   
!     lspace (.true. = real space , .false. = reciprocal space)         
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      LOGICAL lspace 
      REAL u (3), v (3) 
!                                                                       
      INTEGER i 
      REAL w (3) 
!     REAL skalpro 
!                                                                       
      DO i = 1, 3 
      w (i) = v (i) - u (i) 
      ENDDO 
!                                                                       
      IF (lspace) then 
         do_blen = sqrt (skalpro (w, w, cr_gten) ) 
      ELSE 
         do_blen = sqrt (skalpro (w, w, cr_rten) ) 
      ENDIF 
!                                                                       
      END FUNCTION do_blen                          
!*****7*****************************************************************
      REAL function do_bang (lspace, u, v, w) 
!-                                                                      
!     Calculates the angle between vectors v-u and v-w in the space     
!     defined by lspace (.true. = real space ,                          
!       .false. = reciprocal space)                                     
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE errlist_mod 
      USE wink_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      LOGICAL lspace 
      REAL u (3), v (3), w (3) 
!                                                                       
      INTEGER i 
      REAL xx, xy, yy, arg 
      REAL x (3), y (3) 
!     REAL skalpro 
!                                                                       
      do_bang = 0.0 
!                                                                       
      DO i = 1, 3 
      x (i) = v (i) - u (i) 
      y (i) = v (i) - w (i) 
      ENDDO 
!                                                                       
      IF (lspace) then 
         xx = sqrt (skalpro (x, x, cr_gten) ) 
         xy = skalpro (x, y, cr_gten) 
         yy = sqrt (skalpro (y, y, cr_gten) ) 
      ELSE 
         xx = sqrt (skalpro (x, x, cr_rten) ) 
         xy = skalpro (x, y, cr_rten) 
         yy = sqrt (skalpro (y, y, cr_rten) ) 
      ENDIF 
      IF (xx.gt.0.and.yy.gt.0) then 
         arg = xy / (xx * yy) 
         IF (abs (arg) .lt.1.0) then 
            do_bang = acos (arg) / REAL(rad )
         ELSEIF (arg.ge.1.0) then 
            do_bang = 0.0 
         ELSEIF (arg.le. - 1.0) then 
            do_bang = 180.0 
         ENDIF 
      ELSE 
         ier_num = - 32 
         ier_typ = ER_APPL 
         do_bang = 0.0 
      ENDIF 
!                                                                       
      END FUNCTION do_bang                          
!*****7*********************************************************        
      REAL function skalpro (h, k, rten) 
!-                                                                      
!           Calulates the SCALARPRODUCT of two vectors                  
!           1/D**2 = H(I)*K(J)*RTEN(I,J)                                
!-                                                                      
      IMPLICIT none 
!                                                                       
      INTEGER idim 
      PARAMETER (idim = 3) 
      REAL h (idim), k (idim), rten (idim, idim) 
!                                                                       
      INTEGER i, j 
!                                                                       
      skalpro = 0.0 
      DO i = 1, idim 
      DO j = 1, idim 
      skalpro = skalpro + h (i) * k (j) * rten (i, j) 
      ENDDO 
      ENDDO 
!                                                                       
      END FUNCTION skalpro                          
!*****7*****************************************************************
!
SUBROUTINE vekprod (u, v, ww, eps, rten) 
!-                                                                      
!     calculates the VECTORPRODUCT in general triclinic space           
!     with  EPS and RTEN in direct space                                
!     with REPS and GTEN in reciprocal space                            
!+                                                                      
IMPLICIT none 
!                                                                       
INTEGER :: i, j, k, l 
REAL    :: U (3), V (3), WW (3), EPS (3, 3, 3), RTEN (3, 3) 
!                                                                       
DO i = 1, 3 
   ww (i) = 0.0 
   DO j = 1, 3 
      DO k = 1, 3 
         DO l = 1, 3 
            ww (i) = ww (i) + eps (j, k, l) * u (k) * v (l) * rten (j, i) 
         ENDDO 
      ENDDO 
   ENDDO 
ENDDO 
!                                                                       
END SUBROUTINE vekprod                        
!
!*****7*****************************************************************
      SUBROUTINE d2r (line, laenge, lspace) 
!-                                                                      
!     Transforms a vector from direct to reciprocal space or vice versa 
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE trafo_mod
!
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
      USE param_mod 
USE precision_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 3) 
!                                                                       
      CHARACTER(1024) cpara (maxw), line 
      INTEGER lpara (maxw) 
      INTEGER ianz 
      INTEGER i, laenge 
      LOGICAL lspace 
      REAL(KIND=PREC_DP) :: werte (maxw) 
      REAL(KIND=PREC_SP) :: uu, vv 
      REAL(KIND=PREC_SP) :: u (3), v (3) 
!                                                                       
!     REAL skalpro 
!                                                                       
      IF (cr_v.le.0.0) then 
         ier_num = - 35 
         ier_typ = ER_APPL 
         ier_msg (1) = 'A proper unit cell must be defined' 
      ier_msg (2)  = 'for this command to operate       ' 
         RETURN 
      ENDIF 
!                                                                       
      CALL get_params (line, ianz, cpara, lpara, maxw, laenge) 
      IF (ier_num.eq.0) then 
         IF (ianz.eq.3) then 
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
            IF (ier_num.eq.0) then 
               DO i = 1, 3 
               u (i) = werte (i) 
               ENDDO 
               IF (lspace) then 
                  CALL trans (u, cr_gten, v, 3) 
                  uu = sqrt (skalpro (u, u, cr_gten) ) 
                  vv = sqrt (skalpro (v, v, cr_rten) ) 
               ELSE 
                  CALL trans (u, cr_rten, v, 3) 
                  uu = sqrt (skalpro (u, u, cr_rten) ) 
                  vv = sqrt (skalpro (v, v, cr_gten) ) 
               ENDIF 
               IF (abs (uu) .lt.1e-8.or.abs (vv) .lt.1e-8) then 
                  ier_num = - 32 
                  ier_typ = ER_APPL 
                  RETURN 
               ENDIF 
!                                                                       
!     ----Write resulting vector into res_para                          
!                                                                       
               res_para (0) = 6 
               DO i = 1, 3 
               res_para (i) = v (i) 
               ENDDO 
               DO i = 1, 3 
               res_para (i + 3) = v (i) / uu / vv 
               ENDDO 
            ENDIF 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
      ENDIF 
      END SUBROUTINE d2r                            
!*****7*****************************************************************
      SUBROUTINE vprod (line, laenge) 
!-                                                                      
!     Calculates the vector product v X u                               
!     both input vectors and the result can be in direct or             
!     reciprocal spce                                                   
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE trafo_mod
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
      USE param_mod 
USE precision_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 7) 
!                                                                       
      CHARACTER(1024) cpara (maxw), line 
      INTEGER lpara (maxw) 
      INTEGER ianz 
      INTEGER i, laenge 
      REAL(KIND=PREC_DP) :: werte (maxw) 
      REAL(KIND=PREC_SP) :: u (3), v (3), w (3) 
!                                                                       
      IF (cr_v.le.0.0) then 
         ier_num = - 35 
         ier_typ = ER_APPL 
         ier_msg (1) = 'A proper unit cell must be defined' 
      ier_msg (2)  = 'for this command to operate       ' 
         RETURN 
      ENDIF 
!                                                                       
      CALL get_params (line, ianz, cpara, lpara, maxw, laenge) 
      IF (ier_num.eq.0) then 
         IF (ianz.ge.6) then 
!                                                                       
!     ----Set default to direct space                                   
!                                                                       
            IF (ianz.eq.6) then 
               cpara (7) = 'ddd' 
               lpara (7) = 3 
            ENDIF 
            ianz = 6 
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
            IF (ier_num.eq.0) then 
               DO i = 1, 3 
               u (i) = werte (i) 
               v (i) = werte (i + 3) 
               ENDDO 
               IF (cpara (7) (1:1) .eq.'r') then 
!                                                                       
!     ------------Transform vector 1 into direct space                  
!                                                                       
                  CALL trans (u, cr_rten, w, 3) 
                  DO i = 1, 3 
                  u (i) = w (i) 
                  ENDDO 
               ENDIF 
               IF (cpara (7) (2:2) .eq.'r') then 
!                                                                       
!     ------------Transform vector 2 into direct space                  
!                                                                       
                  CALL trans (v, cr_rten, w, 3) 
                  DO i = 1, 3 
                  v (i) = w (i) 
                  ENDDO 
               ENDIF 
!                                                                       
!     ------Calculate vektorprodukt, always in direct space             
!                                                                       
               CALL vekprod (u, v, w, cr_eps, cr_rten) 
!                                                                       
               IF (cpara (7) (3:3) .eq.'r') then 
!                                                                       
!     ------------Transform result vector into reciprocal space         
!                                                                       
                  CALL trans (w, cr_gten, u, 3) 
                  DO i = 1, 3 
                  w (i) = u (i) 
                  ENDDO 
               ENDIF 
!                                                                       
!     ----------Write result vector into res_para                       
!                                                                       
               res_para (0) = 3 
               DO i = 1, 3 
               res_para (i) = w (i) 
               ENDDO 
            ENDIF 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
      ENDIF 
!                                                                       
      END SUBROUTINE vprod                          
!*****7*****************************************************************
      SUBROUTINE do_proj (line, laenge) 
!-                                                                      
!     Calculates the projection of a vector onto another vector and     
!     onto the plane normal to the second vector.                       
!     Both input vectors and the result can be in direct or             
!     reciprocal spce                                                   
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
      USE wink_mod
      USE param_mod 
USE precision_mod
      USE trafo_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 7) 
!                                                                       
      CHARACTER(1024) cpara (maxw), line 
      INTEGER lpara (maxw) 
      INTEGER ianz 
      INTEGER i, laenge 
      LOGICAL lflag 
      REAL(KIND=PREC_DP) :: werte (maxw) 
      REAL(KIND=PREC_SP) :: u (3), v (3), w (3), h (3), x (3) 
      REAL(KIND=PREC_SP) :: vv, uv, uh, hh, uu 
      REAL(KIND=PREC_SP) :: arg 
!                                                                       
!     REAL skalpro 
!                                                                       
      IF (cr_v.le.0.0) then 
         ier_num = - 35 
         ier_typ = ER_APPL 
         ier_msg (1) = 'A proper unit cell must be defined' 
      ier_msg (2)  = 'for this command to operate       ' 
         RETURN 
      ENDIF 
!                                                                       
      lflag = .false. 
!                                                                       
      CALL get_params (line, ianz, cpara, lpara, maxw, laenge) 
      IF (ier_num.eq.0) then 
         IF (ianz.ge.6) then 
!                                                                       
!     ----Set default to direct space                                   
!                                                                       
            IF (ianz.eq.6) then 
               cpara (7) = 'dddd' 
               lpara (7) = 4 
               lflag = .true. 
            ENDIF 
            ianz = 6 
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
            IF (ier_num.eq.0) then 
               DO i = 1, 3 
               u (i) = werte (i) 
               v (i) = werte (i + 3) 
               ENDDO 
               IF (cpara (7) (1:1) .eq.'r') then 
!                                                                       
!     ------------Transform vector 1 into direct space                  
!                                                                       
                  CALL trans (u, cr_rten, w, 3) 
                  DO i = 1, 3 
                  u (i) = w (i) 
                  ENDDO 
                  lflag = .true. 
               ELSEIF (cpara (7) (1:1) .eq.'d') then 
                  lflag = .true. 
               ENDIF 
               IF (cpara (7) (2:2) .eq.'r') then 
!                                                                       
!     ------------Transform vector 2 into direct space                  
!                                                                       
                  CALL trans (v, cr_rten, w, 3) 
                  DO i = 1, 3 
                  v (i) = w (i) 
                  ENDDO 
                  lflag = .true. 
               ELSEIF (cpara (7) (2:2) .eq.'d') then 
                  lflag = .true. 
               ENDIF 
!                                                                       
!     ------Calculate projection onto second vector, always             
!              in direct space                                          
!                                                                       
               uv = skalpro (u, v, cr_gten) 
               uu = skalpro (u, u, cr_gten) 
               vv = skalpro (v, v, cr_gten) 
               IF (vv.eq.0.0) then 
                  ier_num = - 32 
                  ier_typ = ER_APPL 
                  res_para (0) = 0 
                  RETURN 
               ENDIF 
               DO i = 1, 3 
               w (i) = v (i) * uv / vv 
               ENDDO 
!                                                                       
               IF (cpara (7) (3:3) .eq.'r') then 
!                                                                       
!     ------------Transform result vector into reciprocal space         
!                                                                       
                  CALL trans (w, cr_gten, x, 3) 
                  DO i = 1, 3 
                  w (i) = x (i) 
                  ENDDO 
                  lflag = .true. 
               ELSEIF (cpara (7) (3:3) .eq.'d') then 
                  lflag = .true. 
               ENDIF 
!                                                                       
!     ----------Write result vector into res_para                       
!                                                                       
               res_para (0) = 3 
               DO i = 1, 3 
               res_para (i) = w (i) 
               ENDDO 
!                                                                       
!     ------Calculate projection onto plane normal to second vector,    
!           always in direct space                                      
!                                                                       
!                                                                       
!     ------------Transform vector 2 into reciprocal space              
!                                                                       
               CALL trans (v, cr_gten, h, 3) 
               uh = 0.0 
               DO i = 1, 3 
               uh = uh + u (i) * h (i) 
               ENDDO 
               hh = skalpro (h, h, cr_rten) 
               IF (hh.eq.0.0) then 
                  ier_num = - 32 
                  ier_typ = ER_APPL 
                  res_para (0) = 0 
                  RETURN 
               ENDIF 
               DO i = 1, 3 
               w (i) = u (i) - uh * v (i) / hh 
               ENDDO 
!                                                                       
               IF (cpara (7) (4:4) .eq.'r') then 
!                                                                       
!     ------------Transform result vector into reciprocal space         
!                                                                       
                  CALL trans (w, cr_gten, x, 3) 
                  DO i = 1, 3 
                  w (i) = x (i) 
                  ENDDO 
                  lflag = .true. 
               ELSEIF (cpara (7) (4:4) .eq.'d') then 
                  lflag = .true. 
               ENDIF 
!                                                                       
!     ----------Write result vector into res_para                       
!                                                                       
               res_para (0) = 6 
               DO i = 1, 3 
               res_para (i + 3) = w (i) 
               ENDDO 
!                                                                       
!     ------Store angle between vectors in variable res[7]              
!                                                                       
               res_para (0) = 7 
               IF (uu.eq.0.0.or.vv.eq.0.0) then 
                  ier_num = - 32 
                  ier_typ = ER_APPL 
                  res_para (0) = 0 
                  RETURN 
               ENDIF 
               arg = uv / sqrt (uu * vv) 
               IF (abs (arg) .lt.1.0) then 
                  res_para (7) = acos (arg) / REAL(rad )
               ELSEIF (arg.ge.1.0) then 
                  res_para (7) = 0.0 
               ELSEIF (arg.le. - 1.0) then 
                  res_para (7) = 180.0 
               ENDIF 
!                                                                       
               IF (.not.lflag) then 
                  ier_num = - 81 
                  ier_typ = ER_APPL 
      ier_msg (1)  = 'The space flag must be 4 characters long' 
                  ier_msg (2) = 'and consist of ''d'' or ''r'' only.' 
                  res_para (0) = 0 
               ENDIF 
            ENDIF 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
      ENDIF 
!                                                                       
      END SUBROUTINE do_proj                        
END MODULE metric_mod
