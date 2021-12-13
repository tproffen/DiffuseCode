MODULE metric_mod
!
CONTAINS
!*****7*****************************************************************
!
REAL function do_blen (lspace, u, v) 
!-                                                                      
!     Calculates the length of the vector v-u in the space defined by   
!     lspace (.true. = real space , .false. = reciprocal space)         
!+                                                                      
USE discus_config_mod 
USE crystal_mod 
use precision_mod
!
IMPLICIT none 
!                                                                       
LOGICAL , intent(in) :: lspace 
real(kind=PREC_SP), dimension(3), intent(in) ::  u
real(kind=PREC_SP), dimension(3), intent(in) ::  v
!                                                                       
INTEGER :: i 
REAL(kind=PREC_SP), dimension(3) :: w
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
!
!*****7*****************************************************************
!
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
use precision_mod
!
IMPLICIT none 
!                                                                       
LOGICAL , intent(in) :: lspace 
real(kind=PREC_SP), dimension(3), intent(in) ::  u
real(kind=PREC_SP), dimension(3), intent(in) ::  v
real(kind=PREC_SP), dimension(3), intent(in) ::  w
!                                                                       
!     LOGICAL lspace 
!     REAL u (3), v (3), w (3) 
!                                                                       
INTEGER i 
REAL(kind=PREC_SP) ::  xx, xy, yy, arg 
REAL(kind=PREC_SP), dimension(3) ::  x (3), y (3) 
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
!
!*****7*********************************************************        
!
REAL function skalpro (h, k, rten) 
!-                                                                      
!           Calulates the SCALARPRODUCT of two vectors                  
!           1/D**2 = H(I)*K(J)*RTEN(I,J)                                
!-                                                                      
use precision_mod
!
IMPLICIT none 
!                                                                       
real(kind=PREC_SP), dimension(3),   intent(in) :: h
real(kind=PREC_SP), dimension(3),   intent(in) :: k
real(kind=PREC_SP), dimension(3,3), intent(in) :: rten
!
INTEGER , PARAMETER :: idim = 3
!     REAL h (idim), k (idim), rten (idim, idim) 
!                                                                       
INTEGER :: i, j 
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
use precision_mod
!
IMPLICIT none 
!                                                                       
real(kind=PREC_SP), dimension(3),     intent(in)  :: u
real(kind=PREC_SP), dimension(3),     intent(in)  :: v
real(kind=PREC_SP), dimension(3),     intent(out) :: ww
real(kind=PREC_SP), dimension(3,3,3), intent(in)  :: eps
real(kind=PREC_SP), dimension(3,3),   intent(in)  :: rten
!
INTEGER :: i, j, k, l 
!REAL    :: U (3), V (3), WW (3), EPS (3, 3, 3), RTEN (3, 3) 
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
CHARACTER(LEN=*), INTENT(INOUT) :: line 
INTEGER         , INTENT(INOUT) :: laenge 
LOGICAL         , INTENT(IN)    :: lspace 
      CHARACTER(LEN=PREC_STRING) cpara (maxw)
      INTEGER lpara (maxw) 
      INTEGER ianz 
      INTEGER i
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
CHARACTER(LEN=*), INTENT(INOUT) :: line
INTEGER         , INTENT(INOUT) :: laenge 
!
      CHARACTER(LEN=MAX(PREC_STRING,LEN(line))) :: cpara (maxw)
      INTEGER lpara (maxw) 
      INTEGER ianz 
      INTEGER i
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
CHARACTER(LEN=*), INTENT(INOUT) :: line
INTEGER         , INTENT(INOUT) :: laenge 
!
      CHARACTER(LEN=MAX(PREC_STRING,LEN(line))) cpara (maxw)
      INTEGER lpara (maxw) 
      INTEGER ianz 
      INTEGER i
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
!
!*******************************************************************************
!
real function do_area(lspace, u, v)
!-
!  determine the area spanned by vectors u and v
!  area = |u| * |v| * sin(<(u,v))
!+
!
use precision_mod
use trig_degree_mod
!
implicit none
!
logical, intent(in) :: lspace    ! True=direct; false=reciprocal
real(kind=PREC_SP), dimension(3), intent(in) :: u
real(kind=PREC_SP), dimension(3), intent(in) :: v
!
real(kind=PREC_SP), dimension(3), parameter :: NULLV = (/0.0, 0.0, 0.0/)
!
do_area = do_blen(lspace,u, NULLV) * do_blen(lspace, v, NULLV) * &
          sind(do_bang(lspace, u, NULLV, v))
!
end function do_area
!
!*******************************************************************************
!
real function do_volume(lspace, u, v, w)
!-
!  determine the volume spanned by the vectors u,v,w
!+
use precision_mod
use trig_degree_mod
!
implicit none
!
logical, intent(in) :: lspace    ! True=direct; false=reciprocal
real(kind=PREC_SP), dimension(3), intent(in) :: u
real(kind=PREC_SP), dimension(3), intent(in) :: v
real(kind=PREC_SP), dimension(3), intent(in) :: w
!
real(kind=PREC_SP), dimension(3), parameter :: NULLV = (/0.0, 0.0, 0.0/)
real(kind=PREC_SP), dimension(3) :: lengths
real(kind=PREC_SP), dimension(3) :: cosines
!
lengths(1) = do_blen(lspace, u, NULLV)
lengths(2) = do_blen(lspace, v, NULLV)
lengths(3) = do_blen(lspace, w, NULLV)
cosines(1) = cosd(do_bang(lspace, v, NULLV, w))
cosines(2) = cosd(do_bang(lspace, u, NULLV, w))
cosines(3) = cosd(do_bang(lspace, u, NULLV, v))
!
do_volume = lengths(1)*lengths(2)*lengths(3) *                                  &
            sqrt(1.0E0 - cosines(1)**2 - cosines(2)**2 - cosines(3)**2 +        &
                 2.0D0*cosines(1)*cosines(2)*cosines(3)                 )
end function do_volume
!
!*******************************************************************************
!
END MODULE metric_mod
