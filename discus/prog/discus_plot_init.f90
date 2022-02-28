MODULE discus_plot_init_mod
!
CONTAINS
!
!*****7*****************************************************************
!
SUBROUTINE plot_ini_trans (azero,                        &
           pl_tran_g, pl_tran_gi, pl_tran_f, pl_tran_fi, &
           cr_gten, cr_rten, cr_eps)
!-                                                                      
!     Initializes the transformation matrix to cartesian coordinates    
!     with unit cell length 'azero', which should be "1" for XBS and    
!     DIAMOND.                                                          
!+                                                                      
USE metric_mod
USE errlist_mod 
use matrix_mod
use precision_mod
!
IMPLICIT none 
!                                                                       
REAL(kind=PREC_DP),                   INTENT(IN ) :: azero 
REAL(kind=PREC_DP), DIMENSION(4,4)  , INTENT(OUT) :: pl_tran_g
REAL(kind=PREC_DP), DIMENSION(4,4)  , INTENT(OUT) :: pl_tran_gi
REAL(kind=PREC_DP), DIMENSION(4,4)  , INTENT(OUT) :: pl_tran_f
REAL(kind=PREC_DP), DIMENSION(4,4)  , INTENT(OUT) :: pl_tran_fi
REAL(kind=PREC_DP), DIMENSION(3,3)  , INTENT(IN ) :: cr_gten
REAL(kind=PREC_DP), DIMENSION(3,3)  , INTENT(IN ) :: cr_rten
REAL(kind=PREC_DP), DIMENSION(3,3,3), INTENT(IN ) :: cr_eps
!                                                                       
INTEGER               :: i
REAL(kind=PREC_DP)    :: dwert 
REAL(kind=PREC_DP)   , DIMENSION(3) :: u, v ,w !, null 
!                                                                       
!     cartesian b-axis is parallel to b, length = 1                     
!                                                                       
v (1) = 0.0D0 
v (2) = 1.0D0 
v (3) = 0.0D0 
dwert = sqrt (skalpro (v, v, cr_gten) )
v (2) = v (2) / dwert * azero 
pl_tran_g (2, 1) = 0.0D0 
pl_tran_g (2, 2) = v (2) 
pl_tran_g (2, 3) = 0.0D0 
!                                                                       
!     cartesian c-axis is parallel to c*, length = 1                    
!                                                                       
u (1) = 0.0 
u (2) = 0.0 
u (3) = 1.0 
!     CALL trans (u, cr_rten, w, 3) 
w = matmul(cr_rten, u)
dwert = sqrt (skalpro (w, w, cr_gten) )
w (1) = w (1) / dwert * azero 
w (2) = w (2) / dwert * azero 
w (3) = w (3) / dwert * azero 
pl_tran_g (3, 1) = w (1) 
pl_tran_g (3, 2) = w (2) 
pl_tran_g (3, 3) = w (3) 
!                                                                       
!     cartesian a-axis is parallel to vector product (b X c*)           
!                                                                       
call vekprod (v, w, u, cr_eps, cr_rten) 
pl_tran_g (1, 1) = u (1) / azero 
pl_tran_g (1, 2) = u (2) / azero 
pl_tran_g (1, 3) = u (3) / azero 
u (1) = u (1) / azero 
u (2) = u (2) / azero 
u (3) = u (3) / azero 
!                                                                       
!     calculate matrix for atom transformation                          
!                                                                       
DO i = 1, 4 
   pl_tran_g (4, i) = 0.0D0 
   pl_tran_g (i, 4) = 0.0D0 
   pl_tran_f (4, i) = 0.0D0 
   pl_tran_f (i, 4) = 0.0D0 
ENDDO 
pl_tran_f (4, 4) = 1.0D0 
pl_tran_g (4, 4) = 1.0D0 
!                                                                       
!DO i = 1, 4 
!   DO j = 1, 4 
!      pl_tran_gi (i, j) = pl_tran_g (i, j) 
!   ENDDO 
!ENDDO 
!                                                                       
!     CALL invmat4 (pl_tran_gi) 
call matinv(pl_tran_g , pl_tran_gi)
pl_tran_fi = transpose(pl_tran_g)
pl_tran_f  = transpose(pl_tran_gi)
!                                                                       
!     DO i = 1, 3 
!        DO j = 1, 3 
!           pl_tran_f (j, i) = pl_tran_gi (i, j) 
!           pl_tran_fi (j, i) = pl_tran_gi (i, j) 
!        ENDDO 
!     ENDDO 
!     CALL invmat4 (pl_tran_fi) 
!     call matinv(pl_tran_f , pl_tran_fi)
!     write (*,3010) ((pl_tran_g (i,j),j=1,3),i=1,3)                
!     write (*,3020) ((pl_tran_gi(i,j),j=1,3),i=1,3)                
!     write (*,3030) ((pl_tran_f (i,j),j=1,4),i=1,3)                
!     write (*,3040) ((pl_tran_fi(i,j),j=1,4),i=1,3)                
!                                                                       
!3010 FORMAT    (                                                       &
!    &           ' ( a(new) ) = ( ',2(F9.5,','),f9.5,' )   ( a(old) )'/ &
!    &           ' ( b(new) ) = ( ',2(F9.5,','),f9.5,' ) * ( b(old) )'/ &
!    &           ' ( c(new) ) = ( ',2(F9.5,','),f9.5,' )   ( c(old) )'/)
!3020 FORMAT    (                                                       &
!    &           ' ( a(old) ) = ( ',2(F9.5,','),f9.5,' )   ( a(new) )'/ &
!    &           ' ( b(old) ) = ( ',2(F9.5,','),f9.5,' ) * ( b(new) )'/ &
!    &           ' ( c(old) ) = ( ',2(F9.5,','),f9.5,' )   ( c(new) )'/)
!3030 FORMAT    (                                                       &
!    &           ' ( x(new) ) = ( ',2(F9.5,','),f9.5,' )   ( x(old) )', &
!    &           '   (',f9.5,')'/                                       &
!    &           ' ( y(new) ) = ( ',2(F9.5,','),f9.5,' ) * ( y(old) )', &
!    &           ' + (',f9.5,')'/                                       &
!    &           ' ( z(new) ) = ( ',2(F9.5,','),f9.5,' )   ( z(old) )', &
!    &           '   (',f9.5,')'/)                                      
!3040 FORMAT    (                                                       &
!    &           ' ( x(old) ) = ( ',2(F9.5,','),f9.5,' )   ( x(new) )', &
!    &           '   (',f9.5,')'/                                       &
!    &           ' ( y(old) ) = ( ',2(F9.5,','),f9.5,' ) * ( y(new) )', &
!    &           ' + (',f9.5,')'/                                       &
!    &           ' ( z(old) ) = ( ',2(F9.5,','),f9.5,' )   ( z(new) )', &
!    &           '   (',f9.5,')'/)                                      
!                                                                       
END SUBROUTINE plot_ini_trans                 
!
!*****7*****************************************************************
!
SUBROUTINE plot_ini_jmol  (azero,                        &
           pl_tran_g, pl_tran_gi, pl_tran_f, pl_tran_fi, &
           cr_gten, cr_rten, cr_eps)
!-                                                                      
!     Initializes the transformation matrix to cartesian coordinates    
!     with unit cell length 'azero', which should be "1" for jmol
!     JMOL version:
! x: a
! y: c* x a
! z: c*
!+                                                                      
USE metric_mod
USE errlist_mod 
use matrix_mod
use precision_mod
!
IMPLICIT none 
!                                                                       
REAL(kind=PREC_DP),     INTENT(IN ) :: azero 
REAL(kind=PREC_DP), DIMENSION(4,4)  , INTENT(OUT) :: pl_tran_g
REAL(kind=PREC_DP), DIMENSION(4,4)  , INTENT(OUT) :: pl_tran_gi
REAL(kind=PREC_DP), DIMENSION(4,4)  , INTENT(OUT) :: pl_tran_f
REAL(kind=PREC_DP), DIMENSION(4,4)  , INTENT(OUT) :: pl_tran_fi
REAL(kind=PREC_DP), DIMENSION(3,3)  , INTENT(IN ) :: cr_gten
REAL(kind=PREC_DP), DIMENSION(3,3)  , INTENT(IN ) :: cr_rten
REAL(kind=PREC_DP), DIMENSION(3,3,3), INTENT(IN ) :: cr_eps
!                                                                       
INTEGER               :: i
REAL(kind=PREC_DP)    :: dwert 
REAL(kind=PREC_DP)   , DIMENSION(3) :: u, v ,w !, null 
!                                                                       
!     cartesian a-axis is parallel to a, length = 1                     
!                                                                       
u (1) = 1.0D0 
u (2) = 0.0D0 
u (3) = 0.0D0 
dwert = sqrt (skalpro (u, u, cr_gten) )
u (1) = u (1) / dwert * azero 
pl_tran_g (1, 1) = u (1)
pl_tran_g (1, 2) = 0.0D0
pl_tran_g (1, 3) = 0.0D0 
!                                                                       
!     cartesian c-axis is parallel to c*, length = 1                    
!                                                                       
v (1) = 0.0D0 
v (2) = 0.0D0 
v (3) = 1.0D0 
!CALL trans (v, cr_rten, w, 3) 
w = matmul(cr_rten, v)
dwert = sqrt (skalpro (w, w, cr_gten) )
w (1) = w (1) / dwert * azero 
w (2) = w (2) / dwert * azero 
w (3) = w (3) / dwert * azero 
pl_tran_g (3, 1) = w (1) 
pl_tran_g (3, 2) = w (2) 
pl_tran_g (3, 3) = w (3) 
!                                                                       
!     cartesian b-axis is parallel to vector product (c* X a)           
!                                                                       
CALL vekprod (w, u, v, cr_eps, cr_rten) 
pl_tran_g (2, 1) = v (1) / azero 
pl_tran_g (2, 2) = v (2) / azero 
pl_tran_g (2, 3) = v (3) / azero 
v (1) = v (1) / azero 
v (2) = v (2) / azero 
v (3) = v (3) / azero 
!                                                                       
!     calculate matrix for atom transformation                          
!                                                                       
DO i = 1, 4 
   pl_tran_g (4, i) = 0.0D0 
   pl_tran_g (i, 4) = 0.0D0 
   pl_tran_f (4, i) = 0.0D0 
   pl_tran_f (i, 4) = 0.0D0 
ENDDO 
pl_tran_g (4, 4) = 1.0D0 
pl_tran_f (4, 4) = 1.0D0 
!                                                                       
!DO i = 1, 4 
!   DO j = 1, 4 
!      pl_tran_gi (i, j) = pl_tran_g (i, j) 
!   ENDDO 
!ENDDO 
!                                                                       
!ALL invmat4 (pl_tran_gi) 
call matinv(pl_tran_g , pl_tran_gi)
pl_tran_fi = transpose(pl_tran_g)
pl_tran_f  = transpose(pl_tran_gi)
!                                                                       
!DO i = 1, 3 
!   DO j = 1, 3 
!      pl_tran_f (j, i) = pl_tran_gi (i, j) 
!      pl_tran_fi (j, i) = pl_tran_gi (i, j) 
!   ENDDO 
!ENDDO 
!CALL invmat4 (pl_tran_fi) 
!CALL matinv (pl_tran_f , pl_tran_fi) 
!     write (*,3010) ((pl_tran_g (i,j),j=1,3),i=1,3)                
!     write (*,3020) ((pl_tran_gi(i,j),j=1,3),i=1,3)                
!     write (*,3030) ((pl_tran_f (i,j),j=1,4),i=1,3)                
!     write (*,3040) ((pl_tran_fi(i,j),j=1,4),i=1,3)                
!                                                                       
!3010 FORMAT    (                                                       &
!    &           ' ( a(new) ) = ( ',2(F9.5,','),f9.5,' )   ( a(old) )'/ &
!    &           ' ( b(new) ) = ( ',2(F9.5,','),f9.5,' ) * ( b(old) )'/ &
!    &           ' ( c(new) ) = ( ',2(F9.5,','),f9.5,' )   ( c(old) )'/)
!3020 FORMAT    (                                                       &
!    &           ' ( a(old) ) = ( ',2(F9.5,','),f9.5,' )   ( a(new) )'/ &
!    &           ' ( b(old) ) = ( ',2(F9.5,','),f9.5,' ) * ( b(new) )'/ &
!    &           ' ( c(old) ) = ( ',2(F9.5,','),f9.5,' )   ( c(new) )'/)
!3030 FORMAT    (                                                       &
!    &           ' ( x(new) ) = ( ',2(F9.5,','),f9.5,' )   ( x(old) )', &
!    &           '   (',f9.5,')'/                                       &
!    &           ' ( y(new) ) = ( ',2(F9.5,','),f9.5,' ) * ( y(old) )', &
!    &           ' + (',f9.5,')'/                                       &
!    &           ' ( z(new) ) = ( ',2(F9.5,','),f9.5,' )   ( z(old) )', &
!    &           '   (',f9.5,')'/)                                      
!3040 FORMAT    (                                                       &
!    &           ' ( x(old) ) = ( ',2(F9.5,','),f9.5,' )   ( x(new) )', &
!    &           '   (',f9.5,')'/                                       &
!    &           ' ( y(old) ) = ( ',2(F9.5,','),f9.5,' ) * ( y(new) )', &
!    &           ' + (',f9.5,')'/                                       &
!    &           ' ( z(old) ) = ( ',2(F9.5,','),f9.5,' )   ( z(new) )', &
!    &           '   (',f9.5,')'/)                                      
!                                                                       
END SUBROUTINE plot_ini_jmol
!
END MODULE discus_plot_init_mod
