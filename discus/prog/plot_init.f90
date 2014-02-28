MODULE plot_init_mod
!
CONTAINS
!*****7*****************************************************************
      SUBROUTINE plot_ini_trans (azero) 
!-                                                                      
!     Initializes the transformation matrix to cartesian coordinates    
!     with unit cell length 'azero', which should be "1" for XBS and    
!     DIAMOND.                                                          
!+                                                                      
      USE config_mod 
      USE crystal_mod 
      USE metric_mod
      USE plot_mod 
      USE tensors_mod
      USE trafo_mod
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER i, j 
      LOGICAL lspace 
      REAL azero 
      REAL dwert 
      REAL u (3), v (3), w (3), null (3) 
!                                                                       
!     REAL do_blen 
!                                                                       
      DATA null / 0.0, 0.0, 0.0 / 
!                                                                       
      lspace = .true. 
!                                                                       
!     cartesian b-axis is parallel to b, length = 1                     
!                                                                       
      v (1) = 0.0 
      v (2) = 1.0 
      v (3) = 0.0 
      dwert = do_blen (lspace, v, null) 
      v (2) = v (2) / dwert * azero 
      pl_tran_g (2, 1) = 0.0 
      pl_tran_g (2, 2) = v (2) 
      pl_tran_g (2, 3) = 0.0 
!                                                                       
!     cartesian c-axis is parallel to c*, length = 1                    
!                                                                       
      u (1) = 0.0 
      u (2) = 0.0 
      u (3) = 1.0 
      CALL trans (u, cr_rten, w, 3) 
      dwert = do_blen (lspace, w, null) 
      w (1) = w (1) / dwert * azero 
      w (2) = w (2) / dwert * azero 
      w (3) = w (3) / dwert * azero 
      pl_tran_g (3, 1) = w (1) 
      pl_tran_g (3, 2) = w (2) 
      pl_tran_g (3, 3) = w (3) 
!                                                                       
!     cartesian a-axis is parallel to vector product (b X c*)           
!                                                                       
      CALL vekprod (v, w, u, cr_eps, cr_rten) 
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
      pl_tran_g (4, i) = 0.0 
      pl_tran_g (i, 4) = 0.0 
      pl_tran_f (4, i) = 0.0 
      pl_tran_f (i, 4) = 0.0 
      ENDDO 
      pl_tran_f (4, 4) = 1.0 
!                                                                       
      DO i = 1, 4 
      DO j = 1, 4 
      pl_tran_gi (i, j) = pl_tran_g (i, j) 
      ENDDO 
      ENDDO 
!                                                                       
      CALL invmat4 (pl_tran_gi) 
!                                                                       
      DO i = 1, 3 
      DO j = 1, 3 
      pl_tran_f (j, i) = pl_tran_gi (i, j) 
      pl_tran_fi (j, i) = pl_tran_gi (i, j) 
      ENDDO 
      ENDDO 
      CALL invmat4 (pl_tran_fi) 
!DBG      write (*,3010) ((pl_tran_g (i,j),j=1,3),i=1,3)                
!DBG      write (*,3020) ((pl_tran_gi(i,j),j=1,3),i=1,3)                
!DBG      write (*,3030) ((pl_tran_f (i,j),j=1,4),i=1,3)                
!DBG      write (*,3040) ((pl_tran_fi(i,j),j=1,4),i=1,3)                
!                                                                       
! 3010 FORMAT    (                                                       &
!     &           ' ( a(new) ) = ( ',2(F9.5,','),f9.5,' )   ( a(old) )'/ &
!     &           ' ( b(new) ) = ( ',2(F9.5,','),f9.5,' ) * ( b(old) )'/ &
!     &           ' ( c(new) ) = ( ',2(F9.5,','),f9.5,' )   ( c(old) )'/)
! 3020 FORMAT    (                                                       &
!     &           ' ( a(old) ) = ( ',2(F9.5,','),f9.5,' )   ( a(new) )'/ &
!     &           ' ( b(old) ) = ( ',2(F9.5,','),f9.5,' ) * ( b(new) )'/ &
!     &           ' ( c(old) ) = ( ',2(F9.5,','),f9.5,' )   ( c(new) )'/)
! 3030 FORMAT    (                                                       &
!     &           ' ( x(new) ) = ( ',2(F9.5,','),f9.5,' )   ( x(old) )', &
!     &           '   (',f9.5,')'/                                       &
!     &           ' ( y(new) ) = ( ',2(F9.5,','),f9.5,' ) * ( y(old) )', &
!     &           ' + (',f9.5,')'/                                       &
!     &           ' ( z(new) ) = ( ',2(F9.5,','),f9.5,' )   ( z(old) )', &
!     &           '   (',f9.5,')'/)                                      
! 3040 FORMAT    (                                                       &
!     &           ' ( x(old) ) = ( ',2(F9.5,','),f9.5,' )   ( x(new) )', &
!     &           '   (',f9.5,')'/                                       &
!     &           ' ( y(old) ) = ( ',2(F9.5,','),f9.5,' ) * ( y(new) )', &
!     &           ' + (',f9.5,')'/                                       &
!     &           ' ( z(old) ) = ( ',2(F9.5,','),f9.5,' )   ( z(new) )', &
!     &           '   (',f9.5,')'/)                                      
!                                                                       
      END SUBROUTINE plot_ini_trans                 
END MODULE plot_init_mod
