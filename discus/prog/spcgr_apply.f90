MODULE spcgr_apply
!
USE errlist_mod
!
IMPLICIT NONE
!
CONTAINS
!*****7*****************************************************************
      SUBROUTINE get_symmetry_matrices
!-                                                                      
!     Creates all symmetry matrices for the current space group         
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE generate_mod 
      USE gen_add_mod 
      USE sym_add_mod 
      USE unitcell_mod 
      USE wyckoff_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
!                                                                       
      INTEGER igs 
      INTEGER igg 
      INTEGER i, j, k 
!                                                                       
!     Reset all symmetry matrices                                       
!                                                                       
      DO k = 1, SPC_MAX 
      DO i = 1, 4 
      DO j = 1, 4 
      spc_mat (i, j, k) = 0.0 
      ENDDO 
      ENDDO 
      ENDDO 
!                                                                       
!     The first symmetry matrix is the identity matrix                  
!                                                                       
      spc_n = 1 
!                                                                       
      DO i = 1, 4 
      DO j = 1, 4 
      spc_mat (i, j, spc_n) = 0.0 
      ENDDO 
      spc_mat (i, i, spc_n) = 1.0 
      ENDDO 
      spc_spur (spc_n) = 3 
      spc_det (spc_n) = 1 
      CALL get_symmetry_type (SPC_MAX, spc_n, spc_mat, spc_spur,        &
      spc_det, spc_char, spc_xyz)
!                                                                       
!     Loop over all generators                                          
!                                                                       
      DO igs = 1, generspcgr (0, cr_spcgrno) 
      IF (gen_sta.eq.GEN_SYMM) then 
         igg = generspcgr (igs, cr_spcgrno) 
      ELSEIF (gen_sta.eq.GEN_CENTER) then 
         igg = generspcgr_center (igs, cr_spcgrno) 
      ENDIF 
      CALL make_symmetry_matrix (SPC_MAX, spc_n, spc_mat, spc_det,      &
      spc_spur, spc_char, spc_xyz, igg, NG, generators, generpower)     
      IF (ier_num.ne.0) then 
         RETURN 
      ENDIF 
      ENDDO 
!-----      End of loop over all generators                             
!                                                                       
!     Loop over all additional generators                               
!                                                                       
      DO igs = 1, gen_add_n 
      CALL make_symmetry_matrix (SPC_MAX, spc_n, spc_mat, spc_det,      &
      spc_spur, spc_char, spc_xyz, igs, GEN_ADD_MAX, gen_add,           &
      gen_add_power)                                                    
      IF (ier_num.ne.0) then 
         RETURN 
      ENDIF 
      ENDDO 
!-----      End of loop over all additional generators                  
!                                                                       
!     Loop over all additional symmetry operations, just copy these     
!                                                                       
      DO igs = 1, sym_add_n 
      spc_n = spc_n + 1 
      spc_det (spc_n) = 0 
      spc_spur (spc_n) = 0 
      DO i = 1, 4 
      DO j = 1, 4 
      DO k = 1, 4 
      spc_mat (i, j, spc_n) = sym_add (i, j, igs) 
      ENDDO 
      ENDDO 
      ENDDO 
      DO i = 1, 3 
      IF (spc_mat (i, 4, spc_n) .ge.1.0) then 
         spc_mat (i, 4, spc_n) = spc_mat (i, 4, spc_n) - int (spc_mat ( &
         i, 4, spc_n) )                                                 
      ELSEIF (spc_mat (i, 4, spc_n) .lt.0.0) then 
         spc_mat (i, 4, spc_n) = spc_mat (i, 4, spc_n) + int (spc_mat ( &
         i, 4, spc_n) ) + 1                                             
      ENDIF 
      spc_spur (spc_n) = spc_spur (spc_n) + spc_mat (i, i, spc_n) 
      ENDDO 
      spc_det (spc_n) = spc_mat (1, 1, spc_n) * (spc_mat (2, 2, spc_n)  &
      * spc_mat (3, 3, spc_n) - spc_mat (3, 2, spc_n) * spc_mat (2, 3,  &
      spc_n) ) - spc_mat (1, 2, spc_n) * (spc_mat (2, 1, spc_n) *       &
      spc_mat (3, 3, spc_n) - spc_mat (3, 1, spc_n) * spc_mat (2, 3,    &
      spc_n) ) + spc_mat (1, 3, spc_n) * (spc_mat (2, 1, spc_n) *       &
      spc_mat (3, 2, spc_n) - spc_mat (3, 1, spc_n) * spc_mat (2, 2,    &
      spc_n) )                                                          
      CALL get_symmetry_type (SPC_MAX, spc_n, spc_mat, spc_spur,        &
      spc_det, spc_char, spc_xyz)
      ENDDO 
!-----      End of loop over additional symmetry matices                
!                                                                       
      END SUBROUTINE get_symmetry_matrices         
!*****7*****************+***********************************************
      SUBROUTINE make_symmetry_matrix (SPC_MAX, spc_n, spc_mat, spc_det,&
      spc_spur, spc_char, spc_xyz, igg, NG, generators, generpower)     
!-                                                                      
!     Applies the current generator to all existing symmetry matrices   
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER SPC_MAX 
      CHARACTER(65) spc_char (1:SPC_MAX) 
      CHARACTER(87) spc_xyz (1:SPC_MAX) 
      INTEGER spc_n 
      REAL spc_mat (4, 4, 1:SPC_MAX) 
      REAL spc_det (1:SPC_MAX) 
      REAL spc_spur (1:SPC_MAX) 
      INTEGER igg 
      INTEGER NG 
      REAL generators (4, 4, 0:NG) 
      INTEGER generpower (NG) 
!                                                                       
      INTEGER i, j, k 
      INTEGER ipg 
      INTEGER imat 
      INTEGER nmat 
      LOGICAL             :: lexist = .false.
      REAL, PARAMETER     :: EPS    = 1.e-5
      REAL xmat (4, 4) 
      REAL wmat (4, 4) 
!                                                                       
!     For convenience create Identity operator                          
!                                                                       
      DO i = 1, 4 
         DO j = 1, 4 
            xmat (i, j) = 0.0 
         ENDDO 
         xmat (i, i) = 1.0 
      ENDDO 
!                                                                       
!     nmat is the number of symetry operations prior to this generator  
!     if the generator is applied not only to the power of 1, but 2 as  
!     well, it must act on the number of matrices that existed prior    
!     to ist application only.                                          
!                                                                       
      nmat = spc_n 
!                                                                       
!     Loop over all powers of generator igg                             
!                                                                       
      loop_power: DO ipg = 1, generpower (igg) 
!                                                                       
!     --raise power of generator, xmat is a dummy matrix, equal to      
!     --the previous power of the generator                             
!                                                                       
         DO i = 1, 4 
            DO j = 1, 4 
               wmat (i, j) = 0.0 
               DO k = 1, 4 
                  wmat (i, j) = wmat (i, j) + generators (i, k, igg) * xmat (k, j) 
               ENDDO 
            ENDDO 
         ENDDO 
!                                                                       
!     --apply generator to all previous symmetry matrices               
!
         generate: DO imat = 1, nmat 
            spc_n = spc_n + 1 
            spc_det (spc_n) = 0 
            spc_spur (spc_n) = 0 
            spc_mat(:,:,spc_n) = 0
            DO i = 1, 4 
               DO j = 1, 4 
                  DO k = 1, 4 
                     spc_mat (i, j, spc_n) = spc_mat (i, j, spc_n) +  &
                               wmat (i, k) * spc_mat (k, j, imat)                                            
                  ENDDO 
               ENDDO 
            ENDDO 
            DO i = 1, 3 
               IF (spc_mat (i, 4, spc_n) .ge.1.0) then 
                  spc_mat (i, 4, spc_n) = spc_mat (i, 4, spc_n) -     &
                                     int (spc_mat (i, 4, spc_n) )
               ELSEIF (spc_mat (i, 4, spc_n) .lt.0.0) then 
                  spc_mat (i, 4, spc_n) = spc_mat (i, 4, spc_n) +     &
                                     int (spc_mat (i, 4, spc_n) ) + 1
               ENDIF 
               spc_spur (spc_n) = spc_spur (spc_n) + spc_mat (i, i, spc_n) 
            ENDDO 
            spc_det (spc_n) = spc_mat (1, 1, spc_n) * (spc_mat (2, 2, spc_n) * spc_mat (3, 3, spc_n) -      &
                                                       spc_mat (3, 2, spc_n) * spc_mat (2, 3, spc_n)   ) -  &
                              spc_mat (1, 2, spc_n) * (spc_mat (2, 1, spc_n) * spc_mat (3, 3, spc_n) -      &
                                                       spc_mat (3, 1, spc_n) * spc_mat (2, 3, spc_n)   ) +  &
                              spc_mat (1, 3, spc_n) * (spc_mat (2, 1, spc_n) * spc_mat (3, 2, spc_n) -      &
                                                       spc_mat (3, 1, spc_n) * spc_mat (2, 2, spc_n)   )

!        Test if matrix exists already
            is_exist: DO k = 1, spc_n-1
               lexist = .true.
               DO i = 1, 4 
                  DO j = 1, 4 
                     lexist = lexist .and. ABS(spc_mat(i,j,k)-spc_mat(i,j,spc_n)) < EPS
                  ENDDO
               ENDDO
               IF(lexist) THEN
                  spc_n = spc_n -1
                  EXIT loop_power
               ENDIF
            ENDDO is_exist
            CALL get_symmetry_type (SPC_MAX, spc_n, spc_mat, spc_spur,        &
            spc_det, spc_char, spc_xyz)
         ENDDO generate
!                                                                       
!     --Set power of Generator                                          
!                                                                       
         DO i = 1, 4 
            DO j = 1, 4 
               xmat (i, j) = wmat (i, j) 
            ENDDO 
         ENDDO 
      ENDDO loop_power
!                                                                       
      END SUBROUTINE make_symmetry_matrix           
!********************************************************************** 
      SUBROUTINE get_symmetry_type (SPC_MAX, spc_n, spc_mat, spc_spur,  &
      spc_det, spc_char, spc_xyz)
!-                                                                      
!     Determines the xyz triplet, and the letter that describes the     
!     symmetry operation                                                
!+                                                                      
      IMPLICIT none 
!                                                                       
      INTEGER SPC_MAX 
      CHARACTER(65) spc_char (1:SPC_MAX) 
      CHARACTER(87) spc_xyz (1:SPC_MAX) 
      INTEGER spc_n 
      REAL spc_mat (4, 4, 1:SPC_MAX) 
      REAL spc_det (1:SPC_MAX) 
      REAL spc_spur (1:SPC_MAX) 
!                                                                       
      INTEGER i, j, ii, ja, je 
      INTEGER pwr ( - 3:3, - 1:1) 
      INTEGER power 
      CHARACTER(2) typ ( - 3:3, - 1:1) 
      CHARACTER(LEN=1) :: xyz (3) = (/'x','y','z'/)
      CHARACTER(LEN=6) :: vec (-12:12) 
      REAL work (3, 3) 
      REAL add (3) 
      REAL axis (3) 
      REAL posit (3) 
      REAL screw (3) 
      REAL hkl (3) 
      INTEGER           :: mult
      REAL              :: fact
!                                                                       
      DO i = - 3, 3 
         DO j = - 1, 1 
            typ (i, j) = ' ' 
            pwr (i, j) = 0 
         ENDDO 
      ENDDO 
      typ ( 3, 1) = ' 1' 
      typ (-1, 1) = ' 2' 
      typ ( 0, 1) = ' 3' 
      typ ( 1, 1) = ' 4' 
      typ ( 2, 1) = ' 6' 
      typ (-3,-1) = '-1' 
      typ ( 1,-1) = ' m' 
      typ ( 0,-1) = '-3' 
      typ (-1,-1) = '-4' 
      typ (-2,-1) = '-6' 
      pwr ( 3, 1) = 1 
      pwr (-1, 1) = 2 
      pwr ( 0, 1) = 3 
      pwr ( 1, 1) = 4 
      pwr ( 2, 1) = 6 
      pwr (-3,-1) = 1 
      pwr ( 1,-1) = 2 
      pwr ( 0,-1) = 3 
      pwr (-1,-1) = 4 
      pwr (-2,-1) = 6 
!                                                                       
      vec(-12) = '  -1/1' 
      vec(-11) = '-11/12' 
      vec(-10) = '  -5/6' 
      vec( -9) = '  -3/4' 
      vec( -8) = '  -2/3' 
      vec( -7) = ' -7/12' 
      vec( -6) = '  -1/2' 
      vec( -5) = ' -5/12' 
      vec( -4) = '  -1/3' 
      vec( -3) = '  -1/4' 
      vec( -2) = '  -1/6' 
      vec( -1) = ' -1/12' 
      vec ( 0) = '      ' 
      vec ( 1) = ' +1/12' 
      vec ( 2) = '  +1/6' 
      vec ( 3) = '  +1/4' 
      vec ( 4) = '  +1/3' 
      vec ( 5) = ' +5/12' 
      vec ( 6) = '  +1/2' 
      vec ( 7) = ' +7/12' 
      vec ( 8) = '  +2/3' 
      vec ( 9) = '  +3/4' 
      vec (10) = '  +5/6' 
      vec (11) = '+11/12' 
      vec (12) = '  +1/1' 
!                                                                       
      spc_char (spc_n) = ' ' 
      spc_xyz (spc_n) = ' ' 
!                                                                       
      spc_char (spc_n) = typ (NINT(spc_spur (spc_n)), NINT(spc_det (spc_n)) ) 
      power = pwr (NINT(spc_spur (spc_n)), NINT(spc_det (spc_n)) ) 
!                                                                       
      DO i = 1, 3 
         ii = (i - 1) * 29 
         DO j = 1, 3 
         ja = ii + (j - 1) * 7 + 1 
         je = ii + (j - 1) * 7 + 7 
         mult = nint (12.*spc_mat (i, j, spc_n) )
         fact =      (12.*spc_mat (i, j, spc_n) )
         IF( abs(mult-fact) < 0.001 ) THEN
            IF    ( mult == -12 ) THEN
               WRITE(spc_xyz (spc_n) (ja:je), 4000) xyz(j)
            ELSEIF( mult ==   0 ) THEN
               WRITE(spc_xyz (spc_n) (ja:je), 4100) 
            ELSEIF( mult ==  12 ) THEN
               WRITE(spc_xyz (spc_n) (ja:je), 4200) xyz(j)
            ELSE
               WRITE(spc_xyz (spc_n) (ja:je), 4300) mult,xyz(j)
            ENDIF
         ELSE
            WRITE(spc_xyz (spc_n) (ja:je), 5100) fact,xyz(j)
         ENDIF
         ENDDO 
         ja = ii + 22
         je = ii + 29 
         mult = nint (12.*spc_mat (i, 4, spc_n) )
         fact =      (12.*spc_mat (i, 4, spc_n) )
         IF( abs(mult-fact) < 0.001 ) THEN
            IF(mult < -12) THEN
               WRITE(spc_xyz (spc_n) (ja:je), 6000) mult
            ELSEIF(mult >  12) THEN
               WRITE(spc_xyz (spc_n) (ja:je), 6000) mult
            ELSE
               WRITE(spc_xyz (spc_n) (ja:je), 6200) vec(mult)
            ENDIF
         ELSE
            WRITE(spc_xyz (spc_n) (ja:je), 7100) fact
         ENDIF
      ENDDO 
      spc_xyz (spc_n) (86:87) = ' ' 
      i = 87
      call rem_bl(spc_xyz(spc_n),i)
!                                                                       
      DO i = 1, 3 
      DO j = 1, 3 
      work (i, j) = spc_mat (i, j, spc_n) 
      ENDDO 
      add (i) = spc_mat (i, 4, spc_n) 
      ENDDO 
      CALL get_detail (work, add, spc_char (spc_n), power, axis, &
      screw, posit, hkl)                                                
!                                                                       
4000  FORMAT(   '     -',A1)
4100  FORMAT(   '       ')
4200  FORMAT(   '     +',A1)
4300  FORMAT(SP,I3,'/12',A1)
5100  FORMAT(SP,F6.3    ,A1)
6000  FORMAT(SP,I3,'/12',', ')
6200  FORMAT(   A6,      ', ')
7100  FORMAT(SP,F6.3    ,', ')
!
      END SUBROUTINE get_symmetry_type              
!****&******************************************************************
      SUBROUTINE get_detail (work, add, w_char, power, axis, screw, &
      posit, hkl)                                                       
!-                                                                      
!     Determines the local symmetry of the position given in the line   
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE tensors_mod
      USE trafo_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      REAL work (3, 3) 
      REAL add (3) 
      CHARACTER(65) w_char 
      INTEGER power 
      REAL axis (3) 
      REAL screw (3) 
      REAL posit (3) 
      REAL posit1bar (3) 
      REAL hkl (3) 
!                                                                       
      CHARACTER(1) abc (0:3) 
      CHARACTER(3) xyz (3, - 2:2) 
      CHARACTER(3) xxx (3, - 2:2) 
      CHARACTER(3) Oyy (3, - 2:2) 
      CHARACTER(6) ctrans ( - 24:24) 
      INTEGER i, j, k, l 
      INTEGER ia, ie 
      INTEGER ii, jj, kk 
      INTEGER :: nnull = 0
      INTEGER nnullg 
      INTEGER iglide 
      LOGICAL lsearch 
      REAL cp (3, 3) 
      REAL temp (3, 3) 
      REAL t2 (3, 3) 
      REAL sum (3, 3) 
      REAL imat (3, 3) 
      REAL vector (3) 
      REAL p3 (3) 
      REAL hmin 
      REAL det 
      REAL :: factor = 1.0
      REAL scale 
      REAL glide 
      REAL eps 
!                                                                       
      DATA eps / 0.0001 / 
!                                                                       
      abc (0) = 'm' 
      abc (1) = 'a' 
      abc (2) = 'b' 
      abc (3) = 'c' 
      xyz (1, - 2) = '-2x' 
      xyz (1, - 1) = ' -x' 
      xyz (1, 0)  = '  0' 
      xyz (1, 1) = ' +x' 
      xyz (1, 2) = '+2x' 
      xyz (2, - 2) = '-2y' 
      xyz (2, - 1) = ' -y' 
      xyz (2, 0)  = '  0' 
      xyz (2, 1) = ' +y' 
      xyz (2, 2) = '+2y' 
      xyz (3, - 2) = '-2z' 
      xyz (3, - 1) = ' -z' 
      xyz (3, 0)  = '  0' 
      xyz (3, 1) = ' +z' 
      xyz (3, 2) = '+2z' 
      xxx (1, - 2) = '-2x' 
      xxx (1, - 1) = ' -x' 
      xxx (1, 0)  = '  0' 
      xxx (1, 1) = ' +x' 
      xxx (1, 2) = '+2x' 
      xxx (2, - 2) = '-2x' 
      xxx (2, - 1) = ' -x' 
      xxx (2, 0)  = '  0' 
      xxx (2, 1) = ' +x' 
      xxx (2, 2) = '+2x' 
      xxx (3, - 2) = '-2x' 
      xxx (3, - 1) = ' -x' 
      xxx (3, 0)  = '  0' 
      xxx (3, 1) = ' +x' 
      xxx (3, 2) = '+2x' 
      Oyy (1,  - 2)  = '  0' 
      Oyy (1,  - 1)  = '  0' 
      Oyy (1, 0)  = '  0' 
      Oyy (1, 1)  = '  0' 
      Oyy (1, 2)  = '  0' 
      Oyy (2, - 2) = '-2y' 
      Oyy (2, - 1) = ' -y' 
      Oyy (2, 0)  = '  0' 
      Oyy (2, 1) = ' +y' 
      Oyy (2, 2) = '+2y' 
      Oyy (3, - 2) = '-2y' 
      Oyy (3, - 1) = ' -y' 
      Oyy (3, 0)  = '   ' 
      Oyy (3, 1) = ' +y' 
      Oyy (3, 2) = '+2y' 
      ctrans ( - 24) = '-1/1, ' 
      ctrans ( - 23)  = '    , ' 
      ctrans ( - 22)  = '    , ' 
      ctrans ( - 21)  = '    , ' 
      ctrans ( - 20) = '-5/6, ' 
      ctrans ( - 19)  = '    , ' 
      ctrans ( - 18) = '-3/4, ' 
      ctrans ( - 17)  = '    , ' 
      ctrans ( - 16) = '-2/3, ' 
      ctrans ( - 15)  = '    , ' 
      ctrans ( - 14)  = '    , ' 
      ctrans ( - 13)  = '    , ' 
      ctrans ( - 12) = '-1/2, ' 
      ctrans ( - 11)  = '    , ' 
      ctrans ( - 10)  = '    , ' 
      ctrans ( - 9)  = '    , ' 
      ctrans ( - 8) = '-1/3, ' 
      ctrans ( - 7)  = '    , ' 
      ctrans ( - 6) = '-1/4, ' 
      ctrans ( - 5)  = '    , ' 
      ctrans ( - 4) = '-1/6, ' 
      ctrans ( - 3) = '-1/8, ' 
      ctrans ( - 2)  = '    , ' 
      ctrans ( - 1)  = '    , ' 
      ctrans (0)  = '   0, ' 
      ctrans (1)  = '    , ' 
      ctrans (2)  = '    , ' 
      ctrans (3) = '+1/8, ' 
      ctrans (4) = '+1/6, ' 
      ctrans (5)  = '    , ' 
      ctrans (6) = '+1/4, ' 
      ctrans (7)  = '    , ' 
      ctrans (8) = '+1/3, ' 
      ctrans (9)  = '    , ' 
      ctrans (10)  = '    , ' 
      ctrans (11)  = '    , ' 
      ctrans (12) = '+1/2, ' 
      ctrans (13)  = '    , ' 
      ctrans (14)  = '    , ' 
      ctrans (15)  = '    , ' 
      ctrans (16) = '+2/3, ' 
      ctrans (17)  = '    , ' 
      ctrans (18) = '+3/4, ' 
      ctrans (19)  = '    , ' 
      ctrans (20) = '+5/6, ' 
      ctrans (21)  = '    , ' 
      ctrans (22)  = '    , ' 
      ctrans (23)  = '    , ' 
      ctrans (24) = '+1/1, ' 
!                                                                       
      DO i = 1, 3 
      axis (i) = 0.0 
      screw (i) = 0.0 
      posit (i) = 0.0 
      hkl (i) = 0.0 
      ENDDO 
!                                                                       
      IF (w_char (1:2) .eq.' 1') then 
!                                                                       
!     handle (1,w), vector w is a pure centering translation            
         DO i = 1, 3 
         posit (i) = 0.0 
         screw (i) = add (i) 
         ENDDO 
      ELSEIF (w_char (1:2) .eq.'-1') then 
!                                                                       
!     handle (-1,w), vector w/2 is the position of the -1               
         DO i = 1, 3 
         posit1bar (i) = add (i) / 2.0 
         screw (i) = 0.0 
         ENDDO 
      ELSE 
!                                                                       
!     --To detemine the axis, we need to calculate                      
!       (Matrix - Identity)**-1                                         
!                                                                       
         DO i = 1, 3 
         DO j = 1, 3 
         temp (i, j) = work (i, j) 
         sum (i, j) = 0.0 
         ENDDO 
         sum (i, i) = 1.0 
         ENDDO 
!                                                                       
         factor = 1. 
         IF (w_char (1:1) .eq.'-'.or.w_char (2:2) .eq.'m') then 
            factor = - 1. 
         ENDIF 
!                                                                       
         DO i = 1, 3 
         DO j = 1, 3 
         cp (i, j) = work (i, j) 
         ENDDO 
         cp (i, i) = cp (i, i) - factor 
         ENDDO 
!                                                                       
!     --Usually the determinant of (Matrix - Identity) is singular      
!       a 2x2 submatrix can however be solved                           
!                                                                       
         lsearch = .true. 
         l = 4 
         DO while (lsearch.and.l.gt.0) 
         l = l + 1 
         ii = (mod (l - 1 + 1, 3) ) + 1 
         jj = (mod (l - 1 + 2, 3) ) + 1 
         kk = (mod (l - 1 + 3, 3) ) + 1 
         IF (cp (ii, ii) * cp (jj, jj) - cp (ii, jj) * cp (jj, ii)      &
         .ne.0.0) then                                                  
            det = cp (ii, ii) * cp (jj, jj) - cp (ii, jj) * cp (jj, ii) 
            axis (ii) = - (cp (jj, jj) * cp (ii, kk) - cp (ii, jj)      &
            * cp (jj, kk) ) / det                                       
            axis (jj) = - (cp (ii, ii) * cp (jj, kk) - cp (jj, ii)      &
            * cp (ii, kk) ) / det                                       
            axis (kk) = 1.0 
            IF (abs (abs (axis (ii) ) - 0.5) .lt.eps.or.abs (abs (axis (&
            jj) ) - 0.5) .lt.eps) then                                  
               axis (ii) = 2 * axis (ii) 
               axis (jj) = 2 * axis (jj) 
               axis (kk) = 2 * axis (kk) 
            ENDIF 
            lsearch = .false. 
         ENDIF 
         ENDDO 
!                                                                       
!     --Find how many elements are equal to zero, this helps classifying
!                                                                       
         nnull = 0 
         DO i = 1, 3 
         IF (abs (axis (i) ) .lt.eps) then 
            nnull = nnull + 1 
         ENDIF 
         ENDDO 
!                                                                       
!     --Get proper sense of direction                                   
!                                                                       
         IF (nnull.eq.0) then 
            IF (axis (1) * axis (2) * axis (3) .lt.0) then 
               DO i = 1, 3 
               axis (i) = - axis (i) 
               ENDDO 
            ENDIF 
         ELSEIF (nnull.eq.1) then 
            DO i = 1, 3 
            IF (abs (axis (i) ) .lt.eps) then 
               jj = (mod (i - 1 + 1, 3) ) + 1 
               kk = (mod (i - 1 + 2, 3) ) + 1 
               IF (axis (jj) .lt.0.0) then 
                  axis (i) = - axis (i) 
                  axis (jj) = - axis (jj) 
                  axis (kk) = - axis (kk) 
               ENDIF 
            ENDIF 
            ENDDO 
         ENDIF 
      ENDIF 
!                                                                       
!     Determine screw part/glide part, exclude -N axes since t=0 always 
!                                                                       
      IF (w_char (2:2) .ne.'1') then 
         IF (w_char (1:1) .ne.'-') then 
!                                                                       
!     --First add M**(l-1) for l = 1 to power-1                         
            DO l = 2, power 
            DO i = 1, 3 
            DO j = 1, 3 
            sum (i, j) = sum (i, j) + temp (i, j) 
            ENDDO 
            ENDDO 
            IF (l.lt.power) then 
               CALL matmulx (t2, temp, work) 
               DO i = 1, 3 
               DO j = 1, 3 
               temp (i, j) = t2 (i, j) 
               ENDDO 
               ENDDO 
            ENDIF 
            ENDDO 
            DO i = 1, 3 
            screw (i) = 0.0 
            DO j = 1, 3 
            screw (i) = screw (i) + sum (i, j) * add (j) 
            ENDDO 
            screw (i) = screw (i) / power 
            ENDDO 
         ENDIF 
!                                                                       
!     --Subtract glide part from translation part                       
!                                                                       
         DO i = 1, 3 
         vector (i) = add (i) - screw (i) 
         ENDDO 
!                                                                       
!     --Determine Position                                              
!                                                                       
         DO i = 1, 3 
         DO j = 1, 3 
         cp (i, j) = - work (i, j) 
         ENDDO 
         cp (i, i) = 1.0 - work (i, i) 
         ENDDO 
         det = cp (1, 1) * (cp (2, 2) * cp (3, 3) - cp (3, 2) * cp (2,  &
         3) ) - cp (1, 2) * (cp (2, 1) * cp (3, 3) - cp (3, 1) * cp (2, &
         3) ) + cp (1, 3) * (cp (2, 1) * cp (3, 2) - cp (3, 1) * cp (2, &
         2) )                                                           
!                                                                       
         IF (det.ne.0) then 
!     --This allows to calculate the position directly                  
            CALL invmat (imat, cp) 
            DO i = 1, 3 
            posit1bar (i) = 0.0 
            DO j = 1, 3 
            posit1bar (i) = posit1bar (i) + imat (i, j) * vector (j) 
            ENDDO 
            posit (i) = posit1bar (i) 
            ENDDO 
         ELSE 
            DO i = 1, 3 
            DO j = 1, 3 
            temp (i, j) = work (i, j) * factor 
            sum (i, j) = 0.0 
            ENDDO 
            sum (i, i) = power - 1.0 
            ENDDO 
!                                                                       
!     --First add (power-l)*M**(l-1) for l = 1 to power-1               
            DO l = 2, power - 1 
            DO i = 1, 3 
            DO j = 1, 3 
            sum (i, j) = sum (i, j) + (power - l) * temp (i, j) 
            ENDDO 
            ENDDO 
            IF (l.lt.power - 1) then 
               CALL matmulx (t2, temp, work) 
               DO i = 1, 3 
               DO j = 1, 3 
               temp (i, j) = t2 (i, j) 
               ENDDO 
               ENDDO 
            ENDIF 
            ENDDO 
            DO i = 1, 3 
            posit (i) = 0.0 
            DO j = 1, 3 
            posit (i) = posit (i) + sum (i, j) * vector (j) 
            ENDDO 
            posit (i) = posit (i) / power 
            ENDDO 
         ENDIF 
         nnull = 0 
         DO i = 1, 3 
         IF (abs (axis (i) ) .lt.eps) then 
            nnull = nnull + 1 
         ENDIF 
         ENDDO 
      ENDIF 
!                                                                       
!     Cleanup to look like Tables                                       
!                                                                       
      DO i = 1, 3 
      screw (i) = screw (i) - int (screw (i) ) 
      posit (i) = posit (i) - int (posit (i) ) 
      posit1bar (i) = posit1bar (i) - int (posit1bar (i) ) 
      ENDDO 
      IF (w_char (1:2) .eq.' 1') then 
!     --write centerin vector                                           
         IF (abs (screw (1) ) .gt.eps.or.abs (screw (2) )               &
         .gt.eps.or.abs (screw (3) ) .gt.eps) then                      
            WRITE (w_char (4:21), 3100) (ctrans (nint (24 * screw (i) ) &
            ), i = 1, 3)                                                
         ENDIF 
      ELSEIF (w_char (1:2) .eq.'-1') then 
!     --write position of center of symmetry                            
         WRITE (w_char (46:63), 3200) (ctrans (nint (24 * posit1bar (i) &
         ) ), i = 1, 3)                                                 
      ELSEIF (w_char (1:1) .eq.'-') then 
!     --write position of center of symmetry                            
         WRITE (w_char (46:63), 3200) (ctrans (nint (24 * posit1bar (i) &
         ) ), i = 1, 3)                                                 
         w_char (45:45) = ';' 
      ELSEIF (w_char (2:2) .eq.'m') then 
!     --for mirror plane, convert axis to hkl                           
         CALL trans (axis, cr_gten, hkl, 3) 
         hmin = 1.e10 
         DO i = 1, 3 
         IF (abs (hkl (i) ) .gt.eps) then 
            hmin = min (hmin, abs (hkl (i) ) ) 
         ENDIF 
         ENDDO 
         DO i = 1, 3 
         hkl (i) = float (nint (hkl (i) / hmin) ) 
         ENDDO 
      ENDIF 
      IF (w_char (2:2) .ne.'1'.and.w_char (2:2) .ne.'m') then 
!                                                                       
!     --For rotation axis, clean up positions to look like tables       
!                                                                       
         IF (nnull.eq.0) then 
            posit (1) = posit (1) - posit (3) * sign (1., axis (1) )    &
            * sign (1., axis (3) )                                      
            posit (2) = posit (2) - posit (3) * sign (1., axis (2) )    &
            * sign (1., axis (3) )                                      
            posit (3) = 0.0 
         ELSE 
            DO l = 1, 3 
            i = l 
            j = mod (l, 3) + 1 
            k = mod (l + 1, 3) + 1 
            IF (axis (i) .ne.0.and.axis (j) .eq.0.and.axis (k) .eq.0)   &
            then                                                        
               posit (i) = 0.0 
            ELSEIF (axis (i) .eq.axis (j) .and.axis (k) .eq.0) then 
               posit (j) = posit (j) - posit (i) 
               posit (i) = 0.0 
            ELSEIF (axis (i) .gt.0.and.axis (i) .eq. - axis (j)         &
            .and.axis (k) .eq.0) then                                   
               posit (j) = posit (j) + posit (i) 
               posit (i) = 0.0 
            ELSEIF (axis (i) .lt.0.and.axis (i) .eq. - axis (j)         &
            .and.axis (k) .eq.0) then                                   
               posit (i) = posit (i) + posit (j) 
               posit (j) = 0.0 
            ENDIF 
            ENDDO 
         ENDIF 
         DO i = 1, 3 
         screw (i) = screw (i) - int (screw (i) ) 
         posit (i) = posit (i) - int (posit (i) ) 
         ENDDO 
         IF (abs (screw (1) ) .gt.eps.or.abs (screw (2) )               &
         .gt.eps.or.abs (screw (3) ) .gt.eps) then                      
            WRITE (w_char (4:21), 3100) (ctrans (nint (24 * screw (i) ) &
            ), i = 1, 3)                                                
         ENDIF 
         DO l = 1, 3 
         IF (abs (posit (l) ) .gt.eps) then 
            ia = 21 + (l - 1) * 8 + 4 
            ie = 21 + (l - 1) * 8 + 8 
            w_char (ia:ie) = ctrans (nint (24. * posit (l) ) ) 
         ENDIF 
         ENDDO 
         IF (w_char (1:1) .ne.'-') then 
            w_char (45:45) = ' ' 
         ELSE 
            IF (power.gt.2) then 
               w_char (45:45) = ';' 
            ENDIF 
         ENDIF 
      ELSEIF (w_char (2:2) .eq.'m') then 
         DO l = 1, 3 
         i = l 
         j = mod (l, 3) + 1 
         k = mod (l + 1, 3) + 1 
         IF (hkl (i) .ne.0.and.hkl (j) .eq.0.and.hkl (k) .eq.0) then 
            posit (j) = 0.0 
            posit (k) = 0.0 
         ELSEIF (hkl (i) .eq.hkl (j) .and.hkl (k) .eq.0) then 
            posit (i) = posit (j) + posit (i) 
            posit (j) = 0.0 
            posit (k) = 0.0 
         ELSEIF (hkl (i) .gt.0.and.hkl (i) .eq. - hkl (j) .and.hkl (k)  &
         .eq.0) then                                                    
            posit (i) = posit (i) - posit (j) 
            posit (j) = 0.0 
            posit (k) = 0.0 
         ELSEIF (hkl (i) .lt.0.and.hkl (i) .eq. - hkl (j) .and.hkl (k)  &
         .eq.0) then                                                    
            posit (j) = posit (j) - posit (i) 
            posit (i) = 0.0 
            posit (k) = 0.0 
         ENDIF 
         ENDDO 
         DO i = 1, 3 
         screw (i) = screw (i) - int (screw (i) ) 
         posit (i) = posit (i) - int (posit (i) ) 
         ENDDO 
         IF (abs (hkl (1) ) .eq.2..and.abs (hkl (2) ) .eq.1.) then 
            posit (2) = posit (2) + posit (1) * sign (2., hkl (1) )     &
            * sign (1., hkl (2) )                                       
            posit (1) = 0.0 
         ELSEIF (abs (hkl (1) ) .eq.1..and.abs (hkl (2) ) .eq.2.) then 
            posit (1) = posit (1) + posit (2) * sign (2., hkl (2) )     &
            * sign (1., hkl (1) )                                       
            posit (2) = 0.0 
         ENDIF 
         DO l = 1, 3 
         IF (abs (posit (l) ) .gt.eps) then 
            ia = 21 + (l - 1) * 8 + 4 
            ie = 21 + (l - 1) * 8 + 8 
            w_char (ia:ie) = ctrans (nint (24. * posit (l) ) ) 
         ENDIF 
         ENDDO 
         w_char (45:45) = ' ' 
      ENDIF 
!                                                                       
!     Get x,y,z symbol for the axis                                     
!                                                                       
      IF (w_char (2:2) .ne.'1'.and.w_char (2:2) .ne. ('m') ) then 
         DO i = 1, 3 
         ia = 21 + (i - 1) * 8 + 1 
         ie = 21 + (i - 1) * 8 + 3 
         IF (nnull.eq.2) then 
            w_char (ia:ie) = xyz (i, nint (axis (i) ) ) 
         ELSEIF (nnull.eq.1) then 
            IF (axis (1) .ne.0.) then 
               w_char (ia:ie) = xxx (i, nint (axis (i) ) ) 
            ELSE 
               w_char (ia:ie) = Oyy (i, nint (axis (i) ) ) 
            ENDIF 
         ELSEIF (nnull.eq.0) then 
            w_char (ia:ie) = xxx (i, nint (axis (i) ) ) 
         ENDIF 
      IF ( (w_char (ie+1:ie+1) .eq.'+'.or.w_char (ie+1:ie+1) .eq.'-') .a&
     &nd.w_char (ia:ie) .eq.'  0') then                                 
      w_char (ia:ie)  = '  ' 
         ENDIF 
         ENDDO 
         w_char (29:29) = ',' 
         w_char (37:37) = ',' 
      ENDIF 
!                                                                       
!     Get x,y,z symbol for the mirror plane                             
!                                                                       
      IF (w_char (2:2) .eq.'m') then 
         nnull = 0 
         DO i = 1, 3 
         IF (hkl (i) .eq.0.) then 
            nnull = nnull + 1 
         ENDIF 
         ENDDO 
         IF (nnull.eq.2) then 
            DO i = 1, 3 
            ia = 21 + (i - 1) * 8 + 1 
            ie = 21 + (i - 1) * 8 + 3 
            w_char (ia:ie) = xyz (i, nint (1 - hkl (i) ) ) 
            ENDDO 
         ELSEIF (nnull.eq.1) then 
            lsearch = .true. 
            DO l = 3, 5 
            IF (lsearch) then 
               i = mod (l + 2, 3) + 1 
               j = mod (l, 3) + 1 
               k = mod (l + 1, 3) + 1 
               IF (abs (hkl (i) ) .lt.eps) then 
                  ia = 21 + (i - 1) * 8 + 1 
                  ie = 21 + (i - 1) * 8 + 3 
                  w_char (ia:ie) = xyz (i, 1) 
                  jj = 1 
                  scale = - jj * (hkl (j) / hkl (k) ) 
                  kk = - nint (jj * hkl (j) / hkl (k) ) 
                  IF (abs (scale) .lt.1) then 
                     jj = int( jj / abs (scale) )
                     kk = - nint (jj * hkl (j) / hkl (k) ) 
                  ENDIF 
                  IF (abs (hkl (1) ) .lt.eps) then 
                     ia = 21 + (j - 1) * 8 + 1 
                     ie = 21 + (j - 1) * 8 + 3 
                     w_char (ia:ie) = Oyy (j, jj) 
                     ia = 21 + (k - 1) * 8 + 1 
                     ie = 21 + (k - 1) * 8 + 3 
                     w_char (ia:ie) = Oyy (k, kk) 
                  ELSE 
                     ia = 21 + (j - 1) * 8 + 1 
                     ie = 21 + (j - 1) * 8 + 3 
                     w_char (ia:ie) = xxx (j, jj) 
                     ia = 21 + (k - 1) * 8 + 1 
                     ie = 21 + (k - 1) * 8 + 3 
                     w_char (ia:ie) = xxx (k, kk) 
                  ENDIF 
                  lsearch = .false. 
               ENDIF 
            ENDIF 
            ENDDO 
         ENDIF 
         DO i = 1, 3 
         ia = 21 + (i - 1) * 8 + 1 
         ie = 21 + (i - 1) * 8 + 3 
      IF ( (w_char (ie+1:ie+1) .eq.'+'.or.w_char (ie+1:ie+1) .eq.'-') .a&
     &nd.w_char (ia:ie) .eq.'  0') then                                 
      w_char (ia:ie)  = '  ' 
         ENDIF 
         ENDDO 
         w_char (29:29) = ',' 
         w_char (37:37) = ',' 
!                                                                       
!     --get glide plane symbol                                          
!                                                                       
         IF (abs (screw (1) ) .gt.eps.or.abs (screw (2) )               &
         .gt.eps.or.abs (screw (3) ) .gt.eps) then                      
            nnullg = 0 
            iglide = 0 
            DO i = 1, 3 
            IF (screw (i) .eq.0.) then 
               nnullg = nnullg + 1 
            ELSE 
               iglide = i 
            ENDIF 
            ENDDO 
            IF (nnullg.eq.2) then 
               w_char (2:2) = abc (iglide) 
            ELSEIF (nnull.eq.2.and.nnullg.eq.1) then 
               IF (abs (abs (screw (1) ) - abs (screw (2) ) )           &
               .lt.eps.or.abs (abs (screw (1) ) - abs (screw (3) ) )    &
               .lt.eps.or.abs (abs (screw (2) ) - abs (screw (3) ) )    &
               .lt.eps) then                                            
                  glide = max (abs (screw (1) ), abs (screw (2) ),      &
                  abs (screw (3) ) )                                    
                  IF (glide.eq.0.50) then 
                     w_char (2:2) = 'n' 
                  ELSEIF (glide.eq.0.25) then 
                     w_char (2:2) = 'd' 
                  ELSE 
                     w_char (2:2) = 'g' 
                  ENDIF 
               ELSE 
                  w_char (2:2) = 'g' 
               ENDIF 
            ELSEIF (nnull.eq.1.and.nnullg.eq.0) then 
               IF (abs (abs (screw (1) ) - abs (screw (2) ) )           &
               .lt.eps.and.abs (abs (screw (1) ) - abs (screw (3) ) )   &
               .lt.eps) then                                            
                  glide = max (abs (screw (1) ), abs (screw (2) ),      &
                  abs (screw (3) ) )                                    
                  IF (glide.eq.0.50) then 
                     w_char (2:2) = 'n' 
                  ELSEIF (glide.eq.0.25) then 
                     w_char (2:2) = 'd' 
                  ELSE 
                     w_char (2:2) = 'g' 
                  ENDIF 
               ELSE 
                  w_char (2:2) = 'g' 
               ENDIF 
            ELSE 
               w_char (2:2) = 'g' 
            ENDIF 
            WRITE (w_char (4:21), 3100) (ctrans (nint (24 * screw (i) ) &
            ), i = 1, 3)                                                
         ENDIF 
      ENDIF 
!                                                                       
!     get sense of rotation                                             
!                                                                       
      IF (power.gt.2) then 
         vector (1) = 3 
         vector (2) = 7 
         vector (3) = 19 
         DO i = 1, 3 
         p3 (i) = 0.0 
         DO j = 1, 3 
         p3 (i) = p3 (i) + work (i, j) * vector (j) 
         ENDDO 
         ENDDO 
         DO i = 1, 3 
         cp (i, 1) = axis (i) 
         cp (i, 2) = vector (i) 
         cp (i, 3) = p3 (i) 
         ENDDO 
         det = (cp (1, 1) * (cp (2, 2) * cp (3, 3) - cp (3, 2) * cp (2, &
         3) ) - cp (1, 2) * (cp (2, 1) * cp (3, 3) - cp (3, 1) * cp (2, &
         3) ) + cp (1, 3) * (cp (2, 1) * cp (3, 2) - cp (3, 1) * cp (2, &
         2) ) ) * factor                                                
         IF (det.gt.0) then 
            w_char (3:3) = 'P' 
         ELSEIF (det.lt.0) then 
            w_char (3:3) = 'M' 
         ENDIF 
      ENDIF 
!                                                                       
!     -rename a centering vector                                        
!                                                                       
      IF (w_char (1:2) .eq.' 1') then 
         IF (abs (screw (1) ) .gt.eps.or.abs (screw (2) )               &
         .gt.eps.or.abs (screw (2) ) .gt.eps) then                      
            w_char (1:2) = ' t' 
         ENDIF 
      ENDIF 
!                                                                       
!DBG      write(*,*) ' Axis               ',axis                        
!DBG      write(*,*) ' Screw/Glide        ',screw                       
!DBG      write(*,*) ' Position           ',posit                       
!DBG      if(w_char(1:1).eq.'-') then                                   
!DBG      write(*,*) ' Position -1        ',posit1bar                   
!DBG      endif                                                         
!DBG      write(*,*) ' hkl                ',hkl                         
!DBG      write(*,4000)                     no,mod(no-1,48)+1,w_char    
!DBG      write(*,*) '          ',                                      
!DBG     &                      '123456789 123456789 123456789 ',       
!DBG     &                      '123456789 123456789 123456789 '        
!                                                                       
 3100 FORMAT    (' (',a4,',',a4,',',a4,') ') 
 3200 FORMAT    ('  ',a4,',',a4,',',a4,'  ') 
!                                                                       
      END SUBROUTINE get_detail                     
!********************************************************************** 
      SUBROUTINE wyckoff_main (zeile, lp) 
!-                                                                      
!     Determines the local symmetry of the position given in the line   
!+                                                                      
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER FULL, SYMBOL, XYZ, MATRIX 
      PARAMETER (FULL = 0) 
      PARAMETER (SYMBOL = 1) 
      PARAMETER (XYZ = 2) 
      PARAMETER (MATRIX = 3) 
!                                                                       
      CHARACTER(1024) zeile 
      INTEGER lp 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 4) 
!                                                                       
      CHARACTER(1024) cpara (maxw) 
      INTEGER lpara (maxw) 
      INTEGER ianz 
      INTEGER iianz 
      INTEGER mode 
      LOGICAL loutput 
      REAL werte (maxw) 
!                                                                       
      LOGICAL str_comp 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
      IF (ianz.eq.3.or.ianz.eq.4) then 
         iianz = 3 
         CALL ber_params (iianz, cpara, lpara, werte, maxw) 
         IF (ier_num.ne.0) return 
!                                                                       
         loutput = .true. 
         IF (ianz.eq.3.or.str_comp (cpara (4) , 'full', 2, lpara (4) ,  &
         4) ) then                                                      
            mode = FULL 
         ELSEIF (str_comp (cpara (4) , 'symbol', 2, lpara (4) , 6) )    &
         then                                                           
            mode = SYMBOL 
         ELSEIF (str_comp (cpara (4) , 'xyz', 2, lpara (4) , 3) ) then 
            mode = XYZ 
         ELSEIF (str_comp (cpara (4) , 'matrix', 2, lpara (4) , 6) )    &
         then                                                           
            mode = MATRIX 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
            RETURN 
         ENDIF 
         CALL get_wyckoff (werte, loutput, mode) 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
      END SUBROUTINE wyckoff_main                   
!********************************************************************** 
      SUBROUTINE get_wyckoff (vec, loutput, mode) 
!-                                                                      
!     Determines the local symmetry of position xyz within the unit cell
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE wyckoff_mod 
      USE unitcell_mod 
      USE prompt_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      REAL vec (3) 
      LOGICAL loutput 
      INTEGER mode 
!                                                                       
      INTEGER FULL, SYMBOL, XYZ, MATRIX 
      PARAMETER (FULL = 0) 
      PARAMETER (SYMBOL = 1) 
      PARAMETER (XYZ = 2) 
      PARAMETER (MATRIX = 3) 
!                                                                       
      INTEGER i, j 
      INTEGER is 
      INTEGER n_center 
      INTEGER igroup 
      INTEGER :: block = 1
      LOGICAL lident 
      REAL orig (4) 
      REAL copy (4) 
      REAL eps 
!                                                                       
      DATA eps / 0.00001 / 
!                                                                       
      n_center = 1 
      IF (cr_spcgr (1:1) .eq.'P') then 
         n_center = 1 
      ELSEIF (cr_spcgr (1:1) .eq.'A') then 
         n_center = 2 
      ELSEIF (cr_spcgr (1:1) .eq.'B') then 
         n_center = 2 
      ELSEIF (cr_spcgr (1:1) .eq.'C') then 
         n_center = 2 
      ELSEIF (cr_spcgr (1:1) .eq.'I') then 
         n_center = 2 
      ELSEIF (cr_spcgr (1:1) .eq.'F') then 
         n_center = 4 
      ELSEIF (cr_spcgr (1:1) .eq.'R'.and.cr_syst.eq.6) then 
         n_center = 3 
      ENDIF 
      IF (gen_sta.eq.GEN_SYMM) then 
         block = spc_n / n_center 
      ELSEIF (gen_sta.eq.GEN_CENTER) then 
         block = n_center 
      ENDIF 
      IF (loutput) then 
         WRITE (output_io, 900) vec 
      ENDIF 
!                                                                       
      wyc_n = 0 
!                                                                       
!     move position into first unit cell,ia                             
!                                                                       
      CALL firstcell (vec, 3) 
      DO i = 1, 3 
      orig (i) = vec (i) 
      ENDDO 
      orig (4) = 1.0 
!                                                                       
!     apply all symmetry operations to original position                
!                                                                       
      DO is = 1, spc_n 
      IF (gen_sta.eq.GEN_SYMM) then 
         igroup = mod (is - 1, block) + 1 
      ELSEIF (gen_sta.eq.GEN_CENTER) then 
         igroup = (is - 1) / block + 1 
      ENDIF 
      DO i = 1, 4 
      copy (i) = 0.0 
      DO j = 1, 4 
      copy (i) = copy (i) + spc_mat (i, j, is) * orig (j) 
      ENDDO 
      ENDDO 
      CALL firstcell (copy, 4) 
      lident = .true. 
      DO i = 1, 3 
      lident = lident.and.abs (orig (i) - copy (i) ) .lt.eps 
      ENDDO 
      IF (lident) then 
         wyc_n = wyc_n + 1 
         wyc_list (wyc_n) = is 
         IF (loutput) then 
            IF (mode.eq.FULL) then 
               WRITE (output_io, 1000) is, igroup 
               WRITE (output_io, 1100) (spc_mat (1, j, is), j = 1, 4),  &
               spc_char (is), (spc_mat (2, j, is), j = 1, 4), (spc_mat (&
               3, j, is), j = 1, 4), spc_xyz (is)                       
            ELSEIF (mode.eq.SYMBOL) then 
               WRITE (output_io, 3200) is, igroup, spc_char (is) 
            ELSEIF (mode.eq.XYZ) then 
               WRITE (output_io, 4200) is, igroup, spc_xyz (is) 
            ELSEIF (mode.eq.MATRIX) then 
               WRITE (output_io, 5200) is, igroup, (spc_mat (1, j, is), &
               j = 1, 4), (spc_mat (2, j, is), j = 1, 4), (spc_mat (3,  &
               j, is), j = 1, 4)                                        
            ENDIF 
         ENDIF 
      ENDIF 
!                                                                       
      ENDDO 
!                                                                       
      IF (loutput) then 
         WRITE (output_io, 6000) spc_n / wyc_n, wyc_n, spc_n 
      ENDIF 
!                                                                       
  900 FORMAT    (/,' Wyckoff symmetry for position ',3f12.6,/) 
 1000 FORMAT    ('Symmetry No.      [',i3,']  (',i3,')') 
 1100 FORMAT    (  ' ( ',3(f4.1,', '),f8.5,' )','  ',a65,/,             &
     &                    ' ( ',3(f4.1,', '),f8.5,' )',/,               &
     &                    ' ( ',3(f4.1,', '),f8.5,' )','  ',a87,/)      
 3200 FORMAT    ('Symmetry No.      [',i3,']  (',i3,')  ',a65) 
 4200 FORMAT    ('Symmetry No.      [',i3,']  (',i3,')  ',a87) 
 5200 FORMAT    ('Symmetry No.      [',i3,']  (',i3,')  ',              &
     &                    ' ( ',3(f4.1,', '),f8.5,' )',/,               &
     &                32x,' ( ',3(f4.1,', '),f8.5,' )',/,               &
     &                32x,' ( ',3(f4.1,', '),f8.5,' )'   )              
 6000 FORMAT    (/,' Multiplicity;   No of Sym. Op. in Wyckoff group; ',&
     &  '    Highest Multiplicity',/,i8,20x,i8,20x,i8)                  
!                                                                       
      END SUBROUTINE get_wyckoff                    
!********************************************************************** 
      SUBROUTINE symmetry 
!-                                                                      
!     Performs the space group symmetry on the current atom.            
!     cr_natoms             the current atom number                     
!     cr_iscat(cr_natoms)   its scattering type                         
!                           and thus number of chemically identical     
!                           yet symmetrically different atoms           
!                                                                       
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE generate_mod 
      USE gen_add_mod 
      USE molecule_mod 
      USE sym_add_mod 
      USE unitcell_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      REAL eps 
      INTEGER ii, iii, igs, igg 
      INTEGER iiii 
!                                                                       
      DATA eps / 0.00001 / 
!                                                                       
!     ii is the number of the atom the symmetry is to work on           
!                                                                       
      ii = cr_natoms 
!                                                                       
!     Loop over all generators for spacegroup cr_spcgrno                
!                                                                       
      DO igs = 1, generspcgr (0, cr_spcgrno) 
      IF (gen_sta.eq.GEN_SYMM) then 
         igg = generspcgr (igs, cr_spcgrno) 
      ELSEIF (gen_sta.eq.GEN_CENTER) then 
         igg = generspcgr_center (igs, cr_spcgrno) 
      ENDIF 
!                                                                       
!     --iii is the number of the last atom generated by the previous    
!     --generator                                                       
!                                                                       
      iii = cr_natoms 
!                                                                       
      CALL symmetry_gener (NMAX, cr_natoms, cr_pos, cr_iscat, cr_prop,  &
      cr_mole, ii, iii, iii, igg, NG, generators, generpower)                    
!                                                                       
      IF (ier_num.ne.0) then 
         RETURN 
      ENDIF 
      ENDDO 
!                                                                       
!     Loop over all additional generators                               
!                                                                       
      DO igs = 1, gen_add_n 
!                                                                       
!     --iii is the number of the last atom generated by the previous    
!     --generator                                                       
!                                                                       
      iii = cr_natoms 
!                                                                       
      CALL symmetry_gener (NMAX, cr_natoms, cr_pos, cr_iscat, cr_prop,  &
      cr_mole, ii, iii, iii, igs, GEN_ADD_MAX, gen_add, gen_add_power)           
!                                                                       
      IF (ier_num.ne.0) then 
         RETURN 
      ENDIF 
      ENDDO 
!                                                                       
!     Loop over all additional symmetry operations that                 
!       are not generators                                              
!                                                                       
!                                                                       
!     --iiii is the number of the last atom generated by the last       
!     --generator                                                       
!                                                                       
      iiii = cr_natoms 
!                                                                       
      DO igs = 1, sym_add_n 
!                                                                       
!     --iii is the number of the last atom generated by the previous    
!     --symmetry operation                                              
!                                                                       
      iii = cr_natoms 
!                                                                       
      CALL symmetry_gener (NMAX, cr_natoms, cr_pos, cr_iscat, cr_prop,  &
      cr_mole, ii, iii, iiii, igs, SYM_ADD_MAX, sym_add, sym_add_power)          
!                                                                       
      IF (ier_num.ne.0) then 
         RETURN 
      ENDIF 
      ENDDO 
!                                                                       
!     Sort Atoms into molecules                                         
!                                                                       
      IF (mole_l_on) then 
         CALL mole_insert (ii) 
         IF (ier_num.ne.0) then 
            RETURN 
         ENDIF 
!                                                                       
      ENDIF 
!                                                                       
      END SUBROUTINE symmetry                       
!********************************************************************** 
      SUBROUTINE symmetry_gener (NMAX, cr_natoms, cr_pos, cr_iscat,     &
      cr_prop, cr_mole, ii, iii, iiii, igg, NG, generators, generpower)          
!-                                                                      
!     Applies the generator to the current atom                         
!+                                                                      
      USE molecule_mod 
      USE trafo_mod
      IMPLICIT none 
!                                                                       
      INTEGER NMAX 
!
      INTEGER,                       INTENT(INOUT)  :: cr_natoms
      INTEGER, DIMENSION(1:NMAX),    INTENT(INOUT)  :: cr_iscat
      INTEGER, DIMENSION(1:NMAX),    INTENT(INOUT)  :: cr_mole
      INTEGER, DIMENSION(1:NMAX),    INTENT(INOUT)  :: cr_prop
      REAL   , DIMENSION(1:3,1:NMAX),INTENT(INOUT)  :: cr_pos
!                                                                       
      INTEGER NG 
      INTEGER ii, iii, iiii 
      INTEGER igg 
      INTEGER generpower (NG) 
!                                                                       
      REAL generators (4, 4, 0:NG) 
!                                                                       
      INTEGER ia, iaa, ipg 
      INTEGER i, j, k 
      LOGICAL lnew 
      REAL x (4), y (4)
      REAL wmat (4, 4) 
      REAL xmat (4, 4) 
      REAL eps 
      REAL compare (4), previous (4) 
!                                                                       
      DATA eps / 0.00001 / 
!                                                                       
!     For convenience create Identity operator                          
!                                                                       
      DO i = 1, 4 
      DO j = 1, 4 
      xmat (i, j) = 0.0 
      ENDDO 
      xmat (i, i) = 1.0 
      ENDDO 
!                                                                       
!     Loop over all powers of generator igg                             
!                                                                       
      DO ipg = 1, generpower (igg) 
!                                                                       
!     --raise power of generator, xmat is a dummy matrix, equal to      
!     --the previous power of the generator                             
!                                                                       
      DO i = 1, 4 
      DO j = 1, 4 
      wmat (i, j) = 0.0 
      DO k = 1, 4 
      wmat (i, j) = wmat (i, j) + generators (i, k, igg) * xmat (k, j) 
      ENDDO 
      ENDDO 
      ENDDO 
!                                                                       
!     --Multiply all points with current generator                      
!                                                                       
      DO iaa = ii, iiii 
      DO i = 1, 3 
      x (i) = cr_pos (i, iaa) 
      ENDDO 
      x (4) = 1.0 
      CALL trans (x, wmat, y, 4) 
      DO i = 1, 3 
      compare (i) = y (i) 
      ENDDO 
      compare (4) = 1.0 
      CALL firstcell (compare, 4) 
      IF (.not.mole_l_on) then 
!                                                                       
!     ------Transform atom into first unit cell,                        
!           if it is not inside a molecule                              
!                                                                       
         CALL firstcell (y, 4) 
      ENDIF 
      lnew = .true. 
      DO ia = ii, iii 
      DO i = 1, 3 
      previous (i) = cr_pos (i, ia) 
      ENDDO 
      previous (4) = 1.0 
      CALL firstcell (previous, 4) 
!     ------check if atom exists                                        
      IF (abs (compare (1) - previous (1) ) .lt.eps.and.abs (compare (2)&
      - previous (2) ) .lt.eps.and.abs (compare (3) - previous (3) )    &
      .lt.eps) then                                                     
         lnew = .false. 
         GOTO 30 
      ENDIF 
      ENDDO 
   30 CONTINUE 
!                                                                       
!     ------insert atom into crystal                                    
!                                                                       
      IF (lnew) then 
         IF (cr_natoms.lt.nmax) then 
            cr_natoms = cr_natoms + 1 
            cr_pos (1, cr_natoms) = y (1) 
            cr_pos (2, cr_natoms) = y (2) 
            cr_pos (3, cr_natoms) = y (3) 
            cr_iscat (cr_natoms) = cr_iscat (ii) 
            cr_mole (cr_natoms) = cr_mole (ii) 
            cr_prop (cr_natoms) = cr_prop (ii) 
         ELSE 
            ier_num = -10 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
      ENDIF 
      ENDDO 
!                                                                       
!     set xmat to current power of generator                            
!                                                                       
      DO i = 1, 4 
      DO j = 1, 4 
      xmat (i, j) = wmat (i, j) 
      ENDDO 
      ENDDO 
!     --End of loop over all powers                                     
      ENDDO 
      END SUBROUTINE symmetry_gener                 
!********************************************************************** 
      SUBROUTINE mole_insert (ii) 
!-                                                                      
!     Sorts the newly created atoms into the correct molecules.         
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE molecule_mod 
      USE wyckoff_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER ii 
!                                                                       
      INTEGER i, j, l 
      INTEGER m
      INTEGER is 
      INTEGER ifirst 
      INTEGER mole_st 
      INTEGER mole_temp (192) 
!                                                                       
      LOGICAL lsame 
!                                                                       
      REAL vec (3), orig (4), old (3) 
      REAL first (3) 
      REAL eps 
!                                                                       
      DATA eps / 0.00001 / 
!                                                                       
                                                                        
      DO i = 1, 192 
      mole_temp (i) = 0 
      ENDDO 
!                                                                       
!     The first atom of a molecule and its symmetrically equivalent     
!     are each sorted into new molecules                                
!                                                                       
      IF (mole_l_first) then 
         mole_num_act = mole_num_curr - 1 
         DO i = ii, cr_natoms 
         mole_num_act = mole_num_act + 1 
         CALL mole_insert_current (i, mole_num_act) 
         IF (ier_num.ne.0) then 
            RETURN 
         ENDIF 
         ENDDO 
         mole_l_first = .false. 
      ELSE 
!                                                                       
!     --These are atoms further down the line, here we have to check    
!     --into which molecule they belong. The first is inserted into     
!     --the current molecule                                            
!                                                                       
         mole_num_act = mole_num_curr 
         mole_st = mole_len (mole_num_curr) + 1 
!                                                                       
!     --Loop over all secondary atoms                                   
!                                                                       
         DO i = ii, cr_natoms 
!                                                                       
!     ----If not yet in a molecule, add to active molecule              
!                                                                       
         IF (mole_temp (i - ii + 1) .eq.0) then 
            CALL mole_insert_current (i, mole_num_act) 
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
            mole_temp (i - ii + 1) = mole_num_act 
!                                                                       
!     ----Determine coordinates of first atom in active molecule        
!                                                                       
            ifirst = mole_cont (mole_off (mole_num_act) + 1) 
            first (1) = cr_pos (1, ifirst) 
            first (2) = cr_pos (2, ifirst) 
            first (3) = cr_pos (3, ifirst) 
!                                                                       
            CALL get_wyckoff (first, .false., 0) 
!                                                                       
!     ----Create copies of this atom with all Wyckoff symmetry operators
!     ----1.st Wyckoff Symmetry is always identity, thus we             
!       can ignore this                                                 
!                                                                       
            orig (1) = cr_pos (1, i) 
            orig (2) = cr_pos (2, i) 
            orig (3) = cr_pos (3, i) 
            orig (4) = 1.0 
            DO is = 2, wyc_n 
            DO l = 1, 3 
            vec (l) = 0.0 
            DO m = 1, 4 
            vec (l) = vec (l) + spc_mat (l, m, wyc_list (is) ) * orig ( &
            m)                                                          
            ENDDO 
            ENDDO 
            CALL firstcell (vec, 3) 
!                                                                       
!     ------Loop to compare the copies to the remaining original atom   
!                                                                       
            DO j = i + 1, cr_natoms 
            DO l = 1, 3 
            old (l) = cr_pos (l, j) 
            ENDDO 
            CALL firstcell (old, 3) 
            lsame = .true. 
            DO l = 1, 3 
            lsame = lsame.and.abs (old (l) - vec (l) ) .lt.eps 
            ENDDO 
            IF (lsame) then 
!     ----------This is another atom of the same molecule               
               IF (mole_temp (j - ii + 1) .eq.0) then 
                  CALL mole_insert_current (j, mole_num_act) 
                  mole_temp (j - ii + 1) = mole_num_act 
               ENDIF 
            ENDIF 
            ENDDO 
            ENDDO 
            mole_num_act = mole_num_act + 1 
         ENDIF 
         ENDDO 
!     --End of loop over secondary atoms                                
!                                                                       
!     --Makesure that the molecules are not split into different        
!       unit cells                                                      
         CALL first_mole (mole_st) 
      ENDIF 
!                                                                       
      END SUBROUTINE mole_insert                    
!********************************************************************** 
      SUBROUTINE mole_insert_current (iatom, imole) 
!-                                                                      
!     Inserts the last atom into the molecule list as last atom of      
!     the specified molecule.                                           
!+                                                                      
      USE discus_allocate_appl_mod
      USE crystal_mod
      USE molecule_mod 
      USE prop_para_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER, INTENT(IN) :: iatom 
      INTEGER, INTENT(IN) :: imole 
      INTEGER :: n_gene
      INTEGER :: n_symm
      INTEGER :: n_mole
      INTEGER :: n_type
      INTEGER :: n_atom
      LOGICAL :: need_alloc = .false.
!                                                                       
      INTEGER i 
!                                                                       
!     Move the content of all molecules after "imole" one down the list 
!     mole_num_atom : total number of atoms in molecules                
!     mole_num_acur : number of atoms in molecules including "imole"    
!                                                                       
      mole_num_atom = mole_off (mole_num_mole) + mole_len (             &
      mole_num_mole)                                                    
!                                                                       
!     If necessary, create new molecule                                 
!                                                                       
      need_alloc = .false.
      n_gene = MAX( 1, MOLE_MAX_GENE)
      n_symm = MAX( 1, MOLE_MAX_SYMM)
      n_mole =         MOLE_MAX_MOLE
      n_type =         MOLE_MAX_TYPE
      n_atom =         MOLE_MAX_ATOM
      IF (imole > MOLE_MAX_MOLE ) THEN
         n_mole = MOLE_MAX_MOLE + 20
         need_alloc = .true.
      ENDIF
      IF (iatom > MOLE_MAX_ATOM ) THEN
         n_atom = MOLE_MAX_ATOM + 200
         need_alloc = .true.
      ENDIF
      IF ( need_alloc ) THEN
         call alloc_molecule(n_gene, n_symm, n_mole, n_type, n_atom)
      ENDIF
      IF (imole.gt.mole_num_mole) then 
         IF (imole.le.MOLE_MAX_MOLE) then 
            mole_num_mole = mole_num_mole+1 
            mole_len (imole) = 0 
            mole_off (imole) = mole_num_atom 
            mole_type (imole) = mole_num_type 
            mole_char (imole) = mole_char (imole-1) 
            mole_dens (imole) = mole_dens (imole-1) 
!            mole_biso (imole) = mole_biso (imole-1) 
         ELSE 
            ier_num = - 65 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
      ENDIF 
!                                                                       
      mole_num_acur = mole_off (imole) + mole_len (imole) 
!                                                                       
      IF (mole_num_atom + 1.le.MOLE_MAX_ATOM) then 
         DO i = mole_num_atom, mole_num_acur + 1, - 1 
         mole_cont (i + 1) = mole_cont (i) 
         ENDDO 
      ELSE 
         ier_num = - 74 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
!                                                                       
!     move the offset of all molecules after "imole" one down           
!                                                                       
      DO i = imole+1, mole_num_mole 
      mole_off (i) = mole_off (i) + 1 
      ENDDO 
!                                                                       
!     insert atom "iatom" at the end of molecule "imole"                
!                                                                       
      IF (mole_num_acur + 1.le.MOLE_MAX_ATOM) then 
         mole_cont (mole_num_acur + 1) = iatom 
         mole_len (imole) = mole_len (imole) + 1 
         mole_num_acur = mole_num_acur + 1 
      ELSE 
         ier_num = - 74 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
      cr_prop(iatom) = ibset(cr_prop(iatom),PROP_MOLECULE)
      cr_mole(iatom) = imole
!                                                                       
      END SUBROUTINE mole_insert_current            
!********************************************************************** 
      SUBROUTINE first_mole (mole_st) 
!-                                                                      
!     Moves atoms by +- one unit cell to keep molecules concatenated    
!     The algorithm of this subroutine is based on the assumption,      
!     that the first atom is on the highest point of symmetry of the    
!     molecule, i.e. is not copied by any of the symmetry operations    
!     or generators that make up the molecule symmetry. This has the    
!     desired effect that the bond distance of all symmetry             
!     equivalent atoms to this first atom is constant and can be taken  
!     as a reference. If a symmetry operation of the space group        
!     moves an atom out of the current unit cell, one can move this     
!     atom back by integer unit cell vectors until the bond distance    
!     is correct again.                                                 
!     If the molecule does not include an atom on the highest point     
!     of the molecule symmetry, you must insert a "void" on this site.  
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE metric_mod
      USE molecule_mod 
!
      USE prompt_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      LOGICAL lspace 
      PARAMETER (lspace = .true.) 
!                                                                       
      INTEGER mole_st 
      INTEGER i, j, k1, k2, k3 
      INTEGER k1u, k1o, k2u, k2o, k3u, k3o 
      REAL d, dd 
      REAL u (3), v (3) 
      REAL x, y, z 
!                                                                       
!     REAL do_blen 
!                                                                       
!.......calculate metric and reciprocal metric tensor,reciprocal lattice
!       constants and permutation tensors                               
      CALL setup_lattice (cr_a0, cr_ar, cr_eps, cr_gten, cr_reps,       &
      cr_rten, cr_win, cr_wrez, cr_v, cr_vr, .false., cr_gmat, cr_fmat, &
      cr_cartesian)                                                     
!                                                                       
!     Only one atom in the current molecule, return immediately         
!                                                                       
!DBG                                                                    
!DBG      write (output_io,*) 'mole_len(mole_num_curr)',                
!DBG     &                                mole_len(mole_num_curr)       
!DBG      write (output_io,*) 'mole_st                ',mole_st         
                                                                        
      IF (mole_len (mole_num_curr) .eq.1) return 
      IF (mole_st.eq.1) mole_st = 2 
      IF (mole_len (mole_num_curr) .lt.mole_st) return 
!                                                                       
      u (1) = cr_pos (1, mole_cont (mole_off (mole_num_curr) + 1) ) 
      u (2) = cr_pos (2, mole_cont (mole_off (mole_num_curr) + 1) ) 
      u (3) = cr_pos (3, mole_cont (mole_off (mole_num_curr) + 1) ) 
!                                                                       
      v (1) = cr_pos (1, mole_cont (mole_off (mole_num_curr) + mole_st) &
      )                                                                 
      v (2) = cr_pos (2, mole_cont (mole_off (mole_num_curr) + mole_st) &
      )                                                                 
      v (3) = cr_pos (3, mole_cont (mole_off (mole_num_curr) + mole_st) &
      )                                                                 
!                                                                       
      d = do_blen (lspace, u, v) 
!DBG                                                                    
!DBG      write (output_io,5555) u,v,d                                  
!DBG5555      format(3f10.4,2x,3f10.4,2x,f12.4/)                        
!                                                                       
!     Loop over all molecules from current to last                      
!                                                                       
      DO i = mole_num_curr, mole_num_mole 
!DBG                                                                    
!DBG      write (output_io,*) ' Molecule : ',i                          
      u (1) = cr_pos (1, mole_cont (mole_off (i) + 1) ) 
      u (2) = cr_pos (2, mole_cont (mole_off (i) + 1) ) 
      u (3) = cr_pos (3, mole_cont (mole_off (i) + 1) ) 
!                                                                       
!     Loop over all atoms on places mole_st and higher                  
!                                                                       
      DO j = mole_st, mole_len (i) 
!                                                                       
!     ----Determine proper limits for the loop over neighboring         
!          unit cells                                                   
!                                                                       
      x = cr_pos (1, mole_cont (mole_off (i) + j) ) 
      IF ( - 1.0.lt.x.and.x.lt.0.or.x.ge.1.0) then 
         k1u = - int (x) 
         k1o = k1u + 1 
      ELSE 
         k1u = - int (x) - 1 
         k1o = k1u + 1 
      ENDIF 
      y = cr_pos (2, mole_cont (mole_off (i) + j) ) 
      IF ( - 1.0.lt.y.and.y.lt.0.or.y.ge.1.0) then 
         k2u = - int (y) 
         k2o = k2u + 1 
      ELSE 
         k2u = - int (y) - 1 
         k2o = k2u + 1 
      ENDIF 
      z = cr_pos (3, mole_cont (mole_off (i) + j) ) 
      IF ( - 1.0.lt.z.and.z.lt.0.or.z.ge.1.0) then 
         k3u = - int (z) 
         k3o = k3u + 1 
      ELSE 
         k3u = - int (z) - 1 
         k3o = k3u + 1 
      ENDIF 
!DBG                                                                    
!DBG      write (output_io,*) ' x,y,z    ',x,y,z                        
!DBG      write (output_io,*) 'loop over ',k1u,k1o,k2u,k2o,k3u,k3o      
!                                                                       
!     ----perform loop over next unit cells                             
!                                                                       
      DO k1 = 2, - 2, - 1 
      DO k2 = 2, - 2, - 1 
      DO k3 = 2, - 2, - 1 
      v (1) = cr_pos (1, mole_cont (mole_off (i) + j) ) + float (k1) 
      v (2) = cr_pos (2, mole_cont (mole_off (i) + j) ) + float (k2) 
      v (3) = cr_pos (3, mole_cont (mole_off (i) + j) ) + float (k3) 
      dd = do_blen (lspace, u, v) 
!DBG                                                                    
!DBG      write (output_io,5556) u,j,mole_off(i)+j,v,dd                 
!DBG5556      format(3f10.4,2(2x,i2),2x,3f10.4,2x,f12.4)                
      IF (abs (d-dd) .lt.0.01) then 
         GOTO 10 
      ENDIF 
      ENDDO 
      ENDDO 
      ENDDO 
!                                                                       
      ier_num = - 83 
      ier_typ = ER_APPL 
      RETURN 
!                                                                       
   10 CONTINUE 
      cr_pos (1, mole_cont (mole_off (i) + j) ) = v (1) 
      cr_pos (2, mole_cont (mole_off (i) + j) ) = v (2) 
      cr_pos (3, mole_cont (mole_off (i) + j) ) = v (3) 
!DBG                                                                    
!DBG      write (output_io,5557) k1,k2,k3                               
!DBG5557      format('loesung fuer ',3i3)                               
!                                                                       
      ENDDO 
      ENDDO 
!                                                                       
      END SUBROUTINE first_mole                     
!********************************************************************** 
      SUBROUTINE mole_firstcell 
!-                                                                      
!     Moves molecules whose first atom is outside the unit cell into    
!     the first unit cell.                                              
!+                                                                      
      USE discus_config_mod 
      USE molecule_mod 
      USE crystal_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      REAL eps 
      PARAMETER (eps = 1.e-5) 
!                                                                       
      INTEGER i, j, k 
      REAL x, x1 
!                                                                       
!     Loop over all molecules from current to last                      
!                                                                       
      DO j = mole_num_curr, mole_num_mole 
      DO i = 1, 3 
      x = cr_pos (i, mole_cont (mole_off (j) + 1) ) 
      x1 = float (int (x) ) 
      IF (x - x1.lt. - eps) x1 = x1 - 1 
      IF (abs (x1) .gt.eps) then 
!                                                                       
!     ------Loop over all atoms in molecule j                           
!                                                                       
         DO k = 1, mole_len (j) 
         cr_pos (i, mole_cont (mole_off (j) + k) ) = cr_pos (i,         &
         mole_cont (mole_off (j) + k) ) - x1                            
         ENDDO 
      ENDIF 
      ENDDO 
      ENDDO 
      END SUBROUTINE mole_firstcell                 
!********************************************************************** 
      SUBROUTINE firstcell (y, idim) 
!-                                                                      
!     truncates atomic position to fractal position of                  
!     0.0 <= x < 1                                                      
!+                                                                      
      IMPLICIT none 
!                                                                       
      INTEGER idim, i 
      REAL y (idim) 
!                                                                       
      DO i = 1, 3 
      y (i) = y (i) - float (int (y (i) ) ) 
      IF (y (i) .lt.0.0) y (i) = y (i) + 1 
      IF (y (i) .eq.1.0) y (i) = 0.0 
      ENDDO 
      END SUBROUTINE firstcell                      
!*****7**************************************************************** 
      SUBROUTINE setup_lattice (cr_a0, cr_ar, cr_eps, cr_gten, cr_reps, &
      cr_rten, cr_win, cr_wrez, cr_v, cr_vr, lout, cr_gmat, cr_fmat,    &
      cr_cartesian)                                                     
!-                                                                      
!     Updates the crystal lattice and symmetry information              
!+                                                                      
      USE trafo_mod
      IMPLICIT none 
!                                                                       
      INTEGER i 
      LOGICAL lout 
      LOGICAL cr_cartesian 
      REAL cr_a0 (3), cr_ar (3), cr_eps (3, 3, 3), cr_gten (3, 3) 
      REAL cr_reps (3, 3, 3), cr_rten (3, 3), cr_win (3), cr_wrez (3) 
      REAL cr_v, cr_vr, cr_gmat (3, 3), cr_fmat (3, 3) 
      REAL hkl (3) 
      REAL u (3) 
      REAL xc (3) 
      REAL yc (3) 
      REAL zc (3) 
      REAL dist 
!                                                                       
      CALL lattice (cr_a0, cr_ar, cr_eps, cr_gten, cr_reps, cr_rten,    &
      cr_win, cr_wrez, cr_v, cr_vr, lout)                               
      cr_cartesian = cr_a0 (1) .eq.1..and.cr_a0 (2) .eq.1..and.cr_a0 (3)&
      .eq.1..and.cr_win (1) .eq.90..and.cr_win (2) .eq.90..and.cr_win ( &
      3) .eq.90.                                                        
      DO i = 1, 3 
      hkl (i) = 0.0 
      u (i) = 0.0 
      ENDDO 
      hkl (3) = 1.0 
      CALL trafo (hkl, u, xc, yc, zc, cr_gmat, cr_fmat, dist, cr_eps,   &
      cr_gten, cr_rten)                                        
      CALL recip_symm 
!                                                                       
      END SUBROUTINE setup_lattice                  
!*****7**************************************************************** 
      SUBROUTINE recip_symm 
!-                                                                      
!     Creates the symmetry matrices in reciprocal space                 
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE recipro_mod 
      USE rmc_symm_mod
      USE tensors_mod
      IMPLICIT none 
!                                                                       
      INTEGER i, j, k 
      LOGICAL lfriedel_remove 
      LOGICAL lacentric 
      REAL a (3, 3) 
      REAL b (3, 3) 
      REAL h (4, rec_max_sym) 
!                                                                       
!     Set up a general reflection                                       
!                                                                       
      h (1, 1) = 0.1 
      h (2, 1) = 0.2 
      h (3, 1) = 0.3 
      h (4, 1) = 0.0 
      DO i = 1, 4 
      DO j = 1, 4 
      rec_sym (i, j, 1) = 0.0 
      ENDDO 
      rec_sym (i, i, 1) = 1.0 
      ENDDO 
!                                                                       
!     Create all symmetry operations in real space that create          
!     a different vector, ignoring the translational part. The          
!     translational part of each first new matrix is saved, since       
!     this is needed for proper phase assignment.                       
!     lfriedel_remove  Logical variable that signals to remove          
!                      symmetry operations -1                           
!     cr_acentric      logical variable that describes whether the      
!                      space group is acentric or not                   
!                                                                       
      lfriedel_remove = .false. 
      CALL rmc_symmetry (rec_n_sym, h, rec_sym, rec_max_sym,            &
      lfriedel_remove, lacentric)                                       
      cr_acentric = lacentric 
!                                                                       
!     Transform all symmetry operations into reciprocal space           
!                                                                       
      DO k = 1, rec_n_sym 
      DO i = 1, 3 
      DO j = 1, 3 
      a (i, j) = rec_sym (i, j, k) 
      ENDDO 
      ENDDO 
!                                                                       
!     --do transformation q = gSg*                                      
!                                                                       
      CALL matmulx (b, a, cr_rten) 
      CALL matmulx (a, cr_gten, b) 
      DO i = 1, 3 
      DO j = 1, 3 
      rec_sym (i, j, k) = a (i, j) 
      ENDDO 
      ENDDO 
!                                                                       
!DBG      write( *,1000) k                                              
!DBG          do i=1,4                                                  
!DBG            write (*,1010) (rec_sym(i,j,k),j=1,4)                   
!DBG          ENDDO                                                     
!DBG1000      format('**********************'/' Matrix Number ', i4)    
!DBG1010      format(3(f4.1,2x),f7.4)                                   
      ENDDO 
!                                                                       
      END SUBROUTINE recip_symm                     
!*****7*****************************************************************
      SUBROUTINE lattice (a0, ar, eps, gten, reps, rten, win, wrez, vol,&
      vr, lout)                                                         
!+                                                                      
!           Calculates lattice constants, metric and reciprocal metric  
!           tensor, permutation tensors and unit cell volume.           
!     It's done quite some old fashioned way, rather than calculating   
!     the direct metric tensor and its inverse.                         
!-                                                                      
!                                                                       
      USE errlist_mod 
      USE prompt_mod 
      USE wink_mod
!
      IMPLICIT none 
!                                                                       
      INTEGER i, i1, i2, j 
      LOGICAL lout 
      REAL a0 (3), ar (3), eps (3, 3, 3), gten (3, 3) 
      REAL reps (3, 3, 3), rten (3, 3) 
      REAL win (3), wrez (3), vol, vr 
      REAL cosa, cosb, cosg, cos1, cos2, sin1, sin2 
      REAL cosi, sind, cosd, acosd 
!                                                                       
      cosa = cosd (win (1) ) 
      cosb = cosd (win (2) ) 
      cosg = cosd (win (3) ) 
      vol = 1 - cosa * cosa - cosb * cosb - cosg * cosg 
      vol = vol + 2 * cosa * cosb * cosg 
!                                                                       
      IF (vol.gt.0.0) then 
         vol = sqrt (vol) * a0 (1) * a0 (2) * a0 (3) 
         vr = 1. / vol 
!                                                                       
!------ - calculate direct metric tensor                                
!                                                                       
         CALL tensor (gten, a0, win) 
!                                                                       
!------ - calculate reciprocal lattice constants                        
!                                                                       
         DO i = 1, 3 
         i1 = mod (i, 3) + 1 
         i2 = mod (i + 1, 3) + 1 
         ar (i) = a0 (i1) * a0 (i2) * sind (win (i) ) / vol 
         cos1 = cosd (win (i1) ) 
         cos2 = cosd (win (i2) ) 
         cosi = cosd (win (i) ) 
         sin1 = sind (win (i1) ) 
         sin2 = sind (win (i2) ) 
         wrez (i) = acosd ( (cos1 * cos2 - cosi) / (sin1 * sin2) ) 
         ENDDO 
!                                                                       
!------ - calculate reciprocal tensor                                   
!                                                                       
         CALL tensor (rten, ar, wrez) 
!                                                                       
!------ - calculate premutation tensors                                 
!                                                                       
         eps (1, 2, 3) = vol 
         eps (2, 3, 1) = vol 
         eps (3, 1, 2) = vol 
         eps (1, 3, 2) = - vol 
         eps (3, 2, 1) = - vol 
         eps (2, 1, 3) = - vol 
         reps (1, 2, 3) = vr 
         reps (2, 3, 1) = vr 
         reps (3, 1, 2) = vr 
         reps (1, 3, 2) = - vr 
         reps (3, 2, 1) = - vr 
         reps (2, 1, 3) = - vr 
!                                                                       
!------ - output ?                                                      
!                                                                       
         IF (lout) then 
            WRITE (output_io, 2001) (a0 (i), i = 1, 3), (win (i),       &
            i = 1, 3), vol                                              
            WRITE (output_io, 2002) ( (gten (i, j), j = 1, 3), i = 1, 3) 
            WRITE (output_io, 2003) (ar (i), i = 1, 3), (wrez (i),      &
            i = 1, 3), vr                                               
            WRITE (output_io, 2004) ( (rten (i, j), j = 1, 3), i = 1, 3) 
         ENDIF 
!                                                                       
      ELSE 
         ier_num = - 35 
         ier_typ = ER_APPL 
      ENDIF 
!                                                                       
 2001 FORMAT     ( ' Lattice constants :'                               &
     &           ,/,4x,'a',10x,'b',10x,'c', 9x,                         &
     &           'alpha',6x,'beta',7x,'gamma', 6x,'volume',             &
     &           /,6(2X,F9.5),2X,G12.6)                                 
 2002 FORMAT     (/' Metric Tensor     :'/(3(' ',3(2X,F11.5)/))) 
 2003 FORMAT     ( ' Reciprocal Lattice constants :'                    &
     &           ,/,4x,'a*', 9x,'b*', 9x,'c*', 8x,                      &
     &           'alpha*',5x,'beta*',6x,'gamma*', 5x,'volume',          &
     &           /,6(2X,F9.5),2X,G12.6)                                 
 2004 FORMAT     (/' Reciprocal metric tensor     : '/                  &
     &            (3(' ',3(2X,F11.5)/)))                                
      END SUBROUTINE lattice                        
!*****7*****************************************************************
      SUBROUTINE tensor (ten, vec, win) 
!+                                                                      
!     Calculates the metric tensor. Works both for direct and           
!     reciprocal metric tensor.                                         
!-                                                                      
      IMPLICIT none 
!                                                                       
      INTEGER idim 
      PARAMETER (idim = 3) 
!                                                                       
      REAL ten (idim, idim), vec (idim), win (idim) 
      REAL cosd 
      INTEGER i, j 
!                                                                       
      DO i = 1, idim 
      DO j = 1, idim 
      IF (i.ne.j) then 
         ten (i, j) = vec (i) * vec (j) * cosd (win (6 - (i + j) ) ) 
      ELSE 
         ten (i, j) = vec (i) * vec (j) 
      ENDIF 
      ENDDO 
      ENDDO 
      END SUBROUTINE tensor                         
!*****7**************************************************************** 
END MODULE spcgr_apply
