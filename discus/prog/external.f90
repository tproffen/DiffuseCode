MODULE external_four
!
CONTAINS
      SUBROUTINE four_external 
!-                                                                      
!     Calculates the Fourier Transformation using externally            
!     calculated structure factors                                      
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE diffuse_mod 
      USE fourier_sup
      USE external_mod 
      USE metric_mod
      USE molecule_mod 
      USE errlist_mod 
      USE wink_mod
      IMPLICIT none 
!                                                                       
      INTEGER i, j, k, l, m, mm 
      INTEGER ij 
      INTEGER iscat 
      REAL det 
!                                                                       
!     REAL skalpro 
!DBG                                                                    
      exte_orig (1) = 0 
      exte_orig (2) = 0 
      exte_orig (3) = 0 
      exte_scale = 0.01 
!                                                                       
      CALL four_layer 
!                                                                       
!------ zero some arrays                                                
!                                                                       
      DO i = 1, num (1) * num (2) 
      csf (i) = cmplx (0.0, 0.0, KIND=KIND(0.0D0)) 
      acsf (i) = cmplx (0.0, 0.0, KIND=KIND(0.0D0)) 
      dsi (i) = 0.0d0 
      ENDDO 
!                                                                       
!------ preset some tables,                                             
!                                                                       
      CALL four_cexpt 
      CALL four_stltab 
      IF (ier_num.ne.0) return 
      CALL four_formtab 
!                                                                       
!     Loop over all molecules in the crystal                            
!                                                                       
      DO l = 1, mole_num_mole 
!                                                                       
      m = mole_cont (mole_off (l) + 1) 
      DO j = 1, 3 
      DO k = 1, 3 
      mm = mole_cont (mole_off (l) + 1 + k) 
      exte_mat (j, k) = cr_pos (j, mm) - cr_pos (j, m) 
      exte_rmat (k, j) = cr_pos (j, mm) - cr_pos (j, m) 
!DBG_RBN                                                                
!DBG_RBN      THESE TWO LINES ARE NO LONGER NEEDED !?!?!                
!DBG            mm = mole_cont(mole_off(l)+4+k)                         
!DBG            exte_rot(j,k) = cr_pos(j,mm) - cr_pos(j,m)              
      ENDDO 
      exte_mat (j, 4) = cr_pos (j, m) 
      exte_rmat (j, 4) = cr_pos (j, m) 
      ENDDO 
      det = exte_mat (1, 1) * (exte_mat (2, 2) * exte_mat (3, 3)        &
      - exte_mat (2, 3) * exte_mat (3, 2) ) + exte_mat (2, 1) * (       &
      exte_mat (3, 2) * exte_mat (1, 3) - exte_mat (1, 2) * exte_mat (3,&
      3) ) + exte_mat (3, 1) * (exte_mat (1, 2) * exte_mat (2, 3)       &
      - exte_mat (1, 3) * exte_mat (2, 2) )                             
!DBG    write (*,3000) 'Direct space ',((exte_mat  (i,j),j=1,4),i=1,3)  
!DBG    write (*,3000) 'Final Matrix ',((exte_rmat (i,j),j=1,4),i=1,3)  
!WORK          call external_header(mole_type(l))                       
      DO i = 1, num (1) * num (2) 
      tcsf (i) = cmplx (0.0, 0.0, KIND=KIND(0.0D0)) 
      ENDDO 
!                                                                       
!     -- Call the specialised subroutines for standard shapes           
!                                                                       
      IF (mole_char (mole_type (l) ) .eq.MOLE_CUBE) then 
         CALL external_cube (det, l) 
      ELSEIF (mole_char (mole_type (l) ) .eq.MOLE_CYLINDER) then 
         CALL external_cylinder (det, l) 
      ELSEIF (mole_char (mole_type (l) ) .eq.MOLE_SPHERE) then 
         CALL external_sphere (det, l) 
      ELSEIF (mole_char (mole_type (l) ) .eq.MOLE_EDGE) then 
         CALL external_edge (det, l) 
      ENDIF 
!                                                                       
!------ Now we multiply with formfactor                                 
!                                                                       
      iscat = cr_iscat (m) 
      IF (iscat.gt.0) then 
         DO i = 1, num (1) * num (2) 
         csf (i) = csf (i) + tcsf (i) * cfact (istl (i), iscat) 
         ENDDO 
      ENDIF 
      ENDDO 
!                                                                       
!     Finally calclutate DSI                                            
!                                                                       
      DO ij = 1, num (1) * num (2) 
      dsi (ij) = real (csf (ij) * conjg (csf (ij) ) ) 
      ENDDO 
!                                                                       
!DBG3000      format( a30,' : ',4(2x,f9.4)/                             
!DBG     &                2( 30x,' : ',4(2x,f9.4)/))                    
!                                                                       
      END SUBROUTINE four_external                  
!*****7**************************************************************** 
      SUBROUTINE external_cube (det, number) 
!-                                                                      
!     Calculates the Fourier Transformation using externally            
!     calculated structure factors for a cube shaped object             
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE diffuse_mod 
      USE external_mod 
      USE metric_mod
      USE molecule_mod 
      USE errlist_mod 
      USE trig_degree_mod
      USE wink_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER number 
      INTEGER i, j, k 
      INTEGER ii, jj 
      INTEGER ij 
      REAL h (3), hh (3) 
      REAL det 
      REAL phase 
      REAL sf, sfstart 
!                                                                       
!     REAL skalpro 
!     REAL cosd, sind 
!                                                                       
!      -- Loop over all points in reciprocal space                      
!                                                                       
      sfstart = 8. * cr_v * abs (det) * mole_dens (number) 
      ij = 0 
      DO i = 1, num (1) 
      DO j = 1, num (2) 
      ij = ij + 1 
      DO k = 1, 3 
      h (k) = REAL(xm (k) + uin (k) * REAL(i - 1) + vin (k) * REAL(j - 1) )
      ENDDO 
!                                                                       
!     ---- Transform vector                                             
!                                                                       
      phase = 0.0 
      DO ii = 1, 3 
      hh (ii) = 0.0 
      DO jj = 1, 3 
      hh (ii) = hh (ii) + exte_rmat (ii, jj) * h (jj) 
      ENDDO 
      phase = phase+h (ii) * exte_rmat (ii, 4) 
      ENDDO 
      phase = phase * 360.0 
      sf = sfstart 
      DO ii = 1, 3 
      IF (hh (ii) .ne.0) then 
         sf = sf * sin (REAL(zpi) * hh (ii) ) / (REAL(zpi) * hh (ii) ) 
      ENDIF 
      ENDDO 
      tcsf (ij) = tcsf (ij) + cmplx (sf * cosd (phase), sf * sind (     &
      phase) )                                                          
      ENDDO 
      ENDDO 
!                                                                       
      END SUBROUTINE external_cube                  
!*****7**************************************************************** 
      SUBROUTINE external_cylinder (det, number) 
!-                                                                      
!     Calculates the Fourier Transformation using externally            
!     calculated structure factors for a cylindrical object             
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE diffuse_mod 
      USE external_mod 
      USE metric_mod
      USE molecule_mod 
      USE errlist_mod 
      USE trig_degree_mod
      USE wink_mod
      IMPLICIT none 
!                                                                       
      INTEGER number 
      INTEGER i, j, k 
      INTEGER ii, jj 
      INTEGER ij 
      REAL h (3), hh (3), r (3), z (3) 
      REAL det 
      REAL phase 
      REAL sf, dr, dz, sfstart 
      REAL qr, qz 
!                                                                       
      REAL cartesian (3, 3) 
!                                                                       
!     REAL skalpro 
!     REAL cosd, sind 
      REAL bessj1 
!                                                                       
      DATA cartesian / 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 / 
!                                                                       
!      -- Loop over all points in reciprocal space                      
!                                                                       
      sfstart = REAL(zpi) * cr_v * abs (det) * mole_dens (number) 
      ij = 0 
      DO i = 1, num (1) 
      DO j = 1, num (2) 
      ij = ij + 1 
      DO k = 1, 3 
      h (k) = REAL(xm (k) + uin (k) * REAL(i - 1) + vin (k) * REAL(j - 1) )
      ENDDO 
!                                                                       
!     ---- Transform vector                                             
!                                                                       
      phase = 0.0 
      DO ii = 1, 3 
      hh (ii) = 0.0 
      DO jj = 1, 3 
      hh (ii) = hh (ii) + exte_rmat (ii, jj) * h (jj) 
      ENDDO 
      phase = phase+h (ii) * exte_rmat (ii, 4) 
      r (ii) = hh (ii) 
      z (ii) = hh (ii) 
      ENDDO 
      r (3) = 0.0 
      z (1) = 0.0 
      z (2) = 0.0 
      phase = phase * 360.0 
!                                                                       
      dz = sqrt (skalpro (z, z, cartesian) ) 
      dr = sqrt (skalpro (r, r, cartesian) ) 
      qz = REAL(zpi) * dz 
      qr = REAL(zpi) * dr 
      sf = sfstart 
      IF (dr.gt.2e-3) then 
         sf = sf * 2. * bessj1 (qr) / qr 
      ENDIF 
      IF (dz.gt.2e-3) then 
         sf = sf * sin (qz) / qz 
      ENDIF 
      tcsf (ij) = tcsf (ij) + cmplx (sf * cosd (phase), sf * sind (     &
      phase) )                                                          
      ENDDO 
      ENDDO 
!                                                                       
!                                                                       
      END SUBROUTINE external_cylinder              
!*****7**************************************************************** 
      SUBROUTINE external_sphere (det, number) 
!-                                                                      
!     Calculates the Fourier Transformation using externally            
!     calculated structure factors for a spherical object               
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE diffuse_mod 
      USE external_mod 
      USE metric_mod
      USE molecule_mod 
      USE errlist_mod 
      USE trig_degree_mod
      USE wink_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER number 
      INTEGER i, j, k 
      INTEGER ii, jj 
      INTEGER ij 
      REAL h (3), hh (3) 
      REAL det 
      REAL phase 
      REAL ds, sf, sfstart 
      REAL qr 
!                                                                       
      REAL cartesian (3, 3) 
!                                                                       
!     REAL skalpro 
!     REAL cosd, sind 
!                                                                       
      DATA cartesian / 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 / 
!                                                                       
!      -- Loop over all points in reciprocal space                      
!                                                                       
      sfstart = 4. / 3. * REAL(pi) * cr_v * abs (det) * mole_dens (number) 
      ij = 0 
      DO i = 1, num (1) 
      DO j = 1, num (2) 
      ij = ij + 1 
      DO k = 1, 3 
      h (k) = REAL(xm (k) + uin (k) * REAL(i - 1) + vin (k) * REAL(j - 1) )
      ENDDO 
!                                                                       
!     ---- Transform vector                                             
      phase = 0.0 
      DO ii = 1, 3 
      hh (ii) = 0.0 
      DO jj = 1, 3 
      hh (ii) = hh (ii) + exte_rmat (ii, jj) * h (jj) 
      ENDDO 
      phase = phase+h (ii) * exte_rmat (ii, 4) 
      ENDDO 
      phase = phase * 360.0 
!                                                                       
      ds = sqrt (skalpro (hh, hh, cartesian) ) 
      sf = sfstart 
      qr = REAL(zpi) * ds 
      IF (ds.gt.2e-3) then 
         sf = sfstart * 3 * (sin (qr) - qr * cos (qr) ) / (qr) **3 
      ENDIF 
      tcsf (ij) = tcsf (ij) + cmplx (sf * cosd (phase), sf * sind (     &
      phase) )                                                          
      ENDDO 
      ENDDO 
!                                                                       
      END SUBROUTINE external_sphere                
!*****7**************************************************************** 
      SUBROUTINE external_edge (det, number) 
!-                                                                      
!     Calculates the Fourier Transformation using externally            
!     calculated structure factors for an edge shaped object            
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE diffuse_mod 
      USE external_mod 
      USE metric_mod
      USE molecule_mod 
      USE errlist_mod 
      USE prompt_mod 
      USE trig_degree_mod
      USE wink_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER number 
      INTEGER i, j, k 
      INTEGER ii, jj 
      INTEGER ij 
      REAL h (3), hh (3) 
      REAL det 
      REAL phase 
      REAL sf, sfstart 
!DBG_RBN                                                                
      REAL hmin (3) 
!                                                                       
!     REAL skalpro 
!     REAL cosd, sind 
!DBG_RBN                                                                
      hmin (1) = 0.0 
      hmin (2) = 0.0 
      hmin (3) = 0.0 
!                                                                       
!      -- Loop over all points in reciprocal space                      
!                                                                       
      sfstart = 8. * cr_v * abs (det) * mole_dens (number) 
!                                                                       
!-----      --- Find minimum h                                          
!                                                                       
      DO i = 1, 2 
      ij = ij + 1 
      DO k = 1, 3 
      h (k) = REAL(uin (k) * REAL(2 - i) + vin (k) * REAL( - 1 + i) )
      ENDDO 
!                                                                       
!     ---- Transform vector                                             
!                                                                       
      DO ii = 1, 3 
      hh (ii) = 0.0 
      DO jj = 1, 3 
      hh (ii) = hh (ii) + exte_rmat (ii, jj) * h (jj) 
      ENDDO 
      hmin (ii) = max (hmin (ii), abs (hh (ii) ) ) 
      ENDDO 
      ENDDO 
      sfstart = 1.0 
      ij = 0 
      DO i = 1, num (1) 
      DO j = 1, num (2) 
      ij = ij + 1 
      DO k = 1, 3 
      h (k) = REAL(xm (k) + uin (k) * REAL(i - 1) + vin (k) * REAL(j - 1) )
      ENDDO 
!                                                                       
!     ---- Transform vector                                             
!                                                                       
      phase = 0.0 
      DO ii = 1, 3 
      hh (ii) = 0.0 
      DO jj = 1, 3 
      hh (ii) = hh (ii) + exte_rmat (ii, jj) * h (jj) 
      ENDDO 
      phase = phase+h (ii) * exte_rmat (ii, 4) 
      hmin (ii) = min (hmin (ii), abs (hh (ii) ) ) 
      ENDDO 
      phase = phase * 360.0 
      sf = sfstart 
      IF (abs (hh (1) ) .gt.hmin (1) * 3.5.or.abs (hh (2) ) .gt.hmin (2)&
      * 3.5) then                                                       
         sf = 0.0 
      ELSE 
         IF (abs (hh (3) ) .gt.hmin (3) * 1.5) then 
            sf = - sf / (REAL(zpi) * hh (3) ) 
            tcsf (ij) = tcsf (ij) + cmplx (sf * cosd (phase+90),        &
            sf * sind (phase+90) )                                      
            WRITE ( output_io, * ) ' calculating line' 
         ELSE 
            sf = - 0.5 * sf 
            sf = - 1.0 * sf 
            WRITE ( output_io, * ) ' Calculated hh=0' 
            tcsf (ij) = tcsf (ij) + cmplx (sf * cosd (phase), sf * sind &
            (phase) )                                                   
         ENDIF 
      ENDIF 
      ENDDO 
      ENDDO 
!DBG_RBN                                                                
!     WRITE ( output_io , * ) 'hmin ', hmin (1) 
!     WRITE ( output_io , * ) 'hmin ', hmin (2) 
!     WRITE ( output_io , * ) 'hmin ', hmin (3) 
!                                                                       
      END SUBROUTINE external_edge                  
!*****7**************************************************************** 
      SUBROUTINE external_header (exte_num) 
!-                                                                      
!     Opens the external file and reads the header                      
!+                                                                      
      USE discus_config_mod 
      USE external_mod 
      USE molecule_mod 
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER iwr 
      PARAMETER (iwr = 7) 
!                                                                       
      INTEGER(1) header (1024) 
      CHARACTER(80) char_header (EXTE_HLINES) 
!                                                                       
      INTEGER exte_num 
      INTEGER i, j, is
      INTEGER irecl 
!                                                                       
      irecl = 512 
!                                                                       
      CALL no_error 
!                                                                       
      exte_filename = exte_names (exte_num) !(1:exte_length (i) ) 
!                                                                       
      CALL oeffne_external (iwr, exte_filename, irecl) 
!                                                                       
      IF (ier_num.ne.0) then 
         RETURN 
      ENDIF 
!                                                                       
!     Read ASCII Header into integer variable                           
!                                                                       
      DO i = 1, 2 
      READ (iwr, rec = i) (header (j), j = (i - 1) * 512 + 1, (i - 1)   &
      * 512 + 512)                                                      
      ENDDO 
!                                                                       
!     Transform integer Header into Character variable                  
!                                                                       
      DO i = 1, EXTE_HLINES 
      is = (i - 1) * 80 
      DO j = 1, 80 
      char_header (i) (j:j) = char (header (is + j) ) 
      ENDDO 
      ENDDO 
!                                                                       
!     Read character variable corresponding to ASCII Header             
!                                                                       
!     VERSION :                                                         
      READ (char_header (1), 3001) exte_version 
!                                                                       
!     HDRBLKS :                                                         
      READ (char_header (2), 3002) exte_hdrblks 
!                                                                       
!     TITLE :                                                           
      READ (char_header (3), 3003) exte_title (1:72) 
      READ (char_header (4), 3003) exte_title (73:144) 
!                                                                       
!     DIMENS :                                                          
      READ (char_header (5), 3004) exte_nrows, exte_ncols, exte_layer 
!                                                                       
!     NPIXELB:                                                          
      READ (char_header (6), 3005) exte_npixelb 
!                                                                       
!     WORDORD:                                                          
      READ (char_header (7), 3006) exte_wordord 
!                                                                       
!     LONGORD:                                                          
      READ (char_header (8), 3007) exte_longord 
!                                                                       
!     ORIGIN:                                                           
      READ (char_header (9), 3008) exte_orig 
!                                                                       
!     SCALE :                                                           
      READ (char_header (10), 3009) exte_scale 
!                                                                       
 3001 FORMAT    (8x,i10) 
 3002 FORMAT    (8x,i10) 
 3003 FORMAT    (8x,a72) 
 3004 FORMAT    (8x,3i10) 
 3005 FORMAT    (8x,i10) 
 3006 FORMAT    (8x,i10) 
 3007 FORMAT    (8x,i10) 
 3008 FORMAT    (8x,3i10) 
 3009 FORMAT    (8x,f14.7) 
!                                                                       
      END SUBROUTINE external_header                
!*****7**************************************************************** 
      SUBROUTINE oeffne_external (irw, outfile, irecl) 
!-                                                                      
!     Opens the external file for input/output                          
!+                                                                      
      USE errlist_mod 
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER irw, irecl 
      CHARACTER ( * ) outfile 
!                                                                       
      INTEGER ios 
!                                                                       
      ier_num = - 2 
      ier_typ = ER_IO 
!                                                                       
      OPEN (irw, file = outfile, access = 'direct', recl = irecl, form =&
      'unformatted', err = 999, iostat = ios)                           
      ier_num = 0 
      ier_typ = ER_NONE 
  999 CONTINUE 
      IF (ier_num.eq. - 2) then 
         WRITE ( *, 1900) outfile
         WRITE ( output_io, 2000) ios 
      ENDIF 
!                      
 1900 FORMAT    ('Could not open: ', a)
 2000 FORMAT    (' ****SYST****Operating System/Shell Error Number:',i5,&
     &                  '****')                                         
      END SUBROUTINE oeffne_external                
END MODULE external_four
