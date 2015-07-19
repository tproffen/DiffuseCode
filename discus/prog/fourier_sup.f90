MODULE fourier_sup
!
USE errlist_mod 
!
CONTAINS
!*****7*****************************************************************
   SUBROUTINE four_run 
!+                                                                      
!     Main routine to calculate the Fourier transform for the           
!     given plane in reciprocal space. Based on the program             
!     'diffuse' written by B.D. Butler.                                 
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE diffuse_mod 
      USE four_strucf_mod
!                                                                       
      USE prompt_mod 
      IMPLICIT none 
       
!                                                                       
      REAL ss, seknds, dnorm
      INTEGER lbeg (3), csize (3) 
      INTEGER iscat, nlot, ncell, i 
!                                                                       
      ier_num = 0 
!                                                                       
!------ preset some values                                              
!                                                                       
      ss = seknds (0.0) 
      CALL four_layer 
!                                                                       
!------ zero some arrays                                                
!                                                                       
      DO i = 1, num (1) * num (2) * num(3)
         csf (i) = cmplx (0.0, 0.0) 
         acsf (i) = cmplx (0.0, 0.0) 
         dsi (i) = 0.0d0 
      ENDDO 
!                                                                       
!------ preset some tables, calculate average structure                 
!                                                                       
      CALL four_cexpt 
      CALL four_stltab 
      IF (ier_num.ne.0) return 
      CALL four_formtab 
      CALL four_csize (cr_icc, csize, lperiod, ls_xyz) 
      CALL four_aver (ilots, fave) 
!                                                                       
!------ loop over crystal regions (lots)                                
!                                                                       
      loop_lots: DO nlot = 1, nlots 
         DO i = 1, num (1) * num (2) * num(3)
            csf (i) = cmplx (0.0, 0.0) 
         ENDDO 
!                                                                       
         CALL four_ranloc (csize, lbeg) 
         IF (four_log) then 
            IF (ilots.ne.LOT_OFF) then 
               WRITE (output_io, 2000) nlot, nlots, (lbeg (i), i = 1, 3) 
            ELSE 
               WRITE (output_io, 2010) nlot, nlots 
            ENDIF 
         ENDIF 
!                                                                       
!------ - loop over all different atom types                            
!                                                                       
         loop_atoms: DO iscat = 1, cr_nscat 
            CALL four_getatm (iscat, ilots, lbeg, ncell) 
            CALL four_strucf (iscat, .true.) 
!                                                                       
!------ --- Add this part of the structur factor to the total           
!                                                                       
            DO i = 1, num (1) * num (2) *num (3)
               csf (i) = csf (i) + tcsf (i) 
            ENDDO 
            IF (four_log) then 
               WRITE (output_io, 3000) cr_at_lis (iscat), nxat 
            ENDIF 
         ENDDO loop_atoms
!                                                                       
!------ - subtract average structure factor, add intensity              
!                                                                       
         DO i = 1, num (1) * num (2) *num (3)
            csf (i) = csf (i) - acsf (i) 
            dsi (i) = dsi (i) + real (csf (i) * conjg (csf (i) ) ) 
         ENDDO 
      ENDDO loop_lots
!                                                                       
!------ if we had lots, normalise the intensities                       
!                                                                       
      IF (nlots.ne.1) then 
         dnorm = float (cr_icc (1) * cr_icc (2) * cr_icc (3) ) / float (&
         ncell * nlots)                                                 
         dnorm = dnorm**2 
         DO i = 1, num (1) * num (2) *num (3)
         dsi (i) = dnorm * dsi (i) 
         ENDDO 
      ENDIF 
!                                                                       
      CALL four_qinfo 
      ss = seknds (ss) 
      IF (four_log) then 
         WRITE (output_io, 4000) ss 
      ENDIF 
!                                                                       
 2000 FORMAT     (/,' Lot # ',I4,'/',I4,' : starting at unit cell (',   &
     &                       3(I3,1X),')')                              
 2010 FORMAT     (/,' Lot # ',I4,'/',I4,' : complete crystal') 
 3000 FORMAT (  '                   Atom typ = ',A4,7X,'(# ',I9,' )') 
 4000 FORMAT     (/,' Elapsed time    : ',G12.6,' sec') 
      END SUBROUTINE four_run                       
!*****7*****************************************************************
      SUBROUTINE four_aver (lots, ave) 
!+                                                                      
!     This routine calculates the average structure factor <F>.         
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE celltoindex_mod
      USE diffuse_mod 
      USE four_strucf_mod
!     USE modify_mod
      USE random_mod
!                                                                       
      USE prompt_mod 
      IMPLICIT none 
       
!
      INTEGER, INTENT(IN) :: lots
      REAL   , INTENT(IN) :: ave 
!                                                                       
      REAL ran1, norm
      INTEGER isite, iatom, iscat, icell (3) 
      INTEGER scell, ncell, j, ii, jj, kk 
      LOGICAL sel 
!                                                                       
      ave_is_zero: IF (ave.eq.0.0) then 
         RETURN 
      ELSE ave_is_zero  
         IF (four_log) then 
            WRITE (output_io, 1000) 100.0 * ave 
         ENDIF 
         scell = cr_icc (1) * cr_icc (2) * cr_icc (3) 
         ncell = 0 
!                                                                       
!------ - Loop over all unit cells                                      
!                                                                       
         cell_x: DO ii = 1, cr_icc (1) 
         cell_y: DO jj = 1, cr_icc (2) 
         cell_z: DO kk = 1, cr_icc (3) 
         icell (1) = ii 
         icell (2) = jj 
         icell (3) = kk 
!                                                                       
!------ --- get only 'ave'% of those unit cells and compute average     
!------ --- unit cell                                                   
!                                                                       
         sel = (ran1 (idum) .le.ave) 
         IF (sel) then 
            ncell = ncell + 1 
!                                                                       
!------ ----- Loop over all atom types                                  
!                                                                       
            DO iscat = 1, cr_nscat 
               nxat = 0 
!                                                                       
!------ ------- Loop over all sites within unit cell                    
!                                                                       
               DO isite = 1, cr_ncatoms 
                  CALL celltoindex (icell, isite, iatom) 
                  IF (cr_iscat (iatom) .eq.iscat) then 
                     nxat = nxat + 1 
                     DO j = 1, 3 
                     xat (nxat, j) = cr_pos (j, iatom) - float (icell (j)     &
                     - 1) - cr_dim0 (j, 1)                                    
                     ENDDO 
                  ENDIF 
               ENDDO 
               CALL four_strucf (iscat, .true.) 
               DO j = 1, num (1) * num (2) *num (3)
                  acsf (j) = acsf (j) + tcsf (j) 
               ENDDO 
            ENDDO 
         ENDIF 
         ENDDO cell_z
         ENDDO cell_y
         ENDDO cell_x
!                                                                       
!------ - now compute the interference function of the lot shape        
!                                                                       
         CALL four_getav (lots) 
         CALL four_strucf (0, .false.) 
         IF(ncell >0) THEN
            norm = 1.0 / ncell 
         ELSE
            norm = 0.0
            ier_num = +1
            ier_typ = ER_FOUR
            ier_msg(1) = 'Does the crystal consist of just 1 unit cell?'
            ier_msg(2) = 'Increase the percentage for >set aver<'
         ENDIF
         DO j = 1, num (1) * num (2) * num (3)
            acsf (j) = acsf (j) * tcsf (j) * cmplx ( norm, 0.0) 
         ENDDO 
!                                                                       
!------ - write how much of the crystal we actually used                
!                                                                       
         IF (four_log) then 
            WRITE (output_io, 2000) (float (ncell) / float (scell) )    &
            * 100.0                                                     
         ENDIF 
!                                                                       
      ENDIF ave_is_zero
!                                                                       
 1000 FORMAT     (' Calculating <F> with ',F5.1,' % of the crystal ...') 
 2000 FORMAT     (' Used ',F5.1,' % of the crystal to calulate <F> ...') 
      END SUBROUTINE four_aver                      
!*****7*****************************************************************
      SUBROUTINE four_getatm (iscat, lots, lbeg, ncell) 
!+                                                                      
!     This routine creates an atom list of atoms of type 'iscat'        
!     which are within the current lot.                                 
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE celltoindex_mod
      USE diffuse_mod 
!     USE modify_mod
      IMPLICIT none 
!                                                                       
      INTEGER, INTENT(IN)  :: iscat
      INTEGER, INTENT(IN)  :: lots
      INTEGER, DIMENSION(3), INTENT(IN)  :: lbeg
      INTEGER, INTENT(OUT) :: ncell
!                                                                       
      REAL offset (3) 
      REAL x0 (3), xtest (3) 
      INTEGER cr_end 
      INTEGER icell (3), jcell (3) 
      INTEGER i, j, ir, ii, jj, kk, is, ia 
!                                                                       
      nxat = 0 
      ncell = 0 
      cr_end = cr_ncatoms * cr_icc (1) * cr_icc (2) * cr_icc (3) + 1
!                                                                       
!------ No lots required, return all atoms of type 'iscat'              
!                                                                       
      IF (lots.eq.LOT_OFF) then 
         ncell = cr_icc (1) * cr_icc (2) * cr_icc (3) 
         DO i = 1, cr_natoms 
         IF (cr_iscat (i) .eq.iscat) then 
            nxat = nxat + 1 
            DO j = 1, 3 
            xat (nxat, j) = cr_pos (j, i) 
            ENDDO 
         ENDIF 
         ENDDO 
!                                                                       
!------ Box shaped lot                                                  
!                                                                       
      ELSEIF (lots.eq.LOT_BOX) then 
         DO kk = 0, ls_xyz (3) - 1 
         icell (3) = kk + lbeg (3) 
         offset (3) = cr_dim0 (3, 1) + lbeg (3) - 1 
         IF (icell (3) .gt.cr_icc (3) ) then 
            icell (3) = icell (3) - cr_icc (3) 
            offset (3) = offset (3) - float (cr_icc (3) ) 
         ENDIF 
         DO jj = 0, ls_xyz (2) - 1 
         icell (2) = jj + lbeg (2) 
         offset (2) = cr_dim0 (2, 1) + lbeg (2) - 1 
         IF (icell (2) .gt.cr_icc (2) ) then 
            icell (2) = icell (2) - cr_icc (2) 
            offset (2) = offset (2) - float (cr_icc (2) ) 
         ENDIF 
         DO ii = 0, ls_xyz (1) - 1 
         icell (1) = ii + lbeg (1) 
         offset (1) = cr_dim0 (1, 1) + lbeg (1) - 1 
         IF (icell (1) .gt.cr_icc (1) ) then 
            icell (1) = icell (1) - cr_icc (1) 
            offset (1) = offset (1) - float (cr_icc (1) ) 
         ENDIF 
!                                                                       
         ncell = ncell + 1 
         DO is = 1, cr_ncatoms 
         CALL celltoindex (icell, is, ia) 
         IF (cr_iscat (ia) .eq.iscat) then 
            nxat = nxat + 1 
            DO j = 1, 3 
            xat (nxat, j) = cr_pos (j, ia) - offset (j) 
            ENDDO 
         ENDIF 
         ENDDO 
!                                                                       
         DO ir = cr_end, cr_natoms 
         DO j = 1, 3 
         jcell (j) = icell (j) + nint (cr_dim (j, 1) ) - 1 
         ENDDO 
         IF (int (cr_pos (1, ir) ) .eq.jcell (1) .and.int (cr_pos (2,   &
         ir) ) .eq.jcell (2) .and.int (cr_pos (3, ir) ) .eq.jcell (3)   &
         .and.cr_iscat (ir) .eq.iscat) then                             
            nxat = nxat + 1 
            DO j = 1, 3 
            xat (nxat, j) = cr_pos (j, ir) - offset (j) 
            ENDDO 
         ENDIF 
         ENDDO 
!                                                                       
         ENDDO 
         ENDDO 
         ENDDO 
!                                                                       
!------ Ellipsoid shaped lot                                            
!                                                                       
      ELSEIF (lots.eq.LOT_ELI) then 
         DO ii = 1, 3 
         x0 (ii) = float (ls_xyz (ii) ) / 2.0 
         ENDDO 
!                                                                       
         DO kk = 0, ls_xyz (3) - 1 
         icell (3) = kk + lbeg (3) 
         offset (3) = cr_dim0 (3, 1) + lbeg (3) - 1 
         xtest (3) = (float (kk) - x0 (3) + 0.5) **2 / x0 (3) **2 
         IF (icell (3) .gt.cr_icc (3) ) then 
            icell (3) = icell (3) - cr_icc (3) 
            offset (3) = offset (3) - float (cr_icc (3) ) 
         ENDIF 
         DO jj = 0, ls_xyz (2) - 1 
         icell (2) = jj + lbeg (2) 
         offset (2) = cr_dim0 (2, 1) + lbeg (2) - 1 
         xtest (2) = (float (jj) - x0 (2) + 0.5) **2 / x0 (2) **2 
         IF (icell (2) .gt.cr_icc (2) ) then 
            icell (2) = icell (2) - cr_icc (2) 
            offset (2) = offset (2) - float (cr_icc (2) ) 
         ENDIF 
         DO ii = 0, ls_xyz (1) - 1 
         icell (1) = ii + lbeg (1) 
         offset (1) = cr_dim0 (1, 1) + lbeg (1) - 1 
         xtest (1) = (float (ii) - x0 (1) + 0.5) **2 / x0 (1) **2 
         IF (icell (1) .gt.cr_icc (1) ) then 
            icell (1) = icell (1) - cr_icc (1) 
            offset (1) = offset (1) - float (cr_icc (1) ) 
         ENDIF 
!                                                                       
         IF ( (xtest (1) + xtest (2) + xtest (3) ) .le.1.0) then 
            ncell = ncell + 1 
            DO is = 1, cr_ncatoms 
            CALL celltoindex (icell, is, ia) 
            IF (cr_iscat (ia) .eq.iscat) then 
               nxat = nxat + 1 
               DO j = 1, 3 
               xat (nxat, j) = cr_pos (j, ia) - offset (j) 
               ENDDO 
            ENDIF 
            ENDDO 
!                                                                       
            DO ir = cr_end, cr_natoms 
            DO j = 1, 3 
            jcell (j) = icell (j) + nint (cr_dim (j, 1) ) - 1 
            ENDDO 
            IF (int (cr_pos (1, ir) ) .eq.jcell (1) .and.int (cr_pos (2,&
            ir) ) .eq.jcell (2) .and.int (cr_pos (3, ir) ) .eq.jcell (3)&
            .and.cr_iscat (ir) .eq.iscat) then                          
               nxat = nxat + 1 
               DO j = 1, 3 
               xat (nxat, j) = cr_pos (j, ir) - offset (j) 
               ENDDO 
            ENDIF 
            ENDDO 
         ENDIF 
!                                                                       
         ENDDO 
         ENDDO 
         ENDDO 
      ENDIF 
!                                                                       
      END SUBROUTINE four_getatm                    
!*****7*****************************************************************
      SUBROUTINE four_getav (lots) 
!+                                                                      
!     This routine computes the lattice inside the lot                  
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE diffuse_mod 
      IMPLICIT none 
!                                                                       
      INTEGER, INTENT(IN) :: lots
!                                                                       
      REAL xtest (3), x0 (3) 
      INTEGER ii, jj, kk 
!                                                                       
      nxat = 0 
!                                                                       
!------ No lots required, return complete lattice                       
!                                                                       
      IF (lots.eq.LOT_OFF) then 
         DO ii = 0, cr_icc (1) - 1 
         DO jj = 0, cr_icc (2) - 1 
         DO kk = 0, cr_icc (3) - 1 
         nxat = nxat + 1 
         xat (nxat, 1) = float (ii) + cr_dim0 (1, 1) 
         xat (nxat, 2) = float (jj) + cr_dim0 (2, 1) 
         xat (nxat, 3) = float (kk) + cr_dim0 (3, 1) 
         ENDDO 
         ENDDO 
         ENDDO 
!                                                                       
!------ Box shaped lot                                                  
!                                                                       
      ELSEIF (lots.eq.LOT_BOX) then 
         DO kk = 0, ls_xyz (3) - 1 
         DO jj = 0, ls_xyz (2) - 1 
         DO ii = 0, ls_xyz (1) - 1 
         nxat = nxat + 1 
         xat (nxat, 1) = float (ii) 
         xat (nxat, 2) = float (jj) 
         xat (nxat, 3) = float (kk) 
         ENDDO 
         ENDDO 
         ENDDO 
!                                                                       
!------ Ellipsoid shaped lot                                            
!                                                                       
      ELSEIF (lots.eq.LOT_ELI) then 
         DO ii = 1, 3 
         x0 (ii) = float (ls_xyz (ii) ) / 2.0 
         ENDDO 
!                                                                       
         DO kk = 0, ls_xyz (3) - 1 
         xtest (3) = (float (kk) - x0 (3) + 0.5) **2 / x0 (3) **2 
         DO jj = 0, ls_xyz (2) - 1 
         xtest (2) = (float (jj) - x0 (2) + 0.5) **2 / x0 (2) **2 
         DO ii = 0, ls_xyz (1) - 1 
         xtest (1) = (float (ii) - x0 (1) + 0.5) **2 / x0 (1) **2 
         IF ( (xtest (1) + xtest (2) + xtest (3) ) .le.1.0) then 
            nxat = nxat + 1 
            xat (nxat, 1) = float (ii) 
            xat (nxat, 2) = float (jj) 
            xat (nxat, 3) = float (kk) 
         ENDIF 
         ENDDO 
         ENDDO 
         ENDDO 
      ENDIF 
!                                                                       
      END SUBROUTINE four_getav                     
!*****7*****************************************************************
      SUBROUTINE four_ranloc (csize, lbeg) 
!+                                                                      
!     Returns a random cell from within the simulated crystal           
!     which is limited by 'csize'.                                      
!-                                                                      
      USE random_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER, DIMENSION(3), INTENT(IN)  :: csize
      INTEGER, DIMENSION(3), INTENT(OUT) :: lbeg
      REAL ran1 
!                                                                       
      lbeg (1) = int (ran1 (idum) * csize (1) ) + 1 
      lbeg (2) = int (ran1 (idum) * csize (2) ) + 1 
      lbeg (3) = int (ran1 (idum) * csize (3) ) + 1 
!                                                                       
      END SUBROUTINE four_ranloc                    
!*****7*****************************************************************
      SUBROUTINE four_csize (cr_icc, csize, lperiod, ls_xyz) 
!+                                                                      
!     Limits crystal size in case of no periodic boundary               
!     conditions.                                                       
!-                                                                      
      IMPLICIT none 
!                                                                       
      INTEGER, DIMENSION(3), INTENT(IN)  :: cr_icc
      INTEGER, DIMENSION(3), INTENT(OUT) :: csize
      LOGICAL,               INTENT(IN)  :: lperiod 
      INTEGER, DIMENSION(3), INTENT(IN)  :: ls_xyz
!                                                                       
      IF (lperiod) then 
         csize (1) = cr_icc (1) 
         csize (2) = cr_icc (2) 
         csize (3) = cr_icc (3) 
      ELSE 
         csize (1) = cr_icc (1) - ls_xyz (1) 
         csize (2) = cr_icc (2) - ls_xyz (2) 
         csize (3) = cr_icc (3) - ls_xyz (3) 
      ENDIF 
!                                                                       
      END SUBROUTINE four_csize                     
!*****7*****************************************************************
      SUBROUTINE four_layer 
!+                                                                      
!     This routine writes the corners of the plane to be calculated     
!     in the 'diffuse_mod.f90' variables.                                   
!-                                                                      
      USE discus_config_mod 
      USE diffuse_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER i 
!                                                                       
      DO i = 1, 3 
         xm  (i) = eck (i, 1) 
         uin (i) = vi  (i, 1) 
         vin (i) = vi  (i, 2) 
         win (i) = vi  (i, 3) 
      ENDDO 
!                                                                       
      DO i = 1, 3 
         num (i) = inc (i) 
      ENDDO 
!                                                                       
      END SUBROUTINE four_layer                     
!*****7*****************************************************************
      SUBROUTINE four_cexpt 
!+                                                                      
!     This routine initialises the complex exponent table and           
!     is called only at the first Fourier run.                          
!-                                                                      
      USE discus_config_mod 
      USE diffuse_mod 
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      REAL(dp) twopi, xmult, xarg 
      INTEGER i 
!                                                                       
      IF (.not.ffour) then 
         WRITE (output_io, 1000) 
!                                                                       
         twopi = 8.0d0 * datan (1.0d0) 
!                                                                       
         DO i = 0, MASK 
            xmult   = (dble (i) * 1.0d0) / dble (I2PI) 
            xarg    = twopi * xmult 
            cex (i) = cmplx (real( cos (xarg)), real( sin (xarg)) ) 
         ENDDO 
         ffour = .true. 
      ENDIF 
!                                                                       
 1000 FORMAT     (' Computing complex exponent table ...') 
      END SUBROUTINE four_cexpt                     
!*****7*****************************************************************
      SUBROUTINE four_stltab 
!+                                                                      
!     Sets up an integer array containing the corresponding integers    
!     to the formfactor table for each sin(theta)/lambda.               
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE diffuse_mod
      USE quad_mod 
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      REAL q2, h (3) 
      INTEGER i, j, k, l 
!                                                                       
      IF (four_log) then 
         WRITE (output_io, 1000) 
      ENDIF 
!                                                                       
      DO l = 1, num (3) 
      DO j = 1, num (2) 
      DO i = 1, num (1) 
      DO k = 1, 3 
      h (k) = xm (k) + uin (k) * float (i - 1) + vin (k) * float (j - 1) &
                     + win (k) * float (l - 1)
      ENDDO 
      q2 = quad (h, h, cr_rten) / 4.0 
!     k  = (i - 1) * num (2) + j 
      k  = (i - 1) * num (3)* num (2) + (j - 1) * num (3) + l 
!RBN_3D
!write(*,2000) i,j,l, k
!2000 format( 'ijl:k ', 3i4,' : ', i6)
      istl (k) = nint (sqrt (q2) * (1.0 / CFINC) ) 
!                                                                       
      IF (istl (k) .gt.CFPKT) then 
         ier_num = - 3 
         ier_typ = ER_FOUR 
         RETURN 
      ENDIF 
!                                                                       
      ENDDO 
      ENDDO 
      ENDDO 
!                                                                       
 1000 FORMAT     (' Computing sin(theta)/lambda table ...') 
      END SUBROUTINE four_stltab                    
!*****7*****************************************************************
      SUBROUTINE four_formtab 
!+                                                                      
!     This routine sets up the complex formfactor lookup table          
!     for all atom types. The range in sin(theta)/lambda is             
!     0 -> 2 in steps of 0.001. These values can be changed             
!     in the 'diffuse_mod.f90' file.                                        
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE diffuse_mod 
!                                                                       
      USE prompt_mod 
      IMPLICIT none 
       
!                                                                       
      REAL q2, sb, sf, sfp, sfpp 
      INTEGER iq, iscat 
!                                                                       
!     REAL form 
!                                                                       
      IF (four_log) then 
         WRITE (output_io, 1000) 
      ENDIF 
!                                                                       
      DO iscat = 1, cr_nscat 
      DO iq = 0, CFPKT 
      q2 = (float (iq) * CFINC) **2 
      sf = form (iscat, cr_scat, lxray, q2, diff_power) 
!                                                                       
      IF (ano) then 
         sfp = cr_delfr ( (iscat) ) 
         sfpp = cr_delfi ( (iscat) ) 
      ELSE 
         sfp = 0.0 
         sfpp = 0.0 
      ENDIF 
!                                                                       
      IF (ldbw) then 
         sb = exp ( - cr_dw ( (iscat) ) * q2) 
      ELSE 
         sb = 1.0 
      ENDIF 
!                                                                       
      cfact     (iq, iscat) = cmplx (sb * (sf + sfp), sb * sfpp) 
      cfact_pure(iq, iscat) = cmplx (     (sf + sfp),      sfpp) 
      ENDDO 
      ENDDO 
!                                                                       
 1000 FORMAT     (' Computing formfactor lookup table ...') 
      END SUBROUTINE four_formtab                   
!*****7*****************************************************************
      SUBROUTINE four_qinfo 
!+                                                                      
!     Gives information about max/min values for diffuse and            
!     Bragg scattering.                                                 
!-                                                                      
      USE discus_config_mod 
      USE diffuse_mod 
!                                                                     
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
      INTEGER i, j, k, l, nd 
      REAL h (3), dsum, dsum2 
      LOGICAL lbragg 
!                                                                       
      braggmax = - 9999.0 
      braggmin = 1E19 
      diffumax = - 9999.0 
      diffumin = 1E19 
!                                                                       
      nd = 0 
      dsum = 0.0 
      dsum2 = 0.0 
!                                                                       
      DO l = 1, num (3) 
      DO j = 1, num (2) 
      DO i = 1, num (1) 
      lbragg = .true. 
      DO k = 1, 3 
      h (k) = eck (k, 1) + vi (k, 1) * float (i - 1) + &
                           vi (k, 2) * float (j - 1) + &
                           vi (k, 3) * float (l - 1)
      lbragg = lbragg.and.amod (h (k), 1.0) .eq.0.0 
      ENDDO 
!     k = (i - 1) * num (2) + j 
      k = (i - 1) * num (3) * num (2) + (j -1) * num(3) + l 
      IF (lbragg) then 
         braggmax = max (braggmax, dsi (k) ) 
         braggmin = min (braggmin, dsi (k) ) 
      ELSE 
         diffumax = max (diffumax, dsi (k) ) 
         diffumin = min (diffumin, dsi (k) ) 
!                                                                       
         dsum = dsum + dsi (k) 
         dsum2 = dsum2 + dsi (k) **2 
         nd = nd+1 
      ENDIF 
      ENDDO 
      ENDDO 
      ENDDO 
!                                                                       
      IF (nd.ne.0) then 
         diffuave = dsum / float (nd) 
         diffusig = sqrt (dsum2 / float (nd) - diffuave**2) 
      ELSE 
         diffuave = 0.0 
         diffusig = 0.0 
      ENDIF 
!                                                                       
      IF (four_log) then 
         WRITE (output_io, 1000) 
         IF (braggmax.ge.0.0) then 
            WRITE (output_io, 1010) braggmin, braggmax 
         ENDIF 
         WRITE (output_io, 1020) diffumin, diffumax 
         WRITE (output_io, 1030) diffuave, diffusig 
      ENDIF 
!                                                                       
 1000 FORMAT     (/,' ') 
 1010 FORMAT     (  ' Bragg scat.     : ',G12.6,'  -> ',G12.6) 
 1020 FORMAT     (  ' Diffuse scat.   : ',G12.6,'  -> ',G12.6) 
 1030 FORMAT     (  '      Average    : ',G12.6,'  +- ',G12.6) 
      END SUBROUTINE four_qinfo                     
!*****7*****************************************************************
      REAL FUNCTION form (ll, scat, lxray, h2, power) 
!+                                                                      
!       calculates the form factor                                      
!-                                                                      
      USE discus_config_mod
      USE element_data_mod
      IMPLICIT none 
!                                                                       
      INTEGER, INTENT(IN) :: ll
      LOGICAL, INTENT(IN) :: lxray 
      REAL   , DIMENSION(11,0:MAXSCAT), INTENT(INOUT) :: scat ! (11, 0:maxscat) 
      REAL                            , INTENT(IN)    :: h2
      INTEGER, INTENT(IN) :: power
!
      INTEGER   :: i 
!                                                                       
      form = scat (1, ll) 
      IF (lxray) then 
         DO i = 1, power 
         form = form + scat (2 * i, ll) * exp ( - scat (2 * i + 1, ll)  &
         * h2)                                                          
         ENDDO 
      ENDIF 
      END FUNCTION form                             
!*****7*****************************************************************
      SUBROUTINE dlink (ano, lambda, rlambda, diff_radiation, &
                        diff_power) 
!-                                                                      
!     This routine reads wavelength symbols, wavelength values 
!     and atomic form factors from module "element_data_mod"
!     It replaces the old dlink and LAZY_PULVERIX
!+                                                                      
      USE discus_config_mod 
      USE charact_mod
      USE crystal_mod 
      USE element_data_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      LOGICAL             , INTENT(IN)   :: ano
      CHARACTER (LEN = * ), INTENT(IN)   :: lambda 
      REAL                , INTENT(OUT)  :: rlambda
      INTEGER             , INTENT(IN)   :: diff_radiation
      INTEGER             , INTENT(OUT)  :: diff_power
!
      INTEGER , PARAMETER  :: RAD_XRAY = 1
      INTEGER , PARAMETER  :: RAD_NEUT = 2
      INTEGER , PARAMETER  :: RAD_ELEC = 3
!
      CHARACTER (LEN = 4 ) :: element 
      INTEGER    :: i
      INTEGER    :: j
      REAL   , DIMENSION(1:11)  :: temp_scat  ! a1,b1,---a4,b4,c
      REAL   , DIMENSION(1:2)   :: temp_delf  ! delfr, delfi
      REAL                      :: temp_bcoh  ! b_choherent
!
!     CHARACTER(4) line 
!     INTEGER element, symwl, kodlp, i, j 
!     LOGICAL ano, lxray 
!     REAL fa (4, 8), fb (4, 8), fc (8), rlambda 
!     REAL delfr, delfi, fneu 
!     REAL wave 
!                                                                       
      ier_num = -77 
      ier_typ = ER_APPL 
!
      CALL get_wave ( lambda, rlambda, ier_num, ier_typ )
!                                                                       
      IF (ier_num.ne.0) RETURN 
!                                                                       
      any_element: IF (cr_nscat.gt.0) THEN 
         DO i = 1, cr_nscat 
         IF (cr_at_lis (i) /= 'XAXI' .AND. cr_at_lis (i) /= 'YAXI' .and. & 
             cr_at_lis (i) /= 'ZAXI') THEN                  
            IF (cr_scat_int (i) ) then 
               SELECTCASE(diff_radiation)
                  CASE(RAD_NEUT)        !  neutron scattering
                  IF (cr_scat_equ (i) ) then 
                     element =         cr_at_equ (i) 
                  ELSE 
                     element =         cr_at_lis (i) 
                  ENDIF 
!                                                                       
                  CALL symbf ( element, j)
                  IF ( j /= 0 ) THEN 
                     CALL get_scat_neut ( j, temp_bcoh )
                     cr_scat(:,i) = 0.0
                     cr_scat(1,i) = temp_bcoh
                     cr_delfr (i) = 0.0 
                     cr_delfi (i) = 0.0 
                     diff_power   = PER_RAD_POWER(diff_radiation)
                     ier_num = 0
                     ier_typ = ER_NONE 
                  ELSE
                     ier_num = -20 
                     ier_typ = ER_APPL 
                  ENDIF
                  CASE(RAD_XRAY)        !  Xray diffraction
                  IF (cr_scat_equ (i) ) then 
                     element =         cr_at_equ (i) 
                  ELSE 
                     element =         cr_at_lis (i) 
                  ENDIF 
!
                  ier_num = -20 
                  ier_typ = ER_APPL 
!                                                                       
                  CALL symbf ( element, j)
                  IF ( j /= 0 ) THEN 
                     CALL get_scat_xray ( j, temp_scat )
                     cr_scat(:,i) = 0.0
                     cr_scat(:,i) = temp_scat(:)   ! copy temp array into 1st column
                     cr_delfr (i) = 0.0   ! DEVELOPMENT !
                     cr_delfi (i) = 0.0   ! DEVELOPMENT !
                     IF (ano.and.cr_delf_int (i) ) then 
                        CALL get_scat_ano ( j, lambda, temp_delf )
                        cr_delfr (i) = temp_delf(1)
                        cr_delfi (i) = temp_delf(2)
                     ENDIF 
                     diff_power   = PER_RAD_POWER(diff_radiation)
                     ier_num = 0
                  ELSE
                     ier_typ = ER_NONE 
                     ier_typ = ER_APPL 
                  ENDIF
!
                  CASE(RAD_ELEC)        !  Electron diffraction
                  IF (cr_scat_equ (i) ) then 
                     element =         cr_at_equ (i) 
                  ELSE 
                     element =         cr_at_lis (i) 
                  ENDIF 
!
                  ier_num = -20 
                  ier_typ = ER_APPL 
!                                                                       
                  CALL symbf ( element, j)
                  IF ( j /= 0 ) THEN 
                     CALL get_scat_elec ( j, temp_scat)
                     cr_scat(:,i) = 0.0
                     cr_scat(:,i) = temp_scat(:)   ! copy temp array into 1st column
                     cr_delfr (i) = 0.0
                     cr_delfi (i) = 0.0
                     diff_power   = PER_RAD_POWER(diff_radiation)
                     ier_num = 0
                  ELSE
                     ier_typ = ER_NONE 
                     ier_typ = ER_APPL 
                  ENDIF
               END SELECT
!                                                                       
            ENDIF 
         ENDIF 
         ENDDO 
      ELSE any_element
         ier_num = - 21 
         ier_typ = ER_APPL 
      ENDIF any_element
!                                                                       
      END SUBROUTINE dlink                          
!*****7*****************************************************************
      SUBROUTINE calc_000 (rhkl) 
!-                                                                      
!     Calculates the value of F(rh,rk,rl)                               
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE diffuse_mod 
!
      USE param_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      REAL, DIMENSION(3), INTENT(IN) :: rhkl
!
      INTEGER i, j 
      INTEGER shel_inc (3) 
      REAL shel_eck (3, 4) 
      REAL shel_vi (3, 3) 
      COMPLEX shel_acsf 
      REAL shel_dsi 
      COMPLEX shel_tcsf 
!                                                                       
      DO i = 1, 3 
         shel_inc (i) = inc (i) 
      ENDDO 
!
      DO i = 1, 3 
         DO j = 1, 4 
            shel_eck (i, j) = eck (i, j) 
         ENDDO 
         DO j = 1, 3 
            shel_vi (i, j) = vi (i, j) 
         ENDDO 
      ENDDO 
!
      shel_tcsf = csf (1) 
      shel_acsf = acsf (1) 
      shel_dsi  = dsi (1) 
      inc (1) = 1 
      inc (2) = 1 
      inc (3) = 1 
      DO i = 1, 3 
         DO j = 1, 4 
            eck (i, j) = rhkl (i) 
         ENDDO 
         DO j = 1, 3 
            vi (i, j) = 0.0 
         ENDDO 
      ENDDO 
!
      four_log = .false. 
      CALL four_run 
      res_para (1) = real (csf (1) ) 
      res_para (2) = aimag (csf (1) ) 
      res_para (3) = res_para (1) / cr_icc (1) / cr_icc (2) / cr_icc (3) 
      res_para (4) = res_para (2) / cr_icc (1) / cr_icc (2) / cr_icc (3) 
      res_para (0) = 4 
      csf (1) = shel_tcsf 
      acsf(1) = shel_acsf 
      dsi (1) = shel_dsi 
!
      DO i = 1, 3 
         inc (i) = shel_inc (i) 
      ENDDO 
      DO i = 1, 3 
         DO j = 1, 4 
            eck (i, j) = shel_eck (i, j) 
         ENDDO 
         DO j = 1, 3 
            vi (i, j) = shel_vi (i, j) 
         ENDDO 
      ENDDO 
!                                                                       
      END SUBROUTINE calc_000                       
END MODULE fourier_sup
