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
      USE fourier_lmn_mod
!                                                                       
      USE prompt_mod 
      USE precision_mod 
      USE times_mod
      IMPLICIT none 
       
!                                                                       
      REAL ss, seknds
      REAL (KIND=PREC_DP) :: dnorm
      INTEGER lbeg (3), csize (3) 
      INTEGER iscat, nlot, ncell, i 
      INTEGER :: ii          ! Dummy variable
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
         csf (i) = cmplx (0.0D0, 0.0D0, KIND=KIND(0.0D0)) 
         acsf(i) = cmplx (0.0D0, 0.0D0, KIND=KIND(0.0D0)) 
         dsi (i) = 0.0d0 
      ENDDO 
!                                                                       
!------ preset some tables, calculate average structure                 
!                                                                       
      CALL fourier_lmn(eck,vi,inc,lmn,off_shift)
      CALL four_cexpt 
      CALL four_stltab 
      IF (ier_num.ne.0) return 
      CALL four_formtab 
      CALL four_csize (cr_icc, csize, lperiod, ls_xyz) 
      CALL four_aver (ilots, fave, csize) 
!                                                                       
!------ loop over crystal regions (lots)                                
!                                                                       
      loop_lots: DO nlot = 1, nlots 
         DO i = 1, num (1) * num (2) * num(3)
            csf (i) = cmplx (0.0D0, 0.0D0, KIND=KIND(0.0D0)) 
         ENDDO 
!                                                                       
         IF(lot_all) THEN
!           ii      = nlot - 1 
!           lbeg(3) = (ii  )/(cr_icc(1)*cr_icc(2)) + 1
!           ii      = (ii  ) - (lbeg(3)-1)*(cr_icc(1)*cr_icc(2))
!           lbeg(2) = (ii  )/(cr_icc(1)) + 1
!           ii      = (ii  ) - (lbeg(2)-1)*(cr_icc(1))
!           lbeg(1) = MOD(ii,cr_icc(1)) + 1
            ii      = nlot - 1 
            lbeg(3) = (ii  )/(csize(1)*csize(2)) + 1
            ii      = (ii  ) - (lbeg(3)-1)*(csize(1)*csize(2))
            lbeg(2) = (ii  )/(csize(1)) + 1
            ii      = (ii  ) - (lbeg(2)-1)*(csize(1))
            lbeg(1) = MOD(ii,csize(1)) + 1
         ELSE
            CALL four_ranloc (csize, lbeg) 
         ENDIF
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
            IF(ier_num/=0.OR.ier_ctrlc) RETURN      ! An error occured or CTRL-C
         ENDDO loop_atoms
!                                                                       
!------ - subtract average structure factor, add intensity              
!                                                                       
         DO i = 1, num (1) * num (2) *num (3)
            csf (i) = csf (i) - acsf (i) 
            dsi (i) = dsi (i) + DBLE (csf (i) * conjg (csf (i) ) ) 
         ENDDO 
      ENDDO loop_lots
!DO i = 1, num (1) * num (2) *num (3)
!dsi (i) = DBLE (acsf (i) * conjg (acsf (i) ) ) 
!ENDDO 
!                                                                       
!------ if we had lots, normalise the intensities                       
!                                                                       
      IF (nlots.ne.1) then 
         dnorm = DBLE (cr_icc (1) * cr_icc (2) * cr_icc (3) ) / &
                 DBLE (ncell * nlots)                                                 
!         dnorm = dnorm**2    Non-squared gives constant intensity for all sizes and numbers
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
      SUBROUTINE four_aver (lots, ave, csize) 
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
      USE precision_mod
      IMPLICIT none 
       
!
      INTEGER, INTENT(IN) :: lots
      REAL   , INTENT(IN) :: ave 
      INTEGER, DIMENSION(3), INTENT(IN) :: csize
!                                                                       
      REAL ran1
      REAL (KIND=PREC_DP) :: norm
      INTEGER isite, iatom, iscat, icell (3) 
      INTEGER scell, ncell, j, ii, jj, kk 
      LOGICAL, DIMENSION(:,:,:), ALLOCATABLE :: sel_cell
!                                                                       
      ave_is_zero: IF (ave.eq.0.0) then 
         RETURN 
      ELSE ave_is_zero  
         IF (four_log) then 
            WRITE (output_io, 1000) 100.0 * ave 
         ENDIF 
         scell = csize (1) * csize (2) * csize (3) 
         ncell = 0 
         ALLOCATE(sel_cell(csize(1), csize(2), csize(3)))
         sel_cell(:,:,:) = .FALSE.
         DO ii = 1, csize (1) 
            DO jj = 1, csize (2) 
               DO kk = 1, csize (3) 
                  sel_cell(ii,jj,kk) = (ran1 (idum) .le.ave) 
               ENDDO
            ENDDO
         ENDDO
!                                                                       
!------ ----- Loop over all atom types                                  
!                                                                       
            loop_iscat: DO iscat = 1, cr_nscat 
               nxat = 0 
!                                                                       
!------ - Loop over all unit cells                                      
!                                                                       
         cell_x: DO ii = 1, csize (1) 
         cell_y: DO jj = 1, csize (2) 
         cell_z: DO kk = 1, csize (3) 
         icell (1) = ii 
         icell (2) = jj 
         icell (3) = kk 
!                                                                       
!------ --- get only 'ave'% of those unit cells and compute average     
!------ --- unit cell                                                   
!                                                                       
!        sel = (ran1 (idum) .le.ave) 
!        IF (sel) then 
         IF (sel_cell(ii,jj,kk)) then 
            ncell = ncell + 1 
!RBN!                                                                       
!RBN!------ ----- Loop over all atom types                                  
!RBN!                                                                       
!RBN            loop_iscat: DO iscat = 1, cr_nscat 
!RBN               nxat = 0 
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
!RBN               CALL four_strucf_aver (iscat, .true.) 
!RBN               DO j = 1, num (1) * num (2) *num (3)
!RBN                  acsf (j) = acsf (j) + tcsf (j) 
!RBN               ENDDO 
!RBN            ENDDO  loop_iscat
         ENDIF 
         ENDDO cell_z
         ENDDO cell_y
         ENDDO cell_x
               CALL four_strucf (iscat, .true.) 
               DO j = 1, num (1) * num (2) *num (3)
                  acsf (j) = acsf (j) + tcsf (j) 
               ENDDO 
            ENDDO  loop_iscat
         ncell = ncell/cr_nscat
         DEALLOCATE(sel_cell)
!                                                                       
!------ - now compute the interference function of the lot shape        
!                                                                       
         CALL four_getav (lots) 
         CALL four_strucf_aver (0, .false.) 
         IF(ncell >0) THEN
            norm = DBLE(1.0D0 / ncell)
         ELSE
            norm = 0.0D0
            ier_num = +1
            ier_typ = ER_FOUR
            ier_msg(1) = 'Does the crystal consist of just 1 unit cell?'
            ier_msg(2) = 'Increase the percentage for >set aver<'
         ENDIF
         DO j = 1, num (1) * num (2) * num (3)
            acsf (j) = acsf (j) * tcsf (j) * cmplx ( norm, 0.0D0, KIND=KIND(0.0D0))
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
 2000 FORMAT     (' Used ',F5.1,' % of the crystal to calculate <F> ...') 
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
         csize (1) = MAX(1, cr_icc (1) - ls_xyz (1) + 1 )
         csize (2) = MAX(1, cr_icc (2) - ls_xyz (2) + 1 )
         csize (3) = MAX(1, cr_icc (3) - ls_xyz (3) + 1 )
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
      USE precision_mod 
      USE wink_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      REAL(PREC_DP) xmult, xarg 
      INTEGER i 
!                                                                       
      IF (.not.ffour) then 
         WRITE (output_io, 1000) 
!                                                                       
!        zpi = 8.0d0 * datan (1.0d0) 
!                                                                       
         DO i = 0, MASK 
            xmult   = (dble (i) * 1.0d0) / dble (I2PI) 
            xarg    = zpi * xmult 
            cex (i) = CMPLX (DBLE( COS (xarg)), DBLE( SIN (xarg)), KIND=KIND(0.0D0) ) 
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
      h (k) = REAL(xm (k) + uin (k) * REAL (i - 1, KIND=KIND(0.0D0)) &
                          + vin (k) * REAL (j - 1, KIND=KIND(0.0D0)) &
                          + win (k) * REAL (l - 1, KIND=KIND(0.0D0)))
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
      USE precision_mod
      IMPLICIT none 
       
!                                                                       
      REAL                :: q2
      REAL (KIND=PREC_DP) :: sb, sf, sfp, sfpp 
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
      q2 = REAL((REAL (iq, KIND=KIND(0.0D0)) * CFINC) **2 , KIND=KIND(0.0E0))
      sf = DBLE(form (iscat, cr_scat, lxray, q2, diff_power) )
!                                                                       
      IF (ano) then 
         sfp  = DBLE(cr_delfr ( (iscat) ) )
         sfpp = DBLE(cr_delfi ( (iscat) ) )
      ELSE 
         sfp  = 0.0D0 
         sfpp = 0.0D0
      ENDIF 
!                                                                       
      IF (ldbw) then 
         sb = exp ( - DBLE(cr_dw ( (iscat) ) * q2)) * DBLE(cr_occ(iscat))
      ELSE 
         sb = 1.0D0 * DBLE(cr_occ(iscat))
      ENDIF 
!                                                                       
      cfact     (iq, iscat) = cmplx (sb * (sf + sfp), sb * sfpp, KIND=KIND(0.0D0)) 
      cfact_pure(iq, iscat) = cmplx (     (sf + sfp),      sfpp, KIND=KIND(0.0D0)) 
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
         braggmax = max (braggmax, REAL(dsi (k)) ) 
         braggmin = min (braggmin, REAL(dsi (k)) ) 
      ELSE 
         diffumax = max (diffumax, REAL(dsi (k)) ) 
         diffumin = min (diffumin, REAL(dsi (k)) ) 
!                                                                       
         dsum  = dsum  + REAL( dsi (k) )
         dsum2 = dsum2 + REAL( dsi (k) **2) 
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
      USE precision_mod
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
      SUBROUTINE dlink (ano, lambda, rlambda, renergy, l_energy, &
                        diff_radiation, diff_power) 
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
      REAL                , INTENT(INOUT):: renergy
      LOGICAL             , INTENT(IN)   :: l_energy
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
      CALL get_wave ( lambda, rlambda, renergy, l_energy,diff_radiation, &
                      ier_num, ier_typ )
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
      USE precision_mod
!
      IMPLICIT none 
!                                                                       
       
!                                                                       
      REAL, DIMENSION(3), INTENT(IN) :: rhkl
!
      INTEGER i, j 
      INTEGER shel_inc (3) 
      REAL shel_eck (3, 4) 
      REAL shel_vi (3, 3) 
      COMPLEX (KIND=PREC_DP) :: shel_acsf 
      REAL    (KIND=PREC_DP) :: shel_dsi 
      COMPLEX (KIND=PREC_DP) :: shel_tcsf 
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
      res_para (1) =      REAL (csf (1) , KIND=KIND(0.0E0) ) 
      res_para (2) = REAL(AIMAG(csf (1)), KIND=KIND(0.0E0) ) 
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
!
      SUBROUTINE calc_hkl(infile,infile_l, calcfile, calcfile_l,scale,style)
!
      USE crystal_mod 
      USE diffuse_mod 
      USE discus_allocate_appl_mod
      USE get_params_mod
      USE param_mod
!
      IMPLICIT NONE
!
      CHARACTER(LEN=*), INTENT(IN) :: infile
      INTEGER         , INTENT(IN) :: infile_l
      CHARACTER(LEN=*), INTENT(IN) :: calcfile
      INTEGER         , INTENT(IN) :: calcfile_l
      REAL            , INTENT(IN) :: scale
      INTEGER         , INTENT(IN) :: style
!
      INTEGER, PARAMETER :: ird = 54
      INTEGER, PARAMETER :: iwr = 55
      INTEGER, PARAMETER :: HKLF4 = 4
      INTEGER, PARAMETER :: CIF   = 1
!
      CHARACTER(LEN=1024) :: line
      CHARACTER(LEN=1024), DIMENSION(:), ALLOCATABLE   :: ccpara
      INTEGER            , DIMENSION(:), ALLOCATABLE   :: llpara
      INTEGER            :: n_qxy, n_natoms,n_nscat
      INTEGER            :: iostatus
      INTEGER            :: ih,ik,il
      INTEGER            :: ih_min,ik_min,il_min
      INTEGER            :: ih_max,ik_max,il_max
      INTEGER            :: n_refl
      INTEGER            :: indx
      INTEGER            :: startline
      INTEGER            :: j
      INTEGER            :: nentries
      INTEGER            :: length, ianz
      INTEGER            ::   j_h      = 0
      INTEGER            ::   j_k      = 0
      INTEGER            ::   j_l      = 0
      INTEGER            ::   j_iobs   = 0
      INTEGER            ::   j_icalc  = 0
      INTEGER            ::   j_sigi   = 0
      INTEGER            ::   j_fobs   = 0
      INTEGER            ::   j_fcalc  = 0
      INTEGER            ::   j_sigf   = 0
      INTEGER            ::   j_flag   = 0
      REAL               :: rint, sint, wert
      REAL, DIMENSION(7) :: values
      REAL, DIMENSION(3) :: rhkl
!
      n_qxy    = 1
      n_natoms = 1
      n_nscat  = 1
!
      CALL oeffne(ird, infile(1:infile_l),   'old') 
      CALL oeffne(iwr, calcfile(1:calcfile_l), 'unknown') 
      inc(:) = 1
!
      ih_min  = 0
      ih_max  = 0
      ik_min  = 0
      ik_max  = 0
      il_min  = 0
      il_max  = 0
      n_refl  = 0
      startline = 0
      IF(style==HKLF4) THEN
         startline = 1
check:   DO     ! First read, get size and extrema
            READ(ird,'(a)', IOSTAT=iostatus) line
            IF(IS_IOSTAT_END(iostatus)) EXIT check
            IF(line(1:1)=='#' .OR. line(1:1)=='!' .OR. line(1:1)=='d') THEN
               ier_num = -9999
               ier_typ = ER_APPL
               CLOSE(ird) 
               RETURN
            ENDIF
            READ(line,1000, IOSTAT=iostatus) ih,ik,il, rint, sint
            ih_min = MIN(ih_min,ih)
            ih_max = MAX(ih_max,ih)
            ik_min = MIN(ik_min,ik)
            ik_max = MAX(ik_max,ik)
            il_min = MIN(il_min,il)
            il_max = MAX(il_max,il)
            n_refl = n_refl + 1
         ENDDO check
      ELSEIF(style==CIF) THEN
find:    DO
            READ(ird,'(a)', IOSTAT=iostatus) line
            IF(IS_IOSTAT_END(iostatus)) THEN
               ier_num = -9999
               ier_typ = ER_APPL
               CLOSE(ird) 
               RETURN
            ENDIF
            startline = startline + 1
            IF(line(1:14) == '_refln_index_h') EXIT find
         ENDDO find
         nentries = 0
         j_h      = 0
         j_k      = 0
         j_l      = 0
         j_iobs   = 0
         j_icalc  = 0
         j_sigi   = 0
         j_fobs   = 0
         j_fcalc  = 0
         j_sigf   = 0
         j_flag   = 0
         j = 1
entries: DO
            IF(line(1:1)=='#' .or. line == ' ') THEN
            startline = startline + 1
               CYCLE entries
            ENDIF
            IF(line(1:14) == '_refln_index_h')         j_h = j
            IF(line(1:14) == '_refln_index_k')         j_k = j
            IF(line(1:14) == '_refln_index_l')         j_l = j
            IF(line(1:21) == '_refln_F_squared_calc')  j_icalc = j
            IF(line(1:21) == '_refln_F_squared_meas')  j_iobs  = j
            IF(line(1:22) == '_refln_F_squared_sigma') j_sigi  = j
            IF(line(1:22) == '_refln_observed_status') j_flag  = j
            IF(line(1:7) /= '_refln_') EXIT entries
            nentries = nentries + 1
            READ(ird,'(a)', IOSTAT=iostatus) line
            IF(IS_IOSTAT_END(iostatus)) THEN
               ier_num = -9999
               ier_typ = ER_APPL
               CLOSE(ird) 
               RETURN
            ENDIF
            j = j + 1
            startline = startline + 1
         ENDDO entries
         ALLOCATE(ccpara(1:nentries))
         ccpara(:) = ' '
         ALLOCATE(llpara(1:nentries))
         length = LEN_TRIM(line)
check2:  DO
            CALL get_params_blank(line,ianz, ccpara,llpara, nentries, length)
            READ(ccpara(j_h)(1:llpara(j_h)),*) ih 
            READ(ccpara(j_k)(1:llpara(j_k)),*) ik 
            READ(ccpara(j_l)(1:llpara(j_l)),*) il 
            ih_min = MIN(ih_min,ih)
            ih_max = MAX(ih_max,ih)
            ik_min = MIN(ik_min,ik)
            ik_max = MAX(ik_max,ik)
            il_min = MIN(il_min,il)
            il_max = MAX(il_max,il)
            n_refl = n_refl + 1
            READ(ird,'(a)', IOSTAT=iostatus) line
!        READ(ird,*   , IOSTAT=iostatus) ih,ik,il, rint, sint
            IF(IS_IOSTAT_END(iostatus)) EXIT check2
            length = LEN_TRIM(line)
         ENDDO check2
      ELSE
         write(*,*) ' WRONG STYLE ', style
         return
      ENDIF
      vi(:,:) = 0
      vi(1,1) = 1.0
      vi(2,2) = 1.0
      vi(3,3) = 1.0
      inc(1)   = ih_max - ih_min             + 1
      inc(2)   = ik_max - ik_min             + 1
      inc(3)   = il_max - il_min             + 1
      eck(1,1) = ih_min             ! minimum H
      eck(2,1) = ik_min             ! minimum K
      eck(3,1) = il_min             ! minimum L
      eck(1,2) = ih_max             ! maximum H
      eck(2,2) = ik_min             ! minimum K
      eck(3,2) = il_min             ! minimum L
      eck(1,3) = ih_min             ! minimum H
      eck(2,3) = ik_max             ! maximum K
      eck(3,3) = il_min             ! minimum L
      eck(1,4) = ih_min             ! minimum H
      eck(2,4) = ik_min             ! minimum K
      eck(3,4) = il_max             ! maximum L
!
      IF (ier_num == 0) then 
         IF (inc(1) * inc(2) * inc(3) .gt. MAXQXY  .OR.   &
             cr_natoms > DIF_MAXAT                 .OR.   &
             cr_nscat>DIF_MAXSCAT              ) THEN
            n_qxy    = MAX(n_qxy,inc(1) * inc(2)*inc(3),MAXQXY)
            n_natoms = MAX(n_natoms,cr_natoms,DIF_MAXAT)
            n_nscat  = MAX(n_nscat,cr_nscat,DIF_MAXSCAT)
            CALL alloc_diffuse (n_qxy, n_nscat, n_natoms)
            IF (ier_num /= 0) THEN
               RETURN
            ENDIF
         ENDIF
         CALL dlink (ano, lambda, rlambda, renergy, l_energy, &
                     diff_radiation, diff_power) 
         call four_run
         rhkl (1) =  0. 
         rhkl (2) =  0. 
         rhkl (3) =  0. 
         REWIND(ird)
         DO j=1,startline -1 
            READ(ird,'(a)') line
         ENDDO
main:    DO
            IF(style==HKLF4) THEN
               READ(ird,1000, IOSTAT=iostatus) ih,ik,il, rint, sint
            ELSEIF(style==CIF) THEN
               READ(ird,'(a)',IOSTAT=iostatus) line
               length = LEN_TRIM(line)
               CALL get_params_blank(line,ianz, ccpara,llpara, nentries, length)
               DO j=1,nentries-1
                  READ(ccpara(j)(1:llpara(j)),*) values(j)
               ENDDO
               ih = NINT(values(1))
               ik = NINT(values(2))
               il = NINT(values(3))
!              READ(ird,*   , IOSTAT=iostatus) ih,ik,il, rint, sint
            ENDIF
            IF(IS_IOSTAT_END(iostatus)) EXIT main
            indx = (ih-ih_min)*inc(3)*inc(2) + (ik-ik_min)*inc(3) + (il-il_min)  + 1
            wert = REAL(csf(indx)*CONJG(csf(indx)),KIND=KIND(0.0E0))
            sint = SQRT(ABS(wert))
            IF(ih==0 .AND. ik==0 .AND. IL==0) THEN
               WRITE(iwr,1000) ih,ik,il, 0.00, 0.00
            ELSE
               IF(style==HKLF4) THEN
                  WRITE(iwr,1000) ih,ik,il, scale*wert, sqrt(scale)*sint 
               ELSEIF(style==CIF) THEN
                  IF(j_icalc/=0) THEN
                     values(j_icalc) = scale*wert
                  ENDIF
                  WRITE(iwr, 2000) ih,ik,il, (values(j),j=4, nentries-1)
               ENDIF
            ENDIF
         END DO main
      ENDIF 
      IF(ALLOCATED(ccpara))   DEALLOCATE(ccpara)
      IF(ALLOCATED(llpara))   DEALLOCATE(llpara)
1000  FORMAT(3I4,F8.2,F8.2)
2000  FORMAT(3I4,F12.2,F12.2, F12.2)
      CLOSE(ird)
      CLOSE(iwr)
      END SUBROUTINE calc_hkl
!
!*****7*****************************************************************
      SUBROUTINE four_strucf_aver (iscat, lform) 
!+                                                                      
!     Here the complex structure factor of 'nxat' identical atoms       
!     from array 'xat' is computed.                                     
!
!     The phase "iarg0" is calculated via integer math as offset from 
!     phase = 0 at hkl=0.
!-                                                                      
      USE discus_config_mod 
      USE diffuse_mod 
      USE precision_mod
      IMPLICIT none 
!                                                                       
      INTEGER, INTENT(IN) :: iscat 
      LOGICAL, INTENT(IN) :: lform 
!                                                                       
      REAL(PREC_DP)        ::        xincu, xincv , xincw
      REAL(PREC_DP)        ::        oincu, oincv , oincw
      INTEGER (KIND=PREC_INT_LARGE)   :: h, i, ii, j, k, iarg, iarg0, iincu, iincv, iincw
      INTEGER (KIND=PREC_INT_LARGE)   ::                              jincu, jincv, jincw
      INTEGER (KIND=PREC_INT_LARGE), PARAMETER :: shift = -6
!
      INTEGER IAND, ISHFT 
!
!------ zero fourier array                                              
!                                                                       
      tcsf = cmplx (0.0D0, 0.0D0, KIND=KIND(0.0D0)) 
!                                                                       
!------ Loop over all atoms in 'xat'                                    
!                                                                       
      DO k = 1, nxat 
!        xarg0 = xm (1)        * xat(k, 1) + xm (2)         * xat(k, 2) + xm (3)         * xat(k, 3)
         xincu = uin(1)        * xat(k, 1) + uin(2)         * xat(k, 2) + uin(3)         * xat(k, 3)
         xincv = vin(1)        * xat(k, 1) + vin(2)         * xat(k, 2) + vin(3)         * xat(k, 3)
         xincw = win(1)        * xat(k, 1) + win(2)         * xat(k, 2) + win(3)         * xat(k, 3)
         oincu = off_shift(1,1)* xat(k, 1) + off_shift(2,1) * xat(k, 2) + off_shift(3,1) * xat(k, 3)
         oincv = off_shift(1,2)* xat(k, 1) + off_shift(2,2) * xat(k, 2) + off_shift(3,2) * xat(k, 3)
         oincw = off_shift(1,3)* xat(k, 1) + off_shift(2,3) * xat(k, 2) + off_shift(3,3) * xat(k, 3)
!                                                                       
!        iarg0 = nint (64 * I2PI * (xarg0 - int (xarg0) + 0.0d0) ) 
         iincu = nint (64 * I2PI * (xincu - int (xincu) + 0.0d0) ) 
         iincv = nint (64 * I2PI * (xincv - int (xincv) + 0.0d0) ) 
         iincw = nint (64 * I2PI * (xincw - int (xincw) + 0.0d0) ) 
         jincu = nint (64 * I2PI * (oincu - int (oincu) + 0.0d0) ) 
         jincv = nint (64 * I2PI * (oincv - int (oincv) + 0.0d0) ) 
         jincw = nint (64 * I2PI * (oincw - int (oincw) + 0.0d0) ) 
         iarg0 =  lmn(1)*iincu + lmn(2)*iincv + lmn(3)*iincw + &
                  lmn(4)*jincu + lmn(5)*jincv + lmn(6)*jincw
         iarg = iarg0 
!                                                                       
!------ - Loop over all points in Q. 'iadd' is the address of the       
!------ - complex exponent table. 'IADD' divides out the 64 and         
!------ - ISHFT acts as MOD so that the argument stays in the table     
!------ - boundaries.                                                   
!                 iadd      = ISHFT (iarg, - 6) 
!                 iadd      = IAND  (iadd, MASK) 
!                 tcsf (ii) = tcsf (ii) + cex (iadd, MASK) )
!                 iarg      = iarg + iincw
!                                                                       
         ii = 0 
!                                                                       
          DO j = 0, num (1) - 1
             DO i = 0, num (2) - 1
                iarg = iarg0 + iincu*j + iincv*i 
                DO h = 1, num (3) 
                   ii       = ii + 1 
                   tcsf(ii) = tcsf (ii) + cex (IAND  (ISHFT(iarg, shift), MASK) )
                   iarg     = iarg + iincw
                ENDDO 
             ENDDO 
          ENDDO 
      ENDDO 
!
!------ Now we multiply with formfactor                                 
!                                                                       
      IF (lform) then 
         DO  i = 1, num (1) * num (2) * num(3)
!        FORALL( i = 1: num (1) * num (2) * num(3))   !!! DO Loops seem to be faster!
            tcsf (i) = tcsf (i) * cfact (istl (i), iscat) 
!        END FORALL
         END DO
      ENDIF 
!                                                                       
      END SUBROUTINE four_strucf_aver
!
END MODULE fourier_sup
