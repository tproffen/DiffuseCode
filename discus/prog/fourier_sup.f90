MODULE fourier_sup
!
USE errlist_mod 
USE precision_mod
!
COMPLEX(KIND=PREC_DP), DIMENSION(:,:), ALLOCATABLE :: fft_field
COMPLEX(KIND=PREC_DP), DIMENSION(:,:), ALLOCATABLE :: fft_sum  
INTEGER, DIMENSION(2) :: fft_dim      ! Size of fft_field
INTEGER               :: fft_grid     ! No of grid points in a unit cell
!
CONTAINS
!*****7*****************************************************************
!
SUBROUTINE four_run 
!+                                                                      
!     Main routine to calculate the Fourier transform for the           
!     given plane in reciprocal space. Based on the program             
!     'diffuse' written by B.D. Butler.                                 
!-                                                                      
USE discus_config_mod 
USE crystal_mod 
USE diffuse_mod 
USE fourier_conv_mod
USE four_strucf_mod
USE fourier_lmn_mod
!                                                                       
USE prompt_mod 
USE precision_mod 
USE times_mod
USE support_mod
!
IMPLICIT none 
!                                                                       
REAL(KIND=PREC_SP) :: ss!, seknds
REAL(KIND=PREC_DP) :: dnorm
INTEGER :: lbeg (3), csize (3) 
INTEGER :: iscat, nlot, ncell, i 
INTEGER :: ii          ! Dummy variable
!                                                                       
ier_num = 0 
!                                                                       
!------ preset some values                                              
!                                                                       
ss = seknds (0.0) 
!
!                                                                       
!------ preset some values                                              
!                                                                       
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
!write(*,*) ' eck ll ', eck(:,1)
!write(*,*) ' eck lr ', eck(:,2)
!write(*,*) ' eck ul ', eck(:,3)
!write(*,*) ' eck tl ', eck(:,4)
!write(*,*) ' vi abs ', vi(:,1), uin
!write(*,*) ' vi ord ', vi(:,2), vin
!write(*,*) ' vi top ', vi(:,3), win
!write(*,*) ' OFF  x ', off_shift(:,1)
!write(*,*) ' OFF  y ', off_shift(:,2)
!write(*,*) ' OFF  z ', off_shift(:,3)
!write(*,*) ' inc num', inc, num
!write(*,*) ' lmn    ', lmn
IF (ier_num.ne.0) return 
!
CALL four_formtab
!
!----- Test for Friedels law, and reduce calculation time
!
CALL four_test_friedel
CALL four_csize (cr_icc, csize, lperiod, ls_xyz) 
CALL four_aver (ilots, fave, csize) 
!CALL four_fft_prep
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
!call four_strucf_fft(iscat, .TRUE., nlot)
!                                                                       
!------ --- Add this part of the structur factor to the total           
!                                                                       
            DO i = 1, num (1) * num (2) *num (3)
               csf (i) = csf (i) + tcsf (i) 
            ENDDO 
            IF (four_log) then 
               WRITE (output_io, 3000) cr_at_lis (iscat), nxat 
            ENDIF 
            IF(ier_ctrlc) THEN
               ier_num = -14
               ier_typ = ER_COMM
               RETURN
            ENDIF
            IF(ier_num/=0) RETURN      ! An error occured or CTRL-C
         ENDDO loop_atoms
!                                                                       
!------ - subtract average structure factor, add intensity              
!                                                                       
         DO i = 1, num (1) * num (2) *num (3)
            csf (i) = csf (i) - acsf (i) !* fave_sca
            dsi (i) = dsi (i) + DBLE (csf (i) * conjg (csf (i) ) ) 
         ENDDO 
         CALL four_apply_friedel
      ENDDO loop_lots
!DO i = 1, num (1) * num (2) *num (3)
!dsi (i) = DBLE (acsf (i) * conjg (acsf (i) ) ) 
!ENDDO 
!
!----- Reset for Friedels law, and reduce calculation time
!
CALL four_rese_friedel
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
!if(allocated(dsi3d)) deallocate(dsi3d)
!allocate(dsi3d(1:ubound(dsi,1)))
!dsi3d = dsi
CALL four_weight               ! Correct the relative weight of Bragg and diffuse
!
      CALL four_conv           ! Convolute diffraction pattern
!                                                                       
      CALL four_qinfo 
      ss = seknds (ss) 
      IF (four_log) then 
         WRITE (output_io, 4000) ss 
      ENDIF 
!
CALL four_accumulate
!
!write(*,*) ' eck ll ', eck(:,1)
!write(*,*) ' eck lr ', eck(:,2)
!write(*,*) ' eck ul ', eck(:,3)
!write(*,*) ' eck tl ', eck(:,4)
!write(*,*) ' vi abs ', vi(:,1), uin
!write(*,*) ' vi ord ', vi(:,2), vin
!write(*,*) ' vi top ', vi(:,3), win
!write(*,*) ' OFF  x ', off_shift(:,1)
!write(*,*) ' OFF  y ', off_shift(:,2)
!write(*,*) ' OFF  z ', off_shift(:,3)
!write(*,*) ' inc num', inc, num
!write(*,*) ' lmn    ', lmn
!
!CALL four_fft_finalize
!                                                                       
 2000 FORMAT     (/,' Lot # ',I4,'/',I4,' : starting at unit cell (',   &
     &                       3(I3,1X),')')                              
 2010 FORMAT     (/,' Lot # ',I4,'/',I4,' : complete crystal') 
 3000 FORMAT (  '                   Atom typ = ',A4,7X,'(# ',I9,' )') 
 4000 FORMAT     (/,' Elapsed time    : ',G13.6,' sec') 
      END SUBROUTINE four_run                       
!*****7*****************************************************************
!
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
USE lib_random_func
USE prompt_mod 
USE precision_mod
!
IMPLICIT none 
       
!
INTEGER, INTENT(IN) :: lots
REAL(kind=PREC_DP)   , INTENT(IN) :: ave 
INTEGER, DIMENSION(3), INTENT(IN) :: csize
!                                                                       
REAL (KIND=PREC_DP) :: norm
INTEGER :: isite, iatom, iscat, icell (3) 
INTEGER :: scell, ncell, j, ii, jj, kk 
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
                     xat (nxat, j) = cr_pos (j, iatom) - REAL(icell (j)     &
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
!        CALL four_strucf_aver (0, .false.) 
         CALL four_strucf      (0, .false.) 
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
            WRITE (output_io, 2000) (REAL(ncell) / REAL(scell) )    &
            * 100.0                                                     
         ENDIF 
!                                                                       
      ENDIF ave_is_zero
!                                                                       
 1000 FORMAT     (' Calculating <F> with ',F5.1,' % of the crystal ...') 
 2000 FORMAT     (' Used ',F5.1,' % of the crystal to calculate <F> ...') 
      END SUBROUTINE four_aver                      
!
!*****7*****************************************************************
!
subroutine four_weight
!-
!  Correct the weight of Bragg versus diffuse if fave /= 100 and
!  the grid size or lot size does not correspond to 1/cr_icc.
!
!  This is needed only if the full 3D-PDF is calculated.  
!
!  If the volume of a voxel in reciprocal space does not 
!  correspond to the volume of the reciprocal unit cell divided by
!  the number of unit cell dimensions**3 we need to scale
!
!  2D 
!  If the area of the 2D-Pixel does not correspond to the area of
!  the area of the reciprocal unit cell section divided by the 
!  number of unit cell dimensions**2 we need to scale
!
!  1D
!  If the length os a step in reciprocal space does not
!  correspond to the length of the reciprocal unit cell divided by
!  the number of unit cells we need to scale
!
!+
!
use crystal_mod
use diffuse_mod
use metric_mod
!
use lib_errlist_func
use precision_mod
!
implicit none
!
integer               :: isdim      ! Fourier was calculated in 1,2,3 dimensions
integer               :: ncells
integer               :: i, j, k, ii
real(kind=PREC_DP)    :: scalef     ! scale factor to apply to the Bragg reflections
real(kind=PREC_DP), dimension(3  ) :: hkl     ! base vectors to calculate a volume
real(kind=PREC_DP), dimension(3,3) :: bases   ! base vectors to calculate a volume
real(kind=PREC_DP)                 :: voxel
!
scalef = 1.0D0
if_fave:if(fave == 0.000) then
!
   isdim = 3
   IF(num(1)==1) isdim = isdim - 1
   IF(num(2)==1) isdim = isdim - 1
   IF(num(3)==1) isdim = isdim - 1
!write(*,*) ' IN four_weight ISDIM ', num, isdim, dsort
!write(*,*) ' eck ll ', eck(:,1)
!write(*,*) ' eck lr ', eck(:,2)
!write(*,*) ' eck ul ', eck(:,3)
!write(*,*) ' eck tl ', eck(:,4)
!write(*,*) ' vi abs ', vi(:,1), uin
!write(*,*) ' vi ord ', vi(:,2), vin
!write(*,*) ' vi top ', vi(:,3), win
!write(*,*) ' OFF  x ', off_shift(:,1)
!write(*,*) ' OFF  y ', off_shift(:,2)
!write(*,*) ' OFF  z ', off_shift(:,3)
!write(*,*) ' inc num', inc, num
!
   ncells = 1
   if(isdim==1) then              ! 1D correction
      if(num(1)/=1) then
         bases(:,1) = uin
         bases(:,2) = 0.0D0
         bases(:,3) = 0.0D0
         if(bases(1,1)/=0D0) then
             bases(2,2)=1.0D0
             bases(3,3)=1.0D0
         elseif(bases(2,1)/=0D0) then
             bases(3,2)=1.0D0
             bases(1,3)=1.0D0
         elseif(bases(3,1)/=0D0) then
             bases(1,2)=1.0D0
             bases(2,3)=1.0D0
         endif
         if(nlots==1) then
            ncells = cr_icc(1)
         else
            ncells = ls_xyz(1)
         endif
      elseif(num(2)==1) then
         bases(:,1) = 0.0D0
         bases(:,2) = vin
         bases(:,3) = 0.0D0
         if(bases(1,2)/=0D0) then
             bases(2,1)=1.0D0
             bases(3,3)=1.0D0
         elseif(bases(2,2)/=0D0) then
             bases(3,1)=1.0D0
             bases(1,3)=1.0D0
         elseif(bases(3,2)/=0D0) then
             bases(1,1)=1.0D0
             bases(2,3)=1.0D0
         endif
         if(nlots==1) then
            ncells = cr_icc(2)
         else
            ncells = ls_xyz(2)
         endif
      elseif(num(3)==1) then
         bases(:,1) = 0.0D0
         bases(:,2) = 0.0D0
         bases(:,3) = win
         if(bases(1,3)/=0D0) then
             bases(2,1)=1.0D0
             bases(3,2)=1.0D0
         elseif(bases(2,3)/=0D0) then
             bases(3,1)=1.0D0
             bases(1,2)=1.0D0
         elseif(bases(3,3)/=0D0) then
             bases(1,1)=1.0D0
             bases(2,2)=1.0D0
         endif
         if(nlots==1) then
            ncells = cr_icc(3)
         else
            ncells = ls_xyz(3)
         endif
      endif
   elseif(isdim==2) then          ! 2D correction
      if(num(3)==1) then
         bases(:,1) = uin
         bases(:,2) = vin
         bases(:,3) = 1.0D0
         if(bases(1,1)/=0D0 .or. bases(1,2)/=0.0D0) bases(1,3)=0.0D0
         if(bases(2,1)/=0D0 .or. bases(2,2)/=0.0D0) bases(2,3)=0.0D0
         if(bases(3,1)/=0D0 .or. bases(3,2)/=0.0D0) bases(3,3)=0.0D0
         if(nlots==1) then
            ncells = cr_icc(1)*cr_icc(2)
         else
            ncells = ls_xyz(1)*ls_xyz(2)
         endif
      elseif(num(2)==1) then
         bases(:,1) = uin
         bases(:,2) = 1.0D0
         bases(:,3) = win
         if(bases(1,1)/=0D0 .or. bases(1,3)/=0.0D0) bases(1,2)=0.0D0
         if(bases(2,1)/=0D0 .or. bases(2,3)/=0.0D0) bases(2,2)=0.0D0
         if(bases(3,1)/=0D0 .or. bases(3,3)/=0.0D0) bases(3,2)=0.0D0
         if(nlots==1) then
            ncells = cr_icc(1)*cr_icc(3)
         else
            ncells = ls_xyz(1)*ls_xyz(3)
         endif
      elseif(num(3)==1) then
         bases(:,1) = 1.0D0
         bases(:,2) = vin
         bases(:,3) = win
         if(bases(1,2)/=0D0 .or. bases(1,3)/=0.0D0) bases(1,1)=0.0D0
         if(bases(2,2)/=0D0 .or. bases(2,3)/=0.0D0) bases(2,1)=0.0D0
         if(bases(3,2)/=0D0 .or. bases(3,3)/=0.0D0) bases(3,1)=0.0D0
         if(nlots==1) then
            ncells = cr_icc(2)*cr_icc(3)
         else
            ncells = ls_xyz(2)*ls_xyz(3)
         endif
      endif
!write(*,*) ' Base 1 ', bases(:,1)
!write(*,*) ' Base 2 ', bases(:,2)
!write(*,*) ' Base 3 ', bases(:,3)
!   voxel = do_volume(.FALSE., bases(:,1), bases(:,2), bases(:,3))
!   scalef = real(ncells)/cr_vr*voxel
!   write(*,*)  ' Voxel ', voxel
!   write(*,*)  ' Vr/Voxel ', cr_vr/voxel
!   write(*,*)  ' Ratio    ', real(ncells)/cr_vr*voxel, scalef
!
   elseif(isdim==3) then          ! 3D correction
      bases(:,1) = uin
      bases(:,2) = vin
      bases(:,3) = win
      if(nlots==1) then
         ncells = cr_icc(1)*cr_icc(2)*cr_icc(3)
      else
         ncells = ls_xyz(1)*ls_xyz(2)*ls_xyz(3)
   endif
   endif
!
   voxel = do_volume(.FALSE., bases(:,1), bases(:,2), bases(:,3))
   if(voxel>0.0) then
      scalef = real(ncells)/cr_vr*voxel
   else
      call no_error
   endif
endif if_fave
!
if(allocated(dsi3d)) deallocate(dsi3d)
allocate(dsi3d(1:ubound(dsi,1)))
if(fave==0.0 .and. abs(scalef-1.0)>1.0E-5) then
!write(*,*) 'base1  ', bases(:,1)
!write(*,*) 'base2  ', bases(:,2)
!write(*,*) 'base3  ', bases(:,3)
!write(*,*) 'Scalef ', scalef, ncells, cr_vr, voxel
!read(*,*) i
   dsi3d = dsi
   ii= 0
   do k=0,num(3)-1
      do j=0,num(2)-1
         do i=0,num(1)-1
            ii = ii + 1
            hkl(1) = eck(1,1) + vi(1,1)*i + vi(1,2)*j + vi(1,3)*k
            hkl(2) = eck(2,1) + vi(2,1)*i + vi(2,2)*j + vi(2,3)*k
            hkl(3) = eck(3,1) + vi(3,1)*i + vi(3,2)*j + vi(3,3)*k
            if(abs(hkl(1)-nint(hkl(1)))<1.0e-4 .and.                               &
               abs(hkl(2)-nint(hkl(2)))<1.0e-4 .and.                               &
               abs(hkl(3)-nint(hkl(3)))<1.0e-4       ) then
               dsi3d(ii) = dsi(ii)/scalef
            endif
         enddo
      enddo
   enddo
else
   dsi3d = dsi
endif
!
end subroutine four_weight
!
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
            offset (3) = offset (3) - REAL(cr_icc (3) ) 
         ENDIF 
         DO jj = 0, ls_xyz (2) - 1 
         icell (2) = jj + lbeg (2) 
         offset (2) = cr_dim0 (2, 1) + lbeg (2) - 1 
         IF (icell (2) .gt.cr_icc (2) ) then 
            icell (2) = icell (2) - cr_icc (2) 
            offset (2) = offset (2) - REAL(cr_icc (2) ) 
         ENDIF 
         DO ii = 0, ls_xyz (1) - 1 
         icell (1) = ii + lbeg (1) 
         offset (1) = cr_dim0 (1, 1) + lbeg (1) - 1 
         IF (icell (1) .gt.cr_icc (1) ) then 
            icell (1) = icell (1) - cr_icc (1) 
            offset (1) = offset (1) - REAL(cr_icc (1) ) 
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
         x0 (ii) = REAL(ls_xyz (ii) ) / 2.0 
         ENDDO 
!                                                                       
         DO kk = 0, ls_xyz (3) - 1 
         icell (3) = kk + lbeg (3) 
         offset (3) = cr_dim0 (3, 1) + lbeg (3) - 1 
         xtest (3) = (REAL(kk) - x0 (3) + 0.5) **2 / x0 (3) **2 
         IF (icell (3) .gt.cr_icc (3) ) then 
            icell (3) = icell (3) - cr_icc (3) 
            offset (3) = offset (3) - REAL(cr_icc (3) ) 
         ENDIF 
         DO jj = 0, ls_xyz (2) - 1 
         icell (2) = jj + lbeg (2) 
         offset (2) = cr_dim0 (2, 1) + lbeg (2) - 1 
         xtest (2) = (REAL(jj) - x0 (2) + 0.5) **2 / x0 (2) **2 
         IF (icell (2) .gt.cr_icc (2) ) then 
            icell (2) = icell (2) - cr_icc (2) 
            offset (2) = offset (2) - REAL(cr_icc (2) ) 
         ENDIF 
         DO ii = 0, ls_xyz (1) - 1 
         icell (1) = ii + lbeg (1) 
         offset (1) = cr_dim0 (1, 1) + lbeg (1) - 1 
         xtest (1) = (REAL(ii) - x0 (1) + 0.5) **2 / x0 (1) **2 
         IF (icell (1) .gt.cr_icc (1) ) then 
            icell (1) = icell (1) - cr_icc (1) 
            offset (1) = offset (1) - REAL(cr_icc (1) ) 
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
         xat (nxat, 1) = REAL(ii) + cr_dim0 (1, 1) 
         xat (nxat, 2) = REAL(jj) + cr_dim0 (2, 1) 
         xat (nxat, 3) = REAL(kk) + cr_dim0 (3, 1) 
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
         xat (nxat, 1) = REAL(ii) 
         xat (nxat, 2) = REAL(jj) 
         xat (nxat, 3) = REAL(kk) 
         ENDDO 
         ENDDO 
         ENDDO 
!                                                                       
!------ Ellipsoid shaped lot                                            
!                                                                       
      ELSEIF (lots.eq.LOT_ELI) then 
         DO ii = 1, 3 
         x0 (ii) = REAL(ls_xyz (ii) ) / 2.0 
         ENDDO 
!                                                                       
         DO kk = 0, ls_xyz (3) - 1 
         xtest (3) = (REAL(kk) - x0 (3) + 0.5) **2 / x0 (3) **2 
         DO jj = 0, ls_xyz (2) - 1 
         xtest (2) = (REAL(jj) - x0 (2) + 0.5) **2 / x0 (2) **2 
         DO ii = 0, ls_xyz (1) - 1 
         xtest (1) = (REAL(ii) - x0 (1) + 0.5) **2 / x0 (1) **2 
         IF ( (xtest (1) + xtest (2) + xtest (3) ) .le.1.0) then 
            nxat = nxat + 1 
            xat (nxat, 1) = REAL(ii) 
            xat (nxat, 2) = REAL(jj) 
            xat (nxat, 3) = REAL(kk) 
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
USE lib_random_func
      USE random_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER, DIMENSION(3), INTENT(IN)  :: csize
      INTEGER, DIMENSION(3), INTENT(OUT) :: lbeg
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
!
!*****7*****************************************************************
!
SUBROUTINE four_stltab 
!+                                                                      
!     Sets up an integer array containing the corresponding integers    
!     to the formfactor table for each sin(theta)/lambda.               
!-                                                                      
USE discus_config_mod 
USE crystal_mod 
USE diffuse_mod
!USE quad_mod 
use metric_mod ,only: skalpro
USE prompt_mod 
use precision_mod
!                                                                       
IMPLICIT none 
!                                                                       
REAL(kind=PREC_DP) :: q2, h (3) 
INTEGER :: i, j, k, l 
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
!     q2 = quad (h, h, cr_rten) / 4.0 
      q2 = skalpro (h, h, cr_rten) / 4.0 
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
!
IMPLICIT none 
!
!                                                                       
REAL                :: q2
REAL (KIND=PREC_DP) :: sb, sf, sfp, sfpp 
!REAL (KIND=PREC_DP) :: dw
INTEGER iq, iscat 
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
      h (k) = eck (k, 1) + vi (k, 1) * REAL(i - 1) + &
                           vi (k, 2) * REAL(j - 1) + &
                           vi (k, 3) * REAL(l - 1)
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
         diffuave = dsum / REAL(nd) 
         diffusig = sqrt (dsum2 / REAL(nd) - diffuave**2) 
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
 1010 FORMAT     (  ' Bragg scat.     : ',G13.6,'  -> ',G13.6) 
 1020 FORMAT     (  ' Diffuse scat.   : ',G13.6,'  -> ',G13.6) 
 1030 FORMAT     (  '      Average    : ',G13.6,'  +- ',G13.6) 
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
      REAL(kind=PREC_DP)   , DIMENSION(11,0:MAXSCAT), INTENT(INOUT) :: scat ! (11, 0:maxscat) 
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
!
!*****7*****************************************************************
!
SUBROUTINE dlink (ano, lambda, rlambda, renergy, l_energy, &
                  diff_radiation, diff_table, diff_power) 
!-                                                                      
!     This routine reads wavelength symbols, wavelength values 
!     and atomic form factors from module "element_data_mod"
!     It replaces the old dlink and LAZY_PULVERIX
!+                                                                      
USE discus_config_mod 
USE charact_mod
USE crystal_mod 
USE element_data_mod
USE chem_aver_mod
USE param_mod
use precision_mod
!
IMPLICIT none 
!                                                                       
!                                                                       
LOGICAL             , INTENT(IN)   :: ano             ! Anomalous scattering TRUE/FALSE
CHARACTER (LEN = * ), INTENT(IN)   :: lambda          ! Wavelength symbol
REAL(kind=PREC_DP)  , INTENT(OUT)  :: rlambda         ! Wave length value
REAL(kind=PREC_DP)  , INTENT(INOUT):: renergy         ! Wave length energy
LOGICAL             , INTENT(IN)   :: l_energy        ! Wave length specified as energy TRUE/FALSE
INTEGER             , INTENT(IN)   :: diff_radiation  ! xray/neutron/electron
INTEGER             , INTENT(IN)   :: diff_table      ! International / Waasmaier
INTEGER             , INTENT(OUT)  :: diff_power      ! 4 or 5 parameters
!
INTEGER, PARAMETER :: RAD_XRAY = 1
INTEGER, PARAMETER :: RAD_NEUT = 2
INTEGER, PARAMETER :: RAD_ELEC = 3
INTEGER, PARAMETER :: RAD_INTER = 0
LOGICAL, PARAMETER :: LOUT = .FALSE.
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
CALL chem_elem(lout)
!                                                                       
IF (ier_num.ne.0) RETURN 
!                                                                       
any_element: IF (cr_nscat.gt.0) THEN 
   DO i = 1, cr_nscat 
      IF (cr_at_lis (i) /= 'XAXI' .AND. cr_at_lis (i) /= 'YAXI' .and. & 
          cr_at_lis (i) /= 'ZAXI') THEN                  
         IF(res_para(i+1)> 0.0) THEN
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
                     diff_power   = PER_RAD_POWER(diff_radiation, 0)
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
                     if(diff_table==RAD_INTER) then
                        CALL get_scat_xray ( j, temp_scat )
                     else
                        CALL get_scat_xray_waas(j, temp_scat)
                     endif
                     cr_scat(:,i) = 0.0
                     cr_scat(:,i) = temp_scat(:)   ! copy temp array into 1st column
!                    cr_delfr (i) = 0.0   ! DEVELOPMENT !
!                    cr_delfi (i) = 0.0   ! DEVELOPMENT !
                     IF(ano) THEN
                        IF(cr_delf_int (i) ) then 
                           CALL get_scat_ano ( j, lambda, temp_delf )
                           cr_delfr (i) = temp_delf(1)
                           cr_delfi (i) = temp_delf(2)
                        ELSE
                           cr_delfr(i) =   cr_delfr_u(i)
                           cr_delfi(i) =   cr_delfi_u(i)
                        ENDIF 
                     ELSE
                        cr_delfr(i) = 0.0   ! DEVELOPMENT !
                        cr_delfi(i) = 0.0   ! DEVELOPMENT !
                     ENDIF 
                     diff_power   = PER_RAD_POWER(diff_radiation, diff_table)
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
                     diff_power   = PER_RAD_POWER(diff_radiation, 0)
                     ier_num = 0
                  ELSE
                     ier_typ = ER_NONE 
                     ier_typ = ER_APPL 
                  ENDIF
                  j = LEN_TRIM(element)
                  IF(MAXVAL(cr_scat(:,i))==0.0 .AND. &
                     (element(j:j)=='+' .OR. element(j:j)=='-')) THEN
                     ier_num = -173
                     ier_typ = ER_APPL
                     WRITE(ier_msg(1),'(3a,i4,a)') 'Element ',element(1:j), '(',i,')'
                     WRITE(ier_msg(2),'(a)') 'Replace by neutral atom or use in fourier'
                     WRITE(ier_msg(3),'(a)') 'scat <ion_name>, <neutron_name>'
                     RETURN
                  ENDIF
            END SELECT
!                                                                       
            ENDIF 
         ELSE
            cr_scat(:,i) = 0.00
            cr_delfr (i) = 0.0
            cr_delfi (i) = 0.0
         ENDIF
      ENDIF 
   ENDDO 
ELSE any_element
   ier_num = - 21 
   ier_typ = ER_APPL 
ENDIF any_element
!                                                                       
END SUBROUTINE dlink                          
!
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
USE prompt_mod
USE precision_mod
USE support_mod
!
      IMPLICIT NONE
!
      CHARACTER(LEN=*), INTENT(IN) :: infile
      INTEGER         , INTENT(IN) :: infile_l
      CHARACTER(LEN=*), INTENT(IN) :: calcfile
      INTEGER         , INTENT(IN) :: calcfile_l
      REAL(KIND=PREC_DP),INTENT(IN) :: scale
      INTEGER         , INTENT(IN) :: style
!
      INTEGER, PARAMETER :: ird = 54
      INTEGER, PARAMETER :: iwr = 55
      INTEGER, PARAMETER :: HKLF4 = 4
      INTEGER, PARAMETER :: CIF   = 1
!
      CHARACTER(LEN=PREC_STRING) :: line
      CHARACTER(LEN=PREC_STRING), DIMENSION(:), ALLOCATABLE   :: ccpara
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
         WRITE(output_io,*) ' WRONG STYLE ', style
         RETURN
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
                     diff_radiation, diff_table, diff_power) 
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
!
!*******************************************************************************
!
SUBROUTINE four_test_friedel
!
USE diffuse_mod
!
USE precision_mod
!
IMPLICIT NONE
!
REAL(PREC_SP), DIMENSION(3) :: rut   ! right upper top corner
REAL(PREC_SP), DIMENSION(3) :: dia   ! sum of all vi's
REAL(PREC_SP), DIMENSION(3) :: point ! sum of (left lower bottom) and (right upper top)
REAL(PREC_SP)               :: EPS = 1.E-6
!
diff_l_friedel = .FALSE.
IF(.NOT.ano) THEN        ! anomalous scattering is off, test space
   rut(1) = eck(1,1) + vi(1,1)*(inc(1)-1) + vi(1,2)*(inc(2)-1) + vi(1,3)*(inc(3)-1)
   rut(2) = eck(2,1) + vi(2,1)*(inc(1)-1) + vi(2,2)*(inc(2)-1) + vi(2,3)*(inc(3)-1)
   rut(3) = eck(3,1) + vi(3,1)*(inc(1)-1) + vi(3,2)*(inc(2)-1) + vi(3,3)*(inc(3)-1)
!
   diff_l_even(1) = MOD(inc(1),2)==0
   diff_l_even(2) = MOD(inc(2),2)==0
   diff_l_even(3) = MOD(inc(3),2)==0
   point(1) = eck(1,1) + rut(1)
   point(2) = eck(2,1) + rut(2)
   point(3) = eck(3,1) + rut(3)
   dia(1) = vi(1,1) + vi(1,2) + vi(1,3)
   dia(2) = vi(2,1) + vi(2,2) + vi(3,3)
   dia(3) = vi(3,1) + vi(3,2) + vi(3,3)
   IF(.NOT.diff_l_even(1) .AND. .NOT.diff_l_even(2) .AND. .NOT.diff_l_even(3)) THEN   ! All ODD increments
      diff_l_even(0) = .FALSE.    ! All odd
      IF(ABS(point(1))<EPS .AND. ABS(point(3))<EPS .AND.                           &
         ABS(point(2))<EPS                               ) THEN
!                         !left lower bottom and right upper top add up to zero
         diff_eck_u = eck         ! Back up user settings
         diff_vi_u  = vi          ! Back up user settings
         diff_inc_u = inc         ! Back up user settings
         IF(inc(1) > 1) THEN      ! We have an abscissa vector, cut this one
            inc(1) = (inc(1) + 1)/2
            diff_l_friedel = .TRUE.    ! A suitable solution
            diff_idim      = 1
         ELSEIF(inc(1)==1 .AND. inc(2)>1) THEN
            inc(2) = (inc(2) + 1)/2
            diff_l_friedel = .TRUE.    ! A suitable solution
            diff_idim      = 2
         ELSEIF(inc(1)==1 .AND. inc(2)==1 .AND. inc(3)>1) THEN
            inc(3) = (inc(3) + 1)/2
            diff_l_friedel = .TRUE.    ! A suitable solution
            diff_idim      = 3
         ELSE
            diff_l_friedel = .FALSE.   ! Not a suitable solution
            diff_idim      = 0
         ENDIF
      ENDIF
!  ELSEIF(diff_l_even(1) .AND. diff_l_even(2) .AND. diff_l_even(3)) THEN              ! All Even increments
!     diff_l_even(0) = .TRUE.    ! All even
!     IF(ABS(point(1)-dia(1))<EPS .AND. ABS(point(3)-dia(2))<EPS .AND.                           &
!        ABS(point(2)-dia(3))<EPS                                     ) THEN
!                         !left lower bottom and right upper top add up to sum of vi
!        diff_eck_u = eck         ! Back up user settings
!        diff_vi_u  = vi          ! Back up user settings
!        diff_inc_u = inc         ! Back up user settings
!        IF(inc(1) > 1) THEN      ! We have an abscissa vector, cut this one
!           inc(1) = inc(1)/2 + 1
!        ELSEIF(inc(1)==1 .AND. inc(2) > 1) THEN      ! We have an abscissa vector, cut this one
!           inc(2) = inc(2)/2 + 1
!        ENDIF
!     ENDIF
   ELSE
           diff_l_friedel = .FALSE.   ! Not a suitable solution
   ENDIF
!write(*,*) ' llb ', eck(:,1)
!write(*,*) ' rut ', rut(:)
!write(*,*) ' point',point(:)
!write(*,*) ' RES ', diff_l_friedel
!write(*,*) ' INC ', inc(:)
ENDIF 
END SUBROUTINE four_test_friedel
!
!*******************************************************************************
!
SUBROUTINE four_apply_friedel
!-
!  Applies Friedels lay to the Fourier calculation
!+
USE diffuse_mod
!
IMPLICIT NONE
!
INTEGER :: i, j
!
IF(diff_l_friedel) THEN
   j = diff_inc_u(1)*diff_inc_u(2)*diff_inc_u(3) + 1
   DO i = 1, num(1) * num(2) * num(3)
     acsf(j-i) =acsf(i)
      csf(j-i) = csf(i)
      dsi(j-i) = dsi(i)
   ENDDO 
ENDIF
!
END SUBROUTINE four_apply_friedel
!
!*******************************************************************************
!
SUBROUTINE four_rese_friedel
!
! Set the corners, increments back to user values
!
USE diffuse_mod
!
IMPLICIT NONE
!
IF(diff_l_friedel) THEN
   eck = diff_eck_u
   vi  = diff_vi_u
   inc = diff_inc_u
   num = diff_inc_u
ENDIF
!
END SUBROUTINE four_rese_friedel
!
!*******************************************************************************
!
SUBROUTINE four_fft_prep
!
USE crystal_mod
USE diffuse_mod
!
USE precision_mod
!
IMPLICIT NONE
!
fft_grid = 100
IF(ilots==LOT_OFF) THEN
   fft_dim(1) = cr_icc(1)*fft_grid
   fft_dim(2) = cr_icc(2)*fft_grid
ELSE
   fft_dim(1) = ls_xyz(1)*fft_grid
   fft_dim(2) = ls_xyz(2)*fft_grid
ENDIF
ALLOCATE(fft_field(1:fft_dim(1),1:fft_dim(2)))
ALLOCATE(fft_sum  (1:fft_dim(1),1:fft_dim(2)))
fft_field = CMPLX(0.0_PREC_DP, 0.0_PREC_DP)
fft_sum   = CMPLX(0.0_PREC_DP, 0.0_PREC_DP)
!
END SUBROUTINE four_fft_prep
!
!*******************************************************************************
!
SUBROUTINE four_strucf_fft(iscat, lform, nlot)
!
!
!
USE crystal_mod
USE diffuse_mod
!
USE singleton
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: iscat
LOGICAL, INTENT(IN) :: lform
!
INTEGER :: nlot
INTEGER :: loop
INTEGER :: i,j
!
fft_field = CMPLX(0.0_PREC_DP, 0.0_PREC_DP)
!
DO loop=1, nxat
   i = NINT(fft_grid*xat(loop,1))
   j = NINT(fft_grid*xat(loop,2))
   fft_field(i,j) = CMPLX(1.0_PREC_DP, 0.0_PREC_DP)
ENDDO
!
fft_field = fft(fft_field)/REAL(fft_dim(1)*fft_dim(2),KIND=PREC_DP)
fft_sum = fft_sum + fft_field*CONJG(fft_field)
!
END SUBROUTINE four_strucf_fft
!
!*******************************************************************************
!
SUBROUTINE four_fft_finalize
!
USE crystal_mod
!
IMPLICIT NONE
!
IF(ALLOCATED(fft_field)) DEALLOCATE(fft_field)
IF(ALLOCATED(fft_sum  )) DEALLOCATE(fft_sum  )
!
END SUBROUTINE four_fft_finalize
!
!*******************************************************************************
!
SUBROUTINE four_accumulate
!-
!  Accumulates the structure factor / intensity 
!+
USE diffuse_mod
!
USE errlist_mod
!
IMPLICIT NONE
!
IF(four_accum==FOUR_ACCUM_SINGLE) THEN
   IF(four_symm) THEN
      IF(ALLOCATED(csf_sum)) DEALLOCATE(csf_sum)
      IF(ALLOCATED(dsi_sum)) DEALLOCATE(dsi_sum)
      ALLOCATE(csf_sum(LBOUND(csf,1):UBOUND(csf,1)))
      csf_sum = CMPLX(0.0D0, 0.0D0)
      ALLOCATE(dsi_sum(LBOUND(dsi,1):UBOUND(dsi,1)))
      dsi_sum = 0.0D0
      csf_sum = csf
      dsi_sum = dsi
      CALL four_fill_csf                         ! Apply symm, if Set aver was used, fill Bragg
      CALL four_fill_dsi                         ! Apply symm, if Set aver was used, fill Bragg
      csf = csf_sum
      dsi = dsi_sum
   ELSE
      RETURN                                      ! Single diffracation pattern mode no symmetry average
   ENDIF
ELSEIF(four_accum==FOUR_ACCUM_INIT) THEN          ! Initialization, clear *_sum arays
   IF(ALLOCATED(csf_sum)) DEALLOCATE(csf_sum)
   IF(ALLOCATED(dsi_sum)) DEALLOCATE(dsi_sum)
ELSEIF(four_accum==FOUR_ACCUM_ACCUM) THEN         ! Accumulate into current sum
   IF(.NOT.ALLOCATED(csf_sum)) THEN
      ALLOCATE(csf_sum(LBOUND(csf,1):UBOUND(csf,1)))
      csf_sum = CMPLX(0.0D0, 0.0D0)
   ENDIF
   IF(.NOT.ALLOCATED(dsi_sum)) THEN
      ALLOCATE(dsi_sum(LBOUND(dsi,1):UBOUND(dsi,1)))
      dsi_sum = 0.0D0
   ENDIF
   csf_sum = csf_sum + csf
   dsi_sum = dsi_sum + dsi
ELSEIF(four_accum==FOUR_ACCUM_FINISHED) THEN      ! Finished, copy everything into the original
   IF(.NOT.ALLOCATED(csf_sum) .OR. .NOT.ALLOCATED(dsi_sum)) THEN
      ier_num = -175
      ier_typ = ER_APPL
      RETURN
   ENDIF
   CALL four_fill_csf                         ! Apply sym, if Set aver was used, fill Bragg
   CALL four_fill_dsi                         ! Apply sym, if Set aver was used, fill Bragg
   csf = csf_sum
   dsi = dsi_sum
   DEALLOCATE(csf_sum)
   DEALLOCATE(dsi_sum)
ENDIF
!
END SUBROUTINE four_accumulate
!
!*******************************************************************************
!
SUBROUTINE four_fill_csf
!-
!  If set aver was used with faver > 0, we need to fill the Bragg intensities
!+
USE diffuse_mod
!
IMPLICIT NONE
!
COMPLEX(KIND=KIND(0.0D0)), DIMENSION(:,:,:), ALLOCATABLE :: csf_3d
!
INTEGER :: h,k,l
INTEGER :: ii
INTEGER :: np
REAL(KIND=KIND(0.0E0)), DIMENSION(3) :: hkl
REAL(KIND=KIND(0.0E0)), DIMENSION(3) :: step
COMPLEX(KIND=KIND(0.0D0))            :: aaa
!
IF(fave==0.0 .AND. .NOT. four_symm) RETURN          ! nothing to do
!
step(1) = MAX(vi(1,1),vi(1,2), vi(1,3))*0.005      ! Get a step size that
step(2) = MAX(vi(2,1),vi(2,2), vi(2,3))*0.005      ! is used to check if hkl is integer
step(3) = MAX(vi(3,1),vi(3,2), vi(3,3))*0.005
IF(num(1)==1) step(1) = 0.0001
IF(num(2)==1) step(2) = 0.0001
IF(num(3)==1) step(3) = 0.0001
!
ALLOCATE(csf_3d(num(1), num(2), num(3)))
!
ii = 0
      DO h = 1, num (1) 
   DO k = 1, num (2) 
DO l = 1, num (3) 
         ii       = ii + 1 
         csf_3d(h,k,l) = csf_sum(ii)
      ENDDO
   ENDDO
ENDDO
!
IF(four_symm) CALL four_symm_csf(num, csf_3d, eck, vi)
!
!
IF(fave < 0.0) THEN
         DO h = 1, num (1) 
      DO k = 1, num (2) 
   DO l = 1, num (3) 
            hkl(1) = eck(1,1) + (h-1)*vi(1,1) + (k-1)*vi(1,2) + (l-1)*vi(1,3)
            hkl(2) = eck(2,1) + (h-1)*vi(2,1) + (k-1)*vi(2,2) + (l-1)*vi(2,3)
            hkl(3) = eck(3,1) + (h-1)*vi(3,1) + (k-1)*vi(3,2) + (l-1)*vi(3,3)
            IF(ABS(hkl(1)-NINT(hkl(1)))<step(1) .AND.  &
               ABS(hkl(2)-NINT(hkl(2)))<step(2) .AND.  &
               ABS(hkl(3)-NINT(hkl(3)))<step(3)       ) THEN       ! Got integer hkl
               aaa = 0.0
               np  = 0
               IF(h>1) THEN                        ! Add point at h-1
                  aaa = aaa + csf_3d(h-1,k,l)
                  np  = np + 1
               ENDIF
               IF(h<num(1)) THEN                   ! Add point at h+1
                  aaa = aaa + csf_3d(h+1,k,l)
                  np  = np + 1
               ENDIF
               IF(k>1) THEN                        ! Add point at h-1
                  aaa = aaa + csf_3d(h,k-1,l)
                  np  = np + 1
               ENDIF
               IF(k<num(2)) THEN                   ! Add point at h+1
                  aaa = aaa + csf_3d(h,k+1,l)
                  np  = np + 1
               ENDIF
               IF(l>1) THEN                        ! Add point at h-1
                  aaa = aaa + csf_3d(h,k,l-1)
                  np  = np + 1
               ENDIF
               IF(l<num(3)) THEN                   ! Add point at h+1
                  aaa = aaa + csf_3d(h,k,l+1)
                  np  = np + 1
               ENDIF
               IF(np>0) THEN
                  csf_3d(h,k,l) = aaa/REAL(np)
               ENDIF
            ENDIF
         ENDDO
      ENDDO
   ENDDO
ENDIF
!
ii = 0
      DO h = 1, num (1) 
   DO k = 1, num (2) 
DO l = 1, num (3) 
         ii       = ii + 1 
         csf_sum(ii) = csf_3d(h,k,l)
      ENDDO
   ENDDO
ENDDO
!
DEALLOCATE(csf_3d)
!
END SUBROUTINE four_fill_csf
!
!*******************************************************************************
!
SUBROUTINE four_fill_dsi
!-
!  If set aver was used with faver > 0, we need to fill the Bragg intensities
!+
USE diffuse_mod
!
IMPLICIT NONE
!
REAL(KIND=KIND(0.0D0))   , DIMENSION(:,:,:), ALLOCATABLE :: dsi_3d
!
INTEGER :: h,k,l
INTEGER :: ii
INTEGER :: np
REAL, DIMENSION(3) :: hkl
REAL, DIMENSION(3) :: step
REAL(KIND=KIND(0.0D0))               :: aaa
!
IF(fave==0.0 .AND. .NOT. four_symm) RETURN          ! nothing to do
!
step(1) = MAX(vi(1,1),vi(1,2), vi(1,3))*0.005
step(2) = MAX(vi(2,1),vi(2,2), vi(2,3))*0.005
step(3) = MAX(vi(3,1),vi(3,2), vi(3,3))*0.005
IF(num(1)==1) step(1) = 0.0001
IF(num(2)==1) step(2) = 0.0001
IF(num(3)==1) step(3) = 0.0001
!
ALLOCATE(dsi_3d(num(1), num(2), num(3)))
!
ii = 0
      DO h = 1, num (1) 
   DO k = 1, num (2) 
DO l = 1, num (3) 
         ii       = ii + 1 
         dsi_3d(h,k,l) = dsi_sum(ii)
      ENDDO
   ENDDO
ENDDO
!
IF(four_symm) CALL four_symm_dsi(num, dsi_3d, eck, vi)
!
IF(fave < 0.0) THEN
         DO h = 1, num (1) 
      DO k = 1, num (2) 
   DO l = 1, num (3) 
            hkl(1) = eck(1,1) + (h-1)*vi(1,1) + (k-1)*vi(1,2) + (l-1)*vi(1,3)
            hkl(2) = eck(2,1) + (h-1)*vi(2,1) + (k-1)*vi(2,2) + (l-1)*vi(2,3)
            hkl(3) = eck(3,1) + (h-1)*vi(3,1) + (k-1)*vi(3,2) + (l-1)*vi(3,3)
            IF(ABS(hkl(1)-NINT(hkl(1)))<step(1) .AND.  &
               ABS(hkl(2)-NINT(hkl(2)))<step(2) .AND.  &
               ABS(hkl(3)-NINT(hkl(3)))<step(3)       ) THEN       ! Got integer hkl
               aaa = 0.0
               np  = 0
               IF(h>1) THEN                        ! Add point at h-1
                  aaa = aaa + dsi_3d(h-1,k,l)
                  np  = np + 1
               ENDIF
               IF(h<num(1)) THEN                   ! Add point at h+1
                  aaa = aaa + dsi_3d(h+1,k,l)
                  np  = np + 1
               ENDIF
               IF(k>1) THEN                        ! Add point at h-1
                  aaa = aaa + dsi_3d(h,k-1,l)
                     np  = np + 1
               ENDIF
               IF(k<num(2)) THEN                   ! Add point at h+1
                  aaa = aaa + dsi_3d(h,k+1,l)
                  np  = np + 1
               ENDIF
               IF(l>1) THEN                        ! Add point at h-1
                  aaa = aaa + dsi_3d(h,k,l-1)
                  np  = np + 1
               ENDIF
               IF(l<num(3)) THEN                   ! Add point at h+1
                  aaa = aaa + dsi_3d(h,k,l+1)
                  np  = np + 1
               ENDIF
               IF(np>0) THEN
                  dsi_3d(h,k,l) = aaa/REAL(np)
               ENDIF
            ENDIF
         ENDDO
      ENDDO
   ENDDO
ENDIF
!
ii = 0
      DO h = 1, num (1) 
   DO k = 1, num (2) 
DO l = 1, num (3) 
         ii       = ii + 1 
         dsi_sum(ii) = dsi_3d(h,k,l)
      ENDDO
   ENDDO
ENDDO
!
DEALLOCATE(dsi_3d)
!
END SUBROUTINE four_fill_dsi
!
!*******************************************************************************
!
SUBROUTINE four_symm_csf(num, csf_3d, eck, vi)
!-
!   Apply reciprocal space symmetry of point group
!+
USE crystal_mod
use wyckoff_mod
!
USE matrix_mod
use precision_mod
!
IMPLICIT NONE
!
INTEGER                  , DIMENSION(3)                     , INTENT(IN)    :: num
COMPLEX(KIND=KIND(0.0D0)), DIMENSION(num(1), num(2), num(3)), INTENT(INOUT) :: csf_3d
REAL(kind=PREC_DP)       , DIMENSION(1:3, 1:4)           ::  eck
REAL(kind=PREC_DP)       , DIMENSION(1:3, 1:3)           ::  vi
!
INTEGER :: n_center     ! Only use this fraction of symmetry operations
INTEGER :: i
INTEGER :: h,k,l
INTEGER :: h1,k1,l1
INTEGER, DIMENSION(3) :: izero                                    ! Indices of point 0,0,0
REAL(KIND=PREC_DP), DIMENSION(3) :: hkl
REAL(KIND=PREC_DP), DIMENSION(3) :: uvw
COMPLEX(KIND=KIND(0.0D0)), DIMENSION(:,:,:), ALLOCATABLE :: csf_sym
INTEGER                  , DIMENSION(:,:,:), ALLOCATABLE :: weight
REAL(KIND=PREC_DP)       , DIMENSION(3, 3)               :: sym_mat
REAL(KIND=PREC_DP)       , DIMENSION(3, 3)               :: tmp_mat
!
!sym_mat(:,1) = vi(:,1)              ! Find pixels at which hkl = (0,0,0)
!sym_mat(:,2) = vi(:,2)
!sym_mat(:,3) = vi(:,3)
!sym_mat(:,3) = (/ 0.0, 0.0, 1.0/)
!
if(num(1)>1) then
   sym_mat(:,1) = vi(:,1)              ! Find pixels at which hkl = (0,0,0)
else
   sym_mat(:,1) = (/ 1.0, 0.0, 0.0/)
endif
!
if(num(2)>1) then
   sym_mat(:,2) = vi(:,2)              ! Find pixels at which hkl = (0,0,0)
else
   sym_mat(:,2) = (/ 0.0, 1.0, 0.0/)
endif
!
if(num(3)>1) then
   sym_mat(:,3) = vi(:,3)              ! Find pixels at which hkl = (0,0,0)
else
   sym_mat(:,3) = (/ 0.0, 0.0, 1.0/)
endif
!
CALL matinv3(sym_mat, tmp_mat)
!
hkl(:) = -eck(:,1)
uvw = MATMUL(tmp_mat, hkl)
do i=1,3
   if(num(i)>1) then
      izero(i) = NINT(1+uvw(i))            ! Pixel at hkl = (0,0,0)
   else
      izero(i) = 1
   endif
enddo
!
!                                     ! Determine number of primitive operations
!
n_center = 1
IF (cr_spcgr (1:1) .eq.'P') THEN
   n_center = 1
ELSEIF (cr_spcgr (1:1) .eq.'A') THEN  ! Normal space group can be used
   n_center = 2
ELSEIF (cr_spcgr (1:1) .eq.'B') THEN  ! as n_center is identical for
   n_center = 2
ELSEIF (cr_spcgr (1:1) .eq.'C') THEN  ! A B and C, orthorhombic alternative setting
   n_center = 2
ELSEIF (cr_spcgr (1:1) .eq.'I') THEN
   n_center = 2
ELSEIF (cr_spcgr (1:1) .eq.'F') THEN
   n_center = 4
ELSEIF (cr_spcgr (1:1) .eq.'R'.and.cr_syst.eq.6) THEN
   n_center = 3
ENDIF
!
ALLOCATE(csf_sym(num(1), num(2), num(3)))
ALLOCATE(weight (num(1), num(2), num(3)))      ! Weight counts the number that symm adds intensity to a pixel
weight = 1
!
csf_sym = csf_3d                     ! Effectively operation 1
!
!
DO i = 2, 1!spc_n/n_center             ! only apply numbers 2 for P centered part
   sym_mat=spc_mat(1:3,1:3,i)        ! Copy symmetry matrix to local copy
   DO l = 1, num (3) 
      hkl(3) = l - izero(3)
      DO k = 1, num (2) 
         hkl(2) = k - izero(2)
         DO h = 1, num (1) 
            hkl(1) = h - izero(1)
            uvw    = MATMUL(hkl, sym_mat)
            h1 = izero(1) + NINT(uvw(1))
            k1 = izero(2) + NINT(uvw(2))
            l1 = izero(3) + NINT(uvw(3))
            IF(h1>=1 .AND. h1<=num(1) .AND.          &
               k1>=1 .AND. k1<=num(2) .AND.          &
               l1>=1 .AND. l1<=num(3)      ) THEN
               csf_sym(h1,k1,l1) = csf_sym(h1, k1, l1) + csf_3d(h,k,l)
               weight (h1,k1,l1) = weight (h1, k1, l1) + 1
            ENDIF
         ENDDO
      ENDDO
   ENDDO
ENDDO
!
! Apply symmetry weight and divide by number of Symmetry operations
!
csf_3d = csf_sym/weight/REAL(spc_n/n_center)! Copy back into original
!
DEALLOCATE(csf_sym)
DEALLOCATE(weight )
!
END SUBROUTINE four_symm_csf
!
!*******************************************************************************
!
SUBROUTINE four_symm_dsi(num, dsi_3d, eck, vi)
!-
!   Apply reciprocal space symmetry of point group
!+
USE crystal_mod
use wyckoff_mod
!
use precision_mod
USE matrix_mod
!
IMPLICIT NONE
!
INTEGER               , DIMENSION(3)                     , INTENT(IN)    :: num
REAL(KIND=KIND(0.0D0)), DIMENSION(num(1), num(2), num(3)), INTENT(INOUT) :: dsi_3d
REAL(kind=PREC_DP)       , DIMENSION(1:3, 1:4)           ::  eck
REAL(kind=PREC_DP)       , DIMENSION(1:3, 1:3)           ::  vi
!
INTEGER :: n_center     ! Only use this fraction of symmetry operations
INTEGER :: i
INTEGER :: h,k,l
INTEGER :: h1,k1,l1
INTEGER, DIMENSION(3) :: izero                                    ! Indices of point 0,0,0
REAL(KIND=PREC_DP), DIMENSION(3) :: hkl
REAL(KIND=PREC_DP), DIMENSION(3) :: uvw
REAL(KIND=KIND(0.0D0)), DIMENSION(:,:,:), ALLOCATABLE :: dsi_sym
INTEGER               , DIMENSION(:,:,:), ALLOCATABLE :: weight
REAL(KIND=PREC_DP), DIMENSION(3, 3)                   :: sym_mat
REAL(KIND=PREC_DP), DIMENSION(3, 3)                   :: tmp_mat
!
if(num(1)>1) then
   sym_mat(:,1) = vi(:,1)              ! Find pixels at which hkl = (0,0,0)
else
   sym_mat(:,1) = (/ 1000.0, 0.0, 0.0/)
endif
!
if(num(2)>1) then
   sym_mat(:,2) = vi(:,2)              ! Find pixels at which hkl = (0,0,0)
else
   sym_mat(:,2) = (/ 0.0, 1000.0, 0.0/)
endif
!
if(num(3)>1) then
   sym_mat(:,3) = vi(:,3)              ! Find pixels at which hkl = (0,0,0)
else
   sym_mat(:,3) = (/ 0.0, 0.0, 1000.0/)
endif
!
CALL matinv3(sym_mat, tmp_mat)
!
hkl(:) = -eck(:,1)
uvw = MATMUL(tmp_mat, hkl)
do i=1,3
   if(num(i)>1) then
      izero(i) = NINT(1+uvw(i))            ! Pixel at hkl = (0,0,0)
   else
      izero(i) = 1
   endif
enddo
!
!                                     ! Determine number of primitive operations
!write(*,*) ' IZERO ', izero
!
n_center = 1
IF (cr_spcgr (1:1) .eq.'P') THEN
   n_center = 1
ELSEIF (cr_spcgr (1:1) .eq.'A') THEN  ! Normal space group can be used
   n_center = 2
ELSEIF (cr_spcgr (1:1) .eq.'B') THEN  ! as n_center is identical for
   n_center = 2
ELSEIF (cr_spcgr (1:1) .eq.'C') THEN  ! A B and C, orthorhombic alternative setting
   n_center = 2
ELSEIF (cr_spcgr (1:1) .eq.'I') THEN
   n_center = 2
ELSEIF (cr_spcgr (1:1) .eq.'F') THEN
   n_center = 4
ELSEIF (cr_spcgr (1:1) .eq.'R'.and.cr_syst.eq.6) THEN
   n_center = 3
ENDIF
!
ALLOCATE(dsi_sym(num(1), num(2), num(3)))
ALLOCATE(weight (num(1), num(2), num(3)))      ! Weight counts the number that symm adds intensity to a pixel
weight = 1
!
dsi_sym = dsi_3d                     ! Effectively operation 1
!
DO i = 2, spc_n/n_center             ! only apply numbers 2 for P centered part
   sym_mat=spc_mat(1:3,1:3,i)        ! Copy symmetry matrix to local copy
!write(*,'(a, 3F6.1)') ' SYMM ', sym_mat(1,:)
!write(*,'(a, 3F6.1)') ' SYMM ', sym_mat(2,:)
!write(*,'(a, 3F6.1)') ' SYMM ', sym_mat(3,:)
   DO l = 1, num (3) 
      hkl(3) = l - izero(3)
      DO k = 1, num (2) 
         hkl(2) = k - izero(2)
         DO h = 1, num (1) 
            hkl(1) = h - izero(1)
            uvw    = MATMUL(hkl, sym_mat)
            h1 = izero(1) + NINT(uvw(1))
            k1 = izero(2) + NINT(uvw(2))
            l1 = izero(3) + NINT(uvw(3))
            IF(h1>=1 .AND. h1<=num(1) .AND.          &
               k1>=1 .AND. k1<=num(2) .AND.          &
               l1>=1 .AND. l1<=num(3)      ) THEN
               dsi_sym(h1,k1,l1) = dsi_sym(h1, k1, l1) + dsi_3d(h,k,l)
               weight (h1,k1,l1) = weight (h1, k1, l1) + 1
            ENDIF
         ENDDO
      ENDDO
   ENDDO
ENDDO
!
! Apply symmetry weight and divide by number of Symmetry operations
!
dsi_3d = dsi_sym/weight/REAL(spc_n/n_center,kind=PREC_DP)! Copy back into original
!
DEALLOCATE(dsi_sym)
DEALLOCATE(weight )
!
END SUBROUTINE four_symm_dsi
!
!*******************************************************************************
!
END MODULE fourier_sup
