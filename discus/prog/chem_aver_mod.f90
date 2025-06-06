MODULE chem_aver_mod
!
CONTAINS
!*****7*****************************************************************
SUBROUTINE chem_aver (lout, lsite) 
!+                                                                      
!     Calculate average structure and standard deviation                
!-                                                                      
USE discus_config_mod 
USE discus_allocate_appl_mod 
USE crystal_mod 
use dis_estimate_mod
USE atom_env_mod
USE atom_name
USE chem_mod 
USE errlist_mod 
!
use metric_mod
USE param_mod 
use precision_mod
USE prompt_mod 
USE lib_f90_allocate_mod
!
IMPLICIT none 
!                                                                       
LOGICAL, INTENT(IN) :: lout    ! Print output if true
LOGICAL, INTENT(IN) :: lsite   ! Treat different atoms on each site as one 
!                                                                       
!logical, parameter :: lspace =.TRUE.
!
REAL(kind=PREC_DP), DIMENSION(3) ::  p , ez
INTEGER, DIMENSION(3) :: iez
!
INTEGER            :: i, j, k, ii, jj, kk, ia, is, nvalues
!integer            :: i1, i2, i3
INTEGER            :: n_res
LOGICAL            :: flag
!                                                                       
CHARACTER(LEN=9)   :: at_name_i 
!
INTEGER            :: n_atom_cell  ! Dummy for allocation
INTEGER            :: n_max_atom   ! Dummy for allocation
!
INTEGER, DIMENSION(:  ), ALLOCATABLE :: temp_chem_ave_iscat
INTEGER, DIMENSION(:  ), ALLOCATABLE :: temp_chem_ave_anis
REAL(kind=PREC_DP)   , DIMENSION(:,:), ALLOCATABLE :: temp_chem_ave_posit
REAL(kind=PREC_DP)   , DIMENSION(:,:), ALLOCATABLE :: temp_chem_ave_sigma
REAL(kind=PREC_DP)   , DIMENSION(  :), ALLOCATABLE :: temp_chem_ave_bese
REAL(kind=PREC_DP)   , DIMENSION(3  )              :: temp_chem_ave_pos
REAL(kind=PREC_DP)   , DIMENSION(3  )              :: temp_chem_ave_sig
!
!
integer, dimension(3)               :: ncell_out
!integer, dimension(3)               :: iic
real(kind=PREC_DP), dimension(3, 2) :: pdt_dims
integer           , dimension(3)    :: pdt_ilow        ! Unit cell dimensions in periodic
integer           , dimension(3)    :: pdt_ihig        ! low and high inidce
integer                             :: pdt_ncells     ! Number of cells in periodic crystal volume
integer                                           :: pdt_nsite    ! Number of sites in an averaged unit cell
integer                                           :: pdt_asite ! Achieved sites either from set site or find_aver
integer           , dimension(:,  :), allocatable :: pdt_itype    ! Atom types at each site
real(kind=PREC_DP), dimension(:,:)  , allocatable :: pdt_pos  ! Atom positions in average cell
!real(kind=PREC_DP), dimension(3,4*nint(aver)) :: pdt_temp      ! Positions of the average cell
integer           , dimension(:)    , allocatable :: at_site  ! Atom is at this site
REAL(kind=PREC_DP)   , DIMENSION(3  )              :: uvw, shift
!real(kind=PREC_DP)                                 :: dmin        ! minimum distance to average site
!real(kind=PREC_DP)                                 :: dist        ! minimum distance to average site
!real(kind=PREC_DP)                                 :: eps         ! minimum distance to average site
!integer :: istart, ic

real(kind=PREC_DP) :: aver
real(kind=PREC_DP) :: sigma
!
if(cr_icc(1)*cr_icc(2)*cr_icc(3)*cr_ncatoms>cr_natoms) then
   ier_num = -179
   ier_typ= ER_APPL
   ier_msg(1) = 'More atoms in the crystal than aver expects'
   ier_msg(2) = 'Check unit cell, atoms per unit cell, '
   ier_msg(3) = 'and total atom number'
   return
endif
is = 1
!
IF ( CHEM_MAX_AVE_ATOM < MAX(cr_ncatoms, MAXSCAT) .or. &
     CHEM_MAXAT_CELL   <= cr_ncatoms                    ) THEN
   n_atom_cell = MAX(CHEM_MAXAT_CELL,             cr_ncatoms)
   n_max_atom  = MAX(CHEM_MAX_AVE_ATOM, cr_ncatoms, MAXSCAT) + 1
   call alloc_chem_aver ( n_atom_cell, n_max_atom)
ENDIF
IF(.not. lsite) THEN
!  n_atom_cell = MAX(CHEM_MAXAT_CELL, MAXAT_CELL)
   n_atom_cell = MAX(CHEM_MAXAT_CELL, cr_ncatoms)
   IF(ALLOCATED(chem_ave_posit)) DEALLOCATE(chem_ave_posit)
   IF(ALLOCATED(chem_ave_sigma)) DEALLOCATE(chem_ave_sigma)
   ALLOCATE(chem_ave_posit(3,n_atom_cell, MAX(12,cr_nscat)))
   ALLOCATE(chem_ave_sigma(3,n_atom_cell, MAX(12,cr_nscat)))
   chem_ave_posit = 0.0d0
   chem_ave_sigma = 0.0D0
ENDIF
!                                                                       
!------ reset counters                                                  
!                                                                       
IF (cr_ncatoms.gt.CHEM_MAXAT_CELL) then 
   ier_num = - 103 
   ier_typ = ER_APPL 
   ier_msg (1) = 'Adjust the value of the variable' 
   ier_msg (2)  = 'CHEM_MAXAT_CELL in chem_mod.f90 and  ' 
   ier_msg (3)  = 'compile the program             ' 
   RETURN 
ENDIF 
chem_ave_n    = 0    ! (i)   , i=1,cr_ncatoms
chem_ave_bese = 0.0D0  ! (i, k), i=1,cr_ncatoms, k = 1,chem_max_ave_atom
chem_ave_pos  = 0.0D0  ! (j, i), i=1,cr_ncatoms, j = 1,3
chem_ave_sig  = 0.0D0  ! (j, i), i=1,cr_ncatoms, j = 1,3
!
cond_quick: if(chem_quick) then     ! Fast mode
!                                                                       
!------ loop over all unit cells and atoms within unit cell             
!                                                                       
iez(1) = NINT(cr_dim(1,1)) - 1
iez(2) = NINT(cr_dim(2,1)) - 1
iez(3) = NINT(cr_dim(3,1)) - 1
!
loopk: DO k = 1, cr_icc (3) 
   loopj: DO j = 1, cr_icc (2) 
      loopi: DO i = 1, cr_icc (1) 
         ez (1) = REAL(iez(1) + i )
         ez (2) = REAL(iez(2) + j )
         ez (3) = REAL(iez(3) + k )
         loopii: DO ii = 1, cr_ncatoms 
            ia = ( (k - 1) * cr_icc (1) * cr_icc (2) + &
                   (j - 1) * cr_icc (1) + (i - 1) ) * cr_ncatoms + ii                                     
            DO jj = 1, 3 
               p (jj) = cr_pos (jj, ia) - ez (jj) 
               chem_ave_pos (jj, ii) = chem_ave_pos (jj, ii) + p (jj) 
!              chem_ave_sig (jj, ii) = chem_ave_sig (jj, ii) + p (jj)**2 
            ENDDO 
!                                                                       
!------ --- Calculate occupancies ..                                    
!                                                                       
            occup: IF (chem_ave_n (ii) .eq.0) then 
               chem_ave_n (ii) = 1 
               chem_ave_iscat (ii, chem_ave_n (ii) ) = cr_iscat (1, ia) 
               chem_ave_anis  (ii, chem_ave_n (ii) ) = cr_iscat (3, ia) 
               is = 1 
            ELSE  occup
               flag = .true. 
               DO kk = 1, chem_ave_n (ii) 
                  IF (cr_iscat (1,ia) .eq.chem_ave_iscat (ii, kk) ) then 
                     is = kk 
                     flag = .false. 
                  ENDIF 
               ENDDO 
               IF (flag) then 
                  chem_ave_n (ii) = chem_ave_n (ii) + 1 
                  is = chem_ave_n (ii) 
                  IF (chem_ave_n (ii) .gt.CHEM_MAX_AVE_ATOM) then 
                     ier_typ = ER_CHEM 
                     ier_num = - 5 
                     RETURN 
                  ENDIF 
                  chem_ave_iscat (ii, chem_ave_n (ii) ) = cr_iscat (1,ia) 
                  chem_ave_anis  (ii, chem_ave_n (ii) ) = cr_iscat (3,ia) 
               ENDIF 
            ENDIF occup
            IF(.not. lsite) THEN   ! Accumulate individual positions for different atoms
               DO jj = 1, 3 
                  chem_ave_posit(jj,ii,is) = chem_ave_posit(jj,ii,is) + p(jj)
!                 chem_ave_sigma(jj,ii,is) = chem_ave_sigma(jj,ii,is) + p(jj)**2
               ENDDO 
            ENDIF
            chem_ave_bese (ii, is) = chem_ave_bese (ii, is) + 1 
         ENDDO loopii
      ENDDO  loopi
   ENDDO  loopj
ENDDO  loopk
!
! Calculate average positions
!
ia = cr_icc (1) * cr_icc (2) * cr_icc (3) 
IF(lsite) THEN           ! one pos for all types on a single site
   DO ii = 1, cr_ncatoms
      DO jj =1, 3
         chem_ave_pos (jj, ii) = chem_ave_pos (jj, ii)/ REAL(ia)
      ENDDO
   ENDDO
ELSE
   DO i = 1, cr_ncatoms 
      DO k = 1, chem_ave_n (i) 
         DO j = 1, 3 
            chem_ave_posit (j, i, k) = chem_ave_posit (j, i, k) / chem_ave_bese(i,k) 
         ENDDO
      ENDDO
   ENDDO
ENDIF
!
! Now calculate sigmas
!
sloopk: DO k = 1, cr_icc (3) 
   sloopj: DO j = 1, cr_icc (2) 
      sloopi: DO i = 1, cr_icc (1) 
         ez (1) = REAL(iez(1) + i )
         ez (2) = REAL(iez(2) + j )
         ez (3) = REAL(iez(3) + k )
         sloopii: DO ii = 1, cr_ncatoms 
            ia = ( (k - 1) * cr_icc (1) * cr_icc (2) + &
                   (j - 1) * cr_icc (1) + (i - 1) ) * cr_ncatoms + ii                                     
            IF(lsite) THEN                 ! No distinction of atom types
               DO jj = 1, 3 
                  p (jj) = cr_pos (jj, ia) - ez (jj) 
                  chem_ave_sig (jj, ii) = chem_ave_sig (jj, ii) + &
                         (p (jj) - chem_ave_pos (jj, ii))**2
               ENDDO 
            ELSE
               is = 1
               s_site:DO kk = 1, chem_ave_n (ii) 
                  IF (cr_iscat (1,ia) .eq.chem_ave_iscat (ii, kk) ) then 
                     is = kk 
                     EXIT s_site
                  ENDIF
               ENDDO s_site
               DO jj = 1, 3 
                  p (jj) = cr_pos (jj, ia) - ez (jj) 
                  chem_ave_sigma (jj, ii, is) = chem_ave_sigma (jj, ii, is) + &
                         (p (jj) - chem_ave_posit (jj, ii, is))**2
               ENDDO 
            ENDIF
         ENDDO sloopii
      ENDDO  sloopi
   ENDDO  sloopj
ENDDO  sloopk
!
else  cond_quick           ! Exact mode
   call estimate_ncells(ncell_out, pdt_dims, pdt_ilow, pdt_ihig, pdt_ncells)
   if(ier_num/=0) return
   call estimate_ncatom(aver, sigma, pdt_ilow, pdt_ihig, pdt_ncells)
   if(ier_num/=0) return
   chem_ave_n = 1
   allocate(at_site(1:cr_natoms))
   at_site = 0
   call find_average(aver, pdt_dims, pdt_ncells, pdt_nsite, pdt_asite, &
        pdt_itype, pdt_pos, cr_natoms, at_site)
   if(ier_num/=0) return
   chem_ave_pos = pdt_pos
   shift(:) = REAL(-NINT(cr_dim(:,1))+2)    ! Shift to make atom position positive
      chem_ave_n = 0        ! Clear number of atom on a site , reused as number of atom types on the site
!                                                                       
!------ --- Calculate occupancies ..                                    
!                                                                       
   do ia=1, cr_natoms
      ii = at_site(ia)
   occup_exact: IF (chem_ave_n (ii) .eq.0) then 
               chem_ave_n (ii) = 1 
               chem_ave_iscat (ii,               1 ) = cr_iscat (1, ia) 
               chem_ave_anis  (ii,               1 ) = cr_iscat (3, ia) 
               is = 1 
            ELSE  occup_exact
               flag = .true. 
               DO kk = 1, chem_ave_n (ii) 
                  IF (cr_iscat (1,ia) .eq.chem_ave_iscat (ii, kk) ) then 
                     is = kk 
                     flag = .false. 
                  ENDIF 
               ENDDO 
               IF (flag) then 
                  chem_ave_n (ii) = chem_ave_n (ii) + 1 
                  is = chem_ave_n (ii) 
                  IF (chem_ave_n (ii) .gt.CHEM_MAX_AVE_ATOM) then 
                     ier_typ = ER_CHEM 
                     ier_num = - 5 
                     RETURN 
                  ENDIF 
                  chem_ave_iscat (ii, chem_ave_n (ii) ) = cr_iscat (1,ia) 
                  chem_ave_anis  (ii, chem_ave_n (ii) ) = cr_iscat (3,ia) 
               ENDIF 
            ENDIF occup_exact
            IF(.not. lsite) THEN   ! Accumulate individual positions for different atoms
               DO jj = 1, 3 
                  chem_ave_posit(jj,ii,is) = chem_ave_posit(jj,ii,is) + p(jj)
!                 chem_ave_sigma(jj,ii,is) = chem_ave_sigma(jj,ii,is) + p(jj)**2
               ENDDO 
            ENDIF
            chem_ave_bese (ii, is) = chem_ave_bese (ii, is) + 1 
!
!     Calculate sigma
!
         uvw(:) = (cr_pos(:,ia)+shift(:)+ 00.0) - real(int(cr_pos(:,ia)+shift(:)+ 00.0), PREC_DP)   ! Create fractional
         do jj=1, 3
            if(uvw(jj)-chem_ave_pos(jj,ii)<-0.50_PREC_DP) then
               uvw(jj) = uvw(jj) + 1.0_PREC_DP
            elseif(uvw(jj)-chem_ave_pos(jj,ii)>0.50_PREC_DP) then
               uvw(jj) = uvw(jj) - 1.0_PREC_DP
            endif
         enddo
      if(lsite) then                 ! No distinction of atom types
         chem_ave_sig(:, ii) = chem_ave_sig(:, ii) + &
                         (uvw    - chem_ave_pos (:, ii))**2
      else
         is = 1
         s_site_exact:DO kk = 1, chem_ave_n (ii) 
            IF (cr_iscat (1,ia) .eq.chem_ave_iscat (ii, kk) ) then 
               is = kk 
               EXIT s_site_exact
            ENDIF
         ENDDO s_site_exact
            chem_ave_sigma(:, ii, is) = chem_ave_sigma (:, ii, is) + &
                         (uvw    - chem_ave_posit (:, ii, is))**2
      endif
   enddo
endif cond_quick
!
!------ Sort atom types on a given site
!
DO i=1, cr_ncatoms                  ! Loop over all sites
   IF(chem_ave_n(i) > 0) THEN      ! We have atoms on this site
      IF(ALLOCATED(temp_chem_ave_iscat)) DEALLOCATE(temp_chem_ave_iscat)
      IF(ALLOCATED(temp_chem_ave_anis )) DEALLOCATE(temp_chem_ave_anis )
      IF(ALLOCATED(temp_chem_ave_bese )) DEALLOCATE(temp_chem_ave_bese )
      ALLOCATE(temp_chem_ave_iscat(   chem_ave_n(i)))
      ALLOCATE(temp_chem_ave_anis (   chem_ave_n(i)))
      ALLOCATE(temp_chem_ave_bese (   chem_ave_n(i)))
      IF(.NOT. lsite) THEN
         IF(ALLOCATED(temp_chem_ave_posit)) DEALLOCATE(temp_chem_ave_posit)
         IF(ALLOCATED(temp_chem_ave_sigma)) DEALLOCATE(temp_chem_ave_sigma)
         ALLOCATE(temp_chem_ave_posit(3, chem_ave_n(i)))
         ALLOCATE(temp_chem_ave_sigma(3, chem_ave_n(i)))
         temp_chem_ave_posit =  0.0
         temp_chem_ave_sigma =  0.0
      ENDIF
   ENDIF
   temp_chem_ave_iscat =  0
   temp_chem_ave_anis  =  0
   temp_chem_ave_pos   =  0.0
   temp_chem_ave_sig   =  0.0
   temp_chem_ave_bese  =  0.0
   DO k = 1, chem_ave_n (i) 
      j = MINLOC(chem_ave_iscat(i, 1:chem_ave_n (i)),dim=1)
      temp_chem_ave_iscat(  k) = chem_ave_iscat(  i,j)
      temp_chem_ave_anis (  k) = chem_ave_anis (  i,j)
      temp_chem_ave_pos  (:  ) = chem_ave_pos  (:,i  )
      temp_chem_ave_sig  (:  ) = chem_ave_sig  (:,i  )
      temp_chem_ave_bese (  k) = chem_ave_bese (  i,j)
      IF(.NOT. lsite) THEN
         temp_chem_ave_posit(:,k) = chem_ave_posit(:,i,j)
         temp_chem_ave_sigma(:,k) = chem_ave_sigma(:,i,j)
      ENDIF
      chem_ave_iscat(i,j) = 9999
      chem_ave_anis (i,j) = 9999
   ENDDO
   DO k = 1, chem_ave_n (i) 
      chem_ave_iscat(  i,k) = temp_chem_ave_iscat(  k)
      chem_ave_anis (  i,k) = temp_chem_ave_anis (  k)
      chem_ave_pos  (:,i  ) = temp_chem_ave_pos  (:  )
      chem_ave_sig  (:,i  ) = temp_chem_ave_sig  (:  )
      chem_ave_bese (  i,k) = temp_chem_ave_bese (  k)
      IF(.NOT. lsite) THEN
         chem_ave_posit(:,i,k) = temp_chem_ave_posit(:,k)
         chem_ave_sigma(:,i,k) = temp_chem_ave_sigma(:,k)
      ENDIF
   ENDDO
ENDDO
!                                                                       
!------ output of average and sigma                                     
!                                                                       
IF (lout) write (output_io, 1000) 
ia = cr_icc (1) * cr_icc (2) * cr_icc (3) 
IF(lsite) THEN           ! one pos for all types on a single site
   nvalues = 0
   DO i = 1, cr_ncatoms 
      DO j = 1, 3 
         chem_ave_sig (j, i) = chem_ave_sig (j, i) / REAL(ia) !-          &
         IF (chem_ave_sig (j, i) .gt.0.0) then 
            chem_ave_sig (j, i) = sqrt (chem_ave_sig (j, i) ) 
         ELSE 
            chem_ave_sig (j, i) = 0.0 
         ENDIF 
      ENDDO 
      IF (lout) then 
         DO k = 1, chem_ave_n (i) 
            at_name_i = at_name (chem_ave_iscat (i, k) ) 
            WRITE (output_io, 1100) i, at_name_i, (chem_ave_pos (ii, i),   &
            ii = 1, 3), (chem_ave_sig (ii, i), ii = 1, 3), chem_ave_bese ( &
            i, k) / ia                                                     
            nvalues = nvalues + 1
         ENDDO 
      ENDIF 
   ENDDO 
!                                                                       
!------ store results in res_para                                       
!                                                                       
   IF ( (9 * nvalues) .gt.MAXPAR_RES) then 
      n_res = MAX(9 * cr_ncatoms,MAXPAR_RES, MAX_ATOM_ENV)
      CALL alloc_param(n_res)
      MAXPAR_RES = n_res
   ENDIF 
      res_para (0) = 9 * nvalues 
      ii = 0
      DO i = 1, cr_ncatoms 
         DO k = 1, chem_ave_n (i) 
         ii = ii + 1
            res_para ( (ii - 1) * 9     + 1) = i
            res_para ( (ii - 1) * 9     + 2) = chem_ave_iscat (i, k)
         DO j = 1, 3 
            res_para ( (ii - 1) * 9 + j + 2) = chem_ave_pos (j, i) 
         ENDDO 
         DO j = 1, 3 
            res_para ( (ii - 1) * 9 + j + 5) = chem_ave_sig (j, i) 
         ENDDO 
            res_para ( (ii - 1) * 9     + 9) = chem_ave_bese(i, k)/ia
         ENDDO 
      ENDDO 
ELSE
   nvalues = 0
   DO i = 1, cr_ncatoms 
      DO k = 1, chem_ave_n (i) 
         DO j = 1, 3 
            chem_ave_sigma (j, i, k) = chem_ave_sigma (j, i, k) / chem_ave_bese(i,k) !- &
            IF (chem_ave_sigma (j, i, k) .gt.0.0) then 
               chem_ave_sigma (j, i, k) = sqrt (chem_ave_sigma (j, i, k) ) 
            ELSE 
               chem_ave_sigma (j, i, k) = 0.0 
            ENDIF 
         ENDDO 
         IF (lout) then 
           at_name_i = at_name (chem_ave_iscat (i, k) ) 
           WRITE (output_io, 1100) i, at_name_i, (chem_ave_posit (ii, i, k), ii = 1, 3),&
                 (chem_ave_sigma (ii, i, k), ii = 1, 3), chem_ave_bese(i, k) / ia
         ENDIF 
         nvalues = nvalues + 1
      ENDDO 
   ENDDO 
!                                                                       
!------ store results in res_para                                       
!                                                                       
   IF ( (9 * nvalues) .gt.MAXPAR_RES) then 
      n_res = MAX(9 * cr_ncatoms,MAXPAR_RES, MAX_ATOM_ENV)
      CALL alloc_param(n_res)
      MAXPAR_RES = n_res
   ENDIF
      res_para (0) = 9 * nvalues 
      ii = 0
      DO i = 1, cr_ncatoms 
         DO k = 1, chem_ave_n (i) 
         ii = ii + 1
            res_para ( (ii - 1) * 9     + 1) = i
            res_para ( (ii - 1) * 9     + 2) = chem_ave_iscat (i, k)
         DO j = 1, 3 
            res_para ( (ii - 1) * 9 + j + 2) = chem_ave_posit (j, i, k) 
         ENDDO 
         DO j = 1, 3 
            res_para ( (ii - 1) * 9 + j + 5) = chem_ave_sigma (j, i, k) 
         ENDDO 
            res_para ( (ii - 1) * 9     + 9) = chem_ave_bese(i, k)/ia
         ENDDO 
      ENDDO 
ENDIF
!                                                                       
 1000 FORMAT (' Average structure : ',//,                               &
     &        3x,'Site',2x,'atom',11x,'average position',8x,            &
     &        'standard deviation',3x,'occupancy',/,3x,75('-'))         
 1100 FORMAT (3x,i3,2x,a9,2x,3(f7.4,1x),1x,3(f7.4,1x),1x,f7.4) 
      END SUBROUTINE chem_aver                      
!*****7*****************************************************************
SUBROUTINE chem_elem (lout) 
!+                                                                      
!     Show information about elements/rel. amounts within crystal       
!-                                                                      
      USE crystal_mod 
      USE atom_name
USE atom_env_mod
      USE chem_mod 
!                                                                       
      USE errlist_mod 
      USE lib_f90_allocate_mod 
      USE param_mod 
      USE prompt_mod 
      USE wink_mod 
use precision_mod
      IMPLICIT none 
!
LOGICAL, INTENT(IN) :: lout 
!                                                                       
REAL(kind=PREC_DP)             :: proz 
INTEGER          :: i 
INTEGER          :: n_res 
!                                                                       
CHARACTER(LEN=9) :: at_name_i 
!                                                                       
!     Error condition                                                   
!                                                                       
IF (cr_natoms.eq.0) then 
   ier_typ = ER_CHEM 
   ier_num = - 27 
   RETURN 
ENDIF 
!                                                                       
!------ reset counters, ...                                             
!                                                                       
cr_amount(:)    = 0
cr_n_real_atoms = 0
cr_u2aver       = 0.0
!                                                                       
!------ get size of model crystal, rel. amount of elements              
!                                                                       
DO i = 1, cr_natoms 
   cr_amount (cr_iscat (1,i) ) = cr_amount (cr_iscat (1,i) ) + 1 
ENDDO 
!                                                                       
!------ write output                                                    
!                                                                       
IF (lout) write (output_io, 1000) (cr_icc (i), i = 1, 3) 
IF (lout) write (output_io, 1100) cr_natoms, cr_ncatoms, cr_nscat 
IF (i >  MAXPAR_RES) THEN 
   n_res = MAX(cr_nscat, MAXPAR_RES, MAX_ATOM_ENV)
   CALL alloc_param(n_res)
   MAXPAR_RES = n_res
ENDIF
res_para (0) = REAL(cr_nscat) + 1 
DO i = 0, cr_nscat 
   cr_u2aver = cr_u2aver + cr_dw(i) * cr_amount (i)
!   IF (cr_at_lis(cr_iscat (1,i)) /= 'VOID') then
   IF (cr_at_lis(          i ) /= 'VOID') then
      cr_n_real_atoms = cr_n_real_atoms + cr_amount(i)
   ENDIF
   proz = REAL(cr_amount (i) ) / cr_natoms 
   IF (lout) then 
      at_name_i = at_name (i) 
      WRITE (output_io, 1200) at_name_i, proz, cr_amount (i) 
   ENDIF 
   IF (i <= MAXPAR_RES) then 
      res_para (i + 1) = proz 
   ELSE 
      ier_typ = ER_CHEM 
      ier_num = - 2 
   ENDIF 
ENDDO 
cr_u2aver = cr_u2aver/cr_n_real_atoms/8./REAL(pi**2)
!                                                                       
 1000 FORMAT     (' Size of the crystal (unit cells) : ',2(I4,' x '),I4) 
 1100 FORMAT     (' Total number of atoms            : ',I9,/           &
     &                   ' Number of atoms per unit cell    : ',I9,/    &
     &                   ' Number of different atoms        : ',I9,/)   
 1200 FORMAT     ('    Element : ',A9,' rel. abundance : ',F5.3,        &
     &                   '  (',I9,' atoms)')                            
END SUBROUTINE chem_elem                      
!*****7*****************************************************************
SUBROUTINE chem_com (com,lout)
!
! Calculate center of mass for crystal
!
USE crystal_mod
USE prompt_mod
use precision_mod
!
IMPLICIT NONE
!
REAL(kind=PREC_DP), DIMENSION(3), INTENT(OUT) :: com   ! Center of Mass
LOGICAL,            INTENT(IN)  :: lout  ! Screen output flag
!
INTEGER :: i
!
com = 0.0
!
IF(cr_natoms > 0) THEN
   DO i=1,cr_natoms
      com(1) = com(1) + cr_pos(1,i)
      com(2) = com(2) + cr_pos(2,i)
      com(3) = com(3) + cr_pos(3,i)
   ENDDO
   com(1) = com(1)/cr_natoms
   com(2) = com(2)/cr_natoms
   com(3) = com(3)/cr_natoms
!
   IF(lout) THEN
      WRITE(output_io, 1000) com
   ENDIF
ENDIF
!
1000 FORMAT('Center of mass at ',3(2x,F12.3))
!
END SUBROUTINE chem_com 
!
!*******************************************************************************
!
SUBROUTINE get_displacement(line, length)
!
! Calculate the displacement of an atom from its average site
!
!
USE crystal_mod
USE atom_name
USE chem_mod
USE celltoindex_mod
!
USE errlist_mod
USE ber_params_mod
USE get_params_mod
USE precision_mod
USE prompt_mod
USE param_mod
USE precision_mod
USE take_param_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(INOUT) :: line
INTEGER         , INTENT(INOUT) :: length
!
INTEGER, PARAMETER :: MAXW = 4
CHARACTER(LEN=MAX(PREC_STRING,LEN(line))), DIMENSION(MAXW) :: cpara
INTEGER            , DIMENSION(MAXW) :: lpara
REAL(KIND=PREC_DP) , DIMENSION(MAXW) :: werte
!
INTEGER, PARAMETER :: NOPTIONAL = 3
INTEGER, PARAMETER :: O_AVER    = 1
INTEGER, PARAMETER :: O_INDI    = 2
INTEGER, PARAMETER :: O_ECHO    = 3
CHARACTER(LEN=   4), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=MAX(PREC_STRING,LEN(line))), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent!opt. para is present
REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 0 ! Number of values to calculate 
!
DATA oname  / 'aver', 'indi', 'out'   /
DATA loname /  4    ,  4    ,  3      /
!
CHARACTER(LEN=9)   :: at_name_d
!LOGICAL, PARAMETER :: LOLD = .FALSE.
LOGICAL, PARAMETER :: LOUT = .FALSE.
LOGICAL            :: lsite = .TRUE.    ! Average everything onto one site
INTEGER               :: iatom   ! Atom index
INTEGER, DIMENSION(3) :: icell   ! Cell number
INTEGER               :: isite   ! site number
!
INTEGER :: k,kk
INTEGER :: ianz
!
CALL get_params (line, ianz, cpara, lpara, MAXW, length)
IF(ier_num/=0) RETURN
if(ianz<1 .or. ianz>MAXW) then
   ier_num = -1
   ier_typ = ER_COMM
   return
endif
!
opara  =  (/ 'no'    , 'no'   , 'no' /)   ! Always provide fresh default values
lopara =  (/  2,        2     ,  2   /)
owerte =  (/  0.0,      0.0   ,  0.0 /)
!
CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
     oname, loname, opara, lopara, lpresent, owerte)
!
lsite = .TRUE.
IF(opara(O_INDI) == 'yes') lsite = .FALSE.               ! User requested individual sites
!
IF(lsite) THEN                 ! All atoms on a site share their average position
   IF(chem_run_aver .AND. .NOT.opara(O_AVER)=='off') THEN
      CALL chem_aver(LOUT, lsite)
      chem_run_aver = .FALSE.
   ENDIF
ELSE                           ! Individual average positions for a site
   IF(chem_run_aver_ind .AND. .NOT.opara(O_AVER)=='off') THEN
      CALL chem_aver(LOUT, lsite)
      chem_run_aver_ind = .FALSE.
   ENDIF
ENDIF
!
CALL ber_params (ianz, cpara, lpara, werte, MAXW)
IF(ier_num/=0) RETURN
!
iatom = NINT(werte(1))
if(iatom<1 .or. iatom>cr_natoms) then
   ier_num = -137
   ier_typ = ER_APPL
   ier_msg(1) = 'displacement command'
   write(ier_msg(2),'(a, i12)') 'Atom number ', iatom
   return
endif 
CALL indextocell (iatom, icell, isite)
!
IF(opara(O_AVER) == 'yes') CALL chem_aver(LOUT, lsite)   ! User requested fresh aver
!
IF(lsite) THEN
   res_para(1) =  (cr_pos(1,iatom)-chem_ave_pos (1, isite)) -  &
             NINT((cr_pos(1,iatom)-chem_ave_pos (1, isite)))
   res_para(2) =  (cr_pos(2,iatom)-chem_ave_pos (2, isite)) -  &
             NINT((cr_pos(2,iatom)-chem_ave_pos (2, isite)))
   res_para(3) =  (cr_pos(3,iatom)-chem_ave_pos (3, isite)) -  &
             NINT((cr_pos(3,iatom)-chem_ave_pos (3, isite)))
   res_para(0) = 3
ELSE
   kk = -1
   find_type: DO k=1,chem_ave_n(isite)
      IF(cr_iscat(1,iatom)==chem_ave_iscat(isite,k)) THEN
         kk = k
         EXIT find_type
      ENDIF
   ENDDO find_type
   IF(kk>=0) THEN
      res_para(1) =  (cr_pos(1,iatom)-chem_ave_posit(1,isite,kk)) - &
                NINT((cr_pos(1,iatom)-chem_ave_posit(1,isite,kk)))
      res_para(2) =  (cr_pos(2,iatom)-chem_ave_posit(2,isite,kk)) - &
                NINT((cr_pos(2,iatom)-chem_ave_posit(2,isite,kk)))
      res_para(3) =  (cr_pos(3,iatom)-chem_ave_posit(3,isite,kk)) - &
                NINT((cr_pos(3,iatom)-chem_ave_posit(3,isite,kk)))
      res_para(0) = 3
   ENDIF
ENDIF
res_para(0) = 3
IF(opara(O_ECHO)=='yes') THEN
   at_name_d = at_name (cr_iscat (1,iatom) )
   WRITE(output_io, 1000) iatom, at_name_d, res_para(1:3)
1000 FORMAT(i8,1x,a9,3(2x,f12.7))
ENDIF
!
END SUBROUTINE get_displacement
!
!*******************************************************************************
END MODULE chem_aver_mod
