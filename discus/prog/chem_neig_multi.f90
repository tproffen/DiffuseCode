MODULE chem_neig_multi_mod
!-
! Determine neighbors for MMC
!+
CONTAINS
!
!*****7*****************************************************************
!
SUBROUTINE chem_neighbour_multi(jatom, ic, iatom, patom, natom, ncent, maxw, MMC_MAX_CENT)
!+                                                                      
!     Determine neighbours from given atom index 'jatom'.               
!-                                                                      
      USE discus_config_mod 
USE crystal_mod 
USE atom_env_mod 
USE check_bound_mod
USE chem_mod 
USE celltoindex_mod
USE conn_sup_mod
USE do_find_mod
USE metric_mod
USE modify_mod
USE modify_func_mod
!
USE debug_mod 
USE errlist_mod 
USE param_mod 
USE precision_mod
USE prompt_mod 
!                                                                       
IMPLICIT none 
!                                                                       
INTEGER, INTENT(IN)      ::  MAXW    ! Maximum array size
INTEGER, INTENT(IN)      ::  MMC_MAX_CENT    ! Maximum array size
!                                                                       
      INTEGER, INTENT(IN)      ::  jatom   ! Central atom no, get neighbours around jatom
      INTEGER, INTENT(IN)      ::  ic      ! Use correlation definition no. ic
      INTEGER, DIMENSION(0:maxw, MMC_MAX_CENT), INTENT(OUT) :: iatom ! indices of neighbs
      INTEGER, DIMENSION(MMC_MAX_CENT)        , INTENT(OUT) :: natom ! no of neigh
      INTEGER                                  , INTENT(OUT) :: ncent ! no of central atoms
      REAL   , DIMENSION(3, 0:maxw, MMC_MAX_CENT), INTENT(OUT) :: patom ! Coordinates
!
      REAL dist (MMC_MAX_CENT) 
      REAL(KIND=PREC_DP) :: werte (MAX_ATOM_ENV) 
!                                                                       
      REAL u (3), v (3), w (3), uu (3) 
      REAL offset (3)
REAL(KIND=PREC_DP) :: dummy(1) 
      INTEGER jcell (3), icell (3), isite, jsite 
      INTEGER i, j, k, l, ii, iv, katom, ianz 
      INTEGER                    :: is1    ! central atom type
      INTEGER                    :: ino    ! number of connectivity list 
      INTEGER                    :: natoms ! number of atoms in connectivity list 
      INTEGER, DIMENSION(:), ALLOCATABLE :: c_list ! Result of connectivity search
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: c_offs ! Result of connectivity search
      LOGICAL laccept 
      LOGICAL lok 
      LOGICAL ldbg 
!                                                                       
!     LOGICAL atom_allowed 
!     REAL do_bang 
!     REAL do_blen 
!                                                                       
      DO i = 1, MMC_MAX_CENT 
      natom (i) = 0 
      ENDDO 
      dummy = - 1 
      ncent = 1 
!                                                                       
      ldbg = .false. 
!                                                                       
!------ ------------------------------------------------------          
!------ Mode distance                                                   
!------ ------------------------------------------------------          
!                                                                       
      IF (chem_ctyp (ic) .eq.CHEM_DIST) then 
         ncent = 1 
         natom (ncent) = 0 
         DO j = 1, 3 
         u (j) = cr_pos (j, jatom) 
         w (j) = 0.0 
         ENDDO 
         CALL do_find_env (1, dummy, 1, u, chem_rmin (ic), chem_rmax (  &
         ic), chem_quick, chem_period)                                  
         IF (atom_env (0) .ne.0) then 
            DO j = 1, atom_env (0) 
            IF (chem_cang (ic) ) then 
!                                                                       
!------ ------- Check neighbouring angle and symmetry                   
!                                                                       
               DO k = 1, chem_nnei (ic) 
               DO ii = 1, 3 
               v (ii) = atom_pos (ii, j) - cr_pos (ii, jatom) 
               uu (ii) = chem_neig (ii, k, ic) 
               ENDDO 
               IF (abs (do_bang (.true., v, w, uu) )                    &
               .lt.chem_wink_sigma (ic) ) then                          
                  natom (ncent) = natom (ncent) + 1 
                  IF (natom (ncent) .le.maxw) then 
                     iatom (natom (ncent), ncent) = atom_env (j) 
                     patom (1, natom (ncent), ncent) = atom_pos (1, j) 
                     patom (2, natom (ncent), ncent) = atom_pos (2, j) 
                     patom (3, natom (ncent), ncent) = atom_pos (3, j) 
                  ELSE 
                     ier_num = - 23 
                     ier_typ = ER_CHEM 
                     RETURN 
                  ENDIF 
               ENDIF 
               ENDDO 
!                                                                       
!------ ------- Check just distance                                     
!                                                                       
            ELSE 
               natom (ncent) = natom (ncent) + 1 
               IF (natom (ncent) .le.maxw) then 
                  iatom (natom (ncent), ncent) = atom_env (j) 
                  patom (1, natom (ncent), ncent) = atom_pos (1, j) 
                  patom (2, natom (ncent), ncent) = atom_pos (2, j) 
                  patom (3, natom (ncent), ncent) = atom_pos (3, j) 
               ELSE 
                  ier_num = - 23 
                  ier_typ = ER_CHEM 
                  RETURN 
               ENDIF 
            ENDIF 
            ENDDO 
         ENDIF 
         IF (natom (ncent) .gt.0) then 
            iatom (0, ncent) = jatom 
            patom (1, 0, ncent) = cr_pos (1, jatom) 
            patom (2, 0, ncent) = cr_pos (2, jatom) 
            patom (3, 0, ncent) = cr_pos (3, jatom) 
         ENDIF 
!                                                                       
!------ ------------------------------------------------------          
!------ Mode distance-range (New mode for mmc)                          
!------ ------------------------------------------------------          
!                                                                       
      ELSEIF (chem_ctyp (ic) .eq.CHEM_RANGE) then 
         i = 1 
         iv = chem_use_ran (i, ic) 
!                                                                       
!     --Check whether central atom is allowed                           
!                                                                       
         IF (chem_cran_cent (0, iv) .eq. - 1) then 
            laccept = .true. 
         ELSE 
            ianz = chem_cran_cent (0, iv) 
            DO j = 1, chem_cran_cent (0, iv) 
            werte (j) = chem_cran_cent (j, iv) 
!!DBGDISPwrite(*,*) ' neighbours are       ', (iatom(j,ncent),j=1,natom(ncent))
            ENDDO 
            laccept = atom_allowed (jatom, werte, ianz, maxw) 
         ENDIF 
         IF (.not.laccept) return 
!                                                                       
         ncent = 1 
         natom (ncent) = 0 
         DO j = 1, 3 
         u (j) = cr_pos (j, jatom) 
         w (j) = 0.0 
         ENDDO 
         dist (1) = 0.0 
!                                                                       
!     --Check whether neighbour atom is allowed                         
!                                                                       
         IF (chem_cran_neig (0, iv) .eq. - 1) then 
            werte (1) = - 1 
            ianz = 1 
         ELSE 
            ianz = chem_cran_neig (0, iv) 
            DO j = 1, chem_cran_neig (0, iv) 
            werte (j) = chem_cran_neig (j, iv) 
            ENDDO 
         ENDIF 
      CALL do_find_env (1, dummy, 1, u, chem_cran_rmin (iv),  chem_cran_&
     &rmax (iv),  chem_quick, chem_period)                              
         IF (atom_env (0) .ne.0) then 
            DO j = 1, atom_env (0) 
            k = atom_env (j) 
            laccept = atom_allowed (k, werte, ianz, maxw) 
            IF (laccept) then 
!                                                                       
               IF (chem_cran_cang (iv) ) then 
!                                                                       
!------ ------- Check neighbouring angle and symmetry                   
!                                                                       
                  DO k = 1, chem_cran_nuvw (iv) 
                  DO ii = 1, 3 
                  v (ii) = atom_pos (ii, j) - cr_pos (ii, jatom) 
                  uu (ii) = chem_cran_uvw (ii, k, iv) 
                  ENDDO 
                  dummy = do_blen (.true., w, v) 
                  IF (abs (do_bang (.true., v, w, uu) )                 &
                  .lt.chem_cran_wsig (iv) ) then                        
                     natom (ncent) = natom (ncent) + 1 
                     IF (natom (ncent) .le.maxw) then 
                        ii = 1 
                        DO while (dummy(1).gt.dist (ii) .and.ii.lt.natom ( &
                        ncent) )                                        
                        ii = ii + 1 
                        ENDDO 
                        DO l = natom (ncent) - 1, ii, - 1 
                        iatom (l + 1, ncent) = iatom (l, ncent) 
                        patom (1, l + 1, ncent) = patom (1, l, ncent) 
                        patom (2, l + 1, ncent) = patom (2, l, ncent) 
                        patom (3, l + 1, ncent) = patom (3, l, ncent) 
                        dist (l + 1) = dist (l) 
                        ENDDO 
                        iatom (ii, ncent) = atom_env (j) 
                        patom (1, ii, ncent) = atom_pos (1, j) 
                        patom (2, ii, ncent) = atom_pos (2, j) 
                        patom (3, ii, ncent) = atom_pos (3, j) 
                        dist (ii) = dummy (1)
                     ELSE 
                        ier_num = - 23 
                        ier_typ = ER_CHEM 
                        RETURN 
                     ENDIF 
                  ENDIF 
                  ENDDO 
!                                                                       
!------ ------- Check just distance                                     
!                                                                       
               ELSE 
                  natom (ncent) = natom (ncent) + 1 
                  IF (natom (ncent) .le.maxw) then 
                     ii = 1 
                     DO while (dummy(1).gt.dist (ii) .and.ii.lt.natom (    &
                     ncent) )                                           
                     ii = ii + 1 
                     ENDDO 
                     DO l = natom (ncent) - 1, ii, - 1 
                     iatom (l + 1, ncent) = iatom (l, ncent) 
                     patom (1, l + 1, ncent) = patom (1, l, ncent) 
                     patom (2, l + 1, ncent) = patom (2, l, ncent) 
                     patom (3, l + 1, ncent) = patom (3, l, ncent) 
                     dist (l + 1) = dist (l) 
                     ENDDO 
                     iatom (ii, ncent) = atom_env (j) 
                     patom (1, ii, ncent) = atom_pos (1, j) 
                     patom (2, ii, ncent) = atom_pos (2, j) 
                     patom (3, ii, ncent) = atom_pos (3, j) 
                     dist (ii) = dummy(1) 
                  ELSE 
                     ier_num = - 23 
                     ier_typ = ER_CHEM 
                     RETURN 
                  ENDIF 
               ENDIF 
            ENDIF 
            ENDDO 
         ENDIF 
         IF (natom (ncent) .gt.0) then 
            iatom (0, ncent) = jatom 
            patom (1, 0, ncent) = cr_pos (1, jatom) 
            patom (2, 0, ncent) = cr_pos (2, jatom) 
            patom (3, 0, ncent) = cr_pos (3, jatom) 
         ENDIF 
!                                                                       
!     --Only the shortest vectors are to be considered                  
!                                                                       
         IF (chem_cran_short (iv) ) then 
            IF (ldbg) then 
               WRITE ( output_io , * ) ' natom(ncent) ', natom (ncent) 
!      WRITE ( output_io ,  * ) ' nshort       ', chem_cran_nshort (iv) 
            ENDIF 
            DO l = chem_cran_nshort (iv) + 1, natom (ncent) 
            iatom (l, ncent) = 0 
            patom (1, l, ncent) = 0.0 
            patom (2, l, ncent) = 0.0 
            patom (3, l, ncent) = 0.0 
            ENDDO 
            natom (ncent) = min (natom (ncent), chem_cran_nshort (iv) ) 
         ENDIF 
!                                                                       
!------ ------------------------------------------------------          
!------ Mode vector                                                     
!------ ------------------------------------------------------          
!                                                                       
ELSEIF(chem_ctyp(ic) ==  CHEM_VEC) THEN 
   ncent = 1 
   natom (  ncent) = 0 
   iatom (:,ncent) = 0 
   CALL indextocell (jatom, jcell, jsite) 
   DO i = 1, chem_nvec(ic) 
      iv = chem_use_vec(i, ic) 
      IF(jsite == chem_cvec(1, iv) ) THEN     ! Atom jatom is at vector start == central
         iatom(0, ncent) = jatom 
         patom(1, 0, ncent) = cr_pos(1, jatom) 
         patom(2, 0, ncent) = cr_pos(2, jatom) 
         patom(3, 0, ncent) = cr_pos(3, jatom) 
         icell(1) = jcell (1) + chem_cvec (3, iv) 
         icell(2) = jcell (2) + chem_cvec (4, iv) 
         icell(3) = jcell (3) + chem_cvec (5, iv) 
         lok = .false. 
         CALL check_bound(icell, offset, chem_period, lok) 
         IF(lok) THEN 
            isite = chem_cvec(2, iv) 
            CALL celltoindex(icell, isite, katom) 
            natom(ncent) = natom(ncent) + 1 
            iatom(natom (ncent), ncent) = katom 
            patom(1, natom(ncent), ncent) = cr_pos(1, katom) + offset(1)
            patom(2, natom(ncent), ncent) = cr_pos(2, katom) + offset(2)
            patom(3, natom(ncent), ncent) = cr_pos(3, katom) + offset(3)
         ENDIF 
      ELSEIF(jsite == chem_cvec(2, iv) ) THEN ! Atom jatom is at vector end == neighbor
         icell(1) = jcell(1) - chem_cvec(3, iv) 
         icell(2) = jcell(2) - chem_cvec(4, iv) 
         icell(3) = jcell(3) - chem_cvec(5, iv) 
         CALL check_bound(icell, offset, chem_period, lok) 
         IF(lok) THEN 
            natom (ncent) = natom (ncent) + 1 
!
            iatom(natom(ncent), ncent) = jatom 
            patom(1, natom(ncent), ncent) = cr_pos(1, jatom)
            patom(2, natom(ncent), ncent) = cr_pos(2, jatom)
            patom(3, natom(ncent), ncent) = cr_pos(3, jatom)
!
            isite = chem_cvec (1, iv) 
            CALL celltoindex (icell, isite, katom) 
            iatom(0, ncent) = katom
            patom(1, 0, ncent) = cr_pos(1, katom) + offset(1)
            patom(2, 0, ncent) = cr_pos(2, katom) + offset(2)
            patom(3, 0, ncent) = cr_pos(3, katom) + offset(3)
!              iatom (natom (ncent), ncent) = katom 
!              patom(1, natom(ncent), ncent) = cr_pos(1, katom) + offset(1)
!              patom(2, natom(ncent), ncent) = cr_pos(2, katom) + offset(2)
!              patom(3, natom(ncent), ncent) = cr_pos(3, katom) + offset(3)
         ENDIF 
      ENDIF 
   ENDDO 
   IF (natom(ncent) ==  0) THEN 
      ncent = ncent - 1 
      ncent = 0 
   ENDIF 
!                                                                       
!------ ------------------------------------------------------          
!------ Mode connectivity                                                     
!------ ------------------------------------------------------          
!                                                                       
      ELSEIF (chem_ctyp (ic) .eq.CHEM_CON) then 
         ncent = 1 
         natom (ncent) = 0 
!        CALL indextocell (jatom, jcell, jsite) 
         DO i = 1, chem_ncon (ic) 
            iatom (0, ncent)    = jatom             ! Store central atom at entry 0
            patom (1, 0, ncent) = cr_pos (1, jatom) ! Store coord. of central atom
            patom (2, 0, ncent) = cr_pos (2, jatom) 
            patom (3, 0, ncent) = cr_pos (3, jatom) 
            iv = chem_use_con (i, ic)               ! Use connectivity no iv
            IF (cr_iscat(jatom).eq.chem_ccon (1, iv) ) then ! Central has correct type
               is1 = chem_ccon (1, iv)              ! central atom type
               ino = chem_ccon (2, iv)              ! connectivity number
               CALL get_connectivity_list ( jatom, is1, ino, c_list, c_offs, natoms )
               k = natom(ncent)
!write(*,*) ' JATOM ', jatom, natoms, ' :: ',ino, c_list(1:natoms)
               DO j=1,natoms
!if(c_list(j)<1) THEN
!write(*,*) ' BUG   ', j , natoms, ' >> ', j, c_list(1:natoms) , ' << ', c_list(j), ' >>'
!   STOP
!endif
                  iatom(  k+j,ncent) = c_list(j)
                  patom(1,k+j,ncent) = cr_pos(1,c_list(j)) + REAL(c_offs(1,j))
                  patom(2,k+j,ncent) = cr_pos(2,c_list(j)) + REAL(c_offs(2,j))
                  patom(3,k+j,ncent) = cr_pos(3,c_list(j)) + REAL(c_offs(3,j))
               ENDDO 
               natom(ncent) = natom(ncent) + natoms
            ENDIF
         ENDDO 
         IF (natom (ncent) .eq.0) then 
            ncent = ncent - 1 
            ncent = 0 
         ENDIF 
!                                                                       
!------ ------------------------------------------------------          
!------ Mode angular neighbours                                         
!------ ------------------------------------------------------          
!                                                                       
      ELSEIF (chem_ctyp (ic) .eq.CHEM_ANG) then 
         ncent = 0 
         natom (ncent) = 0 
         CALL indextocell (jatom, jcell, jsite) 
         DO i = 1, chem_nwin (ic) 
         ncent = ncent + 1 
         ncent = 1 
         iv = chem_use_win (i, ic) 
         IF (jsite.eq.chem_cwin (1, iv) ) then 
!                                                                       
!     ------Selected atom is center of the angle                        
!                                                                       
            iatom (0, ncent) = jatom 
            patom (1, 0, ncent) = cr_pos (1, jatom) 
            patom (2, 0, ncent) = cr_pos (2, jatom) 
            patom (3, 0, ncent) = cr_pos (3, jatom) 
            icell (1) = jcell (1) + chem_cwin (3, iv) 
            icell (2) = jcell (2) + chem_cwin (4, iv) 
            icell (3) = jcell (3) + chem_cwin (5, iv) 
            CALL check_bound (icell, offset, chem_period, lok) 
            IF (lok) then 
               isite = chem_cwin (2, iv) 
               CALL celltoindex (icell, isite, katom) 
               natom (ncent) = natom (ncent) + 1 
               iatom (natom (ncent), ncent) = katom 
               patom (1, natom (ncent), ncent) = cr_pos (1, katom)      &
               + offset (1)                                             
               patom (2, natom (ncent), ncent) = cr_pos (2, katom)      &
               + offset (2)                                             
               patom (3, natom (ncent), ncent) = cr_pos (3, katom)      &
               + offset (3)                                             
            ENDIF 
            icell (1) = jcell (1) + chem_cwin (7, iv) 
            icell (2) = jcell (2) + chem_cwin (8, iv) 
            icell (3) = jcell (3) + chem_cwin (9, iv) 
            CALL check_bound (icell, offset, chem_period, lok) 
            IF (lok) then 
               isite = chem_cwin (6, iv) 
               CALL celltoindex (icell, isite, katom) 
               natom (ncent) = natom (ncent) + 1 
               iatom (natom (ncent), ncent) = katom 
               patom (1, natom (ncent), ncent) = cr_pos (1, katom)      &
               + offset (1)                                             
               patom (2, natom (ncent), ncent) = cr_pos (2, katom)      &
               + offset (2)                                             
               patom (3, natom (ncent), ncent) = cr_pos (3, katom)      &
               + offset (3)                                             
            ENDIF 
         ELSEIF (jsite.eq.chem_cwin (2, iv) ) then 
!                                                                       
!     ------Selected atom is first neighbour of the central atom        
!                                                                       
            icell (1) = jcell (1) - chem_cwin (3, iv) 
            icell (2) = jcell (2) - chem_cwin (4, iv) 
            icell (3) = jcell (3) - chem_cwin (5, iv) 
            CALL check_bound (icell, offset, chem_period, lok) 
            IF (lok) then 
!     --------The central atom is within the boundary condition         
               isite = chem_cwin (1, iv) 
               CALL celltoindex (icell, isite, katom) 
               iatom (0, ncent) = katom 
               patom (1, 0, ncent) = cr_pos (1, katom) + offset (1) 
               patom (2, 0, ncent) = cr_pos (2, katom) + offset (2) 
               patom (3, 0, ncent) = cr_pos (3, katom) + offset (3) 
               natom (ncent) = natom (ncent) + 1 
               iatom (natom (ncent), ncent) = jatom 
               patom (1, natom (ncent), ncent) = cr_pos (1, jatom) 
               patom (2, natom (ncent), ncent) = cr_pos (2, jatom) 
               patom (3, natom (ncent), ncent) = cr_pos (3, jatom) 
               icell (1) = jcell (1) - chem_cwin (3, iv) + chem_cwin (7,&
               iv)                                                      
               icell (2) = jcell (2) - chem_cwin (4, iv) + chem_cwin (8,&
               iv)                                                      
               icell (3) = jcell (3) - chem_cwin (5, iv) + chem_cwin (9,&
               iv)                                                      
               CALL check_bound (icell, offset, chem_period, lok) 
               IF (lok) then 
                  isite = chem_cwin (6, iv) 
                  CALL celltoindex (icell, isite, katom) 
                  natom (ncent) = natom (ncent) + 1 
                  iatom (natom (ncent), ncent) = katom 
                  patom (1, natom (ncent), ncent) = cr_pos (1, katom)   &
                  + offset (1)                                          
                  patom (2, natom (ncent), ncent) = cr_pos (2, katom)   &
                  + offset (2)                                          
                  patom (3, natom (ncent), ncent) = cr_pos (3, katom)   &
                  + offset (3)                                          
               ENDIF 
            ELSE 
!     --------The central atom is outside the (periodic) boundary       
               natom (ncent) = 0 
            ENDIF 
         ELSEIF (jsite.eq.chem_cwin (6, iv) ) then 
!                                                                       
!     ------Selected atom is second neighbour of the central atom       
!                                                                       
            icell (1) = jcell (1) - chem_cwin (7, iv) 
            icell (2) = jcell (2) - chem_cwin (8, iv) 
            icell (3) = jcell (3) - chem_cwin (9, iv) 
            CALL check_bound (icell, offset, chem_period, lok) 
            IF (lok) then 
               isite = chem_cwin (1, iv) 
               CALL celltoindex (icell, isite, katom) 
               iatom (0, ncent) = katom 
               patom (1, 0, ncent) = cr_pos (1, katom) + offset (1) 
               patom (2, 0, ncent) = cr_pos (2, katom) + offset (2) 
               patom (3, 0, ncent) = cr_pos (3, katom) + offset (3) 
               icell (1) = jcell (1) - chem_cwin (7, iv) + chem_cwin (3,&
               iv)                                                      
               icell (2) = jcell (2) - chem_cwin (8, iv) + chem_cwin (4,&
               iv)                                                      
               icell (3) = jcell (3) - chem_cwin (9, iv) + chem_cwin (5,&
               iv)                                                      
               CALL check_bound (icell, offset, chem_period, lok) 
               IF (lok) then 
                  isite = chem_cwin (2, iv) 
                  CALL celltoindex (icell, isite, katom) 
                  natom (ncent) = natom (ncent) + 1 
                  iatom (natom (ncent), ncent) = katom 
                  patom (1, natom (ncent), ncent) = cr_pos (1, katom)   &
                  + offset (1)                                          
                  patom (2, natom (ncent), ncent) = cr_pos (2, katom)   &
                  + offset (2)                                          
                  patom (3, natom (ncent), ncent) = cr_pos (3, katom)   &
                  + offset (3)                                          
               ENDIF 
               natom (ncent) = natom (ncent) + 1 
               iatom (natom (ncent), ncent) = jatom 
               patom (1, natom (ncent), ncent) = cr_pos (1, jatom) 
               patom (2, natom (ncent), ncent) = cr_pos (2, jatom) 
               patom (3, natom (ncent), ncent) = cr_pos (3, jatom) 
            ELSE 
!     --------The central atom is outside the (periodic) boundary       
               natom (ncent) = 0 
            ENDIF 
         ENDIF 
         ENDDO 
         IF (natom (ncent) .eq.0) then 
            ncent = ncent - 1 
            ncent = 0 
         ENDIF 
!                                                                       
!------ ------------------------------------------------------          
!------ Mode environment                                                
!------ ------------------------------------------------------          
!                                                                       
      ELSEIF (chem_ctyp (ic) .eq.CHEM_ENVIR) then 
         ncent = 0 
         natom (1) = 0 
         DO j = 1, 3 
         u (j) = cr_pos (j, jatom) 
         w (j) = 0.0 
         ENDDO 
!                                                                       
!--------Loop over all environment definitions used by current          
!        neighbor definition                                            
!                                                                       
         DO i = 1, chem_nenv (ic) 
!                                                                       
!----------The selected atom is the central atom of environment         
!          definition                                                   
!                                                                       
         IF (cr_iscat (jatom) .eq.chem_cenv (0, chem_use_env (i, ic) )  &
         .or. - 1.eq.chem_cenv (0, chem_use_env (i, ic) ) ) then        
            DO j = 1, chem_env_neig (chem_use_env (i, ic) ) 
            werte (j) = chem_cenv (j, chem_use_env (i, ic) ) 
            ENDDO 
            ianz = chem_env_neig (chem_use_env (i, ic) ) 
            CALL do_find_env (ianz, werte, maxw, u, chem_rmin_env (ic), &
            chem_rmax_env (ic), chem_quick, chem_period)                
            IF (atom_env (0) .ne.0) then 
!                                                                       
!     --------Neighbours were found, add to list of atoms               
!                                                                       
               ncent = ncent + 1 
               ncent = 1 
               DO j = 1, atom_env (0) 
               IF (natom (ncent) .lt.maxw) then 
                  natom (ncent) = natom (ncent) + 1 
                  iatom (natom (ncent), ncent) = atom_env (j) 
                  patom (1, natom (ncent), ncent) = atom_pos (1, j) 
                  patom (2, natom (ncent), ncent) = atom_pos (2, j) 
                  patom (3, natom (ncent), ncent) = atom_pos (3, j) 
               ELSE 
                  natom (ncent) = 0 
                  ncent = ncent - 1 
                  ier_num = - 23 
                  ier_typ = ER_CHEM 
                  RETURN 
               ENDIF 
               ENDDO 
               iatom (0, ncent) = jatom 
               patom (1, 0, ncent) = cr_pos (1, jatom) 
               patom (2, 0, ncent) = cr_pos (2, jatom) 
               patom (3, 0, ncent) = cr_pos (3, jatom) 
            ENDIF 
         ELSE 
!                                                                       
!     ----The selected atom is a neighbour atom of environment          
!           definition                                                  
!                                                                       
!c           lis_neig = .false.                                         
!c           do j=1,chem_env_neig(chem_use_env(i,ic))                   
!c             if(cr_iscat(jatom).eq.                                   
!C     &                      chem_cenv(j,chem_use_env(i,ic))) then     
!c               lis_neig = .true.                                      
!c             endif                                                    
!c           ENDDO                                                      
!c           if(lis_neig) then                                          
!c             werte(1) = chem_cenv(0,chem_use_env(i,ic))               
!c             ianz     = 1                                             
!c             call do_find_env(ianz,werte,maxw,u,chem_rmin_env(ic),    
!C     &                                  chem_rmax_env(ic),            
!C     &                                  chem_quick,chem_period)       
!                                                                       
!-----      ------- Store the true central atoms in patom               
!                                                                       
!c             do j=1,atom_env(0)                                       
!c               iatom(0,ncent+j)   = atom_env(j)                       
!c               patom(1,0,ncent+j) = atom_pos(1,j)                     
!c               patom(2,0,ncent+j) = atom_pos(2,j)                     
!c               patom(3,0,ncent+j) = atom_pos(3,j)                     
!c             ENDDO                                                    
!c             nnew = atom_env(0)                                       
!c             do k=1,nnew                                              
!c               do j=1,chem_env_neig(chem_use_env(i,ic))               
!c                 werte(j) = chem_cenv(j,chem_use_env(i,ic))           
!c               ENDDO                                                  
!c               do j=1,3                                               
!c                 u(j) = patom(j,0,ncent+1)                            
!c               ENDDO                                                  
!c               ianz = chem_env_neig(chem_use_env(i,ic))               
!c               call do_find_env(ianz,werte,maxw,u,chem_rmin_env(ic),  
!C     &                                    chem_rmax_env(ic),          
!C     &                                    chem_quick,chem_period)     
!c               if (atom_env(0).ne.0) then                             
!c                 ncent = ncent + 1                                    
!c                 ncent =         1                                    
!c                 do j=1,atom_env(0)                                   
!c                   if (natom(ncent).lt.maxw) then                     
!c                     natom(ncent)          = natom(ncent) + 1         
!c                     iatom(natom(ncent),ncent)   = atom_env(j)        
!c                     patom(1,natom(ncent),ncent) = atom_pos(1,j)      
!c                     patom(2,natom(ncent),ncent) = atom_pos(2,j)      
!c                     patom(3,natom(ncent),ncent) = atom_pos(3,j)      
!c                   ELSE                                               
!DBG_SPHERES                                                            
!C  write(*,*) ' NEIGHBOUR, NUMBER OF NEIGHBOURS ',natom(ncent),        
!C  atom_env(0),nnew                                                    
!c                     natom(ncent) = 0                                 
!c                     ncent = ncent - 1                                
!c                     ncent =         0                                
!c                     ier_num = -23                                    
!c                     ier_typ = ER_CHEM                                
!c                     return                                           
!c                   endif                                              
!c                 ENDDO                                                
!c               endif                                                  
!c             ENDDO                                                    
!c           endif                                                      
         ENDIF 
         ENDDO 
      ENDIF 
      ldbg = .false. 
!                                                                       
END SUBROUTINE chem_neighbour_multi           
!
!*****7*****************************************************************
!
END MODULE chem_neig_multi_mod
