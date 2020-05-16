MODULE do_find_mod
!
CONTAINS
!*****7*****************************************************************
      SUBROUTINE do_find (line, laenge) 
!-                                                                      
!     Finds the environment around an atom                              
!+                                                                      
      USE discus_config_mod 
      USE charact_mod 
      USE charact_mod 
      USE crystal_mod 
      USE get_iscat_mod
      USE chem_mod 
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
USE precision_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 200) 
      INTEGER mmaxw 
      PARAMETER (mmaxw = 5) 
!                                                                       
      CHARACTER ( * ) line 
      CHARACTER(LEN=MAX(PREC_STRING,LEN(line))) cpara (maxw) 
      CHARACTER(LEN=MAX(PREC_STRING,LEN(line))) ccpara (mmaxw) 
      INTEGER lpara (maxw) 
      INTEGER llpara (maxw) 
      INTEGER i, ii, ianz, iianz, laenge 
      LOGICAL lnew, fq, fp (3) 
      REAL rmin 
      REAL radius 
      REAL(KIND=PREC_DP) :: werte (maxw) 
      REAL(KIND=PREC_DP) :: wwerte (maxw) 
      REAL x (3) 
!                                                                       
      PARAMETER (lnew = .false.) 
!                                                                       
      LOGICAL str_comp 
!                                                                       
      fp (1) = chem_period (1) 
      fp (2) = chem_period (2) 
      fp (3) = chem_period (3) 
      fq = chem_quick 
!                                                                       
      CALL get_params (line, ianz, cpara, lpara, maxw, laenge) 
      IF (ier_num.eq.0) then 
         IF (str_comp (cpara (1) , 'env', 1, lpara (1) , 3) ) then 
!                                                                       
!     ----Find environment                                              
!                                                                       
            IF (ianz.ge.7) then 
!                                                                       
!     ------copy last five parameters for evaluation                    
!                                                                       
               DO i = ianz - 4, ianz 
               ii = i - (ianz - 5) 
               ccpara (ii) = cpara (i) 
               llpara (ii) = lpara (i) 
               ENDDO 
               iianz = 5 
               CALL ber_params (iianz, ccpara, llpara, wwerte, mmaxw) 
               x (1) = wwerte (1) 
               x (2) = wwerte (2) 
               x (3) = wwerte (3) 
               rmin = wwerte (4) 
               radius = wwerte (5) 
               IF (ier_num.eq.0) then 
!                                                                       
!     -------- shift remaining parameters one left                      
!                                                                       
                  DO i = 2, ianz - 5 
                  cpara (i - 1) = cpara (i) 
                  lpara (i - 1) = lpara (i) 
                  ENDDO 
                  ianz = ianz - 6 
!                                                                       
!     --------Get scattering curves                                     
!                                                                       
                  CALL get_iscat (ianz, cpara, lpara, werte, maxw, lnew) 
                  IF (ier_num.eq.0) then 
                     CALL do_find_env (ianz, werte, maxw, x, rmin,      &
                     radius, fq, fp)                                    
                  ENDIF 
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
         ELSEIF (str_comp (cpara (1) , 'menv', 1, lpara (1) , 4) ) then 
!                                                                       
!     ----Find molecular environment                                    
!                                                                       
            IF (ianz.ge.7) then 
!                                                                       
!     ------copy last three parameters for evaluation                   
!                                                                       
               DO i = ianz - 4, ianz 
               ii = i - (ianz - 5) 
               ccpara (ii) = cpara (i) 
               llpara (ii) = lpara (i) 
               ENDDO 
               iianz = 5 
               CALL ber_params (iianz, ccpara, llpara, wwerte, mmaxw) 
               x (1) = wwerte (1) 
               x (2) = wwerte (2) 
               x (3) = wwerte (3) 
               rmin = wwerte (4) 
               radius = wwerte (5) 
               IF (ier_num.eq.0) then 
!                                                                       
!     -------- shift remaining parameters one left                      
!                                                                       
                  DO i = 2, ianz - 5 
                  cpara (i - 1) = cpara (i) 
                  lpara (i - 1) = lpara (i) 
                  ENDDO 
                  ianz = ianz - 6 
!                                                                       
!     --------Get allowed molecule types                                
!                                                                       
                  IF (str_comp (cpara (1) , 'all', 1, lpara (1) , 3) )  &
                  then                                                  
                     ianz = - 1 
                  ELSE 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  ENDIF 
                  IF (ier_num.eq.0) then 
                     CALL do_find_mol (ianz, werte, maxw, x, rmin, radius)
                  ENDIF
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
      ELSE 
         ier_num = - 10 
         ier_typ = ER_APPL 
      ENDIF 
!                                                                       
      END SUBROUTINE do_find                        
!*****7*****************************************************************
      SUBROUTINE do_find_env (ianz, werte, maxw, x, rmin, rmax, fq, fp) 
!                                                                       
!     This routine finds all atoms around x with a minimal              
!     distance of rmin and a maximum distance of rmax.                  
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE atom_env_mod 
      USE check_bound_mod
      USE check_blen_mod
      USE celltoindex_mod
      USE modify_func_mod
      USE param_mod 
      USE errlist_mod 
      USE sorting_mod 
USE precision_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER ianz, maxw 
      REAL(KIND=PREC_DP) :: werte (maxw) 
      REAL x (3) 
      REAL rmin, rmax 
      LOGICAL fq, fp (3) 
!                                                                       
      INTEGER i, j, k, ii 
      INTEGER :: ix, iy, iz
      INTEGER :: ix1, ix2, iy1,iy2, iz1,iz2
      INTEGER istart (3), iend (3), iii (3), cell (3), iatom 
      REAL offset (3), nooffset (3) 
      LOGICAL ltype 
      INTEGER, DIMENSION(  :), ALLOCATABLE :: tmp_ind
      INTEGER, DIMENSION(  :), ALLOCATABLE :: tmp_env
      REAL   , DIMENSION(:,:), ALLOCATABLE :: tmp_pos
      REAL   , DIMENSION(  :), ALLOCATABLE :: tmp_dis
!                                                                       
      DATA nooffset / 0.0, 0.0, 0.0 / 
!                                                                       
      atom_env (0) = 0 
      res_para (0) = 0 
      IF (fq) then 
!                                                                       
!------ --quick version only looks at neighbouring unit cells           
!                                                                       
         DO i = 1, 3 
         iii (i) = int (x (i) - cr_dim0 (i, 1) ) + 1 
         istart (i) = iii (i) - 1 - int (rmax / cr_a0 (i) )  -1
         iend (i) = iii (i) + 1 + int (rmax / cr_a0 (i) )  + 1
         ENDDO 
!                                                                       
         DO k = istart (3), iend (3) 
         DO j = istart (2), iend (2) 
         DO i = istart (1), iend (1) 
         cell (1) = i 
         cell (2) = j 
         cell (3) = k 
!                                                                       
         CALL check_bound (cell, offset, fp, ltype) 
         IF (ltype) then 
            DO ii = 1, cr_ncatoms 
            CALL celltoindex (cell, ii, iatom) 
              ltype = atom_allowed (iatom, werte, ianz, maxw) .and. &
                      check_select_status(iatom, .true., cr_prop (iatom),  cr_sel_prop)                    
            IF (ltype) then 
               CALL check_blen (x, iatom, rmin, rmax, offset) 
               IF (ier_num.ne.0) return 
            ENDIF 
            ENDDO 
         ENDIF 
!                                                                       
         ENDDO 
         ENDDO 
         ENDDO 
!                                                                       
!     --exact: loop over all atoms in crystal that have been added      
!                                                                       
         DO i = cr_ncatoms * cr_icc (1) * cr_icc (2) * cr_icc (3)       &
         + 1, cr_natoms                                                 
!                                                                       
         IF (fp (1) .or.fp (2) .or.fp (3) ) then 
            ier_num = - 16 
            ier_typ = ER_CHEM 
            ier_msg(1) = 'Number of atoms in crystal is larger than'
            ier_msg(2) = 'atoms_per_unit_cell*number_of_unit_cells'
            ier_msg(3) = 'Use set crystal, 1, 1, 1,n[1] to correct'
            RETURN 
         ENDIF 
!                                                                       
         ltype = atom_allowed (i, werte, ianz, maxw)                &
         .and.check_select_status (i, .true., cr_prop (i), cr_sel_prop)    
         IF (ltype) then 
            CALL check_blen (x, i, rmin, rmax, nooffset) 
            IF (ier_num.ne.0) return 
         ENDIF 
         ENDDO 
      ELSE 
!                                                                       
!     --exact: loop over all atoms in crystal                           
!                                                                       
!!!      IF (fp (1) .or.fp (2) .or.fp (3) ) then 
!!!         ier_num = - 16 
!!!         ier_typ = ER_CHEM 
!!!         ier_msg(1) = 'Chem/exact mode is incompatible with '
!!!!        ier_msg(2) = 'periodic boundary conditions'
!!!         ier_msg(3) = ' '
!!!         RETURN 
!!!      ENDIF 
!                                                                       
!        IF(fp(1)) THEN
!           ix1 = -1
!           ix2 = 1
!        ELSE
!           ix1 = 0
!           ix2 = 0
!        ENDIF
!        IF(fp(2)) THEN
!           iy1 = -1
!           iy2 = 1
!        ELSE
!           iy1 = 0
!           iy2 = 0
!        ENDIF
!        IF(fp(3)) THEN
!           iz1 = -1
!           iz2 = 1
!        ELSE
!           iz1 = 0
!           iz2 = 0
!        ENDIF
         ix1 = 0
         ix2 = 0
         iy1 = 0
         iy2 = 0
         iz1 = 0
         iz2 = 0
         IF(fp(1)) THEN
            IF(x(1)       -cr_dim(1,1) < rmax*1.5) ix1 = -1
            IF(cr_dim(1,2)-x(1)        < rmax*1.5) ix2 =  1
         ENDIF
         IF(fp(2)) THEN
            IF(x(2)       -cr_dim(2,1) < rmax*1.5) iy1 = -1
            IF(cr_dim(2,2)-x(2)        < rmax*1.5) iy2 =  1
         ENDIF
         IF(fp(3)) THEN
            IF(x(3)       -cr_dim(3,1) < rmax*1.5) iz1 = -1
            IF(cr_dim(3,2)-x(3)        < rmax*1.5) iz2 =  1
         ENDIF
   DO i = 1, cr_natoms 
      ltype = atom_allowed (i, werte, ianz, maxw)                    &
             .and.check_select_status (i, .true., cr_prop (i), cr_sel_prop)    
      IF (ltype) then 
         DO ix = ix1, ix2, 1
            offset(1) = ix*cr_icc(1)
            DO iy = iy1, iy2, 1
               offset(2) = iy*cr_icc(2) 
               DO iz = iz1, iz2, 1
                     offset(3) = iz*cr_icc(3)
                  CALL check_blen (x, i, rmin, rmax, offset) 
               ENDDO
            ENDDO
         ENDDO
         IF (ier_num /= 0) RETURN
      ENDIF 
   ENDDO 
ENDIF 
!
!     Sort neighbors according to distance
!
      ALLOCATE(tmp_ind(  1:atom_env(0)))
      ALLOCATE(tmp_env(  1:atom_env(0)))
      ALLOCATE(tmp_pos(3,0:atom_env(0)))
      ALLOCATE(tmp_dis(  1:atom_env(0)))
      tmp_env    = atom_env(1:atom_env(0))
      tmp_ind    = 0
      tmp_pos    = atom_pos
      tmp_dis    = atom_dis(1:atom_env(0))
      CALL indexx(atom_env(0),tmp_dis,tmp_ind)
      DO i=1,atom_env(0)
         atom_env(i)   = tmp_env(  tmp_ind(i))
         atom_pos(:,i) = tmp_pos(:,tmp_ind(i))
         atom_dis(i)   = tmp_dis(  tmp_ind(i))
      ENDDO
      DEALLOCATE(tmp_ind)
      DEALLOCATE(tmp_env)
      DEALLOCATE(tmp_pos)
      DEALLOCATE(tmp_dis)
!                                                                       
      END SUBROUTINE do_find_env                    
!
!*****7*****************************************************************
!
      SUBROUTINE do_find_mol (ianz, werte, maxw, x, rmin, rmax) 
!                                                                       
!     This routine finds all molecules around x with a minimal          
!     distance of rmin and a maximum distance of rmax.                  
!-                                                                      
!                                                                       
      USE discus_config_mod 
      USE check_blen_mod
      USE crystal_mod 
      USE molecule_mod 
      USE mole_env_mod 
      USE errlist_mod 
      USE param_mod 
USE precision_mod
      IMPLICIT none 
       
!                                                                       
      INTEGER ianz, maxw 
      REAL(KIND=PREC_DP) :: werte (maxw) 
      REAL x (3) 
      REAL rmin, rmax 
!                                                                       
      INTEGER i, j 
!                                                                       
      REAL nooffset (3) 
      LOGICAL ltype 
!                                                                       
      DATA nooffset / 0.0, 0.0, 0.0 / 
!                                                                       
      mole_env (0) = 0 
      res_para (0) = 0 
!                                                                       
!     Exact loop over all molecules, no periodic boundary conditions    
!                                                                       
      DO i = 1, mole_num_mole 
      IF (ianz.eq. - 1) then 
         ltype = .true. 
      ELSE 
         ltype = .false. 
         DO j = 1, ianz 
         ltype = i.eq.nint (werte (j) ) 
         ENDDO 
      ENDIF 
      IF (ltype) then 
         CALL check_blen_mol (x, i, rmin, rmax, nooffset) 
      ENDIF 
      IF (ier_num.ne.0) return 
      ENDDO 
!DBG                                                                    
!DBG      do i=1,mole_env(0)                                            
!DBG        write (output_io,*) 'molecule ',mole_env(i)                 
!DBG      ENDDO                                                         
      END SUBROUTINE do_find_mol                    
!
!*****7*****************************************************************
!
END MODULE do_find_mod
