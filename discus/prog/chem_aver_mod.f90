MODULE chem_aver_mod
!
CONTAINS
!*****7*****************************************************************
SUBROUTINE chem_aver (lout, lsite) 
!+                                                                      
!     Calculate average structure and standard deviation                
!-                                                                      
USE config_mod 
USE allocate_appl_mod 
USE crystal_mod 
USE atom_name
USE chem_mod 
USE errlist_mod 
USE param_mod 
USE prompt_mod 
IMPLICIT none 
!                                                                       
LOGICAL, INTENT(IN) :: lout    ! Print output if true
LOGICAL, INTENT(IN) :: lsite   ! Treat different atoms on each site as one 
!                                                                       
REAL, DIMENSION(3) ::  p , ez 
REAL,    DIMENSION(:,:,:), ALLOCATABLE :: chem_ave_posit
REAL,    DIMENSION(:,:,:), ALLOCATABLE :: chem_ave_sigma
!
INTEGER            :: i, j, k, ii, jj, kk, ia, is, nvalues
LOGICAL            :: flag
!                                                                       
CHARACTER(LEN=9)   :: at_name_i 
!
INTEGER            :: n_atom_cell  ! Dummy for allocation
INTEGER            :: n_max_atom   ! Dummy for allocation
!
is = 1
!
IF ( CHEM_MAXAT_CELL   < MAXAT_CELL .or. &
     CHEM_MAX_AVE_ATOM < MAX(cr_ncatoms, MAXSCAT) .or. &
     CHEM_MAXAT_CELL   <= cr_ncatoms                    ) THEN
   n_atom_cell = MAX(CHEM_MAXAT_CELL, MAXAT_CELL, cr_ncatoms)
   n_max_atom  = MAX(CHEM_MAX_AVE_ATOM, cr_ncatoms, MAXSCAT) + 1
   call alloc_chem_aver ( n_atom_cell, n_max_atom)
ENDIF
IF(.not. lsite) THEN
   n_atom_cell = MAX(CHEM_MAXAT_CELL, MAXAT_CELL)
   ALLOCATE(chem_ave_posit(3,n_atom_cell, 10))
   ALLOCATE(chem_ave_sigma(3,n_atom_cell, 10))
   chem_ave_posit = 0.0
   chem_ave_sigma = 0.0
ENDIF
!                                                                       
!------ reset counters                                                  
!                                                                       
IF (cr_ncatoms.gt.CHEM_MAXAT_CELL) then 
   ier_num = - 103 
   ier_typ = ER_APPL 
   ier_msg (1) = 'Adjust the value of the variable' 
   ier_msg (2)  = 'MAXAT_CELL in config_mod.f90 and  ' 
   ier_msg (3)  = 'compile the program             ' 
   RETURN 
ENDIF 
chem_ave_n    = 0    ! (i)   , i=1,cr_ncatoms
chem_ave_bese = 0.0  ! (i, k), i=1,cr_ncatoms, k = 1,chem_max_ave_atom
chem_ave_pos  = 0.0  ! (j, i), i=1,cr_ncatoms, j = 1,3
chem_ave_sig  = 0.0  ! (j, i), i=1,cr_ncatoms, j = 1,3
!     DO i = 1, cr_ncatoms 
!     chem_ave_n (i) = 0 
!     DO k = 1, chem_max_atom 
!     chem_ave_bese (i, k) = 0.0 
!     ENDDO 
!     DO j = 1, 3 
!     chem_ave_pos (j, i) = 0.0 
!     chem_ave_sig (j, i) = 0.0 
!     ENDDO 
!     ENDDO 
!                                                                       
!------ loop over all unit cells ans atoms within unit cell             
!                                                                       
loopk: DO k = 1, cr_icc (3) 
   loopj: DO j = 1, cr_icc (2) 
      loopi: DO i = 1, cr_icc (1) 
         ez (1) = cr_dim0 (1, 1) + float (i - 1) 
         ez (2) = cr_dim0 (2, 1) + float (j - 1) 
         ez (3) = cr_dim0 (3, 1) + float (k - 1) 
         loopii: DO ii = 1, cr_ncatoms 
            ia = ( (k - 1) * cr_icc (1) * cr_icc (2) + &
                   (j - 1) * cr_icc (1) + (i - 1) ) * cr_ncatoms + ii                                     
            DO jj = 1, 3 
               p (jj) = cr_pos (jj, ia) - ez (jj) 
               chem_ave_pos (jj, ii) = chem_ave_pos (jj, ii) + p (jj) 
               chem_ave_sig (jj, ii) = chem_ave_sig (jj, ii) + p (jj)**2 
            ENDDO 
!                                                                       
!------ --- Calculate occupancies ..                                    
!                                                                       
            occup: IF (chem_ave_n (ii) .eq.0) then 
               chem_ave_n (ii) = 1 
               chem_ave_iscat (ii, chem_ave_n (ii) ) = cr_iscat (ia) 
               is = 1 
            ELSE  occup
               flag = .true. 
               DO kk = 1, chem_ave_n (ii) 
                  IF (cr_iscat (ia) .eq.chem_ave_iscat (ii, kk) ) then 
                     is = kk 
                     flag = .false. 
                  ENDIF 
               ENDDO 
               IF (flag) then 
                  chem_ave_n (ii) = chem_ave_n (ii) + 1 
                  is = chem_ave_n (ii) 
                  IF (chem_ave_n (ii) .gt.chem_max_atom) then 
                     ier_typ = ER_CHEM 
                     ier_num = - 5 
                     RETURN 
                  ENDIF 
                  chem_ave_iscat (ii, chem_ave_n (ii) ) = cr_iscat (ia) 
               ENDIF 
            ENDIF occup
            IF(.not. lsite) THEN   ! Accumulate individual positions for different atoms
               DO jj = 1, 3 
                  chem_ave_posit(jj,ii,is) = chem_ave_posit(jj,ii,is) + p(jj)
                  chem_ave_sigma(jj,ii,is) = chem_ave_sigma(jj,ii,is) + p(jj)**2
               ENDDO 
            ENDIF
            chem_ave_bese (ii, is) = chem_ave_bese (ii, is) + 1 
         ENDDO loopii
      ENDDO  loopi
   ENDDO  loopj
ENDDO  loopk
!                                                                       
!------ output of average and sigma                                     
!                                                                       
IF (lout) write (output_io, 1000) 
ia = cr_icc (1) * cr_icc (2) * cr_icc (3) 
IF(lsite) THEN           ! one pos for all types on a single site
   DO i = 1, cr_ncatoms 
      DO j = 1, 3 
         chem_ave_pos (j, i) = chem_ave_pos (j, i) / float (ia) 
         chem_ave_sig (j, i) = chem_ave_sig (j, i) / float (ia) -          &
         chem_ave_pos (j, i) **2                                           
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
         ENDDO 
      ENDIF 
   ENDDO 
!                                                                       
!------ store results in res_para                                       
!                                                                       
   IF ( (6 * cr_ncatoms) .gt.maxpar_res) then 
      ier_typ = ER_CHEM 
      ier_num = - 2 
   ELSE 
      res_para (0) = 6 * cr_ncatoms 
      DO i = 1, cr_ncatoms 
         DO j = 1, 3 
            res_para ( (i - 1) * 6 + j) = chem_ave_pos (j, i) 
         ENDDO 
         DO j = 1, 3 
            res_para ( (i - 1) * 6 + j + 3) = chem_ave_sig (j, i) 
         ENDDO 
      ENDDO 
   ENDIF 
ELSE
   nvalues = 0
   DO i = 1, cr_ncatoms 
      DO k = 1, chem_ave_n (i) 
         DO j = 1, 3 
            chem_ave_posit (j, i, k) = chem_ave_posit (j, i, k) / chem_ave_bese(i,k) 
            chem_ave_sigma (j, i, k) = chem_ave_sigma (j, i, k) / chem_ave_bese(i,k) - &
            chem_ave_posit (j, i, k) **2                                           
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
   IF ( (6 * nvalues) .gt.maxpar_res) then 
      ier_typ = ER_CHEM 
      ier_num = - 2 
   ELSE 
      res_para (0) = 6 * nvalues 
      ii = 0
      DO i = 1, cr_ncatoms 
         DO k = 1, chem_ave_n (i) 
         ii = ii + 1
         DO j = 1, 3 
            res_para ( (ii - 1) * 6 + j)     = chem_ave_posit (j, i, k) 
         ENDDO 
         DO j = 1, 3 
            res_para ( (ii - 1) * 6 + j + 3) = chem_ave_sigma (j, i, k) 
         ENDDO 
         ENDDO 
      ENDDO 
   ENDIF 
ENDIF
IF(.not. lsite) THEN
   DEALLOCATE(chem_ave_posit)
   DEALLOCATE(chem_ave_sigma)
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
      USE config_mod 
      USE crystal_mod 
      USE atom_name
      USE chem_mod 
!                                                                       
      USE errlist_mod 
      USE param_mod 
      USE prompt_mod 
      IMPLICIT none 
       
!                                                                       
      REAL proz 
      INTEGER natom (0:MAXSCAT) 
      INTEGER i 
      LOGICAL lout 
!                                                                       
      CHARACTER(9) at_name_i 
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
      DO i = 0, cr_nscat 
      natom (i) = 0 
      ENDDO 
!                                                                       
!------ get size of model crystal, rel. amount of elements              
!                                                                       
      DO i = 1, cr_natoms 
      natom (cr_iscat (i) ) = natom (cr_iscat (i) ) + 1 
      ENDDO 
!                                                                       
!------ write output                                                    
!                                                                       
      IF (lout) write (output_io, 1000) (cr_icc (i), i = 1, 3) 
      IF (lout) write (output_io, 1100) cr_natoms, cr_ncatoms, cr_nscat 
      res_para (0) = float (cr_nscat) + 1 
      DO i = 0, cr_nscat 
      proz = float (natom (i) ) / cr_natoms 
      IF (lout) then 
         at_name_i = at_name (i) 
         WRITE (output_io, 1200) at_name_i, proz, natom (i) 
      ENDIF 
      IF (i.le.maxpar_res) then 
         res_para (i + 1) = proz 
      ELSE 
         ier_typ = ER_CHEM 
         ier_num = - 2 
      ENDIF 
      ENDDO 
!                                                                       
 1000 FORMAT     (' Size of the crystal (unit cells) : ',2(I4,' x '),I4) 
 1100 FORMAT     (' Total number of atoms            : ',I6,/           &
     &                   ' Number of atoms per unit cell    : ',I6,/    &
     &                   ' Number of different atoms        : ',I6,/)   
 1200 FORMAT     ('    Element : ',A9,' rel. abundance : ',F5.3,        &
     &                   '  (',I6,' atoms)')                            
      END SUBROUTINE chem_elem                      
!*****7*****************************************************************
END MODULE chem_aver_mod
