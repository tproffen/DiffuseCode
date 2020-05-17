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
USE atom_name
USE chem_mod 
USE errlist_mod 
USE param_mod 
USE prompt_mod 
USE lib_f90_allocate_mod
IMPLICIT none 
!                                                                       
LOGICAL, INTENT(IN) :: lout    ! Print output if true
LOGICAL, INTENT(IN) :: lsite   ! Treat different atoms on each site as one 
!                                                                       
REAL, DIMENSION(3) ::  p , ez
INTEGER, DIMENSION(3) :: iez
!REAL,    DIMENSION(:,:,:), ALLOCATABLE :: chem_ave_posit
!REAL,    DIMENSION(:,:,:), ALLOCATABLE :: chem_ave_sigma
!
INTEGER            :: i, j, k, ii, jj, kk, ia, is, nvalues
INTEGER            :: n_res
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
   IF(ALLOCATED(chem_ave_posit)) DEALLOCATE(chem_ave_posit)
   IF(ALLOCATED(chem_ave_sigma)) DEALLOCATE(chem_ave_sigma)
   ALLOCATE(chem_ave_posit(3,n_atom_cell, MAX(12,cr_nscat)))
   ALLOCATE(chem_ave_sigma(3,n_atom_cell, MAX(12,cr_nscat)))
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
                  IF (cr_iscat (ia) .eq.chem_ave_iscat (ii, kk) ) then 
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
!------ output of average and sigma                                     
!                                                                       
IF (lout) write (output_io, 1000) 
ia = cr_icc (1) * cr_icc (2) * cr_icc (3) 
IF(lsite) THEN           ! one pos for all types on a single site
   nvalues = 0
   DO i = 1, cr_ncatoms 
      DO j = 1, 3 
!        chem_ave_pos (j, i) = chem_ave_pos (j, i) / REAL(ia) 
         chem_ave_sig (j, i) = chem_ave_sig (j, i) / REAL(ia) !-          &
!        chem_ave_pos (j, i) **2                                           
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
      n_res = MAX(9 * cr_ncatoms,MAXPAR_RES, CHEM_MAX_NEIG)
      CALL alloc_param(n_res)
      MAXPAR_RES = n_res
!     ier_typ = ER_CHEM 
!     ier_num = - 2 
!  ELSE 
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
!  ENDIF 
ELSE
   nvalues = 0
   DO i = 1, cr_ncatoms 
      DO k = 1, chem_ave_n (i) 
         DO j = 1, 3 
!           chem_ave_posit (j, i, k) = chem_ave_posit (j, i, k) / chem_ave_bese(i,k) 
            chem_ave_sigma (j, i, k) = chem_ave_sigma (j, i, k) / chem_ave_bese(i,k) !- &
!           chem_ave_posit (j, i, k) **2                                           
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
      n_res = MAX(9 * cr_ncatoms,MAXPAR_RES, CHEM_MAX_NEIG)
      CALL alloc_param(n_res)
      MAXPAR_RES = n_res
!     ier_typ = ER_CHEM 
!     ier_num = - 2 
!  ELSE 
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
!  ENDIF 
ENDIF
!IF(.not. lsite) THEN
!   DEALLOCATE(chem_ave_posit)
!   DEALLOCATE(chem_ave_sigma)
!ENDIF
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
      USE discus_config_mod 
      USE crystal_mod 
      USE atom_name
      USE chem_mod 
!                                                                       
      USE errlist_mod 
      USE lib_f90_allocate_mod 
      USE param_mod 
      USE prompt_mod 
      USE wink_mod 
      IMPLICIT none 
!
LOGICAL, INTENT(IN) :: lout 
!                                                                       
REAL             :: proz 
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
   cr_amount (cr_iscat (i) ) = cr_amount (cr_iscat (i) ) + 1 
ENDDO 
!                                                                       
!------ write output                                                    
!                                                                       
IF (lout) write (output_io, 1000) (cr_icc (i), i = 1, 3) 
IF (lout) write (output_io, 1100) cr_natoms, cr_ncatoms, cr_nscat 
IF (i >  MAXPAR_RES) THEN 
   n_res = MAX(cr_nscat, MAXPAR_RES, CHEM_MAX_NEIG)
   CALL alloc_param(n_res)
   MAXPAR_RES = n_res
ENDIF
res_para (0) = REAL(cr_nscat) + 1 
DO i = 0, cr_nscat 
   cr_u2aver = cr_u2aver + cr_dw(i) * cr_amount (i)
!   IF (cr_at_lis(cr_iscat (i)) /= 'VOID') then
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
!
IMPLICIT NONE
!
REAL, DIMENSION(3), INTENT(OUT) :: com   ! Center of Mass
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
      IF(cr_iscat(iatom)==chem_ave_iscat(isite,k)) THEN
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
   at_name_d = at_name (cr_iscat (iatom) )
   WRITE(output_io, 1000) iatom, at_name_d, res_para(1:3)
1000 FORMAT(i8,1x,a9,3(2x,f12.7))
ENDIF
!
END SUBROUTINE get_displacement
!
!*******************************************************************************
END MODULE chem_aver_mod
