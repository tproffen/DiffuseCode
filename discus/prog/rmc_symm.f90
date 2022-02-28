MODULE rmc_symm_mod
!
CONTAINS
!********************************************************************** 
      SUBROUTINE rmc_symmetry (n_hkl, hkl, mat_hkl, max_hkl,            &
      lfriedel_remove, lacentric)                                       
!-                                                                      
!     Performs the space group symmetry on the given [hkl].             
!     Altered from 'subroutine symmetry' in 'structur.f'.               
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE generate_mod 
      USE errlist_mod 
use precision_mod
!
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER, INTENT(OUT) :: n_hkl 
      INTEGER, INTENT(IN)  :: max_hkl
      REAL(kind=PREC_DP)   , INTENT(OUT) :: hkl (4, max_hkl)
      REAL(kind=PREC_DP)   , INTENT(OUT) :: mat_hkl (4, 4, max_hkl) 
      LOGICAL, INTENT(IN)  :: lfriedel_remove 
      LOGICAL, INTENT(OUT) :: lacentric 
!                                                                       
      REAL(kind=PREC_DP) :: x (4), y (4), xmat (4, 4), wmat (4, 4), eps 
      INTEGER i, j, k, l, ii, iii, igs, igg, ipg, ia, iaa 
!                                                                       
!                                                                       
      DATA eps / 0.00001 / 
!                                                                       
      n_hkl = 1 
      ii = n_hkl 
      lacentric = .true. 
!                                                                       
      IF (cr_spcgrno.eq.0) cr_spcgrno = 1 
!                                                                       
!------ Loop over all generators of spacegroup cr_spcgrno               
!                                                                       
      DO igs = 1, generspcgr (0, cr_spcgrno) 
      igg = generspcgr (igs, cr_spcgrno) 
      iii = n_hkl 
      DO i = 1, 4 
      DO j = 1, 4 
      xmat (i, j) = 0.0 
      ENDDO 
      xmat (i, i) = 1.0 
      ENDDO 
!                                                                       
!     --Loop over all powers of generator igg                           
!                                                                       
      DO ipg = 1, generpower (igg) 
      DO i = 1, 4 
      DO j = 1, 4 
      wmat (i, j) = 0.0 
      DO k = 1, 4 
      wmat (i, j) = wmat (i, j) + generators (i, k, igg) * xmat (k, j) 
      ENDDO 
      ENDDO 
      ENDDO 
!                                                                       
!------ ----remove translation part of generator                        
!                                                                       
!         do i=1,3                                                      
!           wmat(i,4)=0.0                                               
!           wmat(4,i)=0.0                                               
!         ENDDO                                                         
!         wmat(4,4)=1.0                                                 
!                                                                       
!     ----Multiply first point with current generator                   
!         Since x(4) = 0.0, removal of translational part is no longer  
!         necessary. The part is indeed needed for proper phase         
!           transferal                                                  
!                                                                       
      DO i = 1, 4 
      x (i) = hkl (i, ii) 
      ENDDO 
      x (4) = 0.0 
!     CALL trans (x, wmat, y, 4) 
      y = matmul(wmat, x)
!                                                                       
!     ----check for existence with all previously generated hkl's '     
!                                                                       
      DO ia = ii, iii 
      IF ( (abs (y (1) - hkl (1, ia) ) .lt.eps.and.abs (y (2) - hkl (2, &
      ia) ) .lt.eps.and.abs (y (3) - hkl (3, ia) ) .lt.eps) ) goto 10   
      IF ( (abs ( - y (1) - hkl (1, ia) ) .lt.eps.and.abs ( - y (2)     &
      - hkl (2, ia) ) .lt.eps.and.abs ( - y (3) - hkl (3, ia) ) .lt.eps)&
      ) then                                                            
         lacentric = .false. 
         IF (lfriedel_remove) then 
            GOTO 10 
         ENDIF 
      ENDIF 
      ENDDO 
!                                                                       
      IF (n_hkl.lt.max_hkl) then 
         n_hkl = n_hkl + 1 
         DO i = 1, 4 
         hkl (i, n_hkl) = y (i) 
         DO j = 1, 4 
         mat_hkl (i, j, n_hkl) = wmat (i, j) 
         ENDDO 
         ENDDO 
      ELSE 
         ier_num = - 13 
         ier_typ = ER_RMC 
         RETURN 
      ENDIF 
!                                                                       
!     ----Multiply all other points with current generator              
!                                                                       
      DO iaa = ii + 1, iii 
      DO i = 1, 4 
      x (i) = hkl (i, iaa) 
      ENDDO 
      x (4) = 0.0 
!     CALL trans (x, wmat, y, 4) 
      y = matmul(wmat, x)
!                                                                       
      IF (n_hkl.lt.max_hkl) then 
         n_hkl = n_hkl + 1 
         DO i = 1, 4 
         hkl (i, n_hkl) = y (i) 
         DO j = 1, 4 
         mat_hkl (i, j, n_hkl) = 0.0 
         DO l = 1, 4 
         mat_hkl (i, j, n_hkl) = mat_hkl (i, j, n_hkl) + mat_hkl (i, l, &
         iaa) * wmat (l, j)                                             
         ENDDO 
         ENDDO 
         ENDDO 
      ELSE 
         ier_num = - 13 
         ier_typ = ER_RMC 
         RETURN 
      ENDIF 
      ENDDO 
!                                                                       
!     --set xmat to current power of generator                          
!                                                                       
      DO i = 1, 4 
      DO j = 1, 4 
      xmat (i, j) = wmat (i, j) 
      ENDDO 
      ENDDO 
      ENDDO 
   10 CONTINUE 
      ENDDO 
      END SUBROUTINE rmc_symmetry                   
!********************************************************************** 
END MODULE rmc_symm_mod
