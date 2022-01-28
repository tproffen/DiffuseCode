MODULE check_blen_mod
!
CONTAINS
!
!*****7*****************************************************************
!
SUBROUTINE check_blen(x, iatom, rmin, rmax, offset) 
!+                                                                      
!     checks if atom 'iatom' is within rmin -> rmax away from           
!     position x(3). All matching atoms are stored in the               
!     arrays 'atom_env' and 'res_para'. The positions of the            
!     atoms are stored in 'atom_pos' to retain information              
!     about possible periodic boundaries.                               
!-                                                                      
USE crystal_mod 
USE atom_env_mod 
USE metric_mod
USE param_mod 
USE errlist_mod 
USE lib_f90_allocate_mod 
use precision_mod
!                                                                       
IMPLICIT none 
!                                                                       
REAL(kind=PREC_DP)              , DIMENSION(3), INTENT(IN) :: x
INTEGER                         , INTENT(IN) :: iatom 
REAL(kind=PREC_DP)              , INTENT(IN) :: rmin
REAL(kind=PREC_DP)              , INTENT(IN) :: rmax
REAL(kind=PREC_DP), DIMENSION(3), INTENT(IN) :: offset
!
INTEGER  :: n_res
!                                                                       
REAL(kind=PREC_DP), DIMENSION(3) :: v
REAL(kind=PREC_DP) :: dist
INTEGER            :: j 
LOGICAL, PARAMETER :: lspace = .TRUE.
!                                                                       
DO j = 1, 3 
   atom_pos (j, 0 ) = x (j) 
ENDDO 
DO j = 1, 3 
   v (j) = cr_pos (j, iatom) + offset (j) 
ENDDO 
dist = do_blen (lspace, x, v) 
IF (dist >= rmin.AND.dist <= rmax) THEN 
   IF (atom_env(0) < MAX_ATOM_ENV) THEN 
      IF(atom_env(0) > MAXPAR_RES) THEN
         n_res = MAX(atom_env(0), NINT(MAXPAR_RES*1.1+10)) 
         CALL alloc_param(n_res)
         MAXPAR_RES = n_res
      ENDIF
      IF (atom_env(0) < MAXPAR_RES) THEN 
         atom_env(0) = atom_env(0) + 1 
         atom_env(atom_env(0) ) = iatom 
         res_para(atom_env(0) ) = dist 
         res_para(0) = REAL(NINT(res_para(0) ) + 1, kind=PREC_DP) 
         DO j = 1, 3 
            atom_pos(j, atom_env(0) ) = v(j) 
         ENDDO 
         atom_dis(atom_env(0) ) = dist
      ELSE 
         ier_num = - 79 
         ier_typ = ER_APPL 
      ENDIF 
   ELSE 
      ier_num = - 45 
      ier_typ = ER_APPL 
   ENDIF 
ENDIF 
!
END SUBROUTINE check_blen                     
!
!*****7*****************************************************************
!
SUBROUTINE check_blen_mol(x, imole, rmin, rmax, offset) 
!+                                                                      
!     checks if molecule 'imole' is within rmin -> rmax awaz from       
!     position x(3). All matching molecules are stored in the           
!     arrays 'mole_env' and 'res_para'. The positions of the            
!     molecules are stored in 'mole_pos' to retain information          
!     about possible periodic boundaries.                               
!-                                                                      
USE crystal_mod 
USE metric_mod
USE molecule_mod 
USE mole_env_mod 
!
USE errlist_mod 
USE param_mod 
use precision_mod
!                                                                       
IMPLICIT none 
!                                                                       
REAL(kind=PREC_DP)              , DIMENSION(3), INTENT(IN) :: x
INTEGER                         , INTENT(IN) :: imole 
REAL(kind=PREC_DP)              , INTENT(IN) :: rmin
REAL(kind=PREC_DP)              , INTENT(IN) :: rmax 
REAL(kind=PREC_DP), DIMENSION(3), INTENT(IN) :: offset
!                                                                       
REAL(kind=PREC_DP) :: v (3), dist !, do_blen 
INTEGER i, j 
LOGICAL, PARAMETER :: lspace =.TRUE.
!                                                                       
!                                                                       
i = mole_cont(mole_off(imole) + 1) 
DO j = 1, 3 
   v(j) = cr_pos(j, i) + offset(j) 
ENDDO 
dist = do_blen(lspace, x, v) 
IF (dist >= rmin.AND.dist <= rmax) THEN 
   IF (mole_env(0) <   MAX_MOLE_ENV) THEN 
      mole_env(0) = mole_env (0) + 1 
      mole_env(mole_env(0) ) = imole 
      res_para(mole_env(0) ) = dist 
      res_para(0) = REAL(NINT(res_para(0) ) + 1, kind=PREC_DP) 
      DO j = 1, 3 
         mole_pos (j, mole_env(0) ) = v (j) 
      ENDDO 
   ELSE 
      ier_num = - 45 
      ier_typ = ER_APPL 
   ENDIF 
ENDIF 
!
END SUBROUTINE check_blen_mol                 
!
!*******************************************************************************
!
END MODULE check_blen_mod
