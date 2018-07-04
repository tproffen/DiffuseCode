MODULE check_blen_mod
!
CONTAINS
!
!*****7*****************************************************************
!
      SUBROUTINE check_blen (x, iatom, rmin, rmax, offset) 
!+                                                                      
!     checks if atom 'iatom' is within rmin -> rmax away from           
!     position x(3). All matching atoms are stored in the               
!     arrays 'atom_env' and 'res_para'. The positions of the            
!     atoms are stored in 'atom_pos' to retain information              
!     about possible periodic boundaries.                               
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE atom_env_mod 
      USE metric_mod
      USE param_mod 
      USE errlist_mod 
      USE lib_f90_allocate_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      REAL x (3), offset (3), rmin, rmax 
      INTEGER iatom 
      INTEGER  :: n_res
!                                                                       
      REAL v (3), dist  !, do_blen 
      INTEGER j 
      LOGICAL lspace 
!                                                                       
      lspace = .true. 
!                                                                       
      DO j = 1, 3 
         atom_pos (j, 0 ) = x (j) 
      ENDDO 
      DO j = 1, 3 
      v (j) = cr_pos (j, iatom) + offset (j) 
      ENDDO 
      dist = do_blen (lspace, x, v) 
      IF (dist.ge.rmin.and.dist.le.rmax) then 
         IF (atom_env (0) .lt.MAX_ATOM_ENV) then 
            IF(atom_env(0) > MAXPAR_RES) THEN
               n_res = MAX(atom_env(0), NINT(MAXPAR_RES*1.1+10),CHEM_MAX_NEIG)
               CALL alloc_param(n_res)
               MAXPAR_RES = n_res
            ENDIF
            IF (atom_env (0) .lt.MAXPAR_RES) then 
               atom_env (0) = atom_env (0) + 1 
               atom_env (atom_env (0) ) = iatom 
               res_para (atom_env (0) ) = dist 
               res_para (0) = float (nint (res_para (0) ) + 1) 
               DO j = 1, 3 
                  atom_pos (j, atom_env (0) ) = v (j) 
               ENDDO 
               atom_dis (   atom_env (0) ) = dist
            ELSE 
               ier_num = - 79 
               ier_typ = ER_APPL 
            ENDIF 
         ELSE 
            ier_num = - 45 
            ier_typ = ER_APPL 
         ENDIF 
      ENDIF 
      END SUBROUTINE check_blen                     
!
!*****7*****************************************************************
!
      SUBROUTINE check_blen_mol (x, imole, rmin, rmax, offset) 
!+                                                                      
!     checks if molecule 'imole' is within rmin -> rmax awaz from       
!     position x(3). All matching molecules are stored in the           
!     arrays 'mole_env' and 'res_para'. The positions of the            
!     molecules are stored in 'mole_pos' to retain information          
!     about possible periodic boundaries.                               
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE metric_mod
      USE molecule_mod 
      USE mole_env_mod 
       
      USE errlist_mod 
      USE param_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      REAL x (3), offset (3), rmin, rmax 
      INTEGER imole 
!                                                                       
      REAL v (3), dist !, do_blen 
      INTEGER i, j 
      LOGICAL lspace 
!                                                                       
      lspace = .true. 
!                                                                       
      i = mole_cont (mole_off (imole) + 1) 
      DO j = 1, 3 
      v (j) = cr_pos (j, i) + offset (j) 
      ENDDO 
      dist = do_blen (lspace, x, v) 
      IF (dist.ge.rmin.and.dist.le.rmax) then 
         IF (mole_env (0) .lt.MAX_MOLE_ENV) then 
            mole_env (0) = mole_env (0) + 1 
            mole_env (mole_env (0) ) = imole 
            res_para (mole_env (0) ) = dist 
            res_para (0) = float (nint (res_para (0) ) + 1) 
            DO j = 1, 3 
            mole_pos (j, mole_env (0) ) = v (j) 
            ENDDO 
         ELSE 
            ier_num = - 45 
            ier_typ = ER_APPL 
         ENDIF 
      ENDIF 
      END SUBROUTINE check_blen_mol                 
END MODULE check_blen_mod
