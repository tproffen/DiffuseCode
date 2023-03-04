module discus_update_mod
!
contains
!
SUBROUTINE discus_ersetz_para (ikl, iklz, string, ll, ww, maxw, ianz) 
!                                                                       
!-                                                                      
!       replaces a substring in an expression by the value of the       
!       appropriate Parameter.                                          
!     This is needed if the parameter is read.                          
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE atom_env_mod 
      USE molecule_mod 
      USE mole_env_mod 
      USE pdf_mod
      USE wyckoff_mod
      USE blanks_mod
      USE errlist_mod 
      USE lib_upd_mod
USE lib_length
USE lib_errlist_func
      USE param_mod 
      USE random_mod 
USE precision_mod
      USE precision_mod
      IMPLICIT none 
!                                                                       
      INTEGER,                    INTENT(IN   ) :: ikl
      INTEGER,                    INTENT(IN   ) :: iklz
      CHARACTER (LEN=*),          INTENT(INOUT) :: string 
      INTEGER,                    INTENT(INOUT) :: ll
      INTEGER,                    INTENT(IN   ) :: maxw
      REAL(KIND=PREC_DP)   , DIMENSION(1:maxw), INTENT(IN   ) :: ww
      INTEGER,                    INTENT(IN   ) :: ianz
!                                                                       
      CHARACTER(LEN=MAX(PREC_STRING,LEN(string))) zeile 
!                                                                       
      INTEGER laenge, ltyp, kpara, kpara2
      INTEGER lcomm 
!
CALL lib_ersetz_para (ikl, iklz, string, ll, ww, maxw, ianz)
IF(ier_num == 0) RETURN
CALL no_error
!                                                                       
      laenge = ll 
      ltyp = 1 
      zeile = ' ' 
      kpara = nint (ww (1) ) 
      kpara2 = 0
      IF (maxw.ge.2) THEN 
         kpara2 = nint (ww (2) ) 
      ENDIF 
!                                                                       
      lcomm = length_com (string, ikl) 
!                                                                       
      IF (lcomm.eq.1) THEN 
!                                                                       
         IF (ikl.gt.lcomm + 1) zeile (1:ikl - lcomm - 1) = string (1:   &
         ikl - lcomm - 1)                                               
!
         IF (string (ikl - 1:ikl - 1) .eq.'x') THEN 
            IF (ianz.eq.1) THEN 
               IF (1.le.kpara.and.kpara.le.NMAX.and.kpara.le.cr_natoms) &
               THEN                                                     
                  WRITE(zeile(ikl - 1:ikl + PREC_WIDTH-2), PREC_F_REAL) cr_pos(1, kpara)
                  zeile (ikl + PREC_MANTIS - lcomm:ikl + PREC_MANTIS - lcomm) = 'd' 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
         ELSEIF (string (ikl - 1:ikl - 1) .eq.'y') THEN 
            IF (ianz.eq.1) THEN 
               IF (1.le.kpara.and.kpara.le.NMAX.and.kpara.le.cr_natoms) &
               THEN                                                     
                  WRITE(zeile(ikl - 1:ikl + PREC_WIDTH-2), PREC_F_REAL) cr_pos(2, kpara)
                  zeile (ikl + PREC_MANTIS - lcomm:ikl + PREC_MANTIS - lcomm) = 'd' 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
         ELSEIF (string (ikl - 1:ikl - 1) .eq.'z') THEN 
            IF (ianz.eq.1) THEN 
               IF (1.le.kpara.and.kpara.le.NMAX.and.kpara.le.cr_natoms) &
               THEN                                                     
                  WRITE(zeile(ikl - 1:ikl + PREC_WIDTH-2), PREC_F_REAL) cr_pos(3, kpara)
                  zeile (ikl + PREC_MANTIS - lcomm:ikl + PREC_MANTIS - lcomm) = 'd' 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
         ELSEIF (string (ikl - 1:ikl - 1) .eq.'m') THEN 
            IF (ianz.eq.1) THEN 
               IF (1.le.kpara.and.kpara.le.NMAX.and.kpara.le.cr_natoms) THEN
                  WRITE(zeile(ikl - 1:ikl + PREC_WIDTH-2), PREC_F_INTE) cr_iscat(kpara)
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
         ELSEIF (string (ikl - 1:ikl - 1) .eq.'b') THEN 
            IF (ianz.eq.1) THEN 
               IF (1.le.kpara.and.kpara.le.cr_nscat) THEN 
                  WRITE(zeile(ikl - 1:ikl + PREC_WIDTH-2), PREC_F_REAL) cr_dw(kpara)
                  zeile (ikl + PREC_MANTIS - lcomm:ikl + PREC_MANTIS - lcomm) = 'd' 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
         ELSEIF (string (ikl - 1:ikl - 1) .eq.'n') THEN 
            IF (ianz.eq.1) THEN 
               IF (kpara.eq.1) THEN 
                  WRITE (zeile (ikl - 1:ikl + PREC_WIDTH-2) , PREC_F_INTE) cr_natoms 
               ELSEIF (kpara.eq.2) THEN 
                  WRITE (zeile (ikl - 1:ikl + PREC_WIDTH-2) , PREC_F_INTE) cr_nscat 
               ELSEIF (kpara.eq.3) THEN 
                  WRITE (zeile (ikl - 1:ikl + PREC_WIDTH-2) , PREC_F_INTE) cr_ncatoms 
               ELSEIF (kpara.eq.4) THEN 
                  WRITE (zeile (ikl - 1:ikl + PREC_WIDTH-2) , PREC_F_INTE)            &
                  mole_num_mole                                         
               ELSEIF (kpara.eq.5) THEN 
                  WRITE (zeile (ikl - 1:ikl + PREC_WIDTH-2) , PREC_F_INTE)            &
                  mole_num_type                                         
               ELSEIF (kpara.eq.6) THEN 
                  WRITE (zeile (ikl - 1:ikl + PREC_WIDTH-2) , PREC_F_INTE)            &
                  mole_num_unit                                         
               ELSEIF (kpara.eq.7) THEN 
                  WRITE (zeile (ikl - 1:ikl + PREC_WIDTH-2) , PREC_F_INTE)            &
                  cr_ncreal
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
         ELSE 
            ier_num = - 2 
            ier_typ = ER_FORT 
         ENDIF 
!                                                                       
      ELSEIF (lcomm.eq.8) THEN 
                                                                        
         IF (string (ikl - 8:ikl - 1) .eq.'mol_cont') THEN 
            IF (ianz.eq.2) THEN 
               IF (ikl.gt.lcomm + 1) zeile (1:ikl - lcomm - 1) = string &
               (1:ikl - lcomm - 1)                                      
               IF (0.lt.kpara.and.kpara.le.mole_num_mole) THEN 
                  IF (0.eq.kpara2) THEN 
                     WRITE (zeile (ikl - 8:ikl + PREC_WIDTH-2) , PREC_F_INTE)         &
                     mole_len (kpara)                                   
                  ELSEIF (0.le.kpara2.and.kpara2.le.mole_len (kpara) )  &
                  THEN                                                  
                     WRITE (zeile (ikl - 8:ikl + PREC_WIDTH-2) , PREC_F_INTE)         &
                     mole_cont (mole_off (kpara) + kpara2)              
                  ELSE 
                     ier_num = - 74 
                     ier_typ = ER_APPL 
                     RETURN 
                  ENDIF 
               ELSE 
                  ier_num = - 63 
                  ier_typ = ER_APPL 
                  RETURN 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
!                                                                       
         ELSEIF (string (ikl - 8:ikl - 1) .eq.'mol_biso') THEN 
            IF (ianz.eq.1) THEN 
               IF (ikl.gt.lcomm + 1) zeile (1:ikl - lcomm - 1) = string &
               (1:ikl - lcomm - 1)                                      
               IF (0.lt.kpara.and.kpara.le.mole_num_type) THEN 
                  WRITE (zeile (ikl - 8:ikl + PREC_WIDTH-2) , PREC_F_REAL)        &
                  mole_biso (kpara)                                     
                  zeile (ikl + PREC_MANTIS - lcomm:ikl + PREC_MANTIS - lcomm) = 'd' 
               ELSE 
                  ier_num = - 64 
                  ier_typ = ER_APPL 
                  RETURN 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
!                                                                       
         ELSEIF (string (ikl - 8:ikl - 1) .eq.'mol_clin') THEN 
            IF (ianz.eq.1) THEN 
               IF (ikl.gt.lcomm + 1) zeile (1:ikl - lcomm - 1) = string &
               (1:ikl - lcomm - 1)                                      
               IF (0.lt.kpara.and.kpara.le.mole_num_type) THEN 
                  WRITE (zeile (ikl - 8:ikl + PREC_WIDTH-2) , PREC_F_REAL)        &
                  mole_clin (kpara)                                     
                  zeile (ikl + PREC_MANTIS - lcomm:ikl + PREC_MANTIS - lcomm) = 'd' 
               ELSE 
                  ier_num = - 64 
                  ier_typ = ER_APPL 
                  RETURN 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
!                                                                       
         ELSEIF (string (ikl - 8:ikl - 1) .eq.'mol_cqua') THEN 
            IF (ianz.eq.1) THEN 
               IF (ikl.gt.lcomm + 1) zeile (1:ikl - lcomm - 1) = string &
               (1:ikl - lcomm - 1)                                      
               IF (0.lt.kpara.and.kpara.le.mole_num_type) THEN 
                  WRITE (zeile (ikl - 8:ikl + PREC_WIDTH-2) , PREC_F_REAL)        &
                  mole_cqua (kpara)                                     
                  zeile (ikl + PREC_MANTIS - lcomm:ikl + PREC_MANTIS - lcomm) = 'd' 
               ELSE 
                  ier_num = - 64 
                  ier_typ = ER_APPL 
                  RETURN 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
!                                                                       
         ELSEIF (string (ikl - 8:ikl - 1) .eq.'mol_dens') THEN 
            IF (ianz.eq.1) THEN 
               IF (ikl.gt.lcomm + 1) zeile (1:ikl - lcomm - 1) = string &
               (1:ikl - lcomm - 1)                                      
               IF (0.lt.kpara.and.kpara.le.mole_num_type) THEN 
                  WRITE (zeile (ikl - 8:ikl + PREC_WIDTH-2) , PREC_F_REAL)        &
                  mole_dens (kpara)                                     
                  zeile (ikl + PREC_MANTIS - lcomm:ikl + PREC_MANTIS - lcomm) = 'd' 
               ELSE 
                  ier_num = - 64 
                  ier_typ = ER_APPL 
                  RETURN 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
!                                                                       
         ELSEIF (string (ikl - 8:ikl - 1) .eq.'mol_type') THEN 
            IF (ianz.eq.1) THEN 
               IF (ikl.gt.lcomm + 1) zeile (1:ikl - lcomm - 1) = string &
               (1:ikl - lcomm - 1)                                      
               IF (0.lt.kpara.and.kpara.le.mole_num_mole) THEN 
                  WRITE (zeile (ikl - 8:ikl + PREC_WIDTH-2) , PREC_F_INTE) mole_type (&
                  kpara)                                                
               ELSE 
                  ier_num = - 64 
                  ier_typ = ER_APPL 
                  RETURN 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
!                                                                       
         ELSEIF (string (ikl - 8:ikl - 1) .eq.'pdf_dens') THEN 
            IF (ianz.eq.1) THEN 
               IF (ikl.gt.lcomm + 1) zeile (1:ikl - lcomm - 1) = string &
               (1:ikl - lcomm - 1)                                      
               WRITE (zeile (ikl - 8:ikl + PREC_WIDTH-2) , PREC_F_REAL) pdf_rho0
               zeile (ikl + PREC_MANTIS - lcomm:ikl + PREC_MANTIS - lcomm) = 'd' 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
!                                                                       
         ELSEIF (string (ikl - 8:ikl - 1) .eq.'pdf_scal') THEN 
            IF (ianz.eq.1) THEN 
               IF (ikl.gt.lcomm + 1) zeile (1:ikl - lcomm - 1) = string &
               (1:ikl - lcomm - 1)                                      
               WRITE (zeile (ikl - 8:ikl + PREC_WIDTH-2) , PREC_F_REAL) pdf_scale
               zeile (ikl + PREC_MANTIS - lcomm:ikl + PREC_MANTIS - lcomm) = 'd' 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
!                                                                       
         ELSEIF (string (ikl - 8:ikl - 1) .eq.'ref_para') THEN 
            IF (ianz.eq.1) THEN 
               IF (ikl.gt.lcomm + 1) zeile (1:ikl - lcomm - 1) = string &
               (1:ikl - lcomm - 1)                                      
               IF (0.lt.kpara.and.kpara.le.MAXPAR_REF   ) THEN 
                  WRITE (zeile (ikl - 8:ikl + PREC_WIDTH-2) , PREC_F_REAL)        &
                  ref_para (kpara)                                     
                  zeile (ikl + PREC_MANTIS - lcomm:ikl + PREC_MANTIS - lcomm) = 'd' 
               ELSE 
                  ier_num = -133 
                  ier_typ = ER_APPL 
                  RETURN 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
         ELSE 
            ier_num = - 2 
            ier_typ = ER_FORT 
         ENDIF 
!                                                                       
      ELSEIF (lcomm.eq.7) THEN 
!                                                                       
         IF (string (ikl - 7:ikl - 1) .eq.'mol_len') THEN 
            IF (ianz.eq.1) THEN 
               IF (ikl.gt.lcomm + 1) zeile (1:ikl - lcomm - 1) = string &
               (1:ikl - lcomm - 1)                                      
               IF (0.lt.kpara.and.kpara.le.mole_num_mole) THEN 
                  WRITE (zeile (ikl - 7:ikl + PREC_WIDTH-2) , PREC_F_INTE) mole_len ( &
                  kpara)                                                
               ELSE 
                  ier_num = - 63 
                  ier_typ = ER_APPL 
                  RETURN 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
         ELSEIF (string (ikl - 7:ikl - 1) .eq.'at_name') THEN 
            IF (ianz.eq.1) THEN 
               IF (ikl.gt.lcomm + 1) zeile (1:ikl - lcomm - 1) = string &
               (1:ikl - lcomm - 1)                                      
               IF (1.le.kpara.and.kpara.le.NMAX.and.kpara.le.cr_natoms) &
               THEN                                                     
                  zeile(ikl-7:ikl+PREC_WIDTH-2) = ''''//cr_at_lis(cr_iscat(kpara))//''''                                                     
               ELSE 
                  ier_num = - 105 
                  ier_typ = ER_APPL 
                  RETURN 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
         ELSEIF (string (ikl - 7:ikl - 1) .eq.'at_type') THEN 
            IF (ianz.eq.1) THEN 
               IF (ikl.gt.lcomm + 1) zeile (1:ikl - lcomm - 1) = string &
               (1:ikl - lcomm - 1)                                      
               IF (0.le.kpara.and.kpara.le.NMAX.and.kpara.le.cr_nscat) THEN
                  zeile(ikl-7:ikl+PREC_WIDTH-2) = ''''//cr_at_lis (kpara)//''''
               ELSE 
                  ier_num = - 122 
                  ier_typ = ER_APPL 
                  RETURN 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
         ELSEIF (string (ikl - 7:ikl - 1) .eq.'in_mole') THEN 
            IF (ianz.eq.1) THEN 
               IF (ikl.gt.lcomm + 1) zeile (1:ikl - lcomm - 1) = string &
               (1:ikl - lcomm - 1)                                      
               IF (1.le.kpara.and.kpara.le.NMAX.and.kpara.le.cr_natoms) &
               THEN                                                     
                  WRITE(zeile(ikl-7:ikl + PREC_WIDTH-2),PREC_F_INTE) cr_mole(kpara)
               ELSE 
                  ier_num = - 105 
                  ier_typ = ER_APPL 
                  RETURN 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
!                                                                       
         ELSE 
            ier_num = - 2 
            ier_typ = ER_FORT 
         ENDIF 
      ELSEIF(lcomm.eq.5) THEN 
         IF(string(ikl - 5:ikl - 1) .eq.'sym_n') THEN 
            IF(ianz.eq.1) THEN 
               WRITE(zeile(ikl-5:ikl+PREC_WIDTH-2),PREC_F_INTE) spc_n
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
!                                                                       
         ELSE 
            ier_num = - 2 
            ier_typ = ER_FORT 
         ENDIF 
!                                                                       
      ELSEIF (lcomm.eq.4) THEN 
!                                                                       
         IF (string (ikl - 4:ikl - 1) .eq.'cdim') THEN 
            IF (ianz.eq.2) THEN 
               IF (ikl.gt.lcomm + 1) zeile (1:ikl - lcomm - 1) = string(1:ikl - lcomm - 1)
      IF (1.le.kpara.and.kpara.le.3.and.1.le.kpara2.and.kpara2.le.2) THEN
                  WRITE (zeile (ikl - 4:ikl + PREC_WIDTH-2) , PREC_F_REAL) cr_dim &
                  (kpara, kpara2)                                       
                  zeile (ikl + PREC_MANTIS - lcomm:ikl + PREC_MANTIS - lcomm) = 'd' 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'prop') THEN 
            IF (ianz.eq.1) THEN 
               IF (1.le.kpara.and.kpara.le.NMAX.and.kpara.le.cr_natoms) THEN
                  WRITE (zeile (ikl - 4:ikl + PREC_WIDTH-2) , PREC_F_INTE) cr_prop (  &
                  kpara)                                                
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'surf') THEN 
            IF (ianz.eq.2) THEN 
               IF (1.le.kpara.and.kpara.le.NMAX.and.kpara.le.cr_natoms .AND. &
                  1<=kpara2 .AND. kpara2 <=3                                 ) THEN
                  WRITE(zeile(ikl - 4:ikl + PREC_WIDTH-2) , PREC_F_INTE) cr_surf(kpara2,kpara)
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'magn') THEN   ! MAGNETIC_WORK
            IF (ianz.eq.2) THEN 
               IF (1.le.kpara.and.kpara.le.NMAX.and.kpara.le.cr_natoms .AND. &
                  0<=kpara2 .AND. kpara2 <=3                                 ) THEN
                  WRITE(zeile(ikl - 4:ikl + PREC_WIDTH-2) , PREC_F_REAL) cr_magn(kpara2,kpara)
                  zeile (ikl + PREC_MANTIS - lcomm:ikl + PREC_MANTIS - lcomm) = 'e' 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'menv') THEN 
            IF (ianz.eq.1) THEN 
               IF (ikl.gt.lcomm + 1) zeile (1:ikl - lcomm - 1) = string(1:ikl - lcomm - 1)
               IF (0.le.kpara.and.kpara.le.MAXPAR_RES) THEN 
               WRITE (zeile (ikl - 4:ikl + PREC_WIDTH-2) , PREC_F_INTE) mole_env (kpara) 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'seed') THEN 
            IF (ianz.eq.1) THEN 
               IF(ikl.gt.lcomm + 1) zeile(1:ikl - lcomm - 1) = string(1:ikl - lcomm - 1)
               WRITE (zeile (ikl - 4:ikl + PREC_WIDTH-2) , PREC_F_INTE) idum
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'rvol') THEN 
            IF (ianz.eq.1) THEN 
               IF(ikl.gt.lcomm + 1) zeile(1:ikl - lcomm - 1) = string(1:ikl - lcomm - 1)
               WRITE (zeile (ikl - 4:ikl + PREC_WIDTH-2) , PREC_F_REAL) cr_vr 
               zeile (ikl + PREC_MANTIS - lcomm:ikl + PREC_MANTIS - lcomm) = 'd' 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
         ELSE 
            ier_num = - 2 
            ier_typ = ER_FORT 
         ENDIF 
      ELSEIF (lcomm.eq.3) THEN 
!                                                                       
         IF (string (ikl - 3:ikl - 1) .eq.'env') THEN 
            IF (ianz.eq.1) THEN 
               IF(ikl.gt.lcomm + 1) zeile(1:ikl - lcomm - 1) = string(1:ikl - lcomm - 1)
               IF (0.le.kpara.and.kpara.le.MAXPAR_RES) THEN 
                  WRITE (zeile (ikl - 3:ikl + PREC_WIDTH-2) , PREC_F_INTE) atom_env (kpara)
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
!                                                                       
!------ - Lattice parameters for current phase lat[n]                   
!                                                                       
         ELSEIF (string (ikl - 3:ikl - 1) .eq.'lat') THEN 
            IF (ianz.eq.1) THEN 
               IF (ikl.gt.lcomm + 1) zeile (1:ikl - lcomm - 1) = string &
               (1:ikl - lcomm - 1)                                      
               IF (kpara.ge.1.and.kpara.le.3) THEN 
                  WRITE(zeile(ikl - 3:ikl + PREC_WIDTH-2), PREC_F_REAL) cr_a0(kpara)
                  zeile (ikl + PREC_MANTIS - lcomm:ikl + PREC_MANTIS - lcomm) = 'd' 
               ELSEIF (kpara.ge.4.and.kpara.le.6) THEN 
                  WRITE(zeile(ikl - 3:ikl + PREC_WIDTH-2), PREC_F_REAL) cr_win(kpara - 3)
                  zeile (ikl + PREC_MANTIS - lcomm:ikl + PREC_MANTIS - lcomm) = 'd' 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
         ELSEIF(string(ikl - 3:ikl - 1)  == 'occ') THEN 
            IF(ianz.eq.1) THEN 
               IF(ikl >  lcomm + 1) zeile(1:ikl - lcomm - 1) = string(1:ikl - lcomm - 1)
               IF(1 <= kpara.AND.kpara <= cr_nscat) THEN 
                  WRITE (zeile (ikl - 3:ikl + PREC_WIDTH-2) , PREC_F_REAL) cr_occ (kpara)
                  zeile (ikl + PREC_MANTIS - lcomm:ikl + PREC_MANTIS - lcomm) = 'd' 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
         ELSEIF (string (ikl - 3:ikl - 1) .eq.'vol') THEN 
            IF (ianz.eq.1) THEN 
               IF(ikl.gt.lcomm + 1) zeile(1:ikl - lcomm - 1) = string(1:ikl - lcomm - 1)
               WRITE (zeile (ikl - 3:ikl + PREC_WIDTH-2) , PREC_F_REAL) cr_v 
               zeile (ikl + PREC_MANTIS - lcomm:ikl + PREC_MANTIS - lcomm) = 'd' 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
         ELSE 
            ier_num = - 2 
            ier_typ = ER_FORT 
         ENDIF 
      ELSE 
         ier_num = - 2 
         ier_typ = ER_FORT 
      ENDIF 
      IF (ier_num.eq.0) THEN 
         ll = laenge+PREC_WIDTH - ltyp - (iklz - ikl + 1) 
         IF(iklz + 1.le.laenge) zeile(ikl + PREC_WIDTH-1:ll) = string(iklz + 1:laenge)
         string = zeile 
      ELSE 
         ll = min (40, laenge) 
         WRITE (ier_msg (1), 8000) string (1:ll) 
      ENDIF 
      ll  = LEN_TRIM(string)
!                                                                       
 8000 FORMAT    (a) 
      END SUBROUTINE discus_ersetz_para                    
!
!*****7*****************************************************************
!
SUBROUTINE discus_upd_para (ctype, ww, maxw, wert, ianz, cstring, substr) 
!-                                                                      
!       updates the parameter specified by ctype, index ww  to the      
!       new value of wert                                               
!+                                                                      
USE discus_config_mod 
USE crystal_mod 
USE atom_env_mod 
USE molecule_mod 
USE do_molecule_alloc
USE prop_para_mod 
USE spcgr_apply, ONLY: setup_lattice
USE surface_mod
!
USE errlist_mod 
USE param_mod 
USE lib_errlist_func
USE lib_upd_mod
USE precision_mod
!
IMPLICIT none 
!                                                                       
CHARACTER (LEN=*),          INTENT(IN) :: ctype 
INTEGER,                    INTENT(IN) :: maxw
INTEGER,                    INTENT(IN) :: ianz 
INTEGER, DIMENSION(1:MAXW), INTENT(IN) :: ww
REAL(KIND=PREC_DP)        , INTENT(IN) :: wert 
CHARACTER (LEN=*),          INTENT(IN) :: cstring
INTEGER, DIMENSION(2), INTENT(IN)    :: substr ! Indices of substring
!
INTEGER :: l
INTEGER :: iwert, owert
!
CALL lib_upd_para (ctype, ww, maxw, wert, ianz, cstring, substr)
IF(ier_num==0 .OR. (ier_num==-40 .AND. ier_typ==ER_FORT)) RETURN
CALL no_error
!                                                                       
      IF (ctype.eq.'x') THEN 
         IF (ianz.eq.1) THEN 
            IF (0.lt.ww (1) .and.ww (1) .le.NMAX.and.ww (1)             &
            .le.cr_natoms) THEN                                         
               cr_pos (1, ww (1) ) = wert 
               l = 1 
               cr_dim (l, 1) = min (cr_dim (l, 1), cr_pos (l, ww (1) ))
               cr_dim (l, 2) = max (cr_dim (l, 2), cr_pos (l, ww (1) ))
            ELSE 
               ier_num = - 8 
               ier_typ = ER_FORT 
            ENDIF 
         ELSE 
            ier_num = - 13 
            ier_typ = ER_FORT 
            RETURN 
         ENDIF 
      ELSEIF (ctype.eq.'y') THEN 
         IF (ianz.eq.1) THEN 
            IF (0.lt.ww (1) .and.ww (1) .le.NMAX.and.ww (1)             &
            .le.cr_natoms) THEN                                         
               cr_pos (2, ww (1) ) = wert 
               l = 2 
               cr_dim (l, 1) = min (cr_dim (l, 1), cr_pos (l, ww (1) ))
               cr_dim (l, 2) = max (cr_dim (l, 2), cr_pos (l, ww (1) ))
            ELSE 
               ier_num = - 8 
               ier_typ = ER_FORT 
            ENDIF 
         ELSE 
            ier_num = - 13 
            ier_typ = ER_FORT 
            RETURN 
         ENDIF 
      ELSEIF (ctype.eq.'z') THEN 
         IF (ianz.eq.1) THEN 
            IF (0.lt.ww (1) .and.ww (1) .le.NMAX.and.ww (1)             &
            .le.cr_natoms) THEN                                         
               cr_pos (3, ww (1) ) = wert 
               l = 3 
               cr_dim (l, 1) = min (cr_dim (l, 1), cr_pos (l, ww (1)))
               cr_dim (l, 2) = max (cr_dim (l, 2), cr_pos (l, ww (1)))
            ELSE 
               ier_num = - 8 
               ier_typ = ER_FORT 
            ENDIF 
         ELSE 
            ier_num = - 13 
            ier_typ = ER_FORT 
            RETURN 
         ENDIF 
      ELSEIF (ctype.eq.'m') THEN 
         IF (ianz.eq.1) THEN 
            IF(0 < ww(1) .AND. ww(1) <= NMAX .AND. ww(1) <= cr_natoms) THEN                                         
               IF(0 <= INT(wert) .AND. INT(wert) <= MAXSCAT .AND. &
                  INT(wert) <= cr_nscat                             ) THEN                                 
                  cr_iscat(ww (1) ) = INT(wert) 
                  IF(cr_iscat(ww(1))>0) THEN
                     cr_prop(ww(1)) = IBSET(cr_prop(ww(1)), PROP_NORMAL)
                  ENDIF 
               ELSE 
                  ier_num = - 97 
                  ier_typ = ER_APPL 
               ENDIF 
            ELSE 
               ier_num = - 8 
               ier_typ = ER_FORT 
            ENDIF 
         ELSE 
            ier_num = - 13 
            ier_typ = ER_FORT 
            RETURN 
         ENDIF 
      ELSEIF (ctype.eq.'b') THEN 
         IF (ianz.eq.1) THEN 
            IF (0.lt.ww (1) .and.ww (1) .le.cr_nscat) THEN 
               cr_dw (ww (1) ) = wert 
            ELSE 
               ier_num = - 8 
               ier_typ = ER_FORT 
            ENDIF 
         ELSE 
            ier_num = - 13 
            ier_typ = ER_FORT 
            RETURN 
         ENDIF 
!                                                                       
!------ Setting lat[n]                                                  
!                                                                       
      ELSEIF (ctype.eq.'lat') THEN 
         IF (ianz.eq.1) THEN 
            IF (ww (1) .ge.1.and.ww (1) .le.3) THEN 
               cr_a0 (ww (1) ) = wert 
               CALL setup_lattice (cr_a0, cr_ar, cr_eps, cr_gten,       &
               cr_reps, cr_rten, cr_win, cr_wrez, cr_v, cr_vr, .false., &
               cr_gmat, cr_fmat, cr_cartesian,                         &
               cr_tran_g, cr_tran_gi, cr_tran_f, cr_tran_fi)
            ELSEIF (ww (1) .ge.4.and.ww (1) .le.6) THEN 
               cr_win (ww (1) - 3) = wert 
               CALL setup_lattice (cr_a0, cr_ar, cr_eps, cr_gten,       &
               cr_reps, cr_rten, cr_win, cr_wrez, cr_v, cr_vr, .false., &
               cr_gmat, cr_fmat, cr_cartesian,                          &
               cr_tran_g, cr_tran_gi, cr_tran_f, cr_tran_fi)
            ELSE 
               ier_num = - 8 
               ier_typ = ER_FORT 
            ENDIF 
         ELSE 
            ier_num = - 13 
            ier_typ = ER_FORT 
            RETURN 
         ENDIF 
      ELSEIF (ctype.eq.'occ') THEN 
         IF (ianz.eq.1) THEN 
            IF (0.lt.ww (1) .and.ww (1) .le.cr_nscat) THEN 
               cr_occ (ww (1) ) = wert 
            ELSE 
               ier_num = - 8 
               ier_typ = ER_FORT 
            ENDIF 
         ELSE 
            ier_num = - 13 
            ier_typ = ER_FORT 
            RETURN 
         ENDIF 
      ELSEIF(ctype == 'surf') THEN 
         IF(ianz == 2) THEN 
            IF (1 <= ww(1) .AND. ww(1) <= cr_natoms .AND.       &
                1 <= ww(2) .AND. ww(2) <=3              ) THEN 
               IF(ABS(NINT(wert))<100) THEN
               cr_surf(ww(2), ww(1)) = NINT(wert)
               IF(cr_surf(1,ww(1))==0 .AND. cr_surf(2,ww(1))==0 .AND. &
                  cr_surf(3,ww(1))==0                                ) THEN
                  cr_surf(0,ww(1)) = 0
                  cr_prop(ww(1)) = IBCLR(cr_prop(ww(1)), PROP_SURFACE_EXT)
               ELSE
                  IF(cr_surf(0,ww(1)) == 0 ) cr_surf(0,ww(1)) = SURF_LOCAL
                  cr_prop(ww(1)) = IBSET(cr_prop(ww(1)), PROP_SURFACE_EXT)
               ENDIF
               ELSE 
                  ier_num = -50
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = -8 
               ier_typ = ER_FORT 
            ENDIF 
         ELSE 
            ier_num = -13 
            ier_typ = ER_FORT 
            RETURN 
         ENDIF 
      ELSEIF(ctype == 'magn') THEN     ! MAGNETIC_WORK
         IF(ianz == 2) THEN 
            IF (1 <= ww(1) .AND. ww(1) <= cr_natoms .AND.       &
                1 <= ww(2) .AND. ww(2) <=3              ) THEN 
               cr_magn(ww(2), ww(1)) = wert
               IF(cr_magn(1,ww(1))==0.0 .AND. cr_magn(2,ww(1))==0.0 .AND. &
                  cr_magn(3,ww(1))==0.0                                ) THEN
                  cr_magn(0,ww(1)) = 0.00
!                 cr_prop(ww(1)) = IBCLR(cr_prop(ww(1)), PROP_SURFACE_EXT)
!              ELSE
!                 IF(cr_magn(0,ww(1)) == 0.0 ) cr_magn(0,ww(1)) = SURF_LOCAL
!                 cr_prop(ww(1)) = IBSET(cr_prop(ww(1)), PROP_SURFACE_EXT)
               ENDIF
            ELSE 
               ier_num = -8 
               ier_typ = ER_FORT 
            ENDIF 
         ELSE 
            ier_num = -13 
            ier_typ = ER_FORT 
            RETURN 
         ENDIF 
      ELSEIF (ctype.eq.'mol_dens') THEN 
         IF (ianz.eq.1) THEN 
            IF (0.le.ww (1) .and.ww (1) .le.mole_num_type) THEN 
               mole_dens (ww (1) ) = wert 
            ELSE 
               ier_num = - 8 
               ier_typ = ER_FORT 
            ENDIF 
         ELSE 
            ier_num = - 13 
            ier_typ = ER_FORT 
            RETURN 
         ENDIF 
      ELSEIF (ctype.eq.'mol_biso') THEN 
         IF (ianz.eq.1) THEN 
            IF (0.le.ww (1) .and.ww (1) .le.mole_num_type) THEN 
               mole_biso (ww (1) ) = wert 
            ELSE 
               ier_num = - 8 
               ier_typ = ER_FORT 
            ENDIF 
         ELSE 
            ier_num = - 13 
            ier_typ = ER_FORT 
            RETURN 
         ENDIF 
      ELSEIF (ctype.eq.'mol_clin') THEN 
         IF (ianz.eq.1) THEN 
            IF (0.le.ww (1) .and.ww (1) .le.mole_num_type) THEN 
               mole_clin (ww (1) ) = wert 
            ELSE 
               ier_num = - 8 
               ier_typ = ER_FORT 
            ENDIF 
         ELSE 
            ier_num = - 13 
            ier_typ = ER_FORT 
            RETURN 
         ENDIF 
      ELSEIF (ctype.eq.'mol_cqua') THEN 
         IF (ianz.eq.1) THEN 
            IF (0.le.ww (1) .and.ww (1) .le.mole_num_type) THEN 
               mole_cqua (ww (1) ) = wert 
            ELSE 
               ier_num = - 8 
               ier_typ = ER_FORT 
            ENDIF 
         ELSE 
            ier_num = - 13 
            ier_typ = ER_FORT 
            RETURN 
         ENDIF 
      ELSEIF (ctype.eq.'mol_type') THEN 
         IF (ianz.eq.1) THEN 
            IF (0 <= ww(1) .AND. ww(1) <=  mole_num_mole) THEN 
               iwert = NINT(wert)
               owert = mole_type(ww(1))
               IF(iwert>mole_num_type) THEN      ! This creates a new type
                  CALL molecule_set_array_size(mole_gene_n, mole_symm_n, mole_num_mole, iwert, mole_num_atom)
                  CALL molecule_copy_prop(owert, iwert)
               ENDIF
               mole_type(ww(1)) = iwert
            ELSE 
               ier_num = - 8 
               ier_typ = ER_FORT 
            ENDIF 
         ELSE 
            ier_num = - 13 
            ier_typ = ER_FORT 
            RETURN 
         ENDIF 
!
      ELSEIF (ctype.eq.'env') THEN 
         IF (ianz.eq.1) THEN 
            IF (0.le.ww (1) .and.ww (1) .le.MAXPAR_RES) THEN 
               atom_env (ww (1) ) = int (wert) 
            ELSE 
               ier_num = - 8 
               ier_typ = ER_FORT 
            ENDIF 
         ELSE 
            ier_num = - 13 
            ier_typ = ER_FORT 
            RETURN 
         ENDIF 
!
      ELSE 
         ier_num = - 2 
         ier_typ = ER_FORT 
         WRITE (ier_msg (1), 8000) ctype 
      ENDIF 
 8000 FORMAT    (a) 
      END SUBROUTINE discus_upd_para                       
!
!*****7***************************************************************  
!
SUBROUTINE discus_calc_intr_spec(string, line, ikl, iklz, ww, laenge, lp)
!-                                                                      
!     These are special intrinsic function for the DISCUS. Any          
!     intrinsic function that references crystallographic values        
!     is found in this subroutine.                                      
!+                                                                      
!
USE discus_config_mod 
USE crystal_mod 
USE atom_env_mod 
USE metric_mod
USE molecule_mod 
!
USE calc_expr_mod
USE errlist_mod 
USE ersetz_mod
USE ber_params_mod
USE get_params_mod
use element_data_mod, only:get_wave
USE lib_length
USE param_mod 
USE do_read_number_mod
USE precision_mod
USE str_comp_mod
use take_param_mod
use trig_degree_mod
use wink_mod
!
IMPLICIT none 
!                                                                       
CHARACTER (LEN=*), INTENT(INOUT) :: string
CHARACTER (LEN=*), INTENT(INOUT) :: line 
INTEGER,           INTENT(IN)    :: ikl
INTEGER,           INTENT(IN)    :: iklz
REAL(KIND=PREC_DP),INTENT(INOUT) :: ww
INTEGER,           INTENT(INOUT) :: laenge
INTEGER,           INTENT(INOUT) :: lp
       
!                                                                       
INTEGER, PARAMETER :: maxw = 9
!                                                                       
character(len=4) :: symbol
CHARACTER(LEN=MAX(PREC_STRING,LEN(string))) cpara (maxw) 
INTEGER lpara (maxw)
INTEGER i, j, k, ianz, lcomm, l 
integer :: lsymbol
LOGICAL :: lspace
logical :: l_energy
integer :: diff_radiation
real(kind=PREC_DP) :: energy
real(kind=PREC_DP) :: rlambda
REAL(kind=PREC_DP), dimension(MAXW) :: werte !(maxw)
real(kind=PREC_DP), dimension(3) :: u
real(kind=PREC_DP), dimension(3) :: v
real(kind=PREC_DP), dimension(3) :: w
real(kind=PREC_DP), dimension(3,3), parameter :: unitmat = reshape((/ 1.0D0, 0.0D0, 0.0D0,  &
                                                                      0.0D0, 1.0D0, 0.0D0,  &
                                                                      0.0D0, 0.0D0, 1.0D0 /), shape(unitmat))
!
integer, parameter :: NOPTIONAL = 4
integer, parameter :: O_THETA   = 1
!integer, parameter :: O_Q       = 2
!integer, parameter :: O_DSTAR   = 3
integer, parameter :: O_LAMBDA  = 4
character(LEN=   6), dimension(NOPTIONAL) :: oname   !Optional parameter names
character(LEN=PREC_STRING), dimension(NOPTIONAL) :: opara   !Optional parameter strings returned
integer            , dimension(NOPTIONAL) :: loname  !Lenght opt. para name
integer            , dimension(NOPTIONAL) :: lopara  !Lenght opt. para name returned
logical            , dimension(NOPTIONAL) :: lpresent!opt. para is present
real(kind=PREC_DP) , dimension(NOPTIONAL) :: owerte   ! Calculated values
integer, parameter                        :: ncalc = 3 ! Number of values to calculate 
!
data oname  / 'theta', 'q    ',  'dstar',  'lambda'   /
data loname /  5     ,  1     ,   5     ,   6         /
opara  =  (/ '0.0000', '0.0000', '0.0000', 'CUA1  ' /)   ! Always provide fresh default values
lopara =  (/  6,        6,        6      ,  4       /)
owerte =  (/  0.0,      0.0,      0.0    ,  1.540510/)
!
!
!                                                                       
lcomm = length_com (string, ikl) 
ier_num = - 1 
ier_typ = ER_FORT 
 werte = 0.0D0
!     DO i = 1, maxw 
!     werte (i) = 0.0 
!     ENDDO 
!                                                                       
cond_lcomm: IF (lcomm.eq.9) THEN 
!                                                                       
!     Calculate average density of a fractional coordinate              
!                                                                       
         IF (string (ikl - 9:ikl - 1) .eq.'aver_dens') THEN 
            CALL get_params (line, ianz, cpara, lpara, 6, lp) 
            IF (ier_num.eq.0) THEN 
               DO i = 1, maxw 
               werte (i) = 0.00001 
               ENDDO 
               IF (ianz.eq.3.or.ianz.eq.4) THEN 
                  DO i = 1, ianz 
                  CALL eval (cpara (i), lpara (i) ) 
                  IF (ier_num.ne.0) THEN 
                     GOTO 999 
                  ENDIF 
                  werte (i) = do_read_number (cpara (i), lpara (i) ) 
                  IF (ier_num.ne.0) THEN 
                     GOTO 999 
                  ENDIF 
                  ENDDO 
!thp               ww = inter_aver_dens(werte,maxw)                     
                  CALL ersetz2 (string, ikl, iklz, ww, 9, laenge) 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
            ENDIF 
         ELSE 
            ier_num = - 3 
            ier_typ = ER_FORT 
         ENDIF 
!                                                                       
ELSEIF (lcomm.eq.8) THEN  cond_lcomm
!                                                                       
         IF (string (ikl - 8:ikl - 1) .eq.'mol_test') THEN 
            CALL get_params (line, ianz, cpara, lpara, 6, lp) 
            IF (ier_num.eq.0) THEN 
               IF (ianz.eq.1) THEN 
                  CALL eval (cpara (1), lpara (1) ) 
                  IF (ier_num.ne.0) THEN 
                     GOTO 999 
                  ENDIF 
                  l = nint (do_read_number (cpara (1), lpara (1) ) ) 
                  IF (ier_num.ne.0) THEN 
                     GOTO 999 
                  ENDIF 
                  res_para (0) = 0 
                  res_para (1) = 0 
                  res_para (2) = 0 
                  ww = 0 
                  DO i = 1, mole_num_mole 
                  DO j = 1, mole_len (i) 
                  IF (mole_cont (mole_off (i) + j) .eq.l) THEN 
                     res_para (0) = 2 
                     res_para (1) = i 
                     res_para (2) = j 
                     ww = REAL(i) 
                     GOTO 8000 
                  ENDIF 
                  ENDDO 
                  ENDDO 
 8000             CONTINUE 
                  CALL ersetz2 (string, ikl, iklz, ww, 8, laenge) 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
            ENDIF 
         ELSE 
            ier_num = - 3 
            ier_typ = ER_FORT 
         ENDIF 
ELSEIF (lcomm.eq.7) THEN  cond_lcomm
!                                                                       
!     Calculate scalar product                                          
!                                                                       
         IF (string (ikl - 7:ikl - 1) .eq.'scalpro') THEN 
            IF (cr_v.le.0) THEN 
               ier_num = - 35 
               ier_typ = ER_APPL 
               ier_msg (1) = 'A proper unit cell must be defined' 
      ier_msg (2)  = 'for this command to operate       ' 
               RETURN 
            ENDIF 
            CALL get_params (line, ianz, cpara, lpara, 7, lp) 
            IF (ier_num.eq.0) THEN 
               IF (ianz.ge.6) THEN 
                  IF (ianz.eq.6) THEN 
                     cpara (7) = 'dd' 
                     lpara (7) = 2 
                  ENDIF 
                  DO i = 1, 6 
                  CALL eval (cpara (i), lpara (i) ) 
                  IF (ier_num.ne.0) THEN 
                     GOTO 999 
                  ENDIF 
                  werte (i) = do_read_number (cpara (i), lpara (i) ) 
                  IF (ier_num.ne.0) THEN 
                     GOTO 999 
                  ENDIF 
                  ENDDO 
                  DO i = 1, 3 
                  u (i) = werte (i) 
                  v (i) = werte (i + 3) 
                  ENDDO 
                  IF (cpara (7) .eq.'rr') THEN 
                     ww = skalpro (u, v, cr_rten) 
                  ELSEIF (cpara (7) .eq.'dd') THEN 
                     ww = skalpro (u, v, cr_gten) 
                  ELSEIF (cpara (7) .eq.'rd'.or.cpara (7) .eq.'dr')     &
                  THEN                                                  
                     ww = skalpro (u, v, unitmat) 
                  ELSE 
                     ier_num = - 81 
                     ier_typ = ER_APPL 
      ier_msg (1)  = 'The space flag must be 2 characters long' 
      ier_msg (2)  = 'and be ''dd'',''rr'',''dr'', or ''rd''  ' 
                  ENDIF 
                  CALL ersetz2 (string, ikl, iklz, ww, 7, laenge) 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
            ENDIF 
         ELSE 
            ier_num = - 3 
            ier_typ = ER_FORT 
         ENDIF 
ELSEIF (lcomm.eq.5) THEN  cond_lcomm
!                                                                       
!     Calculate reciprocal length or d-star value, respectivly          
!                                                                       
   IF (string (ikl - 5:ikl - 1) .eq.'dstar' .or.      &
       string (ikl - 5:ikl - 1) .eq.'qstar'       ) THEN 
      CALL get_params (line, ianz, cpara, lpara, 6, lp) 
      IF (ier_num /= 0) exit cond_lcomm
      call get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                        oname, loname, opara, lopara, lpresent, owerte)
      IF (ier_num /= 0) exit cond_lcomm
      CALL ber_params(ianz, cpara, lpara, werte, maxw)
      IF (ier_num /= 0) exit cond_lcomm
      IF (ianz.eq.3.or.ianz.eq.6) THEN 
         DO i = 1, 3 
            u (i) = werte (i) 
            v (i) = werte (i + 3) 
         ENDDO 
         lspace = .false. 
         ww = do_blen (lspace, u, v) 
      elseif(lpresent(O_THETA) .and. lpresent(O_LAMBDA)) then
         symbol  = opara(O_LAMBDA)(1:lopara(O_LAMBDA))
         lsymbol = lopara(O_LAMBDA)
         cpara(1) = opara(O_LAMBDA)
         lpara(1) = lopara(O_LAMBDA)
         ianz     = 1
         CALL ber_params(ianz, cpara, lpara, werte, maxw)
         IF (ier_num == 0) then
            rlambda = werte(1)
         else
            diff_radiation = 1
            l_energy = .false.
            energy = 0.0
            call get_wave ( symbol , rlambda, energy, l_energy, &
                      diff_radiation,ier_num, ier_typ )
            if(ier_num/=0) exit cond_lcomm
         endif
         ww = 2.0*sind(owerte(O_THETA))/rlambda
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
      if(ier_num==0) then
         if(string (ikl - 5:ikl - 1) .eq.'qstar') ww = ww * zpi
         CALL ersetz2 (string, ikl, iklz, ww, 5, laenge) 
      endif
   ELSE 
      ier_num = - 3 
      ier_typ = ER_FORT 
   ENDIF 
!
ELSEIF (lcomm.eq.4) THEN  cond_lcomm
!                                                                       
!     Calculate a bond angle                                            
!                                                                       
         IF (string (ikl - 4:ikl - 1) .eq.'bang') THEN 
            CALL get_params (line, ianz, cpara, lpara, 9, lp) 
            IF (ier_num.eq.0) THEN 
               IF (ianz.eq.3.or.ianz.eq.4.or.ianz.eq.6.or.ianz.eq.9)    &
               THEN                                                     
!     --------Input are atom numbers                                    
                  IF (str_comp (cpara (1) , 'atom', 1, lpara (1) , 4) ) &
                  THEN                                                  
                     DO i = 2, ianz 
                     CALL eval (cpara (i), lpara (i) ) 
                     IF (ier_num.ne.0) THEN 
                        GOTO 999 
                     ENDIF 
                     werte (i) = do_read_number (cpara (i), lpara (i) ) 
                     IF (ier_num.ne.0) THEN 
                        GOTO 999 
                     ENDIF 
                     ENDDO 
                     i = nint (werte (2) ) 
                     j = nint (werte (3) ) 
                     IF (ianz.eq.3) THEN 
      IF (i.lt.1.or.cr_natoms.lt.i.or.j.lt.1.or.cr_natoms.lt.j) THEN 
                           ier_typ = ER_APPL 
                           ier_num = - 19 
                           RETURN 
                        ELSE 
                           u (1) = cr_pos (1, i) 
                           u (2) = cr_pos (2, i) 
                           u (3) = cr_pos (3, i) 
                           v (1) = 0.0 
                           v (2) = 0.0 
                           v (3) = 0.0 
                           w (1) = cr_pos (1, j) 
                           w (2) = cr_pos (2, j) 
                           w (3) = cr_pos (3, j) 
                        ENDIF 
                     ELSE 
                        k = nint (werte (4) ) 
      IF (i.lt.1.or.cr_natoms.lt.i.or.j.lt.1.or.cr_natoms.lt.j.or.k.lt.1&
     &.or.cr_natoms.lt.k) THEN                                          
                           ier_typ = ER_APPL 
                           ier_num = - 19 
                           RETURN 
                        ELSE 
                           u (1) = cr_pos (1, i) 
                           u (2) = cr_pos (2, i) 
                           u (3) = cr_pos (3, i) 
                           v (1) = cr_pos (1, j) 
                           v (2) = cr_pos (2, j) 
                           v (3) = cr_pos (3, j) 
                           w (1) = cr_pos (1, k) 
                           w (2) = cr_pos (2, k) 
                           w (3) = cr_pos (3, k) 
                        ENDIF 
                     ENDIF 
!     --------Input are atom numbers in an environment
                  ELSEIF (str_comp(cpara(1), 'environment', 1, lpara(1), 11)) THEN
                     DO i = 2, ianz 
                     CALL eval (cpara (i), lpara (i) ) 
                     IF (ier_num.ne.0) THEN 
                        GOTO 999 
                     ENDIF 
                     werte (i) = do_read_number (cpara (i), lpara (i) ) 
                     IF (ier_num.ne.0) THEN 
                        GOTO 999 
                     ENDIF 
                     ENDDO 
                     i = nint (werte (2) ) 
                     j = nint (werte (3) ) 
                     IF (ianz.eq.3) THEN 
                        IF (i.lt.1.or.atom_env(0).lt.i.or.j.lt.1.or. &
                                      atom_env(0).lt.j) THEN 
                           ier_typ = ER_APPL 
                           ier_num = - 19 
                           RETURN 
                        ELSE 
                           u (1) = atom_pos (1, i) 
                           u (2) = atom_pos (2, i) 
                           u (3) = atom_pos (3, i) 
                           v (1) = atom_pos (1, 0)
                           v (2) = atom_pos (2, 0)
                           v (3) = atom_pos (3, 0)
                           w (1) = atom_pos (1, j) 
                           w (2) = atom_pos (2, j) 
                           w (3) = atom_pos (3, j) 
                        ENDIF 
                     ENDIF 
                  ELSE 
!     --------Input are real space coordinates                          
                     if(ianz<6) then
                        ier_num = -6
                        ier_typ = ER_COMM
                        return
                     endif
                     DO i = 1, ianz 
                     CALL eval (cpara (i), lpara (i) ) 
                     IF (ier_num.ne.0) THEN 
                        GOTO 999 
                     ENDIF 
                     werte (i) = do_read_number (cpara (i), lpara (i) ) 
                     IF (ier_num.ne.0) THEN 
                        GOTO 999 
                     ENDIF 
                     ENDDO 
                     IF (ianz.eq.6) THEN 
                        DO i = 4, 6 
                        werte (i + 3) = werte (i) 
                        werte (i) = 0.0 
                        ENDDO 
                     ENDIF 
                     DO i = 1, 3 
                     u (i) = werte (i) 
                     v (i) = werte (i + 3) 
                     w (i) = werte (i + 6) 
                     ENDDO 
                  ENDIF 
                  lspace = .true. 
                  ww = do_bang (lspace, u, v, w) 
                  CALL ersetz2 (string, ikl, iklz, ww, 4, laenge) 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
            ENDIF 
!                                                                       
!     Calculate a bond length                                           
!                                                                       
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'blen') THEN 
            CALL get_params (line, ianz, cpara, lpara, 6, lp) 
            IF (ier_num.eq.0) THEN 
               IF (ianz==2 .or. ianz.eq.3.or.ianz.eq.6) THEN 
!     --------Input are atom numbers                                    
                  IF (str_comp (cpara (1) , 'atom', 1, lpara (1) , 4) ) &
                  THEN                                                  
                     DO i = 2, ianz 
                     CALL eval (cpara (i), lpara (i) ) 
                     IF (ier_num.ne.0) THEN 
                        GOTO 999 
                     ENDIF 
                     werte (i) = do_read_number (cpara (i), lpara (i) ) 
                     IF (ier_num.ne.0) THEN 
                        GOTO 999 
                     ENDIF 
                     ENDDO 
                     i = nint (werte (2) ) 
                     j = nint (werte (3) ) 
      IF (i.lt.1.or.cr_natoms.lt.i.or.j.lt.1.or.cr_natoms.lt.j) THEN 
                        ier_typ = ER_APPL 
                        ier_num = - 19 
                        RETURN 
                     ELSE 
                        u (1) = cr_pos (1, i) 
                        u (2) = cr_pos (2, i) 
                        u (3) = cr_pos (3, i) 
                        v (1) = cr_pos (1, j) 
                        v (2) = cr_pos (2, j) 
                        v (3) = cr_pos (3, j) 
                     ENDIF 
!     --------Input are atom numbers in an environment                                   
                  ELSEIF (str_comp (cpara(1), 'environment', 1, lpara(1), 11)) THEN
                     DO i = 2, ianz 
                     CALL eval (cpara (i), lpara (i) ) 
                     IF (ier_num.ne.0) THEN 
                        GOTO 999 
                     ENDIF 
                     werte (i) = do_read_number (cpara (i), lpara (i) ) 
                     IF (ier_num.ne.0) THEN 
                        GOTO 999 
                     ENDIF 
                     ENDDO 
                     i = nint (werte (2) ) 
                     IF (i.lt.1.or.atom_env(0).lt.i) THEN
                        ier_typ = ER_APPL 
                        ier_num = - 19 
                        RETURN 
                     ELSE 
                        u (1) = atom_pos (1, i) 
                        u (2) = atom_pos (2, i) 
                        u (3) = atom_pos (3, i) 
                        v (1) = atom_pos (1, 0) 
                        v (2) = atom_pos (2, 0) 
                        v (3) = atom_pos (3, 0) 
                     ENDIF 
                  ELSE 
!     --------Input are real space coordinates                          
                     DO i = 1, ianz 
                     CALL eval (cpara (i), lpara (i) ) 
                     IF (ier_num.ne.0) THEN 
                        GOTO 999 
                     ENDIF 
                     werte (i) = do_read_number (cpara (i), lpara (i) ) 
                     IF (ier_num.ne.0) THEN 
                        GOTO 999 
                     ENDIF 
                     ENDDO 
                     IF (ianz.eq.3) THEN 
                        DO i = 4, 6 
                        werte (i) = 0.0 
                        ENDDO 
                     ENDIF 
                     DO i = 1, 3 
                     u (i) = werte (i) 
                     v (i) = werte (i + 3) 
                     ENDDO 
                  ENDIF 
                  lspace = .true. 
                  ww = do_blen (lspace, u, v) 
                  CALL ersetz2 (string, ikl, iklz, ww, 4, laenge) 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
            ENDIF 
!                                                                       
!     Calculate angle in reciprocal space                               
!                                                                       
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'rang') THEN 
            CALL get_params (line, ianz, cpara, lpara, 9, lp) 
            IF (ier_num.eq.0) THEN 
               IF (ianz.eq.6.or.ianz.eq.9) THEN 
                  DO i = 1, ianz 
                  CALL eval (cpara (i), lpara (i) ) 
                  IF (ier_num.ne.0) THEN 
                     GOTO 999 
                  ENDIF 
                  werte (i) = do_read_number (cpara (i), lpara (i) ) 
                  IF (ier_num.ne.0) THEN 
                     GOTO 999 
                  ENDIF 
                  ENDDO 
                  IF (ianz.eq.6) THEN 
                     DO i = 4, 6 
                     werte (i + 3) = werte (i) 
                     werte (i) = 0.0 
                     ENDDO 
                  ENDIF 
                  DO i = 1, 3 
                  u (i) = werte (i) 
                  v (i) = werte (i + 3) 
                  w (i) = werte (i + 6) 
                  ENDDO 
                  lspace = .false. 
                  ww = do_bang (lspace, u, v, w) 
                  CALL ersetz2 (string, ikl, iklz, ww, 4, laenge) 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
            ENDIF 
         ELSE 
            ier_num = - 3 
            ier_typ = ER_FORT 
         ENDIF 
      ELSEIF (lcomm.eq.0) THEN cond_lcomm
         CALL ersetz2 (string, ikl, iklz, ww, 0, laenge) 
      ELSE cond_lcomm
         ier_num = - 3 
         ier_typ = ER_FORT 
      ENDIF cond_lcomm
!                                                                       
  999 CONTINUE 
!                                                                       
IF (ier_num.ne.0) THEN 
   WRITE (ier_msg (1), 9000) string (1:min (40, laenge) ) 
   WRITE (ier_msg (1), 9000) line (1:min (40, laenge) ) 
ENDIF 
!                                                                       
 9000 FORMAT    (a) 
!
END SUBROUTINE discus_calc_intr_spec                 
!
!*****7**************************************************************** 
!
SUBROUTINE discus_calc_intr_log_spec(string, length)
!
! evaluate DISCUS specific logical functions
!
USE errlist_mod 
USE ersetzl_mod
USE variable_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*) , INTENT(INOUT) :: string
INTEGER          , INTENT(INOUT) :: length
!
INTEGER, PARAMETER :: N_LF = 1
!
CHARACTER (LEN=6), DIMENSION(1:N_LF) :: f_names
INTEGER          , DIMENSION(1:N_LF) :: l_names
INTEGER :: j
INTEGER :: ikl, iklz
INTEGER :: ihyp, ihyp2
INTEGER :: n_isexp
LOGICAL :: lres
!
!     LOGICAL :: is_expression
!
DATA f_names /'isprop'/
DATA l_names / 6      /
!INTERFACE
!  LOGICAL FUNCTION is_property(string)
!    CHARACTER(LEN=*), INTENT(IN) :: string
!  END FUNCTION is_property
!END INTERFACE
!
ier_num = 0
ier_typ = ER_NONE
!
DO j=1,N_LF                              ! Loop over all defined functions
   any_isexp: DO                         ! Search unitil no more found
      n_isexp = INDEX(string,f_names(j))
      IF(n_isexp > 0 ) THEN              ! We found a function
         ikl  = n_isexp + INDEX(string(n_isexp+1:length),'(')
         IF(ikl > n_isexp) THEN          ! We found an opening '('
            ihyp = ikl + INDEX(string(ikl+1:length),'''')
            IF(ihyp > ikl) THEN          ! We found an opening ''''
               ihyp2 = ihyp + INDEX(string(ihyp+1:length),'''')
               IF(ihyp2 > ihyp) THEN     ! We found an closing ''''
                  iklz = ihyp2 + INDEX(string(ihyp2+1:length-1),')')
                  IF(iklz > ihyp2) THEN  ! We found an closing ')'
                  SELECT CASE (j)
                  CASE (1)
                     lres = is_property(string(ikl+1:iklz-1))
                  END SELECT
                  CALL ersetzl (string, ikl, iklz, lres, l_names(j), length)
                  ELSE
                     ier_num = -11
                     ier_typ = ER_FORT
                     RETURN
                  ENDIF
               ELSE
                  ier_num = -38
                  ier_typ = ER_FORT
                  RETURN
               ENDIF
            ELSE
               ier_num = -38
               ier_typ = ER_FORT
               RETURN
            ENDIF
         ELSE
            ier_num = -11
            ier_typ = ER_FORT
            RETURN
         ENDIF
      ELSE
         EXIT any_isexp
      ENDIF
   ENDDO any_isexp
ENDDO
!
END SUBROUTINE discus_calc_intr_log_spec
!
!*****7**************************************************************** 
LOGICAL FUNCTION is_property(string)
!
! Evaluate the "isprop" function
!
use crystal_mod
USE prop_para_mod
!
USE errlist_mod
USE ber_params_mod
USE get_params_mod
USE precision_mod
USE take_param_mod
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(IN) :: string
!
INTEGER, PARAMETER :: MAXW      = 3
INTEGER, PARAMETER :: NOPTIONAL = 2
!
CHARACTER(LEN=MAX(PREC_STRING,LEN(string))), DIMENSION(MAXW     ) :: cpara   !parameter strings returned
INTEGER            , DIMENSION(MAXW     ) :: lpara   !Lenght para name
REAL(KIND=PREC_DP) , DIMENSION(MAXW     ) :: werte   !Lenght para name
!
CHARACTER(LEN=   3), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=MAX(PREC_STRING,LEN(string))), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent!opt. para present
REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 0 ! Number of values to calculate 
!
INTEGER :: ianz, iianz
INTEGER :: i
INTEGER :: iatom
INTEGER :: length
LOGICAL :: ltest
!
DATA oname  / 'and', 'or'/
DATA loname /  3,    2   /
!
is_property = .FALSE.
!
opara  =  (/ '      ', '      ' /)   ! Always provide fresh default values
lopara =  (/  0,        0       /)
owerte =  (/  0.0,      0.0     /)
!
length = LEN_TRIM(string)
!
CALL get_params (string, ianz, cpara, lpara, MAXW, length)
IF(ier_num /= 0) RETURN
!
opara  =  (/ '      ', '      ' /)   ! Always provide fresh default values
lopara =  (/  0,        0       /)
owerte =  (/  0.0,      0.0     /)
!
CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
IF(ier_num /=0) RETURN
iianz = 1
CALL ber_params(iianz, cpara, lpara, werte, maxw)
IF(ier_num /=0) RETURN
!
iatom = NINT(werte(1))
CALL del_params(iianz, ianz, cpara, lpara, MAXW)
!
IF(opara(2)/=' ') THEN     ! or is present
   opara(2)  = ' '
   opara(2)(1:lopara(2)-2)  = opara(2)(2:lopara(2)-1)
   lopara(2) = lopara(2) - 2
   ianz = 2
ELSE
   ianz = 1
ENDIF
IF(opara(1)==' ') THEN     ! and is not explicitly given
   opara(1)  = cpara(1)
   lopara(1) = lpara(1)
   IF(lopara(1)==0) lopara(1) = 3
ELSE
   IF(cpara(1)/=' ') THEN   ! 'and:' is given as well as default=> ERROR
      ier_num = -151
      ier_typ = ER_APPL
      RETURN 
   ENDIF
ENDIF
opara(1)  = opara(1)(2:lopara(1)-1)
lopara(1) = lopara(1) - 2
!do i=1,ianz
!  write(*,'(a, i2,a)') ' FINAL ', i, opara(i)(1:lopara(i))
!enddo
!
ltest       = .FALSE.
is_property = .TRUE.
!
IF(    opara(1)=='normal'   ) THEN
   is_property = IBITS(cr_prop(iatom),PROP_NORMAL,1)==1
   ltest = .TRUE.
ELSEIF(opara(1)=='molecule' ) THEN
   is_property = IBITS(cr_prop(iatom),PROP_MOLECULE,1)==1
   ltest = .TRUE.
ELSEIF(opara(1)=='domain'   ) THEN
   is_property = IBITS(cr_prop(iatom),PROP_DOMAIN,1)==1
   ltest = .TRUE.
ELSEIF(opara(1)=='outside'  ) THEN
   is_property = IBITS(cr_prop(iatom),PROP_OUTSIDE,1)==1
   ltest = .TRUE.
ELSEIF(opara(1)=='external' ) THEN
   is_property = IBITS(cr_prop(iatom),PROP_SURFACE_EXT,1)==1
   ltest = .TRUE.
ELSEIF(opara(1)=='internal' ) THEN
   is_property = IBITS(cr_prop(iatom),PROP_SURFACE_INT,1)==1
   ltest = .TRUE.
ELSEIF(opara(1)=='ligand'   ) THEN
   is_property = IBITS(cr_prop(iatom),PROP_LIGAND,1)==1
   ltest = .TRUE.
ELSEIF(opara(1)=='temp'     ) THEN
   is_property = IBITS(cr_prop(iatom),PROP_TEMP  ,1)==1
   ltest = .TRUE.
ELSE
!  Not a single word, test the letters
!write(*,*) ' NOT A WORD ', is_property, MAXPROP
   DO i=1, MAXPROP
!write(*,*) ' LOOP ', i , is_property
      IF(INDEX(opara(1)(1:lopara(1)), c_prop_letter(i:i))>0) THEN  ! Property letter found
         is_property=is_property .AND. IBITS(cr_prop(iatom), i-1, 1)==1  ! Property present
         ltest = .TRUE.
!write(*,*) ' TESTED ',c_prop_letter(i:i), IBITS(cr_prop(iatom), i-1,1)==1, is_property
      ELSEIF(INDEX(opara(1)(1:lopara(1)), c_prop_small (i:i))>0) THEN  ! small Property letter found
         is_property=is_property .AND. IBITS(cr_prop(iatom), i-1, 1)==0  ! property absent
         ltest = .TRUE.
!write(*,*) ' TESTED ',c_prop_small (i:i), IBITS(cr_prop(iatom), i-1,1)==1, is_property
      ENDIF
   ENDDO
ENDIF
IF(ianz==2) THEN    ! 'or:' is specified
!write(*,*) ' OR IS PRESNT '
   IF(    opara(2)=='normal'   ) THEN
      is_property = is_property .OR. IBITS(cr_prop(iatom),PROP_NORMAL,1)==1
      ltest = .TRUE.
   ELSEIF(opara(2)=='molecule' ) THEN
      is_property = is_property .OR. IBITS(cr_prop(iatom),PROP_MOLECULE,1)==1
      ltest = .TRUE.
   ELSEIF(opara(2)=='domain'   ) THEN
      is_property = is_property .OR. IBITS(cr_prop(iatom),PROP_DOMAIN,1)==1
      ltest = .TRUE.
   ELSEIF(opara(2)=='outside'  ) THEN
      is_property = is_property .OR. IBITS(cr_prop(iatom),PROP_OUTSIDE,1)==1
      ltest = .TRUE.
!write(*,*) ' TESTED OUTSIDE  ', IBITS(cr_prop(iatom),PROP_OUTSIDE,1)==1
   ELSEIF(opara(2)=='external' ) THEN
      is_property = is_property .OR. IBITS(cr_prop(iatom),PROP_SURFACE_EXT,1)==1
      ltest = .TRUE.
!write(*,*) ' TESTED EXTERNAL ', IBITS(cr_prop(iatom),PROP_SURFACE_EXT,1)==1
   ELSEIF(opara(2)=='internal' ) THEN
      is_property = is_property .OR. IBITS(cr_prop(iatom),PROP_SURFACE_INT,1)==1
      ltest = .TRUE.
   ELSEIF(opara(2)=='ligand'   ) THEN
      is_property = is_property .OR. IBITS(cr_prop(iatom),PROP_LIGAND,1)==1
      ltest = .TRUE.
   ELSEIF(opara(1)=='temp'     ) THEN
      is_property = is_property .OR. IBITS(cr_prop(iatom),PROP_TEMP,1)==1
      ltest = .TRUE.
   ELSE
!  Not a single word, test the letters
      DO i=1, MAXPROP
         IF(INDEX(opara(1)(1:lopara(1)), c_prop_letter(i:i))>0) THEN  ! Property letter found
            is_property=is_property .OR. IBITS(cr_prop(iatom), i-1, 1)==1
            ltest = .TRUE.
         ELSEIF(INDEX(opara(1)(1:lopara(1)), c_prop_small (i:i))>0) THEN  ! Property letter found
            is_property=is_property .OR. IBITS(cr_prop(iatom), i-1, 1)==0
            ltest = .TRUE.
         ENDIF
      ENDDO
   ENDIF
ENDIF
IF(.NOT.ltest) is_property = .FALSE.
!
END FUNCTION is_property
!
!*****7**************************************************************** 
!
SUBROUTINE discus_validate_var_spec (zeile, lp) 
!-                                                                      
!       checks whether the variable name is legal, DISCUS specific part 
!                                                                       
!       Version : 1.0                                                   
!       Date    : 22 Oct 03                                             
!                                                                       
!       Author  : R.B. Neder  (reinhard.neder@fau.de)    
!+                                                                      
USE reserved_mod
USE discus_config_mod 
USE errlist_mod 
!                                                                       
IMPLICIT none 
!                                                                       
CHARACTER (LEN=*), INTENT(IN) :: zeile 
INTEGER,           INTENT(IN) :: lp 
!
INTEGER  :: i , length
!                                                                       
ier_num = 0 
ier_typ = ER_NONE 
!                                                                       
main: DO i = 1, discus_reserved_n 
!  IF (index (discus_reserved (i), zeile (1:lp) ) .ne.0) THEN 
   length = MAX(LEN_TRIM(discus_reserved(i)), LEN_TRIM(zeile(1:lp)))
   length = MIN(length, LEN(discus_reserved), LEN(zeile))
   IF(discus_reserved (i)(1:length)== zeile(1:length) ) THEN    
      ier_num = - 25 
      ier_typ = ER_FORT 
      EXIT main
   ENDIF 
ENDDO main
!                                                                       
END SUBROUTINE discus_validate_var_spec              
!
!*******************************************************************************
!
SUBROUTINE discus_get_var_type(line,length, var_is_type)
!
! Returns the variable type : INTEGER, REAL, CHARACTER, and Scalar versus field
!
USE constants_mod
USE lib_get_var
USE variable_mod
!
IMPLICIT NONE
!
CHARACTER(LEN=*)     , INTENT(IN)  :: line
INTEGER              , INTENT(IN)  :: length
INTEGER, DIMENSION(3), INTENT(OUT) :: var_is_type
!
INTEGER, PARAMETER :: MAXPAR = 28
CHARACTER(LEN=16), DIMENSION(MAXPAR) :: discus_names
INTEGER          , DIMENSION(MAXPAR) :: discus_type
INTEGER          , DIMENSION(MAXPAR) :: discus_dim
LOGICAL          , DIMENSION(MAXPAR) :: discus_ro 
INTEGER :: i
!
DATA discus_names  &
    /'pdf_scal', 'pdf_dens', 'mol_type', 'mol_dens', 'mol_cont', &
     'mol_cqua', 'mol_clin',                                     &
     'mol_biso', 'mol_len ', 'in_mole ', 'at_type ', 'at_name ', &
     'sym_n   ', 'rvol    ', 'menv    ', 'magn'    , 'cdim    ', 'surf    ', 'vol     ', &
     'occ     ', 'lat     ', 'env     ', 'z       ', 'y       ', &
     'x       ', 'n       ', 'm       ', 'b       '              &
    /
DATA discus_type &
    /  IS_REAL ,   IS_REAL ,   IS_INTE ,   IS_REAL ,   IS_INTE , &
       IS_REAL ,   IS_REAL ,                                     &
       IS_REAL ,   IS_INTE ,   IS_INTE ,   IS_CHAR ,   IS_CHAR , &
       IS_INTE ,   IS_REAL ,   IS_INTE ,   IS_REAL ,   IS_REAL ,   IS_INTE , IS_REAL , &
       IS_REAL ,   IS_REAL ,   IS_INTE ,   IS_REAL ,   IS_REAL , &
       IS_REAL ,   IS_INTE ,   IS_INTE ,   IS_REAL               &
    /
DATA discus_dim  &
    /  IS_VEC  ,   IS_VEC  ,   IS_VEC  ,   IS_VEC  ,   IS_ARR  , &
       IS_VEC  ,   IS_VEC  ,                                     &
       IS_VEC  ,   IS_VEC  ,   IS_VEC  ,   IS_VEC  ,   IS_VEC  , &
       IS_VEC  ,   IS_VEC  ,   IS_VEC  ,   IS_ARR  ,   IS_ARR  ,   IS_ARR  , IS_VEC  , &
       IS_VEC  ,   IS_VEC  ,   IS_VEC  ,   IS_VEC  ,   IS_VEC  , &
       IS_VEC  ,   IS_VEC  ,   IS_VEC  ,   IS_VEC                &
    /
DATA discus_ro  &
    /  .FALSE. ,   .FALSE. ,   .FALSE. ,   .FALSE. ,   .TRUE.  , &
       .FALSE. ,   .FALSE. ,                                     &
       .FALSE. ,   .TRUE.  ,   .TRUE.  ,   .TRUE.  ,   .TRUE.  , &
       .TRUE.  ,   .TRUE.  ,   .TRUE.  ,   .FALSE. ,   .TRUE.  ,   .FALSE. , .TRUE.  , &
       .FALSE. ,   .FALSE. ,   .TRUE.  ,   .FALSE. ,   .FALSE. , &
       .FALSE. ,   .TRUE.  ,   .FALSE. ,   .FALSE.               &
    /
!
var_is_type(:) = IS_UNKNOWN
!
main: DO i=1, MAXPAR
   IF(line(1:length) == discus_names(i)(1:LEN_TRIM(discus_names(i)))) THEN
      var_is_type(1) = discus_type(i)
      var_is_type(2) = discus_dim (i)
      IF(discus_ro(i)) THEN
         var_is_type(3) = IS_READ
      ELSE
         var_is_type(3) = IS_WRITE
      ENDIF
      RETURN
   ENDIF
ENDDO main
!
CALL lib_get_var_type(line, length, var_is_type)
!
!
END SUBROUTINE discus_get_var_type
!
end module discus_update_mod
