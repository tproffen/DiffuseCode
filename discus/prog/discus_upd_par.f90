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
      USE errlist_mod 
      USE param_mod 
      USE random_mod 
      IMPLICIT none 
!                                                                       
      INTEGER,                    INTENT(IN   ) :: ikl
      INTEGER,                    INTENT(IN   ) :: iklz
      CHARACTER (LEN=*),          INTENT(OUT  ) :: string 
      INTEGER,                    INTENT(OUT  ) :: ll
      INTEGER,                    INTENT(IN   ) :: maxw
      REAL   , DIMENSION(1:maxw), INTENT(IN   ) :: ww
      INTEGER,                    INTENT(IN   ) :: ianz
!                                                                       
      CHARACTER(1024) zeile 
!                                                                       
      INTEGER laenge, ltyp, kpara, kpara2
      INTEGER lcomm 
      INTEGER length_com 
!                                                                       
      laenge = ll 
      ltyp = 1 
      zeile = ' ' 
      kpara = nint (ww (1) ) 
      kpara2 = 0
      IF (maxw.ge.2) then 
         kpara2 = nint (ww (2) ) 
      ENDIF 
!                                                                       
      lcomm = length_com (string, ikl) 
!                                                                       
      IF (lcomm.eq.1) then 
!                                                                       
         IF (ikl.gt.lcomm + 1) zeile (1:ikl - lcomm - 1) = string (1:   &
         ikl - lcomm - 1)                                               
         IF (string (ikl - 1:ikl - 1) .eq.'i') then 
            IF (ianz.eq.1) then 
               IF (0.le.kpara.and.kpara.le.MAXPAR) then 
                  WRITE (zeile (ikl - 1:ikl + 13) , '(i15)') inpara (   &
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
         ELSEIF (string (ikl - 1:ikl - 1) .eq.'r') then 
            IF (ianz.eq.1) then 
               IF (0.le.kpara.and.kpara.le.MAXPAR) then 
                  WRITE (zeile (ikl - 1:ikl + 13) , '(e15.8e2)') rpara (&
                  kpara)                                                
                  zeile (ikl + 10:ikl + 10) = 'e' 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
         ELSEIF (string (ikl - 1:ikl - 1) .eq.'x') then 
            IF (ianz.eq.1) then 
               IF (1.le.kpara.and.kpara.le.NMAX.and.kpara.le.cr_natoms) &
               then                                                     
                  WRITE (zeile (ikl - 1:ikl + 13) , '(e15.8e2)') cr_pos &
                  (1, kpara)                                            
                  zeile (ikl + 10:ikl + 10) = 'e' 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
         ELSEIF (string (ikl - 1:ikl - 1) .eq.'y') then 
            IF (ianz.eq.1) then 
               IF (1.le.kpara.and.kpara.le.NMAX.and.kpara.le.cr_natoms) &
               then                                                     
                  WRITE (zeile (ikl - 1:ikl + 13) , '(e15.8e2)') cr_pos &
                  (2, kpara)                                            
                  zeile (ikl + 10:ikl + 10) = 'e' 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
         ELSEIF (string (ikl - 1:ikl - 1) .eq.'z') then 
            IF (ianz.eq.1) then 
               IF (1.le.kpara.and.kpara.le.NMAX.and.kpara.le.cr_natoms) &
               then                                                     
                  WRITE (zeile (ikl - 1:ikl + 13) , '(e15.8e2)') cr_pos &
                  (3, kpara)                                            
                  zeile (ikl + 10:ikl + 10) = 'e' 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
         ELSEIF (string (ikl - 1:ikl - 1) .eq.'m') then 
            IF (ianz.eq.1) then 
               IF (1.le.kpara.and.kpara.le.NMAX.and.kpara.le.cr_natoms) &
               then                                                     
                  WRITE (zeile (ikl - 1:ikl + 13) , '(i15)') cr_iscat ( &
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
         ELSEIF (string (ikl - 1:ikl - 1) .eq.'b') then 
            IF (ianz.eq.1) then 
               IF (1.le.kpara.and.kpara.le.cr_nscat) then 
                  WRITE (zeile (ikl - 1:ikl + 13) , '(e15.8e2)') cr_dw (&
                  kpara)                                                
                  zeile (ikl + 10:ikl + 10) = 'e' 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
         ELSEIF (string (ikl - 1:ikl - 1) .eq.'n') then 
            IF (ianz.eq.1) then 
               IF (kpara.eq.1) then 
                  WRITE (zeile (ikl - 1:ikl + 13) , '(i15)') cr_natoms 
               ELSEIF (kpara.eq.2) then 
                  WRITE (zeile (ikl - 1:ikl + 13) , '(i15)') cr_nscat 
               ELSEIF (kpara.eq.3) then 
                  WRITE (zeile (ikl - 1:ikl + 13) , '(i15)') cr_ncatoms 
               ELSEIF (kpara.eq.4) then 
                  WRITE (zeile (ikl - 1:ikl + 13) , '(i15)')            &
                  mole_num_mole                                         
               ELSEIF (kpara.eq.5) then 
                  WRITE (zeile (ikl - 1:ikl + 13) , '(i15)')            &
                  mole_num_type                                         
               ELSEIF (kpara.eq.6) then 
                  WRITE (zeile (ikl - 1:ikl + 13) , '(i15)')            &
                  mole_num_unit                                         
               ELSEIF (kpara.eq.7) then 
                  WRITE (zeile (ikl - 1:ikl + 13) , '(i15)')            &
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
      ELSEIF (lcomm.eq.8) then 
                                                                        
         IF (string (ikl - 8:ikl - 1) .eq.'mol_cont') then 
            IF (ianz.eq.2) then 
               IF (ikl.gt.lcomm + 1) zeile (1:ikl - lcomm - 1) = string &
               (1:ikl - lcomm - 1)                                      
               IF (0.lt.kpara.and.kpara.le.mole_num_mole) then 
                  IF (0.eq.kpara2) then 
                     WRITE (zeile (ikl - 8:ikl + 13) , '(i15)')         &
                     mole_len (kpara)                                   
                  ELSEIF (0.le.kpara2.and.kpara2.le.mole_len (kpara) )  &
                  then                                                  
                     WRITE (zeile (ikl - 8:ikl + 13) , '(i15)')         &
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
         ELSEIF (string (ikl - 8:ikl - 1) .eq.'mol_biso') then 
            IF (ianz.eq.1) then 
               IF (ikl.gt.lcomm + 1) zeile (1:ikl - lcomm - 1) = string &
               (1:ikl - lcomm - 1)                                      
               IF (0.lt.kpara.and.kpara.le.mole_num_type) then 
                  WRITE (zeile (ikl - 8:ikl + 13) , '(e15.8e2)')        &
                  mole_biso (kpara)                                     
                  zeile (ikl + 3:ikl + 3) = 'e' 
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
         ELSEIF (string (ikl - 8:ikl - 1) .eq.'mol_dens') then 
            IF (ianz.eq.1) then 
               IF (ikl.gt.lcomm + 1) zeile (1:ikl - lcomm - 1) = string &
               (1:ikl - lcomm - 1)                                      
               IF (0.lt.kpara.and.kpara.le.mole_num_type) then 
                  WRITE (zeile (ikl - 8:ikl + 13) , '(e15.8e2)')        &
                  mole_dens (kpara)                                     
                  zeile (ikl + 3:ikl + 3) = 'e' 
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
         ELSEIF (string (ikl - 8:ikl - 1) .eq.'mol_type') then 
            IF (ianz.eq.1) then 
               IF (ikl.gt.lcomm + 1) zeile (1:ikl - lcomm - 1) = string &
               (1:ikl - lcomm - 1)                                      
               IF (0.lt.kpara.and.kpara.le.mole_num_mole) then 
                  WRITE (zeile (ikl - 8:ikl + 13) , '(i15)') mole_type (&
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
         ELSEIF (string (ikl - 8:ikl - 1) .eq.'pdf_dens') then 
            IF (ianz.eq.1) then 
               IF (ikl.gt.lcomm + 1) zeile (1:ikl - lcomm - 1) = string &
               (1:ikl - lcomm - 1)                                      
               WRITE (zeile (ikl - 8:ikl + 13) , '(e15.8e2)') pdf_rho0
               zeile (ikl + 3:ikl + 3) = 'e' 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
!                                                                       
         ELSEIF (string (ikl - 8:ikl - 1) .eq.'pdf_scal') then 
            IF (ianz.eq.1) then 
               IF (ikl.gt.lcomm + 1) zeile (1:ikl - lcomm - 1) = string &
               (1:ikl - lcomm - 1)                                      
               WRITE (zeile (ikl - 8:ikl + 13) , '(e15.8e2)') pdf_scale
               zeile (ikl + 3:ikl + 3) = 'e' 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
!                                                                       
         ELSEIF (string (ikl - 8:ikl - 1) .eq.'ref_para') then 
            IF (ianz.eq.1) then 
               IF (ikl.gt.lcomm + 1) zeile (1:ikl - lcomm - 1) = string &
               (1:ikl - lcomm - 1)                                      
               IF (0.lt.kpara.and.kpara.le.MAXPAR_REF   ) then 
                  WRITE (zeile (ikl - 8:ikl + 13) , '(e15.8e2)')        &
                  ref_para (kpara)                                     
                  zeile (ikl + 3:ikl + 3) = 'e' 
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
      ELSEIF (lcomm.eq.7) then 
!                                                                       
         IF (string (ikl - 7:ikl - 1) .eq.'mol_len') then 
            IF (ianz.eq.1) then 
               IF (ikl.gt.lcomm + 1) zeile (1:ikl - lcomm - 1) = string &
               (1:ikl - lcomm - 1)                                      
               IF (0.lt.kpara.and.kpara.le.mole_num_mole) then 
                  WRITE (zeile (ikl - 7:ikl + 13) , '(i15)') mole_len ( &
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
         ELSEIF (string (ikl - 7:ikl - 1) .eq.'at_name') then 
            IF (ianz.eq.1) then 
               IF (ikl.gt.lcomm + 1) zeile (1:ikl - lcomm - 1) = string &
               (1:ikl - lcomm - 1)                                      
               IF (1.le.kpara.and.kpara.le.NMAX.and.kpara.le.cr_natoms) &
               then                                                     
                  zeile (ikl - 7:ikl + 13) = cr_at_lis (cr_iscat (kpara)&
                  )                                                     
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
         ELSEIF (string (ikl - 7:ikl - 1) .eq.'at_type') then 
            IF (ianz.eq.1) then 
               IF (ikl.gt.lcomm + 1) zeile (1:ikl - lcomm - 1) = string &
               (1:ikl - lcomm - 1)                                      
               IF (0.le.kpara.and.kpara.le.NMAX.and.kpara.le.cr_nscat) THEN
                  zeile (ikl - 7:ikl + 13) = cr_at_lis (kpara)
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
         ELSEIF (string (ikl - 7:ikl - 1) .eq.'in_mole') then 
            IF (ianz.eq.1) then 
               IF (ikl.gt.lcomm + 1) zeile (1:ikl - lcomm - 1) = string &
               (1:ikl - lcomm - 1)                                      
               IF (1.le.kpara.and.kpara.le.NMAX.and.kpara.le.cr_natoms) &
               THEN                                                     
                  WRITE(zeile(ikl-7:ikl + 13),'(i15)') cr_mole(kpara)
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
      ELSEIF(lcomm.eq.5) then 
         IF(string(ikl - 5:ikl - 1) .eq.'sym_n') then 
            IF(ianz.eq.1) then 
               WRITE(zeile(ikl-5:ikl+13),'(i15)') spc_n
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
      ELSEIF (lcomm.eq.4) then 
!                                                                       
         IF (string (ikl - 4:ikl - 1) .eq.'cdim') then 
            IF (ianz.eq.2) then 
               IF (ikl.gt.lcomm + 1) zeile (1:ikl - lcomm - 1) = string &
               (1:ikl - lcomm - 1)                                      
      IF (1.le.kpara.and.kpara.le.3.and.1.le.kpara2.and.kpara2.le.2) the&
     &n                                                                 
                  WRITE (zeile (ikl - 4:ikl + 13) , '(e15.8e2)') cr_dim &
                  (kpara, kpara2)                                       
                  zeile (ikl + 7:ikl + 7) = 'e' 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'prop') then 
            IF (ianz.eq.1) then 
               IF (1.le.kpara.and.kpara.le.NMAX.and.kpara.le.cr_natoms) &
               then                                                     
                  WRITE (zeile (ikl - 4:ikl + 13) , '(i15)') cr_prop (  &
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
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'menv') then 
            IF (ianz.eq.1) then 
               IF (ikl.gt.lcomm + 1) zeile (1:ikl - lcomm - 1) = string &
               (1:ikl - lcomm - 1)                                      
               IF (0.le.kpara.and.kpara.le.MAXPAR_RES) then 
      WRITE (zeile (ikl - 4:ikl + 13) , '(i15    )') mole_env (kpara) 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'seed') then 
            IF (ianz.eq.1) then 
               IF (ikl.gt.lcomm + 1) zeile (1:ikl - lcomm - 1) = string &
               (1:ikl - lcomm - 1)                                      
      WRITE (zeile (ikl - 4:ikl + 13) , '(i15    )') idum
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'rvol') then 
            IF (ianz.eq.1) then 
               IF (ikl.gt.lcomm + 1) zeile (1:ikl - lcomm - 1) = string &
               (1:ikl - lcomm - 1)                                      
               WRITE (zeile (ikl - 4:ikl + 13) , '(e15.8e2)') cr_vr 
               zeile (ikl + 7:ikl + 7) = 'e' 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
         ELSE 
            ier_num = - 2 
            ier_typ = ER_FORT 
         ENDIF 
      ELSEIF (lcomm.eq.3) then 
!                                                                       
         IF (string (ikl - 3:ikl - 1) .eq.'res') then 
            IF (ianz.eq.1) then 
               IF (ikl.gt.lcomm + 1) zeile (1:ikl - lcomm - 1) = string &
               (1:ikl - lcomm - 1)                                      
               IF (0.le.kpara.and.kpara.le.MAXPAR_RES) then 
                  WRITE (zeile (ikl - 3:ikl + 13) , '(e15.8e2)')        &
                  res_para (kpara)                                      
                  zeile (ikl + 8:ikl + 8) = 'e' 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
         ELSEIF (string (ikl - 3:ikl - 1) .eq.'env') then 
            IF (ianz.eq.1) then 
               IF (ikl.gt.lcomm + 1) zeile (1:ikl - lcomm - 1) = string &
               (1:ikl - lcomm - 1)                                      
               IF (0.le.kpara.and.kpara.le.MAXPAR_RES) then 
      WRITE (zeile (ikl - 3:ikl + 13) , '(i15    )') atom_env (kpara) 
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
         ELSEIF (string (ikl - 3:ikl - 1) .eq.'lat') then 
            IF (ianz.eq.1) then 
               IF (ikl.gt.lcomm + 1) zeile (1:ikl - lcomm - 1) = string &
               (1:ikl - lcomm - 1)                                      
               IF (kpara.ge.1.and.kpara.le.3) then 
                  WRITE (zeile (ikl - 3:ikl + 13) , '(e15.8e2)') cr_a0 (&
                  kpara)                                                
                  zeile (ikl + 8:ikl + 8) = 'e' 
               ELSEIF (kpara.ge.4.and.kpara.le.6) then 
                  WRITE (zeile (ikl - 3:ikl + 13) , '(e15.8e2)') cr_win &
                  (kpara - 3)                                           
                  zeile (ikl + 8:ikl + 8) = 'e' 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_FORT 
               ENDIF 
            ELSE 
               ier_num = - 13 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
         ELSEIF (string (ikl - 3:ikl - 1) .eq.'vol') then 
            IF (ianz.eq.1) then 
               IF (ikl.gt.lcomm + 1) zeile (1:ikl - lcomm - 1) = string &
               (1:ikl - lcomm - 1)                                      
               WRITE (zeile (ikl - 3:ikl + 13) , '(e15.8e2)') cr_v 
               zeile (ikl + 8:ikl + 8) = 'e' 
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
      IF (ier_num.eq.0) then 
         ll = laenge+15 - ltyp - (iklz - ikl + 1) 
         IF (iklz + 1.le.laenge) zeile (ikl + 14:ll) = string (iklz + 1:&
         laenge)                                                        
         string = zeile 
      ELSE 
         ll = min (40, laenge) 
         WRITE (ier_msg (1), 8000) string (1:ll) 
      ENDIF 
      CALL rem_bl (string, ll) 
!                                                                       
 8000 FORMAT    (a) 
      END SUBROUTINE discus_ersetz_para                    
!*****7*****************************************************************
      SUBROUTINE discus_upd_para (ctype, ww, maxw, wert, ianz) 
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
      USE errlist_mod 
      USE param_mod 
      IMPLICIT none 
!                                                                       
      CHARACTER (LEN=*),          INTENT(IN) :: ctype 
      INTEGER,                    INTENT(IN) :: maxw
      INTEGER,                    INTENT(IN) :: ianz 
      INTEGER, DIMENSION(1:MAXW), INTENT(IN) :: ww
      REAL   ,                    INTENT(IN) :: wert 
!
      INTEGER :: l
      INTEGER :: iwert, owert
!                                                                       
      IF (ctype.eq.'x') then 
         IF (ianz.eq.1) then 
            IF (0.lt.ww (1) .and.ww (1) .le.NMAX.and.ww (1)             &
            .le.cr_natoms) then                                         
               cr_pos (1, ww (1) ) = wert 
               l = 1 
               cr_dim (l, 1) = amin1 (cr_dim (l, 1), cr_pos (l, ww (1) ))
               cr_dim (l, 2) = amax1 (cr_dim (l, 2), cr_pos (l, ww (1) ))
            ELSE 
               ier_num = - 8 
               ier_typ = ER_FORT 
            ENDIF 
         ELSE 
            ier_num = - 13 
            ier_typ = ER_FORT 
            RETURN 
         ENDIF 
      ELSEIF (ctype.eq.'y') then 
         IF (ianz.eq.1) then 
            IF (0.lt.ww (1) .and.ww (1) .le.NMAX.and.ww (1)             &
            .le.cr_natoms) then                                         
               cr_pos (2, ww (1) ) = wert 
               l = 2 
               cr_dim (l, 1) = amin1 (cr_dim (l, 1), cr_pos (l, ww (1) ))
               cr_dim (l, 2) = amax1 (cr_dim (l, 2), cr_pos (l, ww (1) ))
            ELSE 
               ier_num = - 8 
               ier_typ = ER_FORT 
            ENDIF 
         ELSE 
            ier_num = - 13 
            ier_typ = ER_FORT 
            RETURN 
         ENDIF 
      ELSEIF (ctype.eq.'z') then 
         IF (ianz.eq.1) then 
            IF (0.lt.ww (1) .and.ww (1) .le.NMAX.and.ww (1)             &
            .le.cr_natoms) then                                         
               cr_pos (3, ww (1) ) = wert 
               l = 3 
               cr_dim (l, 1) = amin1 (cr_dim (l, 1), cr_pos (l, ww (1)))
               cr_dim (l, 2) = amax1 (cr_dim (l, 2), cr_pos (l, ww (1)))
            ELSE 
               ier_num = - 8 
               ier_typ = ER_FORT 
            ENDIF 
         ELSE 
            ier_num = - 13 
            ier_typ = ER_FORT 
            RETURN 
         ENDIF 
      ELSEIF (ctype.eq.'m') then 
         IF (ianz.eq.1) then 
            IF (0.lt.ww (1) .and.ww (1) .le.NMAX.and.ww (1)             &
            .le.cr_natoms) then                                         
               IF (0.le.int (wert) .and.int (wert) .le.MAXSCAT.and.int (&
               wert) .le.cr_nscat) then                                 
                  cr_iscat (ww (1) ) = int (wert) 
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
      ELSEIF (ctype.eq.'b') then 
         IF (ianz.eq.1) then 
            IF (0.lt.ww (1) .and.ww (1) .le.cr_nscat) then 
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
      ELSEIF (ctype.eq.'i') then 
         IF (ianz.eq.1) then 
            IF (0.le.ww (1) .and.ww (1) .le.MAXPAR) then 
               inpara (ww (1) ) = int (wert) 
            ELSE 
               ier_num = - 8 
               ier_typ = ER_FORT 
            ENDIF 
         ELSE 
            ier_num = - 13 
            ier_typ = ER_FORT 
            RETURN 
         ENDIF 
      ELSEIF (ctype.eq.'r') then 
         IF (ianz.eq.1) then 
            IF (0.le.ww (1) .and.ww (1) .le.MAXPAR) then 
               rpara (ww (1) ) = wert 
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
!DBG        ELSEIF(ctype.eq.'property') then                            
!DBG        if(ianz.eq.1) then                                          
!DBG          if(0.lt.ww(1) .and. ww(1).le.NMAX  .and.                  
!DBG     &       ww(1).le.cr_natoms   ) then                            
!DBG      if(0 .le.int(wert) .and. int(wert)-1.le.2**(MAXPROP-1)) then  
!DBG             cr_prop (ww(1))=int(wert)                              
!DBG            ELSE                                                    
!DBG              ier_num = -102                                        
!DBG              ier_typ = ER_APPL                                     
!DBG            endif                                                   
!DBG          ELSE                                                      
!DBG            ier_num = -8                                            
!DBG            ier_typ = ER_FORT                                       
!DBG          endif                                                     
!DBG        ELSE                                                        
!DBG          ier_num = -13                                             
!DBG          ier_typ = ER_FORT                                         
!DBG          return                                                    
!DBG        endif                                                       
!                                                                       
!------ Setting lat[n]                                                  
!                                                                       
      ELSEIF (ctype.eq.'lat') then 
         IF (ianz.eq.1) then 
            IF (ww (1) .ge.1.and.ww (1) .le.3) then 
               cr_a0 (ww (1) ) = wert 
               CALL setup_lattice (cr_a0, cr_ar, cr_eps, cr_gten,       &
               cr_reps, cr_rten, cr_win, cr_wrez, cr_v, cr_vr, .false., &
               cr_gmat, cr_fmat, cr_cartesian,                         &
               cr_tran_g, cr_tran_gi, cr_tran_f, cr_tran_fi)
            ELSEIF (ww (1) .ge.4.and.ww (1) .le.6) then 
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
      ELSEIF (ctype.eq.'mol_dens') then 
         IF (ianz.eq.1) then 
            IF (0.le.ww (1) .and.ww (1) .le.mole_num_type) then 
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
      ELSEIF (ctype.eq.'mol_biso') then 
         IF (ianz.eq.1) then 
            IF (0.le.ww (1) .and.ww (1) .le.mole_num_type) then 
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
      ELSEIF (ctype.eq.'mol_type') then 
         IF (ianz.eq.1) then 
            IF (0 <= ww(1) .AND. ww(1) <=  mole_num_mole) then 
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
      ELSEIF (ctype.eq.'res') then 
         IF (ianz.eq.1) then 
            IF (0.le.ww (1) .and.ww (1) .le.MAXPAR_RES) then 
               res_para (ww (1) ) = wert 
            ELSE 
               ier_num = - 8 
               ier_typ = ER_FORT 
            ENDIF 
         ELSE 
            ier_num = - 13 
            ier_typ = ER_FORT 
            RETURN 
         ENDIF 
      ELSEIF (ctype.eq.'env') then 
         IF (ianz.eq.1) then 
            IF (0.le.ww (1) .and.ww (1) .le.MAXPAR_RES) then 
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
      ELSEIF (ctype.eq.'ref_para') THEN
         IF (ianz.eq.1) then 
            IF (0.le.ww (1) .and.ww (1) .le.MAXPAR_REF) then 
               ref_para (ww (1) ) = wert 
            ELSE 
               ier_num = - 8 
               ier_typ = ER_FORT 
            ENDIF 
         ELSE 
            ier_num = - 13 
            ier_typ = ER_FORT 
            RETURN 
         ENDIF 
      ELSE 
         ier_num = - 2 
         ier_typ = ER_FORT 
         WRITE (ier_msg (1), 8000) ctype 
      ENDIF 
 8000 FORMAT    (a) 
      END SUBROUTINE discus_upd_para                       
!*****7***************************************************************  
      SUBROUTINE discus_calc_intr_spec (string, line, ikl, iklz, ww, laenge,   &
      lp)                                                               
!-                                                                      
!     These are special intrinsic function for the DISCUS. Any          
!     intrinsic function that references crystallographic values        
!     is found in this subroutine.                                      
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE atom_env_mod 
      USE metric_mod
      USE molecule_mod 
      USE errlist_mod 
      USE param_mod 
      IMPLICIT none 
!                                                                       
      CHARACTER (LEN=*), INTENT(IN   ) :: string
      CHARACTER (LEN=*), INTENT(INOUT) :: line 
      INTEGER,           INTENT(IN)    :: ikl
      INTEGER,           INTENT(IN)    :: iklz
      INTEGER,           INTENT(IN)    :: laenge
      INTEGER,           INTENT(IN)    :: lp
      REAL   ,           INTENT(OUT)   :: ww
       
!                                                                       
      INTEGER, PARAMETER :: maxw = 9
!                                                                       
      CHARACTER(1024) cpara (maxw) 
      INTEGER lpara (maxw)
      INTEGER i, j, k, ianz, lcomm, l 
      LOGICAL lspace
      REAL werte (maxw), u (3), v (3), w (3) 
      REAL unitmat (3, 3) 
!                                                                       
      INTEGER length_com 
      LOGICAL str_comp 
!     REAL do_blen, do_bang 
!     REAL skalpro 
      REAL do_read_number 
!                                                                       
      DATA unitmat / 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 / 
!                                                                       
      lcomm = length_com (string, ikl) 
      ier_num = - 1 
      ier_typ = ER_FORT 
      DO i = 1, maxw 
      werte (i) = 0.0 
      ENDDO 
!                                                                       
      IF (lcomm.eq.9) then 
!                                                                       
!     Calculate average density of a fractional coordinate              
!                                                                       
         IF (string (ikl - 9:ikl - 1) .eq.'aver_dens') then 
            CALL get_params (line, ianz, cpara, lpara, 6, lp) 
            IF (ier_num.eq.0) then 
               DO i = 1, maxw 
               werte (i) = 0.00001 
               ENDDO 
               IF (ianz.eq.3.or.ianz.eq.4) then 
                  DO i = 1, ianz 
                  CALL eval (cpara (i), lpara (i) ) 
                  IF (ier_num.ne.0) then 
                     GOTO 999 
                  ENDIF 
                  werte (i) = do_read_number (cpara (i), lpara (i) ) 
                  IF (ier_num.ne.0) then 
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
      ELSEIF (lcomm.eq.8) then 
!                                                                       
         IF (string (ikl - 8:ikl - 1) .eq.'mol_test') then 
            CALL get_params (line, ianz, cpara, lpara, 6, lp) 
            IF (ier_num.eq.0) then 
               IF (ianz.eq.1) then 
                  CALL eval (cpara (1), lpara (1) ) 
                  IF (ier_num.ne.0) then 
                     GOTO 999 
                  ENDIF 
                  l = nint (do_read_number (cpara (1), lpara (1) ) ) 
                  IF (ier_num.ne.0) then 
                     GOTO 999 
                  ENDIF 
                  res_para (0) = 0 
                  res_para (1) = 0 
                  res_para (2) = 0 
                  ww = 0 
                  DO i = 1, mole_num_mole 
                  DO j = 1, mole_len (i) 
                  IF (mole_cont (mole_off (i) + j) .eq.l) then 
                     res_para (0) = 2 
                     res_para (1) = i 
                     res_para (2) = j 
                     ww = float (i) 
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
      ELSEIF (lcomm.eq.7) then 
!                                                                       
!     Calculate scalar product                                          
!                                                                       
         IF (string (ikl - 7:ikl - 1) .eq.'scalpro') then 
            IF (cr_v.le.0) then 
               ier_num = - 35 
               ier_typ = ER_APPL 
               ier_msg (1) = 'A proper unit cell must be defined' 
      ier_msg (2)  = 'for this command to operate       ' 
               RETURN 
            ENDIF 
            CALL get_params (line, ianz, cpara, lpara, 7, lp) 
            IF (ier_num.eq.0) then 
               IF (ianz.ge.6) then 
                  IF (ianz.eq.6) then 
                     cpara (7) = 'dd' 
                     lpara (7) = 2 
                  ENDIF 
                  DO i = 1, 6 
                  CALL eval (cpara (i), lpara (i) ) 
                  IF (ier_num.ne.0) then 
                     GOTO 999 
                  ENDIF 
                  werte (i) = do_read_number (cpara (i), lpara (i) ) 
                  IF (ier_num.ne.0) then 
                     GOTO 999 
                  ENDIF 
                  ENDDO 
                  DO i = 1, 3 
                  u (i) = werte (i) 
                  v (i) = werte (i + 3) 
                  ENDDO 
                  IF (cpara (7) .eq.'rr') then 
                     ww = skalpro (u, v, cr_rten) 
                  ELSEIF (cpara (7) .eq.'dd') then 
                     ww = skalpro (u, v, cr_gten) 
                  ELSEIF (cpara (7) .eq.'rd'.or.cpara (7) .eq.'dr')     &
                  then                                                  
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
      ELSEIF (lcomm.eq.5) then 
!                                                                       
!     Calculate reciprocal length or d-star value, respectivly          
!                                                                       
         IF (string (ikl - 5:ikl - 1) .eq.'dstar') then 
            CALL get_params (line, ianz, cpara, lpara, 6, lp) 
            IF (ier_num.eq.0) then 
               IF (ianz.eq.3.or.ianz.eq.6) then 
                  DO i = 1, ianz 
                  CALL eval (cpara (i), lpara (i) ) 
                  IF (ier_num.ne.0) then 
                     GOTO 999 
                  ENDIF 
                  werte (i) = do_read_number (cpara (i), lpara (i) ) 
                  IF (ier_num.ne.0) then 
                     GOTO 999 
                  ENDIF 
                  ENDDO 
                  IF (ianz.eq.3) then 
                     DO i = 4, 6 
                     werte (i) = 0.0 
                     ENDDO 
                  ENDIF 
                  DO i = 1, 3 
                  u (i) = werte (i) 
                  v (i) = werte (i + 3) 
                  ENDDO 
                  lspace = .false. 
                  ww = do_blen (lspace, u, v) 
                  CALL ersetz2 (string, ikl, iklz, ww, 5, laenge) 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
            ENDIF 
         ELSE 
            ier_num = - 3 
            ier_typ = ER_FORT 
         ENDIF 
      ELSEIF (lcomm.eq.4) then 
!                                                                       
!     Calculate a bond angle                                            
!                                                                       
         IF (string (ikl - 4:ikl - 1) .eq.'bang') then 
            CALL get_params (line, ianz, cpara, lpara, 9, lp) 
            IF (ier_num.eq.0) then 
               IF (ianz.eq.3.or.ianz.eq.4.or.ianz.eq.6.or.ianz.eq.9)    &
               then                                                     
!     --------Input are atom numbers                                    
                  IF (str_comp (cpara (1) , 'atom', 1, lpara (1) , 4) ) &
                  then                                                  
                     DO i = 2, ianz 
                     CALL eval (cpara (i), lpara (i) ) 
                     IF (ier_num.ne.0) then 
                        GOTO 999 
                     ENDIF 
                     werte (i) = do_read_number (cpara (i), lpara (i) ) 
                     IF (ier_num.ne.0) then 
                        GOTO 999 
                     ENDIF 
                     ENDDO 
                     i = nint (werte (2) ) 
                     j = nint (werte (3) ) 
                     IF (ianz.eq.3) then 
      IF (i.lt.1.or.cr_natoms.lt.i.or.j.lt.1.or.cr_natoms.lt.j) then 
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
     &.or.cr_natoms.lt.k) then                                          
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
                  ELSEIF (str_comp(cpara(1), 'envi', 1, lpara(1), 4)) THEN
                     DO i = 2, ianz 
                     CALL eval (cpara (i), lpara (i) ) 
                     IF (ier_num.ne.0) then 
                        GOTO 999 
                     ENDIF 
                     werte (i) = do_read_number (cpara (i), lpara (i) ) 
                     IF (ier_num.ne.0) then 
                        GOTO 999 
                     ENDIF 
                     ENDDO 
                     i = nint (werte (2) ) 
                     j = nint (werte (3) ) 
                     IF (ianz.eq.3) then 
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
                     DO i = 1, ianz 
                     CALL eval (cpara (i), lpara (i) ) 
                     IF (ier_num.ne.0) then 
                        GOTO 999 
                     ENDIF 
                     werte (i) = do_read_number (cpara (i), lpara (i) ) 
                     IF (ier_num.ne.0) then 
                        GOTO 999 
                     ENDIF 
                     ENDDO 
                     IF (ianz.eq.6) then 
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
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'blen') then 
            CALL get_params (line, ianz, cpara, lpara, 6, lp) 
            IF (ier_num.eq.0) then 
               IF (ianz==2 .or. ianz.eq.3.or.ianz.eq.6) then 
!     --------Input are atom numbers                                    
                  IF (str_comp (cpara (1) , 'atom', 1, lpara (1) , 4) ) &
                  then                                                  
                     DO i = 2, ianz 
                     CALL eval (cpara (i), lpara (i) ) 
                     IF (ier_num.ne.0) then 
                        GOTO 999 
                     ENDIF 
                     werte (i) = do_read_number (cpara (i), lpara (i) ) 
                     IF (ier_num.ne.0) then 
                        GOTO 999 
                     ENDIF 
                     ENDDO 
                     i = nint (werte (2) ) 
                     j = nint (werte (3) ) 
      IF (i.lt.1.or.cr_natoms.lt.i.or.j.lt.1.or.cr_natoms.lt.j) then 
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
                  ELSEIF (str_comp (cpara(1), 'envi', 1, lpara(1), 4)) THEN
                     DO i = 2, ianz 
                     CALL eval (cpara (i), lpara (i) ) 
                     IF (ier_num.ne.0) then 
                        GOTO 999 
                     ENDIF 
                     werte (i) = do_read_number (cpara (i), lpara (i) ) 
                     IF (ier_num.ne.0) then 
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
                     IF (ier_num.ne.0) then 
                        GOTO 999 
                     ENDIF 
                     werte (i) = do_read_number (cpara (i), lpara (i) ) 
                     IF (ier_num.ne.0) then 
                        GOTO 999 
                     ENDIF 
                     ENDDO 
                     IF (ianz.eq.3) then 
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
         ELSEIF (string (ikl - 4:ikl - 1) .eq.'rang') then 
            CALL get_params (line, ianz, cpara, lpara, 9, lp) 
            IF (ier_num.eq.0) then 
               IF (ianz.eq.6.or.ianz.eq.9) then 
                  DO i = 1, ianz 
                  CALL eval (cpara (i), lpara (i) ) 
                  IF (ier_num.ne.0) then 
                     GOTO 999 
                  ENDIF 
                  werte (i) = do_read_number (cpara (i), lpara (i) ) 
                  IF (ier_num.ne.0) then 
                     GOTO 999 
                  ENDIF 
                  ENDDO 
                  IF (ianz.eq.6) then 
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
      ELSEIF (lcomm.eq.0) then 
         CALL ersetz2 (string, ikl, iklz, ww, 0, laenge) 
      ELSE 
         ier_num = - 3 
         ier_typ = ER_FORT 
      ENDIF 
!                                                                       
  999 CONTINUE 
!                                                                       
      IF (ier_num.ne.0) then 
         WRITE (ier_msg (1), 9000) string (1:min (40, laenge) ) 
         WRITE (ier_msg (1), 9000) line (1:min (40, laenge) ) 
      ENDIF 
!                                                                       
 9000 FORMAT    (a) 
      END SUBROUTINE discus_calc_intr_spec                 
!*****7**************************************************************** 
      SUBROUTINE discus_validate_var_spec (zeile, lp) 
!-                                                                      
!       checks whether the variable name is legal, DISCUS specific part 
!                                                                       
!       Version : 1.0                                                   
!       Date    : 22 Oct 03                                             
!                                                                       
!       Author  : R.B. Neder  (reinhard.neder@fau.de)    
!+                                                                      
      USE discus_config_mod 
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      CHARACTER (LEN=*), INTENT(IN) :: zeile 
      INTEGER,           INTENT(IN) :: lp 
!                                                                       
      INTEGER, PARAMETER :: reserved_n = 33 
                                                                        
      CHARACTER(LEN=12), DIMENSION(1:reserved_n) :: reserved
      INTEGER  :: i 
!                                                                       
      DATA reserved / 'bang', 'blen', 'dstar', 'md_test', 'mol_test',   &
      'rang', 'scalpro', 'x', 'y', 'z', 'm', 'b', 'n', 'cdim', 'lat',   &
      'vol', 'rvol', 'rlat', 'env', 'menv', 'md_num', 'md_cre',         &
      'mc_num', 'md_rad', 'mr_run', 'mc_type', 'mc_orig', 'mc_rad',     &
      'md_next', 'md_dist', 'mol_cont', 'mol_dens', 'mol_len' /         
!                                                                       
      ier_num = 0 
      ier_typ = ER_NONE 
!                                                                       
      DO i = 1, reserved_n 
      IF (index (reserved (i), zeile (1:lp) ) .ne.0) then 
         ier_num = - 25 
         ier_typ = ER_FORT 
      ENDIF 
      ENDDO 
!                                                                       
      END SUBROUTINE discus_validate_var_spec              
