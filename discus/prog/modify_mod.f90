MODULE modify_mod
!+                                                                      
!       Basic atom/molecule manipulation functions and conversion       
!     and other support routines routines.                              
!                                                                       
CONTAINS
!*****7*****************************************************************
      SUBROUTINE do_replace (zeile, lp) 
!-                                                                      
!     Replaces atom(s) or molecule(s) within the crystal ..             
!+                                                                      
      USE discus_config_mod 
      USE discus_allocate_appl_mod
      USE crystal_mod 
      USE celltoindex_mod
      USE get_iscat_mod
      USE modify_func_mod
      USE molecule_mod 
      USE prop_para_mod 
      USE update_cr_dim_mod
      USE ber_params_mod
      USE errlist_mod 
USE lib_errlist_func
USE lib_random_func
      USE get_params_mod
      USE random_mod
USE precision_mod
USE str_comp_mod
!
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 200) 
!                                                                       
CHARACTER(LEN=*), INTENT(INOUT) :: zeile 
INTEGER         , INTENT(INOUT) :: lp 
!
      CHARACTER(LEN=MAX(PREC_STRING,LEN(zeile))) :: cpara (maxw), cc 
      INTEGER lpara (maxw) , ccl
      INTEGER ianz, iianz, jjanz, i 
      INTEGER is1, is2, isite, itype 
      INTEGER ja, jsite, jcell (3) 
      INTEGER            :: n_scat ! dummy for allocation
      LOGICAL lexist, lrepl 
      REAL(KIND=PREC_DP) :: uerte (maxw) 
      REAL(KIND=PREC_DP) :: verte (maxw) 
      REAL(KIND=PREC_DP) :: werte (maxw) 
      REAL prob 
!                                                                       
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
!                                                                       
      IF (ianz.eq.2) then 
         cpara (3) = 'x' 
         lpara (3) = 1 
      ENDIF 
!                                                                       
!------ Here are molecule manipulations                                 
!                                                                       
      IF (str_comp (cpara (3) , 'mol', 2, lpara (3) , 3) ) then 
!                                                                       
!------ - Single molecule will be changed                               
!                                                                       
         IF (ianz.eq.3) then 
            ianz = ianz - 1 
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
            is1 = nint (werte (1) ) 
            is2 = nint (werte (2) ) 
            CALL do_swap_mole (is1, is2, .false.) 
            CALL update_cr_dim 
!                                                                       
!------ - More molecules will be changed                                
!                                                                       
         ELSEIF (ianz.eq.4) then 
            CALL ber_params (2, cpara, lpara, werte, maxw) 
            itype = nint (werte (1) ) 
            is2 = nint (werte (2) ) 
            CALL del_params (3, ianz, cpara, lpara, maxw) 
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
            prob = werte (1) 
!                                                                       
            IF (itype.gt.0.and.itype.le.mole_num_type) then 
               DO ja = 1, mole_num_mole 
               IF (mole_type (ja) .eq.itype.and.ran1 (idum) .le.prob)   &
               then                                                     
                  CALL do_swap_mole (ja, is2, .false.) 
               ENDIF 
               ENDDO 
               CALL update_cr_dim 
            ELSE 
               ier_num = - 64 
               ier_typ = ER_APPL 
            ENDIF 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
!                                                                       
!------ Here are atom manipulations                                     
!                                                                       
      ELSE 
!                                                                       
!------ - Replace all given atoms with given probability                
!                                                                       
         IF (ianz.eq.4) then 
            iianz = 1 
            jjanz = 1 
            CALL get_iscat (iianz, cpara (1), lpara (1), uerte, maxw,   &
            .false.)                                                    
!                                                                       
!     ----The first atom type does not exist, user error                
!                                                                       
            IF (ier_num.ne.0) return 
!                                                                       
!     ----Get second atom                                               
!                                                                       
            CALL get_iscat (jjanz, cpara (2), lpara (2), verte, maxw,   &
            .false.)                                                    
            IF (ier_num.eq. - 27) then 
!                                                                       
!------ ----The second atom type does not exist, create new             
!           scattering curve(s)                                         
!                                                                       
               DO i = 1, iianz 
               IF (cr_nscat + 1 >  MAXSCAT) then 
                  n_scat = cr_nscat + 1
                  call alloc_crystal (n_scat, NMAX)
               ENDIF
               IF (cr_nscat + 1.le.maxscat) then 
                  cr_nscat = cr_nscat + 1 
                  cr_at_lis (cr_nscat) = cpara (2) (1:lpara(2))
!DBG                cr_dw(cr_nscat)     = cr_dw(nint(uerte(i)))         
                  cr_dw (cr_nscat) = 0.05 
                  CALL no_error 
                  verte (i) = REAL(cr_nscat) 
               ELSE 
                  ier_num = - 26 
                  ier_typ = ER_APPL 
                  RETURN 
               ENDIF 
               ENDDO 
            ELSEIF (ier_num.ne.0) then 
!                                                                       
!     ------other errors, return                                        
!                                                                       
               RETURN 
            ENDIF 
            is1 = nint (uerte (1) ) 
            is2 = nint (verte (1) ) 
            IF (str_comp (cpara (3) , 'all', 1, lpara (3) , 3) ) then 
               CALL del_params (3, ianz, cpara, lpara, maxw) 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               isite = - 1 
               prob = werte (1) 
            ELSE 
               CALL del_params (2, ianz, cpara, lpara, maxw) 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               isite = nint (werte (1) ) 
               prob = werte (2) 
            ENDIF 
!                                                                       
      IF (prob.le.0.0.or.prob.gt.1.0.or.isite.gt.cr_ncatoms.or.is1.eq. -&
     & 1.or.is2.eq. - 1) then                                           
               ier_num = - 6 
               ier_typ = ER_COMM 
               RETURN 
            ENDIF 
!                                                                       
            DO ja = 1, cr_natoms 
            IF (isite.eq. - 1) then 
               jsite = - 1 
            ELSE 
               CALL indextocell (ja, jcell, jsite) 
            ENDIF 
            DO i = 1, iianz 
            IF (cr_iscat (ja) .eq.nint (uerte (i) ) .and. (             &
            jsite.eq.isite.or.isite.eq. - 1) .and.ran1 (idum) .le.prob) &
            then                                                        
      IF (check_select_status (ja, .true., cr_prop (ja),  cr_sel_prop) ) THEN
                  cr_iscat (ja) = nint (verte (i) ) 
                  IF (nint (verte (i) ) .gt.0) then 
                     cr_prop (ja) = IBSET (cr_prop (ja), PROP_NORMAL) 
                  ELSE 
                     cr_prop (ja) = IBCLR (cr_prop (ja), PROP_NORMAL) 
                  ENDIF 
                  GOTO 999 
               ENDIF 
            ENDIF 
            ENDDO 
  999       CONTINUE 
            lrepl = .false. 
!           do i=1,iianz                                                
!             lrepl = lrepl .or. cr_iscat(ja).eq.nint(uerte(i))         
!     &                     .and. (jsite.eq.isite .or. isite.eq.-1)     
!           ENDDO                                                       
!           if (cr_iscat(ja).eq.is1               .and.                 
!     &                     (jsite.eq.isite .or. isite.eq.-1) .and.     
!           if (lrepl .and.                                             
!     &         ran1(idum).le.prob                     ) then           
!             cr_iscat(ja) = is2                                        
!           endif                                                       
            ENDDO 
!
!                                                                       
!------ - Replace a single atom                                         
!                                                                       
         ELSEIF (ianz.eq.2.or.ianz.eq.3) then 
            CALL ber_params (1, cpara, lpara, werte, maxw) 
            ja = nint (werte (1) ) 
            iianz = 1 
            CALL get_iscat (iianz, cpara (2), lpara (2), uerte, maxw,   &
            .false.)                                                    
            DO i = 1, iianz
               IF(uerte(i)==-1) THEN
                  ier_num = -6
                  ier_typ = ER_COMM
                  RETURN
               ENDIF
            ENDDO
!                                                                       
            IF (ianz.eq.3.and.ier_num.eq. - 27) then 
               IF (cr_nscat + 1 >  MAXSCAT) then 
                  n_scat = cr_nscat + 1
                  call alloc_crystal (n_scat, NMAX)
               ENDIF
               IF (cr_nscat + 1.le.maxscat) then 
                  cr_nscat = cr_nscat + 1 
                  cr_at_lis (cr_nscat) = cpara (2)  (1:lpara(2))
                  CALL del_params (2, ianz, cpara, lpara, maxw) 
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  cr_dw (cr_nscat) = werte (1) 
                  is1 = cr_nscat 
                  CALL no_error 
               ELSE 
                  ier_num = - 26 
                  ier_typ = ER_APPL 
                  RETURN 
               ENDIF 
            ELSEIF (ianz.eq.3.and.ier_num.eq.0) then 
!                                                                       
!     ------Atom name exists, but since DW was given check this         
!                                                                       
               cc = cpara (2) 
               ccl = lpara (2)
               CALL del_params (2, ianz, cpara, lpara, maxw) 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                                                                        
               lexist = .false. 
               DO i = 1, iianz 
               IF (ABS(cr_dw (NINT(uerte (i) ) ) - werte(1) )< 1.0D-4) THEN 
!                                                                       
!     --------Atom exists with identical DW                             
!                                                                       
                  is1 = nint (uerte (i) ) 
                  lexist = .true. 
               ENDIF 
               ENDDO 
               IF (.not.lexist) then 
                  IF (cr_nscat + 1 >  MAXSCAT) then 
                     n_scat = cr_nscat + 1
                     call alloc_crystal (n_scat, NMAX)
                  ENDIF
                  IF (cr_nscat + 1.le.maxscat) then 
                     cr_nscat = cr_nscat + 1 
                     cr_at_lis (cr_nscat) = cc(1:ccl)
                     cr_dw (cr_nscat) = werte (1) 
                     is1 = cr_nscat 
                     CALL no_error 
                  ELSE 
                     ier_num = - 26 
                     ier_typ = ER_APPL 
                     RETURN 
                  ENDIF 
               ENDIF 
            ELSE 
               is1 = nint (uerte (1) ) 
            ENDIF 
!                                                                       
            IF (is1.eq. - 1.or.ja.le.0.or.ja.gt.cr_natoms) then 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
            IF (ier_num.ne.0) return 
!                                                                       
            cr_iscat (ja) = is1 
            IF (cr_iscat (ja) .gt.0) then 
               cr_prop (ja) = IBSET (cr_prop (ja), PROP_NORMAL) 
            ELSE 
               cr_prop (ja) = IBCLR (cr_prop (ja), PROP_NORMAL) 
            ENDIF 
!                                                                       
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
      ENDIF 
!                                                                       
      END SUBROUTINE do_replace                     
!*****7*****************************************************************
      SUBROUTINE do_swap_mole (idest, isource, lswap) 
!-                                                                      
!     Replaces molecule 'idest' with 'isource'. If 'lswap' is true      
!     both molecule are swapped rather than one replaced by the other.  
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE molecule_mod 
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER idest, isource, ityp, i, k, ii, jj, is, js, i0, j0 
      INTEGER  :: im, jm    ! molecule number 
      INTEGER  :: ip, jp    ! property values
      INTEGER, DIMENSION(0:3) :: iis, jjs ! Surface values
      REAL   , DIMENSION(0:3) :: ris, rjs ! Magnetic moments
      REAL i0pos (3), j0pos (3), ipos (3) 
      LOGICAL lswap 
!                                                                       
      IF (mole_len (idest) .ne.mole_len (isource) ) then 
         ier_num = - 67 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
!                                                                       
      IF (idest.lt.0.or.idest.gt.mole_num_mole.or.isource.lt.0.or.isourc&
     &e.gt.mole_num_mole) then                                          
         ier_num = - 63 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
!                                                                       
!------ Get origin of molecules (atom 1)                                
!                                                                       
      i0 = mole_cont (mole_off (idest) + 1) 
      j0 = mole_cont (mole_off (isource) + 1) 
      DO i = 1, 3 
      i0pos (i) = cr_pos (i, i0) 
      j0pos (i) = cr_pos (i, j0) 
      ENDDO 
!                                                                       
!------ replace/swap atoms                                              
!                                                                       
      DO k = 1, mole_len (idest) 
      ii = mole_cont (mole_off (idest) + k) 
      jj = mole_cont (mole_off (isource) + k) 
      is = cr_iscat (ii) 
      js = cr_iscat (jj) 
      im = cr_mole (ii) 
      jm = cr_mole (jj) 
      iis(:) = cr_surf (:,ii)
      jjs(:) = cr_surf (:,jj)
      ris(:) = cr_magn (:,ii)
      rjs(:) = cr_magn (:,jj)
      ip = cr_prop (ii) 
      jp = cr_prop (jj) 
!                                                                       
      cr_iscat (ii) = js 
!     cr_mole (ii) = jm 
      cr_surf(:,ii)= jjs(:)
      cr_magn(:,ii)= rjs(:)
      cr_prop (ii) = jp 
      DO i = 1, 3 
      ipos (i) = cr_pos (i, ii) 
      cr_pos (i, ii) = cr_pos (i, jj) - j0pos (i) + i0pos (i) 
      ENDDO 
!                                                                       
      IF (lswap) then 
         cr_iscat (jj) = is 
!        cr_mole (jj) = im 
         cr_surf(:,jj)= iis(:)
         cr_magn(:,jj)= ris(:)
         cr_prop (jj) = ip 
         DO i = 1, 3 
         cr_pos (i, jj) = ipos (i) - i0pos (i) + j0pos (i) 
         ENDDO 
      ENDIF 
      ENDDO 
!                                                                       
!------ replace/swap molecule types                                     
!                                                                       
      ityp = mole_type (idest) 
      mole_type (idest) = mole_type (isource) 
      IF (lswap) mole_type (isource) = ityp 
!
!                                                                       
      END SUBROUTINE do_swap_mole                   
!*****7*****************************************************************
      SUBROUTINE do_ins (line, laenge) 
!-                                                                      
!     Inserts a new atom into the structure. The name and position      
!     must be given, optionally the temperature factor will be read.    
!     If an atom of this type exists the temperature factor will be     
!     that of the first type of temperature factor given for this atom  
!     type, ELSE the temperature factor will be 0.                      
!+                                                                      
      USE discus_config_mod 
      USE charact_mod 
      USE discus_allocate_appl_mod
      USE crystal_mod 
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
      USE build_name_mod
USE precision_mod
!
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 15) 
!                                                                       
      CHARACTER (LEN=*), INTENT(IN) :: line
      INTEGER          , INTENT(INOUT) :: laenge
!
      CHARACTER(LEN=4)              :: name 
      CHARACTER(LEN=MAX(PREC_STRING,LEN(line))), DIMENSION(maxw) :: cpara
      INTEGER            , DIMENSION(maxw) :: lpara
      INTEGER i, j, ianz
      INTEGER                              :: new_nmax   = 0
      INTEGER                              :: new_nscat  = 0
      LOGICAL                              :: need_alloc = .false.
      REAL(KIND=PREC_DP) :: werte (maxw) 
!
!     While developing, increment crystal if neede, but keep the check
!
      need_alloc = .false.
      new_nmax   = NMAX
      new_nscat  = MAXSCAT
      IF ( NMAX <= cr_natoms ) THEN
         new_nmax  = max(NMAX+1   , INT(NMAX    * 1.25))
         need_alloc = .true.
      ENDIF
      IF ( MAXSCAT <= cr_nscat ) THEN 
         new_nscat = max(MAXSCAT+1, INT(MAXSCAT * 1.25))
         need_alloc = .true.
      ENDIF
      IF( need_alloc ) THEN
         call alloc_crystal(new_nscat, new_nmax)
         IF ( ier_num /= 0) RETURN
      ENDIF
      IF (cr_natoms.lt.nmax) then 
         CALL get_params (line, ianz, cpara, lpara, maxw, laenge) 
         CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
         IF (ier_num.eq.0) then 
            i = ichar (cpara (1) (1:1) ) 
            IF ( (a.le.i.and.i.le.z) .or. (aa.le.i.and.i.le.zz) ) then 
               name = cpara (1) (1:lpara(1))
               j = 2 
               cpara(1) = '0'
               lpara(1) = 1
            ELSE 
               j = 1 
               name = 'yyyy' 
            ENDIF 
            CALL ber_params(ianz, cpara, lpara, werte, maxw)
!            DO i = j, ianz 
!            lp = lpara (i) 
!            zeile = ' ' 
!            zeile (1:1) = '(' 
!            zeile (2:lp + 1) = cpara (i) (1:lp) 
!            zeile (lp + 2:lp + 2) = ')' 
!            lp = lp + 2 
!            werte (i) = 0 
!            werte (i) = berechne (zeile, lp) 
            IF (ier_num.ne.0) return 
!            ENDDO 
            IF (ianz.ge.4) then 
               CALL do_ins_atom (name, maxw, werte) 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
      ELSE 
         ier_num = -10 
         ier_typ = ER_APPL 
      ENDIF 
      END SUBROUTINE do_ins                         
!*****7*****************************************************************
      SUBROUTINE do_app (line, laenge, lkick) 
!-                                                                      
!     Inserts a new atom into the structure. The name and position      
!     must be given, optionally the temperature factor will be read.    
!     If an atom of this type exists the temperature factor will be     
!     that of the first type of temperature factor given for this atom  
!     type, ELSE the temperature factor will be 0.                      
!                                                                       
!     Appends a new atom                                                
!+                                                                      
      USE discus_config_mod 
      USE charact_mod 
      USE discus_allocate_appl_mod
      USE crystal_mod 
      USE metric_mod
      USE prop_para_mod 
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
USE precision_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 10) 
!                                                                       
      CHARACTER ( * ) line 
      CHARACTER(4) name 
      CHARACTER(LEN=MAX(PREC_STRING,LEN(line))) :: cpara (maxw)
      INTEGER lpara (maxw) 
      INTEGER i, j, k, ianz, laenge
      INTEGER                              :: new_nmax
      INTEGER                              :: new_nscat
      LOGICAL                              :: need_alloc = .false.
      LOGICAL lkick, lspace 
      REAL(KIND=PREC_DP) :: werte (maxw) 
      REAL w (3), v (3) 
!     REAL do_blen 
!                                                                       
      DATA lspace / .true. / 
!                                                                       
!
!     While developing, increment crystal if neede, but keep the check
!
      need_alloc = .false.
      new_nmax   = NMAX
      new_nscat  = MAXSCAT
      IF ( NMAX <= cr_natoms ) then 
         new_nmax  = max(NMAX+1   , INT(NMAX    * 1.25))
         need_alloc = .true.
      ENDIF
      IF ( MAXSCAT <= cr_nscat ) then 
         new_nscat = max(MAXSCAT+1, INT(MAXSCAT * 1.25))
         need_alloc = .true.
      ENDIF
      IF ( need_alloc ) THEN
         call alloc_crystal(new_nscat, new_nmax)
         IF ( ier_num /= 0) RETURN
      ENDIF
!                                                                       
!     If there is still space in the crystal, try to insert the atom    
!                                                                       
      IF (cr_natoms.lt.nmax) then 
         CALL get_params (line, ianz, cpara, lpara, maxw, laenge) 
         IF (ier_num.eq.0) then 
            IF (ianz.ge.8) then 
!                                                                       
!------ ------A minimum of eight parameters is required by 'append'     
!             and 'kick'                                                
!                                                                       
               i = ichar (cpara (1) (1:1) ) 
               IF ( (a.le.i.and.i.le.z) .or. (aa.le.i.and.i.le.zz) )    &
               then                                                     
                  name = cpara (1)(1:lpara(1))
                  j = 2 
                  cpara(1) = '0'
                  lpara(1) = 1
               ELSE 
                  j = 1 
                  name = 'yyyy' 
               ENDIF 
               CALL ber_params(ianz, cpara, lpara, werte, maxw)
!                                                                       
!     ------Evaluate parameters, old form will be replaced by ber_params
!                                                                       
!               DO i = j, ianz 
!               lp = lpara (i) 
!               zeile = ' ' 
!               zeile (1:1) = '(' 
!               zeile (2:lp + 1) = cpara (i) (1:lp) 
!               zeile (lp + 2:lp + 2) = ')' 
!               lp = lp + 2 
!               werte (i) = 0 
!               werte (i) = berechne (zeile, lp) 
               IF (ier_num.ne.0) return 
!               ENDDO 
               j = nint (werte (6) ) 
               k = nint (werte (7) ) 
!                                                                       
!     ------j,k is the range of atoms to compare with the new atom      
!           is this range within the current crystal ?                  
!                                                                       
      IF ( (0.lt.j.and.j.le.cr_natoms.and.0.lt.k.and.k.le.cr_natoms.and.&
     &j.le.k) .or. (j.eq.cr_natoms + 1.and.k.eq.cr_natoms) ) then       
!                                                                       
!     --------Apply default values                                      
!                                                                       
                  IF (ianz.lt.9) werte (9) = werte (8) 
                  IF (ianz.lt.10) werte (10) = werte (9) 
!                                                                       
                  IF(werte(8) .lt.0.0D0) then 
!                                                                       
!     ----------Interpret parameter 8 as bondlength                     
!                                                                       
                     IF (lkick) then 
!                                                                       
!     ------------'kick' mode, remove all atoms from crystal that are   
!                              too close                                
!                                                                       
                        DO i = j, k 
                        w (1) = cr_pos (1, i) 
                        w (2) = cr_pos (2, i) 
                        w (3) = cr_pos (3, i) 
                        v (1) = werte (2) 
                        v (2) = werte (3) 
                        v (3) = werte (4) 
                        IF (do_blen (lspace, w, v) .lt. - werte (8) )   &
                        then                                            
                           cr_iscat (i) = 0 
      cr_prop (i)  = ibclr (cr_prop (i),  PROP_NORMAL) 
                        ENDIF 
                        ENDDO 
                     ELSE 
!                                                                       
!------ ------------'append' mode, insert new atom only if no other     
!                            atom is too close. Vacancies are ignored.  
!                                                                       
                        DO i = j, k 
                        IF (cr_at_lis (cr_iscat (i) ) .ne.'VOID') then 
                           w (1) = cr_pos (1, i) 
                           w (2) = cr_pos (2, i) 
                           w (3) = cr_pos (3, i) 
                           v (1) = werte (2) 
                           v (2) = werte (3) 
                           v (3) = werte (4) 
                           IF (do_blen (lspace, w, v) .lt. - werte (8) )&
                           then                                         
                              GOTO 10 
                           ENDIF 
                        ENDIF 
                        ENDDO 
                     ENDIF 
                     CALL do_ins_atom (name, maxw, werte) 
                  ELSE 
!                                                                       
!     ----------Interpret parameters 8,9,10 as direct lattice units     
!                                                                       
                     IF (lkick) then 
!                                                                       
!     ------------'kick' mode, remove all atoms from crystal that are   
!                              too close                                
!                                                                       
                        DO i = j, k 
                        IF (abs (cr_pos (1, i) - werte (2) ) .lt.werte (&
                        8) .and.abs (cr_pos (2, i) - werte (3) )        &
                        .lt.werte (9) .and.abs (cr_pos (3, i) - werte ( &
                        4) ) .lt.werte (10) ) then                      
!                   cr_pos (1,i)=0.0                                    
!                   cr_pos (2,i)=0.0                                    
!                   cr_pos (3,i)=0.0                                    
                           cr_iscat (i) = 0 
      cr_prop (i)  = ibclr (cr_prop (i),  PROP_NORMAL) 
                        ENDIF 
                        ENDDO 
                     ELSE 
!                                                                       
!     ------------'append' mode, insert new atom only if no other       
!                         atom is too close. Vacancies are ignored.     
!                                                                       
                        DO i = j, k 
                        IF (cr_at_lis (cr_iscat (i) ) .ne.'VOID') then 
                           IF (abs (cr_pos (1, i) - werte (2) )         &
                           .lt.werte (8) .and.abs (cr_pos (2, i)        &
                           - werte (3) ) .lt.werte (9) .and.abs (cr_pos &
                           (3, i) - werte (4) ) .lt.werte (10) ) then   
                              GOTO 10 
                           ENDIF 
                        ENDIF 
                        ENDDO 
                     ENDIF 
                     CALL do_ins_atom (name, maxw, werte) 
                  ENDIF 
               ELSEIF (cr_natoms.eq.0) then 
!                                                                       
!     --------The crystal is empty. Insert atom in any case             
!                                                                       
                  CALL do_ins_atom (name, maxw, werte) 
               ELSE 
!                                                                       
!     --------Wrong range of atoms to compare to the new one            
!                                                                       
                  ier_num = - 2 
                  ier_typ = ER_APPL 
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
   10 CONTINUE 
      END SUBROUTINE do_app                         
!****7******************************************************************
      SUBROUTINE do_ins_atom (name, maxw, werte) 
!-                                                                      
!     Inserts the atom given by name and position in werte into the     
!     structure.                                                        
!+                                                                      
      USE discus_config_mod 
      USE discus_allocate_appl_mod
      USE chem_mod
      USE crystal_mod 
      USE prop_para_mod 
      USE errlist_mod 
USE precision_mod
      USE string_convert_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER, INTENT(IN)     :: maxw
!                                                                       
      CHARACTER (LEN=* )    , INTENT(INOUT) :: name 
      REAL(KIND=PREC_DP) , DIMENSION(maxw), INTENT(IN) :: werte (maxw) 
!
      INTEGER                :: i, l
      INTEGER                :: new_nmax   = 1
      INTEGER                :: new_nscat  = 1
      LOGICAL                :: need_alloc = .false.
      LOGICAL                :: lda 
!                                                                       
!
!     While developing, increment crystal if needed, but keep the check
!
      need_alloc = .false.
      new_nmax   = NMAX
      new_nscat  = MAXSCAT
      IF ( NMAX <= cr_natoms ) then 
         new_nmax  = max(NMAX+1   , INT(NMAX    * 1.25))
         need_alloc = .true.
      ENDIF
      IF ( MAXSCAT <= cr_nscat ) then 
         new_nscat = max(MAXSCAT+1, INT(MAXSCAT * 1.25))
         need_alloc = .true.
      ENDIF
      IF ( need_alloc ) THEN
         call alloc_crystal(new_nscat, new_nmax)
         IF ( ier_num /= 0) RETURN
      ENDIF
      IF (cr_natoms.lt.NMAX) then 
         CALL do_cap (name) 
         i = 0 
         IF (name.eq.'YYYY') then 
            i = nint (werte (1) ) 
            IF (0.le.i.and.i.le.cr_nscat) then 
               ier_num = 0 
               ier_typ = ER_NONE 
               cr_natoms = cr_natoms + 1 
               cr_iscat (cr_natoms) = i 
            ELSE 
               ier_num = - 27 
               ier_typ = ER_APPL 
            ENDIF 
         ELSE 
            lda = name.eq.cr_at_lis (i) .and.ABS(werte(5)- DBLE(cr_dw(i)))<1.D-5
            DO while (.not.lda.and.i.lt.cr_nscat) 
            i = i + 1 
            lda = name.eq.cr_at_lis (i) .and.ABS(werte(5) - DBLE(cr_dw(i))) < 1.D-5
            ENDDO 
            IF (lda) then 
               cr_natoms = cr_natoms + 1 
               cr_iscat (cr_natoms) = i 
            ELSE 
               IF (cr_nscat + 1.le.maxscat) then 
                  cr_natoms = cr_natoms + 1 
                  cr_nscat = cr_nscat + 1 
                  cr_iscat (cr_natoms) = cr_nscat 
                  cr_at_lis (cr_nscat) = name 
                  cr_dw (cr_nscat) = werte (5) 
                  cr_occ(cr_nscat) = 1.0
               ELSE 
                  ier_num = - 26 
                  ier_typ = ER_APPL 
               ENDIF 
            ENDIF 
         ENDIF 
         IF (ier_num.eq.0) then 
            cr_pos (1, cr_natoms) = werte (2) 
            cr_pos (2, cr_natoms) = werte (3) 
            cr_pos (3, cr_natoms) = werte (4) 
            cr_mole (cr_natoms) = 0 
            cr_surf (:,cr_natoms) = 0
            cr_magn (:,cr_natoms) = 0.0
            cr_prop (cr_natoms) = 0 
            cr_prop (cr_natoms)  = ibset (cr_prop (cr_natoms),  PROP_NORMAL) 
            DO l = 1, 3 
               cr_dim(l,1) = amin1(cr_dim(l,1), cr_pos(l,cr_natoms))                                                           
               cr_dim(l,2) = amax1(cr_dim(l,2), cr_pos(l,cr_natoms))                                                           
            ENDDO 
            chem_period(:) = .FALSE.    ! Turn off periodic boundary
            chem_quick     = .FALSE.    ! Turn off quick search mode
         ENDIF 
      ELSE 
         ier_num = - 10 
         ier_typ = ER_APPL 
      ENDIF 
      END SUBROUTINE do_ins_atom                    
!*****7*****************************************************************
SUBROUTINE do_remove (line, ll) 
!-                                                                      
!     Removes a single atom from the structure. line is evaluated to    
!     give the index of the atom.                                       
!+                                                                      
USE discus_config_mod 
USE crystal_mod 
USE molecule_mod 
USE prop_para_mod 
USE ber_params_mod
USE errlist_mod 
USE get_params_mod
USE precision_mod
USE str_comp_mod
IMPLICIT none 
!                                                                       
       
!                                                                       
INTEGER, PARAMETER :: maxw = 6
!                                                                       
CHARACTER( LEN=* ), INTENT(INOUT) :: line 
INTEGER           , INTENT(INOUT) :: ll
!
CHARACTER(LEN=MAX(PREC_STRING,LEN(line))), DIMENSION(MAXW) :: cpara !(maxw) 
INTEGER            , DIMENSION(MAXW) :: lpara !(maxw) 
INTEGER :: i, j 
INTEGER :: istart, iend, ianz 
INTEGER :: tstart, tend 
LOGICAL :: l_mole 
LOGICAL :: l_type 
REAL(KIND=PREC_DP) , DIMENSION(MAXW) :: werte ! (maxw) 
!                                                                       
!                                                                       
l_mole = .false. 
!                                                                       
!     Get parameters                                                    
!                                                                       
CALL get_params (line, ianz, cpara, lpara, maxw, ll) 
IF(ier_num.eq.0) THEN 
!                                                                       
!     --Remove atoms or molecules ?                                     
!                                                                       
   l_mole = .false. 
   tstart = 1 
   tend = mole_num_type 
   IF(ianz >  2 .AND. str_comp(cpara(3), 'molecule', 1, lpara(3), 8)) THEN
      cpara (3) = '0' 
      lpara (3) = 1 
      l_mole = .true. 
   ELSEIF (ianz.gt.1 .AND. str_comp (cpara (2) , 'molecule', 1,     &
      lpara (2) , 8) ) THEN                                          
      DO i = ianz, 3, - 1 
         cpara (i + 1) = cpara (i) 
         lpara (i + 1) = lpara (i) 
      ENDDO 
      cpara (3) = '0' 
      lpara (3) = 1 
      l_mole = .true. 
      cpara (2) = cpara (1) 
      lpara (2) = lpara (1) 
      ianz = ianz + 1 
   ENDIF 
   IF(ianz.gt.4 .AND. str_comp (cpara (4) , 'type', 1, lpara (4) , 8) ) THEN
      cpara (4) = '0' 
      lpara (4) = 1 
      l_type = .true. 
   ENDIF 
!                                                                       
!     --Calculate value of parameters                                   
!                                                                       
   IF (str_comp (cpara (1) , 'all', 1, lpara (1) , 3) ) THEN 
!                                                                       
!     ----Remove all molecules or all atoms                             
!                                                                       
      WRITE (cpara (1), 3000) 1 
      IF (l_mole) THEN 
         WRITE (cpara (2), 3000) mole_num_mole 
      ELSE 
         WRITE (cpara (2), 3000) cr_natoms 
         IF (ianz.eq.1) ianz = 2 
      ENDIF 
      lpara (1) = 11 
      lpara (2) = 11 
   ELSEIF (str_comp (cpara (1) , 'last', 1, lpara (1) , 4) ) THEN 
!                                                                       
!     ----Remove last molecule or last atom                             
!                                                                       
      IF (l_mole) THEN 
         WRITE (cpara (1), 3000) mole_num_mole 
         WRITE (cpara (2), 3000) mole_num_mole 
      ELSE 
         WRITE (cpara (1), 3000) cr_natoms 
         WRITE (cpara (2), 3000) cr_natoms 
         IF (ianz.eq.1) ianz = 2 
      ENDIF 
      lpara (1) = 11 
     lpara (2) = 11 
   ENDIF 
!                                                                       
!     --Select type of molecule                                         
!                                                                       
   IF (str_comp (cpara (5) , 'all', 1, lpara (5) , 3) ) THEN 
!                                                                       
!     ----Remove all types of molecules or all atoms                    
!                                                                       
      WRITE (cpara (5), 3000) 1 
      WRITE (cpara (6), 3000) mole_num_type 
      lpara (5) = 11 
      lpara (6) = 11 
   ELSEIF(str_comp (cpara (5) , 'last', 1, lpara (5) , 4) ) THEN 
!                                                                       
!     ----Remove only the last type of molecule                         
!                                                                       
      WRITE (cpara (5), 3000) mole_num_type 
      WRITE (cpara (6), 3000) mole_num_type 
      lpara (5) = 11 
      lpara (6) = 11 
   ENDIF 
!                                                                       
   CALL ber_params (ianz, cpara, lpara, werte, maxw) 
   IF (ier_num.eq.0) THEN 
      istart = NINT (werte (1) ) 
      IF (ianz.eq.1) THEN 
         iend = istart 
      ELSE 
         iend = nint (werte (2) ) 
      ENDIF 
!                                                                       
!     ----Remove molecules                                              
!                                                                       
      IF (l_mole) THEN            !: Molecules yes/no
         IF(0.lt.istart .AND. istart.le.iend .AND. iend.le.mole_num_mole) THEN
!                                                                       
!     --------Set limits for types that can be removed                  
!                                                                       
               IF (ianz.eq.5) THEN 
                  tstart = INT( werte (5) )
                  tend = INT (werte (5) )
               ELSEIF (ianz.eq.6) THEN 
                  tstart = INT( werte (5) )
                  tend = INT( werte (6) )
               ENDIF 
            IF(0.lt.tstart .AND. tstart.le.tend .AND. tend.le.mole_num_type) THEN 
!                                                                    
!  ----------Loop over all molecules in range                        
!                                                                       
               DO i = istart, iend 
                  IF (tstart.le.mole_type (i) .AND. mole_type (i) .le.tend) THEN
!                                                                       
!     --------------Set molecule type to zero and remove atoms          
!                                                                       
                     mole_type (i) = 0 
                     mole_char (i) = 0 
                     DO j = 1, mole_len (i) 
                        cr_iscat(mole_cont (mole_off (i) + j) ) = 0                                             
!                       cr_mole (mole_cont (mole_off (i) + j) ) = 0
!                       cr_surf (:,mole_cont (mole_off (i) + j) ) = 0
                        cr_prop (mole_cont (mole_off (i) + j) ) =  &
                            IBCLR(cr_prop (mole_cont (mole_off (i) + j) ), PROP_NORMAL)                                    
                     ENDDO 
                  ENDIF 
               ENDDO 
            ELSE 
               ier_num = - 64 
               ier_typ = ER_APPL 
            ENDIF 
         ELSE 
            ier_num = - 63 
            ier_typ = ER_APPL 
         ENDIF 
      ELSE            !: Molecules yes/no
!                                                                       
!     ------If index of atom is within limits, remove atom by           
!           Setting its scattering curve to zero                        
!                                                                       
         IF(0.lt.istart .AND. istart.le.iend .AND. iend.le.cr_natoms) THEN
            DO i = istart, iend 
               cr_iscat (i) = 0 
               cr_prop (i) = IBCLR (cr_prop (i), PROP_NORMAL) 
            ENDDO 
         ELSE 
            ier_num = - 19 
            ier_typ = ER_APPL 
         ENDIF 
      ENDIF            !: Molecules yes/no
   ELSE 
      ier_num = - 6 
      ier_typ = ER_COMM 
   ENDIF 
ENDIF 
!                                                                       
3000 FORMAT(i11) 
!                                                                       
END SUBROUTINE do_remove                      
!
!*****7**************************************************************   
!
SUBROUTINE do_purge (line, length)
!-                                                                      
!     Purges the list of atoms from all deleted atoms                   
!+                                                                      
USE discus_config_mod 
USE conn_mod
USE molecule_mod 
USE update_cr_dim_mod
!
USE errlist_mod 
USE get_params_mod
USE precision_mod
USE string_convert_mod
USE take_param_mod
!
IMPLICIT none 
!
INTEGER, PARAMETER :: MAXW = 1
!
CHARACTER(LEN=*), INTENT(INOUT) :: line
INTEGER            , INTENT(INOUT) :: length
!
INTEGER             :: ianz
!
CHARACTER(LEN=MAX(PREC_STRING,LEN(line))), DIMENSION(maxw) :: cpara
INTEGER            , DIMENSION(maxw) :: lpara
!REAL(KIND=PREC_DP) , DIMENSION(maxw) :: werte
!
INTEGER, PARAMETER :: NOPTIONAL = 1
INTEGER, PARAMETER :: O_TYPE    = 1
CHARACTER(LEN=   4), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=MAX(PREC_STRING,LEN(line))), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent!opt. para is present
REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 0 ! Number of values to calculate 
!
DATA oname  / 'type'   /
DATA loname /  4       /
opara  =  (/ 'no'      /)   ! Always provide fresh default values
lopara =  (/  2        /)
owerte =  (/  0.0      /)
!
CALL get_params (line, ianz, cpara, lpara, maxw, length) 
IF(ier_num /= 0) RETURN
!
CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
IF(ier_num /= 0) RETURN
!
!
IF (mole_num_mole == 0) THEN 
   CALL do_purge_atoms 
ELSEIF (mole_num_mole >  0) THEN 
   CALL do_purge_molecules_new
ENDIF 
!
CALL do_low(opara(O_TYPE))
IF(opara(O_TYPE) == 'yes') THEN
   CALL do_purge_types
ENDIF
CALL update_cr_dim
!
END SUBROUTINE do_purge                       
!
!*****7**************************************************************** 
!
      SUBROUTINE do_purge_atoms 
!-                                                                      
!     Purges the list of atoms from all deleted atoms                   
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod
      USE chem_aver_mod
      USE molecule_mod 
      USE errlist_mod 
      USE param_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER i, ii, idel, ndel 
      LOGICAL lout 
!                                                                       
      DATA lout / .false. / 
!                                                                       
      CALL chem_elem (lout) 
      ndel = nint (res_para (1) * cr_natoms) 
!                                                                       
      IF (ndel.ne.0) then 
         idel = 0 
         DO i = 1, cr_natoms - ndel 
         DO while (i + idel.le.cr_natoms.and.cr_iscat (i + idel) .eq.0) 
         idel = idel + 1 
         ENDDO 
         ii = i + idel 
         cr_pos (1, i) = cr_pos (1, ii) 
         cr_pos (2, i) = cr_pos (2, ii) 
         cr_pos (3, i) = cr_pos (3, ii) 
         cr_iscat (i) = cr_iscat (ii) 
         cr_mole (i) = cr_mole (ii) 
         cr_surf(:,i)= cr_surf(:,ii) 
         cr_magn(:,i)= cr_magn(:,ii) 
         cr_prop (i) = cr_prop (ii) 
         ENDDO 
         cr_natoms = cr_natoms - ndel 
         ndel = 0 
!                                                                       
         CALL do_check_purge (.true.) 
         cr_ncatoms = 1       !Atom number disturbed, 1 atom per unit cell
      ENDIF 
!                                                                       
      END SUBROUTINE do_purge_atoms                 
!*****7**************************************************************** 
      SUBROUTINE do_check_purge (lpurge) 
!-                                                                      
!     This routine checks the neighbor calculation modes. If the        
!     argument is true, order is messed up and we go to exact           
!     modes, otherwise we select speed mode. Any changes are            
!     reported on screen.                                               
!+                                                                      
      USE discus_config_mod 
      USE chem_mod 
      USE pdf_mod 
!
      USE errlist_mod 
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      LOGICAL lpurge 
!                                                                       
      IF (lpurge) THEN 
         WRITE (output_io, 1000) 
!        IF(chem_period(1) .OR. chem_period(2) .OR. chem_period(3) ) THEN
            chem_period (1) = .false. 
            chem_period (2) = .false. 
            chem_period (3) = .false. 
            chem_purge      = .TRUE.      ! Disable use of periodic boundary
            WRITE (output_io, 1100) 'Peridic boundaries DISABLED ..' 
!        ENDIF 
         IF (chem_quick) THEN 
            chem_quick = .false. 
      WRITE (output_io, 1100) 'Chem. neighbor mode set to EXACT ..' 
         ENDIF 
         IF (.not.pdf_lexact) THEN 
            pdf_lexact = .true. 
      WRITE (output_io, 1100) 'PDF calculation mode set to EXACT ..' 
         ENDIF 
      ELSE 
!       To be added ..                                                  
      ENDIF 
!                                                                       
 1000 FORMAT (1x,'Atom order disturbed. Checking calulation modes ..') 
 1100 FORMAT (3x,a) 
!                                                                       
      END SUBROUTINE do_check_purge                 
!*****7**************************************************************** 
      SUBROUTINE do_purge_molecules_new
!-                                                                      
!     Purges the list of atoms from all deleted atoms                   
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE chem_aver_mod
      USE molecule_mod 
      USE errlist_mod 
      USE param_mod 
      IMPLICIT none 
!                                                                       
      INTEGER :: ndel
      INTEGER :: ia                   ! Loop over atoms
      INTEGER :: im                   ! Dummy index for molecule
      INTEGER :: inew                 ! Dummy index for new molecule
      INTEGER :: old_mole_num_mole    ! Original molecule number
      INTEGER :: old_mole_num_type    ! Original molecule type number
      INTEGER :: max_length           ! Original maximum molecule length
      LOGICAL :: lout = .false.       ! no output
!
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: new_mole ! Temporary moleceule properties
      INTEGER, DIMENSION(:  ), ALLOCATABLE :: new_len
      INTEGER, DIMENSION(:  ), ALLOCATABLE :: new_type
      CHARACTER (LEN=200), DIMENSION(:  ), ALLOCATABLE :: new_file
      INTEGER, DIMENSION(:  ), ALLOCATABLE :: new_char
      REAL   , DIMENSION(:  ), ALLOCATABLE :: new_dens
      REAL   , DIMENSION(:  ), ALLOCATABLE :: new_biso
      REAL   , DIMENSION(:  ), ALLOCATABLE :: new_clin
      REAL   , DIMENSION(:  ), ALLOCATABLE :: new_cqua
      REAL   , DIMENSION(:  ), ALLOCATABLE :: new_fuzz
!                                                                       
!                                                                       
      CALL chem_elem (lout) 
      ndel = nint (res_para (1) * cr_natoms) 
!
      old_mole_num_mole = mole_num_mole
      old_mole_num_type = mole_num_type
      max_length        = MAXVAL(mole_len)
!                                                                       
      IF (ndel.ne.0) then 
         CALL do_purge_atoms        ! First remove atoms
!
         ALLOCATE(new_mole(0:old_mole_num_mole,1:max_length))
         ALLOCATE(new_len (0:old_mole_num_mole))
         ALLOCATE(new_type(0:old_mole_num_mole))
         ALLOCATE(new_file(0:old_mole_num_mole))
         ALLOCATE(new_char(0:old_mole_num_mole))
         ALLOCATE(new_dens(0:old_mole_num_mole))
         ALLOCATE(new_biso(0:old_mole_num_type))
         ALLOCATE(new_clin(0:old_mole_num_type))
         ALLOCATE(new_cqua(0:old_mole_num_type))
         ALLOCATE(new_fuzz(0:old_mole_num_mole))
!
         new_mole(:,:) = 0                 ! Initialise all arrays
         new_len (:)   = 0
         new_type(:)   = 0
         new_file(:)   = ' '
         new_char(:)   = 0
         new_dens(:)   = 0
         new_biso(:)   = 0.0
         new_clin(:)   = 0.0
         new_cqua(:)   = 0.0
         new_fuzz(:)   = 0.0
!
         DO ia=1, cr_natoms           ! Loop over all atoms
            IF(cr_mole(ia)/=0) THEN   ! Atom is in a molecule
               im = cr_mole(ia)
               new_len (im)             = new_len (im) + 1 ! increment length
               new_mole(im,new_len(im)) = ia               ! insert atom
            ENDIF
         ENDDO
!
         new_type(0:old_mole_num_mole) = mole_type(0:old_mole_num_mole)         ! Copy all molecule properties
         new_file(0:old_mole_num_mole) = mole_file(0:old_mole_num_mole)
         new_char(0:old_mole_num_mole) = mole_char(0:old_mole_num_mole)
         new_dens(0:old_mole_num_mole) = mole_dens(0:old_mole_num_mole)
         new_biso(0:old_mole_num_type) = mole_biso(0:old_mole_num_type)
         new_clin(0:old_mole_num_type) = mole_clin(0:old_mole_num_type)
         new_cqua(0:old_mole_num_type) = mole_cqua(0:old_mole_num_type)
         new_fuzz(0:old_mole_num_mole) = mole_fuzzy(0:old_mole_num_mole)
!
         mole_len  = 0               ! Clear old molecules
         mole_off  = 0
         mole_type = 0
         mole_cont = 0
         mole_char = 0
         mole_dens = 0
         mole_biso = 0
         mole_clin = 0
         mole_cqua = 0
         mole_file = ' '
!
         inew = 0                    ! No new molecules yet
         DO im=1,mole_num_mole       ! Loop over all old molecules
            IF(new_len(im)>0) THEN   ! This molecule still has atoms
               inew = inew + 1
               mole_len  (inew) = new_len (im)  ! Copy molecule properties
               mole_type (inew) = new_type(im)
               mole_file (inew) = new_file(im)
               mole_char (inew) = new_char(im)
               mole_dens (inew) = new_dens(im)
               mole_biso (mole_type(inew)) = new_biso(new_type(im))
               mole_clin (mole_type(inew)) = new_clin(new_type(im))
               mole_cqua (mole_type(inew)) = new_cqua(new_type(im))
               mole_fuzzy(inew) = new_fuzz(im)
               mole_off  (inew) = mole_off(inew-1) + mole_len(inew-1)
               DO ia=1,new_len(im)   ! Loop over atoms in this molecule
                  mole_cont(mole_off(inew)+ia) = new_mole(im,ia) ! Copy into molecule
                  cr_mole(new_mole(im,ia))     = inew            ! update atom property
               ENDDO
            ENDIF
         ENDDO
         mole_num_mole = inew   ! Update molecule (numbers, types, atom_number)
         mole_num_type = MAXVAL(mole_type)
         mole_num_atom = mole_off(inew) + mole_len(inew)
         mole_num_unit = 1  ! Atom number disturbed, 1 mol per unit cell
!
         DEALLOCATE(new_mole)    ! Free teporary memory
         DEALLOCATE(new_len )
         DEALLOCATE(new_type)
         DEALLOCATE(new_file)
         DEALLOCATE(new_char)
         DEALLOCATE(new_dens)
         DEALLOCATE(new_biso)
         DEALLOCATE(new_clin)
         DEALLOCATE(new_cqua)
         DEALLOCATE(new_fuzz)
!
      ENDIF 
!                                                                       
      END SUBROUTINE do_purge_molecules_new
!*****7**************************************************************** 
      SUBROUTINE do_purge_molecules 
!-                                                                      
!     Purges the list of atoms from all deleted atoms                   
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE chem_aver_mod
      USE molecule_mod 
      USE errlist_mod 
      USE param_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER i, ii, idel, ndel, mdel, k
      INTEGER nmol 
      INTEGER mmol 
      INTEGER iatom 
      LOGICAL lmolecule_deleted 
      LOGICAL lout 
!                                                                       
      DATA lout / .false. / 
!                                                                       
      CALL chem_elem (lout) 
      ndel = nint (res_para (1) * cr_natoms) 
      lmolecule_deleted = .false. 
!                                                                       
      IF (ndel.ne.0) then 
         idel = 0 
         i = 1 
         DO while (i.le.cr_natoms - ndel) 
         IF (cr_iscat (i) .eq.0) then 
!                                                                       
!     ------we have a deleted atom, check for its presence in a         
!     ------ molecule. If found, set reference to atom no to zero.      
!     -------Shift all higher atom numbers one down                     
!                                                                       
            nmol = cr_mole(i)                ! Get molecule no.
            IF(nmol>0) THEN                  ! Atom is in a molecule
            DO nmol = 1, mole_num_mole 
            DO iatom = 1, mole_len (nmol) 
            k = mole_cont (mole_off (nmol) + iatom) 
            IF (i.eq.k) then 
               mole_cont (mole_off (nmol) + iatom) = 0 
               lmolecule_deleted = .true. 
            ELSEIF (i.lt.k) then 
               mole_cont (mole_off (nmol) + iatom) = mole_cont (        &
               mole_off (nmol) + iatom) - 1                             
            ENDIF 
            ENDDO 
            ENDDO 
            ENDIF   ! Atom is in a molecule
!                                                                       
!     ------shift all following atoms one down                          
!                                                                       
            DO ii = i, cr_natoms - idel 
            cr_pos (1, ii) = cr_pos (1, ii + 1) 
            cr_pos (2, ii) = cr_pos (2, ii + 1) 
            cr_pos (3, ii) = cr_pos (3, ii + 1) 
            cr_iscat (ii) = cr_iscat (ii + 1) 
            cr_mole (ii) = cr_mole (ii + 1) 
            cr_surf(:,ii)= cr_surf (:,ii + 1) 
            cr_magn(:,ii)= cr_magn (:,ii + 1) 
            cr_prop (ii) = cr_prop (ii + 1) 
            ENDDO 
            idel = idel + 1 
         ELSE 
            i = i + 1 
         ENDIF 
         ENDDO 
         cr_natoms = cr_natoms - ndel 
         ndel = 0 
!                                                                       
!     --clean up references to atom no. "0" in molecule list            
!       Each time a "0" is found, all following entries are shifted     
!       down by one                                                     
!                                                                       
         IF (lmolecule_deleted) then 
            idel = 0 
            DO nmol = 1, mole_num_mole 
            iatom = 1 
            mdel = 0 
            DO while (iatom.le.mole_len (nmol) - mdel) 
            k = mole_cont (mole_off (nmol) + iatom) 
            DO while (iatom + mdel.le.mole_len (nmol) .and.mole_cont (  &
            mole_off (nmol) + iatom + mdel) .eq.0)                      
            mdel = mdel + 1 
            ENDDO 
            mole_cont (mole_off (nmol) + iatom - idel) = mole_cont (    &
            mole_off (nmol) + iatom + mdel)                             
            iatom = iatom + 1 
            ENDDO 
            mole_len (nmol) = mole_len (nmol) - mdel 
            mole_off (nmol) = mole_off (nmol) - idel 
            idel = idel + mdel 
            ENDDO 
!                                                                       
!     --Check if the complete content of a molecule was purged          
!                                                                       
            nmol = 1 
            DO while (nmol.le.mole_num_mole) 
            IF (mole_len (nmol) .eq.0) then 
               DO mmol = nmol, mole_num_mole-1 
               mole_len (mmol) = mole_len (mmol + 1) 
               mole_off (mmol) = mole_off (mmol + 1) 
               mole_type (mmol) = mole_type (mmol + 1) 
               mole_file (mmol) = mole_file (mmol + 1) 
               mole_char (mmol) = mole_char (mmol + 1) 
               ENDDO 
               mole_num_mole = mole_num_mole-1 
            ELSE 
               nmol = nmol + 1 
            ENDIF 
            ENDDO 
         ENDIF 
         WRITE ( *, 1000) char (7) 
      ENDIF 
!                                                                       
 1000 FORMAT(                                                           &
     & ' ****APPL****All functions that require identical     ****'/    &
     & '         ****number of atoms in ALL unit cells might  ****'/    &
     & '         ****not work correctly after a purge command ****',    &
     & a1)                                                              
!                                                                       
      END SUBROUTINE do_purge_molecules             
!
!*****7**************************************************************** 
!
SUBROUTINE do_purge_types 
!-                                                                      
!     Purges the list of atom typess from all deleted atoms                   
!   If no atoms are present for an atom type this type is purged
!   All atom type move down in the list
!+                                                                      
USE discus_config_mod 
USE crystal_mod
!
IMPLICIT NONE
!
INTEGER, DIMENSION(:), ALLOCATABLE :: n_atom_type
INTEGER, DIMENSION(:), ALLOCATABLE :: new_type
INTEGER :: i, j     ! dummy loop indices
INTEGER :: ndel     ! Number of types to be deleted
!
ALLOCATE(n_atom_type(0:MAXSCAT))
ALLOCATE(new_type(0:MAXSCAT))
!
n_atom_type(:) = 0
!
!  Accumulate number of atoms per type 
!
DO i=1, cr_natoms
   j = cr_iscat(i)
   n_atom_type(j) = n_atom_type(j) + 1
ENDDO
!
! Build lookup table if types are to be deleted
!
ndel = 0
DO i=1, cr_nscat
   IF(n_atom_type(i)==0) THEN
      ndel = ndel + 1
      new_type(i) = 0
   ELSE
      new_type(i) = i - ndel     ! record new location
   ENDIF
ENDDO
!
IF(ndel>0) THEN
   DO i=1, cr_natoms
      cr_iscat(i) = new_type(cr_iscat(i))
   ENDDO
   DO i=1, cr_nscat
      IF(new_type(i) /= 0) THEN
         cr_scat    (:,new_type(i)) = cr_scat    (:,i)
         cr_delfr   (  new_type(i)) = cr_delfr   (  i)
         cr_delfi   (  new_type(i)) = cr_delfi   (  i)
         cr_scat_int(  new_type(i)) = cr_scat_int(  i)
         cr_scat_equ(  new_type(i)) = cr_scat_equ(  i)
         cr_delf_int(  new_type(i)) = cr_delf_int(  i)
         cr_at_lis  (  new_type(i)) = cr_at_lis  (  i)
         cr_at_equ  (  new_type(i)) = cr_at_equ  (  i)
         as_at_lis  (  new_type(i)) = as_at_lis  (  i)
         as_iscat   (  new_type(i)) = as_iscat   (  i)
         as_mole    (  new_type(i)) = as_mole    (  i)
         as_prop    (  new_type(i)) = as_prop    (  i)
         cr_dw      (  new_type(i)) = cr_dw      (  i)
         cr_occ     (  new_type(i)) = cr_occ     (  i)
         as_dw      (  new_type(i)) = as_dw      (  i)
         as_occ     (  new_type(i)) = as_occ     (  i)
         as_pos     (:,new_type(i)) = as_pos     (:,i)
      ENDIF
   ENDDO
   cr_nscat = cr_nscat - ndel
ENDIF
!
DEALLOCATE(n_atom_type)
DEALLOCATE(new_type)
!
END SUBROUTINE do_purge_types
!
!*****7**************************************************************** 
      SUBROUTINE do_copy (line, lp) 
!-                                                                      
!     copies an atom to a new position given in relativ or absolute     
!     fractional coordinates                                            
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
USE precision_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 5) 
!                                                                       
      CHARACTER ( * ) line 
      CHARACTER(LEN=1) :: mode
      CHARACTER(80) zeile 
      CHARACTER(LEN=MAX(PREC_STRING,LEN(line))) :: cpara (maxw) 
      INTEGER lpara (maxw) 
      INTEGER ianz, i, lp, ind 
      REAL(KIND=PREC_DP) :: werte (maxw) 
!                                                                       
      zeile = ' ' 
      CALL get_params (line, ianz, cpara, lpara, maxw, lp) 
      IF (ianz.ne.5) then 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
      IF (ier_num.eq.0) then 
         mode = cpara(1)(1:1)      
         cpara(1) = '0'
         lpara(1) = 1
         CALL ber_params(ianz, cpara, lpara, werte, maxw)
!         DO i = 2, ianz 
!         lp = lpara (i) 
!         zeile = ' ' 
!         zeile (1:1) = '(' 
!         zeile (2:lp + 1) = cpara (i) (1:lp) 
!         zeile (lp + 2:lp + 2) = ')' 
!         lp = lp + 2 
!         werte (i) = berechne (zeile, lp) 
         IF (ier_num.ne.0) return 
!         ENDDO 
         ind = int (werte (2) ) 
         IF (0.lt.ind.and.ind.le.cr_natoms.and.ind.le.NMAX) then 
            zeile = ' ' 
            zeile (1:4) = cr_at_lis (cr_iscat (ind) ) 
            zeile (5:5) = ',' 
            IF (mode            ==  'a') then 
               WRITE (zeile (7:58), 3000) (werte (i), i = 3, 5),        &
               cr_dw (cr_iscat (ind) )                                  
               lp = 58 
            ELSEIF (mode            ==  'r') then 
               WRITE (zeile (7:58), 3000) (werte (i) + cr_pos (i - 2,   &
               ind), i = 3, 5), cr_dw (cr_iscat (ind) )                 
               lp = 58 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
               RETURN 
            ENDIF 
            CALL do_ins (zeile, lp) 
         ELSE 
            ier_num = - 19 
            ier_typ = ER_APPL 
         ENDIF 
      ENDIF 
 3000 FORMAT    (3(e12.5e1,','),e13.5e3) 
      END SUBROUTINE do_copy                        
!*****7*****************************************************************
      SUBROUTINE do_switch (line, lp) 
!-                                                                      
!     Switches two atoms in the structure. The number of the two        
!     atoms must be given. If the last parameter is "mol" then          
!     two molecules are swapped.                                        
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE molecule_mod 
      USE update_cr_dim_mod
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
USE precision_mod
USE str_comp_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 3) 
!                                                                       
      CHARACTER ( * ) line 
      CHARACTER(LEN=MAX(PREC_STRING,LEN(line))) :: cpara (maxw) 
      REAL(KIND=PREC_DP) :: werte (maxw) 
      INTEGER lpara (maxw) 
      INTEGER i, j, ianz, lp, is 
      INTEGER, DIMENSION(0:3) :: iis
      REAL   , DIMENSION(0:3) :: ris
!                                                                       
!                                                                       
      CALL get_params (line, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
!                                                                       
!------ Molecules                                                       
!                                                                       
      IF (ianz.eq.3) then 
         IF (str_comp (cpara (3) , 'mol', 2, lpara (3) , 3) ) then 
            ianz = ianz - 1 
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
            i = nint (werte (1) ) 
            j = nint (werte (2) ) 
      IF (0.lt.i.and.i.le.mole_num_mole.and.0.lt.j.and.j.le.mole_num_mol&
     &e) then                                                           
               CALL do_swap_mole (i, j, .true.) 
               CALL update_cr_dim 
            ELSE 
               ier_num = - 63 
               ier_typ = ER_APPL 
            ENDIF 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
!                                                                       
!------ Atoms                                                           
!                                                                       
      ELSEIF (ianz.eq.2) then 
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         i = nint (werte (1) ) 
         j = nint (werte (2) ) 
         IF (0.lt.i.and.i.le.cr_natoms.and.0.lt.j.and.j.le.cr_natoms)   &
         then                                                           
            is           = cr_iscat (i) 
            cr_iscat (i) = cr_iscat (j) 
            cr_iscat (j) = is 
            is           = cr_prop (i) 
            cr_prop (i)  = cr_prop (j) 
            cr_prop (j)  = is 
            is           = cr_mole (i) 
            cr_mole (i)  = cr_mole (j) 
            cr_mole (j)  = is 
            iis(:)         = cr_surf (:,i)
            cr_surf (:,i)  = cr_surf (:,j)
            cr_surf (:,j)  = iis(:)
            ris(:)         = cr_magn (:,i)
            cr_magn (:,i)  = cr_magn (:,j)
            cr_magn (:,j)  = ris(:)
         ELSE 
            ier_num = - 19 
            ier_typ = ER_APPL 
         ENDIF 
!                                                                       
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
      END SUBROUTINE do_switch                      
!
!*****7*****************************************************************
!
      LOGICAL FUNCTION scat_allowed (is, werte, ianz, maxw) 
!+                                                                      
!     checks if atom i is within the selected atom range in             
!     werte(ianz).                                                      
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE errlist_mod 
USE precision_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER is, ianz, maxw 
      REAL(KIND=PREC_DP) :: werte (maxw) 
!                                                                       
      INTEGER j 
      LOGICAL ltype 
!                                                                       
      IF(NINT(werte(1)) == -1) then 
         ltype = .true. 
      ELSE 
         ltype = .false. 
         DO j = 1, ianz 
         ltype = ltype.or.is.eq.nint (werte (j) ) 
         ENDDO 
      ENDIF 
      scat_allowed = ltype 
      END FUNCTION scat_allowed                     
!*****7*****************************************************************
!
SUBROUTINE atom_select (zeile, lp, lu, lo, latom, &
                        lsite, lus, los,          &
                        sel_atom, lold, lselect,  &
                        ival, repl)                             
!+                                                                      
!     This routine executes the select command                          
!-                                                                      
USE discus_config_mod 
USE crystal_mod 
USE get_iscat_mod
!
USE ber_params_mod
USE errlist_mod 
USE get_params_mod
USE precision_mod
USE take_param_mod
!                                                                       
IMPLICIT none 
!                                                                       
CHARACTER  (LEN=  * ),     INTENT(IN)    :: zeile 
INTEGER,                   INTENT(INOUT) :: lp 
INTEGER,                   INTENT(IN)    :: lu
INTEGER,                   INTENT(IN)    :: lo
LOGICAL, DIMENSION(lu:lo), INTENT(OUT)   :: latom 
INTEGER,                   INTENT(IN)    :: lus
INTEGER,                   INTENT(IN)    :: los
LOGICAL, DIMENSION(lus:los), INTENT(OUT) :: lsite
LOGICAL,                   INTENT(INOUT) :: sel_atom
LOGICAL,                   INTENT(IN)    :: lold
LOGICAL,                   INTENT(IN)    :: lselect 
INTEGER,                   OPTIONAL,  INTENT(IN)    :: ival
INTEGER, DIMENSION(lu:lo), OPTIONAL,  INTENT(OUT)   :: repl (lu:lo)
!
INTEGER  :: maxw
!                                                                       
CHARACTER(LEN=MAX(PREC_STRING,LEN(zeile))), DIMENSION(1:lo+1) :: cpara 
REAL(KIND=PREC_DP) , DIMENSION(1:lo+1) :: werte 
INTEGER            , DIMENSION(1:lo+1) :: lpara 
INTEGER                                :: ianz, i, is
CHARACTER(LEN=PREC_STRING)                    :: line
!
INTEGER, PARAMETER :: NOPTIONAL = 1
INTEGER, PARAMETER :: O_SITE    = 1
CHARACTER(LEN=4), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
CHARACTER(LEN=MAX(PREC_STRING,LEN(zeile))), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent!opt. para is present
REAL(KIND=PREC_DP) , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
INTEGER, PARAMETER                        :: ncalc = 0 ! Number of values to calculate 
!
DATA oname  / 'site'   /
DATA loname /  4       /
opara  =  (/ '0.0000'  /)   ! Always provide fresh default values
lopara =  (/  6        /)
owerte =  (/  0.0      /)
!
!
!                                                                       
maxw = lo+1
CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
IF (ier_num.ne.0) RETURN 
!
CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                  oname, loname, opara, lopara, lpresent, owerte)
!                                                                       
!------ In case we had molecules before, deselect all                   
!                                                                       
IF (.NOT.sel_atom) THEN 
   latom = .false. !  latom (i) = .false. 
   IF(PRESENT(repl)) repl = 0 !  repl (i) = ival 
ENDIF 
!                                                                       
sel_atom = .true. 
!                                                                       
!------ Select/deselect atoms                                           
!                                                                       
CALL get_iscat (ianz, cpara, lpara, werte, maxw, lold) 
IF (ier_num.ne.0) RETURN 
!                                                                       
IF(NINT(werte(1)) == -1) THEN   ! all atoms are selected
   latom  = lselect !  latom (i) = lselect 
   IF(PRESENT(repl)) repl = ival !  repl (i) = ival 
ELSE 
   DO i = 1, ianz 
      is = nint (werte (i) ) 
      IF (is.ge.0.and.is.le.lo) THEN 
         latom (is) = lselect 
         IF(PRESENT(repl)) repl (is) = ival 
      ELSE 
         ier_num = - 27 
         ier_typ = ER_APPL 
      ENDIF 
   ENDDO 
ENDIF 
!
lsite = .TRUE.
!
IF(lpresent(O_SITE)) THEN      !optional parameter 'site:' is present
   IF(opara(O_SITE)=='all') THEN
      lsite = .TRUE.
   ELSE
      line(1:lopara(O_SITE)-2) = opara(O_SITE)(2:lopara(O_SITE)-1)
      lp = lopara(O_SITE)-2
      CALL get_params(line, ianz, cpara, lpara, MAXW, lp)
      IF(ier_num/=0) RETURN
      CALL ber_params(ianz, cpara, lpara, werte, MAXW)
      IF(ier_num/=0) RETURN
!
      lsite = .FALSE.
      DO i = 1, ianz 
         is = NINT(werte(i))
         IF(is >= 0 .AND. is<=cr_ncatoms) THEN 
            lsite(is) = .TRUE.
         ELSE 
            ier_num = -10
            ier_typ = ER_CHEM 
         ENDIF
      ENDDO
   ENDIF
ENDIF
!                                                                       
END SUBROUTINE atom_select                    
!
!*****7*****************************************************************
      SUBROUTINE mole_select (zeile, lp, lu, lo, latom, &
                              sel_atom, lselect,  &
                              ival, repl)                             
!+                                                                      
!     This routine exectues the select command                          
!-                                                                      
      USE discus_config_mod 
      USE molecule_mod 
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
USE precision_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      CHARACTER  (LEN=  * ),     INTENT(IN)    :: zeile 
      INTEGER,                   INTENT(INOUT) :: lp 
      INTEGER,                   INTENT(IN)    :: lu
      INTEGER,                   INTENT(IN)    :: lo
      LOGICAL, DIMENSION(lu:lo), INTENT(OUT)   :: latom 
      LOGICAL,                   INTENT(INOUT) :: sel_atom
      LOGICAL,                   INTENT(IN)    :: lselect 
      INTEGER,                   OPTIONAL,  INTENT(IN)    :: ival
      INTEGER, DIMENSION(lu:lo), OPTIONAL,  INTENT(OUT)   :: repl (lu:lo)
!
      INTEGER                                :: maxw
!                                                                       
      CHARACTER(LEN=MAX(PREC_STRING,LEN(zeile))), DIMENSION(1:lo+1) :: cpara 
      REAL(KIND=PREC_DP) , DIMENSION(1:lo+1) :: werte 
      INTEGER            , DIMENSION(1:lo+1) :: lpara 
!                                                                       
      INTEGER                                :: ianz, i, is
!                                                                       
      maxw = lo+1
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
!                                                                       
!------ In case we had atoms before, deselect all                       
!                                                                       
      IF (sel_atom) THEN 
         latom = .false.  !   latom (i) = .false. 
         IF(PRESENT(repl)) repl = 0   !repl (i) = ival 
      ENDIF 
!                                                                       
      sel_atom = .false. 
!                                                                       
!------ Select/deselect molecules                                       
!                                                                       
      IF (cpara (1) (1:1) .eq.'a') THEN 
         latom = lselect  !   latom (i) = .false. 
         IF(PRESENT(repl)) repl = ival   !repl (i) = ival 
      ELSE 
         CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         DO i = 1, ianz 
            is = nint (werte (i) ) 
            IF (is.gt.0.and.is.le.mole_num_type) THEN 
               latom (is) = lselect 
               IF(PRESENT(repl)) repl (is) = ival 
            ELSE 
               ier_num = - 64 
               ier_typ = ER_APPL 
            ENDIF 
         ENDDO 
      ENDIF 
!                                                                       
      END SUBROUTINE mole_select                    
!
!*****7*****************************************************************
      SUBROUTINE do_change (line, laenge) 
!-                                                                      
!     These subroutine is the main routine for changing various         
!     parameters.                                                       
!                                                                       
!     Version  : 1.0                                                    
!     Date     : 03 Jul 96                                              
!                                                                       
!     Author   : R.B. Neder (reinhard.neder@fau.de)      
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
!                                                                       
      USE param_mod 
      USE prompt_mod 
      USE errlist_mod 
      USE get_params_mod
USE precision_mod
USE str_comp_mod
      IMPLICIT none 
       
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 7) 
!                                                                       
      CHARACTER ( * ) line 
      CHARACTER(LEN=MAX(PREC_STRING,LEN(line))) :: cpara (maxw) 
      INTEGER lpara (maxw) 
      INTEGER ianz 
      INTEGER laenge 
      REAL(KIND=PREC_DP) :: werte (maxw) 
!                                                                       
!                                                                       
      CALL get_params (line, ianz, cpara, lpara, maxw, laenge) 
      IF (ier_num.eq.0) then 
!                                                                       
!     --Interprete first parameter as command                           
!                                                                       
!                                                                       
!     ----Change object properties                                      
!                                                                       
         IF (str_comp (cpara (1) , 'domain', 2, lpara (1) , 6) ) then 
            CALL del_params (1, ianz, cpara, lpara, maxw) 
            CALL change_object (cpara, lpara, werte, ianz, maxw,        &
            - 1)                                                        
         ELSEIF (str_comp (cpara (1) , 'object', 2, lpara (1) , 6) )    &
         then                                                           
            CALL del_params (1, ianz, cpara, lpara, maxw) 
            CALL change_object (cpara, lpara, werte, ianz, maxw, 1) 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
      ENDIF 
!                                                                       
      END SUBROUTINE do_change                      
!*****7*****************************************************************
      SUBROUTINE change_object (cpara, lpara, werte, ianz, maxw, itype) 
!-                                                                      
!     Allows changes of extended objects                                
!                                                                       
!     Version  : 1.0                                                    
!     Date     : 16 Mar 03                                              
!                                                                       
!     Author   : R.B. Neder (reinhard.neder@fau.de)      
!+                                                                      
!                                                                       
      USE discus_config_mod 
      USE crystal_mod 
      USE molecule_mod 
!                                                                       
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
      USE build_name_mod
USE precision_mod
      USE prompt_mod 
USE str_comp_mod
      IMPLICIT none 
       
!                                                                       
      INTEGER maxw 
!                                                                       
      CHARACTER ( * ) cpara (maxw) 
!                                                                       
      INTEGER lpara (maxw) 
      INTEGER ianz 
      INTEGER iianz 
      INTEGER itype 
      INTEGER number 
!                                                                       
      REAL(KIND=PREC_DP) :: werte (maxw) 
!                                                                       
!                                                                       
      IF (str_comp (cpara (1) , 'character', 2, lpara (1) , 6) ) then 
         IF (ianz.eq.3) then 
            CALL del_params (1, ianz, cpara, lpara, maxw) 
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
            iianz = 1 
            CALL ber_params (iianz, cpara, lpara, werte, maxw) 
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
            number = nint (werte (1) ) 
            IF (0.lt.number.and.number.le.mole_num_mole) then 
               IF (itype.gt.0.and.mole_char (number) .ge.0) then 
                  IF (str_comp (cpara (2) , 'atom', 2, lpara (2) , 4) ) &
                  then                                                  
                     mole_char (number) = MOLE_ATOM 
                  ELSEIF (str_comp (cpara (2) , 'cube', 2, lpara (2) ,  &
                  4) ) then                                             
                     mole_char (number) = MOLE_CUBE 
                  ELSEIF (str_comp (cpara (2) , 'cylinder', 2, lpara (2)&
                  , 8) ) then                                           
                     mole_char (number) = MOLE_CYLINDER 
                  ELSEIF (str_comp (cpara (2) , 'sphere', 2, lpara (2) ,&
                  6) ) then                                             
                     mole_char (number) = MOLE_SPHERE 
                  ELSE 
                     ier_num = - 82 
                     ier_typ = ER_APPL 
                  ENDIF 
               ELSEIF (itype.lt.0.and.mole_char (number) .lt.0) then 
                  IF (str_comp (cpara (2) , 'atom', 2, lpara (2) , 4) ) &
                  then                                                  
                     mole_char (number) = MOLE_ATOM 
                  ELSEIF (str_comp (cpara (2) , 'cube', 2, lpara (2) ,  &
                  4) ) then                                             
                     mole_char (number) = MOLE_DOM_CUBE 
                  ELSEIF (str_comp (cpara (2) , 'cylinder', 2, lpara (2)&
                  , 8) ) then                                           
                     mole_char (number) = MOLE_DOM_CYLINDER 
                  ELSEIF (str_comp (cpara (2) , 'sphere', 2, lpara (2) ,&
                  6) ) then                                             
                     mole_char (number) = MOLE_DOM_SPHERE 
                  ELSEIF (str_comp (cpara (2) , 'fuzzy', 2, lpara (2) , &
                  5) ) then                                             
                     mole_char (number) = MOLE_DOM_FUZZY 
                  ELSE 
                     ier_num = - 82 
                     ier_typ = ER_APPL 
                  ENDIF 
               ENDIF 
            ENDIF 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
      ELSEIF (str_comp (cpara (1) , 'density', 2, lpara (1) , 7) ) then 
         IF (ianz.eq.3) then 
            CALL del_params (1, ianz, cpara, lpara, maxw) 
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
            number = nint (werte (1) ) 
            IF (0.lt.number.and.number.le.mole_num_mole) then 
               IF (itype.gt.0.and.mole_char (number) .ge.0) then 
                  mole_fuzzy (number) = werte (2) 
               ENDIF 
            ENDIF 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
      ELSEIF (str_comp (cpara (1) , 'file', 2, lpara (1) , 4) ) then 
         IF (ianz.ge.3) then 
            CALL del_params (1, ianz, cpara, lpara, maxw) 
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
            iianz = 1 
            CALL ber_params (iianz, cpara, lpara, werte, maxw) 
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
            number = nint (werte (1) ) 
            IF (0.lt.number.and.number.le.mole_num_mole) then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
               IF (itype.lt.0.and.mole_char (number) .le.0) then 
                  mole_file (number) = cpara (1) 
               ELSE 
                  ier_num = - 999 
                  ier_typ = ER_APPL 
               ENDIF 
            ENDIF 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
      ELSEIF (str_comp (cpara (1) , 'fuzzy', 2, lpara (1) , 5) ) then 
         IF (ianz.eq.3) then 
            CALL del_params (1, ianz, cpara, lpara, maxw) 
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
            number = nint (werte (1) ) 
            IF (0.lt.number.and.number.le.mole_num_mole) then 
               IF (itype.lt.0.and.mole_char (number) .le.0) then 
                  mole_fuzzy (number) = werte (2) 
               ENDIF 
            ENDIF 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
      ENDIF 
!                                                                       
      END SUBROUTINE change_object                  
!*****7*****************************************************************
      SUBROUTINE copy_mole_char (idest, isource) 
!-                                                                      
!     Copy all (molecule, object, domain properties                     
!                                                                       
!     Version  : 1.0                                                    
!     Date     : 16 Sep 04                                              
!                                                                       
!     Author   : R.B. Neder (reinhard.neder@fau.de)      
!+                                                                      
!                                                                       
      USE discus_config_mod 
      USE molecule_mod 
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER idest 
      INTEGER isource 
!                                                                       
!-----      ----Copy all properties                                     
!                                                                       
      mole_char (idest) = mole_char (isource) 
      mole_dens (idest) = mole_dens (isource) 
!      mole_biso (idest) = mole_biso (isource) 
      mole_file (idest) = mole_file (isource) 
      mole_fuzzy (idest) = mole_fuzzy (isource) 
      mole_type (idest) = mole_type (isource) 
!                                                                       
      END SUBROUTINE copy_mole_char                 
!
!*****7*****************************************************************
!
END MODULE modify_mod
