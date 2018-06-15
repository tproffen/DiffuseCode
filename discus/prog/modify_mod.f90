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
      USE modify_func_mod
      USE molecule_mod 
      USE prop_para_mod 
      USE update_cr_dim_mod
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
      USE random_mod
!
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 200) 
!                                                                       
      CHARACTER ( * ) zeile 
      CHARACTER(1024) cpara (maxw), cc 
      INTEGER lpara (maxw) , ccl
      INTEGER ianz, iianz, jjanz, lp, i 
      INTEGER is1, is2, isite, itype 
      INTEGER ja, jsite, jcell (3) 
      INTEGER            :: n_scat ! dummy for allocation
      LOGICAL lexist, lrepl 
      REAL uerte (maxw) 
      REAL verte (maxw) 
      REAL werte (maxw) 
      REAL prob 
!                                                                       
!     LOGICAL check_select_status 
      LOGICAL str_comp 
      REAL ran1 
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
                  verte (i) = float (cr_nscat) 
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
      IF (check_select_status (.true., cr_prop (ja),  cr_sel_prop) ) the&
     &n                                                                 
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
               IF (cr_dw (nint (uerte (i) ) ) .eq.werte (1) ) then 
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
      ip = cr_prop (ii) 
      jp = cr_prop (jj) 
!                                                                       
      cr_iscat (ii) = js 
      cr_mole (ii) = jm 
      cr_surf(:,ii)= jjs(:)
      cr_prop (ii) = jp 
      DO i = 1, 3 
      ipos (i) = cr_pos (i, ii) 
      cr_pos (i, ii) = cr_pos (i, jj) - j0pos (i) + i0pos (i) 
      ENDDO 
!                                                                       
      IF (lswap) then 
         cr_iscat (jj) = is 
         cr_mole (jj) = im 
         cr_surf(:,jj)= iis(:)
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
      CHARACTER(LEN=1024), DIMENSION(maxw) :: cpara
      INTEGER            , DIMENSION(maxw) :: lpara
      INTEGER i, j, ianz
      INTEGER                              :: new_nmax   = 0
      INTEGER                              :: new_nscat  = 0
      LOGICAL                              :: need_alloc = .false.
      REAL werte (maxw) 
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
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 10) 
!                                                                       
      CHARACTER ( * ) line 
      CHARACTER(4) name 
      CHARACTER(1024) cpara (maxw)
      INTEGER lpara (maxw) 
      INTEGER i, j, k, ianz, laenge
      INTEGER                              :: new_nmax
      INTEGER                              :: new_nscat
      LOGICAL                              :: need_alloc = .false.
      LOGICAL lkick, lspace 
      REAL werte (maxw) 
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
                  IF (werte (8) .lt.0) then 
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
      USE string_convert_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER, INTENT(IN)     :: maxw
!                                                                       
      CHARACTER (LEN=* )    , INTENT(INOUT) :: name 
      REAL , DIMENSION(maxw), INTENT(IN) :: werte (maxw) 
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
            lda = name.eq.cr_at_lis (i) .and.werte (5) .eq.cr_dw (i) 
            DO while (.not.lda.and.i.lt.cr_nscat) 
            i = i + 1 
            lda = name.eq.cr_at_lis (i) .and.werte (5) .eq.cr_dw (i) 
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
IMPLICIT none 
!                                                                       
       
!                                                                       
INTEGER, PARAMETER :: maxw = 6
!                                                                       
CHARACTER( LEN=* ), INTENT(INOUT) :: line 
INTEGER           , INTENT(INOUT) :: ll
!
CHARACTER(LEN=1024), DIMENSION(MAXW) :: cpara !(maxw) 
INTEGER            , DIMENSION(MAXW) :: lpara !(maxw) 
INTEGER :: i, j 
INTEGER :: istart, iend, ianz 
INTEGER :: tstart, tend 
LOGICAL :: l_mole 
LOGICAL :: l_type 
REAL               , DIMENSION(MAXW) :: werte ! (maxw) 
!                                                                       
LOGICAL :: str_comp 
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
!*****7**************************************************************   
      SUBROUTINE do_purge 
!-                                                                      
!     Purges the list of atoms from all deleted atoms                   
!+                                                                      
      USE discus_config_mod 
      USE molecule_mod 
      USE update_cr_dim_mod
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      IF (mole_num_mole == 0) THEN 
         CALL do_purge_atoms 
      ELSEIF (mole_num_mole >  0) THEN 
         CALL do_purge_molecules_new
      ENDIF 
      CALL update_cr_dim
!                                                                       
      END SUBROUTINE do_purge                       
!*****7**************************************************************** 
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
         IF(chem_period(1) .OR. chem_period(2) .OR. chem_period(3) ) THEN
            chem_period (1) = .false. 
            chem_period (2) = .false. 
            chem_period (3) = .false. 
            chem_purge      = .TRUE.      ! Disable use of periodic boundary
            WRITE (output_io, 1100) 'Peridic boundaries DISABLED ..' 
         ENDIF 
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
         ALLOCATE(new_fuzz(0:old_mole_num_mole))
!
         new_mole(:,:) = 0                 ! Initialise all arrays
         new_len (:)   = 0
         new_type(:)   = 0
         new_file(:)   = ' '
         new_char(:)   = 0
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
         new_fuzz(0:old_mole_num_mole) = mole_fuzzy(0:old_mole_num_mole)
!
         mole_len  = 0               ! Clear old molecules
         mole_off  = 0
         mole_type = 0
         mole_cont = 0
         mole_char = 0
         mole_dens = 0
         mole_biso = 0
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
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 5) 
!                                                                       
      CHARACTER ( * ) line 
      CHARACTER(LEN=1) :: mode
      CHARACTER(80) zeile 
      CHARACTER(1024) cpara (maxw) 
      INTEGER lpara (maxw) 
      INTEGER ianz, i, lp, ind 
      REAL werte (maxw) 
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
               WRITE (zeile (7:57), 3000) (werte (i), i = 3, 5),        &
               cr_dw (cr_iscat (ind) )                                  
               lp = 57 
            ELSEIF (mode            ==  'r') then 
               WRITE (zeile (7:57), 3000) (werte (i) + cr_pos (i - 2,   &
               ind), i = 3, 5), cr_dw (cr_iscat (ind) )                 
               lp = 57 
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
 3000 FORMAT    (3(e12.5e1,','),e12.5e3) 
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
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 3) 
!                                                                       
      CHARACTER ( * ) line 
      CHARACTER(1024) cpara (maxw) 
      REAL werte (maxw) 
      INTEGER lpara (maxw) 
      INTEGER i, j, ianz, lp, is 
      INTEGER, DIMENSION(0:3) :: iis
!                                                                       
      LOGICAL str_comp 
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
!*****7*****************************************************************
      SUBROUTINE do_find (line, laenge) 
!-                                                                      
!     Finds the environment around an atom                              
!+                                                                      
      USE discus_config_mod 
      USE charact_mod 
      USE charact_mod 
      USE crystal_mod 
      USE chem_mod 
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 200) 
      INTEGER mmaxw 
      PARAMETER (mmaxw = 5) 
!                                                                       
      CHARACTER ( * ) line 
      CHARACTER(1024) cpara (maxw) 
      CHARACTER(1024) ccpara (mmaxw) 
      INTEGER lpara (maxw) 
      INTEGER llpara (maxw) 
      INTEGER i, ii, ianz, iianz, laenge 
      LOGICAL lnew, fq, fp (3) 
      REAL rmin 
      REAL radius 
      REAL werte (maxw) 
      REAL wwerte (maxw) 
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
      USE celltoindex_mod
      USE modify_func_mod
      USE param_mod 
      USE errlist_mod 
      USE sorting_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER ianz, maxw 
      REAL werte (maxw) 
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
!     LOGICAL atom_allowed 
!     LOGICAL check_select_status 
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
      ltype = atom_allowed (iatom, werte, ianz, maxw) .and.check_select_&
     &status (.true., cr_prop (iatom),  cr_sel_prop)                    
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
         .and.check_select_status (.true., cr_prop (i), cr_sel_prop)    
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
             .and.check_select_status (.true., cr_prop (i), cr_sel_prop)    
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
!*****7*****************************************************************
      SUBROUTINE check_bound (cell, offset, fp, ltype) 
!+                                                                      
!     This routine applies periodic boundaries if 'fp' is TRUE.         
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER cell (3) 
      REAL offset (3) 
      LOGICAL fp (3), ltype, lok (3) 
!                                                                       
      INTEGER i 
      INTEGER j 
!                                                                       
!------ Apply periodic boundaries                                       
!                                                                       
      DO i = 1, 3 
      IF (fp (i) ) then 
         offset (i) = 0.0 
!                                                                       
         IF (cell (i) .lt.1) then 
            IF (cr_icc (i) .eq.1) then 
               offset (i) = float (cell (i) - 1) 
               cell (i) = 1 
            ELSE 
               j = ( - cell (i) / cr_icc (i) + 1) * cr_icc (i) 
               cell (i) = cell (i) + j 
               offset (i) = - float (j) 
            ENDIF 
         ELSEIF (cell (i) .gt.cr_icc (i) ) then 
            IF (cr_icc (i) .eq.1) then 
               offset (i) = float (cell (i) - 1) 
               cell (i) = 1 
            ELSE 
               j = ( (cell (i) - 1) / cr_icc (i) ) * cr_icc (i) 
               cell (i) = cell (i) - j 
               offset (i) = float (j) 
            ENDIF 
         ENDIF 
!                                                                       
         lok (i) = .true. 
!                                                                       
!------ NO periodic boundaries                                          
!                                                                       
      ELSE 
         offset (i) = 0.0 
         lok (i) = (cell (i) .gt.0) .and. (cell (i) .le.cr_icc (i) ) 
      ENDIF 
      ENDDO 
!                                                                       
      ltype = lok (1) .and.lok (2) .and.lok (3) 
!                                                                       
      END SUBROUTINE check_bound                    
!*****7*****************************************************************
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
!*****7*****************************************************************
      LOGICAL FUNCTION scat_allowed (is, werte, ianz, maxw) 
!+                                                                      
!     checks if atom i is within the selected atom range in             
!     werte(ianz).                                                      
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER is, ianz, maxw 
      REAL werte (maxw) 
!                                                                       
      INTEGER j 
      LOGICAL ltype 
!                                                                       
      IF (werte (1) .eq. - 1) then 
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
      SUBROUTINE atom_select (zeile, lp, lu, lo, latom, &
                              sel_atom, lold, lselect,  &
                              ival, repl)                             
!+                                                                      
!     This routine executes the select command                          
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE errlist_mod 
      USE get_params_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      CHARACTER  (LEN=  * ),     INTENT(IN)    :: zeile 
      INTEGER,                   INTENT(INOUT) :: lp 
      INTEGER,                   INTENT(IN)    :: lu
      INTEGER,                   INTENT(IN)    :: lo
      LOGICAL, DIMENSION(lu:lo), INTENT(OUT)   :: latom 
      LOGICAL,                   INTENT(INOUT) :: sel_atom
      LOGICAL,                   INTENT(IN)    :: lold
      LOGICAL,                   INTENT(IN)    :: lselect 
      INTEGER,                   OPTIONAL,  INTENT(IN)    :: ival
      INTEGER, DIMENSION(lu:lo), OPTIONAL,  INTENT(OUT)   :: repl (lu:lo)
!
      INTEGER  :: maxw
!                                                                       
      CHARACTER(LEN=1024), DIMENSION(1:lo+1) :: cpara 
      REAL               , DIMENSION(1:lo+1) :: werte 
      INTEGER            , DIMENSION(1:lo+1) :: lpara 
      INTEGER                                :: ianz, i, is
!                                                                       
      maxw = lo+1
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) RETURN 
!                                                                       
!------ In case we had molecules before, deselect all                   
!                                                                       
      IF (.not.sel_atom) THEN 
         latom = .false. !  latom (i) = .false. 
         IF(PRESENT(repl)) repl = 0 !  repl (i) = ival 
      ENDIF 
!                                                                       
      sel_atom = .true. 
!                                                                       
!------ Select/deselect atoms                                           
!                                                                       
      CALL get_iscat (ianz, cpara, lpara, werte, maxw, lold) 
      IF (ier_num.ne.0) return 
!                                                                       
      IF (werte (1) .eq. - 1) then   ! all atoms are selected
         latom  = lselect !  latom (i) = lselect 
         IF(PRESENT(repl)) repl = ival !  repl (i) = ival 
      ELSE 
         DO i = 1, ianz 
            is = nint (werte (i) ) 
            IF (is.ge.0.and.is.le.cr_nscat) then 
               latom (is) = lselect 
               IF(PRESENT(repl)) repl (is) = ival 
            ELSE 
               ier_num = - 27 
               ier_typ = ER_APPL 
            ENDIF 
         ENDDO 
      ENDIF 
!                                                                       
      END SUBROUTINE atom_select                    
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
      CHARACTER(LEN=1024), DIMENSION(1:lo+1) :: cpara 
      REAL               , DIMENSION(1:lo+1) :: werte 
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
!*****7*****************************************************************
      SUBROUTINE property_select (line, length, sel_mask) 
!-                                                                      
!     Sets the property bits for a select/deselect/mselecet/mdeselect   
!                                                                       
      USE prop_para_mod 
      USE errlist_mod 
      USE get_params_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      CHARACTER ( * ) line 
      INTEGER length 
      INTEGER sel_mask (0:1) 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 25) 
!                                                                       
      CHARACTER(1024) cpara (maxw) 
      INTEGER lpara (maxw) 
      INTEGER ianz 
!                                                                       
      LOGICAL str_comp 
!                                                                       
      CALL get_params (line, ianz, cpara, lpara, maxw, length) 
      IF (ier_num.ne.0) return 
      IF (ianz.lt.2) return 
!                                                                       
      IF (str_comp (cpara (1) , 'ignore', 2, lpara (1) , 6) ) then 
         CALL del_params (1, ianz, cpara, lpara, maxw) 
         CALL property_set (ianz, cpara, lpara, maxw, sel_mask, .false.,&
         .false.)                                                       
      ELSEIF (str_comp (cpara (1) , 'present', 2, lpara (1) , 7) ) then 
         CALL del_params (1, ianz, cpara, lpara, maxw) 
         CALL property_set (ianz, cpara, lpara, maxw, sel_mask, .true., &
         .true.)                                                        
      ELSEIF (str_comp (cpara (1) , 'absent', 2, lpara (1) , 6) ) then 
         CALL del_params (1, ianz, cpara, lpara, maxw) 
         CALL property_set (ianz, cpara, lpara, maxw, sel_mask, .true., &
         .false.)                                                       
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
      END SUBROUTINE property_select                
!*****7*****************************************************************
      SUBROUTINE property_set (ianz, cpara, lpara, maxw, sel_field,     &
      lentry1, lentry2)                                                 
!-                                                                      
!     Sets the property bits for a select/deselect/mselecet/mdeselect   
!                                                                       
      USE prop_para_mod 
      USE modify_func_mod
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER maxw 
      INTEGER ianz 
      CHARACTER ( * ) cpara (maxw) 
      INTEGER lpara (maxw) 
      INTEGER sel_field (0:1) 
      LOGICAL lentry1 
      LOGICAL lentry2 
!                                                                       
      INTEGER i, j 
!                                                                       
      LOGICAL str_comp 
!     INTEGER bit_set 
!                                                                       
      DO i = 1, ianz 
      IF (str_comp (cpara (i) , 'all', 2, lpara (i) , 3) ) then 
         DO j = 0, MAXPROP - 1 
         sel_field (0) = bit_set (sel_field, 0, j, lentry1) 
         sel_field (1) = bit_set (sel_field, 1, j, lentry2) 
         ENDDO 
      ELSEIF (str_comp (cpara (i) , 'normal', 2, lpara (i) , 6) ) then 
         sel_field (0) = bit_set (sel_field, 0, PROP_NORMAL, lentry1) 
         sel_field (1) = bit_set (sel_field, 1, PROP_NORMAL, lentry2) 
      ELSEIF (str_comp (cpara (i) , 'molecule', 2, lpara (i) , 8) )     &
      then                                                              
         sel_field (0) = bit_set (sel_field, 0, PROP_MOLECULE, lentry1) 
         sel_field (1) = bit_set (sel_field, 1, PROP_MOLECULE, lentry2) 
      ELSEIF (str_comp (cpara (i) , 'domain', 2, lpara (i) , 6) ) then 
         sel_field (0) = bit_set (sel_field, 0, PROP_DOMAIN, lentry1) 
         sel_field (1) = bit_set (sel_field, 1, PROP_DOMAIN, lentry2) 
      ELSEIF (str_comp (cpara (i) , 'outside', 2, lpara (i) , 7) ) then 
         sel_field (0) = bit_set (sel_field, 0, PROP_OUTSIDE, lentry1) 
         sel_field (1) = bit_set (sel_field, 1, PROP_OUTSIDE, lentry2) 
      ELSEIF (str_comp (cpara (i) , 'external', 2, lpara (i) , 8) )     &
      then                                                              
         sel_field (0) = bit_set (sel_field, 0, PROP_SURFACE_EXT,       &
         lentry1)                                                       
         sel_field (1) = bit_set (sel_field, 1, PROP_SURFACE_EXT,       &
         lentry2)                                                       
      ELSEIF (str_comp (cpara (i) , 'internal', 2, lpara (i) , 8) )     &
      then                                                              
         sel_field (0) = bit_set (sel_field, 0, PROP_SURFACE_INT,       &
         lentry1)                                                       
         sel_field (1) = bit_set (sel_field, 1, PROP_SURFACE_INT,       &
         lentry2)                                                       
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
         RETURN 
      ENDIF 
      ENDDO 
!                                                                       
      END SUBROUTINE property_set                   
!*****7*****************************************************************
      SUBROUTINE do_find_mol (ianz, werte, maxw, x, rmin, rmax) 
!                                                                       
!     This routine finds all molecules around x with a minimal          
!     distance of rmin and a maximum distance of rmax.                  
!-                                                                      
!                                                                       
      USE discus_config_mod 
      USE crystal_mod 
      USE molecule_mod 
      USE mole_env_mod 
      USE errlist_mod 
      USE param_mod 
      IMPLICIT none 
       
!                                                                       
      INTEGER ianz, maxw 
      REAL werte (maxw) 
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
!*****7*****************************************************************
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
!*****7*****************************************************************
      SUBROUTINE boundary (zeile, lp) 
!+                                                                      
!     This subroutine removes all atomes outside a given boundary.      
!-                                                                      
      USE metric_mod
      USE discus_config_mod 
      USE crystal_mod 
      USE discus_plot_mod
      USE discus_plot_init_mod
      USE point_grp
      USE prop_para_mod 
      USE surface_mod
      USE wyckoff_mod 
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
      USE take_param_mod
!
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER, PARAMETER :: maxw = 9
      INTEGER, PARAMETER :: NOPTIONAL = 6
!                                                                       
      CHARACTER (LEN=* ), INTENT(INOUT) :: zeile 
      INTEGER           , INTENT(INOUT) :: lp 
!                                                                       
!      REAL, PARAMETER :: EPS = 0.000001
      CHARACTER(1024) cpara (maxw) 
      CHARACTER(LEN=   5), DIMENSION(NOPTIONAL) :: oname   !Optional parameter names
      CHARACTER(LEN=1024), DIMENSION(NOPTIONAL) :: opara   !Optional parameter strings returned
      INTEGER            , DIMENSION(NOPTIONAL) :: loname  !Lenght opt. para name
      INTEGER            , DIMENSION(NOPTIONAL) :: lopara  !Lenght opt. para name returned
!      LOGICAL            , DIMENSION(NOPTIONAL) :: lpresent ! Optional parameter is present
      REAL               , DIMENSION(NOPTIONAL) :: owerte   ! Calculated values
      INTEGER, PARAMETER                        :: ncalc = 3 ! Number of values to calculate 
      INTEGER lpara (maxw) 
      INTEGER i, j, k, ianz 
      INTEGER :: special_form
      INTEGER :: special_n = 0
      INTEGER :: nplanes       ! number of planes that at atom is close to
      INTEGER :: iplane        ! Rindex of plane  that at atom is closest to
      LOGICAL lspace 
      LOGICAL linside 
      LOGICAL l_plane 
      LOGICAL l_sphere 
      LOGICAL l_form 
      LOGICAL l_cyl 
      LOGICAL l_ell 
      LOGICAL l_special 
      LOGICAL lwall           ! True if atom is close to cylinder wall
      LOGICAL ltop            ! True if atom is close to cylinder top
!      LOGICAL l_new 
      REAL, DIMENSION(3) :: wall  !local normal at cylinder wall
      REAL, DIMENSION(3) :: top   !local normal at cylinder top
      REAL h (3), d, dstar, radius, height , dshort
      REAL :: hkl(4)!, hklw(4)
      REAL, DIMENSION(3)    :: center       ! center of the shape
      REAL, DIMENSION(3, 2) :: special_hkl
      REAL, DIMENSION(:,:), ALLOCATABLE :: point_hkl ! (3,48)
      REAL, DIMENSION(:,:), ALLOCATABLE, SAVE :: accum_hkl ! (3,48)
      REAL, DIMENSION(:,:), ALLOCATABLE ::  temp_hkl ! (3,48)
      INTEGER               :: point_n
      INTEGER, SAVE         :: accum_n
      REAL, DIMENSION(3)    :: radius_ell
      REAL v (3) 
!      REAL, DIMENSION(3,1) :: col_vec
      REAL null (3) 
      REAL werte (maxw) 
!                                                                       
      LOGICAL str_comp 
!     REAL do_blen 
!                                                                       
      DATA null / 0.0, 0.0, 0.0 / 
      DATA oname  / 'centx', 'centy', 'centz',  'keep ',  'accum', 'exec '/
      DATA loname /  5,       5,       5,        4     ,   5     ,  4    /
      opara  =  (/ '0.0000', '0.0000', '0.0000', 'inside', 'init  ', 'run   ' /)   ! Always provide fresh default values
      lopara =  (/  6,        6,        6      ,  6      ,  4      ,  3 /)
      owerte =  (/  0.0,      0.0,      0.0    ,  0.0    ,  0.0    ,  0.0 /)
!                                                                       
      IF (cr_v.le.0.0) then 
         ier_num = - 35 
         ier_typ = ER_APPL 
         ier_msg (1) = 'A proper unit cell must be defined' 
         ier_msg (2) = 'for this command to operate       ' 
         RETURN 
      ENDIF 
!                                                                       
      l_plane  = .false. 
      l_form   = .false. 
      l_sphere = .false. 
      l_cyl    = .false. 
      l_ell    = .false. 
      linside  = .true. 
      radius   = 0.0
      height   = 0.0
      dstar    = 1.0
      center(:) = 0.0     ! Default center to 0.0, 0.0, 0.0
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF(ier_num /= 0) RETURN
      CALL get_optional(ianz, MAXW, cpara, lpara, NOPTIONAL,  ncalc, &
                        oname, loname, opara, lopara, owerte)
      IF(ier_num /= 0) RETURN
!     Handle optional and global parameters
!     IF (ier_num.eq.0) then 
         center(1:3) = owerte(1:3)   ! As defaults are always provided, we can take it take blindly here
         IF (str_comp (cpara (ianz) , 'outside', 3, lpara (ianz) , 7)) THEN
            linside = .false.
            ianz = ianz - 1
         ELSEIF (str_comp (cpara (ianz) , 'inside', 3, lpara (ianz), 6) ) THEN
            linside = .true.
            ianz = ianz - 1
         ELSE     ! Test optional parameter form
            IF (str_comp (opara(4) , 'outside', 3, lopara(4) , 7)) THEN
               linside = .false.
            ELSEIF (str_comp (opara(4) , 'inside', 3, lopara(4), 6) ) THEN
               linside = .true.
            ENDIF
         ENDIF
         IF (str_comp (cpara (1) , 'hkl',  3, lpara (1) , 3)  .OR.  &
             str_comp (cpara (1) , 'form', 3, lpara (1) , 4) ) THEN 
           l_form = str_comp (cpara (1) , 'form', 3, lpara (1) , 4)
!                                                                       
!     ----Crystal is limited by a plane or by a form of symmetrically equivalent planes
!                                                                       
            lspace = .false. 
            l_special = .FALSE.
            IF(str_comp(cpara(2),'cubeoct',7,lpara(2),7)) THEN
               IF(cr_syst==CR_CUBIC) THEN
                  cpara(5) = cpara(3)
                  lpara(5) = lpara(3)
                  cpara(2) = '1.0'
                  cpara(3) = '0.0'
                  cpara(4) = '0.0'
                  lpara(2) = 3
                  lpara(3) = 3
                  lpara(4) = 3
                  ianz = 5
!                 l_special = .TRUE.
                  special_form = 1
               ELSE
                  ier_num = -142
                  ier_typ = ER_APPL
                  ier_msg(1) = 'Unpredictable forms would result in non-cubic'
                  ier_msg(2) = 'systems. Use explicit forms and combinations'
                  ier_msg(3) = 'in non-cubic systems.'
                  RETURN
               ENDIF
            ELSE
               special_form = 0
            ENDIF
            IF (ianz.ge.4) then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               IF (ier_num.ne.0) return 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               IF (ier_num.ne.0) return 
               DO i = 1, 3 
               h (i) = werte (i) 
               ENDDO 
               dstar = do_blen (lspace, h, null) 
               IF (dstar.le.0.0) then 
                  ier_num = - 32 
                  ier_typ = ER_APPL 
                  RETURN 
               ENDIF 
               IF (ianz.eq.4) then 
                  IF ( werte(4).eq.0) THEN
                     DO i = 1, 3 
                        h (i) = h (i) / dstar *1.0E12
                     ENDDO 
                     dstar = do_blen (lspace, h, null) 
                  ELSE
                     DO i = 1, 3 
                        h (i) = h (i) / dstar / werte (4) 
                     ENDDO 
                     dstar = do_blen (lspace, h, null) 
                  ENDIF 
               ENDIF 
               IF (dstar.le.0.0) then 
                  ier_num = - 32 
                  ier_typ = ER_APPL 
                  RETURN 
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
            l_plane  = .true. 
            l_sphere = .false. 
            l_cyl    = .false. 
            l_ell    = .false. 
            IF(l_form) THEN
               l_plane= .FALSE.
               l_form = .TRUE.
            ENDIF
            IF(special_form == 0) THEN
               special_hkl(:,1) = h(:)
               special_n = 1
            ELSEIF(special_form == 1) THEN            ! CUBEOCTAHEDRON h= 1,0,0
               special_hkl(:,1) = h(:)
               special_hkl(1,2) = h(1) / 2.!/sqrt(3.)
               special_hkl(2,2) = h(1) / 2.!/sqrt(3.)
               special_hkl(3,2) = h(1) / 2.!/sqrt(3.)
               special_n = 2
            ENDIF
!
!           Handle accumulation of flat surfaces
!
            IF(str_comp(opara(5), 'init', 4, lopara(5),4)) THEN
               IF(ALLOCATED(accum_hkl)) DEALLOCATE(accum_hkl)
               accum_n = 0
            ENDIF
            hkl(1:3) = special_hkl(:,1)
            hkl(4)   = 0
            CALL point_init(hkl, point_hkl, point_n)
            IF(l_form) THEN
               DO k = 1,special_n
                  hkl(1:3) = special_hkl(:,k)
                  hkl(4)   = 0
                  CALL point_set(hkl, point_hkl, point_n)
               ENDDO
            ENDIF
            IF(.NOT.ALLOCATED(accum_hkl)) THEN
               ALLOCATE  (accum_hkl(1:3, 1:accum_n + 48))
               accum_n = 0
            ELSE
            IF(accum_n + point_n > UBOUND(accum_hkl,1)) THEN
               IF(ALLOCATED(accum_hkl)) THEN
                  temp_hkl = accum_hkl
                  DEALLOCATE(accum_hkl)
                  ALLOCATE  (accum_hkl(1:3, 1:accum_n + 48))
                  accum_hkl(1:3,1:UBOUND(temp_hkl,2)) = temp_hkl(:,:)
                  DEALLOCATE(temp_hkl)
               ELSE
                  ALLOCATE  (accum_hkl(1:3, 1:accum_n + 48))
               ENDIF
            ENDIF
            ENDIF
            accum_hkl(1:3,accum_n+1:accum_n+point_n) = point_hkl(1:3,1:point_n)
            accum_n = accum_n + point_n
         ELSEIF (str_comp (cpara (1) , 'sphere', 3, lpara (1) , 6) ) THEN
!                                                                       
!     ----Crystal is limited by a sphere                                
!                                                                       
            IF (ianz.eq.2) then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               IF (ier_num.ne.0) return 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               IF (ier_num.ne.0) return 
               radius = werte (1) 
               IF (radius.le.0.0) then 
                  ier_num = - 28 
                  ier_typ = ER_APPL 
                  RETURN 
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
            l_plane  = .false. 
            l_form   = .false. 
            l_sphere = .true. 
            l_cyl    = .false. 
            l_ell    = .false. 
         ELSEIF(str_comp(cpara(1), 'cylinder', 3, lpara(1), 8)) THEN
!                                                                       
!     ----Crystal is limited by a cylinder                              
!                                                                       
            IF (ianz.eq.3) then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               IF (ier_num.ne.0) return 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               IF (ier_num.ne.0) return 
               radius = werte (1) 
               height = werte (2) 
               IF (radius.le.0.0.or.height.le.0.0) then 
                  ier_num = - 28 
                  ier_typ = ER_APPL 
                  RETURN 
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
            l_plane  = .false. 
            l_form   = .false. 
            l_sphere = .false. 
            l_cyl    = .true. 
            l_ell    = .false. 
         ELSEIF(str_comp(cpara(1), 'ellipsoid', 3, lpara(1), 9)) THEN
!                                                                       
!     ----Crystal is limited by a standardized ellipsoid
!                                                                       
            IF (ianz.eq.4) then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               IF (ier_num.ne.0) return 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               IF (ier_num.ne.0) return 
               radius_ell(1:3) = werte(1:3)*0.5  ! User provides diameters
            ELSE 
               ier_num = - 6 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
            l_plane  = .false. 
            l_form   = .false. 
            l_sphere = .false. 
            l_cyl    = .false. 
            l_ell    = .true. 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
            RETURN
         ENDIF 
!
!
         IF (ier_num /= 0) RETURN
            IF ((l_plane .OR. l_form) .AND. str_comp (opara(6) , 'run', 3, lopara(6) , 3)) THEN
form_loop:     DO i = 1, cr_natoms 
                  IF(BTEST(cr_prop(i), PROP_OUTSIDE)) cycle form_loop 
                  IF(linside) THEN
                     DO j=1,accum_n
                         d = 1.0 - (cr_pos (1, i)-center(1)) * accum_hkl (1,j) &
                                 - (cr_pos (2, i)-center(2)) * accum_hkl (2,j) &
                                 - (cr_pos (3, i)-center(3)) * accum_hkl (3,j) 
                         d = d / dstar 
                         h(:) = accum_hkl (:,j)
                         CALL boundarize_atom (center, d, i, linside, SURF_PLANE, h) 
                     ENDDO 
                  ELSE
                     dshort = 1.E8
                     nplanes = 0
                     iplane  = 0 
                     DO j=1,accum_n
                         d = 1.0 - (cr_pos (1, i)-center(1)) * accum_hkl (1,j) &
                                 - (cr_pos (2, i)-center(2)) * accum_hkl (2,j) &
                                 - (cr_pos (3, i)-center(3)) * accum_hkl (3,j)
                         d = d / dstar 
                         IF(ABS(d)<surf_ex_dist(cr_iscat(i) ) ) THEN 
                            nplanes = nplanes + 1
                            iplane  = j
                         ENDIF
                         dshort = MIN(dshort, d)
                     ENDDO 
                     IF(nplanes>2) THEN            ! Atom is at a corner
                        h(:) = NINT(100*cr_pos(:,i))
                        CALL boundarize_atom (center, dshort, i, linside, SURF_CORNER, h) 
                     ELSEIF(nplanes==2) THEN            ! Atom is at an edge 
                        h(:) = NINT(100*cr_pos(:,i))
                        CALL boundarize_atom (center, dshort, i, linside, SURF_EDGE  , h) 
                     ELSEIF(nplanes==1) THEN            ! Atom is at a PLANE 
                        h(:) = accum_hkl (:,iplane)        !WRONG NEEDS WORK
                        CALL boundarize_atom (center, dshort, i, linside, SURF_PLANE, h) 
                     ELSE
                        h(:) = NINT(100*cr_pos(:,i))
                        CALL boundarize_atom (center, dshort, i, linside, SURF_NONE , h) 
                     ENDIF
                  ENDIF
               ENDDO form_loop
               DEALLOCATE  (accum_hkl)      ! reset the accumulation list
               accum_n = 0                  ! reset the accumulation list
            ELSEIF (l_sphere) then 
sphere_loop:   DO i = 1, cr_natoms 
                  IF(BTEST(cr_prop(i), PROP_OUTSIDE)) cycle sphere_loop 
                  v (1) = cr_pos (1, i)-center(1)
                  v (2) = cr_pos (2, i)-center(2)
                  v (3) = cr_pos (3, i)-center(3)
                  d = radius - sqrt (v(1) * v(1) * cr_gten(1, 1)        &
                     +     v(2) * v(2) * cr_gten(2, 2)    &
                     +     v(3) * v(3) * cr_gten(3, 3)    &
                     + 2 * v(1) * v(2) * cr_gten(1, 2)    &
                     + 2 * v(1) * v(3) * cr_gten(1, 3)    &
                     + 2 * v(2) * v(3) * cr_gten(2, 3)    )                                       
                  h(:) = v(:)
                  CALL boundarize_atom (center, d, i, linside, SURF_SPHERE, h)
               ENDDO sphere_loop
            ELSEIF (l_cyl) THEN
cyl_loop:      DO i = 1, cr_natoms 
                  IF(BTEST(cr_prop(i), PROP_OUTSIDE)) cycle cyl_loop 
                  IF(linside) THEN
                  v (1) = cr_pos (1, i)-center(1)
                  v (2) = cr_pos (2, i)-center(2)
                  v (3) = 0.0          !-center(3)
                  d = radius - sqrt(v(1) * v(1) * cr_gten(1, 1)        &
                     +     v(2) * v(2) * cr_gten(2, 2)    &
                     +     v(3) * v(3) * cr_gten(3, 3)    &
                     + 2 * v(1) * v(2) * cr_gten(1, 2)    &
                     + 2 * v(1) * v(3) * cr_gten(1, 3)    &
                     + 2 * v(2) * v(3) * cr_gten(2, 3)    )                                       
                  h(:) = v(:)
                  CALL boundarize_atom (center, d, i, linside, SURF_CYLINDER, h) 
                  IF(BTEST(cr_prop(i), PROP_OUTSIDE)) cycle cyl_loop 
                  v (1) = 0.0          !-center(1)
                  v (2) = 0.0          !-center(2)
                  v (3) = cr_pos (3, i)-center(3)
                  d = height - sqrt(v(1) * v(1) * cr_gten(1, 1)        &
                     +     v(2) * v(2) * cr_gten(2, 2)    & 
                     +     v(3) * v(3) * cr_gten(3, 3)    & 
                     + 2 * v(1) * v(2) * cr_gten(1, 2)    &
                     + 2 * v(1) * v(3) * cr_gten(1, 3)    & 
                     + 2 * v(2) * v(3) * cr_gten(2, 3)    )                                       
                  h(:) = v(:)
                  CALL boundarize_atom (center, d, i, linside, SURF_PLANE, h) 
                  ELSE                             ! 
                     dshort = 1.E8
                     lwall =.FALSE.
                     ltop  =.FALSE.
                     v (1) = cr_pos (1, i)-center(1)
                     v (2) = cr_pos (2, i)-center(2)
                     v (3) = 0.0          -center(3)
                     d = radius - sqrt(v(1) * v(1) * cr_gten(1, 1)        &
                        +     v(2) * v(2) * cr_gten(2, 2)    &
                        +     v(3) * v(3) * cr_gten(3, 3)    &
                        + 2 * v(1) * v(2) * cr_gten(1, 2)    &
                        + 2 * v(1) * v(3) * cr_gten(1, 3)    &
                        + 2 * v(2) * v(3) * cr_gten(2, 3)    )
                     wall(:) = v(:)
                     IF(ABS(d)<surf_ex_dist(cr_iscat(i) ) ) THEN 
                        lwall = .TRUE.               ! Atom is close to wall
                     ENDIF
                     dshort = MIN(dshort, d)
                     IF(BTEST(cr_prop(i), PROP_OUTSIDE)) cycle cyl_loop 
                     v (1) = 0.0          -center(1)
                     v (2) = 0.0          -center(2)
                     v (3) = cr_pos (3, i)-center(3)
                     d = height - sqrt(v(1) * v(1) * cr_gten(1, 1)        &
                        +     v(2) * v(2) * cr_gten(2, 2)    & 
                        +     v(3) * v(3) * cr_gten(3, 3)    & 
                        + 2 * v(1) * v(2) * cr_gten(1, 2)    &
                        + 2 * v(1) * v(3) * cr_gten(1, 3)    & 
                        + 2 * v(2) * v(3) * cr_gten(2, 3)    )
                     top(:) = v(:)
                     IF(ABS(d)<surf_ex_dist(cr_iscat(i) ) ) THEN 
                        ltop  = .TRUE.               ! Atom is close to wall
                     ENDIF
                     dshort = MIN(dshort, d)
                     IF(lwall .AND. .NOT.ltop) THEN   ! Atom is at wall only
                        CALL boundarize_atom (center, dshort, i, linside, SURF_CYLINDER, wall) 
                     ELSEIF(.NOT.lwall .AND. ltop) THEN   ! Atom is at top  only
                        CALL boundarize_atom (center, dshort, i, linside, SURF_PLANE, top) 
                     ELSEIF(lwall .AND. ltop) THEN   ! Atom is at edge
                        h(:) = cr_pos(:,i)
                        CALL boundarize_atom (center, dshort, i, linside, SURF_EDGE , h) 
                     ELSE
                        h(:) = cr_pos(:,i)
                        CALL boundarize_atom (center, dshort, i, linside, SURF_NONE , h) 
                     ENDIF
                  ENDIF
               ENDDO cyl_loop
            ELSEIF (l_ell) then 
               CALL plot_ini_trans (1.0,                          &
                 pl_tran_g, pl_tran_gi, pl_tran_f, pl_tran_fi, &
                 cr_gten, cr_rten, cr_eps)
               radius = (radius_ell(1)*radius_ell(2)*radius_ell(3))**(1./3.)
ell_loop:      DO i = 1, cr_natoms 
                  IF(BTEST(cr_prop(i), PROP_OUTSIDE)) cycle ell_loop 
                  v(:) = cr_pos(:, i)-center(:)
                  v = MATMUL(pl_tran_f(1:3,1:3), v)
                  d = (1. - sqrt((v(1)/radius_ell(1))**2   &
                                +(v(2)/radius_ell(2))**2   &
                                +(v(3)/radius_ell(3))**2   ))* radius
                  h(:) = v(:)
                  CALL boundarize_atom (center, d, i, linside, SURF_SPHERE, h)
               ENDDO ell_loop
            ENDIF 
!        ENDIF 
!     ENDIF 
!                                                                       
      END SUBROUTINE boundary                       
!*****7*****************************************************************
      SUBROUTINE boundarize_atom (center, distance, iatom, linside, &
                                  surface_type, normal) 
!-                                                                      
!     This subroutine sets the boundary property of the atom            
!                                                                       
!     Version  : 1.0                                                    
!     Date     : 07 Sep 10                                              
!                                                                       
!     Author   : R.B. Neder (reinhard.neder@fau.de)      
!+                                                                      
      USE discus_config_mod 
      USE discus_allocate_appl_mod
      USE crystal_mod 
      USE metric_mod
      USE prop_para_mod 
      USE surface_mod 
      USE errlist_mod 
      USE math_sup
      IMPLICIT none 
!                                                                       
       
!                                                                       
      REAL   , DIMENSION(3), INTENT(IN) :: center 
      REAL                 , INTENT(IN) :: distance 
      INTEGER              , INTENT(IN) :: iatom 
      LOGICAL              , INTENT(IN) :: linside 
      INTEGER              , INTENT(IN) :: surface_type
      REAL   , DIMENSION(3), INTENT(IN) :: normal
!
      REAL, DIMENSION(3), PARAMETER :: VNULL = (/ 0.0, 0.0, 0.0 /)
      REAL              , PARAMETER :: TOLERANCE = 5.0   ! Accept a 5 degree tilt for same surface
      LOGICAL           , PARAMETER :: LSPACE = .FALSE.
!
      INTEGER               :: idiv      ! Largest common divisor for the normal
      INTEGER, DIMENSION(3) :: hkl
      REAL   , DIMENSION(3) :: rhkl, u
      REAL                  :: angle
!                                                                       
      IF( cr_nscat > SURF_MAXSCAT) THEN
         CALL alloc_surf ( cr_nscat )
         IF ( ier_num < 0 ) THEN
            RETURN
         ENDIF
      ENDIF
!
      IF((     linside.AND.distance <  0) .OR.  &
         (.not.linside.AND.distance >  0)      ) THEN                            
         cr_iscat (iatom) = 0 
         cr_prop (iatom) = ibclr (cr_prop (iatom), PROP_NORMAL) 
         cr_prop (iatom) = ibset (cr_prop (iatom), PROP_OUTSIDE) 
         IF (abs (distance) .lt.surf_ex_dist (cr_iscat (iatom) ) ) then 
            cr_prop (iatom) = ibset (cr_prop (iatom), PROP_SURFACE_EXT) 
         ENDIF 
         cr_surf (:,iatom) = 0
      ELSE 
         IF (abs (distance) .lt.surf_ex_dist (cr_iscat (iatom) ) ) then 
            cr_prop (iatom) = ibset (cr_prop (iatom), PROP_SURFACE_EXT) 
            IF(cr_surf(0, iatom) == SURF_NONE) THEN  ! Atom was not yet at a surface
               cr_surf(0,   iatom) = surface_type 
               cr_surf(1:3, iatom) = NINT(10.0*normal(:))
            ELSEIF(cr_surf(0, iatom) < SURF_EDGE) THEN    ! Atom was already at a plane, sphere, cylinder wall
               IF(surface_type > SURF_CYLINDER) THEN      ! New location is at least an edge
                  u(:) = cr_pos(:,iatom) - center(:)
                  CALL pos2hkl(u, cr_gten, hkl)
                  cr_surf(0,   iatom) = SURF_CORNER
                  cr_surf(1:3, iatom)  = hkl(1:3)
               ELSEIF(surface_type > SURF_NONE) THEN      ! Should be all other cases
                  rhkl(:) = cr_surf(1:3, iatom)           ! get old normal
                  angle = do_bang(LSPACE, rhkl, VNULL, normal)
                  IF(angle>TOLERANCE) THEN
                      u(:) = cr_pos(:,iatom) - center(:)
                      CALL pos2hkl(u, cr_gten, hkl)
                      cr_surf(0,   iatom) = SURF_EDGE
                      cr_surf(1:3, iatom)  = hkl(1:3)
                  ENDIF
               ENDIF
            ELSE                                          ! Atom was at edge or corner
               u(:) = cr_pos(:,iatom) - center(:)
               CALL pos2hkl(u, cr_gten, hkl)
               cr_surf(0,   iatom) = SURF_CORNER
               cr_surf(1:3, iatom)  = hkl(1:3)
            ENDIF
            IF(.NOT.linside) cr_surf(:, iatom) = -cr_surf(:, iatom)  !invert for inside pointing surface
         ENDIF 
      ENDIF 
      IF(cr_surf(0, iatom) /= 0) THEN
         idiv = gcd(cr_surf(1, iatom), cr_surf(2, iatom), cr_surf(3, iatom))
         IF(idiv /= 0) THEN
             cr_surf(1:3, iatom) = cr_surf(1:3, iatom)/IABS(idiv)
         ENDIF
      ENDIF
!                                                                       
      END SUBROUTINE boundarize_atom                
!
!*****7*****************************************************************
!
      SUBROUTINE pos2hkl(u, gten, hkl)
!
      USE trafo_mod
!
      IMPLICIT NONE
!
      REAL   , DIMENSION(3)  , INTENT(INOUT) :: u
      REAL   , DIMENSION(3,3), INTENT(IN )   :: gten
      INTEGER, DIMENSION(3)  , INTENT(OUT)   :: hkl
!
      REAL   , DIMENSION(3) :: v
      REAL :: uu
!
      uu = SQRT(u(1)**2 + u(2)**2 + u(3)**2)  ! Really rough length
      u(:) = u(:)/uu
      CALL trans (u, gten, v, 3)              ! Transform into reciprocal space
      hkl(1:3) = NINT(10.*v(1:3))             ! Make an approximate integer vector
!
      END SUBROUTINE pos2hkl
!
!*****7*****************************************************************
!
SUBROUTINE do_surface_char(zeile, lp)
!
USE crystal_mod 
USE prop_para_mod
USE param_mod
USE prompt_mod
USE ber_params_mod
USE errlist_mod
USE get_params_mod
IMPLICIT NONE
CHARACTER(LEN=*) , INTENT(INOUT) :: zeile
INTEGER          , INTENT(INOUT) :: lp
!
INTEGER, PARAMETER :: MAXW = 3
CHARACTER(LEN=1024), DIMENSION(1:MAXW) :: cpara
INTEGER            , DIMENSION(1:MAXW) :: lpara
REAL               , DIMENSION(1:MAXW) :: werte
INTEGER               :: surf_char
LOGICAL               :: lshow, lequal
INTEGER                 :: i
INTEGER, DIMENSION(3,6) :: surf_normal
INTEGER, DIMENSION(3)   :: surf_kante
INTEGER, DIMENSION(6)   :: surf_weight
!
INTEGER :: ianz, iatom
!
LOGICAL str_comp
!
cpara(:) = ' '
lpara(:) = 1
CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
IF(ier_num/=0) RETURN
!
! Check if output is desired
!
lshow = .FALSE.
IF (str_comp (cpara (ianz) , 'show', 2, lpara (1) , 4) ) then 
  lshow = .TRUE.
  cpara(ianz) = ' '
  lpara(ianz) = 1
  ianz = ianz - 1
ENDIF
!
! Check if Atoms are restricted to equal or any atom tpye
!
lequal = .TRUE.
IF (str_comp (cpara (2) , 'equal', 2, lpara (1) , 5) ) then 
  lequal = .TRUE.
  cpara(2) = ' '
  lpara(2) = 1
  ianz = ianz - 1
ELSEIF (str_comp (cpara (2) , 'any', 2, lpara (1) , 3) ) then 
  lequal = .FALSE.
  cpara(2) = ' '
  lpara(2) = 1
  ianz = ianz - 1
ENDIF
!
! calculate atom number
!
CALL ber_params (ianz, cpara, lpara, werte, maxw) 
IF(ier_num/=0) RETURN
!
iatom = NINT(werte(1))
CALL surface_character(iatom, 1, cr_natoms, surf_char, surf_normal, surf_kante, surf_weight, lequal)
!
i = 0
res_para(0) = 1
res_para(1) = surf_char
res_para(2:4) = surf_normal(1:3,1)
!IF(surf_char==1) THEN                 ! Planar surface
!   res_para(0) = 4
!   res_para(2:4) = surf_normal(1:3,1)
!ELSEIF(surf_char==2) THEN             ! An edge
!   res_para(0)   = 10
!   res_para(2:4) = surf_normal(1:3,1)
!   res_para(5:7) = surf_normal(1:3,2)
!   res_para(8:10)= surf_kante(1:3)
!ELSEIF(surf_char==3) THEN             ! A corner
!   i = 0
!   fill: DO
!      IF(surf_weight(i+1) == 0) EXIT fill
!      i = i+1
!      res_para((i-1)*3+2:(i-1)*3+4) = surf_normal(1:3,i)
!   ENDDO fill
!   res_para(0) = 3*i+1
!ENDIF
IF(lshow) THEN
   IF(surf_char==0) THEN
      IF(IBITS(cr_prop(iatom),PROP_SURFACE_EXT,1).eq.1 .and.        &  ! real Atom is near surface
         IBITS(cr_prop(iatom),PROP_OUTSIDE    ,1).eq.0       ) THEN    ! real Atom is near surface
         WRITE(output_io, 2000) 
      ELSE
         WRITE(output_io, 2050)
      ENDIF
   ELSEIF(surf_char==1) THEN
      WRITE(output_io, 2100) 'Planar      outside face ',surf_normal(:,1)
   ELSEIF(surf_char==2) THEN
      WRITE(output_io, 2100) 'Spherical   outside face ',surf_normal(:,1)
   ELSEIF(surf_char==3) THEN
      WRITE(output_io, 2100) 'Cylindrical outside face ',surf_normal(:,1)
   ELSEIF(surf_char==4) THEN
      WRITE(output_io, 2100) 'Outside edge,  normal is ',surf_normal(:,1)
   ELSEIF(surf_char==5) THEN
      WRITE(output_io, 2100) 'Outside corner,normal is ',surf_normal(:,1)
   ELSEIF(surf_char==-1) THEN
      WRITE(output_io, 2100) 'Planar       inside face ',surf_normal(:,1)
   ELSEIF(surf_char==-2) THEN
      WRITE(output_io, 2100) 'Spherical    inside face ',surf_normal(:,1)
   ELSEIF(surf_char==-3) THEN
      WRITE(output_io, 2100) 'Cylindrical  inside face ',surf_normal(:,1)
   ELSEIF(surf_char==-4) THEN
      WRITE(output_io, 2100) ' Inside edge,  normal is ',surf_normal(:,1)
   ELSEIF(surf_char==-5) THEN
      WRITE(output_io, 2100) ' Inside corner,normal is ',surf_normal(:,1)
   ENDIF
ENDIF
!
2000 FORMAT(' Surface character could not be determined')
2050 FORMAT(' Atom is not close to a surface')
2100 FORMAT(a,3(I4,2x))
!2200 FORMAT(' Edge   surface,  Edge, dominant normal: ',3(I4,2x), 2x, 3(I4,2x), ':', I4)
!2250 FORMAT('                       secondary normal: ',20x         , 3(I4,2x), ':', I4)
!2300 FORMAT(' Corner surface,  dominant normal: ',3(I4,2x), ':', I4)
!2350 FORMAT('                 secondary normal: ',3(I4,2x), ':', I4)
!3100 FORMAT(' Indented plane,  normal: ',3(I4,2x))
!3200 FORMAT(' Indented edge ,  Edge, dominant normal: ',3(I4,2x), 2x, 3(I4,2x), ':', I4)
!3300 FORMAT(' Indented corner, dominant normal: ',3(I4,2x), ':', I4)
END SUBROUTINE do_surface_char
!
!*****7*****************************************************************
!
SUBROUTINE surface_character(iatom, jstart, jend, surf_char, surf_normal,  &
                             surf_kante, surf_weight, equal)
!
use atom_name
USE atom_env_mod
USE chem_mod
USE crystal_mod
USE metric_mod
USE prop_para_mod
!
USE param_mod
USE prompt_mod
USE math_sup
!
IMPLICIT NONE
INTEGER,                 INTENT(IN)  :: iatom
INTEGER,                 INTENT(IN)  :: jstart
INTEGER,                 INTENT(IN)  :: jend
INTEGER,                 INTENT(OUT) :: surf_char
INTEGER, DIMENSION(3,6), INTENT(OUT) :: surf_normal
INTEGER, DIMENSION(3)  , INTENT(OUT) :: surf_kante
INTEGER, DIMENSION(6)  , INTENT(OUT) :: surf_weight
LOGICAL,                 INTENT(IN)  :: equal
!
!INTEGER, PARAMETER :: SURF_IN_PLANE  = -1
INTEGER, PARAMETER :: SURF_NONE   =  0
INTEGER, PARAMETER :: SURF_PLANE  =  1
INTEGER, PARAMETER :: SURF_EDGE   =  2
INTEGER, PARAMETER :: SURF_CORNER =  3
!
INTEGER, PARAMETER                     :: MAXW = 1
!LOGICAL, PARAMETER :: LNEW = .false.
REAL   , PARAMETER                     :: RADIUS_MIN  = 1.51
REAL   , PARAMETER                     :: RADIUS_STEP = 1.50
REAL   , PARAMETER , DIMENSION(3)      :: NULL        = 0.0
REAL   , PARAMETER                     :: IS_OUTSIDE  = 80.0
REAL   , PARAMETER                     :: IS_PARALLEL = 15.0
!
CHARACTER(LEN=1024), DIMENSION(1:MAXW) :: cpara
INTEGER            , DIMENSION(1:MAXW) :: lpara
REAL               , DIMENSION(1:MAXW) :: werte
!
CHARACTER(LEN=1024)     :: line
INTEGER  , DIMENSION(:), ALLOCATABLE :: neigh
REAL     , DIMENSION(:), ALLOCATABLE :: angles
REAL     , DIMENSION(:), ALLOCATABLE :: sorted
INTEGER  , DIMENSION(:,:), ALLOCATABLE :: surfaces
INTEGER                 :: nsurface
INTEGER                 :: i,j,k,l, ianz
INTEGER                 :: counter
INTEGER                 :: laenge
INTEGER                 :: neigsurf
INTEGER                 :: divisor
LOGICAL  , DIMENSION(3) :: fp
LOGICAL                 :: fq
LOGICAL                 :: lspace
LOGICAL                 :: isfound
LOGICAL                 :: isfirst
REAL                    :: rmin, radius
REAL                    :: alpha, beta
REAL                    :: dstar
REAL                    :: dist, dmin
REAL     , DIMENSION(3) :: x         ! Vector from center to atom
REAL     , DIMENSION(3) :: center    ! Average position of neighboring atoms
INTEGER  , DIMENSION(3) :: rough     ! rough normal 
REAL     , DIMENSION(3) :: u,v,w     ! Vectors from central atom to neighbors
INTEGER  , DIMENSION(3) :: tempsurf  ! Vectors from central atom to neighbors
REAL     , DIMENSION(3) :: realsurf  ! Vectors from central atom to neighbors
!
INTEGER, DIMENSION(0:1) :: temp_sel_prop = 0
!
temp_sel_prop(:) = cr_sel_prop(:)   ! save user settings for property select
cr_sel_prop(:)   = 0                ! Ignore all property selection rules 
surf_char = SURF_NONE
surf_normal = 0
surf_kante  = 0
surf_weight = 0
fp (1) = chem_period (1)
fp (2) = chem_period (2)
fp (3) = chem_period (3)
fq = chem_quick
rmin = 0.5
nsurface = 0
isfirst = .TRUE.
!
IF(IBITS(cr_prop(iatom),PROP_SURFACE_EXT,1).eq.1 .and.        &  ! real Atom is near surface
   IBITS(cr_prop(iatom),PROP_OUTSIDE    ,1).eq.0       ) THEN    ! real Atom is near surface
   IF(cr_surf(0,iatom)/=0) THEN     ! Atom has surface character set
      surf_char = cr_surf(0, iatom)
      surf_normal(:, 1) = cr_surf(1:3, iatom)
      surf_kante(:)     = 0
      surf_weight(:)    = 0
      surf_weight(1)    = 1
      RETURN                        ! Success go back
   ELSE
      RETURN
   ENDIF
ELSE
   RETURN
ENDIF
!OLD CODE TO BE WORKED ON
!
!
IF(IBITS(cr_prop(iatom),PROP_SURFACE_EXT,1).eq.1 .and.        &  ! real Atom is near surface
   IBITS(cr_prop(iatom),PROP_OUTSIDE    ,1).eq.0       ) THEN    ! real Atom is near surface
   cpara(1) = 'all'
   lpara(1) = 3
   ianz     = 1
   x(1)     = cr_pos(1,iatom)              ! Vector from center to atom
   x(2)     = cr_pos(2,iatom)
   x(3)     = cr_pos(3,iatom)
   WRITE(line,'(2(G16.8E3,a1),G16.8E3)') x(1), ',', x(2), ',', x(3)
   lspace = .TRUE.
   laenge = 50
   CALL d2r(line, laenge, lspace)
   rough(1:3) = INT(res_para(1:3))      ! Rough normal 
!
   radius   = RADIUS_MIN
   counter = 0                          ! Prevent infinit loop
   grand: DO                            ! increase search radius if character is not yet found
      counter  = counter + 1
      radius   = radius + RADIUS_STEP
      ianz     = 1
      werte(1) = -1                     ! Find all atom types
!
!     Determine a local center of mass with a larger search radius
!
      center(:) = 0.0
      CALL do_find_env (ianz, werte, maxw, x, rmin, radius+RADIUS_STEP, fq, fp)
      DO i=1, atom_env(0)               ! Pick out surface atom types only
         IF(jstart<=atom_env(i) .AND. atom_env(i)<=jend)  &
               center(:) = center(:) + cr_pos(:,atom_env(i))
      ENDDO
      center(:) = center(:) / atom_env(0)
!
!     Find local environment to construct the surface normal
!
      CALL do_find_env (ianz, werte, maxw, x, rmin, radius, fq, fp)
      ALLOCATE(neigh(0:atom_env(0)))
      neigsurf = 0
      DO i=1, atom_env(0)               ! Pick out surface atom types only
         IF(IBITS(cr_prop(atom_env(i)),PROP_SURFACE_EXT,1).eq.1 .and.        &  ! real Atom is near surface
            IBITS(cr_prop(atom_env(i)),PROP_OUTSIDE    ,1).eq.0       ) THEN    ! real Atom is near surface
            IF(.NOT.equal .OR. (equal .AND. cr_iscat(atom_env(i))==cr_iscat(iatom)) ) THEN
               neigsurf = neigsurf + 1
               neigh(neigsurf) = atom_env(i)
            ENDIF
         ENDIF
      ENDDO
!
      IF(neigsurf >= 3 ) THEN             ! Found three surface atoms as neighbors
!        Determine rough surface as vector from center of mass to actual atom
         WRITE(line,'(2(G16.8E3,a1),G16.8E3)') x(1)-center(1), ',', x(2)-center(2), ',', x(3)-center(3)
         lspace = .TRUE.
         laenge = 50
         CALL d2r(line, laenge, lspace)
         rough(1:3) = INT(10*res_para(1:3))      ! Rough local normal 
!                                         ! Make space for all possible surfaces
         ALLOCATE(surfaces(0:3,neigsurf*neigsurf/2))
         nsurface      = 0
         surfaces(:,:) = 0
         indented: DO i=1, neigsurf-1               ! Loop over all neigbor pairs
            u(:) = cr_pos(:,neigh(i))-x(:)
            inner: DO j=i+1, neigsurf
               v(:) = cr_pos(:,neigh(j))-x(:)
               WRITE(line,2000) u,v
2000 FORMAT(6(G16.8E3,','),'ddr')
               laenge = 105
               CALL vprod(line, laenge)
               tempsurf(:) = NINT(res_para(1:3))  ! Normal to atom triplet
                                                  ! If necessary invert direction
               lspace = .FALSE.
               IF(do_blen(lspace, NULL, FLOAT(tempsurf))>0) THEN
               lspace = .FALSE.
               IF(do_bang(lspace, FLOAT(rough), NULL, FLOAT(tempsurf)) > 90.0) tempsurf(:) = -tempsurf(:)
               WRITE(line,2010) tempsurf(:)
               laenge = 50
               lspace = .FALSE.
               CALL d2r(line, laenge, lspace)     ! Calculate direction in real space
2010 FORMAT(2(G16.8E3,','),G16.8E3) 
               realsurf(:) = res_para(1:3)
!              Test to see if any neighbor is further away than central atom
!
               DO k=1, neigsurf
                  v(:) = cr_pos(:,neigh(k))-x(:)
                  lspace = .TRUE.
                  IF(do_bang(lspace, realsurf, NULL, v)<IS_OUTSIDE) THEN
                     CYCLE inner ! proceed to next pair
                  ENDIF
               ENDDO
!              Not an indented surface, proceed to sort
               isfound = .FALSE.
               lspace = .FALSE.
               DO k=1, nsurface                  ! Loop over all previous surfaces
                  alpha = do_bang(lspace, FLOAT(tempsurf), NULL, FLOAT(surfaces(1:3,k)) )
                  IF(alpha>90.) alpha = 180.-alpha
!
                  IF(alpha<IS_PARALLEL) THEN
                     surfaces(0,k) = surfaces(0,k) + 1   ! increment count of triplets with identical normal
                     isfound = .TRUE.
                  ENDIF
               ENDDO
!
               IF(.NOT.isfound) THEN            ! New surface normal
                  nsurface = nsurface + 1
                  surfaces(0,nsurface) = 1
                  surfaces(1:3,nsurface) = tempsurf
               ENDIF
               ENDIF
            ENDDO inner
         ENDDO indented
!        If only one surface, test if neighbors are on all sides
         IF(nsurface==1) THEN
            ALLOCATE(angles(1:neigsurf+1))
            ALLOCATE(sorted(1:neigsurf+1))
            angles(:) = 0.0
            sorted(:) = 0.0
            WRITE(line,2020) cr_pos(:,neigh(1))-x(:), surfaces(1:3,1)
            laenge = 106
2020 FORMAT(6(G16.8E3,','),'drdd') 
            CALL do_proj(line, laenge)
            u(:) = res_para(4:6)
            angles(1) = 0
            angles(neigsurf+1) = 360.0
            lspace = .TRUE.
            DO k=2,neigsurf
               WRITE(line,2020) cr_pos(:,neigh(k))-x(:), surfaces(1:3,1)
               laenge = 106
               CALL do_proj(line, laenge)
               v(:) = res_para(4:6)
               alpha = do_bang(lspace, u, NULL, v)
               WRITE(line,2000) u,v
               laenge = 105
               CALL vprod(line, laenge)
               w(:) = res_para(1:3)
               IF(do_blen(lspace, NULL, w) > 0.0001) THEN
                  beta = do_bang(lspace, w, NULL, realsurf)
                  IF(beta.gt.90) alpha = 360. - alpha
               ENDIF
               angles(k) = alpha
            ENDDO
            sorted(1) =   0.0
            angles(1) = 400.
            sorted(neigsurf+1) = 360.0
            angles(neigsurf+1) = 400.
            DO k=2,neigsurf
               i=MINLOC(angles,1)
               sorted(k) = angles(i)
               angles(i) = 400.0
            ENDDO
            alpha = 0.0
            DO k=2, neigsurf+1
               beta = sorted(k) - sorted(k-1)
               IF(beta > alpha) alpha = beta
            ENDDO
            DEALLOCATE(angles)
            DEALLOCATE(sorted)
            IF(alpha > 190) THEN       ! Current atom is at the corner of plane!
               nsurface = 3
            ELSEIF(alpha > 170) THEN   ! Current atom is at the edge of plane!
               nsurface = 2
            ENDIF
         ENDIF
         IF(nsurface>2 .AND. isfirst) THEN
            isfirst = .FALSE.
            DEALLOCATE(neigh)
            DEALLOCATE(surfaces)
            CYCLE grand
         ENDIF
         IF(nsurface>0) EXIT grand
         DEALLOCATE(neigh)
         DEALLOCATE(surfaces)
      ELSE                             ! Found less than three neigbors
!         radius = radius + RADIUS_STEP ! Increment search radius
         DEALLOCATE(neigh)
         IF(ALLOCATED(surfaces)) DEALLOCATE(surfaces)
      ENDIF
      IF(counter == 3) THEN 
         EXIT grand   ! Too many trials, give up
      ENDIF
   ENDDO grand
ENDIF
! Normalize
lspace = .false.
DO i=1, MIN(nsurface, UBOUND(surf_weight,1))
   dstar=do_blen(lspace, NULL, FLOAT(surfaces(1:3,i)))
   IF(dstar > 0) THEN
      surfaces(1:3,i) = NINT(FLOAT(surfaces(1:3,i))*10./dstar)
      divisor = IABS(gcd(surfaces(1,i),surfaces(2,i),surfaces(3,i)))
      surfaces(1:3,i) = surfaces(1:3,i)/divisor
   ENDIF
   surf_weight(i) = surfaces(0,i)
ENDDO
!
! might be used to differentiate external from internal surfaces ???
!
!w(1:3)=surfaces(1:3,1)
!beta = do_bang(lspace, w, NULL, x)
!write(*,*) ' normal, coord, angle ', w,x, beta
!
ianz = -1
CALL do_find_env (ianz, werte, maxw, x, rmin, RADIUS+RADIUS_STEP, fq, fp)
dmin = 1
DO i=1, nsurface
   DO k=1, atom_env(0)
   
      dist  = 1.0 - (cr_pos (1, atom_env(k))-x(1)) * surfaces (1,i) &
                  - (cr_pos (2, atom_env(k))-x(2)) * surfaces (2,i)  &
                  - (cr_pos (3, atom_env(k))-x(3)) * surfaces (3,i)
      dmin = MIN(dmin, dist)
   ENDDO
ENDDO
IF(nsurface == 0) THEN
   surf_char = SURF_NONE
ELSEIF(nsurface == 1) THEN
   surf_char        = SURF_PLANE
   surf_normal(:,1) = surfaces(1:3,1)
   surf_kante(:)    = 0
ELSEIF(nsurface == 2) THEN
   surf_char = SURF_EDGE
   i = MAXLOC(surf_weight,1)
   surf_normal(:,1) = surfaces(1:3,i)
   i = 3-i                          !Choose the other normal
   surf_normal(:,2) = surfaces(1:3,i)
   WRITE(line,2000) surfaces(1:3,1),surfaces(1:3,2)
   laenge = 105
   CALL vprod(line, laenge)
   surf_kante(:) = NINT(res_para(1:3)*100)
   dstar=do_blen(lspace, NULL, FLOAT(surf_kante(1:3)))
   surf_kante(1:3) = NINT(FLOAT(surf_kante(1:3))*10./dstar)
ELSEIF(nsurface >  2) THEN
   surf_char = SURF_CORNER
   fill:DO l = 1,6
      i = 1
      k=surfaces(0,1)
      DO j=2,nsurface
         IF(surfaces(0,j) > k ) THEN
            i = j
            k = surfaces(0,j)
         ENDIF
      ENDDO
      IF(k==0) EXIT fill
      surf_normal(:,l) = surfaces(1:3,i)
      surfaces(0,i) = -1
   ENDDO fill
   surf_kante(:)  = 0
ELSE
   surf_char = SURF_NONE
ENDIF
IF(dmin < 0) surf_char = -surf_char   ! Indented surface, return negative character
!
! Divide by greatest common divisor
!
k = 0
DO i=1,3
   k = MAX(k, IABS(surf_kante(i)))
ENDDO
IF(k>0) THEN
   divisor = gcd(surf_kante(1), surf_kante(2),surf_kante(3))
   surf_kante(:) = surf_kante(:)/divisor
ENDIF
IF(ALLOCATED(neigh)) DEALLOCATE(neigh)
IF(ALLOCATED(surfaces)) DEALLOCATE(surfaces)
cr_sel_prop(:) = temp_sel_prop(:)   ! restore user settings for property select
!
END SUBROUTINE surface_character
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
      IMPLICIT none 
       
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 7) 
!                                                                       
      CHARACTER ( * ) line 
      CHARACTER(1024) cpara (maxw) 
      INTEGER lpara (maxw) 
      INTEGER ianz 
      INTEGER laenge 
      REAL werte (maxw) 
!                                                                       
      LOGICAL str_comp 
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
      USE prompt_mod 
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
      REAL werte (maxw) 
!                                                                       
      LOGICAL str_comp 
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
!*****7*****************************************************************
      SUBROUTINE surface_menu 
!-                                                                      
!     Main menu for surface related operations                          
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE molecule_mod 
      USE prop_para_mod 
!
      USE calc_expr_mod
      USE doact_mod 
      USE do_eval_mod
      USE do_wait_mod
      USE errlist_mod 
      USE learn_mod 
      USE class_macro_internal 
      USE prompt_mod 
      USE sup_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxw 
!                                                                       
      CHARACTER(5) befehl 
      CHARACTER(LEN=LEN(prompt)) :: orig_prompt
      CHARACTER(1024) line, zeile
      INTEGER lp, length, lbef 
      INTEGER indxg
      LOGICAL lend
!                                                                       
      INTEGER len_str 
      LOGICAL str_comp 
!                                                                       
      maxw = MAXSCAT
      lend = .false. 
      CALL no_error 
      orig_prompt = prompt
      prompt = prompt (1:len_str (prompt) ) //'/surf' 
!                                                                       
      DO while (.not.lend) 
      CALL get_cmd (line, length, befehl, lbef, zeile, lp, prompt) 
      IF (ier_num.eq.0) then 
         IF (line /= ' '      .and. line(1:1) /= '#' .and. &
             line /= char(13) .and. line(1:1) /= '!'        ) THEN
!                                                                       
!     ----search for "="                                                
!                                                                       
indxg = index (line, '=') 
IF (indxg.ne.0.AND..NOT. (str_comp (befehl, 'echo', 2, lbef, 4) ) &
              .AND..NOT. (str_comp (befehl, 'syst', 2, lbef, 4) )    &
              .AND..NOT. (str_comp (befehl, 'help', 2, lbef, 4) .OR. &
                          str_comp (befehl, '?   ', 2, lbef, 4) )    &
              .AND. INDEX(line,'==') == 0                            ) THEN
!                                                                       
!     ------evaluatean expression and assign the value to a variabble   
!                                                                       
               CALL do_math (line, indxg, length) 
            ELSE 
!                                                                       
!------ ----execute a macro file                                        
!                                                                       
               IF (befehl (1:1) .eq.'@') then 
                  IF (length.ge.2) then 
                     CALL file_kdo (line (2:length), length - 1) 
                  ELSE 
                     ier_num = - 13 
                     ier_typ = ER_MAC 
                  ENDIF 
!                                                                       
!     ----continues a macro 'continue'                                  
!                                                                       
               ELSEIF (str_comp (befehl, 'continue', 2, lbef, 8) ) then 
                  CALL macro_continue (zeile, lp) 
!                                                                       
!------ ----Echo a string, just for interactive check in a macro 'echo' 
!                                                                       
               ELSEIF (str_comp (befehl, 'echo', 2, lbef, 4) ) then 
                  CALL echo (zeile, lp) 
!                                                                       
!      ---Evaluate an expression, just for interactive check 'eval'     
!                                                                       
               ELSEIF (str_comp (befehl, 'eval', 2, lbef, 4) ) then 
                  CALL do_eval (zeile, lp) 
!                                                                       
!     ----exit 'exit'                                                   
!                                                                       
               ELSEIF (str_comp (befehl, 'exit', 2, lbef, 4) ) then 
                  lend = .true. 
!                                                                       
!     ----help 'help','?'                                               
!                                                                       
      ELSEIF (str_comp (befehl, 'help', 2, lbef, 4) .or.str_comp (befehl&
     &, '?   ', 1, lbef, 4) ) then                                      
                  IF (str_comp (zeile, 'errors', 2, lp, 6) ) then 
                     lp = lp + 7 
                     CALL do_hel ('discus '//zeile, lp) 
                  ELSE 
                     lp = lp + 15 
                     CALL do_hel ('discus surface '//zeile, lp) 
                  ENDIF 
!                                                                       
!------- -Operating System Kommandos 'syst'                             
!                                                                       
               ELSEIF (str_comp (befehl, 'syst', 2, lbef, 4) ) then 
                  IF (zeile.ne.' ') then 
                     CALL do_operating (zeile (1:lp), lp) 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
!                                                                       
!------  -----waiting for user input                                    
!                                                                       
               ELSEIF (str_comp (befehl, 'wait', 3, lbef, 4) ) then 
                  CALL do_input (zeile, lp) 
!                                                                       
!     ----run a boundary 'boundary'                                     
!                                                                       
               ELSEIF (str_comp (befehl, 'boundary', 2, lbef, 8) ) then 
                  CALL boundary (zeile, lp) 
!                                                                       
!     ----Test the boundary character of an atom                        
!                                                                       
               ELSEIF (str_comp (befehl, 'character', 2, lbef, 9) ) then 
                  CALL do_surface_char(zeile, lp)
!                                                                       
!     ----define various surface related settings 'set'                 
!                                                                       
               ELSEIF (str_comp (befehl, 'set', 2, lbef, 3) ) then 
                  CALL surf_do_set (zeile, lp) 
!                                                                       
!     ----show current parameters 'show'                                
!                                                                       
               ELSEIF (str_comp (befehl, 'show', 2, lbef, 4) ) then 
                  CALL surf_show 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_COMM 
               ENDIF 
            ENDIF 
         ENDIF 
      ENDIF 
      IF (ier_num.ne.0) THEN 
         CALL errlist 
         IF (ier_sta.ne.ER_S_LIVE) THEN 
            IF (lmakro .OR. lmakro_error) THEN  ! Error within macro or termination errror
               IF(sprompt /= prompt ) THEN
                  ier_num = -10
                  ier_typ = ER_COMM
                  ier_msg(1) = ' Error occured in surface menu'
                  prompt_status = PROMPT_ON 
                  prompt = orig_prompt
                  RETURN
               ELSE
                  IF(lmacro_close) THEN
                     CALL macro_close 
                     prompt_status = PROMPT_ON 
                  ENDIF 
               ENDIF 
            ENDIF 
            IF (lblock) THEN 
               ier_num = - 11 
               ier_typ = ER_COMM 
               prompt_status = PROMPT_ON 
               prompt = orig_prompt
               RETURN 
            ENDIF 
            CALL no_error 
            lmakro_error = .FALSE.
            sprompt = ' '
         ENDIF 
      ENDIF 
      ENDDO 
!
      prompt = orig_prompt
!                                                                       
      END SUBROUTINE surface_menu                   
!*****7*****************************************************************
      SUBROUTINE surf_do_set (zeile, length) 
!+                                                                      
!     This subroutine sets various parameters                           
!-                                                                      
      USE discus_config_mod 
      USE errlist_mod 
      USE get_params_mod
      IMPLICIT none 
!                                                                       
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 20) 
!                                                                       
      CHARACTER ( * ) zeile 
      INTEGER length 
!                                                                       
      CHARACTER(1024) cpara (maxw) 
      INTEGER lpara (maxw) 
      INTEGER ianz 
      REAL werte (maxw) 
!                                                                       
      LOGICAL str_comp 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, length) 
      IF (ier_num.ne.0) return 
      IF (ianz.le.0) return 
!                                                                       
      IF (str_comp (cpara (1) , 'distance', 2, lpara (1) , 2) ) then 
         CALL del_params (1, ianz, cpara, lpara, maxw) 
         CALL surf_set_fuzzy (ianz, cpara, lpara, werte, maxw, 0) 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
      END SUBROUTINE surf_do_set                    
!*****7*****************************************************************
      SUBROUTINE surf_set_fuzzy (ianz, cpara, lpara, werte, maxw, iflag) 
!+                                                                      
!     This subroutine sets the distances between atoms and a surface.   
!-                                                                      
      USE discus_config_mod 
      USE discus_allocate_appl_mod
      USE crystal_mod 
      USE surface_mod 
      USE berechne_mod
      USE errlist_mod 
      USE get_params_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxw 
      INTEGER ianz 
      CHARACTER(1024) cpara (maxw) 
      INTEGER lpara (maxw) 
      REAL werte (maxw) 
      INTEGER iflag 
!                                                                       
      CHARACTER(1024) string 
      INTEGER laenge 
      INTEGER i 
      LOGICAL lold 
      LOGICAL :: linternal = .true.
      LOGICAL :: lexternal = .false.
      REAL distance 
!                                                                       
      LOGICAL str_comp 
!                                                                       
      lold = .false. 
!                                                                       
      IF (ianz.eq.1) then 
         ier_num = - 6 
         ier_typ = ER_COMM 
         RETURN 
      ENDIF 
!                                                                       
      IF( MAXSCAT > SURF_MAXSCAT) THEN
         CALL alloc_surf ( MAXSCAT )
         IF ( ier_num < 0 ) THEN
            RETURN
         ENDIF
      ENDIF
!
!     WRITE ( * , * ) 'IFLAG', iflag 
!     WRITE ( * , * ) 'IANZ ', ianz 
!     DO i = 1, ianz 
!     WRITE ( * , * ) '>', cpara (i) (1:lpara (i) ) , '<' 
!     ENDDO 
      IF (iflag.eq.SURF_SURFACE) THEN 
         IF (str_comp (cpara (1) , 'external', 2, lpara (1) , 8) ) then 
            lexternal = .true. 
            linternal = .false. 
         ELSEIF (str_comp (cpara (1) , 'internal', 2, lpara (1) , 8) )  &
         then                                                           
            lexternal = .false. 
            linternal = .true. 
         ELSEIF (str_comp (cpara (1) , 'all', 2, lpara (1) , 3) ) then 
            lexternal = .true. 
            linternal = .true. 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
            RETURN 
         ENDIF 
         CALL del_params (1, ianz, cpara, lpara, maxw) 
      ELSEIF (IFLAG.eq.SURF_DOMAIN) THEN 
         lexternal = .false. 
         linternal = .true. 
      ENDIF 
!                                                                       
      IF (ianz.gt.1) then 
         IF (str_comp (cpara (1) , 'new', 2, lpara (1) , 3) ) then 
            CALL del_params (1, ianz, cpara, lpara, maxw) 
            IF (ianz.eq.1) then 
               string = '('//cpara (ianz) (1:lpara (ianz) ) //')' 
               laenge = lpara (ianz) + 2 
               distance = berechne (string, laenge) 
               DO i = cr_nscat + 1, SURF_MAXSCAT 
               IF (linternal) then 
                  surf_in_dist (i) = distance 
               ENDIF 
               IF (lexternal) then 
                  surf_ex_dist (i) = distance 
               ENDIF 
               ENDDO 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
               RETURN 
            ENDIF 
         ELSE 
            string = '('//cpara (ianz) (1:lpara (ianz) ) //')' 
            laenge = lpara (ianz) + 2 
            distance = berechne (string, laenge) 
            IF (ier_num.ne.0) return 
            ianz = ianz - 1 
            CALL get_iscat (ianz, cpara, lpara, werte, maxw, lold) 
            IF (ier_num.ne.0) return 
            IF (werte (1) .eq. - 1) then 
               DO i = 0, cr_nscat 
               IF (linternal) then 
                  surf_in_dist (i) = distance 
               ENDIF 
               IF (lexternal) then 
                  surf_ex_dist (i) = distance 
               ENDIF 
               ENDDO 
            ELSE 
               DO i = 1, ianz 
               IF (0.le.nint (werte (i) ) .and.nint (werte (i) )        &
               .le.cr_nscat) then                                       
                  IF (linternal) then 
                     surf_in_dist (nint (werte (i) ) ) = distance 
                  ENDIF 
                  IF (lexternal) then 
                     surf_ex_dist (nint (werte (i) ) ) = distance 
                  ENDIF 
               ELSE 
                  ier_num = - 27 
                  ier_typ = ER_APPL 
                  RETURN 
               ENDIF 
               ENDDO 
            ENDIF 
         ENDIF 
      ELSE 
         IF (str_comp (cpara (1) , 'default', 2, lpara (1) , 7) ) then 
            DO i = 0, SURF_MAXSCAT 
            IF (lexternal) then 
               surf_ex_dist (i) = SURF_DIST_DEF 
            ENDIF 
            IF (linternal) then 
               surf_in_dist (i) = SURF_DIST_DEF 
            ENDIF 
            ENDDO 
         ELSEIF (str_comp (cpara (1) , 'off', 2, lpara (1) , 3) ) then 
            DO i = 0, SURF_MAXSCAT 
            IF (lexternal) then 
               surf_ex_dist (i) = 0.0 
            ENDIF 
            IF (linternal) then 
               surf_in_dist (i) = 0.0 
            ENDIF 
            ENDDO 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
      ENDIF 
!                                                                       
      END SUBROUTINE surf_set_fuzzy                 
!*****7*****************************************************************
      SUBROUTINE surf_show 
!+                                                                      
!     This subroutine shows the surface settings                        
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE surface_mod 
      USE prompt_mod 
      USE errlist_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER i 
!                                                                       
      IF (SURF_MAXSCAT==0) THEN
        WRITE(output_io,*) ' No distances to surfaces have been defined yet'
        WRITE(output_io,*) ' Set distances first'
        RETURN
      ENDIF
      WRITE (output_io, * ) 
      WRITE (output_io, 3000) 
      WRITE (output_io, * ) 
      DO i = 0, cr_nscat 
         WRITE (output_io, 3010) cr_at_lis (i), surf_ex_dist (i) 
      ENDDO 
      IF (cr_nscat.lt.SURF_MAXSCAT) then 
         WRITE (output_io, 3010) 'new', surf_ex_dist (SURF_MAXSCAT) 
      ENDIF 
      WRITE (output_io, * ) 
      WRITE (output_io, 3100) 
      WRITE (output_io, * ) 
      DO i = 0, cr_nscat 
         WRITE (output_io, 3010) cr_at_lis (i), surf_in_dist (i) 
      ENDDO 
      IF (cr_nscat.lt.SURF_MAXSCAT) then 
         WRITE (output_io, 3010) 'new', surf_in_dist (SURF_MAXSCAT) 
      ENDIF 
!                                                                       
 3000 FORMAT (' Distances between atom types and an external surface') 
 3100 FORMAT (' Distances between atom types and an internal surface') 
 3010 FORMAT (' Atom type   ',a4,' at ',f8.3,' Angstroem') 
      END SUBROUTINE surf_show                      
!*****7*****************************************************************
      SUBROUTINE property_menu 
!-                                                                      
!     Main menu for property related operations                         
!+                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE molecule_mod 
      USE prop_para_mod 
!
      USE calc_expr_mod
      USE doact_mod 
      USE do_eval_mod
      USE do_wait_mod
      USE errlist_mod 
      USE learn_mod 
      USE class_macro_internal
      USE prompt_mod 
      USE sup_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxw 
!                                                                       
      CHARACTER(5) befehl 
      CHARACTER(LEN=LEN(prompt)) :: orig_prompt
      CHARACTER(1024) line, zeile
      INTEGER lp, length, lbef 
      INTEGER indxg
      LOGICAL lend
!                                                                       
      INTEGER len_str 
      LOGICAL str_comp 
!                                                                       
      maxw = MAXSCAT
      lend = .false. 
      CALL no_error 
      orig_prompt = prompt
      prompt = prompt (1:len_str (prompt) ) //'/prop' 
!                                                                       
      DO while (.not.lend) 
      CALL get_cmd (line, length, befehl, lbef, zeile, lp, prompt) 
      IF (ier_num.eq.0) then 
         IF (line /= ' '      .and. line(1:1) /= '#' .and. &
             line /= char(13) .and. line(1:1) /= '!'        ) THEN
!                                                                       
!     ----search for "="                                                
!                                                                       
indxg = index (line, '=') 
IF (indxg.ne.0.AND..NOT. (str_comp (befehl, 'echo', 2, lbef, 4) ) &
              .AND..NOT. (str_comp (befehl, 'syst', 2, lbef, 4) )    &
              .AND..NOT. (str_comp (befehl, 'help', 2, lbef, 4) .OR. &
                          str_comp (befehl, '?   ', 2, lbef, 4) )    &
              .AND. INDEX(line,'==') == 0                            ) THEN
!                                                                       
!     ------evaluatean expression and assign the value to a variabble   
!                                                                       
               CALL do_math (line, indxg, length) 
            ELSE 
!                                                                       
!------ ----execute a macro file                                        
!                                                                       
               IF (befehl (1:1) .eq.'@') then 
                  IF (length.ge.2) then 
                     CALL file_kdo (line (2:length), length - 1) 
                  ELSE 
                     ier_num = - 13 
                     ier_typ = ER_MAC 
                  ENDIF 
!                                                                       
!     ----continues a macro 'continue'                                  
!                                                                       
               ELSEIF (str_comp (befehl, 'continue', 2, lbef, 8) ) then 
                  CALL macro_continue (zeile, lp) 
!                                                                       
!------ ----Echo a string, just for interactive check in a macro 'echo' 
!                                                                       
               ELSEIF (str_comp (befehl, 'echo', 2, lbef, 4) ) then 
                  CALL echo (zeile, lp) 
!                                                                       
!      ---Evaluate an expression, just for interactive check 'eval'     
!                                                                       
               ELSEIF (str_comp (befehl, 'eval', 2, lbef, 4) ) then 
                  CALL do_eval (zeile, lp) 
!                                                                       
!     ----exit 'exit'                                                   
!                                                                       
               ELSEIF (str_comp (befehl, 'exit', 2, lbef, 4) ) then 
                  lend = .true. 
!                                                                       
!     ----help 'help','?'                                               
!                                                                       
      ELSEIF (str_comp (befehl, 'help', 2, lbef, 4) .or.str_comp (befehl&
     &, '?   ', 1, lbef, 4) ) then                                      
                  IF (str_comp (zeile, 'errors', 2, lp, 6) ) then 
                     lp = lp + 7 
                     CALL do_hel ('discus '//zeile, lp) 
                  ELSE 
                     lp = lp + 16 
                     CALL do_hel ('discus property '//zeile, lp) 
                  ENDIF 
!                                                                       
!------- -Operating System Kommandos 'syst'                             
!                                                                       
               ELSEIF (str_comp (befehl, 'syst', 2, lbef, 4) ) then 
                  IF (zeile.ne.' ') then 
                     CALL do_operating (zeile (1:lp), lp) 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
!                                                                       
!------  -----waiting for user input                                    
!                                                                       
               ELSEIF (str_comp (befehl, 'wait', 3, lbef, 4) ) then 
                  CALL do_input (zeile, lp) 
!
!     ----Set/clear properties for types or (individual) atoms 'set'         
!
               ELSEIF (str_comp (befehl, 'clear', 2, lbef, 5) ) then 
                  CALL property_set_clr (zeile, lp, .false.) 
!
!     ----Set/clear properties for types or (individual) atoms 'set'         
!
               ELSEIF (str_comp (befehl, 'set', 2, lbef, 3) ) then 
                  CALL property_set_clr (zeile, lp, .true.) 
!                                                                       
!                                                                       
!     ----Define which properties have to be present 'property'         
!                                                                       
               ELSEIF (str_comp (befehl, 'property', 2, lbef, 8) ) then 
                  CALL property_select (zeile, lp, cr_sel_prop) 
!                                                                       
!     ----show current parameters 'show'                                
!                                                                       
               ELSEIF (str_comp (befehl, 'show', 2, lbef, 4) ) then 
                  CALL property_show 
               ELSE 
                  ier_num = - 8 
                  ier_typ = ER_COMM 
               ENDIF 
            ENDIF 
         ENDIF 
      ENDIF 
      IF (ier_num.ne.0) THEN 
         CALL errlist 
         IF (ier_sta.ne.ER_S_LIVE) THEN 
            IF (lmakro .OR. lmakro_error) THEN  ! Error within macro or termination errror
               IF(sprompt /= prompt ) THEN
                  ier_num = -10
                  ier_typ = ER_COMM
                  ier_msg(1) = ' Error occured in property menu'
                  prompt_status = PROMPT_ON 
                  prompt = orig_prompt
                  RETURN
               ELSE
                  IF(lmacro_close) THEN
                     CALL macro_close 
                     prompt_status = PROMPT_ON 
                  ENDIF 
               ENDIF 
            ENDIF 
            IF (lblock) THEN 
               ier_num = - 11 
               ier_typ = ER_COMM 
               prompt_status = PROMPT_ON 
               prompt = orig_prompt
               RETURN 
            ENDIF 
            CALL no_error 
            lmakro_error = .FALSE.
            sprompt = ' '
         ENDIF 
      ENDIF 
      ENDDO 
!
      prompt = orig_prompt
!                                                                       
      END SUBROUTINE property_menu                  
!*****7*****************************************************************
      SUBROUTINE property_set_clr(zeile, lp, set_clr)
!+                                                                      
!     This subroutine sets or clears the property flag of an atom or atom type
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE prop_para_mod 
!
      USE ber_params_mod
      USE errlist_mod 
      USE get_params_mod
      IMPLICIT none 
!
      CHARACTER (LEN=*), INTENT(INOUT) :: zeile
      INTEGER          , INTENT(INOUT) :: lp
      LOGICAL          , INTENT(IN   ) :: set_clr

      INTEGER, PARAMETER                   :: TYPES = 0
      INTEGER, PARAMETER                   :: ATOMS = 1
      INTEGER, PARAMETER                   :: MAXW = 200 
      LOGICAL, PARAMETER                   :: LOLD = .true.
      CHARACTER(LEN=1024), DIMENSION(MAXW) :: cpara
      INTEGER            , DIMENSION(MAXW) :: lpara
      REAL               , DIMENSION(MAXW) :: werte
      INTEGER                              :: ianz
      INTEGER                              :: i,is
      INTEGER                              :: ibit_nr = PROP_NORMAL
      INTEGER                              :: sel_mode
      INTEGER                              :: istart, iend
      LOGICAL, DIMENSION(:), ALLOCATABLE   :: latom
!
      LOGICAL str_comp
!
      CALL get_params (zeile, ianz, cpara, lpara, MAXW, lp) 
      IF (ier_num.ne.0) return 
      IF(str_comp(cpara(1)(1:lpara(1)),'domain', 3, lpara(1), 6)) THEN
         ibit_nr = PROP_DOMAIN
      ELSEIF(str_comp(cpara(1)(1:lpara(1)),'outside', 3, lpara(1), 7)) THEN
         ibit_nr = PROP_OUTSIDE
      ELSEIF(str_comp(cpara(1)(1:lpara(1)),'external', 3, lpara(1), 7)) THEN
         ibit_nr = PROP_SURFACE_EXT
      ELSEIF(str_comp(cpara(1)(1:lpara(1)),'internal', 3, lpara(1), 7)) THEN
         ibit_nr = PROP_SURFACE_INT
      ELSEIF(str_comp(cpara(1)(1:lpara(1)),'ligand', 3, lpara(1), 6)) THEN
         ibit_nr = PROP_LIGAND
      ENDIF
      IF(str_comp(cpara(2)(1:lpara(2)),'types', 3, lpara(2), 5)) THEN
         sel_mode = TYPES
      ELSEIF(str_comp(cpara(2)(1:lpara(2)),'atoms', 3, lpara(2), 5)) THEN
         sel_mode = ATOMS
      ELSE
         ier_num = -6
         ier_typ = ER_COMM
         RETURN
      ENDIF
      CALL del_params (2, ianz, cpara, lpara, MAXW) 
      IF (ier_num.ne.0) return 
      ALLOCATE(latom(0:MAXSCAT))
      typesel: IF(sel_mode == TYPES) THEN
         CALL get_iscat (ianz, cpara, lpara, werte, maxw, lold) 
         IF (ier_num.ne.0) return 
         IF(werte(1)==-1) THEN
            latom = .true.
         ELSE
            latom = .false.
            DO i = 1, ianz 
               is = nint (werte (i) ) 
               IF (is.ge.0.and.is.le.cr_nscat) then 
                  latom (is) = .true. 
               ELSE 
                  ier_num = - 27 
                  ier_typ = ER_APPL 
               ENDIF 
            ENDDO 
         ENDIF
         IF(set_clr) THEN
            DO i=1, cr_natoms
               IF(latom(cr_iscat(i))) THEN
                  cr_prop(i) = IBSET(cr_prop(i),ibit_nr)
               ENDIF
            ENDDO
         ELSE
            DO i=1, cr_natoms
               IF(latom(cr_iscat(i))) THEN
                  cr_prop(i) = IBCLR(cr_prop(i),ibit_nr)
               ENDIF
            ENDDO
         ENDIF
      ELSEIF(sel_mode == ATOMS) THEN typesel
         IF(str_comp(cpara(1)(1:lpara(1)),'all', 3, lpara(1), 3)) THEN
             istart=1
             iend  = cr_natoms
         ELSE
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
            IF (ier_num.ne.0) return 
            IF(ianz==1) THEN
               istart=nint(werte(1))
               iend  =nint(werte(1))
            ELSEIF(ianz==2) THEN
               istart=nint(werte(1))
               iend  =nint(werte(2))
            ELSE
               ier_num = -6
               ier_typ = ER_COMM
               RETURN
            ENDIF
         ENDIF
         IF(set_clr) THEN
            DO i=istart, iend
               cr_prop(i) = IBSET(cr_prop(i),ibit_nr)
            ENDDO
         ELSE
            DO i=istart, iend
               cr_prop(i) = IBCLR(cr_prop(i),ibit_nr)
            ENDDO
         ENDIF
      ENDIF typesel
!
      DEALLOCATE(latom)
!
      END SUBROUTINE property_set_clr
!*****7*****************************************************************
      SUBROUTINE property_show 
!+                                                                      
!     This subroutine shows the property settings                       
!-                                                                      
      USE discus_config_mod 
      USE crystal_mod 
      USE prop_para_mod 
!
      USE errlist_mod 
      USE prompt_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      CHARACTER(32) c_property 
      INTEGER length 
!                                                                       
      WRITE (output_io, 2000) 
      CALL char_prop_2 (c_property, cr_sel_prop (1), cr_sel_prop (0),   &
      length)                                                           
      WRITE (output_io, 2100) c_property (1:length) 
!                                                                       
      CALL surf_show 
!                                                                       
 2000 FORMAT     (/' Property related settings'/) 
                                                                        
 2100 FORMAT    (/' Atom properties              ',/,                   &
     &                   '   N = normal atom            ',/,            &
     &                   '   M = atom in a molecule     ',/,            &
     &                   '   D = atom in a domain       ',/,            &
     &                   '   O = atom outside boundaries',/,            &
     &                   '   E = atom near ext. surface ',/,            &
     &                   '   I = atom near int. surface ',/,            &
     &                   '   L = atom in ligand molecule',/,            &
     &                   '                              : ','NMDOEIL'/, &
     &                   '      absent=- ignored=.      : ',a)          
!                                                                       
      END SUBROUTINE property_show                  
!*****7*****************************************************************
      SUBROUTINE char_prop_1 (c_property, property, length) 
!-                                                                      
!     sets letters for true property bits                               
!     property = 1; ==> Letter, i.e. property is present                
!     property = 0; ==> "-"     i.e. property is absent                 
!+                                                                      
      USE prop_para_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      CHARACTER ( * ) c_property 
      INTEGER property 
      INTEGER length 
!                                                                       
      INTEGER i 
!                                                                       
      DO i = 0, MAXPROP 
      IF (btest (property, i) ) then 
         c_property (i + 1:i + 1) = c_prop_letter (i + 1:i + 1) 
      ELSEIF (.not.btest (property, i) ) then 
         c_property (i + 1:i + 1) = '-' 
      ENDIF 
      ENDDO 
      length = MAXPROP 
      END SUBROUTINE char_prop_1                    
!*****7*****************************************************************
      SUBROUTINE char_prop_2 (c_property, property, mask, length) 
!-                                                                      
!     sets letters for true property bits                               
!     property ; mask = 1; 1 ==> Letter, i.e. property is present       
!     property ; mask = 0; 1 ==> "-" i.e. property must be absent       
!     property ; mask = 0; 0 ==> "." i.e. property is ignored           
!     property ; mask = 1; 0 ==> "." i.e. property is ignored           
!      This last situation should never occur...                        
!+                                                                      
      USE prop_para_mod 
      IMPLICIT none 
!                                                                       
!                                                                       
      CHARACTER ( * ) c_property 
      INTEGER property 
      INTEGER mask 
      INTEGER length 
!                                                                       
      INTEGER i 
!                                                                       
      DO i = 0, MAXPROP 
      IF (btest (property, i) .and.btest (mask, i) ) then 
         c_property (i + 1:i + 1) = c_prop_letter (i + 1:i + 1) 
      ELSEIF (.not.btest (property, i) .and.btest (mask, i) ) then 
         c_property (i + 1:i + 1) = '-' 
      ELSEIF (.not.btest (property, i) .and..not.btest (mask, i) ) then 
         c_property (i + 1:i + 1) = '.' 
      ELSEIF (btest (property, i) .and..not.btest (mask, i) ) then 
         c_property (i + 1:i + 1) = '.' 
      ENDIF 
      ENDDO 
      length = MAXPROP 
      END SUBROUTINE char_prop_2                    
!*****7*****************************************************************
      SUBROUTINE get_iscat (ianz, cpara, lpara, werte, maxw, lnew) 
!-                                                                      
!     Determines the scattering type of the parameter                   
!+                                                                      
      USE discus_config_mod 
      USE charact_mod
      USE crystal_mod 
      USE berechne_mod
      USE do_variable_mod
      USE errlist_mod 
      USE string_convert_mod
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER maxw 
!                                                                       
      CHARACTER ( * ) cpara (maxw) 
      CHARACTER(1024) zeile 
      INTEGER lpara (maxw) 
      INTEGER i, j, l, ianz, jj, jp 
      LOGICAL lnew 
      REAL werte (maxw) 
!                                                                       
!                                                                       
!     ----Select which atoms are included in the wave                   
!                                                                       
!
!     Attempt to replace a string from a variable
!
      DO i=1, ianz
         zeile = ' '
         zeile = cpara(i)(1:lpara(i))
         j     = lpara(1)
         CALL ersetz_variable (zeile, j)
         jj = LEN_TRIM(zeile)
         IF(ier_num == 0 .AND. jj>2 .AND. zeile(1:1)=='''' .AND. &
            zeile(jj:jj) == ''''                                ) THEN
            cpara(i) = ' '
            cpara(i)(1:jj-2) = zeile(2:jj-1)
         ENDIF
      ENDDO
!
      ier_num = 0 
      ier_typ = ER_NONE 
      DO i = 1, maxw 
      werte (i) = 0.0 
      ENDDO 
      jj = 1 
      jp = 0 
      j = 1 
      DO while (j.le.ianz.and.ier_num.eq.0) 
      i = ichar (cpara (j) (1:1) ) 
      IF (cpara (j) .eq.'all') then 
         werte (1) = - 1 
         RETURN 
      ELSEIF ( ( (a.le.i.and.i.le.z) .or. (aa.le.i.and.i.le.zz) ) .and. &
      (index (cpara (j) , '[') .eq.0) ) then                            
         CALL do_cap (cpara (j) ) 
         ier_num = - 27 
         ier_typ = ER_APPL 
         DO i = 0, cr_nscat 
         IF (cpara (j)(1:lpara(j)) .eq.cr_at_lis (i) ) then 
            werte (jj) = i 
            jj = jj + 1 
            jp = jp + 1 
            ier_num = 0 
            ier_typ = ER_NONE 
         ENDIF 
         ENDDO 
         IF (lnew) then 
            IF (j.gt.1) then 
               IF (cr_nscat.lt.MAXSCAT) then 
                  cr_nscat = cr_nscat + 1 
                  cr_at_lis (cr_nscat) = cpara (j) 
                  cr_dw (cr_nscat) = 0.0 
                  werte (jj) = cr_nscat 
                  ier_num = 0 
                  ier_typ = ER_NONE 
               ELSE 
                  ier_num = - 26 
                  ier_typ = ER_APPL 
               ENDIF 
            ENDIF 
         ENDIF 
      ELSE 
         zeile = ' ' 
         l = lpara (j) 
         zeile (1:1) = '(' 
         zeile (2:l + 1) = cpara (j) (1:lpara(j))
         zeile (l + 2:l + 2) = ')' 
         l = l + 2 
         werte (jj) = berechne (zeile, l) 
         IF (ier_num.eq.0) then 
            IF (0.le.nint (werte (jj) ) .and.nint (werte (jj) )         &
            .le.cr_nscat) then                                          
               jj = jj + 1 
               ier_num = 0 
               ier_typ = ER_NONE 
            ELSE 
               ier_num = - 27 
               ier_typ = ER_APPL 
            ENDIF 
         ENDIF 
      ENDIF 
      j = j + 1 
      ENDDO 
      ianz = max (ianz, jp) 
!                                                                       
      END SUBROUTINE get_iscat                      
END MODULE modify_mod
