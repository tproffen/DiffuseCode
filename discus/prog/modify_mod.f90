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
      USE config_mod 
      USE allocate_appl_mod
      USE crystal_mod 
      USE modify_func_mod
      USE molecule_mod 
      USE prop_para_mod 
      USE structur, ONLY: update_cr_dim
      IMPLICIT none 
!                                                                       
       
      include'random.inc' 
      include'errlist.inc' 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 200) 
!                                                                       
      CHARACTER ( * ) zeile 
      CHARACTER(1024) cpara (maxw), cc 
      INTEGER lpara (maxw) 
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
                  cr_at_lis (cr_nscat) = cpara (2) 
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
!                                                                       
            IF (ianz.eq.3.and.ier_num.eq. - 27) then 
               IF (cr_nscat + 1 >  MAXSCAT) then 
                  n_scat = cr_nscat + 1
                  call alloc_crystal (n_scat, NMAX)
               ENDIF
               IF (cr_nscat + 1.le.maxscat) then 
                  cr_nscat = cr_nscat + 1 
                  cr_at_lis (cr_nscat) = cpara (2) 
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
                     cr_at_lis (cr_nscat) = cc 
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
      USE config_mod 
      USE crystal_mod 
      USE molecule_mod 
      IMPLICIT none 
!                                                                       
       
      include'errlist.inc' 
!                                                                       
      INTEGER idest, isource, ityp, i, k, ii, jj, is, js, i0, j0 
      INTEGER ip, jp 
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
      ip = cr_prop (ii) 
      jp = cr_prop (jj) 
!                                                                       
      cr_iscat (ii) = js 
      cr_prop (ii) = jp 
      DO i = 1, 3 
      ipos (i) = cr_pos (i, ii) 
      cr_pos (i, ii) = cr_pos (i, jj) - j0pos (i) + i0pos (i) 
      ENDDO 
!                                                                       
      IF (lswap) then 
         cr_iscat (jj) = is 
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
      USE config_mod 
      USE allocate_appl_mod
      USE crystal_mod 
      IMPLICIT none 
!                                                                       
       
      include'charact.inc' 
      include'errlist.inc' 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 15) 
!                                                                       
      CHARACTER (LEN=*), INTENT(IN) :: line
      INTEGER          , INTENT(IN) :: laenge
!
      CHARACTER(LEN=4)              :: name 
      CHARACTER(LEN=1024), DIMENSION(maxw) :: cpara
      CHARACTER(LEN=1024)                  :: zeile 
      INTEGER            , DIMENSION(maxw) :: lpara
      INTEGER i, j, ianz, lp 
      INTEGER                              :: new_nmax
      INTEGER                              :: new_nscat
      REAL werte (maxw) 
      REAL berechne 
!
!     While developing, increment crystal if neede, but keep the check
!

      IF ( NMAX <= cr_natoms .or. MAXSCAT <= cr_nscat ) then 
         new_nmax  = max(NMAX+1   , INT(NMAX    * 1.25))
         new_nscat = max(MAXSCAT+1, INT(MAXSCAT * 1.25))
         call alloc_crystal(new_nscat, new_nmax)
         IF ( ier_num /= 0) RETURN
      ENDIF
      IF (cr_natoms.lt.nmax) then 
         CALL get_params (line, ianz, cpara, lpara, maxw, laenge) 
         CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
         IF (ier_num.eq.0) then 
            i = ichar (cpara (1) (1:1) ) 
            IF ( (a.le.i.and.i.le.z) .or. (aa.le.i.and.i.le.zz) ) then 
               name = cpara (1) 
               j = 2 
            ELSE 
               j = 1 
               name = 'yyyy' 
            ENDIF 
            DO i = j, ianz 
            lp = lpara (i) 
            zeile = ' ' 
            zeile (1:1) = '(' 
            zeile (2:lp + 1) = cpara (i) (1:lp) 
            zeile (lp + 2:lp + 2) = ')' 
            lp = lp + 2 
            werte (i) = 0 
            werte (i) = berechne (zeile, lp) 
            IF (ier_num.ne.0) return 
            ENDDO 
            IF (ianz.ge.4) then 
               CALL do_ins_atom (name, werte, ianz) 
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
      USE config_mod 
      USE allocate_appl_mod
      USE crystal_mod 
      USE prop_para_mod 
      IMPLICIT none 
!                                                                       
       
      include'charact.inc' 
      include'errlist.inc' 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 10) 
!                                                                       
      CHARACTER ( * ) line 
      CHARACTER(4) name 
      CHARACTER(1024) cpara (maxw), zeile 
      INTEGER lpara (maxw) 
      INTEGER i, j, k, ianz, laenge, lp 
      INTEGER                              :: new_nmax
      INTEGER                              :: new_nscat
      LOGICAL lkick, lspace 
      REAL werte (maxw) 
      REAL w (3), v (3) 
      REAL berechne 
      REAL do_blen 
!                                                                       
      DATA lspace / .true. / 
!                                                                       
!
!     While developing, increment crystal if neede, but keep the check
!

      IF ( NMAX <= cr_natoms .or. MAXSCAT <= cr_nscat ) then 
         new_nmax  = max(NMAX+1   , INT(NMAX    * 1.25))
         new_nscat = max(MAXSCAT+1, INT(MAXSCAT * 1.25))
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
                  name = cpara (1) 
                  j = 2 
               ELSE 
                  j = 1 
                  name = 'yyyy' 
               ENDIF 
!                                                                       
!     ------Evaluate parameters, old form will be replaced by ber_params
!                                                                       
               DO i = j, ianz 
               lp = lpara (i) 
               zeile = ' ' 
               zeile (1:1) = '(' 
               zeile (2:lp + 1) = cpara (i) (1:lp) 
               zeile (lp + 2:lp + 2) = ')' 
               lp = lp + 2 
               werte (i) = 0 
               werte (i) = berechne (zeile, lp) 
               IF (ier_num.ne.0) return 
               ENDDO 
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
                     CALL do_ins_atom (name, werte, ianz) 
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
                     CALL do_ins_atom (name, werte, ianz) 
                  ENDIF 
               ELSEIF (cr_natoms.eq.0) then 
!                                                                       
!     --------The crystal is empty. Insert atom in any case             
!                                                                       
                  CALL do_ins_atom (name, werte, ianz) 
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
      SUBROUTINE do_ins_atom (name, werte, ianz) 
!-                                                                      
!     Inserts the atom given by name and position in werte into the     
!     structure.                                                        
!+                                                                      
      USE config_mod 
      USE allocate_appl_mod
      USE crystal_mod 
      USE prop_para_mod 
      IMPLICIT none 
!                                                                       
      include'errlist.inc' 
!                                                                       
      INTEGER, PARAMETER     :: maxw  = 5
!                                                                       
      CHARACTER (LEN=* )    , INTENT(IN) :: name 
      REAL , DIMENSION(maxw), INTENT(IN) :: werte (maxw) 
      INTEGER               , INTENT(IN) :: ianz 
!
      INTEGER                :: i, l
      INTEGER                :: new_nmax
      INTEGER                :: new_nscat
      LOGICAL                :: lda 
!                                                                       
!
!     While developing, increment crystal if needed, but keep the check
!
      IF ( NMAX <= cr_natoms .or. MAXSCAT <= cr_nscat ) then 
         new_nmax  = max(NMAX+1   , INT(NMAX    * 1.25))
         new_nscat = max(MAXSCAT+1, INT(MAXSCAT * 1.25))
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
            cr_prop (cr_natoms) = 0 
      cr_prop (cr_natoms)  = ibset (cr_prop (cr_natoms),  PROP_NORMAL) 
            DO l = 1, 3 
            cr_dim (l, 1) = amin1 (cr_dim (l, 1), cr_pos (l, cr_natoms) &
            )                                                           
            cr_dim (l, 2) = amax1 (cr_dim (l, 2), cr_pos (l, cr_natoms) &
            )                                                           
            ENDDO 
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
      USE config_mod 
      USE crystal_mod 
      USE molecule_mod 
      USE prop_para_mod 
      IMPLICIT none 
!                                                                       
       
      include'errlist.inc' 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 6) 
!                                                                       
      CHARACTER ( * ) line 
      CHARACTER(1024) cpara (maxw) 
      INTEGER lpara (maxw) 
      INTEGER i, j 
      INTEGER ll, istart, iend, ianz 
      INTEGER tstart, tend 
      LOGICAL l_mole 
      LOGICAL l_type 
      REAL werte (maxw) 
!                                                                       
      LOGICAL str_comp 
!                                                                       
      l_mole = .false. 
!                                                                       
!     Get parameters                                                    
!                                                                       
      CALL get_params (line, ianz, cpara, lpara, maxw, ll) 
      IF (ier_num.eq.0) then 
!                                                                       
!     --Remove atoms or molecules ?                                     
!                                                                       
         l_mole = .false. 
         tstart = 1 
         tend = mole_num_type 
         IF (ianz.gt.2.and.str_comp (cpara (3) , 'molecule', 1, lpara ( &
         3) , 8) ) then                                                 
            cpara (3) = '0' 
            lpara (3) = 1 
            l_mole = .true. 
         ELSEIF (ianz.gt.1.and.str_comp (cpara (2) , 'molecule', 1,     &
         lpara (2) , 8) ) then                                          
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
         IF (ianz.gt.4.and.str_comp (cpara (4) , 'type', 1, lpara (4) , &
         8) ) then                                                      
            cpara (4) = '0' 
            lpara (4) = 1 
            l_type = .true. 
         ENDIF 
!                                                                       
!     --Calculate value of parameters                                   
!                                                                       
         IF (str_comp (cpara (1) , 'all', 1, lpara (1) , 3) ) then 
!                                                                       
!     ----Remove all molecules or all atoms                             
!                                                                       
            WRITE (cpara (1), 3000) 1 
            IF (l_mole) then 
               WRITE (cpara (2), 3000) mole_num_mole 
            ELSE 
               WRITE (cpara (2), 3000) cr_natoms 
               IF (ianz.eq.1) ianz = 2 
            ENDIF 
            lpara (1) = 11 
            lpara (2) = 11 
         ELSEIF (str_comp (cpara (1) , 'last', 1, lpara (1) , 4) ) then 
!                                                                       
!     ----Remove last molecule or last atom                             
!                                                                       
            IF (l_mole) then 
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
         IF (str_comp (cpara (5) , 'all', 1, lpara (5) , 3) ) then 
!                                                                       
!     ----Remove all types of molecules or all atoms                    
!                                                                       
            WRITE (cpara (5), 3000) 1 
            WRITE (cpara (6), 3000) mole_num_type 
            lpara (5) = 11 
            lpara (6) = 11 
         ELSEIF (str_comp (cpara (5) , 'last', 1, lpara (5) , 4) ) then 
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
         IF (ier_num.eq.0) then 
            istart = nint (werte (1) ) 
            IF (ianz.eq.1) then 
               iend = istart 
            ELSE 
               iend = nint (werte (2) ) 
            ENDIF 
!                                                                       
!     ----Remove molecules                                              
!                                                                       
            IF (l_mole) then 
               IF (                                                     &
               0.lt.istart.and.istart.le.iend.and.iend.le.mole_num_mole)&
               then                                                     
!                                                                       
!     --------Set limits for types that can be removed                  
!                                                                       
                  IF (ianz.eq.5) then 
                     tstart = werte (5) 
                     tend = werte (5) 
                  ELSEIF (ianz.eq.6) then 
                     tstart = werte (5) 
                     tend = werte (6) 
                  ENDIF 
      IF (0.lt.tstart.and.tstart.le.tend.and.tend.le.mole_num_type) then 
!                                                                       
!     ----------Loop over all molecules in range                        
!                                                                       
                     DO i = istart, iend 
                     IF (tstart.le.mole_type (i) .and.mole_type (i)     &
                     .le.tend) then                                     
!                                                                       
!     --------------Set molecule type to zero and remove atoms          
!                                                                       
                        mole_type (i) = 0 
                        mole_char (i) = 0 
                        DO j = 1, mole_len (i) 
                        cr_iscat (mole_cont (mole_off (i) + j) )        &
                        = 0                                             
                        cr_prop (mole_cont (mole_off (i) + j) ) = ibclr &
                        (cr_prop (mole_cont (mole_off (i) + j) ),       &
                        PROP_NORMAL)                                    
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
            ELSE 
!                                                                       
!     ------If index of atom is within limits, remove atom by           
!           Setting its scattering curve to zero                        
!                                                                       
               IF (0.lt.istart.and.istart.le.iend.and.iend.le.cr_natoms)&
               then                                                     
                  DO i = istart, iend 
                  cr_iscat (i) = 0 
                  cr_prop (i) = ibclr (cr_prop (i), PROP_NORMAL) 
                  ENDDO 
               ELSE 
                  ier_num = - 19 
                  ier_typ = ER_APPL 
               ENDIF 
            ENDIF 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
      ENDIF 
!                                                                       
 3000 FORMAT    (i11) 
!                                                                       
      END SUBROUTINE do_remove                      
!*****7**************************************************************   
      SUBROUTINE do_purge 
!-                                                                      
!     Purges the list of atoms from all deleted atoms                   
!+                                                                      
      USE config_mod 
      USE molecule_mod 
      IMPLICIT none 
!                                                                       
      include'errlist.inc' 
!                                                                       
      IF (mole_num_mole.eq.0) then 
         CALL do_purge_atoms 
      ELSEIF (mole_num_mole.gt.0) then 
         CALL do_purge_molecules 
      ENDIF 
!                                                                       
      END SUBROUTINE do_purge                       
!*****7**************************************************************** 
      SUBROUTINE do_purge_atoms 
!-                                                                      
!     Purges the list of atoms from all deleted atoms                   
!+                                                                      
      USE config_mod 
      USE crystal_mod
      USE molecule_mod 
      IMPLICIT none 
!                                                                       
      include'errlist.inc' 
      include'param.inc' 
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
         cr_prop (i) = cr_prop (ii) 
         ENDDO 
         cr_natoms = cr_natoms - ndel 
         ndel = 0 
!                                                                       
         CALL do_check_purge (.true.) 
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
      USE config_mod 
      USE chem_mod 
      USE pdf_mod 
      IMPLICIT none 
!                                                                       
      include'errlist.inc' 
      include'prompt.inc' 
!                                                                       
      LOGICAL lpurge 
!                                                                       
      IF (lpurge) then 
         WRITE (output_io, 1000) 
         IF (chem_period (1) .or.chem_period (2) .or.chem_period (3) )  &
         then                                                           
            chem_period (1) = .false. 
            chem_period (2) = .false. 
            chem_period (3) = .false. 
            WRITE (output_io, 1100) 'Peridic boundaries DISABLED ..' 
         ENDIF 
         IF (chem_quick) then 
            chem_quick = .false. 
      WRITE (output_io, 1100) 'Chem. neighbor mode set to EXACT ..' 
         ENDIF 
         IF (.not.pdf_lexact) then 
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
      SUBROUTINE do_purge_molecules 
!-                                                                      
!     Purges the list of atoms from all deleted atoms                   
!+                                                                      
      USE config_mod 
      USE crystal_mod 
      USE molecule_mod 
      IMPLICIT none 
!                                                                       
       
      include'errlist.inc' 
      include'param.inc' 
!                                                                       
      INTEGER i, ii, idel, ndel, mdel 
      INTEGER j, k, l, m 
      INTEGER nmol 
      INTEGER mmol 
      INTEGER ivoid, iatom 
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
!                                                                       
!     ------shift all following atoms one down                          
!                                                                       
            DO ii = i, cr_natoms - idel 
            cr_pos (1, ii) = cr_pos (1, ii + 1) 
            cr_pos (2, ii) = cr_pos (2, ii + 1) 
            cr_pos (3, ii) = cr_pos (3, ii + 1) 
            cr_iscat (ii) = cr_iscat (ii + 1) 
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
      USE config_mod 
      USE crystal_mod 
      IMPLICIT none 
!                                                                       
       
      include'errlist.inc' 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 5) 
!                                                                       
      CHARACTER ( * ) line 
      CHARACTER(80) zeile 
      CHARACTER(1024) cpara (maxw) 
      INTEGER lpara (maxw) 
      INTEGER ianz, i, lp, ind 
      REAL werte (maxw) 
      REAL berechne 
!                                                                       
      zeile = ' ' 
      CALL get_params (line, ianz, cpara, lpara, maxw, lp) 
      IF (ianz.ne.5) then 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
      IF (ier_num.eq.0) then 
         DO i = 2, ianz 
         lp = lpara (i) 
         zeile = ' ' 
         zeile (1:1) = '(' 
         zeile (2:lp + 1) = cpara (i) (1:lp) 
         zeile (lp + 2:lp + 2) = ')' 
         lp = lp + 2 
         werte (i) = berechne (zeile, lp) 
         IF (ier_num.ne.0) return 
         ENDDO 
         ind = int (werte (2) ) 
         IF (0.lt.ind.and.ind.le.cr_natoms.and.ind.le.NMAX) then 
            zeile = ' ' 
            zeile (1:4) = cr_at_lis (cr_iscat (ind) ) 
            zeile (5:5) = ',' 
            IF (cpara (1) (1:1) .eq.'a') then 
               WRITE (zeile (7:57), 3000) (werte (i), i = 3, 5),        &
               cr_dw (cr_iscat (ind) )                                  
               lp = 57 
            ELSEIF (cpara (1) (1:1) .eq.'r') then 
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
      USE config_mod 
      USE crystal_mod 
      USE molecule_mod 
      USE structur, ONLY: update_cr_dim
      IMPLICIT none 
!                                                                       
       
      include'errlist.inc' 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 3) 
!                                                                       
      CHARACTER ( * ) line 
      CHARACTER(1024) cpara (maxw) 
      REAL werte (maxw) 
      INTEGER lpara (maxw) 
      INTEGER i, j, ianz, lp, is 
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
            is = cr_iscat (i) 
            cr_iscat (i) = cr_iscat (j) 
            cr_iscat (j) = is 
            is = cr_prop (i) 
            cr_prop (i) = cr_prop (j) 
            cr_prop (j) = is 
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
      USE config_mod 
      USE crystal_mod 
      USE chem_mod 
      IMPLICIT none 
!                                                                       
      include'charact.inc' 
       
      include'errlist.inc' 
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
               ii = i - (ianz - 4) 
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
                     CALL do_find_mol (ianz, werte, maxw, x, rmin,      &
                     radius, fq, fp)                                    
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
      USE config_mod 
      USE crystal_mod 
      USE atom_env_mod 
      USE modify_func_mod
      IMPLICIT none 
!                                                                       
       
      include'param.inc' 
      include'errlist.inc' 
!                                                                       
      INTEGER ianz, maxw 
      REAL werte (maxw) 
      REAL x (3) 
      REAL rmin, rmax 
      LOGICAL fq, fp (3) 
!                                                                       
      INTEGER i, j, k, ii 
      INTEGER istart (3), iend (3), iii (3), cell (3), iatom 
      REAL offset (3), nooffset (3) 
      LOGICAL ltype 
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
         istart (i) = iii (i) - 1 - int (rmax / cr_a0 (i) ) 
         iend (i) = iii (i) + 1 + int (rmax / cr_a0 (i) ) 
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
            RETURN 
         ENDIF 
!                                                                       
         ltype = atom_allowed (iatom, werte, ianz, maxw)                &
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
         IF (fp (1) .or.fp (2) .or.fp (3) ) then 
            ier_num = - 16 
            ier_typ = ER_CHEM 
            RETURN 
         ENDIF 
!                                                                       
         DO i = 1, cr_natoms 
         ltype = atom_allowed (i, werte, ianz, maxw)                    &
         .and.check_select_status (.true., cr_prop (i), cr_sel_prop)    
         IF (ltype) then 
            CALL check_blen (x, i, rmin, rmax, nooffset) 
            IF (ier_num.ne.0) return 
         ENDIF 
         ENDDO 
      ENDIF 
!                                                                       
      END SUBROUTINE do_find_env                    
!*****7*****************************************************************
      SUBROUTINE check_bound (cell, offset, fp, ltype) 
!+                                                                      
!     This routine applies periodic boundaries if 'fp' is TRUE.         
!-                                                                      
      USE config_mod 
      USE crystal_mod 
      IMPLICIT none 
!                                                                       
       
      include'errlist.inc' 
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
!     checks if atom 'iatom' is within rmin -> rmax awaz from           
!     position x(3). All matching atoms are stored in the               
!     arrays 'atom_env' and 'res_para'. The positions of the            
!     atoms are stored in 'atom_pos' to retain information              
!     about possible periodic boundaries.                               
!-                                                                      
      USE config_mod 
      USE crystal_mod 
      USE atom_env_mod 
      IMPLICIT none 
!                                                                       
       
      include'param.inc' 
      include'errlist.inc' 
!                                                                       
      REAL x (3), offset (3), rmin, rmax 
      INTEGER iatom 
!                                                                       
      REAL v (3), dist, do_blen 
      INTEGER j 
      LOGICAL lspace 
!                                                                       
      lspace = .true. 
!                                                                       
      DO j = 1, 3 
      v (j) = cr_pos (j, iatom) + offset (j) 
      ENDDO 
      dist = do_blen (lspace, x, v) 
      IF (dist.ge.rmin.and.dist.le.rmax) then 
         IF (atom_env (0) .lt.MAX_ATOM_ENV) then 
            IF (atom_env (0) .lt.MAXPAR_RES) then 
               atom_env (0) = atom_env (0) + 1 
               atom_env (atom_env (0) ) = iatom 
               res_para (atom_env (0) ) = dist 
               res_para (0) = float (nint (res_para (0) ) + 1) 
               DO j = 1, 3 
               atom_pos (j, atom_env (0) ) = v (j) 
               ENDDO 
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
      USE config_mod 
      USE crystal_mod 
      IMPLICIT none 
!                                                                       
       
      include'errlist.inc' 
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
      SUBROUTINE celltoindex (icell, isite, iatom) 
!-                                                                      
!       calculates in which unit cell on which site the atom <ia> is    
!+                                                                      
      USE config_mod 
      USE crystal_mod 
      IMPLICIT none 
!                                                                       
       
      include'errlist.inc' 
!                                                                       
      INTEGER iatom, isite, icell (3) 
!                                                                       
      iatom = ( (icell (3) - 1) * cr_icc (1) * cr_icc (2) +       &
                (icell (2) - 1) * cr_icc (1) + (icell (1) - 1) )  &
              * cr_ncatoms + isite        
!                                                                       
      END SUBROUTINE celltoindex                    
!*****7*****************************************************************
      SUBROUTINE indextocell (iatom, icell, isite) 
!-                                                                      
!       calculates in which unit cell on which site the atom <ia> is    
!+                                                                      
      USE config_mod 
      USE crystal_mod 
      IMPLICIT none 
!                                                                       
       
      include'errlist.inc' 
!                                                                       
      INTEGER iatom, isite, icell (3) 
      INTEGER ia 
!                                                                       
      ia = iatom - 1 
!                                                                       
      icell (3) = int (ia / cr_icc (1) / cr_icc (2) / cr_ncatoms) + 1
      ia = ia - (icell (3) - 1) * cr_icc (1) * cr_icc (2) * cr_ncatoms 
      icell (2) = int (ia / cr_icc (1) / cr_ncatoms) + 1 
      ia = ia - (icell (2) - 1) * cr_icc (1) * cr_ncatoms 
      icell (1) = int (ia / cr_ncatoms) + 1 
      isite = ia - (icell (1) - 1) * cr_ncatoms + 1 
!                                                                       
      END SUBROUTINE indextocell                    
!*****7*****************************************************************
      SUBROUTINE atom_select (zeile, lp, lu, lo, latom, &
                              sel_atom, lold, lselect,  &
                              ival, repl)                             
!+                                                                      
!     This routine executes the select command                          
!-                                                                      
      USE config_mod 
      USE crystal_mod 
      IMPLICIT none 
!                                                                       
       
      include'errlist.inc' 
!                                                                       
      CHARACTER  (LEN=  * ),     INTENT(IN)    :: zeile 
      INTEGER,                   INTENT(IN)    :: lp 
      INTEGER,                   INTENT(IN)    :: lu
      INTEGER,                   INTENT(IN)    :: lo
      LOGICAL, DIMENSION(lu:lo), INTENT(OUT)   :: latom 
      LOGICAL,                   INTENT(INOUT) :: sel_atom
      LOGICAL,                   INTENT(IN)    :: lold
      LOGICAL,                   INTENT(IN)    :: lselect 
      INTEGER,                   OPTIONAL,  INTENT(IN)    :: ival
      INTEGER, DIMENSION(lu:lo), OPTIONAL,  INTENT(OUT)   :: repl (lu:lo)
!
      INTEGER  :: sel_mask (0:1) 
      INTEGER  :: maxw
!                                                                       
      CHARACTER(LEN=1024), DIMENSION(1:lo+1) :: cpara 
      REAL               , DIMENSION(1:lo+1) :: werte 
      INTEGER            , DIMENSION(1:lo+1) :: lpara 
      INTEGER                                :: ianz, i, is, j, k
!                                                                       
      LOGICAL str_comp 
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
                              sel_atom, lold, lselect,  &
                              ival, repl)                             
!+                                                                      
!     This routine exectues the select command                          
!-                                                                      
      USE config_mod 
      USE molecule_mod 
      IMPLICIT none 
!                                                                       
      include'errlist.inc' 
!                                                                       
      CHARACTER  (LEN=  * ),     INTENT(IN)    :: zeile 
      INTEGER,                   INTENT(IN)    :: lp 
      INTEGER,                   INTENT(IN)    :: lu
      INTEGER,                   INTENT(IN)    :: lo
      LOGICAL, DIMENSION(lu:lo), INTENT(OUT)   :: latom 
      LOGICAL,                   INTENT(INOUT) :: sel_atom
      LOGICAL,                   INTENT(IN)    :: lold
      LOGICAL,                   INTENT(IN)    :: lselect 
      INTEGER,                   OPTIONAL,  INTENT(IN)    :: ival
      INTEGER, DIMENSION(lu:lo), OPTIONAL,  INTENT(OUT)   :: repl (lu:lo)
!
      INTEGER,             DIMENSION(0:1)    :: sel_mask
      INTEGER                                :: maxw
!                                                                       
      CHARACTER(LEN=1024), DIMENSION(1:lo+1) :: cpara 
      REAL               , DIMENSION(1:lo+1) :: werte 
      INTEGER            , DIMENSION(1:lo+1) :: lpara 
!                                                                       
      INTEGER                                :: ianz, i, is
!                                                                       
      LOGICAL str_comp 
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
      IMPLICIT none 
!                                                                       
      include'errlist.inc' 
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
      IMPLICIT none 
!                                                                       
      include'errlist.inc' 
!                                                                       
      INTEGER maxw 
      INTEGER ianz 
      CHARACTER ( * ) cpara (maxw) 
      INTEGER lpara (maxw) 
      LOGICAL lselect 
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
      SUBROUTINE do_find_mol (ianz, werte, maxw, x, rmin, rmax, fq, fp) 
!                                                                       
!     This routine finds all molecules around x with a minimal          
!     distance of rmin and a maximum distance of rmax.                  
!-                                                                      
!                                                                       
      USE config_mod 
      USE crystal_mod 
      USE molecule_mod 
      USE mole_env_mod 
      IMPLICIT none 
       
      include'param.inc' 
      include'errlist.inc' 
!                                                                       
      INTEGER ianz, maxw 
      REAL werte (maxw) 
      REAL x (3) 
      REAL rmin, rmax 
      LOGICAL fq, fp (3) 
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
      USE config_mod 
      USE crystal_mod 
      USE molecule_mod 
      USE mole_env_mod 
      IMPLICIT none 
!                                                                       
       
      include'param.inc' 
      include'errlist.inc' 
!                                                                       
      REAL x (3), offset (3), rmin, rmax 
      INTEGER imole 
!                                                                       
      REAL v (3), dist, do_blen 
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
      USE config_mod 
      USE crystal_mod 
      IMPLICIT none 
!                                                                       
       
      include'errlist.inc' 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 6) 
!                                                                       
      CHARACTER ( * ) zeile 
      INTEGER lp 
!                                                                       
      CHARACTER(1024) cpara (maxw) 
      INTEGER lpara (maxw) 
      INTEGER i, ianz 
      LOGICAL lspace 
      LOGICAL linside 
      LOGICAL l_plane 
      LOGICAL l_sphere 
      LOGICAL l_cyl 
      REAL h (3), d, dstar, radius, height 
      REAL v (3) 
      REAL null (3) 
      REAL werte (maxw) 
!                                                                       
      LOGICAL str_comp 
      REAL do_blen 
!                                                                       
      DATA null / 0.0, 0.0, 0.0 / 
!                                                                       
      IF (cr_v.le.0.0) then 
         ier_num = - 35 
         ier_typ = ER_APPL 
         ier_msg (1) = 'A proper unit cell must be defined' 
      ier_msg (2)  = 'for this command to operate       ' 
         RETURN 
      ENDIF 
!                                                                       
      linside = .true. 
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.eq.0) then 
         IF (str_comp (cpara (1) , 'hkl', 3, lpara (1) , 3) ) then 
!                                                                       
!     ----Crystal is limited by a plane                                 
!                                                                       
            IF (str_comp (cpara (ianz) , 'outside', 3, lpara (ianz) , 7)&
            ) then                                                      
               linside = .false. 
               ianz = ianz - 1 
            ELSEIF (str_comp (cpara (ianz) , 'inside', 3, lpara (ianz) ,&
            6) ) then                                                   
               linside = .true. 
               ianz = ianz - 1 
            ENDIF 
            lspace = .false. 
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
            ELSE 
               ier_num = - 6 
               ier_typ = ER_FORT 
               RETURN 
            ENDIF 
            l_plane = .true. 
            l_sphere = .false. 
            l_cyl = .false. 
         ELSEIF (str_comp (cpara (1) , 'sphere', 3, lpara (1) , 6) )    &
         then                                                           
!                                                                       
!     ----Crystal is limited by a sphere                                
!                                                                       
            IF (str_comp (cpara (ianz) , 'outside', 3, lpara (ianz) , 7)&
            ) then                                                      
               linside = .false. 
               ianz = ianz - 1 
            ELSEIF (str_comp (cpara (ianz) , 'inside', 3, lpara (ianz) ,&
            6) ) then                                                   
               linside = .true. 
               ianz = ianz - 1 
            ENDIF 
            lspace = .false. 
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
            l_plane = .false. 
            l_sphere = .true. 
            l_cyl = .false. 
         ELSEIF (str_comp (cpara (1) , 'cylinder', 3, lpara (1) , 8) )  &
         then                                                           
!                                                                       
!     ----Crystal is limited by a cylinder                              
!                                                                       
            IF (str_comp (cpara (ianz) , 'outside', 3, lpara (ianz) , 7)&
            ) then                                                      
               linside = .false. 
               ianz = ianz - 1 
            ELSEIF (str_comp (cpara (ianz) , 'inside', 3, lpara (ianz) ,&
            6) ) then                                                   
               linside = .true. 
               ianz = ianz - 1 
            ENDIF 
            lspace = .false. 
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
            l_plane = .false. 
            l_sphere = .false. 
            l_cyl = .true. 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
         ENDIF 
!                                                                       
         IF (ier_num.eq.0) then 
            IF (l_plane) then 
               DO i = 1, cr_natoms 
               d = 1.0 - cr_pos (1, i) * h (1) - cr_pos (2, i) * h (2)  &
               - cr_pos (3, i) * h (3)                                  
               d = d / dstar 
               CALL boundarize_atom (d, i, linside) 
               ENDDO 
            ELSEIF (l_sphere) then 
               DO i = 1, cr_natoms 
               v (1) = cr_pos (1, i) 
               v (2) = cr_pos (2, i) 
               v (3) = cr_pos (3, i) 
               d = radius - sqrt (v (1) * v (1) * cr_gten (1, 1)        &
               + v (2) * v (2) * cr_gten (2, 2) + v (3) * v (3) *       &
               cr_gten (3, 3) + 2 * v (1) * v (2) * cr_gten (1, 2)      &
               + 2 * v (1) * v (3) * cr_gten (1, 3) + 2 * v (2) * v (3) &
               * cr_gten (2, 3) )                                       
               CALL boundarize_atom (d, i, linside) 
               ENDDO 
            ELSEIF (l_cyl) then 
               DO i = 1, cr_natoms 
               v (1) = cr_pos (1, i) 
               v (2) = cr_pos (2, i) 
               v (3) = 0.0 
               d = radius - sqrt (v (1) * v (1) * cr_gten (1, 1)        &
               + v (2) * v (2) * cr_gten (2, 2) + v (3) * v (3) *       &
               cr_gten (3, 3) + 2 * v (1) * v (2) * cr_gten (1, 2)      &
               + 2 * v (1) * v (3) * cr_gten (1, 3) + 2 * v (2) * v (3) &
               * cr_gten (2, 3) )                                       
               CALL boundarize_atom (d, i, linside) 
               v (1) = 0.0 
               v (2) = 0.0 
               v (3) = cr_pos (3, i) 
               d = height - sqrt (v (1) * v (1) * cr_gten (1, 1)        &
               + v (2) * v (2) * cr_gten (2, 2) + v (3) * v (3) *       &
               cr_gten (3, 3) + 2 * v (1) * v (2) * cr_gten (1, 2)      &
               + 2 * v (1) * v (3) * cr_gten (1, 3) + 2 * v (2) * v (3) &
               * cr_gten (2, 3) )                                       
               CALL boundarize_atom (d, i, linside) 
               ENDDO 
            ENDIF 
         ENDIF 
      ENDIF 
!                                                                       
      END SUBROUTINE boundary                       
!*****7*****************************************************************
      SUBROUTINE boundarize_atom (distance, iatom, linside) 
!-                                                                      
!     This subroutine sets the boundary property of the atom            
!                                                                       
!     Version  : 1.0                                                    
!     Date     : 07 Sep 10                                              
!                                                                       
!     Author   : R.B. Neder (reinhard.neder@krist.uni-erlangen.de)      
!+                                                                      
      USE config_mod 
      USE allocate_appl_mod
      USE crystal_mod 
      USE prop_para_mod 
      USE surface_mod 
      IMPLICIT none 
!                                                                       
       
      include'errlist.inc' 
!                                                                       
      REAL distance 
      INTEGER iatom 
      LOGICAL linside 
!                                                                       
      IF( cr_nscat > SURF_MAXSCAT) THEN
         CALL alloc_surf ( cr_nscat )
         IF ( ier_num < 0 ) THEN
            RETURN
         ENDIF
      ENDIF
!
      IF ( (linside.and.distance.lt.0) .or. (                           &
      .not.linside.and.distance.gt.0) ) then                            
         cr_iscat (iatom) = 0 
         cr_prop (iatom) = ibclr (cr_prop (iatom), PROP_NORMAL) 
         cr_prop (iatom) = ibset (cr_prop (iatom), PROP_OUTSIDE) 
         IF (abs (distance) .lt.surf_ex_dist (cr_iscat (iatom) ) ) then 
            cr_prop (iatom) = ibset (cr_prop (iatom), PROP_SURFACE_EXT) 
         ENDIF 
      ELSE 
         IF (abs (distance) .lt.surf_ex_dist (cr_iscat (iatom) ) ) then 
            cr_prop (iatom) = ibset (cr_prop (iatom), PROP_SURFACE_EXT) 
         ENDIF 
      ENDIF 
!                                                                       
      END SUBROUTINE boundarize_atom                
!*****7*****************************************************************
      SUBROUTINE do_change (line, laenge) 
!-                                                                      
!     These subroutine is the main routine for changing various         
!     parameters.                                                       
!                                                                       
!     Version  : 1.0                                                    
!     Date     : 03 Jul 96                                              
!                                                                       
!     Author   : R.B. Neder (reinhard.neder@mail.uni-wuerzburg.de)      
!+                                                                      
      USE config_mod 
      USE crystal_mod 
      IMPLICIT none 
!                                                                       
      include'prompt.inc' 
       
      include'errlist.inc' 
      include'param.inc' 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 7) 
!                                                                       
      CHARACTER ( * ) line 
      CHARACTER(1024) cpara (maxw) 
      INTEGER lpara (maxw) 
      INTEGER ianz 
      INTEGER i, j 
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
!     Author   : R.B. Neder (reinhard.neder@mail.uni-wuerzburg.de)      
!+                                                                      
!                                                                       
      USE config_mod 
      USE crystal_mod 
      USE molecule_mod 
      IMPLICIT none 
!                                                                       
      include'prompt.inc' 
       
      include'errlist.inc' 
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
!     Author   : R.B. Neder (reinhard.neder@mail.uni-wuerzburg.de)      
!+                                                                      
!                                                                       
      USE config_mod 
      USE molecule_mod 
      IMPLICIT none 
!                                                                       
      include'errlist.inc' 
!                                                                       
      INTEGER idest 
      INTEGER isource 
!                                                                       
!-----      ----Copy all properties                                     
!                                                                       
      mole_char (idest) = mole_char (isource) 
      mole_dens (idest) = mole_dens (isource) 
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
      USE config_mod 
      USE crystal_mod 
      USE molecule_mod 
      USE prop_para_mod 
      IMPLICIT none 
!                                                                       
       
      include'doact.inc' 
      include'errlist.inc' 
      include'learn.inc' 
      include'macro.inc' 
      include'prompt.inc' 
!                                                                       
      INTEGER maxw 
      LOGICAL lnew, lold 
!                                                                       
      PARAMETER (lnew = .true., lold = .false.) 
!                                                                       
      CHARACTER(5) befehl 
      CHARACTER(50) prom 
      CHARACTER(1024) line, zeile, cpara (MAXSCAT) 
      INTEGER lpara (MAXSCAT), lp, length, lbef 
      INTEGER indxg, ianz, i 
      INTEGER indxc 
      LOGICAL lend, lspace 
      LOGICAL lselect 
      REAL werte (MAXSCAT) 
!                                                                       
      INTEGER len_str 
      LOGICAL str_comp 
!                                                                       
      maxw = MAXSCAT
      lend = .false. 
      CALL no_error 
!                                                                       
      DO while (.not.lend) 
      prom = prompt (1:len_str (prompt) ) //'/surf' 
      CALL get_cmd (line, length, befehl, lbef, zeile, lp, prom) 
      IF (ier_num.eq.0) then 
         IF (line.ne.' '.and.line (1:1) .ne.'#') then 
!                                                                       
!     ----search for "="                                                
!                                                                       
            indxg = index (line, '=') 
      IF (indxg.ne.0.and..not. (str_comp (befehl, 'echo', 2, lbef, 4) ) &
     &.and..not. (str_comp (befehl, 'syst', 2, lbef, 4) ) .and..not. (st&
     &r_comp (befehl, 'help', 2, lbef, 4) .or.str_comp (befehl, '?   ', &
     &2, lbef, 4) ) ) then                                              
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
      IF (ier_num.ne.0) then 
         CALL errlist 
         IF (ier_sta.ne.ER_S_LIVE) then 
            IF (lmakro) then 
               CALL macro_close 
               prompt_status = PROMPT_ON 
            ENDIF 
            IF (lblock) then 
               ier_num = - 11 
               ier_typ = ER_COMM 
               RETURN 
            ENDIF 
            CALL no_error 
         ENDIF 
      ENDIF 
      ENDDO 
!                                                                       
 9999 CONTINUE 
!                                                                       
      END SUBROUTINE surface_menu                   
!*****7*****************************************************************
      SUBROUTINE surf_do_set (zeile, length) 
!+                                                                      
!     This subroutine sets various parameters                           
!-                                                                      
      USE config_mod 
      IMPLICIT none 
!                                                                       
      include'errlist.inc' 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 20) 
!                                                                       
      CHARACTER ( * ) zeile 
      INTEGER length 
!                                                                       
      CHARACTER(1024) cpara (maxw) 
      CHARACTER(1024) string 
      INTEGER lpara (maxw) 
      INTEGER laenge 
      INTEGER ianz 
      INTEGER i 
      LOGICAL lold 
      LOGICAL linternal 
      LOGICAL lexternal 
      REAL werte (maxw) 
      REAL distance 
!                                                                       
      LOGICAL str_comp 
      REAL berechne 
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
      USE config_mod 
      USE allocate_appl_mod
      USE crystal_mod 
      USE surface_mod 
      IMPLICIT none 
!                                                                       
       
      include'errlist.inc' 
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
      LOGICAL linternal 
      LOGICAL lexternal 
      REAL distance 
!                                                                       
      LOGICAL str_comp 
      REAL berechne 
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
      USE config_mod 
      USE crystal_mod 
      USE surface_mod 
      IMPLICIT none 
!                                                                       
       
      include'errlist.inc' 
      include'prompt.inc' 
!                                                                       
      INTEGER i 
!                                                                       
      IF (SURF_MAXSCAT==0) THEN
        WRITE(*,*) ' No distances to surfaces have been defined yet'
        WRITE(*,*) ' Set distances first'
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
      USE config_mod 
      USE crystal_mod 
      USE molecule_mod 
      USE prop_para_mod 
      IMPLICIT none 
!                                                                       
       
      include'doact.inc' 
      include'errlist.inc' 
      include'learn.inc' 
      include'macro.inc' 
      include'prompt.inc' 
!                                                                       
      INTEGER maxw 
      LOGICAL lnew, lold 
!                                                                       
      PARAMETER (lnew = .true., lold = .false.) 
!                                                                       
      CHARACTER(5) befehl 
      CHARACTER(50) prom 
      CHARACTER(1024) line, zeile
      INTEGER lp, length, lbef 
      INTEGER indxg, ianz, i 
      INTEGER indxc 
      LOGICAL lend, lspace 
      LOGICAL lselect 
!                                                                       
      INTEGER len_str 
      LOGICAL str_comp 
!                                                                       
      maxw = MAXSCAT
      lend = .false. 
      CALL no_error 
!                                                                       
      DO while (.not.lend) 
      prom = prompt (1:len_str (prompt) ) //'/prop' 
      CALL get_cmd (line, length, befehl, lbef, zeile, lp, prom) 
      IF (ier_num.eq.0) then 
         IF (line.ne.' '.and.line (1:1) .ne.'#') then 
!                                                                       
!     ----search for "="                                                
!                                                                       
            indxg = index (line, '=') 
      IF (indxg.ne.0.and..not. (str_comp (befehl, 'echo', 2, lbef, 4) ) &
     &.and..not. (str_comp (befehl, 'syst', 2, lbef, 4) ) .and..not. (st&
     &r_comp (befehl, 'help', 2, lbef, 4) .or.str_comp (befehl, '?   ', &
     &2, lbef, 4) ) ) then                                              
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
      IF (ier_num.ne.0) then 
         CALL errlist 
         IF (ier_sta.ne.ER_S_LIVE) then 
            IF (lmakro) then 
               CALL macro_close 
               prompt_status = PROMPT_ON 
            ENDIF 
            IF (lblock) then 
               ier_num = - 11 
               ier_typ = ER_COMM 
               RETURN 
            ENDIF 
            CALL no_error 
         ENDIF 
      ENDIF 
      ENDDO 
!                                                                       
 9999 CONTINUE 
!                                                                       
      END SUBROUTINE property_menu                  
!*****7*****************************************************************
      SUBROUTINE property_show 
!+                                                                      
!     This subroutine shows the property settings                       
!-                                                                      
      USE config_mod 
      USE crystal_mod 
      USE prop_para_mod 
      IMPLICIT none 
!                                                                       
       
      include'errlist.inc' 
      include'prompt.inc' 
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
     &                   '                              : ','NMDOEI'/,  &
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
      USE config_mod 
      USE crystal_mod 
      IMPLICIT none 
!                                                                       
       
      include'charact.inc' 
      include'errlist.inc' 
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
      REAL berechne 
!                                                                       
!     ----Select which atoms are included in the wave                   
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
         IF (cpara (j) .eq.cr_at_lis (i) ) then 
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
         zeile (2:l + 1) = cpara (j) 
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
