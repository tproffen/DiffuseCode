MODULE structur

IMPLICIT NONE
!
PUBLIC
CONTAINS
!+                                                                      
!     These routines create a new structure. The different options      
!     are (or will be soon) :                                           
!                                                                       
!       - comletely new structure                                       
!       - freestyle edited                                              
!       - generated from plane,space or non-crystallographic group      
!       - read-in structure                                             
!       - user written file                                             
!       - system catalog of standard structures                         
!                                                                       
!*****7*****************************************************************
      SUBROUTINE read_struc 
!-                                                                      
!     Main menu for read command.                                       
!+                                                                      
      USE config_mod 
      USE allocate_appl_mod
      USE crystal_mod 
      USE molecule_mod 
      USE save_mod 
!      USE interface_def
      IMPLICIT none 
!                                                                       
       
      include'doact.inc' 
      include'errlist.inc' 
      include'learn.inc' 
      include'macro.inc' 
      include'prompt.inc' 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 11) 
!                                                                       
      CHARACTER(1024) line, zeile, cpara (maxw) 
      CHARACTER(1024) strucfile 
      CHARACTER (LEN=1024)   :: string
      CHARACTER(50) prom 
      CHARACTER(5) befehl 
      INTEGER lpara (maxw), lp, length 
      INTEGER ce_natoms, lstr, i, j, k, iatom 
      INTEGER ianz, l, n, lbef 
      INTEGER    :: natoms,nscats
      LOGICAL lout 
      REAL hkl (3), u (3), xc (3), yc (3), zc (3), dist 
      REAL werte (maxw) 
!                                                                       
      INTEGER len_str 
      LOGICAL str_comp 
!                                                                       
      CALL no_error 
!                                                                       
      blank = ' ' 
      lout = .true. 
!                                                                       
      lstr = 40 
      cr_icc (1) = 1 
      cr_icc (2) = 1 
      cr_icc (3) = 1 
!                                                                       
      prom = prompt (1:len_str (prompt) ) //'/read' 
      CALL get_cmd (line, length, befehl, lbef, zeile, lp, prom) 
!                                                                       
      IF (ier_num.ne.0) return 
      IF (line (1:1) .eq.' '.or.line (1:1) .eq.'#'.or.line.eq.char (13) &
      ) goto 9999                                                       
!                                                                       
!------ execute a macro file                                            
!                                                                       
      IF (line (1:1) .eq.'@') then 
         IF (length.ge.2) then 
            CALL file_kdo (line (2:length), length - 1) 
         ELSE 
            ier_num = - 13 
            ier_typ = ER_MAC 
         ENDIF 
!                                                                       
!------ Echo a string, just for interactive check in a macro 'echo'     
!                                                                       
      ELSEIF (str_comp (befehl, 'echo', 2, lbef, 4) ) then 
         CALL echo (zeile, lp) 
!                                                                       
!     execute command                                                   
!     help                                                              
!                                                                       
      ELSEIF (str_comp (befehl, 'help', 2, lbef, 4) .or.str_comp (befehl&
     &, '?   ', 1, lbef, 4) ) then                                      
         IF (zeile.eq.' '.or.zeile.eq.char (13) ) then 
            zeile = 'commands' 
            lp = lp + 8 
         ENDIF 
         IF (str_comp (zeile, 'errors', 2, lp, 6) ) then 
            lp = lp + 7 
            CALL do_hel ('discus '//zeile, lp) 
         ELSE 
            lp = lp + 12 
            CALL do_hel ('discus read '//zeile, lp) 
         ENDIF 
!                                                                       
!------  -----waiting for user input                                    
!                                                                       
      ELSEIF (str_comp (befehl, 'wait', 3, lbef, 4) ) then 
         CALL do_input (zeile, lp) 
!                                                                       
      ELSE 
!                                                                       
!     --all other commands                                              
!                                                                       
!     --get command parameters                                          
!                                                                       
         CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
         IF (ier_num.ne.0) then 
            RETURN 
         ENDIF 
         IF (ianz.ge.1) then 
!                                                                       
!     --Build file name                                                 
!                                                                       
            CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
         ENDIF 
!                                                                       
!     --reset epsilon tensors                                           
!                                                                       
         DO i = 1, 3 
         DO j = 1, 3 
         DO k = 1, 3 
         cr_eps (i, j, k) = 0.0 
         cr_reps (i, j, k) = 0.0 
         ENDDO 
         ENDDO 
         ENDDO 
!                                                                       
!     read in old cell, use space group symbol to generate cell         
!     content 'cell'                                                    
!                                                                       
         IF (str_comp (befehl, 'cell', 1, lbef, 4) .or.str_comp (befehl,&
         'lcell', 1, lbef, 5) ) then                                    
            IF (ianz.ge.1) then 
               cr_newtype = str_comp (befehl, 'cell', 1, lbef, 4) 
               CALL rese_cr 
               strucfile = cpara (1) 
               IF (ier_num.eq.0) then 
!                                                                       
!     --------if necessary get crystal size                             
!                                                                       
                  IF (ianz.gt.1) then 
                     cpara (1) = '0.0' 
                     lpara (1) = 3 
                     CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                     IF (ier_num.eq.0) then 
                        DO i = 1, ianz - 1 
                        cr_icc (i) = nint (werte (i + 1) ) 
                        ENDDO 
                     ENDIF 
                  ENDIF 
                  IF (ier_num.eq.0) then 
                     CALL readcell (strucfile) 
!
                     IF (ier_num.eq.0) then 
!                                                                       
!     ----------check whether total number of atoms fits into available 
!     ----------space                                                   
!                                                                       
                        iatom = cr_icc (1) * cr_icc (2) * cr_icc (3)    &
                        * cr_natoms                                     
                        IF (iatom.gt.nmax) then 
                           CALL alloc_crystal ( MAXSCAT, INT(iatom * 1.1))
                           IF (ier_num < 0 ) RETURN
                        ENDIF
!
!                          ier_num = - 10 
!                          ier_typ = ER_APPL 
!                          WRITE (ier_msg (1), 3000) cr_icc (1),        &
!                          cr_icc (2), cr_icc (3)                       
!                          WRITE (ier_msg (2), 3100) cr_natoms 
!                          WRITE (ier_msg (3), 3200) iatom, nmax 
!3000 FORMAT                ('Unit cells   : ',3(i4,2x)) 
!3100 FORMAT                ('Atoms / cell : ',i10) 
!3200 FORMAT                ('Total / max  : ',i10,'/',i10) 
!                          RETURN 
!                       ENDIF 
                        ce_natoms = cr_natoms 
                        cr_ncatoms = cr_natoms 
                        cr_natoms = 0 
                        DO k = 1, cr_icc (3) 
                        DO j = 1, cr_icc (2) 
                        DO i = 1, cr_icc (1) 
                        DO n = 1, ce_natoms 
                        cr_natoms = cr_natoms + 1 
                        cr_iscat (cr_natoms) = cr_iscat (n) 
                        cr_pos (1, cr_natoms) = cr_pos (1, n) + float ( &
                        i - 1)                                          
                        cr_pos (2, cr_natoms) = cr_pos (2, n) + float ( &
                        j - 1)                                          
                        cr_pos (3, cr_natoms) = cr_pos (3, n) + float ( &
                        k - 1)                                          
                        cr_prop (cr_natoms) = cr_prop (n) 
                        ENDDO 
                        ENDDO 
                        ENDDO 
                        ENDDO 
!                                                                       
!     ----------Update crystal dimensions                               
!                                                                       
                        CALL update_cr_dim 
!                                                                       
!     ----------If molecules were read                                  
!                                                                       
                        IF (mole_num_mole.gt.0) then 
                           IF (mole_num_mole * cr_icc (1) * cr_icc (2)  &
                           * cr_icc (3) .le.MOLE_MAX_MOLE) then         
                              mole_num_atom = mole_off (mole_num_mole)  &
                              + mole_len (mole_num_mole)                
                              l = mole_num_mole 
                              mole_num_unit = mole_num_mole 
                              DO i = 2, cr_icc (1) * cr_icc (2) *       &
                              cr_icc (3)                                
                              DO j = 1, mole_num_mole 
                              l = l + 1 
                              mole_len (l) = mole_len (j) 
                              mole_off (l) = mole_off (l -              &
                              mole_num_mole) + mole_num_atom            
                              mole_type (l) = mole_type (j) 
                              mole_char (l) = mole_char (j) 
                              mole_dens (l) = mole_dens (j) 
                              DO k = 1, mole_len (j) 
                              mole_cont (mole_off (l) + k) = mole_cont (&
                              mole_off (j) + k) + (i - 1) * ce_natoms   
                              ENDDO 
                              ENDDO 
                              ENDDO 
                              mole_num_mole = l 
                           ELSE 
                              ier_num = - 65 
                              ier_typ = ER_APPL 
                              RETURN 
                           ENDIF 
                        ENDIF 
!                                                                       
!     ----------Define initial crystal size in fractional coordinates   
!                                                                       
                        DO l = 1, 3 
                        cr_dim0 (l, 1) = float (nint (cr_dim (l, 1) ) ) 
                        cr_dim0 (l, 2) = float (nint (cr_dim (l, 2) ) ) 
                        ENDDO 
                     ENDIF 
                  ENDIF 
               ENDIF 
            ENDIF 
            IF (ier_num.eq.0) then 
!                                                                       
!     ------reset microdomain status                                    
!                                                                       
               CALL do_stack_rese 
            ENDIF 
!                                                                       
!     Free style editing of a structure 'free'                          
!                                                                       
         ELSEIF (str_comp (befehl, 'free', 1, lbef, 4) ) then 
            CALL rese_cr 
            cr_name = 'freely created structure' 
            cr_spcgr (1:1) = 'P' 
            cr_spcgr (2:2) = '1' 
      cr_spcgr (3:16)  = '              ' 
            cr_spcgrno = 1 
            cr_syst = 1 
            CALL get_symmetry_matrices 
            IF (ianz.eq.0) then 
               DO i = 1, 3 
               cr_a0 (i) = 1.0 
               cr_win (i) = 90.0 
               ENDDO 
            ELSEIF (ianz.eq.6) then 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               DO i = 1, 3 
               cr_a0 (i) = werte (i) 
               cr_win (i) = werte (i + 3) 
               ENDDO 
               IF (cr_a0 (1) .le.0.0.or.cr_a0 (2) .le.0.0.or.cr_a0 (3)  &
               .le.0.0.or.cr_win (1) .le.0.0.or.cr_win (2)              &
               .le.0.0.or.cr_win (3) .le.0.0.or.cr_win (1)              &
               .ge.180.0.or.cr_win (2) .ge.180.0.or.cr_win (3)          &
               .ge.180.0) then                                          
                  ier_num = - 93 
                  ier_typ = ER_APPL 
                  ier_msg (1) = 'Error reading unit cell parameters' 
                  RETURN 
               ENDIF 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
               RETURN 
            ENDIF 
            cr_icc (1) = 1 
            cr_icc (2) = 1 
            cr_icc (3) = 1 
            cr_natoms = 0 
            cr_ncatoms = 1 
            cr_nscat = 0 
            as_natoms = 0 
!                                                                       
!     ----reset microdomain status                                      
!                                                                       
            CALL do_stack_rese 
!                                                                       
!     read an old structure 'stru'                                      
!                                                                       
         ELSEIF (str_comp (befehl, 'stru', 1, lbef, 4) ) then 
            IF (ianz.eq.1) then 
               CALL rese_cr 
               sav_r_ncell = .false. 
               strucfile = cpara (1)
               CALL test_file ( strucfile, natoms, nscats, -1 , .false.)
               IF (ier_num /= 0) THEN
                  RETURN
               ENDIF
               IF(natoms > NMAX .or. nscats > MAXSCAT) THEN
                  nscats = MAX(INT(nscats * 1.1), nscats + 2)
                  natoms = MAX(INT(natoms * 1.1), natoms + 10)
                  CALL alloc_crystal (nscats, natoms)
                  IF ( ier_num /= 0 ) RETURN
               ENDIF
!
               CALL readstru (NMAX, MAXSCAT, strucfile, cr_name,        &
               cr_spcgr, cr_a0, cr_win, cr_natoms, cr_nscat, cr_dw,     &
               cr_at_lis, cr_pos, cr_iscat, cr_prop, cr_dim, as_natoms, &
               as_at_lis, as_dw, as_pos, as_iscat, as_prop, sav_ncell,  &
               sav_r_ncell, sav_ncatoms, spcgr_ianz, spcgr_para)        
               IF (ier_num.ne.0) then 
                  RETURN 
               ENDIF 
!                                                                       
!     ------Define initial crystal size in fractional coordinates       
!                                                                       
               DO l = 1, 3 
               cr_dim0 (l, 1) = float (nint (cr_dim (l, 1) ) ) 
               cr_dim0 (l, 2) = float (nint (cr_dim (l, 2) ) ) 
               ENDDO 
!                                                                       
!     ------The crystal size was read from the structure file           
!                                                                       
               IF (sav_r_ncell) then 
                  DO i = 1, 3 
                  cr_icc (i) = sav_ncell (i) 
                  ENDDO 
                  cr_ncatoms = sav_ncatoms 
               ELSE 
!                                                                       
!     ------Define initial crystal size in number of unit cells         
!                                                                       
                  DO i = 1, 3 
                  cr_icc (i) = max (1, int (cr_dim0 (i, 2) - cr_dim0 (i,&
                  1) ) )                                                
                  ENDDO 
!                                                                       
!     ------Define (average) number of atoms per unit cell              
!                                                                       
                  cr_ncatoms = cr_natoms / (cr_icc (1) * cr_icc (2)     &
                  * cr_icc (3) )                                        
               ENDIF 
!                                                                       
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
            IF (ier_num.eq.0) then 
!                                                                       
!     ------reset microdomain status                                    
!                                                                       
               CALL do_stack_rese 
            ENDIF 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
            GOTO 9999 
         ENDIF 
         IF (ier_num.eq.0) then 
            WRITE (output_io, 1000) cr_spcgr, cr_spcgrno 
!.......calculate metric and reciprocal metric tensor,reciprocal lattice
!       constants and permutation tensors                               
            CALL setup_lattice (cr_a0, cr_ar, cr_eps, cr_gten, cr_reps, &
            cr_rten, cr_win, cr_wrez, cr_v, cr_vr, lout, cr_gmat,       &
            cr_fmat, cr_cartesian)                                      
               CALL get_symmetry_matrices 
         ELSE 
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
      ENDIF 
!                                                                       
 9999 CONTINUE 
!                                                                       
 1000 FORMAT    (1x,a16,i5) 
      END SUBROUTINE read_struc                     
!********************************************************************** 
      SUBROUTINE readcell (strucfile) 
!-                                                                      
!           This subroutine reads a unit cell.                          
!+                                                                      
      USE config_mod 
      USE allocate_appl_mod
      USE crystal_mod 
      USE molecule_mod 
      USE save_mod 
      USE wyckoff_mod
      IMPLICIT none 
!                                                                       
       
      include'errlist.inc' 
!                                                                       
      INTEGER ist, maxw 
      PARAMETER (ist = 7, maxw = 5) 
!                                                                       
      CHARACTER ( LEN=* ), INTENT(IN) :: strucfile 
!
      CHARACTER(10) befehl 
      CHARACTER(1024) line, zeile 
      CHARACTER(80) cpara (maxw) 
      INTEGER lpara (maxw) 
      INTEGER property 
      INTEGER i, j, ibl, lbef 
      INTEGER lline 
      INTEGER     :: new_nmax
      INTEGER     :: new_nscat
      INTEGER     :: io_line
      LOGICAL lread, lcell, lout 
      REAL werte (maxw), dw1 
!                                                                       
      INTEGER len_str 
      LOGICAL str_comp 
      LOGICAL :: IS_IOSTAT_END
!                                                                       
      cr_natoms = 0 
      lread     = .true. 
      lcell     = .true. 
      lout      = .false. 
      CALL test_file ( strucfile, new_nmax, new_nscat, -1 , .not.cr_newtype)
      IF (ier_num /= 0) THEN
         CLOSE (ist)
         RETURN
      ENDIF
      CALL oeffne (ist, strucfile, 'old', lread) 
      IF (ier_num /= 0) THEN
         CLOSE (ist)
         RETURN
      ENDIF
!
      DO i = 1, 3 
         cr_dim (i, 1) =  1.e10 
         cr_dim (i, 2) = -1.e10 
      ENDDO 
!                                                                       
!     --Read header of structure file                                   
!                                                                       
         CALL stru_readheader (ist, NMAX, MAXSCAT, lcell, cr_name,      &
         cr_spcgr, cr_at_lis, cr_nscat, cr_dw, cr_a0, cr_win, sav_ncell,&
         sav_r_ncell, sav_ncatoms, spcgr_ianz, spcgr_para)              
      IF (ier_num.ne.0) THEN 
         RETURN 
      ENDIF 
         CALL setup_lattice (cr_a0, cr_ar, cr_eps, cr_gten, cr_reps,    &
         cr_rten, cr_win, cr_wrez, cr_v, cr_vr, lout, cr_gmat, cr_fmat, &
         cr_cartesian)                                                  
!                                                                       
      IF (ier_num /= 0) THEN
         CLOSE (ist)
         RETURN
      ENDIF
!                                                                       
      CALL get_symmetry_matrices 
      IF( NMAX < spc_n*new_nmax) THEN         ! Allocate sufficient atom numbers
         new_nmax = spc_n*new_nmax + 1
         CALL alloc_crystal(MAXSCAT, new_nmax)
        IF ( ier_num /= 0) THEN
            CLOSE (IST)
            RETURN
         ENDIF
      ENDIF
      IF( MAXSCAT < new_nscat) THEN         ! Allocate sufficient atom numbers
         CALL alloc_crystal(new_nscat, NMAX)
        IF ( ier_num /= 0) THEN
            CLOSE (IST)
            RETURN
         ENDIF
      ENDIF
main: DO  ! while (cr_natoms.lt.nmax)  ! end of loop via EOF in input
         ier_num = -49 
         ier_typ = ER_APPL 
         line = ' ' 
!        READ (ist, 2000, end = 2, err = 999) line 
         READ (ist, 2000, IOSTAT=io_line    ) line 
         IF(IS_IOSTAT_END(io_line)) THEN    ! Handle End Of File
            EXIT main
         ELSEIF(io_line /= 0 ) THEN         ! Handle input error
            GOTO 999
         ENDIF
         lline = len_str (line) 
!23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
empty:   IF (line.ne.' '.and.line (1:1) .ne.'#'.and.line.ne.char (13)) THEN
            IF ( NMAX < cr_natoms + spc_n ) THEN ! Allocate sufficient atom numbers
               new_nmax = MAX(NMAX + spc_n + 1, cr_natoms + spc_n+1)
               CALL alloc_crystal(MAXSCAT, new_nmax)
               IF ( ier_num /= 0) THEN
                  CLOSE (IST)
                  RETURN
               ENDIF
               ier_num = -49 
               ier_typ = ER_APPL 
            ENDIF
            IF ( MAXSCAT < cr_nscat + 1 ) THEN ! Allocate sufficient atom types
               new_nscat = MAX(MAXSCAT + 5, INT ( MAXSCAT * 1.025 ))
               CALL alloc_crystal(new_nscat, NMAX)
               IF ( ier_num /= 0) THEN
                  CLOSE (IST)
                  RETURN
               ENDIF
               ier_num = -49 
               ier_typ = ER_APPL 
            ENDIF
               lbef = 10 
               befehl = ' ' 
               ibl = index (line (1:lline) , ' ') 
               IF (ibl.eq.0) THEN 
                  ibl = lline+1 
               ENDIF 
               lbef = min (ibl - 1, lbef) 
               befehl = line (1:lbef) 
typus:         IF (str_comp (befehl, 'molecule', 4, lbef, 8) .or.       &
                   str_comp (befehl, 'domain',   4, lbef, 6) .or.       &
                   str_comp (befehl, 'object',   4, lbef, 6)     ) THEN
!                                                                       
!     ----------Start/End of a molecule                                 
!                                                                       
                  CALL no_error 
                  IF (ibl.le.lline) then 
                     i = lline-ibl 
                     zeile = line (ibl + 1:lline) 
                  ELSE 
                     zeile = ' ' 
                     i = 0 
                  ENDIF 
                  CALL struc_mole_header (zeile, i, .true.) 
                  IF (ier_num.ne.0) THEN
                     CLOSE(IST)
                     RETURN 
                  ENDIF
               ELSE  typus
                  DO j = 1, MAXW 
                  werte (j) = 0.0 
                  ENDDO 
                  werte (5) = 1.0 
                  CALL read_atom_line (line, ibl, lline, as_natoms,     &
                  maxw, werte)                                          
                  IF (ier_num.ne.0.and.ier_num.ne. -49) then 
                     GOTO 999 
                  ENDIF 
                  cr_natoms = cr_natoms + 1 
                  i = cr_natoms 
                  IF (.not.mole_l_on) then 
!                                                                       
!     ------------Transform atom into first unit cell,                  
!                 if it is not inside a molecule                        
!                                                                       
                     CALL firstcell (werte, maxw) 
                  ENDIF 
!                                                                       
                  DO j = 1, 3 
                     cr_pos (j, i) = werte (j) 
                  ENDDO 
                  dw1 = werte (4) 
                  cr_prop (i) = nint (werte (5) ) 
!                                                                       
                  IF (line (1:4) .ne.'    ') then 
                     ibl = ibl - 1 
                     CALL do_cap (line (1:ibl) ) 
!                                                                       
!------ ----------- New option determines whether all atoms in          
!------ ----------- asymmetric unit are considered different atom       
!------ ----------- types ..                                            
!                                                                       
                     IF (.not.cr_newtype) then 
                        DO j = 0, cr_nscat 
                        IF (line (1:ibl) .eq.cr_at_lis (j)              &
                        .and.dw1.eq.cr_dw (j) ) then                    
                           cr_iscat (i) = j 
                           CALL symmetry 
                           IF (ier_num.ne.0) then 
                              CLOSE (IST)
                              RETURN 
                           ENDIF 
                           GOTO 22 
                        ENDIF 
                        ENDDO 
                     ENDIF 
!                                                                       
!------ ----------- end new code                                        
!                                                                       
                     IF (cr_nscat.lt.maxscat) then 
                        cr_nscat = cr_nscat + 1 
                        cr_iscat (i) = cr_nscat 
                        cr_at_lis (cr_nscat) = line (1:ibl) 
                        cr_dw (cr_nscat) = dw1 
!                                                                       
                        as_natoms = as_natoms + 1 
                        as_at_lis (cr_nscat) = cr_at_lis (cr_nscat) 
                        as_iscat (as_natoms) = cr_iscat (i) 
                        as_dw (as_natoms) = cr_dw (cr_nscat) 
                        DO j = 1, 3 
                        as_pos (j, as_natoms) = cr_pos (j, i) 
                        ENDDO 
                        as_prop (as_natoms) = cr_prop (i) 
                        CALL symmetry 
                        IF (ier_num.ne.0) then 
                           CLOSE(IST)
                           RETURN 
                        ENDIF 
                     ELSE 
                        ier_num = -26 
                        ier_typ = ER_APPL 
                        GOTO 2 
                     ENDIF 
   22                CONTINUE 
                  ENDIF 
               ENDIF  typus
            ENDIF empty
      ENDDO main 
!                                                                       
    2    CONTINUE 
         IF (ier_num.eq. -49) then 
            CALL no_error 
!                                                                       
!       move first unit cell into lower left corner of crystal          
!                                                                       
            DO i = 1, cr_natoms 
            DO j = 1, 3 
            cr_pos (j, i) = cr_pos (j, i) - int ( (cr_icc (j) ) / 2) 
!             cr_pos(j,i)=cr_pos(j,i) - int((cr_icc(j)-0.1)/2)          
            cr_dim (j, 1) = amin1 (cr_dim (j, 1), cr_pos (j, i) ) 
            cr_dim (j, 2) = amax1 (cr_dim (j, 2), cr_pos (j, i) ) 
            ENDDO 
            ENDDO 
         ENDIF 
!                                                                       
!     ENDIF 
!                                                                       
  999 CONTINUE 
      CLOSE (ist) 
      IF (ier_num.eq. - 49) then 
         WRITE (ier_msg (1), 3000) as_natoms + 1 
 3000 FORMAT      ('At atom number = ',i8) 
      ENDIF 
!                                                                       
 2000 FORMAT    (a) 
 2010 FORMAT    (a16) 
      END SUBROUTINE readcell                       
!********************************************************************** 
      SUBROUTINE read_atom_line (line, ibl, length, cr_natoms, maxw,    &
      werte)                                                            
!-                                                                      
!     reads a line from the cell file/structure file                    
!+                                                                      
      USE prop_para_mod 
      IMPLICIT none 
!                                                                       
      include'errlist.inc' 
                                                                        
!                                                                       
      CHARACTER ( * ) line 
      INTEGER ibl 
      INTEGER length 
      INTEGER cr_natoms 
      INTEGER maxw 
      REAL werte (maxw) 
!                                                                       
      CHARACTER(1024) cpara (maxw) 
      CHARACTER(1024) string 
      INTEGER lpara (maxw) 
      INTEGER j 
      INTEGER ianz 
!                                                                       
!DBG      write(*,*) 'length ',length                                   
!DBG      write(*,*) '>',line(1:length),'<'                             
!DBG      write(*,*) 'ibl    ',ibl                                      
!DBG      write(*,*) 'maxw   ',maxw                                     
      werte (5) = 1.0 
!                                                                       
!-----      Deal with old four respectively old five column style       
!                                                                       
!     IF (index (line (ibl:length) , ',') .eq.0) then 
         READ (line (ibl:length), *, end = 999, err = 850) (werte (j),  &
         j = 1, 4)                                                      
!                                                                       
!       got four columns, try to read column five                       
!                                                                       
         READ (line (ibl:length), *, end = 800, err = 850) (werte (j),  &
         j = 1, 5)                                                      
  800    CONTINUE 
         CALL no_error 
         GOTO 900 
!     ENDIF 
!                                                                       
  850 CONTINUE
      IF (index (line (ibl:length) , ',') /=  0) then 
      CALL get_params (line (ibl:length), ianz, cpara, lpara, maxw,     &
      length - ibl + 1)                                                 
      IF (ier_num.eq.0) then 
         IF (ianz.eq.4.or.ianz.eq.5) then 
            DO j = 1, ianz 
            string = '(1.0*'//cpara (j) (1:lpara (j) ) //')' 
            cpara (j) = string 
            lpara (j) = lpara (j) + 6 
            ENDDO 
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
            IF (ier_num.ne.0) then 
               ier_msg (1)  = 'Error calculating atom  ' 
               ier_msg (2) = 'coordinates for atom '//line (1:ibl) 
               WRITE (ier_msg (3), 2000) cr_natoms + 1 
               RETURN 
            ENDIF 
            CALL no_error 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
            ier_msg (1) = 'Missing coordinates for ' 
            ier_msg (2) = 'atom '//line (1:ibl) 
            ier_msg (3) = ' ' 
            RETURN 
         ENDIF 
      ELSE 
         ier_msg (1) = 'Error reading parameters for' 
         ier_msg (2) = 'coordinates for atom '//line (1:ibl) 
         WRITE (ier_msg (3), 2000) cr_natoms + 1 
         RETURN 
      ENDIF 
      ELSE 
         ier_msg (1) = 'Error reading parameters for' 
         ier_msg (2) = 'coordinates for atom '//line (1:ibl) 
         WRITE (ier_msg (3), 2000) cr_natoms + 1 
         RETURN 
      ENDIF 
      RETURN 
!                                                                       
  900 CONTINUE 
!                                                                       
!-----      Basic error checks TO FOLLOW                                
!                                                                       
      IF (nint (werte (5) ) < 1              .or.  &
          2**(MAXPROP+1)-1  < nint (werte (5)    )  )THEN
         ier_num = - 102 
         ier_typ = ER_APPL 
      ENDIF 
!                                                                       
  999 CONTINUE 
!                                                                       
 2000 FORMAT    ('Atom Nr. ',i4) 
      END SUBROUTINE read_atom_line                 
!********************************************************************** 
      SUBROUTINE struc_mole_header (zeile, lp, lcell) 
!-                                                                      
!     interprets the 'molecule' lines of a structure file               
!+                                                                      
                                                                        
      USE config_mod 
      USE crystal_mod 
      USE molecule_mod 
      IMPLICIT none 
!                                                                       
       
      include'errlist.inc' 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 21) 
!                                                                       
      CHARACTER ( * ) zeile 
      CHARACTER(1024) cpara (maxw) 
      INTEGER j, ianz 
      INTEGER lpara (maxw), lp 
      REAL werte (maxw) 
      LOGICAL lcell 
!                                                                       
      LOGICAL str_comp 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
                                                                        
      IF (ier_num.eq.0) then 
         IF (ianz.eq.0) then 
!                                                                       
!     --No parameters, start a new Molekule                             
!                                                                       
            IF (mole_num_mole.lt.MOLE_MAX_MOLE) then 
               IF (mole_num_type.lt.MOLE_MAX_TYPE) then 
                  mole_l_on = .true. 
                  mole_l_first = .true. 
                  mole_num_atom = mole_off (mole_num_mole) + mole_len ( &
                  mole_num_mole)                                        
                  mole_num_mole = mole_num_mole+1 
                  mole_num_curr = mole_num_mole 
                  mole_num_type = mole_num_type+1 
                  mole_off (mole_num_mole) = mole_num_atom 
                  mole_type (mole_num_mole) = mole_num_type 
                  mole_gene_n = 0 
                  mole_symm_n = 0 
               ELSE 
                  ier_num = - 66 
                  ier_typ = ER_APPL 
                  RETURN 
               ENDIF 
            ELSE 
               ier_num = - 65 
               ier_typ = ER_APPL 
               RETURN 
            ENDIF 
         ELSE 
!                                                                       
!     --Parameters, interpret parameters                                
!                                                                       
            IF (str_comp (cpara (1) , 'end', 3, lpara (1) , 3) ) then 
!                                                                       
!     ----Turn off molecule                                             
!                                                                       
               IF (lcell) call mole_firstcell 
               mole_l_on = .false. 
!                                                                       
            ELSEIF (str_comp (cpara (1) , 'character', 3, lpara (1) , 9)&
            ) then                                                      
!                                                                       
!     ------Define whether this is a molecule or an object              
!                                                                       
               IF (str_comp (cpara (2) , 'atoms', 2, lpara (2) , 5) )   &
               then                                                     
                  mole_char (mole_num_mole) = MOLE_ATOM 
               ELSEIF (str_comp (cpara (2) , 'cube', 2, lpara (2) , 4) )&
               then                                                     
                  mole_char (mole_num_mole) = MOLE_CUBE 
               ELSEIF (str_comp (cpara (2) , 'cylinder', 2, lpara (2) , &
               8) ) then                                                
                  mole_char (mole_num_mole) = MOLE_CYLINDER 
               ELSEIF (str_comp (cpara (2) , 'sphere', 2, lpara (2) , 6)&
               ) then                                                   
                  mole_char (mole_num_mole) = MOLE_SPHERE 
               ELSEIF (str_comp (cpara (2) , 'edge', 2, lpara (2) , 4) )&
               then                                                     
                  mole_char (mole_num_mole) = MOLE_EDGE 
               ELSEIF (str_comp (cpara (2) , 'domain_cube', 9, lpara (2)&
               , 11) ) then                                             
                  mole_char (mole_num_mole) = MOLE_DOM_CUBE 
               ELSEIF (str_comp (cpara (2) , 'domain_cylinder', 9,      &
               lpara (2) , 15) ) then                                   
                  mole_char (mole_num_mole) = MOLE_DOM_CYLINDER 
               ELSEIF (str_comp (cpara (2) , 'domain_sphere', 9, lpara (&
               2) , 13) ) then                                          
                  mole_char (mole_num_mole) = MOLE_DOM_SPHERE 
               ELSEIF (str_comp (cpara (2) , 'domain_fuzzy', 9, lpara ( &
               2) , 12) ) then                                          
                  mole_char (mole_num_mole) = MOLE_DOM_FUZZY 
               ELSE 
                  ier_num = - 82 
                  ier_typ = ER_APPL 
               ENDIF 
!                                                                       
            ELSEIF (str_comp (cpara (1) , 'file', 3, lpara (1) , 4) )   &
            then                                                        
               mole_file (mole_num_mole) = cpara (2) 
!                                                                       
            ELSEIF (str_comp (cpara (1) , 'density', 3, lpara (1) , 6) )&
            then                                                        
!                                                                       
!     ------Define the scattering density of an object                  
!                                                                       
               cpara (1) = '0' 
               lpara (1) = 1 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               IF (ier_num.eq.0) then 
                  IF (ianz.eq.2) then 
                     mole_dens (mole_num_mole) = werte (2) 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ENDIF 
!                                                                       
            ELSEIF (str_comp (cpara (1) , 'fuzzy', 3, lpara (1) , 5) )  &
            then                                                        
!                                                                       
!     ------Define the minimum distance between atoms in a domain       
!             and the host                                              
!                                                                       
               cpara (1) = '0' 
               lpara (1) = 1 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               IF (ier_num.eq.0) then 
                  IF (ianz.eq.2) then 
                     mole_fuzzy (mole_num_mole) = werte (2) 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ENDIF 
!                                                                       
            ELSEIF (str_comp (cpara (1) , 'generator', 3, lpara (1) , 9)&
            ) then                                                      
!                                                                       
!     ------Define which generators create atoms within the             
!           same molecule                                               
!           Obsolete statement, is done automatically!!! RBN            
!                                                                       
               ier_num = + 2 
               ier_typ = ER_APPL 
!                                                                       
            ELSEIF (str_comp (cpara (1) , 'symmetry', 4, lpara (1) , 8) &
            ) then                                                      
!                                                                       
!     ------Define which symmetries  create atoms within the            
!           same molecule                                               
!           Obsolete statement, is done automatically!!! RBN            
!                                                                       
               ier_num = + 2 
               ier_typ = ER_APPL 
!                                                                       
            ELSEIF (str_comp (cpara (1) , 'type', 3, lpara (1) , 4) )   &
            then                                                        
!                                                                       
!     ------Define the molecule type, if less than current highest      
!     ------type number diminuish type number by one                    
!                                                                       
               cpara (1) = '0' 
               lpara (1) = 1 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               IF (ier_num.eq.0) then 
                  IF (ianz.eq.2) then 
                     IF (nint (werte (2) ) .lt.mole_num_type) then 
                        mole_num_type = mole_num_type-1 
                        mole_type (mole_num_mole) = nint (werte (2) ) 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ENDIF 
!                                                                       
            ELSEIF (str_comp (cpara (1) , 'content', 4, lpara (1) , 7) )&
            then                                                        
!                                                                       
!     ------start reading a molecule content                            
!                                                                       
               cpara (1) = '0' 
               lpara (1) = 1 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               IF (ier_num.eq.0) then 
                  IF (ianz.eq.2.or.ianz.eq.3) then 
                     IF (mole_num_mole.lt.MOLE_MAX_MOLE) then 
                        IF (werte (2) .lt.MOLE_MAX_TYPE) then 
                           IF (mole_l_on) then 
                              mole_type (mole_num_mole) = int (werte (2)&
                              )                                         
                              mole_num_type = max (mole_num_type-1, int &
                              (werte (2) ) )                            
                           ELSE 
                              mole_num_atom = mole_off (mole_num_mole)  &
                              + mole_len (mole_num_mole)                
                              mole_num_mole = mole_num_mole+1 
                              mole_num_curr = mole_num_mole 
                              mole_type (mole_num_mole) = int (werte (2)&
                              )                                         
                              mole_off (mole_num_mole) = mole_num_atom 
                              mole_len (mole_num_mole) = 0 
                              mole_num_type = max (mole_num_type, int ( &
                              werte (2) ) )                             
                              mole_gene_n = 0 
                              mole_symm_n = 0 
                              mole_l_on = .true. 
                           ENDIF 
                        ELSE 
                           ier_num = - 64 
                           ier_typ = ER_APPL 
      ier_msg (1)  = 'First characters of wrong line' 
                           ier_msg (2) = zeile (1:40) 
                           ier_msg (3) = zeile (41:80) 
                        ENDIF 
                     ELSE 
                        ier_num = - 65 
                        ier_typ = ER_APPL 
                        ier_msg (1) = 'First characters of wrong line' 
                        ier_msg (2) = zeile (1:40) 
                        ier_msg (3) = zeile (41:80) 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                     ier_msg (1) = 'First characters of wrong line' 
                     ier_msg (2) = zeile (1:40) 
                     ier_msg (3) = zeile (41:80) 
                  ENDIF 
               ELSE 
                  ier_msg (1) = 'First characters of wrong line' 
                  ier_msg (2) = zeile (1:40) 
                  ier_msg (3) = zeile (41:80) 
               ENDIF 
            ELSEIF (str_comp (cpara (1) , 'atoms', 4, lpara (1) , 5) )  &
            then                                                        
!                                                                       
!     ------read a molecule content                                     
!                                                                       
               IF (mole_l_on) then 
                  cpara (1) = '0' 
                  lpara (1) = 1 
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  IF (ier_num.eq.0) then 
                     DO j = 2, ianz 
                     mole_len (mole_num_mole) = mole_len (mole_num_mole)&
                     + 1                                                
                     mole_cont (mole_off (mole_num_mole) + mole_len (   &
                     mole_num_mole) ) = int (werte (j) )                
                     ENDDO 
                  ELSE 
                     ier_msg (1) = 'First characters of wrong line' 
                     ier_msg (2) = zeile (1:40) 
                     ier_msg (3) = zeile (41:80) 
                  ENDIF 
               ELSE 
                  ier_num = - 65 
                  ier_typ = ER_APPL 
               ENDIF 
            ELSE 
               ier_num = - 84 
               ier_typ = ER_APPL 
               ier_msg (1) = 'First characters of wrong line' 
               ier_msg (2) = zeile (1:40) 
               ier_msg (3) = zeile (41:80) 
            ENDIF 
         ENDIF 
      ENDIF 
!                                                                       
      END SUBROUTINE struc_mole_header              
!********************************************************************** 
      SUBROUTINE readstru (NMAX, MAXSCAT, strucfile, cr_name, cr_spcgr, &
      cr_a0, cr_win, cr_natoms, cr_nscat, cr_dw, cr_at_lis, cr_pos,     &
      cr_iscat, cr_prop, cr_dim, as_natoms, as_at_lis, as_dw, as_pos,   &
      as_iscat, as_prop, sav_ncell, sav_r_ncell, sav_ncatoms,           &
      spcgr_ianz, spcgr_para)                                           
!-                                                                      
!           this subroutine reads an old structur.                      
!+                                                                      
      IMPLICIT none 
!                                                                       
      INTEGER NMAX 
      INTEGER MAXSCAT 
!                                                                       
!     include'crystal3.inc' 
      include'errlist.inc' 
!
      INTEGER,                       INTENT(INOUT)  :: cr_natoms
      INTEGER, DIMENSION(1:NMAX),    INTENT(INOUT)  :: cr_iscat
      INTEGER, DIMENSION(1:NMAX),    INTENT(INOUT)  :: cr_prop
      REAL   , DIMENSION(1:3,1:NMAX),INTENT(INOUT)  :: cr_pos
!                                                                       
      INTEGER sav_ncell (3) 
      INTEGER sav_ncatoms 
      LOGICAL sav_r_ncell 
      LOGICAL lcell 
!                                                                       
      INTEGER ist 
      PARAMETER (ist = 7) 
!                                                                       
      CHARACTER ( * ) strucfile 
      CHARACTER(80) cr_name 
      CHARACTER(16) cr_spcgr 
      CHARACTER(4) cr_at_lis (0:MAXSCAT) 
      CHARACTER(4) as_at_lis (0:MAXSCAT) 
!                                                                       
      INTEGER cr_nscat 
      INTEGER as_natoms 
      INTEGER as_prop (MAXSCAT) 
      INTEGER as_iscat (MAXSCAT) 
!                                                                       
      INTEGER spcgr_ianz 
      INTEGER spcgr_para 
!                                                                       
      REAL cr_a0 (3) 
      REAL cr_win (3) 
      REAL cr_dw (0:MAXSCAT) 
      REAL cr_dim (3, 2) 
      REAL as_pos (3, MAXSCAT) 
      REAL as_dw (0:MAXSCAT) 
!                                                                       
      INTEGER i 
      LOGICAL lread 
!                                                                       
      cr_natoms = 0 
      lread = .true. 
      lcell = .false. 
      CALL oeffne (ist, strucfile, 'old', lread) 
      IF (ier_num.eq.0) then 
         DO i = 1, 3 
         cr_dim (i, 1) = 1.e10 
         cr_dim (i, 2) = - 1.e10 
         ENDDO 
!                                                                       
!     --Read header of structure file                                   
!                                                                       
         CALL stru_readheader (ist, NMAX, MAXSCAT, lcell, cr_name,      &
         cr_spcgr, cr_at_lis, cr_nscat, cr_dw, cr_a0, cr_win, sav_ncell,&
         sav_r_ncell, sav_ncatoms, spcgr_ianz, spcgr_para)              
!                                                                       
         IF (ier_num.eq.0) then 
!                                                                       
            CALL struc_read_atoms (NMAX, MAXSCAT, cr_natoms, cr_nscat,  &
            cr_dw, cr_at_lis, cr_pos, cr_iscat, cr_prop, cr_dim,        &
            as_natoms, as_at_lis, as_dw, as_pos, as_iscat, as_prop)     
         ENDIF 
      ENDIF 
!                                                                       
  999 CONTINUE 
      CLOSE (ist) 
      IF (ier_num.eq. - 49) then 
         WRITE (ier_msg (1), 3000) cr_natoms + 1 
 3000 FORMAT      ('At atom number = ',i8) 
      ENDIF 
!                                                                       
 2000 FORMAT    (a) 
 2010 FORMAT    (a16) 
      END SUBROUTINE readstru                       
!********************************************************************** 
      SUBROUTINE stru_readheader (ist, HD_NMAX, HD_MAXSCAT, lcell, cr_name,   &
      cr_spcgr, cr_at_lis, cr_nscat, cr_dw, cr_a0, cr_win, sav_ncell,   &
      sav_r_ncell, sav_ncatoms, spcgr_ianz, spcgr_para)                 
!-                                                                      
!     This subroutine reads the header of a structure file              
!+                                                                      
      USE gen_add_mod 
      USE sym_add_mod 
      IMPLICIT none 
!                                                                       
      INTEGER HD_NMAX 
      INTEGER HD_MAXSCAT 
!                                                                       
      include'errlist.inc' 
!                                                                       
      CHARACTER(80) cr_name 
      CHARACTER(16) cr_spcgr 
      CHARACTER(4) cr_at_lis (0:HD_MAXSCAT) 
!                                                                       
      INTEGER cr_nscat 
!                                                                       
      INTEGER spcgr_ianz 
      INTEGER spcgr_para 
!                                                                       
      REAL cr_a0 (3) 
      REAL cr_win (3) 
      REAL cr_dw (0:HD_MAXSCAT) 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 13) 
!                                                                       
      CHARACTER(1024) line, cpara (maxw) 
      CHARACTER(1024) zeile 
      CHARACTER(6) befehl 
      INTEGER ist, i, ll, j, islash 
      INTEGER ianz 
!DBG      integer             spcgr_ianz                                
      INTEGER lpara (maxw), lp 
      INTEGER lbef, indxb 
      INTEGER xx_nscat, xx_nadp 
      INTEGER sav_ncell (3) 
      INTEGER sav_ncatoms 
      LOGICAL sav_r_ncell 
      LOGICAL lcell 
      LOGICAL lend 
!DBG      real            spcgr_para                                    
      REAL werte (maxw) 
!                                                                       
      INTEGER len_str 
      LOGICAL str_comp 
!                                                                       
      xx_nscat = 0 
      xx_nadp = 0 
!                                                                       
      ier_num = - 46 
      ier_typ = ER_APPL 
      READ (ist, 2000, end = 999, err = 999) cr_name 
!                                                                       
!     This construction is needed as long a the old cell file format    
!     must be supported.                                                
!                                                                       
!                                                                       
!     The maximum number of significant characters depends on the       
!     length of the character constant befehl.                          
!                                                                       
      lbef = len (befehl) 
      befehl = '    ' 
      indxb = index (cr_name, ' ') 
      lbef = min (indxb - 1, lbef) 
      befehl = cr_name (1:lbef) 
      lbef = len_str (befehl) 
      befehl = cr_name (1:lbef) 
      IF (str_comp (befehl, 'title', 1, lbef, 5) ) then 
!                                                                       
!     Read new header                                                   
!                                                                       
!                                                                       
!     --remove title string from cr_name                                
!                                                                       
         ll = len (cr_name) 
         ll = len_str (cr_name) 
         line = ' ' 
         IF (0.lt.indxb.and.indxb + 1.le.ll) then 
            line (1:ll - indxb) = cr_name (indxb + 1:ll) 
         ELSE 
            line = ' ' 
         ENDIF 
         cr_name = line 
!                                                                       
         CALL no_error 
         DO while (.not.str_comp (befehl, 'atoms', 1, lbef, 5) ) 
         READ (ist, 2000, end = 999, err = 999) line 
         ll = 200 
         ll = len_str (line) 
!                                                                       
!     ----The maximum number of significant characters depends on the   
!     ----length of the character constant befehl.                      
!                                                                       
         lbef = len (befehl) 
      befehl = '    ' 
         indxb = index (line, ' ') 
         lbef = min (indxb - 1, lbef) 
         befehl = line (1:lbef) 
         lbef = len_str (befehl) 
         befehl = line (1:lbef) 
!                                                                       
!     ----command parameters start at the first character following     
!------ ----the blank                                                   
!                                                                       
         zeile = ' ' 
         lp = 0 
         IF (indxb + 1.le.ll) then 
            zeile = line (indxb + 1:ll) 
            lp = ll - indxb 
         ENDIF 
!                                                                       
!     ----Commentary                                                    
!                                                                       
                                                                        
         IF (line.eq.' '.or.line (1:1) .eq.'#'.or.line.eq.char (13) )   &
         then                                                           
            CONTINUE 
!                                                                       
!     ----Space group symbol                                            
!                                                                       
         ELSEIF (str_comp (befehl, 'spcgr', 1, lbef, 5) ) then 
            CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
            IF (ianz.lt.1) then 
               ier_num = - 100 
               ier_typ = ER_APPL 
               RETURN 
            ENDIF 
            islash = index (cpara (1) (1:lpara (1) ) , 'S') 
            DO while (islash.ne.0) 
            cpara (1) (islash:islash) = '/' 
            islash = index (cpara (1) (1:lpara (1) ) , 'S') 
            ENDDO 
            cr_spcgr = cpara (1) 
            spcgr_ianz = ianz - 1 
            ianz = ianz - 1 
            spcgr_para = 1 
            IF (ianz.eq.1) then 
               cpara (1) = cpara (2) 
               lpara (1) = lpara (2) 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               spcgr_para = nint (werte (1) ) 
            ENDIF 
            IF (ier_num.ne.0) then 
               ier_num = - 47 
               ier_typ = ER_APPL 
               RETURN 
            ENDIF 
!                                                                       
!     ----Cell constants                                                
!                                                                       
         ELSEIF (str_comp (befehl, 'cell', 1, lbef, 4) ) then 
            CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
            IF (ier_num.eq.0) then 
               IF (ianz.eq.6) then 
!     --------New style, kommata included                               
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  IF (ier_num.eq.0) then 
                     DO i = 1, 3 
                     cr_a0 (i) = werte (i) 
                     cr_win (i) = werte (i + 3) 
                     ENDDO 
                  ELSE 
                     ier_msg (1) = 'Error reading unit cell parameters' 
                     RETURN 
                  ENDIF 
               ELSE 
!     --------Old style, no kommata included                            
                  ier_num = - 48 
                  ier_typ = ER_APPL 
                  READ (zeile, *, end = 999, err = 999) (cr_a0 (i),     &
                  i = 1, 3), (cr_win (i), i = 1, 3)                     
                  CALL no_error 
               ENDIF 
            ELSE 
               ier_msg (1) = 'Error reading unit cell parameters' 
               RETURN 
            ENDIF 
            IF (cr_a0 (1) .le.0.0.or.cr_a0 (2) .le.0.0.or.cr_a0 (3)     &
            .le.0.0.or.cr_win (1) .le.0.0.or.cr_win (2)                 &
            .le.0.0.or.cr_win (3) .le.0.0.or.cr_win (1)                 &
            .ge.180.0.or.cr_win (2) .ge.180.0.or.cr_win (3) .ge.180.0)  &
            then                                                        
               ier_num = - 93 
               ier_typ = ER_APPL 
               ier_msg (1) = 'Error reading unit cell parameters' 
               RETURN 
            ENDIF 
!                                                                       
!     ----Additional symmetry generators 'generator'                    
!                                                                       
         ELSEIF (str_comp (befehl, 'gene', 1, lbef, 4) ) then 
            IF (gen_add_n.lt.GEN_ADD_MAX) then 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ier_num.eq.0) then 
                  IF (ianz.eq.12.or.ianz.eq.13) then 
                     CALL ber_params (ianz, cpara, lpara, werte,     &
                     maxw)                                           
                     IF (ier_num.eq.0) then 
                        gen_add_n = gen_add_n + 1 
                        DO j = 1, 4 
                        gen_add (1, j, gen_add_n) = werte (j) 
                        gen_add (2, j, gen_add_n) = werte (j + 4) 
                        gen_add (3, j, gen_add_n) = werte (j + 8) 
                        ENDDO 
                        IF (ianz.eq.13) then 
                           gen_add_power (gen_add_n) = nint (werte ( &
                           13) )                                     
                        ELSE 
                           gen_add_power (gen_add_n) = 1 
                        ENDIF 
                     ENDIF 
                  ELSEIF (ianz.eq.1) then 
                     lend = .true. 
                     READ (zeile, *, end = 8000) (werte (j), j = 1,  &
                     13)                                             
                     lend = .false. 
                     gen_add_n = gen_add_n + 1 
                     DO j = 1, 4 
                     gen_add (1, j, gen_add_n) = werte (j) 
                     gen_add (2, j, gen_add_n) = werte (j + 4) 
                     gen_add (3, j, gen_add_n) = werte (j + 8) 
                     ENDDO 
                     gen_add_power (gen_add_n) = nint (werte (13) ) 
 8000                CONTINUE 
                     IF (lend) then 
                        ier_num = - 92 
                        ier_typ = ER_APPL 
                        RETURN 
                     ENDIF 
                  ELSE 
                     ier_num = - 92 
                     ier_typ = ER_APPL 
                  ENDIF 
               ENDIF 
            ELSE 
               ier_num = - 61 
               ier_typ = ER_APPL 
               RETURN 
            ENDIF 
!                                                                       
!     ----Names of atoms to setup specific sequence of scattering curves
!                                                               'scat'  
!                                                                       
         ELSEIF (str_comp (befehl, 'scat', 2, lbef, 4) ) then 
            CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
            IF (ier_num.eq.0) then 
               IF (xx_nscat + ianz.le.HD_MAXSCAT) then 
                  DO i = 1, ianz 
                  CALL do_cap (cpara (i) (1:lpara (i) ) ) 
                  cr_at_lis (xx_nscat + i) = cpara (i) 
                  ENDDO 
                  xx_nscat = xx_nscat + ianz 
                  cr_nscat = max (cr_nscat, xx_nscat) 
               ELSE 
                  ier_num = -26 
                  ier_typ = ER_APPL 
               ENDIF 
            ELSE 
               ier_num = -111
               ier_typ = ER_COMM 
            ENDIF 
!                                                                       
!     ----Displacement parameters to setup specific sequence of         
!                                    scattering curves 'adp'            
!                                                                       
         ELSEIF (str_comp (befehl, 'adp', 2, lbef, 3) ) then 
            CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
            IF (ier_num.eq.0) then 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               IF (ier_num.eq.0) then 
                  IF (xx_nadp + ianz.le.HD_MAXSCAT) then 
                     DO i = 1, ianz 
                     cr_dw (xx_nadp + i) = werte (i) 
                     ENDDO 
                     xx_nadp = xx_nadp + ianz 
                     cr_nscat = max (cr_nscat, xx_nadp) 
                  ELSE 
                     ier_num = - 26 
                     ier_typ = ER_APPL 
                  ENDIF 
               ENDIF 
            ELSE 
               ier_num = -112
               ier_typ = ER_COMM 
            ENDIF 
!                                                                       
!     ----Crystal dimensions and number of atoms per unit cell 'ncell'  
!                                                                       
         ELSEIF (str_comp (befehl, 'ncell', 1, lbef, 5) ) then 
            CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
            IF (ier_num.eq.0) then 
               IF (ianz.eq.4) then 
                  CALL ber_params (ianz, cpara, lpara, werte, maxw) 
                  IF (ier_num.eq.0) then 
                     DO j = 1, 3 
                     sav_ncell (j) = werte (j) 
                     ENDDO 
                     sav_ncatoms = werte (4) 
                     sav_r_ncell = .true. 
                  ENDIF 
               ELSE 
                  ier_num = - 6 
                  ier_typ = ER_COMM 
               ENDIF 
            ENDIF 
!                                                                       
!     ----Additional symmetry operations 'symmetry'                     
!                                                                       
         ELSEIF (str_comp (befehl, 'symm', 2, lbef, 4) ) then 
            IF (sym_add_n.lt.SYM_ADD_MAX) then 
               CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
               IF (ier_num.eq.0) then 
                  IF (ianz.eq.12.or.ianz.eq.13) then 
                     CALL ber_params (ianz, cpara, lpara, werte,     &
                     maxw)                                           
                     IF (ier_num.eq.0) then 
                        sym_add_n = sym_add_n + 1 
                        DO j = 1, 4 
                        sym_add (1, j, sym_add_n) = werte (j) 
                        sym_add (2, j, sym_add_n) = werte (j + 4) 
                        sym_add (3, j, sym_add_n) = werte (j + 8) 
                        ENDDO 
                        IF (ianz.eq.13) then 
                           sym_add_power (sym_add_n) = nint (werte ( &
                           13) )                                     
                        ELSE 
                           sym_add_power (sym_add_n) = 1 
                        ENDIF 
                     ENDIF 
                  ELSE 
                     ier_num = - 6 
                     ier_typ = ER_COMM 
                  ENDIF 
               ENDIF 
            ELSE 
               ier_num = - 62 
               ier_typ = ER_APPL 
               RETURN 
            ENDIF 
!                                                                       
         ELSEIF (str_comp (befehl, 'atoms', 2, lbef, 5) ) then 
            CONTINUE 
         ELSE 
            ier_num = - 89 
            ier_typ = ER_APPL 
         ENDIF 
         IF (ier_num.ne.0) then 
            RETURN 
         ENDIF 
         ENDDO 
!                                                                       
         CALL no_error 
!                                                                       
!     --Determine space group number                                    
!                                                                       
         werte (1) = spcgr_para 
         CALL spcgr_no (spcgr_ianz, maxw, werte) 
      ELSE 
!                                                                       
!     Read old header                                                   
!                                                                       
         ier_num = - 47 
         ier_typ = ER_APPL 
         READ (ist, 2010, end = 999, err = 999) line 
         lp = len_str (line) 
         CALL get_params (line, ianz, cpara, lpara, maxw, lp) 
         cr_spcgr = cpara (1) 
         ianz = ianz - 1 
         IF (ianz.eq.1) then 
            cpara (1) = cpara (2) 
            lpara (1) = lpara (2) 
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
         ENDIF 
         IF (ier_num.eq.0) then 
            ier_num = - 48 
            ier_typ = ER_APPL 
            READ (ist, *, end = 999, err = 999) (cr_a0 (i), i = 1, 3),  &
            (cr_win (i), i = 1, 3)                                      
            CALL no_error 
            CALL spcgr_no (ianz, maxw, werte) 
         ENDIF 
      ENDIF 
!                                                                       
  999 CONTINUE 
!
!                                                                       
 2000 FORMAT    (a) 
 2010 FORMAT    (a16) 
      END SUBROUTINE stru_readheader                
!********************************************************************** 
      SUBROUTINE struc_read_atoms (NMAX, MAXSCAT, cr_natoms, cr_nscat,  &
      cr_dw, cr_at_lis, cr_pos, cr_iscat, cr_prop, cr_dim, as_natoms,   &
      as_at_lis, as_dw, as_pos, as_iscat, as_prop)                      
!-                                                                      
!           This subroutine reads the list of atoms into the            
!       crystal array                                                   
!+                                                                      
      USE molecule_mod 
      IMPLICIT none 
!                                                                       
      INTEGER NMAX 
      INTEGER MAXSCAT 
!                                                                       
!     include'crystal3.inc' 
      include'errlist.inc' 
!
      INTEGER,                       INTENT(INOUT)  :: cr_natoms
      INTEGER, DIMENSION(1:NMAX),    INTENT(INOUT)  :: cr_iscat
      INTEGER, DIMENSION(1:NMAX),    INTENT(INOUT)  :: cr_prop
      REAL   , DIMENSION(1:3,1:NMAX),INTENT(INOUT)  :: cr_pos
!                                                                       
      INTEGER ist, maxw 
      PARAMETER (ist = 7, maxw = 5) 
!                                                                       
      CHARACTER(4) cr_at_lis (0:MAXSCAT) 
      CHARACTER(4) as_at_lis (0:MAXSCAT) 
!                                                                       
      INTEGER cr_nscat 
      INTEGER as_natoms 
      INTEGER as_iscat (MAXSCAT) 
      INTEGER as_prop (MAXSCAT) 
!                                                                       
      REAL cr_a0 (3) 
      REAL cr_win (3) 
      REAL cr_dim (3, 2) 
      REAL cr_dw (0:MAXSCAT) 
      REAL as_pos (3, MAXSCAT) 
      REAL as_dw (0:MAXSCAT) 
!                                                                       
      CHARACTER(10) befehl 
      CHARACTER(80) cpara (maxw) 
      CHARACTER(1024) line, zeile 
      INTEGER lpara (maxw) 
      INTEGER ianz 
      INTEGER i, j, ibl, lbef 
      INTEGER lline 
      REAL werte (maxw), dw1 
!                                                                       
      INTEGER len_str 
      LOGICAL str_comp 
!                                                                       
 1000 CONTINUE 
      ier_num = - 49 
      ier_typ = ER_APPL 
      line = ' ' 
      READ (ist, 2000, end = 2, err = 999) line 
      lline = len_str (line) 
      IF (line.ne.' '.and.line (1:1) .ne.'#'.and.line.ne.char (13) )    &
      then                                                              
         ibl = index (line (1:lline) , ' ') + 1 
         lbef = 10 
         befehl = ' ' 
         ibl = index (line (1:lline) , ' ') 
         IF (ibl.eq.0) then 
            ibl = lline+1 
         ENDIF 
         lbef = min (ibl - 1, lbef) 
         befehl = line (1:lbef) 
         IF (str_comp (befehl, 'molecule', 4, lbef, 8) .or.str_comp (   &
         befehl, 'domain', 4, lbef, 6) .or.str_comp (befehl, 'object',  &
         4, lbef, 6) ) then                                             
!                                                                       
!     ------Start/End of a molecule                                     
!                                                                       
            CALL no_error 
            IF (ibl.le.lline) then 
               zeile = line (ibl:lline) 
               i = lline-ibl + 1 
            ELSE 
               zeile = ' ' 
               i = 0 
            ENDIF 
            CALL struc_mole_header (zeile, i, .false.) 
            IF (ier_num.ne.0) return 
         ELSE 
            CALL read_atom_line (line, ibl, lline, as_natoms, maxw,     &
            werte)                                                      
            IF (ier_num.ne.0.and.ier_num.ne. - 49) then 
               GOTO 999 
            ENDIF 
            IF (cr_natoms.eq.nmax) then 
!                                                                       
!     --------Too many atoms in the structure file                      
!                                                                       
               ier_num = - 10 
               ier_typ = ER_APPL 
               RETURN 
            ENDIF 
            cr_natoms = cr_natoms + 1 
            i = cr_natoms 
            DO j = 1, 3 
            cr_pos (j, i) = werte (j) 
            cr_dim (j, 1) = amin1 (cr_dim (j, 1), cr_pos (j, i) ) 
            cr_dim (j, 2) = amax1 (cr_dim (j, 2), cr_pos (j, i) ) 
            ENDDO 
            dw1 = werte (4) 
            cr_prop (i) = nint (werte (5) ) 
      IF (line (1:4) .ne.'    ') then 
               ibl = ibl - 1 
               CALL do_cap (line (1:ibl) ) 
               DO j = 0, cr_nscat 
               IF (line (1:ibl) .eq.cr_at_lis (j) .and.dw1.eq.cr_dw (j) &
               ) then                                                   
                  cr_iscat (i) = j 
                  GOTO 11 
               ENDIF 
               ENDDO 
               IF (cr_nscat.eq.MAXSCAT) then 
!                                                                       
!     --------  Too many atom types in the structure file               
!                                                                       
                  ier_num = -72 
                  ier_typ = ER_APPL 
                  RETURN 
               ENDIF 
               cr_nscat = cr_nscat + 1 
               cr_iscat (i) = cr_nscat 
               cr_at_lis (cr_nscat) = line (1:ibl) 
               cr_dw (cr_nscat) = dw1 
!                                                                       
               IF (0.0.le.cr_pos (1, i) .and.cr_pos (1, i)              &
               .lt.1.and.0.0.le.cr_pos (2, i) .and.cr_pos (2, i)        &
               .lt.1.and.0.0.le.cr_pos (3, i) .and.cr_pos (3, i) .lt.1) &
               then                                                     
                  as_natoms = as_natoms + 1 
                  as_at_lis (cr_nscat) = cr_at_lis (cr_nscat) 
                  as_iscat (as_natoms) = cr_iscat (i) 
                  as_dw (as_natoms) = cr_dw (cr_nscat) 
                  DO j = 1, 3 
                  as_pos (j, as_natoms) = cr_pos (j, i) 
                  ENDDO 
                  as_prop (as_natoms) = cr_prop (i) 
               ENDIF 
   11          CONTINUE 
!                                                                       
!     --------If we are reading a molecule insert atom into current     
!                                                                       
               IF (mole_l_on) then 
                  CALL mole_insert_current (cr_natoms, mole_num_curr) 
                  IF (ier_num.lt.0.and.ier_num.ne. - 49) then 
                     GOTO 999 
                  ENDIF 
               ENDIF 
            ENDIF 
         ENDIF 
      ENDIF 
      GOTO 1000 
!                                                                       
    2 CONTINUE 
      IF (ier_num.eq. - 49) then 
         CALL no_error 
      ENDIF 
!
!                                                                       
  999 CONTINUE 
      CLOSE (ist) 
!                                                                       
 2000 FORMAT    (a) 
      END SUBROUTINE struc_read_atoms               
!********************************************************************** 
!********************************************************************** 
      SUBROUTINE spcgr_no (ianz, maxw, werte) 
!-                                                                      
!     Interprets the space group symbol. Returns the space group no.    
!+                                                                      
      USE config_mod 
      USE crystal_mod 
      USE spcgr_mod 
      IMPLICIT none 
!                                                                       
       
      include'errlist.inc' 
!                                                                       
      CHARACTER(1024) cpara 
      INTEGER lpara 
      INTEGER ianz, maxw, ii, i 
      REAL werte (maxw) 
      REAL rpara 
!                                                                       
      INTEGER len_str 
!                                                                       
      CALL no_error 
      ii = 1 
      IF (ianz.eq.1) then 
         ii = nint (werte (1) ) 
      ENDIF 
!                                                                       
!     Distinguish between trigonal and rhombohedral settings            
!                                                                       
      IF (cr_spcgr (1:1) .eq.'R') then 
         IF (cr_a0 (1) .eq.cr_a0 (2) .and.cr_win (1) .eq.90..and.cr_win &
         (2) .eq.90.0.and.cr_win (3) .eq.120.0) then                    
            ii = 1 
         ELSEIF (cr_a0 (1) .eq.cr_a0 (2) .and.cr_a0 (1) .eq.cr_a0 (3)   &
         .and.cr_win (1) .eq.cr_win (2) .and.cr_win (1) .eq.cr_win (3) )&
         then                                                           
            ii = 2 
         ELSE 
            ier_num = - 7 
            ier_typ = ER_APPL 
         ENDIF 
      ENDIF 
!                                                                       
      IF (ii.lt.1.or.2.lt.ii) then 
         ier_num = - 7 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
!                                                                       
      cr_spcgrno = 0 
      cr_syst = 0 
      ier_num = - 7 
      ier_typ = ER_APPL 
!                                                                       
      DO i = 1, SPCGR_MAX 
      IF (cr_spcgr.eq.spcgr_name (i) ) then 
         cr_spcgrno = spcgr_num (i, ii) 
         cr_syst = spcgr_syst (cr_spcgrno) 
         CALL no_error 
         GOTO 10 
      ENDIF 
      ENDDO 
!                                                                       
      CALL no_error 
      cpara = cr_spcgr 
      lpara = len_str (cpara) 
      CALL ber_params (1, cpara, lpara, rpara, 1) 
      IF (ier_num.eq.0) then 
         cr_spcgrno = nint (rpara) 
         cr_syst = spcgr_syst (cr_spcgrno) 
         cr_spcgr = spcgr_name (cr_spcgrno) 
      ELSE 
         ier_num = - 7 
         ier_typ = ER_APPL 
      ENDIF 
!                                                                       
   10 CONTINUE 
!                                                                       
      IF (ier_num.eq.0) then 
         ier_num = - 14 
         ier_typ = ER_APPL 
         IF (cr_syst.eq.1) then 
            CALL no_error 
         ELSEIF (cr_syst.eq.2) then 
            IF (cr_win (1) .eq.90.0.and.cr_win (3) .eq.90.0) then 
               CALL no_error 
            ENDIF 
         ELSEIF (cr_syst.eq.3) then 
            IF (cr_win (1) .eq.90.0.and.cr_win (2) .eq.90.0) then 
               CALL no_error 
            ENDIF 
         ELSEIF (cr_syst.eq.4) then 
            IF (cr_win (1) .eq.90.0.and.cr_win (2) .eq.90.0.and.cr_win (&
            3) .eq.90.0) then                                           
               CALL no_error 
            ENDIF 
         ELSEIF (cr_syst.eq.5) then 
            IF (cr_a0 (1) .eq.cr_a0 (2) .and.cr_win (1)                 &
            .eq.90.0.and.cr_win (2) .eq.90.0.and.cr_win (3) .eq.90.0)   &
            then                                                        
               CALL no_error 
            ENDIF 
         ELSEIF (cr_syst.eq.6.or.cr_syst.eq.8) then 
            IF (cr_a0 (1) .eq.cr_a0 (2) .and.cr_win (1)                 &
            .eq.90.0.and.cr_win (2) .eq.90.0.and.cr_win (3) .eq.120.0)  &
            then                                                        
               CALL no_error 
            ENDIF 
         ELSEIF (cr_syst.eq.7) then 
            IF (cr_a0 (1) .eq.cr_a0 (2) .and.cr_a0 (1) .eq.cr_a0 (3)    &
            .and.cr_win (1) .eq.cr_win (2) .and.cr_win (2) .eq.cr_win ( &
            3) ) then                                                   
               CALL no_error 
            ENDIF 
         ELSEIF (cr_syst.eq.9) then 
            IF (cr_a0 (1) .eq.cr_a0 (2) .and.cr_a0 (1) .eq.cr_a0 (3)    &
            .and.cr_win (1) .eq.90.0.and.cr_win (2) .eq.90.0.and.cr_win &
            (3) .eq.90.0) then                                          
               CALL no_error 
            ENDIF 
         ENDIF 
      ENDIF 
      END SUBROUTINE spcgr_no                       
!********************************************************************** 
!*****7*****************************************************************
      SUBROUTINE get_symmetry_matrices 
!-                                                                      
!     Creates all symmetry matrices for the current space group         
!+                                                                      
      USE config_mod 
      USE crystal_mod 
      USE generate_mod 
      USE gen_add_mod 
      USE sym_add_mod 
      USE unitcell_mod 
      USE wyckoff_mod 
      IMPLICIT none 
!                                                                       
       
      include'errlist.inc' 
!                                                                       
!                                                                       
      INTEGER igs 
      INTEGER igg 
      INTEGER i, j, k 
!                                                                       
!     Reset all symmetry matrices                                       
!                                                                       
      DO k = 1, SPC_MAX 
      DO i = 1, 4 
      DO j = 1, 4 
      spc_mat (i, j, k) = 0.0 
      ENDDO 
      ENDDO 
      ENDDO 
!                                                                       
!     The first symmetry matrix is the identity matrix                  
!                                                                       
      spc_n = 1 
!                                                                       
      DO i = 1, 4 
      DO j = 1, 4 
      spc_mat (i, j, spc_n) = 0.0 
      ENDDO 
      spc_mat (i, i, spc_n) = 1.0 
      ENDDO 
      spc_spur (spc_n) = 3 
      spc_det (spc_n) = 1 
      CALL get_symmetry_type (SPC_MAX, spc_n, spc_mat, spc_spur,        &
      spc_det, spc_char, spc_xyz, cr_syst)                              
!                                                                       
!     Loop over all generators                                          
!                                                                       
      DO igs = 1, generspcgr (0, cr_spcgrno) 
      IF (gen_sta.eq.GEN_SYMM) then 
         igg = generspcgr (igs, cr_spcgrno) 
      ELSEIF (gen_sta.eq.GEN_CENTER) then 
         igg = generspcgr_center (igs, cr_spcgrno) 
      ENDIF 
      CALL make_symmetry_matrix (SPC_MAX, spc_n, spc_mat, spc_det,      &
      spc_spur, spc_char, spc_xyz, igg, NG, generators, generpower)     
      IF (ier_num.ne.0) then 
         RETURN 
      ENDIF 
      ENDDO 
!-----      End of loop over all generators                             
!                                                                       
!     Loop over all additional generators                               
!                                                                       
      DO igs = 1, gen_add_n 
      CALL make_symmetry_matrix (SPC_MAX, spc_n, spc_mat, spc_det,      &
      spc_spur, spc_char, spc_xyz, igs, GEN_ADD_MAX, gen_add,           &
      gen_add_power)                                                    
      IF (ier_num.ne.0) then 
         RETURN 
      ENDIF 
      ENDDO 
!-----      End of loop over all additional generators                  
!                                                                       
!     Loop over all additional symmetry operations, just copy these     
!                                                                       
      DO igs = 1, sym_add_n 
      spc_n = spc_n + 1 
      spc_det (spc_n) = 0 
      spc_spur (spc_n) = 0 
      DO i = 1, 4 
      DO j = 1, 4 
      DO k = 1, 4 
      spc_mat (i, j, spc_n) = sym_add (i, j, igs) 
      ENDDO 
      ENDDO 
      ENDDO 
      DO i = 1, 3 
      IF (spc_mat (i, 4, spc_n) .ge.1.0) then 
         spc_mat (i, 4, spc_n) = spc_mat (i, 4, spc_n) - int (spc_mat ( &
         i, 4, spc_n) )                                                 
      ELSEIF (spc_mat (i, 4, spc_n) .lt.0.0) then 
         spc_mat (i, 4, spc_n) = spc_mat (i, 4, spc_n) + int (spc_mat ( &
         i, 4, spc_n) ) + 1                                             
      ENDIF 
      spc_spur (spc_n) = spc_spur (spc_n) + spc_mat (i, i, spc_n) 
      ENDDO 
      spc_det (spc_n) = spc_mat (1, 1, spc_n) * (spc_mat (2, 2, spc_n)  &
      * spc_mat (3, 3, spc_n) - spc_mat (3, 2, spc_n) * spc_mat (2, 3,  &
      spc_n) ) - spc_mat (1, 2, spc_n) * (spc_mat (2, 1, spc_n) *       &
      spc_mat (3, 3, spc_n) - spc_mat (3, 1, spc_n) * spc_mat (2, 3,    &
      spc_n) ) + spc_mat (1, 3, spc_n) * (spc_mat (2, 1, spc_n) *       &
      spc_mat (3, 2, spc_n) - spc_mat (3, 1, spc_n) * spc_mat (2, 2,    &
      spc_n) )                                                          
      CALL get_symmetry_type (SPC_MAX, spc_n, spc_mat, spc_spur,        &
      spc_det, spc_char, spc_xyz, cr_syst)                              
      ENDDO 
!-----      End of loop over additional symmetry matices                
!                                                                       
      END SUBROUTINE get_symmetry_matrices          
!*****7*****************+***********************************************
      SUBROUTINE make_symmetry_matrix (SPC_MAX, spc_n, spc_mat, spc_det,&
      spc_spur, spc_char, spc_xyz, igg, NG, generators, generpower)     
!-                                                                      
!     Applies the current generator to all existing symmetry matrices   
!+                                                                      
      USE config_mod 
      USE crystal_mod 
      IMPLICIT none 
!                                                                       
       
      include'errlist.inc' 
!                                                                       
      INTEGER SPC_MAX 
      CHARACTER(65) spc_char (1:SPC_MAX) 
      CHARACTER(87) spc_xyz (1:SPC_MAX) 
      INTEGER spc_n 
      REAL spc_mat (4, 4, 1:SPC_MAX) 
      REAL spc_det (1:SPC_MAX) 
      REAL spc_spur (1:SPC_MAX) 
      INTEGER igg 
      INTEGER NG 
      REAL generators (4, 4, 0:NG) 
      INTEGER generpower (NG) 
!                                                                       
      INTEGER i, j, k 
      INTEGER ipg 
      INTEGER imat 
      INTEGER nmat 
      REAL xmat (4, 4) 
      REAL wmat (4, 4) 
!                                                                       
!     For convenience create Identity operator                          
!                                                                       
      DO i = 1, 4 
      DO j = 1, 4 
      xmat (i, j) = 0.0 
      ENDDO 
      xmat (i, i) = 1.0 
      ENDDO 
!                                                                       
!     nmat is the number of symetry operations prior to this generator  
!     if the generator is applied not only to the power of 1, but 2 as  
!     well, it must act on the number of matrices that existed prior    
!     to ist application only.                                          
!                                                                       
      nmat = spc_n 
!                                                                       
!     Loop over all powers of generator igg                             
!                                                                       
      DO ipg = 1, generpower (igg) 
!                                                                       
!     --raise power of generator, xmat is a dummy matrix, equal to      
!     --the previous power of the generator                             
!                                                                       
      DO i = 1, 4 
      DO j = 1, 4 
      wmat (i, j) = 0.0 
      DO k = 1, 4 
      wmat (i, j) = wmat (i, j) + generators (i, k, igg) * xmat (k, j) 
      ENDDO 
      ENDDO 
      ENDDO 
!                                                                       
!     --apply generator to all previous symmetry matrices               
!                                                                       
      DO imat = 1, nmat 
      spc_n = spc_n + 1 
      spc_det (spc_n) = 0 
      spc_spur (spc_n) = 0 
      DO i = 1, 4 
      DO j = 1, 4 
      DO k = 1, 4 
      spc_mat (i, j, spc_n) = spc_mat (i, j, spc_n) + wmat (i, k)       &
      * spc_mat (k, j, imat)                                            
      ENDDO 
      ENDDO 
      ENDDO 
      DO i = 1, 3 
      IF (spc_mat (i, 4, spc_n) .ge.1.0) then 
         spc_mat (i, 4, spc_n) = spc_mat (i, 4, spc_n) - int (spc_mat ( &
         i, 4, spc_n) )                                                 
      ELSEIF (spc_mat (i, 4, spc_n) .lt.0.0) then 
         spc_mat (i, 4, spc_n) = spc_mat (i, 4, spc_n) + int (spc_mat ( &
         i, 4, spc_n) ) + 1                                             
      ENDIF 
      spc_spur (spc_n) = spc_spur (spc_n) + spc_mat (i, i, spc_n) 
      ENDDO 
      spc_det (spc_n) = spc_mat (1, 1, spc_n) * (spc_mat (2, 2, spc_n)  &
      * spc_mat (3, 3, spc_n) - spc_mat (3, 2, spc_n) * spc_mat (2, 3,  &
      spc_n) ) - spc_mat (1, 2, spc_n) * (spc_mat (2, 1, spc_n) *       &
      spc_mat (3, 3, spc_n) - spc_mat (3, 1, spc_n) * spc_mat (2, 3,    &
      spc_n) ) + spc_mat (1, 3, spc_n) * (spc_mat (2, 1, spc_n) *       &
      spc_mat (3, 2, spc_n) - spc_mat (3, 1, spc_n) * spc_mat (2, 2,    &
      spc_n) )                                                          
      CALL get_symmetry_type (SPC_MAX, spc_n, spc_mat, spc_spur,        &
      spc_det, spc_char, spc_xyz, cr_syst)                              
      ENDDO 
!                                                                       
!     --Set power of Generator                                          
!                                                                       
      DO i = 1, 4 
      DO j = 1, 4 
      xmat (i, j) = wmat (i, j) 
      ENDDO 
      ENDDO 
      ENDDO 
!                                                                       
      END SUBROUTINE make_symmetry_matrix           
!********************************************************************** 
      SUBROUTINE get_symmetry_type (SPC_MAX, spc_n, spc_mat, spc_spur,  &
      spc_det, spc_char, spc_xyz, cr_syst)                              
!-                                                                      
!     Determines the xyz triplet, and the letter that describes the     
!     symmetry operation                                                
!+                                                                      
      IMPLICIT none 
!                                                                       
      INTEGER SPC_MAX 
      CHARACTER(65) spc_char (1:SPC_MAX) 
      CHARACTER(87) spc_xyz (1:SPC_MAX) 
      INTEGER spc_n 
      REAL spc_mat (4, 4, 1:SPC_MAX) 
      REAL spc_det (1:SPC_MAX) 
      REAL spc_spur (1:SPC_MAX) 
      INTEGER cr_syst 
!                                                                       
      INTEGER i, j, ii, ja, je 
      INTEGER pwr ( - 3:3, - 1:1) 
      INTEGER power 
      CHARACTER(2) typ ( - 3:3, - 1:1) 
      CHARACTER(LEN=1) :: xyz (3) = (/'x','y','z'/)
      CHARACTER(LEN=6) :: vec (-12:12) 
      REAL work (3, 3) 
      REAL add (3) 
      REAL axis (3) 
      REAL posit (3) 
      REAL screw (3) 
      REAL hkl (3) 
      INTEGER           :: mult
      REAL              :: fact
!                                                                       
      DO i = - 3, 3 
         DO j = - 1, 1 
            typ (i, j) = ' ' 
            pwr (i, j) = 0 
         ENDDO 
      ENDDO 
      typ ( 3, 1) = ' 1' 
      typ (-1, 1) = ' 2' 
      typ ( 0, 1) = ' 3' 
      typ ( 1, 1) = ' 4' 
      typ ( 2, 1) = ' 6' 
      typ (-3,-1) = '-1' 
      typ ( 1,-1) = ' m' 
      typ ( 0,-1) = '-3' 
      typ (-1,-1) = '-4' 
      typ (-2,-1) = '-6' 
      pwr ( 3, 1) = 1 
      pwr (-1, 1) = 2 
      pwr ( 0, 1) = 3 
      pwr ( 1, 1) = 4 
      pwr ( 2, 1) = 6 
      pwr (-3,-1) = 1 
      pwr ( 1,-1) = 2 
      pwr ( 0,-1) = 3 
      pwr (-1,-1) = 4 
      pwr (-2,-1) = 6 
!                                                                       
      vec(-12) = '  -1/1' 
      vec(-11) = '-11/12' 
      vec(-10) = '  -5/6' 
      vec( -9) = '  -3/4' 
      vec( -8) = '  -2/3' 
      vec( -7) = ' -7/12' 
      vec( -6) = '  -1/2' 
      vec( -5) = ' -5/12' 
      vec( -4) = '  -1/3' 
      vec( -3) = '  -1/4' 
      vec( -2) = '  -1/6' 
      vec( -1) = ' -1/12' 
      vec ( 0) = '      ' 
      vec ( 1) = ' +1/12' 
      vec ( 2) = '  +1/6' 
      vec ( 3) = '  +1/4' 
      vec ( 4) = '  +1/3' 
      vec ( 5) = ' +5/12' 
      vec ( 6) = '  +1/2' 
      vec ( 7) = ' +7/12' 
      vec ( 8) = '  +2/3' 
      vec ( 9) = '  +3/4' 
      vec (10) = '  +5/6' 
      vec (11) = '+11/12' 
      vec (12) = '  +1/1' 
!                                                                       
      spc_char (spc_n) = ' ' 
      spc_xyz (spc_n) = ' ' 
!                                                                       
      spc_char (spc_n) = typ (NINT(spc_spur (spc_n)), NINT(spc_det (spc_n)) ) 
      power = pwr (NINT(spc_spur (spc_n)), NINT(spc_det (spc_n)) ) 
!                                                                       
      DO i = 1, 3 
         ii = (i - 1) * 29 
         DO j = 1, 3 
         ja = ii + (j - 1) * 7 + 1 
         je = ii + (j - 1) * 7 + 7 
         mult = nint (12.*spc_mat (i, j, spc_n) )
         fact =      (12.*spc_mat (i, j, spc_n) )
         IF( abs(mult-fact) < 0.001 ) THEN
            IF    ( mult == -12 ) THEN
               WRITE(spc_xyz (spc_n) (ja:je), 4000) xyz(j)
            ELSEIF( mult ==   0 ) THEN
               WRITE(spc_xyz (spc_n) (ja:je), 4100) 
            ELSEIF( mult ==  12 ) THEN
               WRITE(spc_xyz (spc_n) (ja:je), 4200) xyz(j)
            ELSE
               WRITE(spc_xyz (spc_n) (ja:je), 4300) mult,xyz(j)
            ENDIF
         ELSE
            WRITE(spc_xyz (spc_n) (ja:je), 5100) fact,xyz(j)
         ENDIF
         ENDDO 
         ja = ii + 22
         je = ii + 29 
         mult = nint (12.*spc_mat (i, 4, spc_n) )
         fact =      (12.*spc_mat (i, 4, spc_n) )
         IF( abs(mult-fact) < 0.001 ) THEN
            IF(mult < -12) THEN
               WRITE(spc_xyz (spc_n) (ja:je), 6000) mult
            ELSEIF(mult >  12) THEN
               WRITE(spc_xyz (spc_n) (ja:je), 6000) mult
            ELSE
               WRITE(spc_xyz (spc_n) (ja:je), 6200) vec(mult)
            ENDIF
         ELSE
            WRITE(spc_xyz (spc_n) (ja:je), 7100) fact
         ENDIF
      ENDDO 
      spc_xyz (spc_n) (86:87) = ' ' 
      i = 87
      call rem_bl(spc_xyz(spc_n),i)
!                                                                       
      DO i = 1, 3 
      DO j = 1, 3 
      work (i, j) = spc_mat (i, j, spc_n) 
      ENDDO 
      add (i) = spc_mat (i, 4, spc_n) 
      ENDDO 
      CALL get_detail (spc_n, work, add, spc_char (spc_n), power, axis, &
      screw, posit, hkl)                                                
!                                                                       
4000  FORMAT(   '     -',A1)
4100  FORMAT(   '       ')
4200  FORMAT(   '     +',A1)
4300  FORMAT(SP,I3,'/12',A1)
5100  FORMAT(SP,F6.3    ,A1)
6000  FORMAT(SP,I3,'/12',', ')
6200  FORMAT(   A6,      ', ')
7100  FORMAT(SP,F6.3    ,', ')
!
      END SUBROUTINE get_symmetry_type              
!****&******************************************************************
      SUBROUTINE get_detail (no, work, add, w_char, power, axis, screw, &
      posit, hkl)                                                       
!-                                                                      
!     Determines the local symmetry of the position given in the line   
!+                                                                      
      USE config_mod 
      USE crystal_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER no 
      REAL work (3, 3) 
      REAL add (3) 
      CHARACTER(65) w_char 
      INTEGER power 
      REAL axis (3) 
      REAL screw (3) 
      REAL posit (3) 
      REAL posit1bar (3) 
      REAL hkl (3) 
!                                                                       
      CHARACTER(1) abc (0:3) 
      CHARACTER(3) xyz (3, - 2:2) 
      CHARACTER(3) xxx (3, - 2:2) 
      CHARACTER(3) Oyy (3, - 2:2) 
      CHARACTER(6) ctrans ( - 24:24) 
      INTEGER i, j, k, l 
      INTEGER ia, ie 
      INTEGER ii, jj, kk 
      INTEGER nnull 
      INTEGER nnullg 
      INTEGER iglide 
      LOGICAL lsearch 
      REAL cp (3, 3) 
      REAL temp (3, 3) 
      REAL t2 (3, 3) 
      REAL sum (3, 3) 
      REAL imat (3, 3) 
      REAL vector (3) 
      REAL p3 (3) 
      REAL hmin 
      REAL det 
      REAL factor 
      REAL scale 
      REAL glide 
      REAL eps 
      REAL sum_a, sum_t, lambda 
!                                                                       
      DATA eps / 0.0001 / 
!                                                                       
      abc (0) = 'm' 
      abc (1) = 'a' 
      abc (2) = 'b' 
      abc (3) = 'c' 
      xyz (1, - 2) = '-2x' 
      xyz (1, - 1) = ' -x' 
      xyz (1, 0)  = '  0' 
      xyz (1, 1) = ' +x' 
      xyz (1, 2) = '+2x' 
      xyz (2, - 2) = '-2y' 
      xyz (2, - 1) = ' -y' 
      xyz (2, 0)  = '  0' 
      xyz (2, 1) = ' +y' 
      xyz (2, 2) = '+2y' 
      xyz (3, - 2) = '-2z' 
      xyz (3, - 1) = ' -z' 
      xyz (3, 0)  = '  0' 
      xyz (3, 1) = ' +z' 
      xyz (3, 2) = '+2z' 
      xxx (1, - 2) = '-2x' 
      xxx (1, - 1) = ' -x' 
      xxx (1, 0)  = '  0' 
      xxx (1, 1) = ' +x' 
      xxx (1, 2) = '+2x' 
      xxx (2, - 2) = '-2x' 
      xxx (2, - 1) = ' -x' 
      xxx (2, 0)  = '  0' 
      xxx (2, 1) = ' +x' 
      xxx (2, 2) = '+2x' 
      xxx (3, - 2) = '-2x' 
      xxx (3, - 1) = ' -x' 
      xxx (3, 0)  = '  0' 
      xxx (3, 1) = ' +x' 
      xxx (3, 2) = '+2x' 
      Oyy (1,  - 2)  = '  0' 
      Oyy (1,  - 1)  = '  0' 
      Oyy (1, 0)  = '  0' 
      Oyy (1, 1)  = '  0' 
      Oyy (1, 2)  = '  0' 
      Oyy (2, - 2) = '-2y' 
      Oyy (2, - 1) = ' -y' 
      Oyy (2, 0)  = '  0' 
      Oyy (2, 1) = ' +y' 
      Oyy (2, 2) = '+2y' 
      Oyy (3, - 2) = '-2y' 
      Oyy (3, - 1) = ' -y' 
      Oyy (3, 0)  = '   ' 
      Oyy (3, 1) = ' +y' 
      Oyy (3, 2) = '+2y' 
      ctrans ( - 24) = '-1/1, ' 
      ctrans ( - 23)  = '    , ' 
      ctrans ( - 22)  = '    , ' 
      ctrans ( - 21)  = '    , ' 
      ctrans ( - 20) = '-5/6, ' 
      ctrans ( - 19)  = '    , ' 
      ctrans ( - 18) = '-3/4, ' 
      ctrans ( - 17)  = '    , ' 
      ctrans ( - 16) = '-2/3, ' 
      ctrans ( - 15)  = '    , ' 
      ctrans ( - 14)  = '    , ' 
      ctrans ( - 13)  = '    , ' 
      ctrans ( - 12) = '-1/2, ' 
      ctrans ( - 11)  = '    , ' 
      ctrans ( - 10)  = '    , ' 
      ctrans ( - 9)  = '    , ' 
      ctrans ( - 8) = '-1/3, ' 
      ctrans ( - 7)  = '    , ' 
      ctrans ( - 6) = '-1/4, ' 
      ctrans ( - 5)  = '    , ' 
      ctrans ( - 4) = '-1/6, ' 
      ctrans ( - 3) = '-1/8, ' 
      ctrans ( - 2)  = '    , ' 
      ctrans ( - 1)  = '    , ' 
      ctrans (0)  = '   0, ' 
      ctrans (1)  = '    , ' 
      ctrans (2)  = '    , ' 
      ctrans (3) = '+1/8, ' 
      ctrans (4) = '+1/6, ' 
      ctrans (5)  = '    , ' 
      ctrans (6) = '+1/4, ' 
      ctrans (7)  = '    , ' 
      ctrans (8) = '+1/3, ' 
      ctrans (9)  = '    , ' 
      ctrans (10)  = '    , ' 
      ctrans (11)  = '    , ' 
      ctrans (12) = '+1/2, ' 
      ctrans (13)  = '    , ' 
      ctrans (14)  = '    , ' 
      ctrans (15)  = '    , ' 
      ctrans (16) = '+2/3, ' 
      ctrans (17)  = '    , ' 
      ctrans (18) = '+3/4, ' 
      ctrans (19)  = '    , ' 
      ctrans (20) = '+5/6, ' 
      ctrans (21)  = '    , ' 
      ctrans (22)  = '    , ' 
      ctrans (23)  = '    , ' 
      ctrans (24) = '+1/1, ' 
!                                                                       
      DO i = 1, 3 
      axis (i) = 0.0 
      screw (i) = 0.0 
      posit (i) = 0.0 
      hkl (i) = 0.0 
      ENDDO 
!                                                                       
      IF (w_char (1:2) .eq.' 1') then 
!                                                                       
!     handle (1,w), vector w is a pure centering translation            
         DO i = 1, 3 
         posit (i) = 0.0 
         screw (i) = add (i) 
         ENDDO 
      ELSEIF (w_char (1:2) .eq.'-1') then 
!                                                                       
!     handle (-1,w), vector w/2 is the position of the -1               
         DO i = 1, 3 
         posit1bar (i) = add (i) / 2.0 
         screw (i) = 0.0 
         ENDDO 
      ELSE 
!                                                                       
!     --To detemine the axis, we need to calculate                      
!       (Matrix - Identity)**-1                                         
!                                                                       
         DO i = 1, 3 
         DO j = 1, 3 
         temp (i, j) = work (i, j) 
         sum (i, j) = 0.0 
         ENDDO 
         sum (i, i) = 1.0 
         ENDDO 
!                                                                       
         factor = 1. 
         IF (w_char (1:1) .eq.'-'.or.w_char (2:2) .eq.'m') then 
            factor = - 1. 
         ENDIF 
!                                                                       
         DO i = 1, 3 
         DO j = 1, 3 
         cp (i, j) = work (i, j) 
         ENDDO 
         cp (i, i) = cp (i, i) - factor 
         ENDDO 
!                                                                       
!     --Usually the determinant of (Matrix - Identity) is singular      
!       a 2x2 submatrix can however be solved                           
!                                                                       
         lsearch = .true. 
         l = 4 
         DO while (lsearch.and.l.gt.0) 
         l = l + 1 
         ii = (mod (l - 1 + 1, 3) ) + 1 
         jj = (mod (l - 1 + 2, 3) ) + 1 
         kk = (mod (l - 1 + 3, 3) ) + 1 
         IF (cp (ii, ii) * cp (jj, jj) - cp (ii, jj) * cp (jj, ii)      &
         .ne.0.0) then                                                  
            det = cp (ii, ii) * cp (jj, jj) - cp (ii, jj) * cp (jj, ii) 
            axis (ii) = - (cp (jj, jj) * cp (ii, kk) - cp (ii, jj)      &
            * cp (jj, kk) ) / det                                       
            axis (jj) = - (cp (ii, ii) * cp (jj, kk) - cp (jj, ii)      &
            * cp (ii, kk) ) / det                                       
            axis (kk) = 1.0 
            IF (abs (abs (axis (ii) ) - 0.5) .lt.eps.or.abs (abs (axis (&
            jj) ) - 0.5) .lt.eps) then                                  
               axis (ii) = 2 * axis (ii) 
               axis (jj) = 2 * axis (jj) 
               axis (kk) = 2 * axis (kk) 
            ENDIF 
            lsearch = .false. 
         ENDIF 
         ENDDO 
!                                                                       
!     --Find how many elements are equal to zero, this helps classifying
!                                                                       
         nnull = 0 
         DO i = 1, 3 
         IF (abs (axis (i) ) .lt.eps) then 
            nnull = nnull + 1 
         ENDIF 
         ENDDO 
!                                                                       
!     --Get proper sense of direction                                   
!                                                                       
         IF (nnull.eq.0) then 
            IF (axis (1) * axis (2) * axis (3) .lt.0) then 
               DO i = 1, 3 
               axis (i) = - axis (i) 
               ENDDO 
            ENDIF 
         ELSEIF (nnull.eq.1) then 
            DO i = 1, 3 
            IF (abs (axis (i) ) .lt.eps) then 
               jj = (mod (i - 1 + 1, 3) ) + 1 
               kk = (mod (i - 1 + 2, 3) ) + 1 
               IF (axis (jj) .lt.0.0) then 
                  axis (i) = - axis (i) 
                  axis (jj) = - axis (jj) 
                  axis (kk) = - axis (kk) 
               ENDIF 
            ENDIF 
            ENDDO 
         ENDIF 
      ENDIF 
!                                                                       
!     Determine screw part/glide part, exclude -N axes since t=0 always 
!                                                                       
      IF (w_char (2:2) .ne.'1') then 
         IF (w_char (1:1) .ne.'-') then 
!                                                                       
!     --First add M**(l-1) for l = 1 to power-1                         
            DO l = 2, power 
            DO i = 1, 3 
            DO j = 1, 3 
            sum (i, j) = sum (i, j) + temp (i, j) 
            ENDDO 
            ENDDO 
            IF (l.lt.power) then 
               CALL matmulx (t2, temp, work) 
               DO i = 1, 3 
               DO j = 1, 3 
               temp (i, j) = t2 (i, j) 
               ENDDO 
               ENDDO 
            ENDIF 
            ENDDO 
            DO i = 1, 3 
            screw (i) = 0.0 
            DO j = 1, 3 
            screw (i) = screw (i) + sum (i, j) * add (j) 
            ENDDO 
            screw (i) = screw (i) / power 
            ENDDO 
         ENDIF 
!                                                                       
!     --Subtract glide part from translation part                       
!                                                                       
         DO i = 1, 3 
         vector (i) = add (i) - screw (i) 
         ENDDO 
!                                                                       
!     --Determine Position                                              
!                                                                       
         DO i = 1, 3 
         DO j = 1, 3 
         cp (i, j) = - work (i, j) 
         ENDDO 
         cp (i, i) = 1.0 - work (i, i) 
         ENDDO 
         det = cp (1, 1) * (cp (2, 2) * cp (3, 3) - cp (3, 2) * cp (2,  &
         3) ) - cp (1, 2) * (cp (2, 1) * cp (3, 3) - cp (3, 1) * cp (2, &
         3) ) + cp (1, 3) * (cp (2, 1) * cp (3, 2) - cp (3, 1) * cp (2, &
         2) )                                                           
!                                                                       
         IF (det.ne.0) then 
!     --This allows to calculate the position directly                  
            CALL invmat (imat, cp) 
            DO i = 1, 3 
            posit1bar (i) = 0.0 
            DO j = 1, 3 
            posit1bar (i) = posit1bar (i) + imat (i, j) * vector (j) 
            ENDDO 
            posit (i) = posit1bar (i) 
            ENDDO 
         ELSE 
            DO i = 1, 3 
            DO j = 1, 3 
            temp (i, j) = work (i, j) * factor 
            sum (i, j) = 0.0 
            ENDDO 
            sum (i, i) = power - 1.0 
            ENDDO 
!                                                                       
!     --First add (power-l)*M**(l-1) for l = 1 to power-1               
            DO l = 2, power - 1 
            DO i = 1, 3 
            DO j = 1, 3 
            sum (i, j) = sum (i, j) + (power - l) * temp (i, j) 
            ENDDO 
            ENDDO 
            IF (l.lt.power - 1) then 
               CALL matmulx (t2, temp, work) 
               DO i = 1, 3 
               DO j = 1, 3 
               temp (i, j) = t2 (i, j) 
               ENDDO 
               ENDDO 
            ENDIF 
            ENDDO 
            DO i = 1, 3 
            posit (i) = 0.0 
            DO j = 1, 3 
            posit (i) = posit (i) + sum (i, j) * vector (j) 
            ENDDO 
            posit (i) = posit (i) / power 
            ENDDO 
         ENDIF 
         nnull = 0 
         DO i = 1, 3 
         IF (abs (axis (i) ) .lt.eps) then 
            nnull = nnull + 1 
         ENDIF 
         ENDDO 
      ENDIF 
!                                                                       
!     Cleanup to look like Tables                                       
!                                                                       
      DO i = 1, 3 
      screw (i) = screw (i) - int (screw (i) ) 
      posit (i) = posit (i) - int (posit (i) ) 
      posit1bar (i) = posit1bar (i) - int (posit1bar (i) ) 
      ENDDO 
      IF (w_char (1:2) .eq.' 1') then 
!     --write centerin vector                                           
         IF (abs (screw (1) ) .gt.eps.or.abs (screw (2) )               &
         .gt.eps.or.abs (screw (3) ) .gt.eps) then                      
            WRITE (w_char (4:21), 3100) (ctrans (nint (24 * screw (i) ) &
            ), i = 1, 3)                                                
         ENDIF 
      ELSEIF (w_char (1:2) .eq.'-1') then 
!     --write position of center of symmetry                            
         WRITE (w_char (46:63), 3200) (ctrans (nint (24 * posit1bar (i) &
         ) ), i = 1, 3)                                                 
      ELSEIF (w_char (1:1) .eq.'-') then 
!     --write position of center of symmetry                            
         WRITE (w_char (46:63), 3200) (ctrans (nint (24 * posit1bar (i) &
         ) ), i = 1, 3)                                                 
         w_char (45:45) = ';' 
      ELSEIF (w_char (2:2) .eq.'m') then 
!     --for mirror plane, convert axis to hkl                           
         CALL trans (axis, cr_gten, hkl, 3) 
         hmin = 1.e10 
         DO i = 1, 3 
         IF (abs (hkl (i) ) .gt.eps) then 
            hmin = min (hmin, abs (hkl (i) ) ) 
         ENDIF 
         ENDDO 
         DO i = 1, 3 
         hkl (i) = float (nint (hkl (i) / hmin) ) 
         ENDDO 
      ENDIF 
      IF (w_char (2:2) .ne.'1'.and.w_char (2:2) .ne.'m') then 
!                                                                       
!     --For rotation axis, clean up positions to look like tables       
!                                                                       
         IF (nnull.eq.0) then 
            posit (1) = posit (1) - posit (3) * sign (1., axis (1) )    &
            * sign (1., axis (3) )                                      
            posit (2) = posit (2) - posit (3) * sign (1., axis (2) )    &
            * sign (1., axis (3) )                                      
            posit (3) = 0.0 
         ELSE 
            DO l = 1, 3 
            i = l 
            j = mod (l, 3) + 1 
            k = mod (l + 1, 3) + 1 
            IF (axis (i) .ne.0.and.axis (j) .eq.0.and.axis (k) .eq.0)   &
            then                                                        
               posit (i) = 0.0 
            ELSEIF (axis (i) .eq.axis (j) .and.axis (k) .eq.0) then 
               posit (j) = posit (j) - posit (i) 
               posit (i) = 0.0 
            ELSEIF (axis (i) .gt.0.and.axis (i) .eq. - axis (j)         &
            .and.axis (k) .eq.0) then                                   
               posit (j) = posit (j) + posit (i) 
               posit (i) = 0.0 
            ELSEIF (axis (i) .lt.0.and.axis (i) .eq. - axis (j)         &
            .and.axis (k) .eq.0) then                                   
               posit (i) = posit (i) + posit (j) 
               posit (j) = 0.0 
            ENDIF 
            ENDDO 
         ENDIF 
         DO i = 1, 3 
         screw (i) = screw (i) - int (screw (i) ) 
         posit (i) = posit (i) - int (posit (i) ) 
         ENDDO 
         IF (abs (screw (1) ) .gt.eps.or.abs (screw (2) )               &
         .gt.eps.or.abs (screw (3) ) .gt.eps) then                      
            WRITE (w_char (4:21), 3100) (ctrans (nint (24 * screw (i) ) &
            ), i = 1, 3)                                                
         ENDIF 
         DO l = 1, 3 
         IF (abs (posit (l) ) .gt.eps) then 
            ia = 21 + (l - 1) * 8 + 4 
            ie = 21 + (l - 1) * 8 + 8 
            w_char (ia:ie) = ctrans (nint (24. * posit (l) ) ) 
         ENDIF 
         ENDDO 
         IF (w_char (1:1) .ne.'-') then 
            w_char (45:45) = ' ' 
         ELSE 
            IF (power.gt.2) then 
               w_char (45:45) = ';' 
            ENDIF 
         ENDIF 
      ELSEIF (w_char (2:2) .eq.'m') then 
         DO l = 1, 3 
         i = l 
         j = mod (l, 3) + 1 
         k = mod (l + 1, 3) + 1 
         IF (hkl (i) .ne.0.and.hkl (j) .eq.0.and.hkl (k) .eq.0) then 
            posit (j) = 0.0 
            posit (k) = 0.0 
         ELSEIF (hkl (i) .eq.hkl (j) .and.hkl (k) .eq.0) then 
            posit (i) = posit (j) + posit (i) 
            posit (j) = 0.0 
            posit (k) = 0.0 
         ELSEIF (hkl (i) .gt.0.and.hkl (i) .eq. - hkl (j) .and.hkl (k)  &
         .eq.0) then                                                    
            posit (i) = posit (i) - posit (j) 
            posit (j) = 0.0 
            posit (k) = 0.0 
         ELSEIF (hkl (i) .lt.0.and.hkl (i) .eq. - hkl (j) .and.hkl (k)  &
         .eq.0) then                                                    
            posit (j) = posit (j) - posit (i) 
            posit (i) = 0.0 
            posit (k) = 0.0 
         ENDIF 
         ENDDO 
         DO i = 1, 3 
         screw (i) = screw (i) - int (screw (i) ) 
         posit (i) = posit (i) - int (posit (i) ) 
         ENDDO 
         IF (abs (hkl (1) ) .eq.2..and.abs (hkl (2) ) .eq.1.) then 
            posit (2) = posit (2) + posit (1) * sign (2., hkl (1) )     &
            * sign (1., hkl (2) )                                       
            posit (1) = 0.0 
         ELSEIF (abs (hkl (1) ) .eq.1..and.abs (hkl (2) ) .eq.2.) then 
            posit (1) = posit (1) + posit (2) * sign (2., hkl (2) )     &
            * sign (1., hkl (1) )                                       
            posit (2) = 0.0 
         ENDIF 
         DO l = 1, 3 
         IF (abs (posit (l) ) .gt.eps) then 
            ia = 21 + (l - 1) * 8 + 4 
            ie = 21 + (l - 1) * 8 + 8 
            w_char (ia:ie) = ctrans (nint (24. * posit (l) ) ) 
         ENDIF 
         ENDDO 
         w_char (45:45) = ' ' 
      ENDIF 
!                                                                       
!     Get x,y,z symbol for the axis                                     
!                                                                       
      IF (w_char (2:2) .ne.'1'.and.w_char (2:2) .ne. ('m') ) then 
         DO i = 1, 3 
         ia = 21 + (i - 1) * 8 + 1 
         ie = 21 + (i - 1) * 8 + 3 
         IF (nnull.eq.2) then 
            w_char (ia:ie) = xyz (i, nint (axis (i) ) ) 
         ELSEIF (nnull.eq.1) then 
            IF (axis (1) .ne.0.) then 
               w_char (ia:ie) = xxx (i, nint (axis (i) ) ) 
            ELSE 
               w_char (ia:ie) = Oyy (i, nint (axis (i) ) ) 
            ENDIF 
         ELSEIF (nnull.eq.0) then 
            w_char (ia:ie) = xxx (i, nint (axis (i) ) ) 
         ENDIF 
      IF ( (w_char (ie+1:ie+1) .eq.'+'.or.w_char (ie+1:ie+1) .eq.'-') .a&
     &nd.w_char (ia:ie) .eq.'  0') then                                 
      w_char (ia:ie)  = '  ' 
         ENDIF 
         ENDDO 
         w_char (29:29) = ',' 
         w_char (37:37) = ',' 
      ENDIF 
!                                                                       
!     Get x,y,z symbol for the mirror plane                             
!                                                                       
      IF (w_char (2:2) .eq.'m') then 
         nnull = 0 
         DO i = 1, 3 
         IF (hkl (i) .eq.0.) then 
            nnull = nnull + 1 
         ENDIF 
         ENDDO 
         IF (nnull.eq.2) then 
            DO i = 1, 3 
            ia = 21 + (i - 1) * 8 + 1 
            ie = 21 + (i - 1) * 8 + 3 
            w_char (ia:ie) = xyz (i, nint (1 - hkl (i) ) ) 
            ENDDO 
         ELSEIF (nnull.eq.1) then 
            lsearch = .true. 
            DO l = 3, 5 
            IF (lsearch) then 
               i = mod (l + 2, 3) + 1 
               j = mod (l, 3) + 1 
               k = mod (l + 1, 3) + 1 
               IF (abs (hkl (i) ) .lt.eps) then 
                  ia = 21 + (i - 1) * 8 + 1 
                  ie = 21 + (i - 1) * 8 + 3 
                  w_char (ia:ie) = xyz (i, 1) 
                  jj = 1 
                  scale = - jj * (hkl (j) / hkl (k) ) 
                  kk = - nint (jj * hkl (j) / hkl (k) ) 
                  IF (abs (scale) .lt.1) then 
                     jj = jj / abs (scale) 
                     kk = - nint (jj * hkl (j) / hkl (k) ) 
                  ENDIF 
                  IF (abs (hkl (1) ) .lt.eps) then 
                     ia = 21 + (j - 1) * 8 + 1 
                     ie = 21 + (j - 1) * 8 + 3 
                     w_char (ia:ie) = Oyy (j, jj) 
                     ia = 21 + (k - 1) * 8 + 1 
                     ie = 21 + (k - 1) * 8 + 3 
                     w_char (ia:ie) = Oyy (k, kk) 
                  ELSE 
                     ia = 21 + (j - 1) * 8 + 1 
                     ie = 21 + (j - 1) * 8 + 3 
                     w_char (ia:ie) = xxx (j, jj) 
                     ia = 21 + (k - 1) * 8 + 1 
                     ie = 21 + (k - 1) * 8 + 3 
                     w_char (ia:ie) = xxx (k, kk) 
                  ENDIF 
                  lsearch = .false. 
               ENDIF 
            ENDIF 
            ENDDO 
         ENDIF 
         DO i = 1, 3 
         ia = 21 + (i - 1) * 8 + 1 
         ie = 21 + (i - 1) * 8 + 3 
      IF ( (w_char (ie+1:ie+1) .eq.'+'.or.w_char (ie+1:ie+1) .eq.'-') .a&
     &nd.w_char (ia:ie) .eq.'  0') then                                 
      w_char (ia:ie)  = '  ' 
         ENDIF 
         ENDDO 
         w_char (29:29) = ',' 
         w_char (37:37) = ',' 
!                                                                       
!     --get glide plane symbol                                          
!                                                                       
         IF (abs (screw (1) ) .gt.eps.or.abs (screw (2) )               &
         .gt.eps.or.abs (screw (3) ) .gt.eps) then                      
            nnullg = 0 
            iglide = 0 
            DO i = 1, 3 
            IF (screw (i) .eq.0.) then 
               nnullg = nnullg + 1 
            ELSE 
               iglide = i 
            ENDIF 
            ENDDO 
            IF (nnullg.eq.2) then 
               w_char (2:2) = abc (iglide) 
            ELSEIF (nnull.eq.2.and.nnullg.eq.1) then 
               IF (abs (abs (screw (1) ) - abs (screw (2) ) )           &
               .lt.eps.or.abs (abs (screw (1) ) - abs (screw (3) ) )    &
               .lt.eps.or.abs (abs (screw (2) ) - abs (screw (3) ) )    &
               .lt.eps) then                                            
                  glide = max (abs (screw (1) ), abs (screw (2) ),      &
                  abs (screw (3) ) )                                    
                  IF (glide.eq.0.50) then 
                     w_char (2:2) = 'n' 
                  ELSEIF (glide.eq.0.25) then 
                     w_char (2:2) = 'd' 
                  ELSE 
                     w_char (2:2) = 'g' 
                  ENDIF 
               ELSE 
                  w_char (2:2) = 'g' 
               ENDIF 
            ELSEIF (nnull.eq.1.and.nnullg.eq.0) then 
               IF (abs (abs (screw (1) ) - abs (screw (2) ) )           &
               .lt.eps.and.abs (abs (screw (1) ) - abs (screw (3) ) )   &
               .lt.eps) then                                            
                  glide = max (abs (screw (1) ), abs (screw (2) ),      &
                  abs (screw (3) ) )                                    
                  IF (glide.eq.0.50) then 
                     w_char (2:2) = 'n' 
                  ELSEIF (glide.eq.0.25) then 
                     w_char (2:2) = 'd' 
                  ELSE 
                     w_char (2:2) = 'g' 
                  ENDIF 
               ELSE 
                  w_char (2:2) = 'g' 
               ENDIF 
            ELSE 
               w_char (2:2) = 'g' 
            ENDIF 
            WRITE (w_char (4:21), 3100) (ctrans (nint (24 * screw (i) ) &
            ), i = 1, 3)                                                
         ENDIF 
      ENDIF 
!                                                                       
!     get sense of rotation                                             
!                                                                       
      IF (power.gt.2) then 
         vector (1) = 3 
         vector (2) = 7 
         vector (3) = 19 
         DO i = 1, 3 
         p3 (i) = 0.0 
         DO j = 1, 3 
         p3 (i) = p3 (i) + work (i, j) * vector (j) 
         ENDDO 
         ENDDO 
         DO i = 1, 3 
         cp (i, 1) = axis (i) 
         cp (i, 2) = vector (i) 
         cp (i, 3) = p3 (i) 
         ENDDO 
         det = (cp (1, 1) * (cp (2, 2) * cp (3, 3) - cp (3, 2) * cp (2, &
         3) ) - cp (1, 2) * (cp (2, 1) * cp (3, 3) - cp (3, 1) * cp (2, &
         3) ) + cp (1, 3) * (cp (2, 1) * cp (3, 2) - cp (3, 1) * cp (2, &
         2) ) ) * factor                                                
         IF (det.gt.0) then 
            w_char (3:3) = 'P' 
         ELSEIF (det.lt.0) then 
            w_char (3:3) = 'M' 
         ENDIF 
      ENDIF 
!                                                                       
!     -rename a centering vector                                        
!                                                                       
      IF (w_char (1:2) .eq.' 1') then 
         IF (abs (screw (1) ) .gt.eps.or.abs (screw (2) )               &
         .gt.eps.or.abs (screw (2) ) .gt.eps) then                      
            w_char (1:2) = ' t' 
         ENDIF 
      ENDIF 
!                                                                       
!DBG      write(*,*) ' Axis               ',axis                        
!DBG      write(*,*) ' Screw/Glide        ',screw                       
!DBG      write(*,*) ' Position           ',posit                       
!DBG      if(w_char(1:1).eq.'-') then                                   
!DBG      write(*,*) ' Position -1        ',posit1bar                   
!DBG      endif                                                         
!DBG      write(*,*) ' hkl                ',hkl                         
!DBG      write(*,4000)                     no,mod(no-1,48)+1,w_char    
!DBG      write(*,*) '          ',                                      
!DBG     &                      '123456789 123456789 123456789 ',       
!DBG     &                      '123456789 123456789 123456789 '        
!                                                                       
 3100 FORMAT    (' (',a4,',',a4,',',a4,') ') 
 3200 FORMAT    ('  ',a4,',',a4,',',a4,'  ') 
 3300 FORMAT    (' ') 
 4000 FORMAT    (i4,2x,i3,2x,a) 
!                                                                       
      END SUBROUTINE get_detail                     
!********************************************************************** 
      SUBROUTINE wyckoff_main (zeile, lp) 
!-                                                                      
!     Determines the local symmetry of the position given in the line   
!+                                                                      
      IMPLICIT none 
!                                                                       
      include'errlist.inc' 
!                                                                       
      INTEGER FULL, SYMBOL, XYZ, MATRIX 
      PARAMETER (FULL = 0) 
      PARAMETER (SYMBOL = 1) 
      PARAMETER (XYZ = 2) 
      PARAMETER (MATRIX = 3) 
!                                                                       
      CHARACTER(1024) zeile 
      INTEGER lp 
!                                                                       
      INTEGER maxw 
      PARAMETER (maxw = 4) 
!                                                                       
      CHARACTER(1024) cpara (maxw) 
      INTEGER lpara (maxw) 
      INTEGER ianz 
      INTEGER iianz 
      INTEGER mode 
      LOGICAL loutput 
      REAL werte (maxw) 
!                                                                       
      LOGICAL str_comp 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
      IF (ier_num.ne.0) return 
      IF (ianz.eq.3.or.ianz.eq.4) then 
         iianz = 3 
         CALL ber_params (iianz, cpara, lpara, werte, maxw) 
         IF (ier_num.ne.0) return 
!                                                                       
         loutput = .true. 
         IF (ianz.eq.3.or.str_comp (cpara (4) , 'full', 2, lpara (4) ,  &
         4) ) then                                                      
            mode = FULL 
         ELSEIF (str_comp (cpara (4) , 'symbol', 2, lpara (4) , 6) )    &
         then                                                           
            mode = SYMBOL 
         ELSEIF (str_comp (cpara (4) , 'xyz', 2, lpara (4) , 3) ) then 
            mode = XYZ 
         ELSEIF (str_comp (cpara (4) , 'matrix', 2, lpara (4) , 6) )    &
         then                                                           
            mode = MATRIX 
         ELSE 
            ier_num = - 6 
            ier_typ = ER_COMM 
            RETURN 
         ENDIF 
         CALL get_wyckoff (werte, loutput, mode) 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
      END SUBROUTINE wyckoff_main                   
!********************************************************************** 
      SUBROUTINE get_wyckoff (vec, loutput, mode) 
!-                                                                      
!     Determines the local symmetry of position xyz within the unit cell
!+                                                                      
      USE config_mod 
      USE crystal_mod 
      USE wyckoff_mod 
      USE unitcell_mod 
      IMPLICIT none 
!                                                                       
       
      include'prompt.inc' 
!                                                                       
      REAL vec (3) 
      LOGICAL loutput 
      INTEGER mode 
!                                                                       
      INTEGER FULL, SYMBOL, XYZ, MATRIX 
      PARAMETER (FULL = 0) 
      PARAMETER (SYMBOL = 1) 
      PARAMETER (XYZ = 2) 
      PARAMETER (MATRIX = 3) 
!                                                                       
      INTEGER i, j 
      INTEGER is 
      INTEGER n_center 
      INTEGER igroup 
      INTEGER block 
      LOGICAL lident 
      REAL orig (4) 
      REAL copy (4) 
      REAL eps 
!                                                                       
      DATA eps / 0.00001 / 
!                                                                       
      n_center = 1 
      IF (cr_spcgr (1:1) .eq.'P') then 
         n_center = 1 
      ELSEIF (cr_spcgr (1:1) .eq.'A') then 
         n_center = 2 
      ELSEIF (cr_spcgr (1:1) .eq.'B') then 
         n_center = 2 
      ELSEIF (cr_spcgr (1:1) .eq.'C') then 
         n_center = 2 
      ELSEIF (cr_spcgr (1:1) .eq.'I') then 
         n_center = 2 
      ELSEIF (cr_spcgr (1:1) .eq.'F') then 
         n_center = 4 
      ELSEIF (cr_spcgr (1:1) .eq.'R'.and.cr_syst.eq.6) then 
         n_center = 3 
      ENDIF 
      IF (gen_sta.eq.GEN_SYMM) then 
         block = spc_n / n_center 
      ELSEIF (gen_sta.eq.GEN_CENTER) then 
         block = n_center 
      ENDIF 
      IF (loutput) then 
         WRITE (output_io, 900) vec 
      ENDIF 
!                                                                       
      wyc_n = 0 
!                                                                       
!     move position into first unit cell,ia                             
!                                                                       
      CALL firstcell (vec, 3) 
      DO i = 1, 3 
      orig (i) = vec (i) 
      ENDDO 
      orig (4) = 1.0 
!                                                                       
!     apply all symmetry operations to original position                
!                                                                       
      DO is = 1, spc_n 
      IF (gen_sta.eq.GEN_SYMM) then 
         igroup = mod (is - 1, block) + 1 
      ELSEIF (gen_sta.eq.GEN_CENTER) then 
         igroup = (is - 1) / block + 1 
      ENDIF 
      DO i = 1, 4 
      copy (i) = 0.0 
      DO j = 1, 4 
      copy (i) = copy (i) + spc_mat (i, j, is) * orig (j) 
      ENDDO 
      ENDDO 
      CALL firstcell (copy, 4) 
      lident = .true. 
      DO i = 1, 3 
      lident = lident.and.abs (orig (i) - copy (i) ) .lt.eps 
      ENDDO 
      IF (lident) then 
         wyc_n = wyc_n + 1 
         wyc_list (wyc_n) = is 
         IF (loutput) then 
            IF (mode.eq.FULL) then 
               WRITE (output_io, 1000) is, igroup 
               WRITE (output_io, 1100) (spc_mat (1, j, is), j = 1, 4),  &
               spc_char (is), (spc_mat (2, j, is), j = 1, 4), (spc_mat (&
               3, j, is), j = 1, 4), spc_xyz (is)                       
            ELSEIF (mode.eq.SYMBOL) then 
               WRITE (output_io, 3200) is, igroup, spc_char (is) 
            ELSEIF (mode.eq.XYZ) then 
               WRITE (output_io, 4200) is, igroup, spc_xyz (is) 
            ELSEIF (mode.eq.MATRIX) then 
               WRITE (output_io, 5200) is, igroup, (spc_mat (1, j, is), &
               j = 1, 4), (spc_mat (2, j, is), j = 1, 4), (spc_mat (3,  &
               j, is), j = 1, 4)                                        
            ENDIF 
         ENDIF 
      ENDIF 
!                                                                       
      ENDDO 
!                                                                       
      IF (loutput) then 
         WRITE (output_io, 6000) spc_n / wyc_n, wyc_n, spc_n 
      ENDIF 
!                                                                       
  900 FORMAT    (/,' Wyckoff symmetry for position ',3f12.6,/) 
 1000 FORMAT    ('Symmetry No.      [',i3,']  (',i3,')') 
 1100 FORMAT    (  ' ( ',3(f4.1,', '),f8.5,' )','  ',a65,/,             &
     &                    ' ( ',3(f4.1,', '),f8.5,' )',/,               &
     &                    ' ( ',3(f4.1,', '),f8.5,' )','  ',a87,/)      
 3200 FORMAT    ('Symmetry No.      [',i3,']  (',i3,')  ',a65) 
 4200 FORMAT    ('Symmetry No.      [',i3,']  (',i3,')  ',a87) 
 5200 FORMAT    ('Symmetry No.      [',i3,']  (',i3,')  ',              &
     &                    ' ( ',3(f4.1,', '),f8.5,' )',/,               &
     &                32x,' ( ',3(f4.1,', '),f8.5,' )',/,               &
     &                32x,' ( ',3(f4.1,', '),f8.5,' )'   )              
 6000 FORMAT    (/,' Multiplicity;   No of Sym. Op. in Wyckoff group; ',&
     &  '    Highest Multiplicity',/,i8,20x,i8,20x,i8)                  
!                                                                       
      END SUBROUTINE get_wyckoff                    
!********************************************************************** 
      SUBROUTINE symmetry 
!-                                                                      
!     Performs the space group symmetry on the current atom.            
!     cr_natoms             the current atom number                     
!     cr_iscat(cr_natoms)   its scattering type                         
!                           and thus number of chemically identical     
!                           yet symmetrically different atoms           
!                                                                       
!+                                                                      
      USE config_mod 
      USE crystal_mod 
      USE generate_mod 
      USE gen_add_mod 
      USE molecule_mod 
      USE sym_add_mod 
      USE unitcell_mod 
      IMPLICIT none 
!                                                                       
       
      include'errlist.inc' 
!                                                                       
      REAL eps 
      INTEGER i, j, ii, iii, igs, igg 
      INTEGER iiii 
      INTEGER mole_i, mole_st 
!                                                                       
!                                                                       
      DATA eps / 0.00001 / 
!                                                                       
!     ii is the number of the atom the symmetry is to work on           
!                                                                       
      ii = cr_natoms 
!                                                                       
!     Loop over all generators for spacegroup cr_spcgrno                
!                                                                       
      DO igs = 1, generspcgr (0, cr_spcgrno) 
      IF (gen_sta.eq.GEN_SYMM) then 
         igg = generspcgr (igs, cr_spcgrno) 
      ELSEIF (gen_sta.eq.GEN_CENTER) then 
         igg = generspcgr_center (igs, cr_spcgrno) 
      ENDIF 
!                                                                       
!     --iii is the number of the last atom generated by the previous    
!     --generator                                                       
!                                                                       
      iii = cr_natoms 
!                                                                       
      CALL symmetry_gener (NMAX, cr_natoms, cr_pos, cr_iscat, cr_prop,  &
      ii, iii, iii, igg, NG, generators, generpower)                    
!                                                                       
      IF (ier_num.ne.0) then 
         RETURN 
      ENDIF 
      ENDDO 
!                                                                       
!     Loop over all additional generators                               
!                                                                       
      DO igs = 1, gen_add_n 
!                                                                       
!     --iii is the number of the last atom generated by the previous    
!     --generator                                                       
!                                                                       
      iii = cr_natoms 
!                                                                       
      CALL symmetry_gener (NMAX, cr_natoms, cr_pos, cr_iscat, cr_prop,  &
      ii, iii, iii, igs, GEN_ADD_MAX, gen_add, gen_add_power)           
!                                                                       
      IF (ier_num.ne.0) then 
         RETURN 
      ENDIF 
      ENDDO 
!                                                                       
!     Loop over all additional symmetry operations that                 
!       are not generators                                              
!                                                                       
!                                                                       
!     --iiii is the number of the last atom generated by the last       
!     --generator                                                       
!                                                                       
      iiii = cr_natoms 
!                                                                       
      DO igs = 1, sym_add_n 
!                                                                       
!     --iii is the number of the last atom generated by the previous    
!     --symmetry operation                                              
!                                                                       
      iii = cr_natoms 
!                                                                       
      CALL symmetry_gener (NMAX, cr_natoms, cr_pos, cr_iscat, cr_prop,  &
      ii, iii, iiii, igs, SYM_ADD_MAX, sym_add, sym_add_power)          
!                                                                       
      IF (ier_num.ne.0) then 
         RETURN 
      ENDIF 
      ENDDO 
!                                                                       
!     Sort Atoms into molecules                                         
!                                                                       
      IF (mole_l_on) then 
         CALL mole_insert (ii) 
         IF (ier_num.ne.0) then 
            RETURN 
         ENDIF 
!                                                                       
      ENDIF 
!                                                                       
      END SUBROUTINE symmetry                       
!********************************************************************** 
      SUBROUTINE symmetry_gener (NMAX, cr_natoms, cr_pos, cr_iscat,     &
      cr_prop, ii, iii, iiii, igg, NG, generators, generpower)          
!-                                                                      
!     Applies the generator to the current atom                         
!+                                                                      
      USE molecule_mod 
      IMPLICIT none 
!                                                                       
      INTEGER NMAX 
!     include'crystal3.inc' 
      include'errlist.inc' 
!
      INTEGER,                       INTENT(INOUT)  :: cr_natoms
      INTEGER, DIMENSION(1:NMAX),    INTENT(INOUT)  :: cr_iscat
      INTEGER, DIMENSION(1:NMAX),    INTENT(INOUT)  :: cr_prop
      REAL   , DIMENSION(1:3,1:NMAX),INTENT(INOUT)  :: cr_pos
!                                                                       
      INTEGER NG 
      INTEGER ii, iii, iiii 
      INTEGER igg 
      INTEGER generpower (NG) 
!                                                                       
      REAL generators (4, 4, 0:NG) 
!                                                                       
      INTEGER ia, iaa, ipg 
      INTEGER i, j, k 
      LOGICAL lnew 
      REAL x (4), y (4), z (4), w (4) 
      REAL wmat (4, 4) 
      REAL xmat (4, 4) 
      REAL eps 
      REAL compare (4), previous (4) 
!                                                                       
      DATA eps / 0.00001 / 
!                                                                       
!     For convenience create Identity operator                          
!                                                                       
      DO i = 1, 4 
      DO j = 1, 4 
      xmat (i, j) = 0.0 
      ENDDO 
      xmat (i, i) = 1.0 
      ENDDO 
!                                                                       
!     Loop over all powers of generator igg                             
!                                                                       
      DO ipg = 1, generpower (igg) 
!                                                                       
!     --raise power of generator, xmat is a dummy matrix, equal to      
!     --the previous power of the generator                             
!                                                                       
      DO i = 1, 4 
      DO j = 1, 4 
      wmat (i, j) = 0.0 
      DO k = 1, 4 
      wmat (i, j) = wmat (i, j) + generators (i, k, igg) * xmat (k, j) 
      ENDDO 
      ENDDO 
      ENDDO 
!                                                                       
!     --Multiply all points with current generator                      
!                                                                       
      DO iaa = ii, iiii 
      DO i = 1, 3 
      x (i) = cr_pos (i, iaa) 
      ENDDO 
      x (4) = 1.0 
      CALL trans (x, wmat, y, 4) 
      DO i = 1, 3 
      compare (i) = y (i) 
      ENDDO 
      compare (4) = 1.0 
      CALL firstcell (compare, 4) 
      IF (.not.mole_l_on) then 
!                                                                       
!     ------Transform atom into first unit cell,                        
!           if it is not inside a molecule                              
!                                                                       
         CALL firstcell (y, 4) 
      ENDIF 
      lnew = .true. 
      DO ia = ii, iii 
      DO i = 1, 3 
      previous (i) = cr_pos (i, ia) 
      ENDDO 
      previous (4) = 1.0 
      CALL firstcell (previous, 4) 
!     ------check if atom exists                                        
      IF (abs (compare (1) - previous (1) ) .lt.eps.and.abs (compare (2)&
      - previous (2) ) .lt.eps.and.abs (compare (3) - previous (3) )    &
      .lt.eps) then                                                     
         lnew = .false. 
         GOTO 30 
      ENDIF 
      ENDDO 
   30 CONTINUE 
!                                                                       
!     ------insert atom into crystal                                    
!                                                                       
      IF (lnew) then 
         IF (cr_natoms.lt.nmax) then 
            cr_natoms = cr_natoms + 1 
            cr_pos (1, cr_natoms) = y (1) 
            cr_pos (2, cr_natoms) = y (2) 
            cr_pos (3, cr_natoms) = y (3) 
            cr_iscat (cr_natoms) = cr_iscat (ii) 
            cr_prop (cr_natoms) = cr_prop (ii) 
         ELSE 
            ier_num = -10 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
      ENDIF 
      ENDDO 
!                                                                       
!     set xmat to current power of generator                            
!                                                                       
      DO i = 1, 4 
      DO j = 1, 4 
      xmat (i, j) = wmat (i, j) 
      ENDDO 
      ENDDO 
!     --End of loop over all powers                                     
      ENDDO 
   10 CONTINUE 
!                                                                       
      END SUBROUTINE symmetry_gener                 
!********************************************************************** 
      SUBROUTINE mole_insert (ii) 
!-                                                                      
!     Sorts the newly created atoms into the correct molecules.         
!+                                                                      
      USE config_mod 
      USE crystal_mod 
      USE molecule_mod 
      USE wyckoff_mod 
      IMPLICIT none 
!                                                                       
       
      include'errlist.inc' 
!                                                                       
      INTEGER ii 
!                                                                       
      INTEGER i, j, l 
      INTEGER m, n 
      INTEGER is 
      INTEGER ifirst 
      INTEGER mole_st 
      INTEGER mole_natoms 
      INTEGER mole_temp (192) 
!                                                                       
      LOGICAL lsame 
!                                                                       
      REAL vec (3), orig (4), old (3) 
      REAL first (3) 
      REAL eps 
!                                                                       
      DATA eps / 0.00001 / 
!                                                                       
                                                                        
      DO i = 1, 192 
      mole_temp (i) = 0 
      ENDDO 
!                                                                       
!     The first atom of a molecule and its symmetrically equivalent     
!     are each sorted into new molecules                                
!                                                                       
      IF (mole_l_first) then 
         mole_num_act = mole_num_curr - 1 
         DO i = ii, cr_natoms 
         mole_num_act = mole_num_act + 1 
         CALL mole_insert_current (i, mole_num_act) 
         IF (ier_num.ne.0) then 
            RETURN 
         ENDIF 
         ENDDO 
         mole_l_first = .false. 
      ELSE 
!                                                                       
!     --These are atoms further down the line, here we have to check    
!     --into which molecule they belong. The first is inserted into     
!     --the current molecule                                            
!                                                                       
         mole_num_act = mole_num_curr 
         mole_st = mole_len (mole_num_curr) + 1 
!                                                                       
!     --Loop over all secondary atoms                                   
!                                                                       
         DO i = ii, cr_natoms 
!                                                                       
!     ----If not yet in a molecule, add to active molecule              
!                                                                       
         IF (mole_temp (i - ii + 1) .eq.0) then 
            CALL mole_insert_current (i, mole_num_act) 
            IF (ier_num.ne.0) then 
               RETURN 
            ENDIF 
            mole_temp (i - ii + 1) = mole_num_act 
!                                                                       
!     ----Determine coordinates of first atom in active molecule        
!                                                                       
            ifirst = mole_cont (mole_off (mole_num_act) + 1) 
            first (1) = cr_pos (1, ifirst) 
            first (2) = cr_pos (2, ifirst) 
            first (3) = cr_pos (3, ifirst) 
!                                                                       
            CALL get_wyckoff (first, .false., 0) 
!                                                                       
!     ----Create copies of this atom with all Wyckoff symmetry operators
!     ----1.st Wyckoff Symmetry is always identity, thus we             
!       can ignore this                                                 
!                                                                       
            orig (1) = cr_pos (1, i) 
            orig (2) = cr_pos (2, i) 
            orig (3) = cr_pos (3, i) 
            orig (4) = 1.0 
            DO is = 2, wyc_n 
            DO l = 1, 3 
            vec (l) = 0.0 
            DO m = 1, 4 
            vec (l) = vec (l) + spc_mat (l, m, wyc_list (is) ) * orig ( &
            m)                                                          
            ENDDO 
            ENDDO 
            CALL firstcell (vec, 3) 
!                                                                       
!     ------Loop to compare the copies to the remaining original atom   
!                                                                       
            DO j = i + 1, cr_natoms 
            DO l = 1, 3 
            old (l) = cr_pos (l, j) 
            ENDDO 
            CALL firstcell (old, 3) 
            lsame = .true. 
            DO l = 1, 3 
            lsame = lsame.and.abs (old (l) - vec (l) ) .lt.eps 
            ENDDO 
            IF (lsame) then 
!     ----------This is another atom of the same molecule               
               IF (mole_temp (j - ii + 1) .eq.0) then 
                  CALL mole_insert_current (j, mole_num_act) 
                  mole_temp (j - ii + 1) = mole_num_act 
               ENDIF 
            ENDIF 
            ENDDO 
            ENDDO 
            mole_num_act = mole_num_act + 1 
         ENDIF 
         ENDDO 
!     --End of loop over secondary atoms                                
!                                                                       
!     --Makesure that the molecules are not split into different        
!       unit cells                                                      
         CALL first_mole (mole_st) 
      ENDIF 
!                                                                       
      END SUBROUTINE mole_insert                    
!********************************************************************** 
      SUBROUTINE mole_insert_current (iatom, imole) 
!-                                                                      
!     Inserts the last atom into the molecule list as last atom of      
!     the specified molecule.                                           
!+                                                                      
      USE molecule_mod 
      IMPLICIT none 
!                                                                       
      include'errlist.inc' 
!                                                                       
      INTEGER iatom 
      INTEGER imole 
!                                                                       
      INTEGER i 
!                                                                       
!     Move the content of all molecules after "imole" one down the list 
!     mole_num_atom : total number of atoms in molecules                
!     mole_num_acur : number of atoms in molecules including "imole"    
!                                                                       
      mole_num_atom = mole_off (mole_num_mole) + mole_len (             &
      mole_num_mole)                                                    
!                                                                       
!     If necessary, create new molecule                                 
!                                                                       
      IF (imole.gt.mole_num_mole) then 
         IF (imole.le.MOLE_MAX_MOLE) then 
            mole_num_mole = mole_num_mole+1 
            mole_len (imole) = 0 
            mole_off (imole) = mole_num_atom 
            mole_type (imole) = mole_num_type 
            mole_char (imole) = mole_char (imole-1) 
            mole_dens (imole) = mole_dens (imole-1) 
         ELSE 
            ier_num = - 65 
            ier_typ = ER_APPL 
            RETURN 
         ENDIF 
      ENDIF 
!                                                                       
      mole_num_acur = mole_off (imole) + mole_len (imole) 
!                                                                       
      IF (mole_num_atom + 1.le.MOLE_MAX_ATOM) then 
         DO i = mole_num_atom, mole_num_acur + 1, - 1 
         mole_cont (i + 1) = mole_cont (i) 
         ENDDO 
      ELSE 
         ier_num = - 74 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
!                                                                       
!     move the offset of all molecules after "imole" one down           
!                                                                       
      DO i = imole+1, mole_num_mole 
      mole_off (i) = mole_off (i) + 1 
      ENDDO 
!                                                                       
!     insert atom "iatom" at the end of molecule "imole"                
!                                                                       
      IF (mole_num_acur + 1.le.MOLE_MAX_ATOM) then 
         mole_cont (mole_num_acur + 1) = iatom 
         mole_len (imole) = mole_len (imole) + 1 
         mole_num_acur = mole_num_acur + 1 
      ELSE 
         ier_num = - 74 
         ier_typ = ER_APPL 
         RETURN 
      ENDIF 
!                                                                       
      END SUBROUTINE mole_insert_current            
!********************************************************************** 
      SUBROUTINE first_mole (mole_st) 
!-                                                                      
!     Moves atoms by +- one unit cell to keep molecules concatenated    
!     The algorithm of this subroutine is based on the assumption,      
!     that the first atom is on the highest point of symmetry of the    
!     molecule, i.e. is not copied by any of the symmetry operations    
!     or generators that make up the molecule symmetry. This has the    
!     desired effect that the bond distance of all symmetry             
!     equivalent atoms to this first atom is constant and can be taken  
!     as a reference. If a symmetry operation of the space group        
!     moves an atom out of the current unit cell, one can move this     
!     atom back by integer unit cell vectors until the bond distance    
!     is correct again.                                                 
!     If the molecule does not include an atom on the highest point     
!     of the molecule symmetry, you must insert a "void" on this site.  
!+                                                                      
      USE config_mod 
      USE crystal_mod 
      USE molecule_mod 
      IMPLICIT none 
!                                                                       
       
      include'errlist.inc' 
      include'prompt.inc' 
!                                                                       
      LOGICAL lspace 
      PARAMETER (lspace = .true.) 
!                                                                       
      INTEGER mole_st 
      INTEGER i, j, k1, k2, k3 
      INTEGER k1u, k1o, k2u, k2o, k3u, k3o 
      REAL d, dd 
      REAL u (3), v (3) 
      REAL x, y, z 
!                                                                       
      REAL do_blen 
!                                                                       
!.......calculate metric and reciprocal metric tensor,reciprocal lattice
!       constants and permutation tensors                               
      CALL setup_lattice (cr_a0, cr_ar, cr_eps, cr_gten, cr_reps,       &
      cr_rten, cr_win, cr_wrez, cr_v, cr_vr, .false., cr_gmat, cr_fmat, &
      cr_cartesian)                                                     
!                                                                       
!     Only one atom in the current molecule, return immediately         
!                                                                       
!DBG                                                                    
!DBG      write (output_io,*) 'mole_len(mole_num_curr)',                
!DBG     &                                mole_len(mole_num_curr)       
!DBG      write (output_io,*) 'mole_st                ',mole_st         
                                                                        
      IF (mole_len (mole_num_curr) .eq.1) return 
      IF (mole_st.eq.1) mole_st = 2 
      IF (mole_len (mole_num_curr) .lt.mole_st) return 
!                                                                       
      u (1) = cr_pos (1, mole_cont (mole_off (mole_num_curr) + 1) ) 
      u (2) = cr_pos (2, mole_cont (mole_off (mole_num_curr) + 1) ) 
      u (3) = cr_pos (3, mole_cont (mole_off (mole_num_curr) + 1) ) 
!                                                                       
      v (1) = cr_pos (1, mole_cont (mole_off (mole_num_curr) + mole_st) &
      )                                                                 
      v (2) = cr_pos (2, mole_cont (mole_off (mole_num_curr) + mole_st) &
      )                                                                 
      v (3) = cr_pos (3, mole_cont (mole_off (mole_num_curr) + mole_st) &
      )                                                                 
!                                                                       
      d = do_blen (lspace, u, v) 
!DBG                                                                    
!DBG      write (output_io,5555) u,v,d                                  
!DBG5555      format(3f10.4,2x,3f10.4,2x,f12.4/)                        
!                                                                       
!     Loop over all molecules from current to last                      
!                                                                       
      DO i = mole_num_curr, mole_num_mole 
!DBG                                                                    
!DBG      write (output_io,*) ' Molecule : ',i                          
      u (1) = cr_pos (1, mole_cont (mole_off (i) + 1) ) 
      u (2) = cr_pos (2, mole_cont (mole_off (i) + 1) ) 
      u (3) = cr_pos (3, mole_cont (mole_off (i) + 1) ) 
!                                                                       
!     Loop over all atoms on places mole_st and higher                  
!                                                                       
      DO j = mole_st, mole_len (i) 
!                                                                       
!     ----Determine proper limits for the loop over neighboring         
!          unit cells                                                   
!                                                                       
      x = cr_pos (1, mole_cont (mole_off (i) + j) ) 
      IF ( - 1.0.lt.x.and.x.lt.0.or.x.ge.1.0) then 
         k1u = - int (x) 
         k1o = k1u + 1 
      ELSE 
         k1u = - int (x) - 1 
         k1o = k1u + 1 
      ENDIF 
      y = cr_pos (2, mole_cont (mole_off (i) + j) ) 
      IF ( - 1.0.lt.y.and.y.lt.0.or.y.ge.1.0) then 
         k2u = - int (y) 
         k2o = k2u + 1 
      ELSE 
         k2u = - int (y) - 1 
         k2o = k2u + 1 
      ENDIF 
      z = cr_pos (3, mole_cont (mole_off (i) + j) ) 
      IF ( - 1.0.lt.z.and.z.lt.0.or.z.ge.1.0) then 
         k3u = - int (z) 
         k3o = k3u + 1 
      ELSE 
         k3u = - int (z) - 1 
         k3o = k3u + 1 
      ENDIF 
!DBG                                                                    
!DBG      write (output_io,*) ' x,y,z    ',x,y,z                        
!DBG      write (output_io,*) 'loop over ',k1u,k1o,k2u,k2o,k3u,k3o      
!                                                                       
!     ----perform loop over next unit cells                             
!                                                                       
      DO k1 = 2, - 2, - 1 
      DO k2 = 2, - 2, - 1 
      DO k3 = 2, - 2, - 1 
      v (1) = cr_pos (1, mole_cont (mole_off (i) + j) ) + float (k1) 
      v (2) = cr_pos (2, mole_cont (mole_off (i) + j) ) + float (k2) 
      v (3) = cr_pos (3, mole_cont (mole_off (i) + j) ) + float (k3) 
      dd = do_blen (lspace, u, v) 
!DBG                                                                    
!DBG      write (output_io,5556) u,j,mole_off(i)+j,v,dd                 
!DBG5556      format(3f10.4,2(2x,i2),2x,3f10.4,2x,f12.4)                
      IF (abs (d-dd) .lt.0.01) then 
         GOTO 10 
      ENDIF 
      ENDDO 
      ENDDO 
      ENDDO 
!                                                                       
      ier_num = - 83 
      ier_typ = ER_APPL 
      RETURN 
!                                                                       
   10 CONTINUE 
      cr_pos (1, mole_cont (mole_off (i) + j) ) = v (1) 
      cr_pos (2, mole_cont (mole_off (i) + j) ) = v (2) 
      cr_pos (3, mole_cont (mole_off (i) + j) ) = v (3) 
!DBG                                                                    
!DBG      write (output_io,5557) k1,k2,k3                               
!DBG5557      format('loesung fuer ',3i3)                               
!                                                                       
      ENDDO 
      ENDDO 
!                                                                       
      END SUBROUTINE first_mole                     
!********************************************************************** 
      SUBROUTINE mole_firstcell 
!-                                                                      
!     Moves molecules whose first atom is outside the unit cell into    
!     the first unit cell.                                              
!+                                                                      
      USE config_mod 
      USE molecule_mod 
      USE crystal_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      REAL eps 
      PARAMETER (eps = 1.e-5) 
!                                                                       
      INTEGER i, j, k 
      REAL x, x1 
!                                                                       
!     Loop over all molecules from current to last                      
!                                                                       
      DO j = mole_num_curr, mole_num_mole 
      DO i = 1, 3 
      x = cr_pos (i, mole_cont (mole_off (j) + 1) ) 
      x1 = float (int (x) ) 
      IF (x - x1.lt. - eps) x1 = x1 - 1 
      IF (abs (x1) .gt.eps) then 
!                                                                       
!     ------Loop over all atoms in molecule j                           
!                                                                       
         DO k = 1, mole_len (j) 
         cr_pos (i, mole_cont (mole_off (j) + k) ) = cr_pos (i,         &
         mole_cont (mole_off (j) + k) ) - x1                            
         ENDDO 
      ENDIF 
      ENDDO 
      ENDDO 
      END SUBROUTINE mole_firstcell                 
!********************************************************************** 
      SUBROUTINE firstcell (y, idim) 
!-                                                                      
!     truncates atomic position to fractal position of                  
!     0.0 <= x < 1                                                      
!+                                                                      
      IMPLICIT none 
!                                                                       
      INTEGER idim, i 
      REAL y (idim) 
!                                                                       
      DO i = 1, 3 
      y (i) = y (i) - float (int (y (i) ) ) 
      IF (y (i) .lt.0.0) y (i) = y (i) + 1 
      IF (y (i) .eq.1.0) y (i) = 0.0 
      ENDDO 
      END SUBROUTINE firstcell                      
!********************************************************************** 
      SUBROUTINE rese_cr 
!                                                                       
!     resets the crystal structure to empty                             
!                                                                       
      USE config_mod 
      USE crystal_mod 
      USE gen_add_mod 
      USE molecule_mod 
      USE sym_add_mod 
      USE save_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER i 
!                                                                       
      cr_natoms = 0 
      as_natoms = 0 
      cr_ncatoms = 1 
      cr_nscat = 0 
      cr_icc       = 1
      cr_cartesian = .false. 
      cr_scat_int (0) = .true. 
      cr_delf_int (0) = .true. 
      cr_scat_equ (0) = .false. 
      cr_sel_prop (0) = 0 
      cr_sel_prop (1) = 0 
      DO i = 1, maxscat 
      cr_at_lis (i) = ' ' 
      cr_at_equ (i) = ' ' 
      as_at_lis (i) = ' ' 
      cr_scat_int (i) = .true. 
      cr_delf_int (i) = .true. 
      cr_scat_equ (i) = .false. 
      ENDDO 
      DO i = 1, 3 
      cr_dim (i, 1) = 0.0 
      cr_dim (i, 2) = 0.0 
      ENDDO 
!                                                                       
      DO i = 0, MOLE_MAX_MOLE 
      mole_len (i) = 0 
      mole_off (i) = 0 
      ENDDO 
      mole_l_on = .false. 
      mole_num_mole = 0 
      mole_num_curr = 0 
      mole_num_act = 0 
      mole_num_type = 0 
      mole_gene_n = 0 
      mole_symm_n = 0 
!                                                                       
      sym_add_n = 0 
      gen_add_n = 0 
!                                                                       
      sav_r_ncell = .false. 
!                                                                       
      END SUBROUTINE rese_cr                        
!********************************************************************** 
      SUBROUTINE recip_symm 
!-                                                                      
!     Creates the symmetry matrices in reciprocal space                 
!+                                                                      
      USE config_mod 
      USE crystal_mod 
      USE recipro_mod 
      IMPLICIT none 
!                                                                       
      INTEGER i, j, k 
      LOGICAL lfriedel_remove 
      LOGICAL lacentric 
      REAL a (3, 3) 
      REAL b (3, 3) 
      REAL h (4, rec_max_sym) 
!                                                                       
!     Set up a general reflection                                       
!                                                                       
      h (1, 1) = 0.1 
      h (2, 1) = 0.2 
      h (3, 1) = 0.3 
      h (4, 1) = 0.0 
      DO i = 1, 4 
      DO j = 1, 4 
      rec_sym (i, j, 1) = 0.0 
      ENDDO 
      rec_sym (i, i, 1) = 1.0 
      ENDDO 
!                                                                       
!     Create all symmetry operations in real space that create          
!     a different vector, ignoring the translational part. The          
!     translational part of each first new matrix is saved, since       
!     this is needed for proper phase assignment.                       
!     lfriedel_remove  Logical variable that signals to remove          
!                      symmetry operations -1                           
!     cr_acentric      logical variable that describes whether the      
!                      space group is acentric or not                   
!                                                                       
      lfriedel_remove = .false. 
      CALL rmc_symmetry (rec_n_sym, h, rec_sym, rec_max_sym,            &
      lfriedel_remove, lacentric)                                       
      cr_acentric = lacentric 
!                                                                       
!     Transform all symmetry operations into reciprocal space           
!                                                                       
      DO k = 1, rec_n_sym 
      DO i = 1, 3 
      DO j = 1, 3 
      a (i, j) = rec_sym (i, j, k) 
      ENDDO 
      ENDDO 
!                                                                       
!     --do transformation q = gSg*                                      
!                                                                       
      CALL matmulx (b, a, cr_rten) 
      CALL matmulx (a, cr_gten, b) 
      DO i = 1, 3 
      DO j = 1, 3 
      rec_sym (i, j, k) = a (i, j) 
      ENDDO 
      ENDDO 
!                                                                       
!DBG      write( *,1000) k                                              
!DBG          do i=1,4                                                  
!DBG            write (*,1010) (rec_sym(i,j,k),j=1,4)                   
!DBG          ENDDO                                                     
!DBG1000      format('**********************'/' Matrix Number ', i4)    
!DBG1010      format(3(f4.1,2x),f7.4)                                   
      ENDDO 
!                                                                       
      END SUBROUTINE recip_symm                     
!*****7**************************************************************** 
      SUBROUTINE update_cr_dim 
!-                                                                      
!     Updates the crystal dimensions to the current values              
!+                                                                      
      USE config_mod 
      USE crystal_mod 
      IMPLICIT none 
!                                                                       
       
!                                                                       
      INTEGER i, j 
!                                                                       
!     Set initial values                                                
!                                                                       
      IF (cr_natoms.gt.0) then 
         DO j = 1, 3 
         cr_dim (j, 1) = cr_pos (j, 1) 
         cr_dim (j, 2) = cr_pos (j, 1) 
         ENDDO 
      ELSE 
         DO j = 1, 3 
         cr_dim (j, 1) = 1.e10 
         cr_dim (j, 1) = - 1.e10 
         ENDDO 
      ENDIF 
!                                                                       
!     Update values from all atoms in crystal                           
!                                                                       
      DO i = 1, cr_natoms 
      DO j = 1, 3 
      cr_dim (j, 1) = amin1 (cr_dim (j, 1), cr_pos (j, i) ) 
      cr_dim (j, 2) = amax1 (cr_dim (j, 2), cr_pos (j, i) ) 
      ENDDO 
      ENDDO 
!                                                                       
      END SUBROUTINE update_cr_dim                  
!*****7**************************************************************** 
      SUBROUTINE setup_lattice (cr_a0, cr_ar, cr_eps, cr_gten, cr_reps, &
      cr_rten, cr_win, cr_wrez, cr_v, cr_vr, lout, cr_gmat, cr_fmat,    &
      cr_cartesian)                                                     
!-                                                                      
!     Updates the crystal lattice and symmetry information              
!+                                                                      
      IMPLICIT none 
!                                                                       
      INTEGER i 
      LOGICAL lout 
      LOGICAL cr_cartesian 
      REAL cr_a0 (3), cr_ar (3), cr_eps (3, 3, 3), cr_gten (3, 3) 
      REAL cr_reps (3, 3, 3), cr_rten (3, 3), cr_win (3), cr_wrez (3) 
      REAL cr_v, cr_vr, cr_gmat (3, 3), cr_fmat (3, 3) 
      REAL hkl (3) 
      REAL u (3) 
      REAL xc (3) 
      REAL yc (3) 
      REAL zc (3) 
      REAL dist 
!                                                                       
      CALL lattice (cr_a0, cr_ar, cr_eps, cr_gten, cr_reps, cr_rten,    &
      cr_win, cr_wrez, cr_v, cr_vr, lout)                               
      cr_cartesian = cr_a0 (1) .eq.1..and.cr_a0 (2) .eq.1..and.cr_a0 (3)&
      .eq.1..and.cr_win (1) .eq.90..and.cr_win (2) .eq.90..and.cr_win ( &
      3) .eq.90.                                                        
      DO i = 1, 3 
      hkl (i) = 0.0 
      u (i) = 0.0 
      ENDDO 
      hkl (3) = 1.0 
      CALL trafo (hkl, u, xc, yc, zc, cr_gmat, cr_fmat, dist, cr_eps,   &
      cr_gten, cr_reps, cr_rten)                                        
      CALL recip_symm 
!                                                                       
      END SUBROUTINE setup_lattice                  
!*****7**************************************************************** 
      SUBROUTINE do_import (zeile, lp) 
!-                                                                      
!     imports a file into discus.cell format                            
!+                                                                      
      IMPLICIT none 
!                                                                       
      include'errlist.inc' 
!                                                                       
      INTEGER maxw 
      PARAMETER (MAXW = 5) 
!                                                                       
      CHARACTER ( * ) zeile 
      INTEGER lp 
!                                                                       
      CHARACTER(1024) cpara (MAXW) 
      INTEGER lpara (MAXW) 
      INTEGER ianz 
!                                                                       
      LOGICAL str_comp 
!                                                                       
      CALL get_params (zeile, ianz, cpara, lpara, MAXW, lp) 
      IF (ier_num.ne.0) then 
         RETURN 
      ENDIF 
!                                                                       
      IF (ianz.ge.1) then 
         IF (str_comp (cpara (1) , 'shelx', 2, lpara (1) , 5) ) then 
            IF (ianz.eq.2) then 
               CALL del_params (1, ianz, cpara, lpara, maxw) 
               IF (ier_num.ne.0) return 
               CALL ins2discus (ianz, cpara, lpara, MAXW) 
            ELSE 
               ier_num = - 6 
               ier_typ = ER_COMM 
            ENDIF 
         ELSE 
            ier_num = - 86 
            ier_typ = ER_APPL 
         ENDIF 
      ELSE 
         ier_num = - 6 
         ier_typ = ER_COMM 
      ENDIF 
!                                                                       
      END SUBROUTINE do_import                      
!*****7**************************************************************** 
      SUBROUTINE ins2discus (ianz, cpara, lpara, MAXW) 
!-                                                                      
!     converts a SHELXL "ins" or "res" file to DISCUS                   
!+                                                                      
      IMPLICIT none 
!                                                                       
      include'errlist.inc' 
!                                                                       
      INTEGER ianz 
      INTEGER MAXW 
      CHARACTER ( * ) cpara (MAXW) 
      INTEGER lpara (MAXW) 
!                                                                       
      INTEGER NFV 
      PARAMETER (NFV = 50) 
!                                                                       
      REAL werte (3) 
!                                                                       
      INTEGER shelx_num 
      PARAMETER (shelx_num = 62) 
      CHARACTER(4) shelx_ign (1:shelx_num) 
      CHARACTER(2) c_atom (20) 
      CHARACTER(4) command 
      CHARACTER(80) line1 
      CHARACTER(80) line2 
      CHARACTER(160) line 
      CHARACTER(1024) infile 
      CHARACTER(1024) ofile 
      INTEGER ird, iwr 
      INTEGER i, j, ii, jj 
      INTEGER ix, iy, iz, idot 
      INTEGER ntyp 
      INTEGER length, length1, length2, lp 
      INTEGER icont 
      INTEGER centering 
      INTEGER ityp 
      INTEGER ifv 
      LOGICAL lread 
      LOGICAL lwrite 
      LOGICAL lmole, lmole_wr 
      LOGICAL lcontinue 
      REAL z, latt (6) 
      REAL xyz (3) 
      REAL sof 
      REAL uiso, uij (6) 
      REAL gen (3, 4) 
      REAL fv (NFV) 
!                                                                       
      INTEGER len_str 
!                                                                       
      DATA shelx_ign / 'ACTA', 'AFIX', 'ANIS', 'BASF', 'BIND', 'BLOC',  &
      'BOND', 'BUMP', 'CGLS', 'CHIV', 'CONF', 'CONN', 'DAMP', 'DANG',   &
      'DEFS', 'DELU', 'DFIX', 'DISP', 'EADP', 'EQIV', 'EXTI', 'EXYZ',   &
      'FEND', 'FLAT', 'FMAP', 'FRAG', 'FREE', 'GRID', 'HFIX', 'HOPE',   &
      'HTAB', 'ISOR', 'L.S.', 'LAUE', 'LIST', 'MERG', 'MORE', 'MOVE',   &
      'MPLA', 'NCSY', 'OMIT', 'PART', 'PLAN', 'REM ', 'RESI', 'RTAB',   &
      'SADI', 'SAME', 'SHEL', 'SIMU', 'SIZE', 'SPEC', 'SUMP', 'STIR',   &
      'SWAT', 'TEMP', 'TIME', 'TWIN', 'UNIT', 'WGHT', 'WPDB', 'ZERR' /  
!                                                                       
      DO i = 1, NFV 
      fv (i) = 0.0 
      ENDDO 
!                                                                       
      lmole = .false. 
      lmole_wr = .true. 
!                                                                       
      CALL do_build_name (ianz, cpara, lpara, werte, maxw, 1) 
      IF (ier_num.ne.0) then 
         RETURN 
      ENDIF 
      infile = cpara (1) 
      i = index (infile, '.') 
      IF (i.eq.0) then 
         infile = cpara (1) (1:lpara (1) ) //'.ins' 
         ofile = cpara (1) (1:lpara (1) ) //'.cell' 
      ELSE 
         ofile = cpara (1) (1:i) //'cell' 
      ENDIF 
      lread = .true. 
      lwrite = .false. 
      ird = 34 
      iwr = 35 
      CALL oeffne (ird, infile, 'old', lread) 
      IF (ier_num.ne.0) then 
         RETURN 
      ENDIF 
      CALL oeffne (iwr, ofile, 'unknown', lwrite) 
      IF (ier_num.ne.0) then 
         RETURN 
      ENDIF 
!                                                                       
      lcontinue = .false. 
      READ (ird, 1000, end = 900, err = 900) line1 
      length1 = len_str (line1) 
      IF (length1.gt.0) then 
         icont = index (line1, '=') 
         IF (icont.gt.0) then 
            READ (ird, 1000, end = 900, err = 900) line2 
            length2 = len_str (line2) 
            line = line1 (1:icont - 1) //' '//line2 (1:length2) 
         ELSE 
            line = line1 
         ENDIF 
      ELSE 
         line = line1 
      ENDIF 
      length = len_str (line) 
      IF (length.gt.0) then 
         command = line (1:4) 
      ELSE 
      command = '    ' 
      ENDIF 
      DO i = 1, shelx_num 
      lcontinue = lcontinue.or.command.eq.shelx_ign (i) 
      ENDDO 
!                                                                       
      DO while (command.ne.'FVAR'.and.command.ne.'MOLE') 
      IF (lcontinue) then 
         CONTINUE 
      ELSEIF (command.eq.'TITL') then 
         WRITE (iwr, 2000) line (6:length) 
         WRITE (iwr, 2100) 
      ELSEIF (command.eq.'CELL') then 
         READ (line (6:length), * ) z, latt 
         WRITE (iwr, 2200) latt 
      ELSEIF (command.eq.'LATT') then 
         READ (line (6:length), * ) centering 
         IF (abs (centering) .eq.1) then 
            CONTINUE 
         ELSEIF (abs (centering) .eq.2) then 
            WRITE (iwr, 2320) 
         ELSEIF (abs (centering) .eq.3) then 
            WRITE (iwr, 2330) 
         ELSEIF (abs (centering) .eq.4) then 
            WRITE (iwr, 2340) 
            WRITE (iwr, 2341) 
         ELSEIF (abs (centering) .eq.5) then 
            WRITE (iwr, 2350) 
         ELSEIF (abs (centering) .eq.6) then 
            WRITE (iwr, 2360) 
         ELSEIF (abs (centering) .eq.7) then 
            WRITE (iwr, 2370) 
         ENDIF 
         IF (centering.gt.0) then 
            WRITE (iwr, 2400) 
         ENDIF 
      ELSEIF (command.eq.'SFAC') then 
         ntyp = 0 
         j = 5 
         DO while (j.lt.length) 
         j = j + 1 
         DO while (j.lt.length.and.line (j:j) .eq.' ') 
         j = j + 1 
         ENDDO 
         IF (j.le.length) then 
            ntyp = ntyp + 1 
            c_atom (ntyp) = ' ' 
            i = 0 
            DO while (j.le.length.and.line (j:j) .ne.' ') 
            i = i + 1 
            c_atom (ntyp) (i:i) = line (j:j) 
            j = j + 1 
            ENDDO 
         ENDIF 
         ENDDO 
         WRITE (iwr, 2500) (c_atom (i) , ',', i = 1, ntyp - 1) , c_atom &
         (ntyp)                                                         
      ELSEIF (command.eq.'SYMM') then 
         lp = length - 5 
         CALL get_params (line (6:length), ianz, cpara, lpara, maxw, lp) 
         IF (ianz.eq.3) then 
            DO i = 1, 3 
            DO jj = 1, 4 
            gen (i, jj) = 0.0 
            ENDDO 
            ix = index (cpara (i) , 'X') 
            IF (ix.gt.0) then 
               gen (i, 1) = 1.0 
               IF (ix.gt.1.and.cpara (i) (ix - 1:ix - 1) .eq.'-') then 
                  gen (i, 1) = - 1.0 
                  cpara (i) (ix - 1:ix - 1) = ' ' 
               ELSEIF (ix.gt.1.and.cpara (i) (ix - 1:ix - 1) .eq.'+')   &
               then                                                     
                  gen (i, 1) = 1.0 
                  cpara (i) (ix - 1:ix - 1) = ' ' 
               ENDIF 
               cpara (i) (ix:ix) = ' ' 
            ENDIF 
            iy = index (cpara (i) , 'Y') 
            IF (iy.gt.0) then 
               gen (i, 2) = 1.0 
               IF (iy.gt.1.and.cpara (i) (iy - 1:iy - 1) .eq.'-') then 
                  gen (i, 2) = - 1.0 
                  cpara (i) (iy - 1:iy - 1) = ' ' 
               ELSEIF (iy.gt.1.and.cpara (i) (iy - 1:iy - 1) .eq.'+')   &
               then                                                     
                  gen (i, 2) = 1.0 
                  cpara (i) (iy - 1:iy - 1) = ' ' 
               ENDIF 
               cpara (i) (iy:iy) = ' ' 
            ENDIF 
            iz = index (cpara (i) , 'Z') 
            IF (iz.gt.0) then 
               gen (i, 3) = 1.0 
               IF (iz.gt.1.and.cpara (i) (iz - 1:iz - 1) .eq.'-') then 
                  gen (i, 3) = - 1.0 
                  cpara (i) (iz - 1:iz - 1) = ' ' 
               ELSEIF (iz.gt.1.and.cpara (i) (iz - 1:iz - 1) .eq.'+')   &
               then                                                     
                  gen (i, 3) = 1.0 
                  cpara (i) (iz - 1:iz - 1) = ' ' 
               ENDIF 
               cpara (i) (iz:iz) = ' ' 
            ENDIF 
            ENDDO 
            DO i = 1, 3 
            idot = index (cpara (i) , '.') 
            IF (idot.eq.0) then 
               cpara (i) (lpara (i) + 1:lpara (i) + 1) = '.' 
               cpara (i) (lpara (i) + 2:lpara (i) + 2) = '0' 
               lpara (i) = lpara (i) + 2 
            ENDIF 
            CALL rem_bl (cpara (i), lpara (i) ) 
            ENDDO 
            CALL ber_params (ianz, cpara, lpara, werte, maxw) 
            gen (1, 4) = werte (1) 
            gen (2, 4) = werte (2) 
            gen (3, 4) = werte (3) 
            WRITE (iwr, 2600) ( (gen (i, j), j = 1, 4), i = 1, 3) 
         ENDIF 
      ENDIF 
!                                                                       
      lcontinue = .false. 
      READ (ird, 1000, end = 900, err = 900) line1 
      length1 = len_str (line1) 
      IF (length1.gt.0) then 
         icont = index (line1, '=') 
         IF (icont.gt.0) then 
            READ (ird, 1000, end = 900, err = 900) line2 
            length2 = len_str (line2) 
            line = line1 (1:icont - 1) //' '//line2 (1:length2) 
         ELSE 
            line = line1 
         ENDIF 
      ELSE 
         line = line1 
      ENDIF 
      length = len_str (line) 
      IF (length.gt.0) then 
         command = line (1:4) 
      ELSE 
      command = '    ' 
      ENDIF 
      DO i = 1, shelx_num 
      lcontinue = lcontinue.or.command.eq.shelx_ign (i) 
      ENDDO 
      ENDDO 
!                                                                       
      WRITE (iwr, 3000) 
!                                                                       
      DO while (command.ne.'HKLF') 
      IF (lcontinue) then 
         CONTINUE 
      ELSEIF (command.eq.'FVAR') then 
         READ (line (6:length), *, end = 800) (fv (i), i = 1, NFV) 
  800    CONTINUE 
      ELSEIF (command.eq.'MOLE') then 
         IF (lmole) then 
            WRITE (iwr, 4000) 'molecule end' 
            lmole_wr = .true. 
         ELSE 
            lmole = .true. 
            lmole_wr = .true. 
         ENDIF 
         CONTINUE 
      ELSE 
         IF (lmole.and.lmole_wr) then 
            WRITE (iwr, 4000) 'molecule' 
            lmole_wr = .false. 
         ENDIF 
         READ (line (6:length), *, end = 850) ityp, xyz, sof, (uij (i), &
         i = 1, 6)                                                      
  850    CONTINUE 
         IF (i.lt.6) then 
            uiso = uij (1) 
         ELSE 
            uiso = (uij (1) + uij (2) + uij (3) ) / 3. 
         ENDIF 
         DO i = 1, 3 
         ifv = nint (xyz (i) / 10.) 
         IF (ifv.eq.1) then 
            xyz (i) = xyz (i) - 10. 
         ELSEIF (ifv.gt.1) then 
            xyz (i) = (xyz (i) - ifv * 10) * fv (ifv) 
         ELSEIF (ifv.lt. - 1) then 
            xyz (i) = (abs (xyz (i) ) + ifv * 10) * (1. - fv (iabs (ifv)&
            ) )                                                         
         ENDIF 
         ENDDO 
!         write(iwr,3100) c_atom(ityp),xyz,float(ityp)                  
         WRITE (iwr, 3100) c_atom (ityp), xyz, uiso 
      ENDIF 
!                                                                       
      lcontinue = .false. 
      READ (ird, 1000, end = 900, err = 900) line1 
      length1 = len_str (line1) 
      IF (length1.gt.0) then 
         icont = index (line1, '=') 
         IF (icont.gt.0) then 
            READ (ird, 1000, end = 900, err = 900) line2 
            length2 = len_str (line2) 
            line = line1 (1:icont - 1) //' '//line2 (1:length2) 
         ELSE 
            line = line1 
         ENDIF 
      ELSE 
         line = line1 
      ENDIF 
      length = len_str (line) 
      IF (length.gt.0) then 
         command = line (1:4) 
      ELSE 
      command = '    ' 
      ENDIF 
      DO i = 1, shelx_num 
      lcontinue = lcontinue.or.command.eq.shelx_ign (i) 
      ENDDO 
      ENDDO 
!                                                                       
  900 CONTINUE 
!                                                                       
      IF (lmole) then 
         WRITE (iwr, 4000) 'molecule end' 
      ENDIF 
!                                                                       
      CLOSE (ird) 
      CLOSE (iwr) 
!                                                                       
 1000 FORMAT    (a) 
 2000 FORMAT    ('title ',a) 
 2100 FORMAT    ('spcgr P1') 
 2200 FORMAT    ('cell ',5(2x,f9.4,','),2x,f9.4) 
 2320 FORMAT    ('gener  1.0, 0.0, 0.0, 0.5,',                          &
     &                     '    0.0, 1.0, 0.0, 0.5,',                   &
     &                     '    0.0, 0.0, 1.0, 0.5   1')                
 2330 FORMAT    ('gener  1.0, 0.0, 0.0, 0.66666667,',                   &
     &                     '    0.0, 1.0, 0.0, 0.33333333,',            &
     &                     '    0.0, 0.0, 1.0, 0.33333333,   2')        
 2340 FORMAT    ('gener  1.0, 0.0, 0.0, 0.0,',                          &
     &                     '    0.0, 1.0, 0.0, 0.5,',                   &
     &                     '    0.0, 0.0, 1.0, 0.5,   1')               
 2341 FORMAT    ('gener  1.0, 0.0, 0.0, 0.5,',                          &
     &                     '    0.0, 1.0, 0.0, 0.0,',                   &
     &                     '    0.0, 0.0, 1.0, 0.5,   1')               
 2350 FORMAT    ('gener  1.0, 0.0, 0.0, 0.0,',                          &
     &                     '    0.0, 1.0, 0.0, 0.5,',                   &
     &                     '    0.0, 0.0, 1.0, 0.5,   1')               
 2360 FORMAT    ('gener  1.0, 0.0, 0.0, 0.5,',                          &
     &                     '    0.0, 1.0, 0.0, 0.0,',                   &
     &                     '    0.0, 0.0, 1.0, 0.5,   1')               
 2370 FORMAT    ('gener  1.0, 0.0, 0.0, 0.5,',                          &
     &                     '    0.0, 1.0, 0.0, 0.5,',                   &
     &                     '    0.0, 0.0, 1.0, 0.0,  1')                
 2400 FORMAT    ('gener -1.0, 0.0, 0.0, 0.0,',                          &
     &                     '    0.0,-1.0, 0.0, 0.0,',                   &
     &                     '    0.0, 0.0,-1.0, 0.0,  1')                
 2500 FORMAT    ('scat ',20(a4,a1)) 
 2600 FORMAT    ('symm ',3(2X,4(1x,f12.8,',')),' 1.') 
!                                                                       
 3000 FORMAT    ('atoms') 
 3100 FORMAT    (a2,2x,4(2x,f9.5)) 
!                                                                       
 4000 FORMAT    (a) 
!                                                                       
      END SUBROUTINE ins2discus                     
!
      SUBROUTINE test_file ( strucfile, natoms, ntypes, init, lcell)
!
!     Determines the number of atoms and atom types in strucfile
!
      IMPLICIT NONE
!
      include'errlist.inc'

!
      CHARACTER (LEN=*), INTENT(IN)    :: strucfile
      INTEGER          , INTENT(INOUT) :: natoms
      INTEGER          , INTENT(INOUT) :: ntypes
      INTEGER          , INTENT(IN)    :: init
      LOGICAL          , INTENT(IN)    :: lcell
!
      INTEGER, PARAMETER                    :: MAXW = 13 
      CHARACTER(LEN=1024), DIMENSION(MAXW)  :: cpara (MAXW) 
      INTEGER            , DIMENSION(MAXW)  :: lpara (MAXW) 
      REAL               , DIMENSION(MAXW)  :: werte (MAXW) 
!
      REAL, PARAMETER                       :: eps = 1e-4
      CHARACTER (LEN=1024)                  :: line
      CHARACTER (LEN=1024)                  :: zeile
      CHARACTER (LEN=   4), DIMENSION(1024), SAVE :: names
      REAL                , DIMENSION(1024), SAVE :: bvals
      INTEGER                               :: ios
      INTEGER                               :: i
      INTEGER                               :: ianz   ! no of arguments
      INTEGER                               :: laenge ! length of input line
      INTEGER                               :: lp     ! length of parameter string
      INTEGER                               :: nscattypes ! no of SCAT arguments 
      INTEGER                               :: nadptypes  ! no of ADP  arguments 
      LOGICAL                               :: new
      REAL                                  :: x,y,z,bval
!
      INTEGER, EXTERNAL :: len_str
      LOGICAL           :: IS_IOSTAT_END
!
      natoms     = 0
      nscattypes = 0
      nadptypes  = 0
      IF ( init == -1 ) then
        names  = ' '
        bvals  = 0.0
        ntypes = 0
      ENDIF
!
      CALL oeffne ( 99, strucfile, 'old', .true. )
      IF ( ier_num /= 0) THEN
          CLOSE ( 99 )
          RETURN
      ENDIF
header: DO
        READ (99,1000, IOSTAT=ios) line
        IF ( ios /= 0 ) THEN
           ier_num = -6
           ier_typ = ER_IO
           CLOSE ( 99 )
           RETURN
        ENDIF
        CALL do_cap (line)
        laenge = len_str(line)
        IF ( laenge .gt. 4 ) then
           zeile = line(5:laenge)
           lp    = laenge - 4
        ELSE
           zeile = ' '
           lp    = 1
        ENDIF
        IF (line(1:4) == 'SCAT' ) then 
            CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
            IF (ier_num.eq.0) then 
               DO i = 1,ianz
                   names(ntypes+i) = cpara(i)(1:lpara(i))
               ENDDO
               nscattypes = ianz
            ELSE
               ier_num = -111
               ier_typ = ER_APPL
               CLOSE(99)
               RETURN
            ENDIF
        ELSEIF (line(1:3) == 'ADP' ) then 
            CALL get_params (zeile, ianz, cpara, lpara, maxw, lp) 
            IF (ier_num.eq.0) then 
               CALL ber_params (ianz, cpara, lpara, werte, maxw) 
               IF (ier_num.eq.0) then 
                  DO i = 1,ianz
                      bvals(ntypes+i) = werte(i)
                  ENDDO
                  nadptypes = ianz
               ELSE
                  ier_num = -111
                  ier_typ = ER_APPL
                  CLOSE(99)
                  RETURN
               ENDIF
            ELSE
               ier_num = -111
               ier_typ = ER_APPL
               CLOSE(99)
               RETURN
            ENDIF
        ENDIF
        IF (line(1:4) == 'ATOM') EXIT header
      ENDDO header
!
      IF (nscattypes /= nadptypes) THEN
         ier_num = -113
         ier_typ = ER_APPL
         CLOSE(99)
         RETURN
      ENDIF
!
      ntypes = MAX(ntypes,nscattypes)
!
main: DO
        READ (99,1000, IOSTAT=ios) line
        IF ( IS_IOSTAT_END(ios) ) EXIT main
        CALL do_cap (line)
        READ (line(5:len_str(line)), *, IOSTAT = ios) x,y,z,bval
isatom: IF ( ios == 0 ) THEN
           natoms = natoms + 1
           new = .true.
types:     DO i=1,ntypes
              IF ( LINE(1:4) == names(i) ) THEN
                 IF ( lcell ) THEN
                    new = .false.
                    EXIT types
                 ELSEIF ( abs(abs(bval)-abs(bvals(i))) < eps ) THEN
                    new = .false.
                    EXIT types
                 ENDIF
              ENDIF
           ENDDO types
           IF ( new ) THEN
              ntypes = ntypes + 1
              names(ntypes) = line(1:4)
              bvals(ntypes) = bval
           ENDIF
        ENDIF isatom
      ENDDO main
!
      CLOSE (99)
!
1000  FORMAT(a)
!
      END SUBROUTINE test_file
END MODULE structur
